"""
Synteny/context analysis around GWAS-significant Panaroo clusters.

Overview:
- Loads Panaroo pangenome outputs (CSV, Rtab, GML) and genome GFFs.
- Resolves GWAS-significant terms to Panaroo cluster IDs.
- For each sample containing a given cluster, locates the corresponding CDS in the GFF
  and extracts a ±window neighborhood to describe local gene context (IDs/products).
- Maps neighbor CDS IDs back to Panaroo clusters to compare observed neighbor clusters
  with Panaroo graph neighbors (from the GML).
- Writes per-cluster:
  * <cluster>_synteny.tsv — row per sample with target and neighbor info and presence flags.
  * <cluster>_signatures.tsv — frequency of product-context signatures across samples.
  * <cluster>_cluster_signatures.tsv — frequency of cluster-context signatures across samples.

Inputs (paths configured below):
- Panaroo CSV: gene_presence_absence_filt_pseudo_length_frag.csv
- Panaroo Rtab: gene_presence_absence_filt_pseudo_length_frag.Rtab
- Panaroo GML: pre_filt_graph.gml
- Panaroo pangenome reference FASTA: pan_genome_reference.fa (loaded but unused here)
- GWAS table: mortality_COGs.txt (expects a column with cluster or variant identifiers)
- GFF folder: GFFs for samples; file basenames must match Panaroo sample column names.

Outputs:
- synteny_context_<GWAS_basename> directory containing TSVs per significant cluster.

Notes:
- The code tries to robustly match GWAS terms to clusters via direct cluster IDs,
  non-unique gene names, and tokenized aliases (handles separators like ~~~, ;, ,).
- GFF parsing uses BCBio.GFF; features considered: CDS, gene, ORF, protein.
- Presence/absence flags combine Rtab presence with GFF localization results.
"""

#%% Load packages and data (Change accordingly)

import pandas as pd
import os
import pickle
import networkx as nx
from Bio import SeqIO
from BCBio import GFF
from collections import Counter, defaultdict

# Set working directory to the GWAS data location
os.chdir('/home/albertotr/OneDrive/Data/Cambridge_Project/GWAS/')

# Configure data folders produced by Panaroo and genome GFFs
panaroo_data = '/home/albertotr/OneDrive/Data/Cambridge_Project/Pangenome_results_HF'
gff_folder = '/home/albertotr/OneDrive/Data/Cambridge_Project/GFFs_genomes/HF_gffs'
dictionary_path = '/home/albertotr/OneDrive/Data/Cambridge_Project/Camille_replacements_foldername.pickle'
# Load dictionary for gene name replacements (if needed)
with open(dictionary_path, 'rb') as f:
    gene_name_dict = pickle.load(f)
# Load Panaroo presence/absence (numeric matrix),
# Panaroo final graph (gene cluster adjacency),
# and pangenome reference (not used in this script but left for completeness).
presence_absence_df = pd.read_table(os.path.join(panaroo_data, 'gene_presence_absence_filt_pseudo_length_frag.Rtab'), sep='\t', header=0)
final_graph_data = nx.read_gml(os.path.join(panaroo_data, 'pre_filt_graph.gml'))
pangenome_reference= SeqIO.parse(os.path.join(panaroo_data, 'pan_genome_reference.fa'), 'fasta')

# GWAS results file to analyze; change to target your association results
GWAS_file = 'mortality_COGs.txt'

# Read and filter GWAS table for significant associations
df = pd.read_table(GWAS_file, sep='\t', header=0)
df_filtered = df[(df['notes'] != 'bad-chisq') & (df['filter-pvalue'] < 0.05)]

# Identify the cluster column in your GWAS table (adjust if needed)
cluster_col = next((c for c in ['variant', 'cluster', 'cluster_id', 'COG', 'COG_id', 'panaroo_cluster']
                    if c in df_filtered.columns), None)
if cluster_col is None:
    raise ValueError("No cluster column found in mortality_COGs.txt. Expected one of: variant, cluster, cluster_id, COG, COG_id, panaroo_cluster")

# Extract significant Panaroo clusters from GWAS
sig_clusters = set(df_filtered[cluster_col].astype(str))

# Map GWAS beta direction to presence/absence (Positive=present, Negative=absent)
cluster_beta_map = {}
cluster_state_map = {}
beta_col = next((c for c in ['beta', 'Beta', 'BETA', 'effect', 'Effect'] if c in df_filtered.columns), None)
if beta_col:
    # If multiple rows per cluster, average beta (sign retained in most cases)
    cluster_beta_map = df_filtered.groupby(cluster_col)[beta_col].mean(numeric_only=True).to_dict()
    cluster_state_map = {cl: ('present' if b >= 0 else 'absent')
                         for cl, b in cluster_beta_map.items()
                         if pd.notnull(b)}
else:
    print("[info] No beta/effect column found; GWAS presence/absence direction will not be annotated.")

#%% Synteny/context analysis around GWAS-significant clusters

OUT_DIR = 'synteny_context_'+ os.path.splitext(os.path.basename(GWAS_file))[0]
os.makedirs(OUT_DIR, exist_ok=True)


def parse_gff_coordinates_bcbiogff(gff_path: str, wanted_ids: set) -> dict:
    """
    Parse a GFF with BCBio to collect coordinates for specific IDs of interest.

    Args:
        gff_path: Path to the GFF file.
        wanted_ids: Set of IDs to match against qualifiers (ID, locus_tag, gene, Parent).

    Returns:
        Dict mapping matched_id -> (start, end, strand, product).

    Notes:
        - Uses 1-based start coordinates and inclusive end.
        - Performs direct ID matching and a substring fallback if direct match fails.
        - Feature types considered: CDS, gene, ORF, protein.
    """
    out = {}
    if not os.path.isfile(gff_path):
        return out
    try:
        with open(gff_path) as handle:
            for record in GFF.parse(handle):
                for feat in record.features:
                    if feat.type not in ('CDS','gene','ORF','protein'):
                        continue
                    start = int(feat.location.start) + 1  # 1-based
                    end = int(feat.location.end)
                    strand = '+' if feat.location.strand == 1 else '-'
                    product = feat.qualifiers.get('product', [''])[0]
                    # Collect candidate IDs from common qualifier keys
                    cand_ids = set()
                    for key in ('ID','locus_tag','gene','Parent'):
                        vals = feat.qualifiers.get(key)
                        if vals:
                            for v in vals:
                                cand_ids.add(str(v))
                    matched = wanted_ids & cand_ids
                    # Also allow substring match if direct match fails (rare)
                    if not matched:
                        for wid in wanted_ids:
                            if any(wid in cid for cid in cand_ids):
                                matched = {wid}
                                break
                    for mid in matched:
                        if mid not in out:  # first occurrence wins
                            out[mid] = (start, end, strand, product)
    except Exception as e:
        print(f"[warn] BCBio parse failed for {gff_path}: {e}")
    return out


def parse_gff_full(gff_path: str):
    """
    Parse all CDS-like features from a GFF into contig-sorted features and an ID lookup.

    Args:
        gff_path: Path to the GFF file.

    Returns:
        contig_index: dict[str, list[dict]] mapping contig -> sorted list of features.
        id_lookup: dict[str, tuple[str, int]] mapping any gene identifier -> (contig, index).

    Feature dict fields:
        - start, end (int): 1-based coordinates.
        - strand (str): '+' or '-'.
        - product (str): feature product/annotation.
        - ids (set[str]): collected identifiers (ID, locus_tag, gene, Parent).

    Notes:
        - Only features of types CDS, gene, ORF, protein are considered.
        - Features are sorted by start then end per contig to enable neighborhood lookup.
    """
    contig_index = defaultdict(list)
    id_lookup = {}
    if not os.path.isfile(gff_path):
        return {}, {}
    try:
        with open(gff_path) as handle:
            for record in GFF.parse(handle):
                contig = str(record.id)
                for feat in record.features:
                    if feat.type not in ('CDS', 'gene', 'ORF', 'protein'):
                        continue
                    start = int(feat.location.start) + 1  # 1-based
                    end = int(feat.location.end)
                    strand = '+' if feat.location.strand == 1 else '-'
                    product = feat.qualifiers.get('product', [''])[0]
                    # Collect candidate IDs from qualifiers
                    cand_ids = set()
                    for key in ('ID', 'locus_tag', 'gene', 'Parent'):
                        vals = feat.qualifiers.get(key)
                        if vals:
                            for v in vals:
                                cand_ids.add(str(v))
                    contig_index[contig].append({
                        'start': start, 'end': end, 'strand': strand,
                        'product': product, 'ids': cand_ids
                    })
            # Sort features and build a fast ID lookup across contigs
            for contig, feats in contig_index.items():
                feats.sort(key=lambda x: (x['start'], x['end']))
                for idx, f in enumerate(feats):
                    for gid in f['ids']:
                        if gid not in id_lookup:
                            id_lookup[gid] = (contig, idx)
    except Exception as e:
        print(f"[warn] GFF full-parse failed for {gff_path}: {e}")
    return contig_index, id_lookup


def neighborhood(contig_index, contig: str, idx: int, window: int = 3):
    """
    Return left/right neighborhood (±window genes) around a feature index on a contig.

    Args:
        contig_index: dict from parse_gff_full containing features per contig.
        contig: contig name.
        idx: index of the target feature within contig_index[contig].
        window: number of genes to include on each side.

    Returns:
        (left_labels, target_label, right_labels) where each label dict includes:
        - id: one identifier string for the feature (or 'NA').
        - product: product annotation.
        - strand: '+' or '-'.
        - start, end: integer coordinates.
    """
    feats = contig_index.get(contig, [])
    left = feats[max(0, idx - window):idx]
    right = feats[idx + 1: idx + 1 + window]
    target = feats[idx] if 0 <= idx < len(feats) else None

    def to_label(f):
        # Use any stable identifier-like string; fall back to product
        if not f:
            return {'id': 'NA', 'product': 'NA', 'strand': 'NA', 'start': None, 'end': None}
        gid = next(iter(f['ids']), 'NA') if f['ids'] else 'NA'
        return {'id': gid, 'product': f['product'], 'strand': f['strand'],
                'start': f['start'], 'end': f['end']}

    left_labels = [to_label(f) for f in left]
    right_labels = [to_label(f) for f in right]
    target_label = to_label(target)
    return left_labels, target_label, right_labels


def first_nonempty_gene_id(cell: str):
    """
    Extract the first gene ID token from a Panaroo CSV cell.

    Args:
        cell: CSV cell content which may contain multi-ID strings separated by ';' or ','.

    Returns:
        The first non-empty token (str) or None if not found.
    """
    if not isinstance(cell, str):
        return None
    # Split on delimiters commonly seen in Panaroo outputs
    for sep in [';', ',']:
        if sep in cell:
            parts = [p.strip() for p in cell.split(sep) if p.strip()]
            return parts[0] if parts else None
    return cell.strip() if cell.strip() else None


def _detect_nonunique_col(pa_csv: pd.DataFrame) -> str | None:
    """
    Detect the column name in Panaroo CSV that holds non-unique gene names.

    Args:
        pa_csv: Panaroo gene_presence_absence CSV dataframe.

    Returns:
        Column name if found, else None. Uses common variants and a fuzzy heuristic.
    """
    # Try common Panaroo names
    candidates = [
        'Non-unique Gene name', 'Non-unique gene name',
        'non_unique_gene_name', 'non-unique gene name'
    ]
    for c in candidates:
        if c in pa_csv.columns:
            return c
    # Fallback: fuzzy search
    for c in pa_csv.columns:
        lc = c.lower()
        if 'non' in lc and 'unique' in lc and 'gene' in lc and 'name' in lc:
            return c
    return None


def _split_multi(cell) -> list[str]:
    """
    Split a cell containing multiple values separated by Panaroo-like delimiters.

    Args:
        cell: Value potentially containing separators like '~~~', ';', ',', '|'.

    Returns:
        List of trimmed tokens.
    """
    if not isinstance(cell, str):
        return []
    # Split on Panaroo-like delimiters and general separators
    seps = ['~~~', ';', ',', '|']
    parts = [cell]
    for sep in seps:
        new_parts = []
        for p in parts:
            new_parts.extend(p.split(sep))
        parts = new_parts
    out = [p.strip() for p in parts if p and p.strip()]
    return out


def _build_gene_name_to_clusters(pa_csv: pd.DataFrame, nonunique_col: str | None) -> dict[str, set]:
    """
    Build a lookup from gene name aliases to Panaroo cluster IDs.

    Args:
        pa_csv: Panaroo CSV dataframe.
        nonunique_col: Column name that contains non-unique gene names.

    Returns:
        dict[name_lower] -> set of Panaroo cluster IDs (strings).
    """
    name_to_clusters = defaultdict(set)
    if nonunique_col and nonunique_col in pa_csv.columns:
        for _, row in pa_csv.iterrows():
            cl = str(row['Gene'])
            names = _split_multi(row[nonunique_col])
            for n in names:
                name_to_clusters[n.lower()].add(cl)
    return name_to_clusters


def _build_sample_gene_to_cluster(pa_csv: pd.DataFrame, sample_cols: list[str]) -> dict[str, dict[str, str]]:
    """
    Build reverse index mapping sample-specific gene IDs to Panaroo clusters.

    Args:
        pa_csv: Panaroo CSV dataframe.
        sample_cols: Columns corresponding to sample names (must match GFF basenames).

    Returns:
        dict[sample] -> dict[gene_id] -> cluster_id
    """
    # sample -> {gene_id -> cluster_id}
    idx = {}
    for _, row in pa_csv.iterrows():
        cl = str(row['Gene'])
        for s in sample_cols:
            cell = row.get(s)
            for gid in _split_multi(cell):
                if not gid:
                    continue
                idx.setdefault(s, {})
                # first mapping wins; Panaroo IDs should be unique within a sample
                idx[s].setdefault(gid, cl)
    return idx


def _resolve_gwas_clusters(sig_terms: set[str], pa_csv: pd.DataFrame) -> list[str]:
    """
    Resolve GWAS terms to actual Panaroo cluster IDs present in the CSV.

    Strategy:
      1) Direct match against pa_csv['Gene'] values (e.g., group_386).
      2) Tokenize GWAS terms and resolve via non-unique gene-name aliases.

    Args:
        sig_terms: Set of terms extracted from the GWAS table.
        pa_csv: Panaroo CSV dataframe.

    Returns:
        Sorted list of resolved cluster IDs (strings). Prints info if none resolved.
    """
    gene_col = pa_csv['Gene'].astype(str)
    all_clusters = set(gene_col)
    nonunique_col = _detect_nonunique_col(pa_csv)
    name_to_clusters = _build_gene_name_to_clusters(pa_csv, nonunique_col)

    resolved = set()
    for term in sig_terms:
        term = str(term).strip()
        if not term:
            continue
        # 1) Direct cluster match (e.g., group_386)
        if term in all_clusters:
            resolved.add(term)
            continue
        # 2) Tokenized name matches from GWAS (handles glpK_1~~~glpK_2~~~glpK, metG)
        tokens = _split_multi(term)
        if not tokens:
            tokens = [term]
        for tok in tokens:
            cand = name_to_clusters.get(tok.lower())
            if cand:
                resolved.update(cand)
    if not resolved:
        print("[info] No GWAS terms resolved to Panaroo clusters. Check GWAS column content and Panaroo CSV.")
    return sorted(resolved)


def _rtab_presence_lookup(rtab: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    """
    Normalize Rtab presence/absence matrix for fast lookup.

    Args:
        rtab: Raw Rtab dataframe loaded from Panaroo.

    Returns:
        (df, sample_cols)
        - df: DataFrame indexed by 'Gene' or first column, values coerced to int (0/1+).
        - sample_cols: list of sample column names.
    """
    # Ensure 'Gene' is the index; detect sample columns
    df = rtab.copy()
    if 'Gene' in df.columns:
        df = df.set_index('Gene')
    else:
        # assume first column is Gene if unnamed
        if df.columns.size > 0 and df.columns[0].lower() != 'gene':
            df = df.set_index(df.columns[0])
    # Normalize numeric to int 0/1+ (any positive -> 1 when interpreted logically)
    df = df.apply(pd.to_numeric, errors='coerce').fillna(0).astype(int)
    return df, list(df.columns)


def analyze_synteny_for_sig_clusters(
    sig_clusters: set,
    window: int = 3,
    include_absent: bool = False,
    include_gwas: bool = False,
    sample_name_map: dict | None = None,
):
    """
    Perform synteny/context analysis around GWAS-significant clusters.

    Steps:
      - Load Panaroo CSV and determine samples overlapping GFF basenames.
      - Build reverse index: sample-specific gene IDs -> cluster IDs.
      - Resolve GWAS terms to Panaroo cluster IDs.
      - Parse GFFs and locate cluster-specific genes per sample.
      - Extract ±window neighborhoods and record IDs/products/cluster mappings.
      - Compare observed neighbor clusters with Panaroo GML graph neighbors.
      - Write per-cluster synteny context and signature summaries.

    Args:
        sig_clusters: Set of terms/IDs from the GWAS, to be resolved to Panaroo clusters.
        window: Number of neighbors to include on each side of the target gene.
        include_absent: If False, exclude rows for samples where the target is absent/missing.
        include_gwas: If False, drop GWAS-related columns from outputs.
        sample_name_map: Optional mapping of sample names to display names for output.
    """
    # Load Panaroo CSV
    pa_csv_path = os.path.join(panaroo_data, 'gene_presence_absence_filt_pseudo_length_frag.csv')
    try:
        pa_csv = pd.read_csv(pa_csv_path)
    except Exception as e:
        raise RuntimeError(f"Could not read {pa_csv_path}: {e}")
    if 'Gene' not in pa_csv.columns:
        raise RuntimeError("Panaroo CSV missing 'Gene' column.")

    # Determine GFFs (sample name assumed to be GFF basename without extension)
    gff_files = {os.path.splitext(f)[0]: os.path.join(gff_folder, f)
                 for f in os.listdir(gff_folder) if f.endswith('.gff')}
    # Intersect samples between CSV and GFFs to ensure we only process samples with GFFs
    sample_cols = [c for c in pa_csv.columns if c in gff_files]
    if not sample_cols:
        print("[warn] No overlap between Panaroo CSV columns and GFF basenames.")
    # Reverse index: sample gene ID -> cluster
    sample_gene_to_cluster = _build_sample_gene_to_cluster(pa_csv, sample_cols)

    # Resolve GWAS terms to clusters
    norm_clusters = _resolve_gwas_clusters(sig_clusters, pa_csv)
    if not norm_clusters:
        return

    # Presence/absence lookup from Rtab
    rtab_df, rtab_samples = _rtab_presence_lookup(presence_absence_df)
    # Restrict to available sample columns if possible; otherwise use all Rtab samples
    rtab_sample_cols = [s for s in rtab_samples if s in sample_cols] or rtab_samples

    # GML graph neighbors: build neighbor sets per cluster node
    gml_neighbors = {}
    try:
        for n in final_graph_data.nodes:
            try:
                gml_neighbors[str(n)] = set(str(nb) for nb in final_graph_data.neighbors(n))
            except Exception:
                pass
    except Exception as e:
        print(f"[warn] Could not derive GML neighbors: {e}")

    # Pre-index all GFFs once for efficiency
    gff_index_cache = {}
    for sample, gff_path in gff_files.items():
        contig_index, id_lookup = parse_gff_full(gff_path)
        gff_index_cache[sample] = (contig_index, id_lookup)

    # Detect non-unique gene-name column for reporting
    nonunique_col = _detect_nonunique_col(pa_csv)

    # Prepare sample display-name mapping (use loaded dictionary if available)
    if sample_name_map is None:
        try:
            sample_name_map = gene_name_dict if isinstance(gene_name_dict, dict) else {}
        except NameError:
            sample_name_map = {}

    # Iterate clusters and collect per-sample contexts
    for cl in norm_clusters:
        # Row for the cluster in Panaroo CSV
        row = pa_csv.loc[pa_csv['Gene'].astype(str) == cl]
        if row.empty:
            # Not in CSV; still write empty with info
            print(f"[info] Cluster {cl} not found in Panaroo CSV; skipping.")
            continue

        contexts = []
        product_sigs = []
        cluster_sigs = []

        for sample in sample_cols:
            # Output display name for figures/TSVs
            sample_disp = sample_name_map.get(sample, sample)
            # Get Panaroo gene ID(s) in this sample for this cluster
            gene_cell = row.iloc[0][sample]
            gene_ids = _split_multi(gene_cell)
            gene_id = gene_ids[0] if gene_ids else None

            # Locate gene in GFF via ID lookup (with substring fallback)
            contig_index, id_lookup = gff_index_cache.get(sample, ({}, {}))
            loc = None
            if gene_id:
                loc = id_lookup.get(gene_id)
                if loc is None:
                    # substring fallback for naming inconsistencies
                    matches = [v for k, v in id_lookup.items() if gene_id in k]
                    loc = matches[0] if matches else None

            # Presence from Rtab (0/1) if available
            present_rtab = None
            if cl in rtab_df.index and sample in rtab_df.columns:
                present_rtab = int(rtab_df.at[cl, sample] > 0)

            # GWAS direction annotations
            gwas_beta = cluster_beta_map.get(cl, None)
            gwas_state = cluster_state_map.get(cl, '')
            if present_rtab is None or not gwas_state:
                matches_dir = ''
            else:
                matches_dir = int((gwas_state == 'present' and present_rtab == 1) or
                                  (gwas_state == 'absent' and present_rtab == 0))

            if loc is None:
                # No GFF hit for this sample: annotate status and presence flags
                status = 'absent_in_gff' if gene_id else 'absent_in_csv_and_gff'
                presence_flag = 'present_rtab_only' if present_rtab == 1 else ('absent_both' if present_rtab == 0 else 'unknown')
                contexts.append({
                    'cluster': cl,
                    'sample': sample_disp,
                    'status': status,
                    'presence_rtab': present_rtab,
                    'presence_flag': presence_flag,
                    'contig': 'NA',
                    'target_id': gene_id or '',
                    'target_product': '',
                    'left_ids': '',
                    'left_products': '',
                    'left_clusters': '',
                    'right_ids': '',
                    'right_products': '',
                    'right_clusters': '',
                    'near_contig_end': '',
                    'gml_neighbors_overlap': '',
                    'non_unique_gene_name': (row.iloc[0][nonunique_col] if nonunique_col else ''),
                    'gwas_beta': gwas_beta,
                    'gwas_presence_state': gwas_state,
                    'matches_gwas_direction': matches_dir
                })
                continue

            # Extract neighborhood around target
            contig, idx_in_contig = loc
            feats = contig_index.get(contig, [])
            left, target, right = neighborhood(contig_index, contig, idx_in_contig, window=window)

            left_ids = [x['id'] for x in left]
            right_ids = [x['id'] for x in right]
            left_prods = [x['product'] for x in left]
            right_prods = [x['product'] for x in right]

            # Map neighbors to Panaroo clusters via reverse index
            smap = sample_gene_to_cluster.get(sample, {})
            left_clusters = [smap.get(gid, '') for gid in left_ids]
            right_clusters = [smap.get(gid, '') for gid in right_ids]

            # Build signatures (products and clusters)
            product_sig = tuple(left_prods + [target['product']] + right_prods)
            cluster_sig = tuple([c for c in left_clusters if c] + [cl] + [c for c in right_clusters if c])
            product_sigs.append(product_sig)
            if cluster_sig:
                cluster_sigs.append(cluster_sig)

            # Presence flag combining GFF and Rtab
            present_gff = 1
            if present_rtab is None:
                presence_flag = 'present_gff_only_unknown_rtab'
            elif present_rtab == 1:
                presence_flag = 'present_both'
            else:
                presence_flag = 'present_gff_only'

            # Contig-end heuristic: fewer neighbors than window implies near end
            near_end = int(len(left) < window or len(right) < window)

            # GML neighbor overlap with observed neighbor clusters
            gml_neigh = gml_neighbors.get(cl, set())
            obs_neigh = set([c for c in left_clusters + right_clusters if c])
            overlap = len(gml_neigh & obs_neigh)
            overlap_str = f"{overlap}/{len(obs_neigh)}" if obs_neigh else "0/0"

            # Collect context row for output
            contexts.append({
                'cluster': cl,
                'sample': sample_disp,
                'status': 'ok',
                'presence_rtab': present_rtab,
                'presence_flag': presence_flag,
                'contig': contig,
                'target_id': target['id'],
                'target_product': target['product'],
                'left_ids': ';'.join(left_ids),
                'left_products': ';'.join(left_prods),
                'left_clusters': ';'.join([c for c in left_clusters if c]),
                'right_ids': ';'.join(right_ids),
                'right_products': ';'.join(right_prods),
                'right_clusters': ';'.join([c for c in right_clusters if c]),
                'near_contig_end': near_end,
                'gml_neighbors_overlap': overlap_str,
                'non_unique_gene_name': (row.iloc[0][nonunique_col] if nonunique_col else ''),
                'gwas_beta': gwas_beta,
                'gwas_presence_state': gwas_state,
                'matches_gwas_direction': matches_dir
            })

        # Optionally remove absent rows to avoid "Absent" annotations in figures
        if not include_absent:
            contexts = [r for r in contexts if r['status'] == 'ok']

        # Derive consensus neighbor clusters from present samples
        present_rows = [r for r in contexts if r['status'] == 'ok']
        neighbor_cluster_counts = Counter()
        for r in present_rows:
            for c_ in (r['left_clusters'].split(';') + r['right_clusters'].split(';')):
                if c_:
                    neighbor_cluster_counts[c_] += 1
        # Expected neighbors: those seen in >=20% of present samples (tunable)
        expected_neighbors = {c for c, cnt in neighbor_cluster_counts.items()
                              if present_rows and (cnt / len(present_rows)) >= 0.2}

        # Classify absence rows (only if we kept them)
        if include_absent:
            for r in contexts:
                if r['status'] == 'ok':
                    r['absence_class'] = ''
                    r['missing_neighbor_fraction'] = ''
                    continue
                # Skip if we have no consensus
                if not expected_neighbors:
                    r['absence_class'] = 'undetermined'
                    r['missing_neighbor_fraction'] = ''
                    continue
                # Presence of expected neighbor clusters in Rtab
                sample = r['sample']
                present_flags = []
                for neigh in expected_neighbors:
                    if neigh in rtab_df.index and sample in rtab_df.columns:
                        present_flags.append(int(rtab_df.at[neigh, sample] > 0))
                if not present_flags:
                    r['absence_class'] = 'undetermined'
                    r['missing_neighbor_fraction'] = ''
                    continue
                frac_missing = 1 - (sum(present_flags) / len(present_flags))
                r['missing_neighbor_fraction'] = round(frac_missing, 3)
                if r['presence_flag'] == 'present_rtab_only':
                    r['absence_class'] = 'annotation_issue'
                elif frac_missing <= 0.2:
                    r['absence_class'] = 'target_only_missing_or_misannotated'
                elif frac_missing >= 0.8:
                    r['absence_class'] = 'region_missing'
                else:
                    r['absence_class'] = 'partial_region_loss'

        # Write per-sample contexts for this cluster
        ctx_df = pd.DataFrame(contexts)

        # Drop GWAS columns if requested
        if not include_gwas:
            ctx_df.drop(columns=['gwas_beta', 'gwas_presence_state', 'matches_gwas_direction'], errors='ignore', inplace=True)

        # Drop presence/absence annotation columns if we excluded absent rows
        if not include_absent:
            ctx_df.drop(columns=['presence_rtab', 'presence_flag', 'absence_class', 'missing_neighbor_fraction'], errors='ignore', inplace=True)
        ctx_path = os.path.join(OUT_DIR, f'{cl}_synteny.tsv')
        ctx_df.to_csv(ctx_path, sep='\t', index=False)

        # Consensus neighbor summary
        cons_rows = []
        for c_, cnt in neighbor_cluster_counts.most_common():
            cons_rows.append({'cluster': cl,
                              'neighbor_cluster': c_,
                              'count_present_samples': cnt,
                              'frequency': round(cnt / len(present_rows), 4)})
        cons_df = pd.DataFrame(cons_rows)
        cons_df.to_csv(os.path.join(OUT_DIR, f'{cl}_consensus_context.tsv'), sep='\t', index=False)

        # Signature frequencies (products)
        sig_counts = Counter(product_sigs)
        sig_rows = [{'cluster': cl, 'count': c, 'signature_products': ';'.join(s)} for s, c in sig_counts.items()]
        sig_df = pd.DataFrame(sig_rows).sort_values(by='count', ascending=False)
        sig_path = os.path.join(OUT_DIR, f'{cl}_signatures.tsv')
        sig_df.to_csv(sig_path, sep='\t', index=False)

        # Signature frequencies (clusters), if any cluster-level signatures recorded
        if cluster_sigs:
            c_counts = Counter(cluster_sigs)
            c_rows = [{'cluster': cl, 'count': c, 'signature_clusters': ';'.join(s)} for s, c in c_counts.items()]
            c_df = pd.DataFrame(c_rows).sort_values(by='count', ascending=False)
            c_path = os.path.join(OUT_DIR, f'{cl}_cluster_signatures.tsv')
            c_df.to_csv(c_path, sep='\t', index=False)

        print(f"[info] Wrote {ctx_path}, {sig_path}" + (f" and {c_path}" if cluster_sigs else ""))

# === Run synteny/context analysis ===
# Use a ±3-gene window by default; absent/GWAS annotations are disabled; sample names remapped if dictionary provided
analyze_synteny_for_sig_clusters(sig_clusters, window=3)

