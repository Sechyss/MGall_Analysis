#%% Load packages and data (Change accordingly)

import pandas as pd
import os
import networkx as nx
from Bio import SeqIO
from BCBio import GFF
from collections import Counter, defaultdict

os.chdir('/home/albertotr/OneDrive/Data/Cambridge_Project/GWAS/')

panaroo_data = '/home/albertotr/OneDrive/Data/Cambridge_Project/Pangenome_results_HF'
gff_folder = '/home/albertotr/OneDrive/Data/Cambridge_Project/GFFs_genomes/HF_gffs'

presence_absence_df = pd.read_table(os.path.join(panaroo_data, 'gene_presence_absence_filt_pseudo_length_frag.Rtab'), sep='\t', header=0)
final_graph_data = nx.read_gml(os.path.join(panaroo_data, 'pre_filt_graph.gml'))
pangenome_reference= SeqIO.parse(os.path.join(panaroo_data, 'pan_genome_reference.fa'), 'fasta')
GWAS_file = 'mortality_COGs.txt'

df = pd.read_table(GWAS_file, sep='\t', header=0)
df_filtered = df[(df['notes'] != 'bad-chisq') & (df['filter-pvalue'] < 0.05)]

# Identify the cluster column in your GWAS table (adjust if needed)
cluster_col = next((c for c in ['variant', 'cluster', 'cluster_id', 'COG', 'COG_id', 'panaroo_cluster']
                    if c in df_filtered.columns), None)
if cluster_col is None:
    raise ValueError("No cluster column found in mortality_COGs.txt. Expected one of: variant, cluster, cluster_id, COG, COG_id, panaroo_cluster")

# Extract significant Panaroo clusters from GWAS
sig_clusters = set(df_filtered[cluster_col].astype(str))

#%% Synteny/context analysis around GWAS-significant clusters

OUT_DIR = 'synteny_context_'+ os.path.splitext(os.path.basename(GWAS_file))[0]
os.makedirs(OUT_DIR, exist_ok=True)


def parse_gff_coordinates_bcbiogff(gff_path: str, wanted_ids: set) -> dict:
    """
    Use BCBio.GFF to parse CDS features.
    Returns dict[matched_id] -> (start,end,strand,product).
    Matches any of feature.qualifiers keys: ID, locus_tag, gene, Parent against wanted_ids.
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
                    # Collect candidate IDs
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
                        if mid not in out:  # first occurrence
                            out[mid] = (start, end, strand, product)
    except Exception as e:
        print(f"[warn] BCBio parse failed for {gff_path}: {e}")
    return out


def parse_gff_full(gff_path: str):
    """
    Parse all CDS-like features from a GFF into:
      - contig_index: dict[contig] -> list of dict(feature fields) sorted by start
      - id_lookup: dict[id_string] -> (contig, index_in_contig_list)
    Candidate IDs include qualifiers: ID, locus_tag, gene, Parent
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
                    # Collect candidate IDs
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
            # sort and build lookup
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
    Return left/right neighborhood (±window genes) around index on a contig.
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
    Panaroo CSV cells often have semi-colon separated gene IDs (or NaN). Take first ID if present.
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
    name_to_clusters = defaultdict(set)
    if nonunique_col and nonunique_col in pa_csv.columns:
        for _, row in pa_csv.iterrows():
            cl = str(row['Gene'])
            names = _split_multi(row[nonunique_col])
            for n in names:
                name_to_clusters[n.lower()].add(cl)
    return name_to_clusters


def _build_sample_gene_to_cluster(pa_csv: pd.DataFrame, sample_cols: list[str]) -> dict[str, dict[str, str]]:
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
    # Resolve GWAS names to Panaroo cluster IDs present in pa_csv['Gene']
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
    # Ensure 'Gene' is the index; detect sample columns
    df = rtab.copy()
    if 'Gene' in df.columns:
        df = df.set_index('Gene')
    else:
        # assume first column is Gene if unnamed
        if df.columns.size > 0 and df.columns[0].lower() != 'gene':
            df = df.set_index(df.columns[0])
    # Normalize numeric
    df = df.apply(pd.to_numeric, errors='coerce').fillna(0).astype(int)
    return df, list(df.columns)


def analyze_synteny_for_sig_clusters(sig_clusters: set, window: int = 3):
    """
    Extended version:
      - robustly resolve GWAS terms to Panaroo clusters (no 'cluster_' prefix assumption)
      - use Rtab presence/absence to mark presence and discordance with GFF
      - map neighbor CDS back to Panaroo clusters
      - compare observed neighbors with Panaroo GML graph neighbors
    """
    # Load Panaroo CSV
    pa_csv_path = os.path.join(panaroo_data, 'gene_presence_absence_filt_pseudo_length_frag.csv')
    try:
        pa_csv = pd.read_csv(pa_csv_path)
    except Exception as e:
        raise RuntimeError(f"Could not read {pa_csv_path}: {e}")
    if 'Gene' not in pa_csv.columns:
        raise RuntimeError("Panaroo CSV missing 'Gene' column.")

    # Determine GFFs
    gff_files = {os.path.splitext(f)[0]: os.path.join(gff_folder, f)
                 for f in os.listdir(gff_folder) if f.endswith('.gff')}
    # Intersect samples between CSV and GFFs
    sample_cols = [c for c in pa_csv.columns if c in gff_files]
    if not sample_cols:
        print("[warn] No overlap between Panaroo CSV columns and GFF basenames.")
    # Reverse index: sample gene ID -> cluster
    sample_gene_to_cluster = _build_sample_gene_to_cluster(pa_csv, sample_cols)

    # Resolve GWAS terms to clusters
    norm_clusters = _resolve_gwas_clusters(sig_clusters, pa_csv)
    if not norm_clusters:
        return

    # Presence/absence lookup
    rtab_df, rtab_samples = _rtab_presence_lookup(presence_absence_df)
    # restrict to available sample columns if possible
    rtab_sample_cols = [s for s in rtab_samples if s in sample_cols] or rtab_samples

    # GML graph neighbors
    gml_neighbors = {}
    try:
        for n in final_graph_data.nodes:
            try:
                gml_neighbors[str(n)] = set(str(nb) for nb in final_graph_data.neighbors(n))
            except Exception:
                pass
    except Exception as e:
        print(f"[warn] Could not derive GML neighbors: {e}")

    # Pre-index all GFFs once
    gff_index_cache = {}
    for sample, gff_path in gff_files.items():
        contig_index, id_lookup = parse_gff_full(gff_path)
        gff_index_cache[sample] = (contig_index, id_lookup)

    # Detect non-unique gene-name column for reporting
    nonunique_col = _detect_nonunique_col(pa_csv)

    for cl in norm_clusters:
        # Row for the cluster
        row = pa_csv.loc[pa_csv['Gene'].astype(str) == cl]
        if row.empty:
            # Not in CSV; still write empty with info
            print(f"[info] Cluster {cl} not found in Panaroo CSV; skipping.")
            continue

        contexts = []
        product_sigs = []
        cluster_sigs = []

        for sample in sample_cols:
            # Get Panaroo gene ID in this sample for this cluster
            gene_cell = row.iloc[0][sample]
            gene_ids = _split_multi(gene_cell)
            gene_id = gene_ids[0] if gene_ids else None

            contig_index, id_lookup = gff_index_cache.get(sample, ({}, {}))
            loc = None
            if gene_id:
                loc = id_lookup.get(gene_id)
                if loc is None:
                    # substring fallback
                    matches = [v for k, v in id_lookup.items() if gene_id in k]
                    loc = matches[0] if matches else None

            present_rtab = None
            if cl in rtab_df.index and sample in rtab_df.columns:
                present_rtab = int(rtab_df.at[cl, sample] > 0)

            if loc is None:
                # No GFF hit for this sample
                status = 'absent_in_gff' if gene_id else 'absent_in_csv_and_gff'
                presence_flag = 'present_rtab_only' if present_rtab == 1 else ('absent_both' if present_rtab == 0 else 'unknown')
                contexts.append({
                    'cluster': cl,
                    'sample': sample,
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
                    'non_unique_gene_name': (row.iloc[0][nonunique_col] if nonunique_col else '')
                })
                continue

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

            # Build signatures
            product_sig = tuple(left_prods + [target['product']] + right_prods)
            cluster_sig = tuple([c for c in left_clusters if c] + [cl] + [c for c in right_clusters if c])
            product_sigs.append(product_sig)
            if cluster_sig:
                cluster_sigs.append(cluster_sig)

            # Presence flag
            present_gff = 1
            if present_rtab is None:
                presence_flag = 'present_gff_only_unknown_rtab'
            elif present_rtab == 1:
                presence_flag = 'present_both'
            else:
                presence_flag = 'present_gff_only'

            # Contig-end heuristic
            near_end = int(len(left) < window or len(right) < window)

            # GML neighbor overlap
            gml_neigh = gml_neighbors.get(cl, set())
            obs_neigh = set([c for c in left_clusters + right_clusters if c])
            overlap = len(gml_neigh & obs_neigh)
            overlap_str = f"{overlap}/{len(obs_neigh)}" if obs_neigh else "0/0"

            contexts.append({
                'cluster': cl,
                'sample': sample,
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
                'non_unique_gene_name': (row.iloc[0][nonunique_col] if nonunique_col else '')
            })

        # Write per-sample contexts
        ctx_df = pd.DataFrame(contexts)
        ctx_path = os.path.join(OUT_DIR, f'{cl}_synteny.tsv')
        ctx_df.to_csv(ctx_path, sep='\t', index=False)

        # Signature frequencies (products)
        sig_counts = Counter(product_sigs)
        sig_rows = [{'cluster': cl, 'count': c, 'signature_products': ';'.join(s)} for s, c in sig_counts.items()]
        sig_df = pd.DataFrame(sig_rows).sort_values(by='count', ascending=False)
        sig_path = os.path.join(OUT_DIR, f'{cl}_signatures.tsv')
        sig_df.to_csv(sig_path, sep='\t', index=False)

        # Signature frequencies (clusters)
        if cluster_sigs:
            c_counts = Counter(cluster_sigs)
            c_rows = [{'cluster': cl, 'count': c, 'signature_clusters': ';'.join(s)} for s, c in c_counts.items()]
            c_df = pd.DataFrame(c_rows).sort_values(by='count', ascending=False)
            c_path = os.path.join(OUT_DIR, f'{cl}_cluster_signatures.tsv')
            c_df.to_csv(c_path, sep='\t', index=False)

        print(f"[info] Wrote {ctx_path}, {sig_path}" + (f" and {c_path}" if cluster_sigs else ""))

# === Run synteny/context analysis ===
# Use a ±3-gene window by default; adjust as needed
analyze_synteny_for_sig_clusters(sig_clusters, window=3)

