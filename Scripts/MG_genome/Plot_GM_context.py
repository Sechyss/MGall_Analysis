import os
import re
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

def _read_cluster_files(ctx_dir: str, cluster: str):
    ctx_path = os.path.join(ctx_dir, f"{cluster}_synteny.tsv")
    prod_sig_path = os.path.join(ctx_dir, f"{cluster}_signatures.tsv")
    clust_sig_path = os.path.join(ctx_dir, f"{cluster}_cluster_signatures.tsv")
    if not os.path.isfile(ctx_path):
        raise FileNotFoundError(f"Missing context file: {ctx_path}")
    ctx_df = pd.read_csv(ctx_path, sep="\t")
    prod_sig = pd.read_csv(prod_sig_path, sep="\t") if os.path.isfile(prod_sig_path) else pd.DataFrame()
    clust_sig = pd.read_csv(clust_sig_path, sep="\t") if os.path.isfile(clust_sig_path) else pd.DataFrame()
    return ctx_df, prod_sig, clust_sig

def _split_semicolon(s):
    if pd.isna(s) or not isinstance(s, str) or s == '':
        return []
    return [x for x in s.split(';') if x != '']

def _auto_window(ctx_df: pd.DataFrame) -> int:
    # Infer window from the longest left/right lists
    lmax = ctx_df['left_products'].dropna().map(lambda x: len(_split_semicolon(x))).max() if 'left_products' in ctx_df else 0
    rmax = ctx_df['right_products'].dropna().map(lambda x: len(_split_semicolon(x))).max() if 'right_products' in ctx_df else 0
    w = int(max(lmax or 0, rmax or 0))
    return max(w, 1)

def _sanitize_token(t):
    if t is None or t == '' or str(t).lower() == 'na':
        return ''
    return str(t)

def _build_tile_matrix(ctx_df: pd.DataFrame, cluster: str, mode: str, window: int):
    """
    Returns:
      samples_order: list[str]
      tokens_matrix: list[list[str]] with columns [-w..-1,0,+1..+w]
      legend_tokens: set[str] used
    mode in {'products','clusters'}
    """
    rows = []
    for _, row in ctx_df.iterrows():
        sample = row['sample']
        left = _split_semicolon(row['left_products' if mode == 'products' else 'left_clusters'])
        right = _split_semicolon(row['right_products' if mode == 'products' else 'right_clusters'])
        # Normalize to fixed window by padding on the left/right
        left = ([''] * (window - len(left))) + left[-window:]
        right = right[:window] + ([''] * (window - len(right)))
        center = _sanitize_token(row['target_product'] if mode == 'products' else cluster)
        tokens = left + [center] + right
        rows.append((sample, tokens, row.get('status',''), row.get('presence_flag','unknown')))
    # Sort rows for readability: status -> presence_flag -> tokens signature
    def sort_key(x):
        sample, tokens, status, pflag = x
        return (0 if status == 'ok' else 1, pflag, tuple(tokens))
    rows.sort(key=sort_key)
    samples_order = [r[0] for r in rows]
    matrix = [r[1] for r in rows]
    legend = set()
    for toks in matrix:
        legend.update([t for t in toks if t])
    return samples_order, matrix, legend

def _palette_for_tokens(tokens: list[str]):
    # Build a color map for up to many tokens; reserve gray for blanks
    uniq = [t for t in tokens if t]
    # Use tab20 then husl to extend
    base = sns.color_palette('tab20', n_colors=min(20, max(1, len(uniq))))
    if len(uniq) > 20:
        base += sns.color_palette('husl', n_colors=len(uniq) - 20)
    color_map = {'': (0.92, 0.92, 0.92)}
    for t, c in zip(uniq, base):
        color_map[t] = c
    return color_map

def _draw_tiles(out_png: str, cluster: str, mode: str, samples: list[str], matrix: list[list[str]], window: int, ctx_df: pd.DataFrame, color_map: dict):
    # Prepare array for imshow: map tokens -> colors
    n_rows = len(samples)
    if n_rows == 0:
        return
    n_cols = 2 * window + 1
    # Build RGB image
    img = np.ones((n_rows, n_cols, 3), dtype=float)
    for i, toks in enumerate(matrix):
        for j, tok in enumerate(toks):
            img[i, j, :] = color_map.get(tok, (0.6, 0.6, 0.6))
    # Figure
    plt.figure(figsize=(max(6, n_cols * 0.5), max(4, n_rows * 0.25)))
    ax = plt.gca()
    ax.imshow(img, aspect='auto', interpolation='nearest')
    # Axes ticks/labels
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(samples, fontsize=7)
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels([f"{i}" for i in range(-window, window + 1)], fontsize=8)
    ax.set_xlabel("Position relative to target (0 = target)")
    ax.set_title(f"{cluster} synteny ({mode})")
    # Grid lines
    ax.set_xticks(np.arange(-.5, n_cols, 1), minor=True)
    ax.set_yticks(np.arange(-.5, n_rows, 1), minor=True)
    ax.grid(which='minor', color='w', linestyle='-', linewidth=0.3, alpha=0.7)
    # Legend (limit tokens)
    used_tokens = []
    for row in matrix:
        for t in row:
            if t and t not in used_tokens:
                used_tokens.append(t)
    max_legend = min(25, len(used_tokens))
    handles = [plt.Line2D([0], [0], marker='s', color='w', markerfacecolor=color_map[t], markersize=8, label=(t[:40] + ('…' if len(t) > 40 else ''))) for t in used_tokens[:max_legend]]
    if handles:
        ax.legend(handles=handles, bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0., title="Tokens")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()

def _draw_signature_bars(out_png: str, cluster: str, prod_sig: pd.DataFrame, clust_sig: pd.DataFrame, top_n=12):
    fig, axes = plt.subplots(1, 2, figsize=(14, 4))
    # Products
    if not prod_sig.empty:
        df = prod_sig.copy().sort_values('count', ascending=False).head(top_n)
        axes[0].barh(range(len(df)), df['count'], color='#4C78A8')
        axes[0].set_yticks(range(len(df)))
        axes[0].set_yticklabels([s[:80] + ('…' if len(s) > 80 else '') for s in df['signature_products']])
        axes[0].invert_yaxis()
        axes[0].set_title('Product signatures')
    else:
        axes[0].axis('off')
        axes[0].set_title('Product signatures (none)')
    # Clusters
    if not clust_sig.empty:
        df = clust_sig.copy().sort_values('count', ascending=False).head(top_n)
        axes[1].barh(range(len(df)), df['count'], color='#F58518')
        axes[1].set_yticks(range(len(df)))
        axes[1].set_yticklabels([s[:80] + ('…' if len(s) > 80 else '') for s in df['signature_clusters']])
        axes[1].invert_yaxis()
        axes[1].set_title('Cluster signatures')
    else:
        axes[1].axis('off')
        axes[1].set_title('Cluster signatures (none)')
    fig.suptitle(f"{cluster} signature frequencies", y=1.02)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200, bbox_inches='tight')
    plt.close()

def _write_status_table(out_tsv: str, ctx_df: pd.DataFrame):
    cols = ['sample', 'status', 'presence_flag', 'gml_neighbors_overlap', 'near_contig_end', 'contig', 'target_id', 'target_product']
    ctx_df[cols].to_csv(out_tsv, sep='\t', index=False)

def list_clusters(ctx_dir: str):
    return sorted(set(re.sub(r'_synteny\.tsv$', '', f) for f in os.listdir(ctx_dir) if f.endswith('_synteny.tsv')))

def main():
    ap = argparse.ArgumentParser(description="Plot synteny context outputs.")
    ap.add_argument('--ctx-dir', required=True, help='Path to synteny_context_ directory')
    ap.add_argument('--out-dir', required=False, default=None, help='Output directory (default: ctx-dir/plots)')
    ap.add_argument('--clusters', nargs='*', default=None, help='Specific clusters to plot (default: all)')
    ap.add_argument('--mode', choices=['both', 'products', 'clusters'], default='both', help='Tile plot mode')
    ap.add_argument('--top-n-sigs', type=int, default=12, help='Top signatures to show in bars')
    args = ap.parse_args()

    ctx_dir = args.ctx_dir
    out_dir = args.out_dir or os.path.join(ctx_dir, 'plots')
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    clusters = args.clusters or list_clusters(ctx_dir)
    if not clusters:
        raise SystemExit(f"No *_synteny.tsv files found under {ctx_dir}")

    for cl in clusters:
        try:
            ctx_df, prod_sig, clust_sig = _read_cluster_files(ctx_dir, cl)
        except Exception as e:
            print(f"[warn] Skipping {cl}: {e}")
            continue

        # Status summary table
        _write_status_table(os.path.join(out_dir, f"{cl}_sample_status.tsv"), ctx_df)

        # Signature bars
        _draw_signature_bars(os.path.join(out_dir, f"{cl}_signature_bars.png"), cl, prod_sig, clust_sig, args.top_n_sigs)

        # Determine window
        window = _auto_window(ctx_df)

        # Tile plots (products/clusters)
        if args.mode in ('both', 'products'):
            samples, matrix, legend = _build_tile_matrix(ctx_df[ctx_df['status'] == 'ok'].copy(), cl, 'products', window)
            cmap = _palette_for_tokens(sorted(list(legend)))
            _draw_tiles(os.path.join(out_dir, f"{cl}_synteny_tiles_products.png"), cl, 'products', samples, matrix, window, ctx_df, cmap)

        if args.mode in ('both', 'clusters'):
            samples, matrix, legend = _build_tile_matrix(ctx_df[ctx_df['status'] == 'ok'].copy(), cl, 'clusters', window)
            cmap = _palette_for_tokens(sorted(list(legend)))
            _draw_tiles(os.path.join(out_dir, f"{cl}_synteny_tiles_clusters.png"), cl, 'clusters', samples, matrix, window, ctx_df, cmap)

        print(f"[info] Wrote plots for {cl} -> {out_dir}")

if __name__ == "__main__":
    main()