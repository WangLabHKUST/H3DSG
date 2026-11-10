# -*- coding: utf-8 -*-
"""Drawing code shared across figures/: the dot matrix and the paired-statistics grid.

It lives here so that the two min-distance dot matrices use one geometry and one
colour/size mapping. They appear side by side in the same figure plate, so their sidebars
(colorbar + size legend) must be the same size -- and because a sidebar is anchored in
ax.transAxes units, **the two axes must have the same absolute dimensions**. Callers
ensure that by arranging cell_w x n_cols and cell_h x n_rows to come out to the same axes
width and height.
"""
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# Colour: the symmetric range of log2(TL / PL)
VMAX_ABS = 1.0
# Size: the bounds on -log10(p); anything beyond is clipped
NLOG10P_MIN = 0.5      # p = 10^-0.5 ~= 0.32
NLOG10P_MAX = 4.0      # p = 1e-4
DOT_SIZE_MIN = 20      # points^2
DOT_SIZE_MAX = 360
SIG_P = 0.05
SIG_EDGE_WIDTH = 0.8

# seaborn's stripplot draws its jitter offsets from numpy's global RNG, so without a seed
# the same script produces a different PNG on every run. Seeding it makes the figures
# byte-reproducible; the value matches the SEED used in pipeline/.
JITTER_SEED = 42


def size_for_p(p, size_scale='radius'):
    """p -> the s of scatter (matplotlib's s is an area, in points^2).

    size_scale='radius' (default): the radius is linear on [sqrt(MIN), sqrt(MAX)], s = r^2.
    size_scale='area': s itself is linear, i.e. area proportional to -log10(p).

    radius is the default because the eye judges diameter rather than area: with a linear
    area the radius grows as a square root, too fast at the start, so a moderately
    significant dot soon looks almost as large as the most significant one. The two
    extremes are the same size under either convention.
    """
    if not np.isfinite(p) or p <= 0:
        val = NLOG10P_MIN
    else:
        val = -np.log10(p)
    val = max(NLOG10P_MIN, min(NLOG10P_MAX, val))
    frac = (val - NLOG10P_MIN) / (NLOG10P_MAX - NLOG10P_MIN)
    if size_scale == 'radius':
        r_min, r_max = np.sqrt(DOT_SIZE_MIN), np.sqrt(DOT_SIZE_MAX)
        return float((r_min + frac * (r_max - r_min)) ** 2)
    return DOT_SIZE_MIN + frac * (DOT_SIZE_MAX - DOT_SIZE_MIN)


def plot_dot_matrix(stats_df, row_order, col_order, out_path, title='',
                    cell_w=0.26, cell_h=0.30,
                    fig_w_min=4.5, fig_h_min=1.0,
                    fig_w_pad=2.1, fig_h_pad=1.45,
                    size_scale='radius', dpi=220,
                    cbar_label='log2(TL / PL)'):
    """Draw a dot matrix: colour = log2_fc, size = -log10(p), black edge where p < 0.05.

    stats_df needs the four columns source / target / log2_fc / mw_p (the output of
    pair_stats).

    Careful: fig_h_min must be pushed below fig_h_pad + n_rows * cell_h, or the height is
    set by fig_h_min, cell_h has no effect at all, and the surplus is spent on row
    spacing -- very noticeable with only three to five rows, where the figure looks tall
    and empty. fig_h_pad reserves a fixed height for the 45-degree x labels and the
    xlabel.

    Careful: with very few rows (2), cell_h must not be too small. The colorbar and the
    size legend are anchored in ax.transAxes units, so a short axes collapses the sidebar
    into a clump with overlapping labels. Measured, the lower bound on cell_h at 2 rows is
    around 0.45.
    """
    sub = stats_df[stats_df['source'].isin(row_order) & stats_df['target'].isin(col_order)]
    n_rows, n_cols = len(row_order), len(col_order)

    fig_w = max(fig_w_min, fig_w_pad + n_cols * cell_w)
    fig_h = max(fig_h_min, fig_h_pad + n_rows * cell_h)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    norm = mcolors.Normalize(vmin=-VMAX_ABS, vmax=VMAX_ABS)
    cmap = plt.get_cmap('RdBu_r')

    xs, ys, sizes, colors, edges, widths = [], [], [], [], [], []
    for _, r in sub.iterrows():
        if not (np.isfinite(r['log2_fc']) and np.isfinite(r['mw_p'])):
            continue
        ys.append(row_order.index(r['source']))
        xs.append(col_order.index(r['target']))
        sizes.append(size_for_p(r['mw_p'], size_scale))
        colors.append(cmap(norm(r['log2_fc'])))
        sig = r['mw_p'] < SIG_P
        edges.append('black' if sig else 'none')
        widths.append(SIG_EDGE_WIDTH if sig else 0.0)

    ax.scatter(xs, ys, s=sizes, c=colors, edgecolors=edges, linewidths=widths)
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(col_order, rotation=45, ha='right', fontsize=8)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(row_order, fontsize=9)
    ax.set_xlim(-0.5, n_cols - 0.5)
    ax.set_ylim(n_rows - 0.5, -0.5)
    ax.set_xlabel('target')
    ax.set_ylabel('source')
    if title:
        ax.set_title(title, fontsize=10)
    ax.grid(False)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)

    fig.patch.set_facecolor('white')
    # tight_layout the main axes first, then anchor the sidebars with insets at a fixed
    # fraction to the right of the axes; calling fig.colorbar(ax=ax) directly would fight
    # ax.legend(bbox_to_anchor=...) over the same right margin.
    plt.tight_layout()

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cax = inset_axes(ax, width='100%', height='100%', loc='lower left',
                     bbox_to_anchor=(1.04, 0.55, 0.04, 0.40),
                     bbox_transform=ax.transAxes, borderpad=0)
    cbar = fig.colorbar(sm, cax=cax, extend='both')
    cbar.set_label(cbar_label, fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    leg_p = [0.5, 0.05, 0.01, 1e-3]
    handles = [plt.scatter([], [], s=size_for_p(p, size_scale), c='grey', edgecolors='none',
                           label=('p=1e-4' if p <= 1e-4 else f'p={p:g}')) for p in leg_p]
    ax.legend(handles=handles, title='-log10(p)',
              bbox_to_anchor=(1.04, 0.45), loc='upper left',
              bbox_transform=ax.transAxes, frameon=False,
              fontsize=7, title_fontsize=8, labelspacing=1.0, borderpad=0.6)

    fig.savefig(out_path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f'saved -> {out_path}')


def pair_box_grid(long_df, pairs, out_path, test_df=None, ncols=4,
                  value_col='log2_d', group_col='group_label',
                  group_order=('PL', 'TL'), group_colors=None,
                  ylabel='log2(avg. min. distance)', dpi=300,
                  panel_w=2.3, panel_h=3.0):
    """One panel per (source, target) pair: a boxplot with jittered points, p value in the
    title.

    pairs    [(source, target, title), ...]; the order is the plotting order
    long_df  one row per sample, needs source / target / value_col / group_col
    test_df  optional, with source / target / mw_p, used to print p in the title; without
             it no p is printed

    Boxplots always use showfliers=False: every panel has jittered points drawn on top, so
    if the boxplot also drew its outliers, the same sample would appear twice (once at the
    outlier position, once at its jittered offset) and read as two samples.
    """
    import pandas as pd
    import seaborn as sns

    group_order = list(group_order)
    if group_colors is None:
        group_colors = {'PL': '#990000', 'TL': '#2b8cbe'}
    if test_df is not None:
        test_lookup = test_df.set_index(['source', 'target'])['mw_p'].to_dict()
    else:
        test_lookup = {}

    np.random.seed(JITTER_SEED)        # see JITTER_SEED above

    n = len(pairs)
    ncols = min(ncols, n)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(panel_w * ncols, panel_h * nrows),
                             sharey=False)
    axes = np.atleast_1d(axes).ravel()

    for ax, (src, tgt, title) in zip(axes, pairs):
        sub = long_df[(long_df['source'] == src) & (long_df['target'] == tgt)]
        sns.boxplot(data=sub, x=group_col, y=value_col, hue=group_col,
                    order=group_order, hue_order=group_order, palette=group_colors,
                    width=0.55, linewidth=1.0, showfliers=False, legend=False, ax=ax,
                    medianprops={'color': 'black', 'linewidth': 1.6})
        sns.stripplot(data=sub, x=group_col, y=value_col, order=group_order,
                      color='black', size=2.6, jitter=0.18, ax=ax)
        p = test_lookup.get((src, tgt))
        ax.set_title(f'{title}\nP = {p:.3g}' if p is not None and np.isfinite(p) else title,
                     fontsize=9)
        ax.set_xlabel('')
        ax.set_ylabel(ylabel if ax is axes[0] else '')
        ax.tick_params(axis='x', labelsize=9)
        ax.tick_params(axis='y', labelsize=8)
        sns.despine(ax=ax, top=True, right=True)

    for ax in axes[n:]:
        ax.axis('off')

    fig.patch.set_facecolor('white')
    plt.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f'saved -> {out_path}')
