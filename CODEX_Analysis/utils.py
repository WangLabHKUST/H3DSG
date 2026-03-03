import pandas as pd
def join_comb(stat_df, celltype_order, class_name='status'):
    
    stat_df[f'i_{class_name}'] = pd.Categorical(stat_df[f'i_{class_name}'], categories=celltype_order,ordered=True)
    stat_df[f'j_{class_name}'] = pd.Categorical(stat_df[f'j_{class_name}'], categories=celltype_order,ordered=True)

    ### i>j
    i_j_mask = stat_df[f'i_{class_name}'] <= stat_df[f'j_{class_name}']
    j_i_mask = stat_df[f'i_{class_name}'] > stat_df[f'j_{class_name}']
    assert (sum(i_j_mask) + sum(j_i_mask)) == stat_df.shape[0]
    stat_df.loc[i_j_mask, 'comb'] = stat_df.loc[i_j_mask, f'i_{class_name}'].astype(str) + '...' + stat_df.loc[i_j_mask, f'j_{class_name}'].astype(str)
    stat_df.loc[j_i_mask, 'comb'] = stat_df.loc[j_i_mask, f'j_{class_name}'].astype(str) + '...' + stat_df.loc[j_i_mask, f'i_{class_name}'].astype(str)

    return stat_df

def add_large_legend(ax=None, 
                     hue_order=None, 
                     palette=None, 
                     data=None, 
                     hue='cell_type',
                     title='Cell Type',
                     marker_size=150,
                     loc='center left',
                     bbox_to_anchor=(1.02, 0.5),
                     ncol=1,
                     fontsize=12,
                     title_fontsize=14,
                     frameon=True,
                     edgecolor='black',
                     linewidth=0.8):
    """Add a standalone legend with larger markers for scatter plots."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import seaborn as sns
    if ax is None:
        ax = plt.gca()
    
    # Determine categories and colors
    if hue_order is None:
        categories = data[hue].unique()
    else:
        categories = hue_order
    
    if palette is None:
        raise ValueError("必须提供 palette 参数（颜色字典）")
    
    # Support both dict palettes and list-like palettes
    if isinstance(palette, dict):
        colors = [palette[cat] for cat in categories]
    else:
        # Seaborn palettes can be list-like or dict-like
        colors = palette if isinstance(palette, list) else list(palette.values())
        colors = colors[:len(categories)]  # Guard against length mismatch
    
    # Build custom legend handles
    handles = [
        Line2D([0], [0], 
               marker='o', 
               color='w', 
               label=cat,
               markerfacecolor=color, 
               markersize=(marker_size ** 0.5),  # Approximate marker area for Line2D
               markeredgecolor=edgecolor,
               markeredgewidth=linewidth)
        for cat, color in zip(categories, colors)
    ]
    
    # Add legend to axis
    ax.legend(handles=handles,
              title=title,
              title_fontsize=title_fontsize,
              fontsize=fontsize,
              loc=loc,
              bbox_to_anchor=bbox_to_anchor,
              ncol=ncol,
              frameon=frameon)
    
    return ax
