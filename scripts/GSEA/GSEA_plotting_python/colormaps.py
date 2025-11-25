"""
Colorblind-safe color definitions for GSEA visualizations.

Provides consistent, accessible color palettes for NES heatmaps
and other GSEA visualizations.
"""

from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt


# Colorblind-safe diverging palette (Blue-White-Orange)
DIVERGING_COLORS = {
    'negative': '#2166AC',  # Blue (downregulated/negative NES)
    'neutral': '#F7F7F7',   # White (neutral)
    'positive': '#B35806',  # Orange (upregulated/positive NES)
}

# Colorblind-safe categorical palette (Okabe-Ito)
CATEGORICAL_COLORS = {
    'blue': '#0072B2',
    'orange': '#D55E00',
    'green': '#009E73',
    'yellow': '#F0E442',
    'sky_blue': '#56B4E9',
    'vermillion': '#E69F00',
    'pink': '#CC79A7',
    'black': '#000000',
    'gray': '#999999',
}


def create_diverging_cmap(name='BlueOrange', n_colors=256):
    """
    Create a Blue-White-Orange colorblind-safe colormap for NES values.

    Parameters
    ----------
    name : str
        Name for the colormap
    n_colors : int
        Number of colors in the colormap

    Returns
    -------
    LinearSegmentedColormap
        Colorblind-safe diverging colormap
    """
    return LinearSegmentedColormap.from_list(
        name,
        [DIVERGING_COLORS['negative'],
         DIVERGING_COLORS['neutral'],
         DIVERGING_COLORS['positive']],
        N=n_colors
    )


def get_colorblind_palette(n_colors=8):
    """
    Get a colorblind-safe categorical palette.

    Parameters
    ----------
    n_colors : int
        Number of colors needed (max 8)

    Returns
    -------
    list
        List of hex color codes
    """
    palette = [
        '#0072B2',  # Blue
        '#D55E00',  # Orange
        '#009E73',  # Green
        '#F0E442',  # Yellow
        '#56B4E9',  # Sky blue
        '#E69F00',  # Vermillion/amber
        '#CC79A7',  # Pink
        '#999999',  # Gray
    ]
    return palette[:n_colors]


def register_custom_cmaps():
    """Register custom colormaps with matplotlib."""
    cmap = create_diverging_cmap('BlueOrange')
    try:
        plt.cm.register_cmap(name='BlueOrange', cmap=cmap)
    except ValueError:
        # Already registered
        pass
