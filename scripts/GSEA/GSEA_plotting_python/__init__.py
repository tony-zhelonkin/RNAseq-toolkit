"""
GSEA Plotting Utilities for Python

Provides colorblind-safe colormaps and reusable heatmap/dotplot rendering
with significance hatching for GSEA visualizations.
"""

from .colormaps import (
    DIVERGING_COLORS,
    create_diverging_cmap,
    get_colorblind_palette
)

from .heatmap_renderer import HeatmapRenderer
from .dotplot_renderer import DotplotRenderer

__all__ = [
    'DIVERGING_COLORS',
    'create_diverging_cmap',
    'get_colorblind_palette',
    'HeatmapRenderer',
    'DotplotRenderer'
]
