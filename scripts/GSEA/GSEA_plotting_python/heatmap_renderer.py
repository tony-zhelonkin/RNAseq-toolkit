"""
Reusable heatmap renderer with significance hatching.

Provides a standardized way to render GSEA heatmaps with:
- Diagonal hatching for non-significant cells (padj >= 0.05), OR
- Values shown for all cells (no hatching, asterisk only for significant)
- Consistent text styling (black, bold for significant)
- Colorblind-safe colormaps
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize

from .colormaps import create_diverging_cmap


class HeatmapRenderer:
    """
    Render heatmaps with diagonal hatching or values-only for non-significant cells.

    Parameters
    ----------
    cmap : Colormap, optional
        Matplotlib colormap to use. Default is Blue-White-Orange.
    vmax : float
        Maximum absolute value for color scaling. Default 3.5.
    hatch_color : str
        Background color for hatched cells. Default '#F0F0F0'.
    hatch_pattern : str
        Matplotlib hatch pattern. Default '///'.
    font_size : int
        Font size for cell annotations. Default 7.
    nonsig_style : str
        Style for non-significant cells: 'hatching' (default) or 'show_values'.
        - 'hatching': Diagonal lines over gray background, no values shown
        - 'show_values': No hatching, show NES values without asterisk
          (color naturally fades as NES approaches 0)
    """

    def __init__(self, cmap=None, vmax=3.5, hatch_color='#F0F0F0',
                 hatch_pattern='///', font_size=7, nonsig_style='hatching'):
        self.cmap = cmap if cmap is not None else create_diverging_cmap()
        self.vmax = vmax
        self.hatch_color = hatch_color
        self.hatch_pattern = hatch_pattern
        self.font_size = font_size
        self.nonsig_style = nonsig_style
        self.norm = Normalize(vmin=-vmax, vmax=vmax)

    def render(self, ax, value_matrix, sig_matrix=None,
               row_labels=None, col_labels=None,
               show_values=True, show_hatching=True):
        """
        Render a complete heatmap with hatching and annotations.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes to render on
        value_matrix : np.ndarray
            2D array of values (e.g., NES scores)
        sig_matrix : np.ndarray, optional
            2D boolean array where True = significant (padj < 0.05)
            If None, all cells are considered significant (no hatching)
        row_labels : list, optional
            Labels for rows (y-axis)
        col_labels : list, optional
            Labels for columns (x-axis)
        show_values : bool
            Whether to show cell values
        show_hatching : bool
            Whether to add hatching for non-significant cells.
            Ignored if nonsig_style='show_values'.

        Returns
        -------
        matplotlib.image.AxesImage
            The heatmap image object
        """
        # Create heatmap
        im = ax.imshow(value_matrix, cmap=self.cmap, norm=self.norm,
                       aspect='auto')

        # Add hatching for non-significant cells (only in 'hatching' mode)
        if show_hatching and sig_matrix is not None and self.nonsig_style == 'hatching':
            self.add_hatching(ax, value_matrix, sig_matrix)

        # Add cell annotations
        if show_values:
            self.annotate_cells(ax, value_matrix, sig_matrix)

        # Set labels
        if row_labels is not None:
            ax.set_yticks(range(len(row_labels)))
            ax.set_yticklabels(row_labels)

        if col_labels is not None:
            ax.set_xticks(range(len(col_labels)))
            ax.set_xticklabels(col_labels, rotation=45, ha='right')

        return im

    def add_hatching(self, ax, value_matrix, sig_matrix):
        """
        Add diagonal hatching for non-significant cells.

        Hatching is added when:
        - Cell value is NA (np.nan)
        - Cell is not significant (sig_matrix[i,j] == False)

        Note: In GSEA exports, NA typically means the pathway was tested
        but did not reach significance threshold (padj >= 0.05).

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes to render on
        value_matrix : np.ndarray
            2D array of values
        sig_matrix : np.ndarray
            2D boolean array where True = significant
        """
        nrows, ncols = value_matrix.shape

        for i in range(nrows):
            for j in range(ncols):
                # Hatch if NA OR not significant
                should_hatch = (np.isnan(value_matrix[i, j]) or
                               (sig_matrix is not None and not sig_matrix[i, j]))

                if should_hatch:
                    rect = Rectangle(
                        (j - 0.5, i - 0.5), 1, 1,
                        fill=True,
                        facecolor=self.hatch_color,
                        hatch=self.hatch_pattern,
                        edgecolor='gray',
                        linewidth=0.5,
                        zorder=2
                    )
                    ax.add_patch(rect)

    def annotate_cells(self, ax, value_matrix, sig_matrix=None, show_missing_as='NT'):
        """
        Add cell annotations based on nonsig_style setting.

        When nonsig_style='hatching' (default):
        - Only significant cells (padj < 0.05) are annotated
        - Hatched cells (non-significant or NA) are left empty
        - Significant cells show: value with asterisk (e.g., "2.1*")

        When nonsig_style='show_values':
        - All cells with values are annotated
        - Significant cells show: value with asterisk, bold (e.g., "2.1*")
        - Non-significant cells show: value without asterisk, normal weight
        - Missing values (NaN) show: 'NT' (not tested) in gray

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes to render on
        value_matrix : np.ndarray
            2D array of values
        sig_matrix : np.ndarray, optional
            2D boolean array where True = significant
        show_missing_as : str, optional
            Text to show for missing (NaN) values. Default 'NT' (not tested).
            Set to None to hide missing values.
        """
        nrows, ncols = value_matrix.shape

        for i in range(nrows):
            for j in range(ncols):
                value = value_matrix[i, j]

                if np.isnan(value):
                    # Handle missing values (NaN)
                    if self.nonsig_style == 'show_values' and show_missing_as:
                        # Show "NT" for not tested
                        ax.text(
                            j, i, show_missing_as,
                            ha='center', va='center',
                            color='#888888',
                            fontweight='normal',
                            fontsize=self.font_size - 1,
                            fontstyle='italic',
                            zorder=3
                        )
                else:
                    is_sig = (sig_matrix is not None and sig_matrix[i, j])

                    if self.nonsig_style == 'show_values':
                        # Show all values: asterisk and bold for significant only
                        if is_sig:
                            text = f'{value:.1f}*'
                            fontweight = 'bold'
                        else:
                            text = f'{value:.1f}'
                            fontweight = 'normal'
                        ax.text(
                            j, i, text,
                            ha='center', va='center',
                            color='black',
                            fontweight=fontweight,
                            fontsize=self.font_size,
                            zorder=3
                        )
                    else:
                        # Default 'hatching' mode: only annotate significant cells
                        if is_sig:
                            text = f'{value:.1f}*'
                            ax.text(
                                j, i, text,
                                ha='center', va='center',
                                color='black',
                                fontweight='bold',
                                fontsize=self.font_size,
                                zorder=3
                            )

    def add_colorbar(self, fig, im, ax=None, label='NES'):
        """
        Add a colorbar to the figure.

        Parameters
        ----------
        fig : matplotlib.figure.Figure
            Figure to add colorbar to
        im : matplotlib.image.AxesImage
            Image object from render()
        ax : matplotlib.axes.Axes, optional
            Axes for colorbar placement
        label : str
            Colorbar label

        Returns
        -------
        matplotlib.colorbar.Colorbar
            The colorbar object
        """
        cbar = fig.colorbar(im, ax=ax, shrink=0.8, aspect=30)
        cbar.set_label(label, fontsize=10)
        return cbar

    def add_legend(self, ax, x=0.02, y=-0.15):
        """
        Add a standard legend explaining hatching/values and significance.

        Legend text adapts based on nonsig_style setting.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes to add legend to
        x, y : float
            Position in axes coordinates
        """
        if self.nonsig_style == 'show_values':
            legend_text = (
                'Bold* = Significant (padj < 0.05)\n'
                'Normal = Not significant'
            )
        else:
            legend_text = (
                'Bold* = Significant (padj < 0.05)\n'
                '/// = Not significant'
            )
        ax.text(
            x, y, legend_text,
            transform=ax.transAxes,
            fontsize=8,
            verticalalignment='top',
            family='monospace',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
        )


def create_significance_matrix(padj_matrix, cutoff=0.05):
    """
    Create a boolean significance matrix from p-adjusted values.

    Parameters
    ----------
    padj_matrix : np.ndarray
        2D array of adjusted p-values
    cutoff : float
        Significance cutoff (default 0.05)

    Returns
    -------
    np.ndarray
        Boolean array where True = significant (padj < cutoff)
    """
    return padj_matrix < cutoff
