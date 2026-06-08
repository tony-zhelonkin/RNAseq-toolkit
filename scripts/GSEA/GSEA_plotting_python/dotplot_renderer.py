"""
Reusable dotplot renderer for GSEA visualizations.

Provides a standardized way to render GSEA dotplots with:
- Dot color representing NES value (diverging colormap)
- Dot size representing significance (-log10(padj))
- Black outline for significant dots (padj < 0.05)
- Colorblind-safe colormaps
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

from .colormaps import create_diverging_cmap


class DotplotRenderer:
    """
    Render dotplots with color for NES and size for significance.

    Parameters
    ----------
    cmap : Colormap, optional
        Matplotlib colormap to use. Default is Blue-White-Orange.
    vmax : float
        Maximum absolute value for color scaling. Default 3.5.
    font_size : int
        Font size for cell annotations. Default 7.
    min_dot_size : float
        Minimum dot size (for non-significant or high padj). Default 10.
    max_dot_size : float
        Maximum dot size (for most significant, low padj). Default 200.
    sig_cutoff : float
        Significance cutoff for black outline. Default 0.05.
    """

    def __init__(self, cmap=None, vmax=3.5, font_size=7,
                 min_dot_size=10, max_dot_size=200, sig_cutoff=0.05):
        self.cmap = cmap if cmap is not None else create_diverging_cmap()
        self.vmax = vmax
        self.font_size = font_size
        self.min_dot_size = min_dot_size
        self.max_dot_size = max_dot_size
        self.sig_cutoff = sig_cutoff
        self.norm = Normalize(vmin=-vmax, vmax=vmax)

    def _calculate_dot_sizes(self, padj_matrix):
        """
        Calculate dot sizes based on -log10(padj).

        Parameters
        ----------
        padj_matrix : np.ndarray
            2D array of adjusted p-values

        Returns
        -------
        np.ndarray
            Dot sizes scaled between min_dot_size and max_dot_size
        """
        # Calculate -log10(padj), handling NaN and 0 values
        with np.errstate(divide='ignore', invalid='ignore'):
            neg_log_padj = -np.log10(padj_matrix)

        # Replace inf (from padj=0) with a large value
        neg_log_padj[np.isinf(neg_log_padj)] = 10

        # For NaN values (missing data), use minimum size
        neg_log_padj[np.isnan(neg_log_padj)] = 0

        # Normalize to [0, 1] range
        # Use -log10(sig_cutoff) as baseline (e.g., -log10(0.05) ≈ 1.3)
        baseline = -np.log10(self.sig_cutoff)

        # Clip values: anything below baseline gets minimum size
        neg_log_padj_clipped = np.clip(neg_log_padj, 0, 10)

        # Scale to dot size range
        # Linear scaling: 0 -> min_size, baseline+ -> max_size
        normalized = neg_log_padj_clipped / 10  # Normalize to [0, 1]
        dot_sizes = self.min_dot_size + normalized * (self.max_dot_size - self.min_dot_size)

        return dot_sizes

    def render(self, ax, value_matrix, padj_matrix,
               row_labels=None, col_labels=None):
        """
        Render a complete dotplot with colored dots.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes to render on
        value_matrix : np.ndarray
            2D array of values (e.g., NES scores)
        padj_matrix : np.ndarray
            2D array of adjusted p-values
        row_labels : list, optional
            Labels for rows (y-axis)
        col_labels : list, optional
            Labels for columns (x-axis)

        Returns
        -------
        matplotlib.collections.PathCollection
            The scatter plot collection
        """
        nrows, ncols = value_matrix.shape

        # Calculate dot sizes
        dot_sizes = self._calculate_dot_sizes(padj_matrix)

        # Prepare data for scatter plot
        x_coords = []
        y_coords = []
        colors = []
        sizes = []
        edgecolors = []
        linewidths = []

        for i in range(nrows):
            for j in range(ncols):
                value = value_matrix[i, j]
                padj = padj_matrix[i, j]

                # Skip NaN values
                if np.isnan(value):
                    continue

                x_coords.append(j)
                y_coords.append(i)

                # Color by NES value
                colors.append(value)

                # Size by significance
                sizes.append(dot_sizes[i, j])

                # Black edge for significant, gray for non-significant
                if not np.isnan(padj) and padj < self.sig_cutoff:
                    edgecolors.append('black')
                    linewidths.append(1.5)
                else:
                    edgecolors.append('gray')
                    linewidths.append(0.5)

        # Create scatter plot
        scatter = ax.scatter(
            x_coords, y_coords,
            c=colors,
            s=sizes,
            cmap=self.cmap,
            norm=self.norm,
            edgecolors=edgecolors,
            linewidths=linewidths,
            alpha=0.9,
            zorder=3
        )

        # Set axis properties
        ax.set_xlim(-0.5, ncols - 0.5)
        ax.set_ylim(-0.5, nrows - 0.5)
        ax.invert_yaxis()  # Match heatmap orientation

        # Add grid for clarity
        ax.set_xticks(range(ncols))
        ax.set_yticks(range(nrows))

        # Set labels
        if row_labels is not None:
            ax.set_yticklabels(row_labels)

        if col_labels is not None:
            ax.set_xticklabels(col_labels)

        return scatter

    def add_colorbar(self, fig, scatter, ax=None, label='NES'):
        """
        Add a colorbar to the figure.

        Parameters
        ----------
        fig : matplotlib.figure.Figure
            Figure to add colorbar to
        scatter : matplotlib.collections.PathCollection
            Scatter plot from render()
        ax : matplotlib.axes.Axes, optional
            Axes for colorbar placement
        label : str
            Colorbar label

        Returns
        -------
        matplotlib.colorbar.Colorbar
            The colorbar object
        """
        cbar = fig.colorbar(scatter, ax=ax, shrink=0.8, aspect=30)
        cbar.set_label(label, fontsize=10)
        return cbar

    def add_size_legend(self, ax, x=1.05, y=0.5):
        """
        Add a legend explaining dot sizes.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes to add legend to
        x, y : float
            Position in axes coordinates
        """
        # Create dummy scatter plots for legend
        padj_values = [0.001, 0.01, 0.05, 0.1]
        with np.errstate(divide='ignore', invalid='ignore'):
            neg_log_padj = [-np.log10(p) for p in padj_values]

        # Calculate corresponding sizes
        baseline = -np.log10(self.sig_cutoff)
        sizes = []
        for nlp in neg_log_padj:
            normalized = np.clip(nlp / 10, 0, 1)
            size = self.min_dot_size + normalized * (self.max_dot_size - self.min_dot_size)
            sizes.append(size)

        # Create legend elements
        legend_elements = []
        for padj, size in zip(padj_values, sizes):
            legend_elements.append(
                plt.scatter([], [], s=size, c='gray', edgecolors='black',
                           linewidths=1, alpha=0.9,
                           label=f'padj = {padj:.3f}')
            )

        ax.legend(handles=legend_elements, title='Dot Size (padj)',
                 loc='center left', bbox_to_anchor=(x, y),
                 frameon=True, framealpha=0.9)

    def add_legend(self, ax, x=0.02, y=-0.15):
        """
        Add a standard legend explaining the dotplot.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes to add legend to
        x, y : float
            Position in axes coordinates
        """
        legend_text = (
            'Dot color = NES value\n'
            'Dot size = -log10(padj)\n'
            'Black edge = Significant (padj < 0.05)'
        )
        ax.text(
            x, y, legend_text,
            transform=ax.transAxes,
            fontsize=8,
            verticalalignment='top',
            family='monospace',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
        )
