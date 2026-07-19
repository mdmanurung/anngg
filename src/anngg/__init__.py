"""anngg: a ggplot2-style plotting layer for scanpy AnnData objects.

The core grammar *is* plotnine. ``gganndata(adata) + aes(...) + geom_*()`` returns
a real :class:`plotnine.ggplot`, with aesthetics resolved from ``AnnData`` through
annplyr. High-level ``plot_*`` helpers reproduce scanpy's core figures, and
``plot_clustermap`` is a separate escape hatch for clustered heatmaps via
PyComplexHeatmap.
"""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

from ._resolve import embedding_coords, gene, obs, obsm
from .clustermap import plot_clustermap
from .grammar import aes, gganndata
from .plots import plot_dotplot, plot_embedding, plot_matrixplot, plot_violin
from .theme import (
    scale_color_celltype,
    scale_color_expression,
    scale_colour_celltype,
    scale_colour_expression,
    scale_fill_celltype,
    scale_fill_expression,
    theme_anngg,
)

try:
    __version__ = version("anngg")
except PackageNotFoundError:  # pragma: no cover - not installed
    __version__ = "0.1.0"

__all__ = [
    "gganndata",
    "aes",
    "gene",
    "obs",
    "obsm",
    "embedding_coords",
    "plot_embedding",
    "plot_dotplot",
    "plot_matrixplot",
    "plot_violin",
    "plot_clustermap",
    "theme_anngg",
    "scale_color_expression",
    "scale_colour_expression",
    "scale_fill_expression",
    "scale_color_celltype",
    "scale_colour_celltype",
    "scale_fill_celltype",
    "__version__",
]
