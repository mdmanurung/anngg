"""Additional marker-summary plots: stacked violin and tracksplot.

Both stack genes as rows (``facet_grid`` over the feature) so many markers can be
compared across groups compactly, mirroring ``sc.pl.stacked_violin`` /
``sc.pl.tracksplot``. Expression is pulled long through annplyr's ``to_tidy``.
"""

from __future__ import annotations

from typing import Sequence

import pandas as pd
import plotnine_extra as pe
from plotnine import (
    aes,
    element_blank,
    element_text,
    facet_grid,
    geom_col,
    geom_violin,
    ggplot,
    labs,
    theme,
)

from ._aggregate import tidy_expression
from ._palette import scale_fill_obs
from .plots import _group_categories, _order_groups, plot_dotplot, plot_matrixplot
from .theme import theme_anngg

__all__ = [
    "plot_stacked_violin",
    "plot_tracksplot",
    "plot_dotplot_grouped",
    "plot_matrixplot_grouped",
]


def _flatten_gene_groups(gene_groups: dict) -> tuple[list, dict]:
    """Flatten a ``{label: [genes]}`` mapping, rejecting a gene listed in two groups."""
    feat_to_group: dict = {}
    genes: list = []
    for label, members in gene_groups.items():
        for g in members:
            if g in feat_to_group:
                raise ValueError(
                    f"gene {g!r} appears in more than one gene group "
                    f"({feat_to_group[g]!r} and {label!r}); each gene must belong to one group."
                )
            feat_to_group[g] = label
            genes.append(g)
    return genes, feat_to_group


def _add_gene_group_facet(plot, gene_groups: dict):
    """Bracket the gene axis into labelled groups via a ``gene_group`` facet column.

    ``gene_groups`` maps a label -> list of genes; the panels appear in dict order.
    Reuses the plot's own aggregated ``.data`` (which carries a ``feature`` column).
    """
    _, feat_to_group = _flatten_gene_groups(gene_groups)
    data = plot.data.copy()
    data["gene_group"] = pd.Categorical(
        data["feature"].astype(str).map(feat_to_group),
        categories=list(gene_groups),
        ordered=True,
    )
    plot.data = data
    return plot + facet_grid(". ~ gene_group", scales="free_x", space="free_x")


def plot_dotplot_grouped(adata, gene_groups: dict, group_by: str, **kwargs):
    """Dotplot with genes bracketed into labelled groups (scanpy ``var_group_labels``)."""
    genes, _ = _flatten_gene_groups(gene_groups)
    return _add_gene_group_facet(plot_dotplot(adata, genes, group_by, **kwargs), gene_groups)


def plot_matrixplot_grouped(adata, gene_groups: dict, group_by: str, **kwargs):
    """Matrixplot with genes bracketed into labelled groups (scanpy ``var_group_labels``)."""
    genes, _ = _flatten_gene_groups(gene_groups)
    return _add_gene_group_facet(plot_matrixplot(adata, genes, group_by, **kwargs), gene_groups)


def plot_stacked_violin(
    adata,
    genes: Sequence[str],
    group_by: str,
    *,
    layer: str | None = None,
    use_raw: bool | None = None,
    scale: str = "width",
    categories_order=None,
):
    """Compact genes-as-rows violin grid across groups (``sc.pl.stacked_violin``)."""
    genes = list(genes)
    tidy = tidy_expression(adata, genes, group_by, layer=layer, use_raw=use_raw)
    tidy = _order_groups(tidy, group_by, categories_order or _group_categories(adata, group_by))
    return (
        ggplot(tidy, aes(group_by, "value", fill=group_by))
        + geom_violin(scale=scale)
        + facet_grid("feature ~ .", scales="free_y")
        + scale_fill_obs(adata, group_by)
        + labs(x="", y="", fill=group_by)
        + theme_anngg()
        + theme(strip_text_y=element_text(angle=0))
        + pe.rotate_x_text(45)
    )


def plot_tracksplot(
    adata,
    genes: Sequence[str],
    group_by: str,
    *,
    layer: str | None = None,
    use_raw: bool | None = None,
    categories_order=None,
):
    """Per-gene expression tracks across cells ordered by group (``sc.pl.tracksplot``)."""
    genes = list(genes)
    tidy = tidy_expression(adata, genes, group_by, layer=layer, use_raw=use_raw)
    tidy = _order_groups(tidy, group_by, categories_order or _group_categories(adata, group_by))

    # order cells by group so each track reads left-to-right by group
    cell = tidy[["obs_name", group_by]].drop_duplicates().sort_values(group_by)
    cell = cell.reset_index(drop=True)
    cell["cell_rank"] = range(len(cell))
    tidy = tidy.merge(cell[["obs_name", "cell_rank"]], on="obs_name")

    return (
        ggplot(tidy, aes("cell_rank", "value", fill=group_by))
        + geom_col(width=1.0)
        + facet_grid("feature ~ .", scales="free_y")
        + scale_fill_obs(adata, group_by)
        + labs(x="cells (ordered by group)", y="", fill=group_by)
        + theme_anngg()
        + theme(
            axis_text_x=element_blank(),
            axis_ticks_major_x=element_blank(),
            strip_text_y=element_text(angle=0),
        )
    )
