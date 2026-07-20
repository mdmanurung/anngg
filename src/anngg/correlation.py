"""Group/sample correlation heatmap from pseudobulk mean-expression profiles.

Inspired by DOTools' ``correlation`` and marsilea's annotated heatmaps: aggregate
mean expression per group (a pseudobulk profile per cell type / sample / batch),
correlate those profiles, and draw the group-by-group correlation as a plotnine
tile plot. Optional hierarchical clustering reorders the axes so related groups
sit together (marsilea/ComplexHeatmap style), and values can be annotated in-cell.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np
import pandas as pd
import plotnine_extra as pe
from plotnine import (
    aes,
    coord_equal,
    element_blank,
    geom_text,
    geom_tile,
    ggplot,
    labs,
    scale_fill_cmap,
    theme,
)

from ._aggregate import aggregate_expression
from .theme import theme_anngg

__all__ = ["plot_correlation"]


def _default_genes(adata) -> list[str]:
    """Highly-variable genes when flagged, else every var (selection, not extraction)."""
    if "highly_variable" in adata.var.columns:
        hv = adata.var_names[np.asarray(adata.var["highly_variable"], dtype=bool)]
        if len(hv) > 1:
            return list(hv)
    return list(adata.var_names)


def _cluster_order(corr: pd.DataFrame) -> list[str]:
    """Leaf order from average-linkage clustering on ``1 - correlation`` distance."""
    from scipy.cluster.hierarchy import leaves_list, linkage
    from scipy.spatial.distance import squareform

    dist = 1.0 - corr.to_numpy()
    np.fill_diagonal(dist, 0.0)
    dist = (dist + dist.T) / 2.0  # enforce symmetry for squareform
    order = leaves_list(linkage(squareform(dist, checks=False), method="average"))
    return [corr.columns[i] for i in order]


def plot_correlation(
    adata,
    group_by: str,
    *,
    genes: Sequence[str] | None = None,
    layer: str | None = None,
    use_raw: bool | None = None,
    method: str = "pearson",
    cluster: bool = True,
    annotate: bool = False,
    cmap: str = "RdBu_r",
):
    """Correlation heatmap between the mean-expression profiles of each group.

    ``genes`` defaults to the highly-variable set (or all genes). ``method`` is any
    :meth:`pandas.DataFrame.corr` method (``"pearson"`` / ``"spearman"`` /
    ``"kendall"``). With ``cluster=True`` the axes are reordered by hierarchical
    clustering; ``annotate=True`` prints each correlation in its cell.
    """
    if genes is None:
        genes = _default_genes(adata)
    genes = list(genes)
    if len(genes) < 2:
        raise ValueError("plot_correlation needs at least two genes to correlate profiles.")

    agg = aggregate_expression(adata, genes, group_by, layer=layer, use_raw=use_raw)
    profiles = agg.pivot(index="feature", columns=group_by, values="mean_expression")
    corr = profiles.corr(method=method)

    order = _cluster_order(corr) if cluster and corr.shape[0] > 2 else list(corr.columns)
    corr = corr.loc[order, order]

    long = (
        corr.rename_axis(index="row", columns="col")
        .reset_index()
        .melt(id_vars="row", var_name="col", value_name="corr")
    )
    long["row"] = pd.Categorical(long["row"], categories=order, ordered=True)
    long["col"] = pd.Categorical(long["col"], categories=order, ordered=True)

    plot = (
        ggplot(long, aes("col", "row", fill="corr"))
        + geom_tile()
        + scale_fill_cmap(cmap_name=cmap, limits=(float(long["corr"].min()), 1.0))
        + coord_equal()
        + labs(x="", y="", fill=f"{method}\ncorrelation")
        + theme_anngg()
        + theme(axis_ticks=element_blank())
        + pe.rotate_x_text(45)
    )
    if annotate:
        plot = plot + geom_text(aes(label="corr"), format_string="{:.2f}", size=7)
    return plot
