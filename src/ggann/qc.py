"""Quality-control plots over ``obs`` metrics and expression.

Early-analysis staples: distributions of per-cell QC metrics, a metric-vs-metric
scatter, and the highest-expressed genes. Metric columns and expression are
pulled through the ``adata.ap`` accessor like everything else.
"""

from __future__ import annotations

import warnings
from typing import Sequence

import numpy as np
import pandas as pd
import plotnine_extra as pe
from plotnine import (
    aes,
    element_blank,
    geom_boxplot,
    geom_violin,
    ggplot,
    labs,
    theme,
)

from ._aggregate import expression_source
from ._palette import scale_fill_obs
from ._resolve import _densify, _raw_wide
from .theme import theme_ggann

__all__ = ["plot_qc_violin", "plot_qc_scatter", "plot_highest_expr_genes"]

# Common QC metric names across scanpy / older workflows.
_DEFAULT_METRICS = [
    "n_genes_by_counts",
    "n_genes",
    "total_counts",
    "n_counts",
    "pct_counts_mt",
    "pct_counts_ribo",
    "percent_mito",
]


def _resolve_metrics(adata, metrics):
    if metrics is not None:
        missing = [m for m in metrics if m not in adata.obs]
        if missing:
            raise KeyError(f"QC metrics not in adata.obs: {missing}")
        return list(metrics)
    found = [m for m in _DEFAULT_METRICS if m in adata.obs]
    if not found:
        raise ValueError(
            "No default QC metrics found in adata.obs; pass metrics=[...] explicitly."
        )
    return found


def plot_qc_violin(
    adata,
    metrics: Sequence[str] | None = None,
    group_by: str | None = None,
    *,
    scale: str = "width",
):
    """Violin distributions of per-cell QC metrics, one facet per metric.

    With ``group_by`` the violins are split by that obs column (stored palette
    reused); without it, one violin per metric.
    """
    metrics = _resolve_metrics(adata, metrics)
    cols = ([group_by] if group_by else []) + metrics
    df = _densify(adata.ap.to_df(obs=cols))
    id_vars = [group_by] if group_by else []
    long = df.melt(id_vars=id_vars, value_vars=metrics, var_name="metric", value_name="value")

    if group_by:
        p = (
            ggplot(long, aes(group_by, "value", fill=group_by))
            + geom_violin(scale=scale)
            + scale_fill_obs(adata, group_by)
        )
    else:
        long["_x"] = "all cells"
        p = ggplot(long, aes("_x", "value")) + geom_violin(scale=scale, fill="#4c72b0")

    p = (
        p
        + pe.facet_wrap("~metric", scales="free_y")
        + labs(x="", y="", fill=group_by or "")
        + theme_ggann()
    )
    if group_by:
        # the fill legend already encodes group; many long cell-type names across
        # shared-x facets are illegible -- drop the redundant x ticks (added last so
        # it wins over theme_ggann's axis_text).
        p = p + theme(axis_text_x=element_blank(), axis_ticks_major_x=element_blank())
    return p


def plot_qc_scatter(adata, x: str, y: str, color: str | None = None, *, size: float = 1.0):
    """Scatter of two obs QC metrics (e.g. total_counts vs pct_counts_mt).

    Thin wrapper over the grammar path; ``color`` may be an obs column or a gene.
    """
    from plotnine import geom_point

    from ._palette import scale_color_obs
    from .grammar import aes as _aes
    from .grammar import gganndata

    mapping = _aes(x, y) if color is None else _aes(x, y, color=color)
    plot = gganndata(adata, mapping) + geom_point(size=size, alpha=0.6)
    if color is not None:
        is_categorical_obs = color in adata.obs and isinstance(
            adata.obs[color].dtype, pd.CategoricalDtype
        )
        if is_categorical_obs:
            plot = plot + scale_color_obs(adata, color)
        elif color in adata.obs or color in adata.var_names or (
            adata.raw is not None and color in adata.raw.var_names
        ):
            # numeric (gene or continuous metric) -> expression colourmap, matching plot_embedding
            from .theme import scale_color_expression

            plot = plot + scale_color_expression()
    return plot


def plot_highest_expr_genes(adata, n: int = 20, *, use_raw: bool = False, layer: str | None = None):
    """Boxplots of the ``n`` genes accounting for the most counts per cell.

    Like ``sc.pl.highest_expr_genes``, this reads ``adata.X`` by default (pass
    ``layer="counts"`` or ``use_raw=True`` to point elsewhere): each cell's values
    are turned into percentages of that cell's total, and genes are ranked by mean
    per-cell percentage. For meaningful results ``adata.X`` should hold counts or
    normalized (not scaled) expression. The matrix is pulled through annplyr.
    """
    kind, lyr = expression_source(adata, layer, use_raw)
    if kind == "raw":
        wide = _raw_wide(adata, list(adata.raw.var_names))
    else:
        wide = _densify(adata.ap.to_df(x=list(adata.var_names), layer=lyr))

    if (wide.to_numpy() < 0).any():
        warnings.warn(
            "plot_highest_expr_genes: the expression matrix has negative values, which "
            "looks like scaled/z-scored data -- 'percent of total counts' will be "
            "meaningless. Pass use_raw=True or layer= to point at counts or "
            "log-normalized values.",
            stacklevel=2,
        )
    n_zero = int((wide.sum(axis=1) == 0).sum())
    if n_zero:
        warnings.warn(
            f"{n_zero} cell(s) with zero total counts excluded from plot_highest_expr_genes.",
            stacklevel=2,
        )
    totals = wide.sum(axis=1).replace(0, np.nan)
    frac = wide.div(totals, axis=0) * 100.0
    top = frac.mean(axis=0).sort_values(ascending=False).head(n).index.tolist()

    long = frac[top].melt(var_name="gene", value_name="percent")
    long["gene"] = pd.Categorical(long["gene"], categories=list(reversed(top)), ordered=True)
    return (
        ggplot(long, aes("gene", "percent"))
        + geom_boxplot(fill="#4c72b0", outlier_alpha=0.2)
        + pe.coord_flip()
        + labs(x="", y="% of total counts per cell")
        + theme_ggann()
    )
