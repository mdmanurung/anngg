"""Per-group expression distribution plots: box, bar and line.

Inspired by DOTools / scplotter's ``FeatureStatPlot`` family. Aesthetics follow
scplotter: box plots overlay the underlying cells as jittered points, bar plots
show the group mean with a summary error bar, and the line plot tracks mean
expression across an ordered variable. All extraction goes through annplyr
(``tidy_expression`` / ``resolve_frame``) and every helper returns a plain
:class:`plotnine.ggplot`, so ``+ scale_* / + theme(...) / + facet_*`` still compose.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np
import pandas as pd
import plotnine_extra as pe
from plotnine import (
    aes,
    geom_boxplot,
    geom_col,
    geom_errorbar,
    geom_jitter,
    geom_line,
    geom_point,
    ggplot,
    labs,
)

from ._aggregate import tidy_expression
from ._palette import scale_color_obs, scale_fill_obs
from ._resolve import plain_name, resolve_frame
from .plots import _feature_facet, _group_categories, _is_numeric, _order_groups
from .theme import theme_anngg

__all__ = ["plot_box", "plot_expression_bar", "plot_expression_line"]


def _summarise(values: pd.core.groupby.SeriesGroupBy, error: str, agg="mean") -> pd.DataFrame:
    """Central tendency (``agg``) and a symmetric error (se / sd / none) per group.

    ``agg`` is any pandas reduction name or callable (``"mean"``, ``"median"``,
    ``np.sum`` ...); the result column is still called ``mean`` for the callers.
    """
    summary = values.agg(mean=agg, sd="std", n="count").reset_index()
    summary["sd"] = summary["sd"].fillna(0.0)
    if error == "se":
        summary["err"] = summary["sd"] / np.sqrt(summary["n"].clip(lower=1))
    elif error == "sd":
        summary["err"] = summary["sd"]
    elif error == "none":
        summary["err"] = 0.0
    else:
        raise ValueError(f"error must be 'se', 'sd' or 'none', got {error!r}")
    summary["ymin"] = summary["mean"] - summary["err"]
    summary["ymax"] = summary["mean"] + summary["err"]
    return summary


def plot_box(
    adata,
    genes: Sequence[str],
    group_by: str,
    *,
    split_by: str | None = None,
    layer: str | None = None,
    use_raw: bool | None = None,
    ncol: int = 1,
    jitter: bool = True,
    jitter_size: float = 0.35,
    jitter_alpha: float = 0.25,
    stats: bool = False,
    categories_order: Sequence[str] | None = None,
):
    """Per-group expression box plots, one facet per gene, with jittered cells overlaid.

    Set ``jitter=False`` for a plain box plot, ``split_by`` for a gene x split facet
    grid, or ``stats=True`` to overlay a group-comparison test via plotnine-extra's
    ``stat_compare_means``.
    """
    genes = list(genes)
    extra = [split_by] if split_by else []
    tidy = tidy_expression(adata, genes, group_by, layer=layer, use_raw=use_raw, extra_obs=extra)
    if categories_order is None:
        categories_order = _group_categories(adata, group_by)
    tidy = _order_groups(tidy, group_by, categories_order)

    plot = ggplot(tidy, aes(group_by, "value", fill=group_by)) + geom_boxplot(
        width=0.7, outlier_alpha=0.0 if jitter else 1.0
    )
    if jitter:
        plot = plot + geom_jitter(width=0.2, height=0.0, size=jitter_size, alpha=jitter_alpha, stroke=0)
    plot = (
        plot
        + _feature_facet(split_by, ncol=ncol)
        + scale_fill_obs(adata, group_by)
        + labs(x="", y="expression", fill=group_by)
        + theme_anngg()
        + pe.rotate_x_text(45)
    )
    if stats:
        plot = plot + pe.stat_compare_means()
    return plot


def plot_expression_bar(
    adata,
    genes: Sequence[str],
    group_by: str,
    *,
    split_by: str | None = None,
    layer: str | None = None,
    use_raw: bool | None = None,
    ncol: int = 1,
    agg="mean",
    error: str = "se",
    categories_order: Sequence[str] | None = None,
):
    """Aggregated expression per group as bars with an error bar, one facet per gene.

    ``agg`` is the bar height (``"mean"`` default, or ``"median"`` / ``np.sum`` / any
    pandas reduction). ``error`` is the summary bar: ``"se"`` (default), ``"sd"`` or
    ``"none"``. ``split_by`` adds a gene x split facet grid. Bars start at zero and
    hide the distribution -- for a distribution-honest view use :func:`plot_box` or
    :func:`anngg.plot_violin`.
    """
    genes = list(genes)
    extra = [split_by] if split_by else []
    tidy = tidy_expression(adata, genes, group_by, layer=layer, use_raw=use_raw, extra_obs=extra)
    if categories_order is None:
        categories_order = _group_categories(adata, group_by)
    tidy = _order_groups(tidy, group_by, categories_order)

    by = [group_by, *extra, "feature"]
    summary = _summarise(tidy.groupby(by, observed=True)["value"], error, agg=agg)
    summary["feature"] = pd.Categorical(summary["feature"], categories=genes, ordered=True)

    ylab = f"{agg} expression" if isinstance(agg, str) else "expression"
    plot = (
        ggplot(summary, aes(group_by, "mean", fill=group_by))
        + geom_col(width=0.7)
        + _feature_facet(split_by, ncol=ncol)
        + scale_fill_obs(adata, group_by)
        + labs(x="", y=ylab, fill=group_by)
        + theme_anngg()
        + pe.rotate_x_text(45)
    )
    if error != "none":
        plot = plot + geom_errorbar(aes(ymin="ymin", ymax="ymax"), width=0.3)
    return plot


def plot_expression_line(
    adata,
    genes: Sequence[str],
    x: str,
    *,
    group_by: str | None = None,
    layer: str | None = None,
    use_raw: bool | None = None,
    ncol: int = 1,
    agg="mean",
    error: str = "se",
    categories_order: Sequence[str] | None = None,
):
    """Aggregated expression trend across an ordered variable ``x``, one facet per gene.

    ``x`` is an obs column (an ordered condition, timepoint or a numeric covariate
    such as pseudotime). With ``group_by`` given, one line is drawn per group and
    coloured by it (e.g. an expression trajectory per cell type across timepoints).
    ``agg`` sets the summarised value (``"mean"`` default, or ``"median"`` etc.);
    ``error`` adds a summary error bar per point (``"se"`` / ``"sd"`` / ``"none"``).
    """
    genes = list(genes)
    xname = plain_name(adata, x)
    gname = plain_name(adata, group_by) if group_by is not None else None

    gene_names = [plain_name(adata, g) for g in genes]
    id_vars = [xname] + ([gname] if gname is not None else [])
    clash = set(gene_names) & set(id_vars)
    if clash:
        raise ValueError(
            f"feature name(s) {sorted(clash)} collide with x/group_by; melt would drop them. "
            f"Rename, or disambiguate with gene()/obs()."
        )

    cols = [x] + ([group_by] if group_by is not None else []) + list(genes)
    frame = resolve_frame(adata, cols, layer=layer, use_raw=use_raw)  # already densified
    long = frame.melt(id_vars=id_vars, value_vars=gene_names, var_name="feature", value_name="value")
    long["feature"] = pd.Categorical(long["feature"], categories=gene_names, ordered=True)

    # Order the x axis: numeric stays numeric; categorical keeps its obs order.
    if not _is_numeric(long[xname]):
        x_cats = _group_categories(adata, x)
        if x_cats is None:
            x_cats = list(pd.unique(long[xname]))
        long[xname] = pd.Categorical(long[xname], categories=x_cats, ordered=True)
    if gname is not None and categories_order is None:
        categories_order = _group_categories(adata, group_by)
    if gname is not None:
        long = _order_groups(long, gname, categories_order)

    summary = _summarise(long.groupby(id_vars + ["feature"], observed=True)["value"], error, agg=agg)
    summary["feature"] = pd.Categorical(summary["feature"], categories=gene_names, ordered=True)

    color = gname if gname is not None else None
    mapping = aes(x=xname, y="mean", color=color, group=color) if color else aes(x=xname, y="mean", group=1)
    ylab = f"{agg} expression" if isinstance(agg, str) else "expression"
    plot = (
        ggplot(summary, mapping)
        + geom_line()
        + geom_point(size=2)
        + pe.facet_wrap("~feature", ncol=ncol, scales="free_y")
        + labs(x=xname, y=ylab, color=gname)
        + theme_anngg()
    )
    if color is not None:
        plot = plot + scale_color_obs(adata, gname)
    if error != "none":
        plot = plot + geom_errorbar(aes(ymin="ymin", ymax="ymax"), width=0.15)
    return plot
