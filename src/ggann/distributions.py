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
    geom_violin,
    ggplot,
    labs,
)

from ._aggregate import tidy_expression
from ._palette import scale_color_obs, scale_fill_obs
from ._resolve import plain_name, resolve_frame
from .plots import (
    _downsample_cells,
    _feature_facet,
    _group_categories,
    _is_numeric,
    _order_groups,
)
from .theme import theme_ggann

__all__ = ["plot_box", "plot_sina", "plot_expression_bar", "plot_expression_line"]


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
    downsample: int | None = None,
):
    """Per-group expression box plots, one facet per gene, with jittered cells overlaid.

    Set ``jitter=False`` for a plain box plot, ``split_by`` for a gene x split facet
    grid, or ``stats=True`` to overlay a group-comparison test via plotnine-extra's
    ``stat_compare_means``. Pass ``downsample=N`` to cap cells per group before the
    (jitter) draw for large data — see :func:`ggann.plot_violin`.

    Note: ``downsample`` subsets the cells the geoms see, so with ``stats=True``
    the p-value is computed on the subsample, and the boxplot's outliers reflect
    only the kept cells. Leave it unset when either must reflect every cell.
    """
    adata = _downsample_cells(adata, group_by, downsample)
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
        + theme_ggann()
        + pe.rotate_x_text(45)
    )
    if stats:
        plot = plot + pe.stat_compare_means()
    return plot


def plot_sina(
    adata,
    genes: Sequence[str],
    group_by: str,
    *,
    split_by: str | None = None,
    layer: str | None = None,
    use_raw: bool | None = None,
    ncol: int = 1,
    size: float = 0.4,
    alpha: float = 0.5,
    violin: bool = True,
    bins: int = 50,
    categories_order: Sequence[str] | None = None,
    downsample: int | None = None,
):
    """Sina / beeswarm of per-group expression, one facet per gene.

    A sina plot spreads each group's cells horizontally in proportion to the local
    density (via plotnine-extra's ``geom_sina``), so it shows every cell like a
    jitter but keeps the shape of a violin. ``violin=True`` draws a faint violin
    behind the points for context. ``downsample=N`` caps cells per group first --
    recommended for large data, since a sina draws one mark per cell.
    """
    genes = list(genes)
    adata = _downsample_cells(adata, group_by, downsample)
    extra = [split_by] if split_by else []
    tidy = tidy_expression(adata, genes, group_by, layer=layer, use_raw=use_raw, extra_obs=extra)
    if categories_order is None:
        categories_order = _group_categories(adata, group_by)
    tidy = _order_groups(tidy, group_by, categories_order)

    # Drive stat_sina by binwidth rather than bins: plotnine 0.15's stat_sina
    # mutates its shared ``bins`` param across panels (the second facet then gets
    # an array where a scalar is expected and crashes). The binwidth branch reads
    # an untouched param each panel, so it is multi-facet safe.
    #
    # Size the binwidth off the NARROWEST gene's range (÷ bins) so every facet
    # gets ~``bins`` bins -- a single global span sized off the widest gene would
    # collapse a narrow-range gene's panel to one or two bins. Floor it by the
    # widest gene's range so a pathologically narrow gene can't drive the widest
    # panel to an unbounded bin count (which would be slow to build).
    ranges = tidy.groupby("feature", observed=True)["value"].agg(lambda s: s.max() - s.min())
    ranges = ranges[ranges > 0]
    if len(ranges):
        binwidth = max(float(ranges.min()) / bins, float(ranges.max()) / 1000)
    else:
        binwidth = 1.0

    plot = ggplot(tidy, aes(group_by, "value", color=group_by))
    if violin:
        plot = plot + geom_violin(
            aes(fill=group_by), color="none", alpha=0.15, scale="width", show_legend=False
        )
    plot = (
        plot
        + pe.geom_sina(size=size, alpha=alpha, stroke=0, binwidth=binwidth)
        + _feature_facet(split_by, ncol=ncol)
        + scale_color_obs(adata, group_by)
        + scale_fill_obs(adata, group_by)
        + labs(x="", y="expression", color=group_by)
        + theme_ggann()
        + pe.rotate_x_text(45)
    )
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
    :func:`ggann.plot_violin`.
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
        + theme_ggann()
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
        + theme_ggann()
    )
    if color is not None:
        plot = plot + scale_color_obs(adata, gname)
    if error != "none":
        plot = plot + geom_errorbar(aes(ymin="ymin", ymax="ymax"), width=0.15)
    return plot
