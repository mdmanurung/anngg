"""High-level plotnine-native plotting helpers mirroring scanpy's core figures.

Each helper follows the same contract: extract a tidy DataFrame with annplyr,
then hand it to plotnine / plotnine-extra. None of them index ``adata`` directly.
All return ordinary :class:`plotnine.ggplot` objects, so they remain fully
composable with ``+ scale_*`` / ``+ theme(...)`` / ``+ facet_*``.
"""

from __future__ import annotations

from typing import Iterable, Sequence

import pandas as pd
import plotnine_extra as pe
from plotnine import (
    aes,
    geom_point,
    geom_tile,
    geom_violin,
    ggplot,
    labs,
    scale_color_cmap,
    scale_fill_cmap,
    scale_size,
)

from ._aggregate import aggregate_expression, tidy_expression
from ._palette import scale_color_obs, scale_fill_obs
from ._resolve import embedding_coords, plain_name, resolve_frame
from .theme import theme_anngg

__all__ = [
    "plot_embedding",
    "plot_features",
    "plot_dotplot",
    "plot_matrixplot",
    "plot_violin",
]


def _is_numeric(series: pd.Series) -> bool:
    return pd.api.types.is_numeric_dtype(series) and not isinstance(
        series.dtype, pd.CategoricalDtype
    )


def _centroid_labels(df: pd.DataFrame, cname: str, xcol: str, ycol: str) -> pd.DataFrame:
    """Median position of each category, for placing a cluster label at its centre."""
    cents = df.groupby(cname, observed=True)[[xcol, ycol]].median().reset_index()
    return cents.rename(columns={cname: "label"})


def plot_embedding(
    adata,
    basis: str = "umap",
    color: str | None = None,
    *,
    split_by: str | None = None,
    layer: str | None = None,
    use_raw: bool | None = None,
    size: float = 1.5,
    alpha: float = 0.9,
    pointdensity: bool | None = None,
    label: bool = False,
    label_size: float = 9,
    low: str = "#d9d9d9",
    high: str = "#2166ac",
):
    """Scatter over an embedding (UMAP/t-SNE/PCA), optionally coloured and split.

    * ``color=None`` -> density-coloured scatter (``geom_pointdensity``), which
      reads well for dense embeddings.
    * categorical ``color`` (e.g. an obs cluster column) -> discrete colours
      (reusing scanpy's stored palette when present).
    * numeric ``color`` (a gene or continuous obs column) -> a gradient.
    * ``split_by`` -> facet the scatter over an obs column (Seurat ``split.by``).
    * ``label=True`` -> for a categorical ``color``, print each category at its
      centroid using repelled (non-overlapping) text, like scplotter's
      ``CellDimPlot`` / Seurat ``label=TRUE``.
    """
    coords = embedding_coords(adata, basis)
    if coords.shape[1] < 2:
        raise ValueError(
            f"Embedding '{basis}' has only {coords.shape[1]} dimension(s); "
            "plot_embedding requires at least 2."
        )
    xcol, ycol = coords.columns[:2]

    split_col = None
    if split_by is not None:
        split_col = resolve_frame(adata, [split_by])[[split_by]]

    def _facet(plot):
        return plot + pe.facet_wrap("~" + split_by) if split_by is not None else plot

    if color is None:
        if pointdensity is None:
            pointdensity = True
        base = coords if split_col is None else coords.join(split_col)
        if pointdensity:
            return _facet(
                ggplot(base, aes(xcol, ycol))
                + pe.geom_pointdensity(size=size, alpha=alpha)
                + labs(color="density")
                + theme_anngg()
            )
        return _facet(pe.DimPlot(base, x=xcol, y=ycol, size=size, alpha=alpha) + theme_anngg())

    # `color` may be a bare name, a prefix string ("gene:CD3D@logcounts") or an accessor
    cname = plain_name(adata, color)
    values = resolve_frame(adata, [color], layer=layer, use_raw=use_raw)
    if cname not in values.columns:
        raise KeyError(f"Could not resolve color={color!r} from obs, genes or obsm.")
    df = coords.join(values[[cname]])
    if split_col is not None and split_by not in df.columns:
        df = df.join(split_col)

    if _is_numeric(df[cname]):
        # Draw low-expression cells first so high-expression cells are not occluded
        # (mirrors scanpy's sc.pl.embedding ordering).
        df = df.sort_values(cname)
        return _facet(
            pe.FeatureDimPlot(
                df, feature=cname, x=xcol, y=ycol, low=low, high=high, size=size, alpha=alpha
            )
            + theme_anngg()
        )
    plot = (
        pe.DimPlot(df, x=xcol, y=ycol, color=cname, size=size, alpha=alpha)
        + scale_color_obs(adata, cname)
        + theme_anngg()
    )
    if label:
        cents = _centroid_labels(df, cname, xcol, ycol)
        plot = plot + pe.geom_text_repel(
            aes(xcol, ycol, label="label"),
            data=cents,
            size=label_size,
            fontweight="bold",
            color="black",
            inherit_aes=False,
        )
    return _facet(plot)


def plot_features(
    adata,
    features: Sequence[str],
    basis: str = "umap",
    *,
    layer: str | None = None,
    use_raw: bool | None = None,
    ncol: int | None = None,
    size: float = 1.2,
    alpha: float = 0.9,
    cmap: str = "magma",
):
    """Multi-gene embedding grid: one faceted panel per feature.

    Panels share a single expression colour scale, which is best for comparing
    magnitudes *across* genes. Because the scale is shared, a low-range gene shown
    next to a high-range one will look faint -- for independent per-gene colour
    bars (like ``sc.pl.umap(color=[...])``), compose separate ``plot_embedding``
    calls with the re-exported ``Wrap`` / ``plot_layout`` instead.
    """
    coords = embedding_coords(adata, basis)
    if coords.shape[1] < 2:
        raise ValueError(f"Embedding '{basis}' has fewer than 2 dimensions.")
    xcol, ycol = coords.columns[:2]

    values = resolve_frame(adata, list(features), layer=layer, use_raw=use_raw)
    feats = list(
        dict.fromkeys(f for f in features if f in values.columns and _is_numeric(values[f]))
    )
    if not feats:
        raise ValueError("plot_features needs at least one numeric feature (gene or metric).")

    df = coords.join(values[feats])
    long = df.melt(
        id_vars=[xcol, ycol], value_vars=feats, var_name="feature", value_name="expression"
    ).sort_values("expression")
    long["feature"] = pd.Categorical(long["feature"], categories=feats, ordered=True)
    return (
        ggplot(long, aes(xcol, ycol, color="expression"))
        + geom_point(size=size, alpha=alpha)
        + pe.facet_wrap("~feature", ncol=ncol)
        + scale_color_cmap(cmap_name=cmap)
        + theme_anngg()
    )


def _group_categories(adata, group_by: str):
    """The obs categorical order for ``group_by`` (e.g. Leiden 0,1,..,10,11), else None.

    annplyr's ``summarize`` / ``to_tidy`` strip the Categorical dtype, so the
    meaningful order must be captured from ``adata.obs`` before aggregation and
    threaded down to :func:`_order_groups`.
    """
    col = adata.obs[group_by]
    if isinstance(col.dtype, pd.CategoricalDtype):
        return list(col.cat.categories)
    return None


def _order_groups(frame: pd.DataFrame, group_by: str, categories_order):
    """Order the group axis, defaulting to the obs categorical order (not alphabetical).

    Groups are nominal (cell identities), so the categorical is left
    ``ordered=False`` -- this keeps discrete-axis placement but lets plotnine pick
    a *qualitative* colour palette when ``group_by`` drives a fill/colour
    aesthetic (an ordered categorical would trigger a misleading sequential ramp).
    """
    if categories_order is not None:
        cats = list(dict.fromkeys(categories_order))  # dedupe, keep order
    else:
        col = frame[group_by]
        if isinstance(col.dtype, pd.CategoricalDtype):
            cats = list(col.cat.categories)
        else:
            cats = sorted(pd.unique(col))
    present = set(pd.unique(pd.Series(frame[group_by]).astype(object)))
    missing = present - set(cats)
    if missing:
        raise ValueError(
            f"categories_order is missing groups present in the data: {sorted(missing)}"
        )
    frame[group_by] = pd.Categorical(frame[group_by], categories=cats, ordered=False)
    frame[group_by] = frame[group_by].cat.remove_unused_categories()
    return frame


def plot_dotplot(
    adata,
    genes: Sequence[str],
    group_by: str,
    *,
    layer: str | None = None,
    use_raw: bool | None = None,
    standard_scale: str | None = None,
    expression_cutoff: float = 0.0,
    cmap: str = "Reds",
    size_range: tuple[float, float] = (0.5, 8.0),
    categories_order: Iterable[str] | None = None,
):
    """Marker dotplot: dot *size* = fraction expressing, *colour* = mean expression.

    Defaults to ``adata.raw`` so values match ``sc.pl.dotplot``.
    """
    agg = aggregate_expression(
        adata,
        genes,
        group_by,
        layer=layer,
        use_raw=use_raw,
        expression_cutoff=expression_cutoff,
        standard_scale=standard_scale,
    )
    if categories_order is None:
        categories_order = _group_categories(adata, group_by)
    agg = _order_groups(agg, group_by, categories_order)
    color_lab = "scaled\nexpression" if standard_scale else "mean\nexpression"
    return (
        ggplot(agg, aes("feature", group_by))
        + geom_point(aes(size="fraction", color="mean_expression"))
        + scale_color_cmap(cmap_name=cmap)
        + scale_size(range=size_range, labels=lambda xs: [f"{x:.0%}" for x in xs])
        + labs(x="", y="", color=color_lab, size="fraction\nexpressing")
        + theme_anngg()
        + pe.rotate_x_text(45)
    )


def plot_matrixplot(
    adata,
    genes: Sequence[str],
    group_by: str,
    *,
    layer: str | None = None,
    use_raw: bool | None = None,
    standard_scale: str | None = None,
    cmap: str = "viridis",
    categories_order: Iterable[str] | None = None,
):
    """Aggregated mean-expression heatmap (genes x groups) as a plotnine tile plot.

    Like ``sc.pl.matrixplot``, this defaults to ``standard_scale=None`` (raw group
    means). Pass ``standard_scale="var"`` to rescale each gene to ``[0, 1]``, which
    keeps cross-gene patterns legible when magnitudes differ widely.
    """
    agg = aggregate_expression(
        adata,
        genes,
        group_by,
        layer=layer,
        use_raw=use_raw,
        standard_scale=standard_scale,
    )
    if categories_order is None:
        categories_order = _group_categories(adata, group_by)
    agg = _order_groups(agg, group_by, categories_order)
    color_lab = "scaled\nexpression" if standard_scale else "mean\nexpression"
    return (
        ggplot(agg, aes("feature", group_by, fill="mean_expression"))
        + geom_tile()
        + scale_fill_cmap(cmap_name=cmap)
        + labs(x="", y="", fill=color_lab)
        + theme_anngg()
        + pe.rotate_x_text(45)
    )


def plot_violin(
    adata,
    genes: Sequence[str],
    group_by: str,
    *,
    layer: str | None = None,
    use_raw: bool | None = None,
    ncol: int = 1,
    scale: str = "width",
    stats: bool = False,
    categories_order: Iterable[str] | None = None,
):
    """Per-group expression distributions, one facet per gene (stacked-violin style).

    Set ``stats=True`` to overlay a group-comparison test via plotnine-extra's
    ``stat_compare_means``.
    """
    genes = list(genes)
    tidy = tidy_expression(adata, genes, group_by, layer=layer, use_raw=use_raw)
    if categories_order is None:
        categories_order = _group_categories(adata, group_by)
    tidy = _order_groups(tidy, group_by, categories_order)
    plot = (
        ggplot(tidy, aes(group_by, "value", fill=group_by))
        + geom_violin(scale=scale)
        + pe.facet_wrap("~feature", ncol=ncol, scales="free_y")
        + scale_fill_obs(adata, group_by)
        + labs(x="", y="expression", fill=group_by)
        + theme_anngg()
        + pe.rotate_x_text(45)
    )
    if stats:
        plot = plot + pe.stat_compare_means()
    return plot
