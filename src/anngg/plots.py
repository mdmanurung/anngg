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

from ._aggregate import aggregate_expression, expression_source
from ._resolve import _densify, embedding_coords, plain_name, resolve_frame
from .theme import theme_anngg

__all__ = ["plot_embedding", "plot_dotplot", "plot_matrixplot", "plot_violin"]


def _is_numeric(series: pd.Series) -> bool:
    return pd.api.types.is_numeric_dtype(series) and not isinstance(
        series.dtype, pd.CategoricalDtype
    )


def plot_embedding(
    adata,
    basis: str = "umap",
    color: str | None = None,
    *,
    layer: str | None = None,
    use_raw: bool | None = None,
    size: float = 1.5,
    alpha: float = 0.9,
    pointdensity: bool | None = None,
    low: str = "#d9d9d9",
    high: str = "#2166ac",
):
    """Scatter over an embedding (UMAP/t-SNE/PCA), optionally coloured.

    * ``color=None`` -> density-coloured scatter (``geom_pointdensity``), which
      reads well for dense embeddings.
    * categorical ``color`` (e.g. an obs cluster column) -> discrete colours.
    * numeric ``color`` (a gene or continuous obs column) -> a gradient.
    """
    coords = embedding_coords(adata, basis)
    if coords.shape[1] < 2:
        raise ValueError(
            f"Embedding '{basis}' has only {coords.shape[1]} dimension(s); "
            "plot_embedding requires at least 2."
        )
    xcol, ycol = coords.columns[:2]

    if color is None:
        if pointdensity is None:
            pointdensity = True
        if pointdensity:
            return (
                ggplot(coords, aes(xcol, ycol))
                + pe.geom_pointdensity(size=size, alpha=alpha)
                + labs(color="density")
                + theme_anngg()
            )
        return pe.DimPlot(coords, x=xcol, y=ycol, size=size, alpha=alpha) + theme_anngg()

    # `color` may be a bare name, a prefix string ("gene:CD3D@logcounts") or an accessor
    cname = plain_name(adata, color)
    values = resolve_frame(adata, [color], layer=layer, use_raw=use_raw)
    if cname not in values.columns:
        raise KeyError(f"Could not resolve color={color!r} from obs, genes or obsm.")
    df = coords.join(values[[cname]])

    if _is_numeric(df[cname]):
        # Draw low-expression cells first so high-expression cells are not occluded
        # (mirrors scanpy's sc.pl.embedding ordering).
        df = df.sort_values(cname)
        return (
            pe.FeatureDimPlot(
                df, feature=cname, x=xcol, y=ycol, low=low, high=high, size=size, alpha=alpha
            )
            + theme_anngg()
        )
    return pe.DimPlot(df, x=xcol, y=ycol, color=cname, size=size, alpha=alpha) + theme_anngg()


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
    categories_order: Iterable[str] | None = None,
):
    """Per-group expression distributions, one facet per gene (stacked-violin style)."""
    genes = list(genes)
    kind, lyr = expression_source(adata, layer, use_raw)
    if kind == "raw":
        tidy = adata.ap.to_tidy(obs=[group_by], raw=genes)
    else:
        tidy = adata.ap.to_tidy(obs=[group_by], x=genes, layer=lyr)
    tidy = _densify(tidy)
    tidy["feature"] = pd.Categorical(tidy["feature"], categories=genes, ordered=True)
    if categories_order is None:
        categories_order = _group_categories(adata, group_by)
    tidy = _order_groups(tidy, group_by, categories_order)
    return (
        ggplot(tidy, aes(group_by, "value", fill=group_by))
        + geom_violin(scale=scale)
        + pe.facet_wrap("~feature", ncol=ncol, scales="free_y")
        + labs(x="", y="expression", fill=group_by)
        + theme_anngg()
        + pe.rotate_x_text(45)
    )
