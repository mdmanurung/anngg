"""Grouped expression aggregation, expressed through annplyr ``summarize``.

Produces the two quantities scanpy's dotplot / matrixplot are built from:

* ``mean_expression`` -- mean of the (raw / X / layer) matrix per group per gene
* ``fraction``        -- fraction of cells with expression above a cutoff

Matching scanpy, aggregation defaults to ``adata.raw`` when present. With
``standard_scale=None`` (the dotplot default) the mean is the raw group mean, so
the numbers reproduce ``sc.pl.dotplot`` / ``sc.pl.matrixplot`` (which also
default to ``standard_scale=None``).
"""

from __future__ import annotations

from typing import Iterable, Literal

import annplyr as ap
import pandas as pd

from ._resolve import _densify, resolve_source

__all__ = ["aggregate_expression", "expression_source", "group_means", "tidy_expression"]

StandardScale = Literal["var", "group", "zscore"]


def tidy_expression(adata, genes, group_by, *, layer=None, use_raw=None):
    """Long per-cell expression ``[obs_name, feature, value, group_by]`` via annplyr.

    Shared by the violin / stacked-violin / tracksplot paths so the source
    picking, densification and feature-ordering live in one place.
    """
    genes = list(genes)
    kind, lyr = resolve_source(adata, layer, use_raw)
    if kind == "raw":
        tidy = adata.ap.to_tidy(obs=[group_by], raw=genes)
    else:
        tidy = adata.ap.to_tidy(obs=[group_by], x=genes, layer=lyr)
    tidy = _densify(tidy)
    tidy["feature"] = pd.Categorical(tidy["feature"], categories=genes, ordered=True)
    return tidy


def expression_source(adata, layer: str | None, use_raw: bool | None) -> tuple[str, str | None]:
    """Decide which matrix to aggregate: ``("raw", None)``, ``("layer", name)`` or ``("x", None)``.

    Thin wrapper over :func:`anngg._resolve.resolve_source` so the aggregation
    path and the grammar path share one source-of-truth (and raise identically
    on ``use_raw=True`` with no ``adata.raw``, or ``layer`` + ``use_raw=True``).
    """
    return resolve_source(adata, layer, use_raw)


def group_means(
    adata, genes: Iterable[str], group_by: str, *, layer=None, use_raw=None
) -> pd.DataFrame:
    """Mean expression per group (index) per gene (columns), mean-only.

    Like :func:`aggregate_expression` but skips the fraction-expressing pass, for
    callers (e.g. :func:`anngg.plot_correlation`) that only need the group means.
    """
    genes = list(genes)
    kind, lyr = expression_source(adata, layer, use_raw)
    mean_expr = {g: ap.mean(ap.col(g)) for g in genes}
    if kind == "raw":
        mean_df = adata.ap.summarize(raw=mean_expr, by=group_by)
    else:
        mean_df = adata.ap.summarize(x=mean_expr, by=group_by, layer=lyr)
    return _densify(mean_df).set_index(group_by)[genes].astype(float)


def _standardize(mean_df: pd.DataFrame, standard_scale: str | None) -> pd.DataFrame:
    """Optionally rescale mean expression across groups/genes.

    ``'var'`` (per-gene 0..1) and ``'group'`` (per-group 0..1) match scanpy's
    ``standard_scale``. ``'zscore'`` (per-gene z-score, population std) is an
    anngg extension not present in scanpy.
    """
    if standard_scale is None:
        return mean_df
    if standard_scale == "var":  # per-gene (column) 0..1
        rng = (mean_df.max() - mean_df.min()).replace(0, 1)
        return (mean_df - mean_df.min()) / rng
    if standard_scale == "group":  # per-group (row) 0..1
        rng = (mean_df.max(axis=1) - mean_df.min(axis=1)).replace(0, 1)
        return mean_df.sub(mean_df.min(axis=1), axis=0).div(rng, axis=0)
    if standard_scale == "zscore":  # per-gene z-score
        std = mean_df.std(ddof=0).replace(0, 1)
        return (mean_df - mean_df.mean()) / std
    raise ValueError(
        f"standard_scale must be None, 'var', 'group' or 'zscore', got {standard_scale!r}"
    )


def aggregate_expression(
    adata,
    genes: Iterable[str],
    group_by: str,
    *,
    layer: str | None = None,
    use_raw: bool | None = None,
    expression_cutoff: float = 0.0,
    standard_scale: StandardScale | None = None,
) -> pd.DataFrame:
    """Return a long DataFrame ``[group_by, feature, mean_expression, fraction]``.

    ``feature`` is an ordered categorical following the order of ``genes`` so
    downstream plots keep the requested gene order. ``standard_scale`` may be
    ``'var'`` / ``'group'`` (scanpy conventions) or ``'zscore'`` (an anngg
    extension); ``None`` leaves the raw group means untouched.
    """
    genes = list(genes)
    kind, lyr = expression_source(adata, layer, use_raw)

    mean_expr = {g: ap.mean(ap.col(g)) for g in genes}
    frac_expr = {g: ap.mean(ap.col(g) > expression_cutoff) for g in genes}

    if kind == "raw":
        mean_df = adata.ap.summarize(raw=mean_expr, by=group_by)
        frac_df = adata.ap.summarize(raw=frac_expr, by=group_by)
    else:
        mean_df = adata.ap.summarize(x=mean_expr, by=group_by, layer=lyr)
        frac_df = adata.ap.summarize(x=frac_expr, by=group_by, layer=lyr)

    mean_df = _densify(mean_df).set_index(group_by)[genes].astype(float)
    frac_df = _densify(frac_df).set_index(group_by)[genes].astype(float)
    mean_df = _standardize(mean_df, standard_scale)

    long = (
        mean_df.reset_index()
        .melt(id_vars=group_by, var_name="feature", value_name="mean_expression")
        .merge(
            frac_df.reset_index().melt(
                id_vars=group_by, var_name="feature", value_name="fraction"
            ),
            on=[group_by, "feature"],
        )
    )
    long["feature"] = pd.Categorical(long["feature"], categories=genes, ordered=True)
    return long
