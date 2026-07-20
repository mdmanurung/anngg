"""Clustered-heatmap escape hatch, delegating to PyComplexHeatmap.

Clustered heatmaps with dendrograms and stacked annotation bars are a grid-based
paradigm that does not fit the grammar of graphics, so they live *outside* the
``gganndata() + geom_*`` path on purpose. This one function is the bridge: it
builds an annplyr-tidied matrix and hands it to
:class:`PyComplexHeatmap.ClusterMapPlotter`, returning that plotter object.

``PyComplexHeatmap`` is an optional dependency; it is imported lazily so that
``import ggann`` never requires it.
"""

from __future__ import annotations

from typing import Iterable, Sequence

import pandas as pd

from ._aggregate import aggregate_expression, expression_source
from ._resolve import _densify

__all__ = ["plot_clustermap"]


def _require_pch():
    try:
        import PyComplexHeatmap as pch  # noqa: N813
    except ImportError as exc:  # pragma: no cover - env-dependent
        raise ImportError(
            "plot_clustermap requires PyComplexHeatmap. "
            "Install it with `pip install ggann[heatmap]` or `pip install PyComplexHeatmap`."
        ) from exc
    return pch


def plot_clustermap(
    adata,
    genes: Sequence[str],
    group_by: str | None = None,
    *,
    layer: str | None = None,
    use_raw: bool | None = None,
    annotations: Iterable[str] | None = None,
    standard_scale: str | None = None,
    z_score: int | None = None,
    row_cluster: bool = True,
    col_cluster: bool = True,
    cmap: str = "viridis",
    show_rownames: bool = True,
    show_colnames: bool = True,
    plot: bool = True,
    **kwargs,
):
    """Clustered heatmap of ``genes`` (rows) across groups or cells (columns).

    * ``group_by`` given -> columns are aggregated group means (like
      ``sc.pl.heatmap`` on pseudobulk). Column annotations default to the group
      identity.
    * ``group_by=None`` -> columns are individual cells; ``annotations`` selects
      ``obs`` columns to show as top annotation bars.

    ``standard_scale`` (ggann-side, per gene/group) and ``z_score``
    (PyComplexHeatmap-side row/column z-score) are mutually exclusive -- applying
    both would normalize twice. Returns the
    :class:`PyComplexHeatmap.ClusterMapPlotter` instance.
    """
    if standard_scale is not None and z_score is not None:
        raise ValueError(
            "standard_scale and z_score are mutually exclusive (both normalize the "
            "matrix); set at most one."
        )
    pch = _require_pch()
    genes = list(genes)

    if group_by is not None:
        agg = aggregate_expression(
            adata, genes, group_by, layer=layer, use_raw=use_raw,
            standard_scale=standard_scale,
        )
        matrix = agg.pivot(index="feature", columns=group_by, values="mean_expression")
        matrix = matrix.loc[genes]
        ann_df = pd.DataFrame({group_by: matrix.columns}, index=matrix.columns)
    else:
        kind, lyr = expression_source(adata, layer, use_raw)
        if kind == "raw":
            wide = adata.ap.to_df(raw=genes)
            wide = wide.rename(columns={f"raw_{g}": g for g in genes})
        else:
            wide = adata.ap.to_df(x=genes, layer=lyr)
        wide = _densify(wide)[genes]
        matrix = wide.T  # genes x cells
        ann_cols = list(annotations) if annotations else []
        ann_df = adata.ap.to_df(obs=ann_cols) if ann_cols else None

    top_annotation = None
    if ann_df is not None and not ann_df.empty:
        # pandas>=3 gives string columns the `str` dtype, which PyComplexHeatmap
        # does not recognise as categorical and then refuses to auto-pick a cmap.
        # Coerce non-numeric annotations to `category` so it colours them discretely.
        # `is_numeric_dtype` is True for bool, so coerce bool columns too -- they
        # are categorical (e.g. is_doublet), not a continuous scale.
        ann_df = ann_df.copy()
        for col in ann_df.columns:
            is_numeric = pd.api.types.is_numeric_dtype(ann_df[col])
            if not is_numeric or pd.api.types.is_bool_dtype(ann_df[col]):
                ann_df[col] = ann_df[col].astype("category")
        top_annotation = pch.HeatmapAnnotation(df=ann_df, axis=1, plot_legend=True)

    return pch.ClusterMapPlotter(
        data=matrix,
        top_annotation=top_annotation,
        z_score=z_score,
        row_cluster=row_cluster,
        col_cluster=col_cluster,
        cmap=cmap,
        show_rownames=show_rownames,
        show_colnames=show_colnames,
        plot=plot,
        **kwargs,
    )
