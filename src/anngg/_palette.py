"""Reuse scanpy's stored category colours so anngg matches the rest of an analysis.

scanpy stores a per-category colour list in ``adata.uns[f"{col}_colors"]``, aligned
to ``adata.obs[col].cat.categories``. When present, anngg picks it up so a cluster
plotted here has the same colours as ``sc.pl.umap`` of the same object.
"""

from __future__ import annotations

import pandas as pd
import plotnine_extra as pe
from plotnine import scale_color_manual, scale_fill_manual

__all__ = ["obs_colors", "scale_color_obs", "scale_colour_obs", "scale_fill_obs"]


def obs_colors(adata, col: str) -> dict | None:
    """Return ``{category: hex}`` from ``adata.uns[f"{col}_colors"]``, or ``None``.

    Falls back to ``None`` when the column is not a categorical or has no stored
    colours (or the colour list is shorter than the categories).
    """
    dtype = adata.obs[col].dtype if col in adata.obs else None
    if not isinstance(dtype, pd.CategoricalDtype):
        return None
    colors = adata.uns.get(f"{col}_colors")
    cats = list(adata.obs[col].cat.categories)
    if colors is None or len(colors) < len(cats):
        return None
    return {cat: str(c) for cat, c in zip(cats, colors)}


def scale_color_obs(adata, col: str, **kwargs):
    """Categorical colour scale using scanpy's stored colours, else a qualitative default."""
    mapping = obs_colors(adata, col)
    if mapping is None:
        return pe.scale_color_hue(**kwargs)
    return scale_color_manual(values=mapping, **kwargs)


# British-spelling alias, mirroring plotnine.
scale_colour_obs = scale_color_obs


def scale_fill_obs(adata, col: str, **kwargs):
    """Categorical fill scale using scanpy's stored colours, else a qualitative default."""
    mapping = obs_colors(adata, col)
    if mapping is None:
        return pe.scale_fill_hue(**kwargs)
    return scale_fill_manual(values=mapping, **kwargs)
