"""Ridgeline (joy) plots of per-group expression distributions.

plotnine has no ridgeline geom, so this builds one: a gaussian KDE per group on a
shared grid, offset vertically by group and drawn as overlapping filled areas
(``geom_ribbon``). One facet per gene; matches scplotter's ``plot_type = "ridge"``.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np
import pandas as pd
import plotnine_extra as pe
from plotnine import (
    aes,
    element_blank,
    geom_line,
    geom_ribbon,
    ggplot,
    labs,
    scale_y_continuous,
    theme,
)

from ._aggregate import tidy_expression
from ._palette import scale_fill_obs
from .plots import _group_categories, _order_groups
from .theme import theme_anngg

__all__ = ["plot_ridge"]


def _kde_curve(values: np.ndarray, grid: np.ndarray) -> np.ndarray:
    """KDE density on ``grid``; zeros for degenerate (constant / tiny) groups."""
    from scipy.stats import gaussian_kde

    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size < 2 or np.ptp(values) == 0:
        return np.zeros_like(grid)
    try:
        dens = gaussian_kde(values)(grid)
    except np.linalg.LinAlgError:  # singular covariance
        return np.zeros_like(grid)
    peak = dens.max()
    return dens / peak if peak > 0 else dens


def plot_ridge(
    adata,
    genes: Sequence[str],
    group_by: str,
    *,
    layer: str | None = None,
    use_raw: bool | None = None,
    ncol: int = 1,
    scale: float = 1.6,
    n_grid: int = 256,
    categories_order: Sequence[str] | None = None,
):
    """Ridgeline plot: one density ridge per group, stacked and overlapping, per gene.

    ``scale`` sets how tall each ridge is relative to its row spacing (>1 overlaps
    neighbours, the classic joyplot look). Groups with fewer than two cells or zero
    variance draw a flat baseline.
    """
    genes = list(genes)
    tidy = tidy_expression(adata, genes, group_by, layer=layer, use_raw=use_raw)
    if categories_order is None:
        categories_order = _group_categories(adata, group_by)
    tidy = _order_groups(tidy, group_by, categories_order)
    order = list(tidy[group_by].cat.categories)
    pos = {g: i for i, g in enumerate(order)}

    rows = []
    for feature, fdf in tidy.groupby("feature", observed=True):
        lo, hi = float(fdf["value"].min()), float(fdf["value"].max())
        if hi <= lo:
            hi = lo + 1.0
        grid = np.linspace(lo, hi, n_grid)
        for g, gdf in fdf.groupby(group_by, observed=True):
            dens = _kde_curve(gdf["value"].to_numpy(), grid) * scale
            base = float(pos[g])
            rows.append(
                pd.DataFrame(
                    {
                        "x": grid,
                        "ymin": base,
                        "ymax": base + dens,
                        group_by: g,
                        "feature": str(feature),
                    }
                )
            )
    long = pd.concat(rows, ignore_index=True)
    long[group_by] = pd.Categorical(long[group_by], categories=order, ordered=True)
    long["feature"] = pd.Categorical(long["feature"], categories=[str(x) for x in genes], ordered=True)

    # Draw top ridges first so lower ones layer IN FRONT and overlap cleanly
    # (a single ribbon layer paints higher rows on top, clipping the peaks below).
    # Each ridge is a filled area with a thin white line tracing only its top edge,
    # not a boxed outline -- so the silhouettes blend instead of looking cut.
    plot = ggplot(long, aes("x"))
    for g in reversed(order):
        gd = long[long[group_by] == g]
        plot = plot + geom_ribbon(
            aes(ymin="ymin", ymax="ymax", fill=group_by), data=gd, alpha=0.9
        )
        plot = plot + geom_line(aes(y="ymax"), data=gd, color="white", size=0.4)
    return (
        plot
        + pe.facet_wrap("~feature", ncol=ncol, scales="free_x")
        + scale_fill_obs(adata, group_by)
        + scale_y_continuous(breaks=list(pos.values()), labels=order)
        + labs(x="expression", y="", fill=group_by)
        + theme_anngg()
        + theme(panel_grid=element_blank())
    )
