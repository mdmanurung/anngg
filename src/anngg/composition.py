"""Cell-type composition / proportion bar plots.

Counts cells per group (optionally within a sample/condition) through annplyr's
``count`` verb, then renders stacked / filled / dodged bars -- the standard
"what fraction of each sample is each cell type" figure that scanpy's basics
don't cover.
"""

from __future__ import annotations

import plotnine_extra as pe
from plotnine import aes, geom_col, ggplot, labs

from ._palette import scale_fill_obs
from .plots import _group_categories, _order_groups
from .theme import theme_anngg

__all__ = ["plot_proportions"]


def plot_proportions(
    adata,
    group_by: str,
    split_by: str | None = None,
    *,
    normalize: bool = True,
    position: str = "stack",
    categories_order=None,
    split_order=None,
):
    """Bar plot of cell-type composition.

    Parameters
    ----------
    group_by
        The categorical whose composition is shown (the bar segments / fill).
    split_by
        Optional obs column for the x-axis (e.g. sample or condition); when
        omitted, one bar per ``group_by`` category.
    normalize
        Plot within-``split_by`` proportions (summing to 1) instead of raw counts.
    position
        ``"stack"`` (default), ``"fill"`` (100 %), or ``"dodge"``.
    """
    by = [split_by, group_by] if split_by else [group_by]
    counts = adata.ap.count(by=by)

    if normalize:
        if split_by:
            denom = counts.groupby(split_by, observed=True)["n"].transform("sum")
        else:
            denom = counts["n"].sum()
        counts["value"] = counts["n"] / denom
        ylab = "proportion of cells"
    else:
        counts["value"] = counts["n"]
        ylab = "number of cells"

    counts = _order_groups(
        counts, group_by, categories_order or _group_categories(adata, group_by)
    )
    if split_by:
        counts = _order_groups(
            counts, split_by, split_order or _group_categories(adata, split_by)
        )

    xvar = split_by if split_by else group_by
    # thin white borders separate the stacked segments cleanly (scplotter-style)
    col_kw = {"color": "white", "size": 0.15} if position != "dodge" else {}
    return (
        ggplot(counts, aes(xvar, "value", fill=group_by))
        + geom_col(position=position, width=0.9, **col_kw)
        + scale_fill_obs(adata, group_by)
        + labs(x="", y=ylab, fill=group_by)
        + theme_anngg()
        + pe.rotate_x_text(45)
    )
