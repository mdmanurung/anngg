"""Default theme and single-cell colour scales for anngg.

These are thin conveniences over plotnine / plotnine-extra so that anngg figures
have a consistent, publication-friendly look while remaining ordinary plotnine
objects the user can override with any ``+ theme(...)`` / ``+ scale_*``.
"""

from __future__ import annotations

import plotnine_extra as pe
from plotnine import element_blank, element_line, element_rect, element_text, theme, theme_bw

__all__ = [
    "theme_anngg",
    "scale_color_expression",
    "scale_colour_expression",
    "scale_fill_expression",
    "scale_color_celltype",
    "scale_colour_celltype",
    "scale_fill_celltype",
]


def theme_anngg(base_size: float = 11, base_family: str | None = None):
    """A clean, boxed theme emulating scplotter / plotthis ``theme_this``.

    The signature scplotter look is a thin black panel border (a full box) with
    no axis lines and no grid -- not the L-shaped axes of ``theme_classic``.
    """
    return theme_bw(base_size=base_size, base_family=base_family) + theme(
        panel_grid=element_blank(),
        panel_border=element_rect(color="black", fill=None, size=1),
        axis_line=element_blank(),
        axis_ticks=element_line(color="black"),
        axis_text=element_text(color="black"),
        legend_key=element_blank(),
        strip_background=element_blank(),
        strip_text=element_text(color="black"),
    )


def scale_color_expression(cmap: str = "Reds", **kwargs):
    """Continuous colour scale for expression (perceptually ordered)."""
    return pe.scale_color_cmap(cmap_name=cmap, **kwargs)


# British-spelling alias, mirroring plotnine.
scale_colour_expression = scale_color_expression


def scale_fill_expression(cmap: str = "viridis", **kwargs):
    """Continuous fill scale for expression heatmaps."""
    return pe.scale_fill_cmap(cmap_name=cmap, **kwargs)


def scale_color_celltype(**kwargs):
    """Categorical colour scale for discrete cell identities."""
    return pe.scale_color_hue(**kwargs)


# British-spelling alias, mirroring plotnine.
scale_colour_celltype = scale_color_celltype


def scale_fill_celltype(**kwargs):
    """Categorical fill scale for discrete cell identities."""
    return pe.scale_fill_hue(**kwargs)
