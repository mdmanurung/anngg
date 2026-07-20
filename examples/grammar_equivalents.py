"""The layer-by-layer plotnine grammar behind each anngg utility helper.

Every plotnine-native ``plot_*`` helper is a thin wrapper: it resolves names
through annplyr into a tidy ``DataFrame``, then stacks plotnine layers. Each
function here rebuilds a helper's figure from that grammar, so you can see there
is no magic -- and drop down to raw grammar whenever a helper does not expose
what you need.

``gganndata(adata, aes(...))`` does the name resolution and returns a real
``plotnine.ggplot`` whose ``.data`` is the tidy frame the layers draw from; from
there it is ordinary plotnine (``+ geom_* + scale_* + theme``).

Escape hatches are intentionally omitted: ``plot_clustermap`` (PyComplexHeatmap)
and ``plot_upset`` (marsilea) are not plotnine and have no grammar equivalent.

Run: ``python examples/grammar_equivalents.py`` -- builds every equivalent and
saves them to ``docs/images/_grammar/``.
"""

from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")

import numpy as np
import plotnine_extra as pe
from plotnine import (
    aes,
    geom_boxplot,
    geom_col,
    geom_errorbar,
    geom_jitter,
    geom_line,
    geom_point,
    geom_tile,
    ggplot,
    guide_legend,
    guides,
    labs,
    scale_color_cmap,
    scale_fill_cmap,
)

import anngg as ag
from anngg import gganndata

GROUP = "bulk_labels"


# --- plot_embedding(color=..., label=True) --------------------------------- #
def embedding_grammar(adata, color=GROUP, basis="umap", label=True):
    coords = ag.embedding_coords(adata, basis)
    xcol, ycol = coords.columns[:2]
    base = gganndata(adata, ag.aes(xcol, ycol, color=color))  # resolves the obs column
    plot = (
        base
        + geom_point(size=1.5, alpha=0.9)
        + ag.scale_color_obs(adata, color)
        + guides(color=guide_legend(override_aes={"size": 4}))
        + ag.theme_anngg()
    )
    if label:
        cents = base.data.groupby(color, observed=True)[[xcol, ycol]].median().reset_index()
        cents = cents.rename(columns={color: "label"})
        plot = plot + pe.geom_label_repel(
            aes(xcol, ycol, label="label"), data=cents, fill="white", color="black",
            inherit_aes=False,
        )
    return plot


# --- plot_density(feature) ------------------------------------------------- #
def density_grammar(adata, feature="CD3D", basis="umap"):
    from pynebulosa import calculate_density

    coords = ag.embedding_coords(adata, basis)
    xcol, ycol = coords.columns[:2]
    # resolve expression through the grammar, then compute the weighted KDE
    resolved = gganndata(adata, ag.aes(xcol, ycol, color=ag.gene(feature))).data
    density = calculate_density(
        resolved[feature].to_numpy(float), resolved[[xcol, ycol]].to_numpy(float)
    )
    df = resolved.assign(density=density).sort_values("density")
    return (
        ggplot(df, aes(xcol, ycol, color="density"))
        + geom_point(size=1.5, alpha=0.9)
        + scale_color_cmap(cmap_name="magma")
        + labs(color="density")
        + ag.theme_anngg()
    )


# --- plot_box(gene, group) ------------------------------------------------- #
def box_grammar(adata, gene="CD3D", group=GROUP):
    # one gene shown; for a facet-per-gene grid melt the genes and add facet_wrap
    return (
        gganndata(adata, ag.aes(group, ag.gene(gene), fill=group))
        + geom_boxplot(width=0.7, outlier_alpha=0.0)
        + geom_jitter(width=0.2, height=0.0, size=0.35, alpha=0.25, stroke=0)
        + ag.scale_fill_obs(adata, group)
        + labs(x="", y="expression")
        + ag.theme_anngg()
        + pe.rotate_x_text(45)
    )


# --- plot_expression_bar(gene, group) -------------------------------------- #
def bar_grammar(adata, gene="CD3D", group=GROUP):
    d = gganndata(adata, ag.aes(group, ag.gene(gene))).data
    s = d.groupby(group, observed=True)[gene].agg(mean="mean", sd="std", n="count").reset_index()
    s["se"] = s["sd"] / np.sqrt(s["n"])
    s["ymin"], s["ymax"] = s["mean"] - s["se"], s["mean"] + s["se"]
    return (
        ggplot(s, aes(group, "mean", fill=group))
        + geom_col(width=0.7)
        + geom_errorbar(aes(ymin="ymin", ymax="ymax"), width=0.3)
        + ag.scale_fill_obs(adata, group)
        + labs(x="", y="mean expression")
        + ag.theme_anngg()
        + pe.rotate_x_text(45)
    )


# --- plot_expression_line(gene, x, group) ---------------------------------- #
def line_grammar(adata, gene="CD3D", x="phase", group=GROUP):
    d = gganndata(adata, ag.aes(x, ag.gene(gene), color=group)).data
    s = d.groupby([x, group], observed=True)[gene].mean().reset_index(name="mean")
    return (
        ggplot(s, aes(x, "mean", color=group, group=group))
        + geom_line()
        + geom_point(size=2)
        + ag.scale_color_obs(adata, group)
        + labs(x=x, y="mean expression")
        + ag.theme_anngg()
    )


# --- plot_correlation(group) ----------------------------------------------- #
def correlation_grammar(adata, group=GROUP, n_genes=100):
    import annplyr as ap

    genes = list(adata.var_names[np.asarray(adata.var["highly_variable"], bool)])[:n_genes]
    means = adata.ap.summarize(x={g: ap.mean(ap.col(g)) for g in genes}, by=group)
    profiles = means.set_index(group)[genes].apply(lambda c: np.asarray(c, float)).T  # genes x groups
    corr = profiles.corr()  # groups x groups
    long = corr.rename_axis("row").reset_index().melt(id_vars="row", var_name="col", value_name="corr")
    return (
        ggplot(long, aes("col", "row", fill="corr"))
        + geom_tile()
        + scale_fill_cmap(cmap_name="RdBu_r")
        + labs(x="", y="", fill="pearson\ncorrelation")
        + ag.theme_anngg()
        + pe.rotate_x_text(45)
    )


# Helper name -> the grammar function that reproduces it.
EQUIVALENTS = {
    "plot_embedding": embedding_grammar,
    "plot_density": density_grammar,
    "plot_box": box_grammar,
    "plot_expression_bar": bar_grammar,
    "plot_expression_line": line_grammar,
    "plot_correlation": correlation_grammar,
}


def main():
    import scanpy as sc

    adata = sc.datasets.pbmc68k_reduced()
    out = os.path.join(os.path.dirname(__file__), "..", "docs", "images", "_grammar")
    os.makedirs(out, exist_ok=True)
    for name, builder in EQUIVALENTS.items():
        plot = builder(adata)
        path = os.path.join(out, f"{name}.png")
        plot.save(path, width=6, height=5, dpi=90, verbose=False)
        print("wrote", os.path.relpath(path, os.path.join(os.path.dirname(__file__), "..")))


if __name__ == "__main__":
    main()
