"""Render one example figure per public plotting function for the API docs.

Writes ``docs/images/api/anngg.<func>.png`` for each helper; the Sphinx extension
``docs/extensions/api_examples.py`` then injects the matching image into that
function's API-reference page. Run locally (all optional deps installed) and
commit the PNGs -- the docs build itself never executes plotting code.

Run: ``python examples/api_examples.py``.
"""

from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import plotnine as p9
import scanpy as sc
from plotnine import geom_point

import anngg as ag
from anngg import aes, gganndata

GROUP = "bulk_labels"
MARKERS = ["CD3D", "CD8A", "NKG7", "GNLY", "MS4A1", "FCGR3A", "CST3"]
GENE_GROUPS = {"T": ["CD3D", "CD8A"], "NK": ["NKG7", "GNLY"], "B": ["MS4A1"], "Mono": ["CST3"]}


def _adata():
    adata = sc.datasets.pbmc68k_reduced()
    sc.tl.rank_genes_groups(adata, GROUP, method="wilcoxon", n_genes=50)
    return adata


def _examples(adata):
    de = ag.rank_genes_df(adata)
    sel = ["CD14+ Monocyte", "CD19+ B", "CD56+ NK", "Dendritic"]
    marker_sets = {g: list(de[de["group"] == g].head(20)["names"]) for g in sel}
    return {
        "anngg.gganndata": lambda: gganndata(adata, aes("UMAP_1", "UMAP_2", color=GROUP))
        + geom_point(size=1.2)
        + ag.theme_anngg(),
        "anngg.plot_embedding": lambda: ag.plot_embedding(adata, "umap", color=GROUP, label=True),
        "anngg.plot_features": lambda: ag.plot_features(adata, MARKERS[:4], basis="umap"),
        "anngg.plot_density": lambda: ag.plot_density(adata, ["CD3D", "NKG7"], joint=True),
        "anngg.plot_dotplot": lambda: ag.plot_dotplot(adata, MARKERS, GROUP),
        "anngg.plot_dotplot_grouped": lambda: ag.plot_dotplot_grouped(adata, GENE_GROUPS, GROUP),
        "anngg.plot_matrixplot": lambda: ag.plot_matrixplot(adata, MARKERS, GROUP, standard_scale="var"),
        "anngg.plot_matrixplot_grouped": lambda: ag.plot_matrixplot_grouped(adata, GENE_GROUPS, GROUP),
        "anngg.plot_violin": lambda: ag.plot_violin(adata, MARKERS[:3], GROUP),
        "anngg.plot_stacked_violin": lambda: ag.plot_stacked_violin(adata, MARKERS, GROUP),
        "anngg.plot_tracksplot": lambda: ag.plot_tracksplot(adata, MARKERS, GROUP),
        "anngg.plot_box": lambda: ag.plot_box(adata, MARKERS[:3], GROUP),
        "anngg.plot_expression_bar": lambda: ag.plot_expression_bar(adata, MARKERS[:3], GROUP),
        "anngg.plot_expression_line": lambda: ag.plot_expression_line(adata, ["CD3D"], x="phase", group_by=GROUP),
        "anngg.plot_proportions": lambda: ag.plot_proportions(adata, GROUP, split_by="phase"),
        "anngg.plot_correlation": lambda: ag.plot_correlation(adata, GROUP),
        "anngg.plot_rank_genes_dotplot": lambda: ag.plot_rank_genes_dotplot(adata, n_genes=3),
        "anngg.plot_rank_genes_matrixplot": lambda: ag.plot_rank_genes_matrixplot(adata, n_genes=3),
        "anngg.plot_volcano": lambda: ag.plot_volcano(adata, group="CD56+ NK"),
        "anngg.plot_qc_violin": lambda: ag.plot_qc_violin(
            adata, metrics=["n_genes", "percent_mito", "n_counts"], group_by=GROUP
        ),
        "anngg.plot_qc_scatter": lambda: ag.plot_qc_scatter(adata, x="n_counts", y="n_genes", color=GROUP),
        "anngg.plot_highest_expr_genes": lambda: ag.plot_highest_expr_genes(adata, n=20),
        "anngg.plot_clustermap": lambda: ag.plot_clustermap(adata, MARKERS, group_by=GROUP),
        "anngg.plot_upset": lambda: ag.plot_upset(marker_sets, min_cardinality=1),
    }


def _save(obj, path):
    """ggplot -> .save; marsilea/PyComplexHeatmap escape hatches -> current figure."""
    if isinstance(obj, p9.ggplot):
        obj.save(path, width=5.5, height=4, dpi=80, verbose=False)
    elif hasattr(obj, "save"):  # marsilea Upset
        obj.save(path, dpi=80)
    else:  # PyComplexHeatmap ClusterMapPlotter rendered onto the current figure
        plt.savefig(path, dpi=80, bbox_inches="tight")
    plt.close("all")


def main():
    adata = _adata()
    out = os.path.join(os.path.dirname(__file__), "..", "docs", "images", "api")
    os.makedirs(out, exist_ok=True)
    for name, builder in _examples(adata).items():
        _save(builder(), os.path.join(out, f"{name}.png"))
        print("wrote", name)


if __name__ == "__main__":
    main()
