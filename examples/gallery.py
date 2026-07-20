"""Generate the README gallery figures into ``docs/images/``.

Run: ``python examples/gallery.py`` (writes committed PNGs referenced by README).
Covers the full anngg surface on ``pbmc68k_reduced``.
"""

from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")

import scanpy as sc

import anngg as ag

HERE = os.path.dirname(__file__)
OUT = os.path.join(HERE, "..", "docs", "images")
os.makedirs(OUT, exist_ok=True)

adata = sc.datasets.pbmc68k_reduced()
adata.obs["condition"] = adata.obs["phase"].astype(str)
group = "bulk_labels"
markers = ["CD3D", "CD8A", "NKG7", "GNLY", "MS4A1", "FCGR3A", "CST3"]
sc.tl.rank_genes_groups(adata, group, method="wilcoxon", n_genes=50)


def save(plot, name, **kw):
    path = os.path.join(OUT, name)
    plot.save(path, verbose=False, **kw)
    print("wrote", os.path.relpath(path, os.path.join(HERE, "..")))


save(ag.plot_embedding(adata, "umap", color=group), "umap_clusters.png", width=7, height=5, dpi=100)
save(ag.plot_features(adata, markers[:4], basis="umap"), "features_grid.png", width=8, height=6, dpi=100)
save(ag.plot_dotplot(adata, markers, group), "dotplot.png", width=7, height=5, dpi=100)
save(ag.plot_rank_genes_dotplot(adata, n_genes=3), "de_dotplot.png", width=9, height=5, dpi=100)
save(ag.plot_volcano(adata, group="CD56+ NK"), "volcano.png", width=6.5, height=5, dpi=100)
save(ag.plot_proportions(adata, group, split_by="condition", position="fill"),
     "proportions.png", width=6, height=5, dpi=100)
save(ag.plot_stacked_violin(adata, markers, group), "stacked_violin.png", width=6, height=8, dpi=100)
save(ag.plot_qc_violin(adata, metrics=["n_genes", "percent_mito", "n_counts"], group_by=group),
     "qc_violin.png", width=7, height=8, dpi=90)

# pyNebulosa-style density + dotools/scplotter-inspired utilities + repel labels
save(ag.plot_embedding(adata, "umap", color=group, label=True),
     "umap_labelled.png", width=7, height=5, dpi=100)
save(ag.plot_density(adata, ["CD3D", "NKG7"], basis="umap", joint=True),
     "density.png", width=10, height=3.2, dpi=100)
save(ag.plot_box(adata, markers[:4], group), "box.png", width=6, height=8, dpi=100)
save(ag.plot_expression_bar(adata, markers[:4], group), "expression_bar.png", width=6, height=8, dpi=100)
save(ag.plot_expression_line(adata, markers[:3], x="phase", group_by=group),
     "expression_line.png", width=7, height=7, dpi=100)
save(ag.plot_correlation(adata, group, cluster=True), "correlation.png", width=6, height=5, dpi=100)

# marsilea UpSet: which top markers are shared across cell types
import matplotlib.pyplot as plt  # noqa: E402

de = ag.rank_genes_df(adata)
sel = ["CD14+ Monocyte", "CD19+ B", "CD56+ NK", "Dendritic"]
marker_sets = {g: list(de[de["group"] == g].head(20)["names"]) for g in sel}
upset = ag.plot_upset(marker_sets, min_cardinality=1)
upset.save(os.path.join(OUT, "upset.png"), dpi=100)
plt.close("all")
print("wrote", os.path.relpath(os.path.join(OUT, "upset.png"), os.path.join(HERE, "..")))

print("\nGallery written to", os.path.relpath(OUT, os.path.join(HERE, "..")))
