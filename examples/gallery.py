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

print("\nGallery written to", os.path.relpath(OUT, os.path.join(HERE, "..")))
