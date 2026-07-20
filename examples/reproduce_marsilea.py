"""Reproduce the scanpy + marsilea how-to figures using ggann's grammar.

Mirrors https://scanpy.readthedocs.io/en/stable/how-to/plotting-with-marsilea.html
but every figure is expressed with the ggplot2 grammar (or, for the clustered
heatmap, the PyComplexHeatmap escape hatch). Run:

    python examples/reproduce_marsilea.py

which writes PNGs to ``examples/_output/``.
"""

from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import scanpy as sc

import ggann as ag

OUT = os.path.join(os.path.dirname(__file__), "_output")
os.makedirs(OUT, exist_ok=True)

adata = sc.datasets.pbmc68k_reduced()
group = "bulk_labels"
markers = ["CD3D", "CD8A", "NKG7", "GNLY", "MS4A1", "FCGR3A", "CST3"]
markers = [g for g in markers if g in adata.raw.var_names]


def save(plot, name, **kw):
    path = os.path.join(OUT, name)
    plot.save(path, verbose=False, **kw)
    print("wrote", path)


# 1. Embedding: grammar of graphics over UMAP, coloured by cluster.
save(
    ag.plot_embedding(adata, basis="umap", color=group),
    "01_umap_clusters.png",
    width=7,
    height=5,
    dpi=120,
)

# 2. Embedding coloured by a marker gene (defaults to adata.raw).
save(
    ag.plot_embedding(adata, basis="umap", color="CD3D"),
    "02_umap_CD3D.png",
    width=6,
    height=5,
    dpi=120,
)

# 3. Marker dotplot (size = fraction expressing, colour = mean expression).
save(
    ag.plot_dotplot(adata, markers, group),
    "03_dotplot.png",
    width=7,
    height=5,
    dpi=120,
)

# 4. Matrixplot: scaled mean-expression heatmap.
save(
    ag.plot_matrixplot(adata, markers, group, standard_scale="var"),
    "04_matrixplot.png",
    width=7,
    height=5,
    dpi=120,
)

# 5. Stacked violins, one facet per gene.
save(
    ag.plot_violin(adata, markers, group),
    "05_violin.png",
    width=6,
    height=12,
    dpi=120,
)

# 6. Escape hatch: clustered heatmap with a cell-type annotation bar.
plt.figure(figsize=(7, 5))
ag.plot_clustermap(adata, markers, group_by=group, standard_scale="var", z_score=None)
plt.savefig(os.path.join(OUT, "06_clustermap.png"), dpi=120, bbox_inches="tight")
plt.close("all")
print("wrote", os.path.join(OUT, "06_clustermap.png"))

print("\nAll figures written to", OUT)
