# anngg

**A ggplot2-style plotting layer for scanpy `AnnData` objects.**

`anngg` gives single-cell users the grammar of graphics over an `AnnData`, the
way `ggplot2` works in R. `gganndata(adata) + aes(...) + geom_*()` returns a real
[`plotnine.ggplot`](https://plotnine.org), and high-level `plot_*` helpers
reproduce scanpy's core figures — every one is a thin stack of grammar layers.

```python
import scanpy as sc
import anngg as ag
from anngg import gganndata, aes
from plotnine import geom_point

adata = sc.datasets.pbmc68k_reduced()

# grammar of graphics, straight over the AnnData
gganndata(adata, aes("UMAP_1", "UMAP_2", color="louvain")) + geom_point()

# or the high-level, scanpy-equivalent helpers
ag.plot_embedding(adata, basis="umap", color="CD3D")
ag.plot_dotplot(adata, ["CD3D", "NKG7", "CST3"], group_by="bulk_labels")
```

## Start here

::::{grid} 1 1 2 2
:gutter: 2

:::{grid-item-card} Install anngg
:link: installation
:link-type: doc

Set up the package and its optional extras.
:::

:::{grid-item-card} Quickstart
:link: quickstart
:link-type: doc

Build your first figures, both ways.
:::

:::{grid-item-card} Gallery
:link: gallery
:link-type: doc

See every figure `anngg` produces on `pbmc68k_reduced`.
:::

:::{grid-item-card} API Reference
:link: api
:link-type: doc

Look up the grammar, helpers, scales and theme.
:::

::::

```{toctree}
:caption: Get Started
:maxdepth: 1
:hidden:

installation
quickstart
gallery
comparisons
```

```{toctree}
:caption: API Reference
:maxdepth: 1
:hidden:

api
```
