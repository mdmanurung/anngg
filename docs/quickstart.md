# Quickstart

Every figure in `anngg` is a real `plotnine.ggplot`. You can build it two ways —
with the low-level grammar or with a high-level helper — and they compose with
plotnine's `+ geom_* / + scale_* / + theme(...)` either way.

```python
import scanpy as sc
import anngg as ag
from anngg import gganndata, aes
from plotnine import geom_point

adata = sc.datasets.pbmc68k_reduced()
```

## The grammar

`gganndata(adata, aes(...))` resolves the names in your `aes()` — obs columns,
genes, embedding coordinates — into a tidy `DataFrame` and returns a `ggplot`:

```python
gganndata(adata, aes("UMAP_1", "UMAP_2", color="louvain")) + geom_point(size=1.5)
```

Force a source with the `gene(...)` / `obs(...)` escapes when a name is
ambiguous, or mix sources across aesthetics:

```python
from anngg import gene
gganndata(adata, aes("UMAP_1", "UMAP_2", color=gene("CD3D"))) + geom_point()
```

## The helpers

The `plot_*` helpers reproduce scanpy's core figures:

```python
ag.plot_embedding(adata, "umap", color="bulk_labels", label=True)
ag.plot_dotplot(adata, ["CD3D", "NKG7", "CST3"], group_by="bulk_labels")
ag.plot_density(adata, ["CD3D", "NKG7"], joint=True)      # needs anngg[density]
ag.plot_box(adata, ["CD3D", "NKG7"], group_by="bulk_labels")
```

Each plotnine-native helper is only a stack of grammar layers — see
[`examples/grammar_equivalents.py`](https://github.com/mdmanurung/anngg/blob/main/examples/grammar_equivalents.py)
for the layer-by-layer twin of every helper.

## Composing panels

`anngg` re-exports plotnine-extra's layout operators, so multi-panel figures
with independent scales are one expression:

```python
from anngg import Wrap
Wrap([ag.plot_embedding(adata, "umap", color=g) for g in ["CD3D", "NKG7"]])
```
