# Quickstart

Every figure in `ggann` is a real `plotnine.ggplot`. You can build it two ways —
with the low-level grammar or with a high-level helper — and they compose with
plotnine's `+ geom_* / + scale_* / + theme(...)` either way.

```python
import scanpy as sc
import ggann as ag
from ggann import gganndata, aes
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
from ggann import gene
gganndata(adata, aes("UMAP_1", "UMAP_2", color=gene("CD3D"))) + geom_point()
```

## The helpers

The `plot_*` helpers reproduce scanpy's core figures:

```python
ag.plot_embedding(adata, "umap", color="bulk_labels", label=True)
ag.plot_dotplot(adata, ["CD3D", "NKG7", "CST3"], group_by="bulk_labels")
ag.plot_density(adata, ["CD3D", "NKG7"], joint=True)      # needs ggann[density]
ag.plot_box(adata, ["CD3D", "NKG7"], group_by="bulk_labels")
```

Each plotnine-native helper is only a stack of grammar layers — see
[`examples/grammar_equivalents.py`](https://github.com/mdmanurung/ggann/blob/main/examples/grammar_equivalents.py)
for the layer-by-layer twin of every helper.

## One consistent look

Call `ag.set_theme()` once to make the ggann theme plotnine's default, so every
figure — even a bare `ggplot(...)` — shares it without `+ theme_ggann()`:

```python
ag.set_theme(base_size=9, family="Arial")   # family is optional; no font is required
gganndata(adata, aes("UMAP_1", "UMAP_2", color="louvain")) + geom_point()  # already themed
```

`ag.sizes` gives a matching font-size scale (`.normal`/`.small`/`.large`/`.title`,
in pt) so text in annotations stays in sync — convert to `geom_text`'s mm unit with
`ag.sizes.geom(...)`. `ag.reset_theme()` restores plotnine's default.

## Composing panels

`ag.compose` assembles a tagged multi-panel figure from a list of plots (each
keeps its own scales), which you then save at an exact physical size:

```python
fig = ag.compose(
    [ag.plot_embedding(adata, "umap", color="bulk_labels"),
     ag.plot_dotplot(adata, markers, "bulk_labels"),
     ag.plot_violin(adata, markers[:1], "bulk_labels"),
     ag.plot_proportions(adata, "bulk_labels", split_by="phase")],
    ncol=2,                 # A/B/C/D panel tags by default
)
fig.save("figure1.pdf", width=180, height=140, units="mm")   # millimetre-exact output
```

`tag_levels=` switches the labels (`"a"`, `"1"`, `"i"`, or `None`), and
`ag.tag_panels(plots)` tags a list you compose by hand with plotnine's `|` / `/`.
