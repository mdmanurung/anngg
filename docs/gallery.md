# Gallery

Every figure below is produced by
[`examples/gallery.py`](https://github.com/mdmanurung/anngg/blob/main/examples/gallery.py)
on `pbmc68k_reduced`. The clean, boxed look emulates
[scplotter](https://pwwang.github.io/scplotter/) / plotthis.

## Embeddings

:::{grid} 1 1 2 2

![UMAP by cluster](images/umap_clusters.png)

![Labelled UMAP](images/umap_labelled.png)
:::

`plot_embedding` — categorical clusters (stored palette) and, with
`label=True`, repelled centroid labels like scplotter's `CellDimPlot`.

![Multi-gene grid](images/features_grid.png)

`plot_features` — one shared-scale panel per gene.

## Gene-weighted density (pyNebulosa)

![Density](images/density.png)

`plot_density` recovers marker signal lost to dropout with a weighted KDE;
`joint=True` adds a co-expression panel.

## Markers & expression

:::{grid} 1 1 2 2

![Dotplot](images/dotplot.png)

![Stacked violin](images/stacked_violin.png)

![Box](images/box.png)

![Expression bar](images/expression_bar.png)
:::

`plot_dotplot`, `plot_stacked_violin`, `plot_box`, `plot_expression_bar`.

## Differential expression

:::{grid} 1 1 2 2

![DE dotplot](images/de_dotplot.png)

![Volcano](images/volcano.png)
:::

`plot_rank_genes_dotplot`, `plot_volcano` (with repelled gene labels).

## Composition, correlation & sets

:::{grid} 1 1 2 2

![Composition](images/proportions.png)

![Correlation](images/correlation.png)

![UpSet](images/upset.png)

![Expression line](images/expression_line.png)
:::

`plot_proportions`, `plot_correlation`, `plot_upset` (marsilea),
`plot_expression_line`.

## Quality control

![QC violins](images/qc_violin.png)

`plot_qc_violin`.
