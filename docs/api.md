# API Reference

```{eval-rst}
.. currentmodule:: anngg
```

## Grammar

The grammar entry point and the aesthetic-resolution escapes.

```{eval-rst}
.. autosummary::
   :toctree: generated
   :nosignatures:

   gganndata
   aes
   gene
   obs
   obsm
   embedding_coords
```

## Embeddings

```{eval-rst}
.. autosummary::
   :toctree: generated
   :nosignatures:

   plot_embedding
   plot_features
   plot_density
```

## Markers & expression summaries

```{eval-rst}
.. autosummary::
   :toctree: generated
   :nosignatures:

   plot_dotplot
   plot_dotplot_grouped
   plot_matrixplot
   plot_matrixplot_grouped
   plot_violin
   plot_ridge
   plot_stacked_violin
   plot_tracksplot
   plot_clustermap
```

## Distributions

```{eval-rst}
.. autosummary::
   :toctree: generated
   :nosignatures:

   plot_box
   plot_expression_bar
   plot_expression_line
```

## Differential expression

```{eval-rst}
.. autosummary::
   :toctree: generated
   :nosignatures:

   rank_genes_df
   plot_rank_genes_dotplot
   plot_rank_genes_matrixplot
   plot_rank_genes_heatmap
   plot_volcano
```

## Composition, correlation & sets

```{eval-rst}
.. autosummary::
   :toctree: generated
   :nosignatures:

   plot_proportions
   plot_correlation
   plot_upset
```

## Pseudobulk & stats

```{eval-rst}
.. autosummary::
   :toctree: generated
   :nosignatures:

   pseudobulk
   stat_compare_means
   stat_pwc
   stat_pvalue_manual
   stat_cor
   stat_regline_equation
   stat_anova_test
   stat_kruskal_test
   stat_central_tendency
   geom_signif
```

## Quality control

```{eval-rst}
.. autosummary::
   :toctree: generated
   :nosignatures:

   plot_qc_violin
   plot_qc_scatter
   plot_highest_expr_genes
```

## Scales, theme & layout

```{eval-rst}
.. autosummary::
   :toctree: generated
   :nosignatures:

   theme_anngg
   set_theme
   reset_theme
   sizes
   scale_color_obs
   scale_fill_obs
   obs_colors
   scale_color_expression
   scale_fill_expression
   scale_color_celltype
   scale_fill_celltype
   geom_text_repel
   geom_label_repel
   compose
   tag_panels
   Wrap
   Stack
   Beside
   plot_layout
   plot_annotation
   ggsave
```
