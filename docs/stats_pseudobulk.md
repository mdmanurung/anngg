# Stats & pseudobulk

## ggpubr-style statistical layers

Because ggann plots are ordinary plotnine objects, the full
[ggpubr](https://rpkgs.datanovia.com/ggpubr/) stat-layer family (shipped by
`plotnine-extra`) composes onto any of them with `+`. ggann re-exports the common
ones so you can reach them straight from `ggann`:

| ggann | ggpubr | what it adds |
|---|---|---|
| `stat_compare_means` | `stat_compare_means` | global or pairwise group-comparison p-values |
| `stat_pwc` | `stat_pwc` | pairwise comparisons with significance brackets |
| `stat_pvalue_manual` | `stat_pvalue_manual` | your own p-values as brackets |
| `stat_cor` | `stat_cor` | correlation coefficient + p on a scatter |
| `stat_regline_equation` | `stat_regline_equation` | regression line equation / R² |
| `stat_anova_test` / `stat_kruskal_test` | same | one-line omnibus tests |
| `stat_central_tendency` | `stat_central_tendency` | a mean / median marker |
| `geom_signif` | (ggsignif) | significance brackets |

```python
import ggann as ag

# correlation annotation on a QC scatter
ag.plot_qc_scatter(adata, x="n_counts", y="n_genes") + ag.stat_cor()

# pairwise significance brackets on a violin
ag.plot_violin(adata, ["CD3D"], "bulk_labels") + ag.stat_pwc()
```

`plot_violin` / `plot_box` also take `stats=True` as a shortcut for
`+ stat_compare_means()`.

## Pseudobulk — then reuse the same grammar

`ag.pseudobulk` aggregates cells to one profile per sample × group (per donor per
cell type) with [decoupler](https://decoupler-py.readthedocs.io/)'s
`pp.pseudobulk` — the same aggregation the
[liana MOFA-cellular](https://liana-py.readthedocs.io/en/latest/notebooks/mofacellular.html)
workflow builds on. It returns a **new AnnData**, so the entire ggann grammar and
every `plot_*` helper work on the pseudobulk object with **no changes** — you just
point them at `pb` instead of `adata`:

```python
import ggann as ag
from ggann import gganndata, aes, gene
from plotnine import geom_boxplot

# cells -> sample x cell-type pseudobulk profiles (needs ggann[pseudobulk])
pb = ag.pseudobulk(adata, sample_col="donor", group_by="cell_type",
                   layer="counts", mode="sum")

# the same grammar, now over pseudobulk samples
gganndata(pb, aes("cell_type", gene("CD3D"), fill="condition")) + geom_boxplot()

# and the same helpers
ag.plot_dotplot(pb, markers, "cell_type", use_raw=False)
ag.plot_correlation(pb, "cell_type")            # sample-profile correlation
ag.plot_violin(pb, markers, "cell_type", use_raw=False) + ag.stat_compare_means()
```

decoupler expects integer counts: pass the counts `layer=` (or `raw=True`), or
`skip_checks=True` to aggregate non-count values. Profiles built from fewer than
`min_cells` cells are dropped.

```{note}
Because pseudobulk is just an AnnData, this composes with everything else:
`plot_correlation` for sample-sample similarity, `stat_compare_means` for
condition tests, `Wrap`/`Stack` for multi-panel figures — the grammar does not
care whether a row is a cell or a pseudobulk sample.
```
