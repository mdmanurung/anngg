"""Plots from differential-expression / marker results (``rank_genes_groups``).

scanpy stores marker results in ``adata.uns['rank_genes_groups']`` and exposes a
tidy view via ``sc.get.rank_genes_groups_df``. These helpers reuse that (never
parsing the recarrays by hand), pick the top markers, and delegate to the
existing ``plot_dotplot`` / ``plot_matrixplot`` -- or, for the volcano, to
plotnine-extra's ``ggvolcano``.
"""

from __future__ import annotations

import plotnine_extra as pe
from plotnine import labs, scale_color_manual

from .plots import plot_dotplot, plot_matrixplot
from .theme import theme_anngg

# Conventional volcano colours: down = blue, non-significant = grey, up = red.
# Keys must match the categories ``plotnine_extra.ggvolcano`` maps to its colour
# aesthetic; unused keys are ignored, so a one-sided result (only "up") is fine.
_VOLCANO_COLORS = {"down": "#3B4CC0", "not significant": "#B8B8B8", "up": "#B40426"}

__all__ = [
    "rank_genes_df",
    "plot_rank_genes_dotplot",
    "plot_rank_genes_matrixplot",
    "plot_rank_genes_heatmap",
    "plot_volcano",
]


def _require_de(adata, key: str):
    if key not in adata.uns:
        raise KeyError(
            f"adata.uns[{key!r}] not found; run sc.tl.rank_genes_groups(adata, ...) first."
        )


def _de_groupby(adata, key: str, groupby):
    if groupby is not None:
        return groupby
    return adata.uns[key]["params"]["groupby"]


def rank_genes_df(adata, group=None, key: str = "rank_genes_groups", **kwargs):
    """Tidy DE table for ``group`` (or all groups), via ``sc.get.rank_genes_groups_df``.

    Columns: ``group, names, scores, logfoldchanges, pvals, pvals_adj``. Extra
    kwargs (``pval_cutoff``, ``log2fc_min``, ...) pass through to scanpy.
    """
    import scanpy as sc

    _require_de(adata, key)
    return sc.get.rank_genes_groups_df(adata, group=group, key=key, **kwargs)


def _top_genes(adata, n_genes: int, key: str) -> list[str]:
    df = rank_genes_df(adata, group=None, key=key)
    top = df.groupby("group", observed=True).head(n_genes)
    return list(dict.fromkeys(top["names"]))  # unique, group-ordered


def plot_rank_genes_dotplot(
    adata, n_genes: int = 5, groupby: str | None = None,
    key: str = "rank_genes_groups", **kwargs,
):
    """Dotplot of the top ``n_genes`` markers per group (``sc.pl.rank_genes_groups_dotplot``)."""
    _require_de(adata, key)
    groupby = _de_groupby(adata, key, groupby)
    genes = _top_genes(adata, n_genes, key)
    return plot_dotplot(adata, genes, groupby, **kwargs)


def plot_rank_genes_matrixplot(
    adata, n_genes: int = 5, groupby: str | None = None,
    key: str = "rank_genes_groups", standard_scale: str | None = "var", **kwargs,
):
    """Matrixplot (group-mean tiles) of the top ``n_genes`` markers per group.

    Matches ``sc.pl.rank_genes_groups_matrixplot`` (group summary), not the
    per-cell ``..._heatmap``.
    """
    _require_de(adata, key)
    groupby = _de_groupby(adata, key, groupby)
    genes = _top_genes(adata, n_genes, key)
    return plot_matrixplot(adata, genes, groupby, standard_scale=standard_scale, **kwargs)


# Back-compat alias (this renders group-mean tiles, i.e. matrixplot semantics).
plot_rank_genes_heatmap = plot_rank_genes_matrixplot


def plot_volcano(
    adata, group: str, key: str = "rank_genes_groups",
    lfc: float = 1.0, padj: float = 0.05, label_top: int = 10, **kwargs,
):
    """Volcano plot (log2FC vs adjusted p-value) for one group's markers.

    Reuses plotnine-extra's ``ggvolcano``; ``lfc``/``padj`` set the fold-change /
    significance cutoffs and ``label_top`` labels the strongest genes. Requires a
    ``rank_genes_groups`` computed with a method that reports p-values and
    log-fold-changes (``wilcoxon`` / ``t-test``, not ``logreg``).
    """
    _require_de(adata, key)
    df = rank_genes_df(adata, group=group, key=key)
    missing = {"logfoldchanges", "pvals_adj"} - set(df.columns)
    if missing:
        raise ValueError(
            f"rank_genes_groups result is missing {sorted(missing)} needed for a volcano; "
            "re-run sc.tl.rank_genes_groups(adata, ..., method='wilcoxon' or 't-test')."
        )
    return (
        pe.ggvolcano(
            df, x="logfoldchanges", y="pvals_adj", label="names",
            p_cutoff=padj, fc_cutoff=lfc, label_top=label_top, **kwargs,
        )
        + scale_color_manual(values=_VOLCANO_COLORS)
        + labs(x="log2 fold change", y="-log10(adjusted p-value)")
        + theme_anngg()
    )
