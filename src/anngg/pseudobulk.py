"""Pseudobulk aggregation, so cell-level AnnData can be summarised to sample level.

Aggregates single cells to one profile per ``sample_col`` x ``group_by`` (e.g. per
donor per cell type) with decoupler's ``pp.pseudobulk`` -- the same aggregation
liana's MOFA-cellular workflow builds on. It returns a **new AnnData**, so every
anngg grammar / ``plot_*`` helper works on the pseudobulk object unchanged:

    pb = ag.pseudobulk(adata, sample_col="donor", group_by="cell_type")
    ag.plot_dotplot(pb, markers, "cell_type", use_raw=False)
    gganndata(pb, aes("cell_type", gene("CD3D"), fill="condition")) + geom_boxplot()

decoupler is an optional dependency (``anngg[pseudobulk]``).
"""

from __future__ import annotations

__all__ = ["pseudobulk"]


def _require_decoupler():
    try:
        import decoupler as dc
    except ImportError as exc:  # pragma: no cover - only without the dep
        raise ImportError(
            "pseudobulk requires decoupler; install with `pip install decoupler` "
            "(bundled in the anngg[pseudobulk] extra)."
        ) from exc
    return dc


def pseudobulk(
    adata,
    sample_col: str,
    group_by: str | None = None,
    *,
    layer: str | None = None,
    raw: bool = False,
    mode: str = "sum",
    min_cells: int = 10,
    skip_checks: bool = False,
):
    """Aggregate cells to a pseudobulk AnnData (one profile per sample x group).

    Parameters
    ----------
    sample_col
        obs column identifying the biological replicate (donor / sample / batch).
    group_by
        obs column to split within each sample (e.g. cell type); ``None`` gives one
        profile per sample.
    layer / raw
        Where the counts live. decoupler expects **integer counts**; pass the counts
        ``layer=`` (or ``raw=True``), or ``skip_checks=True`` to aggregate non-counts.
    mode
        Aggregation: ``"sum"`` (default, for count-based DE) or ``"mean"``.
    min_cells
        Drop pseudobulk profiles built from fewer than this many cells
        (decoupler records the count in ``obs['psbulk_cells']``).

    Returns
    -------
    A pseudobulk :class:`~anndata.AnnData` whose observations are sample x group
    profiles -- ready to hand straight back to any anngg helper or ``gganndata``.
    """
    dc = _require_decoupler()
    pb = dc.pp.pseudobulk(
        adata,
        sample_col=sample_col,
        groups_col=group_by,
        layer=layer,
        raw=raw,
        mode=mode,
        skip_checks=skip_checks,
    )
    if min_cells and "psbulk_cells" in pb.obs.columns:
        pb = pb[pb.obs["psbulk_cells"] >= min_cells].copy()
    return pb
