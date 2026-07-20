"""Tests for plot_ridge and the plot_proportions kind= variants."""

from __future__ import annotations

import numpy as np
import plotnine as p9
import pytest
import scipy.sparse as sp
from anndata import AnnData

import anngg as ag


def _build(plot):
    assert isinstance(plot, p9.ggplot)
    plot._build()
    return plot


def test_plot_ridge(adata, markers, group_key):
    _build(ag.plot_ridge(adata, markers[:2], group_key))
    _build(ag.plot_ridge(adata, markers[:1], group_key))


def test_plot_ridge_degenerate_group(adata, markers, group_key):
    # a constant / tiny group must draw a flat baseline, not crash the KDE
    ad = adata.copy()
    X = ad.raw.X.toarray() if sp.issparse(ad.raw.X) else np.asarray(ad.raw.X)
    smallest = ad.obs[group_key].value_counts().index[-1]
    X[(ad.obs[group_key] == smallest).to_numpy(), :] = 0.0
    ad.raw = AnnData(X, var=ad.raw.var.copy())
    _build(ag.plot_ridge(ad, markers[:1], group_key))


def test_proportions_kinds(adata, group_key):
    _build(ag.plot_proportions(adata, group_key, split_by="phase", kind="area"))
    _build(ag.plot_proportions(adata, group_key, split_by="phase", kind="trend"))
    _build(ag.plot_proportions(adata, group_key, split_by="phase", kind="bar"))


def test_proportions_area_trend_need_split(adata, group_key):
    with pytest.raises(ValueError, match="split_by"):
        ag.plot_proportions(adata, group_key, kind="area")
    with pytest.raises(ValueError, match="split_by"):
        ag.plot_proportions(adata, group_key, kind="trend")


def test_proportions_bad_kind(adata, group_key):
    with pytest.raises(ValueError, match="kind"):
        ag.plot_proportions(adata, group_key, kind="pie")
