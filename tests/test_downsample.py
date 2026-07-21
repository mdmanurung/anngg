"""The ``downsample`` speed lever on the point/violin-heavy plots."""

from __future__ import annotations

import plotnine as p9

import ggann as ag
from ggann.plots import _downsample_cells


def test_downsample_caps_per_group(adata, group_key):
    sub = _downsample_cells(adata, group_key, 20)
    counts = sub.obs[group_key].value_counts()
    assert counts.max() <= 20
    # groups that already had <=20 cells are kept whole
    orig = adata.obs[group_key].value_counts()
    for cat in orig.index:
        assert counts.get(cat, 0) == min(orig[cat], 20)


def test_downsample_caps_total_without_group(adata):
    sub = _downsample_cells(adata, None, 100)
    assert sub.n_obs == 100


def test_downsample_noop_paths(adata, group_key):
    assert _downsample_cells(adata, group_key, None) is adata
    assert _downsample_cells(adata, group_key, adata.n_obs + 1) is adata


def test_downsample_is_deterministic(adata, group_key):
    a = _downsample_cells(adata, group_key, 15)
    b = _downsample_cells(adata, group_key, 15)
    assert list(a.obs_names) == list(b.obs_names)


def test_violin_family_builds_with_downsample(adata, markers, group_key):
    for plot in (
        ag.plot_violin(adata, markers[:2], group_key, downsample=25),
        ag.plot_box(adata, markers[:2], group_key, downsample=25),
        ag.plot_stacked_violin(adata, markers, group_key, downsample=25),
    ):
        assert isinstance(plot, p9.ggplot)
        plot._build()


def test_embedding_downsample_reduces_points(adata, group_key):
    plot = ag.plot_embedding(adata, "umap", color=group_key, downsample=50)
    assert isinstance(plot, p9.ggplot)
    # the coordinate frame handed to plotnine holds only the kept cells
    assert len(plot.data) == 50
    plot._build()


def test_features_builds_with_downsample(adata, markers):
    plot = ag.plot_features(adata, markers[:2], basis="umap", downsample=50)
    assert isinstance(plot, p9.ggplot)
    assert plot.data["expression"].notna().any()
    plot._build()
