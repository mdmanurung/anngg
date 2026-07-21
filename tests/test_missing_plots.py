"""The missing-plots backlog surfaced by the scanpy comparison."""

from __future__ import annotations

import numpy as np
import pandas as pd
import plotnine as p9
import pytest

import ggann as ag


def test_variance_ratio_builds(adata):
    plot = ag.plot_variance_ratio(adata, n_pcs=30)
    assert isinstance(plot, p9.ggplot)
    assert len(plot.data) == 30  # capped at n_pcs
    plot._build()


def test_variance_ratio_missing_pca_raises(adata):
    ad = adata.copy()
    ad.uns.pop("pca", None)
    with pytest.raises(KeyError):
        ag.plot_variance_ratio(ad)


def test_embedding_density_grouped_and_pooled(adata, group_key):
    grouped = ag.plot_embedding_density(adata, "umap", group_key)
    pooled = ag.plot_embedding_density(adata, "umap")
    for plot in (grouped, pooled):
        assert isinstance(plot, p9.ggplot)
        # density is min-max scaled into [0, 1]
        assert plot.data["density"].between(0, 1).all()
        plot._build()


def test_heatmap_builds_and_scales(adata, markers, group_key):
    plain = ag.plot_heatmap(adata, markers, group_key, use_raw=True)
    scaled = ag.plot_heatmap(adata, markers, group_key, use_raw=True, standard_scale="var")
    assert isinstance(plain, p9.ggplot)
    plain._build()
    # each gene scaled independently into [0, 1]
    assert scaled.data["value"].between(0, 1).all()
    scaled._build()


def test_dendrogram_builds_and_autocomputes(adata, group_key):
    ad = adata.copy()  # auto-compute writes to uns; isolate the shared fixture
    ad.uns.pop(f"dendrogram_{group_key}", None)
    plot = ag.plot_dendrogram(ad, group_key)
    assert isinstance(plot, p9.ggplot)
    assert f"dendrogram_{group_key}" in ad.uns  # computed on demand
    plot._build()
    ag.plot_dendrogram(ad, group_key, orientation="left")._build()


def test_dendrogram_bad_orientation(adata, group_key):
    with pytest.raises(ValueError):
        ag.plot_dendrogram(adata, group_key, orientation="sideways")


def test_sina_builds_multi_and_single(adata, markers, group_key):
    multi = ag.plot_sina(adata, markers, group_key, use_raw=True, downsample=100)
    single = ag.plot_sina(adata, markers[:1], group_key, use_raw=True, violin=False)
    for plot in (multi, single):
        assert isinstance(plot, p9.ggplot)
        plot._build()


def _fake_de(n=400, seed=0):
    rng = np.random.RandomState(seed)
    return pd.DataFrame(
        {
            "baseMean": rng.gamma(2.0, 50.0, n),
            "log2FoldChange": rng.normal(0, 2, n),
            "padj": rng.uniform(0, 1, n),
        },
        index=[f"G{i}" for i in range(n)],
    )


def test_ma_builds_and_flags_significant():
    de = _fake_de()
    plot = ag.plot_ma(de, label_top=5)
    assert isinstance(plot, p9.ggplot)
    expected = ((de["baseMean"] > 0) & (de["padj"] < 0.05)).sum()
    assert plot.data["significant"].sum() == expected
    plot._build()


def test_ma_missing_columns_raises():
    with pytest.raises(KeyError):
        ag.plot_ma(pd.DataFrame({"baseMean": [1.0], "log2FoldChange": [0.5]}))
