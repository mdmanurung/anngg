"""Tests for split_by (stat plots) and agg (bar/line)."""

from __future__ import annotations

import numpy as np
import plotnine as p9
import pytest

import anngg as ag


def _build(plot):
    assert isinstance(plot, p9.ggplot)
    plot._build()
    return plot


@pytest.mark.parametrize(
    "call",
    [
        lambda ad, g, mk: ag.plot_violin(ad, mk[:2], g, split_by="phase"),
        lambda ad, g, mk: ag.plot_box(ad, mk[:1], g, split_by="phase"),
        lambda ad, g, mk: ag.plot_dotplot(ad, mk, g, split_by="phase"),
        lambda ad, g, mk: ag.plot_matrixplot(ad, mk, g, split_by="phase", standard_scale="var"),
        lambda ad, g, mk: ag.plot_expression_bar(ad, mk[:2], g, split_by="phase"),
    ],
)
def test_split_by_builds(adata, group_key, markers, call):
    _build(call(adata, group_key, markers))


def test_agg_median_bar_and_line(adata, group_key, markers):
    bar = ag.plot_expression_bar(adata, markers[:2], group_key, agg="median")
    _build(bar)
    assert bar.labels.y == "median expression"
    line = ag.plot_expression_line(adata, markers[:1], x="phase", group_by=group_key, agg="median")
    _build(line)
    assert line.labels.y == "median expression"


def test_agg_callable(adata, group_key, markers):
    _build(ag.plot_expression_bar(adata, markers[:1], group_key, agg=np.sum))


def test_split_by_no_regression(adata, group_key, markers):
    # the default (no split_by) path must be unchanged
    _build(ag.plot_violin(adata, markers, group_key))
    _build(ag.plot_dotplot(adata, markers, group_key))
