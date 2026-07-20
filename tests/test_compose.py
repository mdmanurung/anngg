"""Tests for compose / tag_panels (figure assembly)."""

from __future__ import annotations

import plotnine as p9
import pytest

import ggann as ag
from ggann.layout import _tag_labels


def _panels(adata, group_key, markers, n=4):
    builders = [
        lambda: ag.plot_embedding(adata, "umap", color=group_key),
        lambda: ag.plot_dotplot(adata, markers, group_key),
        lambda: ag.plot_violin(adata, markers[:1], group_key),
        lambda: ag.plot_proportions(adata, group_key, split_by="phase"),
    ]
    return [b() for b in builders[:n]]


def test_tag_labels():
    assert _tag_labels("A", 3) == ["A", "B", "C"]
    assert _tag_labels("a", 2) == ["a", "b"]
    assert _tag_labels("1", 3) == ["1", "2", "3"]
    assert _tag_labels("i", 4) == ["i", "ii", "iii", "iv"]
    with pytest.raises(ValueError, match="tag_levels"):
        _tag_labels("X", 2)


def test_compose_grid_builds(adata, group_key, markers):
    fig = ag.compose(_panels(adata, group_key, markers, 4), ncol=2)
    fig.save  # composition object is saveable
    # rendering the whole composition exercises every panel + the layout
    fig.draw(show=False)


def test_compose_default_and_single(adata, group_key, markers):
    ag.compose(_panels(adata, group_key, markers, 3)).draw(show=False)  # default shape
    ag.compose(_panels(adata, group_key, markers, 1)).draw(show=False)  # single panel


def test_compose_no_tags(adata, group_key, markers):
    ag.compose(_panels(adata, group_key, markers, 2), tag_levels=None).draw(show=False)


def test_tag_panels_returns_tagged_plots(adata, group_key, markers):
    tagged = ag.tag_panels(_panels(adata, group_key, markers, 2), levels="a")
    assert len(tagged) == 2
    assert all(isinstance(p, p9.ggplot) for p in tagged)


def test_compose_empty_raises():
    with pytest.raises(ValueError, match="at least one"):
        ag.compose([])
