"""Every plotnine-native utility helper must have a layer-by-layer grammar twin.

The twins live in ``examples/grammar_equivalents.py``. These tests build each one
so the documented grammar cannot silently drift from the helpers, and assert the
mapping covers every plotnine-native utility (the escape hatches -- plot_clustermap
and plot_upset -- are not plotnine and are intentionally excluded).
"""

from __future__ import annotations

import importlib.util
import os

import plotnine as p9
import pytest

_GE_PATH = os.path.join(os.path.dirname(__file__), "..", "examples", "grammar_equivalents.py")
_spec = importlib.util.spec_from_file_location("grammar_equivalents", _GE_PATH)
ge = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ge)


@pytest.mark.parametrize("name,builder", list(ge.EQUIVALENTS.items()))
def test_grammar_equivalent_builds(adata, name, builder):
    plot = builder(adata)
    assert isinstance(plot, p9.ggplot)
    plot._build()  # force layout so any layer mismatch surfaces here


def test_every_plotnine_utility_has_a_grammar_equivalent():
    expected = {
        "plot_embedding",
        "plot_density",
        "plot_box",
        "plot_expression_bar",
        "plot_expression_line",
        "plot_correlation",
    }
    missing = expected - set(ge.EQUIVALENTS)
    assert not missing, f"utilities lacking a grammar equivalent: {sorted(missing)}"
