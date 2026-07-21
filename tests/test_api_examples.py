"""Guard the per-function API example figures.

Every builder in examples/api_examples.py must produce a figure, so the images
injected into the API-reference pages cannot silently go stale.
"""

from __future__ import annotations

import importlib.util
import os

import plotnine as p9
import pytest

_PATH = os.path.join(os.path.dirname(__file__), "..", "examples", "api_examples.py")
_spec = importlib.util.spec_from_file_location("api_examples", _PATH)
apiex = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(apiex)


# plotnine-extra layout modifiers (not standalone figures)
_NON_FIGURE = {"plot_annotation", "plot_layout"}


@pytest.fixture(scope="module")
def de_adata(adata):
    import scanpy as sc

    # work on a COPY: pbmc68k_reduced ships a logreg rank_genes_groups (no
    # p-values) that another test asserts on; recompute wilcoxon here without
    # mutating the shared session fixture.
    ad = adata.copy()
    sc.tl.rank_genes_groups(ad, "bulk_labels", method="wilcoxon", n_genes=50)
    return ad


def test_every_plot_function_has_an_api_example(de_adata):
    import ggann as ag

    documented = {
        n
        for n in ag.__all__
        if (n.startswith("plot_") or n == "gganndata") and n not in _NON_FIGURE
    }
    covered = {name.split(".", 1)[1] for name in apiex._examples(de_adata)}
    missing = documented - covered
    assert not missing, f"plot functions without an API example: {sorted(missing)}"


def test_api_examples_build(de_adata):
    for name, builder in apiex._examples(de_adata).items():
        obj = builder()
        assert obj is not None, name
        if isinstance(obj, p9.ggplot):
            obj._build()
