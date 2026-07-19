import matplotlib
import pytest

matplotlib.use("Agg")


@pytest.fixture(scope="session")
def adata():
    import scanpy as sc

    return sc.datasets.pbmc68k_reduced()


@pytest.fixture(scope="session")
def markers(adata):
    candidates = ["CD3D", "NKG7", "CST3", "GNLY", "MS4A1", "FCGR3A", "CD8A"]
    return [g for g in candidates if g in adata.raw.var_names][:4]


@pytest.fixture(scope="session")
def group_key():
    return "bulk_labels"
