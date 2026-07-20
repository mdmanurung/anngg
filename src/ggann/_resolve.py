"""Aesthetic resolution: turn names into a tidy per-cell DataFrame via annplyr.

Every plot in ggann is built from a plain :class:`pandas.DataFrame` that is
extracted from the ``AnnData`` *only* through the ``adata.ap`` (annplyr)
accessor -- no direct indexing into ``adata.X`` / ``adata.obs``. This module is
that single extraction layer; the plotting helpers sit on top of it.

Resolution order for a bare name (matching scanpy's precedence, least
surprising first):

1. a column of ``adata.obs``            (per-cell metadata)
2. a gene / feature in ``X`` or a layer (or ``adata.raw`` when ``use_raw``)
3. an embedding coordinate in ``obsm``  (e.g. ``"UMAP_1"``)

When a name is both an obs column and a gene, obs wins and a warning is
emitted. Use the :func:`gene` and :func:`obs` escapes to force a source.
"""

from __future__ import annotations

import re
import warnings
from dataclasses import dataclass
from typing import Iterable

import annplyr as ap
import numpy as np
import pandas as pd

__all__ = ["gene", "obs", "obsm", "embedding_coords"]

_OBSM_TOKEN = re.compile(r"^(?P<basis>.+?)\[(?P<idx>\d+)\]$")


@dataclass(frozen=True)
class Ref:
    """An explicit reference to a data source, produced by :func:`gene`/:func:`obs`."""

    name: str
    source: str  # "gene" | "obs"
    layer: str | None = None
    use_raw: bool | None = None

    def __str__(self) -> str:  # so it can be dropped straight into aes()
        return self.name


def gene(name: str, *, layer: str | None = None, use_raw: bool | None = None) -> Ref:
    """Force ``name`` to resolve as a gene/feature (expression), never obs.

    Optionally pin *this gene's* expression matrix, independent of the plot-wide
    ``layer`` / ``use_raw``. This is how you mix sources in one plot::

        aes(color=gene("CD3D", layer="logcounts"),
            size=gene("CD8A", layer="counts"))

    With neither ``layer`` nor ``use_raw`` given, the gene inherits the plot-wide
    source (which itself defaults to ``adata.raw`` when present).
    """
    return Ref(str(name), "gene", layer=layer, use_raw=use_raw)


def obs(name: str) -> Ref:
    """Force ``name`` to resolve as an ``adata.obs`` column, never a gene."""
    return Ref(str(name), "obs")


@dataclass(frozen=True)
class ObsmRef:
    """A reference to one embedding coordinate, produced by :func:`obsm`."""

    basis: str
    index: int

    def __str__(self) -> str:
        return f"{self.basis}[{self.index}]"


def obsm(basis: str, index: int) -> ObsmRef:
    """Force a specific embedding coordinate, e.g. ``obsm("umap", 0)`` -> ``UMAP_1``."""
    return ObsmRef(str(basis), int(index))


# --------------------------------------------------------------------------- #
# Prefix-string parsing: "obs:phase", "gene:CD3D@logcounts", "obsm:umap[0]"
# --------------------------------------------------------------------------- #
def parse_token(item):
    """Normalise one aes value into a :class:`Ref` / :class:`ObsmRef` / bare string.

    Bare strings (no recognised prefix) are returned unchanged for
    auto-resolution. Recognised strict prefixes:

    * ``obs:<col>``               -- an ``adata.obs`` column
    * ``gene:<name>``             -- expression, using the plot-wide source
    * ``gene:<name>@<layer>``     -- expression from ``adata.layers[<layer>]``
    * ``gene:<name>@raw``/``@X``  -- expression from ``adata.raw`` / ``adata.X``
    * ``obsm:<basis>[<i>]``       -- coordinate ``i`` (0-based) of an embedding

    Non-string, non-Ref values are returned unchanged.
    """
    if isinstance(item, (Ref, ObsmRef)):
        return item
    if not isinstance(item, str) or ":" not in item:
        return item
    prefix, _, body = item.partition(":")
    kind = prefix.strip().lower()
    if kind == "obs":
        return Ref(body.strip(), "obs")
    if kind == "gene":
        name, at, source = body.partition("@")
        name, source = name.strip(), source.strip()
        if not at:
            return Ref(name, "gene")
        low = source.lower()
        if low == "raw":
            return Ref(name, "gene", use_raw=True)
        if low in ("x", ".x"):
            return Ref(name, "gene", use_raw=False)
        return Ref(name, "gene", layer=source)
    if kind == "obsm":
        m = _OBSM_TOKEN.match(body.strip())
        if not m:
            raise ValueError(
                f"obsm token {item!r} must look like 'obsm:umap[0]' (0-based index)."
            )
        return ObsmRef(m.group("basis"), int(m.group("idx")))
    return item  # unrecognised prefix -> leave for plotnine / auto-resolution


def plain_name(adata, item):
    """The plain DataFrame column name a token resolves to (for aes rewriting)."""
    tok = parse_token(item)
    if isinstance(tok, ObsmRef):
        key = embedding_key(adata, tok.basis)
        return f"{_embedding_prefix(key)}_{tok.index + 1}"
    if isinstance(tok, Ref):
        return tok.name
    return tok


# --------------------------------------------------------------------------- #
# Expression source (X / layer / raw) -- the single source of truth
# --------------------------------------------------------------------------- #
def resolve_source(adata, layer: str | None, use_raw: bool | None) -> tuple[str, str | None]:
    """Decide which matrix expression comes from.

    Returns ``("layer", name)``, ``("raw", None)`` or ``("x", None)``. Mirrors
    scanpy: an explicit ``layer`` wins; otherwise ``use_raw`` defaults to ``True``
    when ``adata.raw`` exists. Raises rather than silently substituting a
    different matrix -- both the grammar path and the aggregation path go through
    here, so the two can never disagree.
    """
    if layer is not None:
        if use_raw is True:
            raise ValueError("Cannot specify use_raw=True and a layer at the same time.")
        return "layer", layer
    if use_raw is None:
        use_raw = adata.raw is not None
    if use_raw:
        if adata.raw is None:
            raise ValueError("use_raw=True but adata.raw is None.")
        return "raw", None
    return "x", None


def _gene_universe(adata, kind: str) -> set:
    if kind == "raw":
        return set(adata.raw.var_names)
    return set(adata.var_names)


def _other_gene_universe(adata, kind: str) -> set:
    if kind == "raw":
        return set(adata.var_names)
    return set(adata.raw.var_names) if adata.raw is not None else set()


def _source_label(kind: str, layer: str | None) -> str:
    return {"raw": ".raw", "x": ".X"}.get(kind, f"layer '{layer}'")


# --------------------------------------------------------------------------- #
# Embeddings
# --------------------------------------------------------------------------- #
def _embedding_prefix(key: str) -> str:
    """``X_umap`` -> ``UMAP``, ``X_pca`` -> ``PC`` (Seurat-style coordinate names)."""
    base = key[2:] if key.startswith("X_") else key
    base = base.upper()
    return {"PCA": "PC"}.get(base, base)


def embedding_key(adata, basis: str) -> str:
    """Resolve a user-facing basis name (``"umap"``, ``"X_umap"``, ``"UMAP"``) to an obsm key."""
    candidates = [basis, f"X_{basis}", basis.lower(), f"X_{basis.lower()}"]
    for cand in candidates:
        if cand in adata.obsm:
            return cand
    lower = {k.lower(): k for k in adata.obsm}
    for cand in (basis.lower(), f"x_{basis.lower()}"):
        if cand in lower:
            return lower[cand]
    raise KeyError(
        f"No embedding '{basis}' in adata.obsm (available: {list(adata.obsm)})"
    )


def embedding_coords(adata, basis: str, n: int = 2) -> pd.DataFrame:
    """Return the first ``n`` coordinates of an embedding as a tidy DataFrame.

    Columns are named ``<PREFIX>_<i>`` (e.g. ``UMAP_1``, ``UMAP_2``) to match the
    conventions used by ``plotnine_extra.DimPlot``.
    """
    key = embedding_key(adata, basis)
    df = adata.ap.to_df(obsm={key: ap.everything()}).iloc[:, :n]
    prefix = _embedding_prefix(key)
    df.columns = [f"{prefix}_{i + 1}" for i in range(df.shape[1])]
    return _densify(df)


def _all_embedding_coords(adata, n: int = 2) -> dict[str, tuple[str, int]]:
    """Map every embedding coordinate name -> (obsm key, column index)."""
    out: dict[str, tuple[str, int]] = {}
    for key, arr in adata.obsm.items():
        if getattr(arr, "ndim", 0) != 2:
            continue
        prefix = _embedding_prefix(key)
        for i in range(min(n, arr.shape[1])):
            out.setdefault(f"{prefix}_{i + 1}", (key, i))
    return out


# --------------------------------------------------------------------------- #
# Densification (annplyr returns sparse-backed columns for X/raw)
# --------------------------------------------------------------------------- #
def _densify(df: pd.DataFrame) -> pd.DataFrame:
    for col in df.columns:
        if isinstance(df[col].dtype, pd.SparseDtype):
            df[col] = np.asarray(df[col].sparse.to_dense(), dtype=float)
    return df


def _raw_wide(adata, var_names) -> pd.DataFrame:
    """Wide cells x genes DataFrame from ``adata.raw``, with the ``raw_`` prefix stripped."""
    var_names = list(var_names)
    wide = adata.ap.to_df(raw=var_names).rename(columns={f"raw_{g}": g for g in var_names})
    return _densify(wide)


# --------------------------------------------------------------------------- #
# Frame resolution
# --------------------------------------------------------------------------- #
def resolve_frame(
    adata,
    names: Iterable,
    *,
    layer: str | None = None,
    use_raw: bool | None = None,
    warn_collisions: bool = True,
) -> pd.DataFrame:
    """Build a per-cell DataFrame with one column per requested ``name``.

    ``names`` may contain plain strings or :class:`Ref` objects (from
    :func:`gene`/:func:`obs`). Names that resolve to nothing are silently
    skipped -- they are assumed to be computed aesthetics that plotnine will
    handle from already-present columns. A name that is a gene only in the
    *inactive* matrix (e.g. in ``.X`` while reading from ``.raw``) is skipped
    *with a warning*, so the mistake is visible.
    """
    default_kind, default_layer = resolve_source(adata, layer, use_raw)
    default_set = _gene_universe(adata, default_kind)
    other_set = _other_gene_universe(adata, default_kind)
    obs_cols = set(adata.obs.columns)
    emb = _all_embedding_coords(adata)

    obs_names: list[str] = []
    gene_specs: list[tuple[str, str, str | None]] = []  # (name, kind, layer)
    obsm_specs: list[tuple[str, str, int]] = []  # (col_name, obsm_key, index)
    order: list[str] = []

    def _ref_source(ref: Ref) -> tuple[str, str | None]:
        # an explicit per-gene layer/use_raw fully determines the source;
        # otherwise inherit the plot-wide default.
        if ref.layer is not None:
            return resolve_source(adata, ref.layer, None)
        if ref.use_raw is not None:
            return resolve_source(adata, None, ref.use_raw)
        return default_kind, default_layer

    for raw_item in names:
        tok = parse_token(raw_item)

        if isinstance(tok, ObsmRef):
            key = embedding_key(adata, tok.basis)
            name = f"{_embedding_prefix(key)}_{tok.index + 1}"
            if name in order:
                continue
            obsm_specs.append((name, key, tok.index))
            order.append(name)
            continue

        if isinstance(tok, Ref):
            name = tok.name
            if name in order:
                continue
            if tok.source == "obs":
                _dispatch(name, "obs", obs_cols, obs_names)
            else:  # forced gene
                k, lyr = _ref_source(tok)
                if name not in _gene_universe(adata, k):
                    raise KeyError(f"gene('{name}') not found in {_source_label(k, lyr)}.")
                gene_specs.append((name, k, lyr))
            order.append(name)
            continue

        # bare string -> auto-resolve (obs -> gene -> obsm coordinate)
        name = str(tok)
        if name in order:
            continue
        if name in obs_cols:
            if warn_collisions and name in default_set:
                warnings.warn(
                    f"'{name}' is both an obs column and a gene; using obs. "
                    f"Use 'gene:{name}' to plot expression.",
                    stacklevel=2,
                )
            obs_names.append(name)
        elif name in default_set:
            gene_specs.append((name, default_kind, default_layer))
        elif name in emb:
            key, idx = emb[name]
            obsm_specs.append((name, key, idx))
        elif name in other_set:
            warnings.warn(
                f"'{name}' is a gene in "
                f"{_source_label('raw' if default_kind != 'raw' else 'x', default_layer)} "
                f"but not in the active expression source "
                f"({_source_label(default_kind, default_layer)}); "
                f"pass use_raw={default_kind != 'raw'} to plot it.",
                stacklevel=2,
            )
            continue
        else:
            continue  # computed / literal aesthetic -- leave to plotnine
        order.append(name)

    frames = []
    if obs_names:
        frames.append(_densify(adata.ap.to_df(obs=obs_names)))

    # group genes by their (matrix, layer) so mixed per-aesthetic sources each
    # get one extraction pass
    by_source: dict[tuple[str, str | None], list[str]] = {}
    for name, k, lyr in gene_specs:
        by_source.setdefault((k, lyr), []).append(name)
    for (k, lyr), gnames in by_source.items():
        if k == "raw":
            frames.append(_raw_wide(adata, gnames))
        else:
            kw: dict = {"x": gnames}
            if lyr is not None:
                kw["layer"] = lyr
            frames.append(_densify(adata.ap.to_df(**kw)))

    # group embedding coordinates by obsm key (one extraction, take the columns)
    by_key: dict[str, list[tuple[str, int]]] = {}
    for name, key, idx in obsm_specs:
        by_key.setdefault(key, []).append((name, idx))
    for key, cols in by_key.items():
        coords = embedding_coords(adata, key, n=max(i for _, i in cols) + 1)
        for name, idx in cols:
            frames.append(coords.iloc[:, [idx]].rename(columns={coords.columns[idx]: name}))

    if not frames:
        return pd.DataFrame(index=adata.obs_names)

    out = pd.concat(frames, axis=1)
    ordered = [c for c in order if c in out.columns]
    return out[ordered + [c for c in out.columns if c not in ordered]]


def _dispatch(name, kind, universe, bucket):
    if name not in universe:
        raise KeyError(f"{kind}('{name}') not found in adata.")
    bucket.append(name)
