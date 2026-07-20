# Installation

`anngg` requires Python ≥ 3.12.

```bash
pip install git+https://github.com/mdmanurung/anngg
```

`anngg` builds on [plotnine](https://plotnine.org),
[plotnine-extra](https://github.com/mdmanurung/plotnine-extra) and
[annplyr](https://github.com/mdmanurung/annplyr); the latter two install
automatically from git.

## Optional extras

Some helpers use heavy, single-purpose backends kept out of the core install:

| Extra | Enables | Backend |
|-------|---------|---------|
| `density` | {func}`~anngg.plot_density` | [pyNebulosa](https://github.com/mdmanurung/pyNebulosa) |
| `upset` | {func}`~anngg.plot_upset` | [marsilea](https://marsilea.readthedocs.io/) |
| `heatmap` | {func}`~anngg.plot_clustermap` | [PyComplexHeatmap](https://github.com/DingWB/PyComplexHeatmap) |

```bash
pip install "anngg[density,upset,heatmap] @ git+https://github.com/mdmanurung/anngg"
```

Calling a helper without its extra raises an `ImportError` naming the extra to
install.
