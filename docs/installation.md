# Installation

`ggann` requires Python ≥ 3.12.

```bash
pip install git+https://github.com/mdmanurung/ggann
```

`ggann` builds on [plotnine](https://plotnine.org),
[plotnine-extra](https://github.com/mdmanurung/plotnine-extra) and
[annplyr](https://github.com/mdmanurung/annplyr); the latter two install
automatically from git.

## Optional extras

Some helpers use heavy, single-purpose backends kept out of the core install:

| Extra | Enables | Backend |
|-------|---------|---------|
| `density` | {func}`~ggann.plot_density` | [pyNebulosa](https://github.com/mdmanurung/pyNebulosa) |
| `upset` | {func}`~ggann.plot_upset` | [marsilea](https://marsilea.readthedocs.io/) |
| `heatmap` | {func}`~ggann.plot_clustermap` | [PyComplexHeatmap](https://github.com/DingWB/PyComplexHeatmap) |

```bash
pip install "ggann[density,upset,heatmap] @ git+https://github.com/mdmanurung/ggann"
```

Calling a helper without its extra raises an `ImportError` naming the extra to
install.
