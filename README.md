# alpha_analysis_3

Compute the stellarator alpha-particle loss surrogate f defined in Matt Landreman's
[original note](<docs/references/20261109-01 alpha loss fraction metric.tex>): the fraction
of alpha births whose trapped state is connected to the plasma boundary along contours of
constant second adiabatic invariant J, including the jumps of J at transitions between
trapping classes.

This is the third implementation attempt. It parametrizes trapped states by the ridge
surface of along-field-line maxima of B, which is independent of pitch, so that all
pitch-dependent structure comes from thresholding stored scalar fields. The approach was
prototyped end to end on the five reference equilibria before this repository was created.

## Start here

- [Design and acceptance contract](docs/DESIGN.md)
- [Milestone status](docs/STATUS.md)
- [Instructions for coding agents](AGENTS.md)
- [Milestone skill](.agents/skills/milestone/SKILL.md)
- [Lessons from the earlier attempts](docs/LESSONS.md)
- [Prototype and its evidence](prototype/README.md), [results](docs/validation/prototype/README.md)
- [GitHub and Claude review setup](docs/DEVELOPMENT.md)

Source profile h(ρ) = (1 − ρ²)⁴ with s = ρ², and h = 1 reported alongside. Numerical
accuracy target max(0.02 f, 0.01), an error magnitude. Runtime target under 300 s per
equilibrium including setup and the convergence check, memory under 2 GB. These are
targets, not measured capabilities of the bootstrap.

## Install and check

On the researcher's Mac, use the existing conda interpreter without creating a conda
environment or inheriting its packages:

```bash
/Users/mattland/opt/miniconda3/envs/20220806-03/bin/python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
make check
```

Elsewhere use Python 3.10–3.12 in a clean virtual environment.

```bash
make test       # fast tier: under 120 s; individual tests about 20 s or less
make test-full  # full tier: under 300 s; slow tests about 90 s or less
make lint
make format
make check      # lint + full tier
make smoke      # install and import in a temporary clean environment
```

## Ask an agent to continue

> Use $milestone to implement the next milestone from docs/STATUS.md. Follow
> docs/DESIGN.md, run the real-input checks, open a draft PR, wait for Tests, mark it
> ready, and address the Claude review. Do not merge.
