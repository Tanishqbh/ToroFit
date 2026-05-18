# AI Usage Disclosure

This document discloses the use of AI tools in the development of ToroFit, in the interest of transparency and reproducibility.

## Tools used

- **[Claude Code](https://claude.com/claude-code)** (Anthropic) — an interactive coding assistant powered by Claude. Used during development by the maintainers.

## Scope of AI assistance

AI assistance has contributed to the following areas:

- **Code review and refactoring**: efficiency improvements (vectorisation, removal of redundant computation), correctness fixes (e.g. a bug in `flatten_grid` for the `d = 1` case), and cleanup of unused imports.
- **Packaging metadata**: drafting `pyproject.toml` dependency declarations, the `LICENSE` file, and the `CITATION.cff` file.
- **Documentation**: drafting and editing the `README.md`, `TODO.md`, and inline docstrings.
- **Strategic planning**: outlining the public-release roadmap recorded in `TODO.md`.

## What AI was *not* used for

- The underlying mathematical method (fitting Fourier series on toroidal manifolds via persistent-cohomology circular coordinates) is original work by the authors.
- The core algorithmic design and the choice to integrate with DREiMaC predate and are independent of any AI involvement.

## Responsibility

All AI-generated suggestions were reviewed by the authors before being incorporated. The authors take full responsibility for the correctness, licensing, and content of the code in this repository.

## Why disclose

We disclose AI usage because:

1. It is honest about the tooling involved in producing the codebase.
2. Increasingly, journals (including JOSS) and conferences ask for such disclosures during peer review.
3. Reproducibility benefits when readers know the provenance of design and implementation choices.

This file will be updated if the role of AI in development changes substantively.
