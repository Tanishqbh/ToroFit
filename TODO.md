# ToroFit — Public Alpha Release (Phase 1)

Execution checklist. Design rationale lives in `~/.claude/plans/i-wish-to-market-distributed-aurora.md`.

**Goal:** ship `0.1.0a0` to PyPI so the repo can start accumulating users, stars, and citations during the months before any JOSS submission. Bio integration and JOSS paper are Phases 2 and 3 (sketched at the bottom of the plan file).

## Phase 1

- [x] **1. Fix packaging metadata** — `pyproject.toml`
  Set `version = "0.1.0a0"`, real authors (Tanishq Bhatia + Jose Perea), `dependencies` with lower-bound pins, `[project.optional-dependencies]` for `anndata` / `dev`, classifiers `Intended Audience :: Science/Research` + `Topic :: Scientific/Engineering :: Mathematics`.

- [x] **2. Demote `requirements.txt`**
  Renamed to `requirements-lock.txt`. Dropped editable self-install line. Header note explains it's dev-only.

- [x] **3. Add `LICENSE` file at repo root**
  MIT, © 2026 Tanishq Bhatia, Jose Perea.

- [x] **4. Add `CITATION.cff`**
  Both authors with affiliation Northeastern. ORCID fields are placeholder TODOs (needed before JOSS).

- [x] **5. Rewrite `README.md`**
  Pitch, alpha banner, badges, `pip install torofit` headline, collapsible install-notes section, kept quickstart, "What's coming" roadmap, citation pointer.

- [ ] **6. PyPI publish** — user-driven (needs your PyPI credentials)
  ```bash
  # 1. Verify the name is free
  pip index versions torofit  # should error "no matching distribution"

  # 2. TestPyPI smoke
  twine upload --repository testpypi dist/*
  # in a clean venv:
  pip install --index-url https://test.pypi.org/simple/ \
              --extra-index-url https://pypi.org/simple/ torofit

  # 3. Production
  twine upload dist/*

  # 4. Tag the release in git
  git tag v0.1.0a0 && git push origin v0.1.0a0
  # Then create a matching GitHub Release.
  ```
  Build artifacts are already in `dist/` (`torofit-0.1.0a0.tar.gz` + wheel).

## Verification

- [x] **Existing tests pass** — `pytest` reports 9/9 green after the packaging changes.
- [x] **Fresh build succeeds** — `python -m build` produced `torofit-0.1.0a0.tar.gz` and `torofit-0.1.0a0-py3-none-any.whl`. Wheel METADATA verified: real deps, both authors, LICENSE embedded, README as long description.
- [ ] **Fresh-env install via TestPyPI** — `conda create -n torofit-test python=3.11 && pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ torofit && python -c "import torofit; print(torofit.__version__)"`. Run *after* the TestPyPI upload in step 6.
- [ ] **README renders correctly on GitHub** — eyeball the rendered README and the GitHub Release page after tagging; confirm badges + alpha banner display as expected.

---

## Phase 2 (deferred — bio integration, target `0.2.0`)

Sketched only — full design lives in the plan file.

- AnnData adapter `torofit/anndata.py` exposed as `torofit.an`
- Cell-cycle demo notebook (`docs/notebooks/03_cell_cycle.ipynb`)
- README updated with comp-bio pitch
- GitHub Actions CI (pytest matrix, ruff)
- Tests for the partition pipeline + AnnData round-trip

## Phase 3 (deferred further — JOSS submission)

- `paper.md` + `paper.bib`
- ORCID iDs on every author
- `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, issue/PR templates
- Submit via JOSS editorial workflow (after ~6 months of public usage)
