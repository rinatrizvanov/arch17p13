# Changelog

All notable changes are recorded here. This project follows
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Because every constant in `src/arch17p13/regions.py` changes classifications, any change to
one requires a version bump and an entry here stating what moved and why. See
[`SCOPE.md`](SCOPE.md).

## [0.1.0] — 2026-11-09

First release. Stable within the frozen v0.1 scope: **VCF input only, GRCh38 only, breakpoint
semantics only.**

### Added

- Six-category architecture classification — `Complex`, `Both-single-SV`,
  `Both-independent-SVs`, `TCAB1-only`, `TP53-only`, `Neither` — mirroring the reference
  implementation's ladder, `Complex` precedence included.
- VCF reader using the standard library only, with **no runtime dependencies**. Handles
  symbolic ALTs with `INFO/END` or `INFO/SVLEN`, all four BND bracket orientations,
  `INFO/MATEID` pairing, `FILTER` gating and multi-sample genotype selection.
- `classify_vcf()` and `classify_records()` Python API returning an immutable `Result` with
  the category, a reason string naming the deciding records, the in-locus record count and
  advisory flags.
- `arch17p13 classify` CLI with `--sample`, `--all-filters`, `--explain` and `--json`;
  also runnable as `python -m arch17p13`.
- Stage 1 equivalence suite: 17 cases pinned from the reference implementation at commit
  `a750daf` by `scripts/generate_equivalence.py`. Fifteen agree exactly; the two divergences
  are documented and asserted by name.
- 191 tests at 100% statement and branch coverage, enforced by a coverage gate.
- GitHub Actions CI running pytest, `ruff check`, `ruff format --check` and
  `mypy --strict` on Python 3.10, 3.11, 3.12 and 3.13, plus a job that installs the
  built wheel into a clean environment and smoke-tests the installed console script.
- `py.typed`, so type information is usable by downstream `mypy` under PEP 561.

### Known divergences from the reference implementation

Both are deliberate and tested; full rationale in [`SCOPE.md`](SCOPE.md).

- **Breakpoint semantics, not interval overlap.** The reference's segment path is a forbidden
  variable in the Stage 1 pre-registration, having called 33 of 35 samples `Both-single-SV`
  from one copy-number segment spanning the adjacent coding windows. An SV covering the locus
  without a breakpoint inside it is `Neither` plus a `spans_locus` flag.
- **Contig names are normalised.** `chr17`, `17` and case variants are all accepted. The
  reference compares literally against `chr17`, so it would classify a `17`-spelled VCF as
  entirely off-locus.

### Not included

BED input, GRCh37/hg19 support, liftover, copy-number segment semantics, and any ALT-biology
or gene-dependency inference. Adding any of these is a new minor version with its own scope
note.

Not published to PyPI either; install from the repository with
`pip install git+https://github.com/rinatrizvanov/arch17p13.git`. Packaging for PyPI is its
own piece of work and is not part of v0.1.
