# Changelog

All notable changes are recorded here. This project follows
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Because every constant in `src/arch17p13/regions.py` changes classifications, any change to
one requires a version bump and an entry here stating what moved and why. See
[`SCOPE.md`](SCOPE.md).

## [0.2.0] — 2026-11-11

**Corrects a wrong gene coordinate in 0.1.0. Re-run anything classified with 0.1.0.**

### Fixed

- **`WRAP53_coding` was never lifted from GRCh37.** 0.1.0 carried 7,687,551-7,725,564,
  whose start is not a coordinate but "one base after TP53 ends". TP53 had been lifted
  correctly (+96,682 on both ends); WRAP53 had not. The real GRCh38 extent is
  **7,686,071-7,703,502**, obtained by lifting the Ensembl GRCh37 gene body with the same
  chain file.
  - TP53 and WRAP53 overlap by **1,479 bp** (7,686,071-7,687,550, the top of TP53
    intron 1). Under 0.1.0 they were disjoint, so a single breakpoint could never hit both
    genes: `Both-single-SV` was unreachable from one breakpoint and co-disruption counts
    read zero with every test passing.
  - 0.1.0 also counted 22,062 bp past WRAP53's real end as coding.
- `WRAP53_promoter` moves to 7,681,071-7,686,070 (5 kb upstream of the plus-strand TSS, which
  moved with the gene) and `WRAP53_5kb_flank` to 7,703,503-7,708,502.

### Changed

- **`intergenic_TP53_WRAP53` is removed.** The genes overlap; the span it named does not
  exist. Replaced by **`TP53_WRAP53_overlap`** (7,686,071-7,687,550), which does. Callers
  reading feature names will see the new name.
- Both genes now lie entirely inside the locus window, so `SCOPE.md` §3 no longer records
  WRAP53's tail being gated out — that was an artefact of the unlifted coordinate.
- Equivalence pins regenerated against the corrected reference implementation
  (`9eb36db`), now 18 cases: 16 exact agreements and the same 2 documented divergences.
  Adds `single_breakpoint_in_gene_overlap`, which 0.1.0 could not have produced.

### Added

- `regions.check_genes_overlap`, run at import: raises unless the two coding windows
  intersect. A disjoint table is the failure that hides itself, so it fails loudly instead.
- `regions.OVERLAP_WIDTH`, and tests asserting the overlap directly in both
  `test_regions.py` and `test_intervals.py`, including a breakpoint at 7,686,800 that must
  classify as `Both-single-SV`.
- `tests/data/overlap_single_breakpoint.vcf`.

### Migration

Feature-name consumers: `intergenic_TP53_WRAP53` → `TP53_WRAP53_overlap`. Categories, the
ladder, the CLI and the Python API are unchanged. Anything classified with 0.1.0 should be
re-run: calls driven by a breakpoint in 7,686,071-7,687,550 or above 7,703,502 can move.

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
