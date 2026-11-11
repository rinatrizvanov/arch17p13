# arch17p13

Classify **17p13.1 TP53/WRAP53 structural-variant architecture** from a VCF.

[![CI](https://github.com/rinatrizvanov/arch17p13/actions/workflows/ci.yml/badge.svg)](https://github.com/rinatrizvanov/arch17p13/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12%20%7C%203.13-blue)](pyproject.toml)
[![Coverage](https://img.shields.io/badge/coverage-100%25-brightgreen)](pyproject.toml)
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](LICENSE)

This tool classifies 17p13.1 structural variants into the architecture categories defined in
`label_spec_v0.9` of the reference implementation, operationalizing the TP53/WRAP53
co-disruption described by Keegan et al. 2026. The reference implementation is in
[github.com/rinatrizvanov/osteosarcoma-alt-dependency-research](https://github.com/rinatrizvanov/osteosarcoma-alt-dependency-research).
Testing uses pinned outputs from DepMap inference; these are inferences, not published gold
standard.

**No runtime dependencies.** Python 3.10+, standard library only.

## The locus

TP53 and WRAP53 are convergent overlapping genes: they run head-to-head, share a
bidirectional promoter, and their bodies **intersect over 1,479 bp**. That overlap is the
single thing most likely to trip up a reimplementation, and getting it wrong is silent —
see [Known-wrong in v0.1.0](#known-wrong-in-v010).

```
chr17, GRCh38      locus window  [ 7,560,000 ....................... 7,710,000 )

   TP53  (minus strand)
   7,668,402                                           7,687,550
   +--------------------------------------------------------+
   |                      TP53_coding                        |
   +--------------------------------------------------------+
                                                 +-------------------------+
                                                 |      WRAP53_coding      |
                                                 +-------------------------+
                                                 7,686,071         7,703,502
                                                     WRAP53  (plus strand)

                                                 |<------>|
                                                     (a)
```

**(a) The overlap, chr17:7,686,071–7,687,550.** This is the top of TP53 intron 1, and it
lies inside WRAP53's body. **One breakpoint here disrupts both genes.** That is precisely
what `Both-single-SV` describes, and why it is a separate category from
`Both-independent-SVs` rather than both being folded into one "both" call.

Coordinates are Ensembl gene bodies in GRCh37 lifted with the UCSC hg19-to-hg38 chain,
which is a constant +96,682 across this locus. Both genes sit entirely inside the locus
window.

## Install

Not published to PyPI. Install from the repository:

```bash
pip install git+https://github.com/rinatrizvanov/arch17p13.git
```

Or, for development, from a checkout:

```bash
uv sync
```

## Use

```python
from arch17p13 import classify_vcf

result = classify_vcf("tumour.vcf.gz", sample="TUMOR")
print(result.category)  # 'Both-single-SV'
print(result.reason)  # why, naming the records that decided it
print(result.n_locus_records)  # distinct in-locus SV records
print(result.flags)  # frozenset({'spans_locus'}) when applicable
```

```bash
arch17p13 classify tumour.vcf.gz --sample TUMOR --explain
arch17p13 classify cohort/*.vcf.gz --json
```

`--all-filters` classifies regardless of FILTER; the default keeps `PASS` and `.` only.
Exit status is 0 for any successful classification including `Neither` — a category is a
result, not an error — 1 for an unreadable or malformed file, 2 for a usage error.

## Categories

Tested in this order. The first match wins.

| Category | Condition |
|---|---|
| `Complex` | ≥ 5 distinct in-locus SV records. **Checked first**; no gene call is attempted. |
| `Both-single-SV` | One record hits both `TP53_coding` and `WRAP53_coding`. |
| `Both-independent-SVs` | Both genes hit, by different records. |
| `TCAB1-only` | `WRAP53_coding` hit, `TP53_coding` not. |
| `TP53-only` | `TP53_coding` hit, `WRAP53_coding` not. |
| `Neither` | No in-locus records, or in-locus records hitting neither coding window. |

A record is one SV. Breakend rows linked by `INFO/MATEID` collapse to a single record, so a
reciprocal event touching each gene once is `Both-single-SV`, not `Both-independent-SVs`.

`unknown` exists in the category vocabulary but is never returned; it is there so callers
recording "not assessed" can use the same names.

## Reference Implementation

The categorization logic in this library is a reimplementation of the Stage 1 pre-registered
research code. Reviewers tracing the logic back to its research context should read these
files in `osteosarcoma-alt-dependency-research`:

| Path | What it contains |
|---|---|
| `part1_alt_classification/src/sv_annotator.py:286` | `architecture_call()` — the category decision function this library reimplements |
| `part1_alt_classification/src/sv_annotator.py:76-81` | The `ArchitectureCall` literal: the authoritative category set |
| `part1_alt_classification/src/sv_annotator.py:51-57` | `FEATURES` — GRCh38 feature regions for TP53 and WRAP53 (0-based half-open) |
| `part1_alt_classification/src/sv_annotator.py:68-70` | GRCh37 constants, retained for in-frame classification of PCGP 2014 calls |
| `part1_alt_classification/src/sv_annotator.py:112,163,213,261` | `annotate_breakpoints`, `annotate_segments`, `annotate_fusions`, `annotate_locus_svs` — the ingest paths feeding the call |
| `part1_alt_classification/tests/test_sv_annotator.py` | The original test suite for that logic |
| `part1_alt_classification/src/label_rules.py`, `manifest_builder.py` | How the category was consumed downstream into `master_manifest.csv` |
| `part2_bir_dependency_test/references/` | Frozen read-only copies pinned by `MANIFEST.sha256`, as consumed by the Stage 1 test |
| `part2_bir_dependency_test/docs/analytical_prespec.md` | The Stage 1 pre-registration this categorization was frozen under |

**One biological caveat the reference implementation documents and this library
inherits:** TP53 and WRAP53 share a head-to-head bidirectional promoter and
their gene bodies overlap by roughly 1.5 kb (GRCh37 chr17:7,589,389-7,590,868).
A single breakpoint in that window legitimately flags both genes. This is why
`Both-single-SV` and `Both-independent-SVs` are distinct categories rather than
one "both" category, and any reimplementation that collapses them loses the
distinction the research design depends on.

### Equivalence testing

`tests/test_equivalence.py` compares this library against the reference's category ladder on
17 cases, with expected values pinned by `scripts/generate_equivalence.py` from the reference
at commit `a750daf`. Fifteen agree exactly. Two diverge deliberately, each with a written
explanation and a test asserting it by name:

1. **Breakpoint semantics, not interval overlap.** The reference's segment path is a
   *forbidden variable* in the Stage 1 pre-registration — it called 33 of 35 samples
   `Both-single-SV` off one copy-number segment covering both coding windows. An SV
   covering the locus with no breakpoint inside it is `Neither` here, plus a `spans_locus`
   flag.
2. **Contig names are normalised.** The reference compares literally against `chr17`, so a
   VCF spelling the contig `17` is off-locus to it. Both spellings are common, so this library
   accepts either.

Full rationale for both, and for three further fidelity decisions, is in
[`SCOPE.md`](SCOPE.md).

## Scope

v0.1 is stable within a deliberately narrow scope: **VCF input only, GRCh38 only, breakpoint
semantics only.** BED input, GRCh37, liftover and copy-number segment semantics are out of
scope, and changing any frozen constant requires a version bump. See [`SCOPE.md`](SCOPE.md).

### What this is NOT

- **Not a clinical or diagnostic tool.**
- **Not a telomere-maintenance or ALT classifier.** It reports SV architecture at one locus.
  It infers nothing about ALT status, telomerase activity or gene dependency.
- **Not validated against a published gold standard.** The equivalence pins come from
  DepMap-based inference in the reference project.
- **Not a general VCF library.** It reads enough of the spec to find breakpoints at one 150 kb
  locus.

## Known-wrong in v0.1.0

**v0.1.0 shipped a WRAP53 interval that was never lifted from GRCh37.** It carried
`WRAP53_coding = 7,687,551–7,725,564`, where the start is not a coordinate but "one base
after TP53 ends". TP53 had been lifted correctly (+96,682); WRAP53 had not. The real extent
is 7,686,071–7,703,502.

Two consequences, opposite in sign:

- **The overlap was unrepresentable.** A breakpoint in 7,686,071–7,687,550 returned
  `TP53_coding` and `WRAP53_promoter`, never `WRAP53_coding`. `Both-single-SV` was therefore
  unreachable from a single breakpoint, and any co-disruption count read zero while every
  test stayed green.
- **22,062 bp of phantom gene body.** Breakpoints up to 7,725,564 counted as WRAP53 coding
  hits, past the gene's real end.

If you ran v0.1.0, any `TP53-only` call driven by a breakpoint in 7,686,071–7,687,550 should
have been `Both-single-SV`, and any `TCAB1-only` or both-genes call driven by a breakpoint
above 7,703,502 was spurious. Re-run with v0.2.0.

Guards added so it cannot recur: `regions.check_genes_overlap` raises at import if the two
coding windows stop intersecting, and `tests/test_regions.py::test_coding_windows_overlap`
plus `tests/test_intervals.py::test_one_breakpoint_in_the_overlap_hits_both_genes` assert
the geometry directly. The reference implementation's account is in
[`docs/erratum_wrap53_coordinates.md`](https://github.com/rinatrizvanov/osteosarcoma-alt-dependency-research/blob/main/part1_alt_classification/docs/erratum_wrap53_coordinates.md).

## Development

```bash
uv sync
uv run pytest --cov=arch17p13   # 191 tests, 100% coverage enforced
uv run ruff check .
uv run ruff format --check .
uv run mypy                     # strict, over src/ and tests/
```

All four run in CI on Python 3.10, 3.11, 3.12 and 3.13. A second job builds the wheel,
installs it into a clean environment and runs `scripts/smoke_test.sh` against the installed
console script — the test suite imports from `src/`, so that is the only check that the
packaged artefact actually works.

Regenerating the equivalence pins needs a checkout of the reference repository:

```bash
uv run --with pandas python scripts/generate_equivalence.py \
    --reference-repo ../osteosarcoma-alt-dependency-research
```

See [`CONTRIBUTING.md`](CONTRIBUTING.md).

## Citation

See [`CITATION.cff`](CITATION.cff). Licensed under the MIT Licence — see [`LICENSE`](LICENSE).
