# Scope — v0.1

`arch17p13` v0.1 is **stable within a deliberately narrow scope**. Everything below the
"In scope" heading is frozen: changing any of it changes classifications, so it requires a
version bump and a `CHANGELOG.md` entry stating what moved and why.

## In scope

| | v0.1 |
|---|---|
| Input | VCF 4.x, plain or gzip/bgzip compressed. One sample per call. |
| Assembly | GRCh38 only. |
| Locus | `chr17:7,560,000-7,710,000` (half-open). |
| Semantics | **Breakpoint positions only.** |
| Output | One of six categories, plus a reason string and advisory flags. |

Categories, in the order the ladder tests them:

| Category | Condition |
|---|---|
| `Complex` | ≥ 5 distinct in-locus SV records. Checked first; no gene call is attempted. |
| `Both-single-SV` | One record hits both `TP53_coding` and `WRAP53_coding`. |
| `Both-independent-SVs` | Both genes hit, by different records. |
| `TCAB1-only` | `WRAP53_coding` hit, `TP53_coding` not. |
| `TP53-only` | `TP53_coding` hit, `WRAP53_coding` not. |
| `Neither` | No in-locus records, or in-locus records that hit neither coding window. |

`unknown` is declared in `CATEGORIES` but is never returned. It exists so callers recording
"not assessed" can use the same vocabulary.

## Out of scope for v0.1

- **BED input.** VCF only.
- **GRCh37 / hg19.** No liftover, no in-frame classification of GRCh37 calls.
- **Segment and copy-number overlap semantics.** See the divergence below.
- **Any ALT biology**: no telomere-maintenance inference, no dependency testing, no
  expression stratification. This is a classifier library, nothing more.

## Fidelity to the reference implementation

The reference is `part1_alt_classification/src/sv_annotator.py` in
[osteosarcoma-alt-dependency-research](https://github.com/rinatrizvanov/osteosarcoma-alt-dependency-research).
Five places where this library had to make an explicit choice:

**1. Breakpoint semantics, not interval overlap — a deliberate divergence.**
The reference has two ingest paths. `annotate_breakpoints` asks whether a breakpoint falls
inside a feature; `annotate_segments` asks whether an interval overlaps one. The segment
path produced `sv_architecture_auto`, which the Stage 1 pre-registration lists as a
**forbidden variable**: *"artefactual, 33 of 35 rows reading 'Both-single-SV' from a single
CN segment spanning both WRAP53 and TP53"*
(`autoresearch/iteration_01_tcab1_predicts_alt_status/prespec.md`). The cause is geometric — TP53 and
WRAP53 overlap, so *any* interval covering 7,686,071-7,687,550 overlaps both genes, and
any interval covering both bodies is only ~35 kb wide. Broad 17p loss is near-universal
in osteosarcoma, so overlap semantics would return `Both-single-SV` for almost every real
sample. v0.1 therefore classifies on breakpoints only. An SV that covers the locus without
placing a breakpoint inside it is reported as `Neither` with the `spans_locus` flag, never as
a both-genes call. `tests/test_classify.py::TestSpansLocusDivergence` pins this; if those
tests fail, the divergence has been undone and this section is wrong.

**2. Locus filtering moved into the library.** The reference's `architecture_call` counts
every record handed to it, leaving the locus filter to its caller, while its GRCh37 twin
counts `n_in_window` internally. A VCF is whole-genome, so filtering in the caller would make
every sample `Complex`. This library filters to the locus window before the ladder, matching
the reference's documented intent and its own reason strings ("in locus window").

**3. Both genes lie entirely inside the locus window.** `TP53_coding`
(7,668,402-7,687,550) and `WRAP53_coding` (7,686,071-7,703,502) are both contained in
7,560,000-7,710,000, so no part of either gene is gated out. Until v0.2.0 this section
recorded the opposite — WRAP53 running 15 kb past the window — which was an artefact of an
unlifted coordinate, not a property of the locus. See `CHANGELOG.md`.

**4. No coordinate conversion.** Feature bounds are compared directly against VCF `POS`.
VCF is 1-based and the reference compares its own inputs to the same bounds without
conversion. Introducing a 1-based-to-0-based shift here would move every boundary by one
base relative to the reference.

**5. `Complex` precedence.** Checked before any gene logic, as in the reference. A sample at
or above the threshold never receives a gene-level call even when both coding windows are hit.

**6. Contig names are normalised — a deliberate divergence, found by the equivalence
suite.** The reference compares contig names literally against `chr17`, so a file spelling the
contig `17` is treated as entirely off-locus and classified `Neither`. Both spellings are
common in real VCFs — Ensembl-style references drop the prefix, UCSC-style keep it — and a
standalone tool that silently returns `Neither` for half of them would be a trap. This library
accepts `chr17`, `17`, and case variants, all normalised to `chr17`. The divergence is pinned
by `tests/test_equivalence.py::TestDocumentedDivergences::test_contig_naming_divergence_is_real`
and recorded in the equivalence fixtures, where the reference's `Neither` and this library's
`TP53-only` are both stored for the same input. The reference was not wrong for its own
purpose — its inputs were a single pipeline with one spelling — but the assumption does not
survive being handed an arbitrary VCF.

## Testing posture

Fixture VCFs in `tests/data/` cover one case per category plus the boundary positions where a
call can flip. The Stage 1 equivalence suite pins expected calls derived from the reference
implementation at commit `a750daf`; those pinned values come from DepMap-based inference and
are **inferences, not a published gold standard**.

What the equivalence suite compares, precisely: the reference owns the category ladder, and
VCF parsing has no counterpart there, so the ladder is what gets pinned. Each case in
`scripts/generate_equivalence.py` lists its breakpoints explicitly, written by hand to describe
what the matching fixture VCF contains. Those breakpoints run through the reference's own
`annotate_breakpoints` and `architecture_call`; the test then parses the fixture VCF with this
library and compares categories. Specifying the two sides independently is what makes the test
meaningful — feeding the reference this library's own parser output would have hidden any
parsing bug from both sides.

Of 17 cases, 15 agree with the reference exactly. The two that do not are the divergences in
sections 1 and 6, each carrying a written explanation in the fixture file and a test asserting
it by name. A case is allowed to differ from the reference **only** if it carries that
explanation: `test_undocumented_cases_agree_with_the_reference_exactly` fails otherwise, which
is the guard against redefining behaviour without recording why.

## Change policy

Scope changes need a version bump and a `CHANGELOG.md` entry. Adding BED input, GRCh37
support, or segment semantics is a new minor version with its own scope note, not a patch.
