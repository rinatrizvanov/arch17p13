"""Frozen GRCh38 constants for the 17p13.1 TP53/WRAP53 locus.

Every coordinate and threshold in this module is part of the v0.1 frozen scope.
Changing any value changes classifications, so changes require a version bump
and a ``CHANGELOG.md`` entry (see ``SCOPE.md``).

Coordinate convention
---------------------
Feature bounds are half-open, ``start <= pos < end``, and are compared directly
against the position values supplied by the caller. VCF ``POS`` is 1-based and is
compared without conversion, which is what the reference implementation does with
its own inputs. Introducing a 1-based-to-0-based shift here would silently move
every boundary by one base relative to the reference, so it is deliberately not done.

Why TP53 and WRAP53 cannot be treated as independent intervals
--------------------------------------------------------------
TP53 (minus strand) and WRAP53 (plus strand) are convergent overlapping genes sharing
a bidirectional promoter. Their coding windows **intersect** over 1,479 bp,
chr17:7,686,071-7,687,550, which is the top of TP53 intron 1. Two things follow, and
both matter:

* A *single* breakpoint inside that overlap legitimately disrupts both genes. This is
  what ``Both-single-SV`` exists to describe, and it is why that category is distinct
  from ``Both-independent-SVs``.
* Any interval merely *covering* the region overlaps both genes whether or not it
  breaks either. That is why this library classifies on breakpoint positions rather
  than interval overlap; see ``SCOPE.md``.

v0.1.0 shipped a WRAP53 interval that began one base after TP53 ends, so the overlap
was not representable and a single breakpoint could never hit both genes. See
``CHANGELOG.md`` and the reference implementation's
``docs/erratum_wrap53_coordinates.md``.
"""

from __future__ import annotations

from typing import Final

#: Reference assembly these coordinates are valid for. v0.1 supports this only.
ASSEMBLY: Final = "GRCh38"

#: Contig holding the locus, in normalised form (see :mod:`arch17p13.intervals`).
LOCUS_CHROM: Final = "chr17"

#: Locus window. Breakpoints outside it are ignored entirely.
LOCUS_START: Final = 7_560_000
LOCUS_END: Final = 7_710_000

#: Feature regions as ``(name, chrom, start, end)``, half-open on ``[start, end)``.
#: GRCh38 values: Ensembl gene bodies in GRCh37 lifted with the UCSC hg19-to-hg38
#: chain, which is a constant +96,682 across this locus.
FEATURES: Final[tuple[tuple[str, str, int, int], ...]] = (
    ("TP53_coding", "chr17", 7_668_402, 7_687_550),
    ("TP53_promoter", "chr17", 7_687_551, 7_692_550),  # 5 kb upstream of TSS (minus strand)
    ("TP53_5kb_flank", "chr17", 7_662_402, 7_668_401),  # 5 kb downstream of TES
    ("WRAP53_coding", "chr17", 7_686_071, 7_703_502),
    ("WRAP53_promoter", "chr17", 7_681_071, 7_686_070),  # 5 kb upstream of TSS (plus strand)
    ("WRAP53_5kb_flank", "chr17", 7_703_503, 7_708_502),  # 5 kb downstream of TES
    # The genes overlap, so there is no intergenic span between them to name.
    ("TP53_WRAP53_overlap", "chr17", 7_686_071, 7_687_550),
)

#: Features whose disruption drives the architecture call. Only these two matter to
#: the category ladder; the remaining features are reported but never decisive.
TP53_CODING: Final = "TP53_coding"
WRAP53_CODING: Final = "WRAP53_coding"

#: At or above this many distinct in-locus SV records, the sample is called
#: ``Complex`` and no gene-level call is attempted.
COMPLEX_THRESHOLD: Final = 5

#: Every category this library can emit, in ladder precedence order. ``unknown`` is
#: never returned by classification: it exists for callers recording "not assessed".
CATEGORIES: Final[tuple[str, ...]] = (
    "Complex",
    "Both-single-SV",
    "Both-independent-SVs",
    "TCAB1-only",
    "TP53-only",
    "Neither",
    "unknown",
)


def _coding_bounds(name: str) -> tuple[int, int]:
    """Start and end of a named feature, for the geometry invariant below."""
    return next((start, end) for feature, _chrom, start, end in FEATURES if feature == name)


def check_genes_overlap(tp53: tuple[int, int], wrap53: tuple[int, int]) -> None:
    """Raise unless the two coding windows intersect.

    Not decoration. If TP53 and WRAP53 stop overlapping, a single breakpoint can no
    longer hit both genes, ``Both-single-SV`` becomes unreachable from one breakpoint,
    and every co-disruption count in every downstream analysis silently reads zero
    while the test suite stays green. That is exactly what happened in v0.1.0.
    """
    if min(tp53[1], wrap53[1]) <= max(tp53[0], wrap53[0]):
        raise ValueError(
            f"TP53_coding {tp53[0]:,}-{tp53[1]:,} and WRAP53_coding "
            f"{wrap53[0]:,}-{wrap53[1]:,} do not intersect. They are convergent "
            "overlapping genes and must share at least one base."
        )


check_genes_overlap(_coding_bounds("TP53_coding"), _coding_bounds("WRAP53_coding"))

#: Width of the TP53/WRAP53 overlap, in bases. A single breakpoint in this span hits
#: both genes.
OVERLAP_WIDTH: Final = 1_479
