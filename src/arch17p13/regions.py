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
``TP53_coding`` ends at 7,687,550 and ``WRAP53_coding`` begins at 7,687,551: the two
coding windows are immediately adjacent, sharing a head-to-head bidirectional
promoter. Any interval wider than roughly 57 kb that covers this region therefore
overlaps both genes. That geometry is the reason this library classifies on
breakpoint positions rather than on interval overlap; see ``SCOPE.md``.
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
#: Approximate GRCh38 values derived from Ensembl/GENCODE v45.
FEATURES: Final[tuple[tuple[str, str, int, int], ...]] = (
    ("TP53_coding", "chr17", 7_668_402, 7_687_550),
    ("TP53_promoter", "chr17", 7_687_551, 7_692_550),  # 5 kb upstream of TSS (minus strand)
    ("TP53_5kb_flank", "chr17", 7_662_402, 7_668_401),  # 5 kb downstream of TES
    ("WRAP53_coding", "chr17", 7_687_551, 7_725_564),
    ("WRAP53_promoter", "chr17", 7_682_551, 7_687_550),  # 5 kb upstream of TSS (plus strand)
    ("WRAP53_5kb_flank", "chr17", 7_725_565, 7_730_564),
    ("intergenic_TP53_WRAP53", "chr17", 7_687_000, 7_688_000),
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
