"""Invariants on the frozen constants.

These tests exist to make an accidental coordinate edit fail loudly: every value in
``regions`` changes classifications, so ``SCOPE.md`` requires a version bump to touch
one. A test failing here means either a real mistake or a deliberate scope change
that still needs its CHANGELOG entry.
"""

from __future__ import annotations

from arch17p13.regions import (
    ASSEMBLY,
    CATEGORIES,
    COMPLEX_THRESHOLD,
    FEATURES,
    LOCUS_CHROM,
    LOCUS_END,
    LOCUS_START,
    TP53_CODING,
    WRAP53_CODING,
)


def test_assembly_is_grch38() -> None:
    assert ASSEMBLY == "GRCh38"


def test_locus_window_is_ordered_and_nonempty() -> None:
    assert LOCUS_START < LOCUS_END
    assert LOCUS_CHROM == "chr17"


def test_complex_threshold_matches_reference() -> None:
    assert COMPLEX_THRESHOLD == 5


def test_feature_names_are_unique() -> None:
    names = [name for name, _, _, _ in FEATURES]
    assert len(names) == len(set(names))


def test_every_feature_interval_is_ordered() -> None:
    for name, _chrom, start, end in FEATURES:
        assert start < end, f"{name} has a non-positive width"


def test_both_decisive_features_are_declared() -> None:
    names = {name for name, _, _, _ in FEATURES}
    assert TP53_CODING in names
    assert WRAP53_CODING in names


def test_coding_windows_are_adjacent_not_overlapping() -> None:
    """The geometry behind the segment-path artefact documented in SCOPE.md.

    ``TP53_coding`` ends exactly one base below where ``WRAP53_coding`` begins, so any
    interval covering the junction overlaps both genes. This adjacency is the reason
    v0.1 classifies on breakpoints rather than on interval overlap.
    """
    lookup = {name: (start, end) for name, _, start, end in FEATURES}
    tp53_start, tp53_end = lookup[TP53_CODING]
    wrap53_start, wrap53_end = lookup[WRAP53_CODING]
    assert tp53_end < wrap53_start
    assert wrap53_start - tp53_end == 1
    assert tp53_start < tp53_end < wrap53_start < wrap53_end


def test_wrap53_coding_extends_past_the_locus_window() -> None:
    """A real asymmetry the reference has, preserved deliberately.

    ``WRAP53_coding`` runs to 7,725,564 but the locus window stops at 7,710,000, so a
    breakpoint in the tail of WRAP53 is gated out and reports no features at all. The
    reference behaves the same way because it only annotates in-locus positions.
    """
    lookup = {name: (start, end) for name, _, start, end in FEATURES}
    _, wrap53_end = lookup[WRAP53_CODING]
    assert wrap53_end > LOCUS_END


def test_categories_cover_the_ladder_and_the_sentinel() -> None:
    assert set(CATEGORIES) == {
        "Complex",
        "Both-single-SV",
        "Both-independent-SVs",
        "TCAB1-only",
        "TP53-only",
        "Neither",
        "unknown",
    }
