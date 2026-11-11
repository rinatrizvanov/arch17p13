"""Invariants on the frozen constants.

These tests exist to make an accidental coordinate edit fail loudly: every value in
``regions`` changes classifications, so ``SCOPE.md`` requires a version bump to touch
one. A test failing here means either a real mistake or a deliberate scope change
that still needs its CHANGELOG entry.
"""

from __future__ import annotations

import pytest

from arch17p13.regions import (
    ASSEMBLY,
    CATEGORIES,
    COMPLEX_THRESHOLD,
    FEATURES,
    LOCUS_CHROM,
    LOCUS_END,
    LOCUS_START,
    OVERLAP_WIDTH,
    TP53_CODING,
    WRAP53_CODING,
    check_genes_overlap,
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


def test_coding_windows_overlap() -> None:
    """The central geometric fact of this locus.

    TP53 and WRAP53 are convergent overlapping genes. Their bodies intersect over
    1,479 bp at the top of TP53 intron 1, so one breakpoint there disrupts both. If
    this ever stops being true, ``Both-single-SV`` becomes unreachable from a single
    breakpoint and co-disruption counts silently read zero. v0.1.0 shipped exactly
    that bug; see CHANGELOG.md.
    """
    lookup = {name: (start, end) for name, _, start, end in FEATURES}
    tp53_start, tp53_end = lookup[TP53_CODING]
    wrap53_start, wrap53_end = lookup[WRAP53_CODING]
    assert wrap53_start < tp53_end, "the genes must overlap, not merely abut"
    assert tp53_start < wrap53_start
    assert tp53_end < wrap53_end
    assert min(tp53_end, wrap53_end) - max(tp53_start, wrap53_start) == OVERLAP_WIDTH


def test_overlap_feature_matches_the_intersection() -> None:
    """``TP53_WRAP53_overlap`` must be exactly the intersection, not an approximation."""
    lookup = {name: (start, end) for name, _, start, end in FEATURES}
    tp53 = lookup[TP53_CODING]
    wrap53 = lookup[WRAP53_CODING]
    assert lookup["TP53_WRAP53_overlap"] == (
        max(tp53[0], wrap53[0]),
        min(tp53[1], wrap53[1]),
    )


def test_check_genes_overlap_rejects_disjoint_windows() -> None:
    """The invariant has to actually fire; v0.1.0's table is the regression case."""
    with pytest.raises(ValueError, match="do not intersect"):
        check_genes_overlap((7_668_402, 7_687_550), (7_687_551, 7_725_564))


def test_check_genes_overlap_accepts_the_shipped_table() -> None:
    lookup = {name: (start, end) for name, _, start, end in FEATURES}
    check_genes_overlap(lookup[TP53_CODING], lookup[WRAP53_CODING])


def test_both_genes_fit_inside_the_locus_window() -> None:
    """With WRAP53 lifted, the window contains both genes end to end.

    Under the v0.1.0 table WRAP53 ran 15 kb past the window and its tail was gated
    out. That asymmetry was an artefact of the unlifted coordinate, not a property
    of the locus.
    """
    lookup = {name: (start, end) for name, _, start, end in FEATURES}
    for feature in (TP53_CODING, WRAP53_CODING):
        start, end = lookup[feature]
        assert LOCUS_START <= start < end <= LOCUS_END, f"{feature} escapes the window"


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
