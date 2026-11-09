"""Position-to-feature lookups, including every boundary that changes a call."""

from __future__ import annotations

import pytest

from arch17p13.intervals import features_for_position, normalise_chrom, position_in_locus
from arch17p13.regions import LOCUS_CHROM, LOCUS_END, LOCUS_START


class TestNormaliseChrom:
    @pytest.mark.parametrize("value", ["chr17", "17", "Chr17", "CHR17", " chr17 "])
    def test_locus_contig_spellings(self, value: str) -> None:
        assert normalise_chrom(value) == LOCUS_CHROM

    @pytest.mark.parametrize("value", ["chr1", "1", "chrX", "17_random", "", "chr170"])
    def test_other_contigs_are_none(self, value: str) -> None:
        assert normalise_chrom(value) is None

    def test_none_input(self) -> None:
        assert normalise_chrom(None) is None


class TestPositionInLocus:
    def test_start_is_inclusive(self) -> None:
        assert position_in_locus("chr17", LOCUS_START)

    def test_end_is_exclusive(self) -> None:
        assert not position_in_locus("chr17", LOCUS_END)
        assert position_in_locus("chr17", LOCUS_END - 1)

    def test_before_window(self) -> None:
        assert not position_in_locus("chr17", LOCUS_START - 1)

    def test_wrong_contig(self) -> None:
        assert not position_in_locus("chr1", 7_670_000)


class TestFeaturesForPosition:
    """The boundary table. Each row here is a position where a call can flip."""

    @pytest.mark.parametrize(
        ("pos", "expected"),
        [
            (7_600_000, ()),
            (7_665_000, ("TP53_5kb_flank",)),
            (7_668_402, ("TP53_coding",)),
            (7_670_000, ("TP53_coding",)),
            (7_687_549, ("TP53_coding", "WRAP53_promoter", "intergenic_TP53_WRAP53")),
            (7_687_550, ("intergenic_TP53_WRAP53",)),
            (7_687_551, ("TP53_promoter", "WRAP53_coding", "intergenic_TP53_WRAP53")),
            (7_700_000, ("WRAP53_coding",)),
        ],
    )
    def test_boundary_table(self, pos: int, expected: tuple[str, ...]) -> None:
        assert features_for_position("chr17", pos) == expected

    def test_last_base_of_tp53_coding_is_still_tp53(self) -> None:
        assert "TP53_coding" in features_for_position("chr17", 7_687_549)

    def test_first_base_of_wrap53_coding_is_not_tp53(self) -> None:
        assert "TP53_coding" not in features_for_position("chr17", 7_687_551)

    def test_single_base_gap_between_genes_hits_neither(self) -> None:
        """7,687,550 belongs to no coding window: the one base separating the genes."""
        features = features_for_position("chr17", 7_687_550)
        assert "TP53_coding" not in features
        assert "WRAP53_coding" not in features

    def test_inside_wrap53_but_outside_locus_window_reports_nothing(self) -> None:
        assert position_in_locus("chr17", 7_715_000) is False
        assert features_for_position("chr17", 7_715_000) == ()

    def test_off_contig_reports_nothing(self) -> None:
        assert features_for_position("chr5", 7_670_000) == ()

    def test_results_follow_declaration_order(self) -> None:
        """Feature names come back in FEATURES order, so output is deterministic."""
        assert features_for_position("chr17", 7_687_551) == (
            "TP53_promoter",
            "WRAP53_coding",
            "intergenic_TP53_WRAP53",
        )
