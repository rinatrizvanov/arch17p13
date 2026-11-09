"""The category ladder, its precedence, and the documented divergence."""

from __future__ import annotations

from pathlib import Path

import pytest

from arch17p13 import classify_records, classify_vcf
from arch17p13.classify import SPANS_LOCUS, Result
from arch17p13.regions import COMPLEX_THRESHOLD
from arch17p13.vcf import Breakpoint, SvRecord, read_records
from tests.conftest import write_vcf

GT = ("GT", "0/1")


def _record(record_id: str, *features: str, in_locus: bool = True) -> SvRecord:
    """Build a record that reports exactly ``features``, bypassing VCF parsing."""
    breakpoint_ = Breakpoint(
        chrom="chr17", pos=7_670_000, features=features, in_locus=in_locus, origin="POS"
    )
    return SvRecord(
        record_id=record_id,
        sv_type="BND",
        chrom="chr17",
        pos=7_670_000,
        end=None,
        filter_value="PASS",
        breakpoints=(breakpoint_,),
        spans_locus=False,
    )


class TestCategoriesFromFixtures:
    @pytest.mark.parametrize(
        ("fixture", "expected"),
        [
            ("neither_offlocus.vcf", "Neither"),
            ("neither_in_locus.vcf", "Neither"),
            ("tp53_only.vcf", "TP53-only"),
            ("tcab1_only.vcf", "TCAB1-only"),
            ("both_single_sv.vcf", "Both-single-SV"),
            ("both_independent.vcf", "Both-independent-SVs"),
            ("complex.vcf", "Complex"),
            ("svlen_del.vcf", "Both-single-SV"),
            ("insertion.vcf", "TCAB1-only"),
            ("chrom_naming.vcf", "TP53-only"),
            ("sites_only.vcf", "TP53-only"),
            ("spans_locus.vcf", "Neither"),
            ("wrap53_outside_window.vcf", "Neither"),
        ],
    )
    def test_fixture_category(self, data_dir: Path, fixture: str, expected: str) -> None:
        assert classify_vcf(data_dir / fixture).category == expected

    def test_reciprocal_pair_is_one_sv_not_two(self, data_dir: Path) -> None:
        """The MATEID collapse is what separates Both-single-SV from Both-independent."""
        single = classify_vcf(data_dir / "both_single_sv.vcf")
        independent = classify_vcf(data_dir / "both_independent.vcf")
        assert single.category == "Both-single-SV"
        assert single.n_locus_records == 1
        assert independent.category == "Both-independent-SVs"
        assert independent.n_locus_records == 2


class TestPrecedence:
    def test_complex_outranks_a_both_genes_call(self, data_dir: Path) -> None:
        """complex.vcf hits both coding windows; Complex must still win.

        This mirrors the reference, where the threshold check precedes all gene logic.
        """
        result = classify_vcf(data_dir / "complex.vcf")
        assert result.category == "Complex"
        assert result.n_locus_records >= COMPLEX_THRESHOLD

    def test_threshold_is_inclusive(self) -> None:
        records = [_record(f"r{i}", "TP53_5kb_flank") for i in range(COMPLEX_THRESHOLD)]
        assert classify_records(records).category == "Complex"

    def test_one_below_threshold_is_not_complex(self) -> None:
        records = [_record(f"r{i}", "TP53_5kb_flank") for i in range(COMPLEX_THRESHOLD - 1)]
        assert classify_records(records).category == "Neither"

    def test_both_single_outranks_both_independent(self) -> None:
        """A record hitting both genes wins even when other single-gene records exist."""
        records = [
            _record("pair", "TP53_coding", "WRAP53_coding"),
            _record("solo", "TP53_coding"),
        ]
        assert classify_records(records).category == "Both-single-SV"

    def test_non_coding_features_are_never_decisive(self) -> None:
        records = [_record("r1", "TP53_promoter", "WRAP53_promoter", "intergenic_TP53_WRAP53")]
        assert classify_records(records).category == "Neither"


class TestSpansLocusDivergence:
    """Pins the deliberate divergence from the reference's segment path.

    ``sv_architecture_auto`` in the reference is a forbidden variable: interval-overlap
    semantics called 33 of 35 samples ``Both-single-SV`` because one CN segment spanning
    the adjacent TP53 and WRAP53 windows overlaps both. This library classifies on
    breakpoints, so the same event is ``Neither`` plus a flag. If this test ever fails,
    the divergence has been undone and SCOPE.md is wrong.
    """

    def test_locus_spanning_deletion_is_not_a_both_genes_call(self, data_dir: Path) -> None:
        # The forbidden segment path would call this Both-single-SV, because one
        # interval overlaps both adjacent coding windows. Breakpoint semantics do not.
        result = classify_vcf(data_dir / "spans_locus.vcf")
        assert result.category == "Neither"

    def test_spanning_deletion_sets_the_flag(self, data_dir: Path) -> None:
        result = classify_vcf(data_dir / "spans_locus.vcf")
        assert SPANS_LOCUS in result.flags

    def test_the_flag_is_mentioned_in_the_reason(self, data_dir: Path) -> None:
        result = classify_vcf(data_dir / "spans_locus.vcf")
        assert "span" in result.reason.lower()

    def test_a_focal_event_does_not_set_the_flag(self, data_dir: Path) -> None:
        assert SPANS_LOCUS not in classify_vcf(data_dir / "tp53_only.vcf").flags

    def test_spanning_does_not_suppress_a_real_focal_call(self, tmp_path: Path) -> None:
        """A broad deletion alongside a focal breakend: the focal call still stands."""
        vcf = write_vcf(
            tmp_path / "mixed.vcf",
            [
                (
                    "chr17",
                    "7000000",
                    "broad",
                    "N",
                    "<DEL>",
                    ".",
                    "PASS",
                    "SVTYPE=DEL;END=8000000",
                    *GT,
                ),
                ("chr17", "7670000", "focal", "N", "N[chr5:1[", ".", "PASS", "SVTYPE=BND", *GT),
            ],
        )
        result = classify_vcf(vcf)
        assert result.category == "TP53-only"
        assert SPANS_LOCUS in result.flags
        assert result.n_locus_records == 1


class TestResultShape:
    def test_record_ids_are_sorted_and_deduplicated(self, data_dir: Path) -> None:
        result = classify_vcf(data_dir / "both_independent.vcf")
        assert result.record_ids == ("bnd_a", "bnd_c")

    def test_spanning_only_record_is_not_counted(self, data_dir: Path) -> None:
        result = classify_vcf(data_dir / "spans_locus.vcf")
        assert result.n_locus_records == 0
        assert result.record_ids == ()

    def test_reason_is_always_populated(self, data_dir: Path) -> None:
        # multisample.vcf is excluded on purpose: it refuses to classify without an
        # explicit sample, which TestSampleSelection covers.
        for fixture in sorted(data_dir.glob("*.vcf")):
            if fixture.name == "multisample.vcf":
                continue
            assert classify_vcf(fixture).reason, f"{fixture.name} produced an empty reason"

    def test_reason_names_the_deciding_records(self, data_dir: Path) -> None:
        assert "bnd_a" in classify_vcf(data_dir / "tp53_only.vcf").reason

    def test_as_dict_is_json_serialisable(self, data_dir: Path) -> None:
        import json

        payload = classify_vcf(data_dir / "both_single_sv.vcf").as_dict()
        assert json.loads(json.dumps(payload))["category"] == "Both-single-SV"

    def test_as_dict_sorts_flags(self, data_dir: Path) -> None:
        payload = classify_vcf(data_dir / "spans_locus.vcf").as_dict()
        assert payload["flags"] == [SPANS_LOCUS]

    def test_result_is_immutable(self) -> None:
        result = classify_records([])
        with pytest.raises(AttributeError):
            result.category = "Complex"  # type: ignore[misc]

    def test_empty_input_is_neither(self) -> None:
        result = classify_records([])
        assert result == Result("Neither", "No SVs in locus window", (), 0, frozenset())


class TestOptionsReachClassification:
    def test_pass_only_changes_the_category(self, data_dir: Path) -> None:
        strict = classify_vcf(data_dir / "filtered.vcf")
        lenient = classify_vcf(data_dir / "filtered.vcf", pass_only=False)
        assert strict.category == "TP53-only"
        assert lenient.category == "Both-independent-SVs"

    def test_sample_choice_changes_the_category(self, data_dir: Path) -> None:
        tumour = classify_vcf(data_dir / "multisample.vcf", sample="TUMOR")
        normal = classify_vcf(data_dir / "multisample.vcf", sample="NORMAL")
        assert tumour.category == "TCAB1-only"
        assert normal.category == "TP53-only"

    def test_classify_records_accepts_a_generator(self, data_dir: Path) -> None:
        records = read_records(data_dir / "tp53_only.vcf")
        assert classify_records(iter(records)).category == "TP53-only"
