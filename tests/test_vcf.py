"""VCF parsing: breakpoint derivation, record identity, filtering, errors."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from arch17p13.vcf import VcfError, read_records
from tests.conftest import write_vcf

GT = ("GT", "0/1")


class TestBreakpointDerivation:
    def test_breakend_contributes_pos_and_mate(self, data_dir: Path) -> None:
        (record,) = read_records(data_dir / "tp53_only.vcf")
        assert record.sv_type == "BND"
        assert len(record.breakpoints) == 2
        origins = {bp.origin for bp in record.breakpoints}
        assert origins == {"POS", "MATE"}

    def test_mate_on_another_contig_is_off_locus(self, data_dir: Path) -> None:
        (record,) = read_records(data_dir / "tp53_only.vcf")
        mate = next(bp for bp in record.breakpoints if bp.origin == "MATE")
        assert mate.in_locus is False
        assert mate.features == ()

    def test_deletion_uses_svlen_when_end_absent(self, data_dir: Path) -> None:
        (record,) = read_records(data_dir / "svlen_del.vcf")
        assert record.end == 7_700_000
        assert {bp.origin for bp in record.breakpoints} == {"POS", "END"}

    def test_insertion_has_no_span(self, data_dir: Path) -> None:
        """INS carries SVLEN but describes no interval, so it gets one breakpoint."""
        (record,) = read_records(data_dir / "insertion.vcf")
        assert record.end is None
        assert len(record.breakpoints) == 1

    def test_both_endpoints_of_a_focal_deletion_are_annotated(self, data_dir: Path) -> None:
        (record,) = read_records(data_dir / "svlen_del.vcf")
        assert record.hits_tp53_coding
        assert record.hits_wrap53_coding

    @pytest.mark.parametrize(
        "alt",
        ["N[chr17:7700000[", "N]chr17:7700000]", "[chr17:7700000[N", "]chr17:7700000]N"],
    )
    def test_all_four_bracket_orientations_parse(self, tmp_path: Path, alt: str) -> None:
        vcf = write_vcf(
            tmp_path / "bnd.vcf",
            [("chr17", "7670000", "b1", "N", alt, ".", "PASS", "SVTYPE=BND", *GT)],
        )
        (record,) = read_records(vcf)
        mate = next(bp for bp in record.breakpoints if bp.origin == "MATE")
        assert mate.pos == 7_700_000
        assert mate.in_locus


class TestRecordIdentity:
    def test_mateid_collapses_a_reciprocal_pair(self, data_dir: Path) -> None:
        records = read_records(data_dir / "both_single_sv.vcf")
        assert len({r.record_id for r in records}) == 1

    def test_unlinked_breakends_stay_separate(self, data_dir: Path) -> None:
        records = read_records(data_dir / "both_independent.vcf")
        assert len({r.record_id for r in records}) == 2

    def test_missing_id_is_synthesised_from_line_number(self, tmp_path: Path) -> None:
        vcf = write_vcf(
            tmp_path / "noid.vcf",
            [("chr17", "7670000", ".", "N", "N[chr5:1[", ".", "PASS", "SVTYPE=BND", *GT)],
        )
        (record,) = read_records(vcf)
        assert record.record_id.startswith("__line")

    def test_mateid_list_is_sorted_into_one_id(self, tmp_path: Path) -> None:
        vcf = write_vcf(
            tmp_path / "multimate.vcf",
            [
                (
                    "chr17",
                    "7670000",
                    "zz",
                    "N",
                    "N[chr5:1[",
                    ".",
                    "PASS",
                    "SVTYPE=BND;MATEID=aa,mm",
                    *GT,
                )
            ],
        )
        (record,) = read_records(vcf)
        assert record.record_id == "aa|mm|zz"

    def test_mateid_is_ignored_for_non_breakends(self, tmp_path: Path) -> None:
        vcf = write_vcf(
            tmp_path / "delmate.vcf",
            [
                (
                    "chr17",
                    "7670000",
                    "d1",
                    "N",
                    "<DEL>",
                    ".",
                    "PASS",
                    "SVTYPE=DEL;END=7700000;MATEID=d2",
                    *GT,
                )
            ],
        )
        (record,) = read_records(vcf)
        assert record.record_id == "d1"


class TestLocusRestriction:
    def test_off_locus_records_are_dropped(self, data_dir: Path) -> None:
        assert read_records(data_dir / "neither_offlocus.vcf") == []

    def test_contig_without_chr_prefix_is_accepted(self, data_dir: Path) -> None:
        (record,) = read_records(data_dir / "chrom_naming.vcf")
        assert record.hits_tp53_coding

    def test_breakpoint_past_the_locus_window_is_dropped(self, data_dir: Path) -> None:
        assert read_records(data_dir / "wrap53_outside_window.vcf") == []

    def test_spanning_record_is_kept_and_flagged(self, data_dir: Path) -> None:
        (record,) = read_records(data_dir / "spans_locus.vcf")
        assert record.spans_locus is True
        assert record.in_locus is False
        assert record.features == frozenset()

    def test_focal_record_is_not_flagged_as_spanning(self, data_dir: Path) -> None:
        (record,) = read_records(data_dir / "svlen_del.vcf")
        assert record.spans_locus is False


class TestFilterHandling:
    def test_pass_only_is_the_default(self, data_dir: Path) -> None:
        records = read_records(data_dir / "filtered.vcf")
        assert [r.record_id for r in records] == ["passing"]

    def test_all_filters_can_be_requested(self, data_dir: Path) -> None:
        records = read_records(data_dir / "filtered.vcf", pass_only=False)
        assert {r.record_id for r in records} == {"passing", "lowq"}

    def test_dot_filter_counts_as_passing(self, tmp_path: Path) -> None:
        vcf = write_vcf(
            tmp_path / "dot.vcf",
            [("chr17", "7670000", "d", "N", "N[chr5:1[", ".", ".", "SVTYPE=BND", *GT)],
        )
        assert len(read_records(vcf)) == 1


class TestSampleSelection:
    def test_multisample_requires_a_sample(self, data_dir: Path) -> None:
        with pytest.raises(VcfError, match="pass sample="):
            read_records(data_dir / "multisample.vcf")

    def test_unknown_sample_is_rejected(self, data_dir: Path) -> None:
        with pytest.raises(VcfError, match="not in"):
            read_records(data_dir / "multisample.vcf", sample="MISSING")

    def test_genotype_gates_inclusion(self, data_dir: Path) -> None:
        tumour = read_records(data_dir / "multisample.vcf", sample="TUMOR")
        assert [r.record_id for r in tumour] == ["bnd_t"]

    def test_the_other_sample_sees_the_other_record(self, data_dir: Path) -> None:
        normal = read_records(data_dir / "multisample.vcf", sample="NORMAL")
        assert [r.record_id for r in normal] == ["bnd_n"]

    def test_sites_only_vcf_needs_no_sample(self, data_dir: Path) -> None:
        assert len(read_records(data_dir / "sites_only.vcf")) == 1

    @pytest.mark.parametrize("gt", ["0/0", "./.", ".", "0|0"])
    def test_non_carrier_genotypes_are_skipped(self, tmp_path: Path, gt: str) -> None:
        vcf = write_vcf(
            tmp_path / "gt.vcf",
            [("chr17", "7670000", "g", "N", "N[chr5:1[", ".", "PASS", "SVTYPE=BND", "GT", gt)],
            samples=("ONLY", "OTHER"),
        )
        assert read_records(vcf, sample="ONLY") == []

    @pytest.mark.parametrize("gt", ["0/1", "1|0", "1/1", "0/1:30"])
    def test_carrier_genotypes_are_kept(self, tmp_path: Path, gt: str) -> None:
        vcf = write_vcf(
            tmp_path / "gt.vcf",
            [("chr17", "7670000", "g", "N", "N[chr5:1[", ".", "PASS", "SVTYPE=BND", "GT", gt)],
            samples=("ONLY", "OTHER"),
        )
        assert len(read_records(vcf, sample="ONLY")) == 1


class TestCompressionAndErrors:
    def test_gzipped_vcf_reads_identically(self, data_dir: Path, tmp_path: Path) -> None:
        plain = (data_dir / "tp53_only.vcf").read_bytes()
        gz = tmp_path / "tp53_only.vcf.gz"
        with gzip.open(gz, "wb") as handle:
            handle.write(plain)
        assert read_records(gz) == read_records(data_dir / "tp53_only.vcf")

    def test_missing_file_raises(self, tmp_path: Path) -> None:
        with pytest.raises(VcfError, match="not found"):
            read_records(tmp_path / "absent.vcf")

    def test_truncated_line_raises(self, tmp_path: Path) -> None:
        vcf = write_vcf(tmp_path / "short.vcf", [("chr17", "7670000", "x", "N")])
        with pytest.raises(VcfError, match="at least 8"):
            read_records(vcf)

    def test_non_integer_pos_raises(self, tmp_path: Path) -> None:
        vcf = write_vcf(
            tmp_path / "badpos.vcf",
            [("chr17", "seven", "x", "N", "N[chr5:1[", ".", "PASS", "SVTYPE=BND", *GT)],
        )
        with pytest.raises(VcfError, match="POS is not an integer"):
            read_records(vcf)

    def test_non_integer_end_raises(self, tmp_path: Path) -> None:
        vcf = write_vcf(
            tmp_path / "badend.vcf",
            [("chr17", "7670000", "x", "N", "<DEL>", ".", "PASS", "SVTYPE=DEL;END=lots", *GT)],
        )
        with pytest.raises(VcfError, match="END is not an integer"):
            read_records(vcf)

    def test_non_integer_svlen_raises(self, tmp_path: Path) -> None:
        vcf = write_vcf(
            tmp_path / "badlen.vcf",
            [("chr17", "7670000", "x", "N", "<DEL>", ".", "PASS", "SVTYPE=DEL;SVLEN=big", *GT)],
        )
        with pytest.raises(VcfError, match="SVLEN is not an integer"):
            read_records(vcf)

    def test_svtype_inferred_from_symbolic_alt(self, tmp_path: Path) -> None:
        vcf = write_vcf(
            tmp_path / "nosvtype.vcf",
            [("chr17", "7670000", "x", "N", "<DEL>", ".", "PASS", "END=7700000", *GT)],
        )
        (record,) = read_records(vcf)
        assert record.sv_type == "DEL"

    def test_blank_lines_are_tolerated(self, tmp_path: Path) -> None:
        vcf = tmp_path / "blanks.vcf"
        write_vcf(
            vcf, [("chr17", "7670000", "x", "N", "N[chr5:1[", ".", "PASS", "SVTYPE=BND", *GT)]
        )
        vcf.write_text(vcf.read_text().replace("#CHROM", "\n#CHROM") + "\n\n")
        assert len(read_records(vcf)) == 1


class TestParserEdgeCases:
    """Coverage for the parser's fallbacks, each of which is reachable from real VCFs."""

    def test_locus_breakpoints_property_filters(self, data_dir: Path) -> None:
        (record,) = read_records(data_dir / "tp53_only.vcf")
        assert len(record.breakpoints) == 2
        assert len(record.locus_breakpoints) == 1
        assert record.locus_breakpoints[0].origin == "POS"

    def test_empty_info_column(self, tmp_path: Path) -> None:
        vcf = write_vcf(
            tmp_path / "noinfo.vcf",
            [("chr17", "7670000", "x", "N", "N[chr5:1[", ".", "PASS", ".", *GT)],
        )
        (record,) = read_records(vcf)
        assert record.sv_type == "BND"
        assert record.end is None

    def test_info_with_empty_entries(self, tmp_path: Path) -> None:
        vcf = write_vcf(
            tmp_path / "oddinfo.vcf",
            [
                (
                    "chr17",
                    "7670000",
                    "x",
                    "N",
                    "N[chr5:1[",
                    ".",
                    "PASS",
                    "SVTYPE=BND;;IMPRECISE",
                    *GT,
                )
            ],
        )
        (record,) = read_records(vcf)
        assert record.sv_type == "BND"

    def test_breakend_inferred_from_brackets_without_svtype(self, tmp_path: Path) -> None:
        vcf = write_vcf(
            tmp_path / "nosvtype_bnd.vcf",
            [("chr17", "7670000", "x", "N", "N[chr17:7700000[", ".", "PASS", "IMPRECISE", *GT)],
        )
        (record,) = read_records(vcf)
        assert record.sv_type == "BND"
        assert len(record.breakpoints) == 2

    def test_unrecognised_alt_yields_empty_svtype_and_one_breakpoint(self, tmp_path: Path) -> None:
        vcf = write_vcf(
            tmp_path / "plain.vcf",
            [("chr17", "7670000", "x", "N", "A", ".", "PASS", ".", *GT)],
        )
        (record,) = read_records(vcf)
        assert record.sv_type == ""
        assert len(record.breakpoints) == 1

    def test_breakend_with_no_parsable_mate(self, tmp_path: Path) -> None:
        """SVTYPE=BND but the ALT carries no bracket: POS alone is used."""
        vcf = write_vcf(
            tmp_path / "bndnomate.vcf",
            [("chr17", "7670000", "x", "N", "N", ".", "PASS", "SVTYPE=BND", *GT)],
        )
        (record,) = read_records(vcf)
        assert len(record.breakpoints) == 1

    def test_empty_mateid_falls_back_to_own_id(self, tmp_path: Path) -> None:
        vcf = write_vcf(
            tmp_path / "emptymate.vcf",
            [
                (
                    "chr17",
                    "7670000",
                    "solo",
                    "N",
                    "N[chr5:1[",
                    ".",
                    "PASS",
                    "SVTYPE=BND;MATEID=,",
                    *GT,
                )
            ],
        )
        (record,) = read_records(vcf)
        assert record.record_id == "solo"

    def test_missing_genotype_column_raises(self, tmp_path: Path) -> None:
        """Header promises two samples but the data line carries one."""
        vcf = write_vcf(
            tmp_path / "shortgt.vcf",
            [("chr17", "7670000", "x", "N", "N[chr5:1[", ".", "PASS", "SVTYPE=BND", "GT", "0/1")],
            samples=("FIRST", "SECOND"),
        )
        with pytest.raises(VcfError, match="no genotype column"):
            read_records(vcf, sample="SECOND")

    def test_mateid_pair_keeps_both_rows_for_grouping(self, data_dir: Path) -> None:
        """Both rows survive parsing; collapsing to one SV happens in classification."""
        records = read_records(data_dir / "both_single_sv.vcf")
        assert len(records) == 2
        assert len({r.record_id for r in records}) == 1
