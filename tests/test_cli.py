"""CLI behaviour: output shapes, option plumbing and exit status."""

from __future__ import annotations

import io
import json
from pathlib import Path

import pytest

from arch17p13 import __version__
from arch17p13.cli import main
from tests.conftest import write_vcf

GT = ("GT", "0/1")


def run(*argv: str) -> tuple[int, str, str]:
    """Invoke the CLI with injected streams; returns (status, stdout, stderr)."""
    out, err = io.StringIO(), io.StringIO()
    status = main(list(argv), stdout=out, stderr=err)
    return status, out.getvalue(), err.getvalue()


class TestTextOutput:
    def test_default_is_one_tab_separated_line(self, data_dir: Path) -> None:
        status, out, err = run("classify", str(data_dir / "tp53_only.vcf"))
        assert status == 0
        assert out.strip() == "tp53_only.vcf\tTP53-only"
        assert err == ""

    def test_several_files_give_several_lines(self, data_dir: Path) -> None:
        status, out, _ = run(
            "classify", str(data_dir / "tp53_only.vcf"), str(data_dir / "tcab1_only.vcf")
        )
        assert status == 0
        assert out.strip().splitlines() == [
            "tp53_only.vcf\tTP53-only",
            "tcab1_only.vcf\tTCAB1-only",
        ]

    def test_explain_adds_reason_and_records(self, data_dir: Path) -> None:
        _, out, _ = run("classify", str(data_dir / "both_single_sv.vcf"), "--explain")
        assert "reason:" in out
        assert "records: bnd_a|bnd_b (1 in locus)" in out

    def test_explain_reports_flags(self, data_dir: Path) -> None:
        _, out, _ = run("classify", str(data_dir / "spans_locus.vcf"), "--explain")
        assert "flags: spans_locus" in out

    def test_explain_states_when_no_records_are_in_locus(self, data_dir: Path) -> None:
        _, out, _ = run("classify", str(data_dir / "spans_locus.vcf"), "--explain")
        assert "records: none in locus" in out

    def test_no_flag_line_when_there_are_no_flags(self, data_dir: Path) -> None:
        _, out, _ = run("classify", str(data_dir / "tp53_only.vcf"), "--explain")
        assert "flags:" not in out


class TestJsonOutput:
    def test_single_file_emits_one_object(self, data_dir: Path) -> None:
        _, out, _ = run("classify", str(data_dir / "tp53_only.vcf"), "--json")
        payload = json.loads(out)
        assert isinstance(payload, dict)
        assert payload["category"] == "TP53-only"
        assert payload["file"].endswith("tp53_only.vcf")

    def test_several_files_emit_an_array(self, data_dir: Path) -> None:
        _, out, _ = run(
            "classify", str(data_dir / "tp53_only.vcf"), str(data_dir / "complex.vcf"), "--json"
        )
        payload = json.loads(out)
        assert isinstance(payload, list)
        assert [entry["category"] for entry in payload] == ["TP53-only", "Complex"]

    def test_json_carries_records_and_flags(self, data_dir: Path) -> None:
        _, out, _ = run("classify", str(data_dir / "spans_locus.vcf"), "--json")
        payload = json.loads(out)
        assert payload["flags"] == ["spans_locus"]
        assert payload["n_locus_records"] == 0

    def test_json_is_not_ascii_escaped(self, data_dir: Path) -> None:
        """The reason strings contain an em dash; keep it readable."""
        _, out, _ = run("classify", str(data_dir / "complex.vcf"), "--json")
        assert "\\u2014" not in out
        assert "—" in out


class TestOptionPlumbing:
    def test_sample_selects_the_genotype_column(self, data_dir: Path) -> None:
        _, tumour, _ = run("classify", str(data_dir / "multisample.vcf"), "--sample", "TUMOR")
        _, normal, _ = run("classify", str(data_dir / "multisample.vcf"), "--sample", "NORMAL")
        assert "TCAB1-only" in tumour
        assert "TP53-only" in normal

    def test_all_filters_changes_the_call(self, data_dir: Path) -> None:
        _, strict, _ = run("classify", str(data_dir / "filtered.vcf"))
        _, lenient, _ = run("classify", str(data_dir / "filtered.vcf"), "--all-filters")
        assert "TP53-only" in strict
        assert "Both-independent-SVs" in lenient


class TestExitStatus:
    def test_neither_is_not_an_error(self, data_dir: Path) -> None:
        """A category is a result. Only I/O and parsing problems are failures."""
        status, out, _ = run("classify", str(data_dir / "neither_offlocus.vcf"))
        assert status == 0
        assert "Neither" in out

    def test_missing_file_fails_with_status_one(self, tmp_path: Path) -> None:
        status, out, err = run("classify", str(tmp_path / "absent.vcf"))
        assert status == 1
        assert "not found" in err
        assert out == ""

    def test_multisample_without_sample_fails(self, data_dir: Path) -> None:
        status, _, err = run("classify", str(data_dir / "multisample.vcf"))
        assert status == 1
        assert "pass sample=" in err

    def test_one_bad_file_does_not_stop_the_others(self, data_dir: Path, tmp_path: Path) -> None:
        status, out, err = run(
            "classify", str(tmp_path / "absent.vcf"), str(data_dir / "tp53_only.vcf")
        )
        assert status == 1
        assert "TP53-only" in out
        assert "not found" in err

    def test_malformed_vcf_fails_cleanly(self, tmp_path: Path) -> None:
        vcf = write_vcf(
            tmp_path / "bad.vcf",
            [("chr17", "nope", "x", "N", "N[chr5:1[", ".", "PASS", "SVTYPE=BND", *GT)],
        )
        status, _, err = run("classify", str(vcf))
        assert status == 1
        assert "POS is not an integer" in err


class TestUsage:
    def test_version_flag(self) -> None:
        with pytest.raises(SystemExit) as excinfo:
            run("--version")
        assert excinfo.value.code == 0

    def test_subcommand_is_required(self) -> None:
        with pytest.raises(SystemExit) as excinfo:
            run()
        assert excinfo.value.code == 2

    def test_classify_requires_a_file(self) -> None:
        with pytest.raises(SystemExit) as excinfo:
            run("classify")
        assert excinfo.value.code == 2

    def test_version_matches_the_package(self) -> None:
        assert __version__ == "0.2.0"

    def test_version_has_a_changelog_entry(self) -> None:
        """A release without an entry is the one nobody can reconstruct later."""
        changelog = Path(__file__).resolve().parent.parent / "CHANGELOG.md"
        assert f"## [{__version__}]" in changelog.read_text()


class TestModuleExecution:
    def test_python_m_arch17p13_works(self, data_dir: Path) -> None:
        """``python -m arch17p13`` is a supported entry point alongside the script."""
        import subprocess
        import sys

        completed = subprocess.run(
            [sys.executable, "-m", "arch17p13", "classify", str(data_dir / "tp53_only.vcf")],
            capture_output=True,
            text=True,
            check=False,
        )
        assert completed.returncode == 0
        assert "TP53-only" in completed.stdout
