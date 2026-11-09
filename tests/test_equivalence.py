"""Equivalence against the reference implementation's category ladder.

The pinned values in ``tests/data/stage1_equivalence.json`` were produced by running
the reference ``architecture_call`` at commit ``a750daf``; see
``scripts/generate_equivalence.py`` for how, and ``SCOPE.md`` for what the comparison
does and does not cover. The pins derive from DepMap-based inference and are not a
published gold standard.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from arch17p13 import classify_vcf

_PINS = json.loads((Path(__file__).parent / "data" / "stage1_equivalence.json").read_text())
_CASES: list[dict[str, Any]] = _PINS["cases"]

#: The reference commit the pins are valid for. Changing this needs a regeneration.
EXPECTED_REFERENCE_COMMIT = "a750daf1051be7eae4982cce2e427d4e973bd95f"


def _ids(cases: list[dict[str, Any]]) -> list[str]:
    return [str(case["name"]) for case in cases]


class TestProvenance:
    def test_pins_name_their_reference_commit(self) -> None:
        assert _PINS["provenance"]["reference_commit"] == EXPECTED_REFERENCE_COMMIT

    def test_pins_record_the_module_they_came_from(self) -> None:
        provenance = _PINS["provenance"]
        assert provenance["reference_module"].endswith("sv_annotator.py")
        assert len(provenance["reference_module_sha256"]) == 64

    def test_pins_are_labelled_as_inference_not_gold_standard(self) -> None:
        """Guards the claim the README and SCOPE.md both make about these values."""
        assert "not a published" in _PINS["provenance"]["description"]

    def test_every_case_names_an_existing_fixture(self, data_dir: Path) -> None:
        for case in _CASES:
            assert (data_dir / case["fixture"]).exists(), case["fixture"]

    def test_the_suite_is_not_empty(self) -> None:
        assert len(_CASES) >= 15


class TestLadderEquivalence:
    @pytest.mark.parametrize("case", _CASES, ids=_ids(_CASES))
    def test_library_matches_the_pinned_expectation(
        self, case: dict[str, Any], data_dir: Path
    ) -> None:
        result = classify_vcf(
            data_dir / case["fixture"],
            sample=case["options"]["sample"],
            pass_only=case["options"]["pass_only"],
        )
        assert result.category == case["expected_category"], (
            f"{case['name']}: expected {case['expected_category']}, got {result.category}. "
            f"Reference said {case['reference_category']}."
        )

    @pytest.mark.parametrize(
        "case",
        [c for c in _CASES if not c["divergence"]],
        ids=_ids([c for c in _CASES if not c["divergence"]]),
    )
    def test_undocumented_cases_agree_with_the_reference_exactly(
        self, case: dict[str, Any]
    ) -> None:
        """A case may only differ from the reference if it carries a divergence note.

        This is the guard against silent drift: changing a call without writing down
        why fails here rather than quietly redefining the library's behaviour.
        """
        assert case["expected_category"] == case["reference_category"], (
            f"{case['name']} differs from the reference with no divergence note. "
            "Either restore agreement or document it in SCOPE.md and the generator."
        )

    def test_most_cases_agree_with_the_reference(self) -> None:
        agreeing = [c for c in _CASES if c["expected_category"] == c["reference_category"]]
        assert len(agreeing) >= len(_CASES) - 2


class TestDocumentedDivergences:
    """Each divergence is asserted by name, so undoing one fails loudly."""

    def test_every_divergence_carries_an_explanation(self) -> None:
        for case in _CASES:
            if case["expected_category"] != case["reference_category"]:
                assert case["divergence"], f"{case['name']} diverges with no explanation"

    def test_every_divergence_note_points_at_scope_md(self) -> None:
        for case in _CASES:
            if case["divergence"]:
                assert "SCOPE.md" in case["divergence"], case["name"]

    def test_contig_naming_divergence_is_real(self, data_dir: Path) -> None:
        """The reference compares contigs literally; this library normalises them.

        A VCF spelling the contig ``17`` is entirely off-locus to the reference, which
        returns ``Neither``. Both spellings occur in real VCFs, so normalising is the
        defensible behaviour for a standalone tool — but it is a divergence, and this
        pins it.
        """
        case = next(c for c in _CASES if c["name"] == "contig_without_chr_prefix")
        assert case["reference_category"] == "Neither"
        assert case["expected_category"] == "TP53-only"
        assert classify_vcf(data_dir / case["fixture"]).category == "TP53-only"

    def test_chr_prefixed_equivalent_agrees_with_the_reference(self, data_dir: Path) -> None:
        """The same event spelled ``chr17`` is not a divergence at all."""
        assert classify_vcf(data_dir / "tp53_only.vcf").category == "TP53-only"

    def test_spanning_divergence_is_about_ingest_not_the_ladder(self) -> None:
        """Fed breakpoints, the reference ladder agrees; the divergence is upstream.

        The segment path is what would have called this ``Both-single-SV``, and that
        path is a forbidden variable in the Stage 1 pre-registration.
        """
        case = next(c for c in _CASES if c["name"] == "locus_spanning_deletion")
        assert case["reference_category"] == "Neither"
        assert case["expected_category"] == "Neither"
        assert "forbidden variable" in case["divergence"]
