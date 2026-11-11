#!/usr/bin/env python3
"""Regenerate ``tests/data/stage1_equivalence.json`` from the reference implementation.

Not run by the test suite and not part of the package: it needs a checkout of
`osteosarcoma-alt-dependency-research <https://github.com/rinatrizvanov/osteosarcoma-alt-dependency-research>`_
and that project's pandas environment. The suite consumes the JSON it writes.

    uv run --with pandas python scripts/generate_equivalence.py \
        --reference-repo ../osteosarcoma-alt-dependency-research

What is actually being compared
-------------------------------
The reference owns the *category ladder*; VCF parsing has no counterpart there, so
this pins the ladder. Each case lists its breakpoints explicitly, written by hand to
describe what the matching fixture VCF contains. Those breakpoints go through the
reference's own ``annotate_breakpoints`` and ``architecture_call``, and the resulting
category is pinned. The test then parses the fixture VCF with *this* library and
compares. Because the two sides are specified independently, a parsing bug in this
library shows up as a mismatch instead of being fed to both sides.

Breakpoints are filtered to the locus window before the ladder runs, matching the
decision recorded in ``SCOPE.md`` (the reference leaves that filter to its caller).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

#: The reference commit these pins are valid for.
EXPECTED_COMMIT = "9eb36dbc66b66d6a73a2736b3df0f91b969d1b7f"

#: Path of the reference module inside that repository.
REFERENCE_MODULE = Path("part1_alt_classification/src/sv_annotator.py")

Breakpoints = list[tuple[str, str, int]]  # (record_id, chrom, pos)


def _case(
    name: str,
    fixture: str,
    breakpoints: Breakpoints,
    *,
    sample: str | None = None,
    pass_only: bool = True,
    expected: str | None = None,
    divergence: str | None = None,
) -> dict[str, Any]:
    """One equivalence case. ``expected`` overrides the reference only for a divergence."""
    return {
        "name": name,
        "fixture": fixture,
        "options": {"sample": sample, "pass_only": pass_only},
        "breakpoints": [list(bp) for bp in breakpoints],
        "expected_category_override": expected,
        "divergence": divergence,
    }


# Mate breakpoints on other contigs are included because the reference sees them too:
# they are rows in the annotated frame that contribute no features.
CASES: list[dict[str, Any]] = [
    _case(
        "off_locus_only",
        "neither_offlocus.vcf",
        [
            ("sv1", "chr1", 1_000_000),
            ("sv1", "chr1", 1_050_000),
            ("sv2", "chr17", 2_000_000),
            ("sv2", "chr17", 2_010_000),
        ],
    ),
    _case(
        "tp53_coding_single_breakend",
        "tp53_only.vcf",
        [
            ("bnd_a", "chr17", 7_670_000),
            ("bnd_a", "chr5", 100_000),
        ],
    ),
    _case(
        "wrap53_coding_single_breakend",
        "tcab1_only.vcf",
        [
            ("bnd_a", "chr17", 7_700_000),
            ("bnd_a", "chr5", 100_000),
        ],
    ),
    # One reciprocal event: both rows carry the same record id via MATEID.
    _case(
        "reciprocal_pair_hits_both_genes",
        "both_single_sv.vcf",
        [
            ("bnd_a|bnd_b", "chr17", 7_670_000),
            ("bnd_a|bnd_b", "chr17", 7_700_000),
            ("bnd_a|bnd_b", "chr17", 7_700_000),
            ("bnd_a|bnd_b", "chr17", 7_670_000),
        ],
    ),
    # One breakpoint inside the 1,479 bp TP53/WRAP53 overlap. Unreachable under the
    # pre-erratum coordinate table, where the genes could not intersect at all.
    _case(
        "single_breakpoint_in_gene_overlap",
        "overlap_single_breakpoint.vcf",
        [
            ("bnd_overlap", "chr17", 7_686_800),
            ("bnd_overlap", "chr8", 120_000),
        ],
    ),
    _case(
        "two_events_one_per_gene",
        "both_independent.vcf",
        [
            ("bnd_a", "chr17", 7_670_000),
            ("bnd_a", "chr5", 100_000),
            ("bnd_c", "chr17", 7_700_000),
            ("bnd_c", "chr9", 400_000),
        ],
    ),
    # Five in-locus records, two of which hit coding windows: Complex must still win.
    _case(
        "five_records_trips_complex",
        "complex.vcf",
        [
            ("c1", "chr17", 7_665_000),
            ("c1", "chr5", 100_000),
            ("c2", "chr17", 7_670_000),
            ("c2", "chr5", 200_000),
            ("c3", "chr17", 7_690_000),
            ("c3", "chr5", 300_000),
            ("c4", "chr17", 7_700_000),
            ("c4", "chr5", 400_000),
            ("c5", "chr17", 7_705_000),
            ("c5", "chr5", 500_000),
        ],
    ),
    _case(
        "in_locus_but_no_coding_hit",
        "neither_in_locus.vcf",
        [
            ("n1", "chr17", 7_665_000),
            ("n1", "chr5", 100_000),
            # Post-erratum, 7,687,550 is inside WRAP53. 7,705,000 is in the 3' flank:
            # in the locus, in no coding window, which is what this case is for.
            ("n2", "chr17", 7_705_000),
            ("n2", "chr5", 200_000),
        ],
    ),
    _case(
        "focal_deletion_spanning_both_genes",
        "svlen_del.vcf",
        [
            ("del_focal", "chr17", 7_670_000),
            # 7,705,000 was inside WRAP53 only under the pre-erratum table; the
            # real gene ends at 7,703,502, so the deletion end moves inside it.
            ("del_focal", "chr17", 7_700_000),
        ],
    ),
    _case(
        "insertion_contributes_one_breakpoint",
        "insertion.vcf",
        [
            ("ins_a", "chr17", 7_700_000),
        ],
    ),
    _case(
        "sites_only_vcf",
        "sites_only.vcf",
        [
            ("bnd_a", "chr17", 7_670_000),
            ("bnd_a", "chr5", 100_000),
        ],
    ),
    _case(
        "breakpoint_past_locus_window",
        "wrap53_outside_window.vcf",
        [
            ("bnd_far", "chr17", 7_715_000),
            ("bnd_far", "chr5", 100_000),
        ],
    ),
    _case(
        "filter_default_keeps_pass_only",
        "filtered.vcf",
        [
            ("passing", "chr17", 7_670_000),
            ("passing", "chr5", 200_000),
        ],
    ),
    _case(
        "filter_override_keeps_everything",
        "filtered.vcf",
        [
            ("lowq", "chr17", 7_700_000),
            ("lowq", "chr5", 100_000),
            ("passing", "chr17", 7_670_000),
            ("passing", "chr5", 200_000),
        ],
        pass_only=False,
    ),
    _case(
        "multisample_tumour",
        "multisample.vcf",
        [
            ("bnd_t", "chr17", 7_700_000),
            ("bnd_t", "chr5", 100_000),
        ],
        sample="TUMOR",
    ),
    _case(
        "multisample_normal",
        "multisample.vcf",
        [
            ("bnd_n", "chr17", 7_670_000),
            ("bnd_n", "chr5", 200_000),
        ],
        sample="NORMAL",
    ),
    # --- divergences -------------------------------------------------------------
    _case(
        "locus_spanning_deletion",
        "spans_locus.vcf",
        [("del_broad", "chr17", 7_000_000), ("del_broad", "chr17", 8_000_000)],
        expected="Neither",
        divergence=(
            "Breakpoint semantics agree with the reference here only because this "
            "generator feeds breakpoints. The reference's segment path would call this "
            "Both-single-SV from interval overlap; that path is a forbidden variable in "
            "the Stage 1 pre-registration. This library reports Neither plus the "
            "spans_locus flag. See SCOPE.md section 1."
        ),
    ),
    _case(
        "contig_without_chr_prefix",
        "chrom_naming.vcf",
        [("bnd_a", "17", 7_670_000), ("bnd_a", "5", 100_000)],
        expected="TP53-only",
        divergence=(
            "The reference compares contig names literally against 'chr17', so a VCF "
            "using '17' is treated as entirely off-locus and returns Neither. VCFs in "
            "the wild use both spellings, so this library normalises '17' to 'chr17' "
            "and returns TP53-only. See SCOPE.md section 6."
        ),
    ),
]


def _git(repo: Path, *args: str) -> str:
    return subprocess.run(
        ["git", "-C", str(repo), *args], capture_output=True, text=True, check=True
    ).stdout.strip()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-repo", type=Path, required=True)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).parent.parent / "tests/data/stage1_equivalence.json",
    )
    parser.add_argument(
        "--allow-any-commit",
        action="store_true",
        help="Regenerate against a different reference commit (then update EXPECTED_COMMIT).",
    )
    args = parser.parse_args()

    repo: Path = args.reference_repo.resolve()
    module = repo / REFERENCE_MODULE
    if not module.exists():
        print(f"reference module not found: {module}", file=sys.stderr)
        return 2

    head = _git(repo, "rev-parse", "HEAD")
    if head != EXPECTED_COMMIT and not args.allow_any_commit:
        print(
            f"reference repo is at {head}, expected {EXPECTED_COMMIT}.\n"
            "Check out the expected commit, or pass --allow-any-commit and update "
            "EXPECTED_COMMIT in this script.",
            file=sys.stderr,
        )
        return 2
    if _git(repo, "status", "--porcelain"):
        print("reference repo has uncommitted changes; refusing to pin", file=sys.stderr)
        return 2

    sys.path.insert(0, str(repo / "part1_alt_classification"))

    # Imported here, not at module scope: pandas and the reference module are needed
    # only by this generator and are not dependencies of the package or its test suite.
    import pandas as pd
    from src.sv_annotator import annotate_breakpoints, architecture_call  # type: ignore[import-not-found]

    results = []
    for case in CASES:
        frame = pd.DataFrame(
            [
                {"record_id": rid, "Chrom": chrom, "Pos": pos}
                for rid, chrom, pos in case["breakpoints"]
            ]
        )
        annotated = annotate_breakpoints(
            frame, case["name"], chrom_col="Chrom", pos_col="Pos", record_id_col="record_id"
        )
        # SCOPE.md decision 2: filter to the locus window before the ladder.
        in_locus = annotated[annotated["in_locus"]] if len(annotated) else annotated
        category, reason = architecture_call(in_locus)

        expected = case["expected_category_override"] or category
        results.append(
            {
                "name": case["name"],
                "fixture": case["fixture"],
                "options": case["options"],
                "breakpoints": case["breakpoints"],
                "reference_category": category,
                "reference_reason": reason,
                "expected_category": expected,
                "divergence": case["divergence"],
            }
        )

    payload = {
        "provenance": {
            "description": (
                "Expected architecture calls pinned from the reference implementation's "
                "category ladder. These are DepMap-based inferences, not a published "
                "gold standard."
            ),
            "reference_repository": (
                "https://github.com/rinatrizvanov/osteosarcoma-alt-dependency-research"
            ),
            "reference_commit": head,
            "reference_module": str(REFERENCE_MODULE),
            "reference_module_sha256": hashlib.sha256(module.read_bytes()).hexdigest(),
            "generated_by": "scripts/generate_equivalence.py",
            "generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        },
        "cases": results,
    }
    args.out.write_text(json.dumps(payload, indent=2) + "\n")

    diverged = [c for c in results if c["divergence"]]
    print(f"wrote {args.out} — {len(results)} cases, {len(diverged)} documented divergence(s)")
    for case in results:
        marker = "  !" if case["divergence"] else "   "
        note = ""
        if case["expected_category"] != case["reference_category"]:
            note = f"  (this library: {case['expected_category']})"
        print(f"{marker} {case['name']:40s} reference={case['reference_category']}{note}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
