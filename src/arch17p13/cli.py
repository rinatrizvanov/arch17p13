"""Command-line interface.

    arch17p13 classify sample.vcf.gz
    arch17p13 classify tumour.vcf --sample TUMOR --explain
    arch17p13 classify cohort/*.vcf.gz --json

Exit status is 0 for any successful classification, including ``Neither``: a category
is a result, not an error. Status 1 means a file could not be read or classified, and
2 is reserved for argparse usage errors.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import IO

from . import __version__
from .classify import Result, classify_vcf
from .regions import ASSEMBLY, COMPLEX_THRESHOLD, LOCUS_CHROM, LOCUS_END, LOCUS_START
from .vcf import VcfError

#: Written to stderr once per run when more than one file is classified.
_EPILOG = (
    f"Classifies breakpoints in {LOCUS_CHROM}:{LOCUS_START:,}-{LOCUS_END:,} ({ASSEMBLY}). "
    f"{COMPLEX_THRESHOLD} or more distinct in-locus SV records are reported as Complex. "
    "See SCOPE.md for the frozen v0.1 scope and its divergences from the reference "
    "implementation."
)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="arch17p13",
        description="Classify 17p13.1 TP53/WRAP53 structural-variant architecture from a VCF.",
        epilog=_EPILOG,
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    classify = subparsers.add_parser(
        "classify",
        help="classify one or more VCF files",
        description="Classify one or more VCF files.",
        epilog=_EPILOG,
    )
    classify.add_argument("vcf", nargs="+", type=Path, help="VCF file(s), optionally gzipped")
    classify.add_argument(
        "--sample",
        help="sample column to classify; required for a multi-sample VCF",
    )
    classify.add_argument(
        "--all-filters",
        action="store_true",
        help="classify records regardless of FILTER (default: PASS and '.' only)",
    )
    classify.add_argument(
        "--explain",
        action="store_true",
        help="also print the reason, the contributing record ids and any flags",
    )
    classify.add_argument(
        "--json",
        dest="as_json",
        action="store_true",
        help="emit JSON instead of text; one object per file, or an array for several",
    )
    return parser


def _format_text(path: Path, result: Result, explain: bool) -> str:
    """One line per file, or an indented block when ``--explain`` is given."""
    head = f"{path.name}\t{result.category}"
    if not explain:
        return head
    lines = [head, f"  reason: {result.reason}"]
    if result.record_ids:
        lines.append(
            f"  records: {', '.join(result.record_ids)} ({result.n_locus_records} in locus)"
        )
    else:
        lines.append(f"  records: none in locus ({result.n_locus_records} in locus)")
    if result.flags:
        lines.append(f"  flags: {', '.join(sorted(result.flags))}")
    return "\n".join(lines)


def main(
    argv: Sequence[str] | None = None, stdout: IO[str] | None = None, stderr: IO[str] | None = None
) -> int:
    """Entry point. ``stdout``/``stderr`` are injectable so tests need no capture."""
    out = stdout if stdout is not None else sys.stdout
    err = stderr if stderr is not None else sys.stderr

    args = _build_parser().parse_args(argv)

    payloads: list[dict[str, object]] = []
    text: list[str] = []
    failed = False

    for path in args.vcf:
        try:
            result = classify_vcf(path, sample=args.sample, pass_only=not args.all_filters)
        except VcfError as exc:
            print(f"{path}: {exc}", file=err)
            failed = True
            continue

        if args.as_json:
            payload = result.as_dict()
            payload["file"] = str(path)
            payloads.append(payload)
        else:
            text.append(_format_text(path, result, args.explain))

    if args.as_json:
        # A single file emits one object; several emit an array, so consumers do not
        # have to special-case the common one-file invocation.
        document = payloads[0] if len(payloads) == 1 else payloads
        print(json.dumps(document, indent=2, ensure_ascii=False), file=out)
    elif text:
        print("\n".join(text), file=out)

    return 1 if failed else 0
