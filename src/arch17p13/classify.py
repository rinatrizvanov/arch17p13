"""The architecture category ladder.

Mirrors ``architecture_call`` in the reference implementation
(``part1_alt_classification/src/sv_annotator.py``), including its precedence: the
``Complex`` check runs before any gene-level logic, so a sample with enough SV
records is never given a gene call.

Two faithfulness notes, both recorded in ``SCOPE.md``:

* The reference counts every record it is handed, leaving the locus filter to its
  caller, while its GRCh37 twin counts ``n_in_window``. A VCF is whole-genome, so
  this module filters to the locus window first — matching the documented intent
  and the reason strings the reference emits ("in locus window").
* Records that span the locus without placing a breakpoint inside it are reported
  via the ``spans_locus`` flag and never contribute a gene hit.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

from .regions import COMPLEX_THRESHOLD, TP53_CODING, WRAP53_CODING
from .vcf import SvRecord, read_records

Category = Literal[
    "Complex",
    "Both-single-SV",
    "Both-independent-SVs",
    "TCAB1-only",
    "TP53-only",
    "Neither",
    "unknown",
]

#: Flag set when at least one record covers the locus with no breakpoint inside it.
SPANS_LOCUS = "spans_locus"


@dataclass(frozen=True, slots=True)
class Result:
    """The outcome of classifying one sample.

    Attributes
    ----------
    category
        The architecture call.
    reason
        Human-readable justification, in the spirit of the reference implementation's
        second return value.
    record_ids
        Distinct in-locus SV record identifiers that informed the call, sorted.
    n_locus_records
        How many distinct records placed a breakpoint in the locus window. This is
        the number compared against :data:`COMPLEX_THRESHOLD`.
    flags
        Advisory markers that do not change the category, such as
        :data:`SPANS_LOCUS`.
    """

    category: Category
    reason: str
    record_ids: tuple[str, ...] = ()
    n_locus_records: int = 0
    flags: frozenset[str] = field(default_factory=frozenset)

    def as_dict(self) -> dict[str, object]:
        """A JSON-serialisable view, with ``flags`` sorted for stable output."""
        return {
            "category": self.category,
            "reason": self.reason,
            "record_ids": list(self.record_ids),
            "n_locus_records": self.n_locus_records,
            "flags": sorted(self.flags),
        }


def _features_by_record(records: Iterable[SvRecord]) -> dict[str, frozenset[str]]:
    """Union the features of every row sharing a ``record_id``.

    A MATEID-linked breakend pair contributes two rows with one id; unioning their
    features is what lets a reciprocal event hitting each gene once be recognised as
    a single SV, and so as ``Both-single-SV`` rather than ``Both-independent-SVs``.
    """
    merged: dict[str, frozenset[str]] = {}
    for record in records:
        if not record.in_locus:
            continue
        merged[record.record_id] = merged.get(record.record_id, frozenset()) | record.features
    return merged


def classify_records(records: Iterable[SvRecord]) -> Result:
    """Apply the category ladder to already-parsed records."""
    materialised = list(records)
    spans = any(record.spans_locus for record in materialised)
    flags = frozenset({SPANS_LOCUS}) if spans else frozenset()

    by_record: Mapping[str, frozenset[str]] = _features_by_record(materialised)
    record_ids = tuple(sorted(by_record))
    n_records = len(record_ids)

    if n_records == 0:
        reason = "No SVs in locus window"
        if spans:
            reason += "; one or more SVs span the locus without a breakpoint inside it"
        return Result("Neither", reason, (), 0, flags)

    if n_records >= COMPLEX_THRESHOLD:
        return Result(
            "Complex",
            f"{n_records} distinct SV records in locus — manual review required",
            record_ids,
            n_records,
            flags,
        )

    tp53 = {rid for rid, feats in by_record.items() if TP53_CODING in feats}
    wrap53 = {rid for rid, feats in by_record.items() if WRAP53_CODING in feats}
    both_in_one = tp53 & wrap53

    if both_in_one:
        return Result(
            "Both-single-SV",
            f"Records hitting both WRAP53_coding and TP53_coding: {sorted(both_in_one)}",
            record_ids,
            n_records,
            flags,
        )

    if tp53 and wrap53:
        return Result(
            "Both-independent-SVs",
            f"WRAP53_coding records: {sorted(wrap53)}; TP53_coding records: {sorted(tp53)}",
            record_ids,
            n_records,
            flags,
        )

    if wrap53:
        return Result(
            "TCAB1-only",
            f"WRAP53_coding records: {sorted(wrap53)}",
            record_ids,
            n_records,
            flags,
        )

    if tp53:
        return Result(
            "TP53-only",
            f"TP53_coding records: {sorted(tp53)}",
            record_ids,
            n_records,
            flags,
        )

    return Result(
        "Neither",
        f"{n_records} SV(s) in locus window but none hit TP53_coding or WRAP53_coding — "
        "flag for review",
        record_ids,
        n_records,
        flags,
    )


def classify_vcf(
    path: str | Path,
    *,
    sample: str | None = None,
    pass_only: bool = True,
) -> Result:
    """Read a VCF and classify it. See :func:`arch17p13.vcf.read_records` for options."""
    return classify_records(read_records(path, sample=sample, pass_only=pass_only))
