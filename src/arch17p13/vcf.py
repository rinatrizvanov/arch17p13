"""Minimal VCF reader for structural variants at a single locus.

Deliberately not a general VCF library. It reads enough of the spec to turn SV
records into breakpoint positions at the 17p13.1 locus, using only the standard
library: the locus is 150 kb, so there is nothing to gain from an index, and
avoiding a compiled dependency keeps the package installable anywhere.

Supported record shapes
-----------------------
``<DEL>``, ``<DUP>``, ``<INV>``, ``<CNV>``
    Two breakpoints, at ``POS`` and at the record's end coordinate (``INFO/END``,
    or derived from ``INFO/SVLEN`` when ``END`` is absent).
``<INS>``
    One breakpoint, at ``POS``. An insertion has no span to contribute.
``BND``
    Two breakpoints: ``POS``, and the mate coordinate parsed out of the ALT
    bracket notation. Mates named by ``INFO/MATEID`` collapse to one record, so a
    reciprocal pair counts once both for the ``Both-single-SV`` versus
    ``Both-independent-SVs`` distinction and for the ``Complex`` threshold.

Anything else contributes its ``POS`` alone.
"""

from __future__ import annotations

import gzip
import re
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import IO, Literal

from .intervals import features_for_position, normalise_chrom, position_in_locus
from .regions import LOCUS_END, LOCUS_START, TP53_CODING, WRAP53_CODING

#: ALT bracket notation for a breakend mate, e.g. ``N[chr17:7687600[``.
_BND_ALT = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")

#: Symbolic ALT allele, e.g. ``<DEL>``.
_SYMBOLIC_ALT = re.compile(r"^<([A-Za-z0-9:_]+)>$")

#: SV types that describe a span, and so contribute both endpoints.
_SPANNING_TYPES = frozenset({"DEL", "DUP", "INV", "CNV", "DUP:TANDEM"})

#: FILTER values treated as passing.
_PASSING_FILTERS = frozenset({"PASS", ".", ""})

BreakpointOrigin = Literal["POS", "END", "MATE"]


class VcfError(ValueError):
    """Raised when a VCF cannot be parsed, or when the caller's request is ambiguous."""


@dataclass(frozen=True, slots=True)
class Breakpoint:
    """A single breakpoint position and the locus features it falls in."""

    chrom: str
    pos: int
    features: tuple[str, ...]
    in_locus: bool
    origin: BreakpointOrigin


@dataclass(frozen=True, slots=True)
class SvRecord:
    """One structural variant, reduced to the breakpoints it contributes.

    ``record_id`` is the unit the classifier counts. Two breakend rows joined by
    ``INFO/MATEID`` share one ``record_id``; every other row gets its own.
    """

    record_id: str
    sv_type: str
    chrom: str
    pos: int
    end: int | None
    filter_value: str
    breakpoints: tuple[Breakpoint, ...]
    spans_locus: bool

    @property
    def locus_breakpoints(self) -> tuple[Breakpoint, ...]:
        """Only the breakpoints inside the locus window."""
        return tuple(bp for bp in self.breakpoints if bp.in_locus)

    @property
    def features(self) -> frozenset[str]:
        """Every locus feature touched by any breakpoint of this record."""
        return frozenset(name for bp in self.breakpoints for name in bp.features)

    @property
    def in_locus(self) -> bool:
        """Whether this record places at least one breakpoint in the locus window."""
        return any(bp.in_locus for bp in self.breakpoints)

    @property
    def hits_tp53_coding(self) -> bool:
        return TP53_CODING in self.features

    @property
    def hits_wrap53_coding(self) -> bool:
        return WRAP53_CODING in self.features


def _open_text(path: Path) -> IO[str]:
    """Open a plain or gzip/bgzip-compressed VCF for sequential text reading."""
    with path.open("rb") as probe:
        magic = probe.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("rt", encoding="utf-8")


def _parse_info(field: str) -> dict[str, str]:
    """Parse an INFO column into a dict. Valueless flags map to the empty string."""
    info: dict[str, str] = {}
    if field in (".", ""):
        return info
    for entry in field.split(";"):
        if not entry:
            continue
        key, _, value = entry.partition("=")
        info[key] = value
    return info


def _make_breakpoint(chrom: str, pos: int, origin: BreakpointOrigin) -> Breakpoint:
    canonical = normalise_chrom(chrom)
    in_locus = position_in_locus(chrom, pos)
    return Breakpoint(
        chrom=canonical if canonical is not None else chrom,
        pos=pos,
        features=features_for_position(chrom, pos) if in_locus else (),
        in_locus=in_locus,
        origin=origin,
    )


def _resolve_sv_type(alt: str, info: dict[str, str]) -> str:
    """Determine SVTYPE, preferring INFO and falling back to the symbolic ALT."""
    declared = info.get("SVTYPE", "").strip().upper()
    if declared:
        return declared
    symbolic = _SYMBOLIC_ALT.match(alt.strip())
    if symbolic:
        return symbolic.group(1).upper()
    if _BND_ALT.search(alt):
        return "BND"
    return ""


def _resolve_end(pos: int, sv_type: str, info: dict[str, str], line_no: int) -> int | None:
    """End coordinate for a spanning record, from INFO/END or INFO/SVLEN."""
    raw_end = info.get("END", "").strip()
    if raw_end:
        try:
            return int(raw_end)
        except ValueError as exc:
            raise VcfError(f"line {line_no}: INFO/END is not an integer: {raw_end!r}") from exc
    raw_svlen = info.get("SVLEN", "").split(",")[0].strip()
    if raw_svlen and sv_type in _SPANNING_TYPES:
        try:
            return pos + abs(int(raw_svlen))
        except ValueError as exc:
            raise VcfError(f"line {line_no}: INFO/SVLEN is not an integer: {raw_svlen!r}") from exc
    return None


def _record_id(vcf_id: str, sv_type: str, info: dict[str, str], line_no: int) -> str:
    """Stable identifier, collapsing MATEID-linked breakend pairs onto one id."""
    own = vcf_id.strip()
    if own in (".", ""):
        own = f"__line{line_no}"
    if sv_type != "BND":
        return own
    mate_field = info.get("MATEID", "").strip()
    if not mate_field:
        return own
    mates = [m.strip() for m in mate_field.split(",") if m.strip()]
    if not mates:
        return own
    return "|".join(sorted({own, *mates}))


def _spans_locus(
    pos: int, end: int | None, chrom: str, breakpoints: tuple[Breakpoint, ...]
) -> bool:
    """Whether the record covers the whole locus without a breakpoint inside it.

    This is the shape that interval-overlap semantics would call a both-genes hit and
    that this library deliberately does not: a broad deletion of 17p removes both
    genes, but its breakpoints are nowhere near either one, so it is no evidence of
    the focal co-disruption the categories describe. Reported as a flag instead.
    """
    if end is None or normalise_chrom(chrom) is None:
        return False
    if any(bp.in_locus for bp in breakpoints):
        return False
    low, high = min(pos, end), max(pos, end)
    return low <= LOCUS_START and high >= LOCUS_END


def _genotype_carries_variant(gt: str) -> bool:
    """Whether a GT string indicates the sample carries a non-reference allele."""
    alleles = re.split(r"[/|]", gt.split(":")[0].strip())
    return any(a not in ("0", ".", "") for a in alleles)


def _sample_columns(header: list[str]) -> list[str]:
    """Sample names from a ``#CHROM`` header line, or an empty list for a sites-only VCF."""
    return header[9:] if len(header) > 9 else []


def _iter_data_lines(handle: IO[str]) -> Iterator[tuple[int, str, list[str]]]:
    """Yield ``(line_no, line, header_fields)`` for each data line, tracking the header."""
    header: list[str] = []
    for line_no, raw in enumerate(handle, start=1):
        line = raw.rstrip("\n").rstrip("\r")
        if not line:
            continue
        if line.startswith("##"):
            continue
        if line.startswith("#"):
            header = line.lstrip("#").split("\t")
            continue
        yield line_no, line, header


def read_records(
    path: str | Path,
    *,
    sample: str | None = None,
    pass_only: bool = True,
) -> list[SvRecord]:
    """Read SV records from a VCF, keeping only what the locus classifier needs.

    Parameters
    ----------
    path
        Path to a ``.vcf`` or gzip/bgzip-compressed ``.vcf.gz`` file.
    sample
        Which sample column gates inclusion. Required for a multi-sample VCF, where
        records the sample does not carry are skipped. For a single-sample or
        sites-only VCF this may be left ``None``.
    pass_only
        Keep only records whose FILTER is ``PASS`` or ``.`` (the default). Set
        ``False`` to classify on all records regardless of FILTER.

    Returns
    -------
    Records that place at least one breakpoint in the locus window, plus records
    that span the locus without a breakpoint inside it (carrying ``spans_locus``).
    Everything else is dropped, since it cannot affect the call.

    Raises
    ------
    VcfError
        If the file is malformed, or if it holds several samples and ``sample`` is
        either unset or not among them.
    """
    vcf_path = Path(path)
    if not vcf_path.exists():
        raise VcfError(f"VCF not found: {vcf_path}")

    kept: list[SvRecord] = []

    with _open_text(vcf_path) as handle:
        for line_no, line, header in _iter_data_lines(handle):
            fields = line.split("\t")
            if len(fields) < 8:
                raise VcfError(
                    f"line {line_no}: expected at least 8 tab-separated columns, got {len(fields)}"
                )

            samples = _sample_columns(header)
            if sample is None and len(samples) > 1:
                raise VcfError(
                    f"{vcf_path.name} holds {len(samples)} samples "
                    f"({', '.join(samples)}); pass sample= to choose one"
                )
            if sample is not None and samples:
                if sample not in samples:
                    raise VcfError(
                        f"sample {sample!r} not in {vcf_path.name}; available: {', '.join(samples)}"
                    )
                column = 9 + samples.index(sample)
                if column >= len(fields):
                    raise VcfError(f"line {line_no}: no genotype column for sample {sample!r}")
                if not _genotype_carries_variant(fields[column]):
                    continue

            chrom, raw_pos, vcf_id, _ref, alt, _qual, filter_value, raw_info = fields[:8]
            try:
                pos = int(raw_pos)
            except ValueError as exc:
                raise VcfError(f"line {line_no}: POS is not an integer: {raw_pos!r}") from exc

            if pass_only and filter_value.strip().upper() not in {
                f.upper() for f in _PASSING_FILTERS
            }:
                continue

            info = _parse_info(raw_info)
            sv_type = _resolve_sv_type(alt, info)
            end = _resolve_end(pos, sv_type, info, line_no)

            breakpoints = [_make_breakpoint(chrom, pos, "POS")]
            if sv_type == "BND":
                mate = _BND_ALT.search(alt)
                if mate:
                    breakpoints.append(
                        _make_breakpoint(mate.group(1), int(mate.group(2)), "MATE")
                    )
            elif end is not None and sv_type in _SPANNING_TYPES:
                breakpoints.append(_make_breakpoint(chrom, end, "END"))

            frozen = tuple(breakpoints)
            record = SvRecord(
                record_id=_record_id(vcf_id, sv_type, info, line_no),
                sv_type=sv_type,
                chrom=chrom,
                pos=pos,
                end=end,
                filter_value=filter_value,
                breakpoints=frozen,
                spans_locus=_spans_locus(pos, end, chrom, frozen),
            )

            # Rows that cannot affect the call are dropped here. A MATEID-linked pair
            # keeps both rows: classify_records groups them by record_id and unions
            # their features, so the pair still counts as one SV.
            if record.in_locus or record.spans_locus:
                kept.append(record)

    return kept
