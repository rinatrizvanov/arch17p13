"""Shared paths and helpers for the test suite."""

from __future__ import annotations

from pathlib import Path

import pytest

DATA = Path(__file__).parent / "data"

_HEADER = """##fileformat=VCFv4.3
##reference=GRCh38
##contig=<ID=chr17,length=83257441>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Length of the variant">
##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakend">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
"""


@pytest.fixture
def data_dir() -> Path:
    """Directory holding the committed fixture VCFs."""
    return DATA


def write_vcf(
    path: Path,
    rows: list[tuple[str, ...]],
    samples: tuple[str, ...] = ("SAMPLE1",),
    header: str = _HEADER,
) -> Path:
    """Write a minimal VCF with ``rows`` as tab-separated data lines."""
    columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    if samples:
        columns += ["FORMAT", *samples]
    body = header + "#" + "\t".join(columns) + "\n"
    for row in rows:
        body += "\t".join(row) + "\n"
    path.write_text(body)
    return path
