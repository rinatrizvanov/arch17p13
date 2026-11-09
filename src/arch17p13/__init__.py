"""Classify 17p13.1 TP53/WRAP53 structural-variant architecture from a VCF.

Quickstart
----------
>>> from arch17p13 import classify_vcf
>>> result = classify_vcf("sample.vcf.gz")          # doctest: +SKIP
>>> result.category, result.reason                   # doctest: +SKIP

v0.1 scope: VCF input only, GRCh38 only, breakpoint semantics only. The frozen
scope and its rationale are in ``SCOPE.md``; the research code this reimplements is
listed under "Reference Implementation" in ``README.md``.
"""

from __future__ import annotations

from .classify import SPANS_LOCUS, Category, Result, classify_records, classify_vcf
from .regions import ASSEMBLY, CATEGORIES, COMPLEX_THRESHOLD, FEATURES
from .vcf import Breakpoint, SvRecord, VcfError, read_records

__version__ = "0.1.0"

__all__ = [
    "ASSEMBLY",
    "CATEGORIES",
    "COMPLEX_THRESHOLD",
    "FEATURES",
    "SPANS_LOCUS",
    "Breakpoint",
    "Category",
    "Result",
    "SvRecord",
    "VcfError",
    "__version__",
    "classify_records",
    "classify_vcf",
    "read_records",
]
