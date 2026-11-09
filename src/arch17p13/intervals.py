"""Position-to-feature lookups for the 17p13.1 locus.

This module is the only place that decides whether a coordinate is "in the locus"
and which features it touches. It works on single positions, never on intervals:
interval-overlap semantics are out of scope for v0.1 and the reason is recorded in
``SCOPE.md``.
"""

from __future__ import annotations

from .regions import FEATURES, LOCUS_CHROM, LOCUS_END, LOCUS_START

#: Contig spellings accepted for the locus contig. Anything else is off-locus.
_CHROM_ALIASES = frozenset({"chr17", "17"})


def normalise_chrom(chrom: str | None) -> str | None:
    """Return the canonical contig name, or ``None`` if it is not the locus contig.

    ``chr17`` and ``17`` both normalise to ``chr17``; comparison is case-insensitive
    so ``Chr17`` is accepted too. Every other contig returns ``None``, which callers
    treat as off-locus rather than as an error: a whole-genome VCF is expected to
    contain mostly irrelevant contigs.

    >>> normalise_chrom("17")
    'chr17'
    >>> normalise_chrom("chrX") is None
    True
    """
    if chrom is None:
        return None
    token = chrom.strip().lower()
    if token in _CHROM_ALIASES:
        return LOCUS_CHROM
    return None


def position_in_locus(chrom: str | None, pos: int) -> bool:
    """Whether ``pos`` falls inside the locus window on the locus contig.

    The window is half-open: ``LOCUS_START <= pos < LOCUS_END``.
    """
    return normalise_chrom(chrom) == LOCUS_CHROM and LOCUS_START <= pos < LOCUS_END


def features_for_position(chrom: str | None, pos: int) -> tuple[str, ...]:
    """Names of every feature containing ``pos``, in :data:`FEATURES` order.

    Returns an empty tuple for positions outside the locus window, and may return
    several names because the declared features overlap by design — a breakpoint in
    the intergenic window between the two genes reports both
    ``intergenic_TP53_WRAP53`` and whichever coding or promoter window contains it.

    >>> features_for_position("chr17", 7_670_000)
    ('TP53_coding',)
    >>> features_for_position("chr17", 1_000) == ()
    True
    """
    canonical = normalise_chrom(chrom)
    if canonical is None or not (LOCUS_START <= pos < LOCUS_END):
        return ()
    return tuple(
        name for name, feat_chrom, start, end in FEATURES
        if feat_chrom == canonical and start <= pos < end
    )
