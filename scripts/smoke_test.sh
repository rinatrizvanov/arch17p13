#!/usr/bin/env bash
#
# End-to-end check of an *installed* arch17p13, as opposed to the test suite, which
# imports the package from the source tree. This is the only check that exercises the
# built wheel: the console-script entry point, the package metadata and the fact that
# tests/data/ fixtures are not needed at runtime.
#
# Usage:  scripts/smoke_test.sh [path-to-arch17p13-executable]
#
# With no argument it uses whatever `arch17p13` is on PATH. Run it from anywhere; paths
# are resolved relative to the repository, not the working directory.

set -euo pipefail

CLI="${1:-arch17p13}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA="$REPO_ROOT/tests/data"

expected_version="$(sed -n 's/^__version__ = "\(.*\)"$/\1/p' "$REPO_ROOT/src/arch17p13/__init__.py")"

fail() {
    echo "SMOKE FAILURE: $*" >&2
    exit 1
}

check_eq() {
    local what="$1" expected="$2" actual="$3"
    if [[ "$expected" != "$actual" ]]; then
        fail "$what: expected '$expected', got '$actual'"
    fi
    echo "  ok  $what -> $actual"
}

echo "smoke-testing: $CLI"
command -v "$CLI" >/dev/null 2>&1 || [[ -x "$CLI" ]] || fail "'$CLI' is not executable"

# 1. The installed version is the version in the source tree. Catches a stale wheel.
echo "version"
check_eq "--version" "arch17p13 $expected_version" "$("$CLI" --version)"

# 2. One fixture per category, so a packaging change that broke the region constants or
#    the classification ladder cannot pass.
echo "categories"
while read -r fixture expected; do
    [[ -z "$fixture" ]] && continue
    actual="$("$CLI" classify "$DATA/$fixture" | cut -f2)"
    check_eq "$fixture" "$expected" "$actual"
done <<'CASES'
both_single_sv.vcf	Both-single-SV
both_independent.vcf	Both-independent-SVs
tp53_only.vcf	TP53-only
tcab1_only.vcf	TCAB1-only
complex.vcf	Complex
neither_offlocus.vcf	Neither
chrom_naming.vcf	TP53-only
sites_only.vcf	TP53-only
CASES

# 3. The spans_locus flag survives serialisation. This is the documented divergence from
#    the reference implementation, so it is the one flag worth asserting end to end.
echo "flags"
flags="$("$CLI" classify "$DATA/spans_locus.vcf" --json | tr -d ' \n' | grep -o '"flags":\[[^]]*\]')"
check_eq "spans_locus.vcf flags" '"flags":["spans_locus"]' "$flags"

# 4. Exit statuses, which the CLI contract makes meaningful: a category is a result, not
#    an error; only I/O and usage problems are non-zero.
echo "exit statuses"
status=0; "$CLI" classify "$DATA/neither_offlocus.vcf" >/dev/null || status=$?
check_eq "Neither exits 0" "0" "$status"

status=0; "$CLI" classify "$DATA/does_not_exist.vcf" >/dev/null 2>&1 || status=$?
check_eq "missing file exits 1" "1" "$status"

status=0; "$CLI" classify "$DATA/multisample.vcf" >/dev/null 2>&1 || status=$?
check_eq "ambiguous sample exits 1" "1" "$status"

status=0; "$CLI" --no-such-flag >/dev/null 2>&1 || status=$?
check_eq "usage error exits 2" "2" "$status"

echo "smoke tests passed"
