# Contributing

## Scope first

This project is **deliberately narrow**, and [`SCOPE.md`](SCOPE.md) is the contract. Before
opening a pull request, check whether what you want to add is listed as out of scope for
v0.1 — BED input, GRCh37 support, copy-number segment semantics and anything inferring ALT
biology all are. Those are not rejections in principle; they are new minor versions with
their own scope notes, so please open an issue first rather than a pull request.

Two rules that are not negotiable within a patch release:

1. **Every constant in `src/arch17p13/regions.py` changes classifications.** Touching a
   coordinate, the locus window or `COMPLEX_THRESHOLD` requires a version bump and a
   `CHANGELOG.md` entry stating what moved and why. `tests/test_regions.py` will fail first,
   which is the point.
2. **Divergences from the reference implementation must be written down.** If a change makes
   this library disagree with the pinned reference calls, document it in `SCOPE.md`, add a
   divergence note in `scripts/generate_equivalence.py`, and assert it by name in
   `tests/test_equivalence.py`. The test
   `test_undocumented_cases_agree_with_the_reference_exactly` exists to stop behaviour being
   redefined silently.

## Running the checks

```bash
uv sync
uv run pytest --cov=arch17p13   # coverage gate is 100%, enforced
uv run ruff check .
uv run ruff format .
uv run mypy                     # strict, covers src/ and tests/
```

Tests are plain pytest with no plugins beyond `pytest-cov`. Fixture VCFs live in
`tests/data/`; add one per behaviour rather than overloading an existing file, and prefer a
fixture over constructing records inline when the behaviour involves parsing.

If the suite prints `PytestWarning: (rm_rf) error removing ...` lines, those come from stale
temporary directories left by other projects in the shared temp root, not from this suite.
`uv run pytest --basetemp=.pytest-basetemp` isolates it.

## Regenerating the equivalence pins

Needs a checkout of the reference repository and that project's pandas environment. The
generator refuses to run unless the reference is at the expected commit with a clean tree:

```bash
uv run --with pandas python scripts/generate_equivalence.py \
    --reference-repo ../osteosarcoma-alt-dependency-research
```

To pin against a different reference commit, pass `--allow-any-commit` and update
`EXPECTED_COMMIT` in the generator and `EXPECTED_REFERENCE_COMMIT` in
`tests/test_equivalence.py` in the same change.

## Style

Type hints everywhere; `mypy --strict` must pass on tests as well as `src/`. Docstrings should
say *why* a rule exists when the reason is not obvious from the code — most of the non-trivial
behaviour here follows from the locus geometry or from a decision recorded in `SCOPE.md`, and
a reader should not have to rediscover either.
