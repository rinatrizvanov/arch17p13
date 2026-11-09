"""Allow ``python -m arch17p13`` alongside the installed ``arch17p13`` script."""

from __future__ import annotations

from .cli import main

raise SystemExit(main())
