"""Private input locations for the reorganized collaborator scripts.

This module provides configuration, not analysis. Keep the populated JSON outside
the repository. Importing a plotting script can still load data or create folders;
the repository checks therefore parse source text instead of importing scripts.
"""

from __future__ import annotations

import json
import os
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parent
_MISSING = object()


def _outside_repository(path: Path) -> Path:
    path = path.expanduser().resolve()
    if path == REPOSITORY or REPOSITORY in path.parents:
        raise ValueError("Private configuration, inputs and outputs must be outside the repository")
    return path


def configuration() -> dict:
    filename = os.environ.get("PDAC_PRIVATE_CONFIG")
    if not filename:
        raise RuntimeError("Set PDAC_PRIVATE_CONFIG to a populated private JSON file; see docs/inputs.md")
    with _outside_repository(Path(filename)).open(encoding="utf-8") as handle:
        return json.load(handle)


def private_path(name: str) -> Path:
    value = configuration().get("paths", {}).get(name)
    if not value or value.startswith("/path/to/"):
        raise ValueError(f"Supply a real private location for paths.{name}")
    return _outside_repository(Path(value))


def project_root() -> Path:
    return private_path("project_root")


def figure_root() -> Path:
    return private_path("figures_out")


def cohort(name: str) -> list[str]:
    values = configuration().get("cohorts", {}).get(name)
    if not isinstance(values, list) or not values:
        raise ValueError(f"Populate cohorts.{name} in private configuration, preserving source sample order")
    values = [str(value) for value in values]
    if len(values) != len(set(values)):
        raise ValueError(f"Duplicate identifiers in cohorts.{name}")
    # The supplied file parsers expect numeric identifiers and N/T suffixes.
    if not all(value.isdigit() for value in values):
        raise ValueError(f"cohorts.{name} must contain numeric patient identifiers without tissue suffixes")
    return values


def private_value(name: str, default=_MISSING):
    values = configuration().get("values", {})
    if name in values:
        return values[name]
    if default is not _MISSING:
        return default
    raise ValueError(f"Populate values.{name} in private configuration")


def code_path(relative: str) -> Path:
    result = (REPOSITORY / relative).resolve()
    if REPOSITORY not in result.parents:
        raise ValueError("Code paths must stay inside the repository")
    return result
