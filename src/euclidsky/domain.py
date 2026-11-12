"""The Facetwork domain entry point, kept apart so the tools run without it.

Discovered by the runner through::

    [project.entry-points."facetwork.domains"]
    euclid = "euclidsky.domain:domain"
"""

from __future__ import annotations

from pathlib import Path

from facetwork.domains import DomainPackage

from .handlers import register_handlers

domain = DomainPackage(
    name="euclid",
    ffl_dir=Path(__file__).parent / "ffl",
    register_handlers=register_handlers,
)
