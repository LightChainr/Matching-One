#!/usr/bin/env python3
"""Typed entrypoint retaining the frozen stable axis-pair batch reader."""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

import score_axis_pair_annihilator_stable as stable
import score_axis_pair_annihilator_typed as typed
from wrapping_channels import ObservableDescriptor, map_observable


def calculate_typed(
    root: Path,
    paths: Sequence[Path],
    p_ref: float,
    train_max_L: int,
    q_candidates: Sequence[float],
) -> dict:
    assert ObservableDescriptor is not None and map_observable is not None
    return typed.calculate_typed(
        root, paths, p_ref, train_max_L, q_candidates,
        kernel_key="stable", calculator=stable.calculate,
    )


def main(argv: Sequence[str] | None = None) -> int:
    return typed.main(
        argv,
        kernel_key="stable",
        calculator=stable.calculate,
        reporter=stable.report,
    )


if __name__ == "__main__":
    raise SystemExit(main())
