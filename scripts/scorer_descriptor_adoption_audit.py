#!/usr/bin/env python3
"""Inventory direct and wrapper-mediated ObservableDescriptor adoption."""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
from pathlib import Path
from typing import Mapping


def git_blob_sha(data: bytes) -> str:
    header = f"blob {len(data)}\0".encode("ascii")
    return hashlib.sha1(header + data).hexdigest()


def imports_descriptor_map(path: Path) -> bool:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    imported: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module == "wrapping_channels":
            imported.update(alias.name for alias in node.names)
    return {"ObservableDescriptor", "map_observable"} <= imported


def audit(root: Path, manifest: Mapping[str, object]) -> dict:
    paths = sorted(path.relative_to(root).as_posix() for path in root.glob(str(manifest["corpus_glob"])))
    direct = {path for path in paths if imports_descriptor_map(root / path)}
    wrappers = dict(manifest["typed_wrapper_kernels"])
    standalone = set(manifest["direct_typed_standalone"])

    if set(wrappers) | standalone != direct:
        raise ValueError("manifest direct-typed paths do not equal AST-detected paths")
    missing_kernels = sorted(set(wrappers.values()) - set(paths))
    if missing_kernels:
        raise ValueError("typed wrapper kernels missing from corpus: " + ", ".join(missing_kernels))

    covered_kernels = set(wrappers.values())
    rows = []
    for path in paths:
        if path in direct:
            status = "direct_typed_entrypoint"
        elif path in covered_kernels:
            status = "covered_frozen_kernel"
        else:
            status = "outside_registered_typed_path"
        rows.append(
            {
                "path": path,
                "git_blob_sha": git_blob_sha((root / path).read_bytes()),
                "status": status,
            }
        )

    counts = {
        status: sum(row["status"] == status for row in rows)
        for status in (
            "direct_typed_entrypoint",
            "covered_frozen_kernel",
            "outside_registered_typed_path",
        )
    }
    return {
        "schema": manifest["schema"],
        "corpus_glob": manifest["corpus_glob"],
        "counts": {"total": len(rows), **counts},
        "typed_wrapper_kernels": wrappers,
        "rows": rows,
        "boundary": manifest["boundary"],
    }


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=root)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=root / "analysis/scorer_descriptor_adoption_manifest.json",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    result = audit(args.root, manifest)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
