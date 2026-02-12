#!/usr/bin/env python
"""Collect key TD-DFT metrics from ORCA output files."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


STATE_RE = re.compile(
    r"STATE\s+(\d+).*?([0-9]+\.[0-9]+)\s+eV(?:.*?f[=\s]+([0-9.Ee+-]+))?",
    re.IGNORECASE,
)
ENERGY_RE = re.compile(r"FINAL SINGLE POINT ENERGY\s+(-?[0-9]+\.[0-9]+)")


def parse_orca_output(path: Path) -> dict[str, object]:
    text = path.read_text(encoding="utf-8", errors="ignore")
    lines = text.splitlines()

    terminated = "ORCA TERMINATED NORMALLY" in text
    energies = [float(m.group(1)) for m in ENERGY_RE.finditer(text)]
    final_energy = energies[-1] if energies else float("nan")

    states: list[tuple[int, float, float | None]] = []
    for line in lines:
        if "STATE" not in line or "eV" not in line:
            continue
        m = STATE_RE.search(line)
        if not m:
            continue
        state_idx = int(m.group(1))
        e_ev = float(m.group(2))
        fosc = float(m.group(3)) if m.group(3) is not None else None
        states.append((state_idx, e_ev, fosc))

    s1_e_ev = float("nan")
    s1_fosc = float("nan")
    if states:
        states.sort(key=lambda x: x[0])
        s1_state = states[0]
        s1_e_ev = s1_state[1]
        if s1_state[2] is not None:
            s1_fosc = float(s1_state[2])

    return {
        "terminated_normally": terminated,
        "final_single_point_energy_au": final_energy,
        "n_states_parsed": len(states),
        "s1_energy_ev": s1_e_ev,
        "s1_fosc": s1_fosc,
    }


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--jobs-root",
        default="runs/orca_jobs",
        help="Root folder containing per-candidate ORCA job subfolders",
    )
    p.add_argument(
        "--glob",
        default="**/*.out",
        help="Glob pattern under jobs-root for ORCA output files",
    )
    p.add_argument(
        "--output-csv",
        default="runs/orca_jobs/orca_results_summary.csv",
        help="Summary CSV path",
    )
    return p


def main() -> int:
    args = build_parser().parse_args()
    root = Path(args.jobs_root)
    if not root.exists():
        raise SystemExit(f"jobs-root not found: {root}")

    out_paths = sorted(root.glob(args.glob))
    if not out_paths:
        raise SystemExit(f"No output files found with glob '{args.glob}' in {root}")

    rows = []
    for out_path in out_paths:
        parsed = parse_orca_output(out_path)
        rows.append(
            {
                "candidate_id": out_path.parent.name,
                "stage": out_path.stem,
                "out_file": str(out_path),
                **parsed,
            }
        )

    df = pd.DataFrame(rows).sort_values(["candidate_id", "stage"])
    output_path = Path(args.output_csv)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False, encoding="utf-8-sig")

    ok = int(df["terminated_normally"].sum())
    print(f"Wrote {output_path} with {len(df)} rows | terminated_normally={ok}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
