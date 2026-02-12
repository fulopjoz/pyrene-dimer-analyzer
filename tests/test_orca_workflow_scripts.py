"""Smoke tests for ORCA workflow helper scripts."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]


def _run(cmd: list[str], cwd: Path | None = None) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        cmd,
        cwd=str(cwd or REPO_ROOT),
        check=True,
        text=True,
        capture_output=True,
    )


def test_orca_select_candidates_cli(tmp_path: Path):
    in_csv = tmp_path / "all_conformers.csv"
    out_csv = tmp_path / "orca_candidates.csv"
    df = pd.DataFrame(
        [
            {
                "name": "MolA",
                "conformer_id": 1,
                "classification": "strong_excimer",
                "pi_overlap_pct": 62.0,
                "interplane_distance_A": 3.5,
                "plane_angle_deg": 8.0,
            },
            {
                "name": "MolA",
                "conformer_id": 2,
                "classification": "weak_excimer",
                "pi_overlap_pct": 34.0,
                "interplane_distance_A": 4.3,
                "plane_angle_deg": 28.0,
            },
            {
                "name": "MolB",
                "conformer_id": 3,
                "classification": "monomer",
                "pi_overlap_pct": 5.0,
                "interplane_distance_A": 6.5,
                "plane_angle_deg": 78.0,
            },
        ]
    )
    df.to_csv(in_csv, index=False)

    _run(
        [
            sys.executable,
            "scripts/orca_select_candidates.py",
            "--input-csv",
            str(in_csv),
            "--output-csv",
            str(out_csv),
            "--top-per-molecule",
            "1",
            "--max-total",
            "2",
            "--max-theta",
            "60",
            "--max-distance",
            "4.6",
            "--min-overlap",
            "30",
        ]
    )

    selected = pd.read_csv(out_csv)
    assert len(selected) == 1
    assert int(selected.iloc[0]["conformer_id"]) == 1
    assert "orca_priority_score" in selected.columns
    assert "orca_candidate_rank" in selected.columns


def test_orca_prepare_inputs_and_manifest(tmp_path: Path):
    candidates_csv = tmp_path / "candidates.csv"
    out_dir = tmp_path / "jobs"
    sdf = REPO_ROOT / "tests" / "test_data" / "pyrene_dimer_set_for_MOE.sdf"
    pd.DataFrame(
        [
            {"name": "SmokeA", "conformer_id": 0},
            {"name": "SmokeB", "conformer_id": 1},
        ]
    ).to_csv(candidates_csv, index=False)

    _run(
        [
            sys.executable,
            "scripts/orca_prepare_inputs.py",
            "--candidates-csv",
            str(candidates_csv),
            "--sdf",
            str(sdf),
            "--output-dir",
            str(out_dir),
            "--stages",
            "vertical,s1opt",
            "--nprocs",
            "2",
            "--maxcore",
            "512",
            "--tda",
        ]
    )

    manifest = pd.read_csv(out_dir / "orca_job_manifest.csv")
    assert len(manifest) == 4  # 2 candidates x 2 stages
    assert set(manifest["stage"]) == {"vertical", "s1opt"}

    first_job = Path(manifest.iloc[0]["workdir"])
    assert (first_job / "geom.xyz").exists()
    assert (first_job / "vertical.inp").exists()
    text = (first_job / "s1opt.inp").read_text(encoding="utf-8")
    assert "followiroot true" in text.lower()


def test_orca_collect_results_parser(tmp_path: Path):
    jobs_root = tmp_path / "orca_jobs"
    cand_dir = jobs_root / "MolA_cid00001"
    cand_dir.mkdir(parents=True, exist_ok=True)
    out_file = cand_dir / "vertical.out"
    out_file.write_text(
        "\n".join(
            [
                "Some ORCA output",
                "FINAL SINGLE POINT ENERGY     -1234.567890123",
                "STATE   1: E=   2.3456 eV  f=0.1234",
                "ORCA TERMINATED NORMALLY",
            ]
        ),
        encoding="utf-8",
    )
    out_csv = tmp_path / "orca_results_summary.csv"

    _run(
        [
            sys.executable,
            "scripts/orca_collect_results.py",
            "--jobs-root",
            str(jobs_root),
            "--glob",
            "**/*.out",
            "--output-csv",
            str(out_csv),
        ]
    )

    df = pd.read_csv(out_csv)
    assert len(df) == 1
    assert bool(df.iloc[0]["terminated_normally"])
    assert abs(float(df.iloc[0]["s1_energy_ev"]) - 2.3456) < 1e-6
