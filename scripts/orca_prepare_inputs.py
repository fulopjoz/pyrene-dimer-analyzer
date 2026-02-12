#!/usr/bin/env python
"""Prepare ORCA excited-state job folders from selected conformer candidates."""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Iterable

import pandas as pd
from rdkit import Chem


STAGE_TO_TEMPLATE = {
    "vertical": "vertical_abs.inp.tpl",
    "s1opt": "s1_opt.inp.tpl",
    "s0freq": "s0_opt_freq.inp.tpl",
    "s1freq": "s1_opt_freq.inp.tpl",
    "esd_abs": "esd_abs.inp.tpl",
    "esd_fluor": "esd_fluor.inp.tpl",
}


def _safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text).strip("_") or "candidate"


def _read_template(template_dir: Path, template_name: str) -> str:
    path = template_dir / template_name
    if not path.exists():
        raise FileNotFoundError(f"Template not found: {path}")
    return path.read_text(encoding="utf-8")


def _build_header_line(
    stage: str,
    functional: str,
    basis: str,
    dispersion: str,
    solvent: str,
) -> str:
    keywords: list[str] = [functional, basis, "TightSCF", "RIJCOSX"]
    if dispersion:
        keywords.append(dispersion)
    if solvent:
        keywords.append(f"CPCM({solvent})")

    if stage == "s1opt":
        keywords.append("Opt")
    elif stage == "s0freq":
        keywords.extend(["Opt", "Freq"])
    elif stage == "s1freq":
        keywords.extend(["Opt", "Freq"])
    elif stage == "esd_abs":
        keywords.append("ESD(ABS)")
    elif stage == "esd_fluor":
        keywords.append("ESD(FLUOR)")

    return " ".join(keywords)


def _iter_stages(stages_arg: str) -> list[str]:
    stages = [s.strip() for s in stages_arg.split(",") if s.strip()]
    if not stages:
        raise ValueError("No valid stages requested.")
    unknown = [s for s in stages if s not in STAGE_TO_TEMPLATE]
    if unknown:
        raise ValueError(f"Unknown stage(s): {unknown}")
    return stages


def _write_xyz(mol: Chem.Mol, xyz_path: Path) -> None:
    xyz = Chem.MolToXYZBlock(mol)
    xyz_path.write_text(xyz, encoding="utf-8")


def _load_sdf(sdf_path: Path) -> list[Chem.Mol | None]:
    supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=True)
    if supplier is None:
        raise ValueError(f"Could not read SDF: {sdf_path}")
    return list(supplier)


def _candidate_id(row: pd.Series) -> str:
    name = _safe_name(str(row.get("name", "candidate")))
    conf = int(row["conformer_id"])
    return f"{name}_cid{conf:05d}"


def _render_input(template: str, values: dict[str, str | int | float]) -> str:
    return template.format(**values)


def prepare_jobs(
    candidates: pd.DataFrame,
    mols: list[Chem.Mol | None],
    output_dir: Path,
    template_dir: Path,
    stages: Iterable[str],
    charge: int,
    multiplicity: int,
    nprocs: int,
    maxcore: int,
    functional: str,
    basis: str,
    dispersion: str,
    solvent: str,
    nroots: int,
    iroot: int,
    tda: bool,
    states: str,
    linew: float,
    inlinew: float,
) -> pd.DataFrame:
    manifest_rows: list[dict[str, object]] = []
    stage_list = list(stages)

    for _, row in candidates.iterrows():
        conf_id = int(row["conformer_id"])
        if conf_id < 0 or conf_id >= len(mols):
            raise IndexError(f"conformer_id {conf_id} out of SDF bounds (n={len(mols)})")
        mol = mols[conf_id]
        if mol is None:
            raise ValueError(f"SDF molecule at conformer_id {conf_id} is None")

        cand_id = _candidate_id(row)
        cand_dir = output_dir / cand_id
        cand_dir.mkdir(parents=True, exist_ok=True)

        xyz_name = "geom.xyz"
        xyz_path = cand_dir / xyz_name
        _write_xyz(mol, xyz_path)

        for stage in stage_list:
            template_name = STAGE_TO_TEMPLATE[stage]
            template_text = _read_template(template_dir, template_name)
            inp_name = f"{stage}.inp"
            inp_path = cand_dir / inp_name

            header_line = _build_header_line(
                stage=stage,
                functional=functional,
                basis=basis,
                dispersion=dispersion,
                solvent=solvent,
            )
            values = {
                "header_line": header_line,
                "nprocs": nprocs,
                "maxcore": maxcore,
                "nroots": nroots,
                "iroot": iroot,
                "tda": str(bool(tda)).lower(),
                "charge": charge,
                "mult": multiplicity,
                "xyz_file": xyz_name,
                "states": states,
                "linew": linew,
                "inlinew": inlinew,
                "gshessian_file": "s0freq.hess",
                "eshessian_file": "s1freq.hess",
            }

            inp_path.write_text(_render_input(template_text, values), encoding="utf-8")

            manifest_rows.append(
                {
                    "candidate_id": cand_id,
                    "name": row.get("name", ""),
                    "conformer_id": conf_id,
                    "stage": stage,
                    "input_file": str(inp_path),
                    "workdir": str(cand_dir),
                    "nprocs": nprocs,
                }
            )

    return pd.DataFrame(manifest_rows)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--candidates-csv",
        default="runs/orca_candidates.csv",
        help="Candidate CSV created by orca_select_candidates.py",
    )
    p.add_argument(
        "--sdf",
        default="moe_conformers/cnph_th_cf3_3d_conformers.sdf",
        help="SDF source used in reoptimize workflow",
    )
    p.add_argument(
        "--output-dir",
        default="runs/orca_jobs",
        help="Output directory for ORCA job folders",
    )
    p.add_argument(
        "--template-dir",
        default="orca_templates",
        help="Directory with ORCA input templates",
    )
    p.add_argument(
        "--stages",
        default="vertical,s1opt",
        help="Comma-separated stages "
        "(vertical,s1opt,s0freq,s1freq,esd_abs,esd_fluor)",
    )
    p.add_argument("--charge", type=int, default=0)
    p.add_argument("--multiplicity", type=int, default=1)
    p.add_argument("--nprocs", type=int, default=12)
    p.add_argument("--maxcore", type=int, default=3000, help="MB per core")
    p.add_argument("--functional", default="CAM-B3LYP")
    p.add_argument("--basis", default="def2-SVP")
    p.add_argument("--dispersion", default="D3BJ")
    p.add_argument(
        "--solvent",
        default="toluene",
        help="CPCM solvent name. Use empty string to disable CPCM.",
    )
    p.add_argument("--nroots", type=int, default=10)
    p.add_argument("--iroot", type=int, default=1, help="Target excited root")
    p.add_argument("--tda", action="store_true", help="Enable TDA in TDDFT blocks")
    p.add_argument(
        "--states",
        default="1",
        help="ESD STATES value, e.g. '1' or '1,2,3'",
    )
    p.add_argument("--linew", type=float, default=120.0)
    p.add_argument("--inlinew", type=float, default=800.0)
    return p


def main() -> int:
    args = build_parser().parse_args()

    candidates_path = Path(args.candidates_csv)
    sdf_path = Path(args.sdf)
    out_dir = Path(args.output_dir)
    template_dir = Path(args.template_dir)

    if not candidates_path.exists():
        raise SystemExit(f"Candidate CSV not found: {candidates_path}")
    if not sdf_path.exists():
        raise SystemExit(f"SDF not found: {sdf_path}")

    stages = _iter_stages(args.stages)
    candidates = pd.read_csv(candidates_path)
    if "conformer_id" not in candidates.columns:
        raise SystemExit("Candidate CSV must include conformer_id column")
    if candidates.empty:
        raise SystemExit("Candidate CSV has no rows")

    mols = _load_sdf(sdf_path)
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest = prepare_jobs(
        candidates=candidates,
        mols=mols,
        output_dir=out_dir,
        template_dir=template_dir,
        stages=stages,
        charge=args.charge,
        multiplicity=args.multiplicity,
        nprocs=args.nprocs,
        maxcore=args.maxcore,
        functional=args.functional,
        basis=args.basis,
        dispersion=args.dispersion,
        solvent=args.solvent.strip(),
        nroots=args.nroots,
        iroot=args.iroot,
        tda=bool(args.tda),
        states=args.states,
        linew=args.linew,
        inlinew=args.inlinew,
    )

    manifest_path = out_dir / "orca_job_manifest.csv"
    manifest.to_csv(manifest_path, index=False, encoding="utf-8-sig")
    print(f"Wrote {manifest_path} with {len(manifest)} jobs")
    print(
        "Candidates:",
        manifest["candidate_id"].nunique(),
        "| stages:",
        sorted(manifest["stage"].unique()),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
