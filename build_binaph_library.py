"""Build binaphthalene dimer SMILES library: 4 luminescent R-groups x 16 screening groups.

Luminescent R-groups go at outer ring positions (R_A=1, R_B=12).
Screening groups go at inner ring NaphB position (SCREEN_B=7).

Outputs binaph_dimer_smiles.csv for MOE import and RDKit screening.
"""

import csv
import sys
from rdkit import Chem
from rdkit.Chem import Descriptors

from build_binaph_dimer import build_dimer

NAPH_SMARTS = Chem.MolFromSmarts("c1cccc2ccccc12")

# 4 luminescent R-groups (from user's hand-drawn "Luminiscencni slozka" image)
LUMINESCENT_R = {
    "Pyr": "c1cc2ccc3cccc4ccc(c1)c2c34",           # pyrenyl
    "EtynPyr": "C#Cc1cc2ccc3cccc4ccc(c1)c2c34",    # ethynyl-pyrenyl
    "DCV_Th": "c1cc(sc1)C=C(C#N)C#N",              # dicyanovinyl-thiophene
    "CNPh_Th": "c1cc(sc1)C=C(C#N)c2ccc(C#N)cc2",  # cyano-phenyl-CN-thiophene
}

# 16 screening groups (at MeO position, inner ring NaphB)
SCREENING_GROUPS = {
    "H": "",
    "Me": "C",
    "Et": "CC",
    "nPr": "CCC",
    "iPr": "C(C)C",
    "nBu": "CCCC",
    "tBu": "C(C)(C)C",
    "cHex": "C1CCCCC1",
    "MeO": "OC",
    "OEt": "OCC",
    "F": "F",
    "Cl": "Cl",
    "CF3": "C(F)(F)F",
    "CN": "C#N",
    "NMe2": "N(C)C",
    "Ph": "c1ccccc1",
}

OUTPUT_CSV = "binaph_dimer_smiles.csv"


def main():
    print("=" * 70)
    print("BUILDING BINAPHTHALENE DIMER LIBRARY")
    print(f"  {len(LUMINESCENT_R)} luminescent R-groups x {len(SCREENING_GROUPS)} screening groups")
    print(f"  = {len(LUMINESCENT_R) * len(SCREENING_GROUPS)} total dimers")
    print("=" * 70)

    rows = []
    failures = []

    for r_name, r_smi in LUMINESCENT_R.items():
        print(f"\n--- R-group: {r_name} ---")

        for s_name, s_smi in SCREENING_GROUPS.items():
            name = f"{r_name}_{s_name}"
            mol = build_dimer(r_smi, s_smi, r_name, s_name)

            if mol is None:
                print(f"  {name:20s}: FAILED")
                failures.append(name)
                continue

            smi = Chem.MolToSmiles(mol)
            n_atoms = mol.GetNumAtoms()
            mw = Descriptors.MolWt(mol)
            n_aromatic = sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())
            n_naph = len(mol.GetSubstructMatches(NAPH_SMARTS))

            # Chiral centers
            chiral_info = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
            n_chiral = len(chiral_info)

            # Basic validation
            issues = []
            if n_naph < 2:
                issues.append(f"naph={n_naph}<2")

            status = "OK" if not issues else f"WARN({','.join(issues)})"

            print(f"  {name:20s}: atoms={n_atoms:3d}, MW={mw:7.1f}, "
                  f"arom={n_aromatic:3d}, naph={n_naph}, "
                  f"chiral={n_chiral}  [{status}]")

            rows.append({
                "name": name,
                "r_group": r_name,
                "screen_group": s_name,
                "smiles": smi,
                "n_atoms": n_atoms,
                "mw": round(mw, 1),
                "n_aromatic": n_aromatic,
                "n_naph_matches": n_naph,
                "n_chiral": n_chiral,
                "status": status,
            })

    # Write CSV
    print(f"\n{'=' * 70}")
    print(f"RESULTS: {len(rows)} successful, {len(failures)} failed")
    print(f"{'=' * 70}")

    if rows:
        fieldnames = ["name", "r_group", "screen_group", "smiles", "n_atoms",
                       "mw", "n_aromatic", "n_naph_matches", "n_chiral", "status"]
        with open(OUTPUT_CSV, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        print(f"\nSaved {len(rows)} dimers to {OUTPUT_CSV}")

    if failures:
        print(f"\nFailed: {', '.join(failures)}")

    # Summary by R-group
    print(f"\n{'=' * 70}")
    print("SUMMARY BY R-GROUP")
    print(f"{'=' * 70}")
    for r_name in LUMINESCENT_R:
        r_rows = [r for r in rows if r["r_group"] == r_name]
        if r_rows:
            mw_range = f"{min(r['mw'] for r in r_rows):.0f}-{max(r['mw'] for r in r_rows):.0f}"
            atom_range = f"{min(r['n_atoms'] for r in r_rows)}-{max(r['n_atoms'] for r in r_rows)}"
            print(f"  {r_name:12s}: {len(r_rows):2d} variants, "
                  f"atoms={atom_range}, MW={mw_range}")


if __name__ == "__main__":
    main()
