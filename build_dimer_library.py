"""Build corrected pyrene dimer SMILES for all R-group variants.

Uses RWMol to programmatically construct each dimer with:
- Two pyrene units (peri-fused, 16C each)
- R-groups at positions 6,8 on each pyrene
- Bridge: Pyrene-O-CH(Me)-O-CentralPhenyl-O-CH(Me)-O-Pyrene
- Aminophosphonate on central phenyl
"""

import csv
import sys
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

PYRENE_SMARTS = Chem.MolFromSmarts("c1cc2ccc3cccc4ccc(c1)c2c34")

R_GROUPS = {
    "H": "",
    "Me": "C",
    "Et": "CC",
    "nPr": "CCC",
    "iPr": "C(C)C",
    "nBu": "CCCC",
    "iBu": "CC(C)C",
    "sBu": "C(C)CC",
    "tBu": "C(C)(C)C",
    "nPent": "CCCCC",
    "nHex": "CCCCCC",
    "neoPent": "CC(C)(C)C",
    "cPent": "C1CCCC1",
    "cHex": "C1CCCCC1",
    "Ph": "c1ccccc1",
    "4MePh": "c1ccc(C)cc1",
    "4FPh": "c1ccc(F)cc1",
    "4ClPh": "c1ccc(Cl)cc1",
    "4OMePh": "c1ccc(OC)cc1",
    "4CF3Ph": "c1ccc(C(F)(F)F)cc1",
    "OMe": "OC",
    "OEt": "OCC",
    "NMe2": "N(C)C",
    "F": "F",
    "Cl": "Cl",
    "CF3": "C(F)(F)F",
    "CN": "C#N",
    "CH2OH": "CO",
    "CH2OMe": "COC",
    "Bn": "Cc1ccccc1",
}


def add_r_group_to_rwmol(rw, attach_idx, r_smi):
    """Add an R-group fragment to an RWMol at the specified atom index."""
    r_mol = Chem.MolFromSmiles(r_smi)
    if r_mol is None:
        raise ValueError(f"Invalid R-group SMILES: {r_smi}")
    amap = {}
    for a in r_mol.GetAtoms():
        new_a = Chem.Atom(a.GetAtomicNum())
        new_a.SetIsAromatic(a.GetIsAromatic())
        new_a.SetFormalCharge(a.GetFormalCharge())
        new_a.SetNumExplicitHs(a.GetNumExplicitHs())
        amap[a.GetIdx()] = rw.AddAtom(new_a)
    for b in r_mol.GetBonds():
        rw.AddBond(
            amap[b.GetBeginAtomIdx()],
            amap[b.GetEndAtomIdx()],
            b.GetBondType(),
        )
    rw.AddBond(attach_idx, amap[0], Chem.BondType.SINGLE)
    return amap


def add_pyrene_to_rwmol(rw):
    """Add a pyrene unit to RWMol and return atom index mapping."""
    pyrene = Chem.MolFromSmiles("c1cc2ccc3cccc4ccc(c1)c2c34")
    amap = {}
    for a in pyrene.GetAtoms():
        new_a = Chem.Atom(a.GetAtomicNum())
        new_a.SetIsAromatic(True)
        amap[a.GetIdx()] = rw.AddAtom(new_a)
    for b in pyrene.GetBonds():
        rw.AddBond(
            amap[b.GetBeginAtomIdx()],
            amap[b.GetEndAtomIdx()],
            b.GetBondType(),
        )
    return amap


def build_dimer(r_smi, r_name=""):
    """Build a full pyrene dimer with specified R-groups.

    Structure: Pyrene1(R,R)-O-CH(Me)-O-Phenyl(NH-CH-PO3H2)-O-CH(Me)-O-Pyrene2(R,R)
    R-groups at pyrene positions 6 and 8 (Ring C).
    Bridge attached at pyrene position 0 (Ring A).
    """
    # Start with pyrene 1
    pyrene = Chem.MolFromSmiles("c1cc2ccc3cccc4ccc(c1)c2c34")
    rw = Chem.RWMol(pyrene)

    # Add R-groups to pyrene 1 at positions 6 and 8
    if r_smi:
        for pos in [6, 8]:
            add_r_group_to_rwmol(rw, pos, r_smi)

    # Bridge 1: O-CH(Me)-O from pyrene1 position 0
    o1 = rw.AddAtom(Chem.Atom(8))  # oxygen
    rw.AddBond(0, o1, Chem.BondType.SINGLE)

    c_chiral1 = rw.AddAtom(Chem.Atom(6))  # chiral carbon
    rw.AddBond(o1, c_chiral1, Chem.BondType.SINGLE)

    me1 = rw.AddAtom(Chem.Atom(6))  # methyl
    rw.AddBond(c_chiral1, me1, Chem.BondType.SINGLE)

    o2 = rw.AddAtom(Chem.Atom(8))  # oxygen to phenyl
    rw.AddBond(c_chiral1, o2, Chem.BondType.SINGLE)

    # Central phenyl ring (1,3,5-trisubstituted)
    ph_atoms = []
    for i in range(6):
        a = Chem.Atom(6)
        a.SetIsAromatic(True)
        ph_atoms.append(rw.AddAtom(a))
    for i in range(6):
        rw.AddBond(ph_atoms[i], ph_atoms[(i + 1) % 6], Chem.BondType.AROMATIC)

    # Connect bridge 1 to phenyl position 0
    rw.AddBond(o2, ph_atoms[0], Chem.BondType.SINGLE)

    # Aminophosphonate at phenyl position 2 (1,3,5 pattern)
    c_chiral3 = rw.AddAtom(Chem.Atom(6))  # chiral C
    rw.AddBond(ph_atoms[2], c_chiral3, Chem.BondType.SINGLE)

    n_nh2 = rw.AddAtom(Chem.Atom(7))  # NH2
    rw.AddBond(c_chiral3, n_nh2, Chem.BondType.SINGLE)

    p_atom = rw.AddAtom(Chem.Atom(15))  # P
    rw.AddBond(c_chiral3, p_atom, Chem.BondType.SINGLE)

    o_po1 = rw.AddAtom(Chem.Atom(8))  # P=O
    rw.AddBond(p_atom, o_po1, Chem.BondType.DOUBLE)

    o_po2 = rw.AddAtom(Chem.Atom(8))  # P-OH
    rw.AddBond(p_atom, o_po2, Chem.BondType.SINGLE)

    o_po3 = rw.AddAtom(Chem.Atom(8))  # P-OH
    rw.AddBond(p_atom, o_po3, Chem.BondType.SINGLE)

    # Bridge 2: O-CH(Me)-O from phenyl position 4 (1,3,5 pattern)
    o3 = rw.AddAtom(Chem.Atom(8))
    rw.AddBond(ph_atoms[4], o3, Chem.BondType.SINGLE)

    c_chiral2 = rw.AddAtom(Chem.Atom(6))
    rw.AddBond(o3, c_chiral2, Chem.BondType.SINGLE)

    me2 = rw.AddAtom(Chem.Atom(6))
    rw.AddBond(c_chiral2, me2, Chem.BondType.SINGLE)

    o4 = rw.AddAtom(Chem.Atom(8))
    rw.AddBond(c_chiral2, o4, Chem.BondType.SINGLE)

    # Pyrene 2
    amap2 = add_pyrene_to_rwmol(rw)
    rw.AddBond(o4, amap2[0], Chem.BondType.SINGLE)

    # Add R-groups to pyrene 2 at positions 6 and 8
    if r_smi:
        for pos in [6, 8]:
            add_r_group_to_rwmol(rw, amap2[pos], r_smi)

    # Sanitize
    try:
        mol = rw.GetMol()
        Chem.SanitizeMol(mol)
        return mol
    except Exception as e:
        print(f"  ERROR sanitizing {r_name}: {e}", file=sys.stderr)
        return None


def validate_dimer(mol, r_name):
    """Validate a dimer molecule."""
    if mol is None:
        return False, "None molecule"

    n_atoms = mol.GetNumAtoms()
    mw = Descriptors.MolWt(mol)
    n_arom = sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())

    # Check for pyrene substructures
    matches = mol.GetSubstructMatches(PYRENE_SMARTS)
    n_pyrene = len(matches)

    # Check for chiral centers
    chiral_info = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    n_chiral = len(chiral_info)

    smi = Chem.MolToSmiles(mol)

    valid = n_pyrene >= 2 and n_arom >= 38
    status = "OK" if valid else "WARN"

    return valid, {
        "name": r_name,
        "smiles": smi,
        "n_atoms": n_atoms,
        "mw": round(mw, 1),
        "n_aromatic": n_arom,
        "n_pyrene_matches": n_pyrene,
        "n_chiral": n_chiral,
        "status": status,
    }


def main():
    print("=" * 70)
    print("BUILDING CORRECTED PYRENE DIMER LIBRARY")
    print("=" * 70)

    results = []
    failures = []

    for r_name, r_smi in R_GROUPS.items():
        print(f"\n  Building {r_name} dimer (R = {r_smi if r_smi else 'H'})...")
        mol = build_dimer(r_smi, r_name)

        if mol is None:
            failures.append(r_name)
            print(f"    FAILED: could not build molecule")
            continue

        valid, info = validate_dimer(mol, r_name)

        if not valid:
            failures.append(r_name)
            print(f"    WARN: {info}")
        else:
            print(f"    atoms={info['n_atoms']}, MW={info['mw']}, "
                  f"aromatic={info['n_aromatic']}, pyrenes={info['n_pyrene_matches']}, "
                  f"chiral={info['n_chiral']}")
            results.append(info)

    # Save to CSV
    output_file = "corrected_dimer_smiles.csv"
    with open(output_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "name", "smiles", "n_atoms", "mw", "n_aromatic",
            "n_pyrene_matches", "n_chiral", "status",
        ])
        writer.writeheader()
        writer.writerows(results)

    print(f"\n{'=' * 70}")
    print(f"RESULTS: {len(results)} valid dimers, {len(failures)} failures")
    if failures:
        print(f"Failed: {', '.join(failures)}")
    print(f"Saved to: {output_file}")
    print(f"{'=' * 70}")

    # Also print a summary table
    print(f"\n{'Name':<12} {'Atoms':>5} {'MW':>8} {'Arom':>5} {'Pyr':>4} {'Chi':>4}")
    print("-" * 42)
    for r in results:
        print(f"{r['name']:<12} {r['n_atoms']:>5} {r['mw']:>8.1f} "
              f"{r['n_aromatic']:>5} {r['n_pyrene_matches']:>4} {r['n_chiral']:>4}")


if __name__ == "__main__":
    main()
