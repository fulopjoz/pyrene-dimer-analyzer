"""Fix binaphthalene substituent positions.

Issue: In the old assignment, R and MeO were on opposite ring types
for the upper vs lower naphthalene:
- Upper naph: R on outer ring (correct), bridge on inner ring (correct)
- Lower naph: R on INNER ring (WRONG), MeO on OUTER ring (WRONG)

Fix: Both R groups on outer rings (far from biaryl bond),
both bridge/MeO on inner rings (near biaryl bond).

Binaphthalene SMILES: c1ccc2c(-c3cccc4ccccc34)cccc2c1

Atom numbering:
- Biaryl bond: atom 4 -- atom 5
- Naph A (atoms 0-4, 15-19): connection at atom 4
  - Inner ring (has atom 4): {3, 4, 15, 16, 17, 18}
  - Outer ring: {0, 1, 2, 3, 18, 19}
  - H-bearing inner: 15, 16, 17
  - H-bearing outer: 0, 1, 2, 19
- Naph B (atoms 5-14): connection at atom 5
  - Inner ring (has atom 5): {5, 6, 7, 8, 9, 14}
  - Outer ring: {9, 10, 11, 12, 13, 14}
  - H-bearing inner: 6, 7, 8
  - H-bearing outer: 10, 11, 12, 13
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors
from pathlib import Path

outdir = Path("plots/structure_comparison")
outdir.mkdir(parents=True, exist_ok=True)

BINAPH_SMI = "c1ccc2c(-c3cccc4ccccc34)cccc2c1"


def add_fragment(rw, attach_idx, frag_smi):
    frag = Chem.MolFromSmiles(frag_smi)
    amap = {}
    for a in frag.GetAtoms():
        new_a = Chem.Atom(a.GetAtomicNum())
        new_a.SetIsAromatic(a.GetIsAromatic())
        new_a.SetFormalCharge(a.GetFormalCharge())
        amap[a.GetIdx()] = rw.AddAtom(new_a)
    for b in frag.GetBonds():
        rw.AddBond(amap[b.GetBeginAtomIdx()], amap[b.GetEndAtomIdx()], b.GetBondType())
    rw.AddBond(attach_idx, amap[0], Chem.BondType.SINGLE)
    return amap


def add_binaphthalene(rw):
    binaph = Chem.MolFromSmiles(BINAPH_SMI)
    amap = {}
    for a in binaph.GetAtoms():
        new_a = Chem.Atom(a.GetAtomicNum())
        new_a.SetIsAromatic(True)
        amap[a.GetIdx()] = rw.AddAtom(new_a)
    for b in binaph.GetBonds():
        rw.AddBond(amap[b.GetBeginAtomIdx()], amap[b.GetEndAtomIdx()], b.GetBondType())
    return amap


# ============================================================
# Try multiple position assignments
# ============================================================

# Position sets: (label, R_naphA, bridge_naphA, R_naphB, MeO_naphB)
# NaphA = atoms 0-4,15-19 (connection at 4)
# NaphB = atoms 5-14 (connection at 5)
# User: bridge on upper naph, MeO on lower naph

# The question is which naph is "upper" and which is "lower"
# Both are valid - depends on 3D orientation

POSITION_SETS = {
    "OLD": {
        "desc": "OLD (current): R mixed rings, bridge@7, MeO@1",
        "R_A": 16,   # inner ring of naph A
        "bridge_A": None,
        "MeO_A": 1,  # outer ring of naph A
        "R_B": 12,   # outer ring of naph B
        "bridge_B": 7,  # inner ring of naph B
        "MeO_B": None,
    },
    "FIX1": {
        "desc": "FIX1: Both R on outer rings. NaphA=upper(bridge), NaphB=lower(MeO)",
        "R_A": 1,    # outer ring of naph A
        "bridge_A": 16,  # inner ring of naph A
        "MeO_A": None,
        "R_B": 12,   # outer ring of naph B
        "bridge_B": None,
        "MeO_B": 7,  # inner ring of naph B
    },
    "FIX2": {
        "desc": "FIX2: Both R on outer rings. NaphA=lower(MeO), NaphB=upper(bridge)",
        "R_A": 1,    # outer ring of naph A
        "bridge_A": None,
        "MeO_A": 16,  # inner ring of naph A
        "R_B": 12,   # outer ring of naph B
        "bridge_B": 7,  # inner ring of naph B
        "MeO_B": None,
    },
    "FIX3": {
        "desc": "FIX3: R@19,10 (near junction). NaphA=upper(bridge@15), NaphB=lower(MeO@6)",
        "R_A": 19,   # outer ring, near junction
        "bridge_A": 15,  # inner ring, near biaryl bond
        "MeO_A": None,
        "R_B": 10,   # outer ring, near junction
        "bridge_B": None,
        "MeO_B": 6,  # inner ring, near biaryl bond
    },
}


def build_h_dimer(pos_set):
    """Build H dimer with given position set."""
    rw = Chem.RWMol(Chem.MolFromSmiles(BINAPH_SMI))

    # MeO on naph A or B
    if pos_set.get("MeO_A"):
        o = rw.AddAtom(Chem.Atom(8))
        c = rw.AddAtom(Chem.Atom(6))
        rw.AddBond(pos_set["MeO_A"], o, Chem.BondType.SINGLE)
        rw.AddBond(o, c, Chem.BondType.SINGLE)

    if pos_set.get("MeO_B"):
        o = rw.AddAtom(Chem.Atom(8))
        c = rw.AddAtom(Chem.Atom(6))
        rw.AddBond(pos_set["MeO_B"], o, Chem.BondType.SINGLE)
        rw.AddBond(o, c, Chem.BondType.SINGLE)

    # Bridge from naph A or B: O-CH(Me)-O
    bridge_pos = pos_set.get("bridge_A") or pos_set.get("bridge_B")
    o1 = rw.AddAtom(Chem.Atom(8))
    rw.AddBond(bridge_pos, o1, Chem.BondType.SINGLE)
    c_chiral1 = rw.AddAtom(Chem.Atom(6))
    rw.AddBond(o1, c_chiral1, Chem.BondType.SINGLE)
    me1 = rw.AddAtom(Chem.Atom(6))
    rw.AddBond(c_chiral1, me1, Chem.BondType.SINGLE)
    o2 = rw.AddAtom(Chem.Atom(8))
    rw.AddBond(c_chiral1, o2, Chem.BondType.SINGLE)

    # Central phenyl (1,3,5-trisubstituted)
    ph = []
    for i in range(6):
        a = Chem.Atom(6)
        a.SetIsAromatic(True)
        ph.append(rw.AddAtom(a))
    for i in range(6):
        rw.AddBond(ph[i], ph[(i + 1) % 6], Chem.BondType.AROMATIC)
    rw.AddBond(o2, ph[0], Chem.BondType.SINGLE)

    # Aminophosphonate
    c_ap = rw.AddAtom(Chem.Atom(6))
    rw.AddBond(ph[2], c_ap, Chem.BondType.SINGLE)
    n = rw.AddAtom(Chem.Atom(7))
    rw.AddBond(c_ap, n, Chem.BondType.SINGLE)
    p = rw.AddAtom(Chem.Atom(15))
    rw.AddBond(c_ap, p, Chem.BondType.SINGLE)
    rw.AddBond(p, rw.AddAtom(Chem.Atom(8)), Chem.BondType.DOUBLE)
    rw.AddBond(p, rw.AddAtom(Chem.Atom(8)), Chem.BondType.SINGLE)
    rw.AddBond(p, rw.AddAtom(Chem.Atom(8)), Chem.BondType.SINGLE)

    # Bridge 2: O-CH(Me)-O
    o3 = rw.AddAtom(Chem.Atom(8))
    rw.AddBond(ph[4], o3, Chem.BondType.SINGLE)
    c_chiral2 = rw.AddAtom(Chem.Atom(6))
    rw.AddBond(o3, c_chiral2, Chem.BondType.SINGLE)
    me2 = rw.AddAtom(Chem.Atom(6))
    rw.AddBond(c_chiral2, me2, Chem.BondType.SINGLE)
    o4 = rw.AddAtom(Chem.Atom(8))
    rw.AddBond(c_chiral2, o4, Chem.BondType.SINGLE)

    # Binaphthalene 2
    amap2 = add_binaphthalene(rw)

    # Connect bridge to binaph 2
    bridge2_pos = pos_set.get("bridge_A") or pos_set.get("bridge_B")
    # For binaph 2, bridge connects at same position type as binaph 1
    if pos_set.get("bridge_A"):
        rw.AddBond(o4, amap2[pos_set["bridge_A"]], Chem.BondType.SINGLE)
    else:
        rw.AddBond(o4, amap2[pos_set["bridge_B"]], Chem.BondType.SINGLE)

    # MeO on binaph 2
    if pos_set.get("MeO_A"):
        o = rw.AddAtom(Chem.Atom(8))
        c = rw.AddAtom(Chem.Atom(6))
        rw.AddBond(amap2[pos_set["MeO_A"]], o, Chem.BondType.SINGLE)
        rw.AddBond(o, c, Chem.BondType.SINGLE)
    if pos_set.get("MeO_B"):
        o = rw.AddAtom(Chem.Atom(8))
        c = rw.AddAtom(Chem.Atom(6))
        rw.AddBond(amap2[pos_set["MeO_B"]], o, Chem.BondType.SINGLE)
        rw.AddBond(o, c, Chem.BondType.SINGLE)

    try:
        mol = rw.GetMol()
        Chem.SanitizeMol(mol)
        return mol
    except Exception as e:
        print(f"  ERROR: {e}")
        return None


print("=" * 70)
print("COMPARING POSITION ASSIGNMENTS (H dimer)")
print("=" * 70)

for label, pos_set in POSITION_SETS.items():
    print(f"\n{label}: {pos_set['desc']}")
    mol = build_h_dimer(pos_set)
    if mol:
        smi = Chem.MolToSmiles(mol)
        print(f"  SMILES: {smi}")
        print(f"  Atoms: {mol.GetNumAtoms()}, MW: {Descriptors.MolWt(mol):.1f}")
        AllChem.Compute2DCoords(mol)
        img = Draw.MolToImage(mol, size=(800, 600))
        img.save(str(outdir / f"binaph_H_dimer_{label}.png"))
        print(f"  Saved binaph_H_dimer_{label}.png")
    else:
        print("  FAILED")
