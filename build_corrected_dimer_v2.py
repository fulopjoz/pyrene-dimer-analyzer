"""Build corrected pyrene dimer with MeO on each pyrene unit.

Monomer structure (from hand-drawn image):
- Pyrene core (two naphthalene units joined by peri-bond |)
- Upper naphthalene: R (right, atom 13) + bridge-O (left, atom 8)
- Lower naphthalene: MeO (left, atom 6) + R (right, atom 1)

Pyrene positions:
         10 --- 11
        /         \\
       8(bridge)   13(R)
      /   14---15   \\
     7    |     |    0
      \\   15---14   /
       6(OMe)      1(R)
        \\         /
         4 --- 3
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors
from PIL import Image, ImageDraw
from pathlib import Path

outdir = Path("plots/structure_comparison")
outdir.mkdir(parents=True, exist_ok=True)

PYRENE_SMARTS = Chem.MolFromSmarts("c1cc2ccc3cccc4ccc(c1)c2c34")

# Positions from the hand-drawn image:
R_POS_1 = 13  # upper-right (upper naphthalene)
R_POS_2 = 1   # lower-right (lower naphthalene)
BRIDGE_POS = 8  # upper-left (upper naphthalene)
OME_POS = 6     # lower-left (lower naphthalene)


def add_fragment(rw, attach_idx, frag_smi):
    """Add a SMILES fragment to an RWMol at the specified atom."""
    frag = Chem.MolFromSmiles(frag_smi)
    if frag is None:
        raise ValueError(f"Invalid fragment SMILES: {frag_smi}")
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


def add_pyrene(rw):
    """Add a pyrene unit and return atom index mapping."""
    pyrene = Chem.MolFromSmiles("c1cc2ccc3cccc4ccc(c1)c2c34")
    amap = {}
    for a in pyrene.GetAtoms():
        new_a = Chem.Atom(a.GetAtomicNum())
        new_a.SetIsAromatic(True)
        amap[a.GetIdx()] = rw.AddAtom(new_a)
    for b in pyrene.GetBonds():
        rw.AddBond(amap[b.GetBeginAtomIdx()], amap[b.GetEndAtomIdx()], b.GetBondType())
    return amap


def build_dimer_v2(r_smi, r_name=""):
    """Build corrected pyrene dimer with MeO on each pyrene unit.

    Each pyrene unit has:
    - R at positions 13 and 1 (right side)
    - OMe at position 6 (lower-left)
    - O-bridge at position 8 (upper-left)

    Bridge: Pyrene1-O(8)-CH(Me)-O-CentralPhenyl(aminophosphonate)-O-CH(Me)-O(8)-Pyrene2
    """
    # Start with pyrene 1
    pyrene = Chem.MolFromSmiles("c1cc2ccc3cccc4ccc(c1)c2c34")
    rw = Chem.RWMol(pyrene)

    # Add R-groups at positions 13 and 1 (right side of pyrene)
    if r_smi:
        for pos in [R_POS_1, R_POS_2]:
            add_fragment(rw, pos, r_smi)

    # Add OMe at position 6 (lower-left)
    o_ome = rw.AddAtom(Chem.Atom(8))
    c_ome = rw.AddAtom(Chem.Atom(6))
    rw.AddBond(OME_POS, o_ome, Chem.BondType.SINGLE)
    rw.AddBond(o_ome, c_ome, Chem.BondType.SINGLE)

    # Bridge from position 8 (upper-left): O-CH(Me)-O
    o1 = rw.AddAtom(Chem.Atom(8))
    rw.AddBond(BRIDGE_POS, o1, Chem.BondType.SINGLE)
    c_chiral1 = rw.AddAtom(Chem.Atom(6))
    rw.AddBond(o1, c_chiral1, Chem.BondType.SINGLE)
    me1 = rw.AddAtom(Chem.Atom(6))
    rw.AddBond(c_chiral1, me1, Chem.BondType.SINGLE)
    o2 = rw.AddAtom(Chem.Atom(8))
    rw.AddBond(c_chiral1, o2, Chem.BondType.SINGLE)

    # Central phenyl ring (1,3,5-trisubstituted)
    ph = []
    for i in range(6):
        a = Chem.Atom(6)
        a.SetIsAromatic(True)
        ph.append(rw.AddAtom(a))
    for i in range(6):
        rw.AddBond(ph[i], ph[(i + 1) % 6], Chem.BondType.AROMATIC)
    rw.AddBond(o2, ph[0], Chem.BondType.SINGLE)

    # Aminophosphonate at phenyl position 2
    c_ap = rw.AddAtom(Chem.Atom(6))
    rw.AddBond(ph[2], c_ap, Chem.BondType.SINGLE)
    n = rw.AddAtom(Chem.Atom(7))
    rw.AddBond(c_ap, n, Chem.BondType.SINGLE)
    p = rw.AddAtom(Chem.Atom(15))
    rw.AddBond(c_ap, p, Chem.BondType.SINGLE)
    o_po1 = rw.AddAtom(Chem.Atom(8))
    rw.AddBond(p, o_po1, Chem.BondType.DOUBLE)
    o_po2 = rw.AddAtom(Chem.Atom(8))
    rw.AddBond(p, o_po2, Chem.BondType.SINGLE)
    o_po3 = rw.AddAtom(Chem.Atom(8))
    rw.AddBond(p, o_po3, Chem.BondType.SINGLE)

    # Bridge 2 from phenyl position 4: O-CH(Me)-O
    o3 = rw.AddAtom(Chem.Atom(8))
    rw.AddBond(ph[4], o3, Chem.BondType.SINGLE)
    c_chiral2 = rw.AddAtom(Chem.Atom(6))
    rw.AddBond(o3, c_chiral2, Chem.BondType.SINGLE)
    me2 = rw.AddAtom(Chem.Atom(6))
    rw.AddBond(c_chiral2, me2, Chem.BondType.SINGLE)
    o4 = rw.AddAtom(Chem.Atom(8))
    rw.AddBond(c_chiral2, o4, Chem.BondType.SINGLE)

    # Pyrene 2
    amap2 = add_pyrene(rw)
    rw.AddBond(o4, amap2[BRIDGE_POS], Chem.BondType.SINGLE)

    # Add R-groups to pyrene 2
    if r_smi:
        for pos in [R_POS_1, R_POS_2]:
            add_fragment(rw, amap2[pos], r_smi)

    # Add OMe to pyrene 2
    o_ome2 = rw.AddAtom(Chem.Atom(8))
    c_ome2 = rw.AddAtom(Chem.Atom(6))
    rw.AddBond(amap2[OME_POS], o_ome2, Chem.BondType.SINGLE)
    rw.AddBond(o_ome2, c_ome2, Chem.BondType.SINGLE)

    try:
        mol = rw.GetMol()
        Chem.SanitizeMol(mol)
        return mol
    except Exception as e:
        print(f"  ERROR: {e}")
        return None


# ================================================================
# Build and validate
# ================================================================
print("=" * 70)
print("BUILDING CORRECTED DIMER V2 (with MeO)")
print("=" * 70)

# First, build and render the MONOMER for comparison
print("\n--- Et Monomer (for comparison with hand-drawn image) ---")
rw = Chem.RWMol(Chem.MolFromSmiles("c1cc2ccc3cccc4ccc(c1)c2c34"))
# Et at 13
c1 = rw.AddAtom(Chem.Atom(6)); c2 = rw.AddAtom(Chem.Atom(6))
rw.AddBond(13, c1, Chem.BondType.SINGLE); rw.AddBond(c1, c2, Chem.BondType.SINGLE)
# Et at 1
c3 = rw.AddAtom(Chem.Atom(6)); c4 = rw.AddAtom(Chem.Atom(6))
rw.AddBond(1, c3, Chem.BondType.SINGLE); rw.AddBond(c3, c4, Chem.BondType.SINGLE)
# OMe at 6
o = rw.AddAtom(Chem.Atom(8)); c = rw.AddAtom(Chem.Atom(6))
rw.AddBond(6, o, Chem.BondType.SINGLE); rw.AddBond(o, c, Chem.BondType.SINGLE)
# OH (bridge point) at 8
oh = rw.AddAtom(Chem.Atom(8))
rw.AddBond(8, oh, Chem.BondType.SINGLE)
monomer = rw.GetMol()
Chem.SanitizeMol(monomer)
mono_smi = Chem.MolToSmiles(monomer)
print(f"  Monomer SMILES: {mono_smi}")
AllChem.Compute2DCoords(monomer)
img = Draw.MolToImage(monomer, size=(600, 600))
img.save(str(outdir / "CORRECTED_Et_monomer_v2.png"))
print("  Saved CORRECTED_Et_monomer_v2.png")

# Build dimers for key R-groups
r_groups = {
    "H": "",
    "Me": "C",
    "Et": "CC",
    "iPr": "C(C)C",
    "tBu": "C(C)(C)C",
    "cHex": "C1CCCCC1",
    "Cl": "Cl",
    "CF3": "C(F)(F)F",
}

print("\n--- Building dimers ---")
for name, r_smi in r_groups.items():
    mol = build_dimer_v2(r_smi, name)
    if mol is None:
        print(f"  {name}: FAILED")
        continue

    n_atoms = mol.GetNumAtoms()
    mw = Descriptors.MolWt(mol)
    n_arom = sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())
    pyr = mol.GetSubstructMatches(PYRENE_SMARTS)
    ome = mol.GetSubstructMatches(Chem.MolFromSmarts("[OX2][CH3]"))
    smi = Chem.MolToSmiles(mol)

    print(f"  {name:6s}: atoms={n_atoms}, MW={mw:.1f}, arom={n_arom}, "
          f"pyrene={len(pyr)}, OMe={len(ome)}")

    AllChem.Compute2DCoords(mol)
    img = Draw.MolToImage(mol, size=(800, 600))
    img.save(str(outdir / f"CORRECTED_v2_{name}_dimer.png"))

# Create comparison: old (no OMe) vs new (with OMe) for Et
print("\n--- Side-by-side comparison ---")
old_mol = Chem.MolFromSmiles("CCc1cc(CC)c2ccc3cc(OC(C)Oc4cc(OC(C)Oc5cc6ccc7c(CC)cc(CC)c8ccc(c5)c6c78)cc(C(N)P(=O)(O)O)c4)cc9ccc1c2c93")
new_mol = build_dimer_v2("CC", "Et")

if old_mol and new_mol:
    AllChem.Compute2DCoords(old_mol)
    AllChem.Compute2DCoords(new_mol)

    old_img = Draw.MolToImage(old_mol, size=(500, 400))
    new_img = Draw.MolToImage(new_mol, size=(500, 400))

    comparison = Image.new("RGB", (1000, 480), "white")
    comparison.paste(old_img, (0, 40))
    comparison.paste(new_img, (500, 40))
    draw = ImageDraw.Draw(comparison)
    draw.text((50, 5), "OLD: R@6,8, bridge@0, NO OMe", fill="red")
    draw.text((550, 5), "NEW: R@13,1, bridge@8, OMe@6", fill="green")
    comparison.save(str(outdir / "OLD_vs_NEW_Et_dimer.png"))
    print("  Saved OLD_vs_NEW_Et_dimer.png")

    old_smi = Chem.MolToSmiles(old_mol)
    new_smi = Chem.MolToSmiles(new_mol)
    old_ome = old_mol.GetSubstructMatches(Chem.MolFromSmarts("[OX2][CH3]"))
    new_ome = new_mol.GetSubstructMatches(Chem.MolFromSmarts("[OX2][CH3]"))
    print(f"  OLD: {len(old_ome)} OMe groups")
    print(f"  NEW: {len(new_ome)} OMe groups")
