"""Build 1,1'-binaphthalene dimers with luminescent R-groups and screening substituents.

Structure (FIX1 confirmed):

MONOMER (each half of the dimer):
- 1,1'-binaphthalene core (two naphthalene units connected by single C-C bond)
- NaphA: luminescent R on outer ring (atom 1), bridge-O on inner ring (atom 16)
- NaphB: luminescent R on outer ring (atom 12), screening group on inner ring (atom 7)

DIMER:
  Binaph1(R,screen) - O-CH(Me)-O - CentralPhenyl(aminophosphonate) - O-CH(Me)-O - Binaph2(R,screen)

Atom numbering for 1,1'-binaphthalene (SMILES: c1ccc2c(-c3cccc4ccccc34)cccc2c1):
- Biaryl bond: atom 4 (NaphA) -- atom 5 (NaphB)
- NaphA (atoms 0-4, 15-19): connection at atom 4
  - Inner ring: {3, 4, 15, 16, 17, 18}
  - Outer ring: {0, 1, 2, 3, 18, 19}
- NaphB (atoms 5-14): connection at atom 5
  - Inner ring: {5, 6, 7, 8, 9, 14}
  - Outer ring: {9, 10, 11, 12, 13, 14}

FIX1 position assignments:
- R_A: atom 1 (outer NaphA) - luminescent chromophore
- BRIDGE_A: atom 16 (inner NaphA) - bridge to central phenyl
- R_B: atom 12 (outer NaphB) - luminescent chromophore
- SCREEN_B: atom 7 (inner NaphB) - screening variable (MeO, Et, iPr, etc.)
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors
from pathlib import Path

outdir = Path("plots/structure_comparison")
outdir.mkdir(parents=True, exist_ok=True)

BINAPH_SMI = "c1ccc2c(-c3cccc4ccccc34)cccc2c1"

# FIX1 corrected positions
# NaphA (atoms 0-4, 15-19): R on outer, bridge on inner
R_A = 1          # outer ring NaphA (luminescent R-group)
BRIDGE_A = 16    # inner ring NaphA (bridge to central phenyl)
# NaphB (atoms 5-14): R on outer, screening group on inner
R_B = 12         # outer ring NaphB (luminescent R-group)
SCREEN_B = 7     # inner ring NaphB (screening variable, replaces MeO)


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


def add_binaphthalene(rw):
    """Add a 1,1'-binaphthalene unit and return atom index mapping."""
    binaph = Chem.MolFromSmiles(BINAPH_SMI)
    amap = {}
    for a in binaph.GetAtoms():
        new_a = Chem.Atom(a.GetAtomicNum())
        new_a.SetIsAromatic(True)
        amap[a.GetIdx()] = rw.AddAtom(new_a)
    for b in binaph.GetBonds():
        rw.AddBond(amap[b.GetBeginAtomIdx()], amap[b.GetEndAtomIdx()], b.GetBondType())
    return amap


def build_monomer(r_smi, screen_smi="OC", r_name="", screen_name="MeO"):
    """Build a single binaphthalene monomer with R, screening group, and bridge-OH."""
    rw = Chem.RWMol(Chem.MolFromSmiles(BINAPH_SMI))

    # Luminescent R on NaphA outer ring
    if r_smi:
        add_fragment(rw, R_A, r_smi)

    # Luminescent R on NaphB outer ring
    if r_smi:
        add_fragment(rw, R_B, r_smi)

    # Screening group on NaphB inner ring (replaces hard-coded MeO)
    if screen_smi:
        add_fragment(rw, SCREEN_B, screen_smi)

    # Bridge-OH on NaphA inner ring (attachment point for dimer bridge)
    o_br = rw.AddAtom(Chem.Atom(8))
    rw.AddBond(BRIDGE_A, o_br, Chem.BondType.SINGLE)

    try:
        mol = rw.GetMol()
        Chem.SanitizeMol(mol)
        return mol
    except Exception as e:
        print(f"  ERROR building monomer {r_name}/{screen_name}: {e}")
        return None


def build_dimer(r_smi, screen_smi="OC", r_name="", screen_name="MeO"):
    """Build the full dimer:
    Binaph1(R,screen) - O-CH(Me)-O - Phenyl(aminophosphonate) - O-CH(Me)-O - Binaph2(R,screen)
    """
    # Start with binaphthalene 1
    binaph = Chem.MolFromSmiles(BINAPH_SMI)
    rw = Chem.RWMol(binaph)

    # Luminescent R groups on binaphthalene 1
    if r_smi:
        add_fragment(rw, R_A, r_smi)
        add_fragment(rw, R_B, r_smi)

    # Screening group on binaphthalene 1 (NaphB inner ring)
    if screen_smi:
        add_fragment(rw, SCREEN_B, screen_smi)

    # Bridge from binaphthalene 1 (NaphA inner ring): O-CH(Me)-O
    o1 = rw.AddAtom(Chem.Atom(8))
    rw.AddBond(BRIDGE_A, o1, Chem.BondType.SINGLE)
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

    # Aminophosphonate at phenyl position 3 (1,3,5 pattern)
    c_ap = rw.AddAtom(Chem.Atom(6))  # CH
    rw.AddBond(ph[2], c_ap, Chem.BondType.SINGLE)
    n = rw.AddAtom(Chem.Atom(7))  # NH
    rw.AddBond(c_ap, n, Chem.BondType.SINGLE)
    p = rw.AddAtom(Chem.Atom(15))  # P
    rw.AddBond(c_ap, p, Chem.BondType.SINGLE)
    o_po1 = rw.AddAtom(Chem.Atom(8))  # P=O
    rw.AddBond(p, o_po1, Chem.BondType.DOUBLE)
    o_po2 = rw.AddAtom(Chem.Atom(8))  # P-OH
    rw.AddBond(p, o_po2, Chem.BondType.SINGLE)
    o_po3 = rw.AddAtom(Chem.Atom(8))  # P-OH
    rw.AddBond(p, o_po3, Chem.BondType.SINGLE)

    # Bridge to binaphthalene 2: O-CH(Me)-O
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
    rw.AddBond(o4, amap2[BRIDGE_A], Chem.BondType.SINGLE)

    # Luminescent R groups on binaphthalene 2
    if r_smi:
        add_fragment(rw, amap2[R_A], r_smi)
        add_fragment(rw, amap2[R_B], r_smi)

    # Screening group on binaphthalene 2 (NaphB inner ring)
    if screen_smi:
        add_fragment(rw, amap2[SCREEN_B], screen_smi)

    try:
        mol = rw.GetMol()
        Chem.SanitizeMol(mol)
        return mol
    except Exception as e:
        print(f"  ERROR building dimer {r_name}/{screen_name}: {e}")
        return None


# ================================================================
# Quick validation when run directly
# ================================================================
if __name__ == "__main__":
    NAPH_SMARTS = Chem.MolFromSmarts("c1cccc2ccccc12")

    print("=" * 70)
    print("BUILDING BINAPHTHALENE DIMER (FIX1 positions)")
    print("=" * 70)

    # Test with pyrenyl R-group and MeO screening
    test_r = "c1cc2ccc3cccc4ccc(c1)c2c34"  # pyrenyl
    mol = build_dimer(test_r, "OC", "Pyr", "MeO")
    if mol:
        smi = Chem.MolToSmiles(mol)
        n_atoms = mol.GetNumAtoms()
        mw = Descriptors.MolWt(mol)
        n_arom = sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())
        print(f"  Pyr/MeO: atoms={n_atoms}, MW={mw:.1f}, arom={n_arom}")
        print(f"  SMILES: {smi[:100]}...")
        AllChem.Compute2DCoords(mol)
        img = Draw.MolToImage(mol, size=(900, 700))
        img.save(str(outdir / "BINAPH_dimer_Pyr_MeO.png"))
        print("  Saved BINAPH_dimer_Pyr_MeO.png")
    else:
        print("  FAILED")

    # Test with H R-group and MeO screening (matches FIX1 from fix_positions.py)
    mol_h = build_dimer("", "OC", "H", "MeO")
    if mol_h:
        smi_h = Chem.MolToSmiles(mol_h)
        print(f"\n  H/MeO: atoms={mol_h.GetNumAtoms()}, MW={Descriptors.MolWt(mol_h):.1f}")
        print(f"  SMILES: {smi_h}")
