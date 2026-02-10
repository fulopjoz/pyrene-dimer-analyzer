"""Render pyrene monomer candidates with MeO at various positions.

The hand-drawn image shows:
- H3C-O on upper-left (bridge attachment or MeO)
- MeO on lower-left
- R on upper-right and lower-right

We need to find which atom positions match this layout.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2d
from pathlib import Path
import os

outdir = Path("plots/structure_comparison")
outdir.mkdir(parents=True, exist_ok=True)

# First: render pyrene with atom indices to establish coordinate mapping
pyrene = Chem.MolFromSmiles("c1cc2ccc3cccc4ccc(c1)c2c34")
AllChem.Compute2DCoords(pyrene)

# Get 2D coordinates for each atom
conf = pyrene.GetConformer()
print("Pyrene atom positions (2D coordinates):")
print(f"{'Idx':>3} {'X':>8} {'Y':>8}  Ring membership")
for i in range(pyrene.GetNumAtoms()):
    pos = conf.GetAtomPosition(i)
    # Determine which rings this atom belongs to
    rings = []
    ring_info = pyrene.GetRingInfo()
    for r_idx, ring in enumerate(ring_info.AtomRings()):
        if i in ring:
            rings.append(chr(65 + r_idx))  # A, B, C, D
    ring_str = ",".join(rings)
    degree = pyrene.GetAtomWithIdx(i).GetDegree()
    has_h = "H" if degree == 2 else "junction"
    print(f"  {i:>2}  {pos.x:>8.3f} {pos.y:>8.3f}  Ring {ring_str:>5}  {has_h}")

# Sort H-bearing atoms by position
h_atoms = []
for i in range(pyrene.GetNumAtoms()):
    if pyrene.GetAtomWithIdx(i).GetDegree() == 2:
        pos = conf.GetAtomPosition(i)
        h_atoms.append((i, pos.x, pos.y))

print("\nH-bearing atoms sorted by X (left to right):")
for idx, x, y in sorted(h_atoms, key=lambda t: t[1]):
    side = "LEFT" if x < 0 else "RIGHT"
    vert = "upper" if y > 0 else "lower"
    print(f"  atom {idx:>2}: x={x:>7.3f}, y={y:>7.3f}  ({vert}-{side})")

print("\nH-bearing atoms sorted by Y (bottom to top):")
for idx, x, y in sorted(h_atoms, key=lambda t: t[2]):
    side = "LEFT" if x < 0 else "RIGHT"
    vert = "upper" if y > 0 else "lower"
    print(f"  atom {idx:>2}: x={x:>7.3f}, y={y:>7.3f}  ({vert}-{side})")

# Render pyrene with clear position labels
def render_mol_with_labels(mol, filename, title="", size=(600, 500)):
    drawer = rdMolDraw2d.MolDraw2DCairo(size[0], size[1])
    opts = drawer.drawOptions()
    opts.addAtomIndices = True
    AllChem.Compute2DCoords(mol)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    with open(filename, "wb") as f:
        f.write(drawer.GetDrawingText())

render_mol_with_labels(pyrene, str(outdir / "pyrene_positions.png"))
print(f"\nSaved pyrene_positions.png")

# Now generate candidate monomers
# Based on the image: R on RIGHT side, MeO+bridge-O on LEFT side
# We'll try multiple position combinations

candidates = []

# Identify left-side and right-side atoms
left_atoms = [(i, x, y) for i, x, y in h_atoms if x < 0]
right_atoms = [(i, x, y) for i, x, y in h_atoms if x >= 0]

print(f"\nLeft-side H atoms: {[a[0] for a in sorted(left_atoms, key=lambda t: -t[2])]}")
print(f"Right-side H atoms: {[a[0] for a in sorted(right_atoms, key=lambda t: -t[2])]}")

# Strategy: try several position combinations
# For each: (bridge_O_pos, OMe_pos, R_pos1, R_pos2)
# Image: bridge-O upper-left, MeO lower-left, R upper-right, R lower-right

left_sorted = sorted(left_atoms, key=lambda t: -t[2])  # top to bottom
right_sorted = sorted(right_atoms, key=lambda t: -t[2])  # top to bottom

print(f"\nLeft atoms top-to-bottom: {[(a[0], f'y={a[2]:.2f}') for a in left_sorted]}")
print(f"Right atoms top-to-bottom: {[(a[0], f'y={a[2]:.2f}') for a in right_sorted]}")

# Try different position assignments
position_sets = []

# For bridge-O: try top-left atoms
# For MeO: try bottom-left atoms
# For R: try right-side atoms (top-right and bottom-right)

# Get the top 2 left and top 2 right
if len(left_sorted) >= 2 and len(right_sorted) >= 2:
    # Option 1: bridge at topmost-left, OMe at second-left, R at top-right and second-right
    position_sets.append({
        "label": "v1",
        "bridge_O": left_sorted[0][0],
        "OMe": left_sorted[1][0],
        "R1": right_sorted[0][0],
        "R2": right_sorted[1][0],
    })
    # Option 2: bridge at topmost-left, OMe at third-left, R at top-right and bottom-right
    if len(left_sorted) >= 3:
        position_sets.append({
            "label": "v2",
            "bridge_O": left_sorted[0][0],
            "OMe": left_sorted[2][0],
            "R1": right_sorted[0][0],
            "R2": right_sorted[1][0],
        })
    # Option 3: bridge at second-left, OMe at third-left, R at right pair
    if len(left_sorted) >= 3:
        position_sets.append({
            "label": "v3",
            "bridge_O": left_sorted[1][0],
            "OMe": left_sorted[2][0],
            "R1": right_sorted[0][0],
            "R2": right_sorted[1][0],
        })
    # Option 4: try different right pairs
    if len(right_sorted) >= 3:
        position_sets.append({
            "label": "v4",
            "bridge_O": left_sorted[0][0],
            "OMe": left_sorted[1][0],
            "R1": right_sorted[0][0],
            "R2": right_sorted[2][0],
        })
    # Option 5: use widest-apart right positions
    if len(right_sorted) >= 4:
        position_sets.append({
            "label": "v5",
            "bridge_O": left_sorted[0][0],
            "OMe": left_sorted[1][0],
            "R1": right_sorted[0][0],
            "R2": right_sorted[-1][0],
        })

# Build and render each candidate monomer (with Et as example R-group)
for ps in position_sets:
    label = ps["label"]
    bridge_pos = ps["bridge_O"]
    ome_pos = ps["OMe"]
    r1_pos = ps["R1"]
    r2_pos = ps["R2"]

    print(f"\n--- Candidate {label} ---")
    print(f"  Bridge-O at atom {bridge_pos}")
    print(f"  OMe at atom {ome_pos}")
    print(f"  R (Et) at atoms {r1_pos}, {r2_pos}")

    # Build using RWMol
    rw = Chem.RWMol(Chem.MolFromSmiles("c1cc2ccc3cccc4ccc(c1)c2c34"))

    # Add Et at R1
    c1 = rw.AddAtom(Chem.Atom(6))
    c2 = rw.AddAtom(Chem.Atom(6))
    rw.AddBond(r1_pos, c1, Chem.BondType.SINGLE)
    rw.AddBond(c1, c2, Chem.BondType.SINGLE)

    # Add Et at R2
    c3 = rw.AddAtom(Chem.Atom(6))
    c4 = rw.AddAtom(Chem.Atom(6))
    rw.AddBond(r2_pos, c3, Chem.BondType.SINGLE)
    rw.AddBond(c3, c4, Chem.BondType.SINGLE)

    # Add OMe
    o_me = rw.AddAtom(Chem.Atom(8))
    c_me = rw.AddAtom(Chem.Atom(6))
    rw.AddBond(ome_pos, o_me, Chem.BondType.SINGLE)
    rw.AddBond(o_me, c_me, Chem.BondType.SINGLE)

    # Add OH (bridge attachment point)
    o_bridge = rw.AddAtom(Chem.Atom(8))
    rw.AddBond(bridge_pos, o_bridge, Chem.BondType.SINGLE)

    try:
        mol = rw.GetMol()
        Chem.SanitizeMol(mol)
        smi = Chem.MolToSmiles(mol)
        print(f"  SMILES: {smi}")

        # Render
        AllChem.Compute2DCoords(mol)
        fname = str(outdir / f"monomer_candidate_{label}.png")
        render_mol_with_labels(mol, fname, size=(600, 500))
        print(f"  Saved: {fname}")

    except Exception as e:
        print(f"  ERROR: {e}")

# Also render the CURRENT (wrong) monomer for comparison
print("\n--- Current (WRONG) monomer ---")
print("  R at 6,8; O-bridge at 0; NO OMe")
current = Chem.MolFromSmiles("CCc1cc(CC)c2ccc3cc(O)cc4ccc1c2c43")
AllChem.Compute2DCoords(current)
render_mol_with_labels(current, str(outdir / "monomer_CURRENT_WRONG.png"), size=(600, 500))
print(f"  SMILES: {Chem.MolToSmiles(current)}")
print(f"  Saved monomer_CURRENT_WRONG.png")
