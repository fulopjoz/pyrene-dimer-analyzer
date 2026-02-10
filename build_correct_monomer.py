"""Build the correct pyrene monomer based on the hand-drawn image.

Structure from image (red circle):
- Pyrene = upper naphthalene + lower naphthalene joined by peri-bond
- Upper naphthalene: R (right) + bridge-O (left)
- Lower naphthalene: MeO (left) + R (right)

Pyrene atom positions (RDKit 2D layout):
         10 --- 11          (top vertices)
        /         \
       8           13       (upper inner)
      / \  14--15 / \
     7    |    |    0       (middle)
      \ /  15--14 \ /
       6           1        (lower inner)
        \         /
         4 --- 3            (bottom vertices)

Upper naphthalene H-atoms: 10, 11 (vertices), 8, 13 (inner)
Lower naphthalene H-atoms: 4, 3 (vertices), 6, 1 (inner)
Middle H-atoms: 7, 0

Candidate position sets (all with R-right, MeO/bridge-left):
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from PIL import Image, ImageDraw
from pathlib import Path

outdir = Path("plots/structure_comparison")
outdir.mkdir(parents=True, exist_ok=True)


def build_et_monomer(r1, r2, ome_pos, bridge_pos, label=""):
    """Build Et monomer with specific positions."""
    rw = Chem.RWMol(Chem.MolFromSmiles("c1cc2ccc3cccc4ccc(c1)c2c34"))

    # Et at R1
    c1 = rw.AddAtom(Chem.Atom(6))
    c2 = rw.AddAtom(Chem.Atom(6))
    rw.AddBond(r1, c1, Chem.BondType.SINGLE)
    rw.AddBond(c1, c2, Chem.BondType.SINGLE)

    # Et at R2
    c3 = rw.AddAtom(Chem.Atom(6))
    c4 = rw.AddAtom(Chem.Atom(6))
    rw.AddBond(r2, c3, Chem.BondType.SINGLE)
    rw.AddBond(c3, c4, Chem.BondType.SINGLE)

    # OMe
    o = rw.AddAtom(Chem.Atom(8))
    c = rw.AddAtom(Chem.Atom(6))
    rw.AddBond(ome_pos, o, Chem.BondType.SINGLE)
    rw.AddBond(o, c, Chem.BondType.SINGLE)

    # OH (bridge attachment)
    oh = rw.AddAtom(Chem.Atom(8))
    rw.AddBond(bridge_pos, oh, Chem.BondType.SINGLE)

    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    smi = Chem.MolToSmiles(mol)
    return mol, smi


# Candidate position sets
# (label, R1, R2, OMe_pos, bridge_O_pos)
candidates = [
    # Inner positions (y ≈ ±1.3): R@13,1 on right; OMe@6, bridge@8 on left
    ("A_inner_R13_1_OMe6_Br8", 13, 1, 6, 8),
    # Same but swapped MeO/bridge
    ("B_inner_R13_1_OMe8_Br6", 13, 1, 8, 6),
    # Outer/vertex positions: R@11,3; OMe@4, bridge@10
    ("C_outer_R11_3_OMe4_Br10", 11, 3, 4, 10),
    # Same but swapped
    ("D_outer_R11_3_OMe10_Br4", 11, 3, 10, 4),
    # Mixed: R@13,1; OMe@4, bridge@8
    ("E_mixed_R13_1_OMe4_Br8", 13, 1, 4, 8),
    # Mixed: R@13,3; OMe@6, bridge@8
    ("F_mixed_R13_3_OMe6_Br8", 13, 3, 6, 8),
]

print("=" * 60)
print("CANDIDATE MONOMERS")
print("=" * 60)

mols = []
labels = []

for label, r1, r2, ome, br in candidates:
    mol, smi = build_et_monomer(r1, r2, ome, br, label)
    print(f"\n{label}:")
    print(f"  R at {r1},{r2} (right)")
    print(f"  OMe at {ome} (left)")
    print(f"  bridge-O at {br} (left)")
    print(f"  SMILES: {smi}")

    AllChem.Compute2DCoords(mol)
    mols.append(mol)
    labels.append(label)

    # Render individual
    img = Draw.MolToImage(mol, size=(500, 500))
    img.save(str(outdir / f"monomer_{label}.png"))

# Create comparison grid
print("\n" + "=" * 60)
print("CREATING COMPARISON GRID")
print("=" * 60)

ncols = 3
nrows = 2
cell_w, cell_h = 400, 450
grid = Image.new("RGB", (ncols * cell_w, nrows * cell_h), "white")
draw = ImageDraw.Draw(grid)

for i, (mol, label) in enumerate(zip(mols, labels)):
    row, col = divmod(i, ncols)
    AllChem.Compute2DCoords(mol)
    img = Draw.MolToImage(mol, size=(cell_w, cell_h - 50))
    grid.paste(img, (col * cell_w, row * cell_h))
    # Add label
    short_label = label.split("_", 1)[0] + ": " + label.split("_", 1)[1].replace("_", " ")
    draw.text((col * cell_w + 10, row * cell_h + cell_h - 45),
              short_label, fill="black")

grid.save(str(outdir / "monomer_comparison_grid.png"))
print("Saved monomer_comparison_grid.png")

# Also render the current WRONG monomer for comparison
wrong_mol = Chem.MolFromSmiles("CCc1cc(CC)c2ccc3cc(O)cc4ccc1c2c43")
AllChem.Compute2DCoords(wrong_mol)
img_wrong = Draw.MolToImage(wrong_mol, size=(500, 500))
img_wrong.save(str(outdir / "monomer_CURRENT_no_OMe.png"))
print("Saved monomer_CURRENT_no_OMe.png (for comparison - missing OMe)")

# Print a clear summary
print("\n" + "=" * 60)
print("PLEASE COMPARE WITH YOUR HAND-DRAWN IMAGE")
print("=" * 60)
print("""
Look at monomer_comparison_grid.png and tell me which matches:

A: R at inner-right (13,1), OMe at lower-left (6), bridge at upper-left (8)
B: R at inner-right (13,1), OMe at upper-left (8), bridge at lower-left (6)
C: R at vertex-right (11,3), OMe at bottom-left (4), bridge at top-left (10)
D: R at vertex-right (11,3), OMe at top-left (10), bridge at bottom-left (4)
E: R at 13,1 (inner-right), OMe at bottom-left (4), bridge at upper-left (8)
F: R at 13,3 (right diagonal), OMe at lower-left (6), bridge at upper-left (8)

Or specify the correct positions using the atom numbers from
pyrene_atom_positions.png.
""")
