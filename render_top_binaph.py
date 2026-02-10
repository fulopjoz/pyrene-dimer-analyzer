"""Render top binaphthalene monomer candidates at large size for comparison.

Focus on the most chemically reasonable position combinations.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from PIL import Image, ImageDraw
from pathlib import Path

outdir = Path("plots/structure_comparison")
outdir.mkdir(parents=True, exist_ok=True)

BINAPH_SMI = "c1ccc2c(-c3cccc4ccccc34)cccc2c1"


def build_et_binaph(r_up, bridge, r_lo, ome):
    """Build Et binaphthalene monomer."""
    rw = Chem.RWMol(Chem.MolFromSmiles(BINAPH_SMI))

    # Et at upper R
    c1 = rw.AddAtom(Chem.Atom(6))
    c2 = rw.AddAtom(Chem.Atom(6))
    rw.AddBond(r_up, c1, Chem.BondType.SINGLE)
    rw.AddBond(c1, c2, Chem.BondType.SINGLE)

    # Et at lower R
    c3 = rw.AddAtom(Chem.Atom(6))
    c4 = rw.AddAtom(Chem.Atom(6))
    rw.AddBond(r_lo, c3, Chem.BondType.SINGLE)
    rw.AddBond(c3, c4, Chem.BondType.SINGLE)

    # OMe
    o = rw.AddAtom(Chem.Atom(8))
    c = rw.AddAtom(Chem.Atom(6))
    rw.AddBond(ome, o, Chem.BondType.SINGLE)
    rw.AddBond(o, c, Chem.BondType.SINGLE)

    # Bridge-OH
    ob = rw.AddAtom(Chem.Atom(8))
    rw.AddBond(bridge, ob, Chem.BondType.SINGLE)

    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    return mol


# Key candidates based on chemistry:
# Atom layout (from 2D):
# Upper naph (right): 0,1,2,3,4,15,16,17,18,19. Connect=4
#   H-atoms: 0(far-right-low), 1(far-right-low), 2(near-right-low),
#            15(left-up), 16(left-up), 17(right-up), 19(far-right-mid)
# Lower naph (left): 5,6,7,8,9,10,11,12,13,14. Connect=5
#   H-atoms: 6(right-low), 7(mid-low), 8(left-low),
#            10(far-left-mid), 11(far-left-up), 12(left-up), 13(near-left-up)

# Most likely: R on the outer ring, bridge/OMe on the inner ring near connection
# Or: R adjacent to connection, OMe on far ring

# Top candidates to show:
top = [
    # (label, r_up, bridge, r_lo, ome, description)
    ("A", 19, 15, 6, 13, "R near conn, bridge/OMe near conn (same ring)"),
    ("B", 17, 15, 6, 8, "R near conn(up), OMe mid-ring(lo)"),
    ("C", 0, 15, 6, 10, "R far(up), OMe far(lo)"),
    ("D", 2, 15, 6, 13, "R adj-conn(up), OMe adj-conn(lo)"),
    ("E", 19, 15, 6, 11, "R far-mid(up), OMe far-mid(lo)"),
    ("F", 17, 16, 6, 13, "R&bridge same ring(up), OMe near(lo)"),
    ("G", 0, 16, 6, 10, "R far(up), bridge near(up), OMe far(lo)"),
    ("H", 1, 15, 6, 12, "R far-low(up), OMe far-up(lo)"),
]

mols = []
labels = []
for label, r_up, bridge, r_lo, ome, desc in top:
    try:
        mol = build_et_binaph(r_up, bridge, r_lo, ome)
        smi = Chem.MolToSmiles(mol)
        AllChem.Compute2DCoords(mol)
        mols.append(mol)
        labels.append(f"{label}: R@{r_up},{r_lo} Br@{bridge} OMe@{ome}")
        print(f"  {label}: {desc}")
        print(f"     SMILES: {smi}")

        # Save individual
        img = Draw.MolToImage(mol, size=(600, 600))
        img.save(str(outdir / f"binaph_top_{label}.png"))
    except Exception as e:
        print(f"  {label}: ERROR {e}")

# Create comparison grid 2x4
ncols = 4
nrows = 2
cw, ch = 450, 500
grid = Image.new("RGB", (ncols * cw, nrows * ch), "white")
draw = ImageDraw.Draw(grid)

for i, (mol, label) in enumerate(zip(mols, labels)):
    row, col = divmod(i, ncols)
    AllChem.Compute2DCoords(mol)
    img = Draw.MolToImage(mol, size=(cw, ch - 50))
    grid.paste(img, (col * cw, row * ch))
    draw.text((col * cw + 5, row * ch + ch - 45), label, fill="black")

grid.save(str(outdir / "binaph_top8_comparison.png"))
print(f"\nSaved binaph_top8_comparison.png")

# Also render the numbered binaphthalene rotated vertically for easier comparison
# Use the CoordGen template to orient vertically
binaph = Chem.MolFromSmiles(BINAPH_SMI)
AllChem.Compute2DCoords(binaph)

# Rotate 90 degrees for vertical layout
conf = binaph.GetConformer()
import math
for i in range(binaph.GetNumAtoms()):
    pos = conf.GetAtomPosition(i)
    # Rotate 90 degrees clockwise: (x,y) -> (y, -x)
    conf.SetAtomPosition(i, (pos.y, -pos.x, 0))

for i in range(binaph.GetNumAtoms()):
    binaph.GetAtomWithIdx(i).SetProp("molAtomMapNumber", str(i))

img = Draw.MolToImage(binaph, size=(600, 600))
img.save(str(outdir / "binaph_11_vertical_numbered.png"))
print("Saved binaph_11_vertical_numbered.png (rotated for vertical view)")
