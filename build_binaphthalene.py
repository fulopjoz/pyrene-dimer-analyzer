"""Build the correct aromatic core: two naphthalene units connected by single bond.

User description:
- Two fused benzenes (= naphthalene, upper part)
- Single bond | connecting them
- Two fused benzenes (= naphthalene, lower part)
- Bond from bottom atom of upper naphthalene to top atom of lower naphthalene

This is a binaphthalene system, NOT pyrene.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from PIL import Image, ImageDraw
from pathlib import Path

outdir = Path("plots/structure_comparison")
outdir.mkdir(parents=True, exist_ok=True)

# First, render naphthalene with atom numbers
naph = Chem.MolFromSmiles("c1ccc2ccccc2c1")
AllChem.Compute2DCoords(naph)

# Get naphthalene positions
conf = naph.GetConformer()
print("Naphthalene atom positions:")
for i in range(naph.GetNumAtoms()):
    pos = conf.GetAtomPosition(i)
    deg = naph.GetAtomWithIdx(i).GetDegree()
    h = "H" if deg == 2 else "junction"
    print(f"  atom {i}: ({pos.x:.2f}, {pos.y:.2f})  {h}")

# Render naphthalene with atom map numbers
for i in range(naph.GetNumAtoms()):
    naph.GetAtomWithIdx(i).SetProp("molAtomMapNumber", str(i))
img = Draw.MolToImage(naph, size=(400, 400))
img.save(str(outdir / "naphthalene_numbered.png"))

# Now try different binaphthalene connections
# "from one benzene bottom atom of upper, to upper atom of lower"

# Build several binaphthalene SMILES connecting different positions
binaph_candidates = [
    # label, SMILES
    ("1,1'-binaph", "c1ccc(-c2cccc3ccccc23)c2ccccc12"),
    ("2,2'-binaph", "c1ccc2cc(-c3ccc4ccccc4c3)ccc2c1"),
    ("1,2'-binaph", "c1ccc(-c2ccc3ccccc3c2)c2ccccc12"),
    ("2,1'-binaph", "c1ccc2cc(-c3cccc4ccccc34)ccc2c1"),
    ("8,8'-binaph_a", "c1cccc(-c2cccc3ccccc23)c1c1ccccc1"),  # may not be valid
]

# Also try building with RWMol for more control
def build_binaph_rwmol(upper_pos, lower_pos):
    """Build binaphthalene by connecting two naphthalenes.

    upper_pos: atom index in upper naphthalene to connect FROM
    lower_pos: atom index in lower naphthalene to connect TO
    """
    naph1 = Chem.MolFromSmiles("c1ccc2ccccc2c1")
    rw = Chem.RWMol(naph1)

    # Add second naphthalene
    naph2 = Chem.MolFromSmiles("c1ccc2ccccc2c1")
    amap = {}
    for a in naph2.GetAtoms():
        new_a = Chem.Atom(a.GetAtomicNum())
        new_a.SetIsAromatic(True)
        amap[a.GetIdx()] = rw.AddAtom(new_a)
    for b in naph2.GetBonds():
        rw.AddBond(amap[b.GetBeginAtomIdx()], amap[b.GetEndAtomIdx()], b.GetBondType())

    # Connect by single bond
    rw.AddBond(upper_pos, amap[lower_pos], Chem.BondType.SINGLE)

    try:
        mol = rw.GetMol()
        Chem.SanitizeMol(mol)
        return mol, amap
    except Exception as e:
        return None, None


print("\n" + "=" * 60)
print("BINAPHTHALENE CANDIDATES (RWMol)")
print("=" * 60)

# Get naphthalene bottom/top atoms from coordinates
naph_clean = Chem.MolFromSmiles("c1ccc2ccccc2c1")
AllChem.Compute2DCoords(naph_clean)
nconf = naph_clean.GetConformer()

h_atoms_naph = []
for i in range(naph_clean.GetNumAtoms()):
    if naph_clean.GetAtomWithIdx(i).GetDegree() == 2:
        pos = nconf.GetAtomPosition(i)
        h_atoms_naph.append((i, pos.x, pos.y))

bottom_atoms = sorted(h_atoms_naph, key=lambda t: t[2])[:3]  # lowest y
top_atoms = sorted(h_atoms_naph, key=lambda t: -t[2])[:3]     # highest y

print(f"Top atoms of naphthalene: {[(a[0], f'y={a[2]:.1f}') for a in top_atoms]}")
print(f"Bottom atoms of naphthalene: {[(a[0], f'y={a[2]:.1f}') for a in bottom_atoms]}")

# Try connecting bottom of upper to top of lower
results = []
for bi, bx, by in bottom_atoms:
    for ti, tx, ty in top_atoms:
        mol, amap = build_binaph_rwmol(bi, ti)
        if mol:
            smi = Chem.MolToSmiles(mol)
            n_arom = sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())
            label = f"upper_{bi}_to_lower_{ti}"
            results.append((label, mol, smi, bi, ti))
            print(f"  {label}: {smi} (atoms={mol.GetNumAtoms()}, aromatic={n_arom})")

# Deduplicate by canonical SMILES
seen = set()
unique = []
for label, mol, smi, bi, ti in results:
    if smi not in seen:
        seen.add(smi)
        unique.append((label, mol, smi, bi, ti))

print(f"\n{len(unique)} unique binaphthalene structures")

# Render each
for i, (label, mol, smi, bi, ti) in enumerate(unique):
    AllChem.Compute2DCoords(mol)
    img = Draw.MolToImage(mol, size=(500, 500))
    img.save(str(outdir / f"binaph_{label}.png"))
    print(f"  Saved binaph_{label}.png")

# Also try the SMILES-based candidates
print("\n" + "=" * 60)
print("BINAPHTHALENE CANDIDATES (SMILES)")
print("=" * 60)

for label, smi in binaph_candidates:
    mol = Chem.MolFromSmiles(smi)
    if mol:
        n_atoms = mol.GetNumAtoms()
        n_arom = sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())
        can_smi = Chem.MolToSmiles(mol)
        print(f"  {label}: atoms={n_atoms}, arom={n_arom}, SMILES={can_smi}")
        AllChem.Compute2DCoords(mol)
        img = Draw.MolToImage(mol, size=(500, 500))
        clean_label = label.replace(",", "_").replace("'", "")
        img.save(str(outdir / f"binaph_{clean_label}.png"))
    else:
        print(f"  {label}: INVALID SMILES")

# Create a comparison grid of all unique binaphthalenes
print("\n" + "=" * 60)
print("CREATING COMPARISON GRID")
print("=" * 60)

# Also add SMILES-based that parsed OK
all_mols = []
all_labels = []
for label, mol, smi, bi, ti in unique:
    all_mols.append(mol)
    all_labels.append(f"RW: {label}")

for label, smi in binaph_candidates:
    mol = Chem.MolFromSmiles(smi)
    if mol:
        can_smi = Chem.MolToSmiles(mol)
        # Check if already in list
        if can_smi not in [Chem.MolToSmiles(m) for m in all_mols]:
            all_mols.append(mol)
            all_labels.append(f"SMI: {label}")

ncols = 3
nrows = (len(all_mols) + ncols - 1) // ncols
cell_w, cell_h = 350, 400

grid = Image.new("RGB", (ncols * cell_w, nrows * cell_h), "white")
draw = ImageDraw.Draw(grid)

for i, (mol, label) in enumerate(zip(all_mols, all_labels)):
    row, col = divmod(i, ncols)
    AllChem.Compute2DCoords(mol)
    img = Draw.MolToImage(mol, size=(cell_w, cell_h - 40))
    grid.paste(img, (col * cell_w, row * cell_h))
    draw.text((col * cell_w + 5, row * cell_h + cell_h - 35), label, fill="black")

grid.save(str(outdir / "binaphthalene_comparison.png"))
print(f"Saved binaphthalene_comparison.png ({len(all_mols)} structures)")
