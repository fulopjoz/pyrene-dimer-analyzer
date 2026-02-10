"""Build 1,1'-binaphthalene monomer candidates with R, MeO, bridge-O.

Confirmed core: 1,1'-binaphthalene
SMILES: c1ccc2c(-c3cccc4ccccc34)cccc2c1

From user description:
- Upper naphthalene: R (right) + bridge-O (left)
- Lower naphthalene: MeO (left) + R (right)
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from PIL import Image, ImageDraw
from pathlib import Path

outdir = Path("plots/structure_comparison")
outdir.mkdir(parents=True, exist_ok=True)

# 1,1'-binaphthalene
BINAPH_SMI = "c1ccc2c(-c3cccc4ccccc34)cccc2c1"

# First render with atom numbers
binaph = Chem.MolFromSmiles(BINAPH_SMI)
AllChem.Compute2DCoords(binaph)

# Get atom positions
conf = binaph.GetConformer()
print("1,1'-Binaphthalene atom positions:")
print(f"  Total atoms: {binaph.GetNumAtoms()}")

# Identify upper vs lower naphthalene
# The connecting bond is between two carbon atoms
# Find which atoms are in which naphthalene
naph_smarts = Chem.MolFromSmarts("c1cccc2ccccc12")
matches = binaph.GetSubstructMatches(naph_smarts)
print(f"\n  Naphthalene matches: {len(matches)}")
for i, m in enumerate(matches):
    print(f"    Match {i}: atoms {m}")

# Get 2D positions for all atoms
positions = {}
for i in range(binaph.GetNumAtoms()):
    pos = conf.GetAtomPosition(i)
    deg = binaph.GetAtomWithIdx(i).GetDegree()
    positions[i] = (pos.x, pos.y, deg)
    h_bearing = "H" if deg == 2 else "junction" if deg == 3 else f"deg{deg}"
    print(f"  atom {i:2d}: ({pos.x:>6.2f}, {pos.y:>6.2f})  deg={deg}  {h_bearing}")

# Render with atom map numbers
binaph_labeled = Chem.MolFromSmiles(BINAPH_SMI)
AllChem.Compute2DCoords(binaph_labeled)
for i in range(binaph_labeled.GetNumAtoms()):
    binaph_labeled.GetAtomWithIdx(i).SetProp("molAtomMapNumber", str(i))
img = Draw.MolToImage(binaph_labeled, size=(600, 600))
img.save(str(outdir / "binaph_11_atom_positions.png"))
print("\nSaved binaph_11_atom_positions.png")

# Classify atoms into upper and lower naphthalene
# Use y-coordinates to determine upper/lower
# First, find the two atoms that form the connecting bond (degree 3 atoms
# that are bonded to each other across the two naphthalene units)
connecting_atoms = []
for bond in binaph.GetBonds():
    a1 = bond.GetBeginAtomIdx()
    a2 = bond.GetEndAtomIdx()
    # Both should be degree 3 (junction atoms)
    if (binaph.GetAtomWithIdx(a1).GetDegree() == 3 and
            binaph.GetAtomWithIdx(a2).GetDegree() == 3):
        # Check if they belong to different naphthalene units
        # by checking if removing this bond would disconnect the molecule
        # Simpler: check if they share a ring
        ri = binaph.GetRingInfo()
        shared_ring = False
        for ring in ri.AtomRings():
            if a1 in ring and a2 in ring:
                shared_ring = True
                break
        if not shared_ring:
            connecting_atoms = [a1, a2]
            print(f"\nConnecting bond: atom {a1} -- atom {a2}")

# Now separate into two naphthalene units using BFS from each connecting atom
from collections import deque


def bfs_component(mol, start, exclude_bond_atoms):
    """BFS to find all atoms in one naphthalene unit."""
    visited = {start}
    queue = deque([start])
    while queue:
        atom_idx = queue.popleft()
        for bond in mol.GetAtomWithIdx(atom_idx).GetBonds():
            other = bond.GetOtherAtomIdx(atom_idx)
            if other not in visited and other not in exclude_bond_atoms:
                visited.add(other)
                queue.append(other)
    return visited


if connecting_atoms:
    a1, a2 = connecting_atoms
    unit1 = bfs_component(binaph, a1, {a2})
    unit2 = bfs_component(binaph, a2, {a1})

    # Determine which is upper (higher y) and lower
    y1_avg = sum(positions[i][1] for i in unit1) / len(unit1)
    y2_avg = sum(positions[i][1] for i in unit2) / len(unit2)

    if y1_avg > y2_avg:
        upper_atoms, lower_atoms = unit1, unit2
        upper_connect, lower_connect = a1, a2
    else:
        upper_atoms, lower_atoms = unit2, unit1
        upper_connect, lower_connect = a2, a1

    print(f"\nUpper naphthalene: atoms {sorted(upper_atoms)} (connect at {upper_connect})")
    print(f"Lower naphthalene: atoms {sorted(lower_atoms)} (connect at {lower_connect})")

    # For each naphthalene unit, classify H-bearing atoms by position
    # (left vs right relative to the connecting bond)
    print("\n--- Upper naphthalene H-bearing atoms ---")
    upper_h = []
    for i in sorted(upper_atoms):
        x, y, deg = positions[i]
        if deg == 2:  # H-bearing
            side = "right" if x > positions[upper_connect][0] else "left"
            upper_h.append((i, x, y, side))
            print(f"  atom {i}: ({x:.2f}, {y:.2f}) [{side}]")

    print("\n--- Lower naphthalene H-bearing atoms ---")
    lower_h = []
    for i in sorted(lower_atoms):
        x, y, deg = positions[i]
        if deg == 2:  # H-bearing
            side = "right" if x > positions[lower_connect][0] else "left"
            lower_h.append((i, x, y, side))
            print(f"  atom {i}: ({x:.2f}, {y:.2f}) [{side}]")


# Build monomer candidates
# From user: upper naph has R(right) + bridge-O(left)
#            lower naph has MeO(left) + R(right)
def build_binaph_monomer(r_pos_upper, bridge_pos, r_pos_lower, ome_pos, label=""):
    """Build 1,1'-binaphthalene monomer with Et, OMe, and bridge-OH."""
    rw = Chem.RWMol(Chem.MolFromSmiles(BINAPH_SMI))

    # Et at upper R position
    c1 = rw.AddAtom(Chem.Atom(6))
    c2 = rw.AddAtom(Chem.Atom(6))
    rw.AddBond(r_pos_upper, c1, Chem.BondType.SINGLE)
    rw.AddBond(c1, c2, Chem.BondType.SINGLE)

    # Et at lower R position
    c3 = rw.AddAtom(Chem.Atom(6))
    c4 = rw.AddAtom(Chem.Atom(6))
    rw.AddBond(r_pos_lower, c3, Chem.BondType.SINGLE)
    rw.AddBond(c3, c4, Chem.BondType.SINGLE)

    # OMe at lower left
    o_me = rw.AddAtom(Chem.Atom(8))
    c_me = rw.AddAtom(Chem.Atom(6))
    rw.AddBond(ome_pos, o_me, Chem.BondType.SINGLE)
    rw.AddBond(o_me, c_me, Chem.BondType.SINGLE)

    # Bridge-OH at upper left
    o_br = rw.AddAtom(Chem.Atom(8))
    rw.AddBond(bridge_pos, o_br, Chem.BondType.SINGLE)

    try:
        mol = rw.GetMol()
        Chem.SanitizeMol(mol)
        return mol, Chem.MolToSmiles(mol)
    except Exception as e:
        print(f"  ERROR building {label}: {e}")
        return None, None


print("\n" + "=" * 60)
print("BUILDING MONOMER CANDIDATES")
print("=" * 60)

candidates = []

# Get right and left atoms for each unit
upper_right = [(i, x, y) for i, x, y, s in upper_h if s == "right"]
upper_left = [(i, x, y) for i, x, y, s in upper_h if s == "left"]
lower_right = [(i, x, y) for i, x, y, s in lower_h if s == "right"]
lower_left = [(i, x, y) for i, x, y, s in lower_h if s == "left"]

print(f"\nUpper-right atoms (R): {[a[0] for a in upper_right]}")
print(f"Upper-left atoms (bridge): {[a[0] for a in upper_left]}")
print(f"Lower-right atoms (R): {[a[0] for a in lower_right]}")
print(f"Lower-left atoms (OMe): {[a[0] for a in lower_left]}")

# Build all combinations
for ur_idx, ur_x, ur_y in upper_right:
    for ul_idx, ul_x, ul_y in upper_left:
        for lr_idx, lr_x, lr_y in lower_right:
            for ll_idx, ll_x, ll_y in lower_left:
                label = f"Rup{ur_idx}_Brup{ul_idx}_Rlo{lr_idx}_OMelo{ll_idx}"
                mol, smi = build_binaph_monomer(ur_idx, ul_idx, lr_idx, ll_idx, label)
                if mol:
                    candidates.append((label, mol, smi, ur_idx, ul_idx, lr_idx, ll_idx))

# Deduplicate
seen = set()
unique = []
for c in candidates:
    if c[2] not in seen:
        seen.add(c[2])
        unique.append(c)

print(f"\n{len(candidates)} candidates, {len(unique)} unique")

# Render each unique candidate
for i, (label, mol, smi, ur, ul, lr, ll) in enumerate(unique):
    AllChem.Compute2DCoords(mol)
    img = Draw.MolToImage(mol, size=(500, 500))
    img.save(str(outdir / f"binaph_mono_{i+1}.png"))
    print(f"  #{i+1} {label}: R@up{ur},lo{lr}  bridge@{ul}  OMe@{ll}")
    print(f"       SMILES: {smi}")

# Create comparison grid
if unique:
    ncols = min(4, len(unique))
    nrows = (len(unique) + ncols - 1) // ncols
    cw, ch = 400, 450

    grid = Image.new("RGB", (ncols * cw, nrows * ch), "white")
    draw = ImageDraw.Draw(grid)

    for i, (label, mol, smi, ur, ul, lr, ll) in enumerate(unique):
        row, col = divmod(i, ncols)
        AllChem.Compute2DCoords(mol)
        img = Draw.MolToImage(mol, size=(cw, ch - 50))
        grid.paste(img, (col * cw, row * ch))
        short = f"#{i+1}: R@up{ur},lo{lr} Br@{ul} OMe@{ll}"
        draw.text((col * cw + 5, row * ch + ch - 45), short, fill="black")

    grid.save(str(outdir / "binaph_monomer_grid.png"))
    print(f"\nSaved binaph_monomer_grid.png ({len(unique)} candidates)")
