"""Render candidate pyrene monomers with MeO at different positions.

Since DECIMER OCR can't cleanly parse the hand-drawn image, we'll
generate all reasonable position combinations and render them so the
user can visually compare with their hand-drawn structure.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from PIL import Image, ImageDraw, ImageFont
from pathlib import Path
import itertools

outdir = Path("plots/structure_comparison")
outdir.mkdir(parents=True, exist_ok=True)

pyrene = Chem.MolFromSmiles("c1cc2ccc3cccc4ccc(c1)c2c34")
AllChem.Compute2DCoords(pyrene)

# Get 2D coordinates
conf = pyrene.GetConformer()
h_atoms = []
for i in range(pyrene.GetNumAtoms()):
    if pyrene.GetAtomWithIdx(i).GetDegree() == 2:
        pos = conf.GetAtomPosition(i)
        h_atoms.append((i, pos.x, pos.y))

# Sort by position
left_atoms = sorted([(i, x, y) for i, x, y in h_atoms if x < -0.5], key=lambda t: -t[2])
right_atoms = sorted([(i, x, y) for i, x, y in h_atoms if x > 0.5], key=lambda t: -t[2])
mid_atoms = sorted([(i, x, y) for i, x, y in h_atoms if -0.5 <= x <= 0.5], key=lambda t: -t[2])

print("LEFT atoms (top to bottom):", [(a[0], f"({a[1]:.1f},{a[2]:.1f})") for a in left_atoms])
print("RIGHT atoms (top to bottom):", [(a[0], f"({a[1]:.1f},{a[2]:.1f})") for a in right_atoms])
print("MID atoms:", [(a[0], f"({a[1]:.1f},{a[2]:.1f})") for a in mid_atoms])


def build_monomer_smiles(r_pos1, r_pos2, ome_pos, o_bridge_pos):
    """Build a pyrene monomer with Et at r_pos, OMe, and OH (bridge)."""
    rw = Chem.RWMol(Chem.MolFromSmiles("c1cc2ccc3cccc4ccc(c1)c2c34"))

    # Et at R1
    c1 = rw.AddAtom(Chem.Atom(6))
    c2 = rw.AddAtom(Chem.Atom(6))
    rw.AddBond(r_pos1, c1, Chem.BondType.SINGLE)
    rw.AddBond(c1, c2, Chem.BondType.SINGLE)

    # Et at R2
    c3 = rw.AddAtom(Chem.Atom(6))
    c4 = rw.AddAtom(Chem.Atom(6))
    rw.AddBond(r_pos2, c3, Chem.BondType.SINGLE)
    rw.AddBond(c3, c4, Chem.BondType.SINGLE)

    # OMe
    o_me = rw.AddAtom(Chem.Atom(8))
    c_me = rw.AddAtom(Chem.Atom(6))
    rw.AddBond(ome_pos, o_me, Chem.BondType.SINGLE)
    rw.AddBond(o_me, c_me, Chem.BondType.SINGLE)

    # OH (bridge attachment)
    o_br = rw.AddAtom(Chem.Atom(8))
    rw.AddBond(o_bridge_pos, o_br, Chem.BondType.SINGLE)

    try:
        mol = rw.GetMol()
        Chem.SanitizeMol(mol)
        return Chem.MolToSmiles(mol), mol
    except:
        return None, None


# Generate candidates
# From the hand-drawn image:
# - R (two groups) on one side
# - MeO on the other side
# - O-bridge on the same side as MeO

# We'll try R on right, MeO+bridge on left, and vice versa
candidates = []

# R-group pairs to try (adjacent atoms on same side)
r_pairs_right = []
for i, (idx1, _, _) in enumerate(right_atoms):
    for idx2, _, _ in right_atoms[i+1:]:
        r_pairs_right.append((idx1, idx2))

r_pairs_left = []
for i, (idx1, _, _) in enumerate(left_atoms):
    for idx2, _, _ in left_atoms[i+1:]:
        r_pairs_left.append((idx1, idx2))

# Scenario 1: R on RIGHT, MeO + bridge on LEFT
for r1, r2 in r_pairs_right:
    for i, (ome_idx, _, _) in enumerate(left_atoms):
        for o_br_idx, _, _ in left_atoms:
            if ome_idx == o_br_idx:
                continue
            positions = {r1, r2, ome_idx, o_br_idx}
            if len(positions) == 4:  # all different
                smi, mol = build_monomer_smiles(r1, r2, ome_idx, o_br_idx)
                if smi:
                    candidates.append({
                        "scenario": "R_right_OMeBridge_left",
                        "R": (r1, r2),
                        "OMe": ome_idx,
                        "O_bridge": o_br_idx,
                        "smiles": smi,
                        "mol": mol,
                    })

# Scenario 2: R on LEFT, MeO + bridge on RIGHT
for r1, r2 in r_pairs_left:
    for i, (ome_idx, _, _) in enumerate(right_atoms):
        for o_br_idx, _, _ in right_atoms:
            if ome_idx == o_br_idx:
                continue
            positions = {r1, r2, ome_idx, o_br_idx}
            if len(positions) == 4:
                smi, mol = build_monomer_smiles(r1, r2, ome_idx, o_br_idx)
                if smi:
                    candidates.append({
                        "scenario": "R_left_OMeBridge_right",
                        "R": (r1, r2),
                        "OMe": ome_idx,
                        "O_bridge": o_br_idx,
                        "smiles": smi,
                        "mol": mol,
                    })

print(f"\nGenerated {len(candidates)} candidate monomers")

# Deduplicate by SMILES
seen = set()
unique_candidates = []
for c in candidates:
    if c["smiles"] not in seen:
        seen.add(c["smiles"])
        unique_candidates.append(c)

print(f"After dedup: {len(unique_candidates)} unique candidates")

# Render each candidate
imgs_per_row = 4
total = len(unique_candidates)
nrows = (total + imgs_per_row - 1) // imgs_per_row
img_w, img_h = 350, 350

# Render in batches as grid images
for batch_idx in range(0, total, 16):
    batch = unique_candidates[batch_idx:batch_idx + 16]
    n_in_batch = len(batch)
    ncols = min(4, n_in_batch)
    nrows_b = (n_in_batch + ncols - 1) // ncols

    grid = Image.new("RGB", (ncols * img_w, nrows_b * (img_h + 30)), "white")
    draw = ImageDraw.Draw(grid)

    for i, c in enumerate(batch):
        row, col = divmod(i, ncols)
        mol = c["mol"]
        AllChem.Compute2DCoords(mol)
        img = Draw.MolToImage(mol, size=(img_w, img_h))
        grid.paste(img, (col * img_w, row * (img_h + 30)))

        label = f"R@{c['R']} OMe@{c['OMe']} Obr@{c['O_bridge']}"
        draw.text((col * img_w + 5, row * (img_h + 30) + img_h + 2), label, fill="black")

    fname = str(outdir / f"monomer_candidates_batch{batch_idx // 16 + 1}.png")
    grid.save(fname)
    print(f"Saved {fname} ({n_in_batch} candidates)")

# Also render just the top candidates with clearer labels
# Most likely based on the image: R on right side, OMe + bridge on left
print("\n" + "=" * 60)
print("TOP CANDIDATES (R on right, OMe+bridge on left)")
print("=" * 60)

right_scenario = [c for c in unique_candidates if c["scenario"] == "R_right_OMeBridge_left"]
for i, c in enumerate(right_scenario[:12]):
    print(f"  #{i+1}: R@{c['R']}, OMe@{c['OMe']}, O-bridge@{c['O_bridge']}")
    print(f"       SMILES: {c['smiles']}")

    mol = c["mol"]
    AllChem.Compute2DCoords(mol)
    img = Draw.MolToImage(mol, size=(500, 500))
    fname = str(outdir / f"candidate_R_right_{i+1}.png")
    img.save(fname)

# Also render pyrene with atom numbers for reference
print("\n" + "=" * 60)
print("PYRENE ATOM POSITIONS REFERENCE")
print("=" * 60)

pyrene_labeled = Chem.MolFromSmiles("c1cc2ccc3cccc4ccc(c1)c2c34")
AllChem.Compute2DCoords(pyrene_labeled)

# Create labeled pyrene image
for i in range(pyrene_labeled.GetNumAtoms()):
    pyrene_labeled.GetAtomWithIdx(i).SetProp("molAtomMapNumber", str(i))

img = Draw.MolToImage(pyrene_labeled, size=(500, 500))
img.save(str(outdir / "pyrene_atom_positions.png"))
print("Saved pyrene_atom_positions.png")

# Print position summary
print("\nPosition layout (from 2D coordinates):")
for label, atoms in [("LEFT", left_atoms), ("RIGHT", right_atoms), ("MID", mid_atoms)]:
    for idx, x, y in atoms:
        vert = "upper" if y > 0 else "lower"
        print(f"  Atom {idx:>2} ({label:>5}, {vert:>5}): x={x:>6.2f}, y={y:>6.2f}")
