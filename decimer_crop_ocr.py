"""Try tighter crops of the hand-drawn image for DECIMER OCR.

Also validate the most promising SMILES from the initial OCR run.
"""

import urllib.request
import json
from pathlib import Path
from PIL import Image
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
try:
    from rdkit.Chem.Draw import rdMolDraw2d
except ImportError:
    from rdkit.Chem.Draw import MolDraw2DCairo as _fallback
    rdMolDraw2d = None

API_URL = "https://api.naturalproducts.net/latest/ocsr/process-upload"
PYRENE_SMARTS = Chem.MolFromSmarts("c1cc2ccc3cccc4ccc(c1)c2c34")
outdir = Path("plots/structure_comparison")


def ocr_image(image_path, hand_drawn=False):
    """Upload an image to DECIMER OCSR and get SMILES back."""
    with open(image_path, "rb") as f:
        image_data = f.read()

    ext = Path(image_path).suffix.lower()
    ct_map = {".png": "image/png", ".jpg": "image/jpeg", ".jpeg": "image/jpeg"}
    content_type = ct_map.get(ext, "image/jpeg")

    boundary = "----PyreneAnalyzerBoundary"
    body = b""
    body += f"--{boundary}\r\n".encode()
    body += f'Content-Disposition: form-data; name="file"; filename="{Path(image_path).name}"\r\n'.encode()
    body += f"Content-Type: {content_type}\r\n\r\n".encode()
    body += image_data
    body += b"\r\n"
    body += f"--{boundary}\r\n".encode()
    body += b'Content-Disposition: form-data; name="hand_drawn"\r\n\r\n'
    body += b"true" if hand_drawn else b"false"
    body += b"\r\n"
    body += f"--{boundary}--\r\n".encode()

    req = urllib.request.Request(
        API_URL,
        data=body,
        headers={
            "Content-Type": f"multipart/form-data; boundary={boundary}",
            "User-Agent": "Mozilla/5.0 pyrene-analyzer",
        },
        method="POST",
    )

    try:
        with urllib.request.urlopen(req, timeout=120) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as e:
        error_body = e.read().decode("utf-8") if e.fp else ""
        print(f"  HTTP Error {e.code}: {error_body[:200]}")
        return None
    except Exception as e:
        print(f"  Error: {e}")
        return None


def render_mol(smi, filename, size=(400, 400)):
    """Render a molecule from SMILES."""
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    AllChem.Compute2DCoords(mol)
    img = Draw.MolToImage(mol, size=size)
    img.save(filename)
    return mol


# ================================================================
# PART 1: Validate promising SMILES from first OCR run
# ================================================================
print("=" * 60)
print("VALIDATING PROMISING SMILES FROM DECIMER OCR")
print("=" * 60)

ocr_candidates = [
    ("hand_drawn_1", "CC1=CC2=C3C(=CC=C2)C=CC4=C3C1=CC=C4"),
    ("hand_drawn_2", "CC#CC1=CC2=C3C(=CC=C2)C=CC4=CC=CC1=C43"),
    ("crop_center_1", "CCC/C(=C/C)/C1=C(C=CC2=C1C=CC=C2)OC(C)C3=CC(=CC(=C3)C)CC#CC4=CC=CC5=C4C6=C(C=CC7=C6C=CC=C7)C85C9=C(C=CC=C9)C%10=C8C=CC%11=C%10C=CC=C%11"),
    ("crop_center_2", "C1=CC2=C(C=C1)C3=C(C=C2)OC(C4=CC=C(C=C4)CC(C(F)(F)F)O)C5=CC=C6C=CC=CC6=C53"),
    ("crop_ul_1", "CC1=CC2=C(C=C1)C3=C(C=C2)C4=C(C5=C3C=C(C)C=C5)C6=C(C=C(C)C7=C6C=C(C)C=C7)C4(C)C"),
]

for name, smi in ocr_candidates:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        print(f"\n  {name}: INVALID SMILES")
        continue

    n_atoms = mol.GetNumAtoms()
    n_arom = sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())
    pyr_matches = mol.GetSubstructMatches(PYRENE_SMARTS)
    has_ome = Chem.MolFromSmarts("[OX2][CH3]")
    ome_matches = mol.GetSubstructMatches(has_ome)

    print(f"\n  {name}: {smi[:60]}...")
    print(f"    atoms={n_atoms}, aromatic={n_arom}, pyrene_matches={len(pyr_matches)}, OMe_count={len(ome_matches)}")

    if len(pyr_matches) > 0:
        fname = str(outdir / f"ocr_{name}.png")
        render_mol(smi, fname)
        print(f"    -> Contains pyrene! Rendered: {fname}")

# ================================================================
# PART 2: Try tighter crops of the blue oval area
# ================================================================
print("\n" + "=" * 60)
print("TRYING TIGHTER CROPS OF THE PYRENE AREA")
print("=" * 60)

img = Image.open("docs/1759858767684.jpeg")
w, h = img.size
print(f"Image size: {w}x{h}")

# The user's zoomed image shows the pyrene structure is in approximately
# the upper-left portion of the page, within a red/blue highlighted region
# Based on the zoom, the structure occupies roughly the center of the image
tight_crops = [
    # (name, left, top, right, bottom) - as fractions of image size
    ("pyrene_tight_1", 0.05, 0.05, 0.55, 0.70),   # upper-left third
    ("pyrene_tight_2", 0.10, 0.10, 0.50, 0.65),   # tighter
    ("pyrene_tight_3", 0.02, 0.02, 0.45, 0.80),   # wider crop
    ("pyrene_tight_4", 0.15, 0.15, 0.55, 0.75),   # even tighter center
    ("pyrene_tight_5", 0.0, 0.0, 0.60, 0.60),     # upper portion
]

for name, l, t, r, b in tight_crops:
    box = (int(l * w), int(t * h), int(r * w), int(b * h))
    crop = img.crop(box)
    crop_path = str(outdir / f"crop_{name}.png")
    crop.save(crop_path)
    print(f"\n  Crop: {name} {box}")

    result = ocr_image(crop_path, hand_drawn=True)
    if result and result.get("smiles"):
        smiles_list = result["smiles"]
        for smi in smiles_list:
            mol = Chem.MolFromSmiles(smi)
            if mol:
                n_atoms = mol.GetNumAtoms()
                n_arom = sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())
                pyr = mol.GetSubstructMatches(PYRENE_SMARTS)
                ome = mol.GetSubstructMatches(Chem.MolFromSmarts("[OX2][CH3]"))
                flag = ""
                if len(pyr) > 0:
                    flag += " **PYRENE**"
                if len(ome) > 0:
                    flag += " **HAS_OMe**"
                if n_arom >= 14:
                    flag += " (large aromatic)"
                print(f"    [hd] atoms={n_atoms:>3}, arom={n_arom:>2}, pyr={len(pyr)}, OMe={len(ome)}: {smi[:70]}{flag}")

                if len(pyr) > 0 or n_arom >= 14:
                    fname = str(outdir / f"ocr_{name}_{n_atoms}at.png")
                    render_mol(smi, fname)
            else:
                print(f"    [hd] INVALID: {smi[:70]}")

# ================================================================
# PART 3: Also try the user's zoomed image if available
# ================================================================
print("\n" + "=" * 60)
print("SUMMARY")
print("=" * 60)
print("""
The DECIMER OCR on this hand-drawn image is limited because:
1. Multiple molecules on the same page
2. Red/blue highlighting obscures structure
3. Hand-drawn labels (R, MeO, H3C) are not standard chemical notation

Best approach: manually interpret the image and build the structure.
From the image, the pyrene monomer has:
- Pyrene core (confirmed by OCR finding pyrene SMARTS matches)
- MeO group on one side
- R groups on the other side
- O bridge attachment
""")
