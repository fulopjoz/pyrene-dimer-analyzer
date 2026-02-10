"""Use DECIMER API to OCR chemical structures from the hand-drawn image."""

import urllib.request
import json
import os
from pathlib import Path

API_URL = "https://api.naturalproducts.net/latest/ocsr/process-upload"


def ocr_image(image_path, hand_drawn=False):
    """Upload an image to DECIMER OCSR and get SMILES back."""
    # Read image file
    with open(image_path, "rb") as f:
        image_data = f.read()

    # Determine content type
    ext = Path(image_path).suffix.lower()
    ct_map = {".png": "image/png", ".jpg": "image/jpeg", ".jpeg": "image/jpeg"}
    content_type = ct_map.get(ext, "image/jpeg")

    # Build multipart form data
    boundary = "----PyreneAnalyzerBoundary"
    body = b""

    # File field
    body += f"--{boundary}\r\n".encode()
    body += f'Content-Disposition: form-data; name="file"; filename="{Path(image_path).name}"\r\n'.encode()
    body += f"Content-Type: {content_type}\r\n\r\n".encode()
    body += image_data
    body += b"\r\n"

    # hand_drawn field
    body += f"--{boundary}\r\n".encode()
    body += b'Content-Disposition: form-data; name="hand_drawn"\r\n\r\n'
    body += b"true" if hand_drawn else b"false"
    body += b"\r\n"

    body += f"--{boundary}--\r\n".encode()

    # Make request
    req = urllib.request.Request(
        API_URL,
        data=body,
        headers={
            "Content-Type": f"multipart/form-data; boundary={boundary}",
            "User-Agent": "Mozilla/5.0 pyrene-analyzer",
        },
        method="POST",
    )

    print(f"Uploading {image_path} to DECIMER API (hand_drawn={hand_drawn})...")
    try:
        with urllib.request.urlopen(req, timeout=120) as resp:
            result = json.loads(resp.read().decode("utf-8"))
            return result
    except urllib.error.HTTPError as e:
        error_body = e.read().decode("utf-8") if e.fp else ""
        print(f"HTTP Error {e.code}: {e.reason}")
        print(f"Response: {error_body[:500]}")
        return None
    except Exception as e:
        print(f"Error: {e}")
        return None


if __name__ == "__main__":
    # The hand-drawn image showing the pyrene monomer
    image_path = "docs/1759858767684.jpeg"

    # Try with hand_drawn=True first (better for sketches)
    print("=" * 60)
    print("DECIMER OCR - Hand-drawn mode")
    print("=" * 60)
    result = ocr_image(image_path, hand_drawn=True)
    if result:
        print(f"\nResult: {json.dumps(result, indent=2)}")

    # Also try with hand_drawn=False
    print("\n" + "=" * 60)
    print("DECIMER OCR - Standard mode")
    print("=" * 60)
    result2 = ocr_image(image_path, hand_drawn=False)
    if result2:
        print(f"\nResult: {json.dumps(result2, indent=2)}")

    # If the full image has too much noise, try with cropped versions
    # We can crop the image to just the chemical structure
    try:
        from PIL import Image
        import io

        img = Image.open(image_path)
        w, h = img.size
        print(f"\nOriginal image size: {w}x{h}")

        # The chemical structure is roughly in the center-left area
        # Based on the blue oval in the image
        # Try a few crop regions
        crops = [
            ("full", (0, 0, w, h)),
            ("left_half", (0, 0, w // 2, h)),
            ("center", (w // 6, h // 6, 5 * w // 6, 5 * h // 6)),
            ("upper_left", (0, 0, int(w * 0.6), int(h * 0.7))),
        ]

        for name, box in crops:
            crop = img.crop(box)
            crop_path = f"plots/structure_comparison/crop_{name}.png"
            crop.save(crop_path)
            print(f"\nCrop: {name} ({box})")

            for hd in [True, False]:
                r = ocr_image(crop_path, hand_drawn=hd)
                if r:
                    smiles_list = r.get("smiles", r.get("result", []))
                    if isinstance(smiles_list, list):
                        for s in smiles_list:
                            print(f"  [hd={hd}] {s}")
                    else:
                        print(f"  [hd={hd}] {smiles_list}")

    except ImportError:
        print("\nPillow not installed, skipping crop attempts")
