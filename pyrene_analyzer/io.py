"""
I/O Module
==========

Functions for loading molecular structures and exporting analysis results.

Supported input formats:
- SDF (Structure Data File) - recommended for conformer ensembles
- MOL2 (Tripos Mol2)
- PDB (Protein Data Bank)

Supported output formats:
- CSV (Comma-Separated Values)
- JSON (JavaScript Object Notation)
- Excel (XLSX)
"""

import json
from pathlib import Path
from typing import List, Optional, Tuple, Union

import pandas as pd
from rdkit import Chem


def load_from_sdf(
    filename: Union[str, Path], remove_hydrogens: bool = False
) -> List[Tuple[Chem.Mol, str]]:
    """
    Load molecules with conformers from an SDF file.

    SDF files can contain multiple molecules, each potentially with multiple
    conformers. This function groups conformers by molecule name.

    Args:
        filename: Path to the SDF file.
        remove_hydrogens: If True, remove hydrogen atoms from molecules.

    Returns:
        List of (molecule, name) tuples. Each molecule may contain multiple
        conformers accessible via mol.GetConformer(conf_id).

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If no valid molecules could be loaded.

    Example:
        >>> molecules = load_from_sdf('conformers.sdf')
        >>> for mol, name in molecules:
        ...     print(f"{name}: {mol.GetNumConformers()} conformers")

    Notes:
        - Molecules with the same _Name property are grouped together
        - Energy values are preserved in conformer properties if present
        - Invalid molecules (parsing errors) are silently skipped
    """
    filepath = Path(filename)

    if not filepath.exists():
        raise FileNotFoundError(f"File not found: {filename}")

    supplier = Chem.SDMolSupplier(str(filepath), removeHs=remove_hydrogens)

    # Group molecules by name
    mol_dict = {}
    mol_count = 0

    for mol in supplier:
        if mol is None:
            continue

        mol_count += 1

        # Get molecule name
        if mol.HasProp("_Name") and mol.GetProp("_Name").strip():
            name = mol.GetProp("_Name").strip()
        else:
            name = f"mol_{len(mol_dict)}"

        # Store energy if available
        energy = None
        for prop_name in [
            "energy",
            "Energy",
            "E",
            "ENERGY",
            "mmff94_energy",
            "rel_energy",
            "dE",
            "potential_energy",
        ]:
            if mol.HasProp(prop_name):
                try:
                    energy = float(mol.GetProp(prop_name))
                    break
                except ValueError:
                    pass

        if name not in mol_dict:
            # First occurrence - store the molecule
            mol_dict[name] = mol
            if energy is not None:
                mol.GetConformer().SetDoubleProp("energy", energy)
        else:
            # Add as new conformer to existing molecule
            conf = mol.GetConformer()
            new_conf_id = mol_dict[name].AddConformer(conf, assignId=True)
            if energy is not None:
                mol_dict[name].GetConformer(new_conf_id).SetDoubleProp("energy", energy)

    molecules = [(mol, name) for name, mol in mol_dict.items()]

    if not molecules:
        raise ValueError(f"No valid molecules found in {filename}")

    return molecules


def load_from_mol2(
    filename: Union[str, Path], remove_hydrogens: bool = False
) -> List[Tuple[Chem.Mol, str]]:
    """
    Load a molecule from a MOL2 file.

    MOL2 files typically contain a single molecule with one conformer.

    Args:
        filename: Path to the MOL2 file.
        remove_hydrogens: If True, remove hydrogen atoms.

    Returns:
        List containing a single (molecule, name) tuple.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the molecule could not be parsed.

    Example:
        >>> molecules = load_from_mol2('pyrene_dimer.mol2')
        >>> mol, name = molecules[0]
    """
    filepath = Path(filename)

    if not filepath.exists():
        raise FileNotFoundError(f"File not found: {filename}")

    mol = Chem.MolFromMol2File(str(filepath), removeHs=remove_hydrogens)

    if mol is None:
        raise ValueError(f"Could not parse MOL2 file: {filename}")

    name = filepath.stem

    return [(mol, name)]


def load_from_pdb(
    filename: Union[str, Path], remove_hydrogens: bool = False
) -> List[Tuple[Chem.Mol, str]]:
    """
    Load a molecule from a PDB file.

    PDB files typically contain a single structure. Multi-model PDB files
    are loaded as a single molecule with multiple conformers.

    Args:
        filename: Path to the PDB file.
        remove_hydrogens: If True, remove hydrogen atoms.

    Returns:
        List containing a single (molecule, name) tuple.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the molecule could not be parsed.

    Example:
        >>> molecules = load_from_pdb('structure.pdb')
        >>> mol, name = molecules[0]
    """
    filepath = Path(filename)

    if not filepath.exists():
        raise FileNotFoundError(f"File not found: {filename}")

    mol = Chem.MolFromPDBFile(str(filepath), removeHs=remove_hydrogens)

    if mol is None:
        raise ValueError(f"Could not parse PDB file: {filename}")

    name = filepath.stem

    return [(mol, name)]


def load_molecules(
    filename: Union[str, Path], remove_hydrogens: bool = False
) -> List[Tuple[Chem.Mol, str]]:
    """
    Load molecules from a file, automatically detecting the format.

    Supported formats: SDF, MOL2, MOL, PDB

    Args:
        filename: Path to the molecular structure file.
        remove_hydrogens: If True, remove hydrogen atoms.

    Returns:
        List of (molecule, name) tuples.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the format is not supported or parsing fails.

    Example:
        >>> molecules = load_molecules('conformers.sdf')
        >>> for mol, name in molecules:
        ...     print(f"Loaded {name}")
    """
    filepath = Path(filename)
    suffix = filepath.suffix.lower()

    if suffix == ".sdf":
        return load_from_sdf(filename, remove_hydrogens)
    elif suffix in [".mol2", ".mol"]:
        return load_from_mol2(filename, remove_hydrogens)
    elif suffix == ".pdb":
        return load_from_pdb(filename, remove_hydrogens)
    else:
        raise ValueError(f"Unsupported file format: {suffix}")


def export_to_csv(
    df: pd.DataFrame, filename: Union[str, Path], float_format: str = "%.3f"
) -> None:
    """
    Export analysis results to a CSV file.

    Args:
        df: DataFrame containing analysis results.
        filename: Output CSV file path.
        float_format: Format string for floating point numbers.

    Example:
        >>> export_to_csv(results, 'analysis.csv')
    """
    filepath = Path(filename)
    filepath.parent.mkdir(parents=True, exist_ok=True)

    df.to_csv(filepath, index=False, float_format=float_format)


def export_to_json(
    df: pd.DataFrame,
    filename: Union[str, Path],
    orient: str = "records",
    indent: int = 2,
) -> None:
    """
    Export analysis results to a JSON file.

    Args:
        df: DataFrame containing analysis results.
        filename: Output JSON file path.
        orient: JSON orientation ('records', 'columns', 'index', etc.).
        indent: Indentation level for pretty printing.

    Example:
        >>> export_to_json(results, 'analysis.json')

    Notes:
        The 'records' orientation produces a list of dictionaries,
        one per conformer, which is most suitable for downstream processing.
    """
    filepath = Path(filename)
    filepath.parent.mkdir(parents=True, exist_ok=True)

    # Convert DataFrame to JSON-serializable format
    data = df.to_dict(orient=orient)

    # Handle NaN values
    def clean_nan(obj):
        if isinstance(obj, dict):
            return {k: clean_nan(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [clean_nan(v) for v in obj]
        elif isinstance(obj, float) and pd.isna(obj):
            return None
        return obj

    data = clean_nan(data)

    with open(filepath, "w") as f:
        json.dump(data, f, indent=indent)


def export_to_excel(
    df: pd.DataFrame,
    filename: Union[str, Path],
    sheet_name: str = "Analysis",
    float_format: str = "%.3f",
) -> None:
    """
    Export analysis results to an Excel file.

    Args:
        df: DataFrame containing analysis results.
        filename: Output Excel file path (.xlsx).
        sheet_name: Name of the worksheet.
        float_format: Format string for floating point numbers.

    Example:
        >>> export_to_excel(results, 'analysis.xlsx')

    Notes:
        Requires openpyxl to be installed.
    """
    filepath = Path(filename)
    filepath.parent.mkdir(parents=True, exist_ok=True)

    # Ensure .xlsx extension
    if filepath.suffix.lower() != ".xlsx":
        filepath = filepath.with_suffix(".xlsx")

    df.to_excel(filepath, sheet_name=sheet_name, index=False, float_format=float_format)


def export_results(
    df: pd.DataFrame, filename: Union[str, Path], formats: Optional[List[str]] = None
) -> List[Path]:
    """
    Export analysis results to multiple formats.

    Args:
        df: DataFrame containing analysis results.
        filename: Base output file path (extension will be added/changed).
        formats: List of formats to export ('csv', 'json', 'excel').
                 If None, exports to all formats.

    Returns:
        List of paths to created files.

    Example:
        >>> files = export_results(results, 'analysis', formats=['csv', 'json'])
        >>> print(f"Created: {files}")
    """
    if formats is None:
        formats = ["csv", "json", "excel"]

    filepath = Path(filename)
    base = filepath.parent / filepath.stem

    created_files = []

    if "csv" in formats:
        csv_path = base.with_suffix(".csv")
        export_to_csv(df, csv_path)
        created_files.append(csv_path)

    if "json" in formats:
        json_path = base.with_suffix(".json")
        export_to_json(df, json_path)
        created_files.append(json_path)

    if "excel" in formats:
        excel_path = base.with_suffix(".xlsx")
        export_to_excel(df, excel_path)
        created_files.append(excel_path)

    return created_files
