"""
Aromatic System Registry
========================

Defines AromaticSystem dataclass and a registry of known aromatic
systems with their SMARTS patterns, classification thresholds,
and literature references.

Each aromatic system specifies:
- SMARTS pattern for substructure matching
- Minimum ring atom count for connectivity-based fallback detection
- Classification thresholds for strong/weak excimer vs monomer
- Literature references supporting the threshold values

Usage:
    >>> from pyrene_analyzer.aromatic_systems import get_system, AROMATIC_SYSTEMS
    >>> pyrene = get_system("pyrene")
    >>> print(pyrene.thresholds.strong_distance_range)
    (3.3, 3.7)

    >>> # Register a custom system
    >>> from pyrene_analyzer.aromatic_systems import register_system
    >>> register_system(AromaticSystem(
    ...     name="coronene",
    ...     smarts="...",
    ...     min_ring_atoms=18,
    ...     expected_ring_atoms=24,
    ...     thresholds=ClassificationThresholds(
    ...         strong_angle_max=20.0,
    ...         strong_distance_range=(3.3, 3.7),
    ...         strong_overlap_min=50.0,
    ...         weak_angle_max=60.0,
    ...         weak_distance_max=4.5,
    ...         weak_overlap_min=30.0,
    ...     ),
    ... ))
"""

from dataclasses import dataclass, field
from typing import Dict, Optional, Tuple


@dataclass(frozen=True)
class ClassificationThresholds:
    """
    Literature-backed thresholds for excimer/monomer classification.

    These thresholds define the geometric criteria for classifying
    aromatic dimer conformers as strong excimer, weak excimer, or monomer.

    Attributes:
        strong_angle_max: Maximum plane angle for strong excimer (degrees).
        strong_distance_range: (min, max) interplane distance for strong excimer (A).
        strong_overlap_min: Minimum pi-overlap for strong excimer (%).
        weak_angle_max: Maximum plane angle for weak excimer (degrees).
        weak_distance_max: Maximum interplane distance for weak excimer (A).
        weak_overlap_min: Minimum pi-overlap for weak excimer (%).
        high_angle_warning: Angle above which interplane distance is unreliable (degrees).
    """

    strong_angle_max: float
    strong_distance_range: Tuple[float, float]
    strong_overlap_min: float
    weak_angle_max: float
    weak_distance_max: float
    weak_overlap_min: float
    high_angle_warning: float = 60.0


@dataclass(frozen=True)
class AromaticSystem:
    """
    Definition of an aromatic system with detection pattern and classification criteria.

    Attributes:
        name: Short identifier (e.g., "pyrene", "perylene").
        smarts: SMARTS pattern for substructure matching.
        min_ring_atoms: Minimum atoms in a fused aromatic cluster for
            connectivity-based detection fallback.
        expected_ring_atoms: Typical atom count in the aromatic core.
        thresholds: Classification thresholds for this system.
        references: Literature references supporting the threshold values.
        notes: Additional notes about threshold provenance or limitations.
    """

    name: str
    smarts: str
    min_ring_atoms: int
    expected_ring_atoms: int
    thresholds: ClassificationThresholds
    references: Tuple[str, ...] = field(default_factory=tuple)
    notes: str = ""


# ---------------------------------------------------------------------------
# Global registry
# ---------------------------------------------------------------------------

AROMATIC_SYSTEMS: Dict[str, AromaticSystem] = {}


def register_system(system: AromaticSystem) -> None:
    """
    Register an aromatic system in the global registry.

    Args:
        system: AromaticSystem instance to register.

    Example:
        >>> register_system(AromaticSystem(name="coronene", ...))
    """
    AROMATIC_SYSTEMS[system.name] = system


def get_system(name: str) -> AromaticSystem:
    """
    Look up a registered aromatic system by name.

    Args:
        name: System identifier (case-sensitive).

    Returns:
        The AromaticSystem definition.

    Raises:
        ValueError: If the system name is not registered.
    """
    if name not in AROMATIC_SYSTEMS:
        available = ", ".join(sorted(AROMATIC_SYSTEMS.keys()))
        raise ValueError(
            f"Unknown aromatic system '{name}'. Available: {available}"
        )
    return AROMATIC_SYSTEMS[name]


def list_systems() -> Dict[str, AromaticSystem]:
    """Return a copy of the registry."""
    return dict(AROMATIC_SYSTEMS)


# ---------------------------------------------------------------------------
# Built-in systems
# ---------------------------------------------------------------------------

# Pyrene -- the best-characterized excimer system
register_system(
    AromaticSystem(
        name="pyrene",
        smarts="c1cc2ccc3cccc4ccc(c1)c2c34",
        min_ring_atoms=10,
        expected_ring_atoms=16,
        thresholds=ClassificationThresholds(
            strong_angle_max=20.0,
            strong_distance_range=(3.3, 3.7),
            strong_overlap_min=50.0,
            weak_angle_max=60.0,
            weak_distance_max=4.5,
            weak_overlap_min=30.0,
            high_angle_warning=60.0,
        ),
        references=(
            "Birks & Kazzaz 1968 Proc. R. Soc. A 305, 55 (d=3.34 A)",
            "Ge et al. 2020 J. Mater. Chem. C 8, 10223 (overlap 40-80%)",
            "Basuroy et al. 2021 J. Chem. Phys. 155, 234304 (42% overlap = excimer)",
            "Mazzeo et al. 2024 ChemRxiv (twisted excimer at 50-60 deg)",
            "Marazzi et al. 2024 J. Phys. Chem. Lett. (d=3.24 A TD-DFT)",
        ),
        notes=(
            "Strong overlap threshold set to 50% based on Ge 2020 (40-80% range) "
            "and Basuroy 2021 (42% overlap yields excimer). Previous 70% threshold "
            "was too strict and misclassified real excimers."
        ),
    )
)

# Perylene
register_system(
    AromaticSystem(
        name="perylene",
        smarts="c1ccc2ccc3cccc4ccc(c1)c2c3c45ccccc5",
        min_ring_atoms=14,
        expected_ring_atoms=20,
        thresholds=ClassificationThresholds(
            strong_angle_max=20.0,
            strong_distance_range=(3.4, 3.8),
            strong_overlap_min=50.0,
            weak_angle_max=60.0,
            weak_distance_max=4.5,
            weak_overlap_min=30.0,
            high_angle_warning=60.0,
        ),
        references=(
            "Crystal stacking data: d ~ 3.60 A",
            "Margulies et al. 2017 J. Phys. Chem. A (cofacial perylene dimers)",
        ),
        notes="Distance range widened to 3.4-3.8 A based on perylene crystal data.",
    )
)

# Anthracene
register_system(
    AromaticSystem(
        name="anthracene",
        smarts="c1ccc2cc3ccccc3cc2c1",
        min_ring_atoms=10,
        expected_ring_atoms=14,
        thresholds=ClassificationThresholds(
            strong_angle_max=20.0,
            strong_distance_range=(3.4, 3.93),
            strong_overlap_min=50.0,
            weak_angle_max=60.0,
            weak_distance_max=4.5,
            weak_overlap_min=30.0,
            high_angle_warning=60.0,
        ),
        references=(
            "JACS 2024 Antagonistic effects of distance and overlap (d < 3.93 A)",
            "Schreiber et al. 2011 J. Phys. Chem. A (ab initio aromatic excimers)",
        ),
        notes=(
            "Anthracene excimer formation requires closer approach than pyrene. "
            "Upper distance bound 3.93 A from high-pressure crystal data."
        ),
    )
)

# Naphthalene
register_system(
    AromaticSystem(
        name="naphthalene",
        smarts="c1ccc2ccccc2c1",
        min_ring_atoms=8,
        expected_ring_atoms=10,
        thresholds=ClassificationThresholds(
            strong_angle_max=20.0,
            strong_distance_range=(3.3, 3.7),
            strong_overlap_min=50.0,
            weak_angle_max=60.0,
            weak_distance_max=4.5,
            weak_overlap_min=30.0,
            high_angle_warning=60.0,
        ),
        references=(
            "Schreiber et al. 2011 J. Phys. Chem. A (ab initio aromatic excimers)",
        ),
        notes="Thresholds approximated from pyrene data. Limited naphthalene-specific literature.",
    )
)

# Phenanthrene
register_system(
    AromaticSystem(
        name="phenanthrene",
        smarts="c1ccc2c(c1)cc1ccccc1c2",
        min_ring_atoms=10,
        expected_ring_atoms=14,
        thresholds=ClassificationThresholds(
            strong_angle_max=20.0,
            strong_distance_range=(3.3, 3.7),
            strong_overlap_min=50.0,
            weak_angle_max=60.0,
            weak_distance_max=4.5,
            weak_overlap_min=30.0,
            high_angle_warning=60.0,
        ),
        references=(),
        notes="Thresholds approximated from pyrene data.",
    )
)

# 1,1'-Binaphthalene
# Literature-calibrated thresholds (2026-02-05)
register_system(
    AromaticSystem(
        name="binaphthalene",
        smarts="c1ccc2c(-c3cccc4ccccc34)cccc2c1",
        min_ring_atoms=16,
        expected_ring_atoms=20,
        thresholds=ClassificationThresholds(
            strong_angle_max=20.0,
            strong_distance_range=(3.0, 3.6),  # Updated: DLPNO-CCSD(T)/CBS = 3.43 A
            strong_overlap_min=40.0,
            weak_angle_max=50.0,               # Updated: 60 deg is edge-to-face
            weak_distance_max=4.5,             # Ge et al. (2020): 4.5 A standard
            weak_overlap_min=25.0,
            high_angle_warning=50.0,           # Updated: aligned with weak_angle_max
        ),
        references=(
            "PMC11476719 (2024): DLPNO-CCSD(T)/CBS pyrene dimer d=3.43 A",
            "Ge et al. (2020) J. Mater. Chem. C: pi-overlap 40-80% for excimer",
            "Dai et al. (2024) Molecules 29: twisted excimer 5 kcal/mol more stable",
            "Carter-Fenk & Herbert (2020): steric control of slip-stacking",
            "Ge et al. (2020) J. Mater. Chem. C 8, 10223: 4.5 A weak excimer distance",
        ),
        notes=(
            "Thresholds calibrated from DLPNO-CCSD(T)/CBS pyrene dimer benchmark. "
            "Strong distance 3.0-3.6 A centers on 3.43 A equilibrium. Weak angle "
            "reduced to 50 deg (60 deg represents edge-to-face, not pi-stacking). "
            "Weak distance 4.5 A matches all other systems and Ge et al. (2020) "
            "standard. Binaphthalene has intrinsic dihedral twist (~68 deg in "
            "ground state) due to axial chirality. Overlap thresholds lowered vs "
            "pyrene because twisted geometry achieves less cofacial overlap."
        ),
    )
)
