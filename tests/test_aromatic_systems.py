"""
Tests for the aromatic system registry.
"""

import pytest

from pyrene_analyzer.aromatic_systems import (
    AROMATIC_SYSTEMS,
    AromaticSystem,
    ClassificationThresholds,
    get_system,
    list_systems,
    register_system,
)


class TestClassificationThresholds:
    """Tests for ClassificationThresholds dataclass."""

    def test_create_thresholds(self):
        t = ClassificationThresholds(
            strong_angle_max=20.0,
            strong_distance_range=(3.3, 3.7),
            strong_overlap_min=50.0,
            weak_angle_max=60.0,
            weak_distance_max=4.5,
            weak_overlap_min=30.0,
        )
        assert t.strong_angle_max == 20.0
        assert t.strong_distance_range == (3.3, 3.7)
        assert t.strong_overlap_min == 50.0
        assert t.high_angle_warning == 60.0  # default

    def test_custom_high_angle_warning(self):
        t = ClassificationThresholds(
            strong_angle_max=20.0,
            strong_distance_range=(3.3, 3.7),
            strong_overlap_min=50.0,
            weak_angle_max=60.0,
            weak_distance_max=4.5,
            weak_overlap_min=30.0,
            high_angle_warning=45.0,
        )
        assert t.high_angle_warning == 45.0

    def test_thresholds_are_frozen(self):
        t = ClassificationThresholds(
            strong_angle_max=20.0,
            strong_distance_range=(3.3, 3.7),
            strong_overlap_min=50.0,
            weak_angle_max=60.0,
            weak_distance_max=4.5,
            weak_overlap_min=30.0,
        )
        with pytest.raises(AttributeError):
            t.strong_angle_max = 30.0


class TestAromaticSystem:
    """Tests for AromaticSystem dataclass."""

    def test_create_system(self):
        t = ClassificationThresholds(
            strong_angle_max=20.0,
            strong_distance_range=(3.3, 3.7),
            strong_overlap_min=50.0,
            weak_angle_max=60.0,
            weak_distance_max=4.5,
            weak_overlap_min=30.0,
        )
        system = AromaticSystem(
            name="test",
            smarts="c1ccccc1",
            min_ring_atoms=6,
            expected_ring_atoms=6,
            thresholds=t,
        )
        assert system.name == "test"
        assert system.references == ()
        assert system.notes == ""

    def test_system_is_frozen(self):
        t = ClassificationThresholds(
            strong_angle_max=20.0,
            strong_distance_range=(3.3, 3.7),
            strong_overlap_min=50.0,
            weak_angle_max=60.0,
            weak_distance_max=4.5,
            weak_overlap_min=30.0,
        )
        system = AromaticSystem(
            name="test",
            smarts="c1ccccc1",
            min_ring_atoms=6,
            expected_ring_atoms=6,
            thresholds=t,
        )
        with pytest.raises(AttributeError):
            system.name = "modified"


class TestRegistry:
    """Tests for the aromatic system registry."""

    def test_registry_contains_pyrene(self):
        assert "pyrene" in AROMATIC_SYSTEMS

    def test_registry_contains_all_five_systems(self):
        expected = {"pyrene", "perylene", "anthracene", "naphthalene", "phenanthrene"}
        assert expected.issubset(set(AROMATIC_SYSTEMS.keys()))

    def test_get_system_valid(self):
        system = get_system("pyrene")
        assert system.name == "pyrene"
        assert system.expected_ring_atoms == 16

    def test_get_system_invalid_raises(self):
        with pytest.raises(ValueError, match="Unknown aromatic system"):
            get_system("nonexistent_system")

    def test_get_system_error_lists_available(self):
        with pytest.raises(ValueError, match="pyrene"):
            get_system("nonexistent_system")

    def test_register_custom_system(self):
        t = ClassificationThresholds(
            strong_angle_max=20.0,
            strong_distance_range=(3.3, 3.7),
            strong_overlap_min=50.0,
            weak_angle_max=60.0,
            weak_distance_max=4.5,
            weak_overlap_min=30.0,
        )
        system = AromaticSystem(
            name="_test_custom",
            smarts="c1ccccc1",
            min_ring_atoms=6,
            expected_ring_atoms=6,
            thresholds=t,
        )
        register_system(system)
        assert "_test_custom" in AROMATIC_SYSTEMS
        retrieved = get_system("_test_custom")
        assert retrieved.name == "_test_custom"
        # Clean up
        del AROMATIC_SYSTEMS["_test_custom"]

    def test_list_systems(self):
        systems = list_systems()
        assert isinstance(systems, dict)
        assert "pyrene" in systems
        # Verify it's a copy
        systems["fake"] = None
        assert "fake" not in AROMATIC_SYSTEMS


class TestPyreneThresholds:
    """Verify pyrene thresholds match literature values."""

    def test_pyrene_strong_angle(self):
        t = get_system("pyrene").thresholds
        assert t.strong_angle_max == 20.0

    def test_pyrene_strong_distance(self):
        t = get_system("pyrene").thresholds
        assert t.strong_distance_range == (3.3, 3.7)

    def test_pyrene_strong_overlap_is_50(self):
        """Overlap threshold lowered from 70% to 50% per Ge 2020 / Basuroy 2021."""
        t = get_system("pyrene").thresholds
        assert t.strong_overlap_min == 50.0

    def test_pyrene_weak_thresholds(self):
        t = get_system("pyrene").thresholds
        assert t.weak_angle_max == 60.0
        assert t.weak_distance_max == 4.5
        assert t.weak_overlap_min == 30.0

    def test_pyrene_high_angle_warning(self):
        t = get_system("pyrene").thresholds
        assert t.high_angle_warning == 60.0

    def test_pyrene_smarts(self):
        system = get_system("pyrene")
        assert system.smarts == "c1cc2ccc3cccc4ccc(c1)c2c34"

    def test_pyrene_min_ring_atoms(self):
        system = get_system("pyrene")
        assert system.min_ring_atoms == 10

    def test_pyrene_has_references(self):
        system = get_system("pyrene")
        assert len(system.references) >= 3


class TestOtherSystemThresholds:
    """Verify thresholds for non-pyrene systems."""

    def test_perylene_distance_range(self):
        t = get_system("perylene").thresholds
        assert t.strong_distance_range == (3.4, 3.8)

    def test_anthracene_distance_upper_bound(self):
        t = get_system("anthracene").thresholds
        d_min, d_max = t.strong_distance_range
        assert d_max == pytest.approx(3.93)

    def test_naphthalene_min_ring_atoms(self):
        system = get_system("naphthalene")
        assert system.min_ring_atoms == 8
        assert system.expected_ring_atoms == 10

    def test_phenanthrene_expected_atoms(self):
        system = get_system("phenanthrene")
        assert system.expected_ring_atoms == 14


class TestBinaphthaleneThresholds:
    """Verify binaphthalene thresholds match literature values."""

    def test_binaphthalene_registered(self):
        assert "binaphthalene" in AROMATIC_SYSTEMS

    def test_binaphthalene_smarts(self):
        system = get_system("binaphthalene")
        assert system.smarts == "c1ccc2c(-c3cccc4ccccc34)cccc2c1"

    def test_binaphthalene_ring_atoms(self):
        system = get_system("binaphthalene")
        assert system.min_ring_atoms == 16
        assert system.expected_ring_atoms == 20

    def test_binaphthalene_strong_thresholds(self):
        t = get_system("binaphthalene").thresholds
        assert t.strong_angle_max == 20.0
        assert t.strong_distance_range == (3.0, 3.6)
        assert t.strong_overlap_min == 40.0

    def test_binaphthalene_weak_distance_is_4_5(self):
        """Weak distance 4.5 A per Ge et al. (2020) J. Mater. Chem. C."""
        t = get_system("binaphthalene").thresholds
        assert t.weak_distance_max == 4.5

    def test_binaphthalene_weak_angle(self):
        """Weak angle 50 deg (60 deg is edge-to-face)."""
        t = get_system("binaphthalene").thresholds
        assert t.weak_angle_max == 50.0

    def test_binaphthalene_has_references(self):
        system = get_system("binaphthalene")
        assert len(system.references) >= 4
