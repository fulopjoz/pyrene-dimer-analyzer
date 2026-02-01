"""
Tests for the geometry module.

Tests cover:
- SVD plane fitting
- Plane angle calculations
- Inter-plane distance calculations
- π-overlap calculations
- Centroid distance calculations
- Slip-stack displacement calculations
"""

import numpy as np
import pytest

from pyrene_analyzer.geometry import (
    _create_projection_basis,
    _project_to_2d,
    calculate_centroid_distance,
    calculate_interplane_distance,
    calculate_pi_overlap,
    calculate_plane_angle,
    calculate_slip_stack_displacement,
    calculate_tilt_angle,
    fit_plane_svd,
)


class TestFitPlaneSVD:
    """Tests for fit_plane_svd function."""

    def test_xy_plane(self):
        """Points in XY plane should have normal along Z axis."""
        coords = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [1.0, 1.0, 0.0],
            ]
        )
        centroid, normal = fit_plane_svd(coords)

        # Centroid should be at (0.5, 0.5, 0)
        np.testing.assert_array_almost_equal(centroid, [0.5, 0.5, 0.0])

        # Normal should be along Z axis (either direction)
        assert abs(abs(normal[2]) - 1.0) < 0.01

    def test_xz_plane(self):
        """Points in XZ plane should have normal along Y axis."""
        coords = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0],
                [1.0, 0.0, 1.0],
            ]
        )
        centroid, normal = fit_plane_svd(coords)

        # Normal should be along Y axis
        assert abs(abs(normal[1]) - 1.0) < 0.01

    def test_tilted_plane(self):
        """Test plane tilted at 45 degrees."""
        coords = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 1.0],
                [1.0, 1.0, 1.0],
            ]
        )
        centroid, normal = fit_plane_svd(coords)

        # Normal should be unit vector
        assert abs(np.linalg.norm(normal) - 1.0) < 0.001

    def test_minimum_points(self):
        """Test with minimum 3 points."""
        coords = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
            ]
        )
        centroid, normal = fit_plane_svd(coords)

        # Should work with 3 points
        assert centroid is not None
        assert normal is not None

    def test_too_few_points(self):
        """Should raise error with fewer than 3 points."""
        coords = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
            ]
        )
        with pytest.raises(ValueError):
            fit_plane_svd(coords)

    def test_wrong_dimensions(self):
        """Should raise error with wrong dimensions."""
        coords = np.array(
            [
                [0.0, 0.0],
                [1.0, 0.0],
                [0.0, 1.0],
            ]
        )
        with pytest.raises(ValueError):
            fit_plane_svd(coords)

    def test_normal_is_unit_vector(self):
        """Normal vector should always be unit length."""
        np.random.seed(42)
        for _ in range(10):
            coords = np.random.randn(10, 3)
            _, normal = fit_plane_svd(coords)
            assert abs(np.linalg.norm(normal) - 1.0) < 0.001


class TestCalculatePlaneAngle:
    """Tests for calculate_plane_angle function."""

    def test_parallel_planes(self, sample_coords_parallel):
        """Parallel planes should have angle ~0°."""
        coords1, coords2 = sample_coords_parallel
        angle = calculate_plane_angle(coords1, coords2)
        assert angle < 5.0  # Allow small numerical error

    def test_perpendicular_planes(self, sample_coords_perpendicular):
        """Perpendicular planes should have angle ~90°."""
        coords1, coords2 = sample_coords_perpendicular
        angle = calculate_plane_angle(coords1, coords2)
        assert 85.0 < angle <= 90.0

    def test_45_degree_planes(self, sample_coords_45deg):
        """Planes at 45° should give angle ~45°."""
        coords1, coords2 = sample_coords_45deg
        angle = calculate_plane_angle(coords1, coords2)
        assert 40.0 < angle < 50.0

    def test_angle_range(self):
        """Angle should always be in [0, 90] range."""
        np.random.seed(42)
        for _ in range(20):
            coords1 = np.random.randn(5, 3)
            coords2 = np.random.randn(5, 3)
            angle = calculate_plane_angle(coords1, coords2)
            assert 0 <= angle <= 90

    def test_same_plane(self):
        """Same plane should have angle 0°."""
        coords = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [1.0, 1.0, 0.0],
            ]
        )
        angle = calculate_plane_angle(coords, coords)
        assert angle < 1.0


class TestCalculateInterplaneDistance:
    """Tests for calculate_interplane_distance function."""

    def test_parallel_planes_distance(self, sample_coords_parallel):
        """Parallel planes at z=0 and z=3.5 should have distance 3.5."""
        coords1, coords2 = sample_coords_parallel
        distance = calculate_interplane_distance(coords1, coords2)
        assert abs(distance - 3.5) < 0.1

    def test_same_plane_distance(self):
        """Same plane should have distance 0."""
        coords = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [1.0, 1.0, 0.0],
            ]
        )
        distance = calculate_interplane_distance(coords, coords)
        assert distance < 0.1

    def test_distance_positive(self):
        """Distance should always be positive."""
        np.random.seed(42)
        for _ in range(20):
            coords1 = np.random.randn(5, 3)
            coords2 = np.random.randn(5, 3) + [0, 0, 5]
            distance = calculate_interplane_distance(coords1, coords2)
            assert distance >= 0


class TestCalculateCentroidDistance:
    """Tests for calculate_centroid_distance function."""

    def test_same_centroid(self):
        """Same coordinates should have distance 0."""
        coords = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
            ]
        )
        distance = calculate_centroid_distance(coords, coords)
        assert distance < 0.001

    def test_known_distance(self):
        """Test with known centroid distance."""
        coords1 = np.array(
            [
                [0.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [1.0, 2.0, 0.0],
            ]
        )  # Centroid at (1, 0.67, 0)

        coords2 = np.array(
            [
                [0.0, 0.0, 3.0],
                [2.0, 0.0, 3.0],
                [1.0, 2.0, 3.0],
            ]
        )  # Centroid at (1, 0.67, 3)

        distance = calculate_centroid_distance(coords1, coords2)
        assert abs(distance - 3.0) < 0.1


class TestCalculateSlipStackDisplacement:
    """Tests for calculate_slip_stack_displacement function."""

    def test_perfectly_stacked(self, sample_coords_parallel):
        """Perfectly stacked planes should have slip ~0."""
        coords1, coords2 = sample_coords_parallel
        slip = calculate_slip_stack_displacement(coords1, coords2)
        assert slip < 0.1

    def test_offset_planes(self, sample_coords_offset):
        """Offset planes should have non-zero slip."""
        coords1, coords2 = sample_coords_offset
        slip = calculate_slip_stack_displacement(coords1, coords2)
        assert slip > 1.0  # Significant lateral offset

    def test_slip_positive(self):
        """Slip should always be non-negative."""
        np.random.seed(42)
        for _ in range(20):
            coords1 = np.random.randn(5, 3)
            coords2 = np.random.randn(5, 3)
            slip = calculate_slip_stack_displacement(coords1, coords2)
            assert slip >= 0


class TestCalculatePiOverlap:
    """Tests for calculate_pi_overlap function."""

    def test_full_overlap(self, sample_coords_parallel):
        """Perfectly stacked planes should have high overlap."""
        coords1, coords2 = sample_coords_parallel
        overlap = calculate_pi_overlap(coords1, coords2)
        assert overlap > 80.0  # Should be close to 100%

    def test_no_overlap(self, sample_coords_offset):
        """Offset planes should have lower overlap."""
        coords1, coords2 = sample_coords_offset
        overlap = calculate_pi_overlap(coords1, coords2)
        assert overlap < 50.0  # Reduced overlap due to offset

    def test_overlap_range(self):
        """Overlap should be in [0, 100] range."""
        np.random.seed(42)
        for _ in range(20):
            coords1 = np.random.randn(6, 3)
            coords2 = np.random.randn(6, 3)
            overlap = calculate_pi_overlap(coords1, coords2)
            assert 0 <= overlap <= 100

    def test_without_shapely(self, sample_coords_parallel):
        """Test fallback method without Shapely."""
        coords1, coords2 = sample_coords_parallel
        overlap = calculate_pi_overlap(coords1, coords2, use_shapely=False)
        assert 0 <= overlap <= 100


class TestProjectionFunctions:
    """Tests for projection helper functions."""

    def test_create_projection_basis(self):
        """Test orthonormal basis creation."""
        normal = np.array([0.0, 0.0, 1.0])
        u, v = _create_projection_basis(normal)

        # u and v should be orthogonal to normal
        assert abs(np.dot(u, normal)) < 0.001
        assert abs(np.dot(v, normal)) < 0.001

        # u and v should be orthogonal to each other
        assert abs(np.dot(u, v)) < 0.001

        # u and v should be unit vectors
        assert abs(np.linalg.norm(u) - 1.0) < 0.001
        assert abs(np.linalg.norm(v) - 1.0) < 0.001

    def test_project_to_2d(self):
        """Test 2D projection."""
        coords = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
            ]
        )
        u = np.array([1.0, 0.0, 0.0])
        v = np.array([0.0, 1.0, 0.0])

        proj = _project_to_2d(coords, u, v)

        assert proj.shape == (3, 2)
        np.testing.assert_array_almost_equal(proj[0], [0.0, 0.0])
        np.testing.assert_array_almost_equal(proj[1], [1.0, 0.0])
        np.testing.assert_array_almost_equal(proj[2], [0.0, 1.0])


class TestCalculateTiltAngle:
    """Tests for calculate_tilt_angle function."""

    def test_aligned_planes(self, sample_coords_parallel):
        """Aligned parallel planes should have low tilt."""
        coords1, coords2 = sample_coords_parallel
        tilt = calculate_tilt_angle(coords1, coords2)
        assert tilt < 10.0

    def test_tilt_range(self):
        """Tilt angle should be in [0, 180] range."""
        np.random.seed(42)
        for _ in range(20):
            coords1 = np.random.randn(5, 3)
            coords2 = np.random.randn(5, 3)
            tilt = calculate_tilt_angle(coords1, coords2)
            assert 0 <= tilt <= 180
