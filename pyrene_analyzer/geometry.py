"""
Geometry Module
===============

Functions for calculating geometric properties of pyrene dimer systems.

This module provides core geometric calculations including:
- SVD-based plane fitting
- Plane-plane angle calculation
- Inter-plane distance measurement
- π-π overlap estimation using polygon intersection
- Centroid distance calculation
- Slip-stack displacement analysis

All functions are designed to work with numpy arrays of 3D coordinates.
"""

from typing import Tuple, Optional
import numpy as np
from scipy.spatial import ConvexHull

# Try to import Shapely for accurate polygon intersection
try:
    from shapely.geometry import Polygon
    from shapely.validation import make_valid
    SHAPELY_AVAILABLE = True
except ImportError:
    SHAPELY_AVAILABLE = False


def fit_plane_svd(coords: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Fit a plane to 3D coordinates using Singular Value Decomposition (SVD).
    
    The plane is defined by its centroid and normal vector. SVD is used to find
    the best-fit plane by minimizing the sum of squared distances from points
    to the plane.
    
    Args:
        coords: Nx3 array of atomic coordinates in Angstroms.
        
    Returns:
        Tuple containing:
            - centroid: 3D coordinates of the plane centroid
            - normal: Unit normal vector perpendicular to the plane
            
    Raises:
        ValueError: If coords has fewer than 3 points or wrong dimensions.
        
    Example:
        >>> coords = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]])
        >>> centroid, normal = fit_plane_svd(coords)
        >>> print(f"Normal: {normal}")  # Should be close to [0, 0, 1]
        
    Notes:
        The SVD approach finds the plane that minimizes the sum of squared
        perpendicular distances from all points to the plane. The normal
        vector corresponds to the smallest singular value.
        
    References:
        - Shakarji, C.M. (1998). Least-squares fitting algorithms of the NIST
          algorithm testing system. J. Res. Natl. Inst. Stand. Technol., 103, 633.
    """
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError("coords must be an Nx3 array")
    if coords.shape[0] < 3:
        raise ValueError("At least 3 points are required to fit a plane")
    
    # Calculate centroid
    centroid = coords.mean(axis=0)
    
    # Center the coordinates
    centered = coords - centroid
    
    # SVD decomposition
    # The normal vector is the right singular vector corresponding to
    # the smallest singular value
    _, s, vh = np.linalg.svd(centered)
    
    # The last row of vh is the normal vector
    normal = vh[2, :]
    
    # Ensure unit vector
    norm = np.linalg.norm(normal)
    if norm > 0:
        normal = normal / norm
    
    return centroid, normal


def calculate_plane_angle(
    coords1: np.ndarray,
    coords2: np.ndarray
) -> float:
    """
    Calculate the angle between two planes defined by atom coordinates.
    
    The angle is computed as the angle between the normal vectors of the
    two best-fit planes. The result is always in the range [0, 90] degrees.
    
    Args:
        coords1: Nx3 array of coordinates for the first plane.
        coords2: Mx3 array of coordinates for the second plane.
        
    Returns:
        Angle between planes in degrees (0-90°).
        
    Example:
        >>> # Two parallel planes
        >>> plane1 = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        >>> plane2 = np.array([[0, 0, 3], [1, 0, 3], [0, 1, 3]])
        >>> angle = calculate_plane_angle(plane1, plane2)
        >>> print(f"Angle: {angle:.1f}°")  # Should be ~0°
        
    Notes:
        For pyrene excimer formation:
        - θ < 20°: Strong excimer (parallel stacking)
        - θ = 20-60°: Weak excimer
        - θ > 60°: Monomer emission (perpendicular arrangement)
        
    References:
        - Birks, J.B. (1970). Photophysics of Aromatic Molecules.
    """
    _, normal1 = fit_plane_svd(coords1)
    _, normal2 = fit_plane_svd(coords2)
    
    # Calculate angle between normal vectors
    # Use absolute value since plane orientation is arbitrary
    cos_angle = np.abs(np.dot(normal1, normal2))
    
    # Clip for numerical stability
    cos_angle = np.clip(cos_angle, -1.0, 1.0)
    
    # Convert to degrees
    angle_rad = np.arccos(cos_angle)
    angle_deg = np.degrees(angle_rad)
    
    # Ensure acute angle (0-90°)
    if angle_deg > 90:
        angle_deg = 180 - angle_deg
    
    return float(angle_deg)


def calculate_interplane_distance(
    coords1: np.ndarray,
    coords2: np.ndarray
) -> float:
    """
    Calculate the perpendicular distance between two planes.
    
    The distance is measured as the perpendicular distance from the centroid
    of the second plane to the first plane.
    
    Args:
        coords1: Nx3 array of coordinates for the first plane.
        coords2: Mx3 array of coordinates for the second plane.
        
    Returns:
        Perpendicular distance in Angstroms.
        
    Example:
        >>> plane1 = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        >>> plane2 = np.array([[0, 0, 3.5], [1, 0, 3.5], [0, 1, 3.5]])
        >>> distance = calculate_interplane_distance(plane1, plane2)
        >>> print(f"Distance: {distance:.2f} Å")  # Should be ~3.5 Å
        
    Notes:
        Optimal π-stacking distances for pyrene excimers:
        - 3.3-3.7 Å: Strong excimer formation
        - 3.7-4.5 Å: Weak excimer
        - > 4.5 Å: Monomer emission
        
    References:
        - Stevens, B. (1968). Proc. Royal Soc. A, 305, 55-70.
          (Pyrene crystal excimer: d = 3.34 Å)
    """
    centroid1, normal1 = fit_plane_svd(coords1)
    centroid2, _ = fit_plane_svd(coords2)
    
    # Vector from centroid1 to centroid2
    vec = centroid2 - centroid1
    
    # Perpendicular distance is the projection onto the normal
    distance = np.abs(np.dot(vec, normal1))
    
    return float(distance)


def calculate_centroid_distance(
    coords1: np.ndarray,
    coords2: np.ndarray
) -> float:
    """
    Calculate the Euclidean distance between the centroids of two coordinate sets.
    
    This is a simpler metric than inter-plane distance and represents the
    direct center-to-center distance between two aromatic systems.
    
    Args:
        coords1: Nx3 array of coordinates for the first system.
        coords2: Mx3 array of coordinates for the second system.
        
    Returns:
        Centroid-to-centroid distance in Angstroms.
        
    Example:
        >>> ring1 = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        >>> ring2 = np.array([[5, 0, 0], [6, 0, 0], [5, 1, 0]])
        >>> dist = calculate_centroid_distance(ring1, ring2)
        >>> print(f"Distance: {dist:.2f} Å")
    """
    centroid1 = coords1.mean(axis=0)
    centroid2 = coords2.mean(axis=0)
    
    distance = np.linalg.norm(centroid2 - centroid1)
    
    return float(distance)


def calculate_slip_stack_displacement(
    coords1: np.ndarray,
    coords2: np.ndarray
) -> float:
    """
    Calculate the slip-stack (lateral) displacement between two planes.
    
    This measures how much the two ring systems are offset from perfect
    face-to-face stacking. A value of 0 indicates perfect eclipsed stacking.
    
    Args:
        coords1: Nx3 array of coordinates for the first plane.
        coords2: Mx3 array of coordinates for the second plane.
        
    Returns:
        Lateral displacement in Angstroms.
        
    Example:
        >>> # Perfectly stacked planes
        >>> plane1 = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        >>> plane2 = np.array([[0, 0, 3.5], [1, 0, 3.5], [0, 1, 3.5]])
        >>> slip = calculate_slip_stack_displacement(plane1, plane2)
        >>> print(f"Slip: {slip:.2f} Å")  # Should be ~0 Å
        
    Notes:
        The slip-stack displacement is calculated as the component of the
        centroid-centroid vector that lies in the plane (perpendicular to
        the normal vector).
        
    References:
        - Hunter, C.A. & Sanders, J.K.M. (1990). J. Am. Chem. Soc., 112, 5525.
    """
    centroid1, normal1 = fit_plane_svd(coords1)
    centroid2, _ = fit_plane_svd(coords2)
    
    # Vector from centroid1 to centroid2
    vec = centroid2 - centroid1
    
    # Project onto normal to get perpendicular component
    perpendicular = np.dot(vec, normal1) * normal1
    
    # Lateral component is the remainder
    lateral = vec - perpendicular
    
    return float(np.linalg.norm(lateral))


def _create_projection_basis(normal: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Create an orthonormal basis for projecting 3D points onto a 2D plane.
    
    Args:
        normal: Unit normal vector of the projection plane.
        
    Returns:
        Tuple of two orthonormal vectors (u, v) lying in the plane.
    """
    # Choose an arbitrary vector not parallel to normal
    if abs(normal[0]) < 0.9:
        arbitrary = np.array([1.0, 0.0, 0.0])
    else:
        arbitrary = np.array([0.0, 1.0, 0.0])
    
    # Create orthonormal basis using cross products
    u = np.cross(normal, arbitrary)
    u = u / np.linalg.norm(u)
    v = np.cross(normal, u)
    v = v / np.linalg.norm(v)
    
    return u, v


def _project_to_2d(
    coords: np.ndarray,
    u: np.ndarray,
    v: np.ndarray
) -> np.ndarray:
    """
    Project 3D coordinates onto a 2D plane defined by basis vectors u and v.
    
    Args:
        coords: Nx3 array of 3D coordinates.
        u: First basis vector of the projection plane.
        v: Second basis vector of the projection plane.
        
    Returns:
        Nx2 array of 2D projected coordinates.
    """
    coords_2d = np.zeros((len(coords), 2))
    for i, coord in enumerate(coords):
        coords_2d[i, 0] = np.dot(coord, u)
        coords_2d[i, 1] = np.dot(coord, v)
    return coords_2d


def _calculate_overlap_shapely(
    proj1: np.ndarray,
    proj2: np.ndarray
) -> Tuple[float, float, float]:
    """
    Calculate polygon overlap using Shapely library.
    
    Args:
        proj1: Nx2 array of 2D coordinates for first polygon.
        proj2: Mx2 array of 2D coordinates for second polygon.
        
    Returns:
        Tuple of (overlap_percentage, area1, area2).
    """
    try:
        # Create convex hull polygons
        hull1 = ConvexHull(proj1)
        hull2 = ConvexHull(proj2)
        
        # Get vertices in order
        vertices1 = proj1[hull1.vertices]
        vertices2 = proj2[hull2.vertices]
        
        # Create Shapely polygons
        poly1 = Polygon(vertices1)
        poly2 = Polygon(vertices2)
        
        # Make valid if needed
        if not poly1.is_valid:
            poly1 = make_valid(poly1)
        if not poly2.is_valid:
            poly2 = make_valid(poly2)
        
        # Calculate areas
        area1 = poly1.area
        area2 = poly2.area
        
        # Calculate intersection
        intersection = poly1.intersection(poly2)
        intersection_area = intersection.area
        
        # Calculate overlap as percentage of smaller ring
        min_area = min(area1, area2)
        if min_area > 0:
            overlap_pct = 100.0 * intersection_area / min_area
        else:
            overlap_pct = 0.0
        
        return overlap_pct, area1, area2
        
    except Exception:
        return None, None, None


def _calculate_overlap_convex_hull(
    proj1: np.ndarray,
    proj2: np.ndarray
) -> Tuple[float, float, float]:
    """
    Estimate polygon overlap using convex hull centroid distance method.
    
    This is a fallback method when Shapely is not available.
    
    Args:
        proj1: Nx2 array of 2D coordinates for first polygon.
        proj2: Mx2 array of 2D coordinates for second polygon.
        
    Returns:
        Tuple of (overlap_percentage, area1, area2).
    """
    try:
        hull1 = ConvexHull(proj1)
        hull2 = ConvexHull(proj2)
        area1 = hull1.volume  # In 2D, volume = area
        area2 = hull2.volume
    except Exception:
        return 50.0, None, None  # Default fallback
    
    # Estimate overlap using centroid distance method
    centroid1 = proj1.mean(axis=0)
    centroid2 = proj2.mean(axis=0)
    centroid_dist = np.linalg.norm(centroid1 - centroid2)
    
    # Rough overlap estimation based on centroid distance
    avg_radius = np.sqrt((area1 + area2) / (2 * np.pi))
    if centroid_dist < 2 * avg_radius:
        overlap_pct = 100.0 * (1 - centroid_dist / (2 * avg_radius))
    else:
        overlap_pct = 0.0
    
    return np.clip(overlap_pct, 0, 100), area1, area2


def calculate_pi_overlap(
    coords1: np.ndarray,
    coords2: np.ndarray,
    use_shapely: bool = True
) -> float:
    """
    Calculate the π-π overlap percentage between two aromatic ring systems.
    
    The overlap is calculated by:
    1. Projecting both ring systems onto a common plane (average of the two planes)
    2. Computing the convex hull of each projection
    3. Calculating the intersection area
    4. Returning the overlap as a percentage of the smaller ring area
    
    Args:
        coords1: Nx3 array of coordinates for the first ring system.
        coords2: Mx3 array of coordinates for the second ring system.
        use_shapely: If True and Shapely is available, use accurate polygon
                     intersection. Otherwise, use centroid-based approximation.
        
    Returns:
        Overlap percentage (0-100%).
        
    Example:
        >>> # Two overlapping pyrene rings
        >>> ring1_coords = get_pyrene_coords(mol, conf_id, ring1_atoms)
        >>> ring2_coords = get_pyrene_coords(mol, conf_id, ring2_atoms)
        >>> overlap = calculate_pi_overlap(ring1_coords, ring2_coords)
        >>> print(f"π-overlap: {overlap:.1f}%")
        
    Notes:
        For pyrene excimer formation:
        - > 70%: Strong excimer
        - 30-70%: Weak excimer
        - < 30%: Monomer emission
        
        The Shapely-based calculation provides accurate polygon intersection,
        while the fallback method uses a centroid-distance approximation.
        
    References:
        - Ge, Y. et al. (2020). J. Mater. Chem. C, 8, 10223-10232.
          (π-overlap is more important than distance for excimer formation)
    """
    # Fit planes to both coordinate sets
    _, normal1 = fit_plane_svd(coords1)
    _, normal2 = fit_plane_svd(coords2)
    
    # Use average normal for projection plane
    proj_normal = (normal1 + normal2) / 2
    norm = np.linalg.norm(proj_normal)
    if norm > 0:
        proj_normal = proj_normal / norm
    else:
        proj_normal = normal1  # Fallback if normals are opposite
    
    # Create projection basis
    u, v = _create_projection_basis(proj_normal)
    
    # Project coordinates onto 2D plane
    proj1 = _project_to_2d(coords1, u, v)
    proj2 = _project_to_2d(coords2, u, v)
    
    # Calculate overlap
    if use_shapely and SHAPELY_AVAILABLE:
        result = _calculate_overlap_shapely(proj1, proj2)
        if result[0] is not None:
            return float(np.clip(result[0], 0, 100))
    
    # Fallback to convex hull method
    overlap_pct, _, _ = _calculate_overlap_convex_hull(proj1, proj2)
    return float(np.clip(overlap_pct, 0, 100))


def calculate_tilt_angle(
    coords1: np.ndarray,
    coords2: np.ndarray
) -> float:
    """
    Calculate the tilt angle between two ring systems.
    
    The tilt angle measures the rotation of one ring relative to the other
    around the axis connecting their centroids.
    
    Args:
        coords1: Nx3 array of coordinates for the first ring system.
        coords2: Mx3 array of coordinates for the second ring system.
        
    Returns:
        Tilt angle in degrees (0-180°).
        
    Notes:
        This is different from the plane-plane angle (θ). The tilt angle
        measures rotational alignment, while θ measures angular deviation
        from parallel.
    """
    centroid1, normal1 = fit_plane_svd(coords1)
    centroid2, normal2 = fit_plane_svd(coords2)
    
    # Axis connecting centroids
    axis = centroid2 - centroid1
    axis_norm = np.linalg.norm(axis)
    if axis_norm > 0:
        axis = axis / axis_norm
    else:
        return 0.0
    
    # Project normals onto plane perpendicular to axis
    proj_normal1 = normal1 - np.dot(normal1, axis) * axis
    proj_normal2 = normal2 - np.dot(normal2, axis) * axis
    
    # Normalize
    norm1 = np.linalg.norm(proj_normal1)
    norm2 = np.linalg.norm(proj_normal2)
    
    if norm1 > 0 and norm2 > 0:
        proj_normal1 = proj_normal1 / norm1
        proj_normal2 = proj_normal2 / norm2
        
        cos_angle = np.dot(proj_normal1, proj_normal2)
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        angle = np.degrees(np.arccos(cos_angle))
        return float(angle)
    
    return 0.0
