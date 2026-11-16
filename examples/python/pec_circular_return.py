"""Independent SBR reference for a PEC circular pipe, including its bottom rim.

Unfold the planar PEC reflector at z=0. The native SBR solver searches an
analytical cylinder between Tx=(a,0,d) and image Rx=(a,0,-d), without analytical
return directions or circular-billiard roots. Folding z to |z| restores the
bottom reflection, including simultaneous cylinder/bottom events at the rim.

For the x-polarized e_theta antennas used here, the PEC electric-field image
requires a minus sign on the virtual receive pattern. The cylinder curvature,
neighboring-ray spreading and caustic counts are computed by the native solver.
This construction applies to the planar PEC end reflector used in Section IV.
"""
import numpy as np

ORIENTATION = np.array([[0., 0., 1.], [0., 1., 0.], [-1., 0., 0.]])


def solve_pec_circular_return(rt, depth, radius, offset, frequency,
                              ray_count, wall_limit, max_length, reception_radius):
    if not (depth > 0 and radius > 0 and 0 <= offset < radius):
        raise ValueError('Require d>0 and 0<=a<r.')
    scene = rt.Scene()
    scene.add(rt.Cylinder([0., 0., 0.], [0., 0., 1.], radius,
                          depth + .02, False), rt.Material.pec())
    tx_pattern = rt.Antenna(rt.Isotropic(rt.Polarization.VERTICAL, 1+0j))
    rx_pattern = rt.Antenna(rt.Isotropic(rt.Polarization.VERTICAL, 1+0j))
    tx_pattern.orientation = ORIENTATION
    # A 180-degree z rotation reverses the x-polarized e_theta receive vector.
    # This is the exact PEC image polarity, not a fitted phase correction.
    rx_pattern.orientation = np.diag([-1., -1., 1.]) @ ORIENTATION
    return rt.solve(scene, rt.Transmitter([offset, 0., depth], tx_pattern),
                    rt.Receiver([offset, 0., -depth], rx_pattern),
                    rt.SolverConfig(frequency_hz=frequency, ray_count=ray_count,
                                    max_reflections=wall_limit,
                                    max_distance=max_length+reception_radius,
                                    reception_radius=reception_radius))


def fold_vertices(vertices, tolerance=1e-10):
    """Fold a monotone crossing path and insert its planar reflector hit."""
    vertices = np.asarray(vertices, dtype=float)
    if vertices[0, 2] <= 0 or vertices[-1, 2] >= 0:
        raise ValueError('Unfolded path must cross z=0 from positive to negative.')
    points = [vertices[0].copy()]
    for first, second in zip(vertices[:-1], vertices[1:]):
        if first[2] > tolerance and second[2] < -tolerance:
            fraction = first[2]/(first[2]-second[2])
            points.append(first + fraction*(second-first))
        points.append(second.copy())
    result = np.asarray(points)
    result[:, 2] = np.abs(result[:, 2])
    return result
