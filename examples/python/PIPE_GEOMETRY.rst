Returned rays beyond a straight circular pipe
============================================

Run ``python examples/python/pipe_geometry_extensions.py --output geometry_extensions``
after installing this branch with the visualization dependencies. The example
uses unseeded SBR discovery at 77 GHz in five geometries: a tilted circular pipe
above a horizontal surface, square and elliptical pipes, an upper bore
transition, and a horizontal pipe with a wall port. ``EllipticCylinder`` and
``Annulus`` are exact analytical surfaces; rendering tessellation does not enter
the intersection or spreading calculations. The planar surfaces are PEC in
these geometric examples. Measurement performance for these configurations
has not been evaluated.

Why returns are discrete
------------------------

Fix a finite sequence of smooth specular reflections, excluding grazing hits
and edges. Parameterize the launch direction by two local angular coordinates
``alpha=(alpha_1,alpha_2)``. Near a candidate return, continue the final segment
to a plane through the receiver, perpendicular to that return's final
direction. Let ``F(alpha)`` be the two transverse miss coordinates in this
plane. A returned ray satisfies ``F(alpha)=0``.

If the two-by-two derivative ``J=dF/dalpha`` is nonsingular at a root, the
inverse function theorem makes that root isolated. This is the relevant
noncaustic condition: it applies to each returned path at its endpoint, not
simply to a source position considered independently of the ray. The same
implicit-function argument means that sufficiently small changes of smooth
geometry or source position preserve a nearby isolated return, until a
caustic, grazing event, edge, or visibility boundary is reached.

This argument applies to the smooth portions of all five geometries. It
establishes local isolation of an existing return; it does not guarantee that
every geometry has a detectable liquid-dependent return. A finite path-length
and reflection-order window is used in the examples. No claim of exhaustive
discovery is made for their illustrated subsets. Finite receiver aperture and
measurement resolution broaden ideal point-return directions and delays.

An ellipse provides a useful exception. For semiaxes ``a>b`` and foci
``F1,F2``, a ray from F1 reflects through F2; after the next reflection it
returns to F1. The two ellipse identities ``|P-F1|+|P-F2|=2a`` give a closed
two-reflection transverse length of ``4a`` for a continuum of launch angles.
With an axial planar end reflector this family also has constant 3D length
``sqrt((4a)^2+(2d)^2)``. Its endpoint Jacobian is singular. Merely moving away
from the pipe axis therefore does not exclude caustics: a source at an ellipse
focus is still exceptional. ``test_geometry_extensions`` checks the focal
closure and length against this independent analytical result. The illustrated
elliptical source is away from the foci, and every retained path passes the
Jacobian check.

For a uniform horizontal pipe and a horizontal liquid plane, every surface
normal has zero axial component. Specular reflection therefore conserves the
axial direction component. The net axial displacement is that component
multiplied by the total path length, so a colocated return must have zero
axial component and lie in one cross-sectional plane. In this geometry the
sensing variable is liquid height. It is not the axial distance to an end
reflector. The script checks that the computed paths remain in this plane.

Verification and outputs
------------------------

The script records each path's vertices, reflection sequence, complex field,
receiver miss, specular-law residual, endpoint Jacobian, and singular values.
It requires finite coefficients, receiver miss below 1e-8 m, specular residual
below 1e-8, and smallest Jacobian singular value above 1e-7 m. The figures show
up to three short paths with different reflection sequences in each geometry.
All discovered valid liquid-dependent paths are retained in the JSON files.
The manifest records numerical settings and hashes of the script and binary.

Common sensing chain
--------------------

Panel (f) follows the same chain as the circular-pipe study: returned-ray
model, complex channel impulse response, radar range processing, and distance
or level inference. The range-processing stage uses the transmitted waveform
and receiver processing; it is not restricted to FMCW.

Moving surfaces
---------------

For moving surfaces, a proposed extension uses short measurement bursts,
framewise geometric inference, and robust temporal combination of inferred
levels or likelihoods. Slow
level dynamics can constrain frame-to-frame ambiguities; identifiable surface
tilt can be included as a nuisance parameter. Averaging raw complex fields
can cancel moving echoes. Averaging inferred levels or incoherent power under
approximately constant level is more appropriate, but does not remove a
systematic geometry bias. These are proposed extensions, not simulated
accuracy results.
