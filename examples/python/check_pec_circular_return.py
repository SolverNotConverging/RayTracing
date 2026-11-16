"""Independent checks of PEC reflector unfolding and corner retroreflection."""
from pathlib import Path
import json
import numpy as np
import raytracing as rt
from pec_circular_return import solve_pec_circular_return, fold_vertices, ORIENTATION


def main():
    frequency = 77e9
    antenna = rt.Antenna(rt.Isotropic(rt.Polarization.VERTICAL, 1+0j))
    antenna.orientation = ORIENTATION
    image = rt.Antenna(rt.Isotropic(rt.Polarization.VERTICAL, 1+0j))
    image.orientation = np.diag([-1.,-1.,1.])@ORIENTATION
    mirror = np.diag([1.,1.,-1.])
    for direction in np.random.default_rng(1601).normal(size=(100,3)):
        direction /= np.linalg.norm(direction)
        np.testing.assert_allclose(image.farfield(direction,frequency),
                                   -mirror@antenna.farfield(mirror@direction,frequency),atol=1e-14)
    rows=[]
    for depth,offset,launches in [(2.,.3,4000),(5.,.3,16000),(5.,.1,8000)]:
        radius=.5
        limit=2.5*depth
        order=int(np.floor((1.5*depth+2*offset)/(2*np.sqrt(radius**2-offset**2))))
        scene=rt.Scene()
        scene.add(rt.Cylinder([0,0,(depth+.02)/2],[0,0,1],radius,(depth+.02)/2,False),rt.Material.pec())
        scene.add(rt.Disk([0,0,0],[0,0,1],radius),rt.Material.pec())
        physical=rt.solve(scene,rt.Transmitter([offset,0,depth],antenna),
                          rt.Receiver([offset,0,depth],antenna),rt.SolverConfig(
                              frequency_hz=frequency,ray_count=launches,max_reflections=order+1,
                              max_distance=limit+.48*radius,reception_radius=.48*radius))
        unfolded=solve_pec_circular_return(rt,depth,radius,offset,frequency,launches,order+1,limit,.48*radius)
        old=[p for p in physical.paths if p.field is not None and p.field.path_length<=limit+1e-10
             and sum(h.surface_index==1 for h in p.refinement.geometry.reflections)==1]
        new=[p for p in unfolded.paths if p.field is not None and p.field.path_length<=limit+1e-10]
        errors=[];lengths=[]
        floor=1e-8*max(abs(p.coefficient) for p in old)
        for path in old:
            match=min(new,key=lambda p:np.linalg.norm(path.refinement.launch_direction-p.refinement.launch_direction))
            assert np.linalg.norm(path.refinement.launch_direction-match.refinement.launch_direction)<1e-5
            errors.append(abs(path.coefficient-match.coefficient)/max(abs(path.coefficient),floor))
            lengths.append(abs(path.field.path_length-match.field.path_length))
            assert path.field.caustic_count==match.field.caustic_count
        assert max(errors)<1e-4 and max(lengths)<1e-8
        for x in [-radius,radius]:
            points=fold_vertices([[offset,0,depth],[x,0,0],[offset,0,-depth]])
            np.testing.assert_allclose(points[0],points[-1],atol=1e-14)
            incoming=points[1]-points[0];incoming/=np.linalg.norm(incoming)
            outgoing=incoming.copy();outgoing[[0,2]]*=-1
            np.testing.assert_allclose(outgoing,-incoming,atol=1e-14)
        rows.append(dict(depth_m=depth,offset_m=offset,smooth_paths=len(old),all_paths=len(new),
                         max_relative_complex_error=max(errors),max_length_error_m=max(lengths)))
    print(json.dumps(dict(passed=True,polarization_directions=100,corner_law='passed',cases=rows),indent=2))


if __name__=='__main__':
    main()
