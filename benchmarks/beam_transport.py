"""Analytical validation and convergence of the public 0.5 beam solver.

Run: python benchmarks/beam_transport.py
All CSVs, figures and a machine-readable pass/fail summary go to benchmarks/results/beams.
References are independent normal-incidence slab sums, Gaussian solutions and
paraxial mirror/interface ABCD laws, not calls back into solver coefficients.
"""
from pathlib import Path
import argparse
import csv
import json
import math
import platform
import sys
from datetime import datetime, timezone
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import raytracing as rt


def write_csv(path, rows):
    with path.open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def run(destination):
    if rt.solve.__module__ != 'raytracing.transport':
        raise RuntimeError('This validation requires the current Gaussian-beam transport API')
    destination.mkdir(parents=True,exist_ok=True)
    checks=[]
    def check(name,error,tolerance):
        checks.append(dict(name=name,max_error=float(error),tolerance=tolerance,passed=bool(error<=tolerance)))

    f=10e9; k=2*np.pi*f/rt.C0; w=.1; zr=k*w*w/2
    source=rt.GaussianBeam([0,0,1],[0,0,1],w,launch_distance=-1)
    result=rt.solve(rt.Scene(),source,rt.SolverConfig(frequency_hz=f,max_distance=4))
    e0=np.sqrt(4*rt.ETA0/(np.pi*w*w))
    rows=[]
    for z in np.linspace(.001,3.99,121):
        s=z-1; width=w*np.sqrt(1+(s/zr)**2)
        for radius in (0,.5*width,width,2*width):
            theory=e0/(1+1j*s/zr)*np.exp(1j*k*s)*np.exp(.5j*k*radius**2/(s-1j*zr))
            actual=result.field_at([radius,0,z])[0]
            rows.append(dict(z_m=z,radius_m=radius,actual_real=actual.real,actual_imag=actual.imag,
                             theory_real=theory.real,theory_imag=theory.imag,relative_error=abs(actual-theory)/abs(theory)))
    write_csv(destination/'gaussian_field.csv',rows)
    check('Gaussian complex field across waist',max(r['relative_error'] for r in rows),2e-11)

    # Numerical receiver quadrature must converge to independently normalized
    # source power. The missing tail outside 5w is exp(-50).
    rows=[]
    z=2.3; width=w*np.sqrt(1+((z-1)/zr)**2)
    for count in (9,17,33,65,129,257,513):
        r=np.linspace(0,5*width,count)
        intensity=np.array([np.vdot(v,v).real/(2*rt.ETA0) for v in result.fields_at([[x,0,z] for x in r])])
        integral=np.trapezoid(2*np.pi*r*intensity,r)
        rows.append(dict(radial_samples=count,power_w=integral,theory_w=1.0,absolute_error=abs(integral-1)))
    write_csv(destination/'gaussian_power_convergence.csv',rows)
    check('Gaussian receiver quadrature',rows[-1]['absolute_error'],4e-5)
    check('Gaussian quadrature error reduction',rows[-1]['absolute_error']/rows[0]['absolute_error'],.001)

    # Slab Fabry-Perot amplitude and loss, with both internal reflection branches.
    slab_rows=[]; slab_curves=[]
    thickness=.013
    for epsilon in (4+0j,4+.02j):
        scene=rt.Scene()
        scene.add_volume(rt.Box([0,0,.5+thickness/2],[2,2,thickness/2]),rt.Material.dielectric(epsilon))
        n=np.sqrt(epsilon); p=np.exp(1j*k*n*thickness)
        theory=4*n/(1+n)**2*p/(1-((n-1)/(n+1))**2*p*p)*np.exp(1j*k*(1-thickness))
        errors=[]
        for depth in (2,4,6,8,12,16,24,32):
            r=rt.solve(scene,rt.PlaneWave([0,0,0],[0,0,1]),
                       rt.SolverConfig(frequency_hz=f,max_distance=4,max_interactions=depth,min_power_fraction=0))
            actual=r.field_at([0,0,1])[0]
            error=abs(actual-theory)/abs(theory)
            balance=r.power_balance
            residual=abs(sum(v for key,v in balance.items() if key!='launched')-balance['launched'])/balance['launched']
            slab_rows.append(dict(epsilon_imag=epsilon.imag,max_interactions=depth,relative_error=error,
                                  actual_real=actual.real,actual_imag=actual.imag,theory_real=theory.real,
                                  theory_imag=theory.imag,power_balance_error=residual,
                                  truncated_power=balance.get('interaction_limit',0)))
            errors.append(error)
        check(f'Slab coherent amplitude eps={epsilon}',errors[-1],2e-11)
        check(f'Slab convergence eps={epsilon}',errors[-1]/errors[0],1e-8)
        slab_curves.append((str(epsilon),[2,4,6,8,12,16,24,32],errors))
    write_csv(destination/'slab_convergence.csv',slab_rows)
    check('Slab branch energy accounting',max(r['power_balance_error'] for r in slab_rows),2e-12)

    # Finite Gaussian Fabry-Perot transmission: each round-trip order has
    # a different ABCD propagation length as well as a different carrier phase.
    # Sum 40 analytical orders independently of the solver branch tree.
    finite_slab=[]
    scene=rt.Scene()
    scene.add_volume(rt.Box([0,0,.5+thickness/2],[20,20,thickness/2]),rt.Material.dielectric(4))
    finite=rt.solve(scene,rt.GaussianBeam([0,0,0],[0,0,1],w),
                    rt.SolverConfig(frequency_hz=f,max_distance=4,max_interactions=32,min_power_fraction=0))
    for radius_sample in (0,.03,.08,.15):
        expected=0j
        q0=-1j*zr
        for order in range(40):
            interior=(2*order+1)*thickness
            q=q0+1-thickness+interior/2
            expected += (8/9)*(1/9)**order*q0/q*np.exp(1j*k*(1-thickness+2*interior)) * np.exp(.5j*k*radius_sample**2/q)
        actual=finite.field_at([radius_sample,0,1])[0]
        # The independent series above is normalized to the incident waist field.
        expected *= e0
        finite_slab.append(dict(radius_m=radius_sample,actual_real=actual.real,actual_imag=actual.imag,
                                theory_real=expected.real,theory_imag=expected.imag,
                                relative_error=abs(actual-expected)/abs(expected)))
    write_csv(destination/'gaussian_slab.csv',finite_slab)
    check('Finite Gaussian slab coherent transmission',max(row['relative_error'] for row in finite_slab),2e-11)

    # Separate pruning convergence from interaction-depth convergence.
    scene=rt.Scene(); scene.add_volume(rt.Box([0,0,.5065],[2,2,.0065]),rt.Material.dielectric(4))
    n=2; p=np.exp(1j*k*n*thickness)
    theory=4*n/(1+n)**2*p/(1-((n-1)/(n+1))**2*p*p)*np.exp(1j*k*(1-thickness))
    pruning=[]
    for cutoff in (1e-2,1e-4,1e-6,1e-8,1e-10,1e-14,0):
        r=rt.solve(scene,rt.PlaneWave([0,0,0],[0,0,1]),
                   rt.SolverConfig(frequency_hz=f,max_distance=4,max_interactions=32,min_power_fraction=cutoff))
        pruning.append(dict(power_cutoff=cutoff,relative_error=abs(r.field_at([0,0,1])[0]-theory)/abs(theory),
                            segments=len(r.segments),discarded_power=r.power_balance.get('power_cutoff',0)))
    write_csv(destination/'branch_pruning_convergence.csv',pruning)
    check('Pruning convergence',pruning[-1]['relative_error'],2e-11)

    # Angle sweep: independent lossless textbook Fresnel + Snell and R+T=1.
    rows=[]
    for degrees in (0,10,30,45,math.degrees(math.atan(2)),70,80):
        theta=math.radians(degrees); ci=math.cos(theta); st=math.sin(theta)/2; ct=math.sqrt(1-st*st)
        for pol in ('s','p'):
            d=[math.sin(theta),0,ci]; e=[0,1,0] if pol=='s' else [ci,0,-math.sin(theta)]
            scene=rt.Scene(); scene.add_volume(rt.Box([0,0,1],[20,20,.5]),rt.Material.dielectric(4))
            r=rt.solve(scene,rt.PlaneWave([0,0,0],d,e),
                       rt.SolverConfig(frequency_hz=f,max_distance=5,max_interactions=1,min_power_fraction=0))
            children=[s.state for s in r.segments if len(s.state.history)==1]
            reflected=next(s for s in children if s.history[-1][1]=='reflect')
            transmitted=next(s for s in children if s.history[-1][1]=='transmit')
            expected_r=(ci-2*ct)/(ci+2*ct) if pol=='s' else (2*ci-ct)/(2*ci+ct)
            R=reflected.power/r.launched_power; T=transmitted.power/r.launched_power
            rows.append(dict(degrees=degrees,polarization=pol,R=R,T=T,theory_R=expected_r**2,
                             reflection_error=abs(R-expected_r**2),balance_error=abs(R+T-1),
                             snell_error=abs(transmitted.direction[0]-st)))
    write_csv(destination/'fresnel_angles.csv',rows)
    check('Fresnel angle sweep',max(r['reflection_error'] for r in rows),2e-12)
    check('Snell refraction',max(r['snell_error'] for r in rows),2e-12)
    check('Lossless R+T',max(r['balance_error'] for r in rows),2e-12)

    # Concave spherical mirror: ABCD q_out^-1=q_in^-1-2/R.
    # A PEC spherical surface is an opaque shell, permitting an interior source.
    radius=2; scene=rt.Scene(); scene.add_surface(rt.Sphere([0,0,0],radius),rt.Material.pec())
    source=rt.GaussianBeam([0,0,0],[0,0,1],w)
    r=rt.solve(scene,source,rt.SolverConfig(frequency_hz=f,max_distance=3.99,max_interactions=1))
    reflected=next(s for s in r.segments if s.state.history)
    q0=-1j*zr; qi=q0+radius; qo=1/(1/qi-2/radius)
    rows=[]
    for distance in np.linspace(0,1.98,100):
        from raytracing.beams import propagate
        actual=propagate(reflected.state,distance,rt.VACUUM,f)
        theory=-e0*q0/qi*qo/(qo+distance)*np.exp(1j*k*(radius+distance))
        expected_width=np.sqrt(2/(k*(1/(qo+distance)).imag))
        width=np.sqrt(2/(k*np.linalg.eigvalsh(actual.curvature.imag)[0]))
        rows.append(dict(distance_after_mirror_m=distance,width_m=width,theory_width_m=expected_width,
                         relative_field_error=abs(actual.field[0]-theory)/abs(theory),width_error=abs(width-expected_width)))
    write_csv(destination/'concave_mirror.csv',rows)
    check('Concave mirror complex field and Gouy phase',max(r['relative_field_error'] for r in rows),2e-11)
    check('Concave mirror beam width',max(r['width_error'] for r in rows),2e-12)

    # Oblique spherical reflection: distinct tangential and sagittal powers.
    rows=[]
    for degrees in (0,15,30,45,60):
        theta=math.radians(degrees); point=np.array([radius*math.sin(theta),0,radius*math.cos(theta)])
        origin=point-np.array([0,0,.3])
        r=rt.solve(scene,rt.GaussianBeam(origin,[0,0,1],w),
                   rt.SolverConfig(frequency_hz=f,max_distance=.8,max_interactions=1))
        out=next(s.state for s in r.segments if s.state.history)
        hi=1/(.3-1j*zr)
        actual=np.sort(np.linalg.eigvalsh(out.curvature.real))
        expected=np.sort([hi.real-2/(radius*math.cos(theta)),hi.real-2*math.cos(theta)/radius])
        rows.append(dict(degrees=degrees,curvature_1=actual[0],curvature_2=actual[1],
                         theory_1=expected[0],theory_2=expected[1],max_error=max(abs(actual-expected))))
    write_csv(destination/'astigmatic_mirror.csv',rows)
    check('Oblique spherical mirror astigmatism',max(r['max_error'] for r in rows),2e-11)

    # Finite Gaussian beams through repeated spherical/cylindrical foci.
    # Independent scalar ABCD recurrence in laboratory x/y coordinates:
    # free propagation q -> q+L, mirror 1/q -> 1/q-2/R in
    # each curved plane. A cylinder has no power along its axis (y).
    repeated=[]
    for shape in ('sphere','cylinder'):
        mirror_scene=rt.Scene()
        surface=(rt.Sphere([0,0,0],radius) if shape=='sphere' else
                 rt.Cylinder([0,0,0],[0,1,0],radius,1000,capped=False))
        mirror_scene.add_surface(surface,rt.Material.pec())
        simulation=rt.solve(mirror_scene,rt.GaussianBeam([0,0,0],[0,0,1],w),
                            rt.SolverConfig(frequency_hz=f,max_distance=30.5,max_interactions=8))
        outgoing=sorted((s for s in simulation.segments if s.state.history),
                        key=lambda s: len(s.state.history))
        if len(outgoing)!=8:
            raise RuntimeError(f'{shape}: expected eight reflected Gaussian segments')
        q=np.full(2,-1j*zr); amplitude=complex(e0)
        powers=np.array([2/radius,2/radius if shape=='sphere' else 0])
        for order,segment in enumerate(outgoing,1):
            length=radius if order==1 else 2*radius
            amplitude *= -np.exp(1j*k*length)/np.prod(np.sqrt(1+length/q))
            q=1/(1/(q+length)-powers)
            for distance in (0,.25):
                # Query a traced segment's finite field, including its Gaussian
                # transverse envelope, rather than a native point-ray field.
                qp=q+distance
                expected=amplitude*np.exp(1j*k*distance)/np.prod(np.sqrt(1+distance/q))
                point=segment.state.position+distance*segment.state.direction+np.array([.02,.01,0])
                expected *= np.exp(.5j*k*(.02**2/qp[0]+.01**2/qp[1]))
                actual=segment.field_at(point)[0]
                # Eigenvalue sorting avoids the solver's transverse-frame choice.
                actual_width=np.sort(segment.radii(distance))
                theory_width=np.sort(np.sqrt(2/(k*(1/qp).imag)))
                repeated.append(dict(shape=shape,reflections=order,distance_after_mirror_m=distance,
                                     actual_real=actual.real,actual_imag=actual.imag,
                                     theory_real=expected.real,theory_imag=expected.imag,
                                     relative_field_error=abs(actual-expected)/abs(expected),
                                     width_error=max(abs(actual_width-theory_width))))
    write_csv(destination/'repeated_gaussian_mirrors.csv',repeated)
    check('Repeated spherical/cylindrical Gaussian focus complex field',
          max(row['relative_field_error'] for row in repeated),2e-10)
    check('Repeated spherical/cylindrical Gaussian focus widths',
          max(row['width_error'] for row in repeated),2e-11)

    # Convex spherical dielectric entry: scalar refraction ABCD law.
    # R>0 has its center downstream; q_out^-1=(n1/n2)/q_in
    # -(n2-n1)/(n2*R). Test a finite beam, not a plane wave.
    refraction=[]
    for index in (1.25,1.5,2.0):
        dielectric_scene=rt.Scene()
        dielectric_scene.add_volume(rt.Sphere([0,0,0],radius),rt.Material.dielectric(index**2))
        simulation=rt.solve(dielectric_scene,rt.GaussianBeam([0,0,-radius-.3],[0,0,1],w),
                            rt.SolverConfig(frequency_hz=f,max_distance=.8,max_interactions=1))
        transmitted=next(s for s in simulation.segments
                         if s.state.history and s.state.history[-1][1]=='transmit')
        qi=.3-1j*zr
        qo=1/((1/index)/qi-(index-1)/(index*radius))
        entry=e0*(-1j*zr)/qi*np.exp(1j*k*.3)*2/(1+index)
        for distance in np.linspace(0,.45,20):
            expected=entry*qo/(qo+distance)*np.exp(1j*k*index*distance)
            actual=transmitted.field_at(transmitted.state.position+distance*np.array([0,0,1]))[0]
            expected_width=np.sqrt(2/(k*index*(1/(qo+distance)).imag))
            refraction.append(dict(refractive_index=index,distance_m=distance,
                                   actual_real=actual.real,actual_imag=actual.imag,
                                   theory_real=expected.real,theory_imag=expected.imag,
                                   width_m=transmitted.radii(distance)[0],theory_width_m=expected_width,
                                   relative_field_error=abs(actual-expected)/abs(expected),
                                   width_error=max(abs(transmitted.radii(distance)-expected_width))))
    write_csv(destination/'gaussian_curved_refraction.csv',refraction)
    check('Gaussian curved refraction complex field',max(row['relative_field_error'] for row in refraction),2e-11)
    check('Gaussian curved refraction width',max(row['width_error'] for row in refraction),2e-12)

    # Analytical normal-incidence absorbing-medium attenuation, independent of
    # interface amplitudes. Launch a plane wave inside the volume.
    eps=4+.02j; n=np.sqrt(eps)
    scene=rt.Scene(background=rt.Material.dielectric(eps))
    r=rt.solve(scene,rt.PlaneWave([0,0,0],[0,0,1]),rt.SolverConfig(frequency_hz=f,max_distance=2))
    rows=[]
    for distance in np.linspace(0,1.99,80):
        actual=r.field_at([0,0,distance])[0]; expected=np.exp(1j*k*n*distance)
        rows.append(dict(distance_m=distance,amplitude=abs(actual),theory_amplitude=abs(expected),error=abs(actual-expected)))
    write_csv(destination/'absorption.csv',rows)
    check('Homogeneous absorption and phase',max(r['error'] for r in rows),2e-12)

    # Coherent Gaussian beamlet decomposition, sampled independently at receiver
    # points away from the fitting/launch plane. This tests spatial convergence.
    source=rt.GaussianBeam([0,0,0],[0,0,1],w)
    points=[[x,y,z] for x,y in ((0,0),(.04,.01),(.1,-.03)) for z in (.2,.7,1.3)]
    reference=[]
    for x,y,z in points:
        reference.append(e0/(1+1j*z/zr)*np.exp(1j*k*z)*np.exp(.5j*k*(x*x+y*y)/(z-1j*zr)))
    rows=[]
    for samples in (5,9,13,17,25):
        r=rt.solve(rt.Scene(),source.beamlets(samples),rt.SolverConfig(frequency_hz=f,max_distance=2,min_power_fraction=0))
        actual=r.fields_at(points)[:,0]
        error=np.linalg.norm(actual-reference)/np.linalg.norm(reference)
        rows.append(dict(samples_per_axis=samples,beamlets=samples*samples,relative_complex_field_error=error))
    write_csv(destination/'beamlet_convergence.csv',rows)
    check('Gaussian beamlet spatial convergence',rows[-1]['relative_complex_field_error'],2e-6)
    check('Gaussian beamlet error reduction',rows[-1]['relative_complex_field_error']/rows[0]['relative_complex_field_error'],1e-4)

    fig,axes=plt.subplots(1,3,figsize=(14,4))
    power=list(csv.DictReader((destination/'gaussian_power_convergence.csv').open()))
    axes[0].loglog([int(x['radial_samples']) for x in power],[float(x['absolute_error']) for x in power],'o-')
    axes[0].set(xlabel='Radial receiver samples',ylabel='Absolute power error (W)',title='Gaussian quadrature convergence')
    for label,depths,errors in slab_curves: axes[1].semilogy(depths,errors,'o-',label=label)
    axes[1].set(xlabel='Maximum interface interactions',ylabel='Relative complex field error',title='Coherent slab convergence')
    axes[1].legend(title='Relative permittivity')
    axes[2].semilogy([r['samples_per_axis'] for r in rows],[r['relative_complex_field_error'] for r in rows],'o-')
    axes[2].set(xlabel='Beamlets per transverse axis',ylabel='Relative complex field error',title='Gaussian beamlet convergence')
    for ax in axes: ax.grid(True,which='both',alpha=.3)
    fig.tight_layout(); fig.savefig(destination/'convergence.png',dpi=180); plt.close(fig)
    summary=dict(schema=2,solver_version=rt.__version__,solver_api='raytracing.transport.solve',
                 solver_file=rt.solve.__code__.co_filename,source_types=['GaussianBeam','PlaneWave'],
                 python=sys.version,platform=platform.platform(),
                 utc=datetime.now(timezone.utc).isoformat(),checks=checks,passed=all(x['passed'] for x in checks),
                 assumptions=['paraxial Gaussian beams','homogeneous isotropic regions','weak-loss transport',
                              'opaque PEC/metal','no aperture clipping or edge diffraction'])
    (destination/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    lines=['# Beam transport analytical results','',f'Solver {rt.__version__}; all checks passed: {summary["passed"]}.',
           '', 'API: `raytracing.transport.solve`; finite Gaussian beams and plane-wave interface limits.', '',
           '| Check | Maximum error | Limit | Pass |','| --- | ---: | ---: | --- |']
    lines += [f'| {x["name"]} | {x["max_error"]:.4g} | {x["tolerance"]:.4g} | {x["passed"]} |' for x in checks]
    (destination/'REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    for x in checks: print(f'{"PASS" if x["passed"] else "FAIL"}: {x["name"]}: {x["max_error"]:.4g}')
    if not summary['passed']: raise SystemExit('Analytical benchmark failure; inspect summary.json')


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=Path(__file__).resolve().parent/'results'/'beams')
    run(parser.parse_args().output)
