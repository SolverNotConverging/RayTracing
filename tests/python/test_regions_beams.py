import math
import numpy as np
import pytest
import raytracing as rt


def cfg(**kwargs):
    return rt.SolverConfig(frequency_hz=10e9,max_distance=4,**kwargs)


def wall(z=1, x=0, half=2):
    return rt.Rectangle([x,0,z],[1,0,0],[0,1,0],half,half)


def plane(position=(0,0,0), direction=(0,0,1), pol=(1,0,0)):
    return rt.PlaneWave(position,direction,pol)


def slab(epsilon=4, thickness=.013):
    s = rt.Scene()
    s.add_volume(rt.Box([0,0,.5+thickness/2],[2,2,thickness/2]),rt.Material.dielectric(epsilon))
    return s


def test_material_validation_and_surface_rules():
    with pytest.raises(ValueError): rt.Material.dielectric(4-.1j)
    with pytest.raises(ValueError): rt.Material.dielectric(complex(float('nan'),0))
    with pytest.raises(ValueError): rt.Scene().add_surface(wall(),rt.Material.dielectric(4))
    with pytest.raises(ValueError): rt.Scene().add_volume(wall(),rt.Material.pec())
    with pytest.raises(ValueError): rt.Scene().add_volume(rt.Cylinder([0,0,0],[0,0,1],1,1,False),rt.Material.pec())
    assert rt.Material.dielectric(4+.02j).permittivity(10e9)==4+.02j
    with pytest.raises(ValueError): rt.solve(slab(4+1j),plane(),cfg())


def test_last_volume_wins_and_input_is_copied():
    s=rt.Scene()
    box=rt.Box([0,0,0],[2,2,2])
    a=s.add_volume(box,rt.Material.dielectric(4))
    b=s.add_volume(rt.Sphere([0,0,0],1),rt.Material.dielectric(9))
    box.center[:]=10
    assert s.region_at([0,0,0])==b
    assert s.region_at([1.5,0,0])==a
    assert s.region_at([3,0,0])==-1
    exposed=s.objects
    exposed[0].shape.center[:]=30
    assert s.region_at([1.5,0,0])==a


def test_overlap_priority_and_touching_edge():
    s=rt.Scene()
    s.add_surface(wall(),rt.Material.pec())
    s.add_surface(wall(half=1),rt.Material.metal(1e4))
    event=s.next_event(np.zeros(3),np.array([0.,0,1]),3)
    assert event.object_id==1 and event.kind=='interface'
    event=s.next_event(np.array([1.,0,0]),np.array([0.,0,1]),3)
    assert event.kind=='ambiguous'
    s=rt.Scene()
    s.add_surface(wall(x=-1,half=1),rt.Material.pec())
    s.add_surface(wall(x=1,half=1),rt.Material.pec())
    assert rt.solve(s,plane(),cfg()).terminals[0].reason=='ambiguous'


def test_shared_volume_face_is_one_interface():
    s=rt.Scene()
    s.add_volume(rt.Box([0,0,0],[1,1,.5]),rt.Material.dielectric(4))
    s.add_volume(rt.Box([0,0,1],[1,1,.5]),rt.Material.dielectric(9))
    event=s.next_event(np.zeros(3),np.array([0.,0,1]),2)
    assert event.kind=='interface' and (event.before,event.after)==(0,1)


def test_identical_nested_regions_are_invisible():
    material=rt.Material.dielectric(4)
    s=rt.Scene(background=material)
    s.add_volume(rt.Sphere([0,0,1],.5),material)
    s.add_volume(rt.Box([0,0,1],[.2,.2,.2]),material)
    r=rt.solve(s,plane(),cfg())
    assert len(r.segments)==1 and r.segments[0].length==4
    s.add_surface(wall(1),rt.Material.pec())
    r=rt.solve(s,plane(),cfg(max_interactions=1))
    reflected=next(seg for seg in r.segments if seg.state.history)
    assert reflected.state.history[0][2]==1


def test_opaque_receiver_location_has_no_field():
    s=rt.Scene()
    s.add_volume(rt.Box([.15,0,.5],[.03,.03,.1]),rt.Material.pec())
    r=rt.solve(s,rt.GaussianBeam([0,0,0],[0,0,1],.1),cfg())
    np.testing.assert_allclose(r.field_at([.15,0,.5]),0)


def test_curved_dielectric_is_entered_and_exited():
    for shape in (rt.Sphere([0,0,1],.3),rt.Cylinder([0,0,1],[0,1,0],.3,1,True)):
        s=rt.Scene()
        s.add_volume(shape,rt.Material.dielectric(2.25))
        r=rt.solve(s,rt.GaussianBeam([0,0,0],[0,0,1],.1),cfg(max_interactions=3))
        straight=[segment for segment in r.segments if len(segment.state.history)==2 and
                  all(hit[1]=='transmit' for hit in segment.state.history)]
        assert len(straight)==1 and straight[0].state.region==-1
        assert np.linalg.eigvalsh(straight[0].state.curvature.imag).min()>0


@pytest.mark.parametrize('mode',['surface','entry','exit','inside'])
def test_termination(mode):
    s=rt.Scene()
    if mode=='surface': s.add_termination(wall())
    elif mode=='exit': s.add_termination(rt.Box([0,0,0],[1,1,1]),volume=True,on='exit')
    elif mode=='inside': s.add_termination(rt.Sphere([0,0,0],1),volume=True)
    else: s.add_termination(rt.Box([0,0,1.5],[1,1,.5]),volume=True)
    s.add_surface(wall(2.5),rt.Material.pec())
    r=rt.solve(s,plane(),cfg())
    assert len(r.terminals)==1 and r.terminals[0].reason=='terminated'
    assert len(r.segments)==(0 if mode=='inside' else 1)
    np.testing.assert_allclose(r.field_at([0,0,3]),0)


@pytest.mark.parametrize('volume',[False,True])
@pytest.mark.parametrize('material',[rt.Material.pec(),rt.Material.metal(5.8e7)])
def test_metals_are_opaque(volume,material):
    s=rt.Scene()
    if volume: s.add_volume(rt.Box([0,0,1.1],[2,2,.1]),material)
    else: s.add_surface(wall(),material)
    r=rt.solve(s,plane(),cfg())
    assert len(r.segments)==2
    assert all(step[1]=='reflect' for t in r.terminals for step in t.history)
    np.testing.assert_allclose(r.field_at([0,0,2]),0)


def test_snell_brewster_and_power():
    theta=math.atan(2)
    d=[math.sin(theta),0,math.cos(theta)]
    p=[math.cos(theta),0,-math.sin(theta)]
    r=rt.solve(slab(thickness=.2),plane(direction=d,pol=p),cfg(max_interactions=1,min_power_fraction=0))
    branches=[s for s in r.segments if s.state.history]
    reflected=next(s for s in branches if s.state.history[-1][1]=='reflect')
    transmitted=next(s for s in branches if s.state.history[-1][1]=='transmit')
    assert np.linalg.norm(reflected.state.field)<1e-14
    assert transmitted.state.direction[0]==pytest.approx(math.sin(theta)/2)
    balance=r.power_balance
    assert sum(v for k,v in balance.items() if k!='launched')==pytest.approx(balance['launched'])


def test_total_internal_reflection():
    s=rt.Scene()
    s.add_volume(rt.Box([0,0,0],[5,5,1]),rt.Material.dielectric(4))
    d=[math.sin(.8),0,math.cos(.8)]
    r=rt.solve(s,plane(direction=d,pol=(0,1,0)),cfg(max_interactions=1))
    assert all(step[1]=='reflect' for t in r.terminals for step in t.history)
    assert r.absorbed_power==pytest.approx(0,abs=1e-15)


@pytest.mark.parametrize('epsilon',[4,4+.02j])
def test_coherent_slab_against_closed_form(epsilon):
    thickness=.013
    r=rt.solve(slab(epsilon,thickness),plane(),cfg(max_interactions=32,min_power_fraction=0))
    n=np.sqrt(complex(epsilon)); k=2*np.pi*10e9/rt.C0
    p=np.exp(1j*k*n*thickness)
    expected=(2/(1+n))*(2*n/(n+1))*p/(1-((n-1)/(n+1))**2*p*p)*np.exp(1j*k*(1-thickness))
    np.testing.assert_allclose(r.field_at([0,0,1]),[expected,0,0],atol=1e-12)


def test_focused_gaussian_and_semigroup():
    f=10e9; w=.1; zr=np.pi*w*w*f/rt.C0
    source=rt.GaussianBeam([0,0,1],[0,0,1],w,launch_distance=-1)
    r=rt.solve(rt.Scene(),source,cfg())
    e0=np.sqrt(4*rt.ETA0/(np.pi*w*w))
    for z in (.1,.6,1,1.7,3):
        s=z-1; k=2*np.pi*f/rt.C0
        expected=e0/(1+1j*s/zr)*np.exp(1j*k*s)
        np.testing.assert_allclose(r.field_at([0,0,z]),[expected,0,0],rtol=1e-12,atol=1e-10)
    from raytracing.beams import propagate
    a=r.segments[0].state
    one=propagate(a,2,rt.VACUUM,f)
    many=a
    for _ in range(100): many=propagate(many,.02,rt.VACUUM,f)
    np.testing.assert_allclose(one.field,many.field,rtol=1e-12)


def test_refraction_curvature_normal_incidence():
    r=rt.solve(slab(thickness=.2),rt.GaussianBeam([0,0,0],[0,0,1],.1),cfg(max_interactions=1))
    first=r.segments[0]
    from raytracing.beams import propagate
    hi=propagate(first.state,first.length,rt.VACUUM,10e9).curvature
    transmitted=next(s for s in r.segments if s.state.history and s.state.history[-1][1]=='transmit')
    np.testing.assert_allclose(transmitted.state.curvature,hi/2,atol=1e-12)


def test_result_roundtrip_and_plot(tmp_path):
    r=rt.solve(slab(),rt.GaussianBeam([0,0,0],[0,0,1],.1),cfg(max_interactions=4))
    path=tmp_path/'beam.json'
    rt.save_result(r,path)
    loaded=rt.load_result(path)
    np.testing.assert_allclose(loaded.field_at([.03,.01,1]),r.field_at([.03,.01,1]))
    assert loaded.power_balance==r.power_balance
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    viewer=rt.visualize(loaded,show=False,path=tmp_path/'view.png')
    viewer.close()
    assert (tmp_path/'view.png').stat().st_size>1000


def test_limits_are_diagnostic_and_power_accounted():
    for options,reason in [({'max_interactions':0},'interaction_limit'),
                           ({'max_branches':1},'branch_limit'),
                           ({'min_power_fraction':.2},'power_cutoff')]:
        r=rt.solve(slab(),plane(),cfg(**options))
        assert any(t.reason==reason for t in r.terminals)
        b=r.power_balance
        assert sum(v for k,v in b.items() if k!='launched')==pytest.approx(b['launched'])


def test_gaussian_beamlet_reconstruction():
    source=rt.GaussianBeam([0,0,0],[0,0,1],.1)
    direct=rt.solve(rt.Scene(),source,cfg())
    split=rt.solve(rt.Scene(),source.beamlets(samples=25),cfg(min_power_fraction=0))
    for point in ([0,0,.3],[.03,-.02,.7],[.1,0,1.3]):
        np.testing.assert_allclose(split.field_at(point),direct.field_at(point),rtol=2e-6,atol=1e-6)
    assert any('Coherent source' in w for w in split.warnings)


def test_aperture_source_and_antenna_adapter():
    field=lambda p: np.array([np.exp(-(p[0]**2+p[1]**2)/.1**2),0,0],complex)
    source=rt.ApertureSource([0,0,0],[0,0,1],[.2,.2],9,.05,field)
    r=rt.solve(rt.Scene(),source,cfg(min_power_fraction=0))
    np.testing.assert_allclose(r.field_at([0,0,0]),[1,0,0],atol=1e-8)
    antenna=rt.Antenna(rt.Isotropic())
    source=rt.ApertureSource.from_antenna(antenna,[0,0,-10],[0,0,0],[0,0,1],[.2,.2],5,.1,10e9)
    r=rt.solve(rt.Scene(),source,cfg())
    assert np.isfinite(r.field_at([0,0,.1])).all()
    with pytest.raises(ValueError): rt.solve(rt.Scene(),source,rt.SolverConfig(frequency_hz=20e9))


def test_vtk_pretrace_inspector_never_traces(tmp_path,monkeypatch):
    def forbidden(*args,**kwargs): raise AssertionError('Preview must not trace')
    monkeypatch.setattr(rt,'solve',forbidden)
    import raytracing.transport
    monkeypatch.setattr(raytracing.transport,'solve',forbidden)
    scene=slab()
    tx=rt.GaussianBeam([0,0,.4],[0,0,1],.05,launch_distance=-.4)
    rx=rt.GaussianBeam([.1,0,1],[0,0,-1],.06)
    view=rt.inspect_scene(scene,tx,rx,frequency_hz=77e9,show=False,path=tmp_path/'pretrace.png')
    actors=view.renderer.actors
    assert 'geometry 0' in actors and 'Tx 1 envelope 0' in actors and 'Rx 1 look envelope 0' in actors
    assert (tmp_path/'pretrace.png').stat().st_size>1000
    pixels=np.asarray(view.image,dtype=int)[100:-100,:,:3]
    assert np.count_nonzero(pixels[:,:,2]-pixels[:,:,0]>12)>500  # Material geometry, not only text.
    view.close()
    pattern=rt.AntennaPattern([0,0,1],rt.Antenna(rt.ShortDipole()))
    view=rt.inspect_scene(scene,receivers=pattern,show=False)
    assert 'Rx 1 look pattern' in view.renderer.actors
    view.close()
