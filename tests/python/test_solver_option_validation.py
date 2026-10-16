"""Public-entry validation must preserve the pre-cleanup numerical guards."""
import numpy as np
import pytest
import raytracing as rt
from raytracing import beam


@pytest.mark.parametrize("group,attribute,value", [
    (None,"delay_tolerance",np.inf),
    ("refinement","max_iterations",-1),
    ("refinement","max_backtracks",-1),
    ("refinement","max_direction_step",0),
    ("refinement","max_direction_step",np.inf),
    ("refinement","t_min",-1),
    ("refinement","t_min",np.nan),
    ("spreading","max_step_halvings",-1),
    ("spreading","relative_derivative_tolerance",0),
    ("spreading","relative_derivative_tolerance",np.nan),
    ("spreading","absolute_derivative_tolerance",np.inf),
    ("spreading","minimum_singular_value",0),
    ("spreading","minimum_singular_value_ratio",0),
    ("spreading","minimum_singular_value_ratio",1),
    ("deduplication","angle_tolerance",np.pi+.01),
])
def test_point_solver_rejects_invalid_numerical_controls(group,attribute,value):
    config=rt.SolverConfig(ray_count=8,max_reflections=0,max_distance=2)
    setattr(config if group is None else getattr(config,group),attribute,value)
    with pytest.raises(ValueError):
        rt.solve(rt.Scene(),rt.Transmitter([0,0,0]),rt.Receiver([1,0,0]),config)


@pytest.mark.parametrize("limit", [-.01,.010001,np.inf,np.nan])
def test_beam_solver_enforces_the_supported_weak_loss_range(limit):
    config=beam.SolverConfig(frequency_hz=10e9,max_distance=1,max_loss_ratio=limit)
    with pytest.raises(ValueError):
        beam.solve(beam.Scene(),beam.PlaneWave(),config)


def test_zero_iteration_limits_and_zero_loss_ratio_remain_valid():
    config=rt.SolverConfig(ray_count=8,max_reflections=0,max_distance=2)
    config.refinement.max_iterations=0
    config.refinement.max_backtracks=0
    config.spreading.max_step_halvings=0
    result=rt.solve(rt.Scene(),rt.Transmitter([0,0,0]),rt.Receiver([1,0,0]),config)
    assert len(result.impulse_response.taps)==1
    result=beam.solve(beam.Scene(),beam.PlaneWave(),beam.SolverConfig(
        frequency_hz=10e9,max_distance=1,max_loss_ratio=0))
    assert np.isfinite(result.field_at([0,0,.5])).all()
