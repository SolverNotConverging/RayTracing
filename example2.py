"""Concave opaque PEC mirror: finite focal field and astigmatic beam transport."""
from pathlib import Path
import argparse
import numpy as np
import raytracing as rt


def make_scene():
    scene=rt.Scene()
    # A spherical surface is an opaque shell. Launching inside the shell exposes
    # its concave side without placing the source inside an opaque volume.
    scene.add_surface(rt.Sphere([0,0,0],2),rt.Material.pec())
    beam=rt.GaussianBeam([0,0,0],[0,0,1],.1)
    receiver=rt.GaussianBeam([0,0,1],[0,0,1],.03)
    config=rt.SolverConfig(frequency_hz=77e9,max_interactions=1,max_distance=3.9)
    return scene,beam,receiver,config


def main(show=True):
    output=Path('results/example2')
    output.mkdir(parents=True,exist_ok=True)
    scene,beam,receiver,config=make_scene()
    viewer=rt.inspect_scene(scene,beam,receiver,frequency_hz=config.frequency_hz,path=output/'pretrace.png',show=show)
    viewer.close()
    result=rt.solve(scene,beam,config)
    rt.save_result(result,output/'simulation.json')
    viewer=rt.visualize(result,transmitters=beam,receivers=receiver,path=output/'paths.png',show=show)
    viewer.close()
    rt.plot_field(result,[[x,0,1] for x in np.linspace(-.1,.1,201)],path=output/'focus.png',show=show)
    print(result.power_balance)
    print('Field includes coherent incident and reflected beams.')
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--no-show',action='store_true')
    main(not parser.parse_args().no_show)
