"""Focused Gaussian beam through a lossy dielectric slab."""
from pathlib import Path
import argparse
import numpy as np
import raytracing as rt


def make_scene():
    scene=rt.Scene()
    scene.add_volume(rt.Box([0,0,.51],[.5,.5,.01]),rt.Material.dielectric(4+.02j,name='lossy glass'))
    scene.add_termination(rt.Box([0,0,.8],[2,2,1.8]),volume=True,on='exit')
    beam=rt.GaussianBeam([0,0,.8],[0,0,1],.04,launch_distance=-.8)
    config=rt.SolverConfig(frequency_hz=77e9,max_interactions=20,max_distance=8)
    receiver=rt.GaussianBeam([0,0,1],[0,0,-1],.04)
    return scene,beam,receiver,config


def main(show=True):
    output=Path('results/example1')
    output.mkdir(parents=True,exist_ok=True)
    scene,beam,receiver,config=make_scene()
    viewer=rt.inspect_scene(scene,beam,receiver,frequency_hz=config.frequency_hz,
                            path=output/'pretrace.png',show=show)
    viewer.close()
    result=rt.solve(scene,beam,config)
    rt.save_result(result,output/'simulation.json')
    viewer=rt.visualize(result,transmitters=beam,receivers=receiver,path=output/'paths.png',show=show)
    viewer.close()
    rt.plot_field(result,[[x,0,1] for x in np.linspace(-.15,.15,201)],path=output/'profile.png',show=show)
    print(result.power_balance)
    print('Warnings:',result.warnings)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--no-show',action='store_true')
    main(not parser.parse_args().no_show)
