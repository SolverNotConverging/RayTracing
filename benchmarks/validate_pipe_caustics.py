"""Replay frozen independent analytical circular-pipe solutions using native code.

Run with an installed raytracing wheel: python benchmarks/validate_pipe_caustics.py
"""
from pathlib import Path
import csv
import json
import numpy as np
from raytracing import _core as rt

ROOT=Path(__file__).resolve().parent

def main():
    fixture=json.loads((ROOT/'fixtures/circular_pipe.json').read_text())
    scene=rt.Scene()
    scene.add(rt.Cylinder([0,0,.201],[0,0,1],.0415,.201,False),rt.Material.pec())
    scene.add(rt.Disk([0,0,0],[0,0,1],.0415),rt.Material.pec())
    options=rt.SpreadingOptions(); options.angular_step=1e-5
    rows=[]
    for ray in fixture['rays']:
        source=ray['source']; direction=ray['direction']; sequence=ray['sequence']
        path=rt.refine_path(source,direction,source,scene.surfaces,sequence)
        spread=rt.calculate_spreading(source,direction,source,scene.surfaces,sequence,options)
        assert str(path.status).endswith('CONVERGED')
        assert str(spread.status).endswith('VALID')
        assert spread.caustic_count==ray['mu']
        relative=abs(spread.area_per_solid_angle/ray['area']-1)
        length_error=abs(path.geometry.receiver.path_distance-ray['length'])
        assert relative<1e-6 and length_error<1e-10
        rows.append(dict(offset_m=source[0],family=ray['family'],order=ray['order'],
                         length_m=ray['length'],analytical_mu=ray['mu'],native_mu=spread.caustic_count,
                         relative_area_error=relative,length_error_m=length_error))
    out=ROOT/'results/caustics';out.mkdir(parents=True,exist_ok=True)
    with (out/'pipe_caustics.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    summary=dict(tested_rays=len(rows),families=sorted({r['family'] for r in rows}),
                 caustic_mismatches=0,max_relative_area_error=max(r['relative_area_error'] for r in rows),
                 max_length_error_m=max(r['length_error_m'] for r in rows))
    (out/'pipe_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    matplotlib.rcParams['svg.hashsalt']='raytracing-caustics'
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,1,figsize=(9,7),layout='constrained')
    for family in summary['families']:
        group=[r for r in rows if r['family']==family]
        axes[0].scatter([r['analytical_mu'] for r in group],[r['native_mu'] for r in group],s=18,label=family)
        axes[1].scatter([r['order'] for r in group],[max(r['relative_area_error'],1e-16) for r in group],s=18)
    axes[0].set(xlabel='Analytical conjugate-point count',ylabel='Native conjugate-point count',title='Circular pipe: 634 independent analytical paths')
    axes[0].legend(ncol=3,fontsize=8)
    axes[1].set(xlabel='Wall reflections',ylabel='Relative spreading-area error',yscale='log')
    for ax in axes:ax.grid(alpha=.2)
    for ext in ('png','svg'):
        path=out/f'pipe_caustics.{ext}'
        fig.savefig(path,dpi=200,metadata={'Date':None} if ext=='svg' else None)
        if ext=='svg':path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines())+'\n')
    plt.close(fig)
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
