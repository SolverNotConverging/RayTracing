"""Render native point-ray caustic checks from SpreadingBenchmarks output."""
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['svg.hashsalt']='raytracing-caustics'
import matplotlib.pyplot as plt

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory',type=Path,nargs='?',default=Path('benchmarks/results/caustics'))
    directory=parser.parse_args().directory
    data=np.genfromtxt(directory/'caustics.csv',delimiter=',',names=True,dtype=None,encoding='utf-8')
    fig,axes=plt.subplots(3,2,figsize=(11,9),layout='constrained')
    for col,shape in enumerate(('cylinder','sphere')):
        for angle,color in zip((0,15,30),('tab:blue','tab:orange','tab:green')):
            rows=data[(data['shape']==shape)&(data['angle_deg']==angle)]
            x=rows['distance_m']; label=f'{angle} degrees'
            axes[0,col].plot(x,rows['analytical_factor'],color=color,label=label)
            axes[0,col].plot(x[::5],rows['numerical_factor'][::5],'.',color=color)
            axes[1,col].plot(x,-90*rows['analytical_mu'],color=color)
            axes[1,col].plot(x[::5],-90*rows['numerical_mu'][::5],'.',color=color)
            axes[2,col].semilogy(x,np.maximum(rows['relative_complex_error'],1e-16),color=color)
        axes[0,col].set(title=f'Concave {shape}: lines analytical, dots native',ylabel='Field spreading factor',yscale='log')
        axes[0,col].legend()
        axes[1,col].set(ylabel='Internal caustic phase (degrees)')
        axes[2,col].set(ylabel='Relative complex-field error',xlabel='Distance after reflection (m)')
        for ax in axes[:,col]: ax.grid(alpha=.25)
    for ext in ('png','svg'):
        path=directory/f'caustic_phase_validation.{ext}'
        fig.savefig(path,dpi=200,metadata={'Date':None} if ext=='svg' else None)
        if ext=='svg':path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines())+'\n')
    plt.close(fig)

if __name__=='__main__': main()
