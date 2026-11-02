"""Generate a figure of the *site-site* H2 LJ potential from paper parameters.

Not a simulated trajectory or molecule-molecule potential; four H-site pairs
contribute to the full potential between two rigid H2 molecules.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from build_h2_paper import SIGMA_A, EPS_KCAL_MOL, CUTOFF_A


def lj(r):
    r=np.asarray(r,float)
    return 4*EPS_KCAL_MOL*((SIGMA_A/r)**12-(SIGMA_A/r)**6)


def main():
    r=np.linspace(2.6,CUTOFF_A,650)
    fig,ax=plt.subplots(figsize=(8.2,4.4))
    ax.plot(r,lj(r),linewidth=2.5,color='#2563eb')
    rmin=2**(1/6)*SIGMA_A
    ax.scatter([rmin],[-EPS_KCAL_MOL],color='#ea580c',s=42,zorder=3)
    ax.annotate(f'r_min = {rmin:.3f} Å\nU = -{EPS_KCAL_MOL:.4f} kcal/mol',
        (rmin,-EPS_KCAL_MOL),xytext=(4.4,0.040),
        arrowprops={'arrowstyle':'->','color':'#64748b'},fontsize=9)
    ax.axhline(0,color='#94a3b8',linestyle='--',linewidth=.8)
    ax.axvline(CUTOFF_A,color='#64748b',linestyle=':',linewidth=1.0)
    ax.set_xlim(2.6,12.9);ax.set_ylim(-.025,.065)
    ax.set_xlabel('Separation between H interaction sites (Å)')
    ax.set_ylabel('LJ pair potential (kcal/mol)')
    ax.set_title('Two-site H₂ force field — published H-site Lennard–Jones parameters',pad=13)
    ax.grid(alpha=.18)
    fig.text(.11,.005,'Adam et al. (2024), IJHE 78, Table 1  |  Site–site potential, not molecular COM potential',fontsize=8,color='#64748b')
    fig.tight_layout(rect=(0,.04,1,1))
    d=Path(__file__).parents[2]/'assets'
    d.mkdir(parents=True,exist_ok=True)
    out=d/'h2_paper_lj_potential.png'
    fig.savefig(out,dpi=165,facecolor='white')
    plt.close(fig)
    print('Saved',out)

if __name__=='__main__':main()
