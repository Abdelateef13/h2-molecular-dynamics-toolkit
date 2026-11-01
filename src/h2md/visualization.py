"""3D hydrogen dimers and periodic simulation cell visualization."""
from pathlib import Path

import numpy as np

from .atomistic import AtomicH2Trajectory


def display_positions(tr: AtomicH2Trajectory, frame: int) -> np.ndarray:
    """Wrap the molecular COM, NOT each H separately: preserve visible H--H bonds."""
    com = tr.com[frame]
    return (com % tr.box_length)[:, None, :] + (tr.atoms[frame] - com[:,None,:])


def write_extended_xyz(tr: AtomicH2Trajectory, path, stride: int = 1):
    """Multi-frame extended XYZ; usable in OVITO (H atoms + unit cell)."""
    if stride < 1: raise ValueError('stride >= 1')
    path = Path(path); path.parent.mkdir(exist_ok=True, parents=True)
    L=tr.box_length
    lattice=f'{L:.8f} 0 0 0 {L:.8f} 0 0 0 {L:.8f}'
    with path.open('w',encoding='utf-8') as f:
        for frame in range(0,tr.n_frames,stride):
            x = display_positions(tr, frame)
            f.write(f'{2*tr.n_molecules}\n')
            f.write(f'Lattice="{lattice}" Properties=species:S:1:pos:R:3:molecule:I:1 '
                    f'pbc="T T T" frame={frame} model="illustrative_reduced_units"\n')
            for i in range(tr.n_molecules):
                for j in range(2):
                    a,b,c=x[i,j]
                    f.write(f'H {a:.9f} {b:.9f} {c:.9f} {i+1}\n')
    return path


def plot_configuration(tr: AtomicH2Trajectory, out, frame: int=-1):
    """PNG 3D illustration of *actual* atomic trajectory snapshot."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    from itertools import product, combinations
    x=display_positions(tr,frame)
    L=tr.box_length
    fig=plt.figure(figsize=(9,8),dpi=150)
    ax=fig.add_subplot(111, projection='3d')
    for pair in x:
        ax.plot(pair[:,0],pair[:,1],pair[:,2],lw=1.7,color='#52718f',alpha=.8)
    ax.scatter(x[:,:,0].ravel(),x[:,:,1].ravel(),x[:,:,2].ravel(),
               s=70, color='#78c2ef',edgecolor='#1b4368',linewidth=.6,depthshade=True)
    for s,e in combinations(np.array(list(product([0,L],repeat=3))),2):
        if np.sum(np.abs(s-e)>1e-12)==1:
            ax.plot3D(*zip(s,e),color='#a4b8c8',alpha=.5,lw=1)
    m=.3
    ax.set(xlim=(-m,L+m),ylim=(-m,L+m),zlim=(-m,L+m),
           xlabel='x [σ]',ylabel='y [σ]',zlabel='z [σ]')
    ax.view_init(elev=23,azim=35)
    ax.set_title(f'Atomistic H₂ dimers | {tr.n_molecules} molecules, {tr.n_molecules*2} H atoms\n'
                 'Demonstration potential — reduced units',pad=17,fontsize=14)
    ax.set_box_aspect((1,1,1))
    fig.tight_layout()
    p=Path(out);p.parent.mkdir(exist_ok=True,parents=True)
    fig.savefig(p,bbox_inches='tight',facecolor='white')
    plt.close(fig)
    return p


def write_interactive_html(tr: AtomicH2Trajectory, path, max_frames: int=75):
    """Offline Plotly 3D rotating/zoomable animated atom/bond view."""
    import plotly.graph_objects as go
    from itertools import product, combinations
    L=tr.box_length
    ids=np.unique(np.linspace(0,tr.n_frames-1,min(tr.n_frames,max_frames)).astype(int))

    def traces(frame):
        x=display_positions(tr,frame)
        bonds=[]
        for p in x:
            bonds.extend([p[0],p[1],[np.nan,np.nan,np.nan]])
        b=np.array(bonds)
        return [
            go.Scatter3d(x=x[:,:,0].ravel(),y=x[:,:,1].ravel(),z=x[:,:,2].ravel(),
                mode='markers',name='H atoms',marker=dict(size=6,color='#64b5f6',
                line=dict(width=1,color='#154e75')),
                customdata=np.repeat(np.arange(1,tr.n_molecules+1),2),
                hovertemplate='H atom · molecule %{customdata}<br>x=%{x:.2f}, y=%{y:.2f}, z=%{z:.2f}<extra></extra>'),
            go.Scatter3d(x=b[:,0],y=b[:,1],z=b[:,2],mode='lines',name='H–H bonds',
                line=dict(width=6,color='#8c9eae'),hoverinfo='skip'),
        ]

    frame_objs=[go.Frame(name=str(k),data=traces(k),traces=[0,1]) for k in ids]
    borders=[]
    for s,e in combinations(np.array(list(product([0,L],repeat=3))),2):
        if np.sum(np.abs(s-e)>1e-12)==1:
            borders.extend([s,e,[np.nan,np.nan,np.nan]])
    b=np.array(borders)
    cube=go.Scatter3d(x=b[:,0],y=b[:,1],z=b[:,2],mode='lines',name='Periodic box',
                       line=dict(color='#64748b',width=3),hoverinfo='skip')
    fig=go.Figure(data=traces(int(ids[0]))+[cube],frames=frame_objs)
    fig.update_layout(title=f'H₂ atomistic molecular dynamics · {tr.n_molecules} molecules / {2*tr.n_molecules} H atoms'
                            '<br><sup>Two H atoms per molecule · illustrative reduced-unit potential (not quantitative hydrogen)</sup>',
        scene=dict(xaxis=dict(title='x [σ]',range=[-.3,L+.3]),
                   yaxis=dict(title='y [σ]',range=[-.3,L+.3]),
                   zaxis=dict(title='z [σ]',range=[-.3,L+.3]),aspectmode='cube',bgcolor='white'),
        paper_bgcolor='white',height=720,margin=dict(l=5,r=5,t=80,b=25),
        updatemenus=[dict(type='buttons',showactive=False,x=.02,y=0,
            buttons=[dict(label='▶ Play',method='animate',args=[None,dict(frame=dict(duration=95,redraw=True),fromcurrent=True,transition=dict(duration=0))]),
                     dict(label='⏸ Pause',method='animate',args=[[None],dict(frame=dict(duration=0,redraw=False),mode='immediate')])])],
        sliders=[dict(currentvalue=dict(prefix='Frame: '),pad=dict(t=30),
            steps=[dict(label=str(k),method='animate',args=[[str(k)],dict(mode='immediate',frame=dict(duration=0,redraw=True),transition=dict(duration=0))]) for k in ids])],
        legend=dict(orientation='h',y=1.03,x=0))
    path=Path(path);path.parent.mkdir(exist_ok=True,parents=True)
    fig.write_html(str(path),include_plotlyjs=True,full_html=True)
    return path


def plot_bond_diagnostics(tr: AtomicH2Trajectory, out):
    """Bond-length statistics from explicit H atom trajectories."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    lengths=tr.bond_lengths()
    fig,ax=plt.subplots(1,2,figsize=(10,3.7),dpi=160)
    time=np.arange(tr.n_frames)*tr.frame_dt
    ax[0].plot(time,lengths.mean(axis=1),lw=1.8,color='#1d6384')
    ax[0].fill_between(time,lengths.mean(axis=1)-lengths.std(axis=1),
                        lengths.mean(axis=1)+lengths.std(axis=1),
                        color='#9fcbe4',alpha=.4,label='±1 SD among molecules')
    ax[0].set(xlabel=f'Time [{tr.time_unit}]',ylabel=f'Average H–H length [{tr.length_unit}]',
              title='Bond vibration over time');ax[0].legend(fontsize=8)
    ax[1].hist(lengths.ravel(),bins=30,color='#83b8d5',edgecolor='white')
    ax[1].set(xlabel=f'H–H bond length [{tr.length_unit}]',ylabel='Count',
              title='H–H bond-length distribution')
    fig.suptitle('Explicit two-atom H₂ dynamics · educational potential',fontsize=12)
    fig.tight_layout()
    p=Path(out);p.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(p,bbox_inches='tight');plt.close(fig)
    return p
