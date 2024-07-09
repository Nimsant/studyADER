import matplotlib.pyplot as plt

from PMpolynoms import polynomial_basis

import json
import pandas as pd
import numpy as np

import glob
from icecream import ic
droplist = glob.glob("drop_*")

pallete = plt.get_cmap('Spectral')
fig, ax = plt.subplots()

ax.set_facecolor((0.5, 0.5, 0.5))
Npoints = 400

def plotfile(ifile, filename):
    ax= plt.gca()
    ax.clear()
    print(filename)
    with open(filename) as f:
        s = f.readline()[1:]
        data = json.loads(s)
        df = pd.read_csv(filename, 
                sep='\s+', 
                comment='#', 
                names = [f"u{iq}{ib}" for iq in range(data['NQ']) for ib in range(data['NBx']) ]
                )
    bp = polynomial_basis(data['NBx']-1)['polynomials']

    n_points = max(2,Npoints//df.shape[0])
    for irow, row in df.iterrows():
        x = np.linspace( irow*data['dx'], (irow+1)*data['dx'] )
        xi = np.linspace(0,1)
        y = 0*xi
        for ib in range(data['NBx']):
          y += row[f'u0{ib}']*bp[ib](xi)
        ax.plot(x,y,
                lw = 1 if '-' not in filename else 3,
                color = pallete((ifile+0.5)/len(droplist)),
                label = '_'*(irow!=0)+f"{filename.split('_')[1]}")
    ax.set_title(f"{float(filename.split('_')[1].split('.')[0])*float(data['dt'])}")
    x = np.linspace( 0, data['N']*data['dx'],200 )

    ax.grid()
    ax.legend()
    fig.canvas.draw()

inum = 0
filename = list(sorted(droplist))[inum]
plotfile(inum, filename)

def nextfile(event):
    global inum
    print(event.key)
    if event.key=='right':
        inum = (inum+1)%len(droplist)
    if event.key=='left':
        inum = (inum-1+len(droplist))%len(droplist)
    filename = list(sorted(droplist))[inum]
    plotfile(inum, filename)


cid = fig.canvas.mpl_connect('key_press_event', nextfile)

plt.show()


