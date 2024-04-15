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

Npoints = 400
for ifile, filename in enumerate(droplist):
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
                lw = 1,
                color = pallete((ifile+0.5)/len(droplist)),
                label = '_'*(irow!=0)+f"{filename.split('_')[1]}")
    ax.set_title(f"{filename}")
    x = np.linspace( 0, data['N']*data['dx'],200 )
    phi = x - 0.5*data['T']
    step  = np.heaviside(phi-data['N']/2,.5) 
    step += np.heaviside(-phi,.5) 
    step = 2*step -1
    #ax.plot(x,step,
    #ax.plot(x,np.sin(2*np.pi*(x-0.5*data['T'])/data["N"]),
    ax.plot(x,phi,
            lw=10, alpha=.2,
            label=f"theor_{filename.split('_')[1]}",
            color = pallete((ifile+0.5)/len(droplist)))

ax.grid()
ax.legend()
plt.show()



