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
for ifile, filename in enumerate(sorted(droplist)):
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
    ax.set_title(f"{filename}")
    x = np.linspace( 0, data['N']*data['dx'],200 )

ax.grid()
ax.legend()
plt.show()



