import matplotlib.pyplot as plt
import pandas as pd

what2plot = "Linf"
xaxis = 'N'
xfilter = 'NBx'
filtervalues = list(range(1,9))
Nrow = 4

df = pd.read_csv("run.log",
    comment="#",
    sep = '\s+',
    header=0,
    #skip_blank_lines=True,
    )
df['courant'] = df['dt']/df['dx']
print(df)

def palitra(i):
 return ['r','g','b','y','c','m','k','fuchsia'][i]

fig,axs = plt.subplots( len(filtervalues)//Nrow, Nrow)
#ax.set_title(f"{what2plot}")
for  iM, M in enumerate(filtervalues):
    dfhere = df[df[xfilter]==M]
    if dfhere.shape[0]>0:

        Nxset = [int(i) for i in list(set(dfhere[xaxis]))]
        x = [i for i in Nxset]
        ymin = dfhere[dfhere[xaxis] == min(x)].iloc[0]['Linf']

        ax = axs[iM//Nrow][iM%Nrow]
        ax.set_title(f"{xfilter} = {M}")
        ax.set_xscale('log')
        ax.set_yscale('log')
        y = [i**(-(M))*ymin/(min(x)**(-M)) for i in x]
        ax.plot(x,y,lw=30,color=palitra(iM),alpha=0.1)

        y = [i**(-(1))*ymin/(min(x)**(-1)) for i in x]
        ax.plot(x,y,lw=30,color='grey',alpha=0.1)

        courantset = [i for i in list(set(dfhere['courant']))]
        NBtset = [i for i in list(set(dfhere['NBt']))]

        for icourant, courant in enumerate(courantset):
          for iNBt, NBt in enumerate(NBtset):
            dfhere[dfhere['courant']==courant][dfhere['NBt']==NBt].plot(x=xaxis,y=what2plot,
                    ax=ax,
                    label=f'dx/dt = {courant}, NBt = {NBt}',
                    #color=palitra(icourant),
                    lw =1
                    )
ax.set_ylim(1e-18, 1e1)
plt.show()
