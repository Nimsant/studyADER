import matplotlib.pyplot as plt
import pandas as pd

what2plot = "Linf"
xaxis = 'N'
xfilter = 'NBx'
filtervalues = range(1,5)

df = pd.read_csv("run.log",
    comment="#",
    sep = '\s+',
    header=0,
    #skip_blank_lines=True,
    )
print(df)

def palitra(i):
 return ['r','g','b','y'][i]
Nxset = [int(i) for i in list(set(df[xaxis]))]
x = [i for i in Nxset]
fig,ax = plt.subplots()
ax.set_title(f"{what2plot}")
for  iM, M in enumerate(filtervalues):
    y = [i**(-(M))*9e-1 for i in x]
    ax.plot(x,y,lw=30,color=palitra(iM),alpha=0.1)
    df[df[xfilter]==M].plot(x=xaxis,y=what2plot,
            ax=ax,
            label=f'{xfilter} = {M}',
            color=palitra(iM)
            
            )
plt.xscale('log')
plt.yscale('log')
plt.show()
