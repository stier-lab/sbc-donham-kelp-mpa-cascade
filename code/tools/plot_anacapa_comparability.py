from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path('outputs/anacapa_comparability_20261101');a=pd.read_csv(out/'diagnostic_annual.csv');a=a[(a.scenario=='All mapped stations')&(a.variable=='biomass')]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(1,2,figsize=(11,4.2),sharey=True)
for i,(role,label) in enumerate([('mpa','Protected'),('reference','Reference')]):
 for source,color,marker in [('KFM','#247ba0','o'),('PISCO','#d46b27','^')]:
  d=a[a.Source==source].sort_values('year');ax[i].plot(d.year,d[role],color=color,marker=marker,linewidth=2,markersize=7)
 ax[i].text(.04,.96,chr(65+i),transform=ax[i].transAxes,va='top',fontweight='bold');ax[i].text(.5,.96,label,transform=ax[i].transAxes,va='top',ha='center');ax[i].set_xticks(range(2007,2013));ax[i].tick_params(axis='x',rotation=40);ax[i].set_xlabel('Year');ax[i].set_ylim(-8,410);ax[i].axhline(0,color='#aaa',linewidth=.6);ax[i].grid(axis='y',alpha=.15)
ax[0].set_ylabel('Kelp dry biomass (g m$^{-2}$)');fig.tight_layout();fig.savefig(out/'protected_reference.png',dpi=190);fig.savefig(out/'protected_reference.pdf');plt.close(fig)
# Date-matched year illustration: zones retained as essential identifiers.
p=pd.read_csv(out/'pisco_transects.csv');p=p[(p.year==2009)&p.site.str.contains('LIGHTHOUSE')];fig,ax=plt.subplots(figsize=(7,4.4))
for i,zone in enumerate(['INNER','OUTER']):
 z=p[p.zone==zone].sort_values(['site','transect']);xx=i+pd.Series(range(len(z))).to_numpy()*.045-.11;ax.scatter(xx,z.biomass,c='#d46b27',marker='^',s=55)
ax.scatter([2],[0],color='#247ba0',s=65);ax.set_xticks([0,1,2],['PISCO inner\n4.0–5.2 m','PISCO outer\n7.6–9.8 m','KFM Lighthouse\n~7.6 m']);ax.set_ylabel('Kelp dry biomass (g m$^{-2}$)');ax.axhline(0,color='#aaa',linewidth=.7);ax.grid(axis='y',alpha=.15);fig.tight_layout();fig.savefig(out/'lighthouse_2009.png',dpi=190);fig.savefig(out/'lighthouse_2009.pdf');plt.close(fig)
