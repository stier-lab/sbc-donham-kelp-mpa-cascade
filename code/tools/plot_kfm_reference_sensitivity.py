from pathlib import Path
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[2];O=R/'outputs/kfm_reference_sensitivity_20261103'
plt.rcParams.update({'font.size':11,'font.family':'DejaVu Sans','axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
s=pd.read_csv(O/'sheephead_trends.csv');p=pd.read_csv(O/'pooled_results.csv')
fig,axes=plt.subplots(1,2,figsize=(12,4.1),layout='constrained')
sets=[('Harris Point SMR',[('main','Current · 19 years'),('matched_main_without_wyckoff','Current · same 15 years'),('without_wyckoff','Without Wyckoff · 15 years'),('without_miracle','Without Miracle · 19 years')]),('Santa Barbara Island SMR',[('main','Current'),('without_webster','Without Webster'),('matched_main_without_flagged','Harris omitted · Webster kept'),('without_flagged','Harris and Webster omitted')])]
for ax,(mpa,items),letter in zip(axes,sets,['A','B']):
 for i,(scenario,label) in enumerate(items):
  z=s[s.MPA.eq(mpa)&s.scenario.eq(scenario)].iloc[0];v,lo,hi=z[['estimate_pct','low_pct','high_pct']];color='#53656a' if scenario=='main' else '#267762';ax.errorbar(v,-i,xerr=[[v-lo],[hi-v]],fmt='o',color=color,capsize=3,markersize=6)
 ax.set_yticks(-np.arange(len(items)),[x[1] for x in items]);ax.axvline(0,color='.7',lw=1);ax.set_xlim(-32,35);ax.set_xlabel('Change in density ratio (% per year)');ax.text(0,1.05,letter+'  '+mpa.replace(' SMR',''),transform=ax.transAxes,fontweight='bold');ax.set_ylim(-3.6,.5);ax.grid(axis='x',alpha=.12)
fig.savefig(O/'sheephead_reference_trends.png');fig.savefig(O/'sheephead_reference_trends.pdf');plt.close(fig)
# Three cases: full main, original references with Harris omitted, then Webster omitted too.
sp={'M. franciscanus':'Red urchin','S. purpuratus':'Purple urchin','P. interruptus':'Lobster','S. pulcher':'Sheephead','M. pyrifera':'Kelp'}
q=p[p.Source.eq('All programs')&p.model.eq('post_trend')&p.term.eq('recovery_slope')]
keys=list(q[q.scenario.eq('main')][['Taxa','Resp']].itertuples(index=False,name=None));keys=sorted(keys,key=lambda k:(k[1],sp[k[0]]))
fig,ax=plt.subplots(figsize=(9,6),layout='constrained')
for j,(name,color,dy) in enumerate([('main','#53656a',.19),('matched_main_without_flagged','#be8b37',0),('without_flagged','#267762',-.19)]):
 for i,(taxon,resp) in enumerate(keys):
  z=q[q.scenario.eq(name)&q.Taxa.eq(taxon)&q.Resp.eq(resp)].iloc[0];v,lo,hi=100*np.expm1(z[['estimate','low','high']].astype(float));ax.errorbar(v,-i+dy,xerr=[[v-lo],[hi-v]],fmt=['o','s','^'][j],color=color,capsize=2,markersize=5)
ax.set_yticks(-np.arange(len(keys)),[sp[t]+' · '+('biomass' if r=='Bio' else 'density') for t,r in keys]);ax.axvline(0,color='.7',lw=1);ax.grid(axis='x',alpha=.12);ax.set_xlabel('Change in protected/reference ratio (% per year)');fig.savefig(O/'pooled_reference_trends.png');fig.savefig(O/'pooled_reference_trends.pdf');plt.close(fig)
