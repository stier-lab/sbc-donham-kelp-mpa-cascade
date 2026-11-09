from pathlib import Path
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
out=Path('outputs/model_eligibility_20261101');a=pd.read_csv('outputs/researcher_scope_20261101/annual_input.csv');l=pd.read_csv(out/'model_eligibility.csv');l=l[l.prior_change_fit].sort_values(['Taxa','Resp','Source','MPA']);plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
colors={'Detected baseline trend':'#ba3c31','Nonlinear or level-shift flag':'#b36b00','Unresolved: sparse baseline':'#777777','Unresolved: constant baseline':'#777777','No detected departure; assess precision':'#247ba0'}
short={'M. pyrifera':'Kelp','M. franciscanus':'Red urchin','S. purpuratus':'Purple urchin','S. pulcher':'Sheephead','P. interruptus':'Lobster'}
with PdfPages(out/'baseline_atlas.pdf') as pdf:
 for page,start in enumerate(range(0,len(l),9),1):
  fig,axes=plt.subplots(3,3,figsize=(13,10))
  for ax,(_,r) in zip(axes.flat,l.iloc[start:start+9].iterrows()):
   d=a[(a.series_id==r.series_id)&(a.t<0)].sort_values('year');c=colors[r.status];ax.plot(d.year,d.lnDiff,'o-',color=c,linewidth=.8,markersize=3)
   if np.isfinite(r.slope):ax.plot(d.year,d.lnDiff.mean()+r.slope*(d.year-d.year.mean()),color='black',linestyle='--',linewidth=1)
   ax.text(.02,.97,f'{r.MPA}\n{r.Source} · {short[r.Taxa]} {r.Resp}',transform=ax.transAxes,va='top',fontsize=9,bbox=dict(facecolor='white',alpha=.8,edgecolor='none'));ax.margins(y=.35);ax.tick_params(axis='x',rotation=30);ax.grid(alpha=.15)
  for ax in axes.flat[len(l.iloc[start:start+9]):]:ax.set_visible(False)
  fig.supxlabel('Calendar year');fig.supylabel('Log protected/reference ratio');fig.tight_layout(rect=(.02,.02,1,1));pdf.savefig(fig);fig.savefig(out/f'baseline_page_{page:02d}.png',dpi=110);plt.close(fig)
p=pd.read_csv(out/'pooled.csv');p=p[(p['view']=='All')&(p.metric=='Change')&(p.Taxa=='M. pyrifera')];order=['All_before_after','No_detected_linear_trend','Adequate_no_detected_departure','Tolerance_50','Tolerance_100'];labels=['All before–after comparisons','No detected linear trend','Adequate baseline; no departure detected','1.5-fold stability tolerance','2-fold stability tolerance'];fig,ax=plt.subplots(figsize=(9,4.8))
for i,s in enumerate(order):
 z=p[p.scenario==s]
 if len(z):
  r=z.iloc[0];ax.plot(r.estimate,i,'o',color='#247ba0',markerfacecolor='#247ba0' if np.isfinite(r.low) else 'white',markersize=7)
  if np.isfinite(r.low):ax.plot([r.low,r.high],[i,i],color='#247ba0',linewidth=2)
ax.axvline(0,color='#888',linestyle='--');ax.set_yticks(range(len(order)),labels);ax.invert_yaxis();ax.set_xlabel('Kelp before–after change (log-ratio units)');ax.grid(axis='x',alpha=.15);fig.tight_layout();fig.savefig(out/'kelp_model_eligibility.png',dpi=180);fig.savefig(out/'kelp_model_eligibility.pdf');plt.close(fig)
