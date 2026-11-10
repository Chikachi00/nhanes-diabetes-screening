"""Create Figure 5 from aggregate exports only.
Usage: python figure5.py CODE_ROOT V11_RESULTS OUTPUT_DIRECTORY
"""
from pathlib import Path
import sys,json,os,tempfile
os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'nhanes-diabetes-figure5-cache'))
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from PIL import Image
C,A,O=map(Path,sys.argv[1:4]);O.mkdir(parents=True,exist_ok=True)
s=pd.read_csv(C/'analysis/reviewer_response_v9/results/S21_comparator_summary.csv').set_index('key')
curve=pd.read_csv(A/'Figure5_equation_budget_curve.csv')
bc=pd.read_csv(A/'Figure5_Bang_attainable_points.csv')
blue='#146B8C';orange='#C76529';grey='#4E5963';green='#34765A'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10.5,'axes.labelsize':10.5,'axes.titlesize':11,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
f=plt.figure(figsize=(8.4,6.8),layout=None)
gs=GridSpec(2,2,figure=f,width_ratios=[1.08,1],height_ratios=[1.08,1],left=.085,right=.97,bottom=.19,top=.91,wspace=.36,hspace=.55)
ax=f.add_subplot(gs[:,0])
ax.fill_between(curve.actual_referral_pct,curve.lower95,curve.upper95,color=blue,alpha=.12,lw=0)
ax.plot(curve.actual_referral_pct,curve.sensitivity_pct,color=blue,lw=1.7,label='Equation single cut-off')
b=bc[(bc.referral_pct>=7)&(bc.referral_pct<=65)]
ax.scatter(b.referral_pct,b.sensitivity_pct,marker='s',s=38,facecolors='white',edgecolors=grey,lw=1.3,label='Bang attainable cut-offs',zorder=3)
for r in b.itertuples():
 ax.annotate(f'≥{r.threshold:g}',(r.referral_pct,r.sensitivity_pct),xytext=(-2,-15),textcoords='offset points',fontsize=9,color=grey)
for n,col,mark,off in [('Bang',grey,'s',(-42,-1)),('R0',blue,'o',(5,-12)),('R1',orange,'^',(5,5))]:
 x=s.loc[f'{n}::overall::referral','estimate'];y=s.loc[f'{n}::overall::sensitivity','estimate']
 ax.scatter(x,y,s=52,c=col,marker=mark,zorder=5)
 ax.annotate(n,(x,y),xytext=off,textcoords='offset points',color=col,fontsize=10,fontweight='bold')
ax.axvline(s.loc['Bang::overall::referral','estimate'],color='#ADB5BC',lw=.8,ls=':',zorder=0)
ax.set(xlim=(8,61),ylim=(25,102),xlabel='Adults referred for testing (%)',ylabel='HbA1c-positive adults identified (%)')
ax.set_title('A  Case finding at different budgets',loc='left',pad=14,fontweight='bold')
ax.legend(loc='lower right',frameon=False,fontsize=8.4,handlelength=1.5)
ax.grid(axis='y',color='#E4E7EA',lw=.6)
# B uses separate scales for referrals and positives, both per 1,000 screened.
bgrid=gs[0,1].subgridspec(2,1,hspace=1.3)
rows=[]
for i,(measure,label,xmax,fmt,col) in enumerate([('referrals_per1000','Additional laboratory referrals',75,'.1f',orange),('found_per1000','Additional HbA1c positives',.9,'.2f',green)]):
 a=f.add_subplot(bgrid[i]);r=s.loc[f'R1_minus_R0::overall::{measure}'];est=r.estimate;lo=r.pct_lo;hi=r.pct_hi
 a.errorbar(est,0,xerr=[[est-lo],[hi-est]],fmt='o',color=col,elinewidth=1.7,capsize=4,markersize=6)
 a.set(xlim=(-.025*xmax,xmax),ylim=(-.6,.6),yticks=[])
 a.set_xlabel('Per 1,000 screened',fontsize=9.5,labelpad=1)
 a.set_title(label,loc='left',fontsize=10,pad=2)
 a.text(.98,.84,f'{est:{fmt}} ({lo:{fmt}} to {hi:{fmt}})',transform=a.transAxes,ha='right',va='center',fontsize=9,color=col)
 a.spines['left'].set_visible(False);a.axvline(0,color='#ADB5BC',ls=':',lw=.8)
 a.tick_params(axis='x',labelsize=9)
 rows.append(dict(panel='B',measure=measure,estimate=est,lower95=lo,upper95=hi,interval='percentile'))
box=gs[0,1].get_position(f)
f.text(box.x0,box.y1+.045,'B  R1 adds testing with low yield',fontsize=11,fontweight='bold')
# C positive values denote widening of the high-minus-low access sensitivity gap.
a=f.add_subplot(gs[1,1])
for y,n,col in [(1,'R0',blue),(0,'R1',orange)]:
 r=s.loc[f'{n}_minus_Bang::access_gap::sensitivity'];est=r.estimate;lo=r.normal_lo;hi=r.normal_hi
 a.errorbar(est,y,xerr=[[est-lo],[hi-est]],fmt='o',color=col,capsize=4,markersize=6)
 a.text(0,y+.29,f'{est:+.2f} ({lo:.2f} to {hi:.2f})',fontsize=8.9,ha='center',color=col)
 rows.append(dict(panel='C',strategy=n,measure='change_in_access_gap_pp',estimate=est,lower95=lo,upper95=hi,interval='normal'))
a.axvline(0,color=grey,ls=':',lw=1)
a.set(xlim=(-16,16),ylim=(-.6,1.7),yticks=[0,1],yticklabels=['R1 − Bang','R0 − Bang'],xlabel='Change in sensitivity gap (pp)')
a.tick_params(axis='y',length=0,labelsize=9)
a.set_title('C  Access gap change is uncertain',loc='left',fontweight='bold',pad=13)
a.spines['left'].set_visible(False)
a.text(0,-.53,'← Narrower                         Wider →',ha='center',fontsize=8.5,color=grey)
f.text(.085,.035,'Main self-report Bang implementation, 2013–2018 (n = 11,848).\nComposite-hypertension sensitivity analyses are reported in Supplement S25.',fontsize=9,color=grey)
f.savefig(O/'Figure 5.png',dpi=450,facecolor='white')
f.savefig(O/'Figure 5.pdf',facecolor='white')
f.savefig(O/'Figure 5.svg',facecolor='white')
im=Image.open(O/'Figure 5.png').convert('RGB');im.save(O/'Figure 5.tif',compression='tiff_lzw',dpi=(450,450))
pd.DataFrame(rows).to_csv(A/'Figure5_panels_BC.csv',index=False)
(A/'Figure5_render_metadata.json').write_text(json.dumps(dict(width=im.width,height=im.height,dpi=450,source='aggregate exported results',panels=['budget curve','nested R1-R0 increments','access gap change']),indent=2))
print(O/'Figure 5.tif')
