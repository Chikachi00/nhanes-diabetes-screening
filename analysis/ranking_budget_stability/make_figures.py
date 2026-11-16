"""Scientific figures from aggregate results only.
python make_figures.py RESULTS FIGURES
"""
from pathlib import Path
import os,sys,tempfile,json
os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'nhanes-supplementary-check-mpl'))
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R,O=map(Path,sys.argv[1:3]);O.mkdir(parents=True,exist_ok=True)
S=pd.read_csv(R/'summary_all.csv').set_index('key')
SI=pd.read_csv(R/'budget_simultaneous_intervals.csv').set_index('key')
blue='#156982';orange='#B35F20';grey='#4B535B'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'axes.labelsize':12,
 'axes.titlesize':12,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
scenarios=['common','strict_BP','cycle_H','cycle_I','cycle_J','without_H','without_I','without_J']
labels=['Common sample','Strict BP sample','2013–2014','2015–2016','2017–2018','Without 2013–2014','Without 2015–2016','Without 2017–2018']

def interval(ax,r,y,color=blue):
    ax.plot([r.normal_lower,r.normal_upper],[y,y],color=color,lw=1.8)
    ax.scatter(r.estimate,y,color=color,s=28,zorder=3)
    ax.plot([r.percentile_lower,r.percentile_upper],[y-.12,y-.12],color=orange,lw=1.4)
def save(fig,name):
    fig.savefig(O/(name+'.png'),dpi=450,facecolor='white')
    fig.savefig(O/(name+'.pdf'),facecolor='white')
    fig.savefig(O/(name+'.svg'),facecolor='white')
    plt.close(fig)

f,axes=plt.subplots(1,2,figsize=(9.6,5.7),sharey=True)
f.subplots_adjust(left=.24,right=.975,top=.85,bottom=.30,wspace=.23)
y=np.arange(len(scenarios))[::-1]
for i,(ax,tail,title) in enumerate(zip(axes,['auc::equation_minus_composite','matched5::difference_Bang_composite::sensitivity_pct'],['A  Discrimination','B  Detection at the ≥5 budget'])):
    for yy,sn in zip(y,scenarios): interval(ax,S.loc[sn+'::'+tail],yy)
    ax.axvline(0,color=grey,lw=.9,ls=':');ax.grid(axis='x',color='#E2E5E8',lw=.6)
    ax.set_title(title,loc='left',pad=14,fontweight='bold')
    ax.set_ylim(-.5,len(y)-.3);ax.tick_params(axis='y',length=0)
axes[0].set(yticks=y,yticklabels=labels,xlabel='Equation − expanded Bang AUC',xlim=(-.012,.106))
axes[1].set(xlabel='Sensitivity difference (pp)',xlim=(-13,28))
axes[0].set_xticks([0,.025,.05,.075,.10]);axes[1].set_xticks([-10,0,10,20])
f.text(.24,.045,'Blue: point estimate and normal 95% interval. Orange: percentile 95% interval.\nPaired 1,000-replicate analysis; cycle intervals are pointwise.\nExpanded hypertension definition; activity mapping differs from original Bang.',fontsize=12,color=grey)
save(f,'Figure_A_auc_and_cycle_stability')

f,axes=plt.subplots(1,2,figsize=(9.6,6.8))
f.subplots_adjust(left=.085,right=.98,top=.88,bottom=.22,wspace=.43)
ax=axes[0]
g=[]
for n in range(10,61):
    prefix=f'common::grid::g{n:02d}'
    ref=S.loc[prefix+'::referral_pct'];se=S.loc[prefix+'::sensitivity_pct']
    g.append([ref.estimate,se.estimate,se.percentile_lower,se.percentile_upper])
g=np.array(g)
ax.fill_between(g[:,0],g[:,2],g[:,3],color=blue,alpha=.12,lw=0)
ax.plot(g[:,0],g[:,1],color=blue,lw=1.8,label='Equation (10–60% grid)')
rows=[]
for t in range(-1,11):
    pre=f'common::budget::Bang_composite::k{t}'
    x=S.loc[pre+'::Bang::referral_pct','estimate'];v=S.loc[pre+'::Bang::sensitivity_pct','estimate']
    rows.append((t,x,v))
    ax.scatter(x,v,marker='s',s=30,facecolors='white',edgecolors=grey,zorder=3)
    if t in [4,5,6,7]:ax.annotate(f'≥{t}',(x,v),xytext=(-10,-17),textcoords='offset points',fontsize=12)
ax.scatter([],[],marker='s',facecolors='white',edgecolors=grey,label='Expanded Bang cut-offs')
ax.set(xlim=(-2,102),ylim=(-2,103),xlabel='Adults referred for testing (%)',ylabel='HbA1c-positive adults identified (%)')
ax.set_title('A  Attainable operating points',loc='left',fontweight='bold',pad=16)
ax.legend(loc='lower right',frameon=False,fontsize=12,handlelength=1.2)
ax.grid(axis='y',color='#E2E5E8',lw=.6)
ax=axes[1]
for yy,(t,x,_) in enumerate(rows):
    key=f'common::budget::Bang_composite::k{t}::difference::sensitivity_pct'
    r=S.loc[key];q=SI.loc[key]
    ax.plot([q.simultaneous_lower,q.simultaneous_upper],[yy,yy],color=blue,lw=1.6)
    ax.scatter(r.estimate,yy,color=blue,s=22,zorder=4)
    ax.plot([r.percentile_lower,r.percentile_upper],[yy-.15,yy-.15],color=orange,lw=1.4)
ax.axvline(0,color=grey,ls=':',lw=.9)
ax.set(yticks=range(len(rows)),yticklabels=[f'≥{t} ({x:.1f}%)' for t,x,_ in rows],
       xlabel='Equation − Bang sensitivity (pp)',xlim=(-7,34),ylim=(-.7,len(rows)-.3))
ax.set_title('B  Gains at each Bang budget',loc='left',fontweight='bold',pad=16)
ax.tick_params(axis='y',length=0);ax.grid(axis='x',color='#E2E5E8',lw=.6)
f.text(.085,.055,'Common sample: n = 11,624; 348 HbA1c positives. Bang points are not interpolated.\nPanel A band: pointwise percentile 95% intervals. Panel B blue: simultaneous 95% intervals\nacross 96 declared contrasts; orange: pointwise percentile 95% intervals. Budgets are retrospective.',fontsize=12,color=grey)
save(f,'Figure_B_budget_comparison')
metadata={'dpi':450,'native_width_inches':9.6,'report_width_inches':6.5,
          'minimum_main_native_font_pt':12,'legend_native_font_pt':12,
          'minimum_main_embedded_font_pt':12*6.5/9.6,'minimum_legend_embedded_font_pt':12*6.5/9.6,
          'figures':['Figure_A_auc_and_cycle_stability','Figure_B_budget_comparison'],'participant_data_used':False}
(O/'figure_metadata.json').write_text(json.dumps(metadata,indent=2),encoding='utf-8')
print(json.dumps(metadata,indent=2))
