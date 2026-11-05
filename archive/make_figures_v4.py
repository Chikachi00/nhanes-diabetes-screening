from pathlib import Path
import sys, os, json, hashlib
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'plot_dependencies'))
os.environ['MPLCONFIGDIR']=str(HERE/'qa/mplconfig')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import scienceplots
import numpy as np
import pandas as pd

ROOT=HERE.parents[1]
MAIN=ROOT/'outputs/diabetes_strategy_completion_20261102/results'
CAL=ROOT/'outputs/diabetes_calibration_closeout_20261103/results'
FIG=HERE/'figures';DAT=FIG/'source_data';QA=HERE/'qa'
for p in [FIG,DAT,QA]:p.mkdir(parents=True,exist_ok=True)
plt.style.use(['science','no-latex'])
STYLE={'font.family':'Arial','font.size':9.5,'axes.labelsize':9.5,'axes.titlesize':9.5,
 'xtick.labelsize':9,'ytick.labelsize':9,'legend.fontsize':9,
 'axes.linewidth':.8,'lines.linewidth':1.2,'lines.markersize':5,
 'xtick.direction':'out','ytick.direction':'out','xtick.top':False,'ytick.right':False,
 'axes.spines.top':False,'axes.spines.right':False,'axes.grid':False,
 'figure.facecolor':'white','axes.facecolor':'white','text.color':'#222222',
 'axes.labelcolor':'#222222','xtick.color':'#222222','ytick.color':'#222222',
 'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.dpi':300}
plt.rcParams.update(STYLE)
COL={'Bang':'#62666D','R0':'#0072B2','R1':'#D55E00','R4':'#009E73',
     'high':'#0072B2','low':'#D55E00','original':'#0072B2',
     'intercept_loco':'#CC79A7','joint_loco':'#009E73'}
MM=1/25.4; W=165*MM
audit={'matplotlib':matplotlib.__version__,'scienceplots':scienceplots.__version__ if hasattr(scienceplots,'__version__') else '2.2.2',
       'style':STYLE,'colors':COL,'figure_width_mm':165,'minimum_used_font_pt':8.5,
       'sources':{},'files':[],'scope':'Presentation only. Frozen estimates, CIs, models and thresholds are unchanged.'}
def source(path):
    audit['sources'][str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest()
    return pd.read_csv(path)
def save(fig,name):
    for ext in ['pdf','svg','png']:
        p=FIG/f'{name}.{ext}';fig.savefig(p,facecolor='white')
        audit['files'].append({'path':str(p.relative_to(HERE)),'bytes':p.stat().st_size})
    plt.close(fig)
def panel(ax,letter):
    ax.text(-.03,1.08,letter,transform=ax.transAxes,fontweight='bold',fontsize=10,va='bottom')
flow=source(MAIN/'cohort_flow.csv');flow.to_csv(DAT/'Figure1_cohort_flow.csv',index=False)
fig=plt.figure(figsize=(W,181*MM));ax=fig.add_axes([0,0,1,1]);ax.set_axis_off()
stages=[(0,'Interviewed participants'),(1,'MEC examination; positive weight'),
        (3,'Age ≥20; not pregnant'),(4,'No self-reported diabetes'),
        (6,'Measured HbA1c; access known'),(8,'Common complete-case sample')]
ys=np.linspace(.90,.34,len(stages)); bw=.445;bh=.079
for x,pe,title in [(.035,'dev','Development 2007–2012'),(.52,'eval','Evaluation 2013–2018')]:
    ax.text(x+bw/2,.976,title,ha='center',va='top',fontsize=10,fontweight='bold')
    for j,((order,label),y0) in enumerate(zip(stages,ys)):
        r=flow[(flow.period==pe)&(flow.order==order)].iloc[0]
        box=FancyBboxPatch((x,y0-bh/2),bw,bh,boxstyle='round,pad=0.005,rounding_size=0.005',
             transform=ax.transAxes,edgecolor='#62666D',linewidth=.8,facecolor='white')
        ax.add_patch(box)
        line=label+f'\nN = {int(r.n):,}'
        if order in [4,6,8]:line+=f'; positive = {int(r.observed_positive_n):,}'
        ax.text(x+bw/2,y0,line,ha='center',va='center',fontsize=9,linespacing=1.25)
        if j<len(stages)-1:
            ax.annotate('',xy=(x+bw/2,ys[j+1]+bh/2+.005),xytext=(x+bw/2,y0-bh/2-.005),
                arrowprops={'arrowstyle':'->','color':'#62666D','lw':.9})
ax.text(.5,.242,'Separate evaluation branches (not additive)',ha='center',fontsize=9.5,fontweight='bold')
branches=[('Fasting tests\n2013–2018\nN = 5,659',.035),
 ('Three glycaemic tests\n2013–2016\nN = 3,265',.355),
 ('MI: observed-HbA1c target\n2013–2018\nN = 12,724; positive = 390',.675)]
for text,x in branches:
    box=FancyBboxPatch((x,.105),.29,.105,boxstyle='round,pad=0.005,rounding_size=0.005',
        transform=ax.transAxes,edgecolor='#62666D',linewidth=.8,facecolor='white')
    ax.add_patch(box);ax.text(x+.145,.157,text,ha='center',va='center',fontsize=8.5,linespacing=1.3)
ax.text(.5,.057,'All counts are unweighted. Positive: HbA1c ≥6.5%.',ha='center',fontsize=8.5)
ax.text(.5,.029,'MI excludes 599 participants without measured HbA1c from evaluation.',ha='center',fontsize=8.5)
save(fig,'Figure1_cohort_flow')

e=source(ROOT/'outputs/diabetes_screening_v6/results/evaluation_expanded.csv')
e=e[e.common_complete.astype(str).str.lower().isin(['true','1'])].copy()
assert len(e)==11848 and int(e['case'].sum())==356
dist=[]
for gr in ['high','low']:
    g=e[(e['case']==1)&(e.access==gr)]
    assert len(g)=={'high':251,'low':105}[gr]
    total=g.wt.sum()
    for score in range(-1,11):
        ss=g[g.score==score]
        dist.append({'group':gr,'score':score,'positive_n':len(g),'unweighted_n_at_score':len(ss),
          'weighted_case_share_pct':ss.wt.sum()/total*100})
dist=pd.DataFrame(dist);dist.to_csv(DAT/'Figure2_case_score_distribution.csv',index=False)
fig,axs=plt.subplots(1,2,figsize=(W,85*MM),sharey=True)
fig.subplots_adjust(left=.10,right=.985,bottom=.20,top=.82,wspace=.16)
for ax,gr,letter in zip(axs,['high','low'],['A','B']):
    g=dist[dist.group==gr];n=int(g.positive_n.iloc[0]);s=g.loc[g.score>=5,'weighted_case_share_pct'].sum()
    ax.axvspan(4.5,10.7,color='#F0F1F3',zorder=0)
    ax.bar(g.score,g.weighted_case_share_pct,color=COL[gr],width=.76,edgecolor='white',linewidth=.35)
    ax.axvline(4.5,color='#333333',linestyle='--',linewidth=1)
    ax.text(4.7,33.5,'Score ≥5',fontsize=8.5,va='top')
    ax.set(xlim=(-1.7,10.7),ylim=(0,36),xticks=[-1,1,3,5,7,9],yticks=[0,10,20,30],xlabel='Bang score')
    ax.set_title(('Higher' if gr=='high' else 'Lower')+f' access\nn = {n}; sensitivity = {s:.1f}%',pad=8)
    panel(ax,letter)
axs[0].set_ylabel('Weighted share of positive cases (%)')
save(fig,'Figure2_case_score_distribution')

age=source(ROOT/'outputs/diabetes_age_decomposition_20261103/Figure3_data.csv')
age['comparison']=age['comparison'].replace({'S1_minus_Bang':'R1_minus_Bang'})
age.to_csv(DAT/'Figure3_age_differences.csv',index=False)
fig,axs=plt.subplots(1,2,figsize=(W,94*MM),sharey=True)
fig.subplots_adjust(left=.225,right=.985,bottom=.20,top=.83,wspace=.38)
for ax,measure,letter,lims,ticks in zip(axs,['sensitivity','referral'],['A','B'],[(-10,85),(-15,17)],[[-10,0,20,40,60,80],[-10,0,10]]):
    g=age[age.measure==measure].set_index('age_group').loc[['young','mid','old']]
    y=np.array([2,1,0]);q=g.difference_pp.to_numpy();lo=g.lower95_pp.to_numpy();hi=g.upper95_pp.to_numpy()
    ax.axvline(0,color='#62666D',linestyle='--',linewidth=.9,zorder=0)
    ax.errorbar(q,y,xerr=np.vstack([q-lo,hi-q]),fmt='o',color=COL['R1'],capsize=3,linewidth=1.2,markersize=5.4)
    for qi,yi in zip(q,y):ax.text(qi,yi+.21,f'{qi:+.2f}',ha='center',fontsize=9)
    ax.set(xlim=lims,xticks=ticks,ylim=(-.45,2.48),xlabel='Difference (percentage points)')
    ax.set_title('Sensitivity' if measure=='sensitivity' else 'Referral proportion',pad=9)
    panel(ax,letter)
axs[0].set_yticks([2,1,0],['20–39 years\nn = 43 positive','40–59 years\nn = 145 positive','≥60 years\nn = 168 positive'])
axs[0].tick_params(axis='y',length=0,pad=7)
axs[1].tick_params(axis='y',left=False)
save(fig,'Figure3_age_differences')

curves=source(CAL/'calibration_curves.csv');curves.to_csv(DAT/'Figure4_calibration_curves.csv',index=False)
methods=list(curves.method.unique());assert set(methods)=={'original','intercept_loco','joint_loco'},methods
fig,axs=plt.subplots(2,3,figsize=(W,123*MM))
fig.subplots_adjust(left=.105,right=.985,bottom=.13,top=.91,wspace=.29,hspace=.32)
for j,(me,title,oe,sl) in enumerate([('original','Original',.707,.801),('intercept_loco','Intercept update',1.002,.800),('joint_loco','Intercept + slope',.999,.978)]):
    c=curves[curves.method==me]
    for row in [0,1]:
        ax=axs[row,j];limit=100 if row==0 else 12
        ax.plot([0,limit],[0,limit],color='#62666D',linestyle='--',linewidth=.8)
        ax.fill_between(c.predicted*100,c.lower95*100,c.upper95*100,color=COL[me],alpha=.17,linewidth=0)
        ax.plot(c.predicted*100,c.observed_curve*100,color=COL[me])
        ax.set(xlim=(0,limit),ylim=(0,limit));ax.set_aspect('equal',adjustable='box')
        ticks=[0,50,100] if row==0 else [0,4,8,12]
        ax.set_xticks(ticks);ax.set_yticks(ticks)
        if row==0:
            ax.set_title(title,pad=12)
            ax.text(.05,.94,f'O:E = {oe:.3f}\nSlope = {sl:.3f}',transform=ax.transAxes,va='top',fontsize=8.5,linespacing=1.35)
        ax.text(-.10,1.035,chr(ord('A')+row*3+j),transform=ax.transAxes,fontweight='bold',fontsize=10)
        if j==0:ax.set_ylabel('Observed positivity (%)',labelpad=4)
fig.supxlabel('Predicted probability (%)',fontsize=9.5,y=.025)
save(fig,'Figure4_calibration')

perf=source(MAIN/'main_performance_full_pipeline.csv')
perf['strategy']=perf.strategy.replace({'S0':'R0','S1':'R1','S4':'R4'})
out=perf[perf.group.isin(['overall','old'])&perf.measure.isin(['sensitivity','referral'])&perf.strategy.isin(['Bang','R0','R1','R4'])]
out.to_csv(DAT/'FigureS5_strategy_performance.csv',index=False)
fig,axs=plt.subplots(2,2,figsize=(W,118*MM))
fig.subplots_adjust(left=.10,right=.98,bottom=.11,top=.91,hspace=.45,wspace=.30)
for row,group in enumerate(['overall','old']):
    for j,me in enumerate(['sensitivity','referral']):
        ax=axs[row,j];n='Overall' if group=='overall' else 'Age ≥60 years'
        ax.set_title(n+': '+('sensitivity' if me=='sensitivity' else 'referral'),pad=8)
        for i,name in enumerate(['Bang','R0','R1','R4']):
            r=out[(out.group==group)&(out.measure==me)&(out.strategy==name)].iloc[0]
            ax.errorbar(100*r.estimate,3-i,xerr=np.array([[100*(r.estimate-r.lower95)],[100*(r.upper95-r.estimate)]]),fmt='o',color=COL[name],capsize=2.5,markersize=4.5)
        ax.set(yticks=[3,2,1,0],yticklabels=['Bang','R0','R1','R4'],xlim=(0,105),ylim=(-.4,3.4),xticks=[0,25,50,75,100],xlabel='%')
        panel(ax,chr(ord('A')+row*2+j))
save(fig,'FigureS5_strategy_performance')
(QA/'figure_manifest.json').write_text(json.dumps(audit,indent=2,ensure_ascii=False),encoding='utf-8')
print('Created five consistent figure sets (PDF/SVG/PNG) and exact source CSVs.')
