from pathlib import Path
import sys,os,json,hashlib,shutil
BASE=Path(__file__).resolve().parent
sys.path.insert(0,str(BASE/'glucobench_models_2026-10-14/_deps'))
OUT=BASE/'diabetes_manuscript_revision_20261106_v5'
os.environ['MPLCONFIGDIR']=str(OUT/'qa/mplconfig')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import pandas as pd
import numpy as np
plt.rcParams.update({'font.family':'Arial','font.size':9.5,'axes.labelsize':9.5,'axes.titlesize':9.5,'xtick.labelsize':9,'ytick.labelsize':9,'legend.fontsize':9,'axes.linewidth':.8,'lines.linewidth':1.2,'xtick.direction':'out','ytick.direction':'out','xtick.top':False,'ytick.right':False,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white','axes.facecolor':'white','text.color':'#222222','axes.labelcolor':'#222222','xtick.color':'#222222','ytick.color':'#222222','pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none'})
MM=1/25.4;W=165*MM;FIG=OUT/'figures';DATA=FIG/'source_data';archive=FIG/'v4_archive';archive.mkdir(exist_ok=True)
for name in ['Figure1_cohort_flow','Figure2_case_score_distribution','Figure4_calibration']:
    for ext in ['png','pdf','svg']:shutil.copy2(BASE/'diabetes_manuscript_revision_20261105_v4/figures'/(name+'.'+ext),archive/(name+'.'+ext))
receipt={}
def save(fig,name,number):
    for ext in ['png','pdf','svg']:fig.savefig(FIG/(name+'.'+ext),dpi=390,facecolor='white')
    plt.close(fig)
    shutil.copy2(FIG/(name+'.png'),OUT/f'submission/figures/Figure {number}.png')
    receipt[number]={'name':name,'source_CSV_sha256':hashlib.sha256((DATA/(name+('.csv' if number!=4 else '_curves.csv'))).read_bytes()).hexdigest() if number!=4 else hashlib.sha256((DATA/'Figure4_calibration_curves.csv').read_bytes()).hexdigest(),'raster_dpi':390}
flow=pd.read_csv(DATA/'Figure1_cohort_flow.csv')
fig=plt.figure(figsize=(W,195*MM));ax=fig.add_axes([0,0,1,1]);ax.set_axis_off()
stages=[(0,'Interviewed participants'),(1,'MEC examination; positive weight'),(3,'Age ≥20; not pregnant'),(4,'No self-reported diabetes'),(6,'Measured HbA1c; access known'),(8,'Common complete-case sample')]
ys=np.linspace(.90,.32,len(stages));bh=.094;bw=.445
display=[]
for x,pe,title in [(.035,'dev','Development 2007–2012'),(.52,'eval','Evaluation 2013–2018')]:
    ax.text(x+bw/2,.979,title,ha='center',va='top',fontsize=10,fontweight='bold')
    prev=None
    for j,((order,label),y) in enumerate(zip(stages,ys)):
        r=flow[(flow.period==pe)&(flow.order==order)].iloc[0];n=int(r.n)
        lines=[label,f'N = {n:,}'+(f'; positive = {int(r.observed_positive_n):,}' if order in [4,6,8] else '')]
        ex=prev-n if prev is not None else None
        if ex is not None:lines.append(f'Excluded since previous box: {ex:,}')
        display.append({'period':pe,'order':order,'n':n,'excluded_since_previous_displayed':ex})
        box=FancyBboxPatch((x,y-bh/2),bw,bh,boxstyle='round,pad=0.005,rounding_size=0.005',transform=ax.transAxes,edgecolor='#62666D',linewidth=.8,facecolor='white');ax.add_patch(box)
        ax.text(x+bw/2,y,'\n'.join(lines),ha='center',va='center',fontsize=8.5,linespacing=1.3)
        if j<len(stages)-1:ax.annotate('',xy=(x+bw/2,ys[j+1]+bh/2+.005),xytext=(x+bw/2,y-bh/2-.005),arrowprops={'arrowstyle':'->','color':'#62666D','lw':.9})
        prev=n
ax.text(.5,.236,'Separate evaluation branches (not additive)',ha='center',fontsize=9.5,fontweight='bold')
for t,x in [('Fasting tests\n2013–2018\nN = 5,659',.035),('Three glycaemic tests\n2013–2016\nN = 3,265',.355),('MI: observed-HbA1c target\n2013–2018\nN = 12,724; positive = 390',.675)]:
    ax.add_patch(FancyBboxPatch((x,.105),.29,.105,boxstyle='round,pad=0.005,rounding_size=0.005',transform=ax.transAxes,edgecolor='#62666D',linewidth=.8,facecolor='white'));ax.text(x+.145,.157,t,ha='center',va='center',fontsize=8.5,linespacing=1.3)
ax.text(.5,.064,'Unweighted counts; exclusions compare displayed boxes.',ha='center',fontsize=8.5)
ax.text(.5,.043,'Age and pregnancy eligibility are combined. Positive: HbA1c ≥6.5%.',ha='center',fontsize=8.5)
ax.text(.5,.022,'MI excludes 599 participants without measured HbA1c from evaluation.',ha='center',fontsize=8.5)
pd.DataFrame(display).to_csv(DATA/'Figure1_displayed_exclusions.csv',index=False)
save(fig,'Figure1_cohort_flow',1)

dist=pd.read_csv(DATA/'Figure2_case_score_distribution.csv')
nonzero=dist.loc[dist.weighted_case_share_pct>0,'score'];xmin=int(nonzero.min());xmax=int(nonzero.max())
fig,axs=plt.subplots(1,2,figsize=(W,85*MM),sharey=True);fig.subplots_adjust(left=.10,right=.985,bottom=.20,top=.82,wspace=.16)
for ax,gr,letter,col in zip(axs,['high','low'],['A','B'],['#0072B2','#D55E00']):
    g=dist[dist.group==gr];n=int(g.positive_n.iloc[0]);sens=g.loc[g.score>=5,'weighted_case_share_pct'].sum()
    ax.axvspan(4.5,xmax+.65,color='#F0F1F3',zorder=0);ax.bar(g.score,g.weighted_case_share_pct,color=col,width=.76,edgecolor='white',linewidth=.35)
    ax.axvline(4.5,color='#333333',linestyle='--',linewidth=1);ax.text(4.7,33.5,'Score ≥5',fontsize=8.5,va='top')
    ax.set(xlim=(xmin-.65,xmax+.65),ylim=(0,36),xticks=list(range(xmin,xmax+1,2)),yticks=[0,10,20,30],xlabel='Bang score')
    ax.set_title(('Higher' if gr=='high' else 'Lower')+f' access\nn = {n}; sensitivity = {sens:.1f}%',pad=8)
    ax.text(-.03,1.08,letter,transform=ax.transAxes,fontweight='bold',fontsize=10,va='bottom')
axs[0].set_ylabel('Weighted share of positive cases (%)');save(fig,'Figure2_case_score_distribution',2)

curves=pd.read_csv(DATA/'Figure4_calibration_curves.csv')
fig,axs=plt.subplots(1,3,figsize=(W,85*MM));fig.subplots_adjust(left=.095,right=.985,bottom=.22,top=.82,wspace=.29)
for j,(me,title,col,oe,sl) in enumerate([('original','Original','#0072B2',.707,.801),('intercept_loco','Intercept update','#CC79A7',1.002,.800),('joint_loco','Intercept + slope','#009E73',.999,.978)]):
    ax=axs[j];c=curves[curves.method==me]
    ax.plot([0,12],[0,12],color='#62666D',linestyle='--',linewidth=.8);ax.fill_between(c.predicted*100,c.lower95*100,c.upper95*100,color=col,alpha=.17,linewidth=0);ax.plot(c.predicted*100,c.observed_curve*100,color=col)
    ax.set(xlim=(0,12),ylim=(0,12),xticks=[0,4,8,12],yticks=[0,4,8,12]);ax.set_aspect('equal',adjustable='box');ax.set_title(title,pad=8)
    ax.text(.045,.94,f'O:E = {oe:.3f}\nSlope = {sl:.3f}',transform=ax.transAxes,va='top',fontsize=8,linespacing=1.25)
    ax.text(-.10,1.035,chr(65+j),transform=ax.transAxes,fontweight='bold',fontsize=10)
    if j==0:ax.set_ylabel('Observed positivity (%)')
fig.supxlabel('Predicted probability (%)',fontsize=9.5,y=.08);save(fig,'Figure4_calibration',4)
shutil.copy2(FIG/'Figure3_age_differences.png',OUT/'submission/figures/Figure 3.png')
(OUT/'qa/figure_revision_v5.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2),encoding='utf-8')
print('Figures 1,2,4 redrawn from unchanged CSVs; Figure 3 retained.')
