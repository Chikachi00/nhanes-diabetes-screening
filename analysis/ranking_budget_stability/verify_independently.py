"""Independent Python reconstruction from frozen CSV inputs; no R evaluation code imported.
python verify_independently.py INPUT_ROOT V11_RESULTS RESULTS OUTPUT_JSON
"""
from pathlib import Path
import sys,json,hashlib
import numpy as np
import pandas as pd

V,A,O=map(Path,sys.argv[1:4]); dest=Path(sys.argv[4])
pr=pd.read_csv(V/'reproducibility/main_archive/evaluation_predictions.csv')
fb=pd.read_csv(V/'reproducibility/main_inputs/full_base.csv').set_index('SEQN').loc[pr.SEQN].reset_index()
mp=pd.read_csv(A/'participant_mapping_PRIVATE.csv')
assert np.array_equal(pr.SEQN,mp.SEQN) and np.array_equal(fb.SEQN,pr.SEQN)
assert not pr.SEQN.duplicated().any()
y=(fb.hba1c>=6.5).to_numpy(float); basew=fb.wt.to_numpy(float)
agepoints=np.select([fb.age<40,fb.age<50,fb.age<60],[0,1,2],default=3)
bmi=np.select([fb.bmi<25,fb.bmi<30,fb.bmi<40],[0,1,2],default=3)
wa=fb.waist/2.54
waist=np.where(fb.male==1,np.select([wa<37,wa<40,wa<50],[0,1,2],default=3),
              np.select([wa<31.5,wa<35,wa<49],[0,1,2],default=3))
score=(agepoints+fb.male+fb.famhx+fb.htn+np.maximum(bmi,waist)-(1-fb.inactive)).to_numpy(float)
assert np.array_equal(score>=5,pr.Bang.astype(bool))
scores={'Bang_self':score,'Bang_composite':score-fb.htn.to_numpy()+mp.composite.to_numpy()}
X=np.column_stack([np.ones(len(fb)),fb.age/10-5,10*(fb.waist/fb.height-.55),fb.male,fb.famhx,fb.htn,fb.inactive])
terms=['(Intercept)','a','w','male','famhx','htn','inactive']
locked=pd.read_csv(V/'reproducibility/main_archive/locked_strategy_parameters.csv').set_index('term').value
pars=pd.read_csv(V/'reproducibility/main_archive/bootstrap_parameters.csv')
assert np.array_equal(pars.replicate,np.arange(1,1001))
mult=pd.read_csv(V/'reproducibility/main_inputs/replicate_multipliers.csv').set_index('PSU_key')
keys=fb.cycle+'/'+fb.strata.astype(int).astype(str)+'/'+fb.psu.astype(int).astype(str)
masks={'common':mp.common.to_numpy(bool),'strict_BP':mp.strict_BP.to_numpy(bool)}
for cy in 'HIJ': masks['cycle_'+cy]=masks['common']&(fb.cycle.to_numpy()==cy)
for cy in 'HIJ': masks['without_'+cy]=masks['common']&(fb.cycle.to_numpy()!=cy)
cutoffs=range(int(min(s[masks['common']].min() for s in scores.values())),int(max(s[masks['common']].max() for s in scores.values()))+2)
metricnames=['sensitivity_pct','referral_pct','ppv_pct','found_per1000','referrals_per1000']

def div(a,b): return a/b if b>0 else np.nan
def metrics(r,yy,ww):
    total=ww.sum(); pos=ww@yy; found=ww@(yy*r); ref=ww@r
    return np.array([100*div(found,pos),100*div(ref,total),100*div(found,ref),1000*div(found,total),1000*div(ref,total)])
def auc(s,yy,ww):
    # Group equal scores then accumulate the lower-scoring negative weights.
    k=ww>0;u,inv=np.unique(s[k],return_inverse=True)
    positive=np.bincount(inv,weights=ww[k]*yy[k],minlength=len(u))
    negative=np.bincount(inv,weights=ww[k]*(1-yy[k]),minlength=len(u))
    below=np.cumsum(negative)-negative
    return div(positive@(below+.5*negative),positive.sum()*negative.sum())
def threshold(p,ww,target):
    # Compute attainable blocks with grouping, independently of the R cumulative-end implementation.
    k=ww>0;u,inv=np.unique(p[k],return_inverse=True)
    weights=np.bincount(inv,weights=ww[k],minlength=len(u))[::-1]
    attainable=np.r_[0,np.cumsum(weights)/weights.sum()]
    distance=np.abs(attainable-target);candidates=np.flatnonzero(distance<=distance.min()+1e-14)
    i=candidates[np.argmin(attainable[candidates])]
    return np.inf if i==0 else u[::-1][i-1]
def compute(p,ww,t0):
    out={}
    def add(prefix,v): out.update({prefix+'::'+name:val for name,val in zip(metricnames,v)})
    for sn,k in masks.items():
        yy=y[k];w=ww[k];pp=p[k];ss={n:s[k] for n,s in scores.items()}
        au={n:auc(v,yy,w) for n,v in {'equation':pp,**ss}.items()}
        au.update(equation_minus_self=au['equation']-au['Bang_self'],equation_minus_composite=au['equation']-au['Bang_composite'],composite_minus_self=au['Bang_composite']-au['Bang_self'])
        out.update({sn+'::auc::'+n:v for n,v in au.items()})
        st={'R0':metrics(pp>=t0,yy,w)}
        for n,s in ss.items():
            v=metrics(s>=5,yy,w);st[n]=v
            eq=metrics(pp>=threshold(pp,w,v[1]/100),yy,w);st['Eq_to_'+n]=eq
            add(sn+'::matched5::difference_'+n,eq-v)
            add(sn+'::fixed::R0_minus_'+n,st['R0']-v)
        for n,v in st.items(): add(sn+'::fixed::'+n,v)
        if sn in ['common','strict_BP']:
            for n,s in ss.items():
                for t in cutoffs:
                    v=metrics(s>=t,yy,w);eq=metrics(pp>=threshold(pp,w,v[1]/100),yy,w)
                    prefix=f'{sn}::budget::{n}::k{t}'
                    for entity,value in [('Bang',v),('Equation',eq),('difference',eq-v)]: add(prefix+'::'+entity,value)
        if sn=='common':
            for g in range(10,61): add(sn+f'::grid::g{g:02d}',metrics(pp>=threshold(pp,w,g/100),yy,w))
    return out

# Analytic checks for ties, zero referral, full referral, and unequal sampling weights.
assert abs(auc(np.array([1.,1.,2.]),np.array([1.,0.,1.]),np.ones(3))-.75)<1e-14
assert abs(auc(np.array([1.,1.,2.]),np.array([1.,0.,1.]),np.array([2.,3.,4.]))-5/6)<1e-14
assert np.isinf(threshold(np.array([.5,.5,.1]),np.array([1.,1.,2.]),.25))
assert threshold(np.array([.5,.5,.1]),np.array([1.,1.,2.]),1)==.1
point=pd.read_csv(O/'point_estimates.csv').set_index('key').estimate
reps=pd.read_csv(O/'bootstrap_replicates.csv').set_index('replicate')
S=pd.read_csv(O/'summary_all.csv').set_index('key')
errors=[];nanchecks=0
for b in [0,1,2,100,500,1000]:
    row=locked if b==0 else pars.iloc[b-1]
    w=basew if b==0 else basew*mult.loc[keys,'X'+str(b)].to_numpy(float)
    p=1/(1+np.exp(-(X@row.loc[terms].to_numpy(float))))
    if b==0: assert np.max(np.abs(p-pr.probability))<1e-12
    got=compute(p,w,row.THRESHOLD);expected=point if b==0 else reps.loc[b]
    assert set(got)==set(expected.index)
    for key,value in got.items():
        exp=expected[key]
        if np.isnan(value): assert np.isnan(exp);nanchecks+=1
        else: errors.append(abs(value-exp))
maximum=max(errors);assert maximum<1e-9,maximum
interval_errors=[]
for key,r in S.iterrows():
    vals=reps[key].to_numpy();finite=np.isfinite(vals);assert finite.sum()==r.valid_B
    if finite.any():
        se=vals[finite].std(ddof=1);lo,hi=np.quantile(vals[finite],[.025,.975])
        interval_errors.extend([abs(se-r.se),abs(r.estimate-1.96*se-r.normal_lower),abs(r.estimate+1.96*se-r.normal_upper),abs(lo-r.percentile_lower),abs(hi-r.percentile_upper)])
assert max(interval_errors)<1e-9,max(interval_errors)
sim=pd.read_csv(O/'budget_simultaneous_intervals.csv').set_index('key')
positive=sim.se>1e-12
z=(reps[sim.index[positive]].to_numpy()-sim.loc[positive,'estimate'].to_numpy())/sim.loc[positive,'se'].to_numpy()
mx=np.abs(z).max(axis=1);crit=np.quantile(mx,.95)
simerr=max(np.abs(sim.simultaneous_lower-(sim.estimate-crit*sim.se)).max(),np.abs(sim.simultaneous_upper-(sim.estimate+crit*sim.se)).max())
assert simerr<1e-9 and np.max(np.abs(sim.simultaneous_critical-crit))<1e-9
baseline=pd.read_csv(A/'S25_hypertension_strategy_comparisons.csv').set_index('key')
v11errs=[]
for n in ['Bang_self','Bang_composite','R0']:
    for m in ['sensitivity','referral','ppv','found_per1000','referrals_per1000']:
        nm=m+'_pct' if m in ['sensitivity','referral','ppv'] else m
        v11errs.append(abs(point[f'common::fixed::{n}::{nm}']-baseline.loc[f'{n}::overall::{m}','estimate']))
assert max(v11errs)<1e-9
coverage={sn:{'n':int(k.sum()),'cases':int(y[k].sum())} for sn,k in masks.items()}
result={'verification_status':'VERIFIED','selected_replicates':[1,2,100,500,1000],
        'independent_numeric_checks':len(errors),'expected_NA_checks':nanchecks,'maximum_point_replicate_error':maximum,
        'maximum_interval_error':max(interval_errors),'maximum_simultaneous_interval_error':float(simerr),
        'simultaneous_critical':float(crit),'v11_reproduction_maximum_error':max(v11errs),
        'requested_B':1000,'replicate_rows':len(reps),'metric_columns':len(point),'coverage':coverage,
        'analytic_edge_case_checks':4,'external_validation_added':False,'model_refitted':False,
        'selection_uncertainty_covered':False}
dest.parent.mkdir(parents=True,exist_ok=True)
dest.write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(result,ensure_ascii=False,indent=2))
