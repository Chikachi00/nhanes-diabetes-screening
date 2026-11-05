from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import t
P=Path(__file__).resolve().parent;WS=P.parents[1]
d=pd.read_csv(WS/'outputs/diabetes_full_pipeline_20261029/inputs/evaluation_inputs.csv');f=pd.read_csv(P/'results/evaluation_predictions.csv');assert np.array_equal(d.SEQN,f.SEQN)
fr=pd.read_csv(WS/'outputs/diabetes_full_pipeline_20261029/inputs/full_examined_PSU_frame.csv');keys=lambda z:z.cycle+'/'+z.strata.astype(int).astype(str)+'/'+z.psu.astype(int).astype(str)
ix=pd.Index(keys(fr)).get_indexer(keys(d));assert(ix>=0).all();st=fr.cycle+'/'+fr.strata.astype(int).astype(str);groups=list(fr.groupby(st).indices.values());w=d.wt.to_numpy();y=d['case'].to_numpy()
masks={'overall':np.ones(len(d),bool),'young':d.age_grp.eq('20-39').to_numpy(),'mid':d.age_grp.eq('40-59').to_numpy(),'old':d.age_grp.eq('>=60').to_numpy(),'high':d.access.eq('high').to_numpy(),'low':d.access.eq('low').to_numpy()}
def dfden(den):u=np.unique(ix[den>0]);return len(u)-len(np.unique(st.iloc[u]))
def vu(u):a=np.zeros(len(fr));np.add.at(a,ix,u);return sum(len(g)/(len(g)-1)*np.sum((a[g]-a[g].mean())**2)for g in groups)
items={};meta={}
for strategy in ['Bang','S0','S1','S4']:
 r=f[strategy].to_numpy()
 for g,k in masks.items():
  forms={'sensitivity':(y*r*k,y*k),'specificity':((1-y)*(1-r)*k,(1-y)*k),'ppv':(y*r*k,r*k),'referral':(r*k,k),'found_per1000':(1000*y*r*k,k),'missed_per1000':(1000*y*(1-r)*k,k),'referrals_per1000':(1000*r*k,k)}
  for me,(num,den)in forms.items():
   D=w@den;q=w@num/D;key=f'{strategy}::{g}::{me}';items[key]=(q,w*(num-q*den)/D,dfden(den));meta[key]=(int(k.sum()),int(y[k].sum()))
 a,ua,da=items[f'{strategy}::high::sensitivity'];b,ub,db=items[f'{strategy}::low::sensitivity'];key=f'{strategy}::access_gap::sensitivity';items[key]=(a-b,ua-ub,dfden(y));meta[key]=(len(d),int(y.sum()))
for name,(s1,s0)in {'S0_minus_Bang':('S0','Bang'),'S1_minus_Bang':('S1','Bang'),'S1_minus_S0':('S1','S0'),'S4_minus_S1':('S4','S1'),'S4_minus_Bang':('S4','Bang')}.items():
 for ka in [k for k in items if k.startswith(s1+'::')]:
  kb=s0+'::'+ka.split('::',1)[1];kk=name+'::'+ka.split('::',1)[1];a,ua,dc=items[ka];b,ub,_=items[kb]
  if ka.endswith('::ppv'):
   gr=ka.split('::')[1];dc=dfden((f[s1].to_numpy()+f[s0].to_numpy())*masks[gr])
  items[kk]=(a-b,ua-ub,dc);meta[kk]=meta[ka]
rows=[]
for key,(q,u,dc)in items.items():
 s,g,me=key.split('::');var=vu(u);se=np.sqrt(var);crit=t.ppf(.975,dc);n,c=meta[key];rows.append(dict(strategy=s,group=g,measure=me,n=n,cases=c,estimate=q,variance=var,se=se,domain_df=dc,lower95=q-crit*se,upper95=q+crit*se,analysis='fixed_rule_Taylor'))
out=pd.DataFrame(rows);out.to_csv(P/'results/fixed_rule_taylor.csv',index=False)
orig=pd.read_csv(WS/'outputs/teammate_mi_segment_review_20261101/independent_fixed_strategy_checks.csv');mapping={'Bang':'Bang','S0':'S0_single','S1':'S1_single_agecut','S4':'S4_test_all_60plus'};groupsmap={'overall':'overall','young':'20-39','mid':'40-59','old':'>=60'};errs=[]
for z in out.itertuples():
 if z.strategy not in mapping or z.group not in groupsmap:continue
 a=orig.loc[orig.period.eq('eval')&orig.strategy.eq(mapping[z.strategy])&orig.group.eq(groupsmap[z.group])&orig.measure.eq(z.measure)]
 if len(a):errs.append(max(abs(z.estimate*100-a.estimate.iloc[0]),abs(z.se*100-a.se.iloc[0])))
assert max(errs)<1e-9
pd.DataFrame({'max_source_point_SE_error':[max(errs)],'n_source_checks':[len(errs)]}).to_csv(P/'qa/fixed_Taylor_source_checks.csv',index=False)
print('Fixed-rule Taylor completed:',len(out),'rows; source comparisons:',len(errs))
