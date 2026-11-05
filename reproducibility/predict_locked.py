"""Research-only frozen main-equation arithmetic, with explicit complete inputs.

Usage: python predict_locked.py people.csv predictions.csv
Required columns: age,waist,height,male,famhx,htn,inactive.
Waist and height must use the same unit. Age >=80 uses NHANES's 80 top-code.
This script does not perform diagnosis, imputation, or probability updating.
"""
from pathlib import Path
import csv, math, sys

PARAMETERS = Path(__file__).resolve().parent/'main_archive/locked_strategy_parameters.csv'
with PARAMETERS.open(encoding='utf-8-sig') as f:
    PAR = {r['term']:float(r['value']) for r in csv.DictReader(f)}

def predict(row):
    x = {k:float(row[k]) for k in ['age','waist','height','male','famhx','htn','inactive']}
    if not all(math.isfinite(v) for v in x.values()):
        raise ValueError('All inputs must be observed and finite.')
    if x['age']<20 or x['waist']<=0 or x['height']<=0:
        raise ValueError('Requires age >=20 and positive waist and height in the same unit.')
    if any(x[k] not in (0,1) for k in ['male','famhx','htn','inactive']):
        raise ValueError('Binary indicators must be 0 or 1; unknown is not 0.')
    age=min(x['age'],80)
    features={'a':age/10-5,'w':10*(x['waist']/x['height']-.55),
              **{k:x[k] for k in ['male','famhx','htn','inactive']}}
    eta=PAR['(Intercept)']+sum(PAR[k]*v for k,v in features.items())
    probability=1/(1+math.exp(-eta)) if eta>=0 else math.exp(eta)/(1+math.exp(eta))
    r0=int(probability>=PAR['THRESHOLD'])
    r1=int(probability>=(PAR['THRESHOLD_OLD'] if age>=60 else PAR['THRESHOLD']))
    return {'probability':format(probability,'.17g'),'R0':r0,'R1':r1,'R4':int(age>=60 or r0)}

if __name__=='__main__':
    if len(sys.argv)!=3:
        raise SystemExit(__doc__)
    with open(sys.argv[1],encoding='utf-8-sig',newline='') as f:
        reader=csv.DictReader(f);fields=reader.fieldnames or []
        rows=[dict(r,**predict(r)) for r in reader]
    with open(sys.argv[2],'w',encoding='utf-8-sig',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields+['probability','R0','R1','R4'])
        writer.writeheader();writer.writerows(rows)
