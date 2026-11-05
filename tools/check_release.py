"""Offline integrity, prediction and archived interval checks; no model refit."""
from pathlib import Path
import ast, csv, hashlib, importlib.util, json, math, sys
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
sys.stdout.reconfigure(encoding='utf-8')
def read_csv(path):
    with path.open(encoding='utf-8-sig',newline='') as f: return list(csv.DictReader(f))

manifest=read_csv(ROOT/'SHA256SUMS.csv')
for row in manifest:
    p=ROOT/row['path']
    assert hashlib.sha256(p.read_bytes()).hexdigest()==row['sha256'],row['path']
scripts=list(ROOT.rglob('*.py'))
for p in scripts: ast.parse(p.read_text(encoding='utf-8-sig'),filename=str(p))
provenance=json.loads((ROOT/'SOURCE_PROVENANCE.json').read_text(encoding='utf-8'))
for row in provenance['files']:
    assert row['unchanged'] and hashlib.sha256((ROOT/row['path']).read_bytes()).hexdigest()==row['sha256']

spec=importlib.util.spec_from_file_location('locked',ROOT/'reproducibility/predict_locked.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
synthetic=read_csv(ROOT/'examples/people_synthetic.csv')
pred=[module.predict(r) for r in synthetic]
assert pred[2]==pred[3], 'Age 80 top-code changed'
assert all(p['R4']>=p['R1']>=p['R0'] for p in pred)
for r,p in zip(synthetic,pred):
    a=min(float(r['age']),80)/10-5
    w=10*(float(r['waist'])/float(r['height'])-.55)
    eta=-5.18624341294669+.326735063668154*a+.961234496119346*w
    eta+=.648480842884934*float(r['male'])+.564755920095949*float(r['famhx'])
    eta+=.0220416394625756*float(r['htn'])+.186399703875562*float(r['inactive'])
    assert abs(float(p['probability'])-1/(1+math.exp(-eta)))<1e-15
bad=synthetic[0].copy();bad['famhx']='nan'
try: module.predict(bad)
except ValueError: pass
else: raise AssertionError('Unknown family history accepted')

expected=pd.read_csv(ROOT/'reproducibility/main_archive/interval_summary.csv')
point=pd.read_csv(ROOT/'reproducibility/main_archive/point_metrics.csv').set_index('key').estimate
max_error=0.;cells=0
for analysis,fname in [('full_pipeline','bootstrap_full_metrics.csv'),('fixed_model_bootstrap','bootstrap_fixed_metrics.csv')]:
    z=pd.read_csv(ROOT/'reproducibility/main_archive'/fname)
    assert len(z)==1000 and list(z.replicate)==list(range(1,1001))
    for row in expected[expected.analysis.eq(analysis)].itertuples():
        key='::'.join([row.strategy,row.group,row.measure])
        x=z[key].to_numpy();estimate=float(point[key]);se=np.std(x,ddof=1)
        assert np.isfinite(x).all() and row.valid_B==1000
        pairs=[(estimate,row.estimate),(se,row.bootstrap_se),
               (estimate-1.96*se,row.lower95),(estimate+1.96*se,row.upper95),
               (np.quantile(x,.025),row.percentile_lower95),
               (np.quantile(x,.975),row.percentile_upper95),
               (float(x.mean()),row.bootstrap_mean)]
        max_error=max(max_error,*(abs(a-b) for a,b in pairs));cells+=len(pairs)
assert max_error<1e-10,max_error

raw=read_csv(ROOT/'reproducibility/raw_download_manifest.csv')
assert len(raw)==63 and len({r['file'] for r in raw})==63
assert all(r['CDC_URL'].startswith('https://wwwn.cdc.gov/') and len(r['sha256'])==64 for r in raw)
report=dict(status='PASS',manifest_files=len(manifest),python_scripts_parsed=len(scripts),
            unchanged_archived_files=len(provenance['files']),synthetic_people=len(synthetic),
            archived_interval_rows=len(expected),numeric_comparisons=cells,
            maximum_interval_replay_error=max_error,raw_manifest_files=len(raw),
            r_models_refitted=False,raw_files_downloaded=False)
print(json.dumps(report,ensure_ascii=False,indent=2))
