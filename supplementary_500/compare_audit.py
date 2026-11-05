from pathlib import Path
import csv,json,math,collections
root=Path(__file__).resolve().parent
workspace=root.parents[2]
supplied=workspace/'outputs/diabetes_draft_review_20261104/package/results'

def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))

comparisons=[]
for kind in ['joint_cutoff','access_gap']:
    before=supplied/f'{kind}_v5_20261103'
    after=root/f'{kind}_rerun'
    for path in before.glob('*.csv'):
        new=after/path.name
        if not new.exists():
            comparisons.append({'analysis':kind,'file':path.name,'status':'missing_rerun'});continue
        a,b=read(path),read(new)
        report={'analysis':kind,'file':path.name,'supplied_rows':len(a),'rerun_rows':len(b),'max_absolute_numeric_difference':0.,'text_or_missing_mismatches':0,'numeric_cells_compared':0,'schema_equal':bool(a and b and list(a[0])==list(b[0]))}
        if len(a)!=len(b) or not report['schema_equal']:
            report['status']='shape_mismatch';comparisons.append(report);continue
        for x,y in zip(a,b):
            for k,v in x.items():
                w=y[k]
                try:
                    vf,wf=float(v),float(w)
                    if math.isfinite(vf) and math.isfinite(wf):
                        report['numeric_cells_compared']+=1
                        report['max_absolute_numeric_difference']=max(report['max_absolute_numeric_difference'],abs(vf-wf))
                    elif v!=w:report['text_or_missing_mismatches']+=1
                except (ValueError,TypeError):
                    if v!=w:report['text_or_missing_mismatches']+=1
        report['status']='match_at_1e-8' if report['text_or_missing_mismatches']==0 and report['max_absolute_numeric_difference']<=1e-8 else 'different'
        comparisons.append(report)

rep=root/'access_gap_rerun'
status=read(rep/'replicate_status.csv')
counts=read(rep/'replicate_valid_counts.csv')
mult=read(rep/'replicate_psu_multipliers.csv')
groups=collections.defaultdict(list)
for row in mult:
    groups[(int(row['replicate']),row['stratum_psu'].split()[0])].append(float(row['multiplier']))
bad=[]
for k,values in groups.items():
    nh=len(values)
    if min(values)<0 or abs(sum(values)-nh)>1e-10 or (nh>1 and any(abs(v*(nh-1)/nh-round(v*(nh-1)/nh))>1e-10 for v in values)):bad.append(k)
archive_rows=0
archive_bad_labels=0
archive_replicates=set()
with (rep/'replicate_values.csv').open(encoding='utf-8-sig',newline='') as f:
    for row in csv.DictReader(f):
        archive_rows+=1
        try:
            b=int(row['replicate']);budget=float(row['budget']);index=int(row['budget_index'])
            archive_replicates.add(b)
            if not (1<=b<=500 and 1<=index<=51 and abs(budget-(.25+.005*(index-1)))<1e-10):archive_bad_labels+=1
        except (ValueError,TypeError):archive_bad_labels+=1
summary={'scope':'Original fixed-model conditional bootstrap and original fixed-threshold Taylor analyses. No fitting/selection uncertainty added. No claim of historical blindness.',
 'comparisons':comparisons,
 'replicate_audit':{'requested':500,'completed':len(status),'failed_status_rows':sum(r['status']!='completed' for r in status),'min_valid_across_curve_grid':min(int(r['valid_replicates']) for r in counts),'max_valid_across_curve_grid':max(int(r['valid_replicates']) for r in counts),'cells_with_under_500_valid':sum(int(r['valid_replicates'])<500 for r in counts),'curve_grid_cells':len(counts),'psu_multiplier_rows':len(mult),'replicate_stratum_groups':len(groups),'invalid_multiplier_groups':bad,'replicate_value_rows':archive_rows,'replicate_ids_count':len(archive_replicates),'invalid_archive_labels':archive_bad_labels,'archived':'RDS stores exact 3D replicate array, point estimates, budgets and each pre-draw RNG state; CSV stores each replicate value and PSU multiplier.'},
 'joint_bootstrap':'Not applicable: the original joint_cutoff_v5.R contains no 500-replicate procedure. Its original analyses were rerun without adding one.'}
(root/'verification.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(summary,ensure_ascii=False,indent=2))
