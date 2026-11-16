"""Prepare local v11 references after restoring Data Sheets 2 and 3 together.
python prepare_reference.py RESTORED_ROOT OUTPUT_V11_REFERENCE
Copies a participant-level mapping locally. Do not upload that output directory.
"""
from pathlib import Path
import sys,shutil
V,O=map(Path,sys.argv[1:3]);O.mkdir(parents=True,exist_ok=True)
mapping=V/'analysis/hypertension_budget/participant_mapping.csv'
summary=V/'analysis/hypertension_budget/results/S25_hypertension_strategy_comparisons.csv'
assert mapping.is_file() and summary.is_file(),'Restore v11 Data Sheet 2 and 3 into the same root first'
for src,name in [(mapping,'participant_mapping_PRIVATE.csv'),(summary,summary.name)]:
    target=O/name
    if target.exists(): assert target.read_bytes()==src.read_bytes(),'Refusing to overwrite a different reference file'
    else:shutil.copyfile(src,target)
print('Prepared v11 reference directory locally; participant mapping must remain private to this delivery.')
