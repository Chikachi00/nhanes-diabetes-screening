"""Download the manifest's public NHANES XPT files and check exact hashes.

Usage: python download_nhanes.py OUTPUT_DIRECTORY
Existing matching files are reused. A mismatch is retained with .unmatched
extension and is never accepted as the original analysis input.
"""
from pathlib import Path
import csv, hashlib, sys, urllib.request

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()

def main():
    if len(sys.argv)!=2:raise SystemExit(__doc__)
    out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=True)
    manifest=Path(__file__).resolve().parent/'raw_download_manifest.csv'
    with manifest.open(encoding='utf-8-sig') as f:records=list(csv.DictReader(f))
    for r in records:
        target=out/r['file']
        if target.exists() and sha(target)==r['sha256']:
            print('REUSE',r['file']);continue
        temporary=out/(r['file']+'.download')
        request=urllib.request.Request(r['CDC_URL'],headers={'User-Agent':'NHANES-research-reproduction/1.0'})
        with urllib.request.urlopen(request,timeout=90) as response,temporary.open('wb') as f:
            for chunk in iter(lambda:response.read(1024*1024),b''):f.write(chunk)
        if sha(temporary)!=r['sha256']:
            temporary.replace(out/(r['file']+'.unmatched'))
            raise SystemExit('Hash mismatch for '+r['file']+'; inspect the current CDC release before using it.')
        temporary.replace(target);print('VERIFIED',r['file'])

if __name__=='__main__':main()
