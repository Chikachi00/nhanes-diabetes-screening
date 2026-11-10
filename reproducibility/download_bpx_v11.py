"""Download and hash-check the six public BPX files added for v11.
Usage: python download_bpx_v11.py RAW_DIRECTORY
"""
from pathlib import Path
import csv,hashlib,sys,urllib.request

def main():
    if len(sys.argv)!=2:raise SystemExit(__doc__)
    dest=Path(sys.argv[1]);dest.mkdir(parents=True,exist_ok=True)
    with (Path(__file__).parent/'additional_raw_manifest_v11.csv').open(encoding='utf-8-sig') as f:
        records=list(csv.DictReader(f))
    for r in records:
        p=dest/r['file']
        if p.exists() and hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256']:
            print('VERIFIED',p.name);continue
        req=urllib.request.Request(r['CDC_URL'],headers={'User-Agent':'NHANES-research-reproduction/1.0'})
        with urllib.request.urlopen(req,timeout=60) as response:b=response.read()
        if hashlib.sha256(b).hexdigest()!=r['sha256']:raise SystemExit('Current CDC file differs from archived hash: '+p.name)
        if not b.startswith(b'HEADER RECORD*******LIBRARY HEADER RECORD'):raise SystemExit('Not an XPT file: '+p.name)
        p.write_bytes(b);print('DOWNLOADED AND VERIFIED',p.name)

if __name__=='__main__':main()
