"""Re-render optional production figures from frozen aggregate CSVs, no model fits.
Only standalone A/B/C panel labels and leading panel-letter prefixes are suppressed.
The labelled upload TIFFs remain unchanged.
"""
from pathlib import Path
import os,sys,runpy,shutil,json
C=Path(__file__).resolve().parents[1];R=C;O=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else C/'production_figures';O.mkdir(parents=True,exist_ok=True)
# Install dependencies from requirements.txt before running.
os.environ['MPLCONFIGDIR']=str(O/'_working/mpl')
import matplotlib
matplotlib.use('Agg')
from matplotlib.figure import Figure
from matplotlib.text import Text
from PIL import Image
Q=O/'_working';Q.mkdir(exist_ok=True);P=O
original=Figure.savefig;removed=[]
def save_no_letters(fig,*args,**kwargs):
    for t in fig.findobj(Text):
        s=t.get_text()
        if s in ['A','B','C','D']:removed.append(s);t.set_text('')
        elif len(s)>3 and s[0] in 'ABCD' and s[1:3]=='  ':removed.append(s);t.set_text(s[3:])
    return original(fig,*args,**kwargs)
Figure.savefig=save_no_letters
sys.argv=[str(C/'tools/render_figures.py'),str(Q/'no_letters_base')]
runpy.run_path(sys.argv[0],run_name='__main__')
for number,name in [(1,'Figure1_cohort_flow'),(2,'Figure2_case_score_distribution'),(3,'Figure3_age_differences'),(4,'Figure4_calibration')]:
    source=Q/'no_letters_base/figures'
    for ext in ['png','pdf','svg']:shutil.copy2(source/(name+'.'+ext),P/f'Figure {number}.{ext}')
    im=Image.open(P/f'Figure {number}.png').convert('RGB');im.save(P/f'Figure {number}.tif',compression='tiff_lzw',dpi=(390,390))
sys.argv=[str(C/'analysis/reviewer_response_v11/figure5.py'),str(C),str(C/'analysis/reviewer_response_v11/results'),str(P)]
# The figure writer regenerates only a rendering metadata file in the working copy.
metadata=C/'analysis/reviewer_response_v11/results/Figure5_render_metadata.json';prior=metadata.read_bytes()
runpy.run_path(sys.argv[0],run_name='__main__');metadata.write_bytes(prior)
for number in range(1,6):
    t=P/f'Figure {number}.svg'
    assert not any(x in t.read_text(encoding='utf-8') for x in ['>A</text>','>B</text>','>C</text>'])
metadata={'status':'PASS','purpose':'optional production set, not default review upload','figures':5,'panel_text_suppressed':removed,'source':'Existing plotting code and frozen aggregate CSVs; no model fitting or result recomputation','labelled_main_TIFFs':'unchanged in Figures/'}
(Q/'no_panel_letters_metadata.json').write_text(json.dumps(metadata,indent=2),encoding='utf-8')
(P/'README.txt').write_text('Optional production set without panel letters. Default review uploads are the labelled TIFFs in Figures/. Panel positions: Figure 2 left/right = A/B; Figure 3 left/right = A/B; Figure 4 left/middle/right = A/B/C; Figure 5 left/top-right/bottom-right = A/B/C. Figure 1 is unchanged in content and has no lettered panels. These plots are generated from frozen aggregate datasets; no model was refitted. Do not upload both sets as duplicate figures.\n',encoding='utf-8')
print(json.dumps(metadata,ensure_ascii=False))
