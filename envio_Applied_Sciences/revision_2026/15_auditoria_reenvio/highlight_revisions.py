"""Highlight inserted/replaced words on a layout-identical copy of the clean PDF.

Comparison basis is the preserved local pre-review Applied Sciences PDF. This
does not pretend to show deleted words or to reconstruct journal upload history.
"""
import json,re,hashlib,difflib,subprocess
from pathlib import Path
from collections import defaultdict
import pdfplumber
from pypdf import PdfReader,PdfWriter
from pypdf.annotations import Highlight
from pypdf.generic import ArrayObject,FloatObject,NameObject
ROOT=Path(__file__).resolve().parents[1];OUT=Path(__file__).parent
OLD=ROOT.parent/'fuentes/manuscrito_applied_sciences.pdf'
NEW=ROOT/'08_manuscript/manuscrito_applied_sciences.pdf'
MARKED=NEW.with_name('manuscrito_applied_sciences_CAMBIOS.pdf')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def words(path):
    result=[]
    with pdfplumber.open(path) as pdf:
        for page_no,page in enumerate(pdf.pages):
            for w in page.extract_words(x_tolerance=2,y_tolerance=2,extra_attrs=['size']):
                if w['size']<7.4 or w['top']<48 or w['bottom']>page.height-40:continue
                w.update(page=page_no,height=page.height);result.append(w)
    return result
def norm(w):return re.sub(r'[^\w]+','',w['text']).lower()
def main():
    old=words(OLD);new=words(NEW)
    matcher=difflib.SequenceMatcher(None,[norm(w) for w in old],[norm(w) for w in new],autojunk=False)
    changes=set()
    for tag,i,j,a,b in matcher.get_opcodes():
        if tag in ['insert','replace']:changes.update(range(a,b))
    grouped=defaultdict(list)
    for i in sorted(changes):
        w=new[i];grouped[(w['page'],round(w['top'],1))].append(w)
    writer=PdfWriter(clone_from=NEW);annotations=0
    for (page,top),ww in grouped.items():
        ww.sort(key=lambda w:w['x0']);groups=[]
        for w in ww:
            if not groups or w['x0']-groups[-1][-1]['x1']>10:groups.append([w])
            else:groups[-1].append(w)
        for group in groups:
            x0=min(w['x0'] for w in group);x1=max(w['x1'] for w in group)
            y0=group[0]['height']-max(w['bottom'] for w in group);y1=group[0]['height']-min(w['top'] for w in group)
            quad=ArrayObject([FloatObject(v) for v in [x0,y1,x1,y1,x0,y0,x1,y0]])
            annotation=Highlight(rect=(x0,y0,x1,y1),quad_points=quad,highlight_color='fff2a8',printing=True)
            annotation[NameObject('/CA')]=FloatObject(.4)
            writer.add_annotation(page_number=page,annotation=annotation);annotations+=1
    writer.write(MARKED)
    clean=PdfReader(NEW);marked=PdfReader(MARKED)
    assert len(clean.pages)==len(marked.pages)
    for a,b in zip(clean.pages,marked.pages):assert a.extract_text()==b.extract_text()
    record={'original_basis':str(OLD),'original_sha256':sha(OLD),'clean_sha256':sha(NEW),
      'marked_pdf':str(MARKED),'marked_sha256':sha(MARKED),'pages':len(marked.pages),
      'annotations':annotations,'changed_word_tokens':len(changes),'layout_and_text_identical_to_clean':True,
      'method':'Word-token sequence comparison; yellow marks inserted/replaced words. Deletes not displayed. Hyphenation and equation extraction can produce extra highlights.'}
    (OUT/'highlight_report.json').write_text(json.dumps(record,indent=2))
    dest=ROOT/'13_cierre_editorial/pdf_qa/manuscrito_applied_sciences_CAMBIOS';dest.mkdir(exist_ok=True)
    poppler=Path.home()/'.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/Library/bin/pdftoppm.exe'
    subprocess.run([str(poppler),'-r','110','-png',str(MARKED),str(dest/'page')],check=True)
    print(json.dumps(record,indent=2))
if __name__=='__main__':main()
