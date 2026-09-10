"""Rebuild review redlines with table placeholders and explicit table treatment.
Usage: python build_tracked.py ORIGINAL_BASELINE CURRENT_CANDIDATE OUTPUT_STEM [--cumulative]
Baseline/candidate can be a revision-layout root (05-writing/,03-theory/) or
an archived flat 05-writing baseline with appendix_new.tex in that directory.
"""
from pathlib import Path
import re,sys,subprocess
old=Path(sys.argv[1]).resolve();new=Path(sys.argv[2]).resolve();stem=Path(sys.argv[3]).resolve();cum='--cumulative' in sys.argv

def wroot(p):return p/'05-writing' if (p/'05-writing').is_dir() else p
ow,nw=wroot(old),wroot(new)
def locate(root,name):
 if name.endswith('appendix_new.tex'):
  for q in [root.parent/'03-theory/appendix_new.tex',root/'appendix_new.tex']:
   if q.exists():return q
 q=root/name
 if not q.suffix:q=q.with_suffix('.tex')
 return q
wrapper=(nw/'main_draft.tex').read_text()
def flatten(root):
 def expand(s):
  def inp(m):
   name=m[1]
   if name.startswith('tables/'):
    return '\\TABLEPLACEHOLDER'+Path(name).stem.replace('_','Q')+'ENDPLACEHOLDER'
   p=locate(root,name)
   if not p.exists():raise FileNotFoundError(p)
   return expand(p.read_text())
  return re.sub(r'\\input\{([^}]+)\}',inp,s)
 return expand(wrapper)
stem.parent.mkdir(parents=True,exist_ok=True);oldflat=stem.with_name(stem.name+'_old_flat.tex');newflat=stem.with_name(stem.name+'_new_flat.tex')
oldflat.write_text(flatten(ow));newflat.write_text(flatten(nw))
r=subprocess.run(['latexdiff','--math-markup=whole','--disable-citation-markup',str(oldflat),str(newflat)],check=True,capture_output=True,text=True)
s=r.stdout
stem.with_name(stem.name+'_latexdiff.log').write_text(r.stderr)
def table(m):
 name=m[1].replace('Q','_');cur=(nw/'tables'/f'{name}.tex').read_text()
 if cum and name in ['tab_gauss','tab_map']:
  prior=(ow/'tables'/f'{name}.tex').read_text()
  prior=re.sub(r'\\label\{[^}]+\}','',prior).replace('\\caption{','\\caption*{Before revision: ')
  prior=re.sub(r'(\\begin\{table\}(?:\[[^]]*\])?)',lambda m:m[0]+'\\color{red}',prior)
  prior=prior.replace('\\begin{tabular}','\\begin{adjustbox}{max width=\\textwidth,max totalheight=0.76\\textheight,keepaspectratio}\n\\begin{tabular}').replace('\\end{tabular}','\\end{tabular}\n\\end{adjustbox}')
  return '\\clearpage\n{\\color{red}\n'+prior+'\n\\clearpage}\n{\\color{blue}\n'+cur+'\n\\clearpage}\n'
 if cum and name in ['tab_cahb_transfer','tab_boundary_detail']:
  cur=re.sub(r'(\\begin\{table\}(?:\[[^]]*\])?)',lambda m:m[0]+'\\color{blue}',cur)
  return '{\\color{blue}\n'+cur+'\n}\n'
 return cur
s=re.sub(r'\\TABLEPLACEHOLDER([A-Za-z0-9]+)ENDPLACEHOLDER',table,s)
legend=r'''\begin{center}\Large Colored tracked changes\end{center}
\noindent\textcolor{blue}{Blue underlined text shows additions.}
\textcolor{red}{\sout{Red struck text shows deletions.}}
'''
legend+=('The baseline is the original pre-Codex manuscript preserved on September 10, 2026. Main Tables 1--2 are shown in full before (red) and after (blue); new appendix tables are blue. Other tables retain current layout and unchanged numerical values.\n' if cum else 'The baseline is the clean manuscript immediately before this revision round. Numerical table entries are unchanged.\n')
legend+='Current typography and bibliography are used. Relocations can appear as both deletions and additions.\n\\clearpage\n'
s=s.replace('\\begin{document}','\\hypersetup{hidelinks}\n\\emergencystretch=3em\n\\begin{document}\n\\sloppy\n'+legend,1)
stem.with_suffix('.tex').write_text('\n'.join(line.rstrip() for line in s.splitlines())+'\n')
print(stem.with_suffix('.tex'))
