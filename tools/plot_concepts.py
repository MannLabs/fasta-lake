#!/usr/bin/env python3
"""Regenerate the three conceptual tutorial figures in a fresh directory.

Requires matplotlib. The diagrams are authored worked examples, not empirical data.
Uses openly available DejaVu fonts; the original local house-font export is retained
in the distributed SVG/PDF/PNG assets. Figure geometry and labels are reproducible.
"""
import argparse
import hashlib
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--out', type=Path, required=True)
args = parser.parse_args()
OUT = args.out
OUT.mkdir(parents=True, exist_ok=False)
BLUE = '#0072B2'
ORANGE = '#E69F00'
GRAY = '#7f7f7f'
TEAL = '#009E73'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'svg.fonttype': 'none',
                     'pdf.fonttype': 42, 'font.size': 7.2})

def save(fig, name):
    for extension in ('png', 'svg', 'pdf'):
        fig.savefig(OUT / (name + '.' + extension), dpi=240)
    plt.close(fig)

from matplotlib.patches import FancyBboxPatch
def box(ax, x, y, w, h, text, color=BLUE, fs=8):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.015',
                                fc=color,ec='none',alpha=.10))
    ax.text(x+w/2,y+h/2,text,ha='center',va='center',fontsize=fs,color='#333333')
def arrow(ax,a,b,color=GRAY):
    ax.annotate('',xy=b,xytext=a,arrowprops=dict(arrowstyle='->',lw=1.2,color=color))

fig,ax=plt.subplots(figsize=(7.1,4.0))
fig.subplots_adjust(left=.02,right=.98,bottom=.02,top=.98)
ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
for x,text in [(.02,'Reference lake\nexact sequence identity'),(.37,'Cohort evidence lake\npooled containment'),(.72,'Acquisition draw\nown peptide evidence')]:
    box(ax,x,.74,.25,.17,text)
arrow(ax,(.28,.825),(.355,.825));arrow(ax,(.63,.825),(.705,.825))
box(ax,.02,.40,.40,.19,'De novo predictions\nTop 30% of eligible rows: extraction\nAll eligible rows: razor',TEAL,7.5)
arrow(ax,(.28,.60),(.49,.73),TEAL)
arrow(ax,(.42,.51),(.70,.51),TEAL)
box(ax,.72,.40,.25,.19,'Razor database\nexplicit hash-acc ties')
arrow(ax,(.845,.73),(.845,.60))
box(ax,.72,.07,.25,.17,'SAGE\nidentification and LFQ')
arrow(ax,(.845,.39),(.845,.25))
box(ax,.02,.07,.27,.17,'Matched metaG FASTA\ncomparison search',ORANGE,7.5)
box(ax,.37,.07,.26,.17,'Same mzML\nfor the paired searches',ORANGE,7.5)
ax.plot([.155,.155,.845],[.065,.015,.015],color=ORANGE,lw=1.2)
arrow(ax,(.845,.015),(.845,.065),ORANGE)
arrow(ax,(.64,.155),(.705,.155),ORANGE)
ax.text(.5,.965,'Construction and database searching',ha='center',fontsize=10)
save(fig,'workflow_concept')

fig,axs=plt.subplots(2,1,figsize=(7.1,5.1),constrained_layout=True)
for ax in axs: ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
ax=axs[0]
ax.text(.0,.96,'a',weight='bold',fontsize=10)
for y,seq in [(.72,'MPEPTIDEAK'),(.28,'MPEPTLDEAK')]:
    h=hashlib.sha256(seq.encode()).hexdigest()[:12]
    ax.text(.05,y,seq,fontfamily='monospace',fontsize=13,color=BLUE)
    ax.text(.05,y-.14,'SHA-256 prefix: '+h,fontsize=8,color=GRAY)
    arrow(ax,(.42,y),(.65,.5))
ax.text(.79,.52,'MPEPTLDEAK',ha='center',fontfamily='monospace',fontsize=12,color=TEAL)
ax.text(.79,.32,'I/L-normalized\npeptide matching key',ha='center',fontsize=9)
ax=axs[1]
ax.text(.0,.96,'b',weight='bold',fontsize=10)
headers=['Acquisition','Member-set key','Sample key','Study key']
xs=[.04,.29,.57,.81]
for x,label in zip(xs,headers):ax.text(x,.80,label,fontsize=8.5)
for y,vals in [(.61,['1','A;B','A','B']),(.43,['2','B;C','B','B'])]:
    for x,val in zip(xs,vals):ax.text(x,y,val,fontsize=11,color=BLUE if x==.81 else '#333333')
ax.text(.04,.19,'Observed graph: A → p₁; B → p₁, p₂; C → p₂',fontsize=9)
ax.text(.04,.02,'Greedy cover selects B, explaining both observed peptides.',fontsize=9,color=TEAL)
save(fig,'identity_and_grouping')

# A concrete razor example in which shared-only evidence retains C.
fig,ax=plt.subplots(figsize=(7.1,3.6))
fig.subplots_adjust(left=.04,right=.98,bottom=.05,top=.95)
ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
prots=['A','B','C','D'];peps=['uA','uB','s1','s2','s3']
member=[{0},{1},{0,1},{0,2},{2,3}];assigned=[0,1,0,0,2]
for i,p in enumerate(prots):
    ax.text(.30+i*.15,.89,p,ha='center',fontsize=11)
    ax.text(.30+i*.15,.79,f'U={ [1,1,0,0][i] }, T={ [3,2,2,1][i] }',ha='center',fontsize=8)
for j,p in enumerate(peps):
    y=.66-j*.105;ax.text(.10,y,p,ha='center',va='center',fontsize=10)
    for i in member[j]:
        ax.scatter(.30+i*.15,y,s=90,facecolor=BLUE if i==assigned[j] else 'white',
                    edgecolor=BLUE if i==assigned[j] else GRAY,lw=1)
for i in range(4):
    ax.text(.30+i*.15,.08,'retained' if i<3 else 'omitted',ha='center',fontsize=8.5,
            color=BLUE if i<3 else GRAY)
ax.text(.98,.65,'Filled: assigned\nOpen: also contains',ha='right',fontsize=8)
ax.text(.98,.38,'C survives with\nshared evidence only',ha='right',fontsize=8,color=BLUE)
save(fig,'razor_worked_example')
print('FIGURES '+str(OUT),flush=True)
