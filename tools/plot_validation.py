#!/usr/bin/env python3
"""Regenerate the empirical figure from the distributed Source Data tables."""
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, default=Path(__file__).resolve().parents[1]/'docs/concepts/data')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    blue, orange, gray, green = '#0072B2', '#E69F00', '#7f7f7f', '#009E73'
    plt.rcParams.update({'font.size': 7.2, 'axes.spines.top': False,
                         'axes.spines.right': False, 'svg.fonttype': 'none', 'pdf.fonttype': 42})
    fig, axes = plt.subplots(2, 2, figsize=(7.1,5.4), constrained_layout=True)
    ax = axes[0,0]
    for row, color in zip(rows(args.data/'construction_counts.tsv'), [blue,orange]):
        ax.plot(range(4), [float(row[k]) for k in ('reference','evidence','median_draw','median_razor')],
                'o-', color=color, label=row['cohort'], lw=1, ms=4)
    ax.set_yscale('log'); ax.set_xticks(range(4), ['Reference','Evidence','Sample draw','Razor'])
    ax.set_ylabel('Sequences'); ax.set_title('a  Database reduction', loc='left'); ax.legend(frameon=False)
    ax = axes[0,1]
    pairs = rows(args.data/'hstool_pairwise_jaccard.tsv')
    for x,key,color in [(0,'db_jaccard',gray),(1,'peps_jaccard',blue)]:
        values = np.array([float(row[key]) for row in pairs])
        ax.scatter(x+np.linspace(-.14,.14,len(values)),values,s=10,color=color,alpha=.6)
        ax.plot([x-.23,x+.23],[np.median(values)]*2,color='black',lw=1)
    ax.set_xticks([0,1],['Razor databases','Identified peptides']);ax.set_ylim(0,1)
    ax.set_ylabel('Pairwise Jaccard similarity');ax.set_title('b  HSTOOL technical acquisitions',loc='left')
    ax.text(.03,.94,'45 dependent pairs from 10 acquisitions',transform=ax.transAxes,va='top',fontsize=6.5)
    ax = axes[1,0]
    for row in rows(args.data/'predict_paired_comparison.tsv'):
        values=[int(row['metag_peptides']),int(row['fastalake_peptides'])]
        ax.plot([0,1],values,color=gray,alpha=.6,lw=1)
        ax.scatter([0,1],values,c=[orange,blue],s=15,zorder=3)
        number=row['sample'].split('MicrobPred_')[1].split('_')[0]
        if number in ('10','19'):
            ax.annotate(number,(1,values[1]),xytext=(5,0),textcoords='offset points',fontsize=7)
    ax.set_yscale('log');ax.set_ylim(150,85000);ax.set_xlim(-.15,1.2)
    ax.set_xticks([0,1],['Matched metaG','FastaLake'])
    ax.set_ylabel('Identified canonical peptides');ax.set_title('c  PREDICT paired searches',loc='left')
    ax.text(.03,.94,'One line per acquisition; peptide-q threshold 0.01',transform=ax.transAxes,va='top',fontsize=6.5)
    ax=axes[1,1]
    cv=rows(args.data/'hstool_normalization_sensitivity.tsv')
    for field,label,color in [('raw_cv','Raw LFQ',gray),('median_ratio_normalized_cv','After median-ratio scaling',green)]:
        values=np.sort([float(row[field]) for row in cv])
        ax.plot(values,np.arange(1,len(values)+1)/len(values),color=color,label=label,lw=1.2)
    ax.set_xlabel('Coefficient of variation');ax.set_ylabel('Fraction of complete-case peptides')
    ax.set_title('d  HSTOOL intensity repeatability',loc='left');ax.legend(frameon=False,loc='lower right',fontsize=6.5)
    ax.text(.98,.55,f'{len(cv):,} complete-case peptides\n10 acquisitions; no imputation',transform=ax.transAxes,
            ha='right',va='top',fontsize=6.5)
    for suffix in ('png','svg','pdf'):
        fig.savefig(args.out/('empirical_acceptance.'+suffix),dpi=300,bbox_inches='tight')


if __name__ == '__main__':
    main()
