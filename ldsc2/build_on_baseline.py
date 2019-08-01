#===============================================================================
# build_on_baseline.py
#===============================================================================

"""Generate a set of annotation-specific ld-score files for use with the
baseline model from Finucane et al. 2015
"""




# Imports ======================================================================

import argparse
import funcgenom
import gzip
import os.path
import subprocess

from functools import partial
from multiprocessing import Pool

from ldsc2.env import (
    DIR, ANACONDA_PATH, HAPMAP3_SNPS, PLINKFILES, PLINKFILES_EAS
)




# Functions ====================================================================

def make_annot_file(bed_file, bim_file, annot_file):
    subprocess.run(
        (
            ANACONDA_PATH, os.path.join(DIR, 'ldsc', 'make_annot.py'),
            '--bed-file', bed_file,
            '--bimfile', bim_file,
            '--annot-file', annot_file
        )
    )


def make_annot_file_chrom(chrom, bed_file, bim_prefix, annot_prefix):
    make_annot_file(
        bed_file,
        f'{bim_prefix}.{chrom}.bim',
        f'{annot_prefix}.{chrom}.annot.gz'
    )


def make_annot_files(bed_file, bim_prefix, annot_prefix, processes=1):
    with Pool(processes=processes) as pool:
        pool.map(
            partial(
                make_annot_file_chrom,
                bed_file=bed_file,
                bim_prefix=bim_prefix,
                annot_prefix=annot_prefix
            ),
            range(1, 23)
        )


def compute_ld_scores_chrom(
    chrom,
    annotation,
    output_prefix,
    plink_prefix=os.path.join(PLINKFILES, '1000G.EUR.QC')
):
    subprocess.run(
        (
            ANACONDA_PATH, os.path.join(DIR, 'ldsc', 'ldsc.py'),
            '--l2',
            '--bfile', f'{plink_prefix}.{chrom}',
            '--ld-wind-cm', '1',
            '--annot', f'{output_prefix}.{chrom}.annot.gz',
            '--thin-annot',
            '--out', f'{output_prefix}.{chrom}',
            '--print-snps', f"{os.path.join(HAPMAP3_SNPS, 'hm')}.{chrom}.snp"
        )
    )


def main():
    """main loop"""
    args = parse_arguments()
    make_annot_files(
        bed_file=args.annotation,
        bim_prefix=args.plink_prefix,
        annot_prefix=args.output,
        processes=args.processes
    )
    for chrom in range(1, 23):
        compute_ld_scores_chrom(
            chrom,
            args.annotation,
            args.output
        )


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            'Generate a set of annotation-specific ld-score files for use with '
            'the baseline model from Finucane et al. 2015'
        )
    )
    parser.add_argument(
        'annotation',
        metavar='<path/to/annotations.bed>',
        help='path to .bed file of annotation'
    )
    parser.add_argument(
        'output',
        metavar='<prefix/for/output/files>',
        help='prefix for output files'
    )
    parser.add_argument(
        '--pop',
        choices=('EUR', 'EAS'),
        default='EUR',
        help='reference population (determines default plink files)'
    )
    parser.add_argument(
        '--plink-prefix',
        metavar='<prefix/for/plink/files>',
        help='prefix of user-provided plink files'
    )
    parser.add_argument(
        '--processes',
        metavar='<int>',
        type=int,
        default=1,
        help='number of processes (default: 1)'
    )
    args = parser.parse_args()
    if not args.plink_prefix:
        args.plink_prefix={
            'EUR': os.path.join(PLINKFILES, '1000G.EUR.QC'),
            'EAS': os.path.join(PLINKFILES_EAS, '1000G.EAS.QC')
        }[args.pop]
    return args




# Execute ======================================================================

if __name__ == '__main__':
    main()
