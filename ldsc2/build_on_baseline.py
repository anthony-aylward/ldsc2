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
import sys

from collections import namedtuple
from functools import partial
from multiprocessing import Pool
from pybedtools import BedTool

from ldsc2.env import (
    DIR, ANACONDA_PATH, HAPMAP3_SNPS, PLINKFILES, PLINKFILES_EAS
)

sys.path.append(DIR)
import make_annot



# Functions ====================================================================

def _make_annot_files(bed_file, bim_file, annot_file):
    Args = namedtuple('Args', ('annot_file', 'bimfile'))
    make_annot.make_annot_files(
        Args(annot_file=annot_file, bimfile=bim_file),
        BedTool(bed_file).sort().merge()
    )


def make_annot_files(bed_file, bim_prefix, annot_prefix, processes=1):
    def make_annot_chrom(chrom):
        _make_annot_files(
            bed_file,
            f'{bim_prefix}.{chrom}.bim',
            f'{annot_prefix}.{chrom}.annot.gz'
        )

    with Pool(processes=processes) as pool:
        pool.map(make_annot_chrom, range(1, 23))



def ldsc(args, annotation, chromosome):
    subprocess.run(
        (
            ANACONDA_PATH,
            os.path.join(DIR, 'ldsc.py'),
            '--l2',
            '--bfile', '{}.{}'.format(args.plink_prefix, chromosome),
            '--ld-wind-cm', '1',
            '--annot', '{}.{}.{}.annot.gz'.format(
                args.output,
                annotation,
                chromosome
            ),
            '--out', '{}.{}.{}'.format(
                args.output,
                annotation,
                chromosome
            ),
            '--print-snps', '{}.{}.snp'.format(args.snp_prefix, chromosome)
        )
    )


def main(args):
    """main loop"""
    
    args = parse_arguments()
    for chromosome in range(args.skip_to_chr, 23):
        annotations = construct_annot(args, chromosome)
        for annotation in annotations:
            ldsc(args, annotation, chromosome)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            'Generate a set of annotation-specific ld-score files for use with '
            'the baseline model from Finucane et al. 2015'
        )
    )
    parser.add_argument(
        'blank',
        metavar='<prefix/for/blank.annot.gz>',
        help='prefix of blank .annot.gz files for input'
    )
    parser.add_argument(
        'annotations',
        metavar='<path/to/annotations.bed>',
        help='path to .bed file of annotations'
    )
    parser.add_argument(
        'output',
        metavar='<prefix/for/output/files>',
        help='prefix for output files'
    )
    parser.add_argument(
        '--plink-prefix',
        metavar='<prefix/for/plink/files>',
        default=os.path.join(PLINKFILES, '1000G.EUR.QC'),
        help='prefix of plink files for input'
    )
    parser.add_argument(
        '--snp-prefix',
        metavar='<prefix/for/snp/files>',
        default=os.path.join(PLINKFILES, '1000G_Phase3_plinkfiles'),
        help='prefix of snp files for input'
    )
    parser.add_argument(
        '--processes',
        metavar='<int>',
        type=int,
        default=1,
        help='number of processes [1]'
    )
    args = parser.parse_args()
    if args.processes > 16:
        raise Exception(
            '{} processes, really? Annotating those variants takes a lot of '
            'memory when multiprocessing, and more processes means more memory '
            'consumption. You almost certainly don\'t need more than 16 '
            'processes for this - trust me, it won\'t take THAT long.'
            .format(args.processes)
        )
    return args




# Execute ======================================================================

if __name__ == '__main__':
    main()
