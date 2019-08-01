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
    DIR, ANACONDA_PATH, HAPMAP3_SNPS, PLINKFILES, PLINKFILES_EAS,
    BASELINE, BASELINE_EAS, BLANK, BLANK_EAS
)




# Functions ====================================================================

def write_annot(output_prefix, genome, chrom, annotation):
    with gzip.open(
        '{}.{}.{}.annot.gz'.format(output_prefix, annotation, chrom),
        'wb'
    ) as output_annot:
        output_annot.write(
            (
                '\t'.join(
                    genome.variants_header.tuple + ('ANNOT' + '\n',)
                )
                + '\n'.join(
                    '\t'.join(
                        variant.tuple
                        + (str(int(annotation in variant.annotations)),)
                    )
                    for variant in genome.chromosome[chrom].variants
                )
                + '\n'
            ).encode()
        )


def construct_annot(
    output_prefix,
    chrom,
    annotations,
    blank=BLANK,
    processes=1
):
    with funcgenom.Genome() as genome:
        print(f'loading variants on chromosome {chrom}')
        genome.load_variants(os.path.join(blank, f'blank.{chrom}.annot.gz'))
        genome.sort_variants()
        print(f'loading annotations on chromosome {chrom}')
        genome.load_annotations(annotations)
        annotation_set = set(genome.chromosome[chrom].annotations.keys())
        genome.sort_annotations()
        print(f'annotating variants on chromosome {chrom}')
        genome.annotate_variants(processes=processes)
        print(f'writing output on chromosome {chrom}')
        with Pool(processes=processes) as pool:
            pool.map(
                partial(
                    write_annot,
                    output_prefix=output_prefix,
                    genome=genome,
                    chromosome=chrom
                ),
                annotation_set
            )
        return annotations


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
    for chrom in range(1, 23):
        construct_annot(args.output, chrom, args.annotation)
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
    if args.processes > 16:
        raise Exception(
            f'{args.processes} processes, really? Annotating those variants '
            'takes a lot of memory when multiprocessing, and more processes '
            'means more memory consumption. You almost certainly don\'t need '
            'more than 16 processes for this - trust me, it won\'t take THAT '
            'long.'
        )
    return args




# Execute ======================================================================

if __name__ == '__main__':
    main()
