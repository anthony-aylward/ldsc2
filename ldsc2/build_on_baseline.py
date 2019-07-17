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

from multiprocessing import Pool

from ldsc2.env import DIR, HAPMAP3_SNPS, PLINKFILES, PLINKFILES_EAS



# Functions ====================================================================

def construct_annot(args, chromosome):
    with funcgenom.Genome() as genome:
        print('loading variants on chromosome {}'.format(chromosome))
        genome.load_variants(
            '{}.{}.annot.gz'.format(args.blank, chromosome)
        )
        genome.sort_variants()
        print('loading annotations on chromosome {}'.format(chromosome))
        genome.load_annotations(args.annotations)
        annotations = set(genome.chromosome[chromosome].annotations.keys())
        genome.sort_annotations()
        print('annotating variants on chromosome {}'.format(chromosome))
        genome.annotate_variants(processes=args.processes)
        print('writing output on chromosome {}'.format(chromosome))
        with Pool(processes=args.processes) as pool:
            pool.starmap(
                write_annot,
                ((args, genome, ann, chromosome) for ann in annotations)
            )
        return annotations


def write_annot(args, genome, annotation, chromosome):
    with gzip.open(
        '{}.{}.{}.annot.gz'.format(args.output, annotation, chromosome),
        'w'
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
                    for variant in genome.chromosome[chromosome].variants
                )
                + '\n'
            ).encode()
        )


def ldsc(args, annotation, chromosome):
    subprocess.call(
        (
            ('/home/data/ldsc/ldsc.py',)
        )
        + (
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
        '--skip-to-chr',
        metavar='<int>',
        type=int,
        default=1,
        help='skip to this chromosome [1]'
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
