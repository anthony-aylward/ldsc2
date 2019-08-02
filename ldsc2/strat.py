#===============================================================================
# strat.py
#===============================================================================

"""Run stratified LD-score regression"""




# Imports ======================================================================

import argparse
import itertools
import math
import operator
import os.path
import socket
import subprocess
import sys

from multiprocessing import Pool




# CONSTANTS ====================================================================

HOSTNAME = socket.gethostname()

LD_SCORES_DIR = ((HOSTNAME == 'gatsby.ucsd.edu') * '/home') + '/data/ld-scores'




# Functions ====================================================================

# Low-level functions ----------------------------------------------------------

def munge_sumstats(
    ldsc_dir,
    output_prefix,
    sumstats_file_path,
    snp=None,
    sample_size=None,
    a1=None,
    a2=None,
    maf_min=0.01
):
    subprocess.call(
        ('/home/aaylward/anaconda2/bin/python',) * (HOSTNAME == 'holden')
        + (
            os.path.join(ldsc_dir, 'munge_sumstats.py'),
            '--sumstats', sumstats_file_path,
            '--out', output_prefix,
            '--merge-alleles', os.path.join(LD_SCORES_DIR, 'w_hm3.snplist'),
            '--maf-min', str(maf_min)
        )
        + bool(snp) * ('--snp', str(snp))
        + bool(sample_size) * ('--N', str(sample_size))
        + bool(a1) * bool(a2) * ('--a1', a1, '--a2', a2)
    )


def parse_annotation_list(file_path):
    return tuple(annotation.rstrip() for annotation in open(file_path, 'r'))


def stratified_regression(
    output_prefix,
    sumstats_file_path,
    ld_scores_prefix,
    *annotations,
    intercept=True
):
    subprocess.call(
        ('/home/aaylward/anaconda2/bin/python',) * (HOSTNAME == 'holden')
        + (
            os.path.join(ldsc_dir, 'ldsc.py'),
            '--h2', sumstats_file_path,
            '--w-ld-chr', os.path.join(
                LD_SCORES_DIR,
                'LDSCORE/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.'
            ),
            '--ref-ld-chr', (
                '{},{}/LDSCORE/1000G_EUR_Phase3_baseline/'
                'baseline.'
                .format(
                    '.,'.join(
                        '{}.{}.'.format(ld_scores_prefix, annotation)
                        for annotation in annotations
                    ),
                    LD_SCORES_DIR
                )
            ),
            '--overlap-annot',
            '--frqfile-chr', os.path.join(
                LD_SCORES_DIR,
                'LDSCORE/1000G_Phase3_frq/1000G.EUR.QC.'
            ),
            '--out', output_prefix,
            '--print-coefficients'
        )
        + (not intercept) * ('--no-intercept',)
    )


def generate_statistics(output_prefix, annotations):
    for annotation in annotations:
        with open('{}-{}.results'.format(output_prefix, annotation), 'r') as f:
            coefficient, standard_error, z_score = (
                float(val) for val in f.read().splitlines()[1].split()[-3:]
            )
        yield annotation, coefficient, standard_error, z_score




# Mid-level function definitions -----------------------------------------------

def handle_munging(args):
    """munge input files that need to be munged, and return a list of all
    munged files that will be operated on
    """
    
    if args.munge:
        munge_sumstats_indicator = tuple(
            flag in {'-m', '--munge'}
            for flag in sys.argv
            if flag in {'-s', '--sumstats', '-m', '--munge'}
        )
        output_prefixes_for_munging = tuple(
            itertools.compress(args.output_prefix, munge_sumstats_indicator)
        )
        with Pool(processes=args.processes) as pool:
            pool.starmap(
                munge_sumstats,
                (
                    (
                        args.ldsc_dir,
                        output_prefix,
                        file_to_munge,
                        snp,
                        sample_size,
                        a1,
                        a2,
                        args.maf_min
                    )
                    for file_to_munge, output_prefix, snp, sample_size, a1, a2
                    in itertools.zip_longest(
                        args.munge,
                        output_prefixes_for_munging,
                        args.snp if args.snp else (),
                        args.sample_size if args.sample_size else (),
                        args.a1 if args.a1 else (),
                        args.a2 if args.a2 else ()
                    )
                )
            )
        freshly_munged_files = (
            '{}.sumstats.gz'.format(prefix)
            for prefix in output_prefixes_for_munging
        )
        files_already_munged = (file_path for file_path in args.sumstats)
        return tuple(
            (
                next(freshly_munged_files)
                if boolean
                else next(files_already_munged)
            )
            for boolean in munge_sumstats_indicator
        )
    else:
        return args.sumstats


def perform_stratified_analyses(
    output_prefixes,
    file_paths_ordered_and_munged,
    annotations,
    ld_scores_prefix,
    intercept=True
):
    """perform stratified analyses"""
    
    for output_prefix, sumstats_file_path, annotation in (
        (output_input_pair[0], output_input_pair[1], annotation)
        for output_input_pair, annotation
        in itertools.product(
            zip(output_prefixes, file_paths_ordered_and_munged),
            annotations
        )
    ):
        stratified_regression(
            '{}-{}'.format(output_prefix, annotation),
            sumstats_file_path,
            ld_scores_prefix,
            annotation,
            intercept=intercept
        )


def tabulate_statistics(output_prefix, annotations):
    """Tabulte statistics from the individual annotation results"""

    with open('{}.txt'.format(output_prefix), 'w') as f:
        f.write(
            '\n'.join(
                '\t'.join(str(field) for field in line) for line in (
                    [('annotation', 'coefficient', 'standard_error', 'z_score')]
                    + sorted(
                        generate_statistics(output_prefix, annotations),
                        key=operator.itemgetter(3),
                        reverse=True
                    )
                )
            )
            + '\n'
        )




# High-level functions ---------------------------------------------------------

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Run stratified LD-score regression'
    )
    parser.add_argument(
        '--processes',
        type=int,
        default=1,
        choices=(1, 2),
        help='number of processes'
    )
    parser.add_argument(
        '--munge',
        action='store_true',
        help='munge the input sumstats file'
    )
    parser.add_argument(
        '--snp', help='name of snp column'
    )
    parser.add_argument(
        '-n',
        '--sample-size',
        type=int,
        help='sample size for a sumstats file to be munged'
    )
    parser.add_argument(
        '--a1',
        action='append',
        help='a1 column for a sumstats file to be munged'
    )
    parser.add_argument(
        '--a2',
        action='append',
        help='a2 column for a sumstats file to be munged'
    )
    parser.add_argument(
        '-o',
        '--output-prefix',
        action='append',
        required=True,
        help='output prefix for a file to be processed'
    )
    parser.add_argument(
        '--maf-min',
        type=float,
        default=0.01,
        help='minimum maf to use when munging input'
    )
    
    args = parser.parse_args()
    return args


def main():
    """Perform all main analysis steps"""

    args = parse_arguments()
    file_paths_ordered_and_munged = handle_munging(args)
    annotations = parse_annotation_list(args.annotation_list)
    perform_stratified_analyses(
        args.output_prefix,
        file_paths_ordered_and_munged,
        annotations,
        args.ld_scores,
        intercept=args.intercept
    )
    for output_prefix in args.output_prefix:
        tabulate_statistics(output_prefix, annotations)




# Execute ======================================================================

if __name__ == '__main__':
    main()
