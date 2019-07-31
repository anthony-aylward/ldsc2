#===============================================================================
# download.py
#===============================================================================

# Imports ======================================================================

import gzip
import os
import os.path
import subprocess

from argparse import ArgumentParser
from git import Git
from hashlib import sha256
from urllib.request import urlopen
from shutil import copyfileobj
from tarfile import TarFile
from tempfile import TemporaryDirectory

from ldsc2.env import (
    DIR, ANACONDA_DIR, HAPMAP3_SNPS, PLINKFILES, PLINKFILES_EAS, BASELINE, BASELINE_EAS
)




# Constants ====================================================================

ANACONDA_URL = os.environ.get(
    'LDSC_ANACONDA_URL',
    'https://repo.anaconda.com/archive/Anaconda2-2019.03-Linux-x86_64.sh'
)
ANACONDA_HASH = os.environ.get(
    'LDSC_ANACONDA_HASH',
    'cedfee5b5a3f62fcdac0a1d2d12396d0f232d2213d24d6dc893df5d8e64b8773'
)
CONDA_PATH = os.path.join(ANACONDA_DIR, 'bin', 'conda')
LDSC_GITHUB_REPO = 'https://github.com/bulik/ldsc.git'
PLINKFILES_URL = 'https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_plinkfiles.tgz'
PLINKFILES_EAS_URL = 'https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_EAS_plinkfiles.tgz'
HAPMAP3_SNPS_URL = 'https://data.broadinstitute.org/alkesgroup/LDSCORE/hapmap3_snps.tgz'
BASELINE_URL = 'https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_baseline_v1.1_ldscores.tgz'
BASELINE_EAS_URL = 'https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_EAS_baseline_v1.2_ldscores.tgz'



# Functions ====================================================================

def download_anaconda_install_script(anaconda_install_script_path):
    print(
        'Downloading Anaconda2 install script to '
        f'{anaconda_install_script_path}'
    )
    with urlopen(ANACONDA_URL) as (
        response
    ), open(anaconda_install_script_path, 'wb') as (
        f
    ):
        copyfileobj(response, f)


def check_hash(anaconda_install_script_path):
    print(f'checking hash of {anaconda_install_script_path}')
    with open(anaconda_install_script_path, 'rb') as f:
        if sha256(f.read()).hexdigest() != ANACONDA_HASH:
            raise RuntimeError(f'hash check failed for {ANACONDA_URL}')


def install_anaconda(anaconda_install_script_path):
    input(
        'installing Anaconda2. When prompted, specify the following '
        f'install location:\n{ANACONDA_DIR}\npress ENTER to '
        'continue >>>'
    )
    subprocess.run(('bash', anaconda_install_script_path))


def configure_anaconda():
    print('configuring Anaconda2')
    subprocess.run(
        (
            CONDA_PATH, 'env', 'create',
            '--file', os.path.join(DIR, 'ldsc', 'environment.yml')
        )
    )


def clone_ldsc():
    print(f"cloning the LDSC github repo to {os.path.join(DIR, 'ldsc')}")
    Git(DIR).clone(LDSC_GITHUB_REPO)


def download(
    hapmap3_snps_dir: str = HAPMAP3_SNPS,
    plinkfiles_dir: str = PLINKFILES,
    plinkfiles_eas_dir: str = PLINKFILES_EAS,
    baseline_dir: str = BASELINE,
    baseline_eas_dir: str = BASELINE_EAS,
    quiet: bool = False
):
    for data_url, data_dir in (
        (HAPMAP3_SNPS_URL, hapmap3_snps_dir),
        (PLINKFILES_URL, plinkfiles_dir),
        (PLINKFILES_EAS_URL, plinkfiles_eas_dir),
        (BASELINE_URL, baseline_dir),
        (BASELINE_EAS_URL, baseline_eas_dir)
    ):
        if not quiet:
            print(f'Downloading data to {data_dir}.tgz')
        with urlopen(data_url) as (
            response
        ), open(f'{data_dir}.tgz', 'wb') as (
            f
        ):
            copyfileobj(response, f)
        if not quiet:
            print(f'Extracting data to {data_dir}')
        subprocess.run(
            ('tar', '-C', os.path.dirname(data_dir), '-zxvf', f'{data_dir}.tgz')
        )


def extract_blank_annot():
    baseline_dir = os.path.join(DIR, 'baseline_v1.1')
    baseline_eas_dir = DIR
    blank_dir = os.path.join(DIR, 'blank')
    blank_eas_dir = os.path.join(DIR, 'blank_eas')
    for d in blank_dir, blank_eas_dir:
        if not os.path.isdir(d):
            os.mkdir(d)
    for chrom in range(1, 23):
        with gzip.open(
            os.path.join(baseline_dir, f'baseline.{chrom}.annot.gz'), 'rt'
        ) as f0, gzip.open(
            os.path.join(blank_dir, f'blank.{chrom}.annot.gz'), 'wb'
        ) as f1:
            f1.write(
                ''.join(
                    '\t'.join(line.split()[:4]) + '\n' for line in f0
                ).encode()
            )
        with gzip.open(
            os.path.join(baseline_eas_dir, f'baseline.{chrom}.annot.gz'), 'rt'
        ) as f0, gzip.open(
            os.path.join(blank_eas_dir, f'blank.{chrom}.annot.gz'), 'wb'
        ) as f1:
            f1.write(
                ''.join(
                    '\t'.join(line.split()[:4]) + '\n' for line in f0
                ).encode()
            )


def parse_arguments():
    parser = ArgumentParser(
        description='download components for LDSC'
    )
    parser.add_argument(
        '--ldsc-dir',
        metavar='<path/to/dir/>',
        default=DIR,
        help=f'directory in which to download LDSC data [{DIR}]'
    )
    parser.add_argument(
        '--plinkfiles',
        metavar='<dest/for/plinkfiles/dir>',
        default=PLINKFILES,
        help=(
            'destination for downloaded EUR plink files'
            f'[{PLINKFILES}]'
        )
    )
    parser.add_argument(
        '--plinkfiles-eas',
        metavar='<dest/for/plinkfiles/dir>',
        default=PLINKFILES_EAS,
        help=(
            'destination for downloaded EAS plink files'
            f'[{PLINKFILES_EAS}]'
        )
    )
    parser.add_argument(
        '--baseline',
        metavar='<dest/for/baseline/dir>',
        default=BASELINE,
        help=(
            'destination for downloaded EUR baseline files'
            f'[{BASELINE}]'
        )
    )
    parser.add_argument(
        '--baseline-eas',
        metavar='<dest/for/baselinedir>',
        default=BASELINE_EAS,
        help=(
            'destination for downloaded EAS baseline files'
            f'[{BASELINE_EAS}]'
        )
    )
    parser.add_argument(
        '--hapmap3-snps',
        metavar='<dest/for/hapmap3/dir>',
        default=HAPMAP3_SNPS,
        help=(
            'destination for downloaded SNP files'
            f'[{HAPMAP3_SNPS}]'
        )
    )
    parser.add_argument(
        '--quiet',
        action='store_true',
        help='suppress status updates'
    )
    parser.add_argument(
        '--tmp-dir',
        metavar='<path/to/tmp/dir>',
        help='directory for temporary files'
    )
    return parser.parse_args()


def main():
    args = parse_arguments()
    if os.path.isdir(ANACONDA_DIR):
        use_existing_anaconda_dir = (
            input(
                f'There is already a directory at {ANACONDA_DIR} - is this the '
                'anaconda you wish to use? ([y]/n) >>>'
            ).casefold()
            in {'', 'y', 'yes'}
        )
        if not use_existing_anaconda_dir:
            print(
                'Please change the value of environment variable '
                'LDSC2_ANACONDA_DIR or remove the existing directory at '
                'that location.'
            )
            return
    elif os.path.exists(ANACONDA_DIR):
        raise RuntimeError(
            f'There is a non-directory file at {ANACONDA_DIR}. Please change '
            'the value of environment variable LDSC2_ANACONDA_DIR or '
            'remove the existing file at that location.'
        )
    else:
        use_existing_anaconda_dir = False
    ldsc_dir = os.path.join(DIR, 'ldsc')
    if os.path.isdir(ldsc_dir):
        use_existing_ldsc_dir = (
            input(
                f'There is already a directory at {ldsc_dir} - is this the '
                'LDSC you wish to use? ([y]/n) >>>'
            ).casefold() in {'', 'y', 'yes'}
        )
        if not use_existing_ldsc_dir:
            print(
                'Please change the value of environment variable LDSC2_DIR '
                'or remove the existing directory at that location.'
            )
            return
    elif os.path.exists(ldsc_dir):
        raise RuntimeError(
            f'There is a non-directory file at {ldsc_dir} Please change '
            'the value of environment variable LDSC2_DIR or '
            'remove the existing file at that location.'
        )
    else:
        use_existing_ldsc_dir = False
    if not use_existing_ldsc_dir:
        clone_ldsc()
    if not use_existing_anaconda_dir:
        with TemporaryDirectory(dir=args.tmp_dir) as temp_dir:
            anaconda_install_script_path = os.path.join(
                temp_dir, 'Anaconda2-2019.03-Linux-x86_64.sh'
            )
            download_anaconda_install_script(anaconda_install_script_path)
            check_hash(anaconda_install_script_path)
            install_anaconda(anaconda_install_script_path)
        configure_anaconda()
    download(
        plinkfiles_dir=args.plinkfiles,
        plinkfiles_eas_dir=args.plinkfiles_eas,
        baseline_dir=args.baseline,
        baseline_eas_dir=args.baseline_eas,
        hapmap3_snps_dir=args.hapmap3_snps,
        quiet=args.quiet
    )
    extract_blank_annot()
