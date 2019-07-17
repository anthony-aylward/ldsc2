#===============================================================================
# download.py
#===============================================================================

# Imports ======================================================================

import os
import os.path

from argparse import ArgumentParser
from git import Git
from urllib.request import urlopen
from shutil import copyfileobj
from tarfile import TarFile

from ldsc2.env import DIR, HAPMAP3_SNPS, PLINKFILES, PLINKFILES_EAS




# Constants ====================================================================

LDSC_GITHUB_REPO = 'https://github.com/bulik/ldsc.git'
PLINKFILES_URL = 'https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_plinkfiles.tgz'
PLINKFILES_EAS_URL = 'https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_EAS_plinkfiles.tgz'
HAPMAP3_SNPS_URL = 'https://data.broadinstitute.org/alkesgroup/LDSCORE/hapmap3_snps.tgz'




# Functions ====================================================================

def download(
    ldsc_dir: str = DIR,
    hapmap3_snps_dir: str = HAPMAP3_SNPS,
    plinkfiles_dir: str = PLINKFILES,
    plinkfiles_eas_dir: str = PLINKFILES_EAS,
    quiet: bool = False
):
    if not quiet:
        print(
            f"Cloning LDSC repository to {os.path.join(ldsc_dir, 'LDSC')}"
        )
    Git(ldsc_dir).clone(LDSC_GITHUB_REPO)

    for data_url, data_dir in (
        (HAPMAP3_SNPS_URL, hapmap3_snps_dir),
        (PLINKFILES_URL, plinkfiles_dir),
        (PLINKFILES_EAS_URL, plinkfiles_eas_dir)
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
        with TarFile(f'{data_dir}.tgz') as f:
            f.extractall(data_dir)


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
            'destination for downloaded plink files'
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
    return parser.parse_args()


def main():
    args = parse_arguments()
    download(
        ldsc_dir=args.ldsc_dir,
        plinkfiles_dir=args.plinkfiles,
        plinkfiles_eas_dir=args.plinkfiles,
        hapmap3_snps_dir=args.hapmap3_snps,
        quiet=args.quiet
    )
