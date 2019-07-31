#===============================================================================
# env.py
#===============================================================================

# Imports ======================================================================

import os
import os.path




# Constants ====================================================================

DIR = os.environ.get(
    'LDSC2_DIR',
    os.path.dirname(__file__)
)
ANACONDA_DIR = os.environ.get(
    'LDSC2_ANACONDA_DIR',
    os.path.join(os.path.dirname(__file__), 'anaconda2')
)
ANACONDA_PATH = os.path.join(ANACONDA_DIR, 'envs', 'ldsc', 'bin', 'python')
HAPMAP3_SNPS = os.environ.get(
    'LDSC2_HAPMAP3_SNPS',
    os.path.join(DIR, 'hapmap3_snps')
)
PLINKFILES = os.environ.get(
    'LDSC2_PLINKFILES',
    os.path.join(DIR, '1000G_EUR_Phase3_plink')
)
PLINKFILES_EAS = os.environ.get(
    'LDSC2_PLINKFILES_EAS',
    os.path.join(DIR, '1000G_Phase3_EAS_plinkfiles')
)
BASELINE = os.environ.get(
    'LDSC2_BASELINE',
    os.path.join(DIR, '1000G_Phase3_baseline')
)
BASELINE_EAS = os.environ.get(
    'LDSC2_BASELINE_EAS',
    os.path.join(DIR, '1000G_Phase3_EAS_baseline')
)
BLANK = os.environ.get('LDSC2_BLANK', os.path.join(DIR, 'blank'))
BLANK_EAS = os.environ.get('LDSC2_BLANK_EAS', os.path.join(DIR, 'blank_eas'))