# submission script for local linux box

shopt -s expand_aliases
source /.bashrc
source activate IMPgen1
export OMP_NUM_THREADS=6
IMPRESSION train --prefs settings_train.json --tracetime
