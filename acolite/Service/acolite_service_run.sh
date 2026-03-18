## script to activate conda env and run acolite script
## written by Quinten Vanhellemont, RBINS
## 2026-03-05
## modifications: 2026-03-16 (QV) add path testing for conda env source

## assume miniconda is in the user home
home=`eval echo "~$USER"`
conda_env=acolite

## assume other scripts are in the same path
path=`eval dirname $0`

## source conda
if [ -f "$home/miniconda3/etc/profile.d/conda.sh" ]; then
  source "$home/miniconda3/etc/profile.d/conda.sh"
elif [ -f "/opt/miniconda3/etc/profile.d/conda.sh" ]; then
  source "/opt/miniconda3/etc/profile.d/conda.sh"
else
  echo "Could not source conda.sh"
  echo "Add path to miniconda3/etc/profile.d/conda.sh in acolite_service_run.sh"
  exit 1
fi

## activate conda environment
conda activate $conda_env

## run acolite service
if [ "$#" -eq 1 ]; then
  python $path/acolite_service.py $1 ## pass the argument
else
  python $path/acolite_service.py
fi

## deactivate conda environment
conda deactivate
