# script to run ACOLITE in Docker container
# QV 2021-11-21

## run ACOLITE processing
conda run -n acolite python ./acolite/launch_acolite.py --cli --settings settings

## chown output
chown -R 1000:1000 /output
