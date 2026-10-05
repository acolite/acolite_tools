# script to run ACOLITE in Docker container
# QV 2021-11-21
# QV 2026-10-05 changed to micromamba

## run ACOLITE processing
micromamba run -n acolite python ./acolite/launch_acolite.py --cli --settings settings

## chown output
chown -R 1000:1000 /output
