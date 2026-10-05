## Testing deployment of ACOLITE GitHub code in a docker image
## QV 2021-01-20
## updated 2021-07-15 for generic ACOLITE
## updated 2021-07-26 for retrieving LUTs
## updated 2021-11-21 new run.sh script
## updated 2022-02-22 added new sensors
## updated 2023-10-17 Ubuntu LTS 22.04, added new sensors
## updated 2025-07-24 added conda TOS acceptance, added S2C, WV sensors, added libgdal and updated to last GitHub version
## updated 2026-10-05 moved to micromamba updated to last GitHub version

FROM ubuntu:22.04

# System packages
RUN apt-get update && apt-get install -y curl
RUN apt-get install -y git
RUN apt-get install -y bzip2

## clone acolite
RUN git clone --depth 1 https://github.com/acolite/acolite

# Install micromamba
RUN curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
RUN eval "$(./bin/micromamba shell hook -s posix)"

# Create acolite environment
RUN micromamba create --file acolite/environment.yml
RUN echo "micromamba activate acolite" > ~/.bashrc

## get LUTs
RUN micromamba run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor hyper

RUN micromamba run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor L8_OLI
RUN micromamba run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor L9_OLI
RUN micromamba run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor L5_TM
RUN micromamba run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor L7_ETM

RUN micromamba run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor S2A_MSI
RUN micromamba run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor S2B_MSI
RUN micromamba run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor S2C_MSI

RUN micromamba run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor S3A_OLCI
RUN micromamba run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor S3B_OLCI
RUN micromamba run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor EN1_MERIS

RUN micromamba run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor RapidEye
RUN micromamba run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor SPOT6
RUN micromamba run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor SPOT7
RUN micromamba run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor PHR1A
RUN micromamba run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor PHR1B

RUN micromamba run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor VENµS_VSSC
RUN micromamba run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor WorldView2
RUN micromamba run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor WorldView3
RUN micromamba run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor WorldView4
RUN micromamba run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor WV_LG01
RUN micromamba run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor WV_LG02
RUN micromamba run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor Quickbird2

RUN micromamba run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor PlanetScope_0c
RUN micromamba run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor PlanetScope_0d05
RUN micromamba run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor PlanetScope_0d06
RUN micromamba run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor PlanetScope_0e
RUN micromamba run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor PlanetScope_0f
RUN micromamba run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor PlanetScope_22
RUN micromamba run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor PlanetScope_SD5
RUN micromamba run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor PlanetScope_SD8

RUN micromamba run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor SUOMI-NPP_VIIRS
RUN micromamba run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor JPSS-1_VIIRS
RUN micromamba run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor JPSS-2_VIIRS

# add input, output and settings mount points/files
ADD ./input /input
ADD ./output /output
ADD ./settings settings

## use run script
ADD ./run.sh run.sh
ENTRYPOINT bash run.sh
