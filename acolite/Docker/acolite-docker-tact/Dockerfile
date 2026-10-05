## Testing deployment of ACOLITE GitHub code in a docker image
## QV 2021-01-20
## updated 2021-07-15 for generic ACOLITE
## updated 2021-07-26 for retrieving LUTs
## updated 2021-11-21 libradtran / TACT
## updated 2022-02-22 added L9
## updated 2023-10-17 changed to Ubuntu 22.04 and libRadtran 2.0.5
## updated 2025-07-24 added conda TOS acceptance, added libgdal and updated to last GitHub version, changed to libRadtran history 2.0.5
## updated 2026-10-05 changed from miniconda to micromamba updated to last GitHub version

FROM ubuntu:22.04

# System packages
RUN apt-get update && apt-get install -y curl
RUN apt-get install -y git
RUN apt-get install -y bzip2

## For LibRadtran
RUN apt-get install -y wget  build-essential gcc gfortran libgsl-dev flex
RUN DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata
RUN apt-get install -y netcdf-bin  libnetcdf-dev libnetcdff-dev python3 gawk

## clone acolite
RUN git clone --depth 1 https://github.com/acolite/acolite

## make libradtran directory in ACOLITE external
RUN mkdir acolite/external
WORKDIR acolite/external

# get Libradtran
RUN wget http://www.libradtran.org/download/history/libRadtran-2.0.5.tar.gz
RUN gzip -d libRadtran-2.0.5.tar.gz
RUN tar -xvf libRadtran-2.0.5.tar
WORKDIR /acolite/external/libRadtran-2.0.5

# configure and compile
RUN ./configure
RUN make
RUN make check

## reset workin directory
WORKDIR /

# Install micromamba
RUN curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
RUN eval "$(./bin/micromamba shell hook -s posix)"

# Create acolite environment
RUN micromamba create --file acolite/environment.yml
RUN echo "micromamba activate acolite" > ~/.bashrc

## get LUTs
RUN micromamba run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor L8_OLI
RUN micromamba run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor L9_OLI
RUN micromamba run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor L5_TM
RUN micromamba run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor L7_ETM

# add input, output and settings mount points/files
ADD ./input /input
ADD ./output /output
ADD ./settings settings

## use run script
ADD ./run.sh run.sh
ENTRYPOINT bash run.sh
