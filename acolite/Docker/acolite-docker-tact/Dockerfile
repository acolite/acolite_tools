## Testing deployment of ACOLITE GitHub code in a docker image
## QV 2021-01-20
## updated 2021-07-15 for generic ACOLITE
## updated 2021-07-26 for retrieving LUTs
## updated 2021-11-21 libradtran / TACT

FROM ubuntu:20.04

# System packages 
RUN apt-get update && apt-get install -y curl
RUN apt-get install -y git

## For LibRadtran
RUN apt-get install -y wget  build-essential gcc gfortran libgsl-dev flex
RUN DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata 
RUN apt-get install -y netcdf-bin  libnetcdf-dev libnetcdff-dev python gawk

## clone acolite
RUN git clone --depth 1 https://github.com/acolite/acolite

## make libradtran directory in ACOLITE external
RUN mkdir acolite/external
WORKDIR acolite/external

# get Libradtran
RUN wget http://www.libradtran.org/download/history/libRadtran-2.0.2.tar.gz
RUN gzip -d libRadtran-2.0.2.tar.gz
RUN tar -xvf libRadtran-2.0.2.tar
WORKDIR /acolite/external/libRadtran-2.0.2

# configure and compile
RUN ./configure
RUN make
RUN make check

## reset workin directory
WORKDIR /

# Install miniconda to /miniconda
RUN curl -LO http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b
RUN rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH=/miniconda/bin:${PATH}
RUN conda update -y conda

# Create acolite environment
RUN conda create -n acolite -c conda-forge python=3 numpy matplotlib scipy gdal pyproj astropy cartopy scikit-image pyhdf pyresample netcdf4 requests
RUN conda init bash
RUN echo "conda activate acolite" > ~/.bashrc

## get LUTs
RUN conda run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor L8_OLI
RUN conda run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor L5_TM
RUN conda run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor L7_ETM

# add input, output and settings mount points/files
ADD ./input /input
ADD ./output /output
ADD ./settings settings

## use run script
ADD ./run.sh run.sh
ENTRYPOINT bash run.sh
