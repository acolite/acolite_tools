## Testing deployment of ACOLITE GitHub code in a docker image
## QV 2021-01-20
## updated 2021-07-15 for generic ACOLITE
## updated 2021-07-26 for retrieving LUTs
## updated 2021-11-21 new run.sh script
## updated 2022-02-22 added new sensors
## updated 2023-10-17 Ubuntu LTS 22.04, added new sensors
## updated 2025-07-24 added conda TOS acceptance, added S2C, WV sensors, added libgdal and updated to last GitHub version

FROM ubuntu:22.04

# System packages
RUN apt-get update && apt-get install -y curl
RUN apt-get install -y git

## clone acolite
#RUN git clone --depth 1 https://github.com/acolite/acolite

# Install miniconda to /miniconda
RUN curl -LO http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b
RUN rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH=/miniconda/bin:${PATH}

## accept channel TOS and update
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
RUN conda tos accept --override-channels --channel conda-forge
RUN conda update -y conda

# Create acolite environment
#RUN conda create -n acolite -c conda-forge python=3 numpy matplotlib scipy gdal pyproj astropy cartopy scikit-image pyhdf pyresample netcdf4 requests
#RUN conda create -n acolite -c conda-forge python=3 numpy matplotlib scipy gdal pyproj scikit-image pyhdf pyresample netcdf4 h5py requests pygrib cartopy
RUN conda create -n acolite -c conda-forge python=3 numpy matplotlib scipy gdal libgdal-jp2openjpeg libgdal-netcdf pyproj scikit-image pyhdf pyresample netcdf4 h5py requests pygrib cartopy
RUN conda init bash
RUN echo "conda activate acolite" > ~/.bashrc

## clone acolite
RUN git clone --depth 1 https://github.com/acolite/acolite

## get LUTs
RUN conda run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor hyper

RUN conda run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor L8_OLI
RUN conda run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor L9_OLI
RUN conda run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor L5_TM
RUN conda run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor L7_ETM

RUN conda run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor S2A_MSI
RUN conda run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor S2B_MSI
RUN conda run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor S2C_MSI

RUN conda run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor S3A_OLCI
RUN conda run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor S3B_OLCI
RUN conda run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor EN1_MERIS

RUN conda run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor RapidEye
RUN conda run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor SPOT6
RUN conda run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor SPOT7
RUN conda run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor PHR1A
RUN conda run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor PHR1B

RUN conda run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor VENµS_VSSC
RUN conda run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor WorldView2
RUN conda run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor WorldView3
RUN conda run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor WorldView4
RUN conda run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor WV_LG01
RUN conda run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor WV_LG02
RUN conda run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor Quickbird2

RUN conda run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor PlanetScope_0c
RUN conda run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor PlanetScope_0d05
RUN conda run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor PlanetScope_0d06
RUN conda run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor PlanetScope_0e
RUN conda run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor PlanetScope_0f
RUN conda run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor PlanetScope_22
RUN conda run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor PlanetScope_SD5
RUN conda run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor PlanetScope_SD8

RUN conda run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor SUOMI-NPP_VIIRS
RUN conda run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor JPSS-1_VIIRS
RUN conda run -n acolite python ./acolite/launch_acolite.py --retrieve_luts --sensor JPSS-2_VIIRS

# add input, output and settings mount points/files
ADD ./input /input
ADD ./output /output
ADD ./settings settings

## use run script
ADD ./run.sh run.sh
ENTRYPOINT bash run.sh
