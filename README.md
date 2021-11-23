## ACOLITE Tools
Various scripts and tools that can be used for [ACOLITE](https://github.com/acolite/acolite).

Many of these tools require a clone of the ACOLITE code to be in your python path, for example if you have a git/acolite directory in your $HOME directory:
```
import sys, os
sys.path.append(os.path.expanduser('~/git/acolite'))
```

## ACOLITE
* Docker: Dockerfiles for ACOLITE processing
* satellittdata.no OPeNDAP preprocessing: creates a L1_converted.nc file from a satellittdata.no OPeNDAP URL that can be processed with ACOLITE
