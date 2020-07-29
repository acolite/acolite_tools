## HR-OC preprocessing
Use convert_hroc_merged.py to create a L1_converted.nc file that can be processed with ACOLITE from a merged HR-OC L1 file.

## Scope and limitations
This script converts merged HR-OC L1 files to files that resembles the output of the ACOLITE cropping/merging tool. This file can then be processed by ACOLITE, e.g. by setting the inputfile parameter in your settings file. For the sun and view geometry, full scene interpolated datasets are provided, but these are (not yet) properly supported in ACOLITE processing. A single geometry is stored in the NetCDF global attributes to be used by ACOLITE based on the scene average. **Note: for appropriate processing of these large area merged data, ACOLITE needs to be adapted to take into account the variable geometry.**

## Setup
This requires a clone of the ACOLITE code: `git clone https://github.com/acolite/acolite`, with the path to the base directory to be passed on the command line (--acolite_path), or otherwise be included in your python path, for example, assuming your git clone is in the home directory:
```
import sys, os
sys.path.append(os.path.expanduser("~/git/acolite"))
```

## Examples
CLI: Convert a given scene and output to the $HOME/HROC/Converted directory.
```
python acolite/HROC\ preprocessing/convert_hroc_merged.py --verbosity 1 --acolite_path ~/git/acolite --output ~/HROC/Converted --input ~/HROC/Example/2019-06-02-35T-mosaic-60.nc
```

Python: Import in other script, crop a scene and process immediately with ACOLITE.
```
## Import business
import sys, os
user_home = os.path.expanduser("~")
sys.path.append(user_home+'/git/acolite_tools/acolite/HROC\ preprocessing')
sys.path.append(user_home+'/git/acolite')

## Create cropped and converted input
from convert_hroc_merged import convert_hroc_merged

ifile = " {}/HROC/Example/2019-06-02-35T-mosaic-60.nc".format(user_home)
local_dir = '{}/HROC/Converted'.format(user_home)

ret = convert_hroc_merged(ifile, limit=limit, local_dir=local_dir, verbosity=1)
print('Got scene {}'.format(ret))

## Run acolite processing with some basic settings
import acolite as ac
## set some settings, these are from a paper submitted to Optics Express: Vanhellemont, submitted (May 2020)
settings = {'glint_correction':True,
            'sky_correction_option':'rsky_new',
            'dsf_wave_range':(400,900),
            'ancillary_data':False,
            'l2w_parameters':['t_nechad'],
            'map_l2w' : True}

ac.acolite.acolite_run(inputfile=ret,
                       output='{}/ACOLITE/Output'.format(user_home),
                       settings=settings)
```
