## satellittdata.no OPeNDAP preprocessing
Use convert_opendap_nbs.py to create a L1_converted.nc file that can be processed with ACOLITE from a satellittdata.no OPeNDAP URL. Subsetting is supported using a limit ROI (in geographical coordinates:--limit S,W,N,E) or an image pixel subset (--sub y0,x0,y1,x1).

## Setup
This requires a clone of the ACOLITE code: `git clone https://github.com/acolite/acolite`, with the path to the base directory to be passed on the command line (--acolite_path), or otherwise be included in your python path, for example, assuming your git clone is in the home directory:
```
import sys, os
sys.path.append(os.path.expanduser("~/git/acolite"))
```

## Examples
CLI: Crop a and convert scene over Oslofjord and output to the $HOME/Oslo directory.
```
python acolite/satellittdata.no\ OPeNDAP\ preprocessing/convert_opendap_nbs.py --verbosity 1 —-limit 59.84,10.65,59.91,10.81 --output ~/Oslo --acolite_path ~/git/acolite --input http://nbstds.met.no/thredds/dodsC/NBS/S2B/2020/06/30/S2B_MSIL1C_20200630T124709_N0209_R138_T31XFH_20200630T132039.nc
```

CLI: Download and convert a full scene to the current working directory.
```
python acolite/satellittdata.no\ OPeNDAP\ preprocessing/convert_opendap_nbs.py  'http://nbstds.met.no/thredds/dodsC/NBS/S2B/2020/06/30/S2B_MSIL1C_20200630T124709_N0209_R138_T31XFH_20200630T132039.nc'
```

Python: Import in other script, crop a scene and process immediately with ACOLITE.
```
## Import business
import sys, os
user_home = os.path.expanduser("~")
sys.path.append(user_home+'/git/acolite_tools/acolite/satellittdata.no OPeNDAP preprocessing')
sys.path.append(user_home+'/git/acolite')

## Create cropped and converted input
from convert_opendap_nbs import convert_opendap_nbs

ncf = "http://nbstds.met.no/thredds/dodsC/NBS/S2A/2020/06/26/S2A_MSIL1C_20200626T104031_N0209_R008_T32VNM_20200626T125124.nc"
limit = [59.69, 10.47, 59.91, 10.81]
local_dir = '{}/ACOLITE/Input'.format(user_home)

ret = convert_opendap_nbs(ncf, limit=limit, local_dir=local_dir, verbosity=1)
print('Got scene {}'.format(ret))

## Run acolite processing with some basic settings
import acolite as ac
## set some settings (These are from a paper submitted to Optics Express: Vanhellemont, submitted (May 2020)
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
