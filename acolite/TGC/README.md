## TGC processing
Use tgc.py to run TGC on S2Resampling BEAM-NetCDF file

## Scope and limitations
This script assumed the TOA glint contribution can be estimated from random pixels in B11 by optimising ACOLITE LUTs for a given wind speed and aerosol optical depth. TOA glint is then computed using the band specific geometry and removed from the TOA reflectance.

## Setup
This requires a clone of the "new" generic ACOLITE code: `git clone https://github.com/acolite/acolite --depth=1`, with the path to the base directory to be passed on the command line (--acolite_path), or otherwise be included in your Python path, for example, assuming your git clone is in the home directory:
```
import sys, os
sys.path.append(os.path.expanduser("~/git/acolite"))
```

## Examples
CLI: Convert a given scene and output to the $HOME/TGC/Output directory.
```
python tgc.py --acolite_path ~/git/acolite --output ~/TGC/Output --input ~/TGC/S2RES/S2A_MSIL1C_20190905T113321_N0208_R080_T28RFT_20190905T115747_S2resampling.nc
```

Python:
```
## Import business
import sys, os
user_home = os.path.expanduser("~")
sys.path.append(user_home+'/git/acolite_tools/TGC')
sys.path.append(user_home+'/git/acolite')

from tgc import tgc
inputfile='~/TGC/S2RES/S2A_MSIL1C_20190905T113321_N0208_R080_T28RFT_20190905T115747_S2resampling.nc'
output='~/TGC/Converted'
ret = tgc(ifile, output=output)
```
