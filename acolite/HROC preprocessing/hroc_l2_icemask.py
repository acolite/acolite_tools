## Script to make ice mask for HROC - development
## Written by Quinten Vanhellemont 2020-11-25

## thresholds tested for 5 l2 scenes
#S2B_MSIL1C_20180926T110029_N0206_R094_T30TXR_20180926T152251-hroc-l2.nc
#S2A_MSIL1C_20180607T104021_N0206_R008_T32UMF_20180607T124742-hroc-l2.nc
#S2A_MSIL1C_20180714T135031_N0206_R110_T26WND_20180714T155602-hroc-l2.nc
#S2A_MSIL1C_20180731T102021_N0206_R065_T33VXE_20180731T123616-hroc-l2.nc
#S2B_MSIL1C_20180426T105029_N0206_R051_T31UET_20180426T144550-hroc-l2.nc

## outputs a mask with three products
## ice_mask == 0 -> likely no ice
## ice_mask & 2**0 -> ice detection based on surface level thresholds on blue and nir, misses some ice on land
## ice_mask & 2**1 -> same as above with additional test on vis (443, 492, 560, 665) reflectance max and mean to further mask ice on land
## ice_mask & 2**3 -> same as above with additional test on TOA SWIR 1.6 micron to exclude bright clouds

def hroc_l2_icemask(file, # input file hroc-l2
                    output=None, # output directory for mask file
                    ice_threshold_blue = 0.025, ## blue threshold for ice mask 2**0
                    ice_threshold_nir = 0.015, ## nir threshold for ice mask 2**0
                    vis_mean_threshold = 0.04, ## vis mean threshold for ice mask 2**1
                    vis_max_threshold = 0.05, ## vis max threshold for ice mask 2**1
                    swir_max = 0.12, ## swir threshold for ice mask 2**2
                    add_rgb = False ## whether to add RGB bands (Rrs_a)
                    ):

    from netCDF4 import Dataset
    import os
    import numpy as np

    ## read datasets
    datasets = [ 'Rrs443_a', 'Rrs492_a', 'Rrs560_a',
                 'Rrs665_a', 'Rrs704_a', 'Rrs740_a',
                 'Rrs783_a', 'Rrs833_a', 'Rrs865_a',
                 'B11' ]

    ## read datasets excluding extra NIR
    datasets = [ 'Rrs443_a', 'Rrs492_a', 'Rrs560_a',
                 'Rrs665_a','Rrs865_a', 'B11' ]

    data = {}
    nci = Dataset(file)
    for ds in datasets:
        data[ds] = nci.variables[ds][:]
    nci.close()
    nci = None

    ## blue/nir tests on Rrs
    mask = (data['Rrs443_a']>ice_threshold_blue) & (data['Rrs865_a']>ice_threshold_nir) & (data['Rrs443_a']>data['Rrs865_a'])

    ice_mask = np.zeros((mask.shape[0], mask.shape[1]), dtype=np.int16)
    ice_mask[mask] = np.int16(2**0)

    ## add test on visible bands
    vis = np.dstack((data['Rrs443_a'], data['Rrs492_a'],data['Rrs560_a'], data['Rrs665_a']))
    #nir = np.dstack((data['Rrs704_a'], data['Rrs740_a'], data['Rrs783_a'], data['Rrs833_a'], data['Rrs865_a']))
    mask = mask | (np.nanmean(vis, axis=2)>vis_mean_threshold) & (np.nanmax(vis, axis=2)>vis_max_threshold)
    ice_mask[mask]+= np.int16(2**1)

    ## add swir test
    mask = mask & (data['B11']<swir_max)
    ice_mask[mask]+= np.int16(2**2)

    mask = None

    ## make output filename
    bn = os.path.basename(file)
    idir = os.path.dirname(file)
    odir = output if output is not None else idir
    if not os.path.exists(odir): os.makedirs(odir)
    ofile = '{}/{}'.format(odir, bn.replace('-hroc-l2.nc', '-hroc-ice.nc'))

    if True:
        ## open output file
        nco = Dataset(ofile, 'w')

        ## set up x and y dimensions
        nco.createDimension('x', ice_mask.shape[1])
        nco.createDimension('y', ice_mask.shape[0])

        ## write mask
        var = nco.createVariable('ice_mask',ice_mask.dtype,('y','x'))
        var[:] = ice_mask

        if add_rgb:
            for v in ['Rrs492_a', 'Rrs560_a','Rrs665_a']:
                var = nco.createVariable(v,data[v].dtype,('y','x'), fill_value=np.nan)
                var[:] = data[v]
        nco.close()
        nco = None

    ice_mask = None
    data = None

## simple wrapper
if __name__ == '__main__':
    import argparse
    import sys, os

    parser = argparse.ArgumentParser(description='HR-OC L2 ice-mask')
    parser.add_argument('--input', help='HROC L2 NetCDF file', default=None)
    parser.add_argument('--output', help='Output directory, if not provided will output to input file directory', default=None)

    args, unknown = parser.parse_known_args()

    if args.input is not None:
        ret = hroc_l2_icemask(args.input,output=args.output)
