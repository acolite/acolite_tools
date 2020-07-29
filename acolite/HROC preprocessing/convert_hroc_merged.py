## Script to preprocess merged data from HR-OC project and convert it to a L1 NetCDF file that ACOLITE can process
## Written by Quinten Vanhellemont 2020-07-14

def convert_hroc_merged(ifile,
                        local_dir=None, ## directory to store the converted file *recommended* otherwise wd is used
                        override=True, ## whether to overwrite the ncfile if the path already exists
                        satellite_ = 'Sentinel-2A', ## fix the satellite identifier - to be added in metadata/filename!
                        verbosity=0):

    import acolite as ac
    import numpy as np
    import os
    import datetime, dateutil.parser

    gatts = ac.shared.nc_gatts(ifile)
    datasets = ac.shared.nc_datasets(ifile)

    ## extract info from file name (nothing given in NetCDF header)
    bn, ext = os.path.splitext(os.path.basename(ifile))
    sp = bn.split('-')

    ## for June file
    year, month, day, zone, tile_code, resolution = sp

    ## parse datetime from the input filename
    dt = datetime.datetime(int(year), int(month), int(day))

    ## set up metadata that ACOLITE expects
    metadata = {}
    metadata['TIME'] = dt.isoformat()

    ## fixed from external inputs
    ## needs to be improved
    metadata['SATELLITE'] = satellite_

    if metadata['SATELLITE'] == 'Sentinel-2A':
        metadata['SENSOR'] = 'MSI'
        metadata['SATELLITE_SENSOR'] = 'S2A_MSI'
        metadata['BANDS'] = '1,2,3,4,5,6,7,8,8A,11,12'
        metadata['BANDS'] = '1,2,3,4,5,6,7,8,8A,9,10,11,12'
        metadata['BAND_NAMES'] = 'B1,B2,B3,B4,B5,B6,B7,B8,B8A,B9,B10,B11,B12'
        metadata['BANDS_BESTFIT'] = ['11','12']

    if metadata['SATELLITE'] == 'Sentinel-2B':
        metadata['SENSOR'] = 'MSI'
        metadata['SATELLITE_SENSOR'] = 'S2B_MSI'
        metadata['BANDS'] = '1,2,3,4,5,6,7,8,8A,11,12'
        metadata['BANDS'] = '1,2,3,4,5,6,7,8,8A,9,10,11,12'
        metadata['BAND_NAMES'] = 'B1,B2,B3,B4,B5,B6,B7,B8,B8A,B9,B10,B11,B12'
        metadata['BANDS_BESTFIT'] = ['11','12']

    ## get band averaged wavelengths
    swaves = ac.shared.sensor_wave(metadata['SATELLITE_SENSOR'])
    waves = [int(swaves[band.strip('B')]) for ib, band in enumerate(metadata['BAND_NAMES'].split(','))]
    metadata['WAVES'] = np.array(waves)

    ## datetime info
    dt = dateutil.parser.parse(metadata['TIME'])
    metadata["DOY"] = dt.strftime('%j')
    metadata["SE_DISTANCE"] = ac.shared.distance_se(metadata['DOY'])

    ## some tags that ACOLITE accesses from the metadata
    metadata['SCENE'] = bn
    metadata['TILE_CODE'] = '_'.join([zone, tile_code, resolution])
    metadata['OBASE'] = '{}_{}_{}'.format(metadata['SATELLITE_SENSOR'],dt.strftime('%Y_%m_%d_%H_%M_%S'), metadata['TILE_CODE'])

    ## placeholder for scene projection
    ## can be reversed from the lat, lon and UTM zone?
    pixelsize = (int(resolution), int(resolution))
    xrange = (0,0)
    yrange = (0,0)
    proj4_string = ''

    metadata['xrange'] = xrange
    metadata['yrange'] = yrange
    metadata['proj4_string'] = proj4_string
    metadata['pixel_size'] = pixelsize
    ##

    ## put mean geometry in metadata
    vza_ = np.nanmean(ac.shared.nc_data(ifile, 'view_zenith_mean'))
    vaa_ = np.nanmean(ac.shared.nc_data(ifile, 'view_azimuth_mean'))
    sza_ = np.nanmean(ac.shared.nc_data(ifile, 'sun_zenith'))
    saa_ = np.nanmean(ac.shared.nc_data(ifile, 'sun_azimuth'))

    ## update metadata with sun and view geometry
    metadata['THS'] = sza_
    metadata['THV'] = vza_

    sun_azi = saa_
    view_azi = vaa_

    azi = abs(sun_azi - view_azi)
    while(azi >= 180.):
        azi -= 180.

    metadata['AZI']=abs(azi)
    metadata['ISODATE']=dt.isoformat()
    ## end metadata

    ## required dataset names, and name translation
    ## keep lat lon and save geometry as sza, saa, vza, vaa
    ## ACOLITE expects the wavelength labeled TOA reflectance in the L1 NCDF: rhot_xxx
    required = [{'in':'lon', 'out':'lon'},
                {'in':'lat', 'out':'lat'},
                {'in':'view_zenith_mean', 'out':'vza'},
                {'in':'view_azimuth_mean', 'out':'vaa'},
                {'in':'sun_zenith', 'out':'sza'},
                {'in':'sun_azimuth', 'out':'saa'}]

    for ib, band in enumerate(metadata['BAND_NAMES'].split(',')):
        wk = 'WAVES'
        required.append({'in':band, 'out':'rhot_{}'.format(metadata[wk][ib])})

    ## get output directory
    if local_dir is not None:
        odir = local_dir
    else:
        odir = os.cwd()
    if not os.path.exists(odir): os.makedirs(odir)

    ## output file naming
    ncfile = '{}/{}_L1_converted.nc'.format(odir,metadata['OBASE'])
    if (os.path.exists(ncfile)) & (not override):
        if verbosity > 0: print('File exists {}'.format(ncfile))
        nc.close()
        return(ncfile)

    ## read data and output to new file
    new = True ## creates a new NetCDF file
    for r in required:
        dsi = r['in']
        dso = r['out']
        dat = None

        if dsi not in datasets:
            if verbosity > 0: print('{} not found in {}'.format(dsi, ncf_))
            continue
        else:
            if verbosity > 0: print('Reading {} and writing {}'.format(dsi, dso))

            dat, atti = ac.shared.nc_data(ifile, dsi, attributes=True)
            dat[dat.mask] = np.nan

            atto = {a:atti[a] for a in atti if a not in ['_FillValue']}

            if ('Sentinel-2' in metadata['SATELLITE']) & ('rhot_' in r['out']):
                atto['wavelength'] = float(r['out'].split('_')[1])
                atto['units'] = 1.0

            ac.output.nc_write(ncfile, dso, dat, attributes=metadata, dataset_attributes=atto, new=new)
            new = False ## appends to the NetCDF file

    ## create relative azimuth
    if verbosity > 0: print('Creating raa dataset')
    saa = ac.shared.nc_data(ncfile, 'saa')
    vaa = ac.shared.nc_data(ncfile, 'vaa')
    raa = np.abs(saa-vaa)
    raa[raa>180] -= 180
    ac.output.nc_write(ncfile, 'raa', raa)

    ## return the path to the new file
    return(ncfile)

if __name__ == '__main__':
    import argparse
    import sys, os

    parser = argparse.ArgumentParser(description='HR-OC preprocessor for ACOLITE')
    parser.add_argument('--input', help='NetCDF file containing merged L1C Sentinel-2 data')
    parser.add_argument('--output', help='Output directory, if not provided will output to current working directory', default=None)
    parser.add_argument('--verbosity', help='Verbosity (default=0)', default=0)
    parser.add_argument('--override', help='Overwrite NetCDF files (default=True)', default=True)
    parser.add_argument('--acolite_path', help='Path to ACOLITE code (if not already in Python path) (default=None)', default=None)
    parser.add_argument('--satellite', help='Specifiy satellite (default=Sentinel-2A)', default='Sentinel-2A')

    args, unknown = parser.parse_known_args()

    ## acolite path
    if args.acolite_path is not None:
        sys.path.append(os.path.expanduser(args.acolite_path))
    try:
        import acolite as ac
    except:
        print('Could not import acolite - make sure it is in your Python path, or provide it with --acolite_path')
        exit(1)

    ## parse output directory
    if args.output is not None:
        args.output = os.path.expanduser(args.output)

    ## evaluate override flag
    if type(args.override) == str: args.override = eval(args.override)

    if args.input is not None:

        ## convert the remote file
        ret = convert_hroc_merged(args.input,
                                local_dir=args.output,
                                override=args.override,
                                satellite_=args.satellite,
                                verbosity=int(args.verbosity))
