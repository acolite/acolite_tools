## Script to access OPeNDAP data from satellittdata.no and convert it to a L1 NetCDF file that ACOLITE can process
## Written by Quinten Vanhellemont 2020-07-13

def convert_opendap_nbs(ncf,
                        local_dir=None, ## directory to store the converted file *recommended* otherwise wd is used
                        limit=None, ## limit of ROI [S, W, N, E] to crop data
                        sub=None, ## or subset directly [y0, x0, y1, x1] (counted from top left corner)
                        geometry='mid', ## resolved geometry is currently not supported for NC inputs to ACOLITE
                                        ## here the options are to use "mid" the scene/ROI midpoint geometry (faster)
                                        ## here the options are to use "mean" the scene/ROI average geometry (slower)
                        verbosity=0):

    import netCDF4, os
    import numpy as np
    import dateutil.parser
    from pyproj import Proj
    import acolite as ac

    ## add fillmismatch fix to url
    fm = '[FillMismatch]'
    if fm not in ncf:
        ncf_ = '{}{}'.format(fm, ncf)
    else:
        ncf_ = '{}'.format(ncf)

    ## open NetCDF file
    nc = netCDF4.Dataset(ncf_)

    ## list remote datasets
    datasets = nc.variables.keys()

    ## list remote global attributes
    gatts = {attr : getattr(nc,attr) for attr in nc.ncattrs()}

    ## set up metadata that ACOLITE expects
    metadata = {}
    try:
        metadata['TIME'] = gatts['PRODUCT_START_TIME']
        metadata['SATELLITE'] = gatts['DATATAKE_1_SPACECRAFT_NAME']
    except:
        if verbosity > 0: print('Dataset not recognised')

    ## add BANDS and WAVE data that is normally coming from the metadata scripts
    ## these should be removed from future acolite versions
    ## placeholder there is no L8 data on satellittdata.no
    if metadata['SATELLITE'] == 'LANDSAT_8':
        metadata['SENSOR'] = 'OLI'
        metadata['SATELLITE_SENSOR'] = 'L8_OLI'
        metadata['BANDS_ALL'] = '1,2,3,4,5,6,7,9,10,11'
        metadata['BAND_NAMES'] = '1,2,3,4,5,6,7'
        metadata['BANDS_BESTFIT'] = ['6','7']

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
    if metadata['SENSOR'] ==  'MSI':
        metadata['WAVES'] = np.array(waves)
    else:
        metadata['WAVES_ALL'] = np.array(waves)

    ## datetime info
    dt = dateutil.parser.parse(metadata['TIME'])
    metadata["DOY"] = dt.strftime('%j')
    metadata["SE_DISTANCE"] = ac.shared.distance_se(metadata['DOY'])

    ## some tags that ACOLITE accesses from the metadata
    metadata['SCENE'] = os.path.splitext(gatts['PRODUCT_URI'])[0]
    metadata['TILE_CODE'] = metadata['SCENE'].split('_')[-2]
    metadata['OBASE'] = '{}_{}_{}'.format(metadata['SATELLITE_SENSOR'],dt.strftime('%Y_%m_%d_%H_%M_%S'), metadata['TILE_CODE'])

    ## store the original attributes
    for k in gatts: metadata['{}_{}'.format('ext', k)] = gatts[k]

    ## required dataset names, and name translation
    ## ACOLITE expects the wavelength labeled TOA reflectance in the L1 NCDF: rhot_xxx
    required = [{'in':'lon', 'out':'lon'}, {'in':'lat', 'out':'lat'}]
    for ib, band in enumerate(metadata['BAND_NAMES'].split(',')):
        wk = 'WAVES'
        if wk not in metadata: wk = 'WAVES_ALL'
        required.append({'in':band, 'out':'rhot_{}'.format(metadata[wk][ib])})

    ## get scene projection info
    sdim = nc.variables['lat'].shape
    xscene = int(nc.variables['x'][0].data), int(nc.variables['x'][-1].data)
    yscene = int(nc.variables['y'][0].data), int(nc.variables['y'][-1].data)
    proj4_string = nc.variables['UTM_projection'].proj4_string
    pixelsize_ = (yscene[1]-yscene[0]) / (sdim[0]-1), (xscene[1]-xscene[0]) / (sdim[1]-1)
    pixelsize = abs(pixelsize_[0]), abs(pixelsize_[1])

    ## add scene extent test
    metadata['xrange'] = xscene
    metadata['yrange'] = yscene[1], yscene[0]
    metadata['proj4_string'] = proj4_string
    metadata['pixel_size'] = pixelsize

    p = Proj(proj4_string)

    ## to subset?
    ## if limit given and no sub
    if (sub is None) and (limit is not None):

        if False:
            ## Slow!
            ## we will crop on a subset of the full lat/lon,
            ## so final cropped coordinates are +- stride*10m pixels
            sub_stride = 50

            ## read full datasets
            lat = nc.variables['lat'][0:sdim[0]:sub_stride, 0:sdim[1]:sub_stride]
            lon = nc.variables['lon'][0:sdim[0]:sub_stride, 0:sdim[1]:sub_stride]

            ## get south west corner
            tmp = ((lat - limit[0])**2 + (lon - limit[1])**2)**0.5
            ye, xs = np.where(tmp == tmp.min())
            ye=ye[0]
            xs=xs[0]

            ## get north east corner
            tmp = ((lat - limit[2])**2 + (lon - limit[3])**2)**0.5
            ys, xe = np.where(tmp == tmp.min())
            ys=ys[0]
            xe=xe[0]

            sub = [ys*sub_stride, xs*sub_stride, ye*sub_stride, xe*sub_stride]
            del lat
            del lon
        else:
            ## Fast!
            ## with projection info from the 'UTM_projection' object
            ## does S2 data ever come in Polar Stereographic projection?
            xrange_raw, yrange_raw = p([limit[1],limit[3]],[limit[0],limit[2]])

            xrange_region = [xrange_raw[0] - (xrange_raw[0] % pixelsize[0]),
                             xrange_raw[1] - (xrange_raw[1] % pixelsize[0])]
            yrange_region = [yrange_raw[0] - (yrange_raw[0] % pixelsize[1]),
                             yrange_raw[1] - (yrange_raw[1] % pixelsize[1])]

            xoff = [(i - min(xscene))/pixelsize[0] for i in xrange_region]
            yoff = [(i - min(yscene))/abs(pixelsize[1]) for i in yrange_region]
            yoff = [sdim[1]-yoff[0], sdim[1]-yoff[1]]

            sub = [int(yoff[1]), int(xoff[0]), int(yoff[0]), int(xoff[1])]

            metadata['xrange'] = xrange_region
            metadata['yrange'] = yrange_region

    if sub is not None:
        ## check if the crop is fully outside the scene
        if (sub[0] > sdim[0]) | (sub[2] < 0) |  (sub[1] > sdim[1]) | (sub[3] < 0):
            if verbosity > 0: print('ROI outside of scene')
            return()

        ## partial crops are allowed
        sub[0] = np.max((0, sub[0]))
        sub[1] = np.max((0, sub[1]))
        sub[2] = np.min((sdim[0], sub[2]))
        sub[3] = np.min((sdim[1], sub[3]))
        if verbosity > 0: print('Subsetting scene to {}'.format(sub))
    else:
        if verbosity > 0: print('Processing full scene')

    ## get geometry data
    ## use mid point of scene/region
    if geometry == 'mid':
        if sub is None:
            midy = int(sdim[0]/2)
            midx = int(sdim[1]/2)
        else:
            midy = int((sub[0]+sub[2])/2)
            midx = int((sub[1]+sub[3])/2)

        sza_ = float(nc.variables['sun_zenith'][0, midy, midx].data)
        saa_ = float(nc.variables['sun_azimuth'][0, midy, midx].data)

        vza_ = {b:float(nc.variables['view_zenith_{}'.format(b)][0, midy, midx].data) for b in metadata['BAND_NAMES'].split(',')}
        vaa_ = {b:float(nc.variables['view_azimuth_{}'.format(b)][0, midy, midx].data) for b in metadata['BAND_NAMES'].split(',')}

    ## use mean average of scene/region
    if geometry == 'mean':
        if sub is None:
            sza_ = np.nanmean(nc.variables['sun_zenith'][0, :, :].data)
            saa_ = np.nanmean(nc.variables['sun_azimuth'][0, :, :].data)

            vza_ = {b:np.nanmean(nc.variables['view_zenith_{}'.format(b)][0, :, :].data) for b in metadata['BAND_NAMES'].split(',')}
            vaa_ = {b:np.nanmean(nc.variables['view_azimuth_{}'.format(b)][0, :, :].data) for b in metadata['BAND_NAMES'].split(',')}
        else:
            sza_ = np.nanmean(nc.variables['sun_zenith'][0, sub[0]:sub[2], sub[1]:sub[3]].data)
            saa_ = np.nanmean(nc.variables['sun_azimuth'][0, sub[0]:sub[2], sub[1]:sub[3]].data)

            vza_ = {b:np.nanmean(nc.variables['view_zenith_{}'.format(b)][0, sub[0]:sub[2], sub[1]:sub[3]].data) for b in metadata['BAND_NAMES'].split(',')}
            vaa_ = {b:np.nanmean(nc.variables['view_azimuth_{}'.format(b)][0, sub[0]:sub[2], sub[1]:sub[3]].data) for b in metadata['BAND_NAMES'].split(',')}

    ## update metadata with sun and view geometry
    metadata['THS'] = sza_
    metadata['THV'] = np.nanmean([vza_[b] for b in vza_])

    sun_azi = saa_
    view_azi = np.nanmean([vaa_[b] for b in vaa_])

    azi = abs(sun_azi - view_azi)
    while(azi >= 180.):
        azi -= 180.

    metadata['AZI']=abs(azi)
    metadata['ISODATE']=gatts['PRODUCT_START_TIME']
    ## metadata finished

    ## get output directory
    if local_dir is not None:
        odir = local_dir
    else:
        odir = os.cwd()
    if not os.path.exists(odir): os.makedirs(odir)

    ## output file naming
    ncfile = '{}/{}_L1_converted.nc'.format(odir,metadata['OBASE'])

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
            dim = nc.variables[dsi].shape

            ## read dataset attributes
            atti = {attr : getattr(nc.variables[dsi],attr) for attr in nc.variables[dsi].ncattrs()}
            atto = None

            ## read
            if len(dim) == 2:
                if sub is None:
                    dat = nc.variables[dsi][:,:]
                else:
                    dat = nc.variables[dsi][sub[0]:sub[2], sub[1]:sub[3]]

            ## drop the first dimension for (1, 10980, 10980) data
            elif len(dim) == 3:
                if sub is None:
                    dat = nc.variables[dsi][0, :,:]
                else:
                    dat = nc.variables[dsi][0, sub[0]:sub[2], sub[1]:sub[3]]

            ## apply quantification value to convert to toa reflectance
            if ('Sentinel-2' in metadata['SATELLITE']) & ('rhot_' in r['out']):
                dat = dat.astype(float)/float(gatts['QUANTIFICATION_VALUE'])

                ## keep the attributes we need
                atto = {'units': 1, 'wavelength' : float(atti['wavelength'])}

            ac.output.nc_write(ncfile, dso, dat, attributes=metadata, dataset_attributes=atto, new=new)
            new = False ## appends to the NetCDF file

    ## close NetCDF file
    nc.close()

    ## return the path to the new file
    return(ncfile)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='NBS/OPeNDAP preprocessor for ACOLITE')
    parser.add_argument('--input', help='OPeNDAP URL to L1C Sentinel-2 data')
    parser.add_argument('--output', help='Output directory, if not provided will output to current working directory', default=None)
    parser.add_argument('--limit', help='4 element for cropping ROI in coordinates (default=None)', default=None)
    parser.add_argument('--sub', help='4 element for cropping ROI in pixels (default=None)', default=None)
    parser.add_argument('--geometry', help='Use single geometry for the scene, either "mid" or "mean" for scene/ROI (default=mid)', default="mid")
    parser.add_argument('--verbosity', help='Verbosity (default=0)', default=0)
    parser.add_argument('--acolite_path', help='Path to ACOLITE code (if not already in Python path) (default=None)', default=None)

    args, unknown = parser.parse_known_args()

    ## acolite path
    if args.acolite_path is not None:
        import sys
        sys.path.append(args.acolite_path)
    try:
        import acolite as ac
    except:
        print('Could not import acolite - make sure it is in your Python path, or provide it with --acolite_path')
        exit(1)

    ## parse given limit
    if args.limit is not None:
        sp = args.limit.split(',')
        if len(sp) == 4:
            args.limit = [float(s) for s in sp]
        else:
            args.limit = None

    ## parse given limit
    if args.sub is not None:
        sp = args.sub.split(',')
        if len(sp) == 4:
            args.sub = [int(s) for s in sp]
        else:
            args.sub = None

    if args.input is not None:
        if args.input[-5:] == '.html':
            args.input = args.input[0:-5]

        ## convert the remote file
        ret = convert_opendap_nbs(args.input,
                            local_dir=args.output,
                            limit=args.limit,
                            sub=args.sub,
                            geometry=args.geometry,
                            verbosity=int(args.verbosity))
