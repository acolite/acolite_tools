## tgc.py
## QV 2023-05-23 standalone TGC
##
## requires acolite to be in system path or provided as --acolite_path
## https://github.com/acolite/acolite
## note this is the "new" generic acolite from 2021
## the "old" acolite version has been moved to https://github.com/acolite/acolite_ls2
## package requirements of acolite should work
##
## example: python tgc.py --acolite_path=/storage/git/acolite --input=/storage/Projects/HROC/S2/S2RES/S2A_MSIL1C_20190905T113321_N0208_R080_T28RFT_20190905T115747_S2resampling.nc --output=/storage/Projects/HROC/S2/Standalone/
##
##
## def tgc
## QV 2023-04-27 "Lanzarote" TOA Glint Correction
## based on work done in November 2022-April 2023 for HROC
##
## last updates: 2023-10-30 (QV) added min_pixels keyword

def tgc(ncf, output=None, tgc_ext = 'TGC', method = 'T5', verbosity = 5, override = False,
        lutdw = None, par = 'romix+rsky_t', base_luts = ['ACOLITE-LUT-202110-MOD2'], min_pixels = 500,
        ancillary_data=False, write_rhoi = False, reference_band = '11', glint_threshold = 0.02):
    import os, shutil, time
    import numpy as np
    import acolite as ac
    import scipy.optimize, scipy.ndimage
    import dateutil.parser, datetime

    ## track time
    t0 = time.time()

    ## NDWI settings
    ndwi_green = 560
    ndwi_nir = 1600
    ndwi_threshold_gt = 0.2
    ndwi_threshold_lt = 0.9

    ## seed and fitting settings
    num = 5000
    np.random.seed(2022)
    min_pixels = max((1, min_pixels))

    ## function for fitting wind speed and aot
    def f_fit(x, band, lut):
        wind, aot = x
        romix = lutdw[lut]['rgi'][band]((pressure, lutdw[lut]['ipd'][par],
                                    raa_, vza_, sza, wind, aot))
        return(np.sqrt(np.nanmean(np.square((romix-data)))))

    if not os.path.exists(ncf):
        if verbosity > 0: print('{} does not exist'.format(ncf))
        return

    dn = os.path.dirname(ncf)
    bn = os.path.basename(ncf)
    bn, ext = os.path.splitext(bn)
    if tgc_ext is None: tgc_ext = '{}'.format(method)

    obase = '{bn}_{tgc_ext}.nc'.format(bn=bn, tgc_ext=tgc_ext)
    ofile = '{output}/{obase}'.format(output=output if (output is not None) else dn, obase=obase)

    if os.path.exists(ofile) & (override is False):
        print('{} exists'.format(ofile))
        return()
    print(ofile)

    ## read datasets and gatts
    datasets = ac.shared.nc_datasets(ncf)
    gatts = ac.shared.nc_gatts(ncf)

    ## get date and sensor
    date = gatts['start_date']
    dt = dateutil.parser.parse(date)
    isodate = dt.isoformat()[0:10]
    ftime = dt.hour + dt.minute/60 + dt.second/3600
    sensor = '{}_MSI'.format(os.path.basename(ofile)[0:3])

    ## read sensor rsrd for band names and wavelengths
    ## maybe not needed if we don't generate L1R?
    rsrd = ac.shared.rsr_dict(sensor)[sensor]

    ## read required datasets
    gem = {'data':{}, 'gatts':{}}
    datasets_read = ['B3', 'B8', 'B11',
                     'lat', 'lon', 'view_zenith_mean', 'view_azimuth_mean', 'sun_zenith', 'sun_azimuth',
                     'view_zenith_B11', 'view_azimuth_B11']
    for ds in datasets_read: gem['data'][ds] = ac.shared.nc_data(ncf, ds)

    uoz = 0.3
    uwv = 1.5
    wind = 2
    pressure = 1013.25

    ## try ancillary data
    if ancillary_data:
        clon = np.nanmedian(gem['data']['lon'])
        clat = np.nanmedian(gem['data']['lat'])
        anc = ac.ac.ancillary.get(isodate, clon, clat)

        try:
            uoz = anc['ozone']['interp']/1000. ## convert from MET data
            uwv = anc['p_water']['interp']/10. ## convert from MET data
            wind = ((anc['z_wind']['interp']**2) + (anc['m_wind']['interp']**2))**0.5
            pressure = gem['gatts']['pressure'] = anc['press']['interp']
        except:
            pass
    ## end ancillary data
    print(uoz, uwv, wind, pressure)
    gem['gatts']['wind'] = wind
    gem['gatts']['pressure'] = pressure
    gem['gatts']['sensor'] = sensor

    ## get mean average geometry
    geom_ds = ['sza', 'vza', 'raa', 'pressure', 'wind']
    geom_ds = [ 'view_zenith_mean', 'view_azimuth_mean', 'sun_zenith', 'sun_azimuth', 'pressure', 'wind']
    geom_mean = {k: np.nanmean(gem['data'][k]) if k in gem['data'] else gem['gatts'][k] for k in geom_ds}

    ## get gas transmittance
    tg_dict = ac.ac.gas_transmittance(geom_mean['sun_zenith'], geom_mean['view_zenith_mean'],
                                      uoz=uoz, uwv=uwv, sensor=gem['gatts']['sensor'])
    ## make bands dataset
    gem['bands'] = {}
    for bi, b in enumerate(rsrd['rsr_bands']):
        if b not in gem['bands']:
            gem['bands'][b] = {k:rsrd[k][b] for k in ['wave_mu', 'wave_nm', 'wave_name'] if b in rsrd[k]}
            gem['bands'][b]['rhot_ds'] = 'rhot_{}'.format(gem['bands'][b]['wave_name'])
            gem['bands'][b]['rhos_ds'] = 'rhos_{}'.format(gem['bands'][b]['wave_name'])
            for k in tg_dict:
                if k not in ['wave']: gem['bands'][b][k] = tg_dict[k][b]

    ## read LUT(s)
    if lutdw is None: lutdw = ac.aerlut.import_luts(add_rsky = True, par = par, sensor=gem['gatts']['sensor'], base_luts=base_luts)
    luts = list(lutdw.keys())
    lut = luts[0]

    ## lut index and dimensions
    pid = lutdw[lut]['ipd'][par]
    pressures, pids, raas, vzas, szas, winds, aots = lutdw[lut]['dim']

    ## dimensions and wavelengths
    dims = gem['data']['lon'].shape
    waves = [gem['bands'][b]['wave_nm'] for b in rsrd['rsr_bands']]

    ## simple glint mask
    mask_glint = gem['data']['B11'] > glint_threshold
    mask_nonglint = gem['data']['B11'] <= glint_threshold

    ## compute NDWI
    green = gem['data']['B3']
    nir = gem['data']['B8']
    ndwi = (green-nir)/(green+nir)

    ## compute and erode mask
    ndwi_mask = (ndwi>ndwi_threshold_gt) & (ndwi<ndwi_threshold_lt)
    ## add glint mask
    ndwi_mask[mask_nonglint] = 1
    struct = scipy.ndimage.generate_binary_structure(2, 2)
    ndwi_mask = scipy.ndimage.binary_erosion(ndwi_mask, struct, iterations=5)

    ## find random pixels
    s0 = np.random.choice(np.arange(dims[0]), num, replace=True),\
         np.random.choice(np.arange(dims[1]), num, replace=True)

    ## select pixels not masked
    m = ndwi_mask[s0[0], s0[1]].flatten() * 1.0
    m[m==0] = np.nan
    ms = np.where(ndwi_mask[s0[0], s0[1]].flatten())

    ## test whether enough pixels are available
    if len(ms[0]) < min_pixels:
        if verbosity > 1:
            print('Not enough pixels available ({}/{}) to estimate wind speed and aerosol optical depth.'.format(len(ms[0]), min_pixels))
            print('Not performing TGC for {}.'.format(ncf))
        return(ncf)

    s = s0[0][ms[0]], s0[1][ms[0]]

    ## read reference dataset
    ref_ds = 'B{}'.format(reference_band)
    ref_wv = ref_ds.split('_')[-1]
    data = gem['data'][ref_ds][s[0],s[1]].flatten()
    data /= gem['bands'][reference_band]['tt_gas']

    ## get sun geometry
    sza = gem['data']['sun_zenith'][s[0],s[1]].flatten()
    saa = gem['data']['sun_azimuth'][s[0],s[1]].flatten()
    ## get reference band geometry
    vaa_ = gem['data']['view_azimuth_B{}'.format(reference_band)][s[0],s[1]].flatten()
    raa_ = np.abs(saa-vaa_)
    raa_[raa_>180] -= 360
    raa_ = np.abs(raa_)
    vza_ = gem['data']['view_zenith_B{}'.format(reference_band)][s[0],s[1]].flatten()

    ## fit wind and aot
    ss = scipy.optimize.minimize(f_fit, [0.1, 0.01], args=(reference_band, lut), bounds = [(0.1,20), (0.01, 5)])
    wind_fit, aot_fit = ss.x
    print('Fitted wind={:.2f}m/s aot={:.2f}'.format(wind_fit, aot_fit))

    ## copy inputfile
    if not os.path.exists(os.path.dirname(ofile)):
        os.makedirs(os.path.dirname(ofile))
    shutil.copy(ncf, ofile)

    ## get geometry - mean geom not needed for T5
    sza = gem['data']['sun_zenith']
    saa = gem['data']['sun_azimuth']
    if method == 'T4':
        vza = gem['data']['view_zenith_mean']
        vaa = gem['data']['view_azimuth_mean']
        raa = np.abs(saa-vaa)
        raa[raa>180] -= 360
        raa = np.abs(raa)

    ## run through bands and correct for rhoi
    for band in rsrd['rsr_bands']:
        ds = 'B{}'.format(band)
        print(ds)

        ## correct if high transmittance
        if gem['bands'][band]['tt_gas'] > 0.7:
            tb0 = time.time()
            print(ds, 'computing band specific rhoi for wind={:.2f}m/s aot={:.2f}'.format(wind_fit, aot_fit))
            ## read band data
            d = ac.shared.nc_data(ncf, ds)
            ## band specific view geometry
            vza_ = ac.shared.nc_data(ncf, 'view_zenith_B{}'.format(band))
            vaa_ = ac.shared.nc_data(ncf, 'view_azimuth_B{}'.format(band))
            raa_ = np.abs(saa-vaa_)
            raa_[raa_>180] -= 360
            raa_ = np.abs(raa_)
            ## compute roint
            roint_b = lutdw[lut]['rgi'][band]((pressure, lutdw[lut]['ipd']['rsky_t'],
                                         raa_, vza_, sza, wind_fit, aot_fit))
            ## write rhoi if asked
            if write_rhoi: ac.output.nc_write(ofile, 'rhoi_{}'.format(wv), roint_b, new=new)
            ## remove rhoi from rhot
            d[mask_glint] -= roint_b[mask_glint]

            ## add back in rhoi for mean geometry
            ## not needed for T5
            if method == 'T4':
                print(ds, 'computing mean rhoi for wind={:.2f}m/s aot={:.2f}'.format(wind_fit, aot_fit))
                roint_a = lutdw[lut]['rgi'][band]((pressure, lutdw[lut]['ipd']['rsky_t'],
                                                   raa, vza, sza, wind_fit, aot_fit))
                ## write rhoi if asked
                if write_rhoi: ac.output.nc_write(ofile, 'rhoi_{}_mean'.format(wv), roint_a, new=new)
                ## add mean geometry rhoi
                d[mask_glint] += roint_a[mask_glint]
            ## add zeroes where negative
            d[d<0] = 0

            ## write corrected band
            ac.output.nc_write(ofile, ds, d)
            tb1 = time.time()
            secb = tb1-tb0
            print('band {} took {:.1f} seconds'.format(band, secb))
            #print('band {} took {:.1f} minutes'.format(band, secb/60))

    ## end time
    t1 = time.time()
    sec = t1-t0
    #print('{:.1f} seconds'.format(sec))
    print('Total correction took {:.1f} minutes'.format(sec/60))

    return(ofile)

##
## main fun
if __name__ == '__main__':
    import argparse
    import sys, os

    parser = argparse.ArgumentParser(description='TGC')
    parser.add_argument('--input', help='NetCDF file from S2Resampling processed L1C Sentinel-2 data')
    parser.add_argument('--output', help='Output directory, if not provided will output to same directory as input file', default=None)
    parser.add_argument('--verbosity', help='Verbosity (default=0)', default=0)
    parser.add_argument('--override', help='Overwrite NetCDF files (default=True)', default=True)
    parser.add_argument('--acolite_path', help='Path to ACOLITE code (if not already in Python path) (default=None)', default=None)
    parser.add_argument('--ancillary_data', help='Attempt to use ancillary data for gas transmittances (default=False)', default=False)

    parser.add_argument('--method', help='Method to use, T5 = compute TOA glint per band and remove, T4 = T5, but add band average geometry glint back in (default=T5)', default='T5')
    parser.add_argument('--glint_threshold', help='B11 TOA threshold above which to apply TGC (default=0.01)', default=0.01)
    parser.add_argument('--min_pixels', help='Minimum required pixels to estimate wind speed and aot (default=500)', default=500)

    args, unknown = parser.parse_known_args()

    ## acolite path
    if args.acolite_path is not None:
        sys.path.append(os.path.expanduser(args.acolite_path))
    try:
        import acolite as ac
        print('Imported acolite version {} from {}'.format(ac.version, ac.__path__[0]))
    except:
        print('Could not import acolite - make sure it is in your Python path, or provide it with --acolite_path')
        exit(1)

    ## parse output directory
    if args.output is not None:
        args.output = os.path.expanduser(args.output)

    ## evaluate bool flags
    if type(args.override) == str: args.override = eval(args.override)
    if type(args.ancillary_data) == str: args.ancillary_data = eval(args.ancillary_data)

    if args.input is not None:

        ## convert the remote file
        ret = tgc(args.input, output=args.output, override=args.override, verbosity=int(args.verbosity),
                  ancillary_data=args.ancillary_data, method=args.method, glint_threshold=float(args.glint_threshold),
                  min_pixels=int(args.min_pixels))
