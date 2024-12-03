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
##               2023-11-06 (QV) allow separate estimation and correction
##               2023-11-09 (QV) added grid files and grid interpolation/fill
##               2023-11-13 (QV) update for multi point grid_files, add toa_min
##               2023-11-16 (QV) track band mask (assumed to be add_offset)
##               2024-02-28 (QV) added scaling option, i.e. to compute the per pixel ratio for roint and rhot at reference band
##                               changed minimum aot and wind (to allow valid LUT outputs)
##               2024-04-09 (QV) get platform from gatts if set, added support for aot and wind datasets in inputfile
##               2024-05-22 (QV) replace Sentinel-2 in sensor name with S2
##               2024-05-28 (QV) added aot_min setting via cli
##               2024-11-12 (QV) update for S2C
##               2024-12-03 (QV) added process_high_sza, sza_max parameter

def tgc(ncf, output=None, tgc_ext = 'TGC', method = 'T5', verbosity = 5, override = False,
        estimate = True, estimate_return = False, correct = True,
        wind_input = None, wind_default = 2.0, wind_min = 0.1,
        aot_input = None, aot_default = 0.1, aot_min = 0.01,
        process_high_sza = True, sza_max = 79.999,
        grid_files = None, grid_fill = True, grid_write = False, toa_min = 0.0001,
        lutdw = None, par = 'romix+rsky_t', base_luts = ['ACOLITE-LUT-202110-MOD2'], min_pixels = 500,
        scaling = False,
        ancillary_data=False, write_rhoi = False, reference_band = '11', glint_threshold = 0.02):
    import os, shutil, time, json
    import numpy as np
    import acolite as ac
    import scipy.optimize, scipy.ndimage
    from scipy.interpolate import LinearNDInterpolator
    import dateutil.parser, datetime

    if (estimate is False) & (correct is False):
        print('No estimate or correction required, returning None')
        return

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
    obase_json = '{bn}_{tgc_ext}_parameters.json'.format(bn=bn, tgc_ext=tgc_ext)
    ofile_json = '{output}/{obase}'.format(output=output if (output is not None) else dn, obase=obase_json)

    if not os.path.exists(os.path.dirname(ofile)):
        os.makedirs(os.path.dirname(ofile))

    ## read datasets and gatts
    datasets = ac.shared.nc_datasets(ncf)
    gatts = ac.shared.nc_gatts(ncf)

    ## get date and sensor
    date = gatts['start_date']
    dt = dateutil.parser.parse(date)
    isodate = dt.isoformat()[0:10]
    ftime = dt.hour + dt.minute/60 + dt.second/3600

    ## get platform QV 20240409
    platform = os.path.basename(ofile)[0:3]
    if 'platform' in gatts: platform = gatts['platform']
    sensor = '{}_MSI'.format(platform)
    if 'Sentinel-2' in sensor:
        sensor = sensor.replace('Sentinel-2', 'S2')

    ## get sensor settings (actually only for rsr_version)
    settings = ac.acolite.settings.parse(sensor, merge=False)
    ## read sensor rsrd for band names and wavelengths
    ## maybe not needed if we don't generate L1R?
    if settings['rsr_version'] is None:
        rsr_sensor = '{}'.format(sensor)
    else:
        rsr_sensor = '{}_{}'.format(sensor, settings['rsr_version'])
    rsrd = ac.shared.rsr_dict(rsr_sensor)[rsr_sensor]

    uoz = 0.3
    uwv = 1.5
    wind = 2
    pressure = 1013.25

    ## empty gem dict
    gem = {'data':{}, 'gatts':{}}

    ## read required datasets
    datasets_read = ['lat', 'lon', 'B11',
                     'view_zenith_mean', 'sun_zenith',
                     'view_azimuth_mean', 'sun_azimuth',]
    for ds in datasets_read: gem['data'][ds] = ac.shared.nc_data(ncf, ds)

    if len(np.where(np.isfinite(gem['data']['sun_azimuth']))[0]) == 0:
        print('No valid data in {}'.format(ncf))
        return

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

    if verbosity > 2: print('Used uoz, uwv, pressure: ', uoz, uwv, wind, pressure)

    gem['gatts']['wind'] = wind
    gem['gatts']['pressure'] = pressure
    gem['gatts']['sensor'] = sensor

    ## get mean average geometry
    geom_ds = ['sza', 'vza', 'raa', 'pressure', 'wind']
    geom_ds = [ 'view_zenith_mean', 'view_azimuth_mean', 'sun_zenith', 'sun_azimuth', 'pressure', 'wind']
    geom_mean = {k: np.nanmean(gem['data'][k]) if k in gem['data'] else gem['gatts'][k] for k in geom_ds}

    ## get gas transmittance
    tg_dict = ac.ac.gas_transmittance(geom_mean['sun_zenith'], geom_mean['view_zenith_mean'],
                                          uoz=uoz, uwv=uwv, sensor=rsr_sensor)

    ## make bands dataset
    gem['bands'] = {}
    for bi, b in enumerate(rsrd['rsr_bands']):
        if b not in gem['bands']:
            gem['bands'][b] = {k:rsrd[k][b] for k in ['wave_mu', 'wave_nm', 'wave_name'] if b in rsrd[k]}
            gem['bands'][b]['rhot_ds'] = 'rhot_{}'.format(gem['bands'][b]['wave_name'])
            gem['bands'][b]['rhos_ds'] = 'rhos_{}'.format(gem['bands'][b]['wave_name'])
            for k in tg_dict: ## add gas transmittance to bands dataset
                if k not in ['wave']: gem['bands'][b][k] = tg_dict[k][b]

    ## read LUT(s)
    if lutdw is None: lutdw = ac.aerlut.import_luts(add_rsky = True, par = par, sensor=rsr_sensor, base_luts=base_luts)
    luts = list(lutdw.keys())
    lut = luts[0]
    ## lut index and dimensions
    pid = lutdw[lut]['ipd'][par]
    pressures, pids, raas, vzas, szas, winds, aots = lutdw[lut]['dim']

    ## compute glint mask
    mask_glint = gem['data']['B11'] > glint_threshold

    ## is a grid of aot/wind used?
    grid = False

    ## estimate of wind speed and aot
    if estimate:
        ## read required datasets
        datasets_read = ['B3', 'B8',
                         #'B11',
                         #'lat', 'lon',
                         #'view_zenith_mean', 'sun_zenith',
                         #'view_azimuth_mean', 'sun_azimuth',
                         'view_zenith_B11', 'view_azimuth_B11']
        for ds in datasets_read: gem['data'][ds] = ac.shared.nc_data(ncf, ds)

        gem['gatts']['wind'] = wind
        gem['gatts']['pressure'] = pressure
        gem['gatts']['sensor'] = sensor

        ## dimensions and wavelengths
        dims = gem['data']['lon'].shape
        waves = [gem['bands'][b]['wave_nm'] for b in rsrd['rsr_bands']]

        ## simple glint mask
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
                print('Not performing TGC estimate for {}.'.format(ncf))
            return

        s = s0[0][ms[0]], s0[1][ms[0]]

        ## get average position
        lon_ave = np.nanmean(gem['data']['lon'][s[0],s[1]])
        lat_ave = np.nanmean(gem['data']['lat'][s[0],s[1]])

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
        ss = scipy.optimize.minimize(f_fit, [wind_min, aot_min], args=(reference_band, lut), bounds = [(wind_min,20), (aot_min, 5)])
        wind_fit, aot_fit = ss.x

        ## print fitting time
        t1f = time.time()
        if verbosity > 0: print('Fitting wind speed and aerosol took {:.1f} minutes'.format((t1f-t0)/60))
        if verbosity > 0: print('Fitted wind={:.2f}m/s aot={:.2f} npixels={:.0f}'.format(wind_fit, aot_fit, len(s[0])))

        ## return estimate
        parameters = {'lut': lut, 'wind': wind_fit, 'aot': aot_fit, 'lon': lon_ave, 'lat': lat_ave, 'npixels': len(s[0])}

        if estimate_return:
            with open(ofile_json, 'w') as f:
                json.dump(parameters, f)
            return(parameters)

    else:
        if (correct):
            if (grid_files) is not None:
                if type(grid_files) is not list:
                    grid_files = [grid_files]

                ## load the json data
                res = [json.load(open(f, 'r')) for f in grid_files]
                if verbosity > 0: print('Using wind/aot grid from {} files'.format(len(res)))

                ## extract the grid points
                lons = np.hstack([r['lon'] for r in res]).astype(np.float32).flatten()
                lats = np.hstack([r['lat'] for r in res]).astype(np.float32).flatten()
                aots = np.hstack([r['aot'] for r in res]).astype(np.float32).flatten()
                winds = np.hstack([r['wind'] for r in res]).astype(np.float32).flatten()

                if len(lons) < 3:
                    print('Data read for {} grid points'.format(len(lons)))
                    print('Exiting. Wind/aot grid computation required at least 3 grid points')
                    return

                if verbosity > 0: print('Using {} wind/aot grid points'.format(len(lons)))

                ifun = LinearNDInterpolator((lons, lats), aots)
                aot = ifun((gem['data']['lon'], gem['data']['lat']))
                aot[aot<aot_min] = aot_min
                if grid_fill: aot = ac.shared.fillnan(aot)

                ifun = LinearNDInterpolator((lons, lats), winds)
                wind = ifun((gem['data']['lon'], gem['data']['lat']))
                wind[wind<wind_min] = wind_min
                if grid_fill: wind = ac.shared.fillnan(wind)

                grid = True
            elif ('aot' in datasets) and ('wind' in datasets): ## read anc in dataset  QV 20240409
                if verbosity > 0: print('Reading aot and wind datasets from {}'.format(ncf))
                aot = ac.shared.nc_data(ncf, 'aot')
                aot[aot<aot_min] = aot_min
                wind = ac.shared.nc_data(ncf, 'wind')
                wind[wind<wind_min] = wind_min
                grid = True
            else:
                if (aot_input is None):
                    aot_fit = 1.0 * aot_default
                else:
                    aot_fit = 1.0 * aot_input

                if (wind_input is None):
                    wind_fit = 1.0 * wind_default
                else:
                    wind_fit = 1.0 * wind_input

                if verbosity > 0: print('Using external wind={:.2f}m/s aot={:.2f}'.format(wind_fit, aot_fit))

    ## make output file
    if correct:
        if os.path.exists(ofile) & (override is False):
            print('{} exists'.format(ofile))
            return
        if verbosity > 1: print('Processing to outputfile {}'.format(ofile))

        ## copy inputfile
        shutil.copy(ncf, ofile)

        ## read sun zenith angles
        sza = gem['data']['sun_zenith'] * 1.0
        if not process_high_sza:
            if np.nanmin(sza) > sza_max:
                correct = False
                if verbosity > 1: print('Copied input data to outputfile, but not performing TGC (SZA >= {:.2f})'.format(np.nanmin(sza)))
        else:
            sza[sza > sza_max] = sza_max ## replace with sza_max

    ## perform correction
    if correct:
        if (grid) & (grid_write):
            print('Wrote aot and wind to {}'.format(ofile))
            ac.output.nc_write(ofile, 'aot', aot)
            ac.output.nc_write(ofile, 'wind', wind)
            print('Wrote aot and wind to {}'.format(ofile))

        ## if scaling
        if scaling:
            saa = gem['data']['sun_azimuth'] * 1.0
            #sza = gem['data']['sun_zenith'] * 1.0
            #sza[sza > sza_max] = sza_max

            if 'view_zenith_B{}'.format(reference_band) in gem['data']:
                 vza_ = gem['data']['view_zenith_B{}'.format(reference_band)] * 1.0
                 vaa_ = gem['data']['view_azimuth_B{}'.format(reference_band)] * 1.0
                 ref_data = gem['data']['B{}'.format(reference_band)] * 1.0
            else:
                 vza_ = ac.shared.nc_data(ncf, 'view_zenith_B{}'.format(reference_band))
                 vaa_ = ac.shared.nc_data(ncf, 'view_azimuth_B{}'.format(reference_band))
                 ref_data = ac.shared.nc_data(ncf, 'B{}'.format(reference_band))

            #saa = ac.shared.nc_data(ncf, 'sun_azimuth')
            #sza = ac.shared.nc_data(ncf, 'sun_zenith')
            #vza_ = ac.shared.nc_data(ncf, 'view_zenith_B{}'.format(reference_band))
            #vaa_ = ac.shared.nc_data(ncf, 'view_azimuth_B{}'.format(reference_band))
            ref_data = ac.shared.nc_data(ncf, 'B{}'.format(reference_band))

            raa_ = np.abs(saa-vaa_)
            del saa, vaa_
            raa_[raa_>180] -= 360
            raa_ = np.abs(raa_)

            if grid:
                roint_ref = lutdw[lut]['rgi'][reference_band]((pressure, lutdw[lut]['ipd']['rsky_t'],
                                                 raa_, vza_, sza, wind, aot))
            else:
                print(wind_fit, aot_fit)
                roint_ref = lutdw[lut]['rgi'][reference_band]((pressure, lutdw[lut]['ipd']['rsky_t'],
                                                 raa_, vza_, sza, wind_fit, aot_fit))
            #roint_ratio = ref_data / roint_ref
            roint_ratio = roint_ref / ref_data  # QV 20240409
            roint_ratio[roint_ratio<1] = 1      # QV 20240409

            if write_rhoi:
                ac.output.nc_write(ofile, 'roint_ref', roint_ref)
                ac.output.nc_write(ofile, 'roint_ratio', roint_ratio)
            del roint_ref
        ## end scaling

        ## get geometry - mean geom not needed for T5
        #sza = gem['data']['sun_zenith']
        #sza[sza > sza_max] = sza_max
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
                if grid:
                    print(ds, 'computing band specific rhoi for wind and aot grid')
                else:
                    print(ds, 'computing band specific rhoi for wind={:.2f}m/s aot={:.2f}'.format(wind_fit, aot_fit))

                ## read band data
                d, a = ac.shared.nc_data(ncf, ds, attributes=True)
                if 'add_offset' in a:
                    data_mask = (d == a['add_offset'])
                else:
                    data_mask = d.mask

                ## band specific view geometry
                vza_ = ac.shared.nc_data(ncf, 'view_zenith_B{}'.format(band))
                vaa_ = ac.shared.nc_data(ncf, 'view_azimuth_B{}'.format(band))
                raa_ = np.abs(saa-vaa_)
                raa_[raa_>180] -= 360
                raa_ = np.abs(raa_)
                ## compute roint
                if grid:
                    roint_b = lutdw[lut]['rgi'][band]((pressure, lutdw[lut]['ipd']['rsky_t'],
                                                 raa_, vza_, sza, wind, aot))
                else:
                    roint_b = lutdw[lut]['rgi'][band]((pressure, lutdw[lut]['ipd']['rsky_t'],
                                                 raa_, vza_, sza, wind_fit, aot_fit))

                ## scaling
                #if scaling: roint_b *= roint_ratio
                if scaling: roint_b /= roint_ratio # QV 20240409

                ## write rhoi if asked
                if write_rhoi: ac.output.nc_write(ofile, 'rhoi_{}'.format(band), roint_b)
                ## remove rhoi from rhot
                d[mask_glint] -= roint_b[mask_glint]

                ## add back in rhoi for mean geometry
                ## not needed for T5
                if method == 'T4':
                    if grid:
                        print(ds, 'computing mean rhoi for wind and aot grid')
                        roint_a = lutdw[lut]['rgi'][band]((pressure, lutdw[lut]['ipd']['rsky_t'],
                                                           raa, vza, sza, wind, aot))
                    else:
                        print(ds, 'computing mean rhoi for wind={:.2f}m/s aot={:.2f}'.format(wind_fit, aot_fit))
                        roint_a = lutdw[lut]['rgi'][band]((pressure, lutdw[lut]['ipd']['rsky_t'],
                                                           raa, vza, sza, wind_fit, aot_fit))
                    ## write rhoi if asked
                    if write_rhoi: ac.output.nc_write(ofile, 'rhoi_{}_mean'.format(band), roint_a)
                    ## add mean geometry rhoi
                    d[mask_glint] += roint_a[mask_glint]
                ## fill values less than toa min
                d[d<toa_min] = toa_min

                ## add back the mask
                if 'add_offset' in a:
                    d[data_mask] = a['add_offset']
                else:
                    d.mask = data_mask

                ## write corrected band
                ac.output.nc_write(ofile, ds, d)
                tb1 = time.time()
                secb = tb1-tb0
                print('band {} took {:.1f} seconds'.format(band, secb))

        ## end time
        t1 = time.time()
        sec = t1-t0
        print('Total correction took {:.1f} minutes'.format(sec/60))

        return(ofile)

##
## main fun
if __name__ == '__main__':
    import argparse
    import sys, os
    import numpy as np

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

    parser.add_argument('--write_rhoi', help='Write rhoi data per band (default=False)', default=False)

    parser.add_argument('--estimate', help='Estimate wind speed and aot from scene (default=True)', default=True)
    parser.add_argument('--estimate_return', help='Return only estimate of wind speed and aot from scene (default=False)', default=False)
    parser.add_argument('--correct', help='Perform TGC (default=True)', default=True)
    parser.add_argument('--scaling', help='Scale the estimated roint to the reference band (default=False)', default=False)

    parser.add_argument('--wind', help='Wind speed(s) to use (default=None)', default=None)
    parser.add_argument('--aot', help='Aerosol optical depth(s) to use (default=None)', default=None)
    parser.add_argument('--lat', help='Latitude(s) for wind speed estimate (default=None)', default=None)
    parser.add_argument('--lon', help='Longitude(s) for wind speed estimate (default=None)', default=None)

    parser.add_argument('--grid_files', help='JSON files with wind and aerosol tie points (default=None)', default=None)
    parser.add_argument('--grid_fill', help='Fill grid outside tie points (default=True)', default=True)
    parser.add_argument('--grid_write', help='Write grid in output NetCDF file (default=False)', default=False)

    parser.add_argument('--toa_min', help='Minimum value at TOA after TGC (default=0.0001)', default=0.0001)
    parser.add_argument('--aot_min', help='Minimum aot (default=0.01)', default=0.01)

    parser.add_argument('--process_high_sza', help='Process data for sun zenith angles > 79.999 (default=True)', default=True)

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

    if type(args.estimate) == str: args.estimate = eval(args.estimate)
    if type(args.estimate_return) == str: args.estimate_return = eval(args.estimate_return)
    if type(args.correct) == str: args.correct = eval(args.correct)
    if type(args.scaling) == str: args.scaling = eval(args.scaling)
    if type(args.process_high_sza) == str: args.process_high_sza = eval(args.process_high_sza)

    ## if windspeed and aot are provided
    if type(args.wind) == str:
        args.wind = eval(args.wind)
        if args.wind is not None: args.wind = np.asarray(args.wind).astype(np.float32)
    if type(args.aot) == str:
        args.aot = eval(args.aot)
        if args.aot is not None: args.aot = np.asarray(args.aot).astype(np.float32)
    if type(args.lat) == str:
        args.lat = eval(args.lat)
        if args.lat is not None: args.lat = np.asarray(args.lat).astype(np.float32)
    if type(args.lon) == str:
        args.lon = eval(args.lon)
        if args.lon is not None: args.lon = np.asarray(args.lon).astype(np.float32)

    if type(args.grid_files) == str:
        args.grid_files = args.grid_files.split(',')
    if type(args.grid_fill) == str: args.grid_fill = eval(args.grid_fill)
    if type(args.grid_write) == str: args.grid_write = eval(args.grid_write)
    if type(args.write_rhoi) == str: args.write_rhoi = eval(args.write_rhoi)

    if args.input is not None:
        ret = tgc(args.input, output=args.output, override=args.override, verbosity=int(args.verbosity),
                  ancillary_data=args.ancillary_data, method=args.method, glint_threshold=float(args.glint_threshold),
                  min_pixels=int(args.min_pixels), write_rhoi = args.write_rhoi,
                  estimate = args.estimate, estimate_return = args.estimate_return, correct = args.correct,
                  wind_input = args.wind, aot_input = args.aot,
                  grid_files = args.grid_files, grid_fill = args.grid_fill, grid_write = args.grid_write,
                  toa_min = float(args.toa_min), scaling = args.scaling,
                  process_high_sza = args.process_high_sza, aot_min = float(args.aot_min),
                 )

        if (args.estimate) & (args.estimate_return):
            print(ret)
            sys.exit(ret)
