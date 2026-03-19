## acolite_service
## example NRT processing service
## written by Quinten Vanhellemont, RBINS
## 2026-03-05
## modifications: 2026-03-16 (QV) switch to TOML for configs, allow multiple sites
##                2026-03-17 (QV) changed scene location and deletion
##                2026-03-19 (QV) added passing of settings with acolite key, added rgb_range for ocm

def launch_service(acolite_path = None):
    import sys, os, glob, tomllib
    import datetime, dateutil.parser, shutil
    import numpy as np

    ## get user home
    user_home = os.path.expanduser("~")

    ## find file and load settings
    path = os.path.dirname(__file__)
    if path[-1] == '.': path = path[:-1]

    ## config files
    config_file = '{}/acolite_service_config.toml'.format(path)
    service_config = tomllib.load(open(config_file, 'rb'))
    for k in service_config:
        if type(service_config[k]) is str:
            if service_config[k].startswith('$HOME'):
                service_config[k] = '{}{}'.format(user_home, service_config[k][5:])
            if service_config[k].startswith('$PATH'):
                service_config[k] = '{}{}'.format(path, service_config[k][5:])

    # load second file as site config
    if len(sys.argv) >= 2:
       site_config_file = sys.argv[1]
    else:
       site_config_file = service_config['site_config']

    ## load site configuration
    site_config_dict = tomllib.load(open(site_config_file, 'rb'))
    for site in site_config_dict:
        site_name = site_config_dict[site].get('name')
        for k in site_config_dict[site]:
            if type(site_config_dict[site][k]) is str:
                ## replace $HOME with home directory
                if site_config_dict[site][k].startswith('$HOME'):
                    site_config_dict[site][k] = '{}{}'.format(user_home, site_config_dict[site][k][5:])
                ## replace $PATH with scripts directory
                if site_config_dict[site][k].startswith('$PATH'):
                    site_config_dict[site][k] = '{}{}'.format(path, site_config_dict[site][k][5:])
                ## replace $SITE with site name
                if '$SITE' in site_config_dict[site][k]:
                    site_config_dict[site][k] = site_config_dict[site][k].replace('$SITE', site_name)

    ## if acolite_path not given get it from config
    if (acolite_path is None) & ('acolite_path' in service_config):
        acolite_path = service_config['acolite_path']

    ## assume acolite is in same base path as acolite_tools
    if acolite_path is None:
        acolite_path = path[0:path.find('/acolite_tools')] + '/acolite'

    if acolite_path.startswith('$HOME'):
        acolite_path = '{}{}'.format(user_home, acolite_path[5:])

    sys.path.append(acolite_path)
    import acolite as ac

    ## open base processing settings
    settings_base = tomllib.load(open(service_config['settings_config'], 'rb'))

    ## configure processing settings
    run = True
    settings_sites = {}
    for site in site_config_dict:
        settings_sites[site] = {}
        for processor in settings_base:
            if processor not in site_config_dict[site]['processors']: continue

            ## copy keys from processor base settings
            settings_sites[site][processor] = {k: settings_base[processor][k] for k in settings_base[processor]}
            for k in settings_sites[site][processor]:
               if type(settings_sites[site][processor][k]) == str:
                   if settings_sites[site][processor][k].startswith('@'):
                       if k not in site_config_dict[site]:
                           print('Please configure {} in acolite_service_sites.toml'.format(k))
                           run = False
                       else:
                           settings_sites[site][processor][k] = site_config_dict[site][k]

            ## if acolite key is provided, copy settings
            if 'acolite' in site_config_dict[site]:
                for k in site_config_dict[site]['acolite']:
                    settings_sites[site][processor][k] = site_config_dict[site]['acolite'][k]

            ## replace placeholder values with value from configuration
            #for k in settings_base[processor]:
            #    if type(settings_base[processor][k]) == str:
            #        if settings_base[processor][k].startswith('@'):
            #            if k not in site_config:
            #                print('Please configure {} in acolite_service_config.txt'.format(k))
            #                run = False
            #            else:
            #                settings_base[processor][k] = site_config[k]


    if not run:
        print('Exiting.')
        sys.exit(1)

    ## inputfiles to delete
    inputfiles_to_delete = []

    ## run through sites in configuration
    for site in settings_sites:
        if not site_config_dict[site]['active']: continue

        ## figure out dates to run
        now = datetime.datetime.now()
        if 'datelist' in site_config_dict[site]:
            if type(site_config_dict[site]['datelist']) is not list:
                dates = [site_config_dict[site]['datelist']]
            else:
                dates = [d for d in site_config_dict[site]['datelist']]
        else:
            if 'end_date' not in site_config_dict[site]:
                end = datetime.datetime(now.year, now.month, now.day)
            else:
                end = dateutil.parser.parse(site_config_dict[site]['end_date'])

            if 'start_date' not in site_config_dict[site]:
                cur_date = end - datetime.timedelta(days = site_config_dict[site]['n_days'] - 1)
            else:
                start = dateutil.parser.parse(site_config_dict[site]['start_date'])
                cur_date = datetime.datetime(start.year, start.month, start.day)

            dates = []
            while (cur_date <= end):
                dates.append(cur_date.isoformat()[0:10])
                cur_date += datetime.timedelta(days = 1)
        print('Running for {} day{}, from {} to {}'.format(len(dates), '' if len(dates) == 1 else 's', dates[0], dates[-1]))

        roi = 'POINT ({} {})'.format(site_config_dict[site]['station_lon'], site_config_dict[site]['station_lat'])
        output = site_config_dict[site]['output']

        ## run through dates
        for date in dates:
            diff = now - datetime.datetime(int(date[0:4]), int(date[5:7]), int(date[8:10]))
            print(date, diff)

            ## folder to add data to
            date_folder = '{}'.format(date)
            date_folder = '{}/{}/{}'.format(date[0:4], date[5:7], date[8:10])

            for source in site_config_dict[site]['sources']:
                print('Running {} for {} on {}'.format(source, site, date))
                if source == 'Sentinel-2':
                    ret_ = ac.api.cdse.query(roi = roi, collection = "SENTINEL-2",
                                                           scene = None, start_date = date, end_date = date, attributes = True, verbosity = 5)
                    if ret_ is None:
                        print('Error retrieving Sentinel-2 scenes for {}'.format(date))
                        continue
                    urls, scenes, atts = ret_
                    ## combine scenes per sensor and relative orbit
                    scene_dict = {}
                    for si, sc in enumerate(scenes):
                        orb = sc[0:37]
                        if orb not in scene_dict: scene_dict[orb] = {'urls': [], 'scenes': [], 'atts': []}
                        scene_dict[orb]['urls'].append(urls[si])
                        scene_dict[orb]['scenes'].append(scenes[si])
                        scene_dict[orb]['atts'].append(atts[si])
                elif source == 'Landsat':
                    ret_ = ac.api.earthexplorer.query(roi = roi, landsat_type = 'ot',
                                                      start_date = date, end_date = date, verbosity = 5)
                    if ret_ is None:
                        print('Error retrieving Landsat scenes for {}'.format(date))
                        continue
                    entity_list, identifier_list, dataset_list = ret_

                    ## combine scenes per sensor
                    scene_dict = {}
                    for si, sc in enumerate(entity_list):
                        sp = identifier_list[si].split('_')
                        #orb = '{}_{}'.format(sp[0], sp[3]) ## just sensor and date, WRS could be merged
                        orb = '_'.join((sp[0], sp[1], sp[3])) ## just sensor and date, WRS could be merged
                        if orb not in scene_dict: scene_dict[orb] = {'scenes': [], 'entity_list': [], 'identifier_list': [], 'dataset_list': []}
                        scene_dict[orb]['entity_list'].append(entity_list[si])
                        scene_dict[orb]['identifier_list'].append(identifier_list[si])
                        scene_dict[orb]['dataset_list'].append(dataset_list[si])
                        scene_dict[orb]['scenes'].append(identifier_list[si])

                ## run through orbit combinations
                for orb in scene_dict:
                    ## run processors
                    for processor in settings_sites[site]:
                        ## check if running TACT
                        if (processor == 'TACT') | (settings_sites[site].get('tact_run') is True):
                            if (source not in service_config['tact_sources']):
                                print('Not running TACT for {}'.format(source))
                                continue
                            if (diff.total_seconds() < (86400 * service_config['tact_time_delay'])):
                                print('Skipping {} for {} (tact_time_delay={})'.format(processor, source, service_config['tact_time_delay']))
                                continue
                        settings = None
                        odir = '{}/{}/{}/{}'.format(output, date_folder, orb, processor)
                        print(orb, odir)

                        if (not os.path.exists(odir)) | (site_config_dict[site]['override']):
                            settings = {k: settings_sites[site][processor][k] for k in settings_sites[site][processor]}
                            settings['inputfile'] = []
                            settings['output'] = odir

                            ## download scene if we don't have it
                            for si, scene in enumerate(scene_dict[orb]['scenes']):
                                local_scene = '{}/{}/{}'.format(service_config['l1_scene_directory'], date, scene)

                                if not os.path.exists(local_scene):
                                    print('Downloading {}'.format(scene))
                                    if source == 'Sentinel-2':
                                        local_scene_ = ac.api.cdse.download(scene_dict[orb]['urls'][si],
                                                                            scenes = scene_dict[orb]['scenes'][si],
                                                                            output = '{}/{}'.format(service_config['l1_scene_directory'], date))
                                    elif source == 'Landsat':
                                        local_scene_ = ac.api.earthexplorer.download(scene_dict[orb]['entity_list'][si],
                                                                                     scene_dict[orb]['dataset_list'][si],
                                                                                     scene_dict[orb]['identifier_list'][si],
                                                                                     output = '{}/{}'.format(service_config['l1_scene_directory'], date))

                                    local_scene = local_scene_[0]
                                settings['inputfile'].append(local_scene)
                            settings['inputfile'].sort()

                            ## run processing
                            r = ac.acolite.acolite_run(settings)

                            ## l2r file
                            if 'l2r' in r[0]:
                                ncf = r[0]['l2r'][0]
                                ocm_dataset = 'rhos'
                            elif 'l1r' in r[0]:
                                ncf = r[0]['l1r'][0]
                                ocm_dataset = 'rhot'
                            else:
                                continue


                            rgb_range = [0, 0.15]
                            if 'rgb_min' in settings: rgb_range[0] = settings['rgb_min'][0]
                            if 'rgb_max' in settings: rgb_range[1] = settings['rgb_max'][0]

                            ## cloud masking
                            cm = ac.masking.ocm(ncf, export_maps = True, export_netcdf = True,
                                                rgb_range = rgb_range, dataset = ocm_dataset,
                                                masks = {1: {'name': 'Thick Cloud', 'color': 'red'},
                                                         2: {'name': 'Thin Cloud', 'color': 'orange'},
                                                         3: {'name': 'Cloud Shadow', 'color': 'yellow'}})

                            ## compute mask for current file
                            if site_config_dict[site]['cloud_mask']:
                                dn = os.path.dirname(ncf)
                                if site_config_dict[site].get('mask_target_wkt') is not None:
                                    ## create polygon
                                    poly = ac.shared.polygon_from_wkt(site_config_dict[site]['mask_target_wkt'], file='{}/polygon.geojson'.format(dn))
                                    ## get projection from NetCDF
                                    nc_projection = ac.shared.nc_projection_read(ncf)
                                    nc_dct = ac.shared.nc_projection_dct(nc_projection)
                                    ## create clip mask for polygon
                                    clip_mask = ac.shared.polygon_crop(nc_dct, poly, return_sub = False)
                                    ## compute mask
                                    m = cm[clip_mask == 255].flatten()
                                else:
                                    m = cm.flatten()

                                npix = len(m)
                                nmask = len(np.where(m>0)[0])
                                fmask = nmask/npix
                                ## write mask file
                                mask_file = dn + '/mask.txt'
                                with open(mask_file, 'w', encoding = 'utf-8') as f:
                                    f.write('{}'.format(str(fmask)))

                    ## append inputfile to delete list
                    #if settings is not None:
                    #    for local_scene in settings['inputfile']:
                    #        if (local_scene not in inputfiles_to_delete) & (os.path.exists(local_scene)) & (site_config_dict[site]['delete_l1']):
                    #                inputfiles_to_delete.append(local_scene)

        ## delete inputfiles after processing in case multiple regions use the same files
        #for local_scene in inputfiles_to_delete:
        #    if (os.path.exists(local_scene)):
        #        print('Deleting {}'.format(local_scene))
        #        #shutil.rmtree(local_scene)

    ## clean up l1 scenes - changed logic to store data per date
    ## here delete all scenes for that date if
    if service_config['delete_l1']:
        l1_dirs = glob.glob('{}/*'.format(service_config['l1_scene_directory']))
        #print(service_config['l1_scene_directory'], l1_dirs)
        for l1_dir in l1_dirs:
            l1_paths = glob.glob('{}/*'.format(l1_dir))
            l1_paths.sort()
            nscenes = len(l1_paths)
            l1_date = os.path.basename(l1_dir)
            diff = now - datetime.datetime(int(l1_date[0:4]), int(l1_date[5:7]), int(l1_date[8:10]))
            if (diff.total_seconds() > (86400 * service_config['delete_l1_delay'])):
                if nscenes > 0:
                    print('Deleting {} input scenes for date {} at {}'.format(nscenes, date, l1_dir))
                    for local_scene in l1_paths: shutil.rmtree(local_scene)
                    ## delete directory
                    shutil.rmtree(l1_dir)

if __name__ == '__main__':
    launch_service()
