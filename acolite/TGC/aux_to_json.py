## QV 2023-11-13
## converted to separate function 2024-02-27
##

def aux_to_json(bundle, output=None):
    import os, json
    import acolite as ac

    if output is None:
        output = os.path.dirname(bundle)

    ## parse bundle and get granule
    safe_files = ac.sentinel2.safe_test(bundle)
    granule = safe_files['granules'][0]

    ## get granule auxillary data
    data = ac.sentinel2.auxillary(bundle, granule)

    ## extract grid points
    winds = ((data['u10']['values']**2 + data['v10']['values']**2) ** 0.5).flatten()
    aots = (data['aod550']['values']).flatten()
    lats =  (data['aod550']['latitudes']).flatten()
    lons =  (data['aod550']['longitudes']).flatten()

    ## write parameter JSON file
    parameters = {'wind': [float(v) for v in winds],
                  'lon': [float(v) for v in lons],
                  'lat': [float(v) for v in lats],
                  'aot': [float(v) for v in aots]}
    ofile_json = '{}/{}_grid_anc_TGC_parameters.json'.format(output, os.path.basename(os.path.splitext(bundle)[0]))
    if not os.path.exists(output): os.makedirs(output)
    with open(ofile_json, 'w') as f:
        json.dump(parameters, f)

    return(ofile_json)
