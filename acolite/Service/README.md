## Note
This is an example service for processing Sentinel-2 and Landsat data for a given site (as defined by latitude, longitude and a box size). The acolite_service_run.sh script can be added to a crontab for automated recurring processing, and will activate the named conda environment before running the service script. The default acolite environment can be used, after additionally adding the "omnicloudmask" package if the cloud_mask option is to be used, e.g.

            conda activate acolite
            conda install -c conda-forge omnicloudmask

## Configuration
There are three configuration files that can be edited to configure the service. $HOME and $SITE parameters are replaced by the user home directory and the site name. $PATH is replaced with the acolite_tools/acolite/Service directory location.

acolite_service_config.toml Contains the local path to the main acolite distribution (a clone of https://github.com/acolite/acolite) and the path to the settings and site configurations (defaults are given below). Optionally, the site configuration toml file can be provided as the first argument on the command line.

acolite_service_settings.toml Example base settings for the processing. Multiple processing settings can be provided. In the example DSF, RAdCor, and TACT processing are given. TACT processing is skipped by default for Sentinel-2. Setting values starting with @ are copied from the target settings.

acolite_service_sites.toml Example settings for site and processing configuration. Can be configured either for a date or date range, or counting back a number of days from today. Settings starting with "acolite." will be passed along to ACOLITE in the processing settings, and override those in the acolite_service_settings.toml definition. For example activating the glint correction could be done by "acolite.dsf_residual_glint_correction = true" and disabling negatives masking by "acolite.l2w_mask_negative_rhow = false"
