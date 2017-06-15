import xarray as xr

#Process the original netCDF files downloaded from http://gdo-dcp.ucllnl.org/downscaled_cmip_projections/
#into an ensemble more easily used in R

def process_data(filename, new_filename):
    ds = xr.open_dataset(filename)

    #Make a single time series ensemble
    #of all the different model runs
    ds = ds.mean('projection')

    #Longitude here is from 0 - 360, adjust to -180 - 180
    #Note that if any data was outside the western hemisphere this
    #would have to done differently
    adjusted_longitude = ds['longitude'].values - 360
    ds['longitude'].values = adjusted_longitude

    ds.to_netcdf(new_filename)

process_data('data/cmip5/bcsd5_pr/Extraction_pr.nc', 'data/cmip5/cmip5_precip.nc')
process_data('data/cmip5/bcsd5_temp/Extraction_tas.nc', 'data/cmip5/cmip5_temp.nc')
