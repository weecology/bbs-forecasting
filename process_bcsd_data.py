import xarray as xr

#Process the original netCDF files downloaded from http://gdo-dcp.ucllnl.org/downscaled_cmip_projections/
#into an ensemble more easily used in R

def process_data(filename, new_filename, scale_precip=False):
    ds = xr.open_dataset(filename)

    #Hang on to these
    attributes = ds.attrs

    #Make a single time series ensemble
    #of all the different model runs
    ds = ds.mean('projection')

    #Longitude here is from 0 - 360, adjust to -180 - 180
    #Note that if any data was outside the western hemisphere this
    #would have to done differently
    adjusted_longitude = ds['longitude'].values - 360
    ds['longitude'].values = adjusted_longitude

    #This precip data is in mm/day, convert to total mm/month
    #to match PRISM
    if scale_precip:
        ds*=ds['time.days_in_month']

    ds.attrs = attributes
    ds.to_netcdf(new_filename)

process_data('data/cmip5/bcsd5_pr/Extraction_pr.nc', 'data/cmip5/cmip5_precip.nc', scale_precip=True)
process_data('data/cmip5/bcsd5_temp/Extraction_tas.nc', 'data/cmip5/cmip5_temp.nc')
