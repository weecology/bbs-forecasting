"""Breeding Bird Survey Forecasting Benchmarks"""

import os

import pandas as pd
import numpy as np
import sqlalchemy
from sqlalchemy.engine import reflection

from macroecotools import richness_in_group, abundance_in_group, obs_pred_rsquare

def database_exists(dataset):
    """Check to see if a dataset exists in the database"""
    insp = reflection.Inspector.from_engine(engine)
    return dataset in insp.get_schema_names()

def install_dataset(dataset):
    """Install a dataset using the EcoData Retriever"""
    os.system("retriever install postgres BBS -u postgres")

def get_data(dataset):
    """Get the data needed for the analysis"""
    engine = sqlalchemy.create_engine('postgresql+psycopg2://postgres@localhost/postgres')
    data_path = './data/{}_data.csv'.format(dataset)
    if os.path.exists(data_path):
        return pd.read_csv(data_path)
    else:
        if not database_exists('bbs'):
            install_dataset('bbs')
        #FIXME: This query doesn't currently deal with poorly sampled species (e.g., nocturnal)
        bbs_query = """SELECT (counts.statenum * 1000) + counts.route AS site_id, counts.year,
                       counts.aou AS species_id, counts.speciestotal AS abundance
                       FROM bbs.counts JOIN bbs.weather
                         ON bbs.counts.statenum=bbs.weather.statenum
                         AND bbs.counts.route=bbs.weather.route
                         AND bbs.counts.rpid=bbs.weather.rpid
                         AND bbs.counts.year=bbs.weather.year
                       WHERE bbs.weather.runtype=1 AND bbs.weather.rpid=101;
                    """
        bbs_data = pd.read_sql_query(bbs_query, engine)
        bbs_data.to_csv(data_path, index=False)
        return bbs_data

def filter_timeseries(data, group_cols, date_col, min_years):
    """Filter data to only include time-series with minimum number of years

    Args:
        data: A Pandas data frame one or more date related columns
        group_cols: A list of strings of names of columns to group by when
            counting dates (e.g,. a siteID column)
        date_col: The string name for the column holding the date
        min_years: Minimum number of years to be considered a time-series

    Returns:
        The original data frame filtered to remove time-series < min_years

    """
    #FIXME: Should support filtering to continuous time-series
    return data.groupby(group_cols).filter(lambda x: len(np.unique(x[date_col])) >= min_years)

def benchmark_predictions(time, value, lag=1):
    """Calculate benchmark predictions for a time-series

    Determines the value of the time-series at a given lag and the long-term
    average of the time-series for all points up to that lag to use as
    benchmark predictions for the last time-step.

    Args:
        time: A list like object with a list of times/dates
        value: A list like object with a list of the values at the associated time

    Returns:
        A list including the value at the final time-step, the value at the
        lag, and the average value for all time-steps from the start of the
        time-series up to and including the lag

    """
    assert lag < len(time), "Lag must be less than the length of the time-series"
    ts_data = pd.DataFrame({'time': time, 'value': value})
    ts_data = ts_data.sort('time')
    last_time_val = ts_data['value'].iloc[-1]
    pre_lag_vals = ts_data['value'].iloc[:-lag]
    lag_val = pre_lag_vals.iloc[-1]
    avg_val = np.mean(pre_lag_vals)
    return [last_time_val, lag_val, avg_val]

bbs_data = get_data('bbs')

#Initial Richness Analysis
bbs_comm_timeseries = filter_timeseries(bbs_data, ['site_id'], 'year', 10)
richness = richness_in_group(bbs_comm_timeseries, ['site_id', 'year'], ['species_id'])
richness_by_site = richness.groupby('site_id')

forecast_data = []
for site, site_data in richness_by_site:
    forecast_data.append([site] + benchmark_predictions(site_data['year'], site_data['richness']))

forecast_data = pd.DataFrame(forecast_data, columns=['site', 'last_yr_rich', 'prev_yr_rich', 'avg_rich'])
coefdet_avg_rich = obs_pred_rsquare(forecast_data['last_yr_rich'], forecast_data['avg_rich'])
coefdet_prev_yr_rich = obs_pred_rsquare(forecast_data['last_yr_rich'], forecast_data['prev_yr_rich'])

#Initial Population Abundance Analysis
bbs_pop_timeseries = filter_timeseries(bbs_data, ['site_id', 'species_id'], 'year', 10)
ab_by_pop = bbs_pop_timeseries.groupby(['site_id', 'species_id'])

forecast_data_pop = []
step = 1
for pop, pop_data in ab_by_pop:
    site_id, species_id = pop
    pop_data = pop_data.sort('year')
    last_yr_ab = pop_data['abundance'].iloc[-1]
    other_yrs_ab = pop_data['abundance'].iloc[:-1]
    prev_yr_ab = pop_data.sort('year')['abundance'].iloc[-2]
    avg_ab = np.mean(other_yrs_ab)
    forecast_data_pop.append([site_id, species_id, last_yr_ab, prev_yr_ab, avg_ab])
    step += 1
    if step % 1000 == 0:
        print(step / 200000.0)

forecast_data_pop = pd.DataFrame(forecast_data_pop, columns=['site_id', 'species_id', 'last_yr_ab', 'prev_yr_ab', 'avg_ab'])
coefdet_avg_ab = obs_pred_rsquare(forecast_data_pop['last_yr_ab'], forecast_data_pop['avg_ab'])
coefdet_prev_yr_ab = obs_pred_rsquare(forecast_data_pop['last_yr_ab'], forecast_data_pop['prev_yr_ab'])
