"""Breeding Bird Survey Forecasting Benchmarks"""

import os

import pandas as pd
import numpy as np
import sqlalchemy
from sqlalchemy.engine import reflection

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

def mape(obs, pred):
    """Mean absolute percentage error

    Calculates the mean absolute percentage error = mean(abs(y_obs - y_pred)

    Args:
        obs: A list like object with a list of observed values
        pred: A list like object with a list of predicted values

    Returns:
        float = mean absolute percentage error
    """
    pred = np.array(pred)
    obs = np.array(obs)
    return np.mean(np.abs(pred-obs))