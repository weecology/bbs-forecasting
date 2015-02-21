"""Breeding Bird Survey Forecasting Benchmarks"""

import os

import pandas as pd
import sqlalchemy
from sqlalchemy.engine import reflection

from macroecotools import richness_in_group, abundance_in_group

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
    return data.groupby(group_cols).filter(lambda x: len(np.unique(x[date_col])) >= min_years)


bbs_data = get_data('bbs')
bbs_data_timeseries = filter_timeseries(bbs_data, ['site_id'], 'year', 10)
richness = richness_in_group(bbs_data_timeseries, ['site_id', 'year'], ['species_id'])
