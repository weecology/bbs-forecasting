"""Tests of time-series prediction benchmarking functions"""

from pandas import DataFrame
from pandas.util.testing import assert_frame_equal

from bbs_benchmarks import *

def test_benchmark_predictions():
    time = [1, 2, 3]
    value = [4, 5, 6]
    preds = benchmark_predictions(time, value, lag=1)
    assert preds == [6, 5, 4.5]

def test_filter_timeseries_contiguous():
    data = pd.DataFrame({'site': [1, 1, 1, 1, 2, 2], 'date': [1, 2, 3, 4, 1, 2]})
    filtered = filter_timeseries(data, group_cols='site', date_col='date', min_years=3)
    assert_frame_equal(filtered, pd.DataFrame({'site': [1, 1, 1, 1], 'date': [1, 2, 3, 4]}))

def test_filter_timeseries_noncontiguous_contigtrue():
    data = pd.DataFrame({'site': [1, 1, 1, 1, 2, 2, 2], 'date': [1, 2, 3, 4, 1, 2, 4]})
    filtered = filter_timeseries(data, group_cols='site', date_col='date', min_years=3)
    assert_frame_equal(filtered, pd.DataFrame({'site': [1, 1, 1, 1], 'date': [1, 2, 3, 4]}))

def test_filter_timeseries_noncontiguous_contigfalse():
    data = pd.DataFrame({'site': [1, 1, 1, 1, 2, 2, 2], 'date': [1, 2, 3, 4, 1, 2, 4]})
    filtered = filter_timeseries(data, group_cols='site', date_col='date', min_years=3, contiguous=False)
    assert_frame_equal(filtered, pd.DataFrame({'site': [1, 1, 1, 1, 2, 2, 2], 'date': [1, 2, 3, 4, 1, 2, 4]}))
