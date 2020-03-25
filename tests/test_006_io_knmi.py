# -*- coding: utf-8 -*-
"""
Created on Sat Dec 21 09:24:29 2019

@author: oebbe
"""

from observations.io import io_knmi


def test_download_rd_210():
    knmi_df, variables, stations = io_knmi.download_knmi_data(210, var='RD',
                                                              start='1952',
                                                              end=None,
                                                              interval='daily',
                                                              inseason=False,
                                                              verbose=False)
    return knmi_df, variables, stations


def test_download_ev24_240():
    knmi_df, variables, stations = io_knmi.download_knmi_data(210, var='EV24',
                                                              start='1952',
                                                              end=None,
                                                              interval='daily',
                                                              inseason=False,
                                                              verbose=False)
    return knmi_df, variables, stations


def test_fill_missing_measurements_ev24_278():
    knmi_df, variables, stations = io_knmi.fill_missing_measurements(278,
                                                                     var='EV24',
                                                                     start='1952',
                                                                     end=None,
                                                                     raise_exceptions=False,
                                                                     verbose=False)
    return knmi_df, variables, stations


def test_fill_missing_measurements_rd_278():
    knmi_df, variables, stations = io_knmi.fill_missing_measurements(892,
                                                                     var='RD',
                                                                     start='1952',
                                                                     end=None,
                                                                     raise_exceptions=False,
                                                                     verbose=True)
    return knmi_df, variables, stations
