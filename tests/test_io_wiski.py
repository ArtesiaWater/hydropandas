# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 11:43:27 2019

@author: oebbe
"""
import sys
sys.path.insert(1, "..")
from observations import io_wiski
from observations import observation as obs
import pandas as pd


def test_wiski_download_single():
    # download single file
    
    header, data = io_wiski.read_wiski_file(r".\data\2019-WISKI-test\1016_PBF.csv", 
                                            sep='\s+', headersep=':', 
                                            header_identifier=':', verbose=True,
                                            parse_dates={"datetime": [0,1]},
                                            infer_datetime_format=True,
                                            index_col=["datetime"])
    
    return header, data




def test_wiski_download_single2():
    # download single file
    
    header, data = io_wiski.read_wiski_file(r".\data\2019-WISKI-test\1016_PBF.csv", 
                                            sep='\s+', headersep=':', 
                                            header_identifier=':', verbose=True)
    
    data.index = pd.to_datetime(data['Date']+' '+data['Time'], format='%d-%m-%Y %H:%M:%S')
    data.drop(columns=['Date', 'Time'], inplace=True)
    

    return header, data




def test_wiski_download_single3():
    # download single file
    
    dateparse = lambda x: pd.datetime.strptime(x, '%d-%m-%Y %H:%M:%S')
    
    header, data = io_wiski.read_wiski_file(r".\data\2019-WISKI-test\1016_PBF.csv", 
                                            index_col=['datetime'],
                                            sep='\s+', headersep=':', 
                                            parse_dates={'datetime':[0,1]},
                                            date_parser=dateparse, infer_datetime_format=True,
                                            header_identifier=':', verbose=True)
    
    return header, data
