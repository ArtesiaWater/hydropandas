# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 11:43:27 2019

@author: oebbe
"""
import sys
sys.path.insert(1, "..")
import pandas as pd
from observations import observation as obs
from observations import io_wiski


def test_read_wiski_csv():
    # download single file

    header, data = io_wiski.read_wiski_file(r".\data\2019-WISKI-test\1016_PBF.csv",
                                            sep='\s+', header_sep=':',
                                            header_identifier=':', verbose=True,
                                            parse_dates={"datetime": [0, 1]},
                                            infer_datetime_format=True,
                                            index_col=["datetime"],
                                            translate_dic={'name':'Station Number', 
                                                           'x':'GlobalX',
                                                           'y':'GlobalY'})

    return header, data

#%%
  
def test_read_wiski_zip():
    obs_df = io_wiski.read_wiski_dir(r".\data\2019-WISKI-test\1016_PBF.zip",
                                     ObsClass=obs.GroundwaterObs,
                                     sep='\s+', header_sep=':', 
                                     header_identifier=':', 
                                     parse_dates={"datetime": [0,1]},
                                     index_col=["datetime"], 
                                     translate_dic={'name':'Station Number', 
                                                           'x':'GlobalX',
                                                           'y':'GlobalY'},
                                     verbose=True)
    
    return obs_df