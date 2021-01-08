# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 11:43:27 2019

@author: oebbe
"""
import os

import hydropandas as hpd
from hydropandas.io import io_fews

dirname, fnames = hpd.util.get_files('./tests/data/2019-FEWS-test/WaalenBurg_201810-20190215_prod.zip', 
                                        ".xml")

def test_fews_highmemory():
    
    
    obs_list = io_fews.read_xml_fname(os.path.join(dirname,fnames[0]), 
                                      low_memory=False,
                                      ObsClass=hpd.WaterlvlObs)

    return obs_list


def test_fews_lowmemory():
    obs_list = io_fews.read_xml_fname(os.path.join(dirname,fnames[0]), 
                                      low_memory=True,
                                      ObsClass=hpd.WaterlvlObs)
   
    return obs_list


def test_obscollection_fews_selection():
    obs_list = io_fews.read_xml_filelist(fnames, hpd.WaterlvlObs, dirname)

    
    return obs_list

