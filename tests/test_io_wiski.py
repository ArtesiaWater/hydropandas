# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 11:43:27 2019

@author: oebbe
"""
import sys
sys.path.insert(1, "..")
from observations import io_wiski
from observations import observation as obs


def test_wiski_download_single():
    # download single file
    
    header, data = io_wiski.read_wiski_file(r".\data\2019-WISKI-test\1016_PBF.csv", verbose=True)
    
    return header, data




