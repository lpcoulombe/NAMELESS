#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  4 11:21:34 2026

@author: coulombe
"""


import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from NAMELESS import nameless

import time
import pdb
import yaml



def load_config(path):
    with open(path, "r") as f:
        config = yaml.safe_load(f)
    return config


def main():
    config_path = sys.argv[1]
    config = load_config(config_path)


    # initializing reduction object
    obs = config['observations']
    reduc = nameless.nameless(obs_folder=obs['obs_folder'],red_folder=obs['red_folder'])

    reduc.config = config

    stages = config['stages']
    
    s1_settings = config['s1_settings']
    s2_settings = config['s2_settings']
    
    prep_settings = config['prep_settings']
    cst_settings  = config['cst_settings']
    ext_settings  = config['ext_settings']
    
    reduc.run_reduction(run_s1=stages['run_s1'],run_s2=stages['run_s2'],
                        run_prep=stages['run_prep'],run_cst=stages['run_cst'],
                        run_ext=stages['run_ext'],s1_settings=s1_settings,
                        s2_settings=s2_settings,prep_settings=prep_settings,
                        cst_settings=cst_settings,ext_settings=ext_settings)
    
    
    print('Reduction complete!')

if __name__ == "__main__":
    main()

