#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  4 11:21:34 2026

@author: coulombe
"""


import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from NAMELESS import nameless

import time


# reduction object

reduc = nameless.nameless(obs_folder='WASP_39b_NIRISS')

reduc.run_reduction(run_s1=False,run_s2=False)
