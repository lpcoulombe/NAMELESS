#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 11:42:48 2026

@author: lcoulombe
"""

import os
import urllib.request
from pathlib import Path
from setuptools import setup, find_packages, Command


ROOT_DIR = Path(__file__).resolve().parent

# Define directories

DATA_DIR = ROOT_DIR / "Observations"
AUX_DIR  = ROOT_DIR / "Auxiliary"

for directory in [DATA_DIR, AUX_DIR]:
    directory.mkdir(parents=True, exist_ok=True)

    
# Package installation
setup(
      name="NAMELESS",
      version="0.1.0",
      packages=find_packages(),
      install_requires=[
          "numpy",
          "scipy",
          "matplotlib",
          "astropy",
          "jwst",
          "pastasoss"],
      cmdclass={"download_data": DownloadData,})




