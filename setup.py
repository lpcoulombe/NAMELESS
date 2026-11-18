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


# Files to download
FILES = {
    "opacity_table.dat":
        "https://raw.githubusercontent.com/USERNAME/REPOSITORY/main/my_package_data/opacity_table.dat",

    "chemical_equilibrium.dat":
        "https://raw.githubusercontent.com/USERNAME/REPOSITORY/main/my_package_data/chemical_equilibrium.dat",

    "stellar_spectra.dat":
        "https://raw.githubusercontent.com/USERNAME/REPOSITORY/main/my_package_data/stellar_spectra.dat",
}


def download_files():
    """Download required data files from GitHub."""

    for filename, url in FILES.items():

        destination = AUX_DIR / filename

        if destination.exists():
            print(f"{filename} already exists. Skipping.")
            continue

        print(f"Downloading {filename}...")

        urllib.request.urlretrieve(url, destination)

        print(f"Downloaded: {destination}")


class DownloadData(Command):
    """Custom command to download package data."""

    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        download_files()
    
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
          "jwst"],
      cmdclass={"download_data": DownloadData,})




