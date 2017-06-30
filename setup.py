#!/usr/bin/env python

from setuptools import setup, find_packages
import re

# Imports version number
VERSIONFILE = "ProDuSe/__version.py"
verstrline = open(VERSIONFILE, "rt").read()
verRegex = r"^__version__ = ['\"]([^'\"]*)['\"]"
currentVer = re.search(verRegex, verstrline, re.M)
if currentVer:
    version = currentVer.group(1)
else:
    version = "Unknown"

setup(
    name='ProDuSe',
    version=version,
    description='Process Duplex Sequence Data',
    author='Marco Albuquerque',
    author_email='malbuque@sfu.ca',
    maintainer='Christopher Rushton',
    maintainer_email='ckrushto@sfu.ca',
    include_package_data=True,
    packages=["ProDuSe"],
    url='https://github.com/morinlab/ProDuSe',
    classifiers=[
       "Programming Language :: Python :: 2.7"
       ],
    setup_requires=["numpy"],
    install_requires=[
        "scipy",
        "configparser",
        "configargparse",
        "pysam",
        "numpy",
        "astropy",
	"matplotlib"
        ],
    download_url="https://github.com/morinlab/ProDuSe/dist/ProDuSe-0.2.4.tar.gz",
    license="GNU GPLv3",
    scripts=["bin/produse"],
    data_files = [("ProDuSe", ["LICENSE.txt", "README.md"])])

