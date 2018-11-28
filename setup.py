#!/usr/bin/env python

from setuptools import setup, find_packages
import re
import sys

# Imports version number
VERSIONFILE = "ProDuSe/__version.py"
verstrline = open(VERSIONFILE, "rt").read()
verRegex = r"^__version__ = ['\"]([^'\"]*)['\"]"
currentVer = re.search(verRegex, verstrline, re.M)
if currentVer:
    version = currentVer.group(1)
else:
    version = "Unknown"

# Since scikit-bio hasn't been updated in a while, the version of
# numpy it relies on is not compatible with python 3.7
# To minimize user confusion, check the version number
if sys.version_info >= (3, 7):
    raise NotImplementedError("Python 3.7 and higher are not currently supported")

setup(
    name='ProDuSe',
    version=version,
    description='Process Duplex Sequence Data',
    author='Christopher Rushton',
    author_email='ckrushto@sfu.ca',
    include_package_data=True,
    packages=["ProDuSe"],
    url='https://github.com/morinlab/ProDuSe',
    classifiers=[
       "Programming Language :: Python :: 3"
       ],
    setup_requires=["numpy"],
    install_requires=[  # This ordering is very important!!!!!!
        "fisher",
        "seaborn",
        "scipy",
        "sortedcontainers",
        "configobj",
        "pyfaidx",
        "pysam",
        "packaging",
        "scikit-learn==0.20.1", # To ensure the random forest filter works consistently
        "scikit-bio"
        ],
    download_url="https://github.com/morinlab/ProDuSe/dist/ProDuSe-0.9.4.tar.gz",
    license="GNU GPLv3",
    scripts=["bin/produse"],
    data_files = [("ProDuSe", ["LICENSE.txt", "README.md", "etc/default_filter.pkl"])],
    zip_safe = False
)

