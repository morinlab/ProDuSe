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
        "scikit-learn==0.19.2",
        "scikit-bio"
        ],
    download_url="https://github.com/morinlab/ProDuSe/dist/ProDuSe-0.9.4.tar.gz",
    license="GNU GPLv3",
    scripts=["bin/produse"],
    data_files = [("ProDuSe", ["LICENSE.txt", "README.md", "etc/default_filter.pkl"])],
    zip_safe = False
)

