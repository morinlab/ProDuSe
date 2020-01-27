#!/usr/bin/env python

from setuptools import setup, find_packages
import re
import subprocess
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

# Read in the README for the long description
with open("README.md", "r") as fh:
    long_description = fh.read()

# Due to dependency shenanigans, install all dependencies through pip install of easy_install
# According to everything I have read, they should act the same, but they 100% do not
dependencyList = [
    "seaborn",
    "scipy",
    "fisher",
    "sortedcontainers",
    "configobj",
    "packaging",
    "pyfaidx",
    "pysam!=0.15.3",  # 0.15.3 does not compile
    "scikit-learn==0.22.1",
    "scikit-bio",
    "numpy"
]

# Only run this command if we are installing ProDuSe
if sys.argv[-1] == "install":
    installCom = ["pip", "install"]
    installCom.extend(dependencyList)
    subprocess.check_call(["pip", "install", "numpy"])  # WE ARE INSTALLING NUMPY FIRST BECAUSE ANGRY REASONS
    subprocess.check_call(installCom)

setup(
    name='ProDuSe',
    version=version,
    description='Error supression and variant calling pipeline for Illumina sequencing data',
    long_description=long_description,
    long_description_content_type="text/markdown",
    author='Christopher Rushton',
    author_email='ckrushto@sfu.ca',
    include_package_data=True,
    packages=["ProDuSe"],
    url='https://github.com/morinlab/Dellingr',
    classifiers=[
       "Programming Language :: Python :: 3",
       "Operating System :: Unix",
       "Topic :: Scientific/Engineering :: Bio-Informatics",
       "License :: OSI Approved :: GNU Affero General Public License v3"
       ],
    python_requires='>=3.5',
    install_requires=dependencyList,
    download_url="https://github.com/morinlab/ProDuSe/dist/ProDuSe-0.9.5.tar.gz",
    scripts=["bin/produse"],
    package_data = {"Dellingr": ["LICENSE.txt", "README.md", "etc/default_filter.pkl"]},
    zip_safe = False,
    project_urls={
        "Source": "https://github.com/morinlab/ProDuSe",
        "Documentation": "https://produse.readthedocs.io/en/latest/"
        }
)

