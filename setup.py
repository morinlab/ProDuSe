#!/usr/bin/env python

from distutils.core import setup
from distutils.command.install import INSTALL_SCHEMES

for scheme in INSTALL_SCHEMES.values():
    scheme['data'] = scheme['purelib']

setup(
    name='ProDuSe',
    version="0.1.7",
    description='Process Duplex Sequence Data',
    author='Marco Albuquerque',
    author_email='malbuque@sfu.ca',
    maintainer='Christopher Rushton',
    maintainer_email='ckrushto@sfu.ca',
    packages=["ProDuSe"],
    url='https://github.com/morinlab/ProDuSe',
    classifiers=[
       "Programming Language :: Python :: 2.7"
       ],
    download_url="https://github.com/morinlab/ProDuSe/dist/ProDuse-0.1.7",
    license="GNU GPLv3",
    scripts=["bin/produse"],
    data_files = [("ProDuSe", ["LICENSE.txt", "README.md", "ProDuSe/splitmerge.pl", "ProDuSe/filter_produse.pl"])]
)

