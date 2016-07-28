#!/usr/bin/env python

from distutils.core import setup

setup(
    name='ProDuSe',
    version="0.1.1",
    description='Process Duplex Sequence',
    author='Marco Albuquerque',
    author_email='malbuque@sfu.ca',
    packages=["ProDuSe"],
    url='https://github.com/morinlab/ProDuSe',
    classifiers=[
       "Programming Language :: Python :: 2.7"
       ],
    download_url="https://github.com/morinlab/ProDuSe/tarball/0.1.1",
    license="GNU GPLv3",
    scripts=["bin/produse"],
    data_files = ["LICENSE", "README.md"]
)

