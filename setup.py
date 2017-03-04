#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name='ProDuSe',
    version="0.1.7",
    description='Process Duplex Sequence Data',
    author='Marco Albuquerque',
    author_email='malbuque@sfu.ca',
    maintainer='Christopher Rushton',
    maintainer_email='ckrushto@sfu.ca',
    include_package_data=True,
    # install_requires=["samtools >= 1.0"],
    packages=["ProDuSe"],
    url='https://github.com/morinlab/ProDuSe',
    classifiers=[
       "Programming Language :: Python :: 2.7"
       ],
    setup_requires=["numpy"],
    install_requires=[
        "configparser",
        "numpy",
        "configargparse",
        "pysam",
        "fisher"
        ],
    download_url="https://github.com/morinlab/ProDuSe/dist/ProDuse-0.1.7",
    license="GNU GPLv3",
    scripts=["bin/produse"],
    data_files = [("ProDuSe", ["LICENSE.txt", "README.md", "ProDuSe/splitmerge.pl", "ProDuSe/filter_produse.pl"])]
)

