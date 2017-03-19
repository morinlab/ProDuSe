#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name='ProDuSe',
    version="0.2.1",
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
        "numpy"
        ],
    download_url="https://github.com/morinlab/ProDuSe/dist/ProDuSe-0.18.0.tar.gz",
    license="GNU GPLv3",
    scripts=["bin/produse"],
    data_files = [("ProDuSe", ["LICENSE.txt", "README.md", "ProDuSe/splitmerge.pl", "ProDuSe/filter_produse.pl"])]
)

