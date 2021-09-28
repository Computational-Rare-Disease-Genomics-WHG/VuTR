#!/usr/bin/env python

"""
Setup script for the Python package
- Used for development setup with `pip install --editable .`
- Parsed by conda-build to extract version and metainfo
"""

import setuptools

PKG = 'utr_utils'

setuptools.setup(
    name=PKG,
    # This tag is automatically updated by bump2version
    version='1.0.0',
    description='Utility tools for the UTR application',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url=f'https://github.com/Computational-Rare-Disease-Genomics-WHG/UTR-Visualisation-App/',
    license='MIT',
    packages=[PKG],
    include_package_data=True,
    zip_safe=False,
    scripts=['my-python-project'],
    keywords='bioinformatics',
    classifiers=[
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)
