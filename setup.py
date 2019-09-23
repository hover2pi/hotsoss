#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest', ]

setup(
    author="Joe Filippazzo",
    author_email='jfilippazzo@stsci.edu',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    description="Helpful Organizational Tools for SOSS",
    entry_points={},
    install_requires=['numpy', 'astropy', 'bokeh', 'scipy'],
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='hotsoss',
    name='hotsoss',
    packages=find_packages(include=['hotsoss']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/hover2pi/hotsoss',
    version='0.1.0',
    zip_safe=False,
)
