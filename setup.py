#!/usr/bin/env python

#    setup.py: installation script for the MutateX package
#    Copyright (C) 2015, Matteo Tiberti <matteo.tiberti@gmail.com>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

from setuptools import setup, find_packages

setup(
    name="MutateX",
    version="0.8",
    packages=find_packages(),
    scripts=[ 'bin/ddg2density',
              'bin/ddg2distribution',
              'bin/pdb2labels',
              'bin/ddg2heatmap',
              'bin/ddg2excel',
              #'bin/ddg2csv',
              'bin/ddg2pca',
              'bin/mutatex',
              'bin/ddg2dg',
              'bin/ddg2histo',
              'bin/ddg2logo',
              'bin/ddg2pdb',
              'bin/ddg2summary'],

    install_requires=[ 'biopython',
                       'matplotlib==3.4.3',
		       'numpy==1.21.0',
                       'scipy==1.6.3',
                       'six',
                       'openpyxl',
                       'pyyaml'],
    package_data={
        '': ['LICENSE', '*.md'],
    },
    author="Matteo Tiberti, Thilde Terkelsen, Tycho Canter Cremers",
    author_email="matteo.tiberti@gmail.com",
    description="scripts and utilities to automate FoldX runs",
    keywords="foldx energy mutation",
    url="https://github.com/ELELAB/mutatex",
)

