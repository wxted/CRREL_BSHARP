#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 07:50:40 2017

@author: Ted Letcher
"""
import os
import sys
import setuptools
import socket
import numpy.distutils.core
#Bootstrap a numpy installation before trying to import it.
import imp
try:
    imp.find_module('numpy')
except ImportError:
    import subprocess
    subprocess.call([sys.executable, '-m', 'pip', 'install', 'numpy'])

wrapper = numpy.distutils.core.Extension('ftn_funcs', sources=['src_ftn/fortran_funcs.f95'])

numpy.distutils.core.setup(
      name='bsharp',
      author="Ted Letcher",
      author_email="Theodore.W.Letcher@erdc.dren.mil",
      version='0.2.0',
      ext_modules=[wrapper],
      packages=['bsharp'],
      license='MIT',
      install_requires=['gdal','netCDF4','pygrib']
      )
