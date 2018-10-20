#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from distutils.core import setup, Extension
import numpy as np

ext_modules = [ Extension('louvainLNL', sources = ['louvainLNLmodule.cpp', 'community.cpp', 'graph_binary.cpp']) ]

setup(
        name = 'Louvain_LNL',
        version = '1.0',
        include_dirs = [np.get_include()], #Add Include path of numpy
        ext_modules = ext_modules
      )
