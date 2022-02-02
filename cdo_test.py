#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 17 13:38:02 2021

@author: jing
"""

from cdo import Cdo
from os.path import join

cdo = Cdo()

indir = r'./VCSN'

# sinfon = cdo.sinfon(input=infile)
# cdo.showlevel(input=infile)

# in_files = '' 
# for y in range(1972, 2001):
#     in_files = in_files + join(indir, 'VCSN_Tmax5k_{}.nc '.format(y))

# print(in_files)

# out_file = './VCSN/VCSN_Tmax5k_1972-2000.nc'

# cdo.mergetime(input=in_files, output=out_file, options='-O')

cdo.chname('precipitation_amount,rain', input='./VCSN/drought/VCSN_Rainbc5k_1972-2000.nc', output='./VCSN/VCSN_Rainbc5k_1972-2000_new.nc', options='-O')

# infile = join(indir, r'VCSN_Tmean5k_1972-2020.nc')
# cdo.selyear('1972/2000', input=infile, output=join(indir, 'VCSN_Tmean5k_1972-2000.nc'))