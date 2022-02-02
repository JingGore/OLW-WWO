#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 13:54:35 2021

@author: jing
"""

from cdo import Cdo
from os.path import join

cdo = Cdo()

indir = r'./VCSN/drought'

f_et0 = join(indir, r'VCSN_ET0_1972-2000.nc')
f_kc = join(indir, r'VCSN_kc_maize_1972-2000.nc')
f_et0_kc = join(indir, r'VCSN_ET0_kc_1972-2000.nc')

f_etc = join(indir, r'VCSN_ETc_maize_1972-2000.nc')

# cdo.mul(input = "-et0 " + f_et0 + " -crop_coefficient_maize " +f_kc, output = f_etc)
cdo.merge(input=f_et0+" "+f_kc, output=f_et0_kc, options='-O')
cdo.expr("'etc=et0*kc'",input=f_et0_kc, output=f_etc, options="-f nc")
# cdo.expr("'etc=et0*kc'",input=f_et0+" "+f_kc, output=f_etc, options="-f nc")