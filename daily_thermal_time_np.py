#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 16:20:27 2022

@author: jing
"""

from netCDF4 import Dataset, num2date
import numpy as np
from numpy import ma
from os.path import join

def get_period_fraction(period):
        
    # fraction of 3 hours temperature for the calculation of daily mean temperature based on daily Tmin and Tmax
    p_frac = 0.931 + 0.114 * period - 0.0703 * period**2 + 0.0053 * period**3
    
    return p_frac

def get_thermaltime(temperature, argv):
    '''
    A general fuction to calculate suitability index (0-1) based on fuzzy logic method. 
    array: the environment factor array.  E.g. GDD
    *argv: a list of threshold arguments of the factor to build the fuzzy curve. *Must be sequence of x
    '''
    if len(argv) == 4:
        
        def __EquationStraightLine__(x, x1, y1, x2, y2):
            
            def __ParamsStraightLine(x1, y1, x2, y2):
            
                if x2 - x1  == 0:
                    slope = None
                    intercept = None
                
                else:
                    slope = (y2-y1)/(x2-x1)
                    intercept = y1- slope * x1
                
                return slope, intercept
        
            
            slope, intercept = __ParamsStraightLine(x1, y1, x2, y2)
            
            if slope == None:
                return 0       # here is when line is x = b, so no value of y here. But we set y = 0 is because we use if to represent the unsuitable (suitability = 0)                                                                                                                                                          
            else:
                return slope * x + intercept

            
                       # The input is arrry

        tt = temperature.copy()
        tt = tt.where(tt==0, other=0)
            
        exp_1 = (temperature < argv[0][0])
        exp_2 = ((temperature >= argv[0][0]) & (temperature <= argv[1][0]))
        exp_3 = ((temperature >= argv[1][0]) & (temperature <= argv[2][0]))
        exp_4 = ((temperature >= argv[2][0]) & (temperature <= argv[3][0]))
        exp_5 = (temperature > argv[3][0])
            
        tt = np.where(exp_1, 0, 
                    np.where(exp_2, __EquationStraightLine__(temperature, argv[0][0], argv[0][1], argv[1][0], argv[1][1]),
                            np.where(exp_3, __EquationStraightLine__(temperature, argv[1][0], argv[1][1], argv[2][0], argv[2][1]),
                                     np.where(exp_4, __EquationStraightLine__(temperature, argv[2][0], argv[2][1], argv[3][0], argv[3][1]),
                                              np.where(exp_5, 0, tt)))))
            
        return tt


def main():
    
    # cdo = Cdo()
    no_data = -9999
    
    indir = r'./VCSN'
    outdir = r'./VCSN/drought'
    
    tmin_file = join(indir, 'VCSN_Tmin5k_1972-2000.nc') # 'MinTempCorr_VCSN_BCC-CSM1.1_2006-2010_RCP2.6.nc'
    tmax_file = join(indir, 'VCSN_Tmax5k_1972-2000.nc') # 'MaxTempCorr_VCSN_BCC-CSM1.1_2006-2010_RCP2.6.nc'
    out_file = join(outdir, 'VCSN_TT_1972-2000.nc')
    
    stageTT_dict = {1: 0, 15: 10, 30: 25, 40: 0}
    
    args = [[x, stageTT_dict[x]] for x in stageTT_dict]
    
    # args = []
    # for x in stageTT_dict:
    #     args.append([x, stageTT_dict[x]])
    
    # ds = xr.open_dataset(tmin_file)
    
    
    ncin = Dataset(tmin_file, 'r')
    
    time_var = ncin.variables['time']
    time_data = time_var[:]
    time_units = time_var.units
    tIn = num2date(time_data, time_units)
    
    data = ncin.variables['tmin'][:]
    tmin = ma.masked_values(data, no_data)
    
    ncin = Dataset(tmax_file, 'r')
    data = ncin.variables['tmax'][:]
    tmax = ma.masked_values(data, no_data)
    
    ncin.close()
    
    Absolute_zero = -273.15
    periods = 8
   
    pfs = [get_period_fraction(p) for p in range(1, periods+1)]
    
    tt = np.zeros(tmin.shape)
    nday = tIn.size
    for t in np.arrange(nday):
        tD = tIn[t]
        
        tminD = np.squeeze(tmin[t,:,:]) + Absolute_zero
        tmaxD = np.squeeze(tmax[t,:,:]) + Absolute_zero
        
        t_3hrs = [pf*(tmaxD-tminD) + tminD for pf in pfs]
        tt_3hrs = [get_thermaltime(t_3hr, stageTT_dict) for t_3hr in t_3hrs]
        ttD = sum(tt_3hrs)/periods

        tt[t,:,:] = ma.masked_array(ttD, tmin.mask)
    
if __name__ == '__main__':
    main()    