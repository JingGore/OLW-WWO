# -*- coding: utf-8 -*-
"""
Spyder Editor

Create 8*3hrs based daily mean temperature based on daily t-min and t-max
"""

import os
from os.path import join
import datetime as dt
import numpy as np
import xarray as xr
from cdo import Cdo



def daily_ThermalTime_3hrs(tmin_xarray, tmax_xarray, stageTT_dict, out_file):
    
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
                
            tt = xr.where(exp_1, 0, 
                        xr.where(exp_2, __EquationStraightLine__(temperature, argv[0][0], argv[0][1], argv[1][0], argv[1][1]),
                                xr.where(exp_3, __EquationStraightLine__(temperature, argv[1][0], argv[1][1], argv[2][0], argv[2][1]),
                                         xr.where(exp_4, __EquationStraightLine__(temperature, argv[2][0], argv[2][1], argv[3][0], argv[3][1]),
                                                  xr.where(exp_5, 0, tt)))))
                
            return tt
    
    
    Absolute_zero = -273.15
    periods = 8
   
    pfs = [get_period_fraction(p) for p in range(1, periods+1)]
    
    year_list = list(set(tmin_xarray.time.dt.year.values))
    data_Maarray = tmin_xarray[0].to_masked_array()
    mask_data = data_Maarray.mask
    
    tt = xr.zeros_like(tmin_xarray)
    
    tt.name = 'tt'
    tt.attrs={'long_name': 'Thermal Time', 
              'standard_name': 'tt', 
              'units': 'Degree Days', 
              'description': 'Thermal time over previous 24 hours', 
              'missing': -9999.0}

# =============================================================================
#     for pf in pfs:
#         
#         t_3hrs = pf*(tmax_xarray-tmin_xarray) + tmin_xarray
#         tt_3hrs = get_thermaltime(t_3hrs, stageTT_dict)
#         tt += tt_3hrs
#     
#     tt = tt/periods
#     tt = tt.where(~mask_data, other=np.nan)
#     
#     tt.name = 'tt'
#     tt.attrs={'long_name': 'Thermal Time', 
#               'standard_name': 'tt', 
#               'units': 'Degree Days', 
#               'description': 'Thermal time over previous 24 hours', 
#               'missing': -9999.0}
# =============================================================================
   
    
    year0 = year_list[0]-1
    for tmin, tmax in zip(tmin_xarray, tmax_xarray):
        
        year = int(tmin.time.dt.year.values)
        if year != year0:
            print('Process {} at {}'.format(year, dt.datetime.now()))
            year0 = year
        
        tmin = tmin + Absolute_zero
        tmax = tmax + Absolute_zero
        
        # tt_total = tmin.copy()
        # tt_total = tt_total.where(tt_total==0, other=0)
        
        t_3hrs = [pf*(tmax-tmin) + tmin for pf in pfs]
        tt_3hrs = [get_thermaltime(t_3hr, stageTT_dict) for t_3hr in t_3hrs]
        tt_daily = sum(tt_3hrs)/periods
        
        
# =============================================================================
#         for p in range(1, periods+1):
#             
#             pf = get_period_fraction(p)
#             
#             t_3hrs = pf*(tmax-tmin) + tmin
#             
#             tt_3hrs = get_thermaltime(t_3hrs, stageTT_dict)
#             
#             tt_total += tt_3hrs
#             
#         tt_daily = tt_total / periods
# =============================================================================
        
        tt_daily = tt_daily.where(~mask_data, other=np.nan)
# =============================================================================
#         tt_daily = tt_daily.expand_dims('time')
#         tt_daily.name='thermal_time'
#         
#         if tmin.time.dt.date == tmin_xarray[0].time.dt.date:
#             tt = tt_daily
#         else:
#             tt = xr.merge([tt, tt_daily])
# =============================================================================
        
        tt.loc[tmin.time] = tt_daily
        
        
    tt.to_netcdf(out_file)
        
        

def main():
    
    # cdo = Cdo()
    
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
    
    
    
    with xr.open_dataset(tmin_file) as tmin_ds, xr.open_dataset(tmax_file) as tmax_ds:
        
        daily_ThermalTime_3hrs(tmin_ds.tmin, tmax_ds.tmax, args, out_file)
        
    
    
if __name__ == '__main__':
    main()    