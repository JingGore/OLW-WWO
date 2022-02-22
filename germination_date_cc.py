# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 09:49:34 2020

@author: GuoJ
"""
import numpy as np
from os.path import join
import datetime as dt
# import pandas as pd
from os import walk
import os
import warnings
import xarray as xr


# Use julian day instead of sequential number from 0 to len(annual_file_list)
def get_germination_date(d_xarray, start_date, end_date, day_thre, pcp_thre, gmt_f):
    
    year_list = sorted(list(set(d_xarray.time.dt.year.values)))
    
    for year in year_list:
        print('Process {} at {}'.format(year, dt.datetime.now()))
        
        yearly_pcp = d_xarray.sel(time=slice('{}-{}'.format(year, start_date), '{}-{}'.format(year, end_date)))
        
        reaching_day = yearly_pcp[0].copy()
        reaching_day = reaching_day.where(reaching_day==0, other=0)
        
        
        
        for daily_pcp in yearly_pcp:
            
            day_of_year = int(daily_pcp.time.dt.dayofyear)
            
            period_pcp = yearly_pcp.sel(time=yearly_pcp.time.dt.dayofyear.isin(range(day_of_year, day_of_year+day_thre)))
            
            acc_array = daily_pcp.copy()
            acc_array = acc_array.where(acc_array==0, other=0)
            
            for pcp in period_pcp:
                
                current_day = int(pcp.time.dt.dayofyear)
                acc_array = acc_array + pcp
                
                mask = ((acc_array >= pcp_thre) & (reaching_day == 0))
                reaching_day = xr.where(mask, current_day, reaching_day)
                
        reaching_day = reaching_day.where(reaching_day>0, other=np.nan)
        reaching_day = reaching_day.expand_dims('time')
        # reaching_day.name='sowing_date'
        
        reaching_day.name = 'germination'
        reaching_day.attrs={'long_name': 'Germination Julian Date', 
              'standard_name': 'germination', 
              'units': '', 
              'description': 'Germination Julian Date (from 1 to 365/366)', 
              'missing': -9999.0}
        
        if year == year_list[0]:
            gmt = reaching_day
            
            
            
        else:
            gmt = xr.merge([gmt, reaching_day])
            
    # gmt.name = 'germination'
    # gmt.attrs={'long_name': 'Germination Julian Date', 
    #       'standard_name': 'germination', 
    #       'units': '', 
    #       'description': 'Germination Julian Date (from 1 to 365/366)', 
    #       'missing': -9999.0}
            
        # reaching_day.plot.imshow()
        
    gmt.to_netcdf(gmt_f)
    

def stats_analysis(self, in_dir, year_list):
    
    def get_tif_array_list(in_dir, year_list):

        array_list = []
    
        for (subdirpath, subdirname, filenames) in walk(in_dir):
            for f in filenames:
                # print(f.split('.')[0])
                # print(year_list)
                if (f.split('.')[-1].lower()[:3] == 'tif') and ((f.split('.')[0].split('_')[0] in year_list) or (f.split('.')[0] in year_list)):
                    # print(join(subdirpath, f))
                    array = self.raster.getRasterArray(join(subdirpath, f))
                    # array= np.where(np.isnan(self.ref_array), np.NaN, array)
                    array_list.append(array)
    
        return array_list
    
    def MeanMedianStdOfStack(array_list, vmin=True):

        array_stack = np.dstack(array_list)
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            value_mean = np.nanmean(array_stack, axis=2)
            
            value_std = np.nanstd(array_stack, axis=2)
            
            value_median = np.nanmedian(array_stack, axis=2)
            
            value_min = np.nanmin(array_stack, axis=2)
            
            value_max = np.nanmax(array_stack, axis=2)
            
            value_mean_plus_std = value_mean + value_std
            value_mean_minus_std = value_mean - value_std
        
        return value_mean, value_median, value_std, value_min, value_max, value_mean_plus_std, value_mean_minus_std
    
    array_list = get_tif_array_list(in_dir, year_list)
    
    mean_array, median_arran, std_array, min_array, max_array, mean_plus_std_array, mean_minus_std_array = MeanMedianStdOfStack(array_list)
    
    for array, fname in zip([mean_array, median_arran, std_array, min_array, max_array, mean_plus_std_array, mean_minus_std_array], ['mean', 'median', 'std', 'min', 'max', 'mean_plus_std', 'mean_minus_std']):
        array = np.around(array, 0)
        array = np.where(self.ref_array==self.no_data, self.no_data, array)
        self.raster.array2Raster(array, self.ref_raster, join(in_dir, 'germination_date_{}({}-{}).tif'.format(fname, year_list[0], year_list[-1])))
            
        
        
        
def main():
    
    data_dir = r'K:\climate'
    workspace = r'L:\gisdata\jing\climate\germination_sub-clover'

    
    if not os.path.exists(workspace):
        os.makedirs(workspace, exist_ok = True)
    
    # rcps = ['RCP2.6', 'RCP4.5', 'RCP6.0', 'RCP8.5']
    rcps = ['RCP8.5']
    gcms = ['BCC-CSM1.1', 'CESM1-CAM5', 'GFDL-CM3', 'GISS-EL-R', 'HadGEM2-ES', 'NorESM1-M']
    gcms = ['GISS-EL-R']
    years = ['2006-2010', '2011-2020', '2021-2030', '2031-2040', '2041-2050', '2051-2060', '2061-2070', '2071-2080', '2081-2090', '2091-2100']
    # years = ['1971-1980', '1981-1990', '1991-2000', '2001-2005']

    start_date = '03-01'
    end_date = '04-30'
    day_thre = 7
    pcp_thre = 20

    pcp_key = 'TotalPrecipCorr_VCSN'
    
    for rcp in rcps:
        for gcm in gcms:
            
            d_dir = join(data_dir, rcp, gcm)
            
            gmt_dir = join(workspace, rcp, gcm)
            if not os.path.exists(gmt_dir):
                os.makedirs(gmt_dir, exist_ok = True)
            
# =============================================================================
#             if rcp == 'RCP6.0' and gcm == 'HadGEM2-ES':
#                 years[-1] = '2091-2099'
#             else:
#                 years[-1] = '2091-2100'
# =============================================================================
            
            for y in years:
                pcp_fn = '{}_{}_{}_{}.nc'.format(pcp_key, gcm, y, rcp)
                pcp_f = join(d_dir, pcp_fn)
                
                gmt_file = join(gmt_dir, 'germination_sub-clover_{}_{}_{}.nc'.format(gcm, y, rcp))
                
                print('Processing {}  {}  {}  at  {}.'.format(rcp, gcm, y, dt.datetime.now()))
                
                with xr.open_dataset(pcp_f) as pcp_ds:
                    
                    get_germination_date(pcp_ds.rain, start_date, end_date, day_thre, pcp_thre, gmt_file)

    

    
    # print('Summarising... {} ...'.format(dt.datetime.now()))
    # sow.stats_analysis(processed_dir, year_list)

    print('Finished at {} ...'.format(dt.datetime.now()))
    
    
if __name__ == '__main__':
    main()    
