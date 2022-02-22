# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 14:02:02 2020

Determine date of first trifoliate development

@author: GuoJ
"""
import numpy as np
from os.path import join
# from raster import Raster
import datetime as dt
import os
from os import walk
import xarray as xr

def get_nth_day(year, month, day):
    
    days_in_the_year = (dt.date(year, month, day)-dt.date(year, 1, 1)).days + 1
    
    return days_in_the_year


    
def trifoliate_every_year(tt_xarray, gmt_xarray, start_date, end_date, tt_thre, tt_thre_v6, tri_f):
    
    year_list = sorted(list(set(tt_xarray.time.dt.year.values)))
    
    v3 = xr.zeros_like(gmt_xarray)
    v3.name = 'v3'
    v3.attrs={'long_name': 'First trifoliate Julian Date', 
              'standard_name': 'trifoliate_v3', 
              'units': '', 
              'description': 'First trifoliate Julian Date (from germinations date to the end of June)', 
              'missing': -9999.0}
    
    v6 = xr.zeros_like(gmt_xarray)
    v6.name = 'v6'
    v6.attrs={'long_name': 'Fourth trifoliate Julian Date', 
              'standard_name': 'trifoliate_v6', 
              'units': '', 
              'description': 'Fourth trifoliate Julian Date (from germinations date to the end of June)', 
              'missing': -9999.0}
    
    for year in year_list:
        print('Process {} at {}'.format(year, dt.datetime.now()))
        
        yearly_tt = tt_xarray.sel(time=slice('{}-{}'.format(year, start_date), '{}-{}'.format(year, end_date)))
        gmt = gmt_xarray.sel(time=gmt_xarray.time.dt.year.isin([year]))
        
        reaching_day_v3 = gmt.copy() # the day when accumulated rainfall threshold reached
        reaching_day_v3 = reaching_day_v3.where(reaching_day_v3==0, other=0)
        
        reaching_indicator_v3 = gmt.copy()
        reaching_indicator_v3 = reaching_indicator_v3.where(reaching_indicator_v3==0, other=0) 
        
        reaching_day_v6 = gmt.copy() # the day when accumulated rainfall threshold reached
        reaching_day_v6 = reaching_day_v6.where(reaching_day_v6==0, other=0)
        
        reaching_indicator_v6 = gmt.copy()
        reaching_indicator_v6 = reaching_indicator_v6.where(reaching_indicator_v6==0, other=0) 
        
        tt_accu = gmt.copy()
        tt_accu = tt_accu.where(tt_accu==0, other=0) 
        
        for tt in yearly_tt:
            
            day_of_year = int(tt.time.dt.dayofyear)
            
            mask = (gmt<=day_of_year) 
            tt_accu = xr.where(mask, tt_accu + tt, tt_accu)
            
            # v3
            mask_v3 = (tt_accu>=tt_thre) & (reaching_indicator_v3==0)
            reaching_day_v3 = xr.where(mask_v3, day_of_year, reaching_day_v3)
            reaching_indicator_v3 = xr.where(mask_v3, 1, reaching_indicator_v3)
            
            # v6
            mask_v6 = (tt_accu>=tt_thre_v6) & (reaching_indicator_v6==0)
            reaching_day_v6 = xr.where(mask_v6, day_of_year, reaching_day_v6)
            reaching_indicator_v6 = xr.where(mask_v6, 1, reaching_indicator_v6)
        
        
        reaching_day_v3 = xr.where((reaching_day_v3>0), reaching_day_v3, np.nan)
        v3.loc[gmt.time] = reaching_day_v3
        
        reaching_day_v6 = xr.where((reaching_day_v6>0), reaching_day_v6, np.nan)
        v6.loc[gmt.time] = reaching_day_v6
        # if year == year_list[0]:
        #     v3_6 = xr.merge([reaching_day_v3, reaching_day_v6])

        # else:
        #     v3_6 = xr.merge([v3_6, reaching_day_v3, reaching_day_v6])
    
    v3 = xr.merge([v3, v6])
    v3.to_netcdf(tri_f)
   
        
        
def trifoliate_from_mean_sowing(self, keyword, start_date, end_date, tt_thre, out_dir):
    
    tt_dict = self.__GetFileDictionary__(start_date, end_date, keyword)
    sowing_array = self.ref_array
    sowing_array = sowing_array

    for year in tt_dict:
        annual_file_list = tt_dict[year]
        
        reaching_indicator = np.zeros(self.ref_array.shape)
        reaching_day = np.zeros(self.ref_array.shape) # the day when accumulated rainfall threshold reached
        tt_accu = np.zeros(self.ref_array.shape)
        
        for tt_f in annual_file_list:
            
            tt_date = tt_f.split('\\')[-1].split('.')[0] 
            julian_d = get_nth_day(int(tt_date.split('-')[0]), int(tt_date.split('-')[1]), int(tt_date.split('-')[2]))
            tt = self.raster.getRasterArray(tt_f)  
            
            tt_accu = np.where(np.logical_and(sowing_array<=julian_d, tt_accu<tt_thre), tt_accu + tt, tt_accu)
            
            reaching_day = np.where(np.logical_and(tt_accu>=tt_thre, reaching_indicator==0), julian_d, reaching_day)
            
            reaching_indicator = np.where(np.logical_and(tt_accu>=tt_thre, reaching_indicator==0), 1, reaching_indicator)
            
            
        reaching_day = np.where(reaching_day==0, np.nan, reaching_day)
        reaching_day = np.where(self.ref_array==self.no_data, self.no_data, reaching_day)    
        self.raster.array2Raster(reaching_day, self.ref_raster, join(out_dir, '{}.tif'.format(year)))
        
        print('Processing {} at {} ...'.format(year, dt.datetime.now()))
        
        
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
        self.raster.array2Raster(array, self.ref_raster, join(in_dir, 'trifoliate_date_{}({}-{}).tif'.format(fname, year_list[0], year_list[-1])))

def main():
    
    data_dir = r'L:\gisdata\jing\climate'
    tt_dir = join(data_dir, 'tt_sub-clover')
    gmt_dir = join(data_dir, 'germination_sub-clover')
    workspace = r'L:\gisdata\jing\climate\trifoliate_sub-clover'

    
    if not os.path.exists(workspace):
        os.makedirs(workspace, exist_ok = True)
    
    # rcps = ['RCP2.6', 'RCP4.5', 'RCP6.0', 'RCP8.5']
    rcps = ['RCPpast']
    # gcms = ['BCC-CSM1.1', 'CESM1-CAM5', 'GFDL-CM3', 'GISS-EL-R', 'HadGEM2-ES', 'NorESM1-M']
    gcms = ['CESM1-CAM5', 'GFDL-CM3', 'GISS-EL-R', 'HadGEM2-ES', 'NorESM1-M']
    # years = ['2006-2010', '2011-2020', '2021-2030', '2031-2040', '2041-2050', '2051-2060', '2061-2070', '2071-2080', '2081-2090', '2091-2100']
    years = ['1971-1980', '1981-1990', '1991-2000', '2001-2005']

    start_date = '03-01'
    end_date = '06-30'
    tt_thre = 230
    tt_thre_v6 = 436

    tt_key = 'TT_sub-clover'
    gmt_key = 'germination_sub-clover'
    
    for rcp in rcps:
        for gcm in gcms:
            
            d_tt_dir = join(tt_dir, rcp, gcm)
            d_gmt_dir = join(gmt_dir, rcp, gcm)
            
            tri_dir = join(workspace, rcp, gcm)
            if not os.path.exists(tri_dir):
                os.makedirs(tri_dir, exist_ok = True)
            
# =============================================================================
#             if rcp == 'RCP6.0' and gcm == 'HadGEM2-ES':
#                 years[-1] = '2091-2099'
#             else:
#                 years[-1] = '2091-2100'
# =============================================================================
            
            for y in years:
                tt_fn = '{}_{}_{}_{}.nc'.format(tt_key, gcm, y, rcp)
                tt_f = join(d_tt_dir, tt_fn)
                
                gmt_fn = '{}_{}_{}_{}.nc'.format(gmt_key, gcm, y, rcp)
                gmt_f = join(d_gmt_dir, gmt_fn)
                
                tri_file = join(tri_dir, 'trifoliate_sub-clover_{}_{}_{}.nc'.format(gcm, y, rcp))
                
                print('Processing {}  {}  {}  at  {}.'.format(rcp, gcm, y, dt.datetime.now()))
                
                with xr.open_dataset(tt_f) as tt_ds, xr.open_dataset(gmt_f) as gmt_ds:
                    
                    trifoliate_every_year(tt_ds.tt, gmt_ds.germination, start_date, end_date, tt_thre, tt_thre_v6, tri_file)

    

    
    # print('Summarising... {} ...'.format(dt.datetime.now()))
    # sow.stats_analysis(processed_dir, year_list)

    print('Finished at {} ...'.format(dt.datetime.now()))
    

if __name__ == '__main__':
    main() 