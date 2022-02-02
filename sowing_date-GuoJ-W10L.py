# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 09:49:34 2020

@author: GuoJ
"""
from file_cluster import GetFileCluster
import numpy as np
from os.path import join
import datetime as dt
# import pandas as pd
from os import walk
import os
import warnings

def get_nth_day(year, month, day):
    
    days_in_the_year = (dt.date(year, month, day)-dt.date(year, 1, 1)).days + 1
    
    return days_in_the_year

class FreqAccuPCP(GetFileCluster):
    
# =============================================================================
#     def freq_accu_pcp(self, keyword, start_date, end_date, day_thre, pcp_thre, out_dir):
#         
#         pcp_dict = self.__GetFileDictionary__(start_date, end_date, keyword)
# 
#         for year in pcp_dict:
#             annual_file_list = pcp_dict[year]
#             
#             reaching_day = np.zeros(self.ref_array.shape) # the day when accumulated rainfall threshold reached
#             
#             for i in range(0, len(annual_file_list)):
#                 
#                 # 'days' is the number of accumulated days of rainfall
#                 if i <= len(annual_file_list) - day_thre:
#                     days = day_thre
#                     count_end = 0
#                 else:
#                     days = len(annual_file_list) - i  # when remain days is less than day_thre, just take tha remain days for accumulation
#                     count_end = 1
#                 
#                 acc_array = np.zeros(self.ref_array.shape)
#                 
#                 for day in range(i,i+days):
#                     pcp_array = self.raster.getRasterArray(annual_file_list[day])                      
#                     acc_array = acc_array + pcp_array
#                 
#                 reaching_day = np.where(reaching_day>0, reaching_day, np.where(acc_array>=pcp_thre, i+days, reaching_day))
#                 
#                 if count_end == 1:
#                     break
# 
#             self.raster.array2Raster(reaching_day, self.ref_raster, join(out_dir, '{}.tif'.format(year)))
#             
#             print('Processing {} at {} ...'.format(year, dt.datetime.now()))
# =============================================================================
    
    # Use julian day instead of sequential number from 0 to len(annual_file_list)
    def freq_accu_pcp(self, keyword, start_date, end_date, day_thre, pcp_thre, out_dir):
        
        pcp_dict = self.__GetFileDictionary__(start_date, end_date, keyword)

        for year in pcp_dict:
            annual_file_list = pcp_dict[year]
            
            reaching_day = np.zeros(self.ref_array.shape) # the day when accumulated rainfall threshold reached
            
            first_sowing_date = annual_file_list[0].split('\\')[-1].split('.')[0]
            first_sowing_julian_day = get_nth_day(int(first_sowing_date.split('-')[0]), int(first_sowing_date.split('-')[1]), int(first_sowing_date.split('-')[2]))
            
            last_sowing_date = annual_file_list[-1].split('\\')[-1].split('.')[0]
            last_sowing_julian_day = get_nth_day(int(last_sowing_date.split('-')[0]), int(last_sowing_date.split('-')[1]), int(last_sowing_date.split('-')[2]))
            
            for pcp_f in annual_file_list:
                
                pcp_date = pcp_f.split('\\')[-1].split('.')[0] 
                julian_d = get_nth_day(int(pcp_date.split('-')[0]), int(pcp_date.split('-')[1]), int(pcp_date.split('-')[2]))
                
                # 'days' is the number of accumulated days of rainfall
                if julian_d <= last_sowing_julian_day - day_thre:
                    days = day_thre
                    count_end = 0
                elif julian_d == last_sowing_julian_day:  # the last possible germination date
                     days = 1
                     count_end = 1
                else:
                    days = last_sowing_julian_day - julian_d  # when remain days is less than day_thre, just take tha remain days for accumulation
                    count_end = 0
                
                acc_array = np.zeros(self.ref_array.shape)
                
                for day in range(julian_d-first_sowing_julian_day, julian_d-first_sowing_julian_day+days):
                    current_day = day + first_sowing_julian_day
                    pcp_array = self.raster.getRasterArray(annual_file_list[day])                      
                    acc_array = acc_array + pcp_array
                    reaching_day = np.where(reaching_day>0, reaching_day, np.where(acc_array>=pcp_thre, current_day, reaching_day))
                
                if count_end == 1:
                    break
                
            reaching_day = np.where(reaching_day==0, np.nan, reaching_day)
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
    
    workspace = r'C:\Users\guoj\projects\337003-0037 Regen Hillconutry Lscape 2.2'
    processed_dir = join(workspace, r'processed\sub_clover\germination')
    
    if not os.path.exists(processed_dir):
        os.makedirs(processed_dir, exist_ok = True)
    
    year_list = [str(x) for x in range(1972, 2018)]

    start_date = '0301'
    end_date = '0430'

    keyword = 'Rain5k'
        
    sow = FreqAccuPCP(year_list, join(workspace, r'data\climate_VCSN'))
    
    sow.freq_accu_pcp(keyword, start_date, end_date, 7, 20, processed_dir)
    
    print('Summarising... {} ...'.format(dt.datetime.now()))
    sow.stats_analysis(processed_dir, year_list)

    print('Finished at {} ...'.format(dt.datetime.now()))
    
    
if __name__ == '__main__':
    main()    
