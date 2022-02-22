# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 14:02:02 2020

Determine date of first 50% flowering

@author: GuoJ
"""
import numpy as np
from os.path import join
from raster import Raster
# from file_cluster import GetFileCluster
import datetime as dt
import os
from os import walk
import xarray as xr

def get_nth_day(year, month, day):
    
    days_in_the_year = (dt.date(year, month, day)-dt.date(year, 1, 1)).days + 1
    
    return days_in_the_year

def is_leap_year(year):
    """Determine whether a year is a leap year."""

    return year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)



def flowering_every_year(ref_xarray, tt_xarray, start_date, end_date, tt_thre, flo_f, tt_xarray_2=None):
    
    
    def create_nc(ref_xarray, name, attrs={'long_name': '','standard_name': '', 'units': '', 'description': '', 'missing': -9999.0}):
        
        nc = xr.zeros_like(ref_xarray)
        nc.name = name
        nc.attrs=attrs
        
        return nc
    
    year_list = sorted(list(set(tt_xarray.time.dt.year.values)))
    
    r3_early = create_nc(ref_xarray, 'r3_early', {'long_name': 'First flowering Julian Date of Early cultivar', 
                                            'standard_name': 'flowering_r3_early', 
                                            'units': '', 
                                            'description': 'First flowering (50%) Julian Date (from 21st June) of early cultivar', 
                                            'missing': -9999.0})
    
    r6_early = create_nc(ref_xarray, 'r6_early', {'long_name': 'Pegging Julian Date of Early cultivar', 
                                            'standard_name': 'flowering_r6_early', 
                                            'units': '', 
                                            'description': 'pollinated flower start pegging Julian Date (from 21st June) of early cultivar', 
                                            'missing': -9999.0})
    
    r11_early = create_nc(ref_xarray, 'r11_early', {'long_name': 'Maturity Julian Date of Early cultivar', 
                                            'standard_name': 'flowering_r11_early', 
                                            'units': '', 
                                            'description': 'burr change from green to brown colour and starts to dry off Julian Date (from 21st June) of early cultivar', 
                                            'missing': -9999.0})
    
    r3_late = create_nc(ref_xarray, 'r3_late', {'long_name': 'First flowering Julian Date of Late cultivar', 
                                            'standard_name': 'flowering_r3_late', 
                                            'units': '', 
                                            'description': 'First flowering (50%) Julian Date (from 21st June) of late cultivar', 
                                            'missing': -9999.0})
    
    r6_late = create_nc(ref_xarray, 'r6_late', {'long_name': 'Pegging Julian Date of Late cultivar', 
                                            'standard_name': 'flowering_r6_late', 
                                            'units': '', 
                                            'description': 'pollinated flower start pegging Julian Date (from 21st June) of late cultivar', 
                                            'missing': -9999.0})
    
    r11_late = create_nc(ref_xarray, 'r11_late', {'long_name': 'Maturity Julian Date of Late cultivar', 
                                            'standard_name': 'flowering_r11_late', 
                                            'units': '', 
                                            'description': 'burr change from green to brown colour and starts to dry off Julian Date (from 21st June) of late cultivar', 
                                            'missing': -9999.0})
    
    pheno_dict = {'r3_early': r3_early, 'r6_early': r6_early, 'r11_early': r11_early, 'r3_late': r3_late, 'r6_lat': r6_late, 'r11_late': r11_late}
    
    del r3_early, r6_early, r11_early, r3_late, r6_late, r11_late
    
    
    try:
        if tt_xarray_2.any():
            
            first_year_of_array_2 = tt_xarray_2.sel(time=tt_xarray_2.time.dt.year.isin([list(set(tt_xarray_2.time.dt.year.values))[0]]))
            tt_xarray = xr.concat([tt_xarray, first_year_of_array_2], dim='time')
    
    except:
        
        print('tt_xarray_2 is None')
        year_list = year_list[:-1]
    
    for i, year in enumerate(year_list):
        print('Process {} at {}'.format(year, dt.datetime.now()))
        
        yearly_tt = tt_xarray.sel(time=slice('{}-{}'.format(year, start_date), '{}-{}'.format(year+1, end_date)))
        ref = ref_xarray.sel(time=ref_xarray.time.dt.year.isin([year]))
        
        reaching_day_dict = {'r3_early': '', 'r6_early': '', 'r11_early': '', 'r3_late': '', 'r6_lat': '', 'r11_late': ''}
        reaching_indicator_dict = {'r3_early': '', 'r6_early': '', 'r11_early': '', 'r3_late': '', 'r6_lat': '', 'r11_late': ''}
        
        # initialize array of each pheno stage 
        for pheno in reaching_day_dict:
        
            reaching_day_dict[pheno] = ref.copy() # the day when accumulated tt threshold reached
            reaching_day_dict[pheno] = reaching_day_dict[pheno].where(reaching_day_dict[pheno]==0, other=0)
            
            reaching_indicator_dict[pheno] = ref.copy()
            reaching_indicator_dict[pheno] = reaching_indicator_dict[pheno].where(reaching_indicator_dict[pheno]==0, other=0) 
        
        
        tt_accu = ref.copy()
        tt_accu = tt_accu.where(tt_accu==0, other=0) 
        
        for tt in yearly_tt:
            
            day_of_year = int(tt.time.dt.dayofyear)
            
            tt_accu = tt_accu + tt
            
            
            for pheno in reaching_day_dict:
                mask = (tt_accu>=tt_thre[pheno]) & (reaching_indicator_dict[pheno]==0)
                reaching_day_dict[pheno] = xr.where(mask, day_of_year, reaching_day_dict[pheno])
                reaching_indicator_dict[pheno] = xr.where(mask, 1, reaching_indicator_dict[pheno])
            
            
# =============================================================================
#             # r3 early
#             mask = (tt_accu>=tt_thre['r3_early']) & (reaching_indicator_dict['r3_early']==0)
#             reaching_day_dict['r3_early'] = xr.where(mask, day_of_year, reaching_day_dict['r3_early'])
#             reaching_indicator_dict['r3_early'] = xr.where(mask, 1, reaching_indicator_dict['r3_early'])
#             
#             # r6 early
#             mask = (tt_accu>=tt_thre['r6_early']) & (reaching_indicator_dict['r6_early']==0)
#             reaching_day_dict['r6_early'] = xr.where(mask, day_of_year, reaching_day_dict['r6_early'])
#             reaching_indicator_dict['r6_early'] = xr.where(mask, 1, reaching_indicator_dict['r6_early'])
#             
#             # r11 early
#             mask = (tt_accu>=tt_thre['r11_early']) & (reaching_indicator_dict['r11_early']==0)
#             reaching_day_dict['r11_early'] = xr.where(mask, day_of_year, reaching_day_dict['r11_early'])
#             reaching_indicator_dict['r11_early'] = xr.where(mask, 1, reaching_indicator_dict['r11_early'])
#             
#             # r3 late
#             mask = (tt_accu>=tt_thre['r3_late']) & (reaching_indicator_dict['r3_late']==0)
#             reaching_day_dict['r3_late'] = xr.where(mask, day_of_year, reaching_day_dict['r3_late'])
#             reaching_indicator_dict['r3_late'] = xr.where(mask, 1, reaching_indicator_dict['r3_late'])
#             
#             # r6 late
#             mask = (tt_accu>=tt_thre['r6_late']) & (reaching_indicator_dict['r6_late']==0)
#             reaching_day_dict['r6_late'] = xr.where(mask, day_of_year, reaching_day_dict['r6_late'])
#             reaching_indicator_dict['r6_late'] = xr.where(mask, 1, reaching_indicator_dict['r6_late'])
#             
#             # r11 late
#             mask = (tt_accu>=tt_thre['r11_late']) & (reaching_indicator_dict['r11_late']==0)
#             reaching_day_dict['r11_late'] = xr.where(mask, day_of_year, reaching_day_dict['r11_late'])
#             reaching_indicator_dict['r11_late'] = xr.where(mask, 1, reaching_indicator_dict['r11_late'])
# =============================================================================
            

        
        for pheno in pheno_dict:
            
            reaching_day_dict[pheno] = xr.where((reaching_day_dict[pheno]>0), reaching_day_dict[pheno], np.nan)
            pheno_dict[pheno].loc[ref.time] = reaching_day_dict[pheno]
            
    
    r_all = xr.merge([pheno_dict[p] for p in pheno_dict])
    r_all.to_netcdf(flo_f)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
def flowering(self, keyword, start_date, end_date, tt_thre, out_dir):
    
    tt_dict = self.__GetFileDictionary__(start_date, end_date, keyword)
    sowing_array = self.ref_array
    sowing_array = sowing_array

    for year in tt_dict:
        annual_file_list = tt_dict[year]
        
        reaching_indicator = np.zeros(self.ref_array.shape)
        flowering_day = np.zeros(self.ref_array.shape) # the day when accumulated rainfall threshold reached
        tt_accu = np.zeros(self.ref_array.shape)
        
        for tt_f in annual_file_list:
            
            tt_date = tt_f.split('\\')[-1].split('.')[0] 
            y = int(tt_date.split('-')[0])
            m = int(tt_date.split('-')[1])
            d = int(tt_date.split('-')[2])
            julian_d = get_nth_day(y, m, d)
            
            if y > int(year):  # move to next year, do not recount julian date
                if is_leap_year(int(year)):
                    julian_d = julian_d + 366
                else:
                    julian_d = julian_d + 365
            
            tt = self.raster.getRasterArray(tt_f)  
            
            tt_accu = np.where(tt_accu<tt_thre, tt_accu + tt, tt_accu)
            
            flowering_day = np.where(np.logical_and(tt_accu>=tt_thre, reaching_indicator==0), julian_d, flowering_day)
            
            reaching_indicator = np.where(np.logical_and(tt_accu>=tt_thre, reaching_indicator==0), 1, reaching_indicator)
            
            
        flowering_day = np.where(flowering_day==0, np.nan, flowering_day)
        flowering_day = np.where(self.ref_array==self.no_data, self.no_data, flowering_day)    
        self.raster.array2Raster(flowering_day, self.ref_raster, join(out_dir, '{}.tif'.format(year)))
        
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
        self.raster.array2Raster(array, self.ref_raster, join(in_dir, 'flowering_date_{}({}-{}).tif'.format(fname, year_list[0], year_list[-1])))

def main():
    
    
    data_dir = r'L:\gisdata\jing\climate'
    tt_dir = join(data_dir, 'tt_sub-clover')
    gmt_dir = join(data_dir, 'germination_sub-clover')
    workspace = r'L:\gisdata\jing\climate\flowering_sub-clover'

    
    if not os.path.exists(workspace):
        os.makedirs(workspace, exist_ok = True)
    
    # rcps = ['RCP2.6', 'RCP4.5', 'RCP6.0', 'RCP8.5']
    rcps = ['RCP6.0']
    gcms = ['BCC-CSM1.1', 'CESM1-CAM5', 'GFDL-CM3', 'GISS-EL-R', 'HadGEM2-ES', 'NorESM1-M']
    # gcms = ['BCC-CSM1.1']
    years = ['2006-2010', '2011-2020', '2021-2030', '2031-2040', '2041-2050', '2051-2060', '2061-2070', '2071-2080', '2081-2090', '2091-2100']
    # years = ['1971-1980', '1981-1990', '1991-2000', '2001-2005']

    start_date = '06-21'
    end_date = '06-20'
    tt_thre = {'r3_early': 520, 'r6_early': 690, 'r11_early': 890, 'r3_late': 870, 'r6_lat': 1040, 'r11_late': 1365}

    tt_key = 'TT_sub-clover'
    gmt_key = 'germination_sub-clover'
    
    for rcp in rcps:
        for gcm in gcms:
            
            d_tt_dir = join(tt_dir, rcp, gcm)
            d_gmt_dir = join(gmt_dir, rcp, gcm)
            
            flo_dir = join(workspace, rcp, gcm)
            if not os.path.exists(flo_dir):
                os.makedirs(flo_dir, exist_ok = True)
            
            if rcp == 'RCP6.0' and gcm == 'HadGEM2-ES':
                years[-1] = '2091-2099'
            else:
                years[-1] = '2091-2100'
            
            for i, y in enumerate(years):
                
                tt_fn = '{}_{}_{}_{}.nc'.format(tt_key, gcm, y, rcp)
                tt_f = join(d_tt_dir, tt_fn)
                
                gmt_fn = '{}_{}_{}_{}.nc'.format(gmt_key, gcm, y, rcp)
                gmt_f = join(d_gmt_dir, gmt_fn)
                
                flo_file = join(flo_dir, 'flowering_sub-clover_{}_{}_{}.nc'.format(gcm, y, rcp))
                
                if y != years[-1]:
                    tt_fn_2 = '{}_{}_{}_{}.nc'.format(tt_key, gcm, years[i+1], rcp)
                    tt_f_2 = join(d_tt_dir, tt_fn_2)
                
                else:
                    tt_f_2 = None
                
                
                print('Processing {}  {}  {}  at  {}.'.format(rcp, gcm, y, dt.datetime.now()))
                
                if tt_f_2 is None:
                
                    with xr.open_dataset(tt_f) as tt_ds, xr.open_dataset(gmt_f) as gmt_ds:
                        
                        flowering_every_year(gmt_ds.germination, tt_ds.tt, start_date, end_date, tt_thre, flo_file)

                else:
                    
                    with xr.open_dataset(tt_f) as tt_ds, xr.open_dataset(tt_f_2) as tt_ds_2, xr.open_dataset(gmt_f) as gmt_ds:
                        
                        flowering_every_year(gmt_ds.germination, tt_ds.tt, start_date, end_date, tt_thre, flo_file, tt_ds_2.tt)

    
    # print('Summarising... {} ...'.format(dt.datetime.now()))
    # sow.stats_analysis(processed_dir, year_list)

    print('Finished at {} ...'.format(dt.datetime.now()))
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
# =============================================================================
#     workspace = r'C:\Users\guoj\projects\337003-0037 Regen Hillconutry Lscape 2.2'
#     
#     sowing_dir = join(workspace, r'processed\sub_clover\germination')
#     
#     TT_dir = join(workspace, r'data\climate_VCSN')
#     
#     rst = Raster()
#     
#     year_list = [str(x) for x in range(1972, 2019)]
# 
#     start_date = '0621'
#     end_date = '0620'
# 
#     keyword = 'DailyThermalTime_whiteclover'
#     
#     flr = First50flowering(year_list, TT_dir)
#     
#     tt_thre_r3_early = 520
#     tt_thre_r3_late = 870
#     tt_thre_r6_early = 690
#     tt_thre_r6_late = 1040
#     tt_thre_r11_early = 890
#     tt_thre_r11_late = 1365
#     
#         
#     ref_rst = join(sowing_dir, r'germination_date_mean(1972-2017)_pcp_less_1200.tif')
#     flr.ref_raster = ref_rst
#     ref_array = rst.getRasterArray(ref_rst)
#     flr.ref_array = ref_array
#     
#     processed_dir = join(workspace, r'processed\sub_clover')
#     
#     #------------R3---------------------------
#     early_cultivar_dir = join(processed_dir, r'r3\early_cultivar')
#     if not os.path.exists(early_cultivar_dir):
#         os.makedirs(early_cultivar_dir, exist_ok = True)
#     
#     flr.flowering(keyword, start_date, end_date, tt_thre_r3_early, early_cultivar_dir)
#     
#     print('Summarising... {} ...'.format(dt.datetime.now()))
#     flr.stats_analysis(early_cultivar_dir, year_list[:-1])
#     
#     
#     late_cultivar_dir = join(processed_dir, r'r3\late_cultivar')
#     if not os.path.exists(late_cultivar_dir):
#         os.makedirs(late_cultivar_dir, exist_ok = True)
#     
#     flr.flowering(keyword, start_date, end_date, tt_thre_r3_late, late_cultivar_dir)
#     
#     print('Summarising... {} ...'.format(dt.datetime.now()))
#     flr.stats_analysis(late_cultivar_dir, year_list[:-1])
#     
#     #------------R6---------------------------
#     early_cultivar_dir = join(processed_dir, r'r6\early_cultivar')
#     if not os.path.exists(early_cultivar_dir):
#         os.makedirs(early_cultivar_dir, exist_ok = True)
#     
#     flr.flowering(keyword, start_date, end_date, tt_thre_r6_early, early_cultivar_dir)
#     
#     print('Summarising... {} ...'.format(dt.datetime.now()))
#     flr.stats_analysis(early_cultivar_dir, year_list[:-1])
#     
#     
#     late_cultivar_dir = join(processed_dir, r'r6\late_cultivar')
#     if not os.path.exists(late_cultivar_dir):
#         os.makedirs(late_cultivar_dir, exist_ok = True)
#     
#     flr.flowering(keyword, start_date, end_date, tt_thre_r6_late, late_cultivar_dir)
#     
#     print('Summarising... {} ...'.format(dt.datetime.now()))
#     flr.stats_analysis(late_cultivar_dir, year_list[:-1])
#     
#     #------------R11---------------------------
#     early_cultivar_dir = join(processed_dir, r'r11\early_cultivar')
#     if not os.path.exists(early_cultivar_dir):
#         os.makedirs(early_cultivar_dir, exist_ok = True)
#     
#     flr.flowering(keyword, start_date, end_date, tt_thre_r11_early, early_cultivar_dir)
#     
#     print('Summarising... {} ...'.format(dt.datetime.now()))
#     flr.stats_analysis(early_cultivar_dir, year_list[:-1])
#     
#     
#     late_cultivar_dir = join(processed_dir, r'r11\late_cultivar')
#     if not os.path.exists(late_cultivar_dir):
#         os.makedirs(late_cultivar_dir, exist_ok = True)
#     
#     flr.flowering(keyword, start_date, end_date, tt_thre_r11_late, late_cultivar_dir)
#     
#     print('Summarising... {} ...'.format(dt.datetime.now()))
#     flr.stats_analysis(late_cultivar_dir, year_list[:-1])
# =============================================================================

if __name__ == '__main__':
    main() 