#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 12:26:20 2021

@author: jing
"""
from cdo import Cdo
import datetime as dt
import xarray as xr
import numpy as np
from os.path import join

def get_nth_day(year, month, day):
    
    days_in_the_year = (dt.date(year, month, day)-dt.date(year, 1, 1)).days + 1
    
    return days_in_the_year


def date_initial_stage(xarray_tmean, out_dir):
    
    year_list = list(set(xarray_tmean.time.dt.year.values))
    day_thre = 7
    temp_thre = 273.15 + 13
    
    
    for year in year_list:
        print('Process {} at {}'.format(year, dt.datetime.now()))
        
        year_tmean = xarray_tmean.sel(time=xarray_tmean.time.dt.year.isin([year]))
        
        reaching_day = year_tmean[0].copy()
        reaching_day = reaching_day.where(reaching_day==0, other=0)

        
        for day_tmean in year_tmean:
            
            day_of_year = int(day_tmean.time.dt.dayofyear)
            
            period_tmean = year_tmean.sel(time=year_tmean.time.dt.dayofyear.isin(range(day_of_year, day_of_year+day_thre)))
            
            if period_tmean.sizes['time'] == day_thre:
            
                mean_period_tmean = period_tmean.mean(dim='time', skipna=True, keep_attrs=True)
                mask = ((mean_period_tmean > temp_thre) & (reaching_day == 0))
                reaching_day = xr.where(mask, day_of_year+day_thre, reaching_day)   
                
        reaching_day = reaching_day.where(reaching_day>0, other=np.nan)
        reaching_day = reaching_day.expand_dims('time')
        reaching_day.name='sowing_date'
        
        if year == year_list[0]:
            sowing_day = reaching_day
        else:
            sowing_day = xr.merge([sowing_day, reaching_day])
            
        # reaching_day.plot.imshow()
        
    sowing_day.to_netcdf(join(out_dir, 'VCSN_{}_{}-{}.nc'.format(reaching_day.name, year_list[0], year_list[-1])))
    

def code_pheno_stage(pheno_stage_ds, tt_ds, out_dir):
    
    xarray_tt = tt_ds.tt
    start_date = '09-01'
    end_date = '08-31'
    
    tt_end = 1300
    tt_dev = 0.16 * tt_end
    tt_mid = 0.28 * tt_end + tt_dev
    tt_late = 0.33 * tt_end + tt_mid
    
    year_list = list(set(xarray_tt.time.dt.year.values))
    
    data_Maarray = xarray_tt[0].to_masked_array()
    mask_data = data_Maarray.mask
    
    stege_code = xr.zeros_like(xarray_tt)
    
    stege_code.name = 'psc'
    stege_code.attrs={'long_name': 'Phenology_stage_code', 
              'standard_name': 'psc', 
              'units': '', 
              'description': 'Phenology stage code. 1: Initial; 2: Development; 3: Mid; 4: Late', 
              'missing': -9999.0}
    
    for year in year_list:
        
        print('Process {} at {}'.format(year, dt.datetime.now()))
            
        stage_init = pheno_stage_ds.sowing_date.sel(time=pheno_stage_ds.sowing_date.time.dt.year.isin([year]))[0]
        
        if year == year_list[0]:
            yearly_tt = xarray_tt.sel(time=slice('{}-{}'.format(year, '01-01'), '{}-{}'.format(year+1, end_date)))
        else:
            yearly_tt = xarray_tt.sel(time=slice('{}-{}'.format(year, start_date), '{}-{}'.format(year+1, end_date)))
        
        # stage_init = stage_init.squeeze('time')
        
        accru_tt = stage_init.copy()
        accru_tt = accru_tt.where(accru_tt==0, other=0)
        
        day_dev = stage_init.copy()
        day_dev = day_dev.where(day_dev==0, other=0)
        
        day_mid = stage_init.copy()
        day_mid = day_mid.where(day_mid==0, other=0)
        
        day_late = stage_init.copy()
        day_late = day_late.where(day_late==0, other=0)
        
        day_end = stage_init.copy()
        day_end = day_end.where(day_end==0, other=0)
        
        day_stages = [day_dev, day_mid, day_late, day_end]
        stage_names = ['development_date', 'mid_stage_date','late_stage_date', 'end_stage_date']
        stage_codes = [1,2,3,4]
        
        for daily_tt in yearly_tt: 
            
            daily_stege_code = daily_tt.copy()
            daily_stege_code = daily_stege_code.where(daily_stege_code==0, other=0)
            
            day_of_year = int(daily_tt.time.dt.dayofyear)
            
            y = int(daily_tt.time.dt.year.values)
            if y != year:
                last_date_of_last_year = dt.date.fromisoformat('{}-12-31'.format(year))
                last_day_of_last_year = last_date_of_last_year.timetuple().tm_yday
                day_of_year = last_day_of_last_year + day_of_year


            mask_start_accru_tt = (stage_init < day_of_year)
            
            accru_tt = xr.where(mask_start_accru_tt, accru_tt + daily_tt, accru_tt)
            # accru_tt['time'] = daily_tt.time
            
            mask_before_init = (stage_init > day_of_year)
            mask_init = ((stage_init < day_of_year) & (accru_tt < tt_dev))
            mask_dev = ((stage_init < day_of_year) & (accru_tt >= tt_dev) & (accru_tt < tt_mid))
            mask_mid = ((stage_init < day_of_year) & (accru_tt >= tt_mid) & (accru_tt < tt_late))
            mask_late = ((stage_init < day_of_year) & (accru_tt >= tt_late) & (accru_tt <= tt_end))
            
            daily_stege_code = xr.where(mask_before_init, 0, 
                                        xr.where(mask_init, stage_codes[0],
                                                 xr.where(mask_dev, stage_codes[1],
                                                          xr.where(mask_mid, stage_codes[2],
                                                                   xr.where(mask_late, stage_codes[3], daily_stege_code)
                                                                  )
                                                         )
                                                )
                                       ) 
            
            # daily_stege_code['time'] = daily_tt.time
            
            for i, c in zip(range(0, len(day_stages)), stage_codes):
                mask_stage = ((daily_stege_code == c) & (day_dev == 0))
                day_stages[i] = xr.where(mask_stage, day_of_year, day_stages[i])
            
            daily_stege_code = daily_stege_code.where(~mask_data, other=np.nan)
            # daily_stege_code = daily_stege_code.expand_dims('time')
            # daily_stege_code.name='Phenology_stage_code'
            
            # if daily_tt.time.dt.date == xarray_tt[0].time.dt.date:
            #     stege_code = daily_stege_code
            # else:
            #     stege_code = xr.merge([stege_code, daily_stege_code])
            
            stege_code.loc[daily_tt.time] = daily_stege_code
            
            
        for i, n in zip(range(0, len(day_stages)), stage_names):    
            day_stages[i] = day_stages[i].expand_dims('time')
            day_stages[i] = day_stages[i].where(~mask_data, other=np.nan)
            day_stages[i].name = n 
            # day_stages[i]['time'] = stage_init.time 
            pheno_stage_ds = xr.merge([pheno_stage_ds, day_stages[i]])

    stege_code.to_netcdf(join(out_dir, 'VCSN_{}_{}-{}.nc'.format(daily_stege_code.name, year_list[0], year_list[-1])))
    pheno_stage_ds.to_netcdf(join(out_dir, 'VCSN_growth_stages_date_{}-{}.nc'.format(year_list[0], year_list[-1])))


def k_b_kc_func(grow_stage_ds, pheno_code_ds, out_dir):
    '''
    Calculate the slope (k), and the intercept (b) of kc function
    '''
    pheno_code = pheno_code_ds.psc
    start_date = '09-01'
    end_date = '08-31'
    
    p_code = {'out':0, 'init':1, 'dev':2, 'mid':3, 'late':4}
    
    kc = {'out':0, 'init':0.2, 'dev':1.2, 'mid':1.2, 'end':0.6}
    b_const = {'out':1, 'init':0.2, 'mid':1.2}
    
    year_list = list(set(pheno_code.time.dt.year.values))
    
    data_Maarray = pheno_code[0].to_masked_array()
    mask_data = data_Maarray.mask
    
    kc_array = xr.zeros_like(pheno_code)
    kc_array.name = 'kc'
    kc_array.attrs={'long_name': 'crop coefficient (maize)', 
              'standard_name': 'kc', 
              'units': '', 
              'description': 'crop coefficient (maize) base on maize phenology code', 
              'missing': -9999.0}
    
    for year in year_list:
        
        print('Process {} at {}'.format(year, dt.datetime.now()))

        # sowing_date = grow_stage_ds.sowing_date.sel(time=grow_stage_ds.sowing_date.time.dt.year.isin([year]))[0]
        development_date = grow_stage_ds.development_date.sel(time=grow_stage_ds.development_date.time.dt.year.isin([year]))[0]
        mid_stage_date = grow_stage_ds.mid_stage_date.sel(time=grow_stage_ds.mid_stage_date.time.dt.year.isin([year]))[0]
        late_stage_date = grow_stage_ds.late_stage_date.sel(time=grow_stage_ds.late_stage_date.time.dt.year.isin([year]))[0]
        end_stage_date = grow_stage_ds.end_stage_date.sel(time=grow_stage_ds.end_stage_date.time.dt.year.isin([year]))[0]
        
        if year == year_list[0]:
            yearly_pheno_code = pheno_code.sel(time=slice('{}-{}'.format(year, '01-01'), '{}-{}'.format(year+1, end_date)))
        else:
            yearly_pheno_code = pheno_code.sel(time=slice('{}-{}'.format(year, start_date), '{}-{}'.format(year+1, end_date)))
        
        period_dev2mid = mid_stage_date-development_date
        period_dev2mid_not0 = period_dev2mid.where(period_dev2mid != 0, other=np.nan)
        
        period_late2end = end_stage_date-late_stage_date
        period_late2end_not0 = period_late2end.where(period_late2end != 0, other=np.nan)
        
        
        for daily_pc in yearly_pheno_code: 
            
            daily_k = daily_pc.copy()   # daily_k is the daily slope of kc function
            daily_k = daily_k.where(daily_k==0, other=0)
            
            daily_b = daily_pc.copy()   # daily_b is the daily intercept of kc function
            daily_b = daily_b.where(daily_b==0, other=0)
                       
            mask_out_growth = (daily_pc == p_code['out'])
            mask_init = (daily_pc == p_code['init'])
            mask_dev = (daily_pc == p_code['dev'])
            mask_mid = (daily_pc == p_code['mid'])
            mask_late = (daily_pc == p_code['late'])
            
            daily_k = xr.where((mask_out_growth | mask_init | mask_mid), 0, 
                               xr.where(mask_dev, (kc['dev']-kc['init'])/period_dev2mid_not0,
                                        xr.where(mask_late, (kc['end']-kc['mid'])/period_late2end_not0, daily_k)
                                       )
                              ) 
            
            daily_b = xr.where(mask_out_growth, b_const['out'],
                               xr.where(mask_init, b_const['init'],
                                        xr.where(mask_dev, kc['init'] - daily_k * development_date,
                                                 xr.where(mask_mid, b_const['mid'],
                                                          xr.where(mask_late, kc['mid'] - daily_k * late_stage_date, daily_k)
                                                         )
                                                )
                                       )
                              ) 
            
            day_of_year = int(daily_k.time.dt.dayofyear)
            
            y = int(daily_k.time.dt.year.values)
            if y != year:
                last_date_of_last_year = dt.date.fromisoformat('{}-12-31'.format(year))
                last_day_of_last_year = last_date_of_last_year.timetuple().tm_yday
                day_of_year = last_day_of_last_year + day_of_year
            
            daily_kc = daily_k * day_of_year + daily_b
            
# =============================================================================
#             daily_k['time'] = daily_pc.time
#             daily_k = daily_k.where(~mask_data, other=np.nan)
#             daily_k = daily_k.expand_dims('time')
#             daily_k.name='slope_kc_function'
#             
#             daily_b['time'] = daily_pc.time
#             daily_b = daily_b.where(~mask_data, other=np.nan)
#             daily_b = daily_b.expand_dims('time')
#             daily_b.name='intercept_kc_function'
# =============================================================================
            
            # daily_kc['time'] = daily_pc.time
            daily_kc = daily_kc.where(~mask_data, other=np.nan)
            kc_array.loc[daily_pc.time] = daily_kc
            # daily_kc = daily_kc.expand_dims('time')
            # daily_kc.name='crop_coefficient_maize'
            
            # pheno_code_ds = xr.merge([pheno_code_ds, daily_k, daily_b, daily_kc])
    
    pheno_code_ds = xr.merge([pheno_code_ds, kc_array])
    pheno_code_ds.to_netcdf(join(out_dir, 'VCSN_{}_maize_{}-{}.nc'.format(kc_array.name, year_list[0], year_list[-1])))
    
    
def main():
    
    # cdo = Cdo()
    
    indir = r'./VCSN'
    outdir = r'./VCSN/drought'
    
    f_tmean = join(indir, 'VCSN_Tmean5k_1972-2000.nc')
    f_tt = join(outdir, 'VCSN_TT_1972-2000.nc')
    f_sowing = join(outdir, 'VCSN_sowing_date_1972-2000.nc')
    f_growth_stage_date = join(outdir, 'VCSN_growth_stages_date_1972-2000.nc')
    f_pheno_code = join(outdir, 'VCSN_psc_1972-2000.nc')
    
    
# =============================================================================
#     with xr.open_dataset(f_tmean) as tmean_ds:
#         xarray_tmean = tmean_ds.tmean.sel(time=tmean_ds.tmean.time.dt.month.isin([9,10,11,12]))
#         date_initial_stage(xarray_tmean, outdir)
# =============================================================================
# =============================================================================
#     # step 1 find the initial date (sowing date) based on mean temperature
#     xarray_tmean = cdo.selmonth('9/12', input=f_tmean, returnXArray='tmean')
#     date_initial_stage(xarray_tmean)
# =============================================================================
    
# =============================================================================
#     # step 2 create phenology stage date and code layers
#     
#     with xr.open_dataset(f_tt) as tt_ds, xr.open_dataset(f_sowing) as sowing_ds:
#         code_pheno_stage(sowing_ds, tt_ds, outdir)
# =============================================================================
    
    # step 3 calculate kc
    with xr.open_dataset(f_growth_stage_date) as growth_ds, xr.open_dataset(f_pheno_code) as pheno_code_ds:
        k_b_kc_func(growth_ds, pheno_code_ds, outdir)
    
    
    
if __name__ == '__main__':
    main()