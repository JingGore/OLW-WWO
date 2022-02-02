#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 11:52:06 2021

@author: jing
"""

# from cdo import Cdo
import datetime as dt
import xarray as xr
import numpy as np
import rioxarray as riox
from os.path import join
from cdo import Cdo
import os

def soil_moisture_deficit(pcp_xarray, etc_xarray, taw, raw, year_list, out_f_ped, out_f_smd):
    
    
    
    pcp0 = pcp_xarray[0]
    
    taw = taw.rename({'x': 'longitude','y': 'latitude', 'band': 'time'})
    raw = raw.rename({'x': 'longitude','y': 'latitude', 'band': 'time'})
    temp, taw = xr.align(pcp0, taw, join="override")
    temp, raw = xr.align(pcp0, raw, join="override")
    
    taw = xr.where(((np.isnan(pcp0)) | (taw[0]>10000) | (taw[0]<-9000)), np.nan, taw[0])
    raw = xr.where(((np.isnan(pcp0)) | (raw[0]>10000) | (raw[0]<-9000)), np.nan, raw[0])
    
    swc = taw        #soil water content
    tp = taw - raw   #trigger point
    
    data_Maarray = taw.to_masked_array()
    mask_data = data_Maarray.mask
    
    # aet_list = []
    # ped_list = []
    # smd_list = []
    
    aet_array = xr.zeros_like(pcp_xarray)
    aet_array.name = 'aet'
    aet_array.attrs={'long_name': 'Actural Evapotranspiration', 
              'standard_name': 'aet', 
              'units': 'mm', 
              'description': 'Daily Actural Evapotranspiration', 
              'missing': -9999.0}
    
    ped_array = xr.zeros_like(pcp_xarray)
    ped_array.name = 'ped'
    ped_array.attrs={'long_name': 'Potential Evapotranspiration Deficit', 
              'standard_name': 'ped', 
              'units': 'mm', 
              'description': 'Daily Potential Evapotranspiration Defici. It is the difference between ETc and AET', 
              'missing': -9999.0}
    
    smd_array = xr.zeros_like(pcp_xarray)
    smd_array.name = 'smd'
    smd_array.attrs={'long_name': 'Soil Moisture Deficit', 
              'standard_name': 'smd', 
              'units': 'mm', 
              'description': 'Daily Soil Moisture Deficit. It is the difference between PAW and SWC (soil water content)', 
              'missing': -9999.0}
    
    year = year_list[0] - 1
    for pcp, etc in zip(pcp_xarray, etc_xarray): 
        
        y = int(pcp.time.dt.year.values)
        if y != year:
            print('Process {} at {}'.format(y, dt.datetime.now()))
            year = y
        
        ks = swc/tp
        aet = xr.where(swc>tp, etc, etc*ks)   #actural evapotranspiration
        
        ed = etc - aet  #evapotranspiration deficit
        
        swc = swc + pcp - aet
        
        swc = xr.where(swc>taw, taw, xr.where(swc<0 ,0, swc))
        
        smd = taw - swc

        # aet['time'] = pcp.time
        aet = aet.where(~mask_data, other=np.nan)
        # aet = aet.expand_dims('time')
        # aet.name='AET'
        
        # ed['time'] = pcp.time
        ed = ed.where(~mask_data, other=np.nan)
        # ed = ed.expand_dims('time')
        # ed.name='PED'
        
        # smd['time'] = pcp.time
        smd = smd.where(~mask_data, other=np.nan)
        # smd = smd.expand_dims('time')
        # smd.name='SMD'
        
# =============================================================================
#         if pcp.time.dt.date == pcp0.time.dt.date:
#             et_ds = aet
#             et_ds = xr.merge([et_ds, ed])
#             
#             # aet_ds = aet
#             # ped_ds = ed
#             
#             smd_ds = smd
#         else:
#             # aet_ds = xr.concat([aet_ds, aet], dim="time")
#             # ped_ds = xr.concat([ped_ds, ed], dim="time")
#             # smd_ds = xr.concat([smd_ds, smd], dim="time")
#             et_ds = xr.merge([et_ds, aet, ed])
#             smd_ds = xr.merge([smd_ds, smd])
# =============================================================================
        
        # aet_list.append(aet)
        # ped_list.append(ed)
        # smd_list.append(smd)
        
        aet_array.loc[pcp.time] = aet
        ped_array.loc[pcp.time] = ed
        smd_array.loc[pcp.time] = smd
    
    # aet_ds = xr.concat(aet_list, dim="time")
    # ped_ds = xr.concat(ped_list, dim="time")
    # smd_ds = xr.concat(smd_list, dim="time")
    
    taw['time'] = pcp0.time
    taw = taw.expand_dims('time')
    taw.name='TAW'
    
    raw['time'] = pcp0.time 
    raw = raw.expand_dims('time')
    raw.name='RAW'
    
    smd_ds = xr.merge([smd_array, taw, raw])
    
    ped_ds = xr.merge([ped_array, aet_array])
    
    ped_ds.to_netcdf(out_f_ped)
    smd_ds.to_netcdf(out_f_smd)
        
        
def main():
    
    indir = r'./VCSN/drought'
    
    f_pcp = join(indir, r'VCSN_Rainbc5k_1972-2000.nc')
    f_etc = join(indir, r'VCSN_ETc_maize_1972-2000.nc')
    f_taw = 'PAW.tif'
    f_raw = 'RAW.tif'
    
    with xr.open_dataset(f_pcp) as pcp_ds, xr.open_dataset(f_etc) as etc_ds:
        
        # pcp_ds = xr.open_dataset(f_pcp)
        # etc_ds = xr.open_dataset(f_etc)
        taw = riox.open_rasterio(f_taw)
        raw = riox.open_rasterio(f_raw)
        
        
        year_list = list(set(pcp_ds.time.dt.year.values))
        
        num_group = len(year_list)//10 + 1
        
        group_year_list = [year_list[10*x:10*(x+1)] for x in range(num_group)]
        
        # pcp_xarray = pcp_ds.rain
        # etc_xarray = etc_ds.etc
        
        ped_f_list = ''
        smd_f_list = ''
        
        for y in group_year_list:
            
            pcp_xarray = pcp_ds.rain.sel(time=pcp_ds.rain.time.dt.year.isin(y))
            etc_xarray = etc_ds.etc.sel(time=etc_ds.etc.time.dt.year.isin(y))
            
            out_f_ped = join(indir, r'VCSN_PED_maize_{}-{}.nc'.format(y[0], y[-1]))
            out_f_smd = join(indir, r'VCSN_SMD_maize_{}-{}.nc'.format(y[0], y[-1]))
            
            soil_moisture_deficit(pcp_xarray, etc_xarray, taw, raw, y, out_f_ped, out_f_smd)
    
            ped_f_list = ped_f_list + ' ' + out_f_ped
            smd_f_list = smd_f_list + ' ' + out_f_smd
    
    # del pcp_ds
    # del etc_ds
    
    # print(ped_f_list)
    
    out_merged_ped_f = join(indir, r'VCSN_PED_maize_{}-{}.nc'.format(year_list[0], year_list[-1]))
    out_merged_smd_f = join(indir, r'VCSN_SMD_maize_{}-{}.nc'.format(year_list[0], year_list[-1]))
    
    cdo = Cdo()
    cdo.mergetime(input=ped_f_list, output=out_merged_ped_f, options='-O')
    cdo.mergetime(input=smd_f_list, output=out_merged_smd_f, options='-O')

    for ped_f, smd_f in zip(ped_f_list[1:].split(' '), smd_f_list[1:].split(' ')):
        os.remove(ped_f)
        os.remove(smd_f)


if __name__ == '__main__':
    main()        
        
        
        
        
        
        
        
        
        
        