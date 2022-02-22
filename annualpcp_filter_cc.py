# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 09:21:16 2020

Average annual rainfall

@author: GuoJ
"""
import numpy as np
from os.path import join
from raster import Raster

def main():
    
    workspace = r'C:\Users\guoj\projects\PRJ479 Regen Hillcountry\sub_clover'
    pcp_dir = join(workspace, r'data\cmip5\annual')
    
    processed_dir = join(workspace, r'processed')
    
    rst = Raster()
    year_dict = {'mid': '2031-2060', 'end': '2071-2100'}
    
    rcps = ['RCP8.5']
    gcms = ['BCC-CSM1.1']
    
    pheno = ['emergence', 'v3', 'v6', 'r3_early', 'r3_late', 'r6_early', 'r6_late','r11_early', 'r11_late']
    
    pcp_limit = 1200
    
    for rcp in rcps:
        for gcm in gcms:
            for y in year_dict:
                pcp_file = join(pcp_dir, rcp, gcm, r'annual_mean_pcp_{}.tif'.format(year_dict[y]))
                
                pcp_array = rst.getRasterArray(pcp_file)
                no_data = rst.getNoDataValue(pcp_file)
                
                for p in pheno:
                    
                    d_dir = join(processed_dir, rcp, y, p, gcm)
                
                    for fname in ['mean', 'median', 'std', 'min', 'max', 'mean_plus_std', 'mean_minus_std']:
                
                       rst_array = rst.getRasterArray(join(d_dir, '{}_date_{}({}).tif'.format(p, fname, year_dict[y])))
                       rst_array = np.where(pcp_array>=pcp_limit, no_data, rst_array)
                       
                       rst.array2Raster(rst_array, pcp_file, join(d_dir, '{}_date_{}({})_pcp_less_{}.tif'.format(p, fname, year_dict[y], pcp_limit)))
        

if __name__ == '__main__':
    main() 