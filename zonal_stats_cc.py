# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 13:48:22 2021

@author: GuoJ
"""

from rasterstats import zonal_stats

import pandas as pd

from os.path import join
from os import walk

import datetime as dt


def get_date_from_days(year, days):
    
    days = str(days)
    days.rjust(3 + len(days), '0')
    date = dt.datetime.strptime(str(year) + "-" + days, "%Y-%j").strftime("%m-%d")
    
    return str(date)

def main():
    
    workspace = r'C:\Users\guoj\projects\PRJ479 Regen Hillcountry\sub_clover\processed\RCP8.5\end\summary\BCC-CSM1.1'
    
    zone_shp = r'C:\Users\guoj\projects\common_use\nz-district-2012_wgs84.shp'
    
    year = 1973
    
    for (subdirpath, subdirname, filenames) in walk(workspace):
        for f in filenames:
            if f.split('.')[-1].lower()[:3] == 'tif':
            # if f.split('.')[0].lower()[:4] == 'safe' and f.split('.')[-1].lower()[:3] == 'tif':     
                print(f)
                
                fn = join(subdirpath, f)
                
                zonal_list = zonal_stats(zone_shp, fn, stats="mean", geojson_out=True, included_attributes=['Name_short']) #, copy_properties=True)  
                
                for i, item in enumerate(zonal_list):
                    
                    v = item['properties']
                    
                    if v['mean'] is not None:
                        days = int(v['mean']) + 1
                        if days > 365:
                            days = days - 365
                            
                        date = get_date_from_days(year, days)
                        # date=days
                    else: 
                        date = None
                    
                    if i == 0:
                        df = pd.DataFrame(v, index=[v['dst_id']])
                        df['date'] = '"{}"'.format(date)
                    else:
                        df_1 = pd.DataFrame(v, index=[v['dst_id']])
                        df_1['date'] = '"{}"'.format(date)
                        
                        df = df.append(df_1)
                        
                df.to_csv(join(workspace, '{}.csv'.format(f.split('.')[0])))
        
if __name__ == '__main__':
    
    main()