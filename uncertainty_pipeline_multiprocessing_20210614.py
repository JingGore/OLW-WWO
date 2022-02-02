####Spatial version of WaterYield model###

#####Version 0.0.1#####

###Created  by VVV on 11.09.2015###

###Based on Jonh Dymond code####

###There are some hardcoded paths to files

####Rasters with TAW,RAW, CropCoef, intfrac,intstore should be in one directory

###Rainfalls and PET rasters are in separated directories


import numpy as np
#import numpy.ma as ma
import datetime as dt
import os.path
from os.path import join
#import sys

from os import walk
import random
#sys.path.append('../scripts/')
#sys.path.append('../Plot/')

#import CSVOperation
import ResultAnalysis as ra
import RasterProcessing as RP
import pandas as pd

from raster import Raster
from config import ConfigParameters
import gc
from collections import OrderedDict

#import histogram as Hist
#import numba as nb


def transfer(x):
        return {
            'pa': 'PAST',   '26':'RCP2.6', 
            '45': 'RCP4.5', '60': 'RCP6.0',
            '85': 'RCP8.5',
            'BCC':'BCC-CSM1.1', 'CES':'CESM1-CAM5',
            'GFD':'GFDL-CM3',   'GIS':'GISS-EL-R',
            'HAD':'HadGEM2-ES', 'NOR':'NorESM1-M',
        }.get(x, 'NA')

def datesCal(startDate,endDate):   #function to calculate days between two dates
    dt1=dt.datetime(int(startDate.split('-')[0]), int(startDate.split('-')[1]),int(startDate.split('-')[2]))
    end=dt.datetime(int(endDate.split('-')[0]), int(endDate.split('-')[1]),int(endDate.split('-')[2]))
    step = dt.timedelta(days=1)
    result = []

    while dt1 <= end:
        if dt1.day!=29 or dt1.month!=2:
            result.append(dt1.strftime('%m%d%Y'))
        dt1 += step
    return result

class extractTimeInfo():
    
    def __init__(self, path):
        self.path = path
        self.file = []

        for (dirpath, dirnames, filenames) in walk(path):
            self.file.extend(filenames)
            break
        
    def extractYears(self):
        years = sorted(list(set([y.split('-')[0] for y in self.file])))
        return years
    
    def extractDates(self, year): 
#        dates = sorted(list(d.split('.')[0] for d in self.file if d.split('-')[0] == year and d.split('.')[-1] == 'tif'))
        fstHalfYear = ['07','08','09','10','11','12']
        scdHalfYear = ['01','02','03','04','05','06']
        dates1 = sorted(list(d.split('.')[0] for d in self.file if d.split('-')[0] == str(year) and d.split('.')[-1] == 'tif' and d.split('-')[1] in fstHalfYear))
#        dates2 = sorted(list(d.split('.')[0] for d in self.file if d.split('-')[0] == year and d.split('.')[-1] == 'tiff' and d.split('-')[1] in scdHalfYear))
        dates2 = sorted(list(d.split('.')[0] for d in self.file if d.split('-')[0] == str(int(year)+1) and d.split('.')[-1] == 'tif' and d.split('-')[1] in scdHalfYear))
        if len(dates2) > 0:
            dates = dates1 + dates2
            return dates
        else:
            return []
        
def extractpixellocation(raster, lcraster):
    '''
    Randomly select X number of pixels, X equal to the unipue values of the raster.
    smuc is the list of unique values of smap soil type.
    idxs is the dictionary of coordinates of each selected pixel, including smap soil type, coordinates and land cover type.
    '''
    maskArray = RP.ReadRaster(lcraster, 1, np.int)
    maskoutvalue = [6,16,21]
    
    inArray = RP.ReadRaster(raster, 1, np.int)
    shp = inArray.shape
    smuc = np.unique(inArray).tolist()
    unexpectedvalues = []
    for s in smuc:
        if s < 100:
            unexpectedvalues.append(s)
            smuc.remove(s)
            
    idxs = {}
    selectedsmuc = []
#    selectedluc = {}
    
    while len(selectedsmuc) != len(smuc):
        idx = [random.choice(range(0,shp[0])), random.choice(range(0,shp[1]))]
        if (inArray[idx[0]][idx[1]] not in unexpectedvalues) and (maskArray[idx[0]][idx[1]] not in maskoutvalue):
            if (inArray[idx[0]][idx[1]] in smuc) and (inArray[idx[0]][idx[1]] not in selectedsmuc):
                idx.append(maskArray[idx[0]][idx[1]])
                idxs[inArray[idx[0]][idx[1]]] = idx
                selectedsmuc.append(inArray[idx[0]][idx[1]])
#                selectedluc[maskArray[idx[0]][idx[1]]] = idx

    return selectedsmuc, idxs

#def export2csv(file, head, data):
#    csvw = CSVOperation.CSVWriting()
#    csvw.WriteLines(file, head, data)

def visulization(csvpixel, csvtaw, csvraw, featurename, outpath):
     
    df = pd.read_csv(csvpixel, thousands=',')
    taw = pd.read_csv(csvtaw, thousands=',')
    raw = pd.read_csv(csvraw, thousands=',')
    
    years = sorted(list(set(df.Year)))

    for y in years:
        print(y)
        #Annual distribution
        print("Visualize water demand distribution")
        df1 = df[(df['Year']==y) & (df['Date']=='{}-06-30'.format(y+1))]
        ra.pixelannualwaterdemanddistribution(df1, y, featurename, outpath)
        
        #Daily accumulated
        print("Visualize accumulated water demand")
        df2 = df[(df['Year']==y)]
        ra.accumulateddailywaterdemand(df2, y, featurename, outpath)
        
        # accumulated water demand against taw and raw
        print("Visualize accumulated water demand against taw and raw")
        ytaw = taw[taw['Year']==y]
        yraw = raw[raw['Year']==y]     
        ra.accumulatedwaterdemandagainsttawraw(df1, ytaw, yraw, y, featurename, outpath)
    

def generate_crop_parameters(property_dict, ref_raster):
    rst = Raster()
    rst_array = rst.getRasterArray(ref_raster)
    no_data = rst.getNoDataValue(ref_raster)
    cp_array_dict = OrderedDict()
    
    for p_key in property_dict:
        
        cp_array_dict[p_key] = np.where(rst_array==no_data, no_data, property_dict[p_key])

    return cp_array_dict    



# =============================================================================
# spec = [
#     ('no_data', nb.float32),               # a simple scalar field
#     ('taw', nb.float32[:,:]),          # an array field
#     ('raw', nb.float32[:,:]), 
#     ('tp', nb.float32[:,:]), 
#     ('kc', nb.float32), 
#     ('icf', nb.float32), 
#     ('isc', nb.float32), 
#     ('zeros', nb.float32[:,:]), 
#     ('is_irri', nb.int32), 
# ]
# 
# 
# @nb.jitclass(spec)
# =============================================================================
class WaterBalance(object):
    
    def __init__(self, taw_rst, raw_rst, crop_params_dict, is_irri=False, scenario_id_list=[0,0,0]):  
           
        
        # read the soil and crop properties
        self.rst = Raster()
        
        self.ref_rst = taw_rst
        
        self.no_data = self.rst.getNoDataValue(taw_rst)
        self.taw = self.rst.getRasterArray(taw_rst)
        self.raw = self.rst.getRasterArray(raw_rst)
        self.tp = self.taw - self.raw  # Trigger point
        
        cp_params_array_dict = generate_crop_parameters(crop_params_dict, self.ref_rst)
            
        self.kc = cp_params_array_dict['kc']
        self.icf = cp_params_array_dict['icf']
        self.isc = cp_params_array_dict['isc']
        
        self.taw = np.where(self.taw < 0, self.no_data, self.taw)
        self.raw = np.where(self.raw < 0, self.no_data, self.raw)
        self.kc = np.where(self.kc < 0, self.no_data, self.kc)
        self.icf = np.where(self.icf < 0, self.no_data, self.icf)
        self.isc = np.where(self.isc < 0, self.no_data, self.isc)
        
        #only areas of all input soil property rasters have value will be used in modelling 
        exp = np.logical_and(
                             np.logical_and(
                                            np.logical_and(
                                                           np.logical_and(
                                                                          self.taw != self.no_data, self.isc != self.no_data), 
                                                           self.raw != self.no_data), 
                                            self.kc != self.no_data), 
                             self.icf != self.no_data)
        
        self.ref = np.where(exp, self.taw, self.no_data)
        
        self.taw = np.where(self.ref == self.no_data, self.no_data, self.taw)
        self.raw = np.where(self.ref == self.no_data, self.no_data, self.raw)
        self.raw = np.where(self.raw==0, self.taw*0.4, self.raw)
        self.kc = np.where(self.ref == self.no_data, self.no_data, self.kc)
        self.icf = np.where(self.ref == self.no_data, self.no_data, self.icf)
        self.isc = np.where(self.ref == self.no_data, self.no_data, self.isc)
        self.tp[np.where(self.ref == self.no_data)] = self.no_data
        
        self.zeros=np.zeros(self.ref.shape)
        self.zeros[np.where(self.ref == self.no_data)] = self.no_data
     
        self.is_irri = is_irri
        self.sna_lst = scenario_id_list 
        self.year = 0


    # Define a function to return interception in mm/day......
    def interception(self,rain,intstore,intfrac):
    
        x = intfrac*rain
        x=np.where(x>intstore,intstore,x)
        return x
        
    def daily_water_balance(self, PCP, PET, SWC, TAW, RAW, Kc, ICF, ISC, WRC):
        # the maximum total amount of water that can be intercepted at the day is depends on how much water left on canopy of previous day
        # Interception storage capacity on a perticular day (less water will be intercepted if there is some water remain on canopy)

        ISC = ISC - WRC 
        
        # Newly intercepted water
        NIC = np.where(PCP>0, self.interception(PCP, ISC, ICF), 0) 
        
        # Tatol amount of water on canopy of that day
        IC = NIC + WRC 
    
        # Soil water content of next day (1st time)
        SWC_next_day = SWC + PCP - NIC
        
        #Calculate drainage
        DRG = np.where(SWC_next_day>TAW, SWC_next_day-TAW, 0)
       
        # Soil water content of next day (2nd time)
        SWC_next_day = SWC_next_day - DRG
        
        # Reference crop (grass) evapotranspiration which is actually Potential Evapotranspiration, but using different terminology 
        # The term of ET0 is from FAO, and so is the ETc
        ET0 = PET
        
        # Evaporation from canopy
        EC = np.where(ET0<IC, ET0, IC)
        
        # Water remain on canopy to the next day
        WRC_next_day = IC - EC
        
    
        # Crop evapotranspiration under standard condition (Referenced PET of specific crop)
        ETc = (ET0 - EC) * Kc
        
        # Actual evapotranspiration from crop and soil
        TP = TAW - RAW
        TP[np.where(TAW == self.no_data)] = self.no_data

        Ks = SWC_next_day/TP
        AEc = np.where(SWC_next_day > TP, ETc, ETc*Ks)
        
        # Actual total evapotranspiration
        AET = AEc + EC
        
        # Evapotranspiration deficit (the insufficient part of amount of water that crop would need to achieve full production)
        ED = ETc - AEc
        
        # Soil water content of next day (3Rd time)
        SWC_next_day = np.where(SWC_next_day>AEc, SWC_next_day-AEc, 0)
        
        # Soil moisture deficit
        SMD = TAW - SWC_next_day
        
        return SWC_next_day, IC, DRG, EC, ETc, AEc, AET, ED, SMD, WRC_next_day  
        
    
    
    def CalWatBalance(self, year_list, pcpPath, petPath, outPath, annual_output_names, daily_spatial=True, annual_spatial=True): # y is the iteration number of years used to determine the initial water content of the first day of first year 
        
        
        if (daily_spatial is True) or (annual_spatial is True):
            if not os.path.exists(outPath):
                os.makedirs(outPath, exist_ok = True)
        
        pcptimeExtractor = extractTimeInfo(pcpPath)
        
        last_year = 0
        for year in year_list:
            self.year = int(year)
            # daily climate files in each year
            dates = pcptimeExtractor.extractDates(year)
            # print('Processing year {} at {}.'.format(year, dt.datetime.now()))

            pcp_total = np.zeros(self.ref.shape)
            pet_total = np.zeros(self.ref.shape)
            ped_total = np.zeros(self.ref.shape)
            irr_total = np.zeros(self.ref.shape)
            dt_days_total = np.zeros(self.ref.shape)
        
            # TAW value at the beginning of the model run or the current year is not consecutive to the last year
            if year==year_list[0]:       
                might_reset_swc = True
            elif year-last_year != 1: 
                might_reset_swc = True
            else:
                might_reset_swc = False
            
            last_year = year
            
            for date in dates:
#                print('{} {}'.format(date, dt.datetime.now()))
    #            if '{}'.format(date) == '2004-09-26':
    #                print(date)
                    
                pcpfile = join(pcpPath, "{}.tif".format(date))
                petfile = join(petPath, "{}.tif".format(date))
                pcp = self.rst.getRasterArray(pcpfile)
                pet = self.rst.getRasterArray(petfile)
                #Calculate reference PET (pet_ref) for specific crops based on crop coefficient,
                #the pet_ref will be used as 'PET' for that crop
    
                #### model
                if date == dates[0] and might_reset_swc is True:      # first day
                    swc = np.copy(self.taw)    # intial soil water content
                    wrc = 0               # first day water on canopy
                    swc_next_day, ic, drg, ec, etc, aec, aet, ped, smd, wrc_next_day = self.daily_water_balance(pcp, pet, swc, self.taw, self.raw, self.kc, self.icf, self.isc, wrc)
                else:
                    swc = swc_next_day
                    wrc = wrc_next_day
                    swc_next_day, ic, drg, ec, etc, aec, aet, ped, smd, wrc_next_day = self.daily_water_balance(pcp, pet, swc, self.taw, self.raw, self.kc, self.icf, self.isc, wrc)
                
                
                if self.is_irri is False:
                    # drought days (when swc<raw)
                    dt_days = np.where(swc<self.tp, 1, 0)
                    dt_days_total = dt_days_total + dt_days
                    irr = np.zeros(self.taw.shape)
                else:
                    # Irrigation 
                    irr = np.where(swc<self.tp, 0.9*self.taw-swc, 0)
                    swc_next_day = swc_next_day + irr
                    irr_total = irr_total + irr
                    # here can count the annual irrigation days
                
                # create daily spatial outputs
                if daily_spatial == True:
                    
                    out_dir = join(outPath, 'daily')
                    if not os.path.exists(out_dir):
                        os.makedirs(out_dir, exist_ok = True)
                        
                    out_dir = out_dir.replace('\\','\\\\')
                    
                    out_file = join(out_dir, "PED_{}_s{}_c{}_m{}.tif".format(date,self.sna_lst[0],self.sna_lst[1],self.sna_lst[2]))
                    ped[np.where(self.taw == self.no_data)] = self.no_data
                    self.rst.array2Raster(ped, self.ref_rst, out_file)

                    if self.is_irri is False:
                        out_file = join(out_dir, "DTD_{}_s{}_c{}_m{}.tif".format(date,self.sna_lst[0],self.sna_lst[1],self.sna_lst[2]))
                        dt_days[np.where(self.taw == self.no_data)] = self.no_data
                        self.rst.array2Raster(dt_days, self.ref_rst, out_file)
                    
                    else:
                        out_file = join(out_dir, "IRR{}_s{}_c{}_m{}.tif".format(date,self.sna_lst[0],self.sna_lst[1],self.sna_lst[2]))
                        irr[np.where(self.taw == self.no_data)] = self.no_data
                        self.rst.array2Raster(irr, self.ref_rst, out_file)
                            
                #Annual water balance info (pixel based)
                pcp_total = pcp_total + pcp
                pet_total = pet_total + pet            
                ped_total = ped_total + ped
            
        
            # ***** Model of iteration of one year ends here *****
           
            # Export annual spatial outputs    
            if annual_spatial == True:
                
                out_dir = join(outPath, 'annual')
                if not os.path.exists(out_dir):
                    os.makedirs(out_dir, exist_ok = True)
                
                out_dir = out_dir.replace('\\','\\\\')
                
                pcp_total[np.where(self.taw == self.no_data)] = self.no_data
                RainFile = join(out_dir, "{}_{}_s{}_c{}_m{}.tif".format(year, annual_output_names['pcp'],self.sna_lst[0],self.sna_lst[1],self.sna_lst[2]))
                self.rst.array2Raster(pcp_total, self.ref_rst, RainFile)
                
                pet_total[np.where(self.taw == self.no_data)] = self.no_data
                PETFile = join(out_dir, "{}_{}_s{}_c{}_m{}.tif".format(year, annual_output_names['pet'],self.sna_lst[0],self.sna_lst[1],self.sna_lst[2]))
                self.rst.array2Raster(pet_total, self.ref_rst, PETFile)
                
                ped_total[np.where(self.taw == self.no_data)] = self.no_data
                ped_file = join(out_dir, "{}_{}_s{}_c{}_m{}.tif".format(year, 'PED',self.sna_lst[0],self.sna_lst[1],self.sna_lst[2]))
                self.rst.array2Raster(ped_total, self.ref_rst, ped_file)
                
                if self.is_irri is False:
                    dt_days_total[np.where(self.taw == self.no_data)] = self.no_data
                    dt_days_file = join(out_dir, "{}_{}_s{}_c{}_m{}.tif".format(year, 'DTD',self.sna_lst[0],self.sna_lst[1],self.sna_lst[2]))
                    self.rst.array2Raster(dt_days_total, self.ref_rst, dt_days_file)
                else:
                    irr_total[np.where(self.taw == self.no_data)] = self.no_data
                    irr_file = join(out_dir, "{}_{}_s{}_c{}_m{}.tif".format(year, 'IRR',self.sna_lst[0],self.sna_lst[1],self.sna_lst[2]))
                    self.rst.array2Raster(irr_total, self.ref_rst, irr_file)
    
        # set the initial water content of next year
#        return swc_next_day
        
    
def uncertainty_pipeline(soil_realization):
    
    conf = ConfigParameters(r'config.ini', r'projectConfig')

    out_dir = conf.GetOutputDir()
    
    year_list = conf.GetWorkingYearList()
    
    crop_csv, climate_csv = conf.GetCropClimateCSV()
    
    df_crop = pd.read_csv(crop_csv, thousands=',')
    df_climate = pd.read_csv(climate_csv, thousands=',')
    
    # climate info
    clim_dir = conf.GetClimateCovariateDir()
    pcp_folder, pet_folder = conf.GetClimateParams()
    
    #soil info
    soil_dir = conf.GetSiolCovariateDir()
    taw, raw = conf.GetSiolParams()
    # lc_rst = conf.GetLandcoverRaster()
    
    is_irri = conf.GetIrrigation()

    annual_names = conf.GetOutputParams()
    
    crop_params_dict = OrderedDict([('kc',0),('icf',0),('isc',0)])
    
    if is_irri is True:
        out_dir = '{}_irrigated'.format(out_dir)
    else:
        out_dir = '{}_not_irrigated'.format(out_dir)

        
    out_sub_dir = join(out_dir, 'soil_realization_{}'.format(soil_realization))
    # if not os.path.exists(out_dir):
    #     os.makedirs(out_dir, exist_ok = True)
    
    taw_rst = join(soil_dir, '{}_{}.tif'.format(taw, soil_realization))
    raw_rst = join(soil_dir, '{}_{}.tif'.format(raw, soil_realization))
    
    for crop_index, crop_row in df_crop.iterrows():
        
        crop_params_dict.update([
                                 ('kc', crop_row['kc']),
                                 ('icf',crop_row['icf']),
                                 ('isc',crop_row['isc'])
                                ])

        for climate_index, climate_row in df_climate.iterrows():
            
            if str(climate_row['id'])[0] == '1':  #PAST scenario only
                rcp = climate_row['rcp']
                gcm = climate_row['gcm']
                
                pcpPath = join(clim_dir, rcp, gcm, pcp_folder)
                petPath = join(clim_dir, rcp, gcm, pet_folder)
                
                
                # if either of PCP or PET folder does not exist than do not excute that scenario or climate model
                if os.path.exists(pcpPath) and os.path.exists(petPath):
                    scenario_id_list = [str(soil_realization), str(int(crop_row['sample'])), str(climate_row['id'])]
                    # print('Model start running soil realization {}, crop params {} and climate {} at {}.'.format(soil_realization, crop_row['sample'],climate_row['id'], dt.datetime.now()))
                    wb = WaterBalance(taw_rst, raw_rst, crop_params_dict, is_irri, scenario_id_list)                      
                    wb.CalWatBalance(year_list, pcpPath, petPath, out_sub_dir, annual_names, daily_spatial=False)

                        
                        # print('Water balance info of year {} has been created at {}.'.format(y, dt.datetime.now()))

                        
                else:
                    # pop up some msgs to show difference exists between two climate variables
                    pass
                        
    gc.collect()
         
        
        
        
        
        
        
        
    