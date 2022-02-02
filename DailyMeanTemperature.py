# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 16:00:24 2021

@author: GuoJ
"""
import numpy as np

def get_period_fraction(period):
        
    # fraction of 3 hours temperature for the calculation of daily mean temperature based on daily Tmin and Tmax
    p_frac = 0.931 + 0.114 * period - 0.0703 * period**2 + 0.0053 * period**3
    
    return p_frac

def daily_mean_3hrs(tmin, tmax, period_output=False):
    
    periods = 8
    
    if len(np.shape(tmin)) == 0: # The input is integer
        ttotal = 0
    else:
        ttotal = np.zeros(tmin.shape)  # The input is arrry
    
    
    tmean_3hrs_list = []
    
    for p in range(1, periods+1):
        
        pf = get_period_fraction(p)
        
        tmean_3hrs = pf*(tmax-tmin) + tmin
        
        tmean_3hrs_list.append(tmean_3hrs)
        
        ttotal += tmean_3hrs
        
    tmean_day = ttotal / periods
    
    if period_output == True:
        return tmean_day, tmean_3hrs_list
    else:
        return tmean_day

def daily_mean_avg(tmin, tmax):
    
    return (tmin + tmax) / 2


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

        

        # tt_total = np.zeros(periodic_temperature_array_list[0].shape)
        
        # for array in periodic_temperature_array_list:
            
        #     tt_3hrs_array = np.zeros(array.shape)
            
        
        if len(np.shape(temperature)) == 0: # The input is integer
            
            if temperature < argv[0][0]:
                tt = 0
                
            elif temperature >= argv[0][0] and temperature <= argv[1][0]: 
                tt = __EquationStraightLine__(temperature, argv[0][0], argv[0][1], argv[1][0], argv[1][1])
                
            elif temperature >= argv[1][0] and temperature <= argv[2][0]:
                tt = __EquationStraightLine__(temperature, argv[1][0], argv[1][1], argv[2][0], argv[2][1])
                
            elif temperature >= argv[2][0] and temperature <= argv[3][0]:
                tt = __EquationStraightLine__(temperature, argv[2][0], argv[2][1], argv[3][0], argv[3][1])
                
            elif temperature > argv[3][0]:
                tt = 0
                
        else:                               # The input is arrry

            tt = np.zeros(temperature.shape)
                
            exp_1 = temperature < argv[0][0]
            exp_2 = np.logical_and(temperature >= argv[0][0], temperature <= argv[1][0])
            exp_3 = np.logical_and(temperature >= argv[1][0], temperature <= argv[2][0])
            exp_4 = np.logical_and(temperature >= argv[2][0], temperature <= argv[3][0])
            exp_5 = temperature > argv[3][0]
                
            tt = np.where(exp_1, 0, 
                        np.where(exp_2, __EquationStraightLine__(temperature, argv[0][0], argv[0][1], argv[1][0], argv[1][1]),
                                np.where(exp_3, __EquationStraightLine__(temperature, argv[1][0], argv[1][1], argv[2][0], argv[2][1]),
                                         np.where(exp_4, __EquationStraightLine__(temperature, argv[2][0], argv[2][1], argv[3][0], argv[3][1]),
                                                  np.where(exp_5, 0, tt)))))
            
            # tt_total += tt_3hrs_array
            
            
         
        # tt_daily_mean = tt_total / len(periodic_temperature_array_list)    
            
        return tt
    
    
    else:
                
        print('Error! Incorrect number of argument.')
        return None
    
