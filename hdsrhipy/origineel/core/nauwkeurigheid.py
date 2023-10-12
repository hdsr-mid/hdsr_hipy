# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 09:17:22 2021

@author: HKV lijn in water 2021
Contact: Ruud Hurkmans (hurkmans@hkv.nl)
"""

import os
import sys
import time
import numpy as np
import xarray as xr
import pandas as pd
import rioxarray as rio
import imod
import shutil
from tqdm.auto import tqdm
import logging
logger = logging.getLogger(__name__)
import random

from hdsrhipy import Groundwater
from hdsrhipy import Maatgevend
from hdsrhipy import Runoff

def bootstrap(data, nyear=None, n=3):
    """Functie om willekeurig jaren te samplen volgens het bootstrap-principe"""
    if nyear is None:
        nyear = 8
    if hasattr(data,'time'):
        time = [pd.Timestamp(data['time'].values[i]).year for i in range(len(data['time']))]            
    else:
        time  = [pd.Timestamp(data.iloc[i,0]).year for i in range(len(data.iloc[:,0]))]  
    uyears = list(set(time))   
    samples = []
    for i in range(n):
        datlist = []
        for j in range(nyear):
            yr = random.choice(uyears)
            outinds = [i for i,t in enumerate(time) if t == yr]
            if len(data.shape)==3:
                datlist.append(data[outinds,:,:])
            elif len(data.shape)==2:
                if j==0:
                    outdat = data.iloc[outinds,:]
                else:                    
                    outdat = outdat.append(data.iloc[outinds,:])
            else: 
                raise TypeError('A 2D or 3D array is required.')         
        if hasattr(data,'time'):
            outdat = xr.concat(datlist,'time')
        #else:
        #    outdat.iloc[:,0] = data.iloc[:,0]
        samples.append(outdat)
    return samples 

def sample_nauwkeurigheid(data, bandbreedte, n=3):   
    """Voeg noise toe aan een dataset uit een gegeven bandbreedte"""     
    samples = []
    for i in tqdm(range(n)):
        if len(data.shape)==3:
            noise = bandbreedte[0] + np.random.rand(data.shape[0], data.shape[1], data.shape[2])*(bandbreedte[1]-bandbreedte[0])    
            samples.append(data+noise)
        elif len(data.shape)==2:
            temp = data.iloc[:,1:]
            noise = pd.DataFrame(bandbreedte[0] + np.random.rand(temp.shape[0], temp.shape[1])*(bandbreedte[1]-bandbreedte[0]), columns=temp.columns, index=temp.index)
            sumt = temp.add(noise)
            temp2 = data.copy(deep=True)
            temp2.iloc[:,1:] = sumt
            samples.append(temp2)
        else: 
            raise TypeError('A 2D or 3D array is required.')        
    return samples

def MonteCarlo(variable, samples, bootstrap_n=3, n=3, model_path=None, model_name=None):
    """ Bepaal producten op basis van een aantal samples op basis van nauwkeurigheid en bootstrap. De manier waarop hangt af van het dataformat van de data (raster vs shape)"""
    def makeds(arr, ref):
        ds = ref.copy(deep=True)
        ds.values = arr
        return ds
    
    gw = Groundwater(model_path=model_path, name=model_name)
    mg = Maatgevend()
    r = Runoff()
    
    reslist = []
    for i in tqdm(range(n)):
        dset = samples[np.random.randint(len(samples))]     
        dset2 = bootstrap(dset, n=bootstrap_n)
        for j in dset2:
            if variable == 'seepage':
                reslist.append(gw.seepage_season_means(dataset=j))
                refda = gw.seepage_season_means(dataset=j)[0]           
            elif variable =='gxg':
                j['time'] = pd.date_range(start='2010-01-01', periods=len(j['time']))                
                gw_stats = gw.get_gxg_raster(j)
                gw_stats['gt'] = gw.get_gt_raster(gw_stats)
                reslist.append(gw_stats)
                refda = gw.seepage_season_means(dataset=j)[0]           
            elif variable == 'normative':
                reslist.append(mg.get_q_norm(j))
            elif variable == 'metaswap_mean':
                reslist.append(r.get_season_stat(dataset=j,stat='mean'))
                refda = r.get_season_stat(dataset=j,stat='mean')[0]       
            elif variable == 'metaswap_min':
                reslist.append(r.get_season_stat(dataset=j,stat='min'))
                refda = r.get_season_stat(dataset=j,stat='min')[0]       
            
           
    mins = []
    maxs = []
    if variable =='seepage':
        for s in range(6):
            print(f'Get min/max for season {s+1}')
            # make a list of rasters for this season
            sublist = np.dstack([sub[s] for sub in reslist])
            mins.append(makeds(np.amin(sublist,axis=2), refda)) 
            maxs.append(makeds(np.amax(sublist,axis=2), refda))                      
            
    elif variable=='gxg':    
        sublistghg = np.dstack([sub['ghg'] for sub in reslist])            
        mins.append(makeds(np.amin(sublistghg,axis=2), refda))
        maxs.append(makeds(np.amax(sublistghg,axis=2), refda))            
        
        sublistglg = np.dstack([sub['glg'] for sub in reslist])
        mins.append(makeds(np.amin(sublistglg,axis=2), refda))
        maxs.append(makeds(np.amax(sublistglg,axis=2), refda))    
        
        sublistgvg = np.dstack([sub['gvg'] for sub in reslist])        
        mins.append(makeds(np.amin(sublistgvg,axis=2), refda))
        maxs.append(makeds(np.amax(sublistgvg,axis=2), refda))
        
        sublistgt = np.dstack([sub['gt'] for sub in reslist])        
        mins.append(makeds(np.amin(sublistgt,axis=2), refda))
        maxs.append(makeds(np.amax(sublistgt,axis=2), refda))
    elif variable=='normative':
        mins.append(np.amin(np.dstack([sub['MQAF_LSHA'] for sub in reslist])[0,:,:], axis=1))
        mins.append(np.amin(np.dstack([sub['MQAAN_LSHA'] for sub in reslist])[0,:,:], axis=1))
        maxs.append(np.amax(np.dstack([sub['MQAF_LSHA'] for sub in reslist])[0,:,:], axis=1))
        maxs.append(np.amax(np.dstack([sub['MQAAN_LSHA'] for sub in reslist])[0,:,:], axis=1))                
    elif variable.startswith('metaswap'):
        for s in range(2):
            print(f'Get min/max for season {s+1}')
            # make a list of rasters for this season
            sublist = np.dstack([sub[s] for sub in reslist])
            mins.append(makeds(np.amin(sublist,axis=2), refda)) 
            maxs.append(makeds(np.amax(sublist,axis=2), refda))          
       
    return (mins,maxs)    
