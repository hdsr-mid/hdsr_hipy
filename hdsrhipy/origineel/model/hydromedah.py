# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 09:17:22 2021

@author: HKV lijn in water 2021
Contact: Ruud Hurkmans (hurkmans@hkv.nl)
"""

import imod
from hdsrhipy import Runfile
from hdsrhipy.model import rasters
import hdsrhipy
import sys
import os
import pandas as pd
import geopandas as gpd
from hdsrhipy.model import util
import shutil
import numpy as np
from datetime import datetime as dt
from datetime import timedelta as td
from tqdm.auto import tqdm
from pathlib import Path

class Hydromedah:
    
    """ Methods to set-up and run the Hydromedah model"""
    def __init__(self, hydromedah_path=None, name=None, precip_path=None, evap_path=None):
        
        if hydromedah_path is not None:
            self.data_path = hydromedah_path    
        else:
            return('A path to Hydromedah is required.')
                
        if precip_path is None:
            precip_path = hydromedah_path
        
        if evap_path is None:
            evap_path = hydromedah_path
        
        self.name = name
        self.rf = Runfile(os.path.join(self.data_path,'hdsr_ns.run'), data_path='$DBASE$\\', evap_path=evap_path, precip_path=precip_path)
        
     
    def setup_model(self, start_date=None, end_date=None, resolution=None, metaswap_vars=None, add_surface_water=None, afw_shape = None, firstrun=False):
        """ Functie om Hydromedah run voor te bereiden alle (tijd)-informatie goed te zetten"""
        
        if resolution is None:
            resolution = 250.
        if metaswap_vars is None:
            metaswap_vars = ['ETact','S01','Ssd01', 'qinf']
        
        if firstrun:
            util.prepare_raw_data(self.data_path)                
        
        self.rf.to_version(4.3)
        self.rf.update_metaswap(datadir=self.data_path, start_date=pd.to_datetime(start_date), end_date=pd.to_datetime(end_date), metaswap_vars=metaswap_vars)
        
        # make the calculation grid smaller (from [105000,173000,433000,473000])
        self.rf.data['XMIN'] = 105000.0
        self.rf.data['XMAX'] = 173000.0
        self.rf.data['YMIN'] = 433000.0
        self.rf.data['YMAX'] = 473000.0

        # make the calculation grid finer (from 250 to 100 m)
        self.rf.data['CSIZE'] = resolution
    
        self.rf.change_period(start_date, end_date)
                
        # add surface water from peilgebieden to the runfile
        if add_surface_water:
            gdf_pg = gpd.read_file(os.path.join(self.data_path, afw_shape+'.shp'))
            gdf_pg['CODENR'] = gdf_pg.index + 1                  
            util.add_simgro_surface_water(self.rf, gdf_pg=gdf_pg, run1period=True, datadir=self.data_path)
        
         # save output every timestep (every day)
        self.rf.data['ISAVE'][:] = 1          
        
                
    def run_model(self, model_path=None, use_summerlevel=None, use_winterlevel=None, silent=False):       
        """ Functie om Hydromedah door te rekenen"""
        if model_path is None:
            model_path = self.data_path
        model_path = os.path.join(model_path, 'work', self.name)                
        self.rf.run_imodflow(model_path, self.name, data_path=self.data_path, use_summerlevel=use_summerlevel, use_winterlevel=use_winterlevel, silent=silent)
        

    # inlezen .key-bestand
    def ReadKey(self, filename):
        """ Functie om een SIMGRO key file in te lezen"""
        kf = open(filename+'.key','r')
        lines = kf.readlines()
        kf.close()
        for iline,line in enumerate(lines):
            if line.startswith('*'): continue
            if line.startswith('FORMAT'):
                bytesize = int(line.split('=')[1].strip()[1])
            if line.startswith('PERIOD'):
                period = int(line.split('=')[1].strip())
            if line.startswith('NUMVAR'):
                numvar = int(line.split('=')[1].strip())        
                varlines = lines[iline+3:iline+3+numvar]
                variables = [var.split(' ')[0].strip() for var in varlines]            
            if line.startswith('NUMPTS'):
                numlocs = int(line.split('=')[1].strip())
                loclines = lines[iline+3:iline+3+numlocs]
                locations = [loc.lstrip().split(' ')[0] for loc in loclines]
        return locations, variables, period, bytesize
        
    # inlezen .TIM bestand
    def ReadTim(self,filename):
        """ Functie om een SIMGRO tim file in te lezen"""
        
        tf =open(filename+'.tim','r')
        lines = tf.readlines()
        tf.close()
        times = [(float(line.lstrip().split(' ')[0]), int(line.lstrip().split(' ')[2]))  for line in lines]
        dates = [dt.strptime(str(time[1])+'/01/01 00:00:00','%Y/%m/%d %H:%M:%S') + td(days=time[0]) for time in times]
        dates = [date+td(0,1) if date.minute==59 else date for date in dates]
        return dates
    
    # haal tijdserie op een locatie op
    def GetTimeSeries(self, f,numtim,numlocs,numvar,loc, var ):    
        """ Functie om een tijdserie op te halen uit een binair SIMGRO bestand (oud, V2 wordt gebruikt)"""
        f.seek(loc*numvar*4+var*4,0)
        timeserie = []
        for tim in range(numtim):
            val = float(np.fromfile(f,dtype=np.float32,count=1))
            timeserie.append(val)
            f.seek(numlocs*numvar*4-4,1)       
        return timeserie
    
    def GetTimeSeries_V2(self,f,numtim,numlocs,numvar,loc,var,bytesize,timeindexfilter=None):
        """ Functie om een tijdserie op te halen uit een binair SIMGRO bestand"""
        byte_index = range(var * bytesize, numtim * numlocs * numvar * bytesize, numlocs * numvar * bytesize)
        byte_index = [x + (loc * bytesize * numvar) for x in byte_index]
        # apply timeseries filter
        if timeindexfilter is not None:
            byte_index_filtered = [byte_index[x] for x in timeindexfilter]
            byte_index = list(byte_index_filtered)
        # get the values at the byte indices
        timeserie = [None] * len(byte_index)
        for tim in range(len(byte_index)):
            f.seek(byte_index[tim],0)
            val = float(np.fromfile(f,dtype=np.float32,count=1))
            timeserie[tim] = val
        return timeserie    

    def Read_BDA(self, simgro_path,filename):      
        """ Functie om een binair SIMGRO bestand uit te lezen (oud, niet gebruikt)"""
        # maak lijsten van tijdstappen, variabelen en locaties
        if (filename=='sw_dtsw')|(filename=='sw_dtgw'): 
            var1 = 'Vdsreg'
            var2 = 'Vsues'            
        (locations, period, period, bytesize) =  self.ReadKey(os.path.join(simgro_path,filename))
        (timesteps) = self.ReadTim(os.path.join(simgro_path,filename))
        if timesteps[0].second==59:
           timesteps[0] = timesteps[0]+td(0,1)
        if timesteps[1].second==59:
            timesteps[1] = timesteps[1]+td(0,1)
        tijdstap = np.floor((timesteps[1]-timesteps[0]).total_seconds())
        varind1 = variables.index(var1)
        varind2 = variables.index(var2)

        numpts = len(locations)
        numvar = len(variables)
        numtim = len(timesteps)

        df = pd.DataFrame(np.zeros((numtim-1, numpts))*np.nan)
        df.columns = locations
        df.index = pd.date_range(timesteps[1], timesteps[-1], freq='D')
        df2 = pd.DataFrame(np.zeros((numtim-1, ))*np.nan)
        for iloc,location in tqdm(enumerate(locations),total=numpts):
            locind = locations.index(location)        
            f = open(os.path.join(simgro_path,filename+'.bda'), "rb")            
            timeserie1 = self.GetTimeSeries(f, numtim,numpts,numvar,locind, varind1)[:-1]                           
            timeserie2 = self.GetTimeSeries(f, numtim,numpts,numvar,locind, varind2)[:-1]                           
            df[locations[locind]] = -timeserie1 - timeserie2
            # if iloc==0:
            #     df = pd.DataFrame([dt.strftime(time,"'%Y-%m-%d;%H:%M:%S'") for time in timesteps], columns=['Time'])           
            #     df['Qdssw_'+location] = [f' {tijd:.6E}' for tijd in tijdserie]                    
            # f.close()                  
        return(df)    
    
    def Read_BDA_V2(self,simgro_path,filename,variables, bytesize=None, year_filter=None, locind=None):     
        """ Functie om een binair SIMGRO bestand uit te lezen"""        
        (locations, variabelen, period, key_bytesize) =  self.ReadKey(os.path.join(simgro_path,filename))        
        (timesteps) = self.ReadTim(os.path.join(simgro_path,filename))
        if bytesize is None:
            bytesize = key_bitesize
        tijdstap = np.floor((timesteps[1]-timesteps[0]).total_seconds())
        if tijdstap == 3600.:
            frequency = 'H'
        if tijdstap == 86400.:
            frequency = 'D'
        varind = [variabelen.index(var) for var in variables]
        
        # create neat time series from, with correct interval and proper length
        # takes into accoutn period type (as specified in the key file)
        if period == 1:
            timesteps = [timesteps[0] + td(seconds=tijdstap*x) for x in range(len(timesteps)-1)]
        else:
            timesteps = [timesteps[0] + td(seconds=tijdstap*x) for x in range(len(timesteps))]
        timesteps = list(timesteps)
        
        # get number of locations, variables and timesteps
        numpts = len(locations)
        numvar = len(variabelen)
        numtim = len(timesteps)
        
        # apply year filter
        if year_filter is not None:
            timesteps_years = [dt.strftime(tim, '%Y') for tim in timesteps]
            timefilter_index = [i for i, x in enumerate(timesteps_years) if x == year_filter]
            timesteps = list(map(timesteps.__getitem__, timefilter_index))
        else:
            timefilter_index = None
            
        # loop over all locations
        df = pd.DataFrame(np.zeros((numtim, numpts))*np.nan)
        df.columns = locations
        df.index = pd.date_range(timesteps[0], timesteps[-1], freq=frequency)
        # get time series for a specific location (index)
        if locind is not None:
            locations = [locind]
        for location in tqdm(locations, total=numpts):
            if location == '0':
                # skip swnr 0
                continue
            locind = locations.index(location)
            ts = []
            for vi in varind:
                f = open(os.path.join(simgro_path,filename+'.bda'), "rb")
                ts.append(self.GetTimeSeries_V2(f, numtim, numpts, numvar, locind, vi, bytesize=bytesize, timeindexfilter=timefilter_index))
                f.close()
            df.loc[:,location] = [(-ts[0][i] - ts[1][i])/tijdstap for i in range(len(ts[0]))]
        return(df)

    def read_laterals(self, model_path=None, model_name=None, msw_file=None):
        """ Lateralen uitlezen uit SSIMGRO """
        if msw_file is None: 
            msw_file = 'sw_dtgw'
        if model_path is None:
            model_path = self.data_path
        if model_name is None:
            model_name = self.name
        df = self.Read_BDA_V2(os.path.join(model_path,'work', model_name, 'output','metaswap'),msw_file,['Vdsreg','Vsues'], bytesize=8)             
        return df

    def cleanup(self,model_path=None, name=None, modflow_vars=None, modflow_layers=None, metaswap_files=None):    
        """Opschonen Hydromedah-uitvoer"""
        if model_path is None:
            model_path = self.model_path
        if name is None:
            name = self.name
        if modflow_vars is None:
            print('Are you sure? Specify something to keep otherwise all will be deleted!')
            return
        pad = os.path.join(model_path,'work',name,'output')
        dirlist = os.listdir(pad)
        modflow_vars = modflow_vars + ['mf2005_tmp','metaswap']
        for folder in dirlist:
            if (folder not in modflow_vars):
                if os.path.isdir(os.path.join(pad,folder)):
                    print('Deleting '+folder+'...')
                    shutil.rmtree(os.path.join(pad,folder))
            elif folder == 'metaswap':
                filelist = os.listdir(os.path.join(pad,folder))
                print('Removing large csv files from metaswap...')
                #os.remove(os.path.join(pad,folder,'swmdl1_sw_dtgw.csv'))
                #os.remove(os.path.join(pad,folder,'swmdl1_sw_per.csv'))
                bdalist = [f for f in filelist if f.endswith('bda')]
                for f in bdalist:
                    print('Removing '+f+' (and key/tim) from metaswap...')
                    if f.strip('.bda') not in metaswap_files:                
                        os.remove(os.path.join(pad,folder,f))
                        os.remove(os.path.join(pad,folder,f.strip('.bda')+'.key'))
                        os.remove(os.path.join(pad,folder,f.strip('.bda')+'.tim'))                
            else:
                filelist = os.listdir(os.path.join(pad,folder))
                if filelist[0].endswith('idf'):
                    fends = [f'_l{lay}.idf' for lay in modflow_layers]
                    for layer in modflow_layers:
                        flist = [f for f in filelist if not f[-7:] in fends]
                        print('Removing unnecessary layers from '+folder+'...')
                        [os.remove(os.path.join(pad,folder,f)) for f in flist]
                        