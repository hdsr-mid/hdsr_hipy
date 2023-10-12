# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 09:17:22 2021

@author: HKV lijn in water 2021
Contact: Ruud Hurkmans (hurkmans@hkv.nl)
"""
import sys
import os

package_path = os.path.abspath('../../')
sys.path.append(package_path) 
import pandas as pd
import geopandas as gpd
from hdsrhipy.model import util
from hdsrhipy import Groundwater
import shutil
import numpy as np
from datetime import datetime as dt
from datetime import timedelta as td
from tqdm.auto import tqdm
from pathlib import Path
from rasterstats import zonal_stats
import matplotlib.pyplot as plt

class Maatgevend: 
   
    def __init__(self, name=None, export_path=None, afw_shape=None):        
        
        if export_path is not None:
            self.export_path = Path(export_path) / 'Maatgevend'
            self.export_path.mkdir(parents=True, exist_ok=True)
        
        if afw_shape is None:
            afw_shape = os.path.join(os.path.abspath('.'), '..','data','gis','Afwateringseenheden.shp')       
            
        self.name=name                     
            
        self.afw = gpd.read_file(afw_shape)
        self.afw['CODENR'] = self.afw.index + 1       
                
        
    def subtract_seepage(self, mean_seepage = True, model_path = None, name=None):     
        """ Trek kwel af van de afvoer"""
        if name is None:
            name = self.name
        gw = Groundwater(model_path=model_path, name=name)
        seep = gw.get_seepage()                
        
        self.laterals_nosp = self.laterals.copy(deep=True)
        if mean_seepage:
            # take temporal average from seepage and convert it from mm/d to m3/s , like the laterals
            av_seep  =  seep.mean(dim='time') 
            av_seep.rio.to_raster(self.export_path / 'temp.tif')
            sp_per_afw = zonal_stats( self.afw, self.export_path / 'temp.tif',stats='mean', all_touched=True)        
            seep = [sp_per_afw[int(j)-1]['mean']*self.afw.geometry.area.iloc[int(j)-1] /  (1000. * 86400.) for j in self.laterals.columns]
            for i in tqdm(range(self.laterals_nosp.shape[0]), total=self.laterals_nosp.shape[0]):                 
                # print( 'Processing '+str(i)+' of '+str(self.laterals_nosp.shape[0]))
                self.laterals_nosp.iloc[i,:] = self.laterals.iloc[i,:] - seep 
                # for jj,j in enumerate(self.laterals.columns):
                #     self.laterals_nosp.iloc[i,jj] = self.laterals.iloc[i,jj] - sp_per_afw[int(j)-1]['mean'] * \
                #                                         self.afw.geometry.area.iloc[int(j)-1] / (1000. * 86400.)                        
            self.laterals_nosp.to_csv(self.export_path / str('Laterals_nosp_av_'+self.name+'.csv'), sep=",")
  
        else:
            for i in tqdm(range(self.laterals_nosp.shape[0]), total=self.laterals_nosp.shape[0]):                        
                # take temporal average from seepage and convert it from mm/d to m3/s , like the laterals
                rast = seep.isel({'time':i})
                rast.rio.to_raster(self.export_path / 'temp.tif')     
                sp_per_afw = zonal_stats( self.afw, self.export_path / 'temp.tif',stats='mean', all_touched=True)
                sp = [sp_per_afw[int(j)-1]['mean']*self.afw.geometry.area.iloc[int(j)-1] /  (1000. * 86400.) for j in self.laterals.columns]
                self.laterals_nosp.iloc[i,:] = self.laterals.iloc[i,:] - sp 
                #for jj,j in enumerate(self.laterals.columns):
                #    self.laterals_nosp.iloc[i,jj] = self.laterals.iloc[i,jj] - sp_per_afw[int(j)-1]['mean'] * \
                #                                        self.afw.geometry.area.iloc[int(j)-1] / (1000. * 86400.)                  
                                                        
            self.laterals_nosp.to_csv(self.export_path / str('Laterals_nosp_dyn_'+self.name+'.csv'), sep=",")
        
    def get_laterals(self, seepage_subtracted=True, mean_seepage=True):    
        """ Haal de lateralen op uit de CSV, al dan niet met kwel """
        if seepage_subtracted:
            if mean_seepage:
                self.laterals_nosp = pd.read_csv(self.export_path / str('Laterals_nosp_av_'+self.name+'.csv'), sep=",")
            else:
                self.laterals_nosp = pd.read_csv(self.export_path / str('Laterals_nosp_dyn_'+self.name+'.csv'), sep=",")
            self.laterals_nosp.index = self.laterals_nosp.iloc[:,0] 
            self.laterals_nosp.drop(self.laterals_nosp.columns[0:2], axis=1,inplace=True)            
        else:
            self.laterals = pd.read_csv(self.export_path / str('Laterals_'+self.name+'.csv'), sep=",")
            self.laterals.index = self.laterals.iloc[:,0] 
            self.laterals.drop(self.laterals.columns[0:2], axis=1,inplace=True)
            
    def process_timeseries(self, ts):
        """Bereken de maatgevende afvoer/aanvoer uit een tijdreeks van netto fluxen"""
        times = [pd.Timestamp(ts.index[i]) for i in range(len(ts.index))]
        years = set([pd.Timestamp(ts.index[i]).year for i in range(len(ts.index))])
        timestep = (times[1]-times[0]).total_seconds()/3600.
        ts.iloc[0:3] = np.nan
        ts_af = ts.copy(deep=True)
        ts_af[ts_af < 0.] = 0.                              
        rts1 = ts_af.rolling(int(timestep), closed='both').mean()                    
        rts1[np.isnan(rts1)] = 0. 
        ts_aan = ts.copy(deep=True)                         
        ts_aan[ts_aan > 0. ] = 0
        ts_aan = ts_aan * -1.
        rts10 = ts_aan.rolling(int(10.*timestep), closed='both').mean()
        rts10[np.isnan(rts10)] = 0.  
        mqaf = rts1.sort_values(ascending=True).iloc[int(-1.5*len(years))]        
        mqaan = rts10.sort_values(ascending=True).iloc[int(-0.1*len(years))]
        return [rts1, mqaf, rts10, mqaan]

    def plot_location(self, nr, seepage_subtracted = True):         
        """ plot alle informatie voor een locatie """
        afwid = self.afw.loc[self.afw.CODENR==int(nr),'CODE'].to_string(index=False)
        area = float(self.afw.loc[self.afw.CODENR==int(nr),'geometry'].area)
        m3s_to_mmd = area/(1000.*86400.)
        m3s_to_lsha = area/(1000.*10000.)
        print('Plotting '+afwid)
        ts = self.laterals.loc[:,str(nr)]         
        if seepage_subtracted:
            tsc = self.laterals_nosp.loc[:,str(nr)]            
        fig, axs = plt.subplots(3)
        times = [pd.Timestamp(ts.index[i]) for i in range(len(ts.index))]
        years = set([pd.Timestamp(ts.index[i]).year for i in range(len(ts.index))])
        fig.suptitle('Maatgevende aan/afvoer voor '+afwid)
        fig.set_size_inches(8,8)
        axs[0].plot(times, ts, color='blue',label='Aan/afvoer')        
        if seepage_subtracted:
            axs[0].plot(times, tsc, color='red',label='Aan/afvoer zonder kwel')
        axs[0].plot(times, ts*0., color='black')
        axs[0].set_ylabel('Af/aanvoer [$\mathregular{m^{3} s^{-1}}$]')
        axs[0].legend(ncol=1)
        axs[1].plot(range(len(times)), self.process_timeseries(ts)[0].sort_values(ascending=True)/m3s_to_lsha,label='Afvoer gesorteerd')
        axs[1].plot(len(times)-int(1.5*len(years)), self.process_timeseries(ts)[1]/m3s_to_lsha, 'bo', label='Maatgevende afvoer')
        if seepage_subtracted:
            axs[1].plot(range(len(times)), self.process_timeseries(tsc)[0].sort_values(ascending=True)/m3s_to_lsha,linestyle='dashed', color='blue',label='Afvoer gesorteerd zonder kwel')
            axs[1].plot(len(times)-int(1.5*len(years)), self.process_timeseries(tsc)[1]/m3s_to_lsha, 'bs', label='Maatgevende afvoer zonder kwel')                    
        axs[1].legend(ncol=1)
        axs[1].set_ylabel('Afvoer [$\mathregular{l s^{-1} ha^{-1}}$]')
        axs[2].plot(range(len(times)),self.process_timeseries(ts)[2].sort_values(ascending=True)/m3s_to_lsha,color='red',label='Aanvoer gesorteerd')
        axs[2].plot(len(times)-int(0.1*len(years)), self.process_timeseries(ts)[3]/m3s_to_lsha, 'o',color='red', label='Maatgevende aanvoer')
        if seepage_subtracted:
            axs[2].plot(range(len(times)),self.process_timeseries(tsc)[2].sort_values(ascending=True)/m3s_to_lsha,linestyle='dashed',color='red',label='Aanvoer gesorteerd zonder kwel')
            axs[2].plot(len(times)-int(0.1*len(years)), self.process_timeseries(tsc)[3]/m3s_to_lsha, 's',color='red', label='Maatgevende aanvoer zonder kwel')        
        axs[2].set_ylabel('Aanvoer [$\mathregular{l s^{-1} ha^{-1}}$]')
        axs[2].legend(ncol=1)
        plt.savefig(self.export_path /  str('MG_'+afwid+'_'+self.name+'.jpg'))
                   
    def get_q_norm(self, dataset=None):
        """bereken maatgevende afvoer/aanvoer"""
        
        mg_q = self.afw.copy(deep=True)
        keep = ['CODE','geometry','CODENR']
        for k in self.afw.columns:
            if k not in keep:
                mg_q.drop([k],axis=1,inplace=True) 
        mg_q['MQAF_M3S'] = np.nan
        mg_q['MQAF_LSHA'] = np.nan
        mg_q['MQAAN_M3S'] = np.nan
        mg_q['MQAAN_LSHA'] = np.nan
        mg_q.index = self.afw.index            
        for ind,i in tqdm(self.afw.iterrows(), total=len(self.afw)):
            m3s_to_mmd = i.geometry.area/(1000.*86400.)
            m3s_to_lsha = i.geometry.area/(10000.*1000.)
            if str(i.CODENR) in dataset.columns.to_list():        
                ts = dataset.loc[:,str(i.CODENR)]                        
                (_,mqaf, _, mqaan) = self.process_timeseries(ts)                                
                mg_q.loc[mg_q['CODENR']==i.CODENR, 'MQAF_M3S'] = mqaf
                mg_q.loc[mg_q['CODENR']==i.CODENR, 'MQAF_LSHA'] = mqaf/m3s_to_lsha
            
                mg_q.loc[mg_q['CODENR']==i.CODENR, 'MQAAN_M3S'] = mqaan
                mg_q.loc[mg_q['CODENR']==i.CODENR, 'MQAAN_LSHA'] = mqaan/m3s_to_lsha
        return mg_q
        
    def export_shp(self, dataset, filename):
        """Exporteer een shapefile"""
        filename =self.export_path / str(filename)
        dataset.to_file(filename)
                
    def get_validation(self, afgid=None, afvoer=None, aanvoer=None):                
        """ Vergelijk met WIS data"""
        # read the shapefile of afvoergeibeden and extract the right one
        afg_shape = os.path.join(os.path.abspath('.'), '..','data','gis','Afvoergebieden.shp')                
        afg = gpd.read_file(afg_shape)
        gebied = [i for _,i in afg.iterrows() if str(i.GAFCODE)==afgid]
        
        # make subset of intersecting afwateringseenheden
        sub = self.afw[self.afw.geometry.intersects(gebied[0].geometry)]
        # get the laterals        
        if not hasattr(self, 'laterals'):            
            self.get_laterals(seepage_subtracted=False)
        # sum the flux in m3/s
        num = 0
        for i in sub.CODENR:
            if str(i) in self.laterals.columns:        
                if num==0:
                    netto_hydromedah = self.laterals.loc[:,str(i)]
                else:
                    netto_hydromedah += self.laterals.loc[:,str(i)]
                num+=1
        # make a nice dataframe 
        netto_hydromedah.index = [pd.Timestamp(netto_hydromedah.index[i]) for i in range(len(netto_hydromedah))]
        
        # convert the weird WIS datato floats for this afvoergebied and get the netto flux in dataframe
        af_wis = [float(i) for i in afvoer.loc[:,afgid]]
        aan_wis = [float(i) for i in aanvoer.loc[:,afgid]]
        netto = pd.DataFrame([af_wis[i] - aan_wis[i] for i in range(len(af_wis))], columns=['WIS'], index = afvoer.index)
        # convert to m3/s, like the laterals
        netto['WIS'] = netto['WIS']*float(afg[afg.GAFCODE==afgid].geometry.area)/(1000. * 86400.)
        # merge add the right timesteps
        netto = netto.merge(netto_hydromedah,left_index=True, right_index=True)
        netto.columns = ['WIS','HYDROMEDAH']
        # get the error
        netto['ERROR'] = netto['WIS'] - netto['HYDROMEDAH']
        return netto