# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 09:17:22 2021

@author: HKV lijn in water 2021
Contact: Ruud Hurkmans (hurkmans@hkv.nl)
"""

import sys
import os
from tqdm.auto import tqdm
import pandas as pd
import geopandas as gpd
import imod
import rioxarray as rio
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from hdsrhipy.groundwater.knmistationselector import KnmiStationSelector
#from knmistationselector import KnmiStationSelector
from hdsrhipy.groundwater.timeseries_modelling import *
from hdsrhipy.groundwater.add_stats import *
from hdsrhipy.groundwater.read_csv_timeseries import *

class Groundwater:    
    def __init__(self, model_path=None, name=None, export_path=None):
        if model_path is not None:            
            self.model_path = Path(model_path)
        
        if export_path is not None:            
            self.export_path = Path(export_path) / 'Grondwater'
            self.export_path.mkdir(parents=True, exist_ok=True)
        
        if name is not None:            
            self.name=name                    
            
        
    def get_heads_raster(self, mv_raster=None): 
        """ Haal de stijghoogtes op"""              
        head_folder = self.model_path / 'work' / self.name / 'output' / 'head'
        heads = imod.idf.open(head_folder / "head*_l1.idf")[:,0,:,:]       
                        
        self.heads = heads      
        
    def get_gxg_raster(self, dataset=None, mv_raster=None):
        """Bereken de GxG's"""              
        if mv_raster is None:
            mv_raster = self.model_path / 'work' / self.name / 'maaiveld' / 'maaiveld_m.idf'

        proj4_rd = "+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel +units=m +no_defs +<>"       
        
        if dataset is None:
            dataset = self.heads                
            
        # match resolution and extent of the surface level to the head-rasters
        surflevel = imod.idf.open(mv_raster)
        surflevel = surflevel.rio.write_crs(rio.crs.CRS.from_proj4(proj4_rd))    
        ref = dataset[0,:,:]       
        ref = ref.rio.set_crs(28992, inplace=True)          
        surflevel = surflevel.rio.reproject_match(ref)                   
                
        # calculate GxG's with respect to surface level
        gxg = imod.evaluate.calculate_gxg(dataset)
        
        gxg_mv = (surflevel - gxg)*100.
        return(gxg_mv)        
      
        
    def get_gt_raster(self, gxg=None):
        """ Bepaal de grondwatertrap"""              
        if gxg is None:
            gxg = self.get_gxg()
        
        gxg_mv = gxg
        gt =gxg_mv['glg'] * np.nan
        gt = xr.where((gxg_mv['ghg']<=25)&(gxg_mv['glg']<=50), 1, gt) #Ia
        gt = xr.where((gxg_mv['ghg']>25)&(gxg_mv['glg']<=50), 2, gt) #Ib
        
        gt = xr.where((gxg_mv['ghg']<=25)&(gxg_mv['glg']>50)&(gxg_mv['glg']<=80), 3, gt) #IIa
        gt = xr.where((gxg_mv['ghg']>25)&(gxg_mv['ghg']<=40)&(gxg_mv['glg']>50)&(gxg_mv['glg']<80), 4, gt)#IIb
        gt = xr.where((gxg_mv['ghg']>40)&(gxg_mv['glg']>50)&(gxg_mv['glg']<=80), 5, gt) #IIc
        
        gt = xr.where((gxg_mv['ghg']<=25)&(gxg_mv['glg']>80)&(gxg_mv['glg']<=120), 6, gt) #IIIa
        gt = xr.where((gxg_mv['ghg']>25)&(gxg_mv['ghg']<=40)&(gxg_mv['glg']>80)&(gxg_mv['glg']<=120), 7, gt)#IIIb
        
        gt = xr.where((gxg_mv['ghg']>40)&(gxg_mv['ghg']<=80)&(gxg_mv['glg']>80)&(gxg_mv['glg']<=120), 9, gt) #IVu
        gt = xr.where((gxg_mv['ghg']>80)&(gxg_mv['glg']>80)&(gxg_mv['glg']<=120), 8, gt) #IVc
        
        #gt = xr.where((gxg_mv['ghg']<=25)&(gxg_mv['glg']<=120), 10, gt) # Va
        gt = xr.where((gxg_mv['ghg']<=25)&(gxg_mv['glg']>120)&(gxg_mv['glg']<=180), 11, gt) #Vao
        gt = xr.where((gxg_mv['ghg']<=25)&(gxg_mv['glg']>180), 10, gt) #Vad
        #gt = xr.where((gxg_mv['ghg']>25)&(gxg_mv['ghg']<=40)&(gxg_mv['glg']<=120), 13, gt) #Vb
        gt = xr.where((gxg_mv['ghg']>25)&(gxg_mv['ghg']<=40)&(gxg_mv['glg']>120)&(gxg_mv['glg']<=180), 13, gt) #Vbo
        gt = xr.where((gxg_mv['ghg']>25)&(gxg_mv['ghg']<=40)&(gxg_mv['glg']>180), 12, gt) #Vbd
        
        #gt = xr.where((gxg_mv['ghg']>40)&(gxg_mv['ghg']<=80)&(gxg_mv['glg']<=120), 16, gt) 
        gt = xr.where((gxg_mv['ghg']>40)&(gxg_mv['ghg']<=80)&(gxg_mv['glg']>120)&(gxg_mv['glg']<=180), 19, gt) #VIo
        gt = xr.where((gxg_mv['ghg']>40)&(gxg_mv['ghg']<=80)&(gxg_mv['glg']>180), 14, gt) # VId
        
        #gt = xr.where((gxg_mv['ghg']>80)&(gxg_mv['ghg']<=140)&(gxg_mv['glg']<=120), 19, gt)
        gt = xr.where((gxg_mv['ghg']>80)&(gxg_mv['ghg']<=140)&(gxg_mv['glg']>120)&(gxg_mv['glg']<=180), 18, gt) #VIIo
        gt = xr.where((gxg_mv['ghg']>80)&(gxg_mv['ghg']<=140)&(gxg_mv['glg']>180), 15, gt) #VIId
        
        #gt = xr.where((gxg_mv['ghg']>140)&(gxg_mv['glg']<=120), 22, gt)
        gt = xr.where((gxg_mv['ghg']>140)&(gxg_mv['glg']>120)&(gxg_mv['glg']<=180), 17, gt) #VIIIo
        gt = xr.where((gxg_mv['ghg']>140)&(gxg_mv['glg']>180), 16, gt) #VIIId
        
        return gt
    
    def export_raster(self, dataset, filename):
        """ Exporteer een raster"""              
        dataset.rio.to_raster(self.export_path / filename)        
        
    def get_seepage(self):
        """ Haal de kwel/wegzijging op en zet om naar mm/d"""              
        seep_folder = self.model_path / 'work' / self.name / 'output' / 'bdgflf'
        seepage = imod.idf.open(seep_folder / "bdgflf*_l1.idf")[:,0,:,:]        
        return (1e3*seepage)/float(seepage['dx'])**2
        
    def seepage_season_means(self, dataset=None):
        """ Bereken seizoensgemiddelden"""              
        
        time = [pd.Timestamp(dataset['time'].values[i]).month for i in range(len(dataset['time']))]
        season_means = []
        for s in range(6):
            if s==0: months = [1,2,12]
            if s==1: months = [3,4,5]
            if s==2: months = [6,7,8]
            if s==3: months = [9,10,11]
            if s==4: months = [4,5,6,7,8,9]
            if s==5: months = [10,11,12,1,2,3]
            inds = [it for it,t in enumerate(time) if t in months]
            meanseep = dataset[inds,:,:].mean(dim='time')            
            season_means.append(meanseep)            
        return(season_means)
            
    
    def prep_riverseep(self, offset=None, offset_name=None, seizoen=None, export_path=None): 
        """ Overschrijf waarden in een raster binnen bepaalde shapfile"""              
        if export_path is None:
            export_path = self.export_path
        
        path = os.path.join(os.path.abspath('.'),'..','data','gis')
        peil = imod.idf.open(os.path.join(path, seizoen+'_PEIL_LAAG1_1.IDF'))
        peil = peil.rio.write_crs('epsg:28992', inplace=True)
        
        owshp = gpd.read_file(os.path.join(path,'openwater.shp'))
        owshp = owshp.to_crs('epsg:28992')
            
        clipped = peil.rio.clip(owshp.geometry, owshp.crs, all_touched=False, drop=False)        
        clipped = clipped+offset        
        merged = clipped.fillna(peil)        
        merged.rio.to_raster(os.path.join(export_path,seizoen+'_PEIL_LAAG1_1_'+offset_name+'.tif'))
        imod.idf.write(os.path.join(export_path, seizoen+'_PEIL_LAAG1_1_'+offset_name+'.IDF'), merged, nodata=1e+20, dtype=np.float32)
       
    def get_diffseep(self, hydromedah_path=None, name=None):
        pass
       
    def get_validation_data(self, validation_path=None, catalogue=None, timeseries=None, head_path=None):  
        """ Roep de scripts van Berendrecht Consultancy aan voor de peilbuisanalyse"""
        path = validation_path
        read_csv(path=path, catalogue_file=catalogue, timeseries_file=timeseries)   


        basename = 'catalogus_selected_layer'
        dataset='region'
        subset = 'HDSR'
    
        stationSelector = KnmiStationSelector(path)

        csvFname = os.path.join(path, 'data',basename+ '.csv')
        catalogus = pd.read_csv(csvFname, header=0)
        catalogus.drop(catalogus.columns[0], axis=1, inplace=True)
        catalogus['EIGENAAR'] = catalogus['EIGENAAR'].str.strip()

        nrow = catalogus.shape[0]
        print('Adding KNMI data')
        for idx, row in catalogus.iterrows():
            xyCrd = pd.to_numeric(row[['X_RD_CRD', 'Y_RD_CRD']]).values
            if np.all(np.isfinite(xyCrd)):
                for iType in ('prec', 'evap'):
                    iStation, stnType = stationSelector.getKnmiStation(xyCrd,iType)
                    catalogus.loc[idx, iType.upper() + '_STN']= int(iStation)
                    if iType == 'prec':
                        catalogus.loc[idx, iType.upper() + '_STNTYPE']= stnType
            else:
                catalogus.drop(idx, inplace=True)
        catalogus.to_csv(
            os.path.join(path, 'data', basename + '_knmi.csv'))
        imod.ipf.write(
            os.path.join(path, 'data', basename + '_knmi.ipf'),
            catalogus, indexcolumn=3, assoc_ext='txt', nodata=1e+20)
        
        
        knmi = "station"
        acorr = False
        baseName = "catalogus_selected_layer_knmi.ipf"

  
        ipf_head = Path(path) / 'data' / baseName

        df_merge_P_E_GW, df_GW_info = pre_process_ipf(
            path=path,
            ipf_head=ipf_head,
            knmi=knmi,
            ipf_precipitation=None,
            ipf_evapotranspiration=None,
        )

        # step 2): using create_df_statistics_summary function to create an
        # empty dataframe for the summary of results i.e. df_statistics_summary
        # per location
        df_statistics_summary = create_df_statistics_summary(df_merge_P_E_GW,
                                                             df_GW_info, acorr=acorr)
        
        # step 3): using timeseries_analysis function to perform timeseries analysis
        df_statistics_summary = timeseries_analysis(
            df_statistics_summary, knmi, df_merge_P_E_GW, Path(path), rolling=False,
            acorr=acorr, subset=subset)


        baseName = 'data/catalogus_selected_layer_knmi'        

        ipfPath = os.path.join(path)
        tsPath = os.path.join(path, 'results', 'without rolling', 'ipf')
        outPath = os.path.join(path, 'results', 'with LHM')

        if not os.path.exists(outPath):
            os.makedirs(outPath)

        if not os.path.exists(tsPath):
            os.makedirs(tsPath)

        # read netCDF files as xarray and add to list
        ncData = []        
        ncData.append(imod.idf.open(Path(head_path) / "head*_l1.idf")[:,0,:,:])
                                               
        # read csv file with piezometer data
        csvFile = os.path.join(ipfPath, baseName + '.csv')
        df = pd.read_csv(csvFile, index_col=0)

        # initialize
        for col in outList:
            df[col] = np.NaN

        control = 0

        # define grid properties
        nrows = df.shape[0]
        xcellsize = ncData[0].dx
        ycellsize = ncData[0].dy
        xMin = ncData[0].x.values[0] - 0.5 * xcellsize
        xMax = ncData[0].x.values[-1] + 0.5 * xcellsize
        yMin = ncData[0].y.values[-1] + 0.5 * ycellsize
        yMax = ncData[0].y.values[0] - 0.5 * ycellsize
        gi = {
            'xmin': xMin,
            'xmax': xMax,
            'ymin': yMin,
            'ymax': yMax,
            'xcellsize': xcellsize,
            'ycellsize': ycellsize
            }
        
        # iterate over rows to calculate stats
        # results are stored in temp file after each 100 rows
        # in case the run is interupted unintentionally
        print('Calculating validation statistics')
        for idx, row in tqdm(df.iterrows(), total=len(df)):
            control += 1
            #print ('Processing', control, '/', nrows)
            if np.isfinite(row['LHM LAYER']):
                if row['LHM LAYER'] > 1.0:
                    continue
                if (row['X_RD_CRD']>=xMin
                and row['X_RD_CRD']<=xMax
                and row['Y_RD_CRD']>=yMin
                and row['Y_RD_CRD']<=yMax):
                    df.loc[idx, outList] = addStats(tsPath, outPath, ncData, gi,row['X_RD_CRD'], row['Y_RD_CRD'], row['Buis'], row['LHM LAYER'])
        
            if int(control) == 100*int(control/100):
                #print('Storing results in between')
                df.to_csv(os.path.join(outPath, 'statistics_measured_modeled_temp'
                                       + '.csv'),
                          index=False)
                imod.ipf.write(os.path.join(outPath, 'statistics_measured_modeled_temp'
                                            + '.ipf'),
                               df, indexcolumn=3, assoc_ext='txt', nodata=1e+20)
        
        # write all results to csv and ipf
        df.to_csv(os.path.join(outPath, 'statistics_measured_modeled.csv'),
                  index=False)
        imod.ipf.write(os.path.join(outPath, 'statistics_measured_modeled.ipf'),
                       df, indexcolumn=3, assoc_ext='txt', nodata=1e+20)


    




