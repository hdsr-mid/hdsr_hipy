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
import rasterio
from rasterstats import zonal_stats
from rasterio.fill import fillnodata   
from rasterio.mask import mask
from scipy.ndimage import convolve
from shapely.geometry import mapping
import glob

class Runoff:    
    
    def __init__(self, model_path=None, name=None, export_path=None):
       
        if model_path is not None:
            self.model_path = Path(model_path)
        
        if export_path is not None:                   
            self.export_path = Path(export_path) / 'runoff'
            self.export_path.mkdir(parents=True, exist_ok=True)           
        
        self.name=name    
        
    def get_msw_var(self, variable=None):
        """ Lees data uit een Metaswap uitvoermap"""
        msw_folder = self.model_path / 'work' / self.name / 'output' / 'metaswap' / variable
        var = imod.idf.open(msw_folder / "*_l1.idf")[:,0,:,:]       
        return var
        
    def get_mean(self, dataset):
        """ Gemiddelde van dataset"""
        return dataset.mean(dim='time')
     
    def get_season_stat(self, dataset=None,stat=None):     
        """ Seizoensgemiddelden van dataset"""
        time = [pd.Timestamp(dataset['time'].values[i]).month for i in range(len(dataset['time']))]
        season_means = []
        for s in range(2):
            if s==0: months = [4,5,6,7,8,9]
            if s==1: months = [10,11,12,1,2,3]
            inds = [it for it,t in enumerate(time) if t in months]
            if stat=='mean':
                meanvar = dataset[inds,:,:].mean(dim='time')            
            elif stat=='min':
                meanvar = dataset[inds,:,:].min(dim='time')            
            else:
                print('Only mean and min currently implemented.')
            season_means.append(meanvar)            
        return(season_means)
    
    
    def aggregate_to_shapefile(self, dataset, shapefile=None,output_df = None, outcolid=None):
        """ Rasterinformatie naar vector bestanden"""
        if output_df is None:
            output_df = shapefile
        else:
            outdf = output_df
        dataset.rio.to_raster(os.path.join(self.export_path,'temp.tif'))          
        ds_per_be = zonal_stats(shapefile, os.path.join(self.export_path,'temp.tif'),stats='mean', all_touched=True)        
        outdf[outcolid] = [ds_per_be[i]['mean'] for i in range(len(ds_per_be))]
        return outdf
             
    # functie om percentielen van hoogtemodel te berekenen 
    def compute_zonal_stats(self, fp_vector, fp_raster, stats):
        """ Zonal stats voor AHN"""
        gdf = gpd.read_file(fp_vector)
        
        stats = zonal_stats(
            gdf,
            fp_raster,
            stats=stats
        )    
        
        gdf_sel = gdf[[self.zonalid, 'geometry']]
        df_concat = pd.concat((gdf_sel, pd.DataFrame(stats)), axis=1)
        return df_concat

    #    functie om uitsnede bodemhoogtemodel op te slaan 
    def get_updated_meta(self, f, out_img, out_transform):
        """ Raster meta-informatie"""
        out_meta = f.meta.copy()
        out_meta['transform'] = out_transform
        out_meta['compress'] = 'deflate'
        out_meta['height'] = out_img.shape[1]
        out_meta['width'] = out_img.shape[2]
        #out_meta['crs'] = {'init': 'epsg:28992'}
        out_meta['crs'] = None        
        return out_meta
    
    def ahn_analysis(self, ahndata=None,  bodemeenheden=None, percentage=None, bandbreedte=None, realisation=None):
        """ bepaal voor elke bodemfysische eenheid de statistieken uit een AHN  uitsnede op 0.5m"""
        # Geef id van de zone op waarbinnen de percentages berekend moet worden
        self.zonalid = 'OBJECTID' 
        
        # Geef percentages op
        if percentage is None:            
            percentage = 10
        
        # stats instellingen voor functie percentielen
        stats = ['min', 'percentile_'+str(percentage),'count']
        
        # 3x3 window voor focal mean
        weights = np.ones((3,3))
        
        # Lees water en bodemeenheden bestanden in geopandas
        bodem = gpd.read_file(bodemeenheden)
        water = gpd.read_file(os.path.join(os.path.abspath('.'), '..','data','gis','HSDR_bgt_waterdeel.shp'))
        
        # lees ahn in
        data = rasterio.open(ahndata)
        cellsize = data.transform[0]
        
        # verwijder waterdelen uit bodemeenheden
        bodemeenheden = bodem[bodem['BODEMCODE'] != '|g WATER']
        bodemeenheden = bodemeenheden.reset_index(drop=True)
        
        temp_path = os.path.join(self.export_path, '..', 'temp')
        # proces per bodemeenheid gebied
        for i in range(len(bodemeenheden)): 
            idloc = int(bodemeenheden.loc[i, self.zonalid])
            row = bodemeenheden.iloc[[i]]
            print('Eenheid: '+str(idloc))
            
            geometries = row.geometry.apply(mapping) 
            
            # clip hoogtemodel op geselecteerde bodemeenheid
            cropped_ahn, out_transform = mask(data, geometries, crop=True)
            if bandbreedte is not None:
                noise = bandbreedte[0] + np.random.rand(cropped_ahn.shape[0], cropped_ahn.shape[1],  cropped_ahn.shape[2])*(bandbreedte[1]-bandbreedte[0])    
                cropped_ahn += noise
            out_meta = self.get_updated_meta(data, cropped_ahn, out_transform)
                
            file_cropped_ahn = os.path.join(temp_path, 'cropped_ahn_' + str(idloc) + '.tif')
            with rasterio.open(file_cropped_ahn, 'w', **out_meta) as dest:
                dest.write(cropped_ahn)
            
            # Vull nodata op in geclipt hoogtemdel
            with rasterio.open(file_cropped_ahn) as src:
                profile = src.profile
                arr = src.read(1)
                cropped_ahn = np.resize(cropped_ahn, (cropped_ahn.shape[1], cropped_ahn.shape[2]))
                arr_filled = fillnodata(cropped_ahn, mask=src.read_masks(1), max_search_distance=1000, smoothing_iterations=0)
                
            # In een window van 3x3 wordt een gemiddelde waarde uitgerekend op geclipte hoogtemodel
            ahnfilled = os.path.join(temp_path, 'ahnfilled_' + str(idloc) + '.tif')
            focal_mean = convolve(arr_filled, weights) / np.sum(weights)
            
            # resultaat focal mean opslaan
            with rasterio.open(ahnfilled, 'w', **profile) as dest:
                dest.write_band(1,focal_mean)
           
            # verwijder BGT-water uit zone bodemeenheid
            out_differ = os.path.join(temp_path,'differ_' + str(idloc) + '.shp')
            gdf_differ = gpd.overlay(row, water, how='difference')
            gdf_differ.to_file(out_differ)
              
            # Verkrijg statistiek uit DTM binnen de gedefinieerde zones
            gdf_zonal_stats = self.compute_zonal_stats(out_differ, ahnfilled, stats)
                
            # Verwijder, voeg kolommen toe en bereken per zone een maaiveldcurve: hoogte(m), oppervlak(m2) en volume(m3)
            df1 = pd.DataFrame(gdf_zonal_stats.drop(columns='geometry'))
            df1 = df1.rename(columns={"min": "min_hoogte","count": "oppervlak", "percentile_"+str(percentage): "hoogte"})
            df1["oppervlak"] = df1["oppervlak"]*cellsize*cellsize
            df1["volume"] = abs(df1["hoogte"]*df1["oppervlak"])
                        
            # Resultaat samenvoegen met oorsponkelijke shape en opslaan als shapefile
            resultdef = pd.merge(bodem, df1, on=[self.zonalid, self.zonalid])
            output_file = os.path.join(temp_path,"maaiveldcurve_" + str(idloc) + ".shp")
            resultdef.to_file(output_file)
            
            # opschonen data
            os.remove(file_cropped_ahn)
            os.remove(ahnfilled)
            [os.remove(os.path.join(temp_path,f)) for f in [file for file in os.listdir(temp_path) if file.startswith('differ')]]
         
        # haal maaiveldcurve shapes op
        files = glob.glob(temp_path+'/maaiveldcurve_*.shp')
        
        # maak een lege lijst om dataframes in op te slaan
        li = []
        
        # loop through list of files and read each one into a dataframe and append to list
        for f in files:
            temp_df = gpd.read_file(f)
            li.append(temp_df)
        
        # concate alle dataframes
        df = pd.concat(li, axis=0)
        df.drop([col for col in df.columns if col not in ['OBJECTID','Shape_Area','min_hoogte','hoogte', 'volume','geometry']], axis=1, inplace=True)
        
        [os.remove(os.path.join(temp_path,f)) for f in [file for file in os.listdir(temp_path) if file.startswith('maaiveldcurve')]]

        # schrijf resultaat als shapefile
        if realisation is None:
            df.to_file(os.path.join(self.export_path, 'hdsr_maaiveldcurves.shp'))
        else:
            df.to_file(os.path.join(self.export_path, 'hdsr_maaiveldcurves_'+str(realisation)+'.shp'))
        
        
            