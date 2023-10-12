# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 09:17:22 2021

@author: hurkmans
"""

from datetime import datetime
from pathlib import Path
import requests
import os
import sys
import time
import json
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import pandas as pd
import rioxarray as rio
import imod
import shutil
import logging
logger = logging.getLogger(__name__)

#%$
class Meteorology:
    """ 
    Class to assemble meteorological forcing for the models and perform climate operations
    """
    def __init__(self):
        self.proj = "+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel +units=m +no_defs +<>"

    def write2ascii(self, ds, variable=None, dsvar=None, path=None, timestep=None, scaling_factor=None):
        """
        Write a dataset to ascii-files
        --------------------------------------------------------------------------------------------
        Input: 
            -ds: dataset
            -variable: name to determine file path and name        
            -path: location to put the files
        --------------------------------------------------------------------------------------------
        Output:
            None
        --------------------------------------------------------------------------------------------
        """
        if dsvar is None:
            dsvar = 'variable'
        if scaling_factor is None:
            scaling_factor = 1.0
        forcing_path = path / Path('forcing',variable)        
        forcing_path.mkdir(parents=True, exist_ok=True)        
        ds = ds.where(ds[dsvar] < 1e308)
        ds[dsvar] = ds[dsvar] * scaling_factor
        ds = ds.fillna(-999)
        for i in range(len(ds['time'])):
            arr = ds[dsvar].isel(time=i)
            if timestep=='Hours':
                time='h'
            if timestep=='Days':
                time = 'D'
            time = np.datetime_as_string(arr['time'], unit=time)    
            #time = np.datetime_as_string(arr['time'], unit=time)#[0:10]+np.datetime_as_string(arr['time'], unit='m')[11:13]+np.datetime_as_string(arr['time'], unit='m')[14:16]
            fn = variable+'_'+time+'.asc'    
            imod.rasterio.write( forcing_path / fn, arr, nodata=-999)
         
            
    def from_netcdf(self, input_folder, output_folder, variable, t_unit, dsvar=None, scaling_factor=None):
        if isinstance(input_folder,str):
            input_folder = Path(input_folder)
        file_list = [input_folder/f for f in os.listdir(input_folder) if f.endswith('.nc')]
        for file in file_list:            
            ds = xr.open_dataset(file)
            ds = ds.rio.write_crs(rio.crs.CRS.from_proj4(self.proj))    
            #get a reference raster to match IRC with        
            reff = Path('\\'.join(os.path.abspath(__file__).split('\\')[:-1])) / 'precipitation_meteobase.asc'        
            refrast = imod.rasterio.open(reff)             
            refrast = refrast.rio.set_crs(28992, inplace=True)        
            ds_rd = ds.rio.reproject_match(refrast)                         
            self.write2ascii(ds_rd, variable, dsvar=dsvar, path=output_folder, timestep=t_unit, scaling_factor=scaling_factor)
         
    def select_key(self,dataset_name):
        if dataset_name=='RD85WL':
            return('eyJvcmciOiI1ZTU1NGUxOTI3NGE5NjAwMDEyYTNlYjEiLCJpZCI6ImEyNjkyNGE2YTBiNDQ5N2I5MzQzODkyMjQxOGU5ZDQ2IiwiaCI6Im11cm11cjEyOCJ9')
        elif dataset_name=='EV85WL':
            return('eyJvcmciOiI1ZTU1NGUxOTI3NGE5NjAwMDEyYTNlYjEiLCJpZCI6IjIwMGNjZDBjYWM0ZTQ2NDc5YjgxYTg5Mjk4OGM4ZjZiIiwiaCI6Im11cm11cjEyOCJ9')
        elif dataset_name=='RD85WH':
            return('eyJvcmciOiI1ZTU1NGUxOTI3NGE5NjAwMDEyYTNlYjEiLCJpZCI6IjQwYzJiOGNkOWVmMTQ1ZDhiNjcxMTRmNzkzMDEzMTBjIiwiaCI6Im11cm11cjEyOCJ9')        
        elif dataset_name=='EV85WH':
            return('eyJvcmciOiI1ZTU1NGUxOTI3NGE5NjAwMDEyYTNlYjEiLCJpZCI6ImU1Mjg4NjZhMjMyNzRhMDlhN2M5ZGQzZjIwNGVkYmY4IiwiaCI6Im11cm11cjEyOCJ9')
        elif dataset_name=='RD85GL':
            return('eyJvcmciOiI1ZTU1NGUxOTI3NGE5NjAwMDEyYTNlYjEiLCJpZCI6Ijg4YmI3YWFjMGI0MjRmZWNiY2MyMmI5MTNhMDhmMTVmIiwiaCI6Im11cm11cjEyOCJ9')
        elif dataset_name=='EV85GL':
            return('eyJvcmciOiI1ZTU1NGUxOTI3NGE5NjAwMDEyYTNlYjEiLCJpZCI6IjU4Y2I1MTNmOTM2OTQ4MTI5NGJhYmQ2NGE0YTVlMzcyIiwiaCI6Im11cm11cjEyOCJ9')
        elif dataset_name=='RD85GH':
            return('eyJvcmciOiI1ZTU1NGUxOTI3NGE5NjAwMDEyYTNlYjEiLCJpZCI6IjNlZTM4ZjY5MjIwMzRmMDQ5MzIwMmNjZTExMjJlOTdhIiwiaCI6Im11cm11cjEyOCJ9')
        elif dataset_name=='EV85GH':
            return('eyJvcmciOiI1ZTU1NGUxOTI3NGE5NjAwMDEyYTNlYjEiLCJpZCI6ImM3Njk3MWVlOTMyNjQ0MTJhYzllZTA0ZDU5NjFmZjI4IiwiaCI6Im11cm11cjEyOCJ9')
        
        
    def download_from_KNMI(self, api_key=None, variable=None, dataset_name=None, dataset_version=None, max_keys=None, start_after=None, download_path=None):
        """
        Get data from the KNMI data portal
        --------------------------------------------------------------------------------------------
        Input: 
            -api_key: can be used to get specific key obtained from KNMI. If ommitted, the default open-data key (with limited functionality) is used;
            -dataset_name: name of the dataset as specified in the KNMI data catalogue;
            -dataset_version: version of the dataset as specified in the KNMI data catalogue;
            -max_keys: paramter for the API, default is 10000;
            -start_after: if only part the dataste is needed, start downloading after this filename. Default is dataset start;
            -download_path: location where the results are downloaded; filanames are taken from KNMI.
        --------------------------------------------------------------------------------------------
        Output:
            None
        --------------------------------------------------------------------------------------------
        """
        
        api_url = "https://api.dataplatform.knmi.nl/open-data"
        
        if api_key is None: 
            #api_key = self.select_key(dataset_name)
            # open data key:
            api_key = 'eyJvcmciOiI1ZTU1NGUxOTI3NGE5NjAwMDEyYTNlYjEiLCJpZCI6ImNjOWE2YjM3ZjVhODQwMDZiMWIzZGIzZDRjYzVjODFiIiwiaCI6Im11cm11cjEyOCJ9'
        if max_keys is None:
            max_keys = 10000
        if start_after is None:
            start_after = ' '
        
        
        session = requests.Session()
        session.headers.update({"Authorization": api_key})

        #start_after_filename_prefix = 'INTER_OPER_R___EV24_H__L3__19110515T000000_19110516T000000_0001'
        list_files_response =session.get(f"{api_url}/v1/datasets/{dataset_name}/versions/{dataset_version}/files",
                                           headers={"Authorization": api_key},
                                           params={"maxKeys": max_keys,
                                                   "startAfterFilename": start_after})
        list_files = list_files_response.json()
        dataset_files = list_files.get("files")
        
        temp_path = download_path / 'temp'
        temp_path.mkdir(parents=True, exist_ok=True)
            
        # Retrieve first file in the list files response
        for dataset_file in dataset_files:        
            if not os.path.exists(temp_path / dataset_file['filename']):
                logger.info('Downloading '+dataset_file['filename'])
                filename = dataset_file.get("filename")
                endpoint = f"{api_url}/datasets/{dataset_name}/versions/{dataset_version}/files/{filename}/url"
                #get_file_response = requests.get(endpoint, headers={"Authorization": api_key})
                get_file_response = session.get(endpoint)
                download_url = get_file_response.json().get("temporaryDownloadUrl")
                dataset_file_response = requests.get(download_url)
    
                # Write dataset file to disk
                p = temp_path / filename
                p.write_bytes(dataset_file_response.content)
            else: 
                logger.info(str(temp_path / dataset_file['filename'])+' exists.')
            ds = xr.open_dataset(temp_path / dataset_file['filename'])
            ds = ds.rename({'prediction':'variable'})
            self.write2ascii(ds, variable=variable, path=download_path, timestep='Days')        
            ds.close()
            #os.remove(temp_path / dataset_file['filename'])
            
    def download_from_WIWB(self, credentials=None, datasource=None, variable=None, start=None, end=None, timestep=None, extent=None, download_path=None):
        """
        Get data from the WIWB API. Only grid data can be downloaded (for now).
        --------------------------------------------------------------------------------------------
        Input: 
            -credentials: (username, password)-tuple to gain access to WIWB. Note that also the IP-address should be whitelisted by Hydrologic;
            -datasource: dataset name that is recognized by WIWB;
            -variable: variable string that is part of that dataset
            -start_time: start time (YYYYMMDDHHMMSS)
            -end_time: end time (YYYYMMDDHHMMSS)
            -timestep: for now only 'day' and 'hour' are implemented;
            -extent: extent for which data should be obtained (in RD-coordinates: [xll,yll,xur,yur]). By default the HYDROMEDAH extent is used.           
            -download_path: location for the resulting ascii rasters
        --------------------------------------------------------------------------------------------
        Output:
            None
        --------------------------------------------------------------------------------------------
        """
   
        # API settings
        api_url = 'https://wiwb.hydronet.com/api/'
        if credentials is None:
            return('Credentials are needed, sorry.')            
            
        #if 'Radar' in datasource:            
        url = api_url+'/grids/get'                
        
        if timestep[1]=='D':
            t_unit = 'Days'
            t_val =   int(timestep[0])
        elif timestep[1]=='H':
            t_unit = 'Hours'
            t_val =   int(timestep[0])
        elif timestep[1:]=='min':
            t_unit = 'Minutes'
            t_val = int(timestep[0])
        else:
            logger.error(f'Timestep {timestep} is not valid.')
         
        if extent is None:
            extent = [100000,430000,180000,480000]
         
        if ('International.Radar' in datasource) and (pd.to_datetime(start,format='%Y%m%d%H%M%S') < pd.to_datetime('2019-01-01')):
            logger.error('IRC is not available prior to 2019.')
            sys.exit()
               
        if variable == 'precipitation':
            varcode= 'P'
        elif variable == 'evaporation':
            varcode = 'Evaporation'                        
        else:
            varcode = variable
            
         # request template dictionary                                        
        req_body = {
            "Readers": [{
            "DataSourceCode": "",
            "Settings": {
              "StartDate": "",                  
              "EndDate": "",                  
              "VariableCodes": "",      
              "Interval": { "Type": "", "Value": 0},
              "Extent": {
                "Xll": "",
                "Yll": "",
                "Xur": "",
                "Yur": "",
                "SpatialReference": {"Epsg": 28992}
                }  
              } 
            }],
            "Exporter": { "DataFormatCode": "json" }
        }
   
        #%%
        # update the request
        req_body['Readers'][0].update({'DataSourceCode': datasource})
        req_body['Readers'][0]['Settings'].update({'StartDate': start})            
        req_body['Readers'][0]['Settings'].update({'EndDate': end})
        req_body['Readers'][0]['Settings'].update({'VariableCodes': [varcode]})
        req_body['Readers'][0]['Settings']['Interval'].update({'Type': t_unit, 'Value': t_val})
        req_body['Readers'][0]['Settings']['Extent'].update({'Xll': extent[0], 'Yll': extent[1], 'Xur' : extent[2], 'Yur': extent[3]})
        
        logger.info('Obtaining {datasource} from WIWB...')
        r = requests.post(url, json=req_body, auth=credentials)
        
        if r.status_code != 200:
            logger.error(f"No valid connection to WIWB (error {r.status_code}), aborting download.")
            sys.exit()
            
        resdict = json.loads(r.content)
        
        if '9' in resdict['Meta']["Projections"]:
            proj4 = '+proj=stere +lat_0=90 +lat_ts=60 +lon_0=0 +k=1 +x_0=0 +y_0=0 +a=6378140 +b=6356750'
            res = resdict['Data'][0]['GridDefinition']['CellWidth']*1e3
            xyll = (resdict['Data'][0]['GridDefinition']['Xll']*1e3,resdict['Data'][0]['GridDefinition']['Yll']*1e3)
        else:
            proj4 = resdict["Meta"]['Projections'][list(resdict['Meta']['Projections'].keys())[0]]['ProjectionString']
            res = resdict['Data'][0]['GridDefinition']['CellWidth']
            xyll = (resdict['Data'][0]['GridDefinition']['Xll'],resdict['Data'][0]['GridDefinition']['Yll'])
        ncol=resdict['Data'][0]['GridDefinition']['Columns']
        nrow=resdict['Data'][0]['GridDefinition']['Rows']    
         
        times = pd.date_range(start=pd.to_datetime(start,format='%Y%m%d%H%M%S'),end=pd.to_datetime(end,format='%Y%m%d%H%M%S'), freq=timestep)[1:]    
        dat_arr = np.zeros((len(times), nrow, ncol))*np.nan        
        for i in range(len(times)):
            dat_arr[i,:,:] = np.reshape(resdict['Data'][i]['Data'],(nrow,ncol))          
         
        logger.info('Converting to a dataset in RD-coordinates...')        
        xs = xyll[0] + np.arange(ncol)*res
        ys = xyll[1] + nrow*res -  np.arange(nrow)*res
        ds = xr.Dataset(
                        data_vars=dict(
                            variable=(["time","y", "x"], dat_arr),            
                            ),
                            coords=dict( time=times, x=xs, y=ys)
                            )                
        ds = ds.rio.write_crs(rio.crs.CRS.from_proj4(proj4))    
        #get a reference raster to match IRC with        
        reff = Path('\\'.join(os.path.abspath(__file__).split('\\')[:-1])) / 'precipitation_meteobase.asc'        
        refrast = imod.rasterio.open(reff)             
        refrast = refrast.rio.set_crs(28992, inplace=True)                
        # is not already) and in the process resample to 1000m resolution
        ds_rd = ds.rio.reproject_match(refrast)             
        #ds_rd = ds.rio.reproject(self.proj)
        self.write2ascii(ds_rd, variable, path=download_path, timestep=t_unit)
         
    