# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 09:17:22 2021

@author: HKV lijn in water 2021
Contact: Ruud Hurkmans (hurkmans@hkv.nl)
"""

from pathlib import Path
import hkvsobekpy as hkv
import pandas as pd
import geopandas as gpd

class FlowStats:    
    
    def __init__(self, sobek_path=None, template_path=None, export_path=None):
       
        if sobek_path is not None:
            self.sobek_path = Path(sobek_path)
        
        if export_path is None:                   
            self.export_path = Path(export_path) / 'debietstatistiek' / 'debietstatistiek.shp'
            self.export_path.mkdir(parents=True, exist_ok=True)           
        else:
            self.export_path = export_path

        if template_path is not None:                   
            self.template_path = Path(template_path)             

    # functie om statistieken toe te voegen aan de shapefile
    def overschrijdingsduur(self,x=None, threshold=None):
        return(len(x[x > threshold]))
    def relatieve_overschrijdingsduur(self, x=None, threshold=None):
        return(len(x[x > threshold])/len(x))
    
    def add_statistic(self, output, data, stype,overschrijdingsduur_van_debiet=None, relatieve_overschrijdingsduur_van_debiet=None, debiet_bij_tijdpercentiel=None, name=None):
        if stype=='gemiddelde':
            statcol = data.mean(axis=0)
            statcol.name = name
            output = pd.concat([output, statcol],axis=1, join="inner",ignore_index=False,copy=True)        
        elif stype=='standaard_afwijking':
            statcol = data.std(axis=0)
            statcol.name = name
            output = pd.concat([output, statcol],axis=1, join="inner",ignore_index=False,copy=True)        
        elif stype=='variantie':
            statcol = data.std(axis=0)**2
            statcol.name = name
            output = pd.concat([output, statcol],axis=1, join="inner",ignore_index=False,copy=True)        
        elif stype=='overschrijdingsduur':
            statcol = data.apply(self.overschrijdingsduur, axis=0, threshold=overschrijdingsduur_van_debiet)
            statcol.name = name
            output = pd.concat([output, statcol],axis=1, join="inner",ignore_index=False,copy=True)        
        elif stype=='relatieve_overschrijdingsduur':
            statcol = data.apply(self.relatieve_overschrijdingsduur, axis=0, threshold=relatieve_overschrijdingsduur_van_debiet)
            statcol.name = name
            output = pd.concat([output, statcol],axis=1, join="inner",ignore_index=False,copy=True)        
        elif stype=='percentiel':
            statcol = data.quantile(debiet_bij_tijdpercentiel)
            statcol.name=name
            output = pd.concat([output, statcol],axis=1, join="inner",ignore_index=False,copy=True)        
        else:
            print('Onbekend statistiektype. Kies uit: gemiddelde, standaard_afwijking, variantie, overschrijdingsduur, relatieve_overschrijdingsduur, percentiel')
        return(output)

    def get_flowstats(self, cases2join=None, period=None, overschrijdingsduur_van_debiet=None, relatieve_overschrijdingsduur_van_debiet=None, debiet_bij_tijdpercentiel=None):
        """
        Function to calculate flow stastistics 
        """    
        seg_shape = gpd.read_file(str(self.template_path / 'RchSegments.shp'))
        seg_shape.columns = [col.strip() for col in seg_shape.columns]
        seg_shape.index = seg_shape['ID']
        seg_shape.drop(['NAME', "TYPE", 'PARENTID', 'ID_FROM', 'ID_TO', 'USERID'], axis=1, inplace=True)

        subhislist = []
        for case in cases2join:             
            pad_to_his = str(self.sobek_path / str(case) / 'reachseg.his')
        #     pad_to_his = str(pad_to_sobek / 'reachseg.his')
            reachseg = hkv.read_his.ReadMetadata(pad_to_his)    
            # struc = hkv.read_his.ReadMetadata(str(pad_to_his / 'struc.his'))
            reach_segments = reachseg.DataFrame()
            reach_segments = reach_segments.iloc[32:,:]
            locs= reachseg.GetLocations()
            params= reachseg.GetParameters()
            reach_segments = reach_segments.loc[:,params[params.index('Discharge mean(m³/s)')]]
            reach_segments.head()
            subhislist.append(reach_segments)
        all_data = pd.concat(subhislist)

        sub_data = all_data[all_data.index.strftime('%m%d').isin(period.strftime('%m%d'))]
        output = self.add_statistic(seg_shape, sub_data, 'gemiddelde', name='Av')
        output = self.add_statistic(output, sub_data, 'variantie', name='Var')
        if overschrijdingsduur_van_debiet is not None:
            output = self.add_statistic(output, sub_data, 'overschrijdingsduur', name=overschrijdingsduur_van_debiet['name'], overschrijdingsduur_van_debiet=overschrijdingsduur_van_debiet['value'])
        if relatieve_overschrijdingsduur_van_debiet is not None:
            output = self.add_statistic(output, sub_data, 'relatieve_overschrijdingsduur', name=relatieve_overschrijdingsduur_van_debiet['name'], relatieve_overschrijdingsduur_van_debiet=relatieve_overschrijdingsduur_van_debiet['value'])
        if debiet_bij_tijdpercentiel is not None:
            output = self.add_statistic(output, sub_data, 'percentiel', name=debiet_bij_tijdpercentiel['name'], debiet_bij_tijdpercentiel=debiet_bij_tijdpercentiel['value'])

        output.to_file(str(self.export_path))