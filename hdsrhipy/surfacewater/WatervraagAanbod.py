# -*- coding: utf-8 -*-
"""
Created on Thu Dec 16 15:39:38 2021

@author: Michiel Pezij
mail: pezij@hkv.nl
"""

import os
import copy
import numpy as np
import pandas as pd
import xarray as xr
import geopandas as gpd
from tqdm.auto import tqdm


class WatervraagAanbod():
    
    def __init__(self, fNames=None):
        
        self.invoerdata = dict()
        self.metadata = dict()
        self.resultaten = dict()
        self.intersects = dict()
        self.data = dict()
        self.schaalgebieden = dict()
        self.klimaatscen = dict()
        self.klimaatreeks = dict()
        self.klimaatfactor = dict()
        
        self.invoerdata['top10nl'] = gpd.read_file(fNames['top10nl'])
        self.invoerdata['lsws'] = gpd.read_file(fNames['lsw'])
        self.invoerdata['afwateringseenheden'] = gpd.read_file(fNames['afwateringseenheden'])
        
        # dissolve top10nl
        self.invoerdata['top10nl'] = self.invoerdata['top10nl'].dissolve('typelandgebruik', as_index=False)
        
        # lees mozart data
        self.inlezen_mozart_data(fNames)
        
        # verkrijg kolomnamen
        self.verkrijg_kolomnamen()
        
        # merge koppelcodes top10nl
        self.merge_koppelcodes_top10nl()
        
        print('Data correct ingeladen.')

    def inlezen_mozart_data(self, fNames, schrijf_csv=False, use_dask=False):
        '''
        Functie voor het inlezen van de lswwaterbalans.out-file van Mozart.
        Door de bestandsgrootte wordt Dask hiervoor gebruikt.
        Een optie is om dit bestand weg te schrijven als een csv-bestand.
    
        '''
    
        if use_dask:
    
            import dask.dataframe as dd
    
            # inlezen van Mozart-gegevens
            df = dd.read_csv(fNames['mozart_out'],
                             skiprows=1,
                             delim_whitespace=True)
            
            self.invoerdata['mozart_uitvoer'] = df
    
            if schrijf_csv:
                df.to_csv(r'..\data\mozart\lswwaterbalans.csv',
                          index=False,
                          sep=',',
                          single_file=True)
    
        if use_dask == False:
            # inlezen van Mozart-gegevens
            df = pd.read_csv(fNames['mozart_out'],
                             skiprows=1,
                             delim_whitespace=True)
    
            if schrijf_csv:
                df.to_csv(r'..\data\mozart\lswwaterbalans.csv',
                          index=False,
                          sep=',')
            self.invoerdata['mozart_uitvoer'] = df
            
        self.invoerdata['koppelcodes'] = pd.read_excel(fNames['top10nl_koppeling'])
            
    def verkrijg_kolomnamen(self):
        '''
        functie om kolomnamen van vraag- en aanbodtermen uit dataframe te halen
        '''
        # columns = df_hdsr.columns
        demand_columns = [col for col in self.invoerdata['mozart_uitvoer'] if col.startswith('DEM')]
        alloca_columns = [col for col in self.invoerdata['mozart_uitvoer'] if col.startswith('ALLOC')]
        
        ''' 
        
        de watervraag per lsw bestaat uit:
            DEM_AGRIC: Gevraagd water voor beregening vanuit local surfacewater (m3) 
            DEM_WM: Gevraagd water voor peilbeheer vanuit local surfacewater (m3) 
            DEM_FLUSH: Gevraagd water voor doorspoeling vanuit local surfacewater (m3) 
            DEM_FLUSHRET: Gevraagd water voor doorspoeling (terug) vanuit local surfacewater (m3) 
            DEM_PUBWAT: Gevraagd water voor drinkwaterbereiding vanuit local surfacewater (m3) 
            DEM_INDUS: Gevraagd water voor industrie vanuit local surfacewater (m3)
            DEM_GRHOUSE: Gevraagd water voor glastuinbouw vanuit local surfacewater (m3) 
            DEM_WMTOTAL: Totaal gevraagd water voor peilbeheer, inclusief tekort op de waterbalans 
            DEM_WM_TODW: Gevraagd water voor peilbeheer, inclusief tekort op de waterbalans, 
            aan districtswater (hoofdsysteem)
        
        het wateraanbod per lsw bestaat uit:
            ALLOC_AGRIC: Gealloceerd water voor beregening vanuit local surfacewater (m3) 
            ALLOC_WM: Gealloceerd water voor peilbeheer vanuit local surfacewater (m3) 
            ALLOC_FLUSH: Gealloceerd water voor doorspoeling vanuit local surfacewater (m3) 
            ALLOC_FLUSHR: Gealloceerd water voor doorspoeling (terug) vanuit local surfacewater 
            ALLOC_PUBWAT: Gealloceerd water voor drinkwaterbereiding vanuit local surfacewater 
            ALLOC_INDUS: Gealloceerd water voor industrie vanuit local surfacewater (m3) 
            ALLOC_GRHOUS: Gealloceerd water voor glastuinbouw vanuit local surfacewater (m3)
            ALLOC_WM_DW: Gealloceerd water voor peilbeheer, inclusief tekort op de waterbalans, vanuit districtswater (hoofdsysteem)  
           
        '''        
        
        self.metadata['demand_columns'] = demand_columns
        self.metadata['alloca_columns'] = alloca_columns
        
    def merge_koppelcodes_top10nl(self):
        self.invoerdata['top10nl_gekoppeld'] = self.invoerdata['top10nl'].merge(self.invoerdata['koppelcodes'], on='typelandgebruik')
        
    def selecteer_lsw_data(self, lsw_nr):
        self.invoerdata['lsw'] = self.invoerdata['lsws'][self.invoerdata['lsws']['LSWFINAL']==lsw_nr]
        self.invoerdata['lsw_data'] = self.invoerdata['mozart_uitvoer'][self.invoerdata['mozart_uitvoer']['LSWNR'].isin([lsw_nr])]
        
        # check dat er maar naar 1 lsw wordt gekeken
        assert len(self.invoerdata['lsw']['LSWNR'].unique())==1
        
    def bepaal_intersects(self):
        
        # intersection tussen lsw en schaalgebied
        self.intersects['lsw_schaalgebied'] = self.invoerdata['lsw'].overlay(self.invoerdata['schaalgebied'], how='intersection')


        # intersection tussen lsw's in schaalgebied en top10nl
        self.intersects['lsw_schaalgebied_top10nl'] = self.intersects['lsw_schaalgebied'].overlay(self.invoerdata['top10nl_gekoppeld'], how='intersection')
        
        # intersection tussen lsw en top10nl       
        self.intersects['lsw_top10nl'] = self.invoerdata['lsw'].overlay(self.invoerdata['top10nl_gekoppeld'], how='intersection')
        
        # dissolve intersection tussen lsw en top10nl
        self.intersects['lsw_top10nl_dissolved'] = self.intersects['lsw_top10nl'].dissolve('mozart_koppeling', as_index=False)
        
    def mozart_to_df(self):
        
         # verkrijg start- en einddatum van decades
        t_start = self.invoerdata['lsw_data']['TIMESTART'].astype(int)
        # t_end = self.invoerdata['lsw_data']['TIMEEND'].astype(int)
    
        dt_start = pd.to_datetime(t_start, format='%Y%m%d')
        # dt_end = pd.to_datetime(t_end, format='%Y%m%d')
    
        # selecteer de vraag- en aanbodtermen en maak een dataframe aan
        df_demand = pd.DataFrame(index=dt_start,
                                 data=self.invoerdata['lsw_data'][self.metadata['demand_columns']].values,
                                 columns=self.metadata['demand_columns'])
        df_alloca = pd.DataFrame(index=dt_start,
                                 data=self.invoerdata['lsw_data'][self.metadata['alloca_columns']].values,
                                 columns=self.metadata['alloca_columns'])
    
        # De jaren 2009 en 2010 worden als inspeeljaren gezien bij de validatie
        # van het LHM 4.1 en moeten niet gebruikt worden tijden de analyse!
        # haal ook 2020 en 2021 uit de uitvoer
        df_demand = df_demand.loc['2011':'2020-01-01']
        df_alloca = df_alloca.loc['2011':'2020-01-01']
            
        self.data['watervraag'] = df_demand
        self.data['wateraanbod'] = df_alloca

    def bereken_oppervlak_landgebruik(self):
        self.data['opp_lsw'] = pd.DataFrame(index=self.intersects['lsw_top10nl_dissolved']['mozart_koppeling'],
                                            data=self.intersects['lsw_top10nl_dissolved'].area.values)
        
        self.data['opp_perc'] = self.data['opp_lsw']/self.data['opp_lsw'].sum()
        
    def verkrijg_topnl_schaalgebied(self, code):
        
        self.data['top10nl_schaalgebied'] = self.intersects['lsw_schaalgebied_top10nl'][self.intersects['lsw_schaalgebied_top10nl']['CODE']==code]
        
        self.data['top10nl_peilgebied_dissolved'] = self.data['top10nl_schaalgebied'].dissolve('mozart_koppeling', as_index=False)
        
    def bereken_oppervlak_landgebruik_schaalgebied(self, onzekerheidsfactor):
        self.data['opp_schaalgebied'] = pd.DataFrame(index=self.data['top10nl_peilgebied_dissolved']['mozart_koppeling'],
                                                     data=self.data['top10nl_peilgebied_dissolved'].area.values)
   
        # self.data['opp_schaalgebied'] = self.data['opp_schaalgebied']
        self.data['opp_schaalgebied_perc'] = self.data['opp_schaalgebied']*onzekerheidsfactor/self.data['opp_schaalgebied'].sum()
        
        # zet NaN-waarden om naar 0%
        self.data['opp_schaalgebied_perc'].fillna(0, inplace=True)
        
    def verdeel_watervraag_aanbod(self):
        
        
        self.resultaten['watervraag_aanbod'] = pd.DataFrame(index=self.data['watervraag'].index)
        
        try:
            self.resultaten['watervraag_aanbod']['watervraag_landbouw'] = self.data['watervraag']['DEM_AGRIC']*self.data['opp_schaalgebied_perc'].loc['landbouw'].item()*-1
            self.resultaten['watervraag_aanbod']['wateraanbod_landbouw'] = self.data['wateraanbod']['ALLOC_AGRIC']*self.data['opp_schaalgebied_perc'].loc['landbouw'].item()*-1
            self.resultaten['watervraag_aanbod']['watertekort_landbouw'] = self.resultaten['watervraag_aanbod']['watervraag_landbouw']-self.resultaten['watervraag_aanbod']['wateraanbod_landbouw']
        
        except KeyError:
            self.resultaten['watervraag_aanbod']['watervraag_landbouw'] = 0
            self.resultaten['watervraag_aanbod']['wateraanbod_landbouw'] = 0
            self.resultaten['watervraag_aanbod']['watertekort_landbouw'] = 0
            
        try:
            self.resultaten['watervraag_aanbod']['watervraag_water'] = self.data['watervraag']['DEM_WM']*self.data['opp_schaalgebied_perc'].loc['water'].item()*-1
            self.resultaten['watervraag_aanbod']['wateraanbod_water'] = self.data['wateraanbod']['ALLOC_WM']*self.data['opp_schaalgebied_perc'].loc['water'].item()*-1
            self.resultaten['watervraag_aanbod']['watertekort_water'] = self.resultaten['watervraag_aanbod']['watervraag_landbouw']-self.resultaten['watervraag_aanbod']['wateraanbod_landbouw']
        
        except KeyError:
            self.resultaten['watervraag_aanbod']['watervraag_water'] = 0
            self.resultaten['watervraag_aanbod']['wateraanbod_water'] = 0
            self.resultaten['watervraag_aanbod']['watertekort_water'] = 0
            
    def schaling_naar_csv(self, directory):
        
        for key, data_dict in self.schaalgebieden.items():
            data_dict.to_csv(f'{directory}\mozart_schaalgebied_{key}_geschaald.csv')
            

    def uitvoeren_schaling(self, lsw_nr, onzekerheid_opp):
        
        # selecteer gegevens van specifieke lsw
        self.selecteer_lsw_data(lsw_nr)
        
        # bepaal de intersects tussen de invoergegevens
        self.bepaal_intersects()

        # verkrijg watervraag en wateraanbod voor de lsw
        self.mozart_to_df()

        # bereken oppervlakte van de landgebruiktypen voor de lsw
        self.bereken_oppervlak_landgebruik()

        # verkrijg unieke codes van schaalgebieden in lsw
        self.metadata['lsw_schaalgebied_top10nl'] = self.intersects['lsw_schaalgebied_top10nl']['CODE'].unique()

        # loop over elk schaalgebied in de lsw:
        for code in self.metadata['lsw_schaalgebied_top10nl']:

            # verkrijg top10nl voor schaalgebied
            self.verkrijg_topnl_schaalgebied(code)
            
            # verkrijg oppervlakten en percentages van schaalgebied
            self.bereken_oppervlak_landgebruik_schaalgebied(1+onzekerheid_opp)
            
            # verdeel de watervraag, -aanbod en -tekorten voor de functies indien deze bestaan in schaalgebied
            self.verdeel_watervraag_aanbod()
            
            if code in self.schaalgebieden.keys():
                self.schaalgebieden[code] = self.schaalgebieden[code] + copy.deepcopy(self.resultaten['watervraag_aanbod'])
                
            else:
                self.schaalgebieden[code] = copy.deepcopy(self.resultaten['watervraag_aanbod'])

    # # selecteer gegevens van specifieke lsw
    # test.selecteer_lsw_data(lsw_nr)
    
    # # bepaal de intersects tussen de invoergegevens
    # test.bepaal_intersects()
    
    # # verkrijg watervraag en wateraanbod voor de lsw
    # test.mozart_to_df()
    
    # # bereken oppervlakte van de landgebruiktypen voor de lsw
    # test.bereken_oppervlak_landgebruik()
    
    # # verkrijg unieke codes van schaalgebieden in lsw
    # test.metadata['lsw_schaalgebied_top10nl'] = test.intersects['lsw_schaalgebied_top10nl']['CODE'].unique()

        # # loop over elk schaalgebied in de lsw:
        # for code in test.metadata['lsw_schaalgebied_top10nl'][:]:
            
        #     # verkrijg top10nl voor schaalgebied
        #     test.verkrijg_topnl_schaalgebied(code)
            
        #     # verkrijg oppervlakten en percentages van schaalgebied
        #     test.bereken_oppervlak_landgebruik_schaalgebied(1+onzekerheid_opp)
            
        #     # verdeel de watervraag, -aanbod en -tekorten voor de functies indien deze bestaan in schaalgebied
        #     test.verdeel_watervraag_aanbod()
            
        #     if code in test.schaalgebieden.keys():
        #         test.schaalgebieden[code] = test.schaalgebieden[code] + copy.deepcopy(test.resultaten['watervraag_aanbod'])
                
        #     else:
        #         test.schaalgebieden[code] = copy.deepcopy(test.resultaten['watervraag_aanbod'])

    def mozart_watervraag_wateraanbod(self):
        for key, df_invoer in self.klimaatscen.items():
            df_balans = pd.DataFrame()
            
            # bepaal watervraag per tijdstap (decade)
            df_balans['watervraag'] = df_invoer['demand_agric'] + \
                        df_invoer['demand_flush'] + \
                        df_invoer['demand_pubwat'] + \
                        df_invoer['demand_industry'] + \
                        df_invoer['demand_greenhouse'] + \
                        df_invoer['demand_wmtot']*-1
                        # df_demand['DEM_FLUSHRET'] + \
        
            # bepaald wateraanbod per tijdstap (decade)
            df_balans['wateraanbod'] = df_invoer['alloc_agric'] + \
                        df_invoer['alloc_wm'] + \
                        df_invoer['alloc_flush'] + \
                        df_invoer['alloc_pubwat'] + \
                        df_invoer['alloc_industry'] + \
                        df_invoer['alloc_greenhouse'] #+ \  
        
            self.klimaatreeks[key] = copy.deepcopy(df_balans)
    
    def bepaal_factor_klimaatscenarios(self, directory, scenarios):
        
        for scenario in scenarios:
    
            # definieer map waar de run-resultaten staan
            dir_bp = fr'{directory}\bp2018_LHM_{scenario}'
    
            # verkrijg alle bestandsnamen van de MozartRegio17-netcdfs.
            file_list = []
            for (dirpath, dirnames, filenames) in os.walk(dir_bp):
                for filename in filenames:
                    if filename==f'ZW_LHM_MozartRegio17_{scenario}BP18.nc':
                        file_list.append(os.sep.join([dirpath, filename]))
    
            # selecteer regio
            regio = b'Midden West Nederland - niet extern verzilt'

            # initaliseer lijst voor DataFrames
            list_dfs = []
        
            for file in tqdm(file_list[:]):
                ds = xr.open_dataset(file)
                
                index = np.where(ds.station_names==regio)[0].item()
                
                ds_hdsr = ds[dict(stations=index)].squeeze()
                df_hdsr = ds_hdsr.to_dataframe()
                
                
                df_hdsr.drop(['analysis_time', 'lat', 'lon', 'y', 'x', 'z', 'station_id', 'station_names',
                              'precip','evaporation','drainage_sh','drainage_dp','infiltration_sh',
                              'infiltration_dp','urbanrunoff','upstream','downstream','from_dw',
                              'to_dw','dstorage','balancecheck'],
                             axis=1,
                             inplace=True)
                
                list_dfs.append(copy.deepcopy(df_hdsr))
                
            full_df = pd.concat(list_dfs)
            
            full_df.to_csv(f'volledige_mozart_reeks_HDSR_{scenario}BP18.csv')
            
            self.klimaatscen[scenario] = copy.deepcopy(full_df)
            
            # for key, df in self.klimaatscen.items():
            #     self.klimaatreeks[key] = self.mozart_watervraag_wateraanbod(df)
            # self.klimaatreeks[key] = 
            self.mozart_watervraag_wateraanbod()
            
        self.klimaatfactor['watervraag'] = self.klimaatreeks[scenarios[1]]['watervraag']/self.klimaatreeks[scenarios[0]]['watervraag']
                                                             
        self.klimaatfactor['wateraanbod'] =self.klimaatreeks[scenarios[1]]['wateraanbod']/self.klimaatreeks[scenarios[0]]['wateraanbod']
          
        self.klimaatfactor['watervraag'] = self.klimaatfactor['watervraag'].mean()                                                 
        self.klimaatfactor['wateraanbod'] = self.klimaatfactor['wateraanbod'].mean()                                                 
                                                             
        print(f"Factor watervraag is: {self.klimaatfactor['watervraag']:.2f} [-]")
        print(f"Factor wateraanbod is: {self.klimaatfactor['watervraag']:.2f} [-]")
                                                                           