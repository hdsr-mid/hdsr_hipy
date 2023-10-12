# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 15:08:08 2018

@author: Artesia
"""
#from owslib.wfs import WebFeatureService
import os
import pandas as pd
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
from shapely.strtree import STRtree
from tqdm import tqdm


def get_peilgebieden(fname=None, wfs=None):
    """Get the peilgebieden as a Polygon-GeoDataFrame"""
    if fname is None:
        fname=r'data/downloads/BRPeilgebieden.json'
    # if not os.path.exists(fname):
    #     if wfs is None:
    #         wfs = get_wfs()
    #     response = wfs.getfeature(typename='Datadeler:BRPeilgebieden',outputFormat='json')
    #     out = open(fname, 'wb')
    #     out.write(bytes(response.read(), 'UTF-8'))
    #     out.close()
    # read this data to a geodataframe
    gdf_pg = gpd.read_file(fname)
    # set 'CODE' as the index of the geodataframe
    gdf_pg.set_index('CODE',inplace=True)
    # remove the river Lek
    gdf_pg = gdf_pg.loc[gdf_pg['NAAM']!='Nederrijn en Lek']
    # remove the Amsterdam-Rijnkanaal en Lekkanaal
    gdf_pg = gdf_pg.loc[gdf_pg['NAAM']!='Amsterdam-Rijnkanaal en Lekkanaal']
    # verwijder alle vrije peilen (streefpeil 4, 5 en 99)
    # (hou alleen soort streefpeil 1, 2 of 3 over)
    mask = (gdf_pg.SOORTSTREE==1) | (gdf_pg.SOORTSTREE==2) | (gdf_pg.SOORTSTREE==3)
    gdf_pg = gdf_pg.loc[mask]
    # determine the codenr as an integer
    CODENR = [int(x[2:]) for x in gdf_pg.index]
    # deze moet uniek zijn:
    assert len(CODENR) == len(np.unique(CODENR))
    gdf_pg['CODENR'] = CODENR
    return gdf_pg
    
def get_watervlak(fname=None, wfs=None, datadir=None):
    """Get all the surface water as a Polygon-GeoDataFrame"""
    if fname is None:
        fname = os.path.join(datadir,'downloads/KEUR_2018:Kaart_1A_Watervlak.json').replace(':','_')
    # if not os.path.exists(fname):
    #     if wfs is None:
    #         wfs = get_wfs()
    #     response = wfs.getfeature(typename='KEUR_2018:Kaart_1A_Watervlak',outputFormat='json')
    #     out = open(fname, 'wb')
    #     out.write(bytes(response.read(), 'UTF-8'))
    #     out.close()
    gdf_wv = gpd.read_file(fname)
    return gdf_wv.set_index('id')

def get_waterlijn(fname=None, wfs=None, datadir=None):
    """Get all the surface water as a LineString-GeoDataFrame"""
    # if fname is None:
    #     fname = os.path.join(datadir,'downloads/KEUR_2018:Kaart_1A_Leggervak.json')
    # if not os.path.exists(fname):
    #     if wfs is None:
    #         wfs = get_wfs()
    #     response = wfs.getfeature(typename='KEUR_2018:Kaart_1A_Leggervak',outputFormat='json')
    #     out = open(fname, 'wb')
    #     out.write(bytes(response.read(), 'UTF-8'))
    #     out.close()
    gdf_wl = gpd.read_file(fname)
    return gdf_wl.set_index('id')

def get_waterlijn_met_breedte(fname=None, wfs=None, wl=None, wv=None, datadir=None):
    if fname is None:
        fname = os.path.join(datadir,'downloads/KEUR_2018:Kaart_1A_Leggervak_Breedte.json')
    if os.path.exists(fname):
        wl = wl = gpd.read_file(fname)
    else:
        # stap 1: lees de basisdata
        if wl is None:
            wl = get_waterlijn(wfs=wfs, datadir=datadir)
        if wv is None:
            wv = get_watervlak(wfs=wfs, datadir=datadir)
        # stap 2: versnij de lijnen met de vlakken
        if False:
            # too slow and only introduced in geopandas 0.7
            wl = gpd.clip(wl, wv)
            # bepaal welk polygon bij welke lijn hoort
            wv['AREA'] = wv.area
            wl['LENGTH'] = wl.length
            wl = gpd.sjoin(wl, wv)
        else:
            geometry = wv['geometry']
            for id in wv.index:
                geometry[id].id = id
            s = STRtree(geometry)
            shp_list = []
            for index in tqdm(wl.index):
                iwl = wl.loc[index]
                iwvs = s.query(iwl.geometry)
                for iwv in iwvs:
                    i = iwl.geometry.intersection(iwv)
                    add_wl_to_list(i, iwl, iwv, shp_list)
            wl = gpd.GeoDataFrame(pd.DataFrame(shp_list))
        # stap 4: verwijder de lijnen die niet de langste zijn in een polygoon
        wl['LENGTH'] = wl.length
        wl = wl.sort_values('LENGTH', ascending=False).drop_duplicates(['id_wv'])
        # stap 3: bereken de breedte door het oppervlak van het polygon te delen door de lengte van de lijn
        # oftewel: ken het oppervlak van het polygoon toe aan de langste lijn in het polygon
        wl['AREA'] = wv.area[wl['id_wv']].values
        wl['WIDTH'] = wl['AREA'] / wl['LENGTH']
        # save the result
        wl.to_file(fname, driver='GeoJSON')
    return wl

def add_wl_to_list(i,wl,wv,shp_list):
    """subfunction of get_waterlijn_met_breedte"""
    if not i.is_empty:
        it = i.geometryType()
        if it == 'GeometryCollection':
            for im in i.geoms:
                add_wl_to_list(im, wl, wv, shp_list)
        elif it == 'MultiLineString':
            for im in i.geoms:
                add_wl_to_list(im, wl, wv, shp_list)
        elif it == 'LineString':
            wln = wl.copy()
            wln['geometry'] = i
            wln['id_wv'] = wv.id
            shp_list.append(wln)
        elif it == 'Point':
            # endpoint of the linestring is on the cell-edge
            pass
        elif it == 'MultiPoint':
            # mutiple endpoints of the linestring are on the cell-edge
            pass
        else:
            raise NotImplementedError(
                'geometryType ' + it + ' not yet supprted in add_wl_to_list')
    
#def get_wfs(url='https://geoserver-hdsr.webgispublisher.nl/wfs', version='2.0.0'):
#    """Get an instance of a WebFeatureService (WFS) of the hdsr wfs-server"""
#    return WebFeatureService(url=url, version=version)

def plot_peilgebieden(gdf_pg=None,newfig=True,labels=False,**kwargs):
    """Plot the peilgebieden shape"""
    if newfig:
        f,ax=plt.subplots()
        ax.axis('equal')
    else:
        ax = plt.gca()
        f = ax.figure
    if gdf_pg is None:
        gdf_pg = get_peilgebieden()
    gdf_pg.plot(ax=ax,picker=True,**kwargs)
    if labels:
        for index,row in gdf_pg.iterrows():
            ax.annotate(s=index, xy=row.geometry.centroid.coords[0], ha='center')
    def onpick(event):
        for i in event.ind:
            print(gdf_pg.index[i])
    f.canvas.mpl_connect('pick_event', onpick)
    return f,ax