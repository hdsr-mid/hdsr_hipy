# -*- coding: utf-8 -*-
"""
Created on Mon Nov  8 17:08:51 2021

@author: hurkmans
"""
import os
import sys
package_path = os.path.abspath('../')
sys.path.append(package_path)  
import imod
from hdsrhipy import Runfile
from hdsrhipy.model import rasters
import hdsrhipy
import pandas as pd
import geopandas as gpd
from hdsrhipy.model import util
import shutil
import numpy as np
from datetime import datetime as dt
from datetime import timedelta as td
from tqdm import tqdm

def add_simgro_surface_water(rf, gdf_pg=None, run1period=False, change_conductance=False, datadir=None):
    """Add surface water (peilgebieden) from Simgro to an iMod-simulation"""
        
    print('Running one-timestep simulation to retrieve initial inp-files')
    # first run only the frist period
    rf2 = copy.deepcopy(rf)
    # only run the first timestep
    rf2.change_period(rf2.data['DATE'][0],rf2.data['DATE'][0])
    naam='run1period'
    model_ws = os.path.join(datadir,'work',naam)
    from hdsrhipy.model import hdsr
    if run1period:
        rf2.run_imodflow(model_ws,naam,data_path=datadir)

    # if gdf_pg is None:
    #     # download the layer with peilgebieden
    #     gdf_pg = hdsr.get_peilgebieden()
    
    # generate a flopy model to grid the peilgebieden
    ml = to_flopy(rf2, packages=['dis','lpf'], data_dir=datadir)
    # determine in each cell which CODENR has the largest area
    grid_sw = polygons2grid(gdf_pg, ml, 'CODENR', method='most_common')
    # make a surface_water unit for each remaining peilgebied
    gdf_sw = gdf_pg.set_index('CODENR')
    # only keep surface water units that are inside the model domain
    gdf_sw = gdf_sw.loc[np.unique(grid_sw[~np.isnan(grid_sw)]).astype(int)]
    # calculate area per cell
    cell_area = rf2.data['CSIZE'] * rf2.data['CSIZE']
    for swnr in gdf_sw.index:
        mask = grid_sw == swnr
        # determine the area in the model of each surface water unit
        gdf_sw.at[swnr,'area_grid'] = np.sum(mask)*cell_area
    
    # mana_res.inp
    mana_res = pd.DataFrame({'swnr':gdf_sw.index,'ioma':2,'lvtasm':0.0,'lvtawt':0.0,'dptasu':0.01,
                             'fxsuswsb':0.00346*gdf_sw['area_grid']})
    rf.add_inp(mana_res,'simgro\\mana_res.inp')
    
    # goto_res.inp
    goto_res = pd.DataFrame({'swnr':gdf_sw.index, 'swnrgo':0, 'iofwbk':1}).set_index('swnr',drop=False)
    for swnr in gdf_sw.index:
        if gdf_sw.at[swnr,'SOORTSTREE'] == '1':
            # vastpeil
            goto_res.at[swnr,'lvwrsm'] = gdf_sw.at[swnr,'VASTPEIL']
            goto_res.at[swnr,'lvwrwt'] = gdf_sw.at[swnr,'VASTPEIL']
        elif gdf_sw.at[swnr,'SOORTSTREE'] == '2':
            # flexibelpeil
            goto_res.at[swnr,'lvwrsm'] = gdf_sw.at[swnr,'BOVENPEIL']
            goto_res.at[swnr,'lvwrwt'] = gdf_sw.at[swnr,'BOVENPEIL']
        elif gdf_sw.at[swnr,'SOORTSTREE'] =='3':
            # zp/wp
            goto_res.at[swnr,'lvwrsm'] = gdf_sw.at[swnr,'ZOMERPEIL']
            goto_res.at[swnr,'lvwrwt'] = gdf_sw.at[swnr,'WINTERPEIL']
        else:
            raise(ValueError('Unknown SOORTSTREE: {}'.format(gdf_sw.at[swnr,'SOORTSTREE'])))
        if goto_res.at[swnr,'lvwrsm'] == -999:
            goto_res.at[swnr,'lvwrsm'] = 0.
        if goto_res.at[swnr,'lvwrwt'] == -999:
            goto_res.at[swnr,'lvwrwt'] = 0.
        goto_res.at[swnr,'lvwrlw'] = min(goto_res.at[swnr,'lvwrsm'],goto_res.at[swnr,'lvwrwt'])-0.01
    rf.add_inp(goto_res,'simgro\\goto_res.inp')
    
    # dish_res.inp
    # fswr is 1.5 l/s/ha = 1.5/1000/10000 m3/s/m2
    dish_res = pd.DataFrame({'swnr':gdf_sw.index,'dhwr':0.1,'fswr':1.5/1000/10000*gdf_sw['area_grid'],'swnrgo':0})
    dish_res0 = pd.DataFrame({'swnr':gdf_sw.index,'dhwr':0.0,'fswr':0.0,'swnrgo':0})
    dish_res = dish_res0.append(dish_res, ignore_index=True).sort_values(['swnr','dhwr']).reset_index(drop=True)
    rf.add_inp(dish_res,'simgro\\dish_res.inp')
    
    # swnr_sim.inp
    swnr_sim = pd.DataFrame({'swnr':gdf_sw.index})
    rf.add_inp(swnr_sim,'simgro\\swnr_sim.inp')
    
    # svat2swnr_roff.inp
    # read the existing svat2swnr_roff
    fname = os.path.join(model_ws,rf2.data['OUTPUTDIRECTORY'],'metaswap','svat2swnr_roff.inp')
    svat2swnr_roff = inp.read(fname)
    # read mod2svat
    fname = os.path.join(model_ws,rf2.data['OUTPUTDIRECTORY'],'metaswap','mod2svat.inp')
    svat2mod = inp.read(fname).set_index('svat')

    # get the modflow cell
    modfcell = svat2mod.loc[svat2swnr_roff['svat'],'modfcell']
    # calculate which swnr
    fname = os.path.join(model_ws,rf2.data['OUTPUTDIRECTORY'],'mf2005_tmp',
                     '{}.dxc'.format(os.path.split(model_ws)[-1]))
    lrc = pd.DataFrame(read_dxc(fname), columns=['l','r','c','modfcell'])
    lrc = lrc.set_index('modfcell')
    # from one-based to zero-based
    lrc = lrc-1
    assert np.all(lrc['l']==0)
    swnr = grid_sw[lrc.loc[modfcell,'r'], lrc.loc[modfcell,'c']]
    mask = ~np.isnan(swnr)
    svat2swnr_roff.loc[mask,'swnr'] = swnr[mask].astype(int)
    
    rf.add_inp(svat2swnr_roff,'simgro\\svat2swnr_roff.inp')
        
    # read area_svat.inp
    fname = os.path.join(model_ws,rf2.data['OUTPUTDIRECTORY'],'metaswap','area_svat.inp')
    area_svat = inp.read(fname).set_index('svat',drop=False)
    # generatesvat2swnr_drng.inp
    svat2swnr_drng = pd.DataFrame()
    
    if change_conductance:
        # hoe komen we tot een breedte?
        #wl = hdsr.get_waterlijn(datadir=datadir)
        wl = hdsr.get_waterlijn_met_breedte(datadir=datadir)
        wlm = geodataframe2grid(ml, wl)
        wlm.loc[wlm['IWS_W_WATD']==-99,'IWS_W_WATD'] = np.NaN
        # calculate length again
        wlm['LENGTH'] = wlm.length
        for cat in [1,2,3]:
            mask = wlm['CATEGORIEO'] == cat
            if not mask.any():
                continue
            wl_per_cel = wlm[mask].groupby(['row','col'])
            wl_sum = wl_per_cel.sum()
            wl_mean = wl_per_cel.mean()
            row, col = zip(*wl_per_cel.groups.keys())
            A = rf.data['CSIZE']*rf.data['CSIZE'] # celoppervlak (m2)
            H0 = ml.dis.thickness[0, row, col] # doorstroomde dikte (m)
            kv = ml.lpf.hk[0,row,col]# verticale doorlotendheid (m/d)
            kh = ml.lpf.vka[0,row,col]# horizontale doorlatendheid (m/d)
            c1 = ml.dis.thickness[1,row,col] / ml.lpf.vkcb[0,row,col]# deklaagweerstand (d)
            li = wl_sum['LENGTH'] # lengte van de waterlopen (m)
            Bin = wl_mean['WIDTH'] # bodembreedte (m)
            c0 = 1 # slootbodemweerstand (d)
            #cond = util.lek_conductance(A, H0, kv, kh, c1, li, Bin, c0)
            cond = lek_conductance_array(A, H0, kv, kh, c1, li, Bin, c0)
            wd = wl_mean['IWS_W_WATD'] # waterdiepte
            wd[wd.isna()] = 0
            svat2swnr_drng = cond2svat2swnr_drng(svat2swnr_drng, cond, wd, li, Bin, area_svat,
                                                 svat2mod, grid_sw, goto_res, lrc, rf, ml, sy=cat)
        
    # add rivers
    for lay in [1]:# only the first layer
        for lev in [1,2,3]:
            cond_file = 'OPPERVLAKTEWATER\WINTER\CONDUCTANCE_LAAG{}_{}.IDF'.format(lay,lev)
            botm_file = 'OPPERVLAKTEWATER\WINTER\BODEMHOOGTE_LAAG{}_{}.IDF'.format(lay,lev)
            infi_file = 'OPPERVLAKTEWATER\WINTER\INFFACTOR_LAAG{}_{}.IDF'.format(lay,lev)
            svat2swnr_drng = idf2svat2swnr_drng(svat2swnr_drng,model_ws,cond_file,
                                                botm_file,area_svat,svat2mod,grid_sw,
                                                lrc,rf,ml,'RIV',infi_file=infi_file,sy=lev,
                                                change_conductance=change_conductance, datadir=datadir)
    
    # add drains from 
    cond_file = 'DRAINAGE\CONDUCTANCE_NIVEAU1.IDF'
    botm_file = 'DRAINAGE\DRAINAGE_NIVEAU1.IDF'
    svat2swnr_drng = idf2svat2swnr_drng(svat2swnr_drng,model_ws,cond_file,
                                        botm_file,area_svat,svat2mod,grid_sw,
                                        lrc,rf,ml,'DRN',sy=4,
                                        change_conductance=change_conductance, datadir=datadir)
    rf.add_inp(svat2swnr_drng,'simgro\\svat2swnr_drng.inp')


import geopandas as gpd
resolution = 25.         
hydromedah_path = r'E:\Hydromedah'      
start_date='2010-01-01'
end_date='2010-12-31'
rf = Runfile(os.path.join(hydromedah_path,'hdsr_ns.run'), data_path='$DBASE$\\', evap_path=hydromedah_path, precip_path=hydromedah_path)
rf.to_version(4.3)
rf.update_metaswap(datadir=hydromedah_path, start_date=pd.to_datetime(start_date), end_date=pd.to_datetime(end_date), metaswap_vars=metaswap_vars)


gdf_pg = gpd.read_file(os.path.join(hydromedah_path, 'Afwateringseenheden.shp'))
gdf_pg['CODENR'] = gdf_pg.index + 1                  
add_simgro_surface_water(rf, gdf_pg=gdf_pg, run1period=False, datadir=hydromedah_path)