# -*- coding: utf-8 -*-
"""
Created on Mon Sep  3 12:04:09 2018

@author: ruben
"""
import os
import shutil
from shapely.geometry import Point, Polygon, MultiLineString, MultiPolygon
from shapely.strtree import STRtree
import geopandas as gpd
import pandas as pd
from osgeo import gdal
import hdsrhipy
import numpy as np
import flopy
import zipfile
import time
import imod
import xarray
import copy
from matplotlib.path import Path
import scipy.spatial.qhull as qhull
from hdsrhipy.model import rasters
from hdsrhipy.model import inp
from hdsrhipy.model import hdsr
import tempfile
import rasterio.warp
import warnings

def compress_idfs(pathname,remove=True):
    """Compress all idf-files in the folder pathname to zip-files"""
    if os.path.isdir(pathname):
        idf_files = []
        for root, dirs, files in os.walk(pathname):
            for file in files:
                if file.lower().endswith(".idf"):
                     idf_files.append(os.path.join(root, file))
    elif isinstance(pathname,list):
        idf_files = pathname
    else:
        idf_files=[pathname]
    for file in idf_files:
        with zipfile.ZipFile(file + '.zip', 'w', zipfile.ZIP_DEFLATED) as z:
            z.write(file,os.path.basename(file))
        if remove:
            os.remove(file)

def copy_data_file(src_dir,dst_dir,file):
    """Copy a data-file, and unzip the file if needed"""
    if is_number(file):
        # The file can also be one single number
        return
    else:
        src = os.path.join(src_dir,file)
        dst = os.path.join(dst_dir,file)
        # make sure the path of dst exists
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        if not os.path.exists(src):
            if os.path.exists(src+'.zip'):
                # de data is stored as a zip-file. Extract this file.
                with zipfile.ZipFile(src+'.zip','r') as z:
                    z.extractall(os.path.dirname(dst))
            else:
                # If this is not the case, raise an error
                raise(FileNotFoundError('File not found: {}'.format(src)))
        else:
            # just copy the file
            if src!=dst:
                shutil.copyfile(src,dst)

def is_number(s):
    """Test if a string s represents a number"""
    try:
        float(s)
        return True
    except ValueError:
        return False
    
def geodataframe2grid(mf, gdf_in):
    """Cut a geodataframe gdf_in by the grid of a modflow model mf"""
    geom_col_name = gdf_in._geometry_column_name

    # make a polygon for each of the grid-cells
    # this may result in lines on the edge of the polygon counting double
    grid_polygons = []
    if isinstance(mf.modelgrid,flopy.discretization.structuredgrid.StructuredGrid):
        xv,yv,zv = mf.modelgrid.xyzvertices
        for row in range(mf.modelgrid.nrow):
            for col in range(mf.modelgrid.ncol):
                if False:
                    # the next line is very slow in flopy version 3.2.11
                    vert = mf.modelgrid.get_cell_vertices(row, col)
                else:
                    # so we do it ourselves
                    vert = [(xv[row, col], yv[row, col]),
                            (xv[row, col+1], yv[row, col+1]),
                            (xv[row+1, col+1], yv[row+1, col+1]),
                            (xv[row+1, col], yv[row+1, col]),]
                pol = Polygon(vert)
                pol.row = row
                pol.col = col
                grid_polygons.append(pol)
    elif isinstance(mf.modelgrid,flopy.discretization.vertexgrid.VertexGrid):
        xv,yv,zv = mf.modelgrid.xyzvertices
        for cell2d in range(mf.modelgrid.ncpl):
            if False:
                # the next line is very slow in flopy version 3.2.11
                vert = mf.modelgrid.get_cell_vertices(cell2d)
            else:
                # so we do it ourselves
                vert = list(zip(xv[cell2d],yv[cell2d]))
            pol = Polygon(vert)
            pol.cell2d = cell2d
            grid_polygons.append(pol)
    else:
        raise(Exception('Grid type of {} not supported'.format(type(mf.modelgrid))))

    s = STRtree(grid_polygons)

    shp_list = []
    # cut the lines with the grid
    for index, row in gdf_in.iterrows():
        g = row[geom_col_name]
        gtype = g.geometryType()
        result = s.query(g)
        for r in result:
            i = g.intersection(r)
            _add_shapes_to_list(i,gtype,geom_col_name,row,shp_list,r)
    if len(shp_list)==0:
        warnings.warn('No overlap between model and shape')
        shp_out = gdf_in.loc[[]]
        shp_out['row'] = 0
        shp_out['col'] = 0
    else:
        gdf_out = gpd.GeoDataFrame(pd.DataFrame(shp_list), geometry=geom_col_name)
    return gdf_out

def _add_shapes_to_list(i,geometryType,geom_col_name,row,shp_list,r):
    """subfunction of geodataframe2grid"""
    if geometryType == 'LineString' or geometryType == 'MultiLineString':
        if not i.is_empty:
            it = i.geometryType()
            if it == 'GeometryCollection':
                for im in i.geoms:
                    _add_shapes_to_list(im,geometryType,geom_col_name,row,shp_list,r)
            elif it == 'MultiLineString':
                for im in i.geoms:
                    _add_shapes_to_list(im,geometryType,geom_col_name,row,shp_list,r)
            elif it == 'LineString':
                # TODO: to make sure lines are not counted double
                # do not add the line if the line is on the north or west border of the cell-edge
                rown = row.copy()
                rown[geom_col_name] = i
                if hasattr(r,'cell2d'):
                    rown['cell2d'] = r.cell2d
                else:
                    rown['row'] = r.row
                    rown['col'] = r.col
                shp_list.append(rown)
            elif it == 'Point':
                # endpoint of the linestring is on the cell-edge
                pass
            elif it == 'MultiPoint':
                # mutiple endpoints of the linestring are on the cell-edge
                pass
            else:
                raise NotImplementedError(
                    'geometryType ' + it + ' not yet supprted in geodataframe2grid')
    elif geometryType == 'Polygon' or geometryType == 'MultiPolygon':
        it = i.geometryType()
        if it == 'GeometryCollection':
            for im in i.geoms:
                _add_shapes_to_list(im,geometryType,geom_col_name,row,shp_list,r)
        elif it == 'MultiPolygon':
            for im in i.geoms:
                _add_shapes_to_list(im,geometryType,geom_col_name,row,shp_list,r)
        elif it == 'Polygon':
            rown = row.copy()
            rown[geom_col_name] = i
            if hasattr(r,'cell2d'):
                rown['cell2d'] = r.cell2d
            else:
                rown['row'] = r.row
                rown['col'] = r.col
            shp_list.append(rown)
        elif it == 'Point':
            # endpoint of the polygon is on the cell-edge
            pass
        elif it == 'LineString' or it == 'MultiLineString':
            # one of the edges of the polygon is on a cell-egde
            pass
        else:
            raise NotImplementedError(
                'geometryType ' + it + ' not yet supprted in geodataframe2grid')
    else:
        raise NotImplementedError(
            'geometryType ' + geometryType + ' not yet supprted in geodataframe2grid')

def read_file_from_datadir(fname, datadir='data', verbose=True):
    
    fname_zip = os.path.join(datadir, fname + '.zip')
    if os.path.isfile(fname_zip):
        print('{} exists -> unzip file to temp dir'.format(os.path.join(datadir, fname + '.zip')))
        tempdir = tempfile.TemporaryDirectory().name
        
        zip_ref = zipfile.ZipFile(fname_zip, 'r')
        zip_ref.extractall(tempdir)
        zip_ref.close()
        
        fname_temp = os.path.join(tempdir, os.path.split(fname)[-1])
    else:
        fname_temp = os.path.join(datadir,fname)
    
    return fname_temp

def extract_layer_heights(data_dir, dirname):
    from hdsrhipy.model import rasters
    if not os.path.isdir(os.path.join(data_dir,'extracted',dirname)):
        if not os.path.isdir(os.path.join(data_dir,'extracted')):
            os.mkdir(os.path.join(data_dir,'extracted'))
        # extract the modellagen
        zip_ref = zipfile.ZipFile(os.path.join(data_dir,dirname + '.zip'), 'r')
        zip_ref.extractall(os.path.join(data_dir,'extracted'))
        zip_ref.close()
    if dirname == 'modellagen_25x25m':
        fname25 = os.path.join(data_dir,'extracted','modellagen_25x25m','bot_l8.asc')
        if os.path.exists(fname25):
            return
        # the botm level of layer is missing for 25x25m
        # therefore generate the botm level of layer 8, from 100x100m
        extract_layer_heights(data_dir,'modellagen_100x100m')
        fname100 = os.path.join(data_dir,'extracted','modellagen_100x100m','bot_l8.asc')
        ds100 = gdal.Open(fname100)
        if True:
            # this is needed to take care of NaNs
            x, y = rasters.get_xy(ds100)
            z = rasters.get_values(ds100)
            ds100 = rasters.xyz2gdal_dataset(x, y, z)
        ds25 = rasters.reproject_dataset(ds100, 25, GDALResampleAlg=gdal.GRA_NearestNeighbour)
        if False:
            ds25_0 = gdal.Open(os.path.join(data_dir,'extracted','modellagen_25x25m','top_l8'))
        gdal.Translate(fname25,ds25)

def to_flopy(rf,csize=None,xmin=None,xmax=None,ymin=None,ymax=None,
             packages=['dis', 'lpf', 'bas', 'oc', 'pcg', 'wel', 'drn', 
                       'riv', 'ghb', 'rch'],
             data_dir=None, steady_model=False):
    """
    Transform the iMod runfile to a basic FloPy model.

    This makes that methods developped for Flopy-models can be used
    Like plotting cross-sections or geodataframe2grid
    
    The filenames are specified in the runfile (rf). The data is located in
    the data_dir. The data_dir often contains zipped files. These files are
    unzipped to a temporary folder and read as a modflow model.

    Parameters
    ----------
    rf : RunFile
        an instance of a read Runfile from runfile.py
    csize: float
        Cell size of model. Uses rf.data['CSIZE'] when None
    xmin: float
        Minimum x-coordinate of the model. Uses rf.data['XMIN'] when None
    xmax: float
        Maximum x-coordinate of the model. Uses rf.data['XMAX'] when None
    ymin: float
        Minimum y-coordinate of the model. Uses rf.data['YMIN'] when None
    ymax: float
        Maximum y-coordinate of the model. Uses rf.data['YMAX'] when None

    Returns
    -------
    m : flopy Model

    """
    if csize is None:
        csize = rf.data['CSIZE']
    if xmin is None:
        xmin = rf.data['XMIN']
    if xmax is None:
        xmax = rf.data['XMAX']
    if ymin is None:
        ymin = rf.data['YMIN']
    if ymax is None:
        ymax = rf.data['YMAX']
    if packages is None:
        packages = ['dis']
    m = flopy.modflow.Modflow()
    if 'dis' in packages or 'DIS' in packages:
        nlay = rf.data['NLAY']
        nrow = int((ymax-ymin)/csize)
        ncol = int((xmax-xmin)/csize)
        nper = rf.data['NPER']
        if steady_model:
            nper = 1
        if 'STO' in rf.data.keys():
            steady = [False] * nper
        else:
            steady = True
        if steady_model:
            steady = True
        laycbd = [1]*rf.data['NLAY']
        laycbd[-1]=0
        botm = np.full((nlay+np.array(laycbd).sum(),nrow,ncol),np.NaN)
        if np.mod(csize,100) == 0:
            dirname = 'modellagen_100x100m'
        else:
            dirname = 'modellagen_25x25m'
            #raise(Exception('Make sure csize is a multiple of 100 m. The botm level of layer 8 is missing on a schale of 25x25m'))
            
        extract_layer_heights(data_dir, dirname)
        for l in range(nlay):
            fname = os.path.join(data_dir,'extracted',dirname,'top_l{}.asc'.format(l+1))
            if not os.path.exists(fname):
                fname = fname[:-4]
            if not os.path.exists(fname):
                raise(Exception('Kan bestand {} niet vinden'.format(fname)))
            val = get_grid_values(fname,csize,xmin,xmax,ymin,ymax)
            
            if l==0:
                top = val
            else:
                botm[l*2-1] = val
            fname = os.path.join(data_dir,'extracted',dirname,'bot_l{}.asc'.format(l+1))
            if not os.path.exists(fname):
                fname = fname[:-4]
            if not os.path.exists(fname):
                raise(Exception('Kan bestand {} niet vinden'.format(fname)))
            botm[l*2] = get_grid_values(fname,csize,xmin,xmax,ymin,ymax)
        
        # verwijder laagdiktes van 0.0 meter
        for l in range(len(botm)):
#            if l==2:
#                break
            if l ==0:
                top_l = top
            
            while ((top_l-botm[l])<=0.11).sum() > 0:
                botm[l] = np.where((top_l-botm[l])<=0.11, botm[l]-0.11, botm[l])
            top_l= botm[l].copy()
            
#            if l>0:
#                print(l, (botm[l-1]-botm[l]).min())

        dis = flopy.modflow.ModflowDis(m,nlay=nlay,nrow=nrow,ncol=ncol,
                                       delr=csize,delc=csize,
                                       laycbd=laycbd,top=top,botm=botm,
                                       steady=steady,
                                       xul=xmin,yul=ymax, nper=nper)

        
    if 'lpf' in packages:
        
        #hdry = 3.4028235E+38
        ipakcb = 31
        sf1 = 1e-05
        
        
        laycbd = m.dis.laycbd.array
        hk = np.full((m.nlay,m.nrow,m.ncol),np.NaN)
        for ilay,file in zip(rf.data['KDW']['files']['ILAY'],rf.data['KDW']['files']['FNAME']):
            fname = read_file_from_datadir(file, datadir=data_dir, verbose=True)
            kd = get_grid_values(fname,csize,xmin,xmax,ymin,ymax)
            hk[ilay-1] = kd / m.dis.thickness[ilay-1+m.dis.laycbd[:(ilay-1)].sum()].array
            
        vkcb = np.full((m.nlay,m.nrow,m.ncol),np.NaN)
        for ilay,file in zip(rf.data['VCW']['files']['ILAY'],rf.data['VCW']['files']['FNAME']):
            fname = read_file_from_datadir(file, datadir=data_dir, verbose=True)
            cw = get_grid_values(fname,csize,xmin,xmax,ymin,ymax)
            vkcb[ilay-1] = cw / m.dis.thickness[ilay-1+m.dis.laycbd[:(ilay-1)].sum()].array
        
        
        flopy.modflow.ModflowLpf(m, ipakcb=ipakcb, hk=hk, vka=hk, vkcb=vkcb, 
                                 ss=sf1)


        
    if 'bcf' in packages:
        
        # from ouput\mf2005_tmp
        # anders dan default
        hdry = 3.4028235E+38
        ipakcb = 31
        wetfct = 0.0
        iwetit = 0
        laycon = 0
        sf1 = 1e-05
        
        
        # default
        iwdflg = 0
        ihdwet = 0
        trpy = 1
        
        kd = np.full((m.nlay,m.nrow,m.ncol),np.NaN)
        for ilay,file in zip(rf.data['KDW']['files']['ILAY'],rf.data['KDW']['files']['FNAME']):
            fname = read_file_from_datadir(file, datadir=data_dir, verbose=True)
            kd[ilay-1] = get_grid_values(fname,csize,xmin,xmax,ymin,ymax)
                    
        vcont = np.full((m.nlay-1,m.nrow,m.ncol),np.NaN)
        for ilay,file in zip(rf.data['VCW']['files']['ILAY'],rf.data['VCW']['files']['FNAME']):
            fname = read_file_from_datadir(file, datadir=data_dir, verbose=True)
            cw = get_grid_values(fname,csize,xmin,xmax,ymin,ymax)
            vcont[ilay-1] = 1/cw #/ m.dis.thickness[ilay+m.dis.laycbd[:ilay].sum()].array
        
           
        
        flopy.modflow.ModflowBcf(m, hdry=hdry, ipakcb=ipakcb,
                                 wetfct=wetfct, iwetit=iwetit, 
                                 laycon=laycon, sf1=sf1, tran=kd,
                                 vcont=vcont)
        
    if 'bas' in packages:

        #ibound
        ibound = np.full((m.nlay,m.nrow,m.ncol),np.NaN)
        for ilay, file in zip(rf.data['BND']['files']['ILAY'],rf.data['BND']['files']['FNAME']):
            fname = read_file_from_datadir(file, datadir=data_dir, verbose=True)
            ibound[ilay-1] = get_ibound(fname, rf)

        #strt
        strt = np.full((m.nlay,m.nrow,m.ncol),np.NaN)
        for ilay,file in zip(rf.data['SHD']['files']['ILAY'],rf.data['SHD']['files']['FNAME']):
            fname = read_file_from_datadir(file, datadir=data_dir, verbose=True)
            # use min to make sure there are no NaN's in the result
            strt[ilay-1] = get_grid_values(fname, csize, xmin, xmax, ymin, ymax,
                                           resampling=rasterio.warp.Resampling.min)

        flopy.modflow.ModflowBas(m, ibound=ibound, strt=strt)
        
    # OC package niet in runfile, o.b.v. mf2005_tmp export
    if 'oc' in packages:
        spd = {}
        for i in range(dis.nper):
            spd[(i,0)] = ['save head', 'save budget']
        flopy.modflow.ModflowOc(m, stress_period_data=spd)
    
    if 'pcg' in packages:
        flopy.modflow.ModflowPcg(m, mxiter=rf.data['OUTER'], iter1=rf.data['INNER'], 
                                       hclose=rf.data['HCLOSE'], rclose=rf.data['QCLOSE'],
                                       relax=rf.data['RELAX'], nbpol=1, iprpcg=1, mutpcg=0)
        
    if 'wel' in packages:
        spd = {}
        for sp in rf.data['WEL']['files'].keys():
            lrcd = []
            for ilay,file in zip(rf.data['WEL']['files'][sp]['ILAY'],rf.data['WEL']['files'][sp]['FNAME']):
                fname = read_file_from_datadir(file, datadir=data_dir, verbose=True)
                df = ipf2df(fname, m)
                for rc, Q in zip(df['rc'].values, df['Q'].values):
                    lrcd.append([ilay-1, rc[0], rc[1], Q])

            spd[sp] = lrcd
            
        flopy.modflow.ModflowWel(m, stress_period_data=spd)
        
    if 'drn' in packages:
        spd = idfs2spd(rf, data_dir, 'DRN', m, ipars=[1,0])
        flopy.modflow.ModflowDrn(m, stress_period_data=spd)
    
    if 'riv' in packages:
        #ignore infiltration factor
        spd = idfs2spd(rf, data_dir, 'RIV', m, ipars=[1,0,2])
        flopy.modflow.ModflowRiv(m, stress_period_data=spd)
        
    if 'ghb' in packages:
        spd = idfs2spd(rf, data_dir, 'GHB', m, ipars=[1,0])
        if spd != {}:
            flopy.modflow.ModflowGhb(m, stress_period_data=spd)
        else:
            print('no ghb data')
            
    if 'rch' in packages:
        spd = {}
        for sp in rf.data['RCH']['files'].keys():
            for ilay,file in zip(rf.data['RCH']['files'][sp]['ILAY'],rf.data['RCH']['files'][sp]['FNAME']):
                fname = read_file_from_datadir(file, datadir=data_dir, verbose=True)
                rech = get_grid_values(fname,csize,xmin,xmax,ymin,ymax)
            spd[sp] = rech/1000.
        
        if rf.get_inp_index('mete_svat.inp') is not None:
            pass
        elif rf.get_inp_index('mete_grid.inp') is not None:
            pass
        
        flopy.modflow.ModflowRch(m, rech=spd)

    return m

def idfs2spd(rf, data_dir, pack, ml, ipars=None):
    files = rf.data[pack]['files']
    spd = {}
    npars = rf._get_n(pack)
    if ipars is None:
        ipars = range(npars)
    for sp in files.keys():
        lrc = np.full((0,len(ipars)+3),0)
        nlev = int(files[sp].shape[0]/npars)
        for ilev in range(nlev):
            lrci = None
            for ipar in ipars:
                ifile = ipar*nlev + ilev
                fname = read_file_from_datadir(files[sp].at[ifile,'FNAME'], datadir=data_dir)
                da = imod.idf.open(fname)
                if lrci is None:
                    # calculate in which row and column each value is
                    xmid = da.x.values
                    ymid = da.y.values
                    xmid = xmid-ml.modelgrid.xoffset
                    ymid = ymid-ml.modelgrid.yoffset
                    
                    xedge, yedge = ml.modelgrid.xyedges
                    inx = (xmid>xedge.min()) & (xmid<xedge.max())
                    cols = np.full(xmid.shape,np.NaN)
                    cols[inx] = [np.where(x<xedge[1:])[0][0] for x in xmid[inx]]
                    iny = (ymid>yedge.min()) & (ymid<yedge.max())
                    rows = np.full(ymid.shape,np.NaN)
                    rows[iny] = [np.where(y>yedge[1:])[0][0] for y in ymid[iny]]
                    cols, rows = np.meshgrid(cols,rows)
                    lays = np.full(rows.shape, files[sp].at[ifile, 'ILAY']-1)
                    lrci = np.column_stack((lays.flatten(), rows.flatten(), cols.flatten()))
                lrci = np.column_stack((lrci, da.values.flatten()))
            lrc = np.vstack((lrc, lrci))
        mask = ~np.any(np.isnan(lrc),axis=1)
        spd[sp] = lrc[mask]
    return spd


def ipf2df(fname, m):
    
    with open(fname, 'r') as fo:
        nrows = int(fo.readline().strip())
        nparam = int(fo.readline().strip())
        param_list = []
        for i in range(nparam):
            param_list.append(fo.readline().strip())
            
    df = pd.read_csv(fname, skiprows=nparam+3, names=param_list, sep=' ')
    df['rc'] = df.apply(lambda i: m.sr.get_rc(i[param_list[0]], i[param_list[1]]), axis=1).values
    
    return df


def copy_and_unzip(dirname, data_dir):
    
    
    dir_path = os.path.join(data_dir,'extracted',dirname)
    if not os.path.isdir(dir_path):
        os.mkdir(os.path.join(data_dir,'extracted', dirname))
        src = os.path.join(data_dir,dirname)
        copytree(src, dir_path)
    
    # check if there are any zipfiles in the directory
    if any(fname.endswith('.zip') for fname in os.listdir(dir_path)):
        unzip_files_in_dir(dir_path)
    
def unzip_files_in_dir(dir_path):
    
    for f in os.listdir(dir_path):
            
        # unzip all files in directory
        zip_ref = zipfile.ZipFile(os.path.join(dir_path,f), 'r')
        zip_ref.extractall(dir_path)
        zip_ref.close()
        os.remove(os.path.join(dir_path,f))

def copytree(src, dst, symlinks=False, ignore=None):
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.isdir(s):
            shutil.copytree(s, d, symlinks, ignore)
        else:
            shutil.copy2(s, d)

def get_grid_values(fname,csize,xmin,xmax,ymin,ymax,
                    resampling=rasterio.warp.Resampling.average):
    nrow = int((ymax-ymin)/csize)
    ncol = int((xmax-xmin)/csize)
    val = np.full((nrow,ncol),np.NaN)
    dst_transform = rasterio.transform.from_origin(xmin, ymax, csize, csize)
    if fname.lower().endswith('.idf'):
        fname = fname.replace('\\','/')
        da = imod.idf.open(fname)
        src_transform = imod.util.transform(da)
        source = da.values
        rasterio.warp.reproject(source, val, src_transform=src_transform,
                                dst_transform=dst_transform,
                                src_crs=28992, dst_crs=28992,
                                resampling=resampling, src_nodata=np.NaN)
    else:
        ds = rasterio.open(fname)
        source = rasterio.band(ds,1)
        rasterio.warp.reproject(source, val, dst_transform=dst_transform,
                                src_crs=28992, dst_crs=28992,
                                resampling=resampling)
        if ds.nodata is not None:
            val[val==ds.nodata] = np.NaN
    
    return val
    
def calculate_lek_conductance(ml, wl):
    pass

def extract_changed_files(zipname,pathname,check_time=True,check_size=False,
                          debug=False):
    # Extract each file in a zip-file only when the properties are different
    # With the default arguments this method only checks the modification time
    with zipfile.ZipFile(zipname) as zf:
        infolist = zf.infolist()
        for info in infolist:
            fname = os.path.join(pathname,info.filename)
            extract=False
            if os.path.exists(fname):
                if check_time:
                    tz = time.mktime(info.date_time + (0, 0, -1))
                    tf = os.path.getmtime(fname)
                    if tz!=tf:
                        extract=True
                if check_size:
                    sz = info.file_size
                    sf = os.path.getsize(fname)
                    if sz!=sf:
                        extract=True
            else:
                extract=True
            if extract:
                if debug:
                    print('extracting {}'.format(info.filename))
                zf.extract(info.filename,pathname)
                # set the correct modification time
                # (which is the time of extraction by default)
                tz = time.mktime(info.date_time + (0, 0, -1))
                os.utime(os.path.join(pathname,info.filename), (tz, tz))

def write_idf(x,y,values,fname,make_zip=False):
    a = xarray.DataArray(values,dims=['y','x'],coords=[y,x])
    # make sure the directory existst
    pathname = os.path.dirname(fname)
    if not os.path.isdir(pathname):
        os.makedirs(pathname)
    imod.idf.write(fname,a)
    if make_zip:
        # make a zip of the result, so it requires less space on Git LFS
        compress_idfs(fname)
        
def prepare_raw_data(data_path=None):
    print('Prepare bofek_unsa-files')
    
    # make sure the downloads directory exist
    pathname = os.path.join(data_path,'downloads')
    if not os.path.isdir(pathname):
        os.makedirs(pathname)
    # make sure the extracted directory exist
    pathname = os.path.join(data_path,'extracted')
    if not os.path.isdir(pathname):
        os.makedirs(pathname)
    
    # check for new data (extracted dir is in gitignore)
    extract_changed_files(os.path.join(data_path,'metaswap43','bofek_unsa_1_10.zip'),os.path.join(data_path,'extracted'))
    extract_changed_files(os.path.join(data_path,'metaswap43','bofek_unsa_11_20.zip'),os.path.join(data_path,'extracted'))
    extract_changed_files(os.path.join(data_path,'metaswap43','bofek_unsa_21_30.zip'),os.path.join(data_path,'extracted'))
    extract_changed_files(os.path.join(data_path,'metaswap43','bofek_unsa_31_40.zip'),os.path.join(data_path,'extracted'))
    extract_changed_files(os.path.join(data_path,'metaswap43','bofek_unsa_41_50.zip'),os.path.join(data_path,'extracted'))
    extract_changed_files(os.path.join(data_path,'metaswap43','bofek_unsa_51_60.zip'),os.path.join(data_path,'extracted'))
    extract_changed_files(os.path.join(data_path,'metaswap43','bofek_unsa_61_72.zip'),os.path.join(data_path,'extracted'))
    extract_changed_files(os.path.join(data_path,'bodem','bofek_nl.zip'),os.path.join(data_path,'extracted'))
    
    if False:
        print('Prepare meteobase-files')
        # download and extract data from meteobase, if needed
        for year in range(1990,2018):
            # download one year at a time
            pathname = os.path.join('data','meteobase')
            if not os.path.isdir(pathname):
                os.mkdir(pathname)
            zipname = os.path.join(pathname,'{}.zip'.format(year))
            if not os.path.exists(zipname):
                # download from meteobase
                print('Downloading meteobase data from {}'.format(year))
                bounds = [105000, 433000, 173000, 473000]
                start_date = pd.to_datetime(str(year))
                end_date = pd.to_datetime(str(year+1))
                get_meteo(bounds,start_date,end_date,zipname)
            # extract data from meteobase, when not done allready
            to_path = os.path.join('data\extracted\meteobase',str(year))
            extract_changed_files(zipname,to_path)
            
def polygons2grid(gdf,m,field,method='weighted_average',areacol='area'):
    """Convert polygon-data to a raster grid"""
    # only use polygons that intersect with the model domain
    bounds = get_bounds_polygon(m)
    mask = gdf.intersects(bounds)
    gdf = gdf.loc[mask]
    gdf.geometry = gdf.intersection(bounds)
    if method == 'weighted_average':
        # use a weighted average of each polygon in a cell
        gdf = geodataframe2grid(m,gdf)
        grid=np.full((m.dis.nrow,m.dis.ncol),0.)
        area=np.full((m.dis.nrow,m.dis.ncol),0.)
        for ind,pol in gdf.iterrows():
            grid[pol.row,pol.col] = grid[pol.row,pol.col] + pol.geometry.area * pol[field]
            area[pol.row,pol.col] = area[pol.row,pol.col] + pol.geometry.area
        mask = area>0
        grid[mask] = grid[mask]/area[mask]
        grid[~mask]=np.NaN
    elif method == 'center_grid':
        # the value in each cell is detemined by the polygon in which the cell-center falls
        x = m.sr.xcentergrid.flatten()
        y = m.sr.ycentergrid.flatten()
        grid=np.full((m.dis.nrow,m.dis.ncol),np.NaN)
        for ind,pol in gdf.iterrows():
            mask = [pol.geometry.contains(Point(x[i],y[i])) for i in range(len(x))]
            if np.any(mask):
                grid[np.reshape(mask,(m.nrow,m.ncol))] = pol.VALUE
    elif method == 'most_common':
        # use the value that is most common (has the largest area) in each cell
        gdf = geodataframe2grid(m,gdf)
        gdf[areacol] = gdf.area
        area=np.full((m.dis.nrow,m.dis.ncol),np.NaN)
        grid=np.full((m.dis.nrow,m.dis.ncol),np.NaN)
        val_un = gdf[field].unique()
        for val in val_un:
            mask = gdf[field] == val
            area_per_cel = gdf.loc[mask].groupby(['row','col'])[areacol].sum()
            row = area_per_cel.index.get_level_values(0)
            col = area_per_cel.index.get_level_values(1)
            mask = np.isnan(grid[row,col]) | (area_per_cel>area[row,col])
            area[row[mask],col[mask]] = area_per_cel.loc[mask]
            grid[row[mask],col[mask]] = val
    return grid
            
def get_bounds_polygon(m):
    """make a polygon of the model boundary"""
    extent = m.sr.get_extent()
    return extent_to_polygon(extent)

def extent_to_polygon(extent):
    """Make a Polygon of an extent ([xmin,xmax,ymin,ymax])"""
    bounds = Polygon([(extent[0],extent[2]),
                      (extent[0], extent[3]),
                      (extent[1], extent[3]),
                      (extent[1],extent[2])])
    return bounds

def get_ibound(idf_name,rf):
    """Calculate the ibound from idf_name belonging to the model runfile rf"""
    ibound = get_grid_values(idf_name, rf.data['CSIZE'],
                             rf.data['XMIN'], rf.data['XMAX'],
                             rf.data['YMIN'], rf.data['YMAX'],
                             resampling=rasterio.warp.Resampling.min)
    # also make the border -1
    mask=np.full(ibound.shape,True)
    mask[1:-1,1:-1]=False
    mask = mask & (ibound==1)
    ibound[mask] = -1
    return ibound

def add_simgro_surface_water(rf, gdf_pg=None, run1period=True, change_conductance=True, datadir=None):
    """Add surface water (peilgebieden) from Simgro to an iMod-simulation"""
    print('Running one-timestep simulation to retrieve initial inp-files')
    # first run only the frist period
    rf2 = copy.deepcopy(rf)
    # only run the first timestep
    rf2.change_period(rf2.data['DATE'][0],rf2.data['DATE'][0])
    naam='run1period_'+str(rf.data['CSIZE'])
    model_ws = os.path.join(datadir,'work',naam)
    from hdsrhipy.model import hdsr
    if run1period:
        rf2.run_imodflow(model_ws,naam,data_path=datadir)

    # if gdf_pg is None:
    #     # download the layer with peilgebieden
    #     gdf_pg = hdsr.get_peilgebieden()
    print('Starting flopy operations')
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
    rf.add_inp(mana_res,'simgro/mana_res.inp')
    
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
    rf.add_inp(goto_res,'simgro/goto_res.inp')
    
    # dish_res.inp
    # fswr is 1.5 l/s/ha = 1.5/1000/10000 m3/s/m2
    dish_res = pd.DataFrame({'swnr':gdf_sw.index,'dhwr':0.1,'fswr':1.5/1000/10000*gdf_sw['area_grid'],'swnrgo':0})
    dish_res0 = pd.DataFrame({'swnr':gdf_sw.index,'dhwr':0.0,'fswr':0.0,'swnrgo':0})
    dish_res = dish_res0.append(dish_res, ignore_index=True).sort_values(['swnr','dhwr']).reset_index(drop=True)
    rf.add_inp(dish_res,'simgro/dish_res.inp')
    
    # swnr_sim.inp
    swnr_sim = pd.DataFrame({'swnr':gdf_sw.index})
    rf.add_inp(swnr_sim,'simgro/swnr_sim.inp')
    
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
    
    rf.add_inp(svat2swnr_roff,'simgro/svat2swnr_roff.inp')
        
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
            cond_file = 'OPPERVLAKTEWATER/\WINTER/CONDUCTANCE_LAAG{}_{}.IDF'.format(lay,lev)
            botm_file = 'OPPERVLAKTEWATER/WINTER/BODEMHOOGTE_LAAG{}_{}.IDF'.format(lay,lev)
            infi_file = 'OPPERVLAKTEWATER/WINTER/INFFACTOR_LAAG{}_{}.IDF'.format(lay,lev)
            svat2swnr_drng = idf2svat2swnr_drng(svat2swnr_drng,model_ws,cond_file,
                                                botm_file,area_svat,svat2mod,grid_sw,
                                                lrc,rf,ml,'RIV',infi_file=infi_file,sy=lev,
                                                change_conductance=change_conductance, datadir=datadir)
    
    # add drains from 
    cond_file = 'DRAINAGE/CONDUCTANCE_NIVEAU1.IDF'
    botm_file = 'DRAINAGE/DRAINAGE_NIVEAU1.IDF'
    svat2swnr_drng = idf2svat2swnr_drng(svat2swnr_drng,model_ws,cond_file,
                                        botm_file,area_svat,svat2mod,grid_sw,
                                        lrc,rf,ml,'DRN',sy=4,
                                        change_conductance=change_conductance, datadir=datadir)
    rf.add_inp(svat2swnr_drng,'simgro/svat2swnr_drng.inp')    

    
        
def read_dxc(fname):
    f = open(fname, "r")
    l = f.readline()
    l = f.readline()
    n = int(l.strip())
    lrcs = np.loadtxt(f,dtype=int)
    return lrcs

def xymid_from_attrs(attrs):
    tr = attrs['transform']
    X = np.linspace(tr[2]+tr[0]/2,tr[2]+tr[0]*(attrs['ncol']-0.5),attrs['ncol'])
    Y = np.linspace(tr[5]+tr[4]/2,tr[5]+tr[4]*(attrs['nrow']-0.5),attrs['nrow'])
    return X, Y

def cond2svat2swnr_drng(svat2swnr_drng, cond, wd, li, Bin, area_svat,
                        svat2mod, grid_sw, goto_res, lrc, rf, ml, sy=0):
    row = cond.index.get_level_values('row')
    col = cond.index.get_level_values('col')
    
    # get the svat of the from the row and column
    rc2mod = np.zeros((ml.nrow, ml.ncol), dtype=int)
    rc2mod[lrc['r'],lrc['c']] = lrc.index
    # then select cells that are within a swnr
    # and have a modcell greater than 0 (which means they are coupled to metaswap)
    swnr = grid_sw[row,col]
    mask = ~np.isnan(swnr) & (rc2mod[row,col]>0)
    
    # find the svat within each cell that is not urban area (luk 18)
    area_svat['modfcell'] = svat2mod['modfcell']
    
    # assume there is a non-urban landuse in every cell
    nu = area_svat['luk']!=18
    # if there is only one svat in a modfcell, use this svat
    # even if it is urban area
    for mfc in area_svat['modfcell'].unique():
        inmfc = area_svat['modfcell']==mfc
        if inmfc.sum()==1:
            nu[inmfc] = True
                
    mod2svat = pd.Series(area_svat.loc[nu,'svat'].values,index = area_svat.loc[nu,'modfcell'])
    svat = mod2svat.loc[rc2mod[row[mask],col[mask]]]
    assert not np.any(np.isnan(svat))
    # add data to svat2swnr_drng
    peil = goto_res.loc[swnr[mask],['lvwrsm','lvwrwt']].mean(axis=1)
    # drain depth is surface level - bodem (= peil - waterdiepte)
    dpsw = area_svat.loc[svat,'glk'].values-(peil.values-wd[mask].values)
    
    redr = (Bin[mask] * li[mask]) /cond[mask]
    if sy==1:
        # hoofdwatergangen kunnen infiltreren
        rein = redr.values
    else:
        # andere niet
        rein = 1.0E+05
        
    df = pd.DataFrame({'svat':svat.values,
                       'sy':sy,
                       'dpsw':dpsw,
                       'wisw':Bin[mask].values,
                       'adsw':0.0,
                       'lesw':li[mask].values,
                       'redr':redr.values,
                       'reen':0.0,
                       'rein':rein,
                       'reex':0.0})
    svat2swnr_drng = svat2swnr_drng.append(df)
    
    return svat2swnr_drng
    

def idf2svat2swnr_drng(svat2swnr_drng,model_ws,cond_file,botm_file,area_svat,
                       svat2mod,grid_sw,lrc,rf,ml,pack,infi_file=None,sy=0,
                       change_conductance=False, datadir=None):
    """Replace a surface water idf by items in svat2swnr_drng"""
    cond = imod.idf.open(os.path.join(model_ws,cond_file))
    botm = imod.idf.open(os.path.join(model_ws,botm_file))
    assert np.array_equal(cond.x, botm.x) and np.array_equal(cond.y, botm.y)
    if infi_file:
        infi = imod.idf.open(os.path.join(model_ws,infi_file))
        assert np.array_equal(botm.x, infi.x) and np.array_equal(botm.y, infi.y)
    # determine the x and y locations (centers) of the data
    X = cond.x
    Y = cond.y
    # determine row and column in modflow grid
    # assume the idf and the iMod-grid have no rotation
    row = np.full(cond.shape,ml.dis.nrow) # use nrow as nan
    for ir in range(ml.dis.nrow):
        mask = (Y<ml.sr.ygrid[ir,0]) & (Y>=ml.sr.ygrid[ir+1,0])
        row[mask,:] = ir
    col = np.full(cond.shape,ml.dis.ncol) # use ncol as nan
    for ic in range(ml.dis.ncol):
        mask = (X>ml.sr.xgrid[0,ic]) & (X<=ml.sr.xgrid[0,ic+1])
        col[:,mask] = ic
    # which cells are determined by the surface water in simgro?
    # first select drains that are not NaN and are within the model domain
    mask = ~np.isnan(cond.values) & (row<ml.dis.nrow) & (col<ml.dis.ncol)
    # get the svat of the from the row and column
    rc2mod = np.zeros((ml.nrow, ml.ncol), dtype=int)
    rc2mod[lrc['r'],lrc['c']] = lrc.index
    # then select cells that are within a swnr
    # and have a modcell greater than 0 (which means they are coupled to metaswap)
    mask[mask] = ~np.isnan(grid_sw[row[mask],col[mask]]) & (rc2mod[row[mask],col[mask]]>0)
    if not change_conductance:
        # find the svat within each cell that is not urban area (luk 18)
        area_svat['modfcell'] = svat2mod['modfcell']
        
        # assume there is a non-urban landuse in every cell
        nu = area_svat['luk']!=18
        # if there is only one svat in a modfcell, use this svat
        # even if it is urban area
        for mfc in area_svat['modfcell'].unique():
            inmfc = area_svat['modfcell']==mfc
            if inmfc.sum()==1:
                nu[inmfc] = True
                    
        mod2svat = pd.Series(area_svat.loc[nu,'svat'].values,index = area_svat.loc[nu,'modfcell'])
        svat = mod2svat.loc[rc2mod[row[mask],col[mask]]]
        assert not np.any(np.isnan(svat))
        # add data to svat2swnr_drng
        weerstand = area_svat.loc[svat,'ark']/cond.values[mask]
        if infi_file:
            rein = weerstand/infi.values[mask]
        else:
            rein = 1.0E+05 # hoog indien buisdrain
        # 'swnr':grid_sw[row[mask],col[mask]].astype(int)
        df = pd.DataFrame({'svat':svat.values,
                           'sy':sy,
                           'dpsw':area_svat.loc[svat,'glk'].values-botm.values[mask],
                           'wisw':1.0,
                           'adsw':0.0,
                           'ddsw':weerstand,
                           'redr':weerstand,
                           'reen':0.0,
                           'rein':rein,
                           'reex':0.0})
        svat2swnr_drng = svat2swnr_drng.append(df)
    
    # remove data from idf
    cond = cond.where(~mask, 0.)
    # and write new idf to file
    cond_file_new = os.path.join('extracted',pack,os.path.basename(cond_file))
    imod.idf.write(os.path.join(datadir,cond_file_new),cond)
    # replace IDFs in runfile
    for key in rf.data[pack]['files'].keys():
        mask = rf.data[pack]['files'][key]['FNAME'].str.lower().str.contains(os.path.basename(cond_file.lower()))
        assert np.sum(mask)==1
        rf.data[pack]['files'][key].loc[mask,'FNAME'] = cond_file_new
    return svat2swnr_drng

def get_idf_extent(fname):
    da = imod.idf.open(fname)
    da.plot.imshow
    xstep = (da.x[1]-da.x[0]) / 2
    ystep = (da.y[1]-da.y[0]) / 2
    left, right = da.x[0] - xstep, da.x[-1] + xstep
    bottom, top = da.y[-1] + ystep, da.y[0] - ystep
    return np.array([left, right, bottom, top])
    
def export_ds(ds,fname,driver="AAIGrid"):
    """Export a gdal-dataset to a raster-file.
    
    Parameters
    ----------
    ds : gdal DataSet
        the gdal dataset that contains the raster and coordinate data 
    fname : str
        The filename of the new file
    driver : str
        The driver for the output-file (ascii raster, geoTiff, etc.)
        supported drivers in https://www.gdal.org/frmt_various.html
    
    """
    # supported drivers in https://www.gdal.org/frmt_various.html
    gdal.GetDriverByName(driver).CreateCopy(fname, ds, strict=0)

def export_grid(z,attrs,fname,driver='AAIGrid'):
    """Export a grid z to a raster-file.
    
    Parameters
    ----------
    z : numpy.ndarray
        the raster data
    attrs : OrderedDict
        the grid attributes, which contain the raster coordinates
        the second output from imod.idf.read()
    fname : str
        The filename of the new file
    driver : str
        The driver for the output-file (ascii raster, geoTiff, etc.)
        supported drivers in https://www.gdal.org/frmt_various.html
        
    """
    tr = attrs['transform']
    x = np.linspace(tr[2],tr[2]+tr[0]*attrs['ncol'],attrs['ncol']+1)
    y = np.linspace(tr[5],tr[5]+tr[4]*attrs['nrow'],attrs['nrow']+1)
    ds = rasters.xyz2gdal_dataset(x,y,z)
    export_ds(ds,fname,driver)

def unzip_file(src, dst, force=False, preserve_datetime=False):
    if os.path.exists(dst):
        if not force:
            print("File not unzipped. Destination already exists. Use 'force=True' to unzip.")
            return
    if preserve_datetime:
        zipf = zipfile.ZipFile(src, 'r')
        for f in zipf.infolist():
            zipf.extract(f, path=dst)
            date_time = time.mktime(f.date_time + (0, 0, -1))
            os.utime(os.path.join(dst, f.filename), (date_time, date_time))
        zipf.close()
    else:
        zipf = zipfile.ZipFile(src, 'r')
        zipf.extractall(dst)
        zipf.close()
    return 1


def interp_weights(xy, uv, d=2):
    """Calculate interpolation weights.

    Parameters
    ----------
    xy : np.array
        array containing x-coordinates in first column and y-coordinates
        in second column
    uv : np.array
        array containing coordinates at which interpolation weights should
        be calculated, x-data in first column and y-data in second column
    d : int, optional
        dimension of data? (the default is 2, which works for 2D data)

    Returns
    -------
    vertices: np.array
        array containing interpolation vertices

    weights: np.array
        array containing interpolation weights per point

    Reference
    ---------
    https://stackoverflow.com/questions/20915502/speedup-scipy-griddata-for-multiple-interpolations-between-two-irregular-grids

    """

    tri = qhull.Delaunay(xy)
    simplex = tri.find_simplex(uv)
    vertices = np.take(tri.simplices, simplex, axis=0)
    temp = np.take(tri.transform, simplex, axis=0)
    delta = uv - temp[:, d]
    bary = np.einsum('njk,nk->nj', temp[:, :d, :], delta)
    return vertices, np.hstack((bary, 1 - bary.sum(axis=1, keepdims=True)))


def interpolate(values, vtx, wts):
    """interpolate values at locations defined by vertices and points,
       as calculated by interp_weights function.

    Parameters
    ----------
    values : np.array
        array containing values to interpolate
    vtx : np.array
        array containing interpolation vertices, see interp_weights()
    wts : np.array
        array containing interpolation weights, see interp_weights()
    
    Returns
    -------
    arr: np.array
        array containing interpolated values at locations as given by
        vtx and wts

    Reference
    ---------
    https://stackoverflow.com/questions/20915502/speedup-scipy-griddata-for-multiple-interpolations-between-two-irregular-grids

    """

    return np.einsum('nj,nj->n', np.take(values, vtx), wts)

def load_heads(rf, model_ws, tmin, tmax, kind='head'):
    """Load the heads for all timesteps and layers in a 4d-array
    
    Parameters
    ----------
    rf : RunFile
        an instance of a read Runfile from runfile.py
    model_ws : string
        model working directory
    tmin : datetime
        The start of the desired output
    tmax : datetime
        The end of the desired output
    
    Returns
    -------
    hgr : np.array
        array containing the head values
    attrs : OrderedDict
        the grid attributes, which contain the raster coordinates
    """
    
    nlay = int(rf.data["NLAY"])
    nrow = int((rf.data["YMAX"] - rf.data["YMIN"]) / rf.data["CSIZE"])
    ncol = int((rf.data["XMAX"] - rf.data["XMIN"]) / rf.data["CSIZE"])
    nper = int((tmax - tmin).total_seconds()/(3600*24) + 1)
    hgr = np.empty((nper, nlay, nrow, ncol))
    for i in range(hgr.shape[0]):
        date = (tmin + pd.Timedelta(days=i)).strftime("%Y%m%d")
        for ilay in range(hgr.shape[1]):
            file = '{}_{}_l{}.idf'.format(kind, date, ilay + 1)
            fname = os.path.join(model_ws, rf.data['OUTPUTDIRECTORY'], kind, file)
            head, attrs = imod.idf.read(fname)
            hgr[i, ilay] = head
    return hgr, attrs


def grids_to_lrcd(mask, grid1=None, grid2=None, grid3=None):
    '''convert grids into a layer, row, column, data list
    
    
    CHD: col1 = head at start of stress period
    DRN: col1 = elevation
    GHB: col1 = head
    RIV: col1 = stage
    
    CHD: col2 = head at end of stress period
    DRN: col2 = conductivity
    GHB: col2 = conductivity
    RIV: col2 = conductivity
    
    RIV: col3 = bottom
    
    Parameters
    ----------
    mask : 3d array
        True if cell should be in lrcd
        
        
    
    '''
    lay, row, col = np.where(mask)
    lrcd = []
    for i in range(len(lay)):
        ind = lay[i], row[i], col[i]
        if grid2 is None:
            # Use 1 datalabel
            lrcd.append([lay[i], row[i], col[i], float(grid1[ind])])
        elif grid3 is None:
            # Use 2 datalabels
            lrcd.append([lay[i], row[i], col[i], float(grid1[ind]), float(grid2[ind])])
        else:
            # Use 3 datalabels
            lrcd.append([lay[i], row[i], col[i], float(grid1[ind]), float(grid2[ind]), float(grid3[ind])])

    return lrcd

def get_xy_mid(rf):
    d = rf.data
    x = np.arange(d['XMIN']+d['CSIZE']/2, d['XMAX'], d['CSIZE'])
    y = np.arange(d['YMAX']-d['CSIZE']/2, d['YMIN'], -d['CSIZE'])
    return x, y

def get_output_array(rf, model_ws=None, kind='head', lays=None, dates=None, subdir='', sys=None):
    pathname = os.path.join(model_ws,rf.data['OUTPUTDIRECTORY'],subdir,kind)
    if True:
        if isinstance(lays, int):
            lays = [lays]
        if dates is None:
            dates = rf.data['DATE'][np.where(rf.data['ISAVE'][:rf.data['NPER']])[0]]
        elif isinstance(dates, str):
            dates = [dates]
        if lays is None:
            lays = np.arange(rf.data['NLAY'])+1
        for idate, date in enumerate(dates):
            for ilay, lay in enumerate(lays):
                if sys is None:
                    fname = os.path.join(pathname, '{}_{}_l{}.idf'.format(kind,date,lay))
                    if not os.path.isfile(fname):
                        raise(FileNotFoundError('file not found: {}'.format(fname)))
                    head, attrs = imod.idf.read(fname)
                else:
                    for si, sy in enumerate(sys):
                        fname = os.path.join(pathname, '{}_sys{}_{}_l{}.idf'.format(kind,sy,date,lay))
                        if not os.path.isfile(fname):
                            raise(FileNotFoundError('file not found: {}'.format(fname)))
                        head_sys, attrs = imod.idf.read(fname)
                        if si==0:
                            head = head_sys
                        else:
                            head = head + head_sys
                if idate==0 and ilay==0:
                    heads = np.full((len(dates), len(lays), head.shape[0], head.shape[1]), np.NaN)
                heads[idate, ilay :, :] = head
    else:
        if sys is None:
            path = os.path.join(pathname,'{}_*'.format(kind))
        else:
            path = os.path.join(pathname,'{}_sys{}_*'.format(kind,sys))
        heads = imod.idf.open(path)
        if dates is not None:
            heads = heads.sel(time=pd.to_datetime(dates))
        if lays is not None:
            heads = heads.sel(layer=lays)
        heads = heads.values
    return heads

def get_metaswap_ouput_array(rf, model_ws, kind, lays=1, dates=None, subdir='metaswap'):
    if dates is None:
        dates = get_metaswap_ouput_dates(rf, model_ws)
    bgt = get_output_array(rf, model_ws, kind=kind, subdir=subdir, dates=dates,
                           lays=lays)
    return bgt.squeeze()

def get_metaswap_output_extent(rf, model_ws, kind, date=None, lay=1, subdir='metaswap'):
    if date is None:
        date = get_metaswap_ouput_dates(rf, model_ws)[0]
    fname = os.path.join(model_ws, rf.data['OUTPUTDIRECTORY'],subdir,kind,
                         '{}_{}_l{}.idf'.format(kind, date, lay))
    return get_idf_extent(fname)

def get_metaswap_ouput_dates(rf, model_ws):
    tiop_sim = rf.get_inp('tiop_sim.inp')
    if isinstance(tiop_sim,str):
        tiop_sim = inp.read(tiop_sim, model_ws)
    t = pd.to_datetime([str(y) for y in tiop_sim['iy']]) + pd.to_timedelta(tiop_sim['td']-1, 'd')
    dates = rf.data['DATE'][pd.to_datetime(rf.data['DATE']).isin(t)].values
    return dates

def inpolygon(x,y,polygon,engine='matplotlib'):
    """find out which points defined by x and y are within polygon

    Parameters
    ----------
    x : np.array
        x-coordinates of grid (same shape as y)
    y : np.array
        y-coordinates of grid (same shape as x)
    poolygon : shapely Polygon or MuliPolygon
        the polygon for which you want mask to be True
    engine : str
        Use 'matplotlib' for speed, for all other values it uses shapely

    Returns
    -------
    mask: np.array
        an array of the same shape as x and y: True for points within polygon
    
    """
    shape = x.shape
    points = list(zip(x.flatten(),y.flatten()))
    if engine == 'matplotlib':
        if isinstance(polygon, MultiPolygon):
            mask = np.full((len(points)),False)
            for pol2 in polygon:
                if not isinstance(pol2, Polygon):
                    raise(Exception('{} not supported'.format(type(pol2))))
                if isinstance(pol2.boundary, MultiLineString):
                    xb,yb = pol2.boundary[0].xy
                else:
                    xb,yb = pol2.boundary.xy
                path = Path(list(zip(xb,yb)))
                mask = mask | path.contains_points(points)
        elif isinstance(polygon,Polygon):
            if isinstance(polygon.boundary, MultiLineString):
                xb,yb = polygon.boundary[0].xy
            else:
                xb,yb = polygon.boundary.xy
            path = Path(list(zip(xb,yb)))
            mask = path.contains_points(points)
        else:
            raise(Exception('{} not supported'.format(type(polygon))))
    else:
        mask = [polygon.contains(Point(x,y)) for x,y in points]
        mask = np.array(mask)
    return mask.reshape(shape)

def read_key_header_lines(f, d, n):
    for i in range(n):
        l = f.readline()
        if l.startswith('*'):
            comment = l[1:].strip()
            if len(comment)>0:
                d['COMMENTS'].append(comment)
        else:
            p = l.strip().split('=')
            d[p[0].strip()] = p[1].strip()

def read_key_file(fname_key):
    # read the keyfile
    f = open(fname_key, "r")
    d = {}
    d['COMMENTS'] = []
    read_key_header_lines(f, d, 10)
    d['PERIOD'] = int(d['PERIOD'])
    d['NUMVAR'] = int(d['NUMVAR'])
    l = f.readline()
    var_header = l[1:].strip().split()
    lst = []
    for ivar in range(d['NUMVAR']):
        l = f.readline().strip()
        lst.append(l.split(maxsplit=2))
    d['VAR'] = pd.DataFrame(lst, columns=var_header).set_index('Variable_name')
    
    read_key_header_lines(f, d, 3)
    d['NUMPTS'] = int(d['NUMPTS'])
    l = f.readline()
    pts_header = l[1:].strip().split()
    drop = []
    for col in range(len(pts_header)):
        if '(' in pts_header[col] and ')' in pts_header[col]:
            pts_header[col-1] = '{} {}'.format(pts_header[col-1],pts_header[col])
            drop.append(col)
    for index in sorted(drop, reverse=True):
        del pts_header[index]
    pts_header = [col.strip() for col in pts_header]
    lst = []
    for ipts in range(d['NUMPTS']):
        l = f.readline().strip()
        lst.append(l.split(maxsplit=3))
    d['PTS'] = pd.DataFrame(lst, columns=pts_header, dtype=float).set_index('id')
    d['PTS'].index = d['PTS'].index.astype(int)
    return d

def fast_weighted_average(df,data_col,weight_col,by_col):
    df['_data_times_weight'] = df[data_col]*df[weight_col]
    df['_weight_where_notnull'] = df[weight_col]*pd.notnull(df[data_col])
    g = df.groupby(by_col)
    result = g['_data_times_weight'].sum() / g['_weight_where_notnull'].sum()
    del df['_data_times_weight'], df['_weight_where_notnull']
    return result

def df2gdf(df,xcol='x',ycol='y'):
    """Make a GeoDataFrame from a DataFrame, assuming the geometry are points"""
    gdf = gpd.GeoDataFrame(df.copy(), geometry=[Point((s[xcol],s[ycol])) for i,s in df.iterrows()])
    return gdf

def excel2datetime(excel_datenum, freq="D"):
    """Method to convert excel datetime to pandas timetime objects.

    Parameters
    ----------
    excel_datenum: datetime index
        can be a datetime object or a pandas datetime index.
    freq:

    Returns
    -------
    datetimes: pandas.datetimeindex

    """
    datetimes = pd.to_datetime('1899-12-30') + pd.to_timedelta(excel_datenum, freq)
    return datetimes

def get_layer(dis, i, j, elev):
    """Return the layers for elevations at i, j locations.
    Edited from flopy, takes into account quasi 3d-layers
    Location in the quasi-3d layers belong to the layer above it (like Modpath)

    Parameters
    ----------
    dis : flopy.modflow.ModflowDis object
    i : scaler or sequence
        row index (zero-based)
    j : scaler or sequence
        column index
    elev : scaler or sequence
        elevation (in same units as model)

    Returns
    -------
    k : np.ndarray (1-D) or scalar
        zero-based layer index
    """

    def to_array(arg):
        if not isinstance(arg, np.ndarray):
            return np.array([arg])
        else:
            return arg

    i = to_array(i)
    j = to_array(j)
    elev = to_array(elev)
    botm_ind = np.arange(0, dis.nlay) + np.cumsum(dis.laycbd)
    botms = dis.botm.array[botm_ind, i, j].tolist()
    layers = np.sum(((botms - elev) > 0), axis=0)
    # do not force elevations below model bottom into bottom layer
    # layers[layers > dis.nlay - 1] = dis.nlay - 1
    layers = np.atleast_1d(np.squeeze(layers))
    if len(layers) == 1:
        layers = layers[0]
    return layers

def get_bda_data(rf, model_ws, bda='svat_dtgw.bda', var=None):
    fname_bda = os.path.join(model_ws, rf.data['OUTPUTDIRECTORY'],'metaswap',bda)
    data = np.fromfile(fname_bda, dtype=np.float32)
    
    fname_key = '{}.key'.format(fname_bda[:-4])
    key = read_key_file(fname_key)
    
    nt = int(len(data)/key['NUMVAR']/key['NUMPTS'])
    data = np.reshape(data, (nt, key['NUMPTS'], key['NUMVAR']))
    if var is not None:
        ivar = np.where(key['VAR'].index==var)[0][0]
        data = data[:,:,ivar]
    
    return data, key

def yd_to_t(year, day):
    """Convert year and day in year to datetime."""
    return pd.to_datetime(year, format='%Y')+pd.to_timedelta(day,'d')

def get_bda_time(key, rf, model_ws):
    fname = os.path.join(model_ws, rf.data['OUTPUTDIRECTORY'],'metaswap', key['TIMERFILE'])
    dy = np.loadtxt(fname)
    t = pd.to_datetime(dy[:,1], format='%Y')+pd.to_timedelta(dy[:,0],'d')
    t = yd_to_t(dy[:,1],dy[:,0])
    return t

def replace_with_dict(ar, dic):
    # Extract out keys and values
    k = np.array(list(dic.keys()))
    v = np.array(list(dic.values()))

    # Get argsort indices
    sidx = k.argsort()

    # Drop the magic bomb with searchsorted to get the corresponding
    # places for a in keys (using sorter since a is not necessarily sorted).
    # Then trace it back to original order with indexing into sidx
    # Finally index into values for desired output.
    return v[sidx[np.searchsorted(k,ar,sorter=sidx)]]

def lek_conductance(A, H0, kv, kh, c1, li, Bin, c0):
    """Calculates de lekconductance according to De Lange
    
    Parameters
    ----------
    A : float
        celoppervlak (m2)
    H0 : float
        doorstroomde dikte (m)
    kv : float
        verticale doorlotendheid (m/d)
    kh : float
        horizontale doorlatendheid (m/d)
    c1 : float
        deklaagweerstand (d)
    li : float
        lengte van de waterlopen (m)
    Bin : float
        bodembreedte (m)
    c0 : float
        slootbodemweerstand (d)
    
    Returns
    -------
    float
        Lekconductance (m2/d)
    
    """
    if li>0.001 and Bin>0.001 and A>0.001:
        Bcor = max(Bin, 0.001)
        L = A / li - Bcor
        y = c1 + H0 / kv
        labdaL = np.sqrt(y * kh * H0)
        x = 0.5 * L / labdaL
        FL = x * ctnh(x)
        labdaB = np.sqrt(y * kh * H0 * c0 / (y + c0))
        x = 0.5 * Bcor / labdaB
        FB = x * ctnh(x)
        CL = (c0 + y) * FL + (c0 * L / Bcor) * FB
        CB = (c1 + c0 + H0 / kv) / (CL - c0 * L / Bcor) * CL
        Crad = max(0., L / (np.pi * np.sqrt(kv * kh)) * np.log(4 * H0 / (np.pi * Bcor)))
        pSl = Bcor * li / A
        Wp = 1 / ((1. - pSl) / CL + pSl / CB) + Crad - c1
        return A/Wp
    else:
        return 0.
    
def lek_conductance_array(A, H0, kv, kh, c1, li, Bin, c0):
    """Calculates de lekconductance according to De Lange
    
    Parameters
    ----------
    A : float, numpy.ndarray or pandas.Series
        celoppervlak (m2)
    H0 : float, numpy.ndarray or pandas.Series
        doorstroomde dikte (m)
    kv : float, numpy.ndarray or pandas.Series
        verticale doorlotendheid (m/d)
    kh : float, numpy.ndarray or pandas.Series
        horizontale doorlatendheid (m/d)
    c1 : float, numpy.ndarray or pandas.Series
        deklaagweerstand (d)
    li : numpy.ndarray or pandas.Series
        lengte van de waterlopen (m)
    Bin : float, numpy.ndarray or pandas.Series
        bodembreedte (m)
    c0 : float, numpy.ndarray or pandas.Series
        slootbodemweerstand (d)
    
    Returns
    -------
    numpy.ndarray or pandas.Series
        Lekconductance (m2/d)
    
    """
    
    Bcor = Bin.copy()
    Bcor[Bcor<0.001] = 0.001
    L = A / li - Bcor
    y = c1 + H0 / kv
    labdaL = np.sqrt(y * kh * H0)
    x = 0.5 * L / labdaL
    FL = x * ctnh(x)
    labdaB = np.sqrt(y * kh * H0 * c0 / (y + c0))
    x = 0.5 * Bcor / labdaB
    FB = x * ctnh(x)
    CL = (c0 + y) * FL + (c0 * L / Bcor) * FB
    CB = (c1 + c0 + H0 / kv) / (CL - c0 * L / Bcor) * CL
    Crad = L / (np.pi * np.sqrt(kv * kh)) * np.log(4 * H0 / (np.pi * Bcor))
    Crad[Crad<0.0] = 0.0
    pSl = Bcor * li / A
    Wp = 1 / ((1. - pSl) / CL + pSl / CB) + Crad - c1
    c = A/Wp
    mask = (li<=0.001) | (Bin<=0.001) | (A<=0.001)
    c[mask] = 0.
    return c
              
def ctnh(x):
    return 1 / np.tanh(x)

def get_head_time_series(rf, model_ws, x, y, lay, kind='head'):
    """"Get the series of the head at a certain x and y"""
    fname = os.path.join(model_ws,rf.data['OUTPUTDIRECTORY'],kind,'{}_*'.format(kind))
    head = imod.idf.open(fname)
    return pd.Series(head.sel(layer=lay, x=x, y=y, method='nearest'),
                     index=head.time.values)
