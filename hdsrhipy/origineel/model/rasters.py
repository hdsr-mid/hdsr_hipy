# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 12:31:23 2018

@author: Artesia
"""
from osgeo import osr
from osgeo import gdal
import numpy as np

def xyz2gdal_dataset(x, y, z, projection=28992, no_data_value=None,
                     data_type=gdal.GDT_Float32):
    """
    A method to turn numpy x,y and z values to a gdal-dataset
    x and y should be 1d arrays, describing the cell centers or the boundaries
    y can go up or down, this does not matter
    """
    # we create an in-memory raster
    mem_drv = gdal.GetDriverByName('MEM')
    outDs = mem_drv.Create('', z.shape[1], z.shape[0], 1, data_type)
    outBand = outDs.GetRasterBand(1)
    outBand.WriteArray(z, 0, 0)
    # flush data to disk, set the NoData value and calculate stats
    outBand.FlushCache()
    if no_data_value is not None:
        outBand.SetNoDataValue(no_data_value)
    # georeference the image and set the projection
    if len(x.shape) == 1 and len(y.shape) == 1:
        # x and y are vectors
        if x.shape[0] == z.shape[1] and y.shape[0] == z.shape[0]:
            # x and y are the cell-centers
            dx = np.unique(np.diff(x))
            dy = np.unique(np.diff(y))
            # calculate the cell-borders
            x = np.hstack((x - dx / 2, x[-1] + dx / 2))
            y = np.hstack((y - dy / 2, y[-1] + dy / 2))
            minx = x.min()
            maxy = y.max()
            dy = np.abs(dy)
            dx = np.abs(dx)
        elif x.shape[0] == z.shape[1] + 1 and y.shape[0] == z.shape[0] + 1:
            # x and y are the boundaries
            minx = x.min()
            maxy = y.max()
            dx = np.abs(np.unique(np.diff(x)))
            dy = np.abs(np.unique(np.diff(y)))
        else:
            ValueError('shape of x and y does not correspond to shape of z')
    else:
        ValueError('x and y should be supplied as vectors')
    if dx.size > 1 or dy.size > 1:
        ValueError('gridsize not constant')
    geo_transform = (minx, dx[0], 0, maxy, 0, -dy[0])
    outDs.SetGeoTransform(geo_transform)
    if projection is not None:
        if type(projection) == int:
            sr = osr.SpatialReference()
            sr.ImportFromEPSG(projection)
            projection = sr.ExportToWkt()
        outDs.SetProjection(projection)
    return outDs

def clip_dataset(ds,extent):
    projWin = [extent[0], extent[3], extent[1], extent[2]]
    return gdal.Translate('', ds,  format = 'MEM', projWin = projWin)

def reproject_dataset(ds, pixel_spacing=5000., epsg_from=28992,
                      epsg_to=28992, GDALResampleAlg=gdal.GRA_Average):
    """
    A sample function to reproject and resample a GDAL dataset from within
    Python. The idea here is to reproject from one system to another, as well
    as to change the pixel size. The procedure is slightly long-winded, but
    goes like this:

    1. Set up the two Spatial Reference systems.
    2. Open the original dataset, and get the geotransform
    3. Calculate bounds of new geotransform by projecting the UL corners
    4. Calculate the number of pixels with the new projection & spacing
    5. Create an in-memory raster dataset
    6. Perform the projection
    """
    # Define the UK OSNG, see <http://spatialreference.org/ref/epsg/27700/>
    sr_to = osr.SpatialReference()
    sr_to.ImportFromEPSG(epsg_to)
    sr_from = osr.SpatialReference()
    sr_from.ImportFromEPSG(epsg_from)
    tx = osr.CoordinateTransformation(sr_from, sr_to)
    # Up to here, all  the projection have been defined, as well as a
    # transformation from the from to the  to :)
    # We now open the dataset
    if type(ds) == gdal.Dataset:
        g = ds
    else:
        g = gdal.Open(ds)
    # Get the Geotransform vector
    geo_t = g.GetGeoTransform()
    x_size = g.RasterXSize  # Raster xsize
    y_size = g.RasterYSize  # Raster ysize
    # Work out the boundaries of the new dataset in the target projection
    (ulx, uly, ulz) = tx.TransformPoint(geo_t[0], geo_t[3])
    (lrx, lry, lrz) = tx.TransformPoint(geo_t[0] + geo_t[1] * x_size,
                                        geo_t[3] + geo_t[5] * y_size)
    # See how using 27700 and WGS84 introduces a z-value!
    # Now, we create an in-memory raster
    mem_drv = gdal.GetDriverByName('MEM')
    # The size of the raster is given the new projection and pixel spacing
    # Using the values we calculated above. Also, setting it to store one band
    # and to use Float32 data type.
    dest = mem_drv.Create('', int((lrx - ulx) / pixel_spacing),
                          int((uly - lry) / pixel_spacing), 1,
                          gdal.GDT_Float32)
    # Calculate the new geotransform
    new_geo = (ulx, pixel_spacing, geo_t[2],
               uly, geo_t[4], -pixel_spacing)
    # Set the geotransform
    dest.SetGeoTransform(new_geo)
    dest.SetProjection(sr_to.ExportToWkt())
    # Perform the projection/resampling
    gdal.ReprojectImage(g, dest,
                        sr_from.ExportToWkt(), sr_to.ExportToWkt(),
                        GDALResampleAlg)
    return dest

def get_values(ds):
    """Get the grid-values of the gdal dataset ds"""
    band = ds.GetRasterBand(1)
    H = band.ReadAsArray()
    if H.dtype=='float32':
        H[H == band.GetNoDataValue()] = np.nan
    return H

def get_xy(ds):
    """"Get the x and y coordinates of the borders of the grid of the gdal dataset ds"""
    nY = ds.RasterYSize
    nX = ds.RasterXSize
    gt = ds.GetGeoTransform()
    x = np.linspace(gt[0], gt[0] + nX * gt[1], nX + 1)
    y = np.linspace(gt[3], gt[3] + nY * gt[5], nY + 1)
    return x, y

def get_xy_mid(ds):
    """Get the x and y coordinates of the middle of the grid of the gdal dataset ds"""
    nY = ds.RasterYSize
    nX = ds.RasterXSize
    gt = ds.GetGeoTransform()
    x = np.linspace(gt[0]+gt[1]/2, gt[0] + nX * gt[1] - gt[1]/2, nX)
    y = np.linspace(gt[3]+gt[5]/2, gt[3] + nY * gt[5] - gt[5]/2, nY)
    return x, y