# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 10:41:17 2020

@author: Wilbert Berendrecht
"""

from netCDF4 import Dataset
import numpy as np

class nc(object):
    """nc object Class to read NetCDF file.

       Parameters
       ----------
       fName: string
          filename of netCDF
    """


    def __init__(self, fName, dataLayer):
        """Method to create/initiate a nc object.

        See Also
        --------
        :class:`nc`
        """
        self.fName = fName
        self.ncfile = Dataset(fName, 'r')
        self.gi()
        self.errCnt = 0
        self.ncArr = self.ncfile.variables[dataLayer]


    def gi(self):
        """Method to get geographic information of NetCDF file.
        """
        xvar = self.ncfile.variables['x'][:]
        self.xmin = xvar[0]
        self.xmax = xvar[-1]
        self.xncell = self.ncfile.dimensions['x'].size
        self.xcellsize = (self.xmax - self.xmin) / (self.xncell-1)
        self.xmin = self.xmin - 0.5*self.xcellsize
        self.xmax = self.xmax + 0.5*self.xcellsize


        yvar = self.ncfile.variables['y'][:]
        self.ymin = yvar[-1]
        self.ymax = yvar[0]
        self.yncell = self.ncfile.dimensions['y'].size
        self.ycellsize = (self.ymax - self.ymin) / (self.yncell-1)
        self.ymin = self.ymin - 0.5*self.ycellsize
        self.ymax = self.ymax + 0.5*self.ycellsize

    def xcol(self, x):
        """Method to calculate column number of grid based on x-coordinate.

           Parameters
           ----------
           x: float
              x-coordinate

           Returns
           -------
           col: integer
              column number in grid
        """
        return int((x-self.xmin)/self.xcellsize)

    def yrow(self, y):
        """Method to calculate row number of grid based on y-coordinate.

           Parameters
           ----------
           y: float
              y-coordinate

           Returns
           -------
           row: integer
              row number in grid
        """
        return int((self.ymax-y)/self.ycellsize)

    def getVarGroup(self, group, var):
        """Method to read 1-dimensional variable from group.

           Parameters
           ----------
           group: string
              group name
           var: string
              variable name

           Returns
           -------
           arr: 1-d array
              array of variable values
        """
        if group in self.ncfile.groups:
            if var in self.ncfile.groups[group].variables:
                return self.ncfile.groups[group].variables[var][:]
            else:
                raise Exception('NetCDF groep '+group+
                                ' bevat geen variabele met de naam '+var)
        else:
            raise Exception('NetCDF bevat geen groep met de naam '+group)
        return


    def getVarByAttribute(self, attr):
        """Method to read 1-dimensional variable selected by attribute.

           Parameters
           ----------
           attr: string
              attribute value

           Returns
           -------
           arr: 1-d array
              array of variable values
        """
        sel = self.ncfile.get_variables_by_attributes(type=attr)
        if len(sel)>0:
            if len(sel) == 1:
                return sel[0]
            else:
                raise Exception('NetCDF bevat meerdere variabelen met '
                            + 'attribuut ' + attr)

        else:
            raise Exception('NetCDF bevat geen attribuut met naam ' + attr)
        return


    def getVariable(self, var):
        """Method to read 1-dimensional variable.

           Parameters
           ----------
           group: string
              group name
           var: string
              variable name

           Returns
           -------
           arr: 1-d array
              array of variable values
        """
        if var in self.ncfile.variables:
            return self.ncfile.variables[var][:]
        else:
            raise Exception('NetCDF bevat geen variabele met de naam '+var)
        return


    def getVal(self, xcol, yrow, zlay=None):
        """Method to read value from 2 or 3-d array using
           column and row number.

           Parameters
           ----------
           xcol: integer
              column number
           yrow: integer
              row number
           zlay: integer (optional)
              layer index

           Returns
           -------
           val: float
              value from NetCDF array
        """
        if (xcol>=0 and xcol<self.xncell and yrow>=0 and yrow<self.yncell):
            if len(self.ncArr.shape) == 3:
                if  zlay is None:
                    return self.ncArr[:, yrow, xcol]
                else:
                    return self.ncArr[zlay, yrow, xcol]
            elif len(self.ncArr.shape) == 2:
                return self.ncArr[yrow, xcol]
            else:
                return self.ncArr
        else:
            return np.ma.array(np.NaN, mask=True)

    def getValxy(self, x, y):
        """Method to read value from 2 or 3-d variable using
           x and y coordinates.

           Parameters
           ----------
           x: float
              x coordinate
           y: float
              y coordinate

           Returns
           -------
           val: float
              value from NetCDF array
        """
        if (x>=self.xmin and x<=self.xmax and y>=self.ymin and y<=self.ymax):
            yr = self.yrow(y)
            xc = self.xcol(x)
            return self.getVal(xc, yr)
        else:
            return np.ma.array(np.NaN, mask=True)


    def getValxyz(self, x, y, z):
        """Method to read value from 3-d variable using x, y, z coordinates.

           Parameters
           ----------
           x: float
              x coordinate
           y: float
              y coordinate
           z: int
              z index


           Returns
           -------
           val: float
              value from NetCDF array
        """
        yr = self.yrow(y)
        xc = self.xcol(x)
        return self.getVal(xc, yr, z)

    def getValxyzFast(self, x, y, z):
        """Method to read value from 3-d variable using x, y, z coordinates.

           Parameters
           ----------
           x: float
              x coordinate
           y: float
              y coordinate
           z: int
              z index


           Returns
           -------
           val: float
              value from NetCDF array
        """
        yr = self.yrow(y)
        xc = self.xcol(x)
        return self.ncArr[z, yr, xc]
