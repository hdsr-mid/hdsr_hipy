# -*- coding: utf-8 -*-
"""
This script selects for piezometers in dataset, KNMI stations
for precipitation and evaporation data based on location
and observation period of the piezometers.

Note:
This script is developed based on the data format of input files
used and files generated during validation of LHM 4.1.
Hence, if the user uses different data with different formatting,
all necessary amendments first need to be implemented by the user.

"""

import pastas as ps
from datetime import datetime
import pandas as pd
import numpy as np
from scipy import spatial
import os
import imod

class KnmiStationSelector(object):

    def __init__(self,path):
        """
        KnmiStationSelection class to select KNMI station based on closest
        distance of station to piezometer, requiring complete overlap of series

        Parameters
        ----------
        config : configparser object
            configuration settings
        root : string
            pathname for storing results
        """

        self.root = os.path.join(path, 'knmi')
        if not os.path.exists(self.root):
            os.makedirs(self.root)

        # start- and enddate of period to be considered
        self.startDate = datetime(2011, 1, 1).date()
        self.endDate = datetime(pd.to_datetime('now').year, 7, 31).date()

        # start- and enddate for downloading KNMI data
        self.startKnmi = datetime(1995, 1, 1)
        self.endKnmi = datetime(pd.to_datetime('now').year, 12, 31)

        # minimum length of input data before observation period
        # (required for initialization of time series model)
        self.mininit = 1

        self.knmiMemory = {
            'meteo': {
                'prec': {},
                'evap': {}
                },
            'prec': {}
            }
        self.initialMaxSearch = 1.e12
        self.maxsearch = self.initialMaxSearch

        self.meteoStns = pd.read_csv(os.path.join(self.root, 'meteostns.csv'), header=0)
        self.precStns = pd.read_csv(os.path.join(self.root,'precstns.csv'), header=0)


    def testSeries(self, inSeries):
        """
        Test input series obtained from online KNMI database

        Parameters
        ----------
        inSeries : pandas DataFrame
            input series obtained from KNMI

        Returns
        -------
        lengthtest : boolean
            true if input series overlaps output series
        nyear_init : float
            length of initialization period before output series starts (years)

        """
        lengthtest = False
        nyear_init = 0.
        if inSeries.shape[0] > 2:
            dates = inSeries.index.get_level_values(0)
            startdate = dates[0].date()
            enddate = dates[-1].date()
            if startdate <= self.startDate and enddate >= self.endDate:
                    lengthtest = True
                    nyear_init = (self.startDate - startdate).days / 365.25

        return lengthtest, nyear_init


    def readStnInfo(self, data, xyCrd):
        """
        Extract station information from data and sort by distance
        from well location

        Parameters
        ----------
        data : pandas DataFrame
            table with station information
        xyCrd : list
            xy coordinates of well location

        Returns
        -------
        istn : array
            identifiers for station
        dist : array
            distances to station
        stn : array
            station numbers

        """
        npxy = np.array((data.loc[:, ['X', 'Y']].astype(float).values))
        npstn = np.array((data.iloc[:, 1].astype(int).values))
        tree = spatial.KDTree(npxy)
        istn = tree.query(xyCrd, k=data.shape[0])
        dist = istn[0]
        stn = npstn[istn[1]]

        return istn, dist, stn


    def findStn(self, intype, stn_type, xyCrd):
        """
        Find station that satisfies criteria defined by user

        Parameters
        ----------
        intype : string
            type of input series (prec, evap, evpf)
        stn_type : string
            type of KNMI station (meteo, prec)
        xyCrd : list
            xy coordinates of well location


        Returns
        -------
        iStation : integer
            identification of KNMI station
        iDist : float
            distance of KNMI station from well location
        find_ok : boolean
            True if a station can be found that satisfies criteria

        """

        if stn_type == 'meteo':
            istn, dist, stn = self.readStnInfo(self.meteoStns, xyCrd)
            if intype == 'prec':
                knmiVars = 'RH'
            else:
                knmiVars = 'EV24'
        else:
            istn, dist, stn = self.readStnInfo(self.precStns, xyCrd)
            knmiVars = 'RD'


        # initialize while loop
        iDist = 0.
        lengthtest = False
        inityear = -1.
        i = 0

        find_ok = True

        while ((iDist <= self.maxsearch)
               and (i < len(stn))
               and (not(lengthtest) or inityear < self.mininit)
               ):

            # get data for next nearest station
            iStation = stn[i]
            iDist = dist[i]

            if iDist <= self.maxsearch:

                # test if data is available in memory
                isInMemory = False
                if stn_type == 'meteo':
                    if iStation in self.knmiMemory[stn_type][intype]:
                        inSeries = self.knmiMemory[stn_type][intype][iStation]
                        isInMemory = True
                else:
                    if iStation in self.knmiMemory[stn_type]:
                        inSeries = self.knmiMemory[stn_type][iStation]
                        isInMemory = True

                if not(isInMemory):
                    # download data from internet
                    knmi = ps.read.KnmiStation.download(stns=iStation,
                                                start=self.startKnmi,
                                                end=self.endKnmi,
                                                interval='daily',
                                                vars=knmiVars)
                    inSeries = knmi.data.dropna(subset=[knmiVars])

                    suffix = ''
                    if stn_type == 'meteo':
                        self.knmiMemory[stn_type][intype].update(
                            {iStation: inSeries}
                            )
                        ipfseries = inSeries.reset_index().drop('STN',
                                                                axis='columns')
                        if intype == 'prec':
                            ipfseries = ipfseries.rename(
                                columns={
                                    'YYYYMMDD': 'time',
                                    'RH': 'Calculated'
                                    }
                                )
                            suffix = '_RH'
                        else:
                            ipfseries = ipfseries.rename(
                                columns={
                                    'YYYYMMDD': 'time',
                                    'EV24': 'Calculated'
                                    }
                                )
                    else:
                        self.knmiMemory[stn_type].update({iStation: inSeries})
                        ipfseries = inSeries.reset_index().drop(['STN'], axis='columns')
                        ipfseries = ipfseries.rename(
                            columns={
                                'YYYYMMDD': 'time',
                                'RD': 'Calculated'
                                }
                            )
                        suffix = '_RD'
                    if intype == 'prec':
                        path = os.path.join(self.root, 'precipitation')
                    else:
                        path = os.path.join(self.root, 'evapotranspiration')

                    # save KNMI data to ipf file
                    ipfseries['Calculated'] = (
                        ipfseries['Calculated'] * 1000.
                        ).round(2)
                    path2 = os.path.join(path,
                                          str(iStation) + suffix + '.txt')
                    imod.ipf.write_assoc(path2,
                                          ipfseries,
                                          itype=1,
                                          nodata=-9999.)

                # test for length of series
                lengthtest, inityear = self.testSeries(inSeries)

                i += 1

            else:
                inSeries = None
                find_ok = False

        # if after while-loop the criteria are not met, set find_ok to False
        if inityear >= self.mininit:
            if not(lengthtest):
                find_ok = False
        else:
             find_ok = False

        return iStation, iDist, find_ok


    def getKnmiStation(self, xyCrd, intype):
        """
        Find KNMI station that satisfies criteria

        Parameters
        ----------
        xyCrd : list
            xy coordinates of well location
        intype : string
            'prec' or 'evap'

        Returns
        -------
        iStation : integer
            identification of KNMI station
        stnType : string
            'RH' refers to meteo station - precipitation
            'EV24' refers to meteo station - evapotranspiration
            'RD' refers to precipitation station

        """

        # initialize
        meteo_ok = True
        prec_ok = False

        if intype == 'prec':
            prec_ok = True
            iStation_prec, iDist_prec, prec_ok = self.findStn(intype,
                                                               'prec',
                                                               xyCrd)
        # only search for meteostation until distance is greater than
        # distance of selected precstation
        if prec_ok == True:
            self.maxsearch = iDist_prec

        iStation, iDist, meteo_ok = self.findStn(intype, 'meteo', xyCrd)

        if intype == 'prec':
            # evaluate distances of meteo and prec station and select closest
            use_prec = True
            if prec_ok:
                if meteo_ok:
                    if iDist_prec > iDist:
                        use_prec = False
            else:
                use_prec = False

            if use_prec:
                iDist = iDist_prec
                iStation = iStation_prec
                stnType = 'RD'
            else:
                stnType = 'RH'
        else:
            stnType = 'EV24'

        # reset maxsearch
        self.maxsearch = self.initialMaxSearch

        return int(iStation), stnType



