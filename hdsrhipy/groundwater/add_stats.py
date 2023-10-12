"""
Script to calculate the following validation statistics:
    For observations, model (_MOD), and model error (_ERROR):
    1) GHG: mean highest groundwater level
    2) GLG: mean lowest groundwater level
    3) DYNAMIEK: GHG - GLG
    4) NGXG: number of hydrological years on which GHG/GLG/DYNAMIEK are based
    5) MEAN: average groundwater level

    For residuals and monthly residuals (RESIDUAL_MONTH):
    6) RESIDUAL_ME: mean error of residuals
    7) RESIDUAL_MAE: mean absolute error of residuals
    8) RESIDUAL_RMSE: root mean squared error of residuals
    9) RESNORM_MAE: mean absolute error of normalized residuals
    10) RESNORM_RMSE: root mean squared error of normalized residuals

    For 2018 and 2019:
    11) DIFFHEAD: measured difference between max and min head during drought
    12) GRADHEAD: measured gradient between max and min head during drought
    13) HEADERROR: difference between modeled and measured diffhead
    15) MINTIMEERROR: difference between modeled and measured time of min head
    16) GRADERROR: difference between modeled and measured gradhead
    17) MINDATE: date of measured minimum head

    Gradient errors:
    18) RELERROR_RISE_MONTHLY: relative error for positive monthly head
        gradient
    19) RELERROR_DECLINE_MONTHLY: relative error for negative monthly head
        gradient
    20) RELERROR_RISE_WEEKLY: relative error for positive weekly head gradient
    21) RELERROR_DECLINE_WEEKLY: relative error for negatie weekly head
        gradient
    22) RELERROR_RISE_DAILY: relative error for positive daily head gradient
    23) RELERROR_DECLINE_DAILY: relative error for negative weekly head
        gradient

Note:
This script is developed based on the data format of input files
used and files generated during validation of LHM 4.1.
Hence, if the user uses different data with different formatting,
all necessary amendments first need to be implemented by the user.

"""

import pandas as pd
import numpy as np
import imod
import datetime
from hdsrhipy.groundwater.stats import getGXG #, outlier_test
import os
import xarray as xr
from pathlib import Path
#import configparser
#from readArgs import readArgs


def testDrought(df, year):
    """
    Method to calculate recession statistics
    for drought period during 2018 and 2019

    Parameters
    ----------
    df : pandas DataFrame
       time series data
    year : integer
       year to evaluate (2018 or 2019)

    Returns
    -------
    list of statistics
        diffHeadMeas: measured difference between max and min head
        gradMeas: measured gradient between max and min head
        headError: difference between modeled and measured diffhead
        minTimeError: difference between modeled and measured time of min head
        gradError: difference between modeled and measured gradhead
        minDateMeas: date of measured minimum head

    Notes
    -----
    This method has been developed for the validation of LHM 4.1.
    If the user uses other data or requires other time periods,
    all necessary amendments first need to be implemented by the user.
    """

    # get slice of time series for selected drought period
    sDate = datetime.datetime(year, 4, 1)
    eDate = datetime.datetime(year+1, 1, 1)
    if year == 2018:
        minEndDate = datetime.datetime(year,12,1)
    else:
        minEndDate = datetime.datetime(year,10,1)
    dfYear = df.loc[sDate:eDate]

    # initialize statistics
    gradError = np.NaN
    headError = np.NaN
    diffHeadMeas = np.NaN
    gradMeas = np.NaN
    minTimeError = np.NaN
    minDateMeas = np.NaN

    # period from which maximum head is taken
    maxPeriod = dfYear[:datetime.datetime(year, 6, 1)]
    # period from which minimum head is taken
    minPeriod = dfYear[datetime.datetime(year, 7, 1):]

    # only continue if:
    # - number of observations is more than 126 (70% of 180)
    # - last observation is after minEndDate
    # - maxPeriod and maxPeriod contain observations
    if (dfYear.shape[0] > 126 and dfYear.index[-1] >= minEndDate
        and maxPeriod.shape[0] > 0 and minPeriod.shape[0] > 0):
        # check whether at least 90% of the weeks has one or more observations
        dfWeek = dfYear.asfreq('W').dropna()
        nWeeks = (dfYear.index[-1] - dfYear.index[0]).days/7
        if dfWeek.shape[0] > 0.9*nWeeks:
            # get maximum head
            maxHead = maxPeriod.max()
            # get minimum head
            minHead = minPeriod.min()
            # get date of maximum head
            maxDate = maxPeriod.idxmax()
            # get date of minimum head
            minDate = minPeriod.idxmin()
            # check if heads are rising again after minimum,
            # to make sure that we have found the absolute minimum
            # if so, calculate statistics
            if ((dfYear['head_measured(m)'].values[-1] >
                 minHead['head_measured(m)'] + 0.05)
                and (dfYear['head_modeled(m)'].values[-1] >
                     minHead['head_modeled(m)'] + 0.05)
                and (minDate['head_measured(m)'] < minEndDate)):
                diffHead = maxHead - minHead
                diffTime = (minDate - maxDate).dt.days
                grad = diffHead / diffTime

                # use head_measured instead of head_measured_clean since
                # outlier detection may have resulted in removing one or more
                # extreme low but valid heads
                minTimeError = (minDate['head_modeled(m)']
                             - minDate['head_measured(m)']).days
                headError = round((diffHead['head_modeled(m)']
                             - diffHead['head_measured(m)']), 2)
                gradError = (grad['head_modeled(m)'] /
                             grad['head_measured(m)']) - 1
                diffHeadMeas = diffHead['head_measured(m)']
                gradMeas = grad['head_measured(m)']
                minDateMeas = minDate['head_measured(m)']

    return [diffHeadMeas, gradMeas, headError, minTimeError,
            gradError, minDateMeas]


def getMeanStats(ser):
    """
    Method to calculate statistics for series

    Parameters
    ----------
    ser : pandas Series
       time series to be analyzed

    Returns
    -------
    ret : dictionary
       statistics
    """
    ret = {}
    ret['me'] = np.round(ser.mean(), 2)
    ret['mae'] = np.round(ser.abs().mean(), 2)
    ret['rmse'] = np.round(np.sqrt(ser.pow(2).mean()), 2)
    norm = np.round(ser - ser.mean(), 2)
    ret['mae_norm'] = np.round(norm.abs().mean(), 2)
    ret['rmse_norm'] = np.round(np.sqrt(norm.pow(2).mean()), 2)

    return ret

def relativeError(modeled, measured):
    """
    Method to calculate statistics for series

    Parameters
    ----------
    modeled : pandas Series
       time series of modeled heads
    measured: pandas Series
       time series of measured heads

    Returns
    -------
    ret : list
       relative error for rising and declining heads
    """
    rise = (modeled[modeled > 0].mean()
            / measured[measured > 0].mean()) - 1
    decline = (modeled[modeled < 0].mean()
               / measured[measured < 0].mean()) - 1

    ret = np.round([rise, decline], 2)
    return ret


def getGradient(data):
    """
    Method to calculate rising and declining gradient
    on monthly, weekly, and daily basis over hydrological years
    Apr 2011 - March 2019

    Parameters
    ----------
    data : pandas DataFrame
       time series data of modeled and measured heads

    Returns
    -------
    ret : list
       relative error for rising and declining heads

    Notes
    -----
    This method has been developed for the validation of LHM 4.1.
    If the user uses other data or requires other time periods,
    all necessary amendments first need to be implemented by the user.
    """
    ret = []
    df = data.copy()
    # slice data for required period
    date_start = datetime.datetime(year=2011, month=4, day=1)
    date_end = datetime.datetime(year=2019, month=4, day=1)
    df = df.loc[date_start:date_end, :]
    df['time'] = df.index
    df['dt'] = (df['time'] - df['time'].shift(1)).dt.days
    df['dh_measured(m)'] = df['head_measured(m)'].diff()
    df['dh_modeled(m)'] = df['head_modeled(m)'].diff()
    grad = df[['dh_measured(m)', 'dh_modeled(m)']].div(df['dt'], axis=0)
    gradMonth = grad.resample('M').mean()
    ret.extend(relativeError(gradMonth['dh_modeled(m)'],
                             gradMonth['dh_measured(m)']))
    gradWeek = grad.resample('W').mean()
    ret.extend(relativeError(gradWeek['dh_modeled(m)'],
                             gradWeek['dh_measured(m)']))
    grad = grad.asfreq('D')
    grad = grad.rolling(3, center=True).mean().dropna()
    # only calculate relative gradient error on daily basis if minimum of
    # 2 years of data are available
    if grad.shape[0] > 730:
        ret.extend(relativeError(grad['dh_modeled(m)'],
                                 grad['dh_measured(m)']))
    else:
        ret.extend([np.NaN, np.NaN])

    return ret


def addStats(tsPath, outPath, ncData, gi, x, y, wellId, lay):
    """
    Method to read modeled heads from netCDF files and
    calculate validation statistics.

    Parameters
    ----------
    tsPath : string
        path of time series data (measurements)
    outPath : string
        output directory
    ncData: list of xarray
        list of head data per model layer
    gi : dictionary
        geographical information for grid
    x : float
        x-coordinate of screen
    y : float
        y-coordinate of screen
    wellId : string
        well identification
    lay : integer
        model layer number

    Returns
    -------
    tuple with validation statistics

    Notes
    -----
    The order of the items in the tuple should be identical to those
    in outList as defined in the main program
    """
    txtFile = os.path.join(tsPath, (str(wellId) + '.txt'))
    #print (txtFile)
    if os.path.exists(txtFile):
        # read measured time series data
        df = imod.ipf.read_associated(txtFile, {'delim_whitespace': False})
        dateRange = ncData[int(lay) - 1].time.values
        # read modeled series from xarray
        xcol = int((x-gi['xmin']) / gi['xcellsize'])
        yrow = int((y-gi['ymax']) / gi['ycellsize'])
        values = ncData[int(lay) - 1].isel(x=xcol,y=yrow).values
        modeledSeries = pd.DataFrame(
            values,
            index=dateRange,
            columns=['head_modeled(m)'])
        # merge modeled and measured series by index of measured series
        df = df.merge(modeledSeries, left_on='time', right_index=True,
                      how='left')
        # write merged dataframe to ipf
        if not os.path.exists(os.path.join(outPath, 'ipf')):
            os.makedirs(os.path.join(outPath, 'ipf'))
        outFile = os.path.join(outPath, 'ipf', wellId + '.txt')
        imod.ipf.write_assoc(outFile, df, nodata=-9999.)

        # set index for further processing
        df.set_index('time', inplace=True)

        # if earlier steps do not contain an outlier test, the user
        # may use the code below to do outlier test here, before
        # calculating statistics
        # z_GW = outlier_test(df['head_measured(m)'])
        # threshold = 3
        # df['head_measured_clean(m)'] = df.loc[:, 'head_measured(m)']
        # if z_GW[pd.notnull(z_GW)].shape[0] > 0:
        #     df.loc[(z_GW >= threshold), 'head_measured_clean(m)'] = np.NaN

        retList = []
        # calculate statistics
        gxgMod = getGXG(df['head_modeled(m)'], minyear=3, minval=20,
                        nearest=True, nearlim=4)
        gxgMod['dyn'] = np.round(gxgMod['ghg'] - gxgMod['glg'], 2)
        gxgMod['mean'] = df['head_modeled(m)'].mean()
        gxgClean = getGXG(df['head_measured_clean(m)'], minyear=3,
                          minval=20, nearest=True, nearlim=4)
        gxgClean['dyn'] = np.round(gxgClean['ghg'] - gxgClean['glg'], 2)
        gxgClean['mean'] = df['head_measured_clean(m)'].mean()

        # GHG, GLG, GVG, DYNAMIEK, NGXG, MEAN
        for gxg in ['ghg', 'glg', 'gvg', 'dyn', 'ngxg', 'mean']:
            retList.append(gxgClean[gxg])
        # GHG_MOD, GLG_MOD, GVG_MOD, DYNAMIEK_MOD, NGXG_MOD, MEAN_MOD
        for gxg in ['ghg', 'glg', 'gvg', 'dyn', 'ngxg', 'mean']:
            retList.append(gxgMod[gxg])
        # GHG_ERROR, GLG_ERROR, GVG_ERROR, DYNAMIEK_ERROR, MEAN_ERROR
        for gxg in ['ghg', 'glg', 'gvg', 'dyn', 'mean']:
            retList.append(round((gxgMod[gxg] - gxgClean[gxg]), 2))

        df['residual'] = (df['head_measured_clean(m)']
                          - df['head_modeled(m)'])
        meanStats = getMeanStats(df['residual'])
        # RESIDUAL_ME, RESIDUAL_MAE, RESIDUAL_RMSE, RESNORM_MAE, RESNORM_RMSE
        for item in ['me', 'mae', 'rmse', 'mae_norm', 'rmse_norm']:
            retList.append(meanStats[item])

        df['residual_month'] = df['residual'].resample('M').mean()
        meanStats = getMeanStats(df['residual_month'])
        # RESIDUAL_MONTH_ME, RESIDUAL_MONTH_MAE, RESIDUAL_MONTH_RMSE,
        # RESNORM_MONTH_MAE, RESNORM_MONTH_RMSE
        for item in ['me', 'mae', 'rmse', 'mae_norm', 'rmse_norm']:
            retList.append(meanStats[item])

        for year in [2018, 2019]:
            # DIFFHEAD, GRADHEAD, HEADERROR, MINTIMEERROR, GRADERROR, MINDATE
            retList.extend(testDrought(df[['head_measured(m)',
                                           'head_modeled(m)']],
                                       year))

        # RELERROR_RISE_MONTHLY, RELERROR_DECLINE_MONTHLY,
        # RELERROR_RISE_WEEKLY, RELERROR_DECLINE_WEEKLY,
        # RELERROR_RISE_DAILY, RELERROR_DECLINE_DAILY
        retList.extend(getGradient(df))

    else:
        retList = [np.NaN]*45

    return tuple(retList)


outList = ['GHG', 'GLG', 'GVG', 'DYNAMIEK', 'NGXG', 'MEAN',
           'GHG_MOD', 'GLG_MOD', 'GVG_MOD', 'DYNAMIEK_MOD', 'NGXG_MOD',
           'MEAN_MOD',
           'GHG_ERROR', 'GLG_ERROR', 'GVG_ERROR', 'DYNAMIEK_ERROR',
           'MEAN_ERROR',
           'RESIDUAL_ME', 'RESIDUAL_MAE', 'RESIDUAL_RMSE',
           'RESNORM_MAE', 'RESNORM_RMSE',
           'RESIDUAL_MONTH_ME', 'RESIDUAL_MONTH_MAE', 'RESIDUAL_MONTH_RMSE',
           'RESNORM_MONTH_MAE', 'RESNORM_MONTH_RMSE',
           '2018_DIFFHEAD', '2018_GRADHEAD', '2018_HEADERROR',
           '2018_MINTIMEERROR', '2018_GRADERROR', '2018_MINDATE',
           '2019_DIFFHEAD', '2019_GRADHEAD', '2019_HEADERROR',
           '2019_MINTIMEERROR', '2019_GRADERROR', '2019_MINDATE',
           'RELERROR_RISE_MONTHLY', 'RELERROR_DECLINE_MONTHLY',
           'RELERROR_RISE_WEEKLY', 'RELERROR_DECLINE_WEEKLY',
           'RELERROR_RISE_DAILY', 'RELERROR_DECLINE_DAILY'
           ]
