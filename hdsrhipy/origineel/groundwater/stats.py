"""
Methods to calculate GXG and to detect outliers
"""

import pandas as pd
import numpy as np


def getGXG(ts, minyear=4, minval=20, nearest=False, nearlim=4):
    """
    Calculate GXG (GHG, GVG, GLG) for time series.

    Parameters
    ----------
    ts : pandas Series
        time series
    minyear : integer
        minimum number of years required for calculation of statistics
        (default = 4)
    minval : integer
        minimum number of values within a year required for including
        that year in calculation of statistic
    nearest : boolean
        use nearest value for 14th / 28th in GXG calculation (default: False)
    nearlim : integer
        maximum interval for `nearest` on both sides (default = 4)

    Returns
    -------
    gxgdf : pandas DataFrame
        GHG, GLG, RHG, RLG, and for each statistic the number of years
        on which it is based

    Notes
    -----
    The GHG and GLG (mean highest resp. lowest groundwater level) are
    calculated as follows:

    The time series is resampled to a series with measurements at the
    14th and 28th of each month. If `nearest` is True a bin of +/- `nearlim`
    is used to project the nearest observation to the 14th resp. 28th.
    Else, if `nearest` is False, simply the value at the 14th and 28th
    are taken.
    For each hydrological year (1 April - 31 March) the three highest resp.
    lowest values of the resampled series are averaged.
    The user may define a minimum number of values for each year
    (default `minval` =20). If in one year the number of values is less than
    `minval`, that year will be skipped. The average of the annual values
    over a period of at least `minyear` is used as the GHG resp. GLG.


    """
    if nearest:
        ts_int = ts.asfreq('D').interpolate(method='nearest',
                                            limit=nearlim,
                                            limit_direction='both')
        ts14_28 = ts_int[(ts_int.index.day == 14) |
                         (ts_int.index.day == 28)].dropna()
    else:
        ts14_28 = ts[(ts.index.day == 14) | (ts.index.day == 28)]

    # calculate gvg
    # group series for each year
    grouped = ts14_28.groupby([lambda x: x.year])
    gvgyear = []
    for name, group in grouped:
        grpidx = group.index
        gvg = group[((grpidx.month == 3) |
                     ((grpidx.month == 4) & (grpidx.day == 14)))
                    ].mean()
        try:
            gvgtest = gvg.mean()
        except:
            gvgtest = gvg
        if not np.isnan(gvgtest):
            gvgyear.append(gvg)

    # shift date 3 months to get 'hydrological year' April-March
    tsh = ts.copy()
    tsh.index = tsh.index + pd.DateOffset(months=-3)
    ts14_28.index = ts14_28.index + pd.DateOffset(months=-3)

    # calculate ghg and glg
    # group series for each year
    grouped = ts14_28.groupby([lambda x: x.year])
    glgyear = []
    ghgyear = []
    for name, group in grouped:
        nval = group.shape[0]
        if nval >= minval:
            yearsort = np.sort(group.values, axis=0)
            # get 3 lowest in year
            glgyear.append(yearsort[:3].mean())
            # get 3 highest in year
            ghgyear.append(yearsort[-3:].mean())

    # if enough years, calculate mean over all years
    gxg = {'ghg': np.NaN,
           'glg': np.NaN,
           'gvg': np.NaN,
           'ngxg': 0
           }
    if len(glgyear) >= minyear:
        gxg['ngxg'] = len(glgyear)
        gxg['glg'] = round(np.array(glgyear).mean(), 2)
        gxg['ghg'] = round(np.array(ghgyear).mean(), 2)
    if len(gvgyear) >= minyear:
        gxg['gvg'] = round(np.array(gvgyear).mean(), 2)

    return gxg


def outlier_test(ts):
    """
    Test for outliers based on the modified Z-score test
    (Iglewicz and Hoaglin, 1993).
	In order to account for varying observation frequency,
    the median of the time series is calculated by first
    resampling to the median of each month.

    Parameters
    ----------
    ts : pandas Series with date index
        series with output data

    Returns
    -------
    olv : pandas Series with date index
        modified Z-score


    Notes
    -----
    The modified Z-Score (:math:`M_i`) is computed as

    .. math::
        M_i = \\frac{0.6745(x_i-\\tilde{x})}{MAD}

    where :math:`E(MAD)=0.675\\sigma` for large normal data
    and :math:`\\tilde{x}` is the sample median.

    References
    ----------
    Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect
    and Handle Outliers", The ASQC Basic References in
    Quality Control: Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.
    """

    diff = np.fabs(ts.values - ts.resample('M').median().median())
    mad = np.median(diff)
    if not(mad == 0.):
        modified_z_score = 0.6745 * diff / mad
    else:
        modified_z_score = diff*np.NaN
    olv = pd.Series(modified_z_score,index=ts.index)

    return olv.replace([np.inf, -np.inf], np.nan)
