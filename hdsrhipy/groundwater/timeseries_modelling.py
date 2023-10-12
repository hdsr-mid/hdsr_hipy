"""
This script reads and pre-processes head, precipitation and evapotranspiration
data in ipf format as inputs using iMOD package and performs time series
analysis using PASTAS package and produces the results in png, ipf and xlsx
formats.

Initial version developed by Deltares and updated for LHM 4.1 validation
by adding:
    - updated outlier test
    - step test
    - trend test
    - test for drying up of screen (groundwater level below bottom of screen)
    - added option to read KNMI data by KNMI station;
      requires columns PREC_STN, PREC_STNTYPE, EVAP_STN
      (created with 03_get_knmi.py)

Note:
This script is developed based on the data format of input files
used during validation of LHM 4.1. Hence, if the user uses different
data with different formatting, all necessary amendments
first need to be implemented by the user.

"""

import pandas as pd
import numpy as np
import imod
import datetime
import pastas as ps
from scipy import stats
from tqdm.auto import tqdm
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from sklearn.linear_model import LinearRegression
import warnings
import os
import configparser
from dateutil.relativedelta import relativedelta
import pymannkendall as mk
from scipy.stats import kruskal

warnings.filterwarnings("ignore")


def testDryingUp(ts):
    """
    Test for drying up of monitoring well.

    Parameters
    ----------
    ts : pandas 1d series with date index
        series with time series data

    Returns
    -------
    test : int
       1 if True, 0 if False
    """
    perc = ts.quantile([0.01, 0.10, 0.99])
    prange = min(0.05*(perc[0.99] - perc[0.01]), 0.05)
    test = (perc[0.10] <= perc[0.01] + prange).astype(int)

    return test


def kwTest(ts):
    """
    Use Kruskal-Wallis test to check if the series portrays seasonal effects.

    Parameters
    ----------
    ts : pandas 1d series with date index
        series with time series data

    Returns
    -------
    pval : p-value
        float
    """
    ts.name = "meas"
    freq = pd.infer_freq(ts.index)
    if freq == "3M":
        nper = 4
    else:
        nper = 12
    ts_test = ts.copy()
    seas_slope, seas_intercept = mk.seasonal_sens_slope(ts_test.to_numpy())
    seas_slope = seas_slope / nper
    ts_test = ts_test.to_frame()
    ts_test.loc[:,"t"] = np.arange(len(ts_test))
    slope_b = (ts_test["meas"] - seas_slope * ts_test["t"]).median()
    ts_slope = seas_slope * ts_test["t"] + slope_b
    ts_corr = ts_test["meas"] - ts_slope
    freq = pd.infer_freq(ts.index)
    if freq != "3M":
        ts_corr = ts_corr.resample("3M").median()
    qrt = ts_corr.index.quarter
    stat, pval = kruskal(ts_corr[qrt==1].values,
            ts_corr[qrt==2].values,
            ts_corr[qrt==3].values,
            ts_corr[qrt==4].values)

    return pval


def trendTest(ts):
    """
    Test for linear trend in time series. Time series is resample to
    appropriate frequency preventing too much missing values.
    A seasonality test is used to decide which trend test is applied.
    This method doesn't use a specific test for autocorrelation. Instead it
    simply assumes autocorrelation for monthly series and no autocorrelation
    for quarterly series.

    Parameters
    ----------
    ts : pandas 1d series with date index
        series with time series data

    Returns
    -------
    slope : float
        estimated trend slope
    trend : int
        1 (if trend is present) or 0 (if the trend is absence)
    """
    try:
        # first resample series to appropriate frequency (<30% missing values)
        if ts.shape[0] > 20:
            # resample to monthly values
            tsFiltered = ts.copy()
            tsResample = tsFiltered.resample("M",
                                             label="left",
                                             loffset="14D").median()
            nper = 12
            if (tsResample.count() < 0.7 * tsResample.shape[0]):
                # if series has more than 30% missing values,
                # then change frequency to quarterly
                tsResample = tsFiltered.resample("3M",
                                                 label="left",
                                                 loffset="7W").median()
                nper = 4
                if (tsResample.count() < 0.7 * tsResample.shape[0]):
                    # if series still has more than 30% missing values,
                    # then no trend will be estimated
                    tsResample = None
        else:
            tsResample = None

        if tsResample is not None and tsResample.dropna().shape[0]>20:
            # perform seasonality test
            pval_seasonality = kwTest(tsResample)
            if pval_seasonality > 0.05:
                # no seasonality detected
                if nper == 12:
                    # assume autocorrelation
                    test1= mk.hamed_rao_modification_test(tsResample, lag=3)
                    test2 = mk.yue_wang_modification_test(tsResample, lag=1)
                    slope = test1.slope * nper
                    if (test1.h and test2.h):
                        trend = 1
                    else:
                        trend = 0
                else:
                    # assume no autocorrelation
                    test = mk.original_test(tsResample)
                    slope = test.slope
                    trend = int(test.h)
            else:
                # seasonality detected
                if nper == 12:
                    # assume autocorrelation
                    test = mk.correlated_seasonal_test(tsResample, period=nper)
                else:
                    # assume no autocorrelation
                    test = mk.seasonal_test(tsResample, period=nper)
                slope = test.slope
                trend = int(test.h)
        else:
            slope = np.NaN
            trend = np.NaN
            pval_seasonality = np.NaN
    except:
        slope = np.NaN
        trend = np.NaN
        pval_seasonality = np.NaN

    return slope, trend

def step(ts):
    """
    Test for step in time series

    Parameters
    ----------
    ts : pandas 1d series with date index
        series with output data

    Returns
    -------
    minPval : float
        two-tailed p-value
    maxTtest : float
        t-test statistic
    relStep : float
        difference between means relative to maximum standard deviation
    dateStep : datetime
        date in which step occurs
    """

    startdate = ts.index.get_level_values(0)[0]
    enddate = ts.index.get_level_values(0)[-1]
    tsResample = ts.resample("M").median().dropna()

    minPval = np.nan
    maxTtest = np.nan
    relStep = np.nan
    dateStep = np.nan

    splitStartDate = startdate + relativedelta(years=2)
    splitEndDate = enddate - relativedelta(years=2)
    if splitEndDate > splitStartDate:
        startYear = splitStartDate.year
        startMonth = splitStartDate.month
        endYear = splitEndDate.year
        endMonth = splitEndDate.month
        years = range(startYear, endYear + 1)

        startMonths = range(startMonth, 13)
        endMonths = range(1, endMonth)

        dates = []
        for m in startMonths:
            dates.append([startYear, m])
        for j in years[1:-1]:
            for m in range(1, 13):
                dates.append([j, m])
        for m in endMonths:
            dates.append([endYear, m])
        ttest = []
        pvals = []
        diffmeans = []
        for j in dates:
            splitdate = pd.datetime(int(j[0]), j[1], 1)
            ts1 = tsResample[:splitdate]
            ts2 = tsResample[splitdate:]
            if ts1.shape[0] >= 12 and ts2.shape[0] >= 12:
                res = stats.ttest_ind(ts1, ts2, axis=0, equal_var=True)
                pvals.append(res.pvalue)
                ttest.append(abs(res.statistic))
                diffmeans.append(abs(ts1.mean()-ts2.mean())
                                 / max(ts1.std(), ts2.std()))
            else:
                pvals.append(np.nan)
                ttest.append(np.nan)
                diffmeans.append(np.nan)
        if len(dates) > 0 and np.count_nonzero(~np.isnan(pvals)) > 0:
            idx = np.nanargmin(pvals)
            minPval = pvals[idx]
            maxTtest = ttest[idx]
            relStep = diffmeans[idx]
            dateStep = pd.datetime(dates[idx][0], dates[idx][1], 1).date()

    return minPval, maxTtest, relStep, dateStep

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
    """

    diff = np.fabs(ts.values - ts.resample("M").median().median())
    mad = np.median(diff)
    if not(mad == 0.):
        modified_z_score = 0.6745 * diff / mad
    else:
        modified_z_score = diff*np.NaN
    olv = pd.Series(modified_z_score,index=ts.index)

    return olv.replace([np.inf, -np.inf], np.nan)


def get_knmi_data(knmiData, df, path):
    precDir = path / "knmi" / "precipitation"
    evapDir = path / "knmi" / "evapotranspiration"

    # load precipitation data
    precType = df["PREC_STNTYPE"].values[0]
    stn = int(df["PREC_STN"].values[0])
    fName = os.path.join(precDir, str(stn)+"_"+precType+".txt")
    if (stn in knmiData[precType]):
        df_P = knmiData[precType][stn].copy()
    elif os.path.isfile(fName):
        df_P = imod.ipf.read_associated(fName)
        df_P = (df_P.rename(columns={
                    "time": "Date",
                    "calculated": "Rainfall(mm)"
                })
        )
    else:
        print ("downloading KNMI precipitation data for station", stn)
        knmiDownload = ps.read.KnmiStation.download(
            stns=stn,
            start=datetime.datetime(1995, 1, 1),
            end=datetime.datetime(2021, 1, 1),
            interval="daily",
            vars=precType)
        df_P = knmiDownload.data
        df_P = (df_P.reset_index()
                .rename(columns={
                    "YYYYMMDD": "Date",
                    precType: "Rainfall(mm)"
                })
        )
        df_P["Rainfall(mm)"] = 1000 * df_P["Rainfall(mm)"]
    df_P = df_P[["Date", "Rainfall(mm)"]]
    df_P["Date"] = pd.to_datetime(df_P["Date"].dt.date)
    knmiData[precType][stn] = df_P.copy()

    # load evapotranspiration data
    stn = int(df["EVAP_STN"].values[0])
    fName = os.path.join(evapDir, str(stn)+".txt")
    if (stn in knmiData["EV24"]):
        df_E = knmiData["EV24"][stn].copy()
    elif os.path.isfile(fName):
        df_E = imod.ipf.read_associated(fName)
        df_E = (df_E.rename(columns={
                    "time": "Date",
                    "calculated": "evapotranspiration(mm)"
                })
        )
    else:
        print ("downloading KNMI evapotranspiration data for station",
               stn)
        knmiDownload = ps.read.KnmiStation.download(
            stns=stn,
            start=datetime.datetime(1995, 1, 1),
            end=datetime.datetime(2021, 1, 1),
            interval="daily",
            vars="EV24")
        df_E = knmiDownload.data
        df_E = (df_E
                .reset_index()
                .rename(columns={
                    "YYYYMMDD": "Date",
                    "EV24": "evapotranspiration(mm)"
                })
        )
        df_E["evapotranspiration(mm)"] = (1000 *
            df_E["evapotranspiration(mm)"])

    df_E["Date"] = pd.to_datetime(df_E["Date"].dt.date)
    df_E = df_E[["Date", "evapotranspiration(mm)"]]
    knmiData["EV24"][stn] = df_E.copy()

    return df_P, df_E, knmiData


def pre_process_ipf(path, ipf_head, knmi, ipf_precipitation, ipf_evapotranspiration):
    """
    This function performs the followings:
        1)	reads ipf head data and creates GW info data frame
        2)	reads precipitation data from ipf
        3)	reads evapotranspiration data ipf
        4)	selects/renames columns
        5)	drops duplicates
        6)	resamples head data
        7)	merges precipitation, evapotranspiration and head data
            in one data frame.

    Parameters
    ----------
    path : Path object
        main folder contains input data
    ipf_head: string
        ipf file name for head data
    ipf_precipitation: string
        ipf file name for precipitation data
    ipf_evapotranspiration: string
        ipf file name for evapotranspiration data


    Returns
    -------
    df_merge_P_E_GW: pandas DataFrame
        merged precipitation, evapotranspiration and head data
    df_GW_info: pandas DataFrame
        GW info (e.g. ID,X,Y,etc)

    Note
    ----
    This function is based on the example attached to this script.
    If user is using other ipf files with different formatting,
    this function needs to be amended accordingly.

    """

    # 1) read ipf head data and creates GW info data frame
    print("reading head data : started")
    df_GW = imod.ipf.read(ipf_head, kwargs={"delim_whitespace": False},)
    df_GW_info = (
        df_GW[
            [
                "X_RD_CRD",
                "Y_RD_CRD",
                "Buis",
                "FILTER_TOP_NAP",
                "FILTER_ONDER_NAP",
                "EIGENAAR",
            ]
        ]
        .drop_duplicates(subset="Buis", keep="first")
        .reset_index(drop=True)
        .rename(columns={"Buis": "ID"})
    )

    print("reading head data : done")

    if knmi == "piezo":
        # 2) read precipitation data from ipf
        print("reading precipitation data : started")

        df_P = imod.ipf.read(
            path / "precipitation" / ipf_precipitation,
            kwargs={"delim_whitespace": False},
        )

        print("reading precipitation data : done")

        # 3) read evapotranspiration data from ipf
        print("reading evapotranspiration data : started")

        df_E = imod.ipf.read(
            path / "evapotranspiration" / ipf_evapotranspiration,
            kwargs={"delim_whitespace": False},
        )

        print("reading evapotranspiration data : done")

    if knmi == "piezo":
        # 4-6) Select/rename columns, drop duplicates, resample
        df_P = (
            df_P[["Date", "Calculated", "Buis"]]
            .rename(columns={"Buis": "ID", "Calculated": "Rainfall(mm)"})
            .drop_duplicates(keep=False)
        )
        df_E = (
            df_E[["Date", "Calculated", "Buis"]]
            .rename(columns={"Buis": "ID",
                             "Calculated": "evapotranspiration(mm)"})
            .drop_duplicates(keep=False)
        )
        df_GW = (
            df_GW[["time", "head", "Buis"]]
            .rename(columns={"Buis": "ID", "head": "head_measured(m)",
                             "time": "Date"})
            .drop_duplicates(keep=False)
            .set_index("Date")
            .groupby("ID")
            .resample("d")
            .mean()
            .reset_index()
        )
    else:
        df_GW = (
            df_GW[["time", "head", "Buis",
                   "PREC_STN", "PREC_STNTYPE", "EVAP_STN"]]
            .rename(columns={"Buis": "ID", "head": "head_measured(m)",
                             "time": "Date"})
            )
        stnType = (df_GW[["PREC_STNTYPE", "ID"]].drop_duplicates()
                   .set_index("ID"))
        #drop_duplicates(keep=True)
        df_GW = (df_GW.set_index("Date")
            .groupby("ID")
            .resample("d")
            .mean()
            .reset_index()
        )
        df_GW[["PREC_STN", "EVAP_STN"]] = df_GW[["PREC_STN",
                   "EVAP_STN"]].fillna(method="ffill")
        df_GW = pd.merge(df_GW, stnType, how="left", on=["ID"])

    # 7) merges precipitation, evapotranspiration and head data
    if knmi == "piezo":
        # merging precipitation,evapotranspiration and head data
        df_merge_P_E = pd.merge(df_P, df_E, how="left", on=["Date", "ID"])
        df_merge_P_E_GW = pd.merge(df_merge_P_E, df_GW, how="left",
                                   on=["Date", "ID"])
        print("merging precipitation,evapotranspiration and head data : done")
    else:
        df_merge_P_E_GW = df_GW

    return df_merge_P_E_GW, df_GW_info


def create_df_statistics_summary(df_merge_P_E_GW, df_GW_info, acorr):
    """
    This function creates an empty data frame for the summary of results
    per location including:

    1)	        rmse:root mean squared error of residual
    2)	        rmsn:root mean squared error of noise (innovations)
    3)	        1/(variance of residuals)
    4)	        std of residual
    5)	        number of observation points
    6)	        number of filtered observation points
    7)	        number of years in timeseries model
    8)	        Pearson R2
    9)	        explained variance score
    10)	        evapotranspiration factor
    11)	        stderr of evap factor
    12)	        M0
    13)	        stderr of M0
    14)         noise model decay alpha
    15)	        P50 of simulations
    16)	        P50 of observations
    17)	        height of residual slope
    18)         trend slope
    19)         step test
    20)         test drying up of screen
    21)	        autocorrelation test with lags of 1-365 days (if acorr=True)

    Parameters
    -----------
    df_merge_P_E_GW: pandas DataFrame
        merged precipitation, evapotranspiration and head data
    df_GW_info: pandas DataFrame
        GW info data frame (e.g. ID, X,Y,etc)
    acorr: boolean
        whether or not to include autocorrelation results


    Returns
    -------
    df_statistics_summary: pandas DataFrame
        An empty dataframe for the summary of results. Will be filled
        in the next step by "def time_series_analysis".

    """
    df_statistics_summary = pd.DataFrame(
        data = {
            "ID": df_merge_P_E_GW.ID.unique(),
            "root mean squared error of residual": np.nan,
            "root mean squared error of noise": np.nan,
            "1/Var(residual)": np.nan,
            "std(residual)": np.nan,
            "number of observation points": np.nan,
            "number of filtered observation points": np.nan,
            "number of years in timeseries model": np.nan,
            "Pearson R2": np.nan,
            "explained variance score": np.nan,
            "explained variance score detrended": np.nan,
            "evapotranspiration factor": np.nan,
            "evapotranspiration factor stderr": np.nan,
            "evapotranspiration factor stderr in decimal": np.nan,
            "M0": np.nan,
            "M0 stderr": np.nan,
            "M0 stderr in decimal": np.nan,
            "alpha": np.nan,
            "P50 of simulations": np.nan,
            "P50 of observations": np.nan,
            "height of residual slope": np.nan,
            "trend slope (m/d)": np.nan,
            "trend slope stderr": np.nan,
            "trend slope stderr in decimal": np.nan,
            "P-value t-test for step in mean": np.nan,
            "T-test statistic for step in mean": np.nan,
            "Relative step size": np.nan,
            "Date for step in mean": np.nan,
            "Drying up (1/0)": np.nan
        }
    )
    if acorr:
        df_statistics_summary[[
            "autocorrelation test for lag 1.0",
            "autocorrelation test for lag 2.0",
            "autocorrelation test for lag 3.0",
            "autocorrelation test for lag 4.0",
            "autocorrelation test for lag 5.0",
            "autocorrelation test for lag 6.0",
            "autocorrelation test for lag 7.0",
            "autocorrelation test for lag 8.0",
            "autocorrelation test for lag 9.0",
            "autocorrelation test for lag 10.0",
            "autocorrelation test for lag 12.0",
            "autocorrelation test for lag 13.0",
            "autocorrelation test for lag 14.0",
            "autocorrelation test for lag 30.0",
            "autocorrelation test for lag 61.0",
            "autocorrelation test for lag 90.0",
            "autocorrelation test for lag 120.0",
            "autocorrelation test for lag 150.0",
            "autocorrelation test for lag 180.0",
            "autocorrelation test for lag 210.0",
            "autocorrelation test for lag 240.0",
            "autocorrelation test for lag 270.0",
            "autocorrelation test for lag 300.0",
            "autocorrelation test for lag 330.0",
            "autocorrelation test for lag 365.0"
        ]] = np.nan

    df_statistics_summary = pd.merge(
        df_statistics_summary, df_GW_info, how="left", on="ID"
    )
    df_statistics_summary = df_statistics_summary.set_index("ID")
    return df_statistics_summary


def timeseries_analysis(df_statistics_summary, knmi, df_merge_P_E_GW, path,
                        rolling, acorr, subset):
    """
    This function performs time series analysis as follows:
        1)	reads data per location
        2)	tests for step and removes outliers from head data including
            heads <= -20 m+NAP
        3)	selecting period for training
            Note: The period for training in this script is defined based
            on the validation period of LHM 4.1. If user is using another
            dataset, period for training needs to be amended accordingly.
        4)	selecting period for simulation
            Note: The period for simulation in this script is defined based
            on the validation period of LHM 4.1. If user is using another
            dataset, period for simulation needs to be amended accordingly.
        5)	performs time series analysis using Pastas package
        6)	calculates defined statistics
        7)	addes statistics on df_statistics_summary
        8)	saves results as png
        9)	saves results as ipf
        10)	saves summary of results as xlsx

    Parameters
    ----------
    df_statistics_summary: pandas DataFrame
        Empty data frame for the summary of results. This will be
        created by function of "def create_df_statistics_summary ()",
        hence user does not need to make any action for this parameter.
    df_merge_P_E_GW: pandas DataFrame
        merged precipitation, evapotranspiration and head data
    path: Path object
        folder location for saving the results
    rolling: boolean
        If rolling is defined as “True”, a moving average of 2 and 3 will
        be calculated for precipitation and evapotranspiration time series.
        Time series analysis will then be performed on three sets of time
        series model including:
            i. Time series model with no moving average of precipitation
               and evapotranspiration
            ii. Time series model with moving average of “2” for
                precipitation and evapotranspiration
            iii. Time series model with moving average of “3”
                 for precipitation and evapotranspiration
        Among three time series models listed above the one which has
        the highest explained variance will be selected as the best
        time series model.
    acorr: boolean
        whether or not to include autocorrelation results

    Output
    ------
    png: simulations time series results, and head measurements
        will be saved as png
    ipf: simulations time series results, and head measurements
        will be saved as ipf
    df_statistics_summary: DataFrame containing the summary of results
        for all monitoring locations
    """
    control = 0
    print("Time series analysis : started")

    path = Path(path)
    # Creating folders for saving the results
    if rolling:
        outdir = path / "results" / "with rolling"
        outdir.mkdir(parents=True, exist_ok=True)
    else:
        outdir = path / "results" / "without rolling"
        outdir.mkdir(parents=True, exist_ok=True)

    if knmi == "station":
        # initialize dictionary to store downloaded knmi data
        knmiData = {"EV24": {},
                    "RH": {},
                    "RD": {}
        }
    print('Time series modelling')
    for location in tqdm(df_statistics_summary.index, total=len(df_statistics_summary.index)):
        control = control + 1

        #%% 1) selecting data per location
        df = df_merge_P_E_GW[df_merge_P_E_GW.ID == location]
        if knmi == "station":
            df_P, df_E, knmiData = get_knmi_data(knmiData, df, Path(path))
            # merge precipitation and evapotranspiration data with heads
            df_merge_P_E = pd.merge(df_P, df_E, how="left", on=["Date"])
            df = pd.merge(df.reset_index(), df_merge_P_E, how="right",
                          on=["Date"])
            df = df.drop(["PREC_STN", "EVAP_STN", "PREC_STNTYPE"],
                         axis="columns")

        df = df.sort_values(by=["Date"]).set_index("Date")

        # selecting the period of 2000-2019 (user needs to change
        # this if required)
        # this includes extra years at the start for initializing
        # time series model
        date_start = datetime.datetime(year=2000, month=1, day=1)
        date_end = datetime.datetime(year=2021, month=7, day=31)
        selected_period = pd.date_range(date_start, date_end, freq="d")
        df = pd.DataFrame(df, selected_period)

        #%% 2) testing for step and removing outliers from head data
        head_clean = df[["head_measured(m)"]].dropna()
        number_of_years_in_timeseries_model = (
            head_clean.index.max().year - head_clean.index.min().year + 1
        )
        # remove all values <= -20
        head_clean = head_clean[(head_clean["head_measured(m)"] > -20)]

        # test for drying up of screen
        drying_up = testDryingUp(head_clean["head_measured(m)"])

        # test for step in time series. If significant step, then only
        # use largest section of time series
        pval, tstat, relStep, stepDate = step(head_clean["head_measured(m)"])
        df_statistics_summary.loc[
            location, "P-value t-test for step in mean"
        ] = pval
        df_statistics_summary.loc[
            location, "T-test statistic for step in mean"
        ] = tstat
        df_statistics_summary.loc[
            location, "Relative step size"
        ] = relStep
        df_statistics_summary.loc[
            location, "Date for step in mean"
        ] = stepDate
        # if test shows significant step and relative step size > 3,
        # then only use largest part of series for analysis
        if pval < 0.01 and relStep > 3:
            delta1 = (stepDate - head_clean.index[0].date()).days
            delta2 = (head_clean.index[-1].date() - stepDate).days
            if delta1 >= delta2:
                head_clean = head_clean[:stepDate]
            else:
                head_clean = head_clean[stepDate:]

        # remove outliers using modified z-score test
        z_GW = outlier_test(head_clean["head_measured(m)"])
        threshold = 3
        if z_GW[pd.notnull(z_GW)].shape[0] > 0:
            head_clean = (head_clean)[(z_GW < threshold)]

        if (
            len(head_clean) < 3
            or len(df["Rainfall(mm)"].dropna()) < 3
            or len(df["evapotranspiration(mm)"].dropna()) < 3
        ):
            continue

        head_clean = head_clean.rename(
            columns={"head_measured(m)": "head_measured_clean(m)"}
        )
        df = pd.merge(df, head_clean, how="left", left_index=True,
                      right_index=True)

        #%% 3) selecting period for training
        training_start = df[["head_measured_clean(m)"]].dropna().index.min()
        if ("Provincie Noord-Brabant" in
            df_statistics_summary.loc[location, "EIGENAAR"]):
            training_end = min(
                datetime.datetime(year=2016, month=12, day=31),
                df[["head_measured_clean(m)"]].dropna().index.max(),
            )

        else:
            training_end = df[["head_measured_clean(m)"]].dropna().index.max()
        selected_period_for_training = pd.date_range(
            training_start, training_end, freq="d"
        )
        df_training = pd.DataFrame(df, selected_period_for_training)

        #%% 4) selected period for simulation
        simulation_start = max(
            df["Rainfall(mm)"].dropna().index.min(),
            df["evapotranspiration(mm)"].dropna().index.min(),
        )
        simulation_end = min(
            df["Rainfall(mm)"].dropna().index.max(),
            df["evapotranspiration(mm)"].dropna().index.max(),
        )

        #%% 5) performs time series analysis using Pastas package
        if rolling:
            ml1 = []
            ml2 = []
            ml3 = []
            a = []
            b = []
            c = []
            ml = []
            ml1 = ps.Model(df_training["head_measured_clean(m)"])
            ml2 = ps.Model(df_training["head_measured_clean(m)"])
            ml3 = ps.Model(df_training["head_measured_clean(m)"])
            # adding moving average columns for rainfall and evapotranspiration
            df.loc[:, "Rainfall_rolling_2(mm)"] = (
                df["Rainfall(mm)"].rolling(window=2).mean().round(3)
            )
            df.loc[:, "Rainfall_rolling_3(mm)"] = (
                df["Rainfall(mm)"].rolling(window=3).mean().round(3)
            )
            df.loc[:, "evapotranspiration_rolling_2(mm)"] = (
                df["evapotranspiration(mm)"].rolling(window=2).mean().round(3)
            )
            df.loc[:, "evapotranspiration_rolling_3(mm)"] = (
                df["evapotranspiration(mm)"].rolling(window=3).mean().round(3)
            )

            ts1 = ps.StressModel2(
                [df["Rainfall(mm)"], df["evapotranspiration(mm)"]],
                ps.Gamma,
                name="rainevap",
                settings=("prec", "evap"),
            )
            ts2 = ps.StressModel2(
                [df["Rainfall_rolling_2(mm)"], df["evapotranspiration_rolling_2(mm)"]],
                ps.Gamma,
                name="rainevap",
                settings=("prec", "evap"),
            )
            ts3 = ps.StressModel2(
                [df["Rainfall_rolling_3(mm)"], df["evapotranspiration_rolling_3(mm)"]],
                ps.Gamma,
                name="rainevap",
                settings=("prec", "evap"),
            )

            ml1.add_stressmodel(ts1)
            ml2.add_stressmodel(ts2)
            ml3.add_stressmodel(ts3)

            #  Solve the model
            # set tmin to 1-1-2011 as start for calibration
            ml1.solve(tmin="2000", warmup=1826)
            ml2.solve(tmin="2000", warmup=1826)
            ml3.solve(tmin="2000", warmup=1826)
            a = ml1.stats.evp()
            b = ml2.stats.evp()
            c = ml3.stats.evp()
            ml = []
            # select the best model among rolling 1,2,3 based on highest explained variance score
            if a > b and a > c:
                ml = ml1

            if b > a and b > c:
                ml = ml2

            if c > a and c > b:
                ml = ml3
            if a == b == c:
                ml = ml1
        else:
            ml = []
            ml = ps.Model(df_training["head_measured_clean(m)"])
            ts = ps.StressModel2(
                [df["Rainfall(mm)"], df["evapotranspiration(mm)"]],
                ps.Gamma,
                name="rainevap",
                settings=("prec", "evap"),
            )
            ml.add_stressmodel(ts)
            # set tmin to 1-1-2011 as start for calibration
            ml.solve(tmin="2000", warmup=1826)

        # test for trend in residuals
        slope, trend = trendTest(ml.residuals())
        if trend:
            # add trend model and solve model again
            tm = ps.LinearTrend(
                start=training_start.strftime("%Y-%m-%d"),
                end=training_end.strftime("%Y-%m-%d"),
                name="linear_trend"
                )
            ml.add_stressmodel(tm)
            ml.set_parameter("linear_trend_tstart", vary=False)
            ml.set_parameter("linear_trend_tend", vary=False)
            ml.solve(tmin="2000", warmup=1826)

            detrended = (df_training["head_measured_clean(m)"]
                         - ml.get_contribution("linear_trend")
                         - ml.get_parameters("constant")[0])
            sim = ml.get_contribution("rainevap")
            # calculate evp for detrended series
            evp_detrended = ps.stats.metrics.evp(sim=sim, obs=detrended)
        else:
            evp_detrended = np.NaN

        simulation = ml.simulate(tmin=simulation_start, tmax=simulation_end)
        simulation = simulation.loc[training_start:training_end]

        #%% 6) calculates required statistics

        # autocorrelation_test on noise
        if acorr:
            autocorrelation_test = ps.stats.tests.ljung_box(ml.noise())[1][
                "Accept Ha (alpha=0.05)"
            ].reset_index()
            autocorrelation_test.loc[:, "autocorrelation test for lag"] = [
                "autocorrelation test for lag "
                + str(autocorrelation_test.loc[index, "Lags (Days)"])
                for index in autocorrelation_test.index
            ]
            autocorrelation_test = autocorrelation_test.set_index(
                "autocorrelation test for lag"
            )

        # getting model parameters (M0,evap_factor,trend_slope) and  stderr
        parameters = ml.parameters.loc[:, ["optimal", "stderr", "initial", "vary"]]
        parameters.loc[:, "stderr"] = (
            (parameters.loc[:, "stderr"] / parameters.loc[:, "optimal"])
            .abs()
            .apply("\u00B1{:.2%}".format)
        )

        evap_factor = parameters.loc["rainevap_f", "optimal"]
        evap_factor_stderr = parameters.loc["rainevap_f", "stderr"]
        evap_factor_stderr_decimal = round(
            float(evap_factor_stderr[1:].replace("%", "")) / 100, 3
        )

        M0 = parameters.loc["rainevap_A", "optimal"]
        M0_stderr = parameters.loc["rainevap_A", "stderr"]
        M0_stderr_decimal = round(float(M0_stderr[1:].replace("%", "")) / 100, 3)

        if trend:
            trend_slope = parameters.loc["linear_trend_a", "optimal"]
            trend_slope_stderr = parameters.loc["linear_trend_a", "stderr"]
            trend_slope_stderr_decimal = round(
                float(trend_slope_stderr[1:].replace("%", "")) / 100, 3
            )
        else:
            trend_slope = np.NaN
            trend_slope_stderr = np.NaN
            trend_slope_stderr_decimal = np.NaN


        # getting noise model parameter alpha
        alpha = parameters.loc["noise_alpha", "optimal"]

        # estimating P50_simulated and measured data

        P50_simulated = np.percentile(simulation.values, 50)
        P50_measured = np.percentile(ml.observations().values, 50)

        # check trend in residuals()
        x = np.array([ml.residuals().reset_index().index]).reshape((-1, 1))
        y = np.array(ml.residuals().values)
        reg = LinearRegression().fit(x, y)
        y_min = reg.predict(np.array([x.min()]).reshape((-1, 1)))
        y_max = reg.predict(np.array([x.max()]).reshape((-1, 1)))
        residual_slope_height = float(y_max - y_min)

        # adding the simulation results to df
        df.loc[:, "head_simulated(m)"] = simulation

        df = df.loc["2011-01-01":, ["head_measured(m)",
                                    "head_measured_clean(m)",
                                    "head_simulated(m)",]].round(3)
        df = df.dropna(subset=["head_measured(m)"])

        #%% 7) saves statistics on df_statistics_summary

        df_statistics_summary.loc[
            location, "root mean squared error of residual"
        ] = ml.stats.rmse()
        df_statistics_summary.loc[
            location, "root mean squared error of noise"
        ] = ml.stats.rmsn()
        df_statistics_summary.loc[location, "1/Var(residual)"] = (
            1 / ml.residuals().var()
        )
        df_statistics_summary.loc[
            location, "std(residual)"] = ml.residuals().std()
        df_statistics_summary.loc[
            location, "number of observation points"
            ] = len(df[["head_measured(m)"]].dropna())
        df_statistics_summary.loc[
            location, "number of filtered observation points"
            ] = len(df[["head_measured_clean(m)"]].dropna())
        df_statistics_summary.loc[
            location, "number of years in timeseries model"
            ] = number_of_years_in_timeseries_model
        df_statistics_summary.loc[
            location, "explained variance score"] = ml.stats.evp()
        df_statistics_summary.loc[
            location, "explained variance score detrended"] = evp_detrended
        df_statistics_summary.loc[
            location, "Pearson R2"] = ml.stats.rsq()
        df_statistics_summary.loc[
            location, "evapotranspiration factor"] = evap_factor
        df_statistics_summary.loc[
            location, "evapotranspiration factor stderr"] = evap_factor_stderr
        df_statistics_summary.loc[
            location, "evapotranspiration factor stderr in decimal"
            ] = evap_factor_stderr_decimal
        df_statistics_summary.loc[
            location, "M0"] = M0
        df_statistics_summary.loc[
            location, "M0 stderr"] = M0_stderr
        df_statistics_summary.loc[
            location, "M0 stderr in decimal"] = M0_stderr_decimal
        df_statistics_summary.loc[
            location, "alpha"] = alpha
        df_statistics_summary.loc[
            location, "P50 of simulations"] = P50_simulated
        df_statistics_summary.loc[
            location, "P50 of observations"] = P50_measured
        df_statistics_summary.loc[
            location, "trend slope (m/d)"] = trend_slope
        df_statistics_summary.loc[
            location, "trend slope stderr"] = trend_slope_stderr
        df_statistics_summary.loc[
            location, "trend slope stderr in decimal"
            ] = trend_slope_stderr_decimal
        df_statistics_summary.loc[
            location, "height of residual slope"] = residual_slope_height

        df_statistics_summary.loc[
            location, "drying up of screen (1/0)"] = drying_up

        if acorr:
            for index in autocorrelation_test.index:
                df_statistics_summary.loc[
                    location, index] = autocorrelation_test.loc[
                    index, "Accept Ha (alpha=0.05)"]

        # save results
        #%% 8) saves results as png

        fig = plt.figure()
        plt.style.use("ggplot")
        ax1 = fig.add_subplot()
        ax1.plot(
            df["head_measured(m)"],
            marker="o",
            color="green",
            linestyle="None",
            markersize=1,
            label="measured",
            zorder=0
        )
        ax1.plot(
            df["head_measured_clean(m)"],
            marker="o",
            color="steelblue",
            linestyle="None",
            markersize=1.5,
            label="measured_clean",
            zorder=1
        )
        ax1.plot(df["head_simulated(m)"], color="red", label="simulated")
        ax1.yaxis.set_tick_params(labelsize=7)
        plt.gca().xaxis.set_major_locator(mdates.YearLocator())
        plt.ylabel("head (m)", fontsize=7, fontweight="bold")

        ymin = (
            min(
                df["head_measured_clean(m)"].dropna().values.min(),
                df["head_simulated(m)"].dropna().values.min(),
            )
        )
        ymax = (
            max(
                df["head_measured_clean(m)"].dropna().values.max(),
                df["head_simulated(m)"].dropna().values.max(),
            )
        )
        yrange = ymax - ymin
        ymax = ymax + 0.1*yrange
        ymin = ymin - 0.1*yrange
        plt.ylim(ymin, ymax)
        legend = plt.legend(
            bbox_to_anchor=(1.08, 1), loc="upper left",
            prop={"size": 6}, shadow=True)
        legend.get_frame().set_facecolor("white")
        ax1.set_axisbelow(True)
        ax1.yaxis.grid(color="white", linestyle="dashed")
        ax1.xaxis.grid(color="white", linestyle="dashed")
        ax1.xaxis.set_tick_params(labelsize=7)
        ax1.yaxis.grid(color="white", linestyle="dashed")
        ax1.xaxis.grid(color="white", linestyle="dashed")
        fig.autofmt_xdate()
        plt.suptitle(location, x=0.45, y=1.02, fontsize=8, fontweight="bold")
        plt.tight_layout()
        plt.subplots_adjust(hspace=0.05)
        outdir_png = outdir / "png"
        outdir_png.mkdir(parents=True, exist_ok=True)
        fig.savefig(outdir_png / (str(location) + ".png"),
                    bbox_inches="tight", dpi=150)
        plt.close()
        plt.cla()
        plt.clf()

        #%% 9) saves results as ipf
        df.loc[:, "time"] = df.index
        outdir_ipf = outdir / "ipf"
        outdir_ipf.mkdir(parents=True, exist_ok=True)
        imod.ipf.write_assoc(
            outdir_ipf / (str(location) + ".txt"),
            df,
            itype="timeseries",
            nodata=-9999.0,
        )

        print(str(control) + " out of "
              + str(len(df_merge_P_E_GW.ID.unique())))
        # store statistics in between
        if int(control) == 100*int(control/100):
            print("Storing results in between")
            outdir_statistics = outdir / "statistics"
            outdir_statistics.mkdir(parents=True, exist_ok=True)
            df_statistics_summary.to_excel(outdir_statistics /
                                           "statistics_summary_temp.xlsx")

    #%% 10) save df_statistics_summary
    outdir_statistics = outdir / "statistics"
    outdir_statistics.mkdir(parents=True, exist_ok=True)
    df_statistics_summary.to_excel(outdir_statistics
                                   / "statistics_summary.xlsx")
    print("Time series analysis : Done")

