"""
Script to process data of head observations from regional water managers.
Includes preselection: time series with
    1) less than 3 hydrological years (1 April - 31 March); or
    2) less than 24 measurements in period of 1-4-2018 - 1-10-2018
    are skipped.
The first criterium is used to select only series for which the GXG can
be calculated over at least 3 hydrological years. The second criterium is
used to allow for analysis of recession during drought of 2018.

Note:
This script is developed based on the data format of input files
used and files generated during validation of LHM 4.1.
Hence, if the user uses different data with different formatting,
all necessary amendments first need to be implemented by the user.
"""

import imod
import numpy as np
import os
import pandas as pd
from datetime import datetime
from scipy import stats
from hdsrhipy.groundwater.ModelLayer import ModelLayer
from hdsrhipy.groundwater.stats import *
import warnings
from pathlib import Path

def read_csv(path=None, catalogue_file=None, timeseries_file=None):
    
    dataset = 'regio'
    eigenaar = 'HDSR'
    outRoot = path

    startyear = 2011
    endyear = pd.to_datetime('now').year
    crit = [3, 20] # timeseries is selected if GXG can be calculated
               # based on mimumum of 3 years with 20 GXG observations each

    sel = pd.DataFrame(columns=['beheerder', 'selected'])
    reces = pd.DataFrame(columns=['beheerder', 'selected'])
    desel = pd.DataFrame(columns=['beheerder', 'deselected'])


    outfolder = os.path.join(path, 'data')
    if not(os.path.exists(outfolder)):    
        os.mkdir(outfolder)
    outfile_catalogus = os.path.join(outfolder, 'catalogus_selected.ipf')
    outfile_catalogus_csv = os.path.join(outfolder, 'catalogus_selected.csv')
    path_catalogus = os.path.join(path, catalogue_file)
                                  #'GrondwaterMeetpunten_aanwezig.csv')
    path_timeseries = os.path.join(path, timeseries_file )        
                                  #'GrondwaterMeetpunten_reeksen.csv')

    # read and rework catalogus for further processing
    df_catalogus = pd.read_csv(path_catalogus, index_col=0, header=0)
    df_catalogus.rename(columns={
        'PutFilter': 'Buis',
        'x': 'X_RD_CRD',
        'y': 'Y_RD_CRD',
        'BKf (m NAP)':'FILTER_TOP_NAP',
        'OKf (m NAP)': 'FILTER_ONDER_NAP'
    }, inplace=True)
    # make Buis unique by keeping last row
    df_catalogus = df_catalogus.dropna(subset=['X_RD_CRD', 'Y_RD_CRD'],
                                        how='any')
    df_catalogus.set_index('Buis', inplace=True)
    
    df_catalogus['EIGENAAR'] = eigenaar
    if 'WVP' not in df_catalogus.columns:
        df_catalogus['WVP'] = np.NaN
    df_catalogus['GHG'] = np.NaN
    df_catalogus['GLG'] = np.NaN
    df_catalogus['NGXG'] = np.NaN
    df_catalogus['GVG'] = np.NaN
    df_catalogus['DYNAMIEK'] = np.NaN
    
    
    selected = []
    recession = []
    deselected = []
    
    # read timeseries file filter by filter and write associated .txt file
    # when selection criteria are met.
    # remove filter from catalogus if selection criteria are not met.
    df_timeseries = pd.read_csv(path_timeseries, index_col=0, header=0)
    df_timeseries.index = pd.to_datetime(df_timeseries.index)
    # make index timezone naive
    df_timeseries.index = pd.DatetimeIndex([i.replace(tzinfo=None)
                                            for i in df_timeseries.index])
    # select years in period between startyear and endyear
    df_timeseries = df_timeseries.loc[(
        df_timeseries.index.year >= startyear)
        & (df_timeseries.index.year <= endyear), :]
    
    for buisname in df_catalogus.index:
        if buisname in df_timeseries.columns:
            df = df_timeseries.loc[:, buisname].to_frame()
            # rework df
            df.rename(columns={buisname: 'head'}, inplace=True)
            df['head'] = df['head'].round(3)
            # remove nodata
            df.dropna(subset=['head'], inplace=True)
            # apply simple outlier test before calculating GXG
            with warnings.catch_warnings():
                warnings.simplefilter(action="ignore", category=RuntimeWarning)
                z_GW = np.abs(stats.zscore(df['head']))
                threshold = 3
                df = (df)[(z_GW < threshold)]
            # calculate GXG and preselect series if GXG can be calculated
            gxg = getGXG(df['head'], minyear=crit[0], minval=crit[1],
                          nearest=True, nearlim=4)
            if np.isfinite(gxg['glg']):
                df.index.rename('time', inplace=True)
                df.reset_index(inplace=True)
                selected.append([eigenaar, buisname])
                df_catalogus.loc[buisname, 'GHG'] = gxg['ghg']
                df_catalogus.loc[buisname, 'GLG'] = gxg['glg']
                df_catalogus.loc[buisname, 'NGXG'] = gxg['ngxg']
                df_catalogus.loc[buisname, 'GVG'] = gxg['gvg']
                df_catalogus.loc[buisname, 'DYNAMIEK'] = (gxg['ghg']
                                                          - gxg['glg']
                                                          ).round(2)
                # write to associated file
                path2 = os.path.join(outfolder, buisname + '.txt')
                imod.ipf.write_assoc(path2, df, itype=1, nodata=-9999.)
            elif df[((df.index >= datetime.strptime('1-4-2018',
                                                    '%d-%m-%Y')) &
                        (df.index < datetime.strptime('1-10-2018',
                                                      '%d-%m-%Y'))
                        )
                    ].shape[0] > 24: # preselect series for recession analysis
                df.index.rename('time', inplace=True)
                df.reset_index(inplace=True)
                recession.append([eigenaar, buisname])
                path2 = os.path.join(outfolder, buisname + '.txt')
                imod.ipf.write_assoc(path2, df, itype=1, nodata=-9999.)
            else:
                deselected.append([eigenaar, buisname])
    
    print("# selected = " + str(len(selected)))
    print("# recession = " + str(len(recession)))
    print("# deselected = " + str(len(deselected)))
    
    lst = selected
    lst.extend(recession)
    lst = [x[1] for x in lst]
    
    df_catalogus = df_catalogus[['X_RD_CRD',
                                  'Y_RD_CRD',
                                  'FILTER_TOP_NAP',
                                  'FILTER_ONDER_NAP',
                                  'WVP',
                                  'GHG',
                                  'GLG',
                                  'GVG',
                                  'DYNAMIEK',
                                  'NGXG',
                                  'EIGENAAR']]
    
    
    df_selected = df_catalogus.loc[df_catalogus.index.isin(lst)].reset_index()
    print ('   Total added:', df_selected.shape[0])
    
    sel = sel.append(pd.DataFrame(selected,
                                  columns=['beheerder', 'selected'])
                      )
    reces = reces.append(pd.DataFrame(recession,
                                        columns=['beheerder', 'selected'])
                            )
    desel = desel.append(pd.DataFrame(deselected,
                                      columns=['beheerder', 'deselected'])
    
                          )
    
    sel.to_csv(os.path.join(outfolder, 'selected.csv'), index=False)
    reces.to_csv(os.path.join(outfolder, 'recession.csv'), index=False)
    desel.to_csv(os.path.join(outfolder, 'deselected.csv'), index=False)
    
    df_selected = df_selected[['X_RD_CRD',
                                'Y_RD_CRD',
                                'Buis',
                                'FILTER_TOP_NAP',
                                'FILTER_ONDER_NAP',
                                'WVP',
                                'GHG',
                                'GLG',
                                'GVG',
                                'DYNAMIEK',
                                'NGXG',
                                'EIGENAAR']]
    
    df_selected.to_csv(outfile_catalogus_csv)
    
    # write catalogus as ipf
    imod.ipf.write(outfile_catalogus, df_selected, indexcolumn=3,
                    assoc_ext='txt', nodata=1e+20)
    
    
    baseName = 'catalogus_selected'
    fName = os.path.join(outRoot, 'data', baseName + '.csv')
    data = pd.read_csv(fName, index_col=0, header=0)
    
    if 'WVP' in data.columns:
        data.rename(columns={'WVP': 'LHM LAYER'}, inplace=True)
    else:
        data['LHM LAYER'] = np.NaN
    
    ml = ModelLayer(lhmfolder=os.path.join(path,'lhm'))
    
    layer = (data[['FILTER_TOP_NAP',
                  'FILTER_ONDER_NAP',
                  'X_RD_CRD',
                  'Y_RD_CRD']].astype(float)
                              .apply(lambda x: ml.getLHMLayer(*x), axis=1)
                              .round(0)
                )
    
    data['LHM LAYER'].fillna(layer, inplace=True)
    
    outName = os.path.join(outRoot,'data', baseName + '_layer.csv')
    
    data.to_csv(outName)