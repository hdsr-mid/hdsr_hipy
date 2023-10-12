# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 09:31:17 2018

@author: Artesia
"""
import pandas as pd
from hdsrhipy.model import util
from collections import OrderedDict
import numpy as np
import csv
import os

def write_para_sim(para_sim,fname):
    """
    Write para_sim.inp

    Parameters
    ----------
    para_sim: Pandas Series
        Series containing the para_sim input
    fname: str
        Filename (full path) of the file to write

    Returns
    -------

    """
    with open(fname, "w") as f:
        f.write("*\n")
        f.write("* Simulation period\n")
        f.write("*\n")
        f.write('      IDBG              =       {}\n'.format(para_sim['IDBG']))
        f.write('      IYBG              =       {}\n'.format(para_sim['IYBG']))
        f.write('      IDBGSM            =       {}\n'.format(para_sim['IDBGSM']))
        f.write('      IDEDSM            =       {}\n'.format(para_sim['IDEDSM']))
        f.write("*\n")
        f.write("* Parameters for numerical approximation\n")
        f.write("*\n")
        f.write('      DTGW              =       {}\n'.format(para_sim['DTGW']))
        f.write('      DTSW              =       {}\n'.format(para_sim['DTSW']))
        f.write('      DHMXSW            =       {}\n'.format(para_sim['DHMXSW']))
        if 'unsa_svat_path' in para_sim:
            f.write('      unsa_svat_path    = {}\n'.format(para_sim['unsa_svat_path']))
        f.write("*\n")
        f.write("* Parameters for output\n")
        f.write("*\n")
        f.write('      msg               =       {}\n'.format(para_sim['msg']))
        f.write('      svat_gt           =       {}\n'.format(para_sim['svat_gt']))
        f.write('      svat_per          =       {}\n'.format(para_sim['svat_per']))
        f.write('      svat_per_csv      =       {}\n'.format(para_sim['svat_per_csv']))
        f.write('      svat_dtgw         =       {}\n'.format(para_sim['svat_dtgw']))
        f.write('      svat_dtgw_csv     =       {}\n'.format(para_sim['svat_dtgw_csv']))
        f.write('      svat2gwdtgw       =       {}\n'.format(para_sim['svat2gwdtgw']))
        f.write('      drlink_per        =       {}\n'.format(para_sim['drlink_per']))
        f.write('      sw_per            =       {}\n'.format(para_sim['sw_per']))
        f.write('      sw_per_csv        =       {}\n'.format(para_sim['sw_per_csv']))
        f.write('      sw_dtgw           =       {}\n'.format(para_sim['sw_dtgw']))
        f.write('      sw_dtgw_csv       =       {}\n'.format(para_sim['sw_dtgw_csv']))
        f.write('      sw_hq_dtgw        =       {}\n'.format(para_sim['sw_hq_dtgw']))
        f.write('      sw_dtsw           =       {}\n'.format(para_sim['sw_dtsw']))
        f.write('      sw_hq_dtsw        =       {}\n'.format(para_sim['sw_hq_dtsw']))
        f.write('      ql_svat_per       =       {}\n'.format(para_sim['ql_svat_per']))
        f.write('      ql_sw_per         =       {}\n'.format(para_sim['ql_sw_per']))
        f.write("*\n")
        f.write("* <<end-of-data>>\n")
        
def read_para_sim(fname):
    """
    Read para_sim.inp

    Parameters
    ----------
    fname: str
        Filename (full path) of the file to read

    Returns
    -------
    fname: dict

    """
    with open(fname, "r") as f:
         lines = f.readlines()
    para_sim={}
    for line in lines:
        if not line.startswith('*'):
            parts=line.strip().split('=')
            key = parts[0].strip()
            value = parts[1].split('!')[0]
            if util.is_number(parts[1]):
                para_sim[key] = pd.to_numeric(value)
            else:
                para_sim[key] = value.strip()
    para_sim = pd.Series(para_sim)
    return para_sim

def variable_format(name):
    """
    Determine the format of an inp-file

    Parameters
    ----------
    name: str
        The filename of the file to determine the format for

    Returns
    -------
    fmt: OrderedDict
        Containg the column names as keys and the data-format of the column as values

    """
    # copied formats from SIMGRO V7.1.4 Input and output reference manual
    name = name.lower()
    if name == 'mete_stat.inp':
        fmt = OrderedDict(nmme='I10',lat='F10',alt='F10',zmeasw='F10')
    elif name == 'mete_svat.inp':
        fmt = OrderedDict(td='F15',iy='I5',prec='F10',etref='F10',nmme='I10',
                          tempmn='F10',tempmx='F10',temp='F10',Nrel='F10',
                          rad='F10',hum='F10',wind='F10')
    elif name == 'luse_svat.inp':
        fmt = OrderedDict(lu='I6',luna='A20',vglu='I6',alfafunclu='I6',
                          p1fd='F8',p2fd='F8',p3hfd='F8',p3lfd='F8',p4fd='F8',
                          t3hfd='F8',t3lfd='F8',pbgsplu='F8',frevsplu='F8',
                          gisplu='F8',tigisplu='F8',rpsplu='F8',tdbgsplu='F8',
                          tdedsplu='F8',fecmnlu='F8',albedolu='F8',
                          rscdrylu='F8',rscwetlu='F8',ECmaxlu='F8',
                          ECsloplu='F8')
    elif name == 'fact_svat.inp':
        fmt = OrderedDict(vg='I6',dy='I6',csvg='F8',laivg='F8',vxicvg='F8',
                          faevvg='F8',faeivg='F8',faebsvg='F8',faepdvg='F8',
                          chvg='F8',drpzvg='F8')
    elif name == 'tiop_sim.inp':
        fmt = OrderedDict(td='F15',iy='I6',io='I6',ip='I6')
    elif name == 'area_svat.inp':
        fmt = OrderedDict(svat='I10',ark='F10',glk='F8',tempCbotk='F8',
                          slk='I6',dum='A16',luk='I6',dprzk='F8',nm='I10',
                          cfPm='F8',cfETref='F8')
    elif name == 'infi_svat.inp':
        fmt = OrderedDict(svat='I10',qinfbasic='F8',ctop_down='F8',
                          ctop_up='F8',sc2='F8')
    elif name == 'mete_grid.inp':
        fmt = OrderedDict(td='F',iy='I',precgrid='A',etrefgrid='A',
                          tempmngrid='A',tempmxgrid='A',tempgrid='A',
                          Nrelgrid='A',radgrid='A',humgrid='A',windgrid='A')
    elif name == 'mod2svat.inp':
        fmt = OrderedDict(modfcell='I10',dum='A2',svat='I10',ly='I5')
    elif name == 'svat2swnr_roff.inp':
        fmt = OrderedDict(svat='I10',swnr='I10',vxmu='F8',crunoff='F8',
                          crunon='F8')
    elif name == 'svat2swnr_drng.inp':
        fmt = OrderedDict(svat='I10',sy='I6',dpsw='F8',wisw='F8',
                          adsw='F8',ddsw='F8',lesw='F8',redr='F8',reen='F8',
                          rein='F8',reex='F8',swnr='I10',lvsw='F8')
    elif name == 'mana_res.inp':
        fmt = OrderedDict(swnr='I10',dum='A6',ioma='I6',lvtasm='F8',
                          lvtawt='F8',dum2='A8',dptasu='F8',fxsuswsb='F8',
                          ndta='F10',iotasmnd='I6',iotawtnd='I6',swnrta='I10',
                          iotasmsb='I6',iotawtsb='I6')
    elif name == 'goto_res.inp':
        fmt = OrderedDict(swnr='I10',swnrgo='I10',lvwrsm='F8',lvwrwt='F8',
                          lvwrlw='F8',ndwr='I10',iowrsmnd='I6',iowrwtnd='I6',
                          swnrwr='I10',iowrsmsb='I6', iowrwtsb='I6',
                          iofwbk='I6',glnr='F8',lvcv='F8',alfa='F8',beta='F8')
    elif name == 'dish_res.inp':
        fmt = OrderedDict(swnr='I10',dhwr='F8',fmwr='F8',fswr='F8',
                          swnrgo='I10')
    elif name == 'swnr_sim.inp':
        fmt = OrderedDict(swnr='I10')
    elif name == 'svat2swnr_roff.inp':
        fmt = OrderedDict(svat='I10',swnr='I10',vxmu='F8',crunoff='F8',
                          crunon='F8')
    elif name == 'svat2swnr_drng.inp':
        fmt = OrderedDict(svat='I10',sy='I6',dpsw='F8',wisw='F8',adsw='F8',
                          ddsw='F8',lesw='F8',redr='F8',reen='F8',rein='F8',
                          reex='F8',lvsw='I10',swnr='F8')
    elif name == 'sel_key_svat_dtgw.inp' or name == 'sel_key_svat_per.inp':
        fmt = OrderedDict(key='A10',ioptkey='I10')
    elif name == 'sel_svat_bda.inp' or name == 'sel_svat_dtgw_bda.inp' or name == 'sel_svat_per_bda.inp':
        fmt = OrderedDict(svat='I10')
    else:
        raise(ValueError('Unknown inp-file: {}'.format(name)))
    
    return fmt
        
def read(fname,path=''):
    """
    Read an inp-file.

    Parameters
    ----------
    fname : str
        Name of the file to read (can include the folder)
    path : str
        Name of the path (absolute or relative to SIPS) of fname

    Returns
    -------
    Pandas Series when name is para_sim.inp
    for other inp-files:
    Pandas DataFrame

    """
    
    name = fname.split('\\')[-1]
    fname = fname.replace('\\','/')
    path_fname = os.path.join(path,fname)
    if name=='para_sim.inp':        
        df = read_para_sim(path_fname)
    else:
        fmt_dict = variable_format(name)
        if len(next(iter(fmt_dict.values())))>1:
            # read fixed format
            widths = [int(x[1:]) for x in fmt_dict.values()]
            df = pd.read_fwf(path_fname,widths=widths,names=list(fmt_dict))
        else:
            # read free format
            df = pd.read_csv(path_fname,names=list(fmt_dict))
    df.index.name=fname
    return df

def write(data,fname):
    """
    Write a inp-file

    Parameters
    ----------
    data : Pandas Series or Pandas DataFrame
        Pandas Series when name is para_sim.inp, Pandas DataFrame otherwise
    fname: str
        Full filename of output-file

    Returns
    -------

    """
    # make sure the pathname exists
    fname = fname.replace('/','\\')
    pathname = os.path.dirname(fname)
    if len(pathname)>0 and not os.path.isdir(pathname):
        os.makedirs(pathname)
    name = fname.split('\\')[-1]
    if name=='para_sim.inp':
        write_para_sim(data,fname)
    elif name=='tiop_sim.inp':
        with open(fname,'w') as f:
            for idx,d in data.iterrows():
                f.write(f"{d.td:15.2f}{d.iy:6g}{d.io:6g}\n")
    else:
        fmt_dict = variable_format(name)
        for key in fmt_dict.keys():
            if key not in data.columns:
                print('Adding column {} to {}'.format(key, fname))
                data[key] = np.nan
        # set the columns in the right order
        data = data[list(fmt_dict.keys())]
        
        mask = data.isna()
        if np.any(mask):
            # remove the last columns that only contain NaN's
            last = np.where(~np.all(mask,axis=0))[0][-1]+1
            data = data.iloc[:,:last]
            keys = list(fmt_dict.keys())
            for key in keys[last:]:
                fmt_dict.pop(key)
        
        # check if the format-string contains fieldlength (string is longer than 1 character) 
        if len(next(iter(fmt_dict.values())))>1:
            # when the entire column is NaN, change the format from int to float
            # otherwise writing will result in an error
            for key in fmt_dict.keys():
                if data[key].isna().all():
                    fmt_dict[key] = 'F{}'.format(fmt_dict[key][1:])
            
            # write using fixed format
            if False:
                # let savetxt do the conversion
                # make a string format
                fmt=''
                for val in fmt_dict.values():
                    # transform Simgro-format to numpy string-format
                    fmt = '{}%{}{}'.format(fmt,val[1:],val[0].lower())
                np.savetxt(fname, data.values, fmt=fmt)
            else:
                # convert to string by hand, then write to file
                data_str = data.copy()
                for key in fmt_dict.keys():
                    if fmt_dict[key][0]=='F':
                        n = int(fmt_dict[key][1:])-2
                        data_str[key] = data_str[key].apply(lambda t: '' if pd.isna(t) else str(np.round(t,n)))
                    else:
                        data_str[key] = data_str[key].astype(str)
                fmt=''
                for val in fmt_dict.values():
                    if val[0]=='A':
                        # left align strings
                        fmt = '{}%-{n}.{n}s'.format(fmt,n=val[1:])
                    else:
                        # otherwise, just use right align
                        fmt = '{}%{n}.{n}s'.format(fmt,n=val[1:])
                np.savetxt(fname, data_str.values, fmt=fmt)
        else:
            # write using free format
            data.to_csv(fname,header=False,index=False,quoting=csv.QUOTE_NONNUMERIC)