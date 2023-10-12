# -*- coding: utf-8 -*-
"""
Created on Fri Aug 31 16:54:44 2018

@author: ruben
"""
import pandas as pd
import numpy as np
import os
from hdsrhipy.model import util
from hdsrhipy.model import inp
from hdsrhipy.model import rasters
import shutil
import subprocess as sp
from osgeo import gdal
import imod
from datetime import datetime
from datetime import timedelta

class Runfile():
    """The Runfile class, to read and generate an iMod-runfile."""
    
    def __init__(self, fname, data_path='', version=2.21, precip_path=None, evap_path=None):
        d={}
        with open(fname,'r') as f:
            self.version = version
            
            # read the runfile-data from an existing file
            self.data_path = data_path
            
            if precip_path is None:                
                self.precip_path = os.path.join(data_path)
            else:
                self.precip_path = precip_path
            
            if evap_path is None:                
                self.evap_path = os.path.join(data_path)
            else:
                self.evap_path = evap_path
            
            
            # 1
            # use strip() to remove whitespace characters like `\n` at the end of each line
            d['OUTPUTDIRECTORY'] = self._read_path(f.readline().strip())
            
            # 2
            r = f.readline().replace(',',' ').split()
            d['NLAY'    ] = int(r[0])
            d['MXNLAY'  ] = int(r[1])
            d['NPER'    ] = int(r[2])
            if version<3:
                d['ISS'     ] = int(r[3])
                d['ISCL'    ] = int(r[4])
            else:
                d['SDATE'     ] = int(r[3])
                d['NSCL'    ] = int(r[4])
            d['IFTEST'  ] = int(r[5])
            d['ICONCHK' ] = int(r[6])
            d['IIPF'    ] = int(r[7])
            if version>=3:
                if len(r)>8:
                    d['IUNCONF'] = int(r[8])
                if len(r)>9:
                    d['IFVDL'] = int(r[9])
                if len(r)>10:
                    d['IARMWP'] = int(r[10])
                if len(r)>11:
                    d['IBNDCHK'] = int(r[11])
            for i in range(abs(d['IIPF'])):
                #f.readline()
                raise(NotImplementedError('specify artificial recharge to be read from an IPF file not supported yet'))
            
            # 3
            r = f.readline().replace(',',' ').split()
            d['NMULT'   ] = int(r[0])
            d['IDEBUG'  ] = int(r[1])
            d['IEXPORT'] = int(r[2])
            d['IPOSWEL' ] = int(r[3])
            d['ISCEN'   ] = int(r[4])
            if version>=3:
                if len(r)>5:
                    d['IBDG'] = int(r[5])
                if len(r)>6:
                    d['MINKD'] = float(r[6])
                if len(r)>7:
                    d['MINC'] = float(r[7])
            
            # 4
            r = f.readline().replace(',',' ').split()
            d['OUTER'  ] = int(r[0])
            d['INNER'   ] = int(r[1])
            d['HCLOSE'  ] = float(r[2])
            d['QCLOSE'  ] = float(r[3])
            d['RELAX'   ] = float(r[4])
            if version>=3:
                if len(r)>5:
                    if d['OUTER']>=0:
                        d['NPCOND'] = float(r[5])
                    else:
                        d['PARTOPT'] = int(r[5])
                        if d['PARTOPT']==1:
                            d['LOADFILE'] = self._read_path(f.readline().strip())
                if len(r)>6:
                    if d['OUTER']>=0:
                        d['MAXWBALERROR'] = float(r[6]) 
                    else:
                        d['IDFMERGE '] = int(r[6])
            
            #5
            r = f.readline().replace(',',' ').split()
            if (version<3 or d['NSCL']==1):
                if d['NMULT']<=1:
                    # 5a
                    d['XMIN'    ] = float(r[0])
                    d['YMIN'    ] = float(r[1])
                    d['XMAX'    ] = float(r[2])
                    d['YMAX'    ] = float(r[3])
                    d['CSIZE'   ] = float(r[4])
                    d['BUFFER'  ] = float(r[5])
                else:
                    # 5b
                    d['IACT'    ] = int(r[0])
                    d['XMIN'    ] = float(r[1])
                    d['YMIN'    ] = float(r[2])
                    d['XMAX'    ] = float(r[3])
                    d['YMAX'    ] = float(r[4])
                    d['CSIZE'   ] = float(r[5])
                    d['BUFFER'  ] = float(r[6])
                    if len(r)>7:
                        d['CSUB'] = r[7]
            elif d['NSCL']==2:
                if d['NMULT']<=1:
                    # 5a
                    d['XMIN'    ] = float(r[0])
                    d['YMIN'    ] = float(r[1])
                    d['XMAX'    ] = float(r[2])
                    d['YMAX'    ] = float(r[3])
                    d['CSIZE'   ] = float(r[4])
                    d['MAXCSIZE'] = float(r[5])
                    d['BUFFER'  ] = float(r[6])
                else:
                    # 5b
                    d['IACT'    ] = int(r[0])
                    d['XMIN'    ] = float(r[1])
                    d['YMIN'    ] = float(r[2])
                    d['XMAX'    ] = float(r[3])
                    d['YMAX'    ] = float(r[4])
                    d['CSIZE'   ] = float(r[5])
                    d['MAXCSIZE'] = float(r[6])
                    d['BUFFER'  ] = float(r[7])
                    if len(r)>8:
                        d['CSUB'] = r[8]
            
            # 6
            l = f.readline().strip().lower()
            e = 'Unexpected line in runfile {}: {}'.format(fname,l)
            assert l == 'active modules',e
            
            # 7 t/m 24
            if version<3:
                self.modules=['CAP','BND','SHD','KDW','VCW','STO','PWT','ANI',
                         'HFB','WEL','DRN','RIV','EVT','GHB','RCH','OLF',
                         'CHD','ISG']
                for key in self.modules:
                    r = f.readline().replace(',',' ').split()
                    d[key] = self._read_module(r)
                r = f.readline()
            else:
                # Packages should be specified between brackets: “(“and “)” so iMODFLOW can
                # recognize the keyword. It will be used to compare this with the variable KEY in
                # Data Set 10.
                r = f.readline()
                self.modules = []
                while '(' in r and ')' in r:
                    key = r[r.find("(")+1:r.find(")")].upper()
                    self.modules.append(key)
                    r = r.replace(',',' ').split()
                    d[key] = self._read_module(r)
                    r = f.readline()
                        
            # 25
            d['BNDFILE'] = self._read_path(r.strip())
            
            ## MODULES
            l = f.readline().strip().lower()
            e = 'Unexpected line in runfile {}: {}'.format(fname,l)
            assert l == 'modules for each layer',e
            
            if version<3:
                # Packages
                self._read_layers(f,f.readline(),d,key='CAP')
                self._read_layers(f,f.readline(),d,key='BND')
                self._read_layers(f,f.readline(),d,key='SHD')
                self._read_layers(f,f.readline(),d,key='KDW')
                self._read_layers(f,f.readline(),d,key='VCW')
                self._read_layers(f,f.readline(),d,key='STO')
                self._read_layers(f,f.readline(),d,key='PWT')
                self._read_layers(f,f.readline(),d,key='ANI')
                self._read_layers(f,f.readline(),d,key='HFB')
                
                r = f.readline()
            else:
                r = f.readline()
                while '(' in r and ')' in r:
                    self._read_layers(f,r,d)
                    r = f.readline()
                    
            # MODULES
            l = r.strip().lower()
            e = 'Unexpected line in runfile {}: {}'.format(fname,l)
            assert l == 'packages for each layer and stress-period',e
            
            d['KPER'] = pd.Series(0,range(d['NPER']))
            d['DELT'] = pd.Series(np.NaN,range(d['NPER']))
            d['DATE'] = pd.Series('',range(d['NPER']))
            d['ISAVE'] = pd.Series(0,range(d['NPER']))
            
            r = f.readline()
            for iT in range(d['NPER']):
                r = r.replace(',',' ').split()                
                d['KPER' ][iT] = int(r[0])                
                d['DELT' ][iT] = float(r[1])
                d['DATE' ][iT] = r[2]
                d['ISAVE'][iT] = int(r[3])
                if len(r)>4:
                    if 'ISUMSAVE' not in d:
                        d['ISUMSAVE'] = pd.Series(0,range(d['NPER']))
                    d['ISUMSAVE'][iT] = int(r[4])
                if version<3:
                    self._read_layers(f,f.readline(),d,key='WEL',iT=iT)
                    self._read_layers(f,f.readline(),d,key='DRN',iT=iT)
                    self._read_layers(f,f.readline(),d,key='RIV',iT=iT)
                    self._read_layers(f,f.readline(),d,key='EVT',iT=iT)
                    self._read_layers(f,f.readline(),d,key='GHB',iT=iT)
                    self._read_layers(f,f.readline(),d,key='RCH',iT=iT)
                    self._read_layers(f,f.readline(),d,key='OLF',iT=iT)
                    self._read_layers(f,f.readline(),d,key='CHD',iT=iT)
                    self._read_layers(f,f.readline(),d,key='ISG',iT=iT)
                    r = f.readline()
                else:
                    r = f.readline()
                    while '(' in r and ')' in r:
                        self._read_layers(f,r,d,iT=iT)
                        r = f.readline()
        self.data = d
    
    def _read_module(self,r):
        """Interpret a line in the runfile with the module information."""
        rd = {}
        rd['IPM'] = int(r[0])
        NLSAVE = int(r[1])
        if NLSAVE>0:
            rd['ILSAVE']=[int(l) for l in r[2:2+NLSAVE]]
        return rd
                
    def _read_layers(self,f,r,data,key=None,iT=None):
        """Read the layer information in an iMod RunFile."""
        if key is None:
            key = r[r.find("(")+1:r.find(")")].upper()
        r = r.replace(',',' ').split()
        NFILES=int(r[0])
        if NFILES==-1:
            # use values from previous run, so do nothing (like FloPy)
            pass
        else:
            n = self._get_n(key)
            ls=[]
            for _ in range(NFILES*n):
                r = f.readline().replace(',',' ').split()
                if len(r)==1:
                    ls.append({'FNAME':self._read_path(r[0])})
                elif len(r)==2:
                    ls.append({'ILAY':int(float(r[0])),'FNAME':self._read_path(r[1])})
                elif len(r)==3:
                    ls.append({'FCT':float(r[0]),'IMP':float(r[1]),'FNAME':self._read_path(r[2])})
                else:
                    ls.append({'ILAY':int(float(r[0])),'FCT':float(r[1]),'IMP':float(r[2]),'FNAME':self._read_path(r[3])})
            df = pd.DataFrame(ls)
            if iT is None:
                data[key]['files'] = df
            else:
                if 'files' not in data[key]:
                    data[key]['files'] = {}
                data[key]['files'][iT] = df
    
    def _get_n(self,key):
        """Determine how many types of files determine each package."""
        if key in ['PWT']:
            n=6
        elif key in ['IBS','RIV']:
            n=4
        elif key in ['EVT']:
            n=3
        elif key in ['ANI','DRN','GHB']:
            n=2
        else:
            n=1
        return n
    
    def _read_path(self,r):
        """Read the path of a filename, dropping the data-path-string."""
        if self.data_path:
            return r.replace(self.data_path,'')
        else:
            return r
    
    def write(self,fname,data=None,data_path=None):
        """Write a runfile to file-name fname."""
        if data is None:
            d = self.data
        else:
            d = data
        if data_path is None:
            data_path = self.data_path
        with open(fname, 'w') as f:
            f.write('{}\n'.format(os.path.join(data_path,d['OUTPUTDIRECTORY'])))
            if self.version<3:
                f.write('{},{},{},{},{},{},{},{}\n'.format(d['NLAY'],d['MXNLAY'],d['NPER'],d['ISS'],d['ISCL'],d['IFTEST'],d['ICONCHK'],d['IIPF']))
            else:
                f.write('{},{},{},{},{},{},{},{}'.format(d['NLAY'],d['MXNLAY'],d['NPER'],d['SDATE'],d['NSCL'],d['IFTEST'],d['ICONCHK'],d['IIPF']))
                if 'IUNCONF' in d:
                    f.write(',{}'.format(d['IUNCONF']))
                    if 'IFVDL' in d:
                        f.write(',{}'.format(d['IFVDL']))
                        if 'IARMWP' in d:
                            f.write(',{}'.format(d['IARMWP']))
                f.write('\n')
            assert d['IIPF']==0, NotImplementedError('IIPF>0 niet ondersteund')
            f.write('{},{},{},{},{}'.format(d['NMULT'],d['IDEBUG'],d['IEXPORT'],d['IPOSWEL'],d['ISCEN']))
            if self.version>3:
                if 'IBDG' in d:
                    f.write(',{}'.format(d['IBDG']))
                    if 'MINKD' in d:
                        f.write(',{}'.format(d['MINKD']))
                        if 'MINC' in d:
                            f.write(',{}'.format(d['MINC']))
            f.write('\n')
            f.write('{},{},{},{},{}'.format(d['OUTER'],d['INNER'],d['HCLOSE'],d['QCLOSE'],d['RELAX']))
            if d['OUTER']>=0:
                if 'NPCOND' in d:
                    f.write(',{}'.format(d['NPCOND']))
                    if 'MAXWBALERROR' in d:
                        f.write(',{}'.format(d['MAXWBALERROR']))
            else:
                if 'PARTOPT' in d:
                    f.write(',{}'.format(d['PARTOPT']))
                    assert d['PARTOPT']!=1, NotImplementedError('PARTOPT==1 niet ondersteund')
                    if 'IDFMERGE' in d:
                        f.write(',{}'.format(d['IDFMERGE']))
            f.write('\n')
            if self.version<3 or d['NSCL']==1:
                if d['NMULT']<=1:
                    f.write('{},{},{},{},{},{}'.format(d['XMIN'],d['YMIN'],d['XMAX'],d['YMAX'],d['CSIZE'],d['BUFFER']))
                else:
                    f.write('{},{},{},{},{},{},{}'.format(d['IACT'],d['XMIN'],d['YMIN'],d['XMAX'],d['YMAX'],d['CSIZE'],d['BUFFER']))
            elif d['NSCL']==2:
                if d['NMULT']<=1:
                    f.write('{},{},{},{},{},{},{}'.format(d['XMIN'],d['YMIN'],d['XMAX'],d['YMAX'],d['CSIZE'],d['MAXCSIZE'],d['BUFFER']))
                else:
                    f.write('{},{},{},{},{},{},{}'.format(d['IACT'],d['XMIN'],d['YMIN'],d['XMAX'],d['YMAX'],d['CSIZE'],d['MAXCSIZE'],d['BUFFER']))
            if self.version<3 or d['NSCL']==1 or d['NSCL']==2:
                if 'CSUB' in d:
                    f.write(',{}'.format(d['CSUB']))
                f.write('\n')
            f.write('{}\n'.format('ACTIVE MODULES'))
            for key in self.modules:
                if 'ILSAVE' in d[key]:
                    f.write('{},{},{} ({})\n'.format(d[key]['IPM'],len(d[key]['ILSAVE']),
                            ",".join([str(x) for x in d[key]['ILSAVE']]),key))
                else:
                    f.write('{},{} ({})\n'.format(d[key]['IPM'],0,key))
                        
            f.write('{}\n'.format(os.path.join(data_path,d['BNDFILE'])))
            f.write('{}\n'.format('MODULES FOR EACH LAYER'))
            for key in self.modules:
                if not isinstance(d[key]['files'],dict):
                    self._write_layers(f,d,key,data_path=data_path)
            
            f.write('{}\n'.format('PACKAGES FOR EACH LAYER AND STRESS-PERIOD'))
            for iT in range(d['NPER']):
                f.write('{},{},{},{}'.format(d['KPER'][iT],d['DELT'][iT],d['DATE'][iT],d['ISAVE'][iT]))
                if 'ISUMSAVE' in d:
                    f.write(',{}'.format(d['ISUMSAVE'][iT]))
                f.write('\n')
                for key in self.modules:
                    if isinstance(d[key]['files'],dict):
                        self._write_layers(f,d,key,data_path=data_path,iT=iT)

    def _write_layers(self,f,d,key,data_path='',iT=None):
        """Write the input of the layers from one (timestep of a) package."""
        if iT is None:
            df = d[key]['files']
        else:
            if iT in d[key]['files']:
                df = d[key]['files'][iT]
            else:
                # use input from previous timestep
                f.write('{},({})\n'.format(-1,key))
                return
        # write the header
        n = self._get_n(key)
        if self.version<3:
            f.write('{},{}\n'.format(int(df.shape[0]/n),key))
        else:
            f.write('{},({})\n'.format(int(df.shape[0]/n),key))
        # write the layer-information
        for i,r in df.iterrows():
            if isinstance(r.FNAME,pd.Series) or isinstance(r.FNAME,pd.DataFrame):
                # this probably is para_sim.inp (Series)
                # or one of the other inp-files (DataFrame)
                fname = os.path.join(data_path,r.FNAME.index.name)
                #inp.write(r.FNAME,fname)
            elif util.is_number(r.FNAME):
                fname = r.FNAME
            else:
                fname = os.path.join(data_path,r.FNAME)
            if 'FCT' in r.index and 'IMP' in r.index and not np.isnan(r['FCT']) and not np.isnan(r['IMP']):
                if 'ILAY' in r.index and not np.isnan(r['ILAY']):
                    f.write('{},{},{},{}\n'.format(r.ILAY,r.FCT,r.IMP,fname))
                else:
                    f.write('{},{},{}\n'.format(r.FCT,r.IMP,fname))
            else:
                if 'ILAY' in r.index and not np.isnan(r['ILAY']):
                    f.write('{},{}\n'.format(r.ILAY,fname))
                else:
                    f.write('{}\n'.format(fname))
            
    def run_imodflow(self,model_ws,name,data_path=None,exefile=None,silent=False, use_summerlevel=None, use_winterlevel=None):  
        """Run imodflow in the directory model_ws."""
        # prepare the model directory
        if not os.path.isdir(model_ws):
            os.makedirs(model_ws)
        
        if True:
            # remove the output-directory
            OUTPUTDIRECTORY = os.path.join(model_ws,self.data['OUTPUTDIRECTORY'])
            if os.path.isdir(OUTPUTDIRECTORY):
                shutil.rmtree(OUTPUTDIRECTORY)
            
        if exefile is None:
            if self.version == 4.21:
                shutil.copyfile(os.path.join(data_path,'tools\\fmpich2.dll'),os.path.join(model_ws,'fmpich2.dll'))
                shutil.copyfile(os.path.join(data_path,'tools\\I_accepted_v4_2.txt'),os.path.join(model_ws,'I_accepted_v4_2.txt'))
                shutil.copyfile(os.path.join(data_path,'tools\\mpich2mpi.dll'),os.path.join(model_ws,'mpich2mpi.dll'))
                shutil.copyfile(os.path.join(data_path,'tools\\mpich2nemesis.dll'),os.path.join(model_ws,'mpich2nemesis.dll'))
                exefile = os.path.join(data_path,'tools\\IMODFLOW_V4_2_1_METASWAP_SVN1233_X64R.EXE')
            elif self.version == 4.3:
                shutil.copyfile(os.path.join(data_path,'tools\\fmpich2.dll'),os.path.join(model_ws,'fmpich2.dll'))
                shutil.copyfile(os.path.join(data_path,'tools\\I_accepted_v4_2.txt'),os.path.join(model_ws,'I_accepted_v4_2.txt'))
                shutil.copyfile(os.path.join(data_path,'tools\\mpich2mpi.dll'),os.path.join(model_ws,'mpich2mpi.dll'))
                shutil.copyfile(os.path.join(data_path,'tools\\mpich2nemesis.dll'),os.path.join(model_ws,'mpich2nemesis.dll'))
                exefile = os.path.join(data_path,'tools\\Driver_MetaSWAP_SVN1341_iMOD43_SVN1958.exe')
            elif self.version == 2.21:
                exefile = os.path.join(data_path,'tools\\imodflow_v2.2.1_metaswap.exe')
            else:
                raise(ValueError('Version {} not supported'.format(self.version)))
            
        # copy the iModflow-executable
        file = exefile.split('\\')[-1]
        exefile_new = os.path.join(model_ws,file)
        shutil.copyfile(exefile,exefile_new)

        # write the runfile
        print('Write iModflow Runfile')
        runfile=os.path.join(model_ws,name+'.run')
        if self.version>3:
            # use relative paths
            self.write(runfile,data_path='.\\')
        else:
            # relative paths are not allowed for older versions of iModflow
            self.write(runfile,data_path=model_ws)

        
        # copy or extract the data-files
        if data_path is not None:
            # copy extra files
            print('Copying recharge-files')
            # use files from MeteoBase
            index = self.get_inp_index('mete_grid.inp')
            if index:
                mete_grid = self.get_inp('mete_grid.inp')
                if isinstance(mete_grid,str):
                    mete_grid = inp.read(os.path.join(model_ws,mete_grid))
                for file in mete_grid['precgrid']:
                    util.copy_data_file(data_path,model_ws,file)
                for file in mete_grid['etrefgrid']:
                    util.copy_data_file(data_path,model_ws,file)
                if self.version>3:
                    # use relative paths
                    mete_grid.precgrid = [os.path.join('.\\',x) for x in mete_grid.precgrid]
                    mete_grid.etrefgrid = [os.path.join('.\\',x) for x in mete_grid.etrefgrid]
                else:
                    # relative paths are not allowed for older versions of iModflow
                    mete_grid.precgrid = [os.path.join(model_ws,x) for x in mete_grid.precgrid]
                    mete_grid.etrefgrid = [os.path.join(model_ws,x) for x in mete_grid.etrefgrid]
                #inp.write(mg,os.path.join(model_ws,fname))
                self.update_inp(mete_grid,'mete_grid.inp')
            
            print('Copying model-files. This can take a while.')
            util.copy_data_file(data_path,model_ws,self.data['BNDFILE'])
            for key in self.modules:
                if self.data[key]['IPM']:
                    files = self.data[key]['files']
                    if isinstance(files,dict):
                        # multiple timesteps
                        for iT in files.keys():
                            if not files[iT].empty:
                                for file in files[iT]['FNAME']:
                                    util.copy_data_file(data_path,model_ws,file)
                    else:
                        if not files.empty:
                            for file in files['FNAME']:
                                if isinstance(file,pd.Series) or isinstance(file,pd.DataFrame):
                                    # this probably is para_sim.inp (Series)
                                    # or one of the other inp-files (DataFrame)
                                    fname = os.path.join(model_ws,file.index.name)  
                                    fname = fname.replace('/','\\')
                                    print(fname)
                                    inp.write(file,fname)
                                else:
                                    util.copy_data_file(data_path,model_ws,file)
        if use_summerlevel is not None:            
            shutil.copy(use_summerlevel, os.path.join(model_ws,'OPPERVLAKTEWATER','ZOMER','PEIL_LAAG1_1.IDF'))
            shutil.copy(use_winterlevel, os.path.join(model_ws,'OPPERVLAKTEWATER','WINTER','PEIL_LAAG1_1.IDF'))
            
        # run executable and display output in Console
        print('Run iModflow')
        batch_file = os.path.join(model_ws,name+'.bat')
        with open(batch_file, 'w') as file:
            file.writelines('"{}" "{}"'.format(exefile_new.split('\\')[-1],runfile.split('\\')[-1]))
        proc = sp.Popen(batch_file,stdout=sp.PIPE, stderr=sp.STDOUT, cwd=model_ws)
        while True:
            line = proc.stdout.readline()
            c = line.decode('utf-8')
            if c != '':
                c = c.rstrip('\r\n')
                if not silent:
                    print('{}'.format(c))
            else:
                break
            
    def to_version(self,version):
        """Change the version of a Runfile."""
        from_version = self.version
        d = self.data
        if from_version<3 and version>=3:
            d['SDATE'] = 0
            d['NSCL'] = 1
            d.pop('ISS')
            d.pop('ISCL')
            
        if from_version<4.1 and version>=4.1:
            # ISG definitions now must include FCT and IMP parameter specifications.
            for iT in d['ISG']['files']:
                if 'FCT' not in d['ISG']['files'][iT].columns:
                    d['ISG']['files'][iT]['FCT'] = 1.0
                if 'IMP' not in d['ISG']['files'][iT].columns:
                    d['ISG']['files'][iT]['IMP'] = 1.0
        self.version = version
        
    def update_metaswap(self, mete_grid=True, land_use=True, svat_dtgw=0, svat_per=0, datadir=None, start_date=None, end_date=None, metaswap_vars=None):
        """Update the metaswap to new soil types and grids for meteo."""
        # add storage cofficient
        self.data['STO']['IPM']=1
        self.data['STO']['files'] = pd.DataFrame([{'ILAY':ilay,'FCT':1.0,'IMP':0.0,'FNAME':10**-5} for ilay in range(1,self.data['NLAY']+1)])
        
        self.datadir = datadir
        self.start_date = start_date
        self.end_date = end_date
        
        # add extra variables for metaswap
        ind = np.hstack((np.arange(9),8,9,8,8,8,8,8,8,8,8,8,8,8,np.arange(10,19)))
        df = self.data['CAP']['files'].iloc[ind].copy()
        df.reset_index(inplace=True,drop=True)
        df.loc[8,'FNAME']=100 # ARC mm/d
        df.loc[11,'FNAME']=0.05 # PDU m
        df.loc[12,'FNAME']=0.2 # PDR m
        df.loc[13,'FNAME']=0 # OFU day
        df.loc[14,'FNAME']=0 # OFR day
        df.loc[15,'FNAME']=0 # ONU day
        df.loc[16,'FNAME']=0 # ONR day
        df.loc[17,'FNAME']=1.0 # QIU m/d
        df.loc[18,'FNAME']=1.0 # QIR m/d
        df.loc[19,'FNAME']=-9999 # PWT m+MSL
        df.loc[20,'FNAME']=1.0 # SFC -
        df.loc[21,'FNAME']=1.0 # CFC -
        
        
        fnames = self.get_inp_fnames(df['FNAME'])
        if land_use:
            # replace some files by newer ones from WUR
            df.loc[fnames=='metaswap\\fact_svat.inp','FNAME']=r'metaswap43\\fact_svat.inp'
            df.loc[fnames=='metaswap\\luse_svat.inp','FNAME']=r'metaswap43\\luse_svat.inp'
            # update landuse
            fname_old = 'metaswap\\landgebruik.idf'
            fname = util.read_file_from_datadir(fname_old, datadir=datadir)
            a = imod.idf.open(fname)
            val = a.values.astype(int)
            koppel={1:1,2:2,3:3,4:4,5:5,6:6,8:8,9:9,10:10,11:11,12:12,16:16,18:18,
                    19:18,20:11,21:12,22:11,23:1,24:15,25:1,26:1,35:15,36:14,37:14,
                    38:1,41:13,42:13,43:11,45:1,46:15}
            assert np.all(np.unique(val) == np.array(list(koppel.keys())))
            a.values = util.replace_with_dict(val, koppel)
            fname_new = r'extracted\\landgebruik.idf'
            imod.idf.write(os.path.join(datadir,fname_new), a)
            df.loc[fnames==fname_old,'FNAME']=fname_new
        else:
            fname_old = 'metaswap\\fact_svat.inp'
            fact_svat = inp.read(os.path.join(datadir,fname_old))
            fact_svat.loc[fact_svat['faepdvg']<0.01,'faepdvg'] = 0.01
            fname_new = 'extracted\\fact_svat.inp'
            inp.write(fact_svat, os.path.join(datadir,fname_new))
            df.loc[fnames==fname_old,'FNAME']=fname_new
            
        df.loc[fnames=='metaswap\\unsa_svat.bda','FNAME']='metaswap43\\init_svat.inp'
        
        # set metaswap files in runfile
        self.data['CAP']['files']=df
        
        # add unsa_svat_path to para_sim.inp
        para_sim = self.get_inp('para_sim.inp')
        if isinstance(para_sim, str):
            para_sim = inp.read(para_sim,datadir)
        unsa_svat_path = os.path.join(datadir,'extracted\\bofek_unsa/')
        para_sim['unsa_svat_path']='"{}"'.format(unsa_svat_path)
        para_sim['svat_dtgw'] = svat_dtgw
        para_sim['svat_per'] = svat_per
        self.update_inp(para_sim,'para_sim.inp','extracted\\para_sim.inp')
        
        # update sel_key_svat_per.inp
        sel_key_svat_per = self.get_inp('sel_key_svat_per.inp')
        if isinstance(sel_key_svat_per, str):
            sel_key_svat_per = inp.read(sel_key_svat_per, datadir)
        for var in metaswap_vars:
            sel_key_svat_per.loc[sel_key_svat_per.key==var,'ioptkey'] = 1
        self.update_inp(sel_key_svat_per, 'sel_key_svat_per.inp', 'metaswap/sel_key_svat_per.inp')
            
        # update tiop_sim
        tiop_sim = pd.DataFrame(np.zeros((len(pd.date_range(start_date, end_date)),4)), columns=['td','iy','io','id'], index = range(len(pd.date_range(start_date, end_date))))
        for dd,d in enumerate(pd.date_range(start_date,end_date)):
            tiop_sim.at[dd,'td'] =f"{(d - pd.Timestamp(str(d.year)+'-01-01')).total_seconds()/(3600.*24.):15.2f}"
            tiop_sim.at[dd,'iy'] =f"{d.year:6g}"
            tiop_sim.at[dd,'io'] = f"{7:6g}"
            tiop_sim.at[dd,'ip'] = f" "            
        self.update_inp(tiop_sim,'tiop_sim.inp','metaswap\\tiop_sim.inp')                               
        
        
        if mete_grid:            
            # replace meteo-stations with meteo-grids from MeteoBase
            self.add_mete_grid()
        
    def add_mete_grid(self):
        """Add mete-grid to runfile."""
        #if pd.to_datetime(self.data['DATE'][0])<pd.to_datetime('1990'):
            # first remove the year 1989, as we do not have grids for 1989
            #self.change_period('1990',self.data['DATE'].iloc[-1])
        self.change_period(self.start_date, self.end_date)
        # generate our own mete_grid from MeteoBase
        mete_grid = self.get_mete_grid()
        # replace mete_svat.inp by mete_grid.inp (and maybe remove meteostation.idf?)
        self.update_inp(mete_grid,r'metaswap\mete_svat.inp','extracted\mete_grid.inp')
    
    def get_mete_grid(self,years=None):
        leap = np.arange(1904,2100,4)
        """Get the mete-grid file with links to all the meteo-grids."""
        if years is None:
            dates = pd.to_datetime(self.data['DATE'])
            years = pd.Index(dates).year.unique()
        pathname = os.path.join(self.datadir, 'forcing')    
        preclist = [os.path.join(os.path.join(self.precip_path, 'forcing','precipitation'), f) for f in os.listdir(os.path.join(self.precip_path, 'forcing','precipitation'), )]
        evaplist = [os.path.join(os.path.join(self.evap_path, 'forcing','evaporation'), f) for f in os.listdir(os.path.join(self.evap_path, 'forcing','evaporation'), )]
        if os.path.isdir(pathname):
            fname = os.path.join(pathname,'mete_grid.inp')                      
            mete_grid = []
            for year in years:               
                mg = inp.read(fname)     
                # if year == years[-1]:
                #     mg = mg.iloc[0:366,:]
                #     mg.iy = np.repeat(year,366)
                #     mg.td = np.arange(366.)                    
                # if year in leap:
                #     mg = mg.iloc[0:366,:]
                #     mg.iy = np.repeat(year,366)
                #     mg.td = np.arange(366.)                    
                # else:
                ndays = len([p for p in preclist if pd.to_datetime(os.path.splitext(os.path.basename(p))[0].split('_')[1],format='%Y%m%d').year == year])                    
                mg = mg.iloc[0:ndays,:]
                mg.iy = np.repeat(year,ndays)
                mg.td = np.arange(float(ndays))
                mg.precgrid = [p for p in preclist if pd.to_datetime(os.path.splitext(os.path.basename(p))[0].split('_')[1],format='%Y%m%d').year == year]
                mg.etrefgrid = [e for e in evaplist if pd.to_datetime(os.path.splitext(os.path.basename(e))[0].split('_')[2],format='%Y%m%d').year == year]                
                mete_grid.append(mg)
        else:
            raise(ValueError(f'{pathname} does not exist.'))
        mete_grid = pd.concat(mete_grid).reset_index(drop=True)
        return mete_grid
    
    def change_period(self, start_date, end_date):
        """
        Change the simulation period of a RunFile.
        
        When the end_date is after the last timestep, the data of the last
        timestep is copied to all periods after the last timestep. Only the
        meteo-data from grid-files is updated for the new period.
    
        Parameters
        ----------
        self : RunFile
            an instance of a read Runfile from runfile.py
        start_date : str or datetime
            The first day of the new simulation period
        end_date : str or datetime
            The last day of the new simulation period
    
        """
        start_date = pd.to_datetime(start_date)
        end_date = pd.to_datetime(end_date)
        dates = pd.to_datetime(self.data['DATE'])
        if start_date<dates.iloc[0]:
            raise(ValueError('start date is before start of model input'))
        if end_date>dates.iloc[-1]:
            # raise(ValueError('end date is after end of model input'))
            print('End date is after end of model input. Copying last timestep of input (except for mete_grid.inp)')
            dates_extra = pd.date_range(dates.iloc[-1]+pd.to_timedelta(1,'d'),end_date)
            index = np.arange(self.data['NPER'],self.data['NPER']+len(dates_extra))
            dates = dates.append(pd.Series(dates_extra,index=index))
            self.data['NPER'] = self.data['NPER'] + len(dates_extra)
            self.data['KPER'] = self.data['KPER'].append(pd.Series(index+1,index=index))
            self.data['DELT'] = self.data['DELT'].append(pd.Series(1.0,index=index))
            date = [d.strftime('%Y%m%d') for d in dates_extra]
            self.data['DATE'] = self.data['DATE'].append(pd.Series(date,index=index))
            self.data['ISAVE'] = self.data['ISAVE'].append(pd.Series(1,index=index))
            
            # TODO: for all the packages, copy the last full year off the original input
            # now the last timestep of the original input is used for all new timesteps
            if self.get_inp_index('mete_grid.inp') is not None:
                # add meteorological data
                years = pd.Index(dates).year.unique()
                # get the full meteo-data for the simulation years
                # later this will be trucated to the exact simulation period
                mete_grid = self.get_mete_grid(years)
                self.update_inp(mete_grid, 'mete_grid.inp')
            
            
        # inT = np.where((dates>=start_date) & (dates<=end_date))[0]
        inT = dates.index[(dates>=start_date) & (dates<=end_date)]
        self.data['NPER'] = len(inT)
        self.data['KPER'] = (self.data['KPER'][inT]-inT[0]).reset_index(drop=True)
        self.data['DELT'] = self.data['DELT'][inT].reset_index(drop=True)
        self.data['DATE'] = self.data['DATE'][inT].reset_index(drop=True)
        self.data['ISAVE'] = self.data['ISAVE'][inT].reset_index(drop=True)
        for module in self.modules:
            files = self.data[module]['files']
            # a dictionary is used for data that can change per timestep:
            if isinstance(files,dict):
                # find out which timestep contains the first files
                if inT[0] not in files:
                    for iT in range(inT[0],-1,-1):
                        if iT in files:
                            first_files = files[iT]
                            break
                
                # remove or reindex files
                for iT in list(files):
                    if iT<inT[0] or iT>inT[-1]:
                        # remove the files
                        files.pop(iT)
                    elif inT[0]!=0:
                        # change the timestep
                        files[iT-inT[0]] = files.pop(iT)
                        
                # and set the first set of files
                if 0 not in files:
                    files[0] = first_files
        for fname in self.get_inp_fnames():
            # if 'mete_svat.inp' in fname:
            #     mete_svat_org = self.get_inp('mete_svat')
            #     if isinstance(mete_svat_org,str):
            #         mete_svat_org = inp.read(os.path.join(self.datadir,mete_svat_org))
                
            #     t_org = util.yd_to_t(mete_svat_org['iy'], mete_svat_org['td'])
            #     if not np.any(t_org<=start_date) or not np.any(t_org>=end_date):
            #         print('Downloading KNMI-data')
                    
            #         import pastas as pt
            #         # download precipitation
            #         stns = mete_svat_org['nmme'].unique()
            #         one_day = pd.to_timedelta(1,'d')
            #         knmi_p = pt.read.KnmiStation.download(start=start_date+one_day,
            #                                               end=end_date+one_day,
            #                                               vars='RD', 
            #                                               stns=stns)
                    
            #         # download evaporation for de Bilt
            #         knmi_e = pt.read.KnmiStation.download(start=start_date,
            #                                               end=end_date,
            #                                               vars='EV24',
            #                                               stns=260)
                    
            #         # set 9 and 1 hour to 0
            #         knmi_p.data.index = knmi_p.data.index.normalize()
            #         knmi_e.data.index=knmi_e.data.index.normalize()
                    
            #         # Set the time to the start of the interval (which is what SIMGRO needs)
            #         knmi_p.data.index = knmi_p.data.index - one_day
            #         knmi_e.data.index = knmi_e.data.index - one_day
                    
            #         if set(knmi_p.data.index) != set(knmi_e.data.index):
            #             raise(Exception('De tijdstippen van de neerslag- en de verdmapingsdata komen niet overeen'))
            #         RD = knmi_p.data.pivot(columns='STN',values='RD')
            #         if set(stns) != set(RD.columns):
            #             raise(Exception('Het is niet gelukt voor alle stations data te downloaden'))
            #         if np.any(RD.isna()):
            #             raise(Exception('Er zitten NaNs in de neerslag-data'))
            #         if np.any(knmi_e.data['EV24'].isna()):
            #             raise(Exception('Er zitten NaNs in de verdampings-data'))
            #         td = ((RD.index-util.yd_to_t(RD.index.year,0))/one_day).repeat(len(stns))
            #         iy = RD.index.year.repeat(len(stns))
            #         prec = np.ravel(RD)*1000
            #         etref = knmi_e.data.loc[RD.index.repeat(len(stns)),'EV24']*1000
            #         nmme = list(RD.columns)*RD.shape[0]
            #         mete_svat = pd.DataFrame({'td':td, 'iy':iy, 'prec':prec, 'etref':etref, 'nmme':nmme})
                    
            #         self.update_inp(mete_svat,fname,'extracted\\mete_svat.inp')
            if 'mete_grid.inp' in fname:
                mete_grid = self.get_inp(fname)
                if isinstance(mete_grid,str):
                    mete_grid = inp.read(os.path.join(self.datadir,fname))
                t = pd.to_datetime(mete_grid['iy'].astype(str)) + pd.to_timedelta(mete_grid['td'],unit='d')
                mask = (t>=start_date) & (t<=end_date+pd.to_timedelta(1,'d'))
                mete_grid = mete_grid.loc[mask]
                self.update_inp(mete_grid,fname,'extracted\\mete_grid.inp')
            if 'para_sim.inp' in fname:
                para_sim = self.get_inp(fname)
                if isinstance(para_sim,str):
                    para_sim = inp.read(os.path.join(self.datadir,fname))
                para_sim['IYBG'] = start_date.year
                dt = start_date-pd.to_datetime(str(start_date.year))
                para_sim['IDBG'] = dt.total_seconds()/86400
                self.update_inp(para_sim,fname,'extracted\\para_sim.inp')
    
    def get_inp_fnames(self, FNAME=None):
        """Get a Pandas Series with all the filenames (or constants) of the meteswap-input."""
        if FNAME is None:
            FNAME = self.data['CAP']['files']['FNAME']
        FNAME = FNAME.copy()
        for ind,fname in FNAME.iteritems():
            if isinstance(fname,str):
                pass
            elif isinstance(fname,float) or isinstance(fname,int):
                FNAME.loc[ind] = str(fname)
            elif isinstance(fname,pd.DataFrame) or isinstance(fname,pd.Series):
                FNAME.loc[ind] = fname.index.name
            #elif isinstance(fname,dict):
            #    FNAME.loc[ind] = fname['fname']
            else:
                raise(TypeError(type(fname)))
        return FNAME
    
    def get_inp_index(self,name):
        fnames = self.get_inp_fnames()
        #ind = np.where(FNAME.astype(str).str.contains(name))[0]
        ind = np.where([name in fname for fname in fnames])[0]
        if len(ind)==1:
            return fnames.index[ind[0]]
        else:
            return None
        
    def get_inp(self,name):
        return self.data['CAP']['files']['FNAME'].loc[self.get_inp_index(name)]

    def update_inp(self,inp_data,fname_old,fname_new=None):
        if False:
            if fname_new is None:
                fname_old = fname_new
            # write file and place filename in dataframe
            inp.write(inp_data,os.path.join('data',fname_new))
            mask = self.data['CAP']['files']['FNAME']==fname_old
            self.data['CAP']['files'].at[mask,'FNAME'] = fname_new
        else:
            # place data directly in dataframe and write when running the model
            # does not work yet
            if fname_new is None:
                fname_new = self.get_inp(fname_old)
                if not isinstance(fname_new,str):
                    fname_new = fname_new.index.name
            inp_data.index.name = fname_new
            ind = self.get_inp_index(fname_old)
            self.data['CAP']['files'].at[ind,'FNAME'] = inp_data
    
    def add_inp(self,inp_data,fname):
        if True:
            # write file and place filename in dataframe
            inp.write(inp_data,os.path.join(self.datadir,fname))
            self.data['CAP']['files'] = self.data['CAP']['files'].append({'FNAME':fname},ignore_index=True)
        else:
            # place data directly in dataframe and write when running the model
            # does not work yet
            inp_data.index.name = fname
            self.data['CAP']['files'] = self.data['CAP']['files'].append({'FNAME':inp_data},ignore_index=True)
            
    def update_surface_level(self,data_dir,tif_name,idf_name):
        ds = gdal.Open(os.path.join(data_dir,tif_name))
        ds = rasters.clip_dataset(ds,self.get_extent())
        ds2 = rasters.reproject_dataset(ds,pixel_spacing=self.data['CSIZE'],
                                        GDALResampleAlg=gdal.GRA_Med)
        x,y=rasters.get_xy_mid(ds2)
        mv = rasters.get_values(ds2)
        util.write_idf(x,y,mv,os.path.join(data_dir,idf_name))
        self.data['CAP']['files'].at[5,'FNAME'] = idf_name
    
    def update_land_use(self,data_dir,tif_name,idf_name_lgn,idf_name_uba,
                        wss2lgn=r'data\landgebruik\lgqgis2moz.csv',
                        wss2uba=r'data\landgebruik\lgqgis2verhard.csv',
                        figures=False):
        """Update the landuse with data from the waterschadeschatter."""
        ds = gdal.Open(os.path.join(data_dir,tif_name))
        ds = rasters.clip_dataset(ds,self.get_extent())
        # zet de dataset om naar de benodigde grid-grootte (voor de IDF)
        ds2 = rasters.reproject_dataset(ds,pixel_spacing=self.data['CSIZE'],
                                        GDALResampleAlg=gdal.GRA_Med)
        # read the values
        lg_org = rasters.get_values(ds2)
        lg_org = lg_org.astype(int)
        
        # read the lookup-table of the land-use
        l2m = pd.read_csv(wss2lgn,index_col='LG',names=['LG','MOZ','DESC'],header=0)
        lgn = lg_org.copy()
        for LG in np.unique(lg_org):
            mask = lg_org == LG
            lgn[mask] = l2m.at[LG,'MOZ']
        
        # read the lookup-table of the urban area
        l2v = pd.read_csv(wss2uba,index_col='LG',names=['LG','PV','MOZ','DESC'],header=0)
        uba = lg_org.copy()
        for LG in np.unique(lg_org):
            mask = lg_org == LG
            uba[mask] = l2v.at[LG,'PV']
        
        # write to idf's
        x,y=rasters.get_xy_mid(ds2)
        util.write_idf(x,y,lgn,os.path.join(data_dir,idf_name_lgn))
        self.data['CAP']['files'].at[1,'FNAME'] = idf_name_lgn
        util.write_idf(x,y,uba,os.path.join(data_dir,idf_name_uba))
        if self.version<3:
            self.data['CAP']['files'].at[9,'FNAME'] = idf_name_lgn
        else:
            self.data['CAP']['files'].at[10,'FNAME'] = idf_name_lgn
        
        if figures:
            # get the coordinates of the edges
            X,Y = rasters.get_xy(ds2)
            
            f,ax = plot.model_map(figsize=(10,5.65))
            pcm = ax.pcolormesh(X,Y,lg_org)
            plot.colorbar_outside(pcm,ax=ax)
            f.tight_layout(pad=0.0)
            f.savefig(r'figures\landgebruik_origineel')
            
            f,ax = plot.model_map(figsize=(10,5.7))
            pcm = ax.pcolormesh(X,Y,lgn)
            plot.colorbar_outside(pcm,ax=ax)
            f.tight_layout(pad=0.0)
            f.savefig(r'figures\landgebruik_omgezet')
            
            f,ax = plot.model_map(figsize=(10,5.7))
            pcm = ax.pcolormesh(X,Y,uba)
            plot.colorbar_outside(pcm,ax=ax)
            f.tight_layout(pad=0.0)
            f.savefig(r'figures\oppervlakte_stedelijk_gebied')
            
    def update_soil_physical_unit(self, data_dir, idf_name, figures=False):
        """Update the soil physical unit with data from Bofek."""
        util.extract_changed_files(os.path.join(data_dir,r'bodem\bofek_nl.zip'),
                                   os.path.join(data_dir,r'data\extracted'))

        # read the data
        fname = os.path.join(data_dir,r'extracted\bofek_nl.asc')
        ds = gdal.Open(fname)
        ds = rasters.clip_dataset(ds,self.get_extent())
        val_org = rasters.get_values(ds)
        x,y = rasters.get_xy_mid(ds)
        
        # read koppeltabel bofek2bodem
        fname = os.path.join(data_dir,r'bodem\bofek2bodem.csv')
        b2b = pd.read_csv(fname,index_col='BOFEK',squeeze=True)
        
        if False:
            # make a selection (is allready done by clip_dataset)
            inX = (x>self.data['XMIN']) & (x<self.data['XMAX'])
            inY = (y>self.data['YMIN']) & (y<self.data['YMAX'])
            x=x[inX]
            y=y[inY]
            val_org = val_org[np.ix_(inY,inX)]
        
        # calulate bodemfysische_eenheid from bofek
        val_new = val_org.copy()
        for bofek in np.unique(val_org):
            mask = val_org == bofek
            val_new[mask] = b2b[bofek]
        
        # write as idf
        fname = os.path.join(data_dir,idf_name)
        util.write_idf(x,y,val_new,fname)
        self.data['CAP']['files'].at[3,'FNAME'] = idf_name
        
        if figures:
            f,ax = plot.model_map(figsize=(10,5.7))
            pl = val_org.astype(float)
            pl[pl==-9999]=np.NaN
            pcm = ax.pcolormesh(x,y,pl)
            plot.colorbar_outside(pcm,ax=ax)
            f.tight_layout(pad=0.0)
            f.savefig(r'figures\bofek')
            
            f,ax = plot.model_map(figsize=(10,5.7))
            pcm = ax.pcolormesh(x,y,val_new)
            plot.colorbar_outside(pcm,ax=ax)
            f.tight_layout(pad=0.0)
            f.savefig(r'figures\bodemfysische_eenheid')
        
    def get_extent(self):
        return [self.data['XMIN'],self.data['XMAX'],self.data['YMIN'],self.data['YMAX']]
        
def test():
    """Test some iModflow-runfiles."""
    run = Runfile(r'data\AZURE_1.0.3_2001-2005.run',version=4.0)
    run = Runfile(r'A:\model\Ov_bereg\basis.run',version=4.0)
    run = Runfile('data\hdsr_ns.run',data_path='$DBASE$\\',version=2.21)

    