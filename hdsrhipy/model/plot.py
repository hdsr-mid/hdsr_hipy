# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 14:40:29 2018

@author: Artesia
"""
import numpy as np
import os
import matplotlib.pyplot as plt
import imod
#import util
import flopy
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import MultipleLocator
import pandas as pd
#import inp
import matplotlib
import warnings
import xarray

def head_map(rf,model_ws=None,date=None,lay=None,newfig=True):
    """
    Generate a figure with a map of the head

    Parameters
    ----------
    rf : RunFile
        an instance of a read Runfile from runfile.py
    model_ws : str
        The model working directory
    date : str
        A string that represents the date in iModflow-format (yyyymmdd)
    lay : int
        The layer number (one-based) that you want to plot
    newfig : bool
        True to generate a new figure

    Returns
    -------
    f : matplotlib Figure
    ax : matplotlib Axes

    """
    f,ax = get_figure(newfig)
    
    if lay is None:
        lay=1
    if isinstance(rf, flopy.modflow.Modflow):
        ml = rf
        if model_ws is None:
            model_ws = ml.model_ws
        fname = os.path.join(model_ws, '{}.hds'.format(ml.name))
        hds = flopy.utils.binaryfile.HeadFile(fname)
        if date is not None:
            # t = ml.dis.start_datetime + hds.get_times()
            raise(Exception('date not yet supported for flopy-model'))
        head = hds.get_data()[lay-1]
        head[head==ml.bas6.hnoflo] = np.NaN
        extent = ml.modelgrid.extent
        im = plt.imshow(head,extent=extent)
        if newfig:
            ax.axis(extent)
    else:
        # read the head of layer 1 of the last timestep with output
        if date is None:
            # the last timestep with output
            iT = np.where(rf.data['ISAVE'][:rf.data['NPER']])[0][-1]
            date = rf.data['DATE'][iT]
        file = 'head_{}_l{}.idf'.format(date,lay)
        fname = os.path.join(model_ws,rf.data['OUTPUTDIRECTORY'],'head',file)
        head = imod.idf.open(fname)
        im = xarray.plot.imshow(head[0,0], add_colorbar=False)
        if newfig:
            ax.axis('equal')
    if newfig:
        colorbar_outside(im)
        f.tight_layout(pad=0.0)
    return f,ax

def head_diff_map(rf1,rf2,model_ws1,model_ws2,date=None,lay=None,newfig=True,
                  threshold=0.01):
    """
    Generate a figure with a map of the difference in head between two
    simulations. The figure show the head in scenario 2 minus the head in 
    scenario 1.

    Parameters
    ----------
    rf1 : RunFile
        an instance of a read Runfile from runfile.py of sceanrio 1 (the reference)
    rf2 : RunFile
        an instance of a read Runfile from runfile.py of sceanrio 2
    model_ws1 : str
        The model working directory of sceanrio 1 (the reference)
    model_ws2 : str
        The model working directory of scenario 2
    date : str
        A string that represents the date in iModflow-format (yyyymmdd)
    lay : int
        The layer number (one-based) that you want to plot
    newfig : bool
        True to generate a new figure
    threshold : float
        Make absolute values smaller than threshold transparent (NaN)

    Returns
    -------
    f : matplotlib Figure
    ax : matplotlib Axes

    """
    f,ax = get_figure(newfig)
    if date is None:
        # the last timestep with output
        iT = np.where(rf1.data['ISAVE'][:rf1.data['NPER']])[0][-1]
        date = rf1.data['DATE'][iT]
    if lay is None:
        lay=1
    # read the head of layer 1 of the last timestep with output
    file = 'head_{}_l{}.idf'.format(date,lay)
    fname = os.path.join(model_ws1,rf1.data['OUTPUTDIRECTORY'],'head',file)
    head1 = imod.idf.open(fname)
    fname = os.path.join(model_ws2,rf2.data['OUTPUTDIRECTORY'],'head',file)
    head2 = imod.idf.open(fname)
    assert np.array_equal(head1.x, head2.x) and np.array_equal(head1.y, head2.y), 'Dimensions do not agree'
    dh = head2 - head1
    dh = dh.where(np.abs(dh)>=threshold,np.NaN)
    im = dh[0,0].plot.imshow(add_colorbar=False)
    if newfig:
        ax.set_title('dh')
        ax.axis('equal')
        colorbar_outside(im)
        f.tight_layout(pad=0.0)
    return f,ax
    
    
def head_cross_section(rf,model_ws,date=None,line=None,newfig=True):
    """
    Generate a figure with a cross-section of the head

    Parameters
    ----------
    rf : RunFile
        an instance of a read Runfile from runfile.py
    model_ws : str
        The model working directory
    date : str
        A string that represents the date in iModflow-format (yyyymmdd)
    line : list
        list of vertices (x,y) along which the corss-section is drawn
    newfig : bool
        True to generate a new figure

    Returns
    -------
    f : matplotlib Figure
    ax : matplotlib Axes

    """
    f,ax = get_figure(newfig)
    if date is None:
        # the last timestep with output
        iT = np.where(rf.data['ISAVE'][:rf.data['NPER']])[0][-1]
        date = rf.data['DATE'][iT]
    # make a flopy model, to plot a cross-section
    m = util.to_flopy(rf, packages=['dis'])
    file = 'head_{}_*.idf'.format(date)
    path = os.path.join(model_ws,rf.data['OUTPUTDIRECTORY'],'head',file)
    head = imod.idf.open(path)[0].values   
    if line is None:
        #line = {'row':int(m.dis.nrow/2)}
        ymean = (rf.data['YMIN'] + rf.data['YMAX'])/2 + rf.data['CSIZE']/2
        line=[(rf.data['XMIN'],ymean),(rf.data['XMAX'],ymean)]
    cs=flopy.plot.PlotCrossSection(model=m,line={'line':line})
    cs.plot_grid()
    pcm = cs.plot_array(head,head=head)
    if newfig:
        plt.title('Run succeeded!\n{}, {}'.format(date,str(line).replace('{','').replace('}','')))
        colorbar_outside(pcm)
        f.tight_layout(pad=0.0)
    return f,ax

def head_time_series(rf,model_ws,x=None,y=None,lays=1,newfig=True, **kwargs):
    """
    Generate a figure with a time series of the head at a certain x and y

    Parameters
    ----------
    rf : RunFile
        an instance of a read Runfile from runfile.py
    model_ws : str
        The model working directory
    x : float
        The x-coordinate you want to plot the head at
    y : float
        The y-coordinate you want to plot the head at
    lays : int of list
        The layer(s) you want to plot the head of
    newfig : bool
        True to generate a new figure

    Returns
    -------
    f : matplotlib Figure
    ax : matplotlib Axes

    """
    f,ax = get_figure(newfig)
    if x is None:
        x = np.mean([rf.data['XMIN'],rf.data['XMAX']])
    if y is None:
        y = np.mean([rf.data['YMIN'],rf.data['YMAX']])
    if not isinstance(lays,list):
        lays=[lays]
    for lay in lays:
        s = util.get_head_time_series(rf,model_ws,x,y,lay)
        add_legend=False
        if 'label' in kwargs:
            label = kwargs.pop('label')
        else:
            if len(lays)>1:
                add_legend=True
                label = 'Layer {}'.format(lay)
            else:
                label = ''
        s.plot(ax=ax,label=label,**kwargs)
        if add_legend:
            ax.legend()
        if newfig:
            plt.title('Head at x={} and y={}'.format(x,y))
            
    return f,ax

def get_figure(newfig):
    if newfig:
        f,ax = plt.subplots()
    else:
        ax = plt.gca()
        f = ax.figure
    return f,ax
	
def colorbar_inside(mappable=None, ax=None, width=0.2, height="90%", loc=5, **kw):
    """"Add a colorbar inside an ax"""
    if ax is None:
        ax = plt.gca()
    cax = inset_axes(ax, width=width, height=height, loc=loc)
    cb = plt.colorbar(mappable, cax=cax, ax=ax, **kw)
    if loc == 1 or loc == 4 or loc == 5:
        cax.yaxis.tick_left()
        cax.yaxis.set_label_position("left")
    return cb

def colorbar_outside(mappable=None, ax=None,**kwargs):
    """"Add a colorbar just outside an ax, tighter than the normal colorbar"""
    ca = plt.gca()
    if ax is None:
        ax = ca
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.05)
    plt.sca(ca)
    cb = plt.colorbar(mappable, cax = cax, **kwargs)
    return cb

def model_map(figsize=(10,6),extent=[105000.,173000.,433000.,473000.],base=10000,
              nrows=1, ncols=1, **kwargs):
    """"Make a map of the HYDROMEDAH model area, to plot properties in"""
    f,axes = plt.subplots(figsize=figsize, nrows=nrows, ncols=ncols, **kwargs)
    if nrows==1 and ncols==1:
        axes = np.array([axes])
    for ax in axes.flat:
        ax.axis('scaled')
        ax.axis(extent)
        ax.xaxis.set_major_locator(MultipleLocator(base))
        ax.yaxis.set_major_locator(MultipleLocator(base))
        rotate_yticklabels(ax)
        ax.grid()
    f.tight_layout(pad=0.0)
    
    if nrows==1 and ncols==1:
        axes = axes[0]
    return f,axes

def title_inside(title, ax=None, x=0.5, y=0.98,
                 horizontalalignment='center', verticalalignment='top', **kwargs):
    """"Place a title inside a matplotlib axes, at the top"""
    if ax is None:
        ax = plt.gca()
    return ax.text(x,y,title, horizontalalignment=horizontalalignment,
                   verticalalignment=verticalalignment, transform=ax.transAxes,
                   **kwargs)

def bda(var, rf, model_ws, it='mean', bda='svat_dtgw.bda', ax=None):
    
    data, key = util.get_bda_data(rf, model_ws, bda, var=var)
    

    if it=='mean':
        data = data.mean(axis=0)
    else:
        data = data[it]
    return bda_array(data, key, rf, model_ws, ax=ax)

def bda_array(data, key, rf, model_ws, ax=None, **kwargs):
    s = pd.Series(data, index=key['PTS'].index, name='values')
    
    # read mod2svat
    fname = os.path.join(model_ws,rf.data['OUTPUTDIRECTORY'],'metaswap',
                         'mod2svat.inp')
    mod2svat = inp.read(fname).set_index('svat')
    
    df = pd.concat((s,key['PTS']['area (m2)'],mod2svat['modfcell']),axis=1)
    
    s = util.fast_weighted_average(df,'values','area (m2)','modfcell')
    
    ncol = int((rf.data['XMAX']-rf.data['XMIN']) / rf.data['CSIZE'])
    nrow = int((rf.data['YMAX']-rf.data['YMIN']) / rf.data['CSIZE'])
    array = np.full((nrow, ncol), np.NaN)
    
    fname = os.path.join(model_ws,rf.data['OUTPUTDIRECTORY'],'mf2005_tmp',
                         '{}.dxc'.format(os.path.split(model_ws)[-1]))
    lrcs = util.read_dxc(fname)
    
    array[lrcs[:,1]-1, lrcs[:,2]-1] = s.loc[lrcs[:,3]]
    if ax is None:
        f,ax = model_map();
    im = plt.imshow(array, extent=rf.get_extent(), **kwargs)
    return im

def RdGrBl():
    """A colormap that goes from dark-red, to light-green, to dark blue
    
    This colorbar can be used to plot values that are good in the middle (light
    green), bad when they are high (dark blue) and bad when they are low (dark
    red) """
    #colors = ['r','g','b']
    colors = [(1, 0, 0), (0.5, 1, 0.5), (0, 0, 1)]
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('RdGrBl', colors)
    return cmap
    
def gdf_column(gdf, column, vmin=None, vmax=None, norm=None, ncolors=None, boundaries=None, markersize=300,
               cmap='viridis', log=False, text=True, fontsize=6, fmt='{:0.2f}',nan_color=(1., 1., 1., 1.),
               colorbar_kw={}, ax=None):
    if ax is None:
        ax = plt.gca()
    if not isinstance(column, str):
        gdf = gdf.copy()
        gdf['column'] = column
        column = 'column'
    if nan_color is None:
        gdf = gdf.loc[~np.isnan(gdf[column])]
    if isinstance(cmap, str):
        cmap = matplotlib.cm.get_cmap(cmap)
    if ncolors is not None and boundaries is not None:
        warnings.warn('ncolors is ignored when specifying boundaries')
    if boundaries is not None and (vmin is not None or vmax is not None):
        warnings.warn('vmin and vmax are ignored when specifying boundaries')
    if norm is None:
        if vmin is None:
            vmin = gdf[column].min()
        if vmax is None:
            vmax = gdf[column].max()
        if (ncolors is None) and (boundaries is None):
            if log:
                norm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax)
            else:
                norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
        if (ncolors is None) and (boundaries is None):
            if log:
                norm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax)
            else:
                norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
        else:
            if boundaries is None:
                if log:
                    boundaries = np.logspace(
                        np.log10(vmin), np.log10(vmax), ncolors + 1)
                else:
                    boundaries = np.linspace(vmin, vmax, ncolors + 1)
            norm = matplotlib.colors.BoundaryNorm(boundaries, cmap.N)
    else:
        if ncolors is not None:
            warnings.warn('ncolors is ignored when specifying norm')
        if boundaries is not None:
            warnings.warn('boundaries is ignored when specifying norm')
        if vmin is not None:
            warnings.warn('vmin is ignored when specifying norm')
        if vmax is not None:
            warnings.warn('vmax is ignored when specifying norm')
        if log:
            warnings.warn('log is ignored when specifying norm')
            log = False
    if (colorbar_kw is not None) and log and ('format' not in colorbar_kw):
        colorbar_kw['format'] = matplotlib.ticker.ScalarFormatter()
    if (colorbar_kw is not None) and (log == False) and ('spacing' not in colorbar_kw):
        colorbar_kw['spacing'] = 'proportional'

    facecolor = []
    for x in gdf[column]:
        if np.isnan(x):
            facecolor.append(nan_color)
        else:
            facecolor.append(cmap(norm(x)))
    facecolor = pd.Series(facecolor, index=gdf.index)
    
    gdf.plot(facecolor=facecolor, edgecolor='k',
             ax=ax, markersize=markersize, picker=True)
    if text:
        for name in gdf.index:
            if not np.isnan(gdf.at[name, column]):
                s = fmt.format(gdf.at[name, column])
                color = 'k'
                if np.mean(facecolor[name][:3]) < 0.4:
                    color = 'w'
                ax.text(gdf.at[name, 'geometry'].x, gdf.at[name, 'geometry'].y, s, fontsize=fontsize,
                        color=color, ha='center', va='center', clip_on=True)
    if colorbar_kw is not None:
        # add a colorbar inside the axes, first generate dummy data with the same norm and cmap
        sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        sm._A = []
        colorbar_inside(sm, ax=ax, **colorbar_kw)
        
    return ax

def rotate_yticklabels(ax):
    """Rotate the labels on the y-axis 90 degrees to save space""" 
    yticklabels = ax.yaxis.get_ticklabels()
    plt.setp(yticklabels, rotation=90, verticalalignment='center')