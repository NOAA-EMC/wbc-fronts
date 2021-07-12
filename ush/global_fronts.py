#!/bin/env python
"""
Global Fronts Locator and Hausdorff Metrics
Version 1.1
Todd Spindler
IMSG@NCEP/EMC

Changes
25 July 2019 -- ported to Phase 3 (Mars)
1 June 2020 -- radio killed the video star
24 Aug 2020 -- ported to Hera
13 May 2021 -- refactored code to use xarray instead of netCDF4
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()
import numpy as np
import numpy.ma as ma
import matplotlib.image as image
import matplotlib.colors as colors
from matplotlib.dates import MonthLocator, DateFormatter, DayLocator
from matplotlib.dates import WeekdayLocator, MO, TU, WE, TH, FR, SA, SU
from cartopy import crs
import xarray as xr
from datetime import datetime, timedelta
import seawater
import subprocess
import sqlite3
import pandas as pd
from multiprocessing import Pool
import io
import os, sys
import warnings
from read_navy import read_navo, read_navoceano  # local module
import hausdorff as haus                         # local module

warnings.filterwarnings("ignore")

# global parameters 
imageDir='/work/noaa/marine/erob/soca-diag-eric/' 
dbFile='/home/erobin/fronts/fix/global_fronts.db'
modelDir=[]
# set to old matplotlib defaults
plt.style.use('classic')

# task settings
maxjobs=11            # Number of parallel tasks (python multiprocessing)
WANT_POOL=False        # True will turn on parallel processing
UPDATE_DB=False        # Update the accumulated stats database
WANT_STATS_PLOTS=True # Read the stats database and plot timeseries 
QUANTILE=0.75         # setting for Quantile Hausdorff

# scale the colormap for nonsymmetric colorbars
class MidpointNormalize(colors.Normalize):
     def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
          self.midpoint = midpoint
          colors.Normalize.__init__(self, vmin, vmax, clip)
          
     def __call__(self, value, clip=None):
          result, is_scalar = self.process_value(value)
          (vmin,), _ = self.process_value(self.vmin)
          (vmax,), _ = self.process_value(self.vmax)
          resdat = np.asarray(result.data)
          result = np.ma.array(resdat, mask=result.mask, copy=False)
          x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
          res = np.interp(result, x, y)
          result = np.ma.array(res, mask=result.mask, copy=False)
          if is_scalar:
               result = result[0]
          return result

#------------------------
# initialize the database
#------------------------
def init_db(dbfile,table):
    
        
    # connect to db file
    conn = sqlite3.connect(dbfile,detect_types=sqlite3.PARSE_DECLTYPES,timeout=30.0)
    c = conn.cursor()    
    
    # Create table if needed
    c.execute(
        'CREATE TABLE IF NOT EXISTS '+table+' (date timestamp,'\
                                           + ' forecast int,'\
                                           + ' haus_navy real,'\
                                           + ' modhaus_navy real,'\
                                           + ' quanthaus_navy real,'\
                                           + ' unique(date,forecast))')
                                  
    # Save (commit) the changes
    conn.commit()
                 
    # close the connection
    conn.close()

#-----------------------------------------------
# update the database                           
#-----------------------------------------------
def update_db(region,thedate,fcst,haus_navy,modhaus_navy,quanthaus_navy):
     
    init_db(dbFile,region['db'])    

    # connect to db file
    conn = sqlite3.connect(dbFile,detect_types=sqlite3.PARSE_DECLTYPES,timeout=30.0)
    c = conn.cursor()    

    # Larger example that inserts many records at a time
    update = (thedate,fcst,haus_navy,modhaus_navy,quanthaus_navy)
    
    c.execute('REPLACE INTO '+region['db']+' VALUES (?,?,?,?,?)', update)
    
    # Save (commit) the changes
    conn.commit()
                 
    # close the connection
    conn.close()    
    
#-------------------------------------
# load all of the regions of interest 
#-------------------------------------
def load_regions():
    region={} 
    
    region['gulfstream']={}
    region['gulfstream']['name']='Gulf Stream North Wall'
    region['gulfstream']['source']='NAVO'
    region['gulfstream']['db']='gulfstream'
    region['gulfstream']['lims']=[-85, -40, 20, 50]  # (minlon,maxlon,minlat,maxlat)
    region['gulfstream']['loc']=[12,400,12]          # T,depth,znum
    
    region['loopcurrent']={}
    region['loopcurrent']['name']='Loop Current'
    region['loopcurrent']['source']='NAVOCEANO'
    region['loopcurrent']['db']='loopcurrent'
    region['loopcurrent']['lims']=[-100, -75, 15, 35]
    region['loopcurrent']['loc']=[12,400,12]

    region['kuroshio']={}
    region['kuroshio']['name']='N. Wall Kuroshio'
    region['kuroshio']['source']='NAVOCEANO'
    region['kuroshio']['db']='kuroshio'
    region['kuroshio']['lims']=[115, 162, 25, 45]
    region['kuroshio']['loc']=[16,200,9]
    
    return region

#-----------------------------------------------
# read the global rtofs file                 	
#-----------------------------------------------
def read_model(fcst,validDate,region):

    theDate=validDate-timedelta(fcst/24)  # run date    
    zdepth=region['loc'][2]
    
    if fcst == 0:
        fcst_str='n024'
    else:
        fcst_str='f{:03n}'.format(fcst)

    if not os.path.exists(modelDir):
        print('missing model data for',str(fcst),'hour fcst','('+validDate.strftime('%Y-%m-%d')+')')
        print('looking for',modelDir)
        model=None
    else:
        model={'lon':[],'lat':[]}

#        # open the datasets
 #       nc_sal=xr.open_dataset(f'{modelDir}/rtofs_glo_3dz_{fcst_str}_daily_3zsio.nc',decode_times=True)
  #      nc_tmp=xr.open_dataset(f'{modelDir}/rtofs_glo_3dz_{fcst_str}_daily_3ztio.nc',decode_times=True)
   #     nc_ssh=xr.open_dataset(f'{modelDir}/rtofs_glo_2ds_{fcst_str}_diag.nc',decode_times=True)
    #    nc_u=xr.open_dataset(f'{modelDir}/rtofs_glo_3dz_{fcst_str}_daily_3zuio.nc',decode_times=True)
     #   nc_v=xr.open_dataset(f'{modelDir}/rtofs_glo_3dz_{fcst_str}_daily_3zvio.nc',decode_times=True)
        nc_sal=xr.open_dataset("/work/noaa/marine/erob/soca-diag-eric/2020_used/{}/wtxtbul/rtofs_glo_3dz_n024_daily_3zsio.nc".format(validDate.strftime('%Y%m%d')), decode_times=True)
        nc_tmp=xr.open_dataset("/work/noaa/marine/erob/soca-diag-eric/2020_used/{}/wtxtbul/rtofs_glo_3dz_n024_daily_3ztio.nc".format(validDate.strftime('%Y%m%d')),decode_times=True)
        nc_ssh=xr.open_dataset('/work/noaa/marine/erob/soca-diag-eric/2020_used/{}/wtxtbul/rtofs_glo_2ds_n024_daily_diag.nc'.format(validDate.strftime('%Y%m%d')),decode_times=True)
        nc_u=xr.open_dataset('/work/noaa/marine/erob/soca-diag-eric/2020_used/{}/wtxtbul/rtofs_glo_3dz_n024_daily_3zuio.nc'.format(validDate.strftime('%Y%m%d')) ,decode_times=True)
        nc_v=xr.open_dataset('/work/noaa/marine/erob/soca-diag-eric/2020_used/{}/wtxtbul/rtofs_glo_3dz_n024_daily_3zvio.nc'.format(validDate.strftime('%Y%m%d')),decode_times=True)
        # drop singleton dimensions
        nc_sal=nc_sal.squeeze()
        nc_tmp=nc_tmp.squeeze()
        nc_ssh=nc_ssh.squeeze()
        nc_u=nc_u.squeeze()
        nc_v=nc_v.squeeze()

        [minlon,maxlon,minlat,maxlat]=region['lims']
        
        salin=nc_sal['salinity'][zdepth,].where( \
            (nc_sal.Longitude>=minlon%360) & (nc_sal.Longitude<=maxlon%360) & \
            (nc_sal.Latitude>=minlat) & (nc_sal.Latitude<=maxlat),drop=True)
        temp=nc_tmp['temperature'][zdepth,].where( \
            (nc_tmp.Longitude>=minlon%360) & (nc_tmp.Longitude<=maxlon%360) & \
            (nc_tmp.Latitude>=minlat) & (nc_tmp.Latitude<=maxlat),drop=True)
        sst=nc_tmp['temperature'][0,].where( \
            (nc_tmp.Longitude>=minlon%360) & (nc_tmp.Longitude<=maxlon%360) & \
            (nc_tmp.Latitude>=minlat) & (nc_tmp.Latitude<=maxlat),drop=True)
        ssh=nc_ssh['ssh'].where( \
            (nc_ssh.Longitude>=minlon%360) & (nc_ssh.Longitude<=maxlon%360) & \
            (nc_ssh.Latitude>=minlat) & (nc_ssh.Latitude<=maxlat),drop=True)
        u=nc_u['u'][0,].where( \
            (nc_u.Longitude>=minlon%360) & (nc_u.Longitude<=maxlon%360) & \
            (nc_u.Latitude>=minlat) & (nc_u.Latitude<=maxlat),drop=True)
        v=nc_v['v'][0,].where( \
            (nc_v.Longitude>=minlon%360) & (nc_v.Longitude<=maxlon%360) & \
            (nc_v.Latitude>=minlat) & (nc_v.Latitude<=maxlat),drop=True)
        
        MT=pd.to_datetime(nc_tmp['MT'].values)

        model['lat']=temp['Latitude'].values
        model['lon']=temp['Longitude'].values
                
        model['temp']=seawater.temp(salin.to_masked_array(),temp.to_masked_array(),seawater.pres(region['loc'][1],model['lat']),2000)
        model['sst']=sst.to_masked_array()
        model['ssh']=ssh.to_masked_array()
        model['u']=u.to_masked_array()
        model['v']=v.to_masked_array()
        model['current']=(u.to_masked_array()**2+v.to_masked_array()**2)**0.5
        model['vdate']=validDate
        model['rundate']=theDate
        model['fcst']=fcst
        
    return model
#--------------------------------------------------------------------------
# calculate the Hausdorff metrics           	                           
#   this computes three hausdorff metrics, original, modified and quantile 
#--------------------------------------------------------------------------
def hausdorff_metrics(model,navy,region):
    
    # generate contour for hausdorff and plotting
    navy=navy[region['name'].upper()][model['vdate']]
    lims=region['lims']
        
    #clip to the navy frontal region
    clipmask = ((model['lat']<min(navy['lat'])-1) |
                (model['lat']>max(navy['lat'])+1) |
                (model['lon']<min(navy['lon'])%360-1) |
                (model['lon']>max(navy['lon'])%360+1))

    # Pull out the temperature contour at the specified depth and clip to the navy frontal region
    model['temp'].mask=ma.mask_or(model['temp'].mask,clipmask)    
    CS=plt.contour(model['lon'],model['lat'],model['temp'],[region['loc'][2]])
    
    # sort and select the longest contour segment.  This removes closed eddies.
    [segs]=CS.allsegs
    seglen=np.array([s.shape[0] for s in segs])
    i=np.argsort(seglen).tolist()
    i.reverse()
    model_front=segs[i[0]]  # the longest segment is (hopefully) the main front

    haus_navy=[]
    modhaus_navy=[]
    quanthaus_navy=[]
    navy_front=np.stack((navy['lon']%360,navy['lat'])).T
    haus_navy=haus.hausdorff(model_front,navy_front)
    modhaus_navy=haus.mod_hausdorff(model_front,navy_front)
    quanthaus_navy=haus.quantile_hausdorff(model_front,navy_front,QUANTILE)

    if haus_navy==[]:
        haus_navy=np.nan
    if modhaus_navy==[]:
        modhaus_navy=np.nan
    if quanthaus_navy==[]:
        quanthaus_navy=np.nan
            
    return haus_navy, modhaus_navy, quanthaus_navy
#-----------------------------------------------
# create the maps                           	
#-----------------------------------------------
def plot_map(model,param,navy,region):

    #navy=navy[region['name'].upper()][model['vdate']]
    reg=region['db']    
    lims=region['lims']
    source=region['source']

    # generate contour for hausdorff and plotting
    haus_navy,modhaus_navy,quanthaus_navy=hausdorff_metrics(model,navy,region)
    
    data=model[param]
        
    fig=plt.figure(dpi=150)
    ax=plt.axes(projection=crs.Mercator())
    ax.set_extent(np.array(lims)%360,crs=crs.PlateCarree())
    
    if param=='ssh':
        pcm='bwr'
        rcm='yellow'
        ncm='black'
        ocm='white'
    elif param=='current':
        pcm='viridis'
        rcm='red'
        ncm='black'
        ocm='white'
    else:
        pcm='nipy_spectral'
        rcm='white'
        ncm='black'
        ocm='white'
        
    #m.pcolormesh(x,y,data,cmap=plt.cm.gist_rainbow_r)
    if param=='ssh':
        plt.contourf(model['lon'],model['lat'],data,30,cmap=pcm,
            norm=MidpointNormalize(midpoint=0.,
            vmin=data.min(),
            vmax=data.max()),alpha=0.4,
            transform=crs.PlateCarree())
    else:
        plt.contourf(model['lon'],model['lat'],data,30,cmap=pcm,alpha=0.4,
            transform=crs.PlateCarree())
    cbar=plt.colorbar()
    cbar.ax.tick_params(labelsize='x-small')
    if param=='current':
        rskip=int(np.floor(model['lon'].shape[0]/60))
        cskip=int(np.floor(model['lon'].shape[1]/60))
        plt.quiver(model['lon'][::rskip,::cskip],model['lat'][::rskip,::cskip],
                 model['u'][::rskip,::cskip],model['v'][::rskip,::cskip],
                 units='inches',scale_units='inches',width=0.01,scale=6,
                 color='white',transform=crs.PlateCarree())
    
    #x2,y2=m(navy[region['name'].upper()][model['vdate']]['lon']%360,
    #      navy[region['name'].upper()][model['vdate']]['lat'])
    x2=navy[region['name'].upper()][model['vdate']]['lon']%360
    y2=navy[region['name'].upper()][model['vdate']]['lat']
    if source=='NAVOCEANO':
        plt.plot(x2,y2,'-',color=ncm,label=source,linewidth=2,
        zorder=1,transform=crs.PlateCarree())
    else:
        for key in list(navy.keys()):
            #x2,y2=m(navy[key][model['vdate']]['lon']%360,navy[key][model['vdate']]['lat'])
            x2=navy[key][model['vdate']]['lon']%360
            y2=navy[key][model['vdate']]['lat']
            if key.find('NORTH')>=0:
                plt.plot(x2,y2,'-',color=ncm,label='NAVO',linewidth=2,zorder=1,
                    transform=crs.PlateCarree())
            else:
                plt.plot(x2,y2,'-',color=ncm,linewidth=2,zorder=1,
                    transform=crs.PlateCarree())
                
    CS=plt.contour(model['lon'],model['lat'],model['temp'],[region['loc'][2]],
                 colors=rcm,linestyles='-',linewidths=2,zorder=10,
                 transform=crs.PlateCarree())
    CS.collections[0].set_label('RTOFS')
        
    plt.legend(loc='upper left',fontsize='xx-small',facecolor='lightgrey')
    ax.coastlines()
    gl=ax.gridlines(draw_labels=True)
    gl.xlabels_top=False
    gl.ylabels_right=False
        
    plt.title('Global RTOFS '+region['name']+' Location\n'+ \
        "{:03n}".format(model['fcst'])+'H fcst valid '+model['vdate'].strftime('%B %d, %Y')+ \
        ' model run date '+model['rundate'].strftime('%B %d, %Y')+'\n'+ \
        str(region['loc'][0])+'$^\circ$C Isotherm and '+ \
        str(region['loc'][1])+' m with '+ \
        param.upper(),fontsize='small') 
            
    txt=plt.text(0.50,0.95, \
    	# "{:>16}{:5.2f}".format('Hausdorff = ',haus_navy)+' '+source+'\n'+ \
        f'     Modified Hausdorff: {modhaus_navy:5.2f} {source}\n'+\
        f'{int(QUANTILE*100)}th Quantile Hausdorff: {quanthaus_navy:5.2f} {source}\n',
        fontsize='x-small',horizontalalignment='center',verticalalignment='top',
        multialignment='left',transform=ax.transAxes,color='black',
        bbox=dict(facecolor='white'))
    txt.set_multialignment('right')
    txt.set_family('monospace')

    # add some branding and dates
#    noaa_logo=image.imread('work/noaa/marine/lliu/tools/fronts/logo/NOAA_logo.png')
   # nws_logo=image.imread('work/noaa/marine/lliu/tools/fronts/logo/NWS_logo.png')
 #   fig.figimage(noaa_logo,
  #    yo=fig.get_figheight()*fig.dpi-noaa_logo.shape[0])
   # fig.figimage(nws_logo,
    #  xo=fig.get_figwidth()*fig.get_dpi()-nws_logo.shape[1],
     # yo=fig.get_figheight()*fig.get_dpi()-nws_logo.shape[0])
    plt.annotate('NCEP/EMC Verification Post Processing Product Generation Branch',
      xy=(0.01,0.01),xycoords='figure fraction',
      horizontalalignment='left',fontsize='x-small')
    plt.annotate(f'{datetime.now():%d %b %Y} $on Hera$',
      xy=(0.99,0.01),xycoords='figure fraction',
      horizontalalignment='right',fontsize='x-small')    
    
    plt.savefig(imageDir+'/'+model['vdate'].strftime('%Y%m%d')+'/'+region['db']+'_location_'+param+'_'+"{:03n}".format(model['fcst'])+'.png',dpi=fig.dpi)
    print('saved plot to'.format((imageDir+'/'+model['vdate'].strftime('%Y%m%d')+'/'+region['db']+'_location_'+param+'_'+"{:03n}".format(model['fcst'])+'.png')))
    plt.close()
    
    return

#-----------------------------------------------
# extract stats from sqlite db and plot       	
#-----------------------------------------------
def plot_stats(region):
        
    # set up date formats for plot
    loc = WeekdayLocator(byweekday=MO)
    dateFmt=DateFormatter("%b %Y")

    # get all dates
    dbfile='/home/erobin/fronts/fix/global_fronts.db'
    conn = sqlite3.connect(dbfile,detect_types=sqlite3.PARSE_DECLTYPES,timeout=30.0)
    df=pd.read_sql_query('SELECT date, forecast FROM gulfstream '+\
                         'UNION '+\
                         'SELECT date, forecast FROM loopcurrent '+\
                         'UNION '+\
                         'SELECT date, forecast FROM kuroshio '+\
                         'ORDER BY date', conn)
    conn.close()
    
    # extract missing dates and print out blockDates.js file
    df=df[df.forecast==0]  # only nowcast for this
    df['date'] = pd.to_datetime(df.date)
    df.set_index('date',inplace=True)
    df=df.resample('1D').mean()  # get rid of dupes (what are they doing here?)
    missing_dates=df[np.isnan(df.forecast)].index.strftime('%Y-%m-%d')
    df=df.dropna()
    with open(imageDir+'/blockDates.js','w') as f:
        f.write('var minDate=new Date("{}");\n'.format(df.index.strftime('%Y/%m/%d')[0]))
        f.write('var maxDate=new Date("{}");\n'.format(df.index.strftime('%Y/%m/%d')[-1]))
        f.write('var blockDates=new Object();\n')
        f.write('blockDates={};\n'.format(missing_dates.format()))

    conn = sqlite3.connect(dbfile,detect_types=sqlite3.PARSE_DECLTYPES,timeout=30.0)
    df=pd.read_sql_query(
    	'SELECT DISTINCT date, forecast, modhaus_navy, quanthaus_navy FROM '
        +region['db']+' ORDER BY date AND forecast',conn)                
    conn.close()

    df=df[df.date>=pd.Timestamp.now()-pd.Timedelta(weeks=12)]
    df.set_index('date',inplace=True)
    
    for fcst in df.forecast.unique():
        fig=plt.figure(dpi=150)
        ax=plt.axes()
        df.groupby('forecast').quanthaus_navy.plot(ax=ax,x='date',fontsize='x-small',color='green',alpha=0.5,label='_nolegend_')
        df[df.forecast==fcst].quanthaus_navy.plot(ax=ax,x='date',fontsize='x-small',color='red',linewidth=2,label=str(fcst)+' hr fcst')            
        df[df.forecast==0].quanthaus_navy.plot(ax=ax,x='date',fontsize='x-small',color='blue',linewidth=2,label='nowcast')
        plt.legend(fontsize='xx-small')
        plt.grid(axis='both',which='both')
        ax.xaxis.set_major_formatter(dateFmt)
        ax.xaxis.set_major_locator(MonthLocator())
        ax.xaxis.set_minor_locator(DayLocator(bymonthday=15))
        #fig.autofmt_xdate()
        plt.xlabel('')        
        plt.title(f'{int(QUANTILE*100)}th Quantile Hausdorff for {region["name"]} {fcst:03n}Z forecast',fontsize='small')
        
        # add some branding and dates
        #add_mmab_logos()
 #       noaa_logo=image.imread('work/noaa/marine/lliu/tools/fronts/logo/NOAA_logo.png')
#        nws_logo=image.imread('work/noaa/marine/lliu/tools/fronts/logo/NWS_logo.png')
  #      fig.figimage(noaa_logo,
   #         yo=fig.get_figheight()*fig.dpi-noaa_logo.shape[0])
    #    fig.figimage(nws_logo,
     #       xo=fig.get_figwidth()*fig.dpi-nws_logo.shape[1],
      #      yo=fig.get_figheight()*fig.dpi-nws_logo.shape[0])
        plt.annotate('NCEP/EMC/Verification Post Processing Product Generation Branch',
            xy=(0.01,0.01),xycoords='figure fraction',
            horizontalalignment='left',fontsize='x-small')
        plt.annotate(datetime.now().strftime('%d %b %Y'),
            xy=(0.99,0.01),xycoords='figure fraction',
            horizontalalignment='right',fontsize='x-small')
            
        if not os.path.isdir(imageDir+'/stats'):
            os.makedirs(imageDir+'/stats')
        plt.savefig(imageDir+'/stats/'+region['db']+'_hausdorff_'+"{:03n}".format(fcst)+'.png',dpi=fig.dpi)
        plt.close()
        
#--------------
# do the thing 
#--------------
def process_region(region,fcst,navy,theDate):
    
    # load RTOFS
    model=read_model(fcst,theDate,region)
    print('finished reading model')
    # compute the metrics
    haus_navy,modhaus_navy,quanthaus_navy=hausdorff_metrics(model,navy,region)
    print('finished hausdorff metrics')
    # update the dbase
    if UPDATE_DB:
        update_db(region,model['vdate'],fcst,haus_navy,modhaus_navy,quanthaus_navy)
    print('finished updating DB')
    # create imagedir by date if needed
    if not os.path.isdir(imageDir+'/'+model['vdate'].strftime('%Y%m%d')):
        os.makedirs(imageDir+'/'+model['vdate'].strftime('%Y%m%d'))
    if not os.path.isdir(imageDir+'/stats'):
        os.makedirs(imageDir+'/stats')
    print('starting plots')
    # make pretty pictures
    plot_map(model,'sst',navy,region)
    plot_map(model,'ssh',navy,region)
    plot_map(model,'current',navy,region),
    return

#------------------------------------
# start of main routine              
#------------------------------------
if __name__ == '__main__':
    
    if len(sys.argv)==1:
        n=datetime.now()-timedelta(1)  #yesterday
        theDate=datetime(n.year,n.month,n.day)
    else:
        theDate=datetime.strptime(sys.argv[1],'%Y%m%d')
        
    print('Starting WBC Fronts at',datetime.now(),'for',theDate)

    # set model directory RTOFS (global)
    modelDir='/work/noaa/marine/erob/soca-diag-eric/2020_used/'+theDate.strftime('%Y%m%d')+'/wtxtbul'
        
    regions=load_regions()
    
    # process navy frontal messages
    navo=read_navo(theDate, f'{theDate:%Y}')
    navoceano=read_navoceano(f'{theDate:%Y}')

    if WANT_POOL:
        pool=Pool(processes=maxjobs)
            
    FoundDate=False
    for (reg,region) in list(regions.items()):
        print('region',reg)
        # select the correct Navy source for frontal data
        if reg=='gulfstream':
            navy=navo
        elif reg=='loopcurrent':
            navy=navoceano
        elif reg=='kuroshio':
            navy=navoceano
        elif reg=='azores':
            navy=navoceano

        if theDate not in navy[region['name'].upper()]:
            print(theDate.strftime('%Y%m%d'),'not found -- skipping',reg)
            continue
        else:
            FoundDate=True
        
        # set your forecast range here
        #for fcst in range(0,193,24):
        for fcst in range(0,24,24):  # for testing, only does nowcast
            print('fcst',fcst)
            if WANT_POOL:
                result=pool.apply_async(process_region,(region,fcst,navy,theDate))
            else:
                process_region(region,fcst,navy,theDate)
                
    if WANT_POOL:
        #finish up pool
        print('closing pool')
        pool.close()
        print('waiting for workers to exit')
        pool.join()

    # create stats plots only if some new data has been found
    if WANT_STATS_PLOTS:
        for (reg,region) in list(regions.items()):
            plot_stats(region)
                
