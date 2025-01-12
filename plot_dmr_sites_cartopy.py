#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 11:53:58 2023

@author: JiM
plots site location for particular Maine DMR participants 
Derived from my plt_emolt_sites.py code on github (rather than Basemap)
Started in 2025 where I needed to read "emolt_sites_extra.csv" lookup table 
"""

minyrs=0
maxyrs=55
crit_years=.1*365*24 # actually the "numpts" needed (typically 0.1 years)
region='Maine' # other options are "all", "Maine", "Mass_DMF_sites", "Downeast", etc
border=0.01 # border in lat/lon coordinates to extend boundary
title='Maine_DMR_ventless_participants_2021-2024'
plt_filename=title+'_map.png'

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
#import shapely.geometry as sgeom
import cartopy
import cartopy.crs as ccrs
import pandas as pd
import numpy as np
import matplotlib as mpl
from adjustText import adjust_text
import conversions

def getsite_latlon(site): # returns lat lon given site code from lookup table
    df=pd.read_csv('emolt_sites_extra.csv')
    if site=='ALL':
        sites=np.unique(df['SITE'].str[:2])# returns initials for each participant
        for i in range(len(sites)):
            if i==0:
                df1=df[df['SITE'].str[:2]==sites[i]]
            else:
                df1=pd.concat([df1,df[df['SITE'].str[:2]==sites[i]]])
        lats=df1['LAT'].values
        lons=df1['LON'].values
    else:
        df1=df[df['SITE']==site]
        lats=df1['LAT'].values[0]
        lons=df1['LON'].values[0]
    return lats,lons

# Get sites not yet loaded into emolt_sites_extra.csv
columns=['SITE','LAT','LON','dep','fn','ln','ventless']
df21=pd.read_csv('2021VTSTemp_header_2021.dat',header=None,names=columns)
df22=pd.read_csv('2022VTSTemp_with_hr_header_2022.dat',header=None,names=columns)
df23=pd.read_csv('2023VTSTemp_header_dd.dat',header=None,names=columns)
#for k in range(len(df24)):
#    [df24['LAT'][k],df24['LON'][k]]=conversions.dm2dd(df24['LAT'].values[k],df24['LON'].values[k])
df24=pd.read_csv('2024VTSTemp_header_dd.dat',header=None,names=columns)
df2=pd.concat([df21,df22,df23,df24])
lats=df2['LAT'].values
lons=df2['LON'].values
'''
df=pd.read_csv('emolt_sites_extra.csv') # has numyears, syr, eyr, and numpts columns derived from "fix_emolt_sites.py"
df.rename(columns={'LAT_DDMM': 'LAT','LON_DDMM': 'LON'}, inplace=True)
if len(region)>0: # get only specified region
    if region=='Maine':
        df=df[(df['LAT']>42.) & (df['LON']<-66.)]
    elif region=='Downeast':
            df=df[(df['LAT']>44.) & (df['LON']<-66.)]    
lats=df['LAT'].values
lons=df['LON'].values
'''
ax = plt.axes(projection=cartopy.crs.PlateCarree())
ax.add_feature(cartopy.feature.LAND,color='lightgray')
ax.add_feature(cartopy.feature.OCEAN,color='white')
ax.add_feature(cartopy.feature.COASTLINE)
ax.add_feature(cartopy.feature.BORDERS, linestyle=':')
gls=ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,color="None")
gls.top_labels=False   # suppress top labels
gls.right_labels=False
ax.set_xlim([-71., -67.])
ax.set_ylim([43., 45.])
#xticks = np.arange(-71, -70.4, .2)
#yticks = np.arange(41.4, 42.5, .2)
#plt.yticks(ticks=yticks)
#plt.xticks(ticks=xticks)
ax.set_title(title)
#plt.suptitle(' spanning '+'%0.f' % min(numyears)+'-'+'%0.f' % max(numyears)+' years',size=10)
#sc=plt.scatter(x=lons,y=lats,s=numyears*2,c='red',alpha=1)#,transform=ccrs.PlateCarree())
fishers=list(np.unique(df2.fn))
texts = []
for k in range(len(fishers)):
    df3=df2[df2['fn'].values==fishers[k]]
    sc=plt.scatter(x=df3['LON'].values,y=df3['LAT'].values,s=20,c='C'+str(k),alpha=1)
    if df3['fn'].values[0] in ['Jordan','Sam','Joe','Travis','Nate','Brandon','Jim']:
        osx=.2
    else:
        osx=.4
    if df3['fn'].values[0] in ['Travis']:#,'Brian']:
        osy=.2
    elif df3['fn'].values[0] in ['Justin']:
        osy=-.2
    else:
        osy=.0
    texts.append(ax.text(df3['LON'].values[0]+osx,df3['LAT'].values[0]+osy,df3['fn'].values[0], color='C'+str(k), 
                                  ha='center',va='bottom', fontweight='bold',fontsize=8))
    '''
    for i, point in df2.iterrows():
            texts.append(ax.text(df2['LON'].values[i],df2['LAT'].values[i],df2['fn'].values[i], color="red", 
                                  ha='center',va='bottom', fontweight='bold',fontsize=8))#, path_effects=[pe.withStroke(linewidth=3,foreground="white")]))
    '''
adjust_text(texts)

plt.show()
plt.savefig(plt_filename)