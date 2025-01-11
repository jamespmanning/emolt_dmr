#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 29 06:54:19 2024

@author: user
"""

# routine to read in MaineMDR's "consolidated" temperature data and export as .dat file as follows:
# emolt_site,sn,Ps,datetime,yd,degF,salt,depth(fth),depth_ngdc(fth)  where depth_ngdc is derived from NGDC 
# MS10,5651,01,2013-06-21 08:07:00,    171.34,     43.81,99.999,26.
#
# Jim Manning March 2017
# modified in Dec 2020 to work in Python 3, deal with 2017-2020 data
# this code might be followed by "plt_dmr_sites.py" and "plt_dmr_temps.py"
# Mofified 21 Dec 2020 to read in another file "07_19Sites.csv" to find the serial number
# Modified 22 Dec 2020 to deal with different header input file format in 2020
# Modified 23 Dec 2020 to despike time series relative to a 30pt rolling medium in "clean_time_series" function
# Mofified 29 Dec 2022 to process 2021 data
# Mofified 19 Apr 2023 to process 2022 data
import pandas as pd
from conversions import dd2dm,c2f,m2fth
from datetime import datetime
from datetime import timedelta as td
import dateutil.parser as dparser
import numpy as np
import netCDF4
import os
from matplotlib import pyplot as plt
import warnings
#from pandas.core.common import SettingWithCopyWarning
#warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)
##### FUNCTIONS
def write_html_file(Sc,year,fisherman,sn):
    # writes an html file for this case (or appends to existing for this fishermen)
    #file="/net/pubweb_html/epd/ocean/MainPage/lob/"+fisherman+".html"
    #file=webdir+fisherman+".html"
    file=webdir+Sc[0:2].lower()+".html"
    ex=os.path.exists(file)
    if ex==True:
        with open(file, 'r+') as fd:
            contents = fd.readlines()
            #insert_string='<tr><td></td><td></td><td><a href="'+fisherman+'_'+str(year)+'.png">'+str(year)+'</td></tr>'
            insert_string='<tr><td>'+Sc+'</td><td>'+sn+'</td><td><a href=o'+Sc.lower()+'01.png>'+str(year)+'</td></tr>'
            for index, line in enumerate(contents):
                if '</table>' in line:
                    contents.insert(index -1, insert_string)
                    break
            fd.seek(0)
            fd.writelines(contents)
            fd.close()
    else:
        #fin = open("/net/pubweb_html/epd/ocean/MainPage/lob/template.html", "rt")
        #fout = open("/net/pubweb_html/epd/ocean/MainPage/lob/"+fisherman+'.html', "wt")
        fin = open(webdir+"template.html", "rt")
        #fout = open(webdir+fisherman+'.html', "wt")
        fout = open(webdir+Sc[0:2].lower()+'.html', "wt")
        for line in fin:
            if 'AA' in line:
                fout.write(line.replace('AA', fisherman))
            elif 'XX01' in line:
                #fout.write(line.replace('XX01', Sc).replace('ZZZ', str(int(sn))).replace('first setting','<a href="./fsrs/'+Sc+'_'+fisherman+'.png">'+str(year)))
                fout.write(line.replace('XX01', Sc).replace('ZZZ', str(int(sn))).replace('first setting','<a href=o'+Sc.lower()+'01.png>'+str(year)))
            elif 'next setting?' in line:
                fout.write(line.replace('next setting?', '<a href=o'+Sc.lower()+'.png>'+str(year)))
            else:
                fout.write(line)
        fin.close()
        fout.close()
####################################
def clean_time_series(FFR,threshold,max_bad_expected):
      # function to clean time series record of "threshold" times the standard deviations 
      # various methods of despiking are coded but, in Oct 2020, JiM hardcoded the "rolling_median" method
      # borrowed this from "fsrs2emolt.py".
      # Input "FFR" is the raw dataframe that has one column called "temp" and does a 30 pt rolling window check
      # where it returns the same dataframe with a new columns called "temp_despiked" and "temp_roll_med"
      # Modified by JiM in Jan 2021 to require "max_bad_expected" which might be, for example, the # of hauls that might cause a spike in the time series
      FFR['temp_roll_med']=FFR['temp'].rolling(window=30,center=True).median().fillna(method='bfill').fillna(method='ffill')
      difference=np.abs(FFR['temp']-FFR['temp_roll_med'])
      outlier_idx=difference > threshold
      FFR['temp_despiked']=FFR['temp']
      FFR['temp_despiked'][outlier_idx]=np.nan
      num_spikes=FFR.temp_despiked.isnull().sum()
      if num_spikes>max_bad_expected: # likely not to happen in one season so redo with larger threshold
          print(str(num_spikes)+' spikes removed in 1st pass so threshold doubled')
          threshold_new=threshold*2
          outlier_idx=difference > threshold_new
          FFR['temp_despiked']=FFR['temp']
          FFR['temp_despiked'][outlier_idx]=np.nan
          num_spikes=FFR.temp_despiked.isnull().sum()
          print(str(num_spikes)+' spikes removed in 2nd pass')
      if num_spikes>0:
          print(str(num_spikes)+' spikes removed')
      return FFR,num_spikes
def parse(datet,hr):
     dd=dparser.parse(datet)+td(hr/24)
     return dd
def get_depth(loni,lati,mindist_allowed):
    # routine to get depth (meters) using vol1 from NGDC
    url='https://www.ngdc.noaa.gov/thredds/dodsC/crm/crm_vol1.nc'
    nc = netCDF4.Dataset(url).variables 
    lon=nc['x'][:]
    lat=nc['y'][:]
    xi,yi,min_dist= nearlonlat_zl(lon,lat,loni,lati) 
    if min_dist>mindist_allowed:
      depth=np.nan
    else:
      depth=nc['z'][yi,xi].data
    return float(depth)#,min_dist
def nearlonlat_zl(lon,lat,lonp,latp): # needed for the next function get_FVCOM_bottom_temp 
    """ 
    used in "get_depth"
    """ 
    # approximation for small distance 
    cp=np.cos(latp*np.pi/180.) 
    dx=(lon-lonp)*cp
    dy=lat-latp 
    xi=np.argmin(abs(dx)) 
    yi=np.argmin(abs(dy))
    min_dist=111*np.sqrt(dx[xi]**2+dy[yi]**2)
    return xi,yi,min_dist
#### HARDCODES##############################
dirout=''
webdir='web/'
year=2022 #processing one year at a time
max_bad_expected=5 # expected number of spikes due to hauling
mindist_allowed=0.4# minimum distance from nearest NGDC depth in km
threshold=3. # number of standard deviations to define a spike in the "clean_time_series" function (typically 3.0)

#### MAIN CODE #######
if year==2016:# The following three lines were the input in the MArch 2017 code
    dirout='/net/data5/jmanning/tidbit/dmr/'
    in_outputfile='/net/data5/jmanning/tidbit/dmr/2016VTSConsolidatedTemps'
    input_header_filename='/net/data5/jmanning/dmr/2016_VTS_Tidbits.csv'# site information with columns: Fisherman, Onset SN, Date Initialized, Start Recording Date,Start Recording Time,Recording Frequency,Downloaded,Deployment Date,Deployment Time,Retrieval Date,Retrieval Time,Site,Mean Latitude,Mean Longitude,Mean Depth (fm),Depth Stratum (fm),Comments
elif (year==2017) | (year==2018) | (year==2019):
    in_outputfile='VTStemps_locations07_19BPJ120820'# csv file sent from Blaise on 8 Dec 2020 that already has lat/lon
    input_header_filename='VTSsites07_19 BPJ_020420.xlsx'
    file_with_sn='07_19Sites.csv'
elif year==2020:
    in_outputfile='2020VTSTemp_BPJ_locations'
    input_header_filename='VTS Tidbit Info 2020.csv'
elif year==2021:
    in_outputfile='2021VTSTemp' # contains all info needed
elif year==2022:
    in_outputfile='2022VTSTemp_with_hr' # contains all info needed    
################################################

# read input header files
# df now has ['NOAA_SITE', 'LAT', 'LON', 'DEPTH', 'FISHER', 'VTS_SITE', 'Probe #','YEAR']
# but, in 2020, the info available are ['Mean Latitude','Mean Longitude','Mean Depth (fm)','Fisherman','Onset SN','Deployment Date']
if (year==2016) | (year==2020):#rename columns to be consistent with the 2017-2019 format
    df=pd.read_csv(input_header_filename)
    df.rename(columns = {'Site':'VTS_SITE','Mean Latitude': 'LAT','Mean Longitude': 'LON','Mean Depth (fm)':'DEPTH','Fisherman':'FISHER','Onset SN':'Probe #'},inplace = True)
elif (year==2021) | (year==2022):
    df=pd.read_csv(in_outputfile+'.csv')
    df = df.drop_duplicates('site')# reduces to keep only info needed in header
    df.rename(columns = {'site':'VTS_SITE','mean_lat': 'LAT','mean_lon': 'LON','mean_depth':'DEPTH','captain':'FISHER','serial_number':'Probe #'},inplace=True)#,'Onset SN':'Probe #'},inplace = True)
    df['LAT'] = df.groupby('VTS_SITE')['LAT'].transform('mean')
    df['LON'] = df.groupby('VTS_SITE')['LON'].transform('mean')
    #df['Probe #']=999999
else:
    df=pd.read_excel(input_header_filename) # reads in site/header information with columns: Fisherman, Onset SN, Date Initialized, Start Recording Date,Start Recording Time,Recording Frequency,Downloaded,Deployment Date,Deployment Time,Retrieval Date,Retrieval Time,Site,Mean Latitude,Mean Longitude,Mean Depth (fm),Depth Stratum (fm),Comments


if (year>=2017) and (year<=2019):# read input header files
    df=df[df['YEAR']==year]# delimits header to one year
    dfsn=pd.read_csv(file_with_sn)# has columns "NOAA_SITE,	LAT	,LON,	DEPTH,	FISHER,	VTS_SITE,	Probe #,	YEAR"
    dfsn=dfsn[dfsn['YEAR']==year]# delimits header to one year

# read data file
print('reading '+in_outputfile+' ...')
#dfd=pd.read_csv(in_outputfile+'.csv',parse_dates={'datet'},index_col='datet',date_parser=parse)
dfd=pd.read_csv(in_outputfile+'.csv',parse_dates=[['date','hr']])
dfd['datet']=pd.to_datetime(dfd['date_hr'], format='%m/%d/%Y %H')
dfd=dfd.set_index('datet')
#dfd['datet']=dfd['datet'].dt.datetime
# dfd now has 'ID', 'site', 'year', 'month', 'day', 'yrdy', 'temp', 'lat','lon' with datetime index
if (year>=2017) and (year<=2019):
    dfd=dfd[dfd['year']==year]# data header to one year

# open output files
fout=open(in_outputfile+'_fixed_'+str(year)+'.dat','w')
fheader=open(in_outputfile+'_header_'+str(year)+'.dat','w')
#scodes=['TJ03','TJ04','ID02','SM01','SH01','ID13','SM02','SH02','TO01','BB01','YN01','UP01','TL04','YN02','UP02','TL05']
scodes={'Dustin Delano':'UD','Peter Miller':'ZM','Trevor Jessiman':'RJ','Michael Dawson':'ID','Josh Miller':'SM',\
        'Sam Hyler':'SH','Travis Otis':'TO','BillyBob':'BB','Ryder Noyes':'YN','Justin Papkee':'UP','Terry Lagasse':'TL',\
        'Brian Tripp':'IT','Jordan Drouin':'OD','Joseph Locurto':'JL','Ladd Olson':'LO','Elijah Joyce':'EJ','Joshua Joyce':'JJ',\
            'William Olson':'WO','Tony Bray':'NB','John Mccarthy':'HM'}
if year==2017:    
    inc_start={'Dustin Delano':1,'Peter Miller':1,'Trevor Jessiman':3,'Mike Dawson':14,'Josh Miller':3,\
       'Sam Hyler':3,'Travis Otis':2,'BillyBob':2,'Ryder Noyes':3,'Justin Papkee':1,'Terry Lagasse':6,
       'Brian Tripp':1,'Jordan Drouin':9}# incremental site code for this individual determined by lookin at their sites from previous years using "select first)name,last_name from emoltdbs.emolt_site where site like 'XX%';"
elif year==2018:    
    inc_start={'Dustin Delano':3,'Peter Miller':3,\
        'Travis Otis':2,'BillyBob':2,'Ryder Noyes':3,'Justin Papkee':2,'Terry Lagasse':8,
        'Brian Tripp':3,'Jordan Drouin':9}# incremental site code for this individual determined by lookin at their sites from previous years using "select first)name,last_name from emoltdbs.emolt_site where site like 'XX%';"
elif year==2019:    
    inc_start={'Dustin Delano':5,'Peter Miller':5,'Trevor Jessiman':5,'Mike Dawson':14,\
        'Sam Hyler':7,'Travis Otis':4,'BillyBob':2,'Ryder Noyes':5,'Justin Papkee':4,'Terry Lagasse':10,
        'Brian Tripp':5,'Jordan Drouin':11,'Joe Locurto':1}# incremental site code for this individual determined by lookin at their sites from previous years using "select first)name,last_name from emoltdbs.emolt_site where site like 'XX%';"
elif year==2020:    
    inc_start={'Mike Dawson':16,'Travis Otis':6,'Ryder Noyes':7,'Justin Papkee':5,'Terry Lagasse':12,
        'Jordan Drouin':13,'Joe Locurto':3,'Elijah Joyce':1,'Ladd Olson':1}#r this individual determined by lookin at their sites from previous years using "select first)name,last_name from emoltdbs.emolt_site where site like 'XX%';"
elif year==2021:    
    inc_start={'Jordan Drouin':15, 'Joe Locurto':5, 'Joshua Joyce ':1, 'Tony Bray':1,
           'William Olson':1, 'Mike Dawson':18, 'Terry Lagasse':14,
           'Justin Papkee':7, 'Ryder Noyes':9, 'Ladd Olson':3}#r this individual determined by lookin at their sites from previous years using "select first)name,last_name from emoltdbs.emolt_site where site like 'XX%';"
elif year==2022:    
    inc_start={'Jordan Drouin':17, 'Joseph Locurto':7, 'Joshua Joyce':3,'John Mccarthy':1,
           'Michael Dawson':20, 'Terry Lagasse':16,'Ryder Noyes':11,'Travis Otis':8}#r this individual determined by lookin at their sites from previous years using "select first)name,last_name from emoltdbs.emolt_site where site like 'XX%';"
    
# loop through each line in the header file w/fishermen's name  and generate lines in the output header file
fn,ln=[],[]

df['LON'].astype('float64')# convertingstring to float
df['LAT'].astype('float64')
if df['LON'].values[0]>0: # found that some longitudes were documented as positive
    df['LON']=-1*df['LON']
    
dep_ngdc,diffdep=[],[]
for k in range(df.VTS_SITE.count()):
  if np.isnan(df['LAT'].values[k]): # needed this in the case of one record in 2018 where Justin Papkee had no lat/lon
        continue
  if year>=2021: # position already in dd format
      lat=df['LAT'].values[k]
      lon=df['LON'].values[k]
  else:
      [lat,lon]=dd2dm(df['LAT'].values[k],df['LON'].values[k])
  if df['FISHER'].values[k]=='Jordan Drouin': # needed special case here to change threshold
      threshold=1.0
      print('changed threshold to 1.0 in Jordan case')
  else:
      threshold=3.0
  if df['FISHER'].values[k]=='BillyBob': 
     fn='Billy';ln='Bob'
  else:
     [fn,ln]=df['FISHER'].values[k].split()# gets first_name and last_name
  inc=inc_start[df['FISHER'].values[k]]+len(list(np.where(df['FISHER'].values[0:k]==df['FISHER'].values[k]))[0])# increments the emolt_site code 
  print(scodes[df['FISHER'].values[k]],inc,fn,ln,'%0.2f'%lat,'%0.2f'%lon,df['DEPTH'].values[k]) # writing these out in case we already have a site for him
  #Sc=raw_input('4-digit eMOLT sitecode? ')
  Sc=scodes[df['FISHER'].values[k]]+str(inc).zfill(2) # gets site code from the dictionary hardcoded above and adds a two digit consecutive number w/leading zero
  #Ps=raw_input('Probe setting?')
  Ps='01'
  site=str(int(df['VTS_SITE'].values[k])) # this is DMR site Ventless Trap site code
  if (df['Probe #'].values[k]=='nan') | (np.isnan(df['Probe #'].values[k])):
      #SN=dfsn['Probe #'][np.where(dfsn['VTS_SITE']==int(site))[0]].values[0] # this seemed to work in 2018 but not 2019????
      SN=dfsn['Probe #'].loc[dfsn['VTS_SITE']==int(site)].values[0]
      if type(SN)==float:
          SN='nan' # needs to be a string when writing it to the file
  else:
      SN=str(int(df['Probe #'].values[k]))
  sal='99.999'
  dep=df['DEPTH'].values[k]
  
  if year<2021:
      depth_ngdc=-1*m2fth(get_depth(df['LON'].values[k],df['LAT'].values[k],mindist_allowed))
      print(dep,'%0.1f'%depth_ngdc)
      dep_ngdc.append(depth_ngdc) # save this for plotting later
      if (~np.isnan(depth_ngdc)) and (~np.isinf(depth_ngdc)):
          diffdep.append(np.abs((depth_ngdc-dep)/depth_ngdc))
      else:
          diffdep.append(np.nan)
      if np.isnan(dep): # in this where the header file did not have depth, get it from NGDC given position info
          dep=depth_ngdc 
  
  fheader.write(Sc+","+'%0.2f'%lat+","+'%0.2f'%lon+","+'%0.1f'%dep+","+fn+","+ln+","+site+"\n")

  # select the data for this site, clean it, and plot it
  dfds=dfd.loc[dfd['site']==int(site)]
  dfds,num_spikes=clean_time_series(dfds,threshold,max_bad_expected)
  
  fig, ax = plt.subplots()
  dfds['temp'].plot()
  dfds['temp_despiked'].plot()
  plt.title(fn+' '+ln+' '+Sc+' (VTS='+site+') in '+'%0.1f'%dep+' fathoms with '+str(num_spikes)+' removed.')
  plt.xlabel(year)
  plt.ylabel('degF')
  fig.autofmt_xdate()
  plt.show()
  fig.savefig(webdir+'o'+Sc.lower()+Ps+'.png')
  
  
  # now loop through data and write out to .dat
  for j in range(len(dfds)):
      dtime=dfds.index[j]
      #dtime=datetime.utcfromtimestamp(dfds['datet'].values[j].tolist()/1e9)
      #dtime=dfds['datet'][j]
      #yd=(dtime.replace(tzinfo=None)-datetime(dtime.year,1,1,0,0,0).replace(tzinfo=None)).total_seconds()/86400.
      yd=(dtime-datetime(dtime.year,1,1,0,0,0).replace(tzinfo=None)).total_seconds()/86400.
      fout.write(Sc+","+SN+","+Ps+","+str(dtime)+","+'%0.4f'%(yd)+","+'%0.2f'%dfds['temp_despiked'].values[j]+","+sal+","+'%0.1f'%dep+"\n")
      #fout.write(Sc+","+SN[-5:-1]+","+Ps+","+str(dtime)+","+'%0.4f'%(yd)+","+'%0.2f'%c2f(dfd['temp'][j])[0]+","+sal+","+'%0.1f'%dep+"\n"

  # finally, output the records in the html file for this fisherman
  write_html_file(Sc,year,df['FISHER'].values[k],SN)  # commented out for 2021 case since I didn't havce a "template.html" file handy
  
fout.close()
fheader.close()
if (year==2017) | (year==2020): # these are years we can compare logged depth with NGDC
    fig, ax = plt.subplots()
    if np.inf in diffdep:
        diffdep[diffdep.index(np.inf)]=np.nan
    indices=[i for i, x in enumerate(dep_ngdc) if x >6] # get all the cases where the depth > 6 fths
    diffdep_deep=[diffdep[i] for i in indices]
    print('mean diffdep = ',str(np.nanmean(diffdep)),' max =',str(max(diffdep)))
    print('mean diffdep_deep = ',str(np.nanmean(diffdep_deep)),' max =',str(max(diffdep_deep)))
    #print(diffdep)
    plt.plot(df['DEPTH'].values[0:df.VTS_SITE.count()],dep_ngdc,'r*',markersize=40)
    plt.title('mean percent difference in logged depth vs in NGDC depth in > 6fth = '+"{:#.1f}".format(np.abs(100*np.nanmean(diffdep_deep))))
    plt.ylabel('Fathoms')
    fig.savefig(webdir+'depth_vs_ngdc_'+str(year)+'.png')