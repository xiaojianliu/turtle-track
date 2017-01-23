# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 12:50:38 2017

@author: xiaojian
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import  Axes3D
from matplotlib import  cm
from matplotlib.ticker import LinearLocator,FormatStrFormatter
import netCDF4
import datetime as dt
import pandas as pd
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import csv
import math
from matplotlib.path import Path
from Back_forecast_function import get_fvcom
import sys
from scipy import  interpolate
from matplotlib.path import Path
from dateutil.parser import parse
from math import radians, cos, sin, atan, sqrt  
from matplotlib.dates import date2num,num2date
from sympy import * 
def sh_bindata(x, y, z, xbins, ybins):
    """
    Bin irregularly spaced data on a rectangular grid.

    """
    ix=np.digitize(x,xbins)
    iy=np.digitize(y,ybins)
    xb=0.5*(xbins[:-1]+xbins[1:]) # bin x centers
    yb=0.5*(ybins[:-1]+ybins[1:]) # bin y centers
    zb_mean=np.empty((len(xbins)-1,len(ybins)-1),dtype=z.dtype)
    zb_median=np.empty((len(xbins)-1,len(ybins)-1),dtype=z.dtype)
    zb_std=np.empty((len(xbins)-1,len(ybins)-1),dtype=z.dtype)
    zb_num=np.zeros((len(xbins)-1,len(ybins)-1),dtype=int)    
    for iix in range(1,len(xbins)):
        for iiy in range(1,len(ybins)):
#            k=np.where((ix==iix) and (iy==iiy)) # wrong syntax
            k,=np.where((ix==iix) & (iy==iiy))
            zb_mean[iix-1,iiy-1]=np.mean(z[k])
            zb_median[iix-1,iiy-1]=np.median(z[k])
            zb_std[iix-1,iiy-1]=np.std(z[k])
            zb_num[iix-1,iiy-1]=len(z[k])
            
    return xb,yb,zb_mean,zb_median,zb_std,zb_num
def nearest_point( lon, lat, lons, lats, length):  #0.3/5==0.06
    '''Find the nearest point to (lon,lat) from (lons,lats),
       return the nearest-point (lon,lat)
       author: Bingwei'''
    p = Path.circle((lon,lat),radius=length)
    #numpy.vstack(tup):Stack arrays in sequence vertically
    points = np.vstack((lons.flatten(),lats.flatten())).T  
    
    insidep = []
    #collect the points included in Path.
    for i in xrange(len(points)):
        if p.contains_point(points[i]):# .contains_point return 0 or 1
            insidep.append(points[i])  
    # if insidep is null, there is no point in the path.
    if not insidep:
        print 'There is no model-point near the given-point.'
        raise Exception()
    #calculate the distance of every points in insidep to (lon,lat)
    distancelist = []
    for i in insidep:
        ss=math.sqrt((lon-i[0])**2+(lat-i[1])**2)
        distancelist.append(ss)
    # find index of the min-distance
    mindex = np.argmin(distancelist)
    # location the point
    lonp = insidep[mindex][0]; latp = insidep[mindex][1]
        
    return lonp,latp
def calculate_point( lon, lat, lons, lats,xixi):  #0.3/5==0.06
        '''Find the point to (lon,lat) from (lons,lats),
           return the nearest-point (lon,lat)
           author: xiaojian'''
        x=[]
        y=[]        
        for i in np.arange(len(lons)):            
            x.append((lon-lons[i])*(111.111*np.cos(lat*np.pi/180)))
            y.append((lat-lats[i])*111.111)
        #print x
        s=[]
        
        for a in np.arange(len(x)):
            s.append(np.sqrt(x[a]*x[a]+y[a]*y[a]))
            
        ss1=[]
        ii1=[]
        ss2=[]
        ii2=[]
        print 'min(abs(np.array(s)))',min(abs(np.array(s)))
        print 'max(abs(np.array(s)))',max(abs(np.array(s)))
        
        
        for a in np.arange(len(x)):
            if (abs(s[a])-xixi)<0.01 and lons[a]>lon:
                ss1.append(s[a])
                ii1.append(a)
        #print 'ss1',ss1
        if ss1==[]:
            ii1=[]
            for a in np.arange(len(x)):
                if (abs(s[a])-xixi)<0.01 and lats[a]>lat:
                    ss1.append(s[a])
                    ii1.append(a)
            for a in np.arange(len(x)):
                if (abs(s[a])-xixi)<0.01 and lats[a]<lat:
                    ss2.append(s[a])
                    ii2.append(a)
        else:
            for a in np.arange(len(x)):
                if (abs(s[a])-xixi)<0.01 and lons[a]<lon:
                    ss2.append(s[a])
                    ii2.append(a)
        
        #print 'ss1',ss1
        #print 'ss2',ss2
        c1=np.argmin(abs(np.array(ss1)-xixi))
        c2=np.argmin(abs(np.array(ss2)-xixi))
        #print 'c1,c2',c1,c2
        return lons[ii1[c1]],lats[ii1[c1]],lons[ii2[c2]],lats[ii2[c2]]
######## Hard codes ##########
FNCL='necscoast_worldvec.dat'
CL=np.genfromtxt(FNCL,names=['lon','lat'])
tg2014 = np.genfromtxt('julu2015.csv',dtype=None,names=['id','day','time','lat','lon'],delimiter=',')
'''
plt.figure()
plt.scatter(tg2014['lon'],tg2014['lat'],marker='o',color='red',s=5) 
plt. plot(CL['lon'],CL['lat'])
plt.axis([-70.7,-69.95,41.63,42.12])
plt.show()
'''
lon=tg2014['lon']#[-70.011,-70.087]
lat=tg2014['lat']#[41.859,41.767]
dd=0.00003969#0.000324648324
dd=dd*3
day=tg2014['day']
time=tg2014['time']
endtimes=[]
days=7
wind=0
wind_get_type='FVCOM'
Model='30yr'
for a in np.arange(len(day)):
    if day[a][-6:-3]=='Dec':
        #print 1
        endtimes.append(dt.datetime(2015,12,int(day[a][:-7]),int(time[a][:-3]),int(time[a][-2:])))
    if day[a][-6:-3]=='Nov':
        #print 1
        endtimes.append(dt.datetime(2015,11,int(day[a][:-7]),int(time[a][:-3]),int(time[a][-2:])))
    if day[a][-6:-3]=='Jan':
        #print 1
        endtimes.append(dt.datetime(2015,1,int(day[a][:-7]),int(time[a][:-3]),int(time[a][-2:])))
    if day[a][-6:-3]=='Oct':
        #print 1
        endtimes.append(dt.datetime(2015,10,int(day[a][:-7]),int(time[a][:-3]),int(time[a][-2:])))
starttimes=[]
a=[]
b=[]
xx=[]
yy=[]
wei_lon=[]
wei_lat=[]
po=dict(lon=[],lat=[])
plt.figure() 
plt.scatter(lon,lat,marker='o',color='red',s=10)
for a in np.arange(len(endtimes)):
    
    endtimes[a]=endtimes[a]+timedelta(hours=5)
    end=endtimes[a]
    st=endtimes[a]-timedelta(hours=days*24)
    starttimes.append(endtimes[a]-timedelta(hours=days*24))
    lon1,lat1=nearest_point(lon[a],lat[a],CL['lon'],CL['lat'],0.2)
    po['lon'].append(lon1)
    po['lat'].append(lat1)
    a1=(lat1-lat[a])/(lon1-lon[a])
    b1=lat1-a1*lon1
    #a.append(a1)
    #b.append(b1)

    x=symbols('x')
    x12=solve(Eq((1+a1*a1)*x**2+(-2*lon1+2*a1*b1-2*a1*lat1)*x+lon1*lon1+b1*b1+lat1*lat1-2*b1*lat1,dd),x)
    y12=[]
    for aa in np.arange(len(x12)):
        y12.append(x12[aa]*a1+b1)
    
    xx.append(x12)
    yy.append(y12)
    get_obj =  get_fvcom(Model)
    url_fvcom = get_obj.get_url(st,end)                
    b_points = get_obj.get_data(url_fvcom)
    back=dict(lon=[],lat=[],time=[],deep=[]) 
    if a==3 or a==32 or a==40 or a==41 or a==42:
        
        back,windspeed= get_obj.get_track(float(lon[a]),float(lat[a]),-1,st,wind,wind_get_type,0)        
        plt.scatter(lon[a],lat[a],marker='o',color='red',s=10)
        #plt.scatter(float(x12[0]),float(y12[0]),marker='o',color='yellow',s=10)
        plt.plot(back['lon'],back['lat'],'p-')
    if a==9 or a==11 or a==27 or a==20 or a==49:
        continue
        
    elif (x12[0]-lon1>0 and y12[0]-lat1>0 and lon1-lon[a]>0 and lat1-lat[a]>0) or (x12[0]-lon1>0 and y12[0]-lat1<0 and lon1-lon[a]>0 and lat1-lat[a]<0) or (x12[0]-lon1<0 and y12[0]-lat1<0 and lon1-lon[a]<0 and lat1-lat[a]<0) or (x12[0]-lon1<0 and y12[0]-lat1>0 and lon1-lon[a]<0 and lat1-lat[a]>0):
        print 'a,a0',a
        plt.scatter(float(x12[0]),float(y12[0]),marker='o',color='yellow',s=10)
        
        try:   
            
            back,windspeed= get_obj.get_track(float(x12[0]),float(y12[0]),-1,st,wind,wind_get_type,0)        
            plt.scatter(lon[a],lat[a],marker='o',color='red',s=10)
            plt.scatter(float(x12[0]),float(y12[0]),marker='o',color='yellow',s=10)
            plt.plot(back['lon'],back['lat'],'p-')
        except:
            continue
    else:
        print 'a,a1',a
        plt.scatter(float(x12[1]),float(y12[1]),marker='o',color='yellow',s=10)
        
        try:   
            #back=dict(lon=[],lat=[],time=[],deep=[]) 
            back,windspeed= get_obj.get_track(float(x12[1]),float(y12[1]),-1,st,wind,wind_get_type,0)        
            plt.scatter(lon[a],lat[a],marker='o',color='red',s=10)
            plt.scatter(float(x12[1]),float(y12[1]),marker='o',color='yellow',s=10)
            plt.plot(back['lon'],back['lat'],'g-')
        except:
            continue
    plt.scatter(back['lon'][-1],back['lat'][-1],color='green',s=10)
    wei_lon.append(back['lon'][-1])
    wei_lat.append(back['lat'][-1]) 
plt.plot(CL['lon'],CL['lat'],'b-')
plt.axis([-70.75,-69.8,41.63,42.12])
np.save('lon2015',wei_lon)
np.save('lat2015',wei_lat)
