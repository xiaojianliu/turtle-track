# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 14:11:54 2016

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
import numpy as np
import csv
import math
from matplotlib.path import Path
from Back_forecast_function import get_fvcom
import sys
from scipy import  interpolate
import datetime as dt
from matplotlib.path import Path
import netCDF4
from dateutil.parser import parse
import numpy as np
import math
import pandas as pd
from datetime import datetime, timedelta
from math import radians, cos, sin, atan, sqrt  
from matplotlib.dates import date2num,num2date
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
lon=[-70.011,-70.087]
lat=[41.859,41.767]
x=0.015
x1=0.0115
w=0.6
days=7
endtimes =[dt.datetime(2014,11,21,12,0,0),dt.datetime(2014,11,22,12,0,0)]
starttimes=[dt.datetime(2014,11,21,12,0,0)-timedelta(hours=days*24),dt.datetime(2014,11,22,12,0,0)-timedelta(hours=days*24)]
wind=0
wind_get_type='FVCOM'
FNCL='necscoast_worldvec.dat'
CL=np.genfromtxt(FNCL,names=['lon','lat'])
data = np.genfromtxt('sea.csv',dtype=None,names=['x','y','h'],delimiter=',')
cl1=dict(lon=[],lat=[])
for a in np.arange(len(CL['lon'])):
    if CL['lon'][a]>lon[0]-x and CL['lon'][a]<lon[0]+x and CL['lat'][a]<lat[0]+x and CL['lat'][a]>lat[0]-x:
        cl1['lon'].append(CL['lon'][a])
        cl1['lat'].append(CL['lat'][a])

xx=max(cl1['lon'])-min(cl1['lon'])
yy=max(cl1['lat'])-min(cl1['lat'])
print xx,yy
if xx>yy:
    xxx=np.arange(min(cl1['lon']),max(cl1['lon']),0.00001)
    lo=interpolate.interp1d(cl1['lon'],cl1['lat'],kind='cubic') 
    yyy=lo(xxx)
else:
    yyy=np.arange(min(cl1['lat']),max(cl1['lat']),0.00001)
    lo1=interpolate.interp1d(cl1['lat'],cl1['lon'],kind='cubic')
    xxx=lo1(yyy)
print len(cl1['lon']),min(cl1['lon']),max(cl1['lon'])
print len(cl1['lat']),min(cl1['lat']),max(cl1['lat'])
cl1['lon']=[]
cl1['lat']=[]
cl1['lon']=xxx
cl1['lat']=yyy
print len(xxx),min(xxx),max(xxx)
print len(yyy),min(yyy),max(yyy)

zj=[]
for a in np.arange(len(cl1['lon'])):
    zj.append((cl1['lon'][a]-lon[0])*(cl1['lon'][a]-lon[0])+(cl1['lat'][a]-lat[0])*(cl1['lat'][a]-lat[0]))
zz=np.argmin(np.array(zj))
lon1=cl1['lon'][zz]
lat1=cl1['lat'][zz]

cl2=dict(lon=[],lat=[])
for a in np.arange(len(CL['lon'])):
    if CL['lon'][a]>lon[1]-x1 and CL['lon'][a]<lon[1]+x1 and CL['lat'][a]<lat[1]+x1 and CL['lat'][a]>lat[1]-x1:
        cl2['lon'].append(CL['lon'][a])
        cl2['lat'].append(CL['lat'][a])

xx=max(cl2['lon'])-min(cl2['lon'])
yy=max(cl2['lat'])-min(cl2['lat'])
print xx,yy
if xx>yy:
    xxx=np.arange(min(cl2['lon']),max(cl2['lon']),0.00001)
    lo=interpolate.interp1d(cl2['lon'],cl2['lat'],kind='cubic') 
    yyy=lo(xxx)
else:
    yyy=np.arange(min(cl2['lat']),max(cl2['lat']),0.00001)
    lo1=interpolate.interp1d(cl2['lat'],cl2['lon'],kind='cubic')
    xxx=lo1(yyy)
print len(cl2['lon']),min(cl2['lon']),max(cl2['lon'])
print len(cl2['lat']),min(cl2['lat']),max(cl2['lat'])
cl2['lon']=[]
cl2['lat']=[]
cl2['lon']=xxx
cl2['lat']=yyy
print len(xxx),min(xxx),max(xxx)
print len(yyy),min(yyy),max(yyy)

zj1=[]
for a in np.arange(len(cl2['lon'])):
    zj1.append((cl2['lon'][a]-lon[1])*(cl2['lon'][a]-lon[1])+(cl2['lat'][a]-lat[1])*(cl2['lat'][a]-lat[1]))
zz1=np.argmin(np.array(zj1))
lon2=cl2['lon'][zz1]
lat2=cl2['lat'][zz1]
x=[]
y=[]
h=[]
x=data['x']
y=data['y']
h=data['h']

lon11,lat11,lon12,lat12=calculate_point(lon[0],lat[0],cl1['lon'],cl1['lat'],w)
lloo1=[]
lloo2=[]
llaa1=[]
llaa2=[]
lloo1.append(lon11)
lloo1.append(lon12)
llaa1.append(lat11)
llaa1.append(lat12)
print 'lon11,lat11,lon12,lat12',lon11,lat11,lon12,lat12
lon21,lat21,lon22,lat22=calculate_point(lon[1],lat[1],cl2['lon'],cl2['lat'],w)
lloo2.append(lon21)
lloo2.append(lon22)
llaa2.append(lat21)
llaa2.append(lat22)
llo=[]
lla=[]
llo.append(lloo1)
llo.append(lloo2)
lla.append(llaa1)
lla.append(llaa2)
print 'lon21,lat21,lon22,lat22',lon21,lat21,lon22,lat22
plt.figure()
plt.plot(CL['lon'],CL['lat'])
plt.scatter(lon11,lat11,s=5,color='green')
plt.scatter(lon12,lat12,s=5,color='green')
plt.scatter(lon21,lat21,s=5,color='green')
plt.scatter(lon22,lat22,s=5,color='green')
plt.scatter(lon[0],lat[0],s=5,color='red')
plt.scatter(lon[1],lat[1],s=5,color='red')
#plt.scatter(lon1,lat1,s=5,color='yellow')
#plt.scatter(lon2,lat2,s=5,color='yellow')
plt.plot(cl1['lon'],cl1['lat'])
plt.plot(cl2['lon'],cl2['lat'])
plt.axis([-70.75,-70,41.63,42.12])
for a in np.arange(len(lon)):
    m=(lla[a][1]-lla[a][0])/(llo[a][1]-llo[a][0])
    b=lla[a][0]-m*llo[a][0]
    lox=np.arange(min(llo[a]),max(llo[a]),0.000000001)
    loy=[]
    for v in np.arange(len(lox)):
        loy.append(lox[v]*m+b)
    s=[]
    for v in np.arange(len(lox)):
        s.append(np.sqrt((lox[v]-lon[a])*(lox[v]-lon[a])+(loy[v]-lat[a])*(loy[v]-lat[a])))
    suo=np.argmin(s)
    xxin=lox[suo]
    yxin=loy[suo]
    plt.scatter(xxin,yxin,s=10,color='pink')
    m1=(yxin-lat[a])/(xxin-lon[a])
    b1=yxin-m1*xxin
    if xxin<lon[a]:
        lox1=np.arange(xxin-0.02,xxin,0.00000001)
    else:
        lox1=np.arange(xxin,xxin+0.02,0.00000001)
    loy1=[]
    for v in np.arange(len(lox1)):
        loy1.append(lox1[v]*m1+b1)
    #plt.scatter(lox1,loy1,s=5,color='blue')
    s1=[]
    for v in np.arange(len(lox1)):
        s1.append(np.sqrt((lox1[v]-xxin)*111.111*np.cos(yxin*np.pi/180)*111.111*np.cos(yxin*np.pi/180)*(lox1[v]-xxin)+(loy1[v]-yxin)*111.111*111.111*(loy1[v]-yxin)))
    suo1=np.argmin(abs(np.array(s1)-0.7))
    xxin1=lox1[suo1]
    yxin1=loy1[suo1]
    plt.scatter(xxin1,yxin1,s=5,color='yellow')
    Model='30yr'
    back=dict(lon=[],lat=[],time=[],deep=[]) 
    get_obj =  get_fvcom(Model)
    url_fvcom = get_obj.get_url(starttimes[a],endtimes[a])                
    b_points = get_obj.get_data(url_fvcom)
    try:
                    
        back,windspeed= get_obj.get_track(xxin1,yxin1,-1,starttimes[a],wind,wind_get_type,0)
                
    except:
        continue
    plt.plot(back['lon'],back['lat'],'r-')
    plt.text(back['lon'][0]-0.4,back['lat'][0],'endtime=%s'%(str(endtimes[a])))
    plt.text(back['lon'][-1]-0.4,back['lat'][-1],'starttime=%s'%(str(endtimes[a]-timedelta(hours=len(back['lon'])))))
lo1=dict(lon=[],lat=[])
lo2=dict(lon=[],lat=[])
plt.savefig('aaaaaaa',dpi=800)
plt.show()


