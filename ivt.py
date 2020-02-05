# -*- coding: utf-8 -*-
"""
Created on Tue May 21 16:20:53 2019

@author: CoeFamily
"""

import numpy as np
from matplotlib import pyplot
from matplotlib.cm import get_cmap
from cartopy import crs
from cartopy.feature import NaturalEarthFeature
from netCDF4 import Dataset
from wrf import getvar, to_np, get_cartopy, latlon_coords,  interplevel, vinterp
import metpy.calc as mpcalc
from metpy.units import units
import cartopy.feature as cfeature
import os

filenamesgfs = os.listdir("/media/mariofire/4TBExternal/newrf/netcdfdata/gfs_c3c6")
filenamesgfs.sort()

c3gfs = Dataset("/media/mariofire/4TBExternal/newrf/netcdfdata/gfs_c3c6/"+filenamesgfs[0])
c6gfs= Dataset("/media/mariofire/4TBExternal/newrf/netcdfdata/gfs_c3c6/"+filenamesgfs[1])
c3ivt = c3gfs.variables['ivt']
c6ivt = c6gfs.variables['ivt']

ivt = np.zeros((71,99))

#C3dates
#Import the Data files
c3ncfile = [0] * 24
c6ncfile = [0] * 20
#C3dates
c3ncfile[0] = Dataset("H:/reg/wrfout_d01_2006-10-26")
c3ncfile[1] = Dataset("H:/reg/wrfout_d02_2006-10-26")
c3ncfile[2] = Dataset("H:/reg/wrfout_d01_2006-10-27")
c3ncfile[3] = Dataset("H:/reg/wrfout_d02_2006-10-27")

c3ncfile[4] = Dataset("H:/reg/wrfout_d01_2008-04-27")
c3ncfile[5] = Dataset("H:/reg/wrfout_d02_2008-04-27")
c3ncfile[6] = Dataset("H:/reg/wrfout_d01_2008-04-28")
c3ncfile[7] = Dataset("H:/reg/wrfout_d02_2008-04-28")

c3ncfile[8] = Dataset("H:/reg/wrfout_d01_2008-12-10")
c3ncfile[9] = Dataset("H:/reg/wrfout_d02_2008-12-10")
c3ncfile[10] = Dataset("H:/reg/wrfout_d01_2008-12-11")
c3ncfile[11] = Dataset("H:/reg/wrfout_d02_2008-12-11")

c3ncfile[12] = Dataset("H:/reg/wrfout_d01_2002-05-12")
c3ncfile[13] = Dataset("H:/reg/wrfout_d02_2002-05-12")
c3ncfile[14] = Dataset("H:/reg/wrfout_d01_2002-05-13")
c3ncfile[15] = Dataset("H:/reg/wrfout_d02_2002-05-13")

c3ncfile[16] = Dataset("H:/reg/wrfout_d01_2002-10-15")
c3ncfile[17] = Dataset("H:/reg/wrfout_d02_2002-10-15")
c3ncfile[18] = Dataset("H:/reg/wrfout_d01_2002-10-16")
c3ncfile[19] = Dataset("H:/reg/wrfout_d02_2002-10-16")

c3ncfile[20] = Dataset("H:/reg/wrfout_d01_2003-10-27")
c3ncfile[21] = Dataset("H:/reg/wrfout_d02_2003-10-27")
c3ncfile[22] = Dataset("H:/reg/wrfout_d01_2003-10-28")
c3ncfile[23] = Dataset("H:/reg/wrfout_d02_2003-10-28")

#C6dates
c6ncfile[0] = Dataset("H:/reg/wrfout_d01_2004-04-11")
c6ncfile[1] = Dataset("H:/reg/wrfout_d02_2004-04-11")
c6ncfile[2] = Dataset("H:/reg/wrfout_d01_2004-04-12")
c6ncfile[3] = Dataset("H:/reg/wrfout_d02_2004-04-12")

c6ncfile[4] = Dataset("H:/reg/wrfout_d01_2005-11-27")
c6ncfile[5] = Dataset("H:/reg/wrfout_d02_2005-11-27")
c6ncfile[6] = Dataset("H:/reg/wrfout_d01_2005-11-28")
c6ncfile[7] = Dataset("H:/reg/wrfout_d02_2005-11-28")

c6ncfile[8] = Dataset("H:/reg/wrfout_d01_2007-10-21")
c6ncfile[9] = Dataset("H:/reg/wrfout_d02_2007-10-21")
c6ncfile[10] = Dataset("H:/reg/wrfout_d01_2007-10-22")
c6ncfile[11] = Dataset("H:/reg/wrfout_d02_2007-10-22")

c6ncfile[12] = Dataset("H:/reg/wrfout_d01_2008-03-06")
c6ncfile[13] = Dataset("H:/reg/wrfout_d02_2008-03-06")
c6ncfile[14] = Dataset("H:/reg/wrfout_d01_2008-03-07")
c6ncfile[15] = Dataset("H:/reg/wrfout_d02_2008-03-07")

c6ncfile[16] = Dataset("H:/reg/wrfout_d01_2002-12-31")
c6ncfile[17] = Dataset("H:/reg/wrfout_d02_2002-12-31")
c6ncfile[18] = Dataset("H:/reg/wrfout_d01_2003-01-01")
c6ncfile[19] = Dataset("H:/reg/wrfout_d02_2003-01-01")

ii = 0
while ii < 24:
    for k in range(16,25):

        #Use the mixing ratio to compute the specific humidity at each level
        #Compute every 50hPa from 1000hPa to 100 hPa
        w = getvar(c3ncfile[ii], "QVAPOR", timeidx=k)
        p = getvar(c3ncfile[ii],"pressure", timeidx=k)
        #w = vinterp(ncfile,h,'p',[1000,950,900,850,800,750,700,650,600,550,500,450,400,350,300,250,200,150,100])
        
        #Compute the specific humidity using q = (w/w+1)
        q = w / (1+w)
        
        #Get the vector wind V at each level by interpolating the windspeed to each level
        uv = getvar(c3ncfile[ii], 'uvmet', timeidx=k)
        #vspd = vinterp(ncfile,vsp,'p',[1000,950,900,850,800,750,700,650,600,550,500,450,400,350,300,250,200,150,100])
        #Find the IVT using the equation IVT = 1/g * sum(wspd*q*p)/number of levels
        #Make sure to remove all NaN values from the array
        #Create the pressure array based on the size of the other arrays
        
        parray = np.zeros((np.size(q,0), np.size(q,1), np.size(q,2)))
        i = 0
        j = 0
        k = 0
        pnew = 1000
        
        while i < np.size(parray,0):
            while j < np.size(parray,1):
                while k < np.size(parray,2):
                    parray[i,j,k] = pnew
                    k = k + 1
                
                k = 0
                j = j + 1
                
            i = i + 1
            j = 0
            k = 0
            pnew = pnew - 50
        
        #Method 1: Do layer from 100 hPa to 1000 hPa
        newarray = np.sqrt(((1/9.81)*(q[17,:,:]*uv[0,17,:,:]*p[17,:,:] - q[0,:,:] * uv[0,0,:,:]*p[0,:,:]))**2 + ((1/9.81)*(q[17,:,:]*uv[1,17,:,:]*p[17,:,:] - q[0,:,:]*uv[1,0,:,:]*p[0,:,:]))**2)
        ivt = ivt + newarray
    ii = ii + 4
ivt = ivt / 6
#ivt = ivt[:,:] - c3ivt[:,:]
levels = np.arange(0,300,25)
projec = getvar(c3ncfile[0],"pressure")
cart_proj = get_cartopy(projec)
lats, lons = latlon_coords(p)
fig = pyplot.figure(figsize=(10,7.5))
geo_axes = pyplot.axes(projection=cart_proj)
geo_axes.set_extent([-95, -65, 30, 45])
states = NaturalEarthFeature(category = 'cultural',
                       scale = '50m',
                       facecolor = 'none',
                       name = 'admin_1_states_provinces_shp')
geo_axes.add_feature(states,linewidth=1.,edgecolor="black")
geo_axes.coastlines('50m',linewidth=0.8)
geo_axes.gridlines(color="black", linestyle="dotted")
pyplot.contourf(to_np(lons), to_np(lats),
                             to_np(ivt),levels,
                     cmap=get_cmap('jet'),
                             transform=crs.PlateCarree())
pyplot.colorbar(ax=geo_axes, shrink=.86)
pyplot.savefig('c3_day2_dim1_reg_ivt_anom.png')
#ii = 1
#while ii < 16:
#    for k in range(16,25):
#        #Use QVapor for the water vapor mixing ratio
#        #Use the mixing ratio to compute the specific humidity at each level
#        #Compute every 50hPa from 1000hPa to 100 hPa
#        w = getvar(c3ncfile[ii], "QVAPOR", timeidx=k)
#        p = getvar(c3ncfile[ii],"pressure", timeidx=k)
#        #w = vinterp(ncfile,h,'p',[1000,950,900,850,800,750,700,650,600,550,500,450,400,350,300,250,200,150,100])
#        
#        #Compute the specific humidity using q = (w/w+1)
#        q = w / (1+w)
#        
#        #Get the vector wind V at each level by interpolating the windspeed to each level
#        usp = getvar(c3ncfile[ii], 'ua', timeidx=k)
#        #uspd = vinterp(ncfile,usp,'p',[1000,950,900,850,800,750,700,650,600,550,500,450,400,350,300,250,200,150,100])
#        
#        vsp = getvar(c3ncfile[ii], 'va', timeidx=k)
#        
#        wsp = getvar(c3ncfile[ii], 'wspd', timeidx=k)
#        #vspd = vinterp(ncfile,vsp,'p',[1000,950,900,850,800,750,700,650,600,550,500,450,400,350,300,250,200,150,100])
#        #Find the IVT using the equation IVT = 1/g * sum(wspd*q*p)/number of levels
#        #Make sure to remove all NaN values from the array
#        #Create the pressure array based on the size of the other arrays
#        
#        parray = np.zeros((np.size(q,0), np.size(q,1), np.size(q,2)))
#        i = 0
#        j = 0
#        k = 0
#        pnew = 1000
#        
#        while i < np.size(parray,0):
#            while j < np.size(parray,1):
#                while k < np.size(parray,2):
#                    parray[i,j,k] = pnew
#                    k = k + 1
#                
#                k = 0
#                j = j + 1
#                
#            i = i + 1
#            j = 0
#            k = 0
#            pnew = pnew - 50
#        
#        #Method 1: Do layer from 100 hPa to 1000 hPa
#        newarray = (wsp[0,:,:]*q[0,:,:])*p[0,:,:] - (wsp[18,:,:]*q[18,:,:]*p[18,:,:])
#        ivt = ivt + newarray / 9.81
#    ii = ii + 4
#ivt = ivt / 4
#levels = np.arange(0,250,25)
#projec = getvar(c3ncfile[0],"pressure")
#cart_proj = get_cartopy(projec)
#lats, lons = latlon_coords(p)
#fig = pyplot.figure(figsize=(10,7.5))
#geo_axes = pyplot.axes(projection=cart_proj)
#geo_axes.set_extent([-95, -65, 30, 45])
#states = NaturalEarthFeature(category = 'cultural',
#                       scale = '50m',
#                       facecolor = 'none',
#                       name = 'admin_1_states_provinces_shp')
#geo_axes.add_feature(states,linewidth=1.,edgecolor="black")
#geo_axes.coastlines('50m',linewidth=0.8)
#geo_axes.gridlines(color="black", linestyle="dotted")
#pyplot.contourf(to_np(lons), to_np(lats),
#                             to_np(ivt),levels,
#                     cmap=get_cmap('jet'),
#                             transform=crs.PlateCarree())
#pyplot.colorbar(ax=geo_axes, shrink=.86)
#pyplot.savefig('c3_day2_dim2_ivt.png')
ivt = np.zeros((71,99))
ii = 2
while ii < 24:
    for k in range(8,17):
        #Use QVapor for the water vapor mixing ratio
        #Use the mixing ratio to compute the specific humidity at each level
        #Compute every 50hPa from 1000hPa to 100 hPa
        w = getvar(c3ncfile[ii], "QVAPOR", timeidx=k)
        p = getvar(c3ncfile[ii],"pressure", timeidx=k)
        #w = vinterp(ncfile,h,'p',[1000,950,900,850,800,750,700,650,600,550,500,450,400,350,300,250,200,150,100])
        
        #Compute the specific humidity using q = (w/w+1)
        q = w / (1+w)
        
        #Get the vector wind V at each level by interpolating the windspeed to each level
        uv = getvar(c3ncfile[ii], 'uvmet', timeidx=k)
        #vspd = vinterp(ncfile,vsp,'p',[1000,950,900,850,800,750,700,650,600,550,500,450,400,350,300,250,200,150,100])
        #Find the IVT using the equation IVT = 1/g * sum(wspd*q*p)/number of levels
        #Make sure to remove all NaN values from the array
        #Create the pressure array based on the size of the other arrays
        
        parray = np.zeros((np.size(q,0), np.size(q,1), np.size(q,2)))
        i = 0
        j = 0
        k = 0
        pnew = 1000
        
        while i < np.size(parray,0):
            while j < np.size(parray,1):
                while k < np.size(parray,2):
                    parray[i,j,k] = pnew
                    k = k + 1
                
                k = 0
                j = j + 1
                
            i = i + 1
            j = 0
            k = 0
            pnew = pnew - 50
        
        #Method 1: Do layer from 100 hPa to 1000 hPa
        newarray = np.sqrt(((1/9.81)*(q[17,:,:]*uv[0,17,:,:]*p[17,:,:] - q[0,:,:] * uv[0,0,:,:]*p[0,:,:]))**2 + ((1/9.81)*(q[17,:,:]*uv[1,17,:,:]*p[17,:,:] - q[0,:,:]*uv[1,0,:,:]*p[0,:,:]))**2)
        ivt = ivt + newarray
    ii = ii + 4
ivt = ivt / 6
#ivt = ivt[:,:] - c3ivt[:,:]
levels = np.arange(0,300,25)
projec = getvar(c3ncfile[0],"pressure")
cart_proj = get_cartopy(projec)
lats, lons = latlon_coords(p)
fig = pyplot.figure(figsize=(10,7.5))
geo_axes = pyplot.axes(projection=cart_proj)
geo_axes.set_extent([-95, -65, 30, 45])
states = NaturalEarthFeature(category = 'cultural',
                       scale = '50m',
                       facecolor = 'none',
                       name = 'admin_1_states_provinces_shp')
geo_axes.add_feature(states,linewidth=1.,edgecolor="black")
geo_axes.coastlines('50m',linewidth=0.8)
geo_axes.gridlines(color="black", linestyle="dotted")
pyplot.contourf(to_np(lons), to_np(lats),
                             to_np(ivt),levels,
                     cmap=get_cmap('jet'),
                             transform=crs.PlateCarree())
pyplot.colorbar(ax=geo_axes, shrink=.86)
pyplot.savefig('c3_day1_dim1_reg_ivt_anom.png')
#Import the Data files
c3ncfile = [0] * 24
c6ncfile = [0] * 20
#C3dates
c3ncfile[0] = Dataset("H:/reg/wrfout_d01_2006-10-26")
c3ncfile[1] = Dataset("H:/reg/wrfout_d02_2006-10-26")
c3ncfile[2] = Dataset("H:/reg/wrfout_d01_2006-10-27")
c3ncfile[3] = Dataset("H:/reg/wrfout_d02_2006-10-27")

c3ncfile[4] = Dataset("H:/reg/wrfout_d01_2008-04-27")
c3ncfile[5] = Dataset("H:/reg/wrfout_d02_2008-04-27")
c3ncfile[6] = Dataset("H:/reg/wrfout_d01_2008-04-28")
c3ncfile[7] = Dataset("H:/reg/wrfout_d02_2008-04-28")

c3ncfile[8] = Dataset("H:/reg/wrfout_d01_2008-12-10")
c3ncfile[9] = Dataset("H:/reg/wrfout_d02_2008-12-10")
c3ncfile[10] = Dataset("H:/reg/wrfout_d01_2008-12-11")
c3ncfile[11] = Dataset("H:/reg/wrfout_d02_2008-12-11")

c3ncfile[12] = Dataset("H:/reg/wrfout_d01_2002-05-12")
c3ncfile[13] = Dataset("H:/reg/wrfout_d02_2002-05-12")
c3ncfile[14] = Dataset("H:/reg/wrfout_d01_2002-05-13")
c3ncfile[15] = Dataset("H:/reg/wrfout_d02_2002-05-13")

c3ncfile[16] = Dataset("H:/reg/wrfout_d01_2002-10-15")
c3ncfile[17] = Dataset("H:/reg/wrfout_d02_2002-10-15")
c3ncfile[18] = Dataset("H:/reg/wrfout_d01_2002-10-16")
c3ncfile[19] = Dataset("H:/reg/wrfout_d02_2002-10-16")

c3ncfile[20] = Dataset("H:/reg/wrfout_d01_2003-10-27")
c3ncfile[21] = Dataset("H:/reg/wrfout_d02_2003-10-27")
c3ncfile[22] = Dataset("H:/reg/wrfout_d01_2003-10-28")
c3ncfile[23] = Dataset("H:/reg/wrfout_d02_2003-10-28")

#C6dates
c6ncfile[0] = Dataset("H:/reg/wrfout_d01_2004-04-11")
c6ncfile[1] = Dataset("H:/reg/wrfout_d02_2004-04-11")
c6ncfile[2] = Dataset("H:/reg/wrfout_d01_2004-04-12")
c6ncfile[3] = Dataset("H:/reg/wrfout_d02_2004-04-12")

c6ncfile[4] = Dataset("H:/reg/wrfout_d01_2005-11-27")
c6ncfile[5] = Dataset("H:/reg/wrfout_d02_2005-11-27")
c6ncfile[6] = Dataset("H:/reg/wrfout_d01_2005-11-28")
c6ncfile[7] = Dataset("H:/reg/wrfout_d02_2005-11-28")

c6ncfile[8] = Dataset("H:/reg/wrfout_d01_2007-10-21")
c6ncfile[9] = Dataset("H:/reg/wrfout_d02_2007-10-21")
c6ncfile[10] = Dataset("H:/reg/wrfout_d01_2007-10-22")
c6ncfile[11] = Dataset("H:/reg/wrfout_d02_2007-10-22")

c6ncfile[12] = Dataset("H:/reg/wrfout_d01_2008-03-06")
c6ncfile[13] = Dataset("H:/reg/wrfout_d02_2008-03-06")
c6ncfile[14] = Dataset("H:/reg/wrfout_d01_2008-03-07")
c6ncfile[15] = Dataset("H:/reg/wrfout_d02_2008-03-07")

c6ncfile[16] = Dataset("H:/reg/wrfout_d01_2002-12-31")
c6ncfile[17] = Dataset("H:/reg/wrfout_d02_2002-12-31")
c6ncfile[18] = Dataset("H:/reg/wrfout_d01_2003-01-01")
c6ncfile[19] = Dataset("H:/reg/wrfout_d02_2003-01-01")

ii = 0
while ii < 24:
    for k in range(16,25):
        #Use QVapor for the water vapor mixing ratio
        #Use the mixing ratio to compute the specific humidity at each level
        #Compute every 50hPa from 1000hPa to 100 hPa
        w = getvar(c3ncfile[ii], "QVAPOR", timeidx=k)
        p = getvar(c3ncfile[ii],"pressure", timeidx=k)
        #w = vinterp(ncfile,h,'p',[1000,950,900,850,800,750,700,650,600,550,500,450,400,350,300,250,200,150,100])
        
        #Compute the specific humidity using q = (w/w+1)
        q = w / (1+w)
        
        #Get the vector wind V at each level by interpolating the windspeed to each level
        #Get the vector wind V at each level by interpolating the windspeed to each level
        uv = getvar(c3ncfile[ii], 'uvmet', timeidx=k)
        #vspd = vinterp(ncfile,vsp,'p',[1000,950,900,850,800,750,700,650,600,550,500,450,400,350,300,250,200,150,100])
        #Find the IVT using the equation IVT = 1/g * sum(wspd*q*p)/number of levels
        #Make sure to remove all NaN values from the array
        #Create the pressure array based on the size of the other arrays
        
        parray = np.zeros((np.size(q,0), np.size(q,1), np.size(q,2)))
        i = 0
        j = 0
        k = 0
        pnew = 1000
        
        while i < np.size(parray,0):
            while j < np.size(parray,1):
                while k < np.size(parray,2):
                    parray[i,j,k] = pnew
                    k = k + 1
                
                k = 0
                j = j + 1
                
            i = i + 1
            j = 0
            k = 0
            pnew = pnew - 50
        
        #Method 1: Do layer from 100 hPa to 1000 hPa
        newarray = np.sqrt(((1/9.81)*(q[17,:,:]*uv[0,17,:,:]*300 - q[0,:,:] * uv[0,0,:,:]*1000))**2 + ((1/9.81)*(q[17,:,:]*uv[1,17,:,:]*300 - q[0,:,:]*uv[1,0,:,:]*1000))**2)
        ivt = ivt + newarray
    ii = ii + 4
ivt = ivt / 6
#ivt = ivt[:,:] - c3ivt[:,:]
levels = np.arange(0,300,25)
projec = getvar(c3ncfile[0],"pressure")
cart_proj = get_cartopy(projec)
lats, lons = latlon_coords(p)
fig = pyplot.figure(figsize=(10,7.5))
geo_axes = pyplot.axes(projection=cart_proj)
geo_axes.set_extent([-95, -65, 30, 45])
states = NaturalEarthFeature(category = 'cultural',
                       scale = '50m',
                       facecolor = 'none',
                       name = 'admin_1_states_provinces_shp')
geo_axes.add_feature(states,linewidth=1.,edgecolor="black")
geo_axes.coastlines('50m',linewidth=0.8)
geo_axes.gridlines(color="black", linestyle="dotted")
pyplot.contourf(to_np(lons), to_np(lats),
                             to_np(ivt),levels,
                     cmap=get_cmap('jet'),
                             transform=crs.PlateCarree())
pyplot.colorbar(ax=geo_axes, shrink=.86)
pyplot.savefig('c3_day2_dim1_reg_ivt_anom.png')
#ii = 1
#while ii < 16:
#    for k in range(16,25):
#        #Use QVapor for the water vapor mixing ratio
#        #Use the mixing ratio to compute the specific humidity at each level
#        #Compute every 50hPa from 1000hPa to 100 hPa
#        w = getvar(c3ncfile[ii], "QVAPOR", timeidx=k)
#        p = getvar(c3ncfile[ii],"pressure", timeidx=k)
#        #w = vinterp(ncfile,h,'p',[1000,950,900,850,800,750,700,650,600,550,500,450,400,350,300,250,200,150,100])
#        
#        #Compute the specific humidity using q = (w/w+1)
#        q = w / (1+w)
#        
#        #Get the vector wind V at each level by interpolating the windspeed to each level
#        usp = getvar(c3ncfile[ii], 'ua', timeidx=k)
#        #uspd = vinterp(ncfile,usp,'p',[1000,950,900,850,800,750,700,650,600,550,500,450,400,350,300,250,200,150,100])
#        
#        vsp = getvar(c3ncfile[ii], 'va', timeidx=k)
#        
#        wsp = getvar(c3ncfile[ii], 'wspd', timeidx=k)
#        #vspd = vinterp(ncfile,vsp,'p',[1000,950,900,850,800,750,700,650,600,550,500,450,400,350,300,250,200,150,100])
#        #Find the IVT using the equation IVT = 1/g * sum(wspd*q*p)/number of levels
#        #Make sure to remove all NaN values from the array
#        #Create the pressure array based on the size of the other arrays
#        
#        parray = np.zeros((np.size(q,0), np.size(q,1), np.size(q,2)))
#        i = 0
#        j = 0
#        k = 0
#        pnew = 1000
#        
#        while i < np.size(parray,0):
#            while j < np.size(parray,1):
#                while k < np.size(parray,2):
#                    parray[i,j,k] = pnew
#                    k = k + 1
#                
#                k = 0
#                j = j + 1
#                
#            i = i + 1-200,200
#            j = 0
#            k = 0
#            pnew = pnew - 50
#        
#        #Method 1: Do layer from 100 hPa to 1000 hPa
#        newarray = (wsp[0,:,:]*q[0,:,:])*p[0,:,:] - (wsp[18,:,:]*q[18,:,:]*p[18,:,:])
#        ivt = ivt + newarray / 9.81
#    ii = ii + 4
#ivt = ivt / 4
#levels = np.arange(0,250,25)
#projec = getvar(c3ncfile[0],"pressure")
#cart_proj = get_cartopy(projec)
#lats, lons = latlon_coords(p)
#fig = pyplot.figure(figsize=(10,7.5))
#geo_axes = pyplot.axes(projection=cart_proj)
#geo_axes.set_extent([-95, -65, 30, 45])
#states = NaturalEarthFeature(category = 'cultural',
#                       scale = '50m',
#                       facecolor = 'none',
#                       name = 'admin_1_states_provinces_shp')
#geo_axes.add_feature(states,linewidth=1.,edgecolor="black")
#geo_axes.coastlines('50m',linewidth=0.8)
#geo_axes.gridlines(color="black", linestyle="dotted")
#pyplot.contourf(to_np(lons), to_np(lats),
#                             to_np(ivt),levels,
#                     cmap=get_cmap('jet'),
#                             transform=crs.PlateCarree())
#pyplot.colorbar(ax=geo_axes, shrink=.86)
#pyplot.savefig('c3_day2_dim2_ivt.png')
ivt = np.zeros((71,99))
ii = 2
while ii < 24:
    for k in range(8,17):
        #Use QVapor for the water vapor mixing ratio
        #Use the mixing ratio to compute the specific humidity at each level
        #Compute every 50hPa from 1000hPa to 100 hPa
        w = getvar(c3ncfile[ii], "QVAPOR", timeidx=k)
        p = getvar(c3ncfile[ii],"pressure", timeidx=k)
        #w = vinterp(ncfile,h,'p',[1000,950,900,850,800,750,700,650,600,550,500,450,400,350,300,250,200,150,100])
        
        #Compute the specific humidity using q = (w/w+1)
        q = w / (1+w)
        
        #Get the vector wind V at each level by interpolating the windspeed to each level
        uv = getvar(c3ncfile[ii], 'uvmet', timeidx=k)
        u = getvar(c3ncfile[ii], 'ua', timeidx=k)
        v = getvar(c3ncfile[ii], 'va', timeidx=k)
        #vspd = vinterp(ncfile,vsp,'p',[1000,950,900,850,800,750,700,650,600,550,500,450,400,350,300,250,200,150,100])
        #Find the IVT using the equation IVT = 1/g * sum(wspd*q*p)/number of levels
        #Make sure to remove all NaN values from the array
        #Create the pressure array based on the size of the other arrays
        
        parray = np.zeros((np.size(q,0), np.size(q,1), np.size(q,2)))
        i = 0
        j = 0
        k = 0
        pnew = 1000
        
        while i < np.size(parray,0):
            while j < np.size(parray,1):
                while k < np.size(parray,2):
                    parray[i,j,k] = pnew
                    k = k + 1
                
                k = 0
                j = j + 1
                
            i = i + 1
            j = 0
            k = 0
            pnew = pnew - 50
        
        #Method 1: Do layer from 300 hPa to 1000 hPa
        newarray = np.sqrt(((1/9.81)*(q[17,:,:]*u[17,:,:]*300 - q[0,:,:] * u[0,:,:]*1000))**2 + ((1/9.81)*(q[17,:,:]*v[17,:,:]*300 - q[0,:,:]*v[0,:,:]*1000))**2)
        ivt = ivt + newarray
    ii = ii + 4
ivt = ivt / 6
#ivt = ivt[:,:] - c3ivt[:,:]
levels = np.arange(0,300,25)
projec = getvar(c3ncfile[0],"pressure")
cart_proj = get_cartopy(projec)
lats, lons = latlon_coords(p)
fig = pyplot.figure(figsize=(10,7.5))
geo_axes = pyplot.axes(projection=cart_proj)
geo_axes.set_extent([-95, -65, 30, 45])
states = NaturalEarthFeature(category = 'cultural',
                       scale = '50m',
                       facecolor = 'none',
                       name = 'admin_1_states_provinces_shp')
geo_axes.add_feature(states,linewidth=1.,edgecolor="black")
geo_axes.coastlines('50m',linewidth=0.8)
geo_axes.gridlines(color="black", linestyle="dotted")
pyplot.contourf(to_np(lons), to_np(lats),
                             to_np(ivt),levels,
                     cmap=get_cmap('jet'),
                             transform=crs.PlateCarree())
pyplot.colorbar(ax=geo_axes, shrink=.86)
pyplot.savefig('c3_day1_dim1_reg_ivt_anom.png')