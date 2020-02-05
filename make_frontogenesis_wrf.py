#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 10:01:33 2019

@author: mariofire
"""

import numpy as np
from matplotlib import pyplot
from matplotlib.cm import get_cmap
from cartopy import crs
from cartopy.feature import NaturalEarthFeature
from netCDF4 import Dataset
from wrf import getvar, to_np, get_cartopy, latlon_coords,  interplevel
import metpy.calc as mpcalc
from metpy.units import units
import cartopy.feature as cfeature
import os
import xarray as xr

#pyplot.switch_backend('agg')

#C3dates
#Import the Data files
c3ncfile = [0] * 40
c6ncfile = [0] * 40
c3nncfile = [0] * 40
c6nncfile = [0] * 40
#C3dates
c3ncfile[0] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2006-10-26")
c3ncfile[1] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2006-10-26")
c3ncfile[2] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2006-10-27")
c3ncfile[3] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2006-10-27")

c3ncfile[4] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2008-04-27")
c3ncfile[5] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2008-04-27")
c3ncfile[6] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2008-04-28")
c3ncfile[7] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2008-04-28")

c3ncfile[8] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2008-12-10")
c3ncfile[9] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2008-12-10")
c3ncfile[10] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2008-12-11")
c3ncfile[11] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2008-12-11")

c3ncfile[12] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2002-05-12")
c3ncfile[13] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2002-05-12")
c3ncfile[14] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2002-05-13")
c3ncfile[15] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2002-05-13")

c3ncfile[16] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2002-10-15")
c3ncfile[17] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2002-10-15")
c3ncfile[18] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2002-10-16")
c3ncfile[19] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2002-10-16")

c3ncfile[20] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2003-10-27")
c3ncfile[21] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2003-10-27")
c3ncfile[22] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2003-10-28")
c3ncfile[23] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2003-10-28")

c3ncfile[24] = ""
c3ncfile[25] = ""
c3ncfile[26] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2000-03-11")
c3ncfile[27] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2000-03-11")

c3ncfile[28] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2000-04-02")
c3ncfile[29] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2000-04-02")
c3ncfile[30] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2000-04-03")
c3ncfile[31] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2000-04-03")

c3ncfile[32] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2001-08-18")
c3ncfile[33] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2001-08-18")
c3ncfile[34] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2001-08-19")
c3ncfile[35] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2001-08-19")

c3ncfile[36] = ""
c3ncfile[37] = ""
c3ncfile[38] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2006-06-26")
c3ncfile[39] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2006-06-26")

#C3dates nodh
c3nncfile[0] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d01_2006-10-26")
c3nncfile[1] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d02_2006-10-26")
c3nncfile[2] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d01_2006-10-27")
c3nncfile[3] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d02_2006-10-27")

c3nncfile[4] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d01_2008-04-27")
c3nncfile[5] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d02_2008-04-27")
c3nncfile[6] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d01_2008-04-28")
c3nncfile[7] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d02_2008-04-28")

c3nncfile[8] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d01_2008-12-10")
c3nncfile[9] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d02_2008-12-10")
c3nncfile[10] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d01_2008-12-11")
c3nncfile[11] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d02_2008-12-11")

c3nncfile[12] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d01_2002-05-12")
c3nncfile[13] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d02_2002-05-12")
c3nncfile[14] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d01_2002-05-13")
c3nncfile[15] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d02_2002-05-13")

c3nncfile[16] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d01_2002-10-15")
c3nncfile[17] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d02_2002-10-15")
c3nncfile[18] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d01_2002-10-16")
c3nncfile[19] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d02_2002-10-16")

c3nncfile[20] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d01_2003-10-27")
c3nncfile[21] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d02_2003-10-27")
c3nncfile[22] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d01_2003-10-28")
c3nncfile[23] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d02_2003-10-28")

c3nncfile[24] = ""
c3nncfile[25] = ""
c3nncfile[26] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2000-03-11")
c3nncfile[27] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2000-03-11")

c3nncfile[28] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2000-04-02")
c3nncfile[29] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2000-04-02")
c3nncfile[30] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2000-04-03")
c3nncfile[31] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2000-04-03")

c3nncfile[32] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2001-08-18")
c3nncfile[33] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2001-08-18")
c3nncfile[34] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2001-08-19")
c3nncfile[35] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2001-08-19")

c3nncfile[36] = ""
c3nncfile[37] = ""
c3nncfile[38] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2006-06-26")
c3nncfile[39] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2006-06-26")

#C6dates
c6ncfile[0] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2004-04-11")
c6ncfile[1] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2004-04-11")
c6ncfile[2] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2004-04-12")
c6ncfile[3] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2004-04-12")

c6ncfile[4] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2005-11-27")
c6ncfile[5] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2005-11-27")
c6ncfile[6] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2005-11-28")
c6ncfile[7] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2005-11-28")

c6ncfile[8] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2007-10-21")
c6ncfile[9] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2007-10-21")
c6ncfile[10] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2007-10-22")
c6ncfile[11] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2007-10-22")

c6ncfile[12] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2008-03-06")
c6ncfile[13] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2008-03-06")
c6ncfile[14] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2008-03-07")
c6ncfile[15] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2008-03-07")

c6ncfile[16] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2002-12-31")
c6ncfile[17] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2002-12-31")
c6ncfile[18] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2003-01-01")
c6ncfile[19] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2003-01-01")

c6ncfile[20] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2000-12-15")
c6ncfile[21] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2000-12-15")
c6ncfile[22] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2000-12-16")
c6ncfile[23] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2000-12-16")

c6ncfile[24] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2001-03-28")
c6ncfile[25] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2001-03-28")
c6ncfile[26] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2001-03-29")
c6ncfile[27] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2001-03-29")

c6ncfile[28] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2000-06-13")
c6ncfile[29] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2000-06-13")
c6ncfile[30] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2000-06-14")
c6ncfile[31] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2000-06-14")

c6ncfile[32] = ""
c6ncfile[33] = ""
c6ncfile[34] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2002-11-16")
c6ncfile[35] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2002-11-16")

c6ncfile[36] = ""
c6ncfile[37] = ""
c6ncfile[38] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2006-01-17")
c6ncfile[39] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2006-01-17")

#C6datesnodh
c6nncfile[0] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d01_2004-04-11")
c6nncfile[1] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d02_2004-04-11")
c6nncfile[2] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d01_2004-04-12")
c6nncfile[3] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d02_2004-04-12")

c6nncfile[4] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d01_2005-11-27")
c6nncfile[5] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d02_2005-11-27")
c6nncfile[6] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d01_2005-11-28")
c6nncfile[7] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d02_2005-11-28")

c6nncfile[8] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d01_2007-10-21")
c6nncfile[9] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d02_2007-10-21")
c6nncfile[10] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d01_2007-10-22")
c6nncfile[11] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d02_2007-10-22")

c6nncfile[12] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d01_2008-03-06")
c6nncfile[13] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d02_2008-03-06")
c6nncfile[14] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d01_2008-03-07")
c6nncfile[15] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d02_2008-03-07")

c6nncfile[16] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d01_2002-12-31")
c6nncfile[17] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d02_2002-12-31")
c6nncfile[18] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d01_2003-01-01")
c6nncfile[19] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d02_2003-01-01")

c6nncfile[20] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d01_2000-12-15")
c6nncfile[21] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d02_2000-12-15")
c6nncfile[22] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d01_2000-12-16")
c6nncfile[23] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d02_2000-12-16")

c6nncfile[24] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d01_2001-03-28")
c6nncfile[25] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d02_2001-03-28")
c6nncfile[26] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d01_2001-03-29")
c6nncfile[27] = Dataset("H:/newrf/wrf_files/nodh/wrfout_d02_2001-03-29")

c6nncfile[28] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2000-06-13")
c6nncfile[29] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2000-06-13")
c6nncfile[30] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2000-06-14")
c6nncfile[31] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2000-06-14")

c6nncfile[32] = ""
c6nncfile[33] = ""
c6nncfile[34] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2002-11-16")
c6nncfile[35] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2002-11-16")

c6nncfile[36] = ""
c6nncfile[37] = ""
c6nncfile[38] = Dataset("H:/newrf/wrf_files/reg/wrfout_d01_2006-01-17")
c6nncfile[39] = Dataset("H:/newrf/wrf_files/reg/wrfout_d02_2006-01-17")


#Set the background dimensions of the plot for the outer domain
def plot_background(ax):
    ax.set_extent([-95, -65, 30, 45])
    ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.5)
    ax.add_feature(cfeature.STATES, linewidth=0.5)
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    return ax

i = 2
frontot = np.zeros((71,99))
htot = np.zeros((71,99))
ttot = np.zeros((71,99))
frontotc3_reg = np.zeros((71,99))
frontotc3_nodh = np.zeros((71,99))
frontotc6_reg = np.zeros((71,99))
frontotc6_nodh = np.zeros((71,99))
hc3_reg = np.zeros((71,99))
hc3_nodh = np.zeros((71,99))
hc6_reg = np.zeros((71,99))
hc6_nodh = np.zeros((71,99))
tc3_reg = np.zeros((71,99))
tc3_nodh = np.zeros((71,99))
tc6_reg = np.zeros((71,99))
tc6_nodh = np.zeros((71,99))
count = 0
while(i < 40):
    modulo = i % 4
    if(modulo == 0):
        day = 2
        dim = 1
    if(modulo == 2):
        day = 1 
        dim = 1
    if(day == 1 and dim == 1):
        time = 8
        et = 17
    elif( day == 1 and dim == 2):
        time = 24
        et = 49
    elif( day == 2 and dim == 1):
        time = 15
        et = 25
    else:
        time = 48
        et = 73
    st=time
    nc1 = c6ncfile[i]
    date = nc1.START_DATE
    yr = date[0:4]
    dy = date[5:7]
    mh = date[8:10]
    date = yr + dy+mh
    print(date)
    while( time < et):
        level = 850 * units.hPa
        h = getvar(nc1, "z", timeidx=time, units = 'dm')
        u = getvar(nc1, "ua", timeidx=time)
        v = getvar(nc1, "va", timeidx=time)
        temp = getvar(nc1, "tc", timeidx=time)
        p = getvar(nc1,"pressure", timeidx=time)
        hght_850 = interplevel(h,p,850)
        tmpk_850 = interplevel(temp,p,850)
        tmpk_850.metpy.convert_units('kelvin')
        uwnd_850 = interplevel(u,p,850)
        vwnd_850 = interplevel(v,p,850)
        uwnd_850 = uwnd_850.values *units('m/s')
        vwnd_850 = vwnd_850.values * units('m/s')
        lats, lons = latlon_coords(hght_850)
        # Calculate potential temperature for frontogenesis calculation
        thta_850 = mpcalc.potential_temperature(level, tmpk_850)

        dx, dy = mpcalc.lat_lon_grid_deltas(lons, lats)

        fronto_850 = mpcalc.frontogenesis(thta_850, uwnd_850, vwnd_850, dx, dy, dim_order='yx')

        frontot = frontot + fronto_850
        htot = htot + hght_850
        tmpk_850.metpy.convert_units('degC')
        ttot = ttot + tmpk_850
        time = time + 1
    frontotc6_reg[:,:] = frontotc6_reg + frontot[:,:] *1000*100*3600*3
    hc6_reg[:,:] = hc6_reg[:,:] + htot[:,:] 
    tc6_reg[:,:] = tc6_reg[:,:] + ttot[:,:] 

    frontot = np.zeros((71,99))
    htot = np.zeros((71,99))
    ttot = np.zeros((71,99))
    time = st
    nc1 = c6nncfile[i]
    date = nc1.START_DATE
    yr = date[0:4]
    dy = date[5:7]
    mh = date[8:10]
    date = yr + dy+mh
    print(date)
    while( time < et):
        level = 850 * units.hPa
        h = getvar(nc1, "z", timeidx=time, units = 'dm')
        u = getvar(nc1, "ua", timeidx=time)
        v = getvar(nc1, "va", timeidx=time)
        temp = getvar(nc1, "tc", timeidx=time)
        p = getvar(nc1,"pressure", timeidx=time)
        hght_850 = interplevel(h,p,850)
        tmpk_850 = interplevel(temp,p,850)
        uwnd_850 = interplevel(u,p,850)
        vwnd_850 = interplevel(v,p,850)
        uwnd_850 = uwnd_850.values *units('m/s')
        vwnd_850 = vwnd_850.values * units('m/s')
        lats, lons = latlon_coords(hght_850)
        tmpk_850.metpy.convert_units('kelvin')
        # Calculate potential temperature for frontogenesis calculation
        thta_850 = mpcalc.potential_temperature(level, tmpk_850)

        dx, dy = mpcalc.lat_lon_grid_deltas(lons, lats)

        fronto_850 = mpcalc.frontogenesis(thta_850, uwnd_850, vwnd_850, dx, dy, dim_order='yx')

        frontot = frontot + fronto_850
        htot = htot + hght_850
        ttot = ttot + tmpk_850 - 273.15
        time = time + 1
    frontotc6_nodh[:,:] = frontotc6_nodh[:,:] + frontot[:,:] *1000*100*3600*3
    hc6_nodh[:,:] = hc6_nodh[:,:]+htot[:,:] 
    tc6_nodh[:,:] = tc6_nodh[:,:]+ttot[:,:] 
    frontot = np.zeros((71,99))
    htot = np.zeros((71,99))
    ttot = np.zeros((71,99))
    count = count + 1
    i = i + 4

frontot = np.zeros((71,99))
htot = np.zeros((71,99))
ttot = np.zeros((71,99))  
count = 0
i = 2
while(i < 40):
    modulo = i % 4
    if(modulo == 0):
        day = 2
        dim = 1
    if(modulo == 2):
        day = 1 
        dim = 1
    if(day == 1 and dim == 1):
        time = 8
        et = 17
    elif( day == 1 and dim == 2):
        time = 24
        et = 49
    elif( day == 2 and dim == 1):
        time = 15
        et = 25
    else:
        time = 48
        et = 73
    st=time
    nc1 = c3ncfile[i]
    date = nc1.START_DATE
    yr = date[0:4]
    dy = date[5:7]
    mh = date[8:10]
    date = yr + dy+mh
    print(date)
    while( time < et):
        level = 850 * units.hPa
        h = getvar(nc1, "z", timeidx=time, units = 'dm')
        u = getvar(nc1, "ua", timeidx=time)
        v = getvar(nc1, "va", timeidx=time)
        temp = getvar(nc1, "tc", timeidx=time)
        p = getvar(nc1,"pressure", timeidx=time)
        hght_850 = interplevel(h,p,850)
        tmpk_850 = interplevel(temp,p,850)
        uwnd_850 = interplevel(u,p,850)
        vwnd_850 = interplevel(v,p,850)
        uwnd_850 = uwnd_850.values *units('m/s')
        vwnd_850 = vwnd_850.values * units('m/s')
        lats, lons = latlon_coords(hght_850)
        tmpk_850.metpy.convert_units('kelvin')
        # Calculate potential temperature for frontogenesis calculation
        thta_850 = mpcalc.potential_temperature(level, tmpk_850)

        dx, dy = mpcalc.lat_lon_grid_deltas(lons, lats)

        fronto_850 = mpcalc.frontogenesis(thta_850, uwnd_850, vwnd_850, dx, dy, dim_order='yx')

        frontot = frontot + fronto_850
        htot = htot + hght_850
        ttot = ttot + tmpk_850 - 273.15
        time = time + 1
    frontotc3_reg = frontotc3_reg+frontot *1000*100*3600*3
    hc3_reg = hc3_reg+htot 
    tc3_reg = tc3_reg+ttot 

    frontot = np.zeros((71,99))
    htot = np.zeros((71,99))
    ttot = np.zeros((71,99))
    time = st
    nc1 = c3nncfile[i]
    date = nc1.START_DATE
    yr = date[0:4]
    dy = date[5:7]
    mh = date[8:10]
    date = yr + dy+mh
    print(date)
    while( time < et):
        level = 850 * units.hPa
        h = getvar(nc1, "z", timeidx=time, units = 'dm')
        u = getvar(nc1, "ua", timeidx=time)
        v = getvar(nc1, "va", timeidx=time)
        temp = getvar(nc1, "tc", timeidx=time)
        p = getvar(nc1,"pressure", timeidx=time)
        hght_850 = interplevel(h,p,850)
        tmpk_850 = interplevel(temp,p,850)
        uwnd_850 = interplevel(u,p,850)
        vwnd_850 = interplevel(v,p,850)
        uwnd_850 = uwnd_850.values *units('m/s')
        vwnd_850 = vwnd_850.values * units('m/s')
        lats, lons = latlon_coords(hght_850)
        # Calculate potential temperature for frontogenesis calculation
        tmpk_850.metpy.convert_units('kelvin')
        thta_850 = mpcalc.potential_temperature(level, tmpk_850)

        dx, dy = mpcalc.lat_lon_grid_deltas(lons, lats)

        fronto_850 = mpcalc.frontogenesis(thta_850, uwnd_850, vwnd_850, dx, dy, dim_order='yx')

        frontot = frontot + fronto_850
        htot = htot + hght_850
        ttot = ttot + tmpk_850 - 273.15
        time = time + 1
    frontotc3_nodh[:,:] = frontotc3_nodh[:,:]+frontot[:,:] *1000*100*3600*3
    hc3_nodh[:,:] = hc3_nodh[:,:]+htot[:,:] 
    tc3_nodh[:,:] = tc3_nodh[:,:]+ttot[:,:] 
    frontot = np.zeros((71,99))
    htot = np.zeros((71,99))
    ttot = np.zeros((71,99))
    count = count + 1
    i = i + 4

#frontotc3_reg2 = frontotc3_reg.mean(axis=0)
#frontotc3_nodh2 = frontotc3_nodh.mean(axis=0)
#frontotc6_reg2 = frontotc6_reg.mean(axis=0)
#frontotc6_nodh2 = frontotc6_nodh.mean(axis=0)
frontotc3_reg = frontotc3_reg / 10
frontotc3_nodh = frontotc3_nodh / 10
frontotc6_reg = frontotc6_reg / 10
frontotc6_nodh = frontotc6_nodh / 10
hc3_reg = hc3_reg / 10
hc3_nodh = hc3_nodh / 10
hc6_reg = hc6_reg / 10
hc6_nodh = hc6_nodh / 10
tc3_reg = tc3_reg / 10
tc3_nodh = tc3_nodh / 10
tc6_reg = tc6_reg / 10
tc6_nodh = tc6_nodh / 10
#hc3_reg2 = hc3_reg.mean(axis=0)
#hc3_nodh2 = hc3_nodh.mean(axis=0)
#hc6_reg2 = hc6_reg.mean(axis=0)
#hc6_nodh2 = hc6_nodh.mean(axis=0)
#tc3_reg2 = tc3_reg.mean(axis=0)
#tc3_nodh2 = tc3_nodh.mean(axis=0)
#tc6_reg2 = tc6_reg.mean(axis=0)
#tc6_nodh2 = tc6_nodh.mean(axis=0)

filenames = os.listdir("H:/ERA5/era5_front")
filenames.sort()
i = 0

front_overall = np.zeros((20,60,120))
h850_overall = np.zeros((20,60,120))
t850_overall = np.zeros((20,60,120))
while i < len(filenames):
    ds = xr.open_dataset('H:/ERA5/era5_front/' + filenames[i])
    
    #Calculate the potential temperature using metpy 
    thta_850 = mpcalc.potential_temperature(850 * units.hPa, ds.t)
    thta_850 = xr.DataArray(thta_850).mean('dim_0')
    
    #Convert lat/lon data to a grid spacing
    dx, dy = mpcalc.lat_lon_grid_deltas(ds.coords['longitude'], ds.coords['latitude'])
    
    #Get mean for the whole day for the wind components
    u = ds.u.mean('time')
    v = ds.v.mean('time')
    
    #Add back missing units to the DataArrays
    u.attrs['units'] = 'm/s'
    v.attrs['units'] = 'm/s'
    thta_850.attrs['units'] = 'kelvin'
    
    #Calculate frontogenesis using metpy
    fronto_850 = mpcalc.frontogenesis(thta_850,u,v, dx, dy, dim_order='yx')
    
    fronto_850 = fronto_850 * 1000 * 100 * 3600 * 3
    front_new = fronto_850[180:240,1060:1180]
    front_overall[i,:,:] = front_overall[i,:,:] + front_new[:,:]
    h850_overall[i,:,:] = h850_overall[i,:,:] + (ds.z.mean('time')[180:240,1060:1180]/9.81)
    t850_overall[i,:,:] = t850_overall[i,:,:] + (ds.t.mean('time')[180:240,1060:1180] - 273.15)
    del ds
    i = i + 1

dates = ['20061028','20000312',
         '20000404','20000615',
         '20001217','20010330',
         '20010820','20020514',
         '20021017','20021117',
         '20030102','20031029',
         '20040413','20051129',
         '20060118','20060627',
         '20071023','20080308',
         '20080429','20081212']

c3 = [0,1,2,6,7,8,11,15,18,19]
c6 = [3,4,5,9,10,12,13,14,16,17]

c3_dates = [15,18,19,6,7,10]
c6_dates = [3,4,11,12,16,17,10]
ds = xr.open_mfdataset('H:/ERA5/era5_front/' + filenames[i-1])
#c3_dates = [0,1,5,6,7,10,14,15,18,19]
#c6_dates = [2,3,4,8,9,11,12,13,16,17] 
lons2 = ds.coords['longitude']
lats2 = ds.coords['latitude']

lons2 = lons2 - 360

cart_proj = crs.PlateCarree()
fig, axarr = pyplot.subplots(nrows=2, ncols=3, figsize=(12, 8), constrained_layout=True,
                  subplot_kw={'projection': cart_proj})

axlist = axarr.flatten()
for ax in axlist:
    #D1
        plot_background(ax)
#ivt = ivt[:,:] - c3ivt[:,:]
levels = np.arange(-5,5,0.5)
levels2 = np.arange(100,200,8)
levels3 = np.arange(-20,20,4)


ax0 = axlist[0].contourf(to_np(lons2[1060:1180]), np.flipud(to_np(lats2[180:240])),
                             np.flipud(to_np(np.mean(front_overall[c3_dates,:,:],axis=0))),levels,
                     cmap=get_cmap('seismic'),extend='both',
                             transform=crs.PlateCarree())
csf = axlist[0].contour(to_np(lons2[1060:1180]), np.flipud(to_np(lats2[180:240])),
                             np.flipud(to_np(np.mean(h850_overall[c3_dates,:,:],axis=0)))/10,levels2,colors='black',
                             transform=crs.PlateCarree())
cs = axlist[0].contour(to_np(lons2[1060:1180]), np.flipud(to_np(lats2[180:240])),
                             np.flipud(to_np(np.mean(t850_overall[c3_dates,:,:],axis=0))),levels3,colors='blue',transform=crs.PlateCarree())
axlist[0].clabel(csf, fmt='%d')
axlist[0].clabel(cs, fmt='%d')
fig.colorbar(ax0,ax=axlist[0], shrink=.86)
titlename = "C3 Reanalysis"
axlist[0].set_title(titlename)

ax1 = axlist[1].contourf(to_np(lons), to_np(lats),
                             to_np(frontotc3_reg/10),levels,
                     cmap=get_cmap('seismic'),extend='both',
                             transform=crs.PlateCarree())
csf = axlist[1].contour(to_np(lons),to_np(lats),to_np(hc3_reg/10),levels2,colors='black',transform=crs.PlateCarree())
cs = axlist[1].contour(to_np(lons),to_np(lats),to_np(tc3_reg/10),levels3,colors='blue',transform=crs.PlateCarree())
axlist[1].clabel(csf, fmt='%d')
axlist[1].clabel(cs, fmt='%d')
fig.colorbar(ax1,ax=axlist[1], shrink=.86)
titlename = "C3 Control"
axlist[1].set_title(titlename)

ax2 = axlist[2].contourf(to_np(lons), to_np(lats),
                             to_np(frontotc3_nodh/10),levels,
                     cmap=get_cmap('seismic'),extend='both',
                             transform=crs.PlateCarree())
csf = axlist[2].contour(to_np(lons),to_np(lats),to_np(hc3_nodh/10),levels2,colors='black',transform=crs.PlateCarree())
cs = axlist[2].contour(to_np(lons),to_np(lats),to_np(tc3_nodh/10),levels3,colors='blue',transform=crs.PlateCarree())
axlist[2].clabel(csf, fmt='%d')
axlist[2].clabel(cs, fmt='%d')
fig.colorbar(ax2,ax=axlist[2], shrink=.86)
titlename = "C3 No Diabatic Heating"
axlist[2].set_title(titlename)

ax33 = axlist[3].contourf(to_np(lons2[1060:1180]), np.flipud(to_np(lats2[180:240])),
                             np.flipud(to_np(np.mean(front_overall[c6_dates,:,:],axis=0))),levels,
                     cmap=get_cmap('seismic'),extend='both',
                             transform=crs.PlateCarree())
csf = axlist[3].contour(to_np(lons2[1060:1180]), np.flipud(to_np(lats2[180:240])),
                             np.flipud(to_np(np.mean(h850_overall[c6_dates,:,:],axis=0)))/10,levels2,colors='black',
                             transform=crs.PlateCarree())
cs = axlist[3].contour(to_np(lons2[1060:1180]), np.flipud(to_np(lats2[180:240])),
                             np.flipud(to_np(np.mean(t850_overall[c6_dates,:,:],axis=0))),levels3,colors='blue',transform=crs.PlateCarree())
axlist[3].clabel(csf, fmt='%d')
axlist[3].clabel(cs, fmt='%d')
fig.colorbar(ax33,ax=axlist[3], shrink=.86)
titlename = "C6 Reanalysis"
axlist[3].set_title(titlename)

ax3 = axlist[4].contourf(to_np(lons), to_np(lats),
                             to_np(frontotc6_reg/10),levels,
                     cmap=get_cmap('seismic'),extend='both',
                             transform=crs.PlateCarree())
csf = axlist[4].contour(to_np(lons),to_np(lats),to_np(hc6_reg/10),levels2,colors='black',transform=crs.PlateCarree())
cs = axlist[4].contour(to_np(lons),to_np(lats),to_np(tc6_reg/10),levels3,colors='blue',transform=crs.PlateCarree())
axlist[4].clabel(csf, fmt='%d')
axlist[4].clabel(cs, fmt='%d')
fig.colorbar(ax3,ax=axlist[4], shrink=.86)
titlename = "C6 Control"
axlist[4].set_title(titlename)

ax4 = axlist[5].contourf(to_np(lons), to_np(lats),
                             to_np(frontotc6_nodh/10),levels,
                     cmap=get_cmap('seismic'),extend='both',
                             transform=crs.PlateCarree())
csf = axlist[5].contour(to_np(lons),to_np(lats),to_np(hc6_nodh/10),levels2,colors='black',transform=crs.PlateCarree())
cs = axlist[5].contour(to_np(lons),to_np(lats),to_np(tc6_nodh/10),levels3,colors='blue',transform=crs.PlateCarree())
axlist[5].clabel(csf, fmt='%d')
axlist[5].clabel(cs, fmt='%d')
fig.colorbar(ax4,ax=axlist[5], shrink=.86)
titlename = "C6 No Diabatic Heating"
axlist[5].set_title(titlename)


fig.suptitle('Frontogenesis (shaded), H850 (black), and T850 (Blue)')
pyplot.savefig('D1_1_frontogenesis.png',bbox_inches='tight')