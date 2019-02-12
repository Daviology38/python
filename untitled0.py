
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 09:35:10 2019

@author: CoeFamily
"""

#!/usr/bin/env python

# Working script to generate maps from wrfout netCDF files
# using matplot lib with basemap
# Basemap coding from David John Gagne II
# Written by Luke Madaus for use with operational WRF domains

import sys,getopt
from netCDF4 import Dataset
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import os
from datetime import datetime, timedelta
import coltbls as coltbls
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid import make_axes_locatable
import matplotlib.axes as maxes
import glob
import numpy as np
from matplotlib.cm import get_cmap
from matplotlib.colors import from_levels_and_colors
from cartopy import crs
from cartopy.feature import NaturalEarthFeature
from wrf import getvar, to_np, get_cartopy, latlon_coords, cartopy_xlim, cartopy_ylim, interplevel, smooth2d


#import calc_wrf_severe as severe

# Set the default domain to be d02
dom = 'd01'
var = 'snow'
export_flag = 1
filename = '../wrfout_' + dom
runtotal = 0
runtotal1 = 0
restart_time = 0

# Set up a command-line argument structure to allow
# for command-line changes of variables.
# f --> the name of the domain we want to use
(opts,args)=getopt.getopt(sys.argv[1:],'f:v:r:e')
for o,a in opts:
	if o=="-f":
		filename = a
	if o=="-v":
		var = str(a)
	if o=="-e":
		export_flag = 1	
	if o=="-r":
		restart_time = int(a)
	

# Skip is the length between outputs
skip =0.5 

# Directory to move images to (if requested)
outdir = 'C:/Users/Coefamily/'

nc = Dataset("WRFV3/run/wrfout_d03_2019-02-12_00:00:00","r")


# Grab these variables for now
temps =  nc.variables['T2']
u_wind_ms = nc.variables['U10']
v_wind_ms = nc.variables['V10']
psfc = nc.variables['PSFC']
T = nc.variables['T']
times = nc.variables['Times']


# Thin factor is used for thinning out wind barbs
thin = 5 




# BEGIN ACTUAL PROCESSING HERE
# x_dim and y_dim are the x and y dimensions of the model
# domain in gridpoints
x_dim = len(nc.dimensions['west_east'])
y_dim = len(nc.dimensions['south_north'])

# Get the grid spacing
dx = float(nc.DX)
dy = float(nc.DY)

width_meters = dx * (x_dim - 1)
height_meters = dy * (y_dim - 1)

cen_lat = float(nc.CEN_LAT)
cen_lon = float(nc.CEN_LON)
truelat1 = float(nc.TRUELAT1)
truelat2 = float(nc.TRUELAT2)
standlon = float(nc.STAND_LON)
truelat1 = float(60)
truelat2 = float(10)





# Draw the base map behind it with the lats and
# lons calculated earlier
m = Basemap(resolution='i',projection='lcc',\
    width=width_meters,height=height_meters,\
    lat_0=cen_lat,lon_0=cen_lon,lat_1=truelat1,\
    lat_2=truelat2)

# This sets the standard grid point structure at full resolution
x,y = m(nc.variables['XLONG'][0],nc.variables['XLAT'][0])
xlat,ylons = m(-97.37,35.16)
# This sets a thinn-ed out grid point structure for plotting
# wind barbs at the interval specified in "thin"
x_th,y_th = m(nc.variables['XLONG'][0,::thin,::thin],\
	nc.variables['XLAT'][0,::thin,::thin])


# Set universal figure margins
width = 10
height = 8

plt.figure(figsize=(width,height))
plt.rc("figure.subplot", left = .001)
plt.rc("figure.subplot", right = .999)
plt.rc("figure.subplot", bottom = .001)
plt.rc("figure.subplot", top = .999)





def timestring(wrftime,curtime):
    curtime_str = '%02.0f' % curtime
    wrfdt = datetime.strptime(wrftime,'%Y-%m-%d_%H:%M:%S')
    outtime = '%sZ F%s' % (wrfdt.strftime('%a %Y%m%d/%H%M'),curtime_str)
    return outtime

def drawmap(DATA,TITLESTRING,PROD,UNITS):
    F = plt.gcf()  # Gets the current figure

    m.drawstates(color='k', linewidth=1.25)
    m.drawcoastlines(color='k')
    m.drawcountries(color='k', linewidth=1.25)
    m.scatter(xlat,ylons)
	#m.readshapefile(shapefile='/data/geog/shapefiles/fe_2007_40_county.shp',name='COUNTY',drawbounds='True')
	#m.readshapefile(shapefile='/data/geog/shapefiles/fe_2007_48_county.shp',name='COUNTY',drawbounds='True')
	#plt.suptitle('%s' % UNITS, fontsize = 11, x = 0.08, y = 0.105)
    plt.title('UML WRF-ARW %s (%s)   Valid: %s' % (TITLESTRING, UNITS, curtimestring), \
		fontsize=11,bbox=dict(facecolor='white', alpha=0.65),\
		x=0.5,y=.95,weight = 'demibold',style='oblique', \
		stretch='normal', family='sans-serif')

    # Code to make the colorbar outside of the main axis, on the bottom, and lined up
    ax = plt.gca()  # Gets the current axes
    divider = make_axes_locatable(ax)  # Lets us move axes around
    #cax = divider.append_axes("bottom", size="2%",pad=-0.02,axes_class=maxes.Axes) # Adds an axis for the colorbar
    #F.add_axes(cax)  # Adds the new axis to the figure as the current working axis
    bar = plt.colorbar(DATA,orientation='horizontal',format='%4.2f',extend='both') # Plots colorbar in new axis 
    bar.ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(base=1.0)) # Make the colorbars numbers nice
    bar.update_ticks()

    file_id = '%s_%s_f%02d' % (dom, PROD, time+restart_time)
    filename = '%s.png' % (file_id)
	
    plt.savefig(filename,bbox_inches='tight') # Saves the figure with small margins
    plt.close()

    #if export_flag == 1:
    # Convert the figure to a gif file
    os.system('convert -render -flatten %s %s.gif' % (filename, file_id))
    #os.system('rm -f %s' % filename)

def kuchera():
	#load the variables and interpolate the temperature between levels from surface to 500 mb
	tk = getvar(nc,"tk", timeidx = time)
	p = getvar(nc,"pressure", timeidx=time)
	tk500 = interplevel(tk,p,500)
	tk600 = interplevel(tk,p,600)
	tk700 = interplevel(tk,p,700)
	tk850 = interplevel(tk,p,850)
	tk925 = interplevel(tk,p,925)
	t2 = getvar(nc, "T2", timeidx=time)
	overallarray = np.zeros((len(tk500[:,0]), len(tk500[0,:])))
	#find the max temp at each gridpoint
	for i in range(len(tk500[:,0])):
		for j in range(len(tk500[0,:])):
			tmp = np.zeros(6)
			tmp[0] = tk500[i,j]
			tmp[1] = tk600[i,j]
			tmp[2] = tk700[i,j]
			tmp[3] = tk850[i,j]
			tmp[4] = tk925[i,j]
			tmp[5] = t2[i,j]
			x = tmp.max()
			if( x > 271.16):
				x = 12 + (2.*(271.16-x))
			else:
				x = 12 + (1. *(271.16-x))
			overallarray[i,j] = x
	return overallarray


def twometertemp():
	t2 = getvar(nc, "T2", timeidx=time)
	cart_proj = get_cartopy(t2)
	lats, lons = latlon_coords(t2)
	fig = plt.figure(figsize=(10,7.5))
	geo_axes = plt.axes(projection=cart_proj)
	states = NaturalEarthFeature(category = 'cultural',
				    scale = '50m', 
				    facecolor = 'none',
				    name = 'admin_1_states_provinces_shp')
	geo_axes.add_feature(states,linewidth=1.,edgecolor="black")
	geo_axes.coastlines('50m',linewidth=0.8)
	geo_axes.gridlines(color="black", linestyle="dotted")
	t2_levels = T_LEVS = range(-10,100,2)
	# Manually setting the t2 RGB colors (normalized to 1)
	t2_rgb = np.array([[181,82,0], [181,82,0],
                  [198,107,8], [206,107,8],
                  [231,140,8], [239,156,8],
                  [247,173,24], [255,189,41],
                  [255,212,49], [255,222,66],
                  [255,239,90], [247,255,123],
                  [214,255,132], [181,231,148],
                  [156,222,156], [132,222,132],
                  [112,222,112], [82,222,82],
                  [57,222,57], [33,222,33],
                  [8,206,8], [0,165,0],
                  [0,140,0]]) / 255.0
    
	#t2_cmap, t2_norm = from_levels_and_colors(t2_levels, t2_rgb, extend="both")
	plt.contour(to_np(lons), to_np(lats),
		      to_np((t2-273.15)*(9/5)+32.), levels = t2_levels,
		      colors="black",
		      transform = crs.PlateCarree())	
	plt.contourf(to_np(lons), to_np(lats),
		      to_np((t2-273.15)*(9/5)+32.), levels = t2_levels,
		      cmap=coltbls.sftemp(),
		      extend="both",
		      transform = crs.PlateCarree())
	bar = plt.colorbar(orientation='horizontal',format='%4.2f',extend='both',shrink=.86)
        plt.title('UWML WRF-ARW 2-Meter Temperature Valid: %s' % (curtimestring), \
		fontsize=11,bbox=dict(facecolor='white', alpha=0.65),\
		x=0.5,y=.95,weight = 'demibold',style='oblique', \
		stretch='normal', family='sans-serif')

    	# Code to make the colorbar outside of the main axis, on the bottom, and lined up
    	ax = plt.gca()  # Gets the current axes
    	divider = make_axes_locatable(ax)  # Lets us move axes around
    	#cax = divider.append_axes("bottom", size="2%",pad=-0.02,axes_class=maxes.Axes) # Adds an axis for the colorbar
    	#F.add_axes(cax)  # Adds the new axis to the figure as the current working axis 
    	bar.ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(base=1.0)) # Make the colorbars numbers nice
    	bar.update_ticks()
	filename = "2metertemp_f%02d" % (time+restart_time)
	plt.savefig(filename,bbox_inches='tight') # Saves the figure with small margins



def plot_sim_reflect():
    print("    SIM REFLECTIVITY")
    # Set Figure Size (1000 x 800)
    plt.figure(figsize=(width,height),frameon=False)
    QR = nc.variables['QRAIN']
    try:
        QS = nc.variables['QSNOW']
    except:
        QS = np.zeros(np.shape(QR))
    # Define 'constant' densities (kg m-3)
    rhor = 1000
    rhos = 100
    rhog = 400
    rhoi = 917

    # Define "fixed intercepts" (m-4)
    Norain = 8.0E6
    #Nosnow = 2.0E7
    Nosnow = 2.0E6*np.exp(-0.12 * (temps[time]-273))
    Nograu = 4.0E6


    # First, find the density at the first sigma level
    # above the surface
    density = np.divide(psfc[time],(287.0 * temps[time]))
    #print "Rho: ", np.mean(density)
    Qra = QR[time,1]
    Qsn = QS[time,1]
    Qra = np.nan_to_num(Qra)
    Qsn = np.nan_to_num(Qsn)

    # Calculate slope factor lambda
    lambr = np.divide((3.14159 * Norain * rhor), np.multiply(density, Qra))
    lambr = lambr ** 0.25

    #lambs = np.divide((3.14159 * Nosnow * rhoi), np.multiply(density, Qsn))
    #lambs = lambs ** 0.25
    lambs = np.exp(-0.0536 * (temps[time] - 273))
	
    # Calculate equivalent reflectivity factor
    Zer = (720.0 * Norain * (lambr ** -7.0)) * 1E18
    Zes = (0.224 * 720.0 * Nosnow * (lambr ** -7.0) * (rhos/rhoi) ** 2) * 1E18
    Zes_int = np.divide((lambs * Qsn * density), Nosnow)
    Zes = ((0.224 * 720 * 1E18) / (3.14159 * rhor) ** 2) * Zes_int ** 2 



    Ze = np.add(Zer, Zes)
    #Ze = Zer
    # Convert to dBZ
    
    dBZ = 10 * np.log10(Ze)	
    dBZ = np.nan_to_num(dBZ)
    units = 'dBZe'
    print "      MAX: ", np.max(dBZ)
    # Now plot
    REF_LEVELS = range(5,90,5)
    SREFLECT=plt.contourf(x,y,dBZ,REF_LEVELS,cmap='jet')
    #SREFLECT=plt.contourf(x,y,dBZ)
    title = 'Simulated Surface Reflectivity'
    prodid = 'sref'

    drawmap(SREFLECT, title, prodid, units) 	

def plot_comp_reflect():
    print("    COMP REFLECTIVITY")
    # Set Figure Size (1000 x 800)
    plt.figure(figsize=(width,height),frameon=False)
    QR = nc.variables['QRAIN']
    try:
        QS = nc.variables['QSNOW']
    except:
        QS = np.zeros(np.shape(QR))

    # Define 'constant' densities (kg m-3)
    rhor = 1000
    rhos = 100
    rhog = 400
    rhoi = 917

    # Define "fixed intercepts" (m-4)
    Norain = 8.0E6
    #Nosnow = 2.0E7
    Nosnow = 2.0E6*np.exp(-0.12 * (temps[time]-273))
    Nograu = 4.0E6


    # First, find the density at the first sigma level
    # above the surface
    density = np.divide(psfc[time],(287.0 * temps[time]))
    #print "Rho: ", np.mean(density)
    Qra_all = QR[time]
    Qsn_all = QS[time]

    for j in range(len(Qra_all[1,:,1])):
		curcol_r = []
		curcol_s = []
		for i in range(len(Qra_all[1,1,:])):
				maxrval = np.max(Qra_all[:,j,i])
				maxsval = np.max(Qsn_all[:,j,i])
				curcol_r.append(maxrval)		
				curcol_s.append(maxsval)
		np_curcol_r = np.array(curcol_r)
		np_curcol_s = np.array(curcol_s)
		if j == 0:
			Qra = np_curcol_r
			Qsn = np_curcol_s
		else:
			Qra = np.row_stack((Qra, np_curcol_r))
			Qsn = np.row_stack((Qsn, np_curcol_s))

    #print "Qra shp: ", np.shape(Qra)
    #print "Den shp: ", np.shape(density)
	


    # Calculate slope factor lambda
    lambr = np.divide((3.14159 * Norain * rhor), np.multiply(density, Qra))
    lambr = lambr ** 0.25

    #lambs = np.divide((3.14159 * Nosnow * rhoi), np.multiply(density, Qsn))
    #lambs = lambs ** 0.25
    lambs = np.exp(-0.0536 * (temps[time] - 273))
	
    # Calculate equivalent reflectivity factor
    Zer = (720.0 * Norain * (lambr ** -7.0)) * 1E18
    Zes = (0.224 * 720.0 * Nosnow * (lambr ** -7.0) * (rhos/rhoi) ** 2) * 1E18
    Zes_int = np.divide((lambs * Qsn * density), Nosnow)
    Zes = ((0.224 * 720 * 1E18) / (3.14159 * rhor) ** 2) * Zes_int ** 2 



    Ze = np.add(Zer, Zes)
    #Ze = Zer
    # Convert to dBZ
    dBZ = 10 * np.log10(Ze)	
    dBZ = np.nan_to_num(dBZ)
    units = 'dBZe'
    print "      MAX: ", np.max(dBZ)
    # Now plot
    REF_LEVELS = range(5,90,5)
    CREFLECT=plt.contourf(x,y,dBZ,REF_LEVELS,cmap='jet')
    #SREFLECT=plt.contourf(x,y,dBZ)
    title = 'Simulated Composite Reflectivity'
    prodid = 'cref'

    drawmap(CREFLECT, title, prodid, units) 	

def plot_precip(runtotal):
    print("    PRECIP")
    # Set Figure Size (1000 x 800)
    plt.figure(figsize=(width,height),frameon=False)
    rainc =  nc.variables['RAINC']
    rainnc = nc.variables['RAINNC']

    # First, find out if this is first time or not
    # Based on skip.  This should be total from each output time
    if time == 0:
		prev_total = rainc[time] + rainnc[time]
    else:
		prev_total = rainc[time-1] + rainnc[time-1]
    total_accum = rainc[time] + rainnc[time]
    precip_tend = total_accum
    #runtotal = runtotal + total_accum
	
    # Convert from mm to in
    precip_tend = total_accum * .0393700787
    runtotal = runtotal + precip_tend
    units = 'in'
    #PCP_LEVELS = [0.01,0.03,0.05,0.10,0.15,0.20,0.25,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.00,1.25,1.50,1.75,2.00,2.50]
    PRECIP=plt.contourf(x,y,runtotal,cmap='jet')
    #plt.jet()
    title = 'Total Precip'
    prodid = 'precip'

    drawmap(PRECIP, title, prodid, units) 	
    


def plot_pwat():
    print("    PRECIP. WATER")
    # Set Figure Size (1000 x 800)
    plt.figure(figsize=(width,height),frameon=False)
    g = 9.81
    P = nc.variables['P']
    PB = nc.variables['PB']
    Qv = nc.variables['QVAPOR']
    
    # First we need an array of pressures
    Pr = P[time] + PB[time]
    print "SHAPE Pr: ", np.shape(Pr)
    
    Qvap = Qv[time,:,:,:]
    print "SHAPE Qvap: ", np.shape(Qvap)

    # Now go through each point
    for j in range(len(Pr[0])):
		currow_pwat = []
		for i in range(len(Pr[0,0])):
			curcol_pwat = []		
			for k in range(len(Pr)-1):
				curdp = (Pr[k,j,i] - Pr[k+1,j,i]) * 0.01
				curpwat = curdp * Qvap[k,j,i]
				curcol_pwat.append(curpwat)		
			np_curcol_pwat = np.array(curcol_pwat)
			point_pwat = np.sum(curcol_pwat)
			currow_pwat.append(point_pwat)		
		
		np_currow_pwat = np.array(currow_pwat)	
		if j == 0:
			total_pwat = np_currow_pwat
		else:
			total_pwat = np.row_stack((np_currow_pwat, total_pwat))

		
    pwat = np.divide(total_pwat,g)
    print "Len j: ", len(Pr)
    print "Len i: ", len(Pr[0])

    print "SHAPE: ", np.shape(pwat)
    print "MAX: ", np.max(pwat)	
    PCP_LEVELS = [0.01,0.03,0.05,0.10,0.15,0.20,0.25,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.00,1.25,1.50,1.75,2.00,2.50]
    PWAT_LEVS = [0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.1,2.3,2.4,2.5]

    PWAT=plt.contourf(x,y,pwat,PCP_LEVELS,cmap='jet')
    units = 'mm'
    title = 'Precipitable Water (in)'
    prodid = 'pwat'

    drawmap(PWAT, title, prodid, units) 	
    


def plot_snowfall(runtotal1):
    print("    SNOWFALL")
    # Set Figure Size (1000 x 800)
    plt.figure(figsize=(width,height),frameon=False)
    
    snownc = nc.variables['SNOWNC']
    snowh = nc.variables['SNOWH']
    graupel = nc.variables['GRAUPELNC']
    # First, find out if this is first time or not
	
    #if time == 0:
    #	prev_total = snownc[time]
    #else:
    #	prev_total = snownc[time-1]
    snow_accum = np.multiply(snownc[time]-graupel[time], .0393701)
    #precip_tend = total_accum  - prev_total
    snowlevs = np.arange(0.,25.,1)

    ratio = kuchera()
    runtotal1 = runtotal1 + snow_accum
    runtotal1 = np.multiply(runtotal1, ratio)
    snow_prev = np.multiply(snowh[time], 39.37)
	
    #SNOW_LEVS = range(1,20,1)
    #SNOWP_LEVS = [0.25,0.5,0.75,1,1.5,2,2.5,3,4,5,6,8,10,12,14,16,18]

    #plt.cool()
    #plt.cool()
    SNOWF=plt.contourf(x,y,runtotal1, levels=snowlevs, extend='max',cmap=coltbls.snow2())
    #SNOWP=plt.contour(x,y,snow_prev,SNOWP_LEVS, extend='max', cmap=coltbls.grays(), linewidth = 0.75)
    title = 'Total Accum. Snowfall (in.)'
    prodid = 'snow'
    units = 'in'

    drawmap(SNOWF, title, prodid, units) 	

def plot_precip_type():
    print("    PRECIP TYPE")
    # Set Figure Size (1000 x 800)
    plt.figure(figsize=(width,height),frameon=False)

    sr = nc.variables['SR']
    tsk = nc.variables['TSK']
    rainc =  nc.variables['RAINC']
    rainnc = nc.variables['RAINNC']


    type_pct = sr[time]
    if time == 0:
		prev_total = rainc[time] + rainnc[time]
    else:
		prev_total = rainc[time-1] + rainnc[time-1]
    total_accum = rainc[time] + rainnc[time]
    precip_tend = total_accum  - prev_total
	
    snow_precip = []
    mix_precip = []
    rain_precip = []
    for j in range(len(precip_tend)):
		cur_col_rain = []
		cur_col_mix = []
		cur_col_snow = []
		for i in range(len(precip_tend[0])):
			if (0.20 < type_pct[j,i] < 0.90):
				cur_col_mix.append(precip_tend[j,i])
				cur_col_snow.append(0.)
				cur_col_rain.append(0.)	
			elif (type_pct[j,i] >= 0.90):
				cur_col_mix.append(0.)
				cur_col_snow.append(precip_tend[j,i])
				cur_col_rain.append(0.)	
				#print type_pct[j,i]	

			else:
				cur_col_mix.append(0.)
				cur_col_snow.append(0.)
				cur_col_rain.append(precip_tend[j,i])
		snow_precip.append(cur_col_snow)
		mix_precip.append(cur_col_mix)
		rain_precip.append(cur_col_rain)		
    #print snow_precip
    #raw_input()
    PCP_LEVELS = [0.01,0.03,0.05,0.10,0.15,0.20,0.25,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.00,1.25,1.50,1.75,2.00,2.50]

    MIXT=plt.contourf(x,y,mix_precip,PCP_LEVELS,extend='max',cmap=coltbls.mixprecip1())
    RAINT=plt.contourf(x,y,rain_precip,PCP_LEVELS,extend='max',cmap=coltbls.rain1())
    SNOWT=plt.contourf(x,y,snow_precip,PCP_LEVELS,extend='max',cmap=coltbls.snow1())

    ftemps = (9./5.)*(temps[time]-273) + 32
    ftsk = (9./5.)*(tsk[time] - 273) + 32
    T=plt.contour(x,y,ftemps,[32],colors='red',linestyles='solid')
    TS=plt.contour(x,y,ftsk,[32],colors='purple',linestyles='dashdot')

    #RAINT=plt.contourf(x,y,type_pct,extend='max',cmap=matplotlib.cm.copper)

    title = 'Frozen Precipitation'
    prodid = 'ptype'
    units = 'in'

    drawmap(RAINT, title, prodid, 'in') 	





def plot_swdown():
    print("    SWDOWN")
    # Set Figure Size (1000 x 800)
    plt.figure(figsize=(width,height),frameon=False)
    swdown = nc.variables['SWDOWN']

	
    SWDOWN=plt.contourf(x,y,swdown[time],cmap=matplotlib.cm.bone_r)
    title = 'Shortwave Radiation at Sfc'
    prodid = 'swdown'
    units = 'W/m\xb2'

    drawmap(SWDOWN, title, prodid, units) 	



def plot_olr():
    print("    OLR/IRSAT")
    # Set Figure Size (1000 x 800)
    plt.figure(figsize=(width,height),frameon=False)
    olr = nc.variables['OLR']
	
    sbc = .000000056704	
    ir_T = ((olr[time] / sbc) ** (0.25)) - 273

    IR_LEVS = [-100,-90,-85,-80,-75,-70,-65,-60,-55,-50,-46,-42,-38,-36,\
		-34,-32,-30,-28,-26,-24,-22,-20,-18,-16,-14,-12,-10,\
		-8,-6,-4,-2,0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,\
		30,32,34,36,38,40,42,44,46,48,50,52]

    #OLR=plt.contourf(x,y,ir_T,IR_LEVS,cmap=matplotlib.cm.spectral_r)
    OLR=plt.contourf(x,y,ir_T,IR_LEVS,cmap=coltbls.irsat())
    #OLR=plt.contourf(x,y,olr[time],range(95,350,5),cmap=coltbls.bw_irsat())
    title = 'TOA Inferred Temperature'
    prodid = 'olr'
    units = u"\u00B0" + "C"	

    drawmap(OLR, title, prodid, units) 	



def plot_surface():
    print("    SURFACE")
    # Set Figure Size (1000 x 800)
    plt.figure(figsize=(width,height),frameon=False)


    # Convert Surface Pressure to Mean Sea Level Pressure	
    stemps = temps[time]+6.5*nc.variables['HGT'][time]/1000.
    mslp = nc.variables['PSFC'][time]*np.exp(9.81/(287.0*stemps)*nc.variables['HGT'][time])*0.01 + (6.7 * nc.variables['HGT'][time] / 1000)

    # Convert Celsius Temps to Fahrenheit
    ftemps = (9./5.)*(temps[time]-273) + 32


    T_LEVS = range(-10,125,5)

    # Contour and fill the temperature
    T=plt.contourf(x,y,ftemps,T_LEVS,cmap=coltbls.sftemp())

    # Contour the pressure
    P=plt.contour(x,y,mslp,V=2,colors='k',linewidths=1.5)
    plt.clabel(P,inline=1,fontsize=8,fmt='%1.0f',inline_spacing=1)

    #plt.clabel(T,inline=1,fontsize=10)

    # Convert winds from m/s to kts and then draw barbs	
    u_wind_kts = u_wind_ms[time] * 1.94384449
    v_wind_kts = v_wind_ms[time] * 1.94384449
    plt.barbs(x_th,y_th,u_wind_kts[::thin,::thin],\
		v_wind_kts[::thin,::thin], length=5,\
		sizes={'spacing':0.2},pivot='middle')

    title = 'Sfc Temp, MSLP (mb), 10m Wind (kts)'
    prodid = 'pmsl'
    units = u"\u00B0" + "F"	

    drawmap(T, title, prodid, units)

def plot_sfwind():
    print("    10M WIND")
    # Set Figure Size (1000 x 800)
    plt.figure(figsize=(width,height),frameon=False)

    # Convert winds from m/s to kts and then draw barbs	
    u_wind_kts = u_wind_ms[time] * 1.94384449
    v_wind_kts = v_wind_ms[time] * 1.94384449
    windmag = np.power(np.power(u_wind_kts,2)+np.power(v_wind_kts,2), 0.5)
    WIND_LEVS = range(10,46,2)
    W=plt.contourf(x,y,windmag,WIND_LEVS,extend='max')

    plt.barbs(x_th,y_th,u_wind_kts[::thin,::thin],\
		v_wind_kts[::thin,::thin], length=5,\
		sizes={'spacing':0.2},pivot='middle')



    # Convert Surface Pressure to Mean Sea Level Pressure	
    stemps = temps[time]+6.5*nc.variables['HGT'][time]/1000.
    mslp = nc.variables['PSFC'][time]*np.exp(9.81/(287.0*stemps)*nc.variables['HGT'][time])*0.01 + (6.7 * nc.variables['HGT'][time] / 1000)

    # Contour the pressure
    #PLEVS = range(900,1050,5)
    #P=plt.contour(x,y,mslp,PLEVS,V=2,colors='k',linewidths=1.5)
    #plt.clabel(P,inline=1,fontsize=8,fmt='%1.0f',inline_spacing=1)



    title = 'Sfc MSLP (mb), 10m Wind (kts)'
    prodid = 'wind'
    units = "kts"	

    drawmap(W, title, prodid, units)

def plot_dwp():
    print "    DEWPOINT"
    # Set Figure Size (1000 x 800)
    plt.figure(figsize=(width,height),frameon=False)
    qhum = nc.variables['Q2']
	
	
    # Convert Surface Pressure to Mean Sea Level Pressure
    stemps = temps[time]+6.5*nc.variables['HGT'][time]/1000.
    mslp = psfc[time]*np.exp(9.81/(287.0*stemps)*nc.variables['HGT'][time])*0.01
    # Find saturation vapor pressure
    es = 6.112 * np.exp(17.67 * temps[time]/(temps[time] + 243.5))
    w = qhum[time]/(1-qhum[time])
    e = (w * psfc[time] / (.622 + w)) / 100
    Td_C = (243.5 * np.log(e/6.112))/(17.67-np.log(e/6.112))
    Td_F = (Td_C * 9 / 5) + 32


    DP_LEVS = range(-10,85,1)
    DP_CLEVS = range(40,90,10)

    # Contour and fill the dewpoint temperature		
    Td=plt.contourf(x,y,Td_F,DP_LEVS,cmap=coltbls.dewpoint1(),extend='min')
    Td_lev = plt.contour(x,y,Td_F,DP_CLEVS,colors='k',linewidths=.5)
    plt.clabel(Td_lev,inline=1,fontsize=7,fmt='%1.0f',inline_spacing=1)

    # Contour the pressure
    # P=plt.contour(x,y,mslp,V=2,colors='k',linewidths=1.5)
    # plt.clabel(P,inline=1,fontsize=8,fmt='%1.0f',inline_spacing=1)
    
    #plt.clabel(T,inline=1,fontsize=10)

    # Convert winds from m/s to kts and then draw barbs	
    u_wind_kts = u_wind_ms[time] * 1.94384449
    v_wind_kts = v_wind_ms[time] * 1.94384449
    plt.barbs(x_th,y_th,u_wind_kts[::thin,::thin],\
		v_wind_kts[::thin,::thin], length=5,\
		sizes={'spacing':0.2},pivot='middle')

    title = 'Surface Dwp, 10m Wind (kts)'
    prodid = 'dewp'
    units = u"\u00B0" + "F"	

    drawmap(Td, title, prodid, units)


def plot_zlcl():
    print "    ZLCL"
    # Set Figure Size (1000 x 800)
    plt.figure(figsize=(width,height),frameon=False)
    qhum = nc.variables['Q2']

    # Find saturation vapor pressure
    es = 6.112 * np.exp(17.67 * temps[time]/(temps[time] + 243.5))
    w = qhum[time]/(1-qhum[time])
    e = (w * psfc[time] / (.622 + w)) / 100
    Td_C = (243.5 * np.log(e/6.112))/(17.67-np.log(e/6.112))
    Td_F = (Td_C * 9 / 5) + 32

    # Calculate the LCL height
    z_lcl = 125.0 * np.subtract((temps[time]-273),Td_C)

    # Contour and fill the dewpoint temperature		
    ZLCL=plt.contourf(x,y,z_lcl)

    title = 'LCL Height'
    prodid = 'zlcl'
    units = 'm'

    drawmap(ZLCL, title, prodid)
		
def plot_thte():
    """Plot surface theta-e map"""
    print "    THETA-E"
    plt.figure(figsize=(width,height),frameon=False)
    qhum = nc.variables['Q2']
    
    thte = (temps[time] + qhum[time] * 2500000.0/1004.0) * (100000/psfc[time]) ** (287.0/1004.0) 	
    THTE_LEVS = range(270,360,5)
    THTE = plt.contourf(x,y,thte,THTE_LEVS,cmap=coltbls.thetae(),extend='max')
    
    u_wind_kts = u_wind_ms[time] * 1.94384449
    v_wind_kts = v_wind_ms[time] * 1.94384449
    plt.barbs(x_th,y_th,u_wind_kts[::thin,::thin],\
        v_wind_kts[::thin,::thin], length=5,\
        sizes={'spacing':0.2},pivot='middle')

    title = 'Theta-e, 10 m Wind (kt)'
    prodid = 'thte'
    units = 'K'

    drawmap(THTE, title, prodid, units)

def plot_cape():
    print("    CAPE")
    # Set Figure Size (1000 x 800)
    plt.figure(figsize=(width,height),frameon=False)
    P = nc.variables['P']
    PB = nc.variables['PB']
    Qv = nc.variables['QVAPOR']

    # Need pressures, temps and mixing ratios
    PR = P[time] + PB[time]
    W = Qv[time]/(1-Qv[time])
    TH = np.add(T[time],290.)
    T_K = np.multiply(TH, np.power(np.divide(PR,1000.),(287.04/1004.)))
    PR_h = PR / 100.
	
    print "PR: ", np.shape(PR)
    print "W: ", np.shape(W)
    print "T_K: ", np.shape(T_K)


    for j in range(len(T_K[1,:,1])):
		curcol_c = []
		for i in range(len(T_K[1,1,:])):
				sparms = severe.CAPESOUND(PR_h[:,j,i],T_K[:,j,i],W[:,j,i])
				curcol_c.append(sparms[1])		
		np_curcol_c = np.array(curcol_c)
		if j == 0:
			cape = np_curcol_c
		else:
			cape = np.row_stack((cape, np_curcol_c))

    print "CAPE: ", np.shape(cape)
    # Now plot
    CAPE_LEVS = range(500,6000,250)
    SCAPE=plt.contourf(x,y,cape,CAPE_LEVS)
    #SREFLECT=plt.contourf(x,y,dBZ)
    title = 'SBCAPE (J/kg)'
    prodid = 'cape'
    units = 'J/Kg'

    drawmap(SCAPE, title, prodid, units) 	

def plot_srhel():
    print("    SR Helicity")
    # Set Figure Size (1000 x 800)
    plt.figure(figsize=(width,height),frameon=False)
    P = nc.variables['P']
    PB = nc.variables['PB']
    UU = nc.variables['U']
    VV = nc.variables['V']
    PH = nc.variables['PH']
    PHB = nc.variables['PHB']
    
    # Need pressures, temps and mixing ratios
    PR = P[time] + PB[time]
    PHT = np.add(PH[time],PHB[time])
    ZH = np.divide(PHT, 9.81)
    U = UU[time]
    V = VV[time]
    

    for j in range(len(U[1,:,1])):
		curcol_c = []
		curcol_Umo = []
		curcol_Vmo = []
		for i in range(len(V[1,1,:])):
				sparms = severe.SRHEL_CALC(U[:,j,i], V[:,j,i], ZH[:,j,i], PR[:,j,i])
				curcol_c.append(sparms[0])		
				curcol_Umo.append(sparms[1])
				curcol_Vmo.append(sparms[2])
		np_curcol_c = np.array(curcol_c)
		np_curcol_Umo = np.array(curcol_Umo)
		np_curcol_Vmo = np.array(curcol_Vmo)

		if j == 0:
			srhel = np_curcol_c
			U_srm = np_curcol_Umo
			V_srm = np_curcol_Vmo
		else:
			srhel = np.row_stack((srhel, np_curcol_c))
			U_srm = np.row_stack((U_srm, np_curcol_Umo))
			V_srm = np.row_stack((V_srm, np_curcol_Vmo))

    #print "       SRHEL: ", np.shape(srhel)

    # Now plot
    SRHEL_LEVS = range(50,800,50)
    srhel = np.nan_to_num(srhel)
    SRHEL=plt.contourf(x,y,srhel,SRHEL_LEVS)

    u_mo_kts = U_srm * 1.94384449
    v_mo_kts = V_srm * 1.94384449
    plt.barbs(x_th,y_th,u_mo_kts[::thin,::thin],\
		v_mo_kts[::thin,::thin], length=5,\
		sizes={'spacing':0.2},pivot='middle')
    title = '0-3 km SRHelicity, Storm Motion (kt)'
    prodid = 'hlcy'
    units = "m" + u'\u00B2' + '/s' + u'\u00B2'

    drawmap(SRHEL, title, prodid, units) 	

def plot_plcl():
    print("    LCL Pressure")
    # Set Figure Size (1000 x 800)
    plt.figure(figsize=(width,height),frameon=False)
    P = nc.variables['P']
    PB = nc.variables['PB']
    Qv = nc.variables['QVAPOR']

    # Need pressures, temps and mixing ratios
    PR = P[time][0] + PB[time][0]
    W = Qv[time][0]/(1-Qv[time][0])
    TH = np.add(T[time][0],290.)
    T_K = np.multiply(TH, np.power(np.divide(PR,1000.),(287.04/1004.)))
    PR_h = PR / 100.


    plcl = severe.PLCL_CALC(PR_h, T_K, W)

    plcl = np.nan_to_num(plcl)

    PLCL=plt.contourf(x,y,plcl)

    title = 'Surface-Based LCL Pressure'
    prodid = 'plcl'
    units = "hPa"

    drawmap(PLCL, title, prodid, units) 	



# Check to see if we are exporting
if export_flag == 1:
    dom = 'wrf'
# Begin looping through times
for time in range(len(temps[:,0,0])):
    print 'Plotting time ',time*skip+restart_time

    curtimestring = timestring(''.join(times[time]),time*skip+restart_time)
	
    if var == 'temp':
		plot_surface()
    elif var == 'precip':
		plot_precip(runtotal)
    elif var == 'dwp':
		plot_dwp()
    elif var == 'zlcl':
		plot_zlcl()
    elif var == 'thte':
		plot_thte()
    elif var == 'snow':
		plot_snowfall(runtotal)
    elif var == 'swdown':
		plot_swdown()
    elif var == 'olr':
		plot_olr()
    elif var == 'ptype':
		plot_precip_type()
    elif var == 'sref':
		plot_sim_reflect()
    elif var == 'cref':
		plot_comp_reflect()
    elif var == 'pwat':
		plot_pwat()
    elif var == 'cape':
		plot_cape()
    elif var == 'srhel':
		plot_srhel()
    elif var == 'wind':
		plot_sfwind()
    elif var == 'temp2':
		twometertemp()
    else:
        plot_surface()
        plot_precip(runtotal)
        #plot_dwp()
        #plot_comp_reflect()
        #plot_sfwind()
        plot_sim_reflect()
	plot_comp_reflect()
        #plot_thte()
        #plot_cape()
        plot_snowfall(runtotal1)
        plot_precip_type()
        #plot_swdown()
        #plot_olr()
	twometertemp()

os.system("convert -delay 40 -loop 0 wrf_snow*.png snow.gif")
os.system("convert -delay 40 -loop 0 wrf_precip*.png precip.gif")
os.system("convert -delay 40 -loop 0 wrf_pmsl*.png pmsl.gif")
os.system("convert -delay 40 -loop 0 wrf_sref*.png sref.gif")
os.system("convert -delay 40 -loop 0 wrf_cref*.png cref.gif")
os.system("convert -delay 40 -loop 0 wrf_ptype*.png ptype.gif")
os.system("convert -delay 40 -loop 0 2metertemp*.png 2meter.gif")

