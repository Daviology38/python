import numpy as np
from matplotlib import pyplot
from matplotlib.cm import get_cmap
from matplotlib.colors import from_levels_and_colors
from cartopy import crs
from cartopy.feature import NaturalEarthFeature
from netCDF4 import Dataset
from wrf import getvar, to_np, get_cartopy, latlon_coords, cartopy_xlim, cartopy_ylim, interplevel, smooth2d

pyplot.switch_backend('agg')

#List of potential variables in wrf file
variables = ['UV wind components (m/s)', 'W z-wind component (m/s)', ' H Geopotential Height (m2/s2)','T2 2-meter temp (K)', 'PSFC Surface Pressure (Pa)','UV10 UV 10-meter wind (m/s)','H_DIABATIC Microphysics Latent Heat (K/s)','RTHCUTEN Theta tendency (PA K /s)', 'RQVCUTEN coupled Q_V Tendency ( (Pa*kg) / (kg*s) )', 'Rain (mm)','tc (C)', 'tk (K)', 'height model height (km)','pressure (hPa)', 'td Dew point (C)','rh (%)']

#Import the Data file
ncfile = Dataset("wrfout_d01_2019-01-25_00:00:00")

#Ask if want to make a composite comparison
comp = raw_input("Do you want to make a height compost    y/n    ")

if( comp.lower() == "y"):
	time = 0
	endtime = 17
	while time <= endtime:
		year = raw_input("Enter a year  yyyy       ")
		m = raw_input("Enter the month  mm    ")
		d1 = raw_input("Enter the first day    dd       ")
		d2 = raw_input("Enter the second day    dd      ")
		dim = raw_input("Which dimension?  Enter d01 or d02        ")
	
		ncfile2 = "wrfout_" + dim + "_" + year + "-" + m + "-" + d1 + "_00:00:00"
		ncfile3 = "wrfout_" + dim + "_" + year + "-" + m + "-" + d2 + "_00:00:00"
		lev = 500
		h = getvar(ncfile2, "z", timeidx=time, units="dm")
        	p = getvar(ncfile2,"pressure", timeidx=time)
        	ht2 = interplevel(h,p,lev)
		h = getvar(ncfile3, "z", timeidx=time, units="dm")
        	p = getvar(ncfile3,"pressure", timeidx=time)
        	ht3 = interplevel(h,p,lev)
		htnew = ht3 - ht2 
        	cart_proj = get_cartopy(htnew)
        	lats, lons = latlon_coords(htnew)
        	fig = pyplot.figure(figsize=(10,7.5))
        	geo_axes = pyplot.axes(projection=cart_proj)
        	# Make the contour lines and fill them.
                contours = pyplot.contour(to_np(lons), to_np(lats), 
                               to_np(htnew), 15, colors="black",
                               transform=crs.PlateCarree())
                pyplot.clabel(contours, inline=1, fontsize=10, fmt="%i")
                pyplot.contourf(to_np(lons), to_np(lats), 
                               to_np(htnew), 15, 
                               cmap=get_cmap('jet'),
                               transform=crs.PlateCarree())
                pyplot.colorbar(ax=geo_axes, shrink=.86)
                name = "Heightanom" + str(lev) + str(time) + ".png"
                pyplot.savefig(name)
		time = time + 1

#Ask user for variable

var = raw_input("Which variable to plot? type w to get list of variables    ")

if( var.lower() == "w"):
	print variables 
	var = raw_input("Which variable to plot?    ")

#Ask for input level to analyze at
lev = raw_input("Pick a level to plot at (hPa): 1000, 950, 900, 850, 800, 750, 700, 650, 600, 550, 500, 450, 400, 350, 300, 250, 200, 150, 100    ")

#Enter model end run time here
endtime = 16

#start time at 0
time = 0

#time loop
while time <= endtime:

#Get data to plot
	if(var.lower() == "t2"):
		t2 = getvar(ncfile, "T2", timeidx=time)
		cart_proj = get_cartopy(t2)
		lats, lons = latlon_coords(t2)
		fig = pyplot.figure(figsize=(10,7.5))
		geo_axes = pyplot.axes(projection=cart_proj)
		states = NaturalEarthFeature(category = 'cultural',
				     scale = '50m', 
				     facecolor = 'none',
				     name = 'admin_1_states_provinces_shp')
		geo_axes.add_feature(states,linewidth=1.,edgecolor="black")
		geo_axes.coastlines('50m',linewidth=0.8)
		geo_axes.gridlines(color="black", linestyle="dotted")
		t2_levels = np.arange(int((t2.min()-273.15)*(9/5)+32.), int((t2.max()-273.15)*(9/5)+32.), 5.)
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
		pyplot.contour(to_np(lons), to_np(lats),
		       to_np((t2-273.15)*(9/5)+32.), levels = t2_levels,
		       colors="black",
		       transform = crs.PlateCarree())	
		pyplot.contourf(to_np(lons), to_np(lats),
		       to_np((t2-273.15)*(9/5)+32.), levels = t2_levels,
		       cmap="jet",
		       extend="both",
		       transform = crs.PlateCarree())
		pyplot.colorbar(ax=geo_axes,shrink=.86)
        	name = "t2" + str(time) + ".png"
		titlename = "2 Meter Temperature at " + str(time * 3) + "z"
		pyplot.title(titlename)
		pyplot.savefig(name)  
		time = time + 1

	if(var.lower() == "psfc"):
		slp = getvar(ncfile, "slp", timeidx=time)
		smooth_slp = smooth2d(slp,3)
		cart_proj = get_cartopy(slp)
		lats, lons = latlon_coords(slp)
		fig = pyplot.figure(figsize=(10,7.5))
		geo_axes = pyplot.axes(projection=cart_proj)
		states = NaturalEarthFeature(category = 'cultural',
				     scale = '50m', 
				     facecolor = 'none',
				     name = 'admin_1_states_provinces_shp')
		geo_axes.add_feature(states,linewidth=.5)
		geo_axes.coastlines('50m',linewidth=0.8)
		# Set the contour levels so that all plots match
		#levels = np.arange(960.,1040.,2.5)

		# Make the contour lines and fill them.
		contours = pyplot.contour(to_np(lons), to_np(lats), 
               		       to_np(smooth_slp), 15, colors="black",
                               transform=crs.PlateCarree())
		pyplot.clabel(contours, inline=1, fontsize=10, fmt="%i")
		pyplot.contourf(to_np(lons), to_np(lats), 
                               to_np(smooth_slp), 15, 
			       cmap=get_cmap('jet'),
                               transform=crs.PlateCarree())                              
		pyplot.colorbar(ax=geo_axes, shrink=.86)
		name = "SLP" + str(time) + ".png"
		pyplot.savefig(name)  
		time = time + 1
	
	if(var.lower() == "h"):
		h = getvar(ncfile, "z", timeidx=time, units="dm")
		p = getvar(ncfile,"pressure", timeidx=time)
		ht = interplevel(h,p,lev)
		cart_proj = get_cartopy(ht)
		lats, lons = latlon_coords(ht)
		fig = pyplot.figure(figsize=(10,7.5))
		geo_axes = pyplot.axes(projection=cart_proj)
		states = NaturalEarthFeature(category = 'cultural',
				     scale = '50m', 
				     facecolor = 'none',
				     name = 'admin_1_states_provinces_shp')
		geo_axes.add_feature(states,linewidth=.5)
		geo_axes.coastlines('50m',linewidth=0.8)

		# Make the contour lines and fill them.
		contours = pyplot.contour(to_np(lons), to_np(lats), 
               		       to_np(ht), 15, colors="black",
                               transform=crs.PlateCarree())
		pyplot.clabel(contours, inline=1, fontsize=10, fmt="%i")
		pyplot.contourf(to_np(lons), to_np(lats), 
                               to_np(ht), 15, 
			       cmap=get_cmap('jet'),
                               transform=crs.PlateCarree())                              
		pyplot.colorbar(ax=geo_axes, shrink=.86)
		name = "Height" + str(lev) + str(time) + ".png"
		pyplot.savefig(name)		
		time = time + 1

	if(var.lower() == "uv"):
		u = getvar(ncfile, "ua", timeidx=time, units="kt")
		v = getvar(ncfile, "va", timeidx=time, units="kt")
		wspd = getvar(ncfile, "wspd_wdir", timeidx=time, units = "kt")[0,:]
		z = getvar(ncfile, "z", timeidx=time, units="dm")
		p = getvar(ncfile,"pressure",timeidx=time)
		ht = interplevel(z,p,lev)
		ulev = interplevel(u, p, lev)
		vlev = interplevel(v, p, lev)
		wspdlev = interplevel(wspd, p, lev)
		ht = interplevel(z, p, lev)
		cart_proj = get_cartopy(ht)
		lats, lons = latlon_coords(ht)
		fig = pyplot.figure(figsize=(10,7.5))
		geo_axes = pyplot.axes(projection=cart_proj)
		states = NaturalEarthFeature(category = 'cultural',
				     scale = '50m', 
				     facecolor = 'none',
				     name = 'admin_1_states_provinces_shp')
		geo_axes.add_feature(states,linewidth=.5)
		geo_axes.coastlines('50m',linewidth=0.8)

		# Make the contour lines and fill them.
		contours = pyplot.contour(to_np(lons), to_np(lats), 
               		       to_np(ht), 15, colors="black",
                               transform=crs.PlateCarree())
		pyplot.clabel(contours, inline=1, fontsize=10, fmt="%i")
		levels = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 					     90, 100, 110]    
		wspd_contours = pyplot.contourf(to_np(lons), 
                             to_np(lats), 
                             to_np(wspdlev), 
                             levels=levels,
                             cmap=get_cmap("rainbow"),
                             transform=crs.PlateCarree())

		pyplot.colorbar(wspd_contours, ax=geo_axes, orientation="horizontal", 
             		     pad=.05, shrink=.75)

		
		thin = [int(x/10.) for x in lons.shape]
		pyplot.barbs(to_np(lons[::thin[0], ::thin[1]]), 
          		to_np(lats[::thin[0], ::thin[1]]), 
          		to_np(ulev[::thin[0], ::thin[1]]),
          		to_np(vlev[::thin[0], ::thin[1]]), 
          		length=6,
          		transform=crs.PlateCarree())                      
		name = "Winds" + str(lev) + str(time) + ".png"
		pyplot.savefig(name)		
		time = time + 1
		pyplot.close()

	if(var.lower() == "uv10"):
		slp = getvar(ncfile, "slp", timeidx=time)
		u_sfc = getvar(ncfile, "ua", timeidx = time, units="kt")[0,:]
		v_sfc = getvar(ncfile, "va", timeidx = time, units="kt")[0,:]
		smooth_slp = smooth2d(slp,3)
		cart_proj = get_cartopy(slp)
		lats, lons = latlon_coords(slp)
		fig = pyplot.figure(figsize=(10,7.5))
		geo_axes = pyplot.axes(projection=cart_proj)
		states = NaturalEarthFeature(category = 'cultural',
				     scale = '50m', 
				     facecolor = 'none',
				     name = 'admin_1_states_provinces_shp')
		geo_axes.add_feature(states,linewidth=.5)
		geo_axes.coastlines('50m',linewidth=0.8)
		# Set the contour levels so that all plots match
		#levels = np.arange(960.,1040.,2.5)

		# Make the contour lines and fill them.
		contours = pyplot.contour(to_np(lons), to_np(lats), 
               		       to_np(smooth_slp), 15, colors="black",
                               transform=crs.PlateCarree())
		pyplot.clabel(contours, inline=1, fontsize=10, fmt="%i")
		# Plot the wind barbs, but only plot ~10 barbs in each direction.
		thin = [int(x/10.) for x in lons.shape]
		pyplot.barbs(to_np(lons[::thin[0], ::thin[1]]), 
             			to_np(lats[::thin[0], ::thin[1]]), 
             			to_np(u_sfc[::thin[0], ::thin[1]]), 
             			to_np(v_sfc[::thin[0], ::thin[1]]),
             			transform=crs.PlateCarree())                           
		pyplot.colorbar(ax=geo_axes, shrink=.86)
		name = "UV10m" + str(time) + ".png"
		pyplot.savefig(name)  
		time = time + 1

	if(var.lower() == "h_diabatic"):
		h = getvar(ncfile, "H_DIABATIC", timeidx=time)
		p = getvar(ncfile,"pressure", timeidx=time)
		hd = interplevel(h,p,lev)
		cart_proj = get_cartopy(hd)
		lats, lons = latlon_coords(hd)
		fig = pyplot.figure(figsize=(10,7.5))
		geo_axes = pyplot.axes(projection=cart_proj)
		states = NaturalEarthFeature(category = 'cultural',
				     scale = '50m', 
				     facecolor = 'none',
				     name = 'admin_1_states_provinces_shp')
		geo_axes.add_feature(states,linewidth=.5)
		geo_axes.coastlines('50m',linewidth=0.8)

		# Make the contour lines and fill them.
		
		pyplot.contourf(to_np(lons), to_np(lats), 
                               to_np(hd), 15, 
			       cmap=get_cmap('jet'),
                               transform=crs.PlateCarree())                              
		pyplot.colorbar(ax=geo_axes, shrink=.86)
		name = "Diabatic Heating" + str(lev) + str(time) + ".png"
		pyplot.savefig(name)		
		time = time + 1


