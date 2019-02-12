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
ncfile = Dataset("WRFV3/run/wrfout_d02_2019-02-11_00:00:00")

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

	#2 meter temperature
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

	#surface pressure
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
		geo_axes.add_feature(states,linewidth=1.,edgecolor="black")
		geo_axes.coastlines('50m',linewidth=0.8)
		geo_axes.gridlines(color="black", linestyle="dotted")
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
		titlename = "Sea Level Pressure at " + str(time * 3) + "z"
		pyplot.title(titlename)
		pyplot.savefig(name)  
		time = time + 1

	#H500 with filled contours
	if(var.lower() == "h"):
		h = getvar(ncfile, "z", timeidx=time, units="dm")
		p = getvar(ncfile,"pressure", timeidx=time)
		ht = interplevel(h,p,500)
		cart_proj = get_cartopy(ht)
		lats, lons = latlon_coords(ht)
		fig = pyplot.figure(figsize=(10,7.5))
		geo_axes = pyplot.axes(projection=cart_proj)
		states = NaturalEarthFeature(category = 'cultural',
				     scale = '50m', 
				     facecolor = 'none',
				     name = 'admin_1_states_provinces_shp')
		geo_axes.add_feature(states,linewidth=1.,edgecolor="black")
		geo_axes.coastlines('50m',linewidth=0.8)
		geo_axes.gridlines(color="black", linestyle="dotted")
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
		titlename = "500 hPa Heights at " + str(time * 3) + "z"
		pyplot.title(titlename)
		pyplot.savefig(name)		
		time = time + 1

	#850 heights and winds
	if(var.lower() == "uv8h"):
		u = getvar(ncfile, "ua", timeidx=time, units="kt")
		v = getvar(ncfile, "va", timeidx=time, units="kt")
		wspd = getvar(ncfile, "wspd_wdir", timeidx=time, units = "kt")[0,:]
		z = getvar(ncfile, "z", timeidx=time, units="dm")
		p = getvar(ncfile,"pressure",timeidx=time)
		ht = interplevel(z,p,850)
		ulev = interplevel(u, p, 850)
		vlev = interplevel(v, p, 850)
		wspdlev = interplevel(wspd, p, 850)
		cart_proj = get_cartopy(ht)
		lats, lons = latlon_coords(ht)
		fig = pyplot.figure(figsize=(10,7.5))
		geo_axes = pyplot.axes(projection=cart_proj)
		states = NaturalEarthFeature(category = 'cultural',
				     scale = '50m', 
				     facecolor = 'none',
				     name = 'admin_1_states_provinces_shp')
		geo_axes.add_feature(states,linewidth=1.,edgecolor="black")
		geo_axes.coastlines('50m',linewidth=0.8)
		geo_axes.gridlines(color="black", linestyle="dotted")
		
		# Make the contour lines and fill them.
		contours = pyplot.contour(to_np(lons), to_np(lats), 
               		       to_np(ht), 15, colors="black",
                               transform=crs.PlateCarree())
		pyplot.clabel(contours, inline=1, fontsize=10, fmt="%i") 
		contours = pyplot.contour(to_np(lons), to_np(lats), 
               		       to_np(ht), 15, colors="black",
                               transform=crs.PlateCarree())
		pyplot.clabel(contours, inline=1, fontsize=10, fmt="%i")

		
		thin = [int(x/10.) for x in lons.shape]
		pyplot.barbs(to_np(lons[::thin[0], ::thin[1]]), 
          		to_np(lats[::thin[0], ::thin[1]]), 
          		to_np(ulev[::thin[0], ::thin[1]]),
          		to_np(vlev[::thin[0], ::thin[1]]), 
          		length=6,
          		transform=crs.PlateCarree())                      
		name = "Winds_H" + str(lev) + str(time) + ".png"
		titlename = "Winds and Heights at 850 hPa at " + str(time * 3) + "z"
		pyplot.title(titlename)
		pyplot.savefig(name)		
		time = time + 1
		pyplot.close()

	#850 heights, temperatures and winds
	if(var.lower() == "uv8"):
		u = getvar(ncfile, "ua", timeidx=time, units="kt")
		v = getvar(ncfile, "va", timeidx=time, units="kt")
		wspd = getvar(ncfile, "wspd_wdir", timeidx=time, units = "kt")[0,:]
		z = getvar(ncfile, "z", timeidx=time, units="dm")
		p = getvar(ncfile,"pressure",timeidx=time)
		temp = getvar(ncfile, "temp", timeidx=time, units = "degF")
		ht = interplevel(z,p,850)
		ulev = interplevel(u, p, 850)
		vlev = interplevel(v, p, 850)
		wspdlev = interplevel(wspd, p, 850)
		tlev = interplevel(temp, p, 850)
		cart_proj = get_cartopy(ht)
		lats, lons = latlon_coords(ht)
		fig = pyplot.figure(figsize=(10,7.5))
		geo_axes = pyplot.axes(projection=cart_proj)
		states = NaturalEarthFeature(category = 'cultural',
				     scale = '50m', 
				     facecolor = 'none',
				     name = 'admin_1_states_provinces_shp')
		geo_axes.add_feature(states,linewidth=1.,edgecolor="black")
		geo_axes.coastlines('50m',linewidth=0.8)
		geo_axes.gridlines(color="black", linestyle="dotted")
		tlevs = np.arange(int(tlev.min()), int(tlev.max()), 5.)
		# Make the contour lines and fill them.
		contours = pyplot.contour(to_np(lons), to_np(lats), 
               		       to_np(ht), 15, colors="black",
                               transform=crs.PlateCarree())
		pyplot.clabel(contours, inline=1, fontsize=10, fmt="%i")
		pyplot.contourf(to_np(lons), to_np(lats),
		       to_np(tlev), levels = tlevs,
		       cmap="jet",
		       extend="both",
		       transform = crs.PlateCarree())
		pyplot.colorbar(ax=geo_axes,shrink=.86)

		
		thin = [int(x/10.) for x in lons.shape]
		pyplot.barbs(to_np(lons[::thin[0], ::thin[1]]), 
          		to_np(lats[::thin[0], ::thin[1]]), 
          		to_np(ulev[::thin[0], ::thin[1]]),
          		to_np(vlev[::thin[0], ::thin[1]]), 
          		length=6,
          		transform=crs.PlateCarree())                      
		name = "Winds_T_H" + "850" + str(time) + ".png"
		titlename = "Winds, Temperature, and Heights at 850 hPa at " + str(time * 3) + "z"
		pyplot.title(titlename)
		pyplot.savefig(name)		
		time = time + 1
		pyplot.close()

	#700 heights and winds
	if(var.lower() == "uv7h"):
		u = getvar(ncfile, "ua", timeidx=time, units="kt")
		v = getvar(ncfile, "va", timeidx=time, units="kt")
		wspd = getvar(ncfile, "wspd_wdir", timeidx=time, units = "kt")[0,:]
		z = getvar(ncfile, "z", timeidx=time, units="dm")
		p = getvar(ncfile,"pressure",timeidx=time)
		ht = interplevel(z,p,700)
		ulev = interplevel(u, p, 700)
		vlev = interplevel(v, p, 700)
		wspdlev = interplevel(wspd, p, 700)
		cart_proj = get_cartopy(ht)
		lats, lons = latlon_coords(ht)
		fig = pyplot.figure(figsize=(10,7.5))
		geo_axes = pyplot.axes(projection=cart_proj)
		states = NaturalEarthFeature(category = 'cultural',
				     scale = '50m', 
				     facecolor = 'none',
				     name = 'admin_1_states_provinces_shp')
		geo_axes.add_feature(states,linewidth=1.,edgecolor="black")
		geo_axes.coastlines('50m',linewidth=0.8)
		geo_axes.gridlines(color="black", linestyle="dotted")
		# Make the contour lines and fill them.
		contours = pyplot.contour(to_np(lons), to_np(lats), 
               		       to_np(ht), 15, colors="black",
                               transform=crs.PlateCarree())
		pyplot.clabel(contours, inline=1, fontsize=10, fmt="%i")
		levels = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 					     90, 100, 110]    
		contours = pyplot.contour(to_np(lons), to_np(lats), 
               		       to_np(ht), 15, colors="black",
                               transform=crs.PlateCarree())
		pyplot.clabel(contours, inline=1, fontsize=10, fmt="%i")

		
		thin = [int(x/10.) for x in lons.shape]
		pyplot.barbs(to_np(lons[::thin[0], ::thin[1]]), 
          		to_np(lats[::thin[0], ::thin[1]]), 
          		to_np(ulev[::thin[0], ::thin[1]]),
          		to_np(vlev[::thin[0], ::thin[1]]), 
          		length=6,
          		transform=crs.PlateCarree())                      
		name = "Winds_H" + "700" + str(time) + ".png"
		titlename = "Winds and Heights at 700 hPa at " + str(time * 3) + "z"
		pyplot.title(titlename)
		pyplot.savefig(name)		
		time = time + 1
		pyplot.close()

	#700 heights, temperatures and winds
	if(var.lower() == "uv7"):
		u = getvar(ncfile, "ua", timeidx=time, units="kt")
		v = getvar(ncfile, "va", timeidx=time, units="kt")
		wspd = getvar(ncfile, "wspd_wdir", timeidx=time, units = "kt")[0,:]
		z = getvar(ncfile, "z", timeidx=time, units="dm")
		p = getvar(ncfile,"pressure",timeidx=time)
		temp = getvar(ncfile, "temp", timeidx=time, units = "degF")
		ht = interplevel(z,p,700)
		ulev = interplevel(u, p, 700)
		vlev = interplevel(v, p, 700)
		wspdlev = interplevel(wspd, p, 700)
		tlev = interplevel(temp, p, 700)
		cart_proj = get_cartopy(ht)
		lats, lons = latlon_coords(ht)
		fig = pyplot.figure(figsize=(10,7.5))
		geo_axes = pyplot.axes(projection=cart_proj)
		states = NaturalEarthFeature(category = 'cultural',
				     scale = '50m', 
				     facecolor = 'none',
				     name = 'admin_1_states_provinces_shp')
		geo_axes.add_feature(states,linewidth=1.,edgecolor="black")
		geo_axes.coastlines('50m',linewidth=0.8)
		geo_axes.gridlines(color="black", linestyle="dotted")
		tlevs = np.arange(int(tlev.min()), int(tlev.max()), 5.)
		# Make the contour lines and fill them.
		contours = pyplot.contour(to_np(lons), to_np(lats), 
               		       to_np(ht), 15, colors="black",
                               transform=crs.PlateCarree())
		pyplot.clabel(contours, inline=1, fontsize=10, fmt="%i")
		pyplot.contourf(to_np(lons), to_np(lats),
		       to_np(tlev), levels = tlevs,
		       cmap="jet",
		       extend="both",
		       transform = crs.PlateCarree())
		pyplot.colorbar(ax=geo_axes,shrink=.86)

		
		thin = [int(x/10.) for x in lons.shape]
		pyplot.barbs(to_np(lons[::thin[0], ::thin[1]]), 
          		to_np(lats[::thin[0], ::thin[1]]), 
          		to_np(ulev[::thin[0], ::thin[1]]),
          		to_np(vlev[::thin[0], ::thin[1]]), 
          		length=6,
          		transform=crs.PlateCarree())                      
		name = "Winds_T_H" + "700" + str(time) + ".png"
		titlename = "Winds, Temperature, and Heights at 700 hPa at " + str(time * 3) + "z"
		pyplot.title(titlename)
		pyplot.savefig(name)		
		time = time + 1
		pyplot.close()

#700 heights, RH, and winds
	if(var.lower() == "uv7rh"):
		u = getvar(ncfile, "ua", timeidx=time, units="kt")
		v = getvar(ncfile, "va", timeidx=time, units="kt")
		wspd = getvar(ncfile, "wspd_wdir", timeidx=time, units = "kt")[0,:]
		z = getvar(ncfile, "z", timeidx=time, units="dm")
		p = getvar(ncfile,"pressure",timeidx=time)
		rh = getvar(ncfile, "rh", timeidx=time)
		ht = interplevel(z,p,700)
		ulev = interplevel(u, p, 700)
		vlev = interplevel(v, p, 700)
		wspdlev = interplevel(wspd, p, 700)
		rhlev = interplevel(rh, p, 700)
		cart_proj = get_cartopy(ht)
		lats, lons = latlon_coords(ht)
		fig = pyplot.figure(figsize=(10,7.5))
		geo_axes = pyplot.axes(projection=cart_proj)
		states = NaturalEarthFeature(category = 'cultural',
				     scale = '50m', 
				     facecolor = 'none',
				     name = 'admin_1_states_provinces_shp')
		geo_axes.add_feature(states,linewidth=1.,edgecolor="black")
		geo_axes.coastlines('50m',linewidth=0.8)
		geo_axes.gridlines(color="black", linestyle="dotted")
		
		# Make the contour lines and fill them.
		contours = pyplot.contour(to_np(lons), to_np(lats), 
               		       to_np(ht), 15, colors="black",
                               transform=crs.PlateCarree())
		pyplot.clabel(contours, inline=1, fontsize=10, fmt="%i")

		levels = [50,55,60,65,70,75,80,85,90,95,100]
		pyplot.contourf(to_np(lons), to_np(lats),
		       to_np(rhlev), levels = levels,
		       cmap="Oranges",
		       extend="both",
		       transform = crs.PlateCarree())
		pyplot.colorbar(ax=geo_axes,shrink=.86)

		
		thin = [int(x/10.) for x in lons.shape]
		pyplot.barbs(to_np(lons[::thin[0], ::thin[1]]), 
          		to_np(lats[::thin[0], ::thin[1]]), 
          		to_np(ulev[::thin[0], ::thin[1]]),
          		to_np(vlev[::thin[0], ::thin[1]]), 
          		length=6,
          		transform=crs.PlateCarree())                      
		name = "Winds_RH_H" + "700" + str(time) + ".png"
		titlename = "Winds, RH, and Heights at 700 hPa at " + str(time * 3) + "z"
		pyplot.title(titlename)
		pyplot.savefig(name)		
		time = time + 1
		pyplot.close()


	#10 Meter Wind with SLP
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
		geo_axes.gridlines(color="black", linestyle="dotted")
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
		titlename = "10 meter wind and Sea Level Pressure at " + str(time * 3) + "z"
		pyplot.title(titlename)
		pyplot.savefig(name)  
		time = time + 1
	
	#250 hPa Jet Stream
	if(var.lower() == "uv250"):
		u = getvar(ncfile, "ua", timeidx=time, units="kt")
		v = getvar(ncfile, "va", timeidx=time, units="kt")
		wspd = getvar(ncfile, "wspd_wdir", timeidx=time, units = "kt")[0,:]
		z = getvar(ncfile, "z", timeidx=time, units="dm")
		p = getvar(ncfile,"pressure",timeidx=time)
		ht = interplevel(z,p,250)
		ulev = interplevel(u, p, 250)
		vlev = interplevel(v, p, 250)
		wspdlev = interplevel(wspd, p, 250)
		ht = interplevel(z, p, 250)
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
		geo_axes.gridlines(color="black", linestyle="dotted")

		# Make the contour lines and fill them.
		contours = pyplot.contour(to_np(lons), to_np(lats), 
               		       to_np(ht), 15, colors="black",
                               transform=crs.PlateCarree())
		pyplot.clabel(contours, inline=1, fontsize=10, fmt="%i")
		levels = [60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230]    
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
		name = "Winds250" + str(lev) + str(time) + ".png"
		titlename = "Winds and Heights at 250 hPa at " + str(time * 3) + "z"
		pyplot.title(titlename)
		pyplot.savefig(name)		
		time = time + 1
		pyplot.close()

	#SLP and DBZ
	if(var.lower() == "dbz"):
		slp = getvar(ncfile, "slp", timeidx=time)
		dbz = getvar(ncfile, "dbz", timeidx=time)
		p = getvar(ncfile, "pressure", timeidx=time)
		dbzlev = interplevel(dbz,p,1000)
		print dbzlev
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
		geo_axes.gridlines(color="black", linestyle="dotted")
		# Set the contour levels so that all plots match
		#levels = np.arange(960.,1040.,2.5)
		levels = [5 + 5*n for n in range(15)]
		# Make the contour lines and fill them.
		contours = pyplot.contour(to_np(lons), to_np(lats), 
               		       to_np(smooth_slp), 15, colors="black",
                               transform=crs.PlateCarree())
		pyplot.clabel(contours, inline=1, fontsize=10, fmt="%i")
		dbz_contours = pyplot.contourf(to_np(lons), to_np(lats), to_np(dbzlev), levels = levels, cmap = get_cmap("jet"), transform=crs.PlateCarree())
		pyplot.colorbar(dbz_contours, ax=geo_axes)
		name = "DBZ" + str(time) + ".png"
		titlename = "Reflectivity and Sea Level Pressure at " + str(time * 3) + "z"
		pyplot.title(titlename)
		pyplot.savefig(name)  
		time = time + 1
	

	
