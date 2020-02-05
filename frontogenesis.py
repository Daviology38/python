import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import metpy.calc as mpcalc
from metpy.units import units
import numpy as np
from netCDF4 import Dataset
from wrf import getvar, to_np, get_cartopy, latlon_coords, cartopy_xlim, cartopy_ylim, interplevel, smooth2d

ncfile = Dataset("/media/mariofire/4TBExternal/newrf/wrf_files/reg/wrfout_d01_2006-10-26")

level = 850 * units.hPa
h = getvar(ncfile, "z", timeidx=0, units = 'dm')
u = getvar(ncfile, "ua", timeidx=0)
v = getvar(ncfile, "va", timeidx=0)
temp = getvar(ncfile, "tc", timeidx=0)
p = getvar(ncfile,"pressure", timeidx=0)
hght_850 = interplevel(h,p,850)
tmpk_850 = interplevel(temp,p,850)
uwnd_850 = interplevel(u,p,850)
vwnd_850 = interplevel(v,p,850)
lats, lons = latlon_coords(hght_850)
uwnd_850 = uwnd_850[:,:]
vwnd_850 = vwnd_850[:,:]
# Calculate potential temperature for frontogenesis calculation
thta_850 = mpcalc.potential_temperature(level, tmpk_850)

dx, dy = mpcalc.lat_lon_grid_deltas(lons, lats)

fronto_850 = mpcalc.frontogenesis(thta_850, uwnd_850.values*units('m/s'), vwnd_850.values*units('m/s'), dx, dy, dim_order='yx')

# A conversion factor to get frontogensis units of K per 100 km per 3 h
convert_to_per_100km_3h = 1000*100*3600*3

# Set map projection
mapcrs = ccrs.LambertConformal(central_longitude=-100, central_latitude=35,
                               standard_parallels=(30, 60))

# Set projection of the data (GFS is lat/lon)
datacrs = ccrs.PlateCarree()

# Start figure and limit the graphical area extent
fig = plt.figure(1, figsize=(14, 12))
ax = plt.subplot(111, projection=mapcrs)
ax.set_extent([-90, -60, 30, 50], ccrs.PlateCarree())

# Add map features of Coastlines and States
ax.add_feature(cfeature.COASTLINE.with_scale('50m'))
ax.add_feature(cfeature.STATES.with_scale('50m'))

# Plot 850-hPa Frontogenesis
clevs_tmpc = np.arange(-40, 41, 2)
cf = ax.contourf(lons, lats, fronto_850*convert_to_per_100km_3h, np.arange(-8, 8.5, 0.5),
                 cmap=plt.cm.bwr, extend='both', transform=datacrs)
cb = plt.colorbar(cf, orientation='horizontal', pad=0, aspect=50, extendrect=True)
cb.set_label('Frontogenesis K / 100 km / 3 h')

# Plot 850-hPa Temperature in Celsius
csf = ax.contour(lons, lats, tmpk_850, clevs_tmpc, colors='grey',
                 linestyles='dashed', transform=datacrs)
plt.clabel(csf, fmt='%d')

# Plot 850-hPa Geopotential Heights
clevs_850_hght = np.arange(0, 8000, 30)
cs = ax.contour(lons, lats, hght_850, clevs_850_hght, colors='black', transform=datacrs)
plt.clabel(cs, fmt='%d')

# Plot 850-hPa Wind Barbs only plotting every fifth barb
#wind_slice = (slice(None, None, 5), slice(None, None, 5))
#ax.barbs(lons[wind_slice[0]], lats[wind_slice[1]],
         #uwnd_850[wind_slice] * units('kt'), vwnd_850[wind_slice] * units('kt'),
         #color='black', transform=datacrs)

# Plot some titles
plt.title('GFS 850-hPa Geopotential Heights (m), Temp (C), and Winds', loc='left')
plt.title('Valid Time:', loc='right')

plt.show()


