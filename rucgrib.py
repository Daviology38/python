#! /usr/local/anaconda2/bin/python
import pygrib
import math
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import shiftgrid
import numpy as np
 
plt.figure(figsize=(12,8))
 
grib = 'ruc2_252_20060206_1900_000.grb' # Set the file name of your input GRIB file
grbs = pygrib.open(grib)
 
grb = grbs.select(name = '10 metre U wind component')[0]
data = grb.values
datatmp = grb.values

grb2 = grbs.select(name ='10 metre V wind component')[0]
data2 = grb.values
data2tmp = grb.values
 

lat, lon = grb.latlons()

m = Basemap(projection='lcc', lon_0=-72, lat_0=45, llcrnrlon = -80, llcrnrlat = 40, urcrnrlat = 50, urcrnrlon = -65)

nxv = 41; nyv = 41

udat, vdat, xv, yv = m.rotate_vector( data[36:,:], data2[36:,:], lon[36:,:], lat[36:,:], returnxy=True)

x, y = m(lon,lat)
 

i = 0
j = 0
rotcon = 0.422618
lonxx = -95.0

while j < 301:
	while i < 225:
		angle2 = rotcon * (lon[i][j] - lonxx) * 0.017453
		sinx2 = math.sin(math.radians(angle2))
		cosx2 = math.cos(math.radians(angle2))
		data[i][j] = cosx2*datatmp[i][j] + sinx2*data2tmp[i][j]
		data2[i][j] = -1*sinx2*datatmp[i][j] + cosx2*data2tmp[i][j]
		i = i + 1
	j = j + 1
	i = 0


cs = m.quiver(x[::3,::3], y[::3,::3], data[::3,::3], data2[::3,::3],units='width')
 
m.drawcoastlines()
m.drawmapboundary()
m.drawparallels(np.arange(-90.,120.,30.),labels=[1,0,0,0])
m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1])
 

plt.title('CAMS AOD forecast') # Set the name of the variable to plot
plt.savefig(grib+'.png') # Set the output file name
