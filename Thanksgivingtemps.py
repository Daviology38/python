from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset 
import os


#Check to see if out directory exists and if not, make it
outdir = "Thanksgivingtemps/"
if not os.path.exists(outdir):
    os.makedirs(outdir)


#open the lat lon files
    
latlon1 = Dataset("G:/MERRA/2meter_temp/01011980.SUB.nc4", mode='r')

#put the lat and lon data into arrays based on which grid type it is
    
lats = np.squeeze(latlon1.variables['lat'][:])
lons = np.squeeze(latlon1.variables['lon'][:])

# mesh the grids together for plotting purposes
    
LON,LAT = np.meshgrid(lons,lats)
    
#close the lat lon files
    
latlon1.close()

#Load in the dates of Thanksgiving for the past 30+ years
with open('C:/Users/CoeFamily/Documents/MATLAB/thanksgivingdatessince1980.txt') as f:
    dates = f.readlines()

# you may also want to remove whitespace characters like `\n` at the end of each line
dates = [x.strip() for x in dates] 
year = 1980
Toverall = np.zeros((lats.shape[0],lons.shape[0]))
for d in range(len(dates)):
    
    #open the datafile for the date
    filename = "G:/MERRA/2meter_temp/" + str(dates[d]) + ".SUB.nc4"
    tmp = Dataset(filename, mode="r")
    T = np.zeros((lats.shape[0],lons.shape[0]))
    count = 0
    for q in range(0,23):
      for mm in range (lats.shape[0]):
	   for nn in range( lons.shape[0]):
            if(tmp.variables['T2M'][q,mm,nn] > T[mm][nn]):
             T[mm][nn] =  tmp.variables['T2M'][q,mm,nn]
            
       
    
    #Average the temperature out and then convert from K to C        
    T = (((T) - 273.15) * (9/5)) + 32
    for mm in range(lats.shape[0]):
        for nn in range(lons.shape[0]):
            Toverall[mm][nn] = Toverall[mm][nn] + T[mm][nn]
    #map = Basemap(projection='merc', lat_0 = 43.1, lon_0 = -71.725,
    #        resolution = 'h', area_thresh = 0.1,
    #        llcrnrlon=-73.72, llcrnrlat=41.12,
    #        urcrnrlon=-69.73, urcrnrlat=45.08)
    
    
    #map.drawcoastlines()
    ## map.fillcontinents(color='coral',lake_color='aqua')
    #map.drawcountries()
    
    #map.drawstates()
    #map.drawmapboundary()
    #map.drawmapboundary(fill_color='#46bcec')
    #x,y = map(LON,LAT) 
    #clevs = np.arange(-10,60,2)  
    #cnplot = map.contourf(x,y,T[:][:],clevs,cmap=plt.cm.jet)

    #-- add colorbar
    #cbar = map.colorbar(cnplot,location='bottom',pad="10%")      #-- pad: distance between map and colorbar
    #cbar.set_label('deg F')
    
    #Change the title here tomatch the month we are doing this for
    #title = 'Thanksgiving Mean Temp'
    
    #title1 = outdir + '/' + 'Thanksgiving' + str(year)+'.png'
    
    #plt.title(title)
    
    #plt.savefig(title1, bbox_inches='tight')
    
    #plt.gcf().clear()
    
    year = year + 1
    
    

Toverall = Toverall / (len(dates))
map = Basemap(projection='merc', lat_0 = 43.1, lon_0 = -71.725,
            resolution = 'h', area_thresh = 0.1,
            llcrnrlon=-73.72, llcrnrlat=41.12,
            urcrnrlon=-69.73, urcrnrlat=45.08)
    
    
map.drawcoastlines()
# map.fillcontinents(color='coral',lake_color='aqua')
map.drawcountries()
    
map.drawstates()
map.drawmapboundary()
map.drawmapboundary(fill_color='#46bcec')
x,y = map(LON,LAT) 
clevs = np.arange(-10,60,2)  
cnplot = map.contourf(x,y,Toverall[:][:],clevs,cmap=plt.cm.jet)

#-- add colorbar
cbar = map.colorbar(cnplot,location='bottom',pad="10%")      #-- pad: distance between map and colorbar
cbar.set_label('deg F')

#Change the title here tomatch the month we are doing this for
title = 'Thanksgiving Mean Temp'
    
title1 = outdir + '/' + 'Thanksgivingoverall.png'
    
plt.title(title)
    
plt.savefig(title1, bbox_inches='tight')
    
plt.gcf().clear()