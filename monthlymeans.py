#This program loads all of the data, splits it by pattern type and then prints pics for each

#Import a bunch of stuff
import numpy as np
import matplotlib as plt
from netCDF4 import Dataset 
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import os
import datetime
import pickle
import csv
import scipy as sp
import scipy.ndimage as ndimage
import scipy.io as sio
from sklearn.metrics import mean_squared_error
from math import sqrt

#Check to see if out directory exists and if not, make it
outdir = "Allyear/pix"
if not os.path.exists(outdir):
    os.makedirs(outdir)

missing = 1 * 10**15

#Load in the Matlab File
mat_contents = sio.loadmat('C:/Users/CoeFamily/Documents/MATLAB/monthlydata.mat')


#Put Variables into an array (change the 3 letter month to the month you need)
h500 = mat_contents['octh500']
mslp = mat_contents['octmslp']
u850 = mat_contents['octu']
v850 = mat_contents['octv']

#Place data into a new array (need to subtract 1 from the CI due to index differences)
dset1 = h500[:,:,:]

#open the lat lon files    
latlon1 = Dataset("G:/Merra_new/19790101.nc4", mode='r')
    
#put the lat and lon data into arrays based on which grid type it is    
lats = np.squeeze(latlon1.variables['YDim'][:])
lons = np.squeeze(latlon1.variables['XDim'][:])

# mesh the grids together for plotting purposes
LON,LAT = np.meshgrid(lons,lats)

#close the lat lon files    
latlon1.close()

i = 0
year = 1979
#Plot the images for each of the clusters  
while i <= 31:    
    m = Basemap(projection='merc',
        resolution='h',area_thresh = 1000.0,lat_0=0,lon_0=100.,     	llcrnrlon=-90,
        llcrnrlat=30,urcrnrlon=-60,urcrnrlat=50)
    m.drawcoastlines()
    m.drawcountries()
    m.drawstates()
    
    #Average the data based on the CI Value
    H = np.zeros((lats.shape[0], lons.shape[0]))
    S = np.zeros((lats.shape[0], lons.shape[0]))
    U = np.zeros((lats.shape[0], lons.shape[0]))
    V = np.zeros((lats.shape[0], lons.shape[0]))
    count  = 0
    for q in range(0,31):
	for mm in range (lats.shape[0]):
		for nn in range( lons.shape[0]):
			H[mm][nn] = H[mm][nn]+ h500[q +(i*31)][mm][nn]
			S[mm][nn] = S[mm][nn] + mslp[q+(i*31)][mm][nn]
			U[mm][nn] = U[mm][nn] + u850[q+(i*31)][mm][nn]
			V[mm][nn] = V[mm][nn] + v850[q+(i*31)][mm][nn]
        count = count + 1

    for mm in range(lats.shape[0]):
	for nn in range(lons.shape[0]):
		H[mm][nn] = H[mm][nn] / (count*10)
		S[mm][nn] = S[mm][nn] / (count*100)
		U[mm][nn] = U[mm][nn] / count
		V[mm][nn] = V[mm][nn] / count    

    x,y = m(LON,LAT)  
    Z_500 = ndimage.gaussian_filter(S[:][:], sigma=0, order=0)        
    cd = m.contour(x,y,H[:][:],colors='blue',)
    plt.clabel(cd, inline=True, fmt='%1.0f', fontsize=10, colors='black')
    m.quiver(x[:,:],y[:,:],U[:][:],V[:][:],scale=None, scale_units='inches')
    cs = m.contour(x,y,Z_500, colors='black',)
    plt.clabel(cs, inline=True, fmt='%1.0f', fontsize=12, colors='k')
    m.drawparallels(np.arange(-80.,81.,10.), labels=[1,0,0,0],fontsize=10)
    m.drawmeridians(np.arange(-180.,181.,10.), labels=[0,0,0,1], fontsize = 10)
    
    
   
    #Change the title here tomatch the month we are doing this for
    title = 'October Mean'
    
    title1 = outdir + '/' + 'Oct' + str(year)+'.png'
    
    plt.title(title)
    
    plt.savefig(title1, bbox_inches='tight')
    
    plt.gcf().clear()

    i = i + 1
    year = year + 1

