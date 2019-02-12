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
outdir = "C:/Users/CoeFamily/Documents/MATLAB/SON_new/"
if not os.path.exists(outdir):
    os.makedirs(outdir)

missing = 1 * 10**15

#Load in the Matlab File
mat_contents = sio.loadmat('C:/Users/CoeFamily/Documents/MATLAB/SON_new/SONdailymean_95/CI_results.mat')

#Put Variables into an array
K = mat_contents['K']
ci = 9

#Ask user which CI value to make pics for
#ci = int(raw_input("Enter CI value (Between 1-25)"))

#Place data into a new array (need to subtract 1 from the CI due to index differences)
dset1 = K[:,ci-1]

#open the lat lon files    
latlon1 = Dataset("G:/Merra_new/19790101.nc4", mode='r')
    
#put the lat and lon data into arrays based on which grid type it is    
lats = np.squeeze(latlon1.variables['YDim'][:])
lons = np.squeeze(latlon1.variables['XDim'][:])

# mesh the grids together for plotting purposes
LON,LAT = np.meshgrid(lons,lats)

#close the lat lon files    
latlon1.close()

#preallocate all of our multi-dimensional lists
h500 = np.zeros((dset1.shape[0],lats.shape[0],lons.shape[0]))
SLP = np.zeros((dset1.shape[0],lats.shape[0],lons.shape[0]))
u = np.zeros((dset1.shape[0],lats.shape[0],lons.shape[0]))
v = np.zeros((dset1.shape[0],lats.shape[0],lons.shape[0]))
h500t = np.zeros((lats.shape[0],lons.shape[0]))
SLPt = np.zeros((lats.shape[0],lons.shape[0])) 
ut = np.zeros((lats.shape[0],lons.shape[0]))
vt = np.zeros((lats.shape[0],lons.shape[0]))

#enter the start and end date for the loop to put all of the data into an array where ( row equals the date, columns equal the lats and pages equal the lons)
tnum = '1979'
enum = '2011'
m = 1
d = 1
year = 1979
p = 0
  
while tnum != enum:
    if m < 10:
        strm = '0' + str(m)
    else:
        strm = str(m)     
    if d < 10:
        strd = '0' + str(d)
    else:
        strd = str(d)
    
    stry = str(year)
    
    str1 = stry + strm + strd
    if(strm == '09' or strm == '10' or strm == '11'):
        filename = 'G:/Merra_new/'+ str1 +'.nc4'
        filename1 = 'G:/Merra_new/' + str1 + '.nc4'

            
        tmp = Dataset(filename1, mode='r')   
        h500m = tmp.variables['H'][4,4,:,:]                          
        SLPm = tmp.variables['SLP'][4,:,:]
        um = tmp.variables['U'][4,2,:,:]
        vm = tmp.variables['V'][4,2,:,:]
	
        for mm in range(lats.shape[0]):
            for nn in range(lons.shape[0]):
                h500[p][mm][nn] = h500m[mm][nn]       
                SLP[p][mm][nn] = SLPm[mm][nn]
                u[p][mm][nn] = um[mm][nn]
                v[p][mm][nn] = vm[mm][nn]
                if(u[p][mm][nn] == missing):
                    u[p][mm][nn] = 0
                if(v[p][mm][nn] == missing):
                    v[p][mm][nn] = 0    
                if(SLP[p][mm][nn] == missing):
                    SLP[p][mm][nn] = 0
                if(h500[p][mm][nn] == missing):
                    h500[p][mm][nn] = 0
        p = p + 1
    #modulo operator to check for a leap year and add the extra date if needed        
    modulo = year % 4        
    
    if(m == 1 and d < 31):
        d = d + 1
    elif( m == 1 and d == 31):
        m = m + 1
        d = 1
    elif( m == 2 and d < 28):
        d = d + 1
    elif( m == 2 and d == 28 and modulo == 0):
        d = d + 1
    elif( m == 2 and d == 28 and modulo != 0):
        m = m + 1
        d = 1
    elif( m == 2 and d > 28):
        m = 3
        d = 1
    elif( m == 3 and d < 31):
        d = d + 1
    elif( m == 3 and d == 31):
        m = m + 1
        d = 1
    elif( m == 4 and d < 30):
        d = d + 1
    elif( m == 4 and d == 30):
        m = m + 1
        d = 1
    elif( m == 5 and d < 31):
        d = d + 1
    elif( m == 5 and d == 31):
        m = m + 1 
        d = 1
    elif( m == 6 and d < 30):
        d = d + 1
    elif( m == 6 and d == 30):
        d = 1
        m = m + 1
    elif( m == 7 and d < 31):
        d = d + 1
    elif( m == 7 and d == 31):
        d = 1
        m = m + 1
    elif( m == 8 and d < 31):
        d = d + 1
    elif( m == 8 and d == 31):
        d = 1 
        m = m + 1
    elif( m == 9 and d < 30):
        d = d + 1
    elif( m == 9 and d == 30):
        d = 1
        m = m + 1
    elif( m == 10 and d < 31):
        d = d + 1
    elif( m == 10 and d == 31):
        d = 1 
        m = m + 1
    elif( m == 11 and d < 30):
        d = d + 1
    elif( m == 11 and d == 30):
        d = 1
        m = m + 1
    elif( m == 12 and d < 31):
        d = d + 1
    else:
        m = 1
        d = 1
        year = year + 1
        tnum = str(year)
        
i = 1 
#Plot the images for each of the clusters  
while i <= ci:    
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
    for q in range(dset1.shape[0]):
	if(dset1[q] == i):
		for mm in range (lats.shape[0]):
			for nn in range( lons.shape[0]):
				H[mm][nn] = H[mm][nn]+ h500[q][mm][nn]
				S[mm][nn] = S[mm][nn] + SLP[q][mm][nn]
				U[mm][nn] = U[mm][nn] + u[q][mm][nn]
				V[mm][nn] = V[mm][nn] + v[q][mm][nn]
                count = count + 1

    for mm in range(lats.shape[0]):
	for nn in range(lons.shape[0]):
		H[mm][nn] = H[mm][nn] / (count*10)
		S[mm][nn] = S[mm][nn] / (count*100)
		U[mm][nn] = U[mm][nn] / count
		V[mm][nn] = V[mm][nn] / count    

    x,y = m(LON,LAT)  
    Z_500 = ndimage.gaussian_filter(S[:][:], sigma=5, order=0)        
    cd = m.contour(x,y,H[:][:],colors='blue',)
    plt.clabel(cd, inline=True, fmt='%1.0f', fontsize=10, colors='black')
    m.quiver(x[:,:],y[:,:],U[:][:],V[:][:],scale=None, scale_units='inches')
    cs = m.contour(x,y,Z_500, colors='black',)
    plt.clabel(cs, inline=True, fmt='%1.0f', fontsize=12, colors='k')
    m.drawparallels(np.arange(-80.,81.,10.), labels=[1,0,0,0],fontsize=10)
    m.drawmeridians(np.arange(-180.,181.,10.), labels=[0,0,0,1], fontsize = 10)
    
    
   
    
    title = 'Cluster ' + str(i)
    
    title1 = outdir + '/_' + str(i)
    
    plt.title(title)
    
    plt.savefig(title1, bbox_inches='tight')
    
    plt.gcf().clear()

    i = i + 1





