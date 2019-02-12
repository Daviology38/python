# -*- coding: utf-8 -*-
"""
Created on Sat Feb 09 11:46:48 2019

@author: CoeFamily
"""

import numpy as np
import scipy.io as sio
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import os
import scipy.ndimage as ndimage
import scipy.io as sio
import math
from sklearn.metrics import mean_squared_error



#Function to find the 2d correlation coefficient
def mean2(x):
    y = np.sum(x) / np.size(x)
    return y

def corr2(a,b):
    a = a - mean2(a)
    b = b - mean2(b)
    
    r = (a*b).sum() / math.sqrt((a*a).sum() * (b*b).sum())
    return r






#Make directory to save data files to if it doesn't already exist
outdir = "kmeans_daily/cluster/"
if not os.path.exists(outdir):
    os.makedirs(outdir)

#pick number of clusters to analyze
clustnum = 9

#load in the data
#Define our missing data
missing = 1 * 10**15

#open the overall data file
#Load in the Matlab Files
mat_contents = sio.loadmat('C:/Users/CoeFamily/Documents/MATLAB/monthlydata.mat')
#mat_contents2 = sio.loadmat('C:/Users/CoeFamily/Documents/MATLAB/CI_results_07.mat')
mat_contents2 = sio.loadmat('C:/Users/CoeFamily/Documents/MATLAB/SONdailymean_95/CI_results.mat')

#Put Variables into an array (change the 3 letter month to the month you need)
h500oct = mat_contents['octh500']
h500sep = mat_contents['seph500']
h500nov = mat_contents['novh500']
K = mat_contents2['K']
K = K[:,8]

#combine the h500 arrays
h500 = np.zeros((2912,h500oct.shape[1],h500oct.shape[2]))
i = 0
y = 0
countsep = 0
countoct = 0
countnov = 0

while i <=2911:
    if y <=29:
        h500[i][:][:] = h500sep[countsep][:][:]
        y = y + 1
        countsep = countsep + 1
    elif (y > 29 and y <=60):
        h500[i][:][:] = h500oct[countoct][:][:]
        y = y + 1
        countoct = countoct + 1
    else:
        h500[i][:][:] = h500nov[countnov][:][:]
        y = y + 1
        countnov = countnov + 1
        
    if y == 91:
        y = 0
        
    i = i + 1    
    
#Create the anomaly fields for each day
#Start by using the seasonal mean of H500

h500a = h500.mean(axis=0)   

#Next separate the data by WT

i = 0
wt1 = np.zeros((18,25))
wt2 = np.zeros((18,25))
wt3 = np.zeros((18,25))
wt4 = np.zeros((18,25))
wt5 = np.zeros((18,25))
wt6 = np.zeros((18,25))
wt7 = np.zeros((18,25))
wt8 = np.zeros((18,25))
wt9 = np.zeros((18,25))
count1 = 0
count2 = 0
count3 = 0
count4 = 0
count5 = 0
count6 = 0
count7 = 0
count8 = 0
count9 = 0

while i <= 2911:
    if(K[i] == 1):
        wt1[:][:] = wt1[:][:] + h500[i][:][:]
        count1 = count1 + 1
    elif(K[i] == 2):
        wt2[:][:] = wt2[:][:] + h500[i][:][:]
        count2 = count2 + 1
    elif(K[i] == 3):
        wt3[:][:] = wt3[:][:] + h500[i][:][:]
        count3 = count3 + 1
    elif(K[i] == 4):
        wt4[:][:] = wt4[:][:] + h500[i][:][:]
        count4 = count4 + 1
    elif(K[i] == 5):
        wt5[:][:] = wt5[:][:] + h500[i][:][:]
        count5 = count5 + 1
    elif(K[i] == 6):
        wt6[:][:] = wt6[:][:] + h500[i][:][:]
        count6 = count6 + 1
#    else:
#        wt7[:][:] = wt7[:][:] + h500[i][:][:]
#        count7 = count7 + 1
    elif(K[i] == 7):
        wt7[:][:] = wt7[:][:] + h500[i][:][:]
        count7 = count7 + 1
    elif(K[i] == 8):
        wt8[:][:] = wt8[:][:] + h500[i][:][:]
        count8 = count8 + 1
    else:
        wt9[:][:] = wt9[:][:] + h500[i][:][:]
        count9 = count9 + 1
    i = i + 1
    
wt1 = wt1 / count1
wt2 = wt2 / count2 
wt3 = wt3 / count3
wt4 = wt4 / count4
wt5 = wt5 / count5
wt6 = wt6 / count6
wt7 = wt7 / count7
wt8 = wt8 / count8
wt9 = wt9 / count9

#make the anomaly arrays
i = 0
h500aoverall = np.zeros((2912,18,25))
while i <= 2911:
    h500aoverall[i][:][:] = h500[i][:][:] - h500a[:][:]
    i = i + 1
wt1a = wt1[:][:] - h500a[:][:]
wt2a = wt2[:][:] - h500a[:][:]
wt3a = wt3[:][:] - h500a[:][:]
wt4a = wt4[:][:] - h500a[:][:]
wt5a = wt5[:][:] - h500a[:][:]
wt6a = wt6[:][:] - h500a[:][:]
wt7a = wt7[:][:] - h500a[:][:]  
wt8a = wt8[:][:] - h500a[:][:]
wt9a = wt9[:][:] - h500a[:][:]  

#Now compute the correlation, rmse and bias for each day of the SON season
corrs = np.zeros((2911,9))
rmse = np.zeros((2911,9))
bias = np.zeros((2911,9))

i = 0
count1 = 0
count2 = 0 
count3 = 0
count4 = 0
count5 = 0
count6 = 0
count7 = 0
count8 = 0
count9 = 0
while i <= 2911:
    if( K[i] == 1):
        corrs[count1][0] = corr2(wt1a,h500aoverall[i][:][:])
        rmse[count1][0] = math.sqrt(mean_squared_error(wt1a,h500aoverall[i][:][:]))
        bias[count1][0] = (wt1a[:][:] - h500aoverall[i][:][:]).mean()
        count1 = count1 + 1
    elif( K[i] == 2):
        corrs[count2][1] = corr2(wt2a, h500aoverall[i][:][:])
        rmse[count2][1] = math.sqrt(mean_squared_error(wt2a,h500aoverall[i][:][:]))
        bias[count2][1] = (wt2a[:][:] - h500aoverall[i][:][:]).mean()
        count2 = count2 + 1
    elif( K[i] == 3):
        corrs[count3][2] = corr2(wt3a, h500aoverall[i][:][:])
        rmse[count3][2] = math.sqrt(mean_squared_error(wt3a,h500aoverall[i][:][:]))
        bias[count3][2] = (wt3a[:][:] - h500aoverall[i][:][:]).mean()
        count3 = count3 + 1
    elif( K[i] == 4):
        corrs[count4][3] = corr2(wt4a, h500aoverall[i][:][:])
        rmse[count4][3] = math.sqrt(mean_squared_error(wt4a,h500aoverall[i][:][:]))
        bias[count4][3] = (wt4a[:][:] - h500aoverall[i][:][:]).mean()
        count4 = count4 + 1
    elif( K[i] == 5):
        corrs[count5][4] = corr2(wt5a, h500aoverall[i][:][:])
        rmse[count5][4] = math.sqrt(mean_squared_error(wt5a,h500aoverall[i][:][:]))
        bias[count5][4] = (wt5a[:][:] - h500aoverall[i][:][:]).mean()
        count5 = count5 + 1
    elif( K[i] == 6):
        corrs[count6][5] = corr2(wt6a, h500aoverall[i][:][:])
        rmse[count6][5] = math.sqrt(mean_squared_error(wt6a,h500aoverall[i][:][:]))
        bias[count6][5] = (wt6a[:][:] - h500aoverall[i][:][:]).mean()
        count6 = count6 + 1
#    else:
#        corrs[count7][6] = corr2(wt7a, h500aoverall[i][:][:])
#        rmse[count7][6] = math.sqrt(mean_squared_error(wt7a,h500aoverall[i][:][:]))
#        bias[count7][6] = (wt7a[:][:] - h500aoverall[i][:][:]).mean()
#        count7 = count7 + 1
    elif( K[i] == 7):
        corrs[count7][6] = corr2(wt7a, h500aoverall[i][:][:])
        rmse[count7][6] = math.sqrt(mean_squared_error(wt7a,h500aoverall[i][:][:]))
        bias[count7][6] = (wt7a[:][:] - h500aoverall[i][:][:]).mean()
        count7 = count7 + 1
    elif( K[i] == 8):
        corrs[count8][7] = corr2(wt8a, h500aoverall[i][:][:])
        rmse[count8][7] = math.sqrt(mean_squared_error(wt8a,h500aoverall[i][:][:]))
        bias[count8][7] = (wt8a[:][:] - h500aoverall[i][:][:]).mean()
        count8 = count8 + 1
    else:
        corrs[count9][8] = corr2(wt9a, h500aoverall[i][:][:])
        rmse[count9][8] = math.sqrt(mean_squared_error(wt9a,h500aoverall[i][:][:]))
        bias[count9][8] = (wt9a[:][:] - h500aoverall[i][:][:]).mean()
        count9 = count9 + 1

    i = i + 1

#create the histogram for correlation coefficient
binscorr = np.arange(-1.,1.,.1)
i = 0
fig = plt.figure(figsize=(16,10))
fig.subplots_adjust(hspace=0.4, wspace=0.4)
while i <= clustnum -1:
    tmp = corrs[:,i]
    tmp = tmp[tmp != 0]
    hist, bins = np.histogram(tmp, bins = binscorr)
    width = 0.7 * (bins[1]-bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    fig.add_subplot(3,3,i+1)    
    plt.bar(center, hist, align='center', width=width)
    plt.ylim(0,175,25)
    plt.xlim(-1,1,.1)
    if( i !=clustnum):
        name = "WT" + str(i+1)
    else:
        name = "WT" + str(i)
    plt.title(name)
    i = i + 1
fig.suptitle('Correlation Coefficients')
plt.savefig(outdir+'correlation.png', bbox_inches='tight')


#create the histogram for rmse
binsrmse = np.arange(0,250,20)
i = 0
fig = plt.figure(figsize=(16,10))
fig.subplots_adjust(hspace=0.4, wspace=0.4)
while i <= clustnum-1:
    tmp = rmse[:,i]
    tmp = tmp[tmp != 0]
    hist, bins = np.histogram(tmp, bins = binsrmse)
    width = 0.7 * (bins[1]-bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    fig.add_subplot(3,3,i+1)    
    plt.bar(center, hist, align='center', width=width)
    plt.ylim(0,250,25)
    plt.xlim(0,250,25)
    if( i !=clustnum):
        name = "WT" + str(i+1)
    else:
        name = "WT" + str(i)
    plt.title(name)
    i = i + 1
fig.suptitle('RMSE')
plt.savefig(outdir+'rmse.png', bbox_inches='tight')

#create the histogram for bias
binsbias = np.arange(-200,200,20)
i = 0
fig = plt.figure(figsize=(16,10))
fig.subplots_adjust(hspace=0.4, wspace=0.4)
while i <= clustnum-1:
    tmp = bias[:,i]
    tmp = tmp[tmp != 0]
    hist, bins = np.histogram(tmp, bins = binsbias)
    width = 0.7 * (bins[1]-bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    fig.add_subplot(3,3,i+1)    
    plt.bar(center, hist, align='center', width=width)
    plt.ylim(0,120,10)
    plt.xlim(-200,200,50)
    if( i !=clustnum):
        name = "WT" + str(i+1)
    else:
        name = "WT" + str(i)
    plt.title(name)
    i = i + 1
fig.suptitle('Bias')
plt.savefig(outdir+'bias.png', bbox_inches='tight')