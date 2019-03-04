import xarray as xr
import dask
from dask.diagnostics import ProgressBar
import pandas as pd
import numpy as np
import scipy
from scipy import signal
import matplotlib.pyplot as plt
from scipy.stats import skew
from scipy.stats import kurtosis
from scipy.stats import ttest_1samp
import cartopy.crs as ccrs
import datetime
import matplotlib.lines as mlines
import numpy.ma as ma
ocean=xr.open_dataset('thetao_.nc')
ocean_lon = ocean.sel(lat=slice(-10,0)).mean(['lat'])


def isolevel(lev1, lev2, t1, t2, t=293.15):
    return (lev1 + (t-t1)*((lev2-lev1)/(t2-t1)))



"""
for i, l in ocean_lon[i]
    iso=isolevel(ocean_lon[i, l].lev, ocean_lon[i, l+1].lev, ocean_lon[i, l], ocean_lon[i, l+1])
print iso
"""
#take vector at time=0, lat = -5 lon = 75
OV = ocean.thetao[:,:,:,:]
#print OV
#plt.figure(1)
#plt.plot(OV.values, OV.lev)
#plt.gca().invert_yaxis()
#plt.show()

#this is wrong, as the isolevel is not always between these levels
iso= isolevel(OV[:,12].lev, OV[:,13].lev, OV[:,12], OV[:,13])
plt.figure(1)
plt.pcolormesh(iso[0], vmin=90, vmax=300)
plt.colorbar()

OV=ocean.thetao[:,:,5,75]
print OV 

OV=ocean.thetao[:,:,5,75]




#create empty array to output data into, this will be good for one point
level=xr.DataArray(np.zeros((len(ocean.time))), dims=['time'], coords={'time': ocean.time})
level[:]=np.nan
#this loop goes through every timestep for this exact point and interpolates the 20 degree isotherm level.
for i, t in enumerate (ocean.thetao.time):
    print (i, t.values)
    OV=ocean.thetao[i,:,5,75]
    T=np.sort(OV.values)
    lev=np.array(OV.lev)[np.argsort(OV.values)]
    level[i] = float(np.interp(293.15, T, lev))

print level
#look at changes through time of isotherm level at this point
plt.figure(2)
plt.plot(level)
#plt.show()

#now try for longtudinal profile of the level



level=xr.DataArray(np.zeros((len(ocean.time), len(ocean.lon), len(ocean.lat))), dims=['time','lon','lat'], coords={'time': ocean.time, 'lon': ocean.lon, 'lat': ocean.lat})
level[:]=np.nan

"""
for i, t, in enumerate (ocean.thetao.time):
    print (i, t.values)
    OV=ocean.thetao[i,:,:,:]
    T=np.sort(OV.values, axis=1)
    lev=np.array(OV.lev)[np.argsort(OV.values,axis=1)]
    level[i] = float(np.interp(293.15, T, lev))

print level
"""
with ProgressBar():
    for j, t, in enumerate (ocean.lon):
        for k, t in enumerate (ocean.lat):
            for i, t, in enumerate (ocean.thetao.time):               
                OV=ocean.thetao[i,:,k,j]
                T=np.sort(OV.values)
                lev=np.array(OV.lev)[np.argsort(OV.values)]
                level[i,j,k] = float(np.interp(293.15, T, lev))

print level

plt.figure(3)
plt.plot(level.mean('time'))
plt.show()
"""

 


ocean_lon=ocean_lon.thetao[0]
T=np.sort(ocean_lon.values, axis=1)
print T
lev=np.array(ocean_lon.lev)[np.argsort(ocean_lon.values, axis = 1)]
print lev
print np.interp(293.15, T, lev)
"""
