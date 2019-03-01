import xarray as xr
import dask
import pandas as pd
import numpy as np
import scipy
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import datetime
import matplotlib.lines as mlines
import numpy.ma as ma
ocean=xr.open_dataset('thetao_.nc')
ocean_lon = ocean.sel(lat=slice(-10,0)).mean(['lat'])
#note this dataset is was remapped and latitude limted in CDO 
#it was also merged with all other times. 
"""
this is shell script used to create the set
module load cdo
cdo -mergetime /g/data1/rr3/publications/CMIP5/output1/CSIRO-BOM/ACCESS1-0/historical/mon/ocean/Omon/r1i1p1/latest/thetao/*.nc /short/e14/sm2435/thetao.nc
cdo -sellevidx,1/24 -sellonlatbox,0,360,-10,10 -remapbil,global_1 /short/e14/sm2435/thetao.nc /short/e14/sm2435/thetao_.nc
"""
#find where the 20 degree isotherm lies at a set point in lon and lat
#create empty array to output data into, this will be good for one point
level=xr.DataArray(np.zeros((len(ocean.time))), dims=['time'], coords={'time': ocean.time})
level[:]=np.nan
#this loop goes through every timestep for this exact point( 75E, 5S) and interpolates the 20 degree isotherm level.
for i, t in enumerate (ocean.thetao.time):
    print (i, t.values)
    OV=ocean.thetao[i,:,5,75]#at a specific lat lon
    T=np.sort(OV.values)#sorts the array so it is increasing (neded for interp function)
    lev=np.array(OV.lev)[np.argsort(OV.values)]#indexes the level based on the above sorting
    level[i] = float(np.interp(293.15, T, lev))#interpolates level and outputs into array created above

print level
#look at changes through time of isotherm level at this point
plt.figure(1)
plt.plot(level)
plt.show()

#now  do this loop over every lat lon to get surface of the level of 20 deg isotherm
#create new array with time lon lat for outputs
level=xr.DataArray(np.zeros((len(ocean.time), len(ocean.lon), len(ocean.lat))), dims=['time','lon','lat'], coords={'time': ocean.time, 'lon': ocean.lon, 'lat': ocean.lat})
level[:]=np.nan

#run the same loop but at all points
for i, t, in enumerate (ocean.thetao.time):
    print (i, t.values)
    OV=ocean.thetao[i,:,:,:]
    T=np.sort(OV.values)
    lev=np.array(OV.lev)[np.argsort(OV.values)]
    level[i] = float(np.interp(293.15, T, lev))


 
