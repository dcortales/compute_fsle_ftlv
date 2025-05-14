# FSLE computation DUACS 1/16

import os,itertools,sys
import numpy as np
import xarray as xr
from skimage.feature import peak_local_max
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon 

# FSLE and FTLV fortran modules
import ftlv_fsle_pdua_iv as lvddn # Fortran module lagrangian vorticity

# Data import 
import netCDF4 as nc
import math
import scipy.io as sio

# Velocity function
from prepvel import * # Vorticity and velocity loading

# Create .nc file
from netCDF4 import Dataset as datas
import h5py

# Velocity field files ====================================================

# Define date_list velocities
dates_vel       = []
files_vel       = []

path_vel        = r'/path_velocities/'
lst             = os.listdir(path_vel)
lst.sort()
for filename in lst:
        if filename.endswith(".nc"):
            dates_vel.append(filename[24:24+8]) # dates within the filename at the position 24:24+8 with the form YYYYMMDD.
            files_vel.append(path_vel+'/'+filename)

# Define grid for lyapunov exponents:
## Regular grid with velocity field resolution ##
fnv                  = xr.open_dataset(files_vel[tt2])
u,v                  = fnv.variables['ugos'].T,fnv.variables['vgos'].T
ulon,vlat            = fnv.variables['longitude'].data,fnv.variables['latitude'].data

min_lon, max_lon     = np.floor(np.min(ulon[:,0])), np.floor(np.max(ulon[:,0]))
min_lat, max_lat     = np.floor(np.min(vlat[:,0])), np.floor(np.max(vlat[:,0]))
ndim                 = 1 # Ratio between velocity resolution and lyapunov exponent resolution
dx, dy               = (ulon[0,1]-ulon[0,0])/ndim, (vlat[1,0]-vlat[0,0])/ndim

lon, lat             = np.arange(min_lon, max_lon+dx,dx), np.arange(min_lat, max_lat+dy,dy)

numdays = 90 # Lagrangian integration interval
print(str(numdays)+"-day integration period")

# Daily snapshot ========================================================
for tt in np.arange(0,len(dates_vel)):
    print(dates_vel[tt])
    
    # Velocity peparation ---------------------------------------------------
    [ug,vg,mascara,ulon,vlat]               = pathprep_fsle(dates_vel,files_vel,tt,numdays)  
    dimlon, dimlat, dimtim                  = ug.shape[0], ug.shape[1], ug.shape[2]     # Velocity field dimensions
    
    ug, vg                                  = np.asfortranarray(ug), np.asfortranarray(vg)
    dimlyx, dimlyy, dimlyt                  = lon.shape[0],lat.shape[1],numdays*24+1    # Lagrangian variables field dimensions
    
    # FLSE computation:    
    exponentelyapunov                       = lvddn.lyapunov_module.compute_fsle(dimlon,dimlat,dimtim,ug,vg,ulon,vlat,masknew,lon,lat,dimlyx,dimlyy,dimlyt,4) 

    # Create .nc file with daily FSLE ------------------------------------------    
    print('Save daily FSLE field')
    file_out                = r'/path_out/FSLE_%s_1_24.nc'%(dates_vel[tt])
    ncfile                  = datas(file_out,mode='w',format='NETCDF4')
        
    lond, latd = ncfile.createDimension('x',dimlyx), ncfile.createDimension('y',dimlyy)
    
    lat_file                = ncfile.createVariable('lat',np.float64,('x','y',))
    lat_file.units          = 'degrees_north'
    lat_file.long_name      = 'latitude'
    
    lon_file                = ncfile.createVariable('lon',np.float64,('x','y',))
    lon_file.units          = 'degrees_east'
    lon_file.long_name      = 'longitude'
    
    fsle_file               = ncfile.createVariable('fsle',np.float64,('x','y',))
    fsle_file.units         = 's-1'
    fsle_file.long_name     = 'finite-size lyapunov exponents'
    
    lat_file[:,:], lon_file[:,:]    = lat, lon
    fsle_file[:,:]                  = exponentelyapunov
    
    ncfile.close()
    

