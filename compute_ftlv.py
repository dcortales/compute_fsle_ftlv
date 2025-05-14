# No-slip condition at the boundaries -> Difference between 0 and nan in the boundaries

import os
import numpy as np
import xarray as xr
import time
from netCDF4 import Dataset as datas
import ftlv_fsle_computation as lvddn       # Fortran module lagrangian vorticity
from prepvel import pathprep                # Vorticity and velocity loading

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

# Lagrangian grid:
lon, lat               = np.arange(min_lon, max_lon+dx,dx), np.arange(min_lat, max_lat+dy,dy)
numdays                = 40 # Lagrangian integration interval
dimlyx, dimlyy, dimlyt = lon.shape[0],lat.shape[1],numdays*24+1    # Lagrangian variables field dimensions

print(str(numdays)+"-day integration period")

# FTLV LOOP ================================================================
for tt in np.arange(0,len(dates_vel)-numdays):
        print(dates_vel[tt])

        [ug,vg,zeta,mascara,ulon,vlat]  = pathprep(dates_vel,files_vel,tt,numdays)      # Velocity fields numdays days before tt
        dimlon, dimlat, dimtim          = ug.shape[0], ug.shape[1], ug.shape[2]         # Velocity field dimensions
        ug, vg                          = np.asfortranarray(ug), np.asfortranarray(vg)
        
        LV                              = lvddn.lyapunov_module.compute_ftlv(dimlon,dimlat,dimtim,ug,vg,ulon,vlat,mascara,XX,YY,dimlyx,dimlyy) 
        
        # Create .nc file with daily FTLV ------------------------------------------    
        print('Save daily FTLV field')
        file_out                = r'/path_out/ftlv_%s_1_24_RK_%s.nc'%(dates_vel[tt],numdays)
        ncfile                  = datas(file_out,mode='w',format='NETCDF4')
            
        lond, latd = ncfile.createDimension('x',dimlyx), ncfile.createDimension('y',dimlyy)
        
        lat_file                = ncfile.createVariable('lat',np.float64,('x','y',))
        lat_file.units          = 'degrees_north'
        lat_file.long_name      = 'latitude'
        
        lon_file                = ncfile.createVariable('lon',np.float64,('x','y',))
        lon_file.units          = 'degrees_east'
        lon_file.long_name      = 'longitude'
        
        ftlv_file               = ncfile.createVariable('ftlv',np.float64,('x','y',))
        ftlv_file.units         = 's-1'
        ftlv_file.long_name     = 'finite-time lagrangian vorticity'
        
        lat_file[:,:], lon_file[:,:]    = lat, lon
        ftlv_file[:,:]                  = LV
        
        ncfile.close()
        
        print(file_out)
        print('Saved daily FTLV field')
