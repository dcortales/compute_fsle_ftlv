# No-slip condition at the boundaries -> Difference between 0 and nan in the boundaries

import os
import numpy as np
import xarray as xr
import time
from netCDF4 import Dataset as datas
import ftlv_fsle_computation as lvddn       # Fortran module lagrangian vorticity
from prepvel import pathprep             # Vorticity and velocity loading

# Velocity field files ====================================================

anno            = np.arange(2022,2023)
month           = np.arange(1,13)

# DUACS 1/4

# Define date_list velocities

dates_vel       = []
files_vel       = []

for an in anno:
    for m in month:
        path_vel        = r'/media/dcortes/EXTERNAL_DRIVE_1/Diego/DUACS16/SEALEVEL_EUR_PHY_L4_MY_008_068/%s/%02d'%(an,m)
        lst             = os.listdir(path_vel)
        lst.sort()
        for filename in lst:
            if filename.endswith(".nc"):
                dates_vel.append(filename[24:24+8])
                files_vel.append(path_vel+'/'+filename)

# Load FTLV file for lyapunov grid ========================================

path_param      = r'/media/dcortes/EXTERNAL_DRIVE_1/Diego/FTLV_20160401_1_24.nc'
#path_param      = r'E:/Data/4D_MEDSEA_SSH/nc/FTLV/2016/1_24/FTLV_20160401_1_24.nc'
filedat         = xr.open_dataset(path_param)
XX,YY           = filedat.variables['lon'].T.data,filedat.variables['lat'].T.data


for numdays in np.arange(10,)

# Lagrangian computation grid:
dimlyx, dimlyy, dimlyt      = XX.shape[0],XX.shape[1],numdays*24+1          # Lagrangian variables field dimensions

# FTLV LOOP ================================================================
for tt in np.arange(0,9313):
        print(dates_vel[tt])
        
        start_time                      = time.time()
        [ug,vg,zeta,mascara,ulon,vlat]  = pathprep(dates_vel,files_vel,tt,numdays)      # Fields n days before ftlv
        dimlon, dimlat, dimtim          = ug.shape[0], ug.shape[1], ug.shape[2]         # Velocity field dimensions
        ug, vg                          = np.asfortranarray(ug), np.asfortranarray(vg)
                            
        print('- Computing FTLV ...')
        
        LV                              = lvddn.lyapunov_module.compute_ftlv(dimlon,dimlat,dimtim,ug,vg,ulon,vlat,mascara,XX,YY,dimlyx,dimlyy) 

        end_time                        = time.time()
        elapsed_time                    = end_time - start_time
        
        print("- Computation time: ", round(elapsed_time,3),"s") 
        
        # Save .nc file:
        #file_out                        = r'/media/dcortes/EXTERNAL_DRIVE_1/Diego/DUACS16/nc/ftlv_atlas/%s/duacs16_ftlv_%s.nc'%(dates_vel[tt][0:4],dates_vel[tt])
        #file_out                        = r'/media/dcortes/EXTERNAL_DRIVE_1/Diego/DUACS16/nc/calypso2018_case_study/duacs16_ftlv_%s_RK_%s.nc'%(dates_vel[tt],numdays)
        file_out                        = r'/media/dcortes/EXTERNAL_DRIVE_1/Diego/DUACS16/nc/coherency_study/numdays/duacs16_ftlv_%s_RK_%s.nc'%(numdays,dates_vel[tt],numdays)
        ncfile                          = datas(file_out,mode='w',format='NETCDF4')
        
        # Create dimensions:
        lond, latd                      = ncfile.createDimension('x',dimlyx), ncfile.createDimension('y',dimlyy)
        
        lat_file                        = ncfile.createVariable('lat',np.float64,('x','y',))
        lat_file.units                  = 'degrees_north'
        lat_file.long_name              = 'latitude'
        
        lon_file                        = ncfile.createVariable('lon',np.float64,('x','y',))
        lon_file.units                  = 'degrees_east'
        lon_file.long_name              = 'longitude'
        
        ftlv_file                       = ncfile.createVariable('ftlv',np.float64,('x','y',))
        ftlv_file.units                 = 's-1'
        ftlv_file.long_name             = 'finite-time lagrangian vorticity (velocity boundary = NaN)'
        
        lat_file[:,:], lon_file[:,:]    = YY, XX
        ftlv_file[:,:]                  = LV
        
        ncfile.close()
        
        print(file_out)
        print('Saved daily FTLV field')