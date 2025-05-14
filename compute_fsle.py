# FSLE computation DUACS 1/16

import os,itertools,sys
import numpy as np
import xarray as xr
from skimage.feature import peak_local_max
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon 

# Parameters and subfunctions:
from floater import rclv
from subfunctions_for_RCLV_atlas_ed import *     # Subfunctions
sys.path.append('../')
from config import *      # Parameters

import copy

import time
from netCDF4 import Dataset as datas    # Note: python is case-sensitive!

# FTLV and LAVD fortran modules
import ftlv_fsle_pdua_iv as lvddn # Fortran module lagrangian vorticity
import ftlv_fsle_pdua_iv_ekst as lvddnes # Fortran module lagrangian vorticity

# Data import 
import netCDF4 as nc
import math
import scipy.io as sio
import pickle

# Plot
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from prepvel import * # Vorticity and velocity loading
from coherency_function import test_coherency, test_coherency_low

# Create .nc file
from netCDF4 import Dataset as datas
import h5py

# Velocity field files ====================================================

anno            = np.arange(2016,2023)

# Define date_list velocities
dates_vel       = []
files_vel       = []

path_vel        = r'/media/dcortes/EXTERNAL_DRIVE_1/Diego/DUACS/'
#path_vel        = r'E:/Data/4D_MEDSEA_SSH/'
lst             = os.listdir(path_vel)
lst.sort()
for filename in lst:
        if filename.endswith(".nc"):
            dates_vel.append(filename[24:24+8])
            files_vel.append(path_vel+'/'+filename)

# Load FTLV file for lyapunov grid ========================================

path_param      = r'/media/dcortes/EXTERNAL_DRIVE_1/Diego/DUACS/nc/2016/1_24/FSLE_20160410_1_24.nc'
filedat         = xr.open_dataset(path_param)
lon,lat         = filedat.variables['lon'].T.data,filedat.variables['lat'].T.data

# Parameters:
dx, dy          = lon[0,1]-lon[0,0], lat[1,0]-lat[0,0]

path_param      = r'/media/dcortes/EXTERNAL_DRIVE_1/Diego//dt_europe_allsat_ekst_l4_20160101.nc'
filedat         = xr.open_dataset(path_param)
ut,vt           = filedat.variables['vt'].T.data,filedat.variables['ut'].T.data
lont,latt       = filedat.variables['longitude'].T.data,filedat.variables['latitude'].T.data
# Define mask
maskland,maskocean = 0,1
mascarat = np.zeros([len(lont),len(latt)])    

ut       = np.array(ut)
vt       = np.array(vt)
for li in np.arange(0,len(lont)):
    for lj in np.arange(0,len(latt)):
        if ut[li,lj] == -999 or ut[li,lj] == -999:
            mascarat[li,lj] = maskland
        else:
            mascarat[li,lj] = maskocean

# Parameters:
dx, dy          = lon[0,1]-lon[0,0], lat[1,0]-lat[0,0]

masknew = copy.deepcopy(mascara)
masknew[192:192+344,82:82+128] = mascara[192:192+344,82:82+128]+mascarat
masknew[masknew == 1] = 0


numdays = 90
print(str(numdays)+"-day integration period")

# Daily snapshot ========================================================
for tt in np.arange(92,93+365):
    print(dates_vel[tt])

    start = time.time()
    
    # Velocity peparation ---------------------------------------------------
    [ug,vg,mascara,ulon,vlat]               = pathprep_fsle(dates_vel,files_vel,tt,numdays)  
    dimlon, dimlat, dimtim                  = ug.shape[0], ug.shape[1], ug.shape[2]     # Velocity field dimensions
    
    masknew                                 = copy.deepcopy(mascara)
    masknew[192:192+344,82:82+128]          = mascara[192:192+344,82:82+128]+mascarat
    masknew[masknew == 1]                   = 0
    ug[masknew == 0], vg[masknew == 0]      = -999.999, -999.999
    
    ug, vg                                  = np.asfortranarray(ug), np.asfortranarray(vg)
    dimlyx, dimlyy, dimlyt                  = lon.shape[0],lat.shape[1],numdays*24+1    # Lagrangian variables field dimensions
    exponentelyapunov                       = lvddn.lyapunov_module.compute_fsle(dimlon,dimlat,dimtim,ug,vg,ulon,vlat,masknew,lon,lat,dimlyx,dimlyy,dimlyt,4) 
    
    end = time.time()
    # print('Computation time: ' + str(round(end - start,3)) + 's')

    # print('Save daily FSLE field')
    # file_out = r'/media/dcortes/EXTERNAL_DRIVE_1/Diego/DUACS/FSLEmask/%s/FSLE_%s_1_24.nc'%(dates_vel[tt][0:4],dates_vel[tt])
    # ncfile = datas(file_out,mode='w',format='NETCDF4')
    # print(file_out)
    
    # # Create dimensions:
    # lond, latd = ncfile.createDimension('x',dimlyx), ncfile.createDimension('y',dimlyy)
    
    # lat_file                = ncfile.createVariable('lat',np.float64,('x','y',))
    # lat_file.units          = 'degrees_north'
    # lat_file.long_name      = 'latitude'
    
    # lon_file                = ncfile.createVariable('lon',np.float64,('x','y',))
    # lon_file.units          = 'degrees_east'
    # lon_file.long_name      = 'longitude'
    
    # fsle_file               = ncfile.createVariable('fsle',np.float64,('x','y',))
    # fsle_file.units         = 's-1'
    # fsle_file.long_name     = 'finite-size lyapunov exponents'
    
    # lat_file[:,:], lon_file[:,:]    = lat, lon
    # fsle_file[:,:]                  = exponentelyapunov
    
    # ncfile.close()
    
    
    
## DUACS + Ek + St ========================================================
# Velocity field files ====================================================

anno            = np.arange(2016,2023)

# Define date_list velocities
dates_vel       = []
files_vel       = []

path_vel        = r'/media/dcortes/Elements SE/Data/AGESC_Med/daily_av/'
#path_vel        = r'E:/Data/4D_MEDSEA_SSH/'
lst             = os.listdir(path_vel)
lst.sort()
for filename in lst:
        if filename.endswith(".nc"):
            dates_vel.append(filename[25:25+11])
            files_vel.append(path_vel+'/'+filename)

# Load FTLV file for lyapunov grid ========================================

path_param      = r'/media/dcortes/Elements SE/Data/DUACS/nc/2016/1_24/FSLE_2016_1_24/FSLE_20160410_1_24.nc'
filedat         = xr.open_dataset(path_param)
lon,lat         = filedat.variables['lon'].T.data,filedat.variables['lat'].T.data

# Parameters:
dx, dy          = lon[0,1]-lon[0,0], lat[1,0]-lat[0,0]

numdays = 90
print(str(numdays)+"-day integration period")

# Daily snapshot ========================================================
for tt in np.arange(92*4+365+7,93*4+365*4,4):
    print(dates_vel[tt])

    start = time.time()
    
    # Velocity peparation ---------------------------------------------------
    [ug,vg,mascara,ulon,vlat]               = pathprep_fsle_ekst(dates_vel,files_vel,tt,numdays)  
    dimlon, dimlat, dimtim                  = ug.shape[0], ug.shape[1], ug.shape[2]     # Velocity field dimensions
        
    ug, vg                                  = np.asfortranarray(ug), np.asfortranarray(vg)
    dimlyx, dimlyy, dimlyt                  = lon.shape[0],lat.shape[1],numdays*24+1    # Lagrangian variables field dimensions
    exponentelyapunov                       = lvddnes.lyapunov_module.compute_fsle(dimlon,dimlat,dimtim,ug,vg,ulon,vlat,mascara,lon,lat,dimlyx,dimlyy,dimlyt,4) 
    
    end = time.time()
    print('Computation time: ' + str(round(end - start,3)) + 's')

    print('Save daily FSLE field')
    file_out = r'/media/dcortes/EXTERNAL_DRIVE_1/Diego/DUACS+/FSLEmask/%s/FSLE_%s_1_24.nc'%(dates_vel[tt][0:4],dates_vel[tt])
    ncfile = datas(file_out,mode='w',format='NETCDF4')
    print(file_out)
    
    # Create dimensions:
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

exponentelyapunov[exponentelyapunov==-2000] = np.nan
bounds = np.arange(0,0.5,0.01)   
plt.figure(figsize=(20,15))
ax0 = plt.axes(projection=ccrs.PlateCarree())
land = cfeature.NaturalEarthFeature(
    'physical', 'land', scale='50m', edgecolor='none', facecolor=cfeature.COLORS['land'])
ax0.add_feature(land, facecolor='black')
plt.contourf(lon,lat,abs(exponentelyapunov),bounds,cmap='Greys',transform=ccrs.PlateCarree())


fsle_file = '/media/dcortes/Elements SE/Data/4D_MEDSEA_SSH/nc/2016/1_24/FSLE_2016_1_24/FSLE_20160402_1_24.nc'
filedat         = xr.open_dataset(fsle_file)
lon,lat         = filedat.variables['lon'].T.data,filedat.variables['lat'].T.data
fsle = np.squeeze(filedat.variables['fsle'].T.data)
fsle[fsle==-1000] = 0
fsle[fsle==-2000] = 0

bounds = np.arange(0,0.5,0.001)  
plt.figure(figsize=(20,15))
ax0 = plt.axes(projection=ccrs.PlateCarree())
land = cfeature.NaturalEarthFeature(
    'physical', 'land', scale='50m', edgecolor='none', facecolor=cfeature.COLORS['land'])
ax0.add_feature(land, facecolor='black')
plt.contourf(lon,lat,abs(fsle),bounds,cmap='Greys',transform=ccrs.PlateCarree())


lonm,latm = np.meshgrid(ulon,vlat)
lontm,lattm = np.meshgrid(lont,latt)
plt.figure(figsize=(50,50))
plt.scatter(lonm,latm,s=10,color='black')
plt.scatter(lontm,lattm,s=5,color='red')

plt.figure(figsize=(50,50))
plt.contour(lonm,latm,mascara.T,[0,1])
plt.scatter(lonm,latm,s=10,c=mascara.T)
plt.contour(lontm,lattm,mascarat.T,[0,1],edgecolor = 'black')
plt.scatter(lontm,lattm,s=5,color='red')


plt.figure(figsize=(50,50))
#plt.contour(lonm,latm,masknew.T,[0,1,2])
plt.scatter(lonm,latm,s=10,c=masknew.T)
plt.contour(lonm,latm,mascara.T,[0,1],edgecolor='red')
plt.contour(lontm,lattm,mascarat.T,[0,1],edgecolor = 'black')

plt.figure(figsize=(50,50))
plt.scatter(lontm,lattm,s=10,c=mascarat.T)
plt.contour(lontm,lattm,mascarat.T,[0,1],edgecolor = 'black')

masknew = copy.deepcopy(mascara)
masknew[192:192+344,82:82+128] = mascara[192:192+344,82:82+128]+mascarat
masknew[masknew == 1] = 0


plt.figure()
plt.contourf(ug[:,:,0]*mascara)
plt.colorbar()

lonm, latm = np.meshgrid(ulon,vlat)
plt.figure(figsize=(50,50))
plt.scatter(lonm,latm,s=10,c=mascara)

