
import numpy as np
import xarray as xr

# Function velocity interpolation ==========================================

def velinterp(lon0,lat0,u,v,xlon,ylat,mascara):
    # INTERPOLATION FUNCTION U-COMPONENT
    
    # Parameters:
    Radio       = 6371000
    maskocean   = 1
    maskland    = 0

    # Grid:
    dx      = xlon[1]-xlon[0] 					                                # u-component velocity field grid
    dy 	    = ylat[1]-ylat[0]					                                # v-component velocity field grid
    
    xred 	= ((lon0-xlon[0])/dx)+1
    yred	= ((lat0-ylat[0])/dy)+1
    i       = int(xred)                                                         # x-position in lon/lat grid
    j       = int(yred)                                                         # y-position in lon/lat grid

    
    dxred	= 1
    dyred	= 1  
    tx      = (xred-i)/dxred                                                    # x-scale factor
    ux      = (yred-j)/dyred                                                    # y-scale factor
    
    # Velocity:
    v1      = u[i,j]
    v2	    = u[i+1,j]
    v3	    = u[i+1,j+1]
    v4	    = u[i,j+1]
        
    # Interpolation
    norm    = 0
    w1      = 0
		
    if int(mascara[i,j]) is maskocean:
        norm    = norm + (1 - tx) * (1 - ux)
        w1      = w1 + (1 - tx) * (1 - ux ) * v1
		
    if int(mascara[i+1,j]) is maskocean:
        norm	= norm + tx * (1 - ux)
        w1	    = w1 + tx * (1 - ux) * v2

    if int(mascara[i+1,j+1]) is maskocean:
        norm	= norm + tx * ux
        w1	    = w1 + tx * ux * v3

    if int(mascara[i,j+1]) is maskocean:
        norm	= norm + (1 - tx) * ux
        w1      = w1 + (1 - tx) * ux * v4

    if norm > 1e-5:
        w1      = w1 / norm
    else:
        w1      = 0

    w       = w1 * 86400 / (Radio * np.pi * np.cos((np.pi / 180) * lat0) / 180)     # Degrees/day
	
	# INTERPOLATION FUNCTION V-COMPONENT
    l       = int(xred)
    m       = int(yred)

    u1      = v[l,m]
    u2      = v[l+1,m]
    u3      = v[l+1,m+1]
    u4      = v[l,m+1]	

    ty      = (xred-l)/dxred
    uy      = (yred-m)/dyred

    norm    = 0
    g1      = 0
	
    if int(mascara[l,m]) is maskocean:
        norm    = norm + (1 - ty) * (1 - uy)
        g1      = g1 + (1 - ty) * (1 - uy) * u1
    
    if int(mascara[l+1,m]) is maskocean:
        norm    = norm + ty * (1 - uy)
        g1      = g1 + ty * (1 - uy) * u2
    
    if int(mascara[l+1,m+1]) is maskocean:
        norm    = norm + ty * uy
        g1      = g1 + ty * uy * u3
       
    if int(mascara[l,m+1]) is maskocean:
        norm    = norm + (1 - ty) * uy
        g1      = g1 + (1 - ty) * uy * u4

    if norm > 1e-5:
        g1      = g1/norm
    else:
        g1      = 0

    g 	= g1 * 86400 / (Radio * np.pi / 180)
    
    return w,g

# Funcion trajectory ======================================================

def pathprep(dates_vel,files_vel,tt,numdays):
    # Loading velocity loop %% 40 days
    
    # Parameters:
    Radio       = 6371000
    maskocean   = 1
    maskland    = 0
    
    ii = 0
    for filename in files_vel:
        ii += 1
        fn = filename.find(dates_vel[tt]+'_')
        if fn != -1:
            ind_ini = ii-numdays-1
            ind_fin = ii-1
            
    # Load velicity field during 40 days defore snapshot
    it = 0
    print('- Loading velocity field ...')
    for tt2 in np.arange(ind_ini,ind_fin+1):
        it += 1
        fnv = xr.open_dataset(files_vel[tt2])
        u,v = fnv.variables['ugos'].T,fnv.variables['vgos'].T
        ulon,vlat = fnv.variables['longitude'].data,fnv.variables['latitude'].data
        if it == 1:
            ug = u
            vg = v
        ug = np.append(ug,u,axis=2)
        vg = np.append(vg,v,axis=2)
        
    # vorticity computation
    
    #print("Vorticity computation")
    du_y    = np.zeros(ug.shape)
    dv_x    = np.zeros(ug.shape)    
    
    for li in np.arange(0,len(ulon)):
        for lj in np.arange(0,len(vlat)-2):
            du_y[li,lj,:] = ug[li,lj+2,:]*(180/(Radio*np.pi*np.cos(np.pi*vlat[lj+2]/180)))-ug[li,lj,:]*(180/(Radio*np.pi*np.cos(np.pi*vlat[lj]/180)))
                
    for li in np.arange(0,len(ulon)-2):
        for lj in np.arange(0,len(vlat)):
            dv_x[li,lj,:] = vg[li+2,lj,:]*(180/(Radio*np.pi)) - vg[li,lj,:]*(180/(Radio*np.pi))

    dx      = ulon[1]-ulon[0]
    dy      = vlat[1]-vlat[0]    

    dvdx    = dv_x/(2*dx)	
    dudy    = du_y/(2*dy)	
		
    zeta    = dvdx - dudy

    # Define mask
    mascara = np.zeros([len(ulon),len(vlat)])    

    u       = np.array(u)
    v       = np.array(v)
    for li in np.arange(0,len(ulon)):
        for lj in np.arange(0,len(vlat)):
            if np.isnan(zeta[li,lj,0]) == 1 or np.isnan(u[li,lj,0]) == 1 or np.isnan(v[li,lj,0]) == 1:
                mascara[li,lj] = maskland
            else:
                mascara[li,lj] = maskocean
                
    #print('Applying Atlantic and Black Sea masks ...')

    for li in np.arange(0,len(ulon)):
        for lj in np.arange(0,len(vlat)):
            if vlat[lj] >= 41 and ulon[li] <= 0:
                ug[li,lj,:]     = -999.999;
                vg[li,lj,:]     = -999.999;
                mascara[li,lj]  = maskland;
            elif vlat[lj] >= 40 and ulon[li] >= 26.77:
                ug[li,lj,:]     = -999.999;
                vg[li,lj,:]     = -999.999;  
                mascara[li,lj]  = maskland;
       
    ug      = np.where(np.isnan(ug), -999.999, ug)
    vg      = np.where(np.isnan(vg), -999.999, vg)
    zeta    = np.where(np.isnan(zeta), -999.999, zeta)
    
    return ug,vg,zeta,mascara,ulon,vlat

def pathprep_fsle(dates_vel,files_vel,tt,numdays):
    # Loading velocity loop %% 40 days
    
    # Parameters:
    Radio       = 6371000
    maskocean   = 1
    maskland    = 0
    
    ii = 0
    for filename in files_vel:
        ii += 1
        fn = filename.find(dates_vel[tt]+'_')
        if fn != -1:
            ind_ini = ii-numdays-1
            ind_fin = ii-1
            
    # Load velicity field during 40 days defore snapshot
    it = 0
    print('- Loading velocity field ...')
    for tt2 in np.arange(ind_ini,ind_fin+1):
        it += 1
        fnv = xr.open_dataset(files_vel[tt2])
        u,v = fnv.variables['ugos'].T,fnv.variables['vgos'].T
        ulon,vlat = fnv.variables['longitude'].data,fnv.variables['latitude'].data
        if it == 1:
            ug = u
            vg = v
        ug = np.append(ug,u,axis=2)
        vg = np.append(vg,v,axis=2)
        

    # Define mask
    mascara = np.zeros([len(ulon),len(vlat)])    

    u       = np.array(u)
    v       = np.array(v)
    for li in np.arange(0,len(ulon)):
        for lj in np.arange(0,len(vlat)):
            if np.isnan(u[li,lj,0]) == 1 or np.isnan(v[li,lj,0]) == 1:
                mascara[li,lj] = maskland
            else:
                mascara[li,lj] = maskocean
                
    #print('Applying Atlantic and Black Sea masks ...')

    for li in np.arange(0,len(ulon)):
        for lj in np.arange(0,len(vlat)):
            if vlat[lj] >= 41 and ulon[li] <= 0:
                ug[li,lj,:]     = -999.999;
                vg[li,lj,:]     = -999.999;
                mascara[li,lj]  = maskland;
            elif vlat[lj] >= 40 and ulon[li] >= 26.77:
                ug[li,lj,:]     = -999.999;
                vg[li,lj,:]     = -999.999;  
                mascara[li,lj]  = maskland;
       
    ug      = np.where(np.isnan(ug), -999.999, ug)
    vg      = np.where(np.isnan(vg), -999.999, vg)
    
    return ug,vg,mascara,ulon,vlat

def pathprep_fsle_ekst(dates_vel,files_vel,tt,numdays):
    # Loading velocity loop %% 40 days
    
    # Parameters:
    Radio       = 6371000
    maskocean   = 1
    maskland    = 0
    
    ii = 0
    for filename in files_vel:
        ii += 1
        fn = filename.find(dates_vel[tt]+'.nc')
        if fn != -1:
            ind_ini = ii-numdays-1
            ind_fin = ii-1
            
    # Load velicity field during 40 days defore snapshot
    it = 0
    print('- Loading velocity field ...')
    for tt2 in np.arange(ind_ini,ind_fin+1):
        it += 1
        fnv = xr.open_dataset(files_vel[tt2])
        u,v = np.expand_dims(fnv.variables['ut'].T,axis=2),np.expand_dims(fnv.variables['vt'].T,axis=2)
        ulon,vlat = fnv.variables['longitude'].data,fnv.variables['latitude'].data
        if it == 1:
            ug = u
            vg = v
        ug = np.append(ug,u,axis=2)
        vg = np.append(vg,v,axis=2)
        

    # Define mask
    mascara = np.zeros([len(ulon),len(vlat)])    

    u       = np.array(u)
    v       = np.array(v)
    for li in np.arange(0,len(ulon)):
        for lj in np.arange(0,len(vlat)):
            if u[li,lj] == -999 or v[li,lj] == -999:
                mascara[li,lj] = maskland
            else:
                mascara[li,lj] = maskocean
                
    #print('Applying Atlantic and Black Sea masks ...')

    for li in np.arange(0,len(ulon)):
        for lj in np.arange(0,len(vlat)):
            if vlat[lj] >= 41 and ulon[li] <= 0:
                ug[li,lj,:]     = -999;
                vg[li,lj,:]     = -999;
                mascara[li,lj]  = maskland;
            elif vlat[lj] >= 40 and ulon[li] >= 26.77:
                ug[li,lj,:]     = -999;
                vg[li,lj,:]     = -999;  
                mascara[li,lj]  = maskland;
       
    ug      = np.where(np.isnan(ug), -999, ug)
    vg      = np.where(np.isnan(vg), -999, vg)
    
    return ug,vg,mascara,ulon,vlat
	