module lyapunov_module
implicit real*4(a-h,o-z)

real(8)                         :: icontrol, xinicial
real(8)                         :: dx, dy, dt
real(8)                         :: h, dlx, dly, delta0, t
real(8)                         :: xk1, xk12, xk13, xk14, yk2, yk22, yk23, yk24
real(8)                         :: tiempo

integer                         :: i, j, it, itt, inidia, li, lj, jsteps, ndaysmax, nstep0, idia, ndays, nparametroh
integer                         :: mascoceano, masctierra

! Arrays
integer, parameter              :: nvec = 4                                              
real(8), dimension(nvec)        :: dist, xvini, yvini, xv, yv, distinicial                
	
! Constants
real(8), parameter              :: pi = 3.14159265359d0                                  !Pi
real(8), parameter              :: Radio = 6371000.0                                     !Earth radius
	
! Common block
common/variables/a,b,zeta,mascara,xlon,ylat,dy,dx,dt,h

contains

  subroutine compute_lagtraj(nx, ny, nt, a, b, xlon, ylat, mascara, xl, yl, mx, my, mt, xp, yp, zp, zint)
implicit real*8(a-h,o-z)

integer, intent(in)     :: nx, ny, nt
real(8)                 :: x, y, newx, newy, t, rx, ry
real(4)                 :: icontrolo
real(8), allocatable    :: zeta(:,:,:)
integer, intent(in)     :: mx, my, mt
real(8), intent(in)     :: a(:,:,:), b(:,:,:)
real(8), intent(in)     :: mascara(:,:)
real(8), intent(in)     :: xl(:,:), yl(:,:)
real(8), intent(in)     :: xlon(:), ylat(:)
real(8), intent(out)    :: xp(mx,my,mt), yp(mx,my,mt), zp(mx,my,mt), zint(mx,my)

allocate(zeta(nx,ny,nt))

     
! 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 			  		LOAD GRID DATA
! 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	

! 	SPATIAL PARAMETERS:
	
        dx              = (xlon(2) - xlon(1))        ! u-component velocity field grid
        dy              = (ylat(2) - ylat(1))        ! v-component velocity field grid
        
! 	TIME PARAMETERS

        dt              = 1.0000                     ! 1 day
        nstep0          = nt-2                       ! tiempo inicial (ej: 365+125: 5Mayo 2015)
        delta_hours     = 1
        deltat          = delta_hours/24             ! paso tiempo integracion Runge-Kutta
        ndays           = nt-2                       ! intervalo temporal integracion Runge-Kutta
        nparametroh     = floor(abs(1/deltat))
        ndaysmax        = ndays*nparametroh          ! intervalo temporal total de integracion Runge-Kutta en pasos
        h               = -deltat                    ! Runge-Kutta interval
        
! 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 				        INITIAL STEP
! 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        jsteps          = nstep0 +1 
        tiempo          = 0.0000
        
        do lj=1,my
               
        do li=1,mx
        
        ! Initial parameters
        
        icontrolo    = 0                                 ! Define control values
        x            = xl(li,lj)                         ! Initial longitude       
        y            = yl(li,lj)                         ! Initual latitude
        rx           = real(nint(x * 10.0d0**6)) / 10.0d0**6
        ry           = real(nint(y * 10.0d0**6)) / 10.0d0**6
        xp(li,lj,1)  = rx
        yp(li,lj,1)  = ry
        
        zeta_interp1 = 0.0000                            ! Initial vorticity
        
        ! Vorticity over the trajectory
        
        t            = jsteps + float(itt - 1) * 0
        zp(li,lj,1)  = zeta_interp1 + qq(t,rx,ry,nx,ny,nt,zeta,mascara,xlon,ylat)
        
!       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 				  RUNGE-KUTTA EVOLUTION
! 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
        do itt = 1,ndaysmax
			
                t       = jsteps + float(itt - 1) * h
        
                rx      = real(nint(x * 10.0d0**6)) / 10.0d0**6
                ry      = real(nint(y * 10.0d0**6)) / 10.0d0**6
                
                xk1     = w(t,rx,ry,nx,ny,nt,a,mascara,xlon,ylat)
                yk2     = g(t,rx,ry,nx,ny,nt,b,mascara,xlon,ylat)
                xk12    = w(t+0.5*h,rx+0.5*h*xk1,ry+0.5*h*yk2,nx,ny,nt,a,mascara,xlon,ylat)
                yk22    = g(t+0.5*h,rx+0.5*h*xk1,ry+0.5*h*yk2,nx,ny,nt,b,mascara,xlon,ylat)
                xk13    = w(t+0.5*h,rx+0.5*h*xk12,ry+0.5*h*yk22,nx,ny,nt,a,mascara,xlon,ylat)
                yk23    = g(t+0.5*h,rx+0.5*h*xk12,ry+0.5*h*yk22,nx,ny,nt,b,mascara,xlon,ylat)
                xk14    = w(t+h,x+h*xk13,y+h*yk23,nx,ny,nt,a,mascara,xlon,ylat)
                yk24    = g(t+h,x+h*xk13,y+h*yk23,nx,ny,nt,b,mascara,xlon,ylat)
                
        !       Remove beaching particles
					
                if ((xk1==0.0000).and.(yk2==0.0000)) then
                        icontrolo   = 1
                        zint(li,lj) = 0
                        goto 1000 ! break loop
                end if
                
                newx = xk1 + 2.0 * xk12 + 2.0 * xk13 + xk14
                newy = yk2 + 2.0 * yk22 + 2.0 * yk23 + yk24
        
                x    = rx + (h / 6.0) * newx
                y    = ry + (h / 6.0) * newy
			
        !       Break the loop if x > xlon or y > ylat
        
                if (x.ge.xlon(nx)) then 
                        icontrolo   = 1
                        zint(li,lj) = 0.0000
                        goto 1000 ! break loop
                end if
 			
                if (y.ge.ylat(ny)) then
                        icontrolo   = 1
                        zint(li,lj) = 0.0000
                        goto 1000 ! break loop
                end if
                
        !       Interpolation vorticity over trajectory
	
                zeta_interp1 = zeta_interp1 + qq(t,rx,ry,nx,ny,nt,zeta,mascara,xlon,ylat)
                tiempo       = dfloat(itt) * abs(h)
			
        !       Entire trajectories
        	
                xp(li,lj,1+itt) = x
                yp(li,lj,1+itt) = y
                zp(li,lj,1+itt) = qq(t,rx,ry,nx,ny,nt,zeta,mascara,xlon,ylat) !zeta_interp1
        
        !End loop RK evolution
        end do
        
!       Finite-Time Lagrangian Vorticity:        
        zint(li,lj) = (1 / tiempo) * zeta_interp1 * (dt * abs(h))

        1000 continue
        
!       End spatial loop
        end do
        end do
        
end subroutine compute_lagtraj

subroutine compute_ftlv(nx, ny, nt, a, b, xlon, ylat, mascara, xl, yl, mx, my, zint)

implicit real*8(a-h,o-z)

real(8), intent(in)     :: a(:,:,:), b(:,:,:)                           ! Velocities: zonal and meridional components
integer, intent(in)     :: nx, ny, nt                                   ! zonal, meridional and temporal velocity dimensions	
real(8), intent(in)     :: xlon(:), ylat(:)                             ! velocity longitude / latitude grid	
real(8), allocatable    :: zeta(:,:,:)                                  ! Relative vorticity	
real(8), intent(out)    :: zint(mx,my)                                  ! finite-time lagrangian vorticity
integer, intent(in)     :: mx,my                                        ! ftlv grid dimensions
real(8), intent(in)     :: xl(:,:), yl(:,:)                             ! ftlv longitude / latitude grid	
real(8)                 :: x, y, newx, newy, t, rx, ry                  ! Trajectory position variables
real(4)                 :: icontrolo                                    ! Control value
real(8), intent(in)     :: mascara(:,:)                                 ! Ocean mask  

real(8), parameter      :: pi=3.14159265359d0                           ! Pi
real(8), parameter      :: Radio=6371000.0                              ! Earth radius

allocate(zeta(nx,ny,nt))

! 	Compute z for the entire period

        zeta            = vort(a,b,nx,ny,nt,mascara,xlon,ylat)

! 	SPATIAL PARAMETERS:

        dx              = (xlon(2)-xlon(1))                             ! u-component velocity field grid
        dy              = (ylat(2)-ylat(1))                             ! v-component velocity field grid
        dt              = 1.0000                                        ! 1 day 
	
! 	TIME PARAMETERS

        tiempo         = 0.0000
        nstep0         = nt - 2                                         ! initial time
        delta_hours    = 1
        deltat         = delta_hours / 24                               ! paso tiempo integracion Runge-Kutta
        ndays          = nt - 2                                         ! intervalo temporal integracion Runge-Kutta
        nparametroh    = floor(abs(1 / deltat))
        ndaysmax       = ndays * nparametroh                            ! intervalo temporal total de integracion Runge-Kutta en pasos
	
! 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 					SPATIAL LOOP
! 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        jsteps         = nstep0 +1
        
        do lj = 1,my
        do li = 1,mx
        
! 	 Initial values:
        x            = xl(li,lj)                                        ! Initial longitude       
        y            = yl(li,lj)                                        ! Initial latitude
        zeta_interp1 = 0.0000                                           ! Initial vorticity
        
! 	Parameters:
        icontrolo    = 0                                                ! Define control value
        t0           = dfloat(jsteps)                                   ! Initial time
        h            = -deltat                                          ! Runge-Kutta interval
        
! 	Particule trajectory starting point (round 6 decimals)
        rx           = real(nint(x * 10.0d0**6)) / 10.0d0**6
        ry           = real(nint(y * 10.0d0**6)) / 10.0d0**6
  
        t            = jsteps + float(itt - 1) * 0
        
!       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 				  RUNGE-KUTTA EVOLUTION
! 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
        do itt = 1,ndaysmax
			
                t = jsteps + float(itt - 1) * h
        
                rx = real(nint(x * 10.0d0**6)) / 10.0d0**6
                ry = real(nint(y * 10.0d0**6)) / 10.0d0**6
                
                xk1 = w(t,rx,ry,nx,ny,nt,a,mascara,xlon,ylat)
                yk2 = g(t,rx,ry,nx,ny,nt,b,mascara,xlon,ylat)
                xk12= w(t+0.5*h,rx+0.5*h*xk1,ry+0.5*h*yk2,nx,ny,nt,a,mascara,xlon,ylat)
                yk22= g(t+0.5*h,rx+0.5*h*xk1,ry+0.5*h*yk2,nx,ny,nt,b,mascara,xlon,ylat)
                xk13= w(t+0.5*h,rx+0.5*h*xk12,ry+0.5*h*yk22,nx,ny,nt,a,mascara,xlon,ylat)
                yk23= g(t+0.5*h,rx+0.5*h*xk12,ry+0.5*h*yk22,nx,ny,nt,b,mascara,xlon,ylat)
                xk14= w(t+h,x+h*xk13,y+h*yk23,nx,ny,nt,a,mascara,xlon,ylat)
                yk24= g(t+h,x+h*xk13,y+h*yk23,nx,ny,nt,b,mascara,xlon,ylat)
                
        !       Remove beaching particles
					
                if ((xk1==0.0000).and.(yk2==0.0000)) then
                        icontrolo          = 1
                        zint(li,lj) = 0
                        goto 1000 ! break loop
                end if
        
        ! 	Update trajectory
        
                newx    = xk1 + 2.0 * xk12 + 2.0 * xk13 + xk14
                newy    = yk2 + 2.0 * yk22 + 2.0 * yk23 + yk24
        
                x       = rx + (h / 6.0) * newx
                y       = ry + (h / 6.0) * newy
			
        !       Break the loop if x > xlon or y > ylat
        
                if (x.ge.xlon(nx)) then 
                        icontrolo   = 1
                        zint(li,lj) = 0.0000
                        goto 1000 ! break loop
                end if
 			
                if (y.ge.ylat(ny)) then
                        icontrolo   = 1
                        zint(li,lj) = 0.0000
                        goto 1000 ! break loop
                end if
                
        !       Interpolation vorticity over trajectory
	
                zeta_interp1 = zeta_interp1 + qq(t,rx,ry,nx,ny,nt,zeta,mascara,xlon,ylat)
                tiempo       = dfloat(itt) * abs(h)
        
        !       End loop RK evolution
        
        end do
        
!       Finite-Time Lagrangian Vorticity:        

        zint(li,lj) = (1 / tiempo) * zeta_interp1 * (dt * abs(h))

        1000 continue
        
!       End spatial loop
        end do
        end do
        
end subroutine compute_ftlv

subroutine compute_fsle(nx, ny, nt,a,b,xlon,ylat,mascara,xl,yl,mx,my,mt,ndiv,exponentelyapunov)
implicit real*8(a-h,o-z)

integer, intent(in)     :: nx, ny, nt, ndiv
real(8)                 :: x, y, newx, newy, t, rx, ry,maska
real(4)                 :: icontrolo
integer, intent(in)     :: mx,my,mt
real(8), intent(in)     :: a(:,:,:), b(:,:,:)
real(8)                 :: aa, ba
real(8), intent(in)     :: mascara(:,:)
real(8), intent(in)     :: xl(:,:), yl(:,:)
real(8), intent(in)     :: xlon(:), ylat(:)
real(8), intent(out)    :: exponentelyapunov(mx,my)

real(8), parameter      :: pi=3.14159265359d0 !Pi
real(8), parameter      :: Radio=6371000.0 !Earth radius
integer, parameter      :: nvec = 4

!       DUMMY PARAMETERS:
        aa = a(1,1,1)
        ba = b(1,1,1)
        maska = mascara(1,1)
        x = xl(1,1)
        y = yl(1,1)
        mta = mt

newx = x
newy = y
t = 1
rx = x
ry = y
rx = xlon(1)
ry = ylat(1)
        
! 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 			  		LOAD GRID DATA
! 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	

! 	SPATIAL PARAMETERS:
	
        dx = (xlon(2)-xlon(1)) ! u-component velocity field grid
        dy = (ylat(2)-ylat(1)) ! v-component velocity field grid
        dt = 1.0000 ! 1 day

!	LYAPUNOV GRID SEPARATION
	
        dlx = dx/dfloat(ndiv)
        dly = dy/dfloat(ndiv)

!	LYAPUNOV CONDITIONS
        deg = 1
        delta0 = deg/24
        deltaf=delta0*30
	
! 	TIME PARAMETERS
	
        nstep0         = nt-2                         ! tiempo inicial (ej: 365+125: 5Mayo 2015)
	
        delta_hours    = 1
        deltat         = delta_hours/24             ! paso tiempo integracion Runge-Kutta
        ndays          = nt-2                         ! intervalo temporal integracion Runge-Kutta
        nparametroh    = floor(abs(1/deltat))
        ndaysmax       = ndays*nparametroh          ! intervalo temporal total de integracion Runge-Kutta en pasos
	
!	MASCARA PARAMETERS:

        masctierra = 0
        mascoceano = 1
        
!       DUMMY PARAMETERS:
        aa = a(1,1,1)
        ba = b(1,1,1)

! 	LYAPUNOV COORDENADADAS:

        jsteps       = nstep0 +1
        
! 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 				        INITIAL STEP
! 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
        do lj=1,my
               
        do li=1,mx
        
        icontrolo    = 0                                 ! Define control value
	
        x            = xl(li,lj)                                ! Initial longitude       
        y            = yl(li,lj)                                ! Initual latitude
		
        t0           = dfloat(jsteps)                    ! Initial time

        h            = -deltat                           ! Runge-Kutta interval
        
        ! Particule trajectory starting point (round 6 decimals)
        
        rx = real(nint(x * 10.0d0**6)) / 10.0d0**6
        ry = real(nint(y * 10.0d0**6)) / 10.0d0**6

        t = jsteps+float(itt-1)*0
	
	! Neighbours initial position
        do j=1,nvec
                vecino   = float(j)/float(nvec)
                xvini(j) = x+cos(2.0000*pi*vecino)*dlx
                yvini(j) = y+sin(2.0000*pi*vecino)*dly
        end do
	
        do j = 1,nvec
                xinicial       = dcos(y*pi/180)*dcos(yvini(j)*pi/180)*dcos((xvini(j)-x)*pi/180)&
                &+dsin(y*pi/180)*dsin(yvini(j)*pi/180)
                distinicial(j) = abs(acos(xinicial)*180.0000/pi)
        end do

        do j = 1,nvec
                xv(j) = xvini(j)
                yv(j) = yvini(j)
        end do       
        
!       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 				  RUNGE-KUTTA EVOLUTION
! 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
        do itt = 1,ndaysmax
			
                t = jsteps+float(itt-1)*h
        
                rx = real(nint(x * 10.0d0**6)) / 10.0d0**6
                ry = real(nint(y * 10.0d0**6)) / 10.0d0**6
                
                xk1 = w(t,rx,ry,nx,ny,nt,a,mascara,xlon,ylat)
                yk2 = g(t,rx,ry,nx,ny,nt,b,mascara,xlon,ylat)
                xk12= w(t+0.5*h,rx+0.5*h*xk1,ry+0.5*h*yk2,nx,ny,nt,a,mascara,xlon,ylat)
                yk22= g(t+0.5*h,rx+0.5*h*xk1,ry+0.5*h*yk2,nx,ny,nt,b,mascara,xlon,ylat)
                xk13= w(t+0.5*h,rx+0.5*h*xk12,ry+0.5*h*yk22,nx,ny,nt,a,mascara,xlon,ylat)
                yk23= g(t+0.5*h,rx+0.5*h*xk12,ry+0.5*h*yk22,nx,ny,nt,b,mascara,xlon,ylat)
                xk14= w(t+h,x+h*xk13,y+h*yk23,nx,ny,nt,a,mascara,xlon,ylat)
                yk24= g(t+h,x+h*xk13,y+h*yk23,nx,ny,nt,b,mascara,xlon,ylat)
                
        !       Remove beaching particles
					
                if ((xk1==0.0000).and.(yk2==0.0000)) then
                        icontrolo   = 1
                        exponentelyapunov(li,lj) = -2000
                        goto 1000 ! break loop
                end if
                
                newx = xk1+2.0*xk12+2.0*xk13+xk14
                newy = yk2+2.0*yk22+2.0*yk23+yk24
        
                x = rx+(h/6.0)*newx
                y = ry+(h/6.0)*newy
               	
        !       Break the loop if x > xlon or y > ylat
        
                if (x.ge.xlon(nx)) then 
                        icontrolo   = 1
                        exponentelyapunov(li,lj) = -1000
                        goto 1000 ! break loop
                end if
 			
                if (y.ge.ylat(ny)) then
                        icontrolo   = 1
                        exponentelyapunov(li,lj) = -1000
                        goto 1000 ! break loop
                end if
                
! 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 			  NEAREST NEIGHBOURS RUNGE-KUTTA EVOLUTION
! 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%					

                do jj=1,nvec!# neighbours

                        xk1= w(t,xv(jj),yv(jj),nx,ny,nt,a,mascara,xlon,ylat)
                        yk2= g(t,xv(jj),yv(jj),nx,ny,nt,b,mascara,xlon,ylat)
                        xk12= w(t+0.5*h,xv(jj)+0.5*h*xk1,yv(jj)+0.5*h*yk2,nx,ny,nt,a,mascara,xlon,ylat)
                        yk22= g(t+0.5*h,xv(jj)+0.5*h*xk1,yv(jj)+0.5*h*yk2,nx,ny,nt,b,mascara,xlon,ylat)
                        xk13= w(t+0.5*h,xv(jj)+0.5*h*xk12,yv(jj)+0.5*h*yk22,nx,ny,nt,a,mascara,xlon,ylat)
                        yk23= g(t+0.5*h,xv(jj)+0.5*h*xk12,yv(jj)+0.5*h*yk22,nx,ny,nt,b,mascara,xlon,ylat)
                        xk14= w(t+h,xv(jj)+h*xk13,yv(jj)+h*yk23,nx,ny,nt,a,mascara,xlon,ylat)
                        yk24= g(t+h,xv(jj)+h*xk13,yv(jj)+h*yk23,nx,ny,nt,b,mascara,xlon,ylat)
				
!Remove beaching particles

                        if ((xk1==0.0000).and.(yk2==0.0000)) then
                                icontrolo= 1
                                exponentelyapunov(li,lj)= -2000
                                goto 1000 !break loop
                        end if
	
                        xv(jj)= xv(jj)+(h/6.0)*(xk1+2.0*xk12+2.0*xk13+xk14)
                        yv(jj)= yv(jj)+(h/6.0)*(yk2+2.0*yk22+2.0*yk23+yk24)
               end do
	
!DISTANCE MAIN PARTICLE AND NEIGHBOURS

                        do j=1,nvec
                                xinterm= cos(y*pi/180.0000)*cos(yv(j)*pi/180.0000)*cos((xv(j)-x)*pi/180.0000)&
                                &+sin(y*pi/180.0000)*sin(yv(j)*pi/180.0000)
                                dist(j)= abs(acos(xinterm)*180.0000/pi)
                        end do
deltafcm=deltaf

! Distancia maxima entre vecinos
distanciamaxima= dist(1)
nvecino = 1
do j=1,nvec
if (distanciamaxima.lt.dist(j)) then
distanciamaxima=dist(j)
nvecino = j
end if
end do

distanciainicial= abs(distinicial(nvecino))

if (distanciamaxima.ge.deltafcm) then

tiempo= float(itt)*abs(h)
exponentelyapunov(li,lj)= 1.0000/tiempo*log(distanciamaxima/distanciainicial)
icontrolo= 1
goto 1000 !break loop
end if

        !End loop RK evolution
        end do

        1000 continue      
        
        end do
        end do
        
end subroutine compute_fsle


! 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 			   INTERPOLATION FUNCTION U-COMPONENT
! 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	

function w(t, x, y, nx, ny, nt, a, mascara, xlon, ylat) result(res)
	
implicit real*8(b-h,o-z)
integer, intent(in)             :: nx, ny, nt
real(8), intent(in)             :: t, x, y
real(8), intent(in)             :: a(nx,ny,nt)
real(8), intent(in)             :: mascara(nx,ny)
real(8), intent(in)             :: xlon(nx), ylat(ny)
real(8)                         :: norm
real(8)                         :: res, xred, yred
integer                         :: i, j, it

double precision, parameter     :: pi=3.14159265359d0                           ! Pi
double precision, parameter     :: Radio=6371000.0                              ! Earth radius
real(8)                         :: dx, dy, dt, masctierra, mascoceano

        masctierra = 0
        mascoceano = 1

! 	Grid parameters

        dx      = xlon(2)-xlon(1) ! u-component velocity field grid
        dy      = ylat(2)-ylat(1) ! v-component velocity field grid
        dt      = 1.0000          ! 1 day

        xred    = ((x-xlon(1))/dx)+1.0000
        yred    = ((y-ylat(1))/dy)+1.0000

!       Integrers

        i       = int(xred)
        j       = int(yred)
        it      = int(t)

!       Interpolation in time T

        v1      = a(i,j,it)
        v2      = a(i+1,j,it)
        v3      = a(i+1,j+1,it)
        v4      = a(i,j+1,it)
   
        dxred   = 1.0000
        dyred   = 1.0000
	
        tx      = (xred - dfloat(i)) / dxred
        ux      = (yred - dfloat(j)) / dyred
	
        norm    = 0.0000
        w1      = 0.0000
			
        if (mascara(i,j).gt.0.5) then
                norm    = norm + (1.0 - tx) * (1.0 - ux)
                w1      = w1 + (1.0 - tx) * (1.0 - ux) * v1
        end if
		
        if (mascara(i+1,j).gt.0.5) then
                norm    = norm + tx * (1.0 - ux)
                w1      = w1 + tx * (1.0 - ux) * v2
        end if

        if (mascara(i+1,j+1).gt.0.5) then
                norm    = norm + tx * ux
                w1      = w1 + tx * ux * v3
        end if

        if (mascara(i,j+1).gt.0.5) then
                norm    = norm + (1.0 - tx) * ux
                w1      = w1 + (1.0 - tx) * ux * v4
        end if

        if (norm>1e-5) then
                w1      = w1 / norm
        else
                w1      = 0.0000
        end if

!       Interpolation in time T+1

        norm    = 0.0000
        w2      = 0.0000

        v12     = a(i,j,it+1)
        v22     = a(i+1,j,it+1)
        v32     = a(i+1,j+1,it+1)
        v42     = a(i,j+1,it+1)

        tx2     = (xred - dfloat(i)) / dxred
        ux2     = (yred - dfloat(j)) / dyred
	
        if (mascara(i,j).gt.0.5) then
                norm    = norm + (1.0 - tx2) * (1.0 - ux2)
                w2      = w2 + (1.0 - tx2) * (1.0 - ux2) * v12
        end if

        if (mascara(i+1,j).gt.0.5) then
                norm    = norm + tx2 * (1.0 - ux2)
                w2      = w2 + tx2 * (1.0 - ux2) * v22
        end if
		
        if (mascara(i+1,j+1).gt.0.5) then
                norm    = norm + tx2 * ux2
                w2      = w2 + tx2 * ux2 * v32
        end if

        if (mascara(i,j+1).gt.0.5) then
                norm    = norm + (1.0 - tx2) * ux2
                w2      = w2 + (1.0 - tx2) * ux2 * v42
        end if

        if (norm > 1e-5) then
                w2      = w2 / norm
        else
                w2      = 0.0000
        end if

!       Interpolation T - T+1

        px      = (t - dfloat(it)) / dt
        res     = ((1.0000 - px) * w1 + px * w2) * 86400 / (Radio * pi * cos((pi / 180.0) * y) / 180.0000)               ! Degrees/day
	
end function w
	
	
! 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 			   INTERPOLATION FUNCTION V-COMPONENT
! 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function g(t,x,y,nx,ny,nt,b,mascara,xlon,ylat) result(res)

implicit real*8(a-h,o-z)
integer, intent(in):: nx, ny, nt
real(8), intent(in):: t,x,y
real(8), intent(in):: b(nx,ny,nt)
real(8), intent(in):: mascara(nx,ny)
real(8), intent(in):: xlon(nx), ylat(ny)
real(8):: norm
real(8):: res, xred, yred
integer:: it, l, m

double precision, parameter :: pi=3.14159265359d0 !Pi
double precision, parameter :: Radio=6371000.0 !Earth radius
real(8) :: dx, dy, dt, masctierra, mascoceano
	
	
masctierra = 0
mascoceano = 1

dx = xlon(2)-xlon(1) ! u-component velocity field grid
dy = ylat(2)-ylat(1) ! v-component velocity field grid
dt = 1.0000 ! 1 day

xred= ((x-xlon(1))/dx)+1.0000
yred= ((y-ylat(1))/dy)+1.0000

dxred= 1.0000
dyred= 1.0000

l= int(xred)
m= int(yred)
it= int(t)

!Interpolation in time T

u1= b(l,m,it)
u2= b(l+1,m,it)
u3= b(l+1,m+1,it)
u4= b(l,m+1,it)

ty= (xred-dfloat(l))/dxred
uy= (yred-dfloat(m))/dyred

	
norm= 0.0000
g1= 0.0000
	
!if (m.lt.ny) then
if (mascara(l,m).gt.0.5) then
norm = norm+(1.0-ty)*(1.0-uy)
g1= g1+(1.0-ty)*(1.0-uy)*u1
end if

if (mascara(l+1,m).gt.0.5) then
norm= norm+ty*(1.0-uy)
g1= g1+ty*(1.0-uy)*u2
end if
		
if (mascara(l+1,m+1).gt.0.5) then
norm= norm+ty*uy
g1= g1+ty*uy*u3
end if

if (mascara(l,m+1).gt.0.5) then
norm= norm+(1.0-ty)*uy
g1= g1+(1.0-ty)*uy*u4
end if


if (norm>1e-5) then
g1=g1/norm
else
g1=0.0000
end if

!Interpolation in time T+1

u12= b(l,m,it+1)
u22= b(l+1,m,it+1)
u32= b(l+1,m+1,it+1)
u42= b(l,m+1,it+1)

ty2= (xred-dfloat(l))/dxred
uy2= (yred-dfloat(m))/dyred

norm= 0.0000
g2 = 0.0000
		
if (mascara(l,m).gt.0.5) then
norm= norm+(1.0-ty2)*(1.0-uy2)
g2= g2+(1.0-ty2)*(1.0-uy2)*u12
end if
	      
if (mascara(l+1,m).gt.0.5) then
norm= norm+ty2*(1.0-uy2)
g2= g2+ty2*(1.0-uy2)*u22
end if

if (mascara(l+1,m+1).gt.0.5) then
norm= norm+ty2*uy2
g2= g2+ty2*uy2*u32
end if

if (mascara(l,m+1).gt.0.5) then
norm= norm+(1.0-ty2)*uy2
g2= g2+(1.0-ty2)*uy2*u42
end if


if (norm>1e-5) then
g2=g2/norm
else
g2=0.0000
end if
	
! Interpolation T - T+1

py = (t-dfloat(it))/dt
res= ((1.0000-py)*g1+py*g2)*86400.00/(Radio*pi/180.0000)
	
	
end function g
	

! 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 			   INTERPOLATION FUNCTION SCALAR
! 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	

function qq(t,x,y,nx,ny,nt,zeta,mascara,xlon,ylat) result(res)

implicit real*8(a-h,o-z)
integer, intent(in):: nx, ny, nt
real(8), intent(in):: t,x,y
real(8), intent(in):: zeta(nx,ny,nt)
real(8), intent(in):: mascara(nx,ny)
real(8), intent(in):: xlon(nx), ylat(ny)
real(8):: norm
real(8):: res, xred, yred
integer:: i, j, it, mascoceano, masctierra

real(8), parameter :: pi=3.14159265359d0 !Pi
real(8), parameter :: Radio=6371000.0 !Earth radius
real(8) :: dx, dy, dt, dxred, dyred

masctierra = 0
mascoceano = 1
controlo = 0

dx = xlon(2)-xlon(1) ! u-component velocity field grid
dy = ylat(2)-ylat(1) ! v-component velocity field grid
dt = 1.0000 ! 1 day

xred= ((x-xlon(1))/dx)+1.0000
yred= ((y-ylat(1))/dy)+1.0000

i= int(xred)
j= int(yred)
it= int(t)
res=0

if (it==0) then
       icontrolo=1
        !exponentelyapunov(li,lj)=0.0000
       goto 1000! break loop
       end if
	
! Interpolation in time T

var1= zeta(i,j,it)
var2= zeta(i+1,j,it)
var3= zeta(i+1,j+1,it)
var4= zeta(i,j+1,it)
   
dxred= 1.0000
dyred= 1.0000
	
tx= (xred-dfloat(i))/dxred
ux= (yred-dfloat(j))/dyred

norm= 0.0000
s1= 0.0000
	
if (mascara(i,j).gt.0.5) then
norm=norm+(1.0-tx)*(1.0-ux)
s1= s1+(1.0-tx)*(1.0-ux)*var1
end if
		
if (mascara(i+1,j).gt.0.5) then
norm=norm+tx*(1.0-ux)
s1= s1+tx*(1.0-ux)*var2
end if

if (mascara(i+1,j+1).gt.0.5) then
norm=norm+tx*ux
s1= s1+tx*ux*var3
end if

if (mascara(i,j+1).gt.0.5) then
norm=norm+(1.0-tx)*ux
s1= s1+(1.0-tx)*ux*var4
end if

if (norm>1e-5) then
s1=s1/norm
else
s1=0.0000
end if

! Interpolation in time T+1

norm= 0.0000
s2= 0.0000
	
if (it.lt.nt) then
var12= zeta(i,j,it+1)
var22= zeta(i+1,j,it+1)
var32= zeta(i+1,j+1,it+1)
var42= zeta(i,j+1,it+1)

tx2= (xred-dfloat(i))/dxred
ux2= (yred-dfloat(j))/dyred
	
if (mascara(i,j).gt.0.5) then
norm=norm+(1.0-tx2)*(1.0-ux2)
s2= s2+(1.0-tx2)*(1.0-ux2)*var12
end if

if (mascara(i+1,j).gt.0.5) then
norm=norm+tx2*(1.0-ux2)
s2= s2+tx2*(1.0-ux2)*var22
end if
		
if (mascara(i+1,j+1).gt.0.5) then
norm=norm+tx2*ux2
s2= s2+tx2*ux2*var32
end if
 
if (mascara(i,j+1).gt.0.5) then
norm=norm+(1.0-tx2)*ux2
s2= s2+(1.0-tx2)*ux2*var42
end if
end if

if (norm>1e-5) then
s2=s2/norm
else
s2=0.0000
end if

! Interpolation T - T+1
px= (t-dfloat(it))/dt
res= ((1.0000-px)*s1+px*s2)

1000 continue
	
end function qq


function vort(a, b, nx, ny, nt, mascara, xlon, ylat) result(res)

implicit real*8(a-h,o-z)
integer, intent(in)             :: nx, ny, nt
real(8), intent(in)             :: a(nx,ny,nt), b(nx,ny,nt)
real(8), intent(in)             :: mascara(nx,ny)
real(8), intent(in)             :: xlon(nx), ylat(ny)
real(8),dimension(nx,ny,nt)     :: res
integer                         :: li, lj, it

real(8), parameter              :: pi = 3.14159265359d0                                                ! Pi
real(8), parameter              :: Radio = 6371000.0                                                   ! Earth radius
real(8)                         :: dx, dy

real(8), allocatable            :: zeta(:,:,:), du_y(:,:,:), dv_x(:,:,:), dvdx(:,:,:), dudy(:,:,:)

allocate(zeta(nx,ny,nt))
allocate(du_y(nx,ny,nt))
allocate(dv_x(nx,ny,nt))
allocate(dvdx(nx,ny,nt))
allocate(dudy(nx,ny,nt))

masctierra = 0
mascoceano = 1


dx = (xlon(2)-xlon(1)) ! u-component velocity field grid
dy = (ylat(2)-ylat(1)) ! v-component velocity field grid

do li=1,nx-1
        do lj = 1,ny-1
                du_y(li,lj,:) = 0.0000
                dv_x(li,lj,:) = 0.0000
        end do
end do
	
do li=1,nx
        do lj = 1,ny-2
                if (mascara(li,lj)==mascoceano) then
                        if (mascara(li+1,lj)==mascoceano) then
                                if (mascara(li+1,lj+1)==mascoceano) then
                                        if (mascara(li,lj+1)==mascoceano) then
                                                du_y(li,lj,:) = a(li,lj+2,:)*(float(180)/(Radio*pi*cos(pi*ylat(lj+2)/float(180))))& 
                                                & - a(li,lj,:)*(float(180)/(Radio*pi*cos(pi*ylat(lj)/float(180))))
                                        end if
                                end if
                        end if
                end if
        end do
end do
	
do li=1,nx-2
        do lj = 1,ny
                if (mascara(li,lj)==mascoceano) then
                        if (mascara(li+1,lj)==mascoceano) then
                                if (mascara(li+1,lj+1)==mascoceano) then
                                        if (mascara(li,lj+1)==mascoceano) then
                                                dv_x(li,lj,:) = b(li+2,lj,:)*(180/(Radio*pi))- b(li,lj,:)*(180/(Radio*pi))
                                        end if
                                end if
                        end if
                end if
        end do 
end do
	
do i=1,nx
        do j = 1,ny
                do it = 1,nt
                        dvdx(i,j,it) = 0.0000
                        dudy(i,j,it) = 0.0000
                end do
        end do
end do

do i=1,nx-2
        do j = 2,ny-1
                do it = 1,nt
                        dvdx(i+1,j,it) = dv_x(i,j,it) / (2 * dx)
                end do
        end do
end do
	
do i=2,nx-1
        do j = 1,ny-2
                do it = 1,nt
                        dudy(i,j+1,it) = du_y(i,j,it) / (2 * dy)
                end do
        end do
end do
	
res = dvdx - dudy

end function vort

end module lyapunov_module
	
