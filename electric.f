      PROGRAM MANYDROP 
c
c  Iterative Code for the solution to the time-evolving drop problem
c  in unbounded two-dimensional Stokes flow.
c
c  Runge Kutta is the time-stepping method
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c
      parameter (ndmax=3000, kmax = 150, nmax = ndmax*kmax)
c
c     NDMAX is the max number of points allowed on each hole boundary
c     KMAX is the maximum number of holes.
c     NMAX is the maximum dimension of the unknown density vector.
c          It is twice the size of the maximum number of points, since
c          the density is a complex quantity.
c
c  Initial Curve Description Variables
c
      dimension ai(kmax), bi(kmax), thetax(kmax), ncyc(kmax), 
     *          area0(kmax), area(kmax), nd(kmax)
      dimension alph(nmax)
      complex*16 zk(kmax)
c
c  Geometry arrays
c
      dimension rkappa(nmax), dsdth(nmax), arcl0(kmax),
     *          arcl(kmax)
      complex*16 z(ndmax,kmax), dz(nmax), u(nmax), zbond
c
c  Theta-L arrays
c
      dimension beta(nmax)
c
c  FFT arrays
c
      parameter (ifwork = 10*ndmax+50, nzwork=32*nmax+100)
      dimension amode(nmax), fmode(nmax)
      complex*16 za(nmax), zai(nmax), wsave(ifwork,kmax), zwork(nzwork)
      complex*16 z2(nmax)
c
c  ARRAYS FOR RUNGE KUTTA
c
      parameter (neqmax = 2*nmax)
      dimension ybar(neqmax), phi(4*neqmax), y(neqmax), dy(neqmax)
c
c  TIME ARRAYS
c
      parameter (ntmax = 100000)
      dimension tarray(ntmax), tarea(ntmax)  
c      
      common /problem/ zk, nd, k, nbk, neq
      common /density/ u
      common /flowinfo/ cap, zbond, E0
      common /wallinfo / xleft, slength, iwall  
      common /doublend/ icut       
c
      REAL*4 TIMEP(2), ETIME
      character*1 iresp
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0, 1.d0)
c
c
c  Read in initial data 
c
         call READINI (k, nd, nbk, neq, ai, bi, ncyc, zk, thetax, dt,  
     1                 tfinal, iplot, ipost, nplot, cap, zbond, iwall,
     2                 E0)
c
c  Construct initial bubble shape
c
         idump = 0
         call PRIN2 ('ai = *', ai, k)
         call PRIN2 ('bi = *', bi, k)
         if (idump.eq.0) then
           call GEOINIT (k, nd, nbk, ai, bi, thetax, ncyc, zk, 
     *                   z, dsdth, rkappa, iplot, area0, ndmax)
            call CLUSTPARAM (k, nd, nbk, dsdth, arcl0, beta, 
     *                       za, zai, wsave, ifwork)
            call GEOSHIFT (k, nd, nbk, ai, bi, thetax, ncyc, zk,  
     *                     beta, z, rkappa, iplot, ndmax)
ccc            call DUMPDATA (k, nd, nbk, neq, area0, z, ndmax)
            istart = 0
            cmin = 0.d0
            ntime = 0
            do kbod = 1, k
ccc               call PRINF ('KBOD = *', kbod, 1)
               ist = istart+1
               rkmax = 0.d0
               do i = 1, nd(kbod)
                  alph(i) = 2.d0*pi*(i-1.d0)/nd(kbod)
                  rkmax = max(rkmax, dabs(rkappa(istart+i)))
               end do
               call PRIN2 ('    INITIAL KAPPA MAX = *', rkmax, 1)  
               ds = arcl0(kbod)/(2.d0*pi)
               call PRIN2 ('    INITIAL Salpha = *', ds, 1)
               call RSPLOT (alph, rkappa(ist), nd(kbod), 1, 21)
               call RSCPLOT (z(1,kbod), nd(kbod), 1, iplot)
ccc               call RSZ_MOVIE (z(1,kbod),0.d0, nd(kbod), kbod, k,  
ccc     1                         ntime, iplot)
               coeff = dt*nd(kbod)/arcl0(kbod)
               cmin = max(cmin,coeff)
               istart = istart+nd(kbod)
            end do
           else
            call READDATA (k, nd, nbk, neq, area0, arcl0, z, wsave, 
     *                     ndmax,ifwork, dt, tstart)
            tfinal = tstart + tfinal
            call PRIN2 ('  dt = *', dt, 1)
            call PRIn2 ('  tstart = *', tstart, 1)
            call PRIn2 ('  tfinal = *', tfinal, 1)
            call PRIN2 ('  area0 = *', area0, k)
            call PRIN2 ('  arcl0 = *', arcl0, k)
            cmin = 0.d0
            do kbod = 1, k
               call RSCPLOT (z(1,kbod), nd(kbod), kbod, iplot)
               coeff = dt*nd(kbod)/arcl0(kbod)
               cmin = max(cmin,coeff)
               istart = istart+nd(kbod)
            end do
         end if
         call PRIN2 (' cmin = *', cmin, 1)
c
c  Plotting data
c
         dtplot = tfinal/nplot
         tplot = dtplot
         eps = 1.d-6         
c
c  Set input data for IVP solver
c
         tbeg = etime(timep)
 200     call INITCOND (k, ndmax, nd, nbk, neq, z, ybar)
c  
c  Time integration:  Loops until time > tfinal
c
         areatot0 = 0.d0
         do kbod = 1, k
            areatot0 = areatot0 + area0(kbod)
         end do
         call PRIN2 ('  TOTAL INITIAL AREA = *', areatot0, 1)
         idoub = 0
         time = tstart 
  100    time = time + dt
         if (time.gt.tfinal) then
            dt = tfinal - time + dt
            time = tfinal 
         end if
         ntime = ntime + 1
            call PRINd2 ('+++++++++++++++++++++  TIME = *',time,1)
            call PRIN2 ('         DT = *',dt,1)
            call PRINf ('    ND = *',nd,k)
            call PRINF (' nplot = *', nplot, 1)
c
            call RK2 (k, neq, time, dt, phi, ybar, y, dy)
c 
c  Now compute shape
c
            icut = 0
            istart = 0
            istart2 = 0
            areatot = 0.d0
            akapmax = 0.d0
            errsmax = 0.d0
            do kbod = 1, k
               n = nd(kbod)/2
ccc               call PRINF ('kbod = *', kbod, 1)
               ist = istart+1
               ist2 = istart2+1
               call EXTRACTSOL (nd(kbod), neq, area(kbod), arcl(kbod),  
     *                          z(1,kbod), dz(ist), ybar(ist2), za, 
     *                          zai, wsave(1,kbod), errds, rkmax, 
     *                          rkappa)
               areatot = areatot + area(kbod)
ccc               call PRIN2 ('  ARC LENGTH = *', arcl, 1)
ccc               call PRIN2 ('     ERR DS = *', errds, 1)
ccc               call PRIN2 ('  KAPPA MAX = *', rkmax, 1)
               akapmax = max(akapmax,rkmax)
               errsmax = max(errsmax,errds)
               call GETFOURIER (nd(kbod), amode, fmode, z(1,kbod), za,
     *                           wsave(1,kbod), idoub)
               if (idoub.eq.1) then
                   call RSLOGPLOT(amode,fmode,n,1,50) 
ccc                   stop
               end if
               if (idoub.eq.1) then
                  icut = 1
                  call PRINF ('  DOUBLING ND ON KBOD = *', kbod, 1)
                  call DOUBLE (nd(kbod), z(1,kbod), z2, za, 
     *                         wsave(1,kbod), zwork, nzwork, ndmax,
     *                         istop)
                  n = nd(kbod)/2
                  if (istop.eq.1) then
                     nd(kbod) = nd(kbod)/2
                     n = nd(kbod)/2
                     call GETFOURIER (nd(kbod), amode, fmode, 
     *                           z(1,kbod), za, wsave(1,kbod), idoub)
ccc                     call RSLOGPLOT (amode, fmode, n, 1, 50)
                     stop
                  end if
                  call PRINF ('  new nd = *', nd(kbod), 1)
                  if (mod(ntime,nplot).eq.0) then
                     call RSCPLOT (z(1,kbod), nd(kbod), 1, iplot)
ccc                     call RSZ_MOVIE (z(1,kbod),time, nd(kbod), kbod, k,  
ccc     1                               ntime, iplot)
                  end if
                  istart = istart+nd(kbod)/2
                  istart2 = istart2+nd(kbod)
                 else
                  if (mod(ntime,nplot).eq.0) then
                     call RSCPLOT (z(1,kbod), nd(kbod), 1, iplot)
                     call RSLOGPLOT (amode, fmode, n, 1, 50)
ccc                     call RSZ_MOVIE (z(1,kbod),time, nd(kbod), kbod, k,  
ccc     1                               ntime, iplot)
                  end if
                  istart = istart+nd(kbod)
                  istart2 = istart2+2*nd(kbod)
               end if
            end do
ccc
c   compute new time step
c
            if (icut.eq.1) then
               call PRIn2 ('  arcl0 = *', arcl0, k)
               call PRIn2 ('  arcl = *', arcl, k)
               dtmin = 1.d10
               nbod_max = 0
               do kbod = 1, k
                  dtkbod = cmin*arcl(kbod)/nd(kbod)
                  dtmin = min(dtmin,dtkbod)
                  nbod_max = max(nbod_max,nd(kbod))
               end do
               nplot = 2*nplot
               dt = dtmin            
               call PRIN2 ('  DOUBLED ND, NEW DT = *', dt, 1)
               call INITCOND (k, ndmax, nd, nbk, neq, z, ybar)
              else
               dtmin = 1.d10
               do kbod = 1, k
                  dtkbod = cmin*arcl(kbod)/nd(kbod)
                  dtmin = min(dtmin,dtkbod)
               end do
               dt = dtmin            
            end if
            erra = dabs(areatot0-areatot)/areatot0
ccc            call PRIN2 ('  KAPPA MAX = *', akapmax, 1)
ccc            call PRIN2 ('  ERR DS MAX = *', errsmax, 1)
            call PRIN2 ('  TOTAL AREA = *', areatot, 1)
            call PRIN2 ('  ERROR IN TOTAL AREA = *', erra, 1)
         if (time.lt.tfinal-1.d-8) goto 100
         tend = etime(timep)
         call PRIND2 ('TOTAL TIME IN IVP = *',tend-tbeg,1)
         istart = 0
         do kbod = 1, k
            ist = istart+1
ccc            call RSZ_MOVIE (z(1,kbod),time, nd(kbod), kbod, k,  
ccc     1                               ntime, iplot)
            call RSCPLOT (z(1,kbod), nd(kbod), 1, iplot)
            call GETFOURIER (nd(kbod), amode, fmode, z(1,kbod), za, 
     *                       wsave(1,kbod), idoub)
            call RSLOGPLOT (amode, fmode, n, 1, 50)
            istart = istart+nd(kbod)
         end do
         close(21)
         call DUMPDATA (k, nd, nbk, neq, area, arcl, z, ndmax,
     *                  dt,time)
c
      STOP
      END
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine READINI (k, nd, nbk, neq, ai, bi,ncyc, zk, thetax, dt,  
     *                    tfinal, iplot, ipost, nplot, cap, zbond, 
     *                    iwall, E0)
c---------------
c
      implicit real*8 (a-h,o-z)
      dimension ai(*), bi(*), thetax(*), ncyc(*), nd(*)
      complex*16 zk(*), zbond
c
c     ND        = number of points on each hole boundary
c     K         =  number of hole boundaries
c     DT, NTIME = time stepping parameters
c     AI, BI    = axis of ellipse or starfish (BI = eps)
c     ZK        = hole centres
c     NCYC      = number of arms in starfish
c                 if NCYC = 0, the hole is an ellipse
c     THETAX    = orientation of ellipse (rotation angle measured
c                 from positive x-axis).
c     RADIUS    = a maximal radius for each particle.
c
         pi = 4.d0*datan(1.d0)
c         
         call PRINI (6, 13)
         OPEN(2, FILE='bubbin',
     1        status='OLD')
         READ (2,*) k
         CALL PRINF(' K = *',K,1)
         read (2,*) (nd(kk), kk=1,k)
         nbk = nd(1)
         CALL PRINF(' ND = *',ND,k)
         read (2,*) dt,tfinal,nplot
         call PRINF (' NPLOT = *',nplot, 1)
         call PRIN2 (' DTIME = *',dt,1)
         call PRIN2 (' TFINAL = *',tfinal,1)
         read (2,*) cap, bondx, bondy, E0
         zbond = dcmplx(bondx, bondy)
         call PRIN2 (' cap = *', cap, 1)
         call PRIN2 (' zbond = *', zbond, 2)
         nbk = 0
         do kbod = 1,k
            read (2,*) ai(kbod),bi(kbod),thetax(kbod),cx,cy,ncyc(kbod)
            zk(kbod) = dcmplx(cx,cy)
            nbk = nbk + nd(kbod)
         end do
         CALL PRINF(' NBK = *',NBK,1)
         neq = 2*nbk
         if (cdabs(zbond).gt.0.d0) then
            iwall = 1
            call PRINF ('  WALL BOUNDED *', 1, 1)
           else
            iwall = 0
         end if
         close (2)
c
         iplot = 10
         open (unit = iplot, file = 'bubble.m')
         ipost = 11
         open (unit = ipost, file = 'time.data')
         open (unit=21, file = 'kappa.m')
         open (unit=50, file = 'fourier.m')
c 
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine GEOINIT (k, nd, nbk, ai, bi, thetax, ncyc, zk,
     *                    z, dsdth, rkappa, iplot, area0, nmax)
c
c---------------
c
c  Provides initial curve description for resampler
c  The number of points on this curve is some factor of nd
c
      implicit real*8 (a-h, o-z)
      dimension ai(k), bi(k), thetax(k), ncyc(k), area0(k)
      dimension dsdth(nbk), rkappa(nbk), nd(k)
      complex*16 z(nmax,k), zk(k), eye
c
         eye = dcmplx(0.d0, 1.d0)
         pi = 4.d0*datan(1.d0)
c
c
         istart = 0
         do kbod = 1,k
            dth = 2.d0*pi/nd(kbod)
            if (ncyc(kbod).eq.0) then
               area0(kbod) = pi*ai(kbod)*bi(kbod)
              else
               area0(kbod) = pi*(ai(kbod)**2 - ncyc(kbod)*bi(kbod)**2)
            end if
            call PRINF ('+++++  NBOD = *',kbod,1)
            call PRIN2 ('  INITIAL AREA = *',area0(kbod),1)
            do i= 1, nd(kbod)
c
c     Calculate geometry variables at each subinterval endpoint
c
               index = istart + i
               th = dth*(i-1.d0)
               call INITPARAM (th, ai(kbod), bi(kbod), thetax(kbod),
     *                         ncyc(kbod), xp, yp, xdot, ydot,
     *                         xddot, yddot)
c
               z(i,kbod) = zk(kbod) + dcmplx(xp,yp)
               dsdth(index) = dsqrt((xdot)**2 + (ydot)**2)
               rkappa(istart+i) = (xdot*yddot-ydot*xddot)/
     *                            dsdth(index)**3
            end do
ccc            call RSCPLOT (z(1,kbod),nd,1,iplot)
            istart = istart+nd(kbod)
         end do
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine GEOSHIFT (k, nd, nbk, ai, bi, thetax, ncyc, zk, beta,
     *                     z, rkappa, iplot, nmax)
c
c---------------
c
c  Given the equispaced parametrization imbedded in beta, compute
c  the geometry data at equispaced points in arclength
c
      implicit real*8 (a-h, o-z)

      dimension ai(k), bi(k), thetax(k), ncyc(k)
      dimension beta(nbk), rkappa(nbk), nd(k)
      complex*16 zk(k) 
      complex*16 z(nmax,k)
c
         istart = 0
         do kbod = 1, k
            do i = 1, nd(kbod)
               th = beta(istart+i)
               call INITPARAM (th, ai(kbod), bi(kbod), thetax(kbod),
     *                         ncyc(kbod), xp, yp, xdot, ydot,
     *                         xddot, yddot)
c
               dsdt = dsqrt(xdot**2 + ydot**2)
               rkappa(istart+i) = (xdot*yddot-ydot*xddot)/dsdt**3
               z(i,kbod) = zk(kbod) + dcmplx(xp, yp)
            end do
ccc            call RSCPLOT (z(1,kbod), nd(kbod), kbod+1, iplot)
            istart = istart+nd(kbod)
         end do
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine INITPARAM (th, ai, bi, thetax, ncyc, xp, yp, xdot,
     *                      ydot, xddot, yddot)
c---------------
c
c  The initial curve parameterization is given in term of theta
c  This routine returns curve points and derivative data given
c  theta
c
      implicit real*8 (a-h, o-z)
c
         if (ncyc.eq.0) then
            cs = dcos(th-thetax)
            sn = dsin(th-thetax)
            a2 = (ai)**2
            b2 = (bi)**2
            rnorm = dsqrt(b2*cs**2 + a2*sn**2)
            radius = ai*bi/rnorm
            rdot = -ai*bi*cs*sn*(-b2+a2)/rnorm**3
            rddot =  -ai*bi*(2.d0*a2*b2*cs**2*sn**2
     *                     +a2**2*cs**4 + a2*b2-a2**2+b2**2*cs**4
     *                     -2.d0*b2**2*cs**2)/rnorm**5
            xp = radius*dcos(th)
            yp = radius*dsin(th)
c
            cs = dcos(th)
            sn = dsin(th)
            xdot = rdot*cs - radius*sn
            ydot = rdot*sn + radius*cs
            xddot = rddot*cs - 2.d0*rdot*sn - radius*cs
            yddot = rddot*sn + 2.d0*rdot*cs - radius*sn
           elseif (ncyc.lt.0) then
            R = ai
            anu = bi
            den = (1.d0-2.d0*anu*dcos(2.d0*th)+anu**2)
     *           *dsqrt(1.d0+anu**2)
            ddot = 4.d0*anu*dsin(2.d0*th)*dsqrt(1.d0+anu**2)
            dddot = 8.d0*anu*dcos(2.d0*th)*dsqrt(1.d0+anu**2)
            xp = (1.d0-anu**2)*(1.d0-anu)*R*dsqrt(2.d0)*dcos(th)/den
            yp = (1.d0-anu**2)*(1.d0+anu)*R*dsqrt(2.d0)*dsin(th)/den
            xdot = -(1.d0-anu**2)*(1.d0-anu)*R*dsqrt(2.d0)*dsin(th)/den
     *             -xp*ddot/den            
            xddot = -xp 
     *              +(1.d0-anu**2)*(1.d0-anu)*R*dsqrt(2.d0)
     *                            *dsin(th)*ddot/den**2
     *              - xdot*ddot/den - xp*dddot/den 
     *              + xp*ddot**2/den**2
            ydot = (1.d0-anu**2)*(1.d0+anu)*R*dsqrt(2.d0)*dcos(th)/den
     *             -yp*ddot/den
            yddot = -yp 
     *              -(1.d0-anu**2)*(1.d0+anu)*R*dsqrt(2.d0)
     *                            *dcos(th)*ddot/den**2
     *              - ydot*ddot/den - yp*dddot/den 
     *              + yp*ddot**2/den**2 
           else
            eps = bi
            cs = dcos(th)
            sn = dsin(th)
            snn = dsin(NCYC*th)
            csn = dcos(NCYC*th)
            xp = ai*cs + bi*csn
            yp = ai*sn - bi*snn
            xdot = -ai*sn - NCYC*bi*snn
            xddot = -ai*cs - NCYC**2*bi*csn
            ydot = ai*cs - NCYC*bi*csn
            yddot = -ai*sn + NCYC**2*bi*snn
         end if
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine INITCOND (k, nmax, nd, nbk, neq, z, y)
c
c---------------
c Sets initial conditions for z
c
      implicit real*8 (a-h, o-z)
      complex*16 z(nmax,k)
      dimension y(neq), nd(k)
c
         istart = 0
         nbk = 0
         do kbod = 1, k
            do i = 1, nd(kbod)
               y(istart+2*i-1) = dreal(z(i,kbod))
               y(istart+2*i) = dimag(z(i,kbod))
            end do
            istart = istart + 2*nd(kbod)
            nbk = nbk + nd(kbod)
         end do
         neq = 2*nbk
c
      return
      end        
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine GETFOURIER (nd, amode, fmode, z, za, wsave, idoub)
c---------------
C  Get Fourier coefficients of theta and plot them
c
      implicit real*8 (a-h, o-z)
      dimension amode(*), fmode(*)
      complex*16 za(nd), wsave(*), z(*)
c
         n = nd/2
         pi = 4.d0*datan(1.d0)
c
            do i = 1, nd
               alpha = 2.d0*pi*(i-1)/nd
               za(i) = z(i) 
            end do
            call DCFFTF (nd, za, wsave)
            do i = 1, nd
               za(i) = za(i)/nd
            end do
c                  
            amode(1) = 0.d0
            fmode(1) = cdabs(za(1))
            do i = 1, n-1
               amode(i+1) = float(i)
               fmode(i+1) = cdabs(za(i+1)) + cdabs(za(nd-i+1))
ccc               if (fmode(i+1).le.0.d0) fmode(i+1) = fmode(i)
            end do
c
c  Doubling test
c
         tol = 1.d-11
         idoub = 0
         ns = 3*n/4
         do i = n-10,n-1
            if (fmode(i+1).gt.tol) idoub = 1
         end do  
ccc         call PRINF ('  IDOUB = *', idoub, 1)          
c
      return
      end         
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine EXTRACTSOL (nd, neq, area, arcl, z, dz, rku, za, zai, 
     *                       wsave, errds, rkmax, rkappa)
c
c---------------
c Extracts geometry data from RK solver arrays
c
      implicit real*8 (a-h, o-z)
      dimension rku(neq), rkappa(nd)
      complex*16 z(nd), za(nd), zai(nd), wsave(*), dz(nd)
c
         n = nd/2
         pi = 4.d0*datan(1.d0)      
c
         do i = 1, nd
            z(i) = dcmplx(rku(2*i-1),rku(2*i))
         end do
c
c  Now compute tha from za
         call FDIFFF (z, dz, nd, wsave) 
         smin = 1.d10
         smax = 0.d0
         do i = 1, nd
            xdot = dreal(dz(i))
            ydot = dimag(dz(i))
            dsdt = dsqrt(xdot**2 + ydot**2)
            dz(i) = -dz(i)
            dsda = cdabs(dz(i))
            smin = min(smin, dsda)
            smax = max(smax, dsda)
         end do
         call FDIFFF (dz, za, nd, wsave)
c
c  Compute curvature and area
c
         rkmax = 0.d0
         do i = 1, nd
            xdot = dreal(dz(i))
            ydot = dimag(dz(i))
            xddot = dreal(za(i))
            yddot = dimag(za(i))
            dsdt = dsqrt(xdot**2 + ydot**2)
            rkappa(i) = dabs((xdot*yddot-ydot*xddot)/dsdt**3)
            rkmax = max(rkmax,dabs(rkappa(i)))
            za(i) = dimag(z(i))*dreal(dz(i))
            zai(i) = cdabs(dz(i))
         end do
         call DCFFTF (nd, za, wsave)
         call DCFFTF (nd, zai, wsave)
         area = cdabs(2.d0*pi*za(1)/nd)
         arcl = cdabs(2.d0*pi*zai(1)/nd)
ccc         call PRIN2 ('  dsda min = *', smin, 1)
ccc         call PRIn2 ('  dsda max = *', smax, 1)
ccc         call PRIN2 ('  KAPPA MAX = *', rkmax, 1)
         errds = (smax-smin)/smax
c
      return
      end        
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine EXTRACTDAT (k, nd, nbk, x, y, z, zk, dz, tha, dsda, 
     *                       rkappa, ckern, rku, za, zai, wsave, 
     *                       ifmax, zc, xleft, slength, zbond, icut) 
c
c---------------
c Extracts geometry data from RK solver arrays
c
      implicit real*8 (a-h, o-z)
      dimension nd(k)
      dimension tha(*), rkappa(*), x(*), y(*), rku(*), dsda(*)
      complex*16 z(*), dz(*), ckern(*), zk(k), eye, zc(k)
      complex*16 za(*), zai(*), wsave(ifmax,k), zbond
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0,1.d0) 
c              
         do i = 1, nbk
            z(i) = dcmplx(rku(2*i-1),rku(2*i))
         end do
c
         istart = 0
         xmin = 1.d10
         xmax = -1.d10
         ymin = 1.d10
         ymax = -1.d10
         do kbod = 1, k
            call DCFFTI (nd(kbod), wsave(1,kbod))
            ist = istart+1
            do i = 1, nd(kbod)
               x(istart+i) = dreal(z(istart+i))
               y(istart+i) = dimag(z(istart+i))
               xmin = min(xmin, x(istart+i))
               xmax = max(xmax, x(istart+i))
               ymin = min(ymin, y(istart+i))
               ymax = max(ymax, y(istart+i))
               za(i) = z(istart+i)
            end do
c
c  Calculate geometry variables
c         
            call DCFFTF (nd(kbod), za(1), wsave(1,kbod))
            do i = 1, nd(kbod)
               za(i) = za(i)/nd(kbod)
            end do
            zk(kbod) = za(1)
            call FDIFFF (z(istart+1), za, nd(kbod), wsave(1,kbod)) 
            call FDIFFF (za, zai, nd(kbod), wsave(1,kbod)) 
            dmin = 1.d10
            dmax = 0.d0
            err = 0.d0
            do i = 1, nd(kbod)
               xdot = dreal(za(i))
               ydot = dimag(za(i))
               xddot = dreal(zai(i))
               yddot = dimag(zai(i))
               dsdt = dsqrt(xdot**2 + ydot**2)
               rkappa(istart+i) = -(xdot*yddot-ydot*xddot)/dsdt**3
               ckern(istart+i) = -dcmplx(xddot,yddot)
     *                            /dcmplx(xdot,ydot)
               tha(istart+i) = -rkappa(istart+i)*dsdt
               dz(istart+i) = -dcmplx(xdot, ydot)
               dsda(istart+i) = cdabs(dz(istart+i))
               za(i) = dimag(z(istart+i))*dreal(dz(istart+i))
            end do
            call DCFFTF (nd(kbod), za, wsave(1,kbod))
            area = 2.d0*pi*dreal(za(1)/nd(kbod))
ccc            call PRIN2 ('  area = *', area, 1)
            zc(kbod) = eye*zbond*area/(8.d0*pi)
ccc            call PRIN2 ('   arg = *', arg(istart+1), nd(kbod))
            istart = istart+nd(kbod)
         end do
ccc         call PRIN2 (' zk = *', zk, 2*k)
ccc         call PRIN2 (' zc = *', zc, 2*k)
ccc         call PRIN2 (' zbond = *', zbond, 2)
         xleft = xmin
         slength = max(xmax-xmin, ymax)
ccc         call PRIN2 ('  xleft = *', xleft, 1)
ccc         call PRIN2 ('  slength = *', slength, 1)
c
      return
      end        
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine DOUBLE (nd, z, z2, za, wsave, zwork, nzwork, ndmax, 
     *                   istop)
c
c---------------
c Extracts geometry data from RK solver arrays
c
      implicit real*8 (a-h, o-z)
      complex*16 z(2*nd), z2(2*nd)
      complex*16 za(2*nd), wsave(*), zwork(nzwork)
c
         pi = 4.d0*datan(1.d0) 
c
c  Double the number of points, i.e. pad the spectrum for z
c
         ND2 = 2*nd
         NDM1= nd-1
         ND2M1= ND2-1
         CALL FTRPINC(zwork,nzwork,IP2,IPT,NDM1,ND2M1)
         ndo = nd
         nd = 2*nd
         istop = 0
         if (nd.gt.ndmax) then
            write (6,*) 'NEW ND TOO BIG IN DOUBLE'
            write (6,*) 'Increase NDMAX, ndmax = ', ndmax
            istop = 1
           else
            n = nd/2
c
            call DCFFTI (nd, wsave)
            CALL FINTERC(z,za,NDM1,ND2M1,zwork,nzwork,IP2,IPT)
            do i = 1, nd
               z(i) = za(i)
            end do
         end if
c
      return
      end        
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine RK2 (k, neq, time, dt, phi, ybar, y, dy)
c
c---------------
c Do Runge-Kutta step (Modified Euler)
c
      implicit real*8 (a-h, o-z)
      dimension ybar(neq), y(neq), dy(neq), phi(neq,0:3)
c
c  Compute f0, f1 and put in phi array
c
         t = time
         call F (k, neq, t, ybar, dy)
         do i = 1, neq
            phi(i,0) = dt*dy(i)
            y(i) = ybar(i) + phi(i,0)
         end do
         t = time+dt
         call F (k, neq, t, y, dy)
         do i = 1, neq
            phi(i,1) = dt*dy(i)
            ybar(i) = ybar(i) + 0.5d0*phi(i,0) + 0.5d0*phi(i,1)
         end do
c
      return
      end 
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine INTERFACENEW (nd, k, nbk, h, x, y, dz, u, Un, rkappa, 
     *                         ckern, wsave, zk, ifwork, zc, iwall,
     *                         xleft, slength, cap, dsda, zb) 
c
c---------------
      implicit real*8 (a-h,o-z)
      dimension nd(k), Un(nbk), rkappa(nbk), h(k), x(nbk), y(nbk),
     *          dsda(nbk)
      complex*16 u(nbk), ckern(nbk), zk(k)
      complex*16 zn, zrhs, eye, dz(nbk), zc(k),
     *           wsave(ifwork,k), zdiff, zck, zdis, z
c
      parameter (nsp = 500000, ncwork = 500000, npmax = 500000)
      dimension ier(5), inform(6), iwork(ncwork)
      complex*16 zwork(nsp), udot(npmax), gradw(npmax)
      complex*16 zvel, zpsi, zcorr, zself, arg(npmax), arga(npmax),
     *           zkern, gradwr, zphi, zphip, z2pii, zb(k)
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0,1.d0)
         z2pii = 1.d0/(2.d0*pi*eye)
c
c  Before first call to PVINTEV, first interpolate x and y to twice
c  the number of points
c
         istart = 0
         do nbod = 1,k
            h(nbod) = 2.d0*pi/nd(nbod)
            zb(nbod) = 0.d0
            call FDIFFF (u(istart+1), udot(istart+1), nd(nbod), 
     *                   wsave(1,nbod))
            istart = istart+nd(nbod)
         end do
c         
         LEVMAX = 7
         TOL = 1.D-13
         NTARG = 0
         iflag1 = 0
         if (iwall.eq.1) then
            iflag2 = 1
           else
            iflag2 = 0
         end if
         call PRINI (0,0)
         call CADAPMOD (iflag1,iflag2,levmax,x,y,dz,1,k,nd,u,h,
     1                 iwork,ncwork,zwork,nsp,xtarg,ytarg,ntarg,gradw,
     2                 phipr,phippr,tol,inform,ier,xleft,slength)
         call PRINI (6,13)
         if (ier(1).ne.0) then
            call PRINF ('  ERROR IN CADAPMOD, IER = *', ier, 5)
            call PRINF ('  INFORM = *', inform, 6)
            stop
         end if
c
c  Add on corrections
c
         istart = 0
         flux = 0.d0
         do kbod = 1, k
            umax = 0.d0
            do i = 1,nd(kbod)
               zn = -eye*dz(istart+i)/dsda(istart+i)
               z = dcmplx(x(istart+i), y(istart+i))
               zcorr = 0.5d0*h(kbod)*2.d0*dreal(ckern(istart+i))
     *                 *u(istart+i) - 2.d0*h(kbod)*udot(istart+i)
               zcorr = zcorr*z2pii
               zself = dconjg(u(istart+i))*h(kbod)*
     *                 cdabs(dz(istart+i))*0.5d0*rkappa(istart+i)
     *                 *zn*zn/pi
               gradw(istart+i) = gradw(istart+i) + zself + zcorr
               zck = 0.d0
               do nbod = 1, k
                  zdis = dcmplx(x(istart+i),y(istart+i)) - zk(nbod)
                  zck = zck + 2.d0*zc(nbod)*dlog(cdabs(zdis)) 
     *            + dconjg(zc(nbod))*zdis/dconjg(zdis)
               end do
               zphi = 0.d0
               zphip = 0.d0
               zpsi = -eye*0.5d0*cap*z
               if (iwall.eq.1) 
     *            call REFSING (z,k,zk,zc,zb,zphi,zphip,zpsi)
               gradwr = zphi + z*dconjg(zphip) + dconjg(zpsi)
               gradw(istart+i) = gradw(istart+i) + gradwr + zck
               zvel = -eye*( gradw(istart+i) )
               Un(istart+i) = dreal(zvel*dconjg(zn))
               flux = flux + h(kbod)*Un(istart+i)*dsda(istart+i)
               umax = max(umax, dabs(Un(istart+i)))
            end do
ccc            call PRIN2 ('  umax = *', umax, 1)
            istart = istart+nd(kbod)
         end do
         call PRIN2 ('  FLUX = *', flux, 1)
c
c  Output data to be read in by Stokes solver
c
ccc         open (unit=59, file='stokes.data')
ccc         write (59,*) k, (nd(kk),kk=1,k), nbk
ccc         write (59, *) iwall
ccc         write (59,'(e22.16,2x,e22.16)') xleft, slength
ccc         do i = 1,k
ccc            write (59,'(e22.16,2x,e22.16)') zk(i)
ccc         end do
ccc         do i = 1,k
ccc            write (59,'(e22.16,2x,e22.16)') zc(i)
ccc         end do
ccc         do i = 1, nbk
ccc            write (59,'(e22.16,2x,e22.16)') x(i), y(i)
ccc         end do
ccc         do i = 1, nbk
ccc            write (59,'(e22.16,2x,e22.16)') dz(i)
ccc         end do
ccc         do i = 1, nbk
ccc            write (59,'(e22.16,2x,e22.16)') rkappa(i)
ccc         end do
ccc         do i = 1, nbk
ccc            write (59,'(e22.16,2x,e22.16)') gradw(i)
ccc         end do 
ccc         close(59)  
ccc         stop               
c
      return
      end 
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine WALLCHECK (nd, k, nbk, h, x, y, dz, u) 
c
c---------------
      implicit real*8 (a-h,o-z)
      dimension nd(k), h(k), x(nbk), y(nbk)
      complex*16 u(nbk), dz(nbk)
c
      parameter (nsp = 5000, ncwork = 5000, npmax = 5000)
      dimension ier(5), inform(6), iwork(ncwork)
      dimension xtarg(npmax), ytarg(npmax)
      complex*16 zwork(nsp), gradw(npmax)
      complex*16 phipr(npmax), phippr(npmax)
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0,1.d0)
c
         do kbod = 1, k
            h(kbod) = 2.d0*pi/nd(kbod)
         end do
         ntarg = 15
         dx = 4.d0/(ntarg)
         do i = 1, ntarg
            ytarg(i) = 0.d0
            xtarg(i) = -2.d0 + (i-0.5d0)*dx
         end do 
c         
         LEVMAX = 4
         TOL = 1.D-14
         iflag1 = 1
         iflag2 = 1
         XLEFT = -2.d0
         SLENGTH = 4.d0
         call PRINI (0,0)
         call CADAPMOD (iflag1, iflag2, levmax, x, y, dz, 1, k, nd, u,
     *                  h, iwork, ncwork, zwork, nsp, xtarg, ytarg, 
     *                  ntarg, gradw, phipr, phippr, tol, inform, ier, 
     *                  xleft, slength)         
         call PRINI (6,13)
         if (ier(1).ne.0) then
            call PRINF ('  ERROR IN CADAPMOD, IER = *', ier, 5)
            call PRINF ('  INFORM = *', inform, 6)
            stop
         end if
         call PRINF (' ntarg = *', ntarg, 1)
         call PRIN2 ('  gradw on wall = *', gradw, 30)
c
      return
      end 
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine F (kk, neqm,time,rku,rkup)
c---------------
c
      implicit real*8 (a-h,o-z)
      dimension rku(*), rkup(*)
c
      parameter (ndmax=3000, kmax = 150, nmax = ndmax*kmax)
c
c  Geometry arrays
c
      dimension rkappa(nmax), x(nmax), y(nmax), dsda(nmax), arg(nmax)
      dimension nd(kmax), h(kmax)
      complex*16 zki, zko, zk(kmax), z(nmax), dz(nmax), ckern(nmax) 
c
c  Theta-L arrays
c
      dimension theta(nmax), tha(nmax), Un(nmax), T(nmax)
c
c Laplace arrays
      dimension ak(kmax), ulap(nmax+kmax)
      dimension phi(nmax), psi(nmax), dphi_dn(nmax)
c
c  Velocity and integral equation arrays
c
      dimension pk(kmax)
      complex*16 u(nmax), zvel(nmax), za(kmax), zb(kmax), zc(kmax)
c
c  FFT and Work arrays
c
      parameter (ifwork = 10*ndmax+50, nzwork=16*nmax+100)
      parameter (nbox=2**12, niwork=5*nbox+3*nmax+10)
      dimension iwork(niwork)
      complex*16 zf(nmax), zfi(nmax), wsave(ifwork,kmax), eye, zbond
      dimension amode(nmax), fmode(nmax)
c
c  Matrix equation variables for GMRES
c
c  MAXL is the maximum nubmer of GMRES iterations performed
c       before restarting.
c  LRWORK is the dimension of a real workspace needed by DGMRES.
c  LIWORK is the dimension of an integer workspace needed by DGMRES.
c  GMWORK and IWORK are work arrays used by DGMRES
c
      parameter (maxl = 100, ncmax = 2*nmax+kmax)
      parameter (lrwork=10+ncmax*(maxl+6)+maxl*(maxl+3), liwork=30)
      dimension gmwork(lrwork), igwork(liwork)
      dimension rhs(ncmax), soln(ncmax)
c
c  Force and Torque arrays
c
      dimension xforce(kmax), yforce(kmax), torque(kmax), rw(100000)
      dimension x2(2*nmax), y2(2*nmax)
      complex*16 zw(500000), z2(500000), w(500000)     
c
      common /problem/ zk, nd, k, nbk, neq
      common /inteqn/ dsda, rkappa, z, dz
      common /density/ u
      common /flowinfo/ cap, zbond, E0
      common /wallinfo / xleft, slength, iwall         
      common /doublend/ icut       
      REAL*4 TIMEP(2), ETIME
c
         pi = 4.d0*datan(1.d0)  
         eye = dcmplx(0.d0,1.d0)    
         call PRINI (6, 13)
         k = kk
         nbk = neqm/2

c
c  Solve the integral equation and compute the velocity on the boundary
c          
         call EXTRACTDAT (k, nd, nbk, x, y, z, zk, dz, tha, dsda, 
     *                    rkappa, ckern, rku, zf, zfi, wsave, ifwork, 
     *                    zc, xleft, slength, zbond, icut)
c
c Solve Laplace's equation 1st
         call GETBCs (k, nd, nbk, zk, z, rhs)
ccc         call PRIN2 (' rhs  = *', rhs,nbk)
         igwork(1) = maxl
         call SOLVE_LAPL (nd, k, nbk, rhs, soln, ulap, ak, 
     1                    gmwork, lrwork, igwork, liwork, dsda)
         call FASMVP_LAP (k, nd, nbk, zk, z, dz, rkappa, dsda, 
     1                       ak, soln, rhs)
ccc         call prin2 (' residual after solve = *', rhs, nbk)
ccc       stop
         call DIRICH_NEU_MAP (k, nd, nbk, z, ulap, dz, zk, ak, phi, 
     1                        psi, dphi_dn, dsda)
c
c Stokes Solve     
         call SETVELF (k, nd, nbk, z, zk, dz, arg, zb, zc, zf, zfi, 
     *                 wsave, ifwork, rhs, iwall, zbond, cap, dsda, 
     2                 dphi_dn, E0)
         igwork(1) = maxl
         call SOLVEI (nd, k, nbk, rhs, soln, pk, pinf, za, zk, z, dz, 
     *                u, gmwork, lrwork, igwork, liwork)
         call PRIN2 ('  pk = *', pk, k)
         call PRIN2 ('  pinf = *', pinf, 1)
         call INTERFACENEW (nd, k, nbk, h, x, y, dz, u, Un, rkappa, 
     *                      ckern, wsave, zk, ifwork, zc, iwall, xleft,
     *                      slength, cap, dsda, zb) 
ccc            call WALLCHECK (nd, k, nbk, h, x, y, dz, u) 
ccc         call PRIn2 ('  Un 1 = *', Un, nd(1))
ccc         call PRIn2 ('  Un 2 = *', Un(1+nd(1)), nd(2))
ccc         stop
c
c  Now, compute the time derivative for the fourier coefficients of
c  theta_alpha
c      
         istart = 0
         istarto = 0
         do kbod = 1, k
            n = nd(kbod)/2
            ist = istart+1
            call EQUIT (nd(kbod), tha(ist), Un(ist), T, zf, zfi, 
     *                  wsave(1,kbod))
            do i = 1, nd(kbod)
               zf(i) = -Un(istart+i)*eye*dz(istart+i) 
     *               - T(i)*dz(istart+i)
               zf(i) = zf(i)/dsda(istart+i)
            end do
            call PRINI (6,13)
c         
            call DCFFTF (nd(kbod), zf, wsave(1,kbod))
            do i = 1, nd(kbod)
               zf(i) = zf(i)/nd(kbod)
            end do
ccc            call PRIN2 ('  fft of rhs = *', zf, 2*nd(kbod))
            call FILTER (nd(kbod), zf)
            call DCFFTB (nd(kbod), zf, wsave(1,kbod))
            do i = 1, nd(kbod)
               rkup(istarto+2*i-1) = dreal(zf(i))
               rkup(istarto+2*i) = dimag(zf(i))
            end do
            istart = istart+nd(kbod)
            istarto = istarto+2*nd(kbod)
         end do
ccc         tend = etime(timep)
ccc         call PRIN2 ('  time in F = *', tend-tbeg, 1)
ccc         stop
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine PCALC (nd, k, nbk, z, dz, u, pinf)
c
c---------------
c
      implicit real*8 (a-h, o-z)
      dimension nd(k)
      complex*16 z(nbk), dz(nbk), u(nbk), zphip, om, eye, zp, zdis
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0,1.d0)
c
   1     write (6,*) 'ENTER Z '
          read (5,*) x,y
         zp = dcmplx(x,y)
         zphip = 0.d0
         istart = 0
         do kbod = 1, k
            h = 2.d0*pi/nd(kbod)
            do i = 1, nd(kbod)   
               om = -eye*u(istart+i)
               zdis = z(istart+i) - zp
               zphip = zphip + h*dz(istart+i)*om/zdis**2
            end do
            istart = istart + nd(kbod)
         end do 
         zphip = zphip/(2.d0*pi*eye)
         pressure = -4.d0*dimag(zphip) + pinf
         call PRInd2 ('  pressure = *', pressure, 2)
         write (6,*) 'continue, yes=1'
          read (5,*) iresp
          if (iresp.eq.1) goto 1
c
      return
      end                             
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine FILTER (nd, cn)
c
c---------------
c
      implicit real*8 (a-h, o-z)
      complex*16 cn(nd)
         n = nd/2
         eps = 1.d-13
         iord = 25
c
         do i = 0, n
            ratio = float(i)/float(n)
            ffilt = dexp(-10.d0*ratio**iord)
            cn(i+1) = ffilt*cn(i+1)
            cr = dreal(cn(i+1))
            ci = dimag(cn(i+1))
            if (dabs(dreal(cn(i+1))).lt.eps) 
     *         cn(i+1) = dcmplx(0.d0, ci)
            if (dabs(dimag(cn(i+1))).lt.eps) 
     *         cn(i+1) = dcmplx(dreal(cn(i+1)), 0.d0)
         end do
         cn(n+1) = 0.d0 
         do i = 1, n-1
            ratio = float(i)/float(n)
            ffilt = dexp(-10.d0*ratio**iord)
            cn(nd-i+1) = ffilt*cn(nd-i+1)
            cr = dreal(cn(nd-i+1))
            ci = dimag(cn(nd-i+1))
            if (dabs(dreal(cn(nd-i+1))).lt.eps) 
     *         cn(nd-i+1) = dcmplx(0.d0, ci)
            if (dabs(dimag(cn(nd-i+1))).lt.eps) 
     *         cn(nd-i+1) = dcmplx(dreal(cn(nd-i+1)), 0.d0)
         end do
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine CLUSTPARAM (k, nd, nbk, dsdth, arcl, beta, 
     *                       csa, csai, wsave, ifwork)
c
c---------------
      implicit real*8 (a-h, o-z)
      dimension dsdth(nbk), beta(nbk), arcl(k), nd(k)
      complex*16 csa(*), csai(*), wsave(ifwork,k), eye, FOUREVAL
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0, 1.d0)
c
         istart = 0
         do kbod = 1, k
            call DCFFTI (nd(kbod), wsave(1,kbod))
            n = nd(kbod)/2
            do i = 1, nd(kbod)
               csa(i) = dsdth(istart+i)
            end do
            call DCFFTF (nd(kbod), csa, wsave(1,kbod))
            do i = 1, nd(kbod)
               csa(i) = csa(i)/nd(kbod)
            end do
            csa(n+1) = 0.d0
            arcl(kbod) = 2.d0*pi*dreal(csa(1))
c
c  Compute fourier coefficients of integral
c
            call FOURINT (nd(kbod), csa, csai)
c
c  Determine constant of integration so that first point on
c  closed curve corresponds to s = 0
c
            csai(1) = 0.d0
            fk = dreal(FOUREVAL (nd(kbod), csai, 0.d0))
            csai(1) = -fk
ccc            call PRINd2 (' csai(1) after fourint = *', csai(1), 2)
ccc            call PRINd2 ('  CKI = *', csai, 2*nd)
c
c  Now use Newton's method to compute mapping to equi-arclength
c  distribution
c
            tol = 1.d-14
            beta0 = 0.d0
            do i = 1, nd(kbod)
               th = 2.d0*pi*(i-1)/nd(kbod)
               do niter = 1, 100
                  fk = dreal(FOUREVAL (nd(kbod), csai, beta0))
                  fk = fk + csa(1)*beta0 
     *                - (i-1.d0)*arcl(kbod)/nd(kbod)
                  fpk = dreal(FOUREVAL (nd(kbod), csa, beta0))
                  del = -fk/fpk
                  beta0 = beta0 + del
                  iter = niter
                  if (dabs(del).lt.tol) goto 100
               end do
               call PRINF('  NEWTON FAILED TO CONVERGE *',1,1)
               stop
 100           continue
               beta(istart+i) = beta0
            end do
ccc            call PRIN2 ('  BETA IN CLUST = *', beta(istart+1), nd)
c
c  Compute arc length
c            
            do i = 1, nd(kbod)
               csa(i) = dsdth(istart+i)
            end do
            call DCFFTF (nd(kbod), csa, wsave(1,kbod))
            call PRIN2 ('  INITIAL ARC LENGTH = *', arcl(kbod), 1)
            istart = istart + nd(kbod)
         end do
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine EQUIT (nd, tha, Un, T, zf, zfi, wsave)
c
c---------------
c
c  Computes T(alpha, t) for equispaced parametrization
c
      implicit real*8 (a-h, o-z)
      dimension tha(nd), Un(nd), T(nd)
      complex*16 zf(nd), zfi(nd), wsave(*)
c      
         pi = 4.d0*datan(1.d0)
         dalph = 2.d0*pi/nd
         n = nd/2
c
         do i = 1, nd
            zf(i) = - tha(i)*Un(i)
         end do
         call DCFFTF (nd, zf, wsave)
         do i = 1, nd
            zf(i) = zf(i)/nd
         end do
         zf(n+1) = 0.d0
         dadt = dreal(zf(1))
         call FOURINT (nd, zf, zfi)
         call DCFFTB (nd, zfi, wsave)
         do i = 1, nd
            alpha = dalph*(i-1.d0)
            T(i) = -dreal(zf(1)*alpha + zfi(i))
     *             + alpha*dadt
         end do
c         
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine FOURINT (nd, zf, zfi)
c
c---------------
c  Given the fourier coefficients of a function zf, compute the
c  fourier coefficients of it's integral and store in zfi
c  NOTE:  a constant term must be added in later
c
      implicit real*8 (a-h, o-z)
      complex*16 zf(nd), zfi(nd), eye
c
         eye = dcmplx(0.d0, 1.d0)
         pi = 4.d0*datan(1.d0)
c
         n = nd/2
         zfi(1) = 0.d0
         do i = 1, n-1
            zfi(i+1) = zf(i+1)/(i*eye)
         end do
         zfi(n+1) = 0.d0
         do i = 1, n-1
            zfi(nd-i+1) = -zf(nd-i+1)/(i*eye)
         end do
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      complex*16 function FOUREVAL (nd, ck, theta)
c
c--------------
      implicit real*8 (a-h, o-z)
      complex*16 ck(nd), zsum, eye
c
         eye = dcmplx(0.d0,1.d0)
c
         n = nd/2
         zsum = ck(1)
         do i = 1, n-1
            zsum = zsum + ck(i+1)*cdexp(eye*i*theta)
         end do
         do i = 1, n-1
            zsum = zsum + ck(nd-i+1)*cdexp(-i*eye*theta)
         end do
         FOUREVAL = zsum
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine DUMPDATA (k, nd, nbk, neq, area0, arcl0, z, nmax,
     *                     dt, time)
c
c---------------
c
      implicit real*8 (a-h, o-z)
      dimension area0(k), arcl0(k), nd(k)
      complex*16 z(nmax,k)
c
         open (unit=56, file = 'geometry.data')
         write (56, *) k, nbk, neq
         write (56, *) (nd(kk),kk=1,k)
         write (56,'(e22.16,2x,e22.16)') dt, time
         do kbod = 1, k
            write (56,'(e22.16,2x,e22.16)') area0(kbod)
         end do
         do kbod = 1, k
            write (56,'(e22.16,2x,e22.16)') arcl0(kbod)
         end do
         do kbod = 1,k
            do i = 1, nd(kbod)
               write (56,'(e22.16,2x,e22.16)') z(i,kbod)
            end do
         end do
         close (56)
c
      return
      end         
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine READDATA (k, nd, nbk, neq, area0, arcl0, z, wsave, 
     *                     nmax, ifwork, dt, tstart)
c
c---------------
c
      implicit real*8 (a-h, o-z)
      dimension area0(k), arcl0(k), nd(k)
      complex*16 z(nmax,k), wsave(ifwork,k)
c
         open (unit=56, file = 'geometry.data')
         read (56, *) k, nbk, neq
         read (56, *) (nd(kk),kk=1,k)
         read (56, *) dt, tstart
         do kbod = 1, k
            read (56,*) area0(kbod)
            call DCFFTI (nd(kbod),wsave(1,kbod))
         end do
         do kbod = 1, k
            read (56,*) arcl0(kbod)
         end do
         do kbod = 1, k
            do i = 1, nd(kbod)
               read (56,*) x, y
               z(i,kbod) = dcmplx(x,y)
            end do
         end do
         close (56)
c
      return
      end         
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine CALCTORQUE (k, nd, nbk, z, zk, dz, rhs, zf, zfi, 
     *                       wsave, za)
c
c---------------
c
      implicit real*8 (a-h, o-z)
      dimension rhs(2*nbk), nd(k)
      complex*16 z(nbk), dz(nbk), zk(k), eye, zstress, ztorque, zdis
      complex*16 zf(*), zfi(*), wsave(*), zforce, za(k)
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0,1.d0)
         h = 2.d0*pi/nd(1)
         call PRINF ('  got in calctorque *', 1, 1)
c
         istart = 0
         istartc = 0
         do kbod = 1, k
            ztorque = 0.d0
            do i = 1, nd(kbod)
               zstress = dcmplx(rhs(istartc+2*i-1),rhs(istartc+2*i))
               zf(i) = zstress
            end do
            call PRIN2 ('  zstress in calctorque = *', zf, 2*nd(1))
            call DCFFTF (nd(1), zf, wsave)
            do i = 1, nd(1)
               zf(i) = zf(i)/nd(1)
            end do
            call PRIN2 ('  fft of zstress = *', zf, 2*nd(1))
            stop
            call FDIFFF (zf, zfi, nd, wsave)
            do i = 1, nd(kbod)
               zdis = z(istart+i) - zk(kbod)
               ds = cdabs(dz(istart+i))
               zforce = -2.d0*zfi(i)
ccc               call PRIN2 ('  zforce in calctorque = *', zforce, 2)
               ztorque = ztorque - zdis*dconjg(zforce)*h
            end do
            torque = dimag(ztorque)
            call PRIN2 ('  torque = *', torque, 1)
            istart = istart + nd(kbod)
            istartc = istartc + 2*nd(kbod)
         end do
c
      return
      end         
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine CHECKDENS (k, nd, nbk, h, z, dz, dsdth, rkappa, u,
     *                      za, zb, zk, pinf)
c
c---------------
c
      implicit real*8 (a-h, o-z)
      parameter (maxwrk = 900000, maxiwrk = 400000)
      dimension x(maxwrk), y(maxwrk), dsdth(nbk), rkappa(nbk), h(k)
      dimension  nd(k)
      complex*16 dz(nbk), zko, zki, u(nbk)
      complex*16 stress(nbk) 
      complex*16 z(nbk), zu, zn, zx, zself, eye, zb0, zz, zw
      complex*16 zf, zfp, zg, zc, zd, zpress, zdis
      complex*16 za(k), zb(k), zk(k)
c
c  WORK ARRAYS
c
      dimension ier(5),inform(6),iwork(maxiwrk), yy(100000)
      complex*16 cwork(maxwrk)
c
c  Dummy variables
c
      dimension XTARG(1), YTARG(1)
      complex*16 fp(1), fpp(1)
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0,1.d0)  
c
         do kbod = 1, k
            h(kbod) = 2.d0*pi/nd(kbod)
         end do
         do i = 1, nbk
            x(i) = dreal(z(i))
            y(i) = dimag(z(i))
         end do
         call PRIN2 ('pinf in checkedens = *', pinf, 1)
c         
         iflag1 = 0
         iflag2 = 0
         levmax = 5
         lenw = maxwrk
         leniw = maxiwrk
         tol = 1.d-12
         call PRINI (0,0)
         call CADAP(iflag1,iflag2,levmax,x,y,dz,1,k,nd,u, h,
     1                  iwork,leniw,cwork,lenw,xtarg,ytarg,ntarg,
     2                  stress,fp,fpp,tol,inform,ier,xleft,slength)
         call PRINI (0,13)
         if (ier(1).ne.0) then
            write (6,*) '  ****   ERROR IN CADAP   ****'
            call PRINF (' IER = *',IER,5)
            call PRINF (' INFORM = *',INFORM,6)
            write (6,*) '  IER(1) = ',ier(1)
            write (6,*) '  INFORM = ',(inform(ii),ii=1,6)
            stop
         end if
c
c  Add on the identity, the self interactions, and the singular 
c  contributions
c
         istart = 0
         istartc = 0
         do kbod = 1,k
            do i = 1,nd(kbod)
               zu = u(istart+i)
               zn = -eye*dz(istart+i)/cdabs(dz(istart+i))
               qwi = h(kbod)
c               
               zf = 0.5d0*zu - za(kbod)
               zfp = 0.d0 
               zg = 0.5d0*dconjg(zu)
               do nbod = 1,k
                  zg = zg + zb(nbod)/(z(istart+i) - zk(nbod))
               end do
               zself = zu*qwi*dsdth(istart+i)
     *                   *0.5d0*rkappa(istart+i)/pi
               zself = zself + dconjg(zu)*qwi*dsdth(istart+i)*0.5d0
     *                          *rkappa(istart+i)*zn*zn/pi
               stress(istart+i) = stress(istart+i) + zself 
     *                    + zf + z(i)*dconjg(zfp) + dconjg(zg)
     *                    + 0.5d0*pinf*z(istart+i)
            end do
            istart = istart+nd(kbod)
            istartc = istartc + 2*nd(kbod)
         end do
         call PRIN2 ('stress on 1 = *', stress, 2*nd(1))
         call PRIN2 ('stress on 2 = *', stress(1+nd(1)), 2*nd(1))
         stop
c
      return
      end                                  
