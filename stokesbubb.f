c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine SETVELF (k, nd, nbk, z, zk, dz, arg, zb, zc, zf, zfi, 
     *                    wsave, ifwork, rhs, iwall, zbond, cap, dsda,
     2                    dphi_dn, E0)
c---------------
c
c  The right hand side is given by integrating components of the stress
c
c     i PHI + Z CONJG(i PHI') + CONJG(i PSI) =  -sigma * n
c
c  where sigma = surface tension
c
      implicit real*8 (a-h,o-z)
      dimension rhs(*), nd(k), arg(nbk), dsda(nbk), dphi_dn(nbk)
      complex*16 z(nbk), dz(nbk), zc(k), zf(nbk), zfi(nbk), 
     *           wsave(ifwork,k), zk(k), zb(k), zf1
      complex*16 zn, eye, zvel, zgrav, zck, zrhs, zdis, ck
      complex*16 zphi, zphip, zpsi, zckr, zbond, zstr
c
         eye = dcmplx(0.d0,1.d0)
         pi = 4.d0*datan(1.d0)
c
ccc         call PRIN2 ('  zc in setvel = *', zc, 2*k)
         do kbod = 1, k
            zb(kbod) = 0.d0
         end do
c
         istart = 0
         istart2 = 0
         do kbod = 1, k
c
c  Contribution to rhs from bond numbers and reflected C_k
c         
            do i = 1, nd(kbod)
               zfi(i) = zbond*dconjg(z(istart+i))*dz(istart+i) 
     *                  + dconjg(zbond)*z(istart+i)*dz(istart+i)
            end do
            call DCFFTF (nd(kbod), zfi, wsave(1,kbod))
            do i = 1, nd(kbod)
               zfi(i) = zfi(i)/nd(kbod)
            end do
            call FOURINT (nd(kbod), zfi, zf)
            call DCFFTB (nd(kbod), zf, wsave(1,kbod)) 
c
            do i = 1, nd(kbod)
               alpha = 2.d0*pi*(i-1.d0)/nd(kbod)
               zn = -0.5d0*dz(istart+i)/dsda(istart+i)
               zgrav = -eye*0.25d0*(alpha*zfi(1) + zf(i))
               rhs(istart2+2*i-1) = dreal(zn + zgrav)
               rhs(istart2+2*i) = dimag(zn + zgrav)
               if (iwall.eq.1) then
                  zphi = 0.d0
                  zphip = 0.d0
                  zpsi = 0.d0
                  call REFSING (z(istart+i),k,zk,zc,zb,zphi,zphip,zpsi)
                  zckr = zphi-z(istart+i)*dconjg(zphip)-dconjg(zpsi)
                  rhs(istart2+2*i-1) = rhs(istart2+2*i-1) - dreal(zckr)
                  rhs(istart2+2*i) = rhs(istart2+2*i) - dimag(zckr)
                 else
                  zphi = 0.d0
                  zphip = 0.d0
                  zpsi = -eye*0.5d0*cap*z(istart+i)
                  zstr = zphi-z(istart+i)*dconjg(zphip)-dconjg(zpsi)
               end if
               rhs(istart2+2*i-1) = rhs(istart2+2*i-1) - dreal(zstr)
               rhs(istart2+2*i) = rhs(istart2+2*i) - dimag(zstr)
            end do
c            
c  Add on contributions from C_k
c
            do nbod = 1, k
               do i = 1, nd(kbod)
                  zf(i) = (z(istart+i) - zk(nbod))
     *                 /cdabs(z(istart+i)-zk(nbod))
               end do
               call FDIFFF (zf, zfi, nd(kbod), wsave(1,kbod))
               do i = 1, nd(kbod)
                  zf(i) = -eye*zfi(i)/zf(i)
               end do
               call DCFFTF (nd(kbod), zf, wsave(1,kbod))
               do i = 1, nd(kbod)
                  zf(i) = zf(i)/nd(kbod)
               end do
               zf1 =  dreal(zf(1))
               call FOURINT (nd(kbod), zf, zfi)
               call DCFFTB (nd(kbod), zfi, wsave(1,kbod))
               rad = cdabs(z(istart+1)-zk(nbod))
               arg0 = dasin(dimag(z(istart+1)-zk(nbod))/rad)
               if (dreal(z(istart+1)-zk(nbod)).lt.0.d0)
     *            arg0 = pi - arg0  
               do i = 1, nd(kbod)
                  alpha = 2.d0*pi*(i-1.d0)/nd(kbod)
                  arg(i) = arg0 + alpha*zf1 + dreal(zfi(i)-zfi(1))
                  zck = 2.d0*eye*zc(nbod)*arg(i) 
     *                - dconjg(zc(nbod))*cdexp(2.d0*eye*arg(i))
                  zrhs = -zck   
                  rhs(istart2+2*i-1) = rhs(istart2+2*i-1) + dreal(zrhs)
                  rhs(istart2+2*i) = rhs(istart2+2*i) + dimag(zrhs)
               end do
            end do
c               
            istart = istart + nd(kbod)
            istart2 = istart2 + 2*nd(kbod)
         end do
c
c  Contribution to rhs from Electricity
c      
         istart = 0
         istart2 = 0
         do kbod = 1, k
            do i = 1, nd(kbod)
               zfi(i) = dz(istart+i)*(dphi_dn(istart+i))**2
            end do
            call DCFFTI (nd(kbod), wsave(1,kbod))
            call DCFFTF (nd(kbod), zfi, wsave(1,kbod))
            do i = 1, nd(kbod)
               zfi(i) = zfi(i)/nd(kbod)
            end do
            zf1 = zfi(1)
            call FOURINT (nd(kbod), zfi, zf)
            call DCFFTB (nd(kbod), zf, wsave(1,kbod))
            do i = 1, nd(kbod)
               zstr = alpha*zf1 + zf(i)-zf(1)
               zstr = 0.25d0*eye*E0*zstr
               rhs(istart2+2*i-1) = rhs(istart2+2*i-1) + dreal(zstr)
               rhs(istart2+2*i) = rhs(istart2+2*i) + dimag(zstr)
            end do
            istart = istart + nd(kbod)
            istart2 = istart2 + 2*nd(kbod)
         end do
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine SETVELFO (k, nd, nbk, z, zk, dz, arg, zb, zc, zf, zfi, 
     *                    wsave, ifwork, rhs, iwall, zbond, cap, dsda)
c---------------
c
c  The right hand side is given by integrating components of the stress
c
c     i PHI + Z CONJG(i PHI') + CONJG(i PSI) =  -sigma * n
c
c  where sigma = surface tension
c
      implicit real*8 (a-h,o-z)
      dimension rhs(*), nd(k), arg(nbk)
      complex*16 z(nbk), dz(nbk), zc(k), zf(nbk), zfi(nbk), 
     *           wsave(ifwork,k), zk(k), zb(k)
      complex*16 zn, eye, zvel, zgrav, zck, zrhs, zdis, ck
      complex*16 zphi, zphip, zpsi, zckr, zbond
c
         eye = dcmplx(0.d0,1.d0)
         pi = 4.d0*datan(1.d0)
c
ccc         call PRIN2 ('  zc in setvel = *', zc, 2*k)
c
         do i = 1, 2*nbk
            rhs(i) = 0.d0
         end do         
c
         istart = 0
         istart2 = 0
         do kbod = 1, k
            zb(kbod) = 0.d0
c
c  Contribution to rhs from bond numbers
c         
            do nbod = 1, k
               do i = 1, nd(kbod)
                  zf(i) = (z(istart+i) - zk(nbod))
     *                 /cdabs(z(istart+i)-zk(nbod))
               end do
               call FDIFFF (zf, zfi, nd(kbod), wsave(1,kbod))
               do i = 1, nd(kbod)
                  zf(i) = -eye*zfi(i)/zf(i)
               end do
               call DCFFTF (nd(kbod), zf, wsave(1,kbod))
               do i = 1, nd(kbod)
                  zf(i) = zf(i)/nd(kbod)
               end do
               zf1 =  dreal(zf(1))
               call FOURINT (nd(kbod), zf, zfi)
               call DCFFTB (nd(kbod), zfi, wsave(1,kbod))
               rad = cdabs(z(istart+1)-zk(nbod))
               arg0 = dasin(dimag(z(istart+1)-zk(nbod))/rad)
               if (dreal(z(istart+1)-zk(nbod)).lt.0.d0)
     *            arg0 = pi - arg0  
ccc               call PRIN2 ('  arg0 = *', arg0, 1)             
               do i = 1, nd(kbod)
                  alpha = 2.d0*pi*(i-1.d0)/nd(kbod)
                  arg(i) = arg0 + alpha*zf1 + dreal(zfi(i)-zfi(1))
                  zck = 2.d0*eye*zc(nbod)*arg(i) 
     *                - dconjg(zc(nbod))*cdexp(2.d0*eye*arg(i))
                  zrhs = -zck   
                  rhs(istart2+2*i-1) = rhs(istart2+2*i-1) + dreal(zrhs)
                  rhs(istart2+2*i) = rhs(istart2+2*i) + dimag(zrhs)
               end do
ccc               call PRIN2 ('  arg = *', arg, nd(kbod))
            end do
c               
            do i = 1, nd(kbod)
ccc               zfi(i) = 2.d0*eye*dimag(z(istart+i))*dreal(dz(istart+i)) 
ccc     *                - 2.d0*dimag(z(istart+i))*dimag(dz(istart+i))
               zfi(i) = zbond*dconjg(z(istart+i))*dz(istart+i) 
     *                  + dconjg(zbond)*z(istart+i)*dz(istart+i)
            end do
            call DCFFTF (nd(kbod), zfi, wsave(1,kbod))
            do i = 1, nd(kbod)
               zfi(i) = zfi(i)/nd(kbod)
            end do
            call FOURINT (nd(kbod), zfi, zf)
            call DCFFTB (nd(kbod), zf, wsave(1,kbod)) 
c
            do i = 1, nd(kbod)
               alpha = 2.d0*pi*(i-1.d0)/nd(kbod)
               zn = -0.5d0*dz(istart+i)/cdabs(dz(istart+i))
ccc               zgrav = 0.25d0*bond*(alpha*zfi(1) + zf(i))
               zgrav = -eye*0.25d0*(alpha*zfi(1) + zf(i))
               zrhs = zgrav   
               rhs(istart2+2*i-1) = rhs(istart2+2*i-1) +
     *                    dreal(zn + zrhs)
               rhs(istart2+2*i) = rhs(istart2+2*i) +
     *                    dimag(zn + zrhs)
            end do
            istart = istart + nd(kbod)
            istart2 = istart2 + 2*nd(kbod)
         end do
c  if iwall = 1, add on reflected ck's
c
ccc         call PRINF (' iwall in setvelf = *', iwall, 1)
         if (iwall.eq.1) then
            do i = 1, nbk
               zphi = 0.d0
               zphip = 0.d0
               zpsi = 0.d0
               call REFSING (z(i),k,zk,zc,zb,zphi,zphip,zpsi)
               zckr = zphi - z(i)*dconjg(zphip) - dconjg(zpsi)
ccc               call PRIn2 ('  zckr = *', zckr, 2)
               rhs(2*i-1) = rhs(2*i-1) - dreal(zckr)
               rhs(2*i) = rhs(2*i) - dimag(zckr)
            end do
         end if
         call PRINd2 ('rhs = *', rhs, 2*nd(1))
ccc         stop
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine SETVEL (k, nd, nbk, z, zk, dz, arg, zb, zc, zf, zfi, 
     *                   wsave, ifwork, rhs, iwall, bond, cap)
c---------------
c
c  The right hand side is given by integrating components of the stress
c
c     i PHI + Z CONJG(i PHI') + CONJG(i PSI) =  -sigma * n
c
c  where sigma = surface tension
c
      implicit real*8 (a-h,o-z)
      dimension rhs(*), nd(k), arg(nbk)
      complex*16 z(nbk), dz(nbk), zc(k), zf(nbk), zfi(nbk), 
     *           wsave(ifwork,k), zk(k), zb(k)
      complex*16 zn, eye, zstr
c
         eye = dcmplx(0.d0,1.d0)
         pi = 4.d0*datan(1.d0)
c
         istart = 0
         istart2 = 0
         do kbod = 1, k
            zc(kbod) = 0.d0
            do i = 1, nd(kbod)
               zn = -0.5d0*dz(istart+i)/cdabs(dz(istart+i))
               zstr = 0.5d0*cap*dconjg(z(istart+i))
               rhs(istart2+2*i-1) = dreal(zn - zstr) 
               rhs(istart2+2*i) = dimag(zn - zstr) 
            end do
            istart = istart + nd(kbod)
            istart2 = istart2 + 2*nd(kbod)
         end do
c         
c
      return
      end
c
c
c---------------
      subroutine SOLVEI (nd, k, nbk, rhs, soln, pk, pinf, za, zk, z,  
     *                   dz, u, rwork, lrwork, iwork, liwork)
c---------------
c
      implicit real*8 (a-h,o-z)
      external MATVEC, MSOLVE
c
c  Geometry
c
      dimension nd(k), pk(k)
c
c  System
c
      dimension soln(*), rhs(*)
c
c  DGMRES work arrays
c
      dimension rwork(lrwork),iwork(liwork)
c
c  Solution
c
      complex*16 u(nbk), za(k), z(nbk), dz(nbk), zb(k), zk(k)
      complex*16 eye, zx, zdis
c
c  Timings
c
      real*4 timep(2), etime
c
         eye = dcmplx(0.d0,1.d0)
         pi = 4.d0*datan(1.d0)
c
c  solve linear system using GMRES.
c
c
c     parameters for DGMRES
c
         itol = 0
         tol = 1.0d-10
         isym = 0
         do i=2,liwork
            iwork(i) = 0
         enddo
c
c  Preconditioner flag
c
         iwork(4) = 0
c
c  Restart flag
c  
         iwork(5) = -1      
c
c     provide initial guess soln
c
         norder = 2*nbk 
         do i=1,norder
            soln(i) = rhs(i)
         enddo
ccc         call PRIN2 ('  rhs = 1 *', rhs, 2*nd(1))
ccc         call PRIN2 ('  rhs = 2 *', rhs(1+2*nd(1)), 2*nd(1))
c
         t0 = etime(timep)
         call DGMRES (norder, rhs, soln, nelt, ia, ja, a, isym,
     1               MATVEC, MSOLVE, itol, tol, itmax, iter, err,  
     1               ierr, 6, sb, sx, rwork, lrwork, iwork, 
     1               liwork, rw, iw)
         call PRINI (6,13)
         call PRINF ('  # GMRES ITERATIONS = *',iter,1)
         if (ierr.gt.2) then
            call PRINF ('  SOMETHING WRONG IN GMRES, IERR = *',ierr,1)
            call PRINF ('  iwork = *',iwork,10)
            stop
           elseif (ierr.ge.0) then
            t1 = etime(timep)
            tsec = t1 - t0
            call PRIN2 (' time in Gmres = *', tsec, 1)
c
c  unpack RHS into U
c
            do i = 1,nbk
               u(i) = dcmplx(soln(2*i-1),soln(2*i))
            end do   
ccc            call PRIN2 ('   SOLUTION 1 = *',u,2*nd(1))   
ccc            call PRIN2 ('   SOLUTION 2 = *',u(1+nd(1)),2*nd(1)) 
c
c  compute pk and pinf
c  
         istart = 0
         do kbod = 1, k
            pk(kbod) = 0.d0
            h = 2.d0*pi/nd(kbod)
            do i = 1, nd(kbod)
               zdis = z(istart+i) - zk(kbod)
	       zx = dreal(u(istart+i)*dz(istart+i)/zdis**2)
               pk(kbod) = pk(kbod) + 2.d0*h*zx/pi
            end do
            istart = istart+nd(kbod)
         end do
         pinf = -pk(1)
         do kbod = 1, k
            pk(kbod) = pk(kbod) + pinf
         end do
c
c  calculate za, zb
c
            istart = 0
            circ = 0.d0
            do kbod = 1, k
               za(kbod) = 0.d0
               h = 2.d0*pi/nd(kbod)
               do i = 1,nd(kbod)
                  zn = -eye*dz(istart+i)
                  za(kbod) = za(kbod) + h*u(istart+i)
               end do
               istart = istart + nd(kbod)
            end do
            call prin2 (' A_k = *', za(1), 2*k)
            call prin2 (' p_k = *', pk(1), k)
            call PRIN2 (' pinf = *', pinf, 1)
         end if
c
      return
      end
c
c
c---------------
      subroutine MATVEC (N, XX, YY, NELT, IA, JA, A, ISYM)
c---------------
c
c  Required by DGMRES with this precise calling sequence.
c  We ignore most of the parameters except N, XX and YY
c  and call the fast multipole routine FASMVP to do the actual
c  work.
c
      implicit double precision (a-h,o-z)
      dimension xx(n), yy(n)
      parameter (ndmax=3000, kmax = 150, nmax = ndmax*kmax)
c      
      common /problem/ zk, nd, k, nbk, neq
      common /inteqn/ dsdth, rkappa, z, dz
c
      dimension dsdth(nmax), rkappa(nmax), h(kmax), nd(kmax), pk(kmax)
      complex*16 zk(kmax), dz(nmax), z(nmax), za(kmax)
c
c  local work arrays
c
      real*4 timep(2), etime
      dimension x(nmax), y(nmax)
      complex*16 dens(nmax), stress(nmax)
c
         pi = 4.d0*datan(1.d0)      
c
         t0 = etime(timep)
         do nbod = 1,k
            h(nbod) = 2.d0*pi/nd(nbod)
         end do
         do i = 1, nbk
            x(i) = dreal(z(i))
            y(i) = dimag(z(i))
         end do
c         
ccc         call PRIn2 ('  xx = *', xx, n)
         call FASMVP (nd, nbk, k, x, y, dz, dsdth, rkappa, 
     *                zk, xx, yy, h, za, pk, stress, dens)
ccc         call MSOLVE(N, yy, xx, NELT, IA, JA, A, ISYM, RWORK, IWORK)
ccc         call PRIn2 ('  xx = *', xx, n)
         t1 = etime(timep)
         tsec = t1 - t0
ccc         WRITE(13,*) 'TIME IN SECONDS FOR MATVEC = ',TSEC
ccc         WRITE(6,*) 'TIME IN SECONDS FOR MATVEC = ',TSEC
c
      RETURN
      END
c
c
c---------------
      subroutine MSOLVE(N, R, U, NELT, IA, JA, A, ISYM, RWORK, IWORK)
c---------------
c
c     Another routine required by DGMRES. It allows the use of a
c     preconditioner.
c
      implicit double precision (A-H,O-Z)
      parameter (ndmax=3000, kmax = 150, nmax = ndmax*kmax)
c      
      common /problem/ zk, nd, k, nbk, neq
      common /inteqn/ dsdth, rkappa, z, dz
c
      dimension dsdth(nmax), rkappa(nmax), h(kmax), bk(kmax)
      dimension nd(kmax)
      complex*16 zk(kmax), dz(nmax), z(nmax), zko, zki, eye
      complex*16 zc, zd, zn, zg, zu, zr 
c
      dimension r(n), u(n)
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0, 1.d0)
c
      RETURN
      END
c
c
c********1*********2*********3*********4*********5*********6*********7**
      subroutine FASMVP (nd, nbk, k, x, y, dz, dsdth, rkappa, zk, 
     *                   u, w, h, za, pk, stress, dens)
c---------------
c
c  f + z conjg(f') + conjg(g) = stress
c
      implicit real*8 (a-h,o-z)
      parameter (maxwrk = 1000000, maxiwrk = 1000000)
      dimension x(nbk), y(nbk), dsdth(nbk), rkappa(nbk), h(k)
      dimension u(*), w(*), nd(k), pk(k)
      complex*16 za(k), zb(k), zk(k), dz(nbk), zko, zki
      complex*16 stress(nbk), dens(nbk), zcent
      complex*16 z, zu, zn, zx, zself, eye, zb0, zz, zw
      complex*16 zf, zfp, zg, zc, zd, zpress, zdis, zkern
c
c  WORK ARRAYS
c
      dimension ier(5),inform(6),iwork(maxiwrk)
      complex*16 cwork(maxwrk)
c
c  Dummy variables
c
      dimension XTARG(1), YTARG(1)
      complex*16 fp(1), fpp(1)
      common /wallinfo / xleft, slength, iwall         
c
         pi = 4.d0*datan(1.d0)
         eye = dcmplx(0.d0,1.d0)
c
c  Complexify density data and extract zb
c
         do i = 1,nbk
            dens(i) = dcmplx(u(2*i-1),u(2*i))
         end do
c
c  Compute pk, za and pinf
c
         istart = 0
         do kbod = 1, k
            za(kbod) = 0.d0
            pk(kbod) = 0.d0
            do i = 1, nd(kbod)
               z = dcmplx(x(istart+i),y(istart+i))
               zdis = z - zk(kbod)
	       zx = dreal(dens(istart+i)*dz(istart+i)/zdis**2)
               pk(kbod) = pk(kbod) + 2.d0*h(kbod)*dreal(zx)/pi
               za(kbod) = za(kbod) + h(kbod)*dens(istart+i)
            end do
            istart = istart+nd(kbod)
         end do
         pinf = -pk(1)
         do kbod = 1, k
            pk(kbod) = pk(kbod) + pinf
         end do
ccc         call PRIn2 ('  pk = *', pk, k)
c
c  Call fast multipole to compute GRAD W at the points on the boundary
c  This computation does not include the self-interactions and 
c  contributions from the sources.
c
         iflag1 = 0
         if (iwall.eq.1) then
            iflag2 = 1
ccc            call PRIN2 ('  xleft = *', xleft, 1)
ccc            call PRIn2 ('  slength = *', slength, 1)
           else
            iflag2 = 0
         end if
         levmax = 7
         lenw = maxwrk
         leniw = maxiwrk
         tol = 1.d-13
         call PRINI (0,0)
         call CADAPDR(iflag1,iflag2,levmax,x,y,dz,1,k,nd,dens, h,
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
               zu = dens(istart+i)
               z = dcmplx(x(istart+i),y(istart+i))
               zn = -eye*dz(istart+i)/dsdth(istart+i)
               qwi = h(kbod)
c               
               zf = 0.5d0*zu + za(kbod)
               zfp = 0.d0 
               zg = 0.5d0*dconjg(zu)
c
               zself = zu*qwi*dsdth(istart+i)
     *                   *0.5d0*rkappa(istart+i)/pi
               zself = zself - dconjg(zu)*qwi*dsdth(istart+i)*0.5d0
     *                          *rkappa(istart+i)*zn*zn/pi
               stress(istart+i) = stress(istart+i) + zself 
     *                    - 0.5d0*eye*(pinf-pk(kbod))*z
               w(istartc+2*i-1) = dreal(dens(istart+i) + za(kbod) 
     *                                  + stress(istart+i))
               w(istartc+2*i) = dimag(dens(istart+i) + za(kbod) 
     *                                  + stress(istart+i))
            end do
            istart = istart+nd(kbod)
            istartc = istartc + 2*nd(kbod)
         end do
c
      return
      end
