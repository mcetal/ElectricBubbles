c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine GETBCs (k, nd, nbk, zk, z, rhs)
c---------------
c
      implicit real*8 (a-h,o-z)
      dimension rhs(*), nd(k)
      complex*16 zk(k), z(nbk), delz
c
         pi = 4.d0*datan(1.d0)
c         
         do i = 1, nbk
	    rhs(i) = dreal(z(i))
ccc            rhs(i) = dreal(1.d0/z(i)-zk(1))
	    rhs(i) = -2.d0*rhs(i)
         end do
         do kbod = 1, k
            rhs(nbk+kbod) = 0
ccc            u(nbk+kbod) = 0
         end do
c 
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine CHECK (nd, k, nbk, z, dz, zk, dsdth, u, ak)
c---------------
c
      implicit real*8 (a-h,o-z)
      dimension u(nbk), nd(k), ak(k), dsdth(nbk)
      complex*16 z(nbk), zk(k), dz(nbk), zstar, zexact, zn, delz, eye,
     1           z_cauchy
c
         pi = 4.d0*datan(1.d0)
	 eye = dcmplx(0.d0,1.d0)
	 h = 2.d0*pi/nd(1)
c  
         zstar = dcmplx(3.d0,0.5d0)
         delz = zstar
	 exact = dreal(delz)/(cdabs(delz))**2
c 
         sum = 0.d0
	 do i = 1, nbk
	    sum = sum + 0.5d0*h*u(i)*dsdth(i)/pi
	 end do
	 call PRIN2 (' far field constant = *', sum, 1)
ccc         answer = sum
         answer = 0.d0
         z_cauchy = 0.d0
         istart = 0
         do kbod = 1, k
            h = 2.d0*pi/nd(kbod)
            do j = 1, nd(kbod)
               delz = z(istart+j) - zstar
               zn = -eye*dz(istart+j)/dsdth(istart+j)
	       akern = dreal(delz*dconjg(zn))/(cdabs(delz))**2
	       answer = answer + 0.5d0*h*akern*u(istart+j)
     1                                *dsdth(istart+j)/pi
	       z_cauchy = z_cauchy + u(istart+j)*dz(istart+j)*h/delz
	    end do
	    istart = istart + nd(kbod)
         end do
         z_cauchy = -z_cauchy/(2.d0*pi*eye)
         call PRIN2 (' z_cauchy = *', z_cauchy, 2)
         call PRIN2 (' exact answer = *', exact, 1)
         call PRIN2 (' calc answer = *', -answer, 1)
c 
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine DIRECT_MATVEC (k, nd, nbk, zk, z, dz, rkappa, dsdth, 
     1                          ak, u, w)
c---------------
c
      implicit real*8 (a-h,o-z)
      dimension rkappa(nbk), nd(k), dsdth(nbk), u(nbk), ak(k), w(nbk)
      complex*16 zk(k), z(nbk), dz(nbk), eye, delz, zn, zcauchy
c
         pi = 4.d0*datan(1.d0)
	 eye = dcmplx(0.d0,1.d0)
c
c  extract off strength of logarithmic singularities
	 do kbod = 1, k
	    ak(kbod) = u(nbk+kbod)
	 end do
	 sum = 0.d0
	 istart = 0
	 do kbod = 1, k
	    h = 2.d0*pi/nd(kbod)
	    do i = 1, nd(kbod)
	       sum = sum + 0.5d0*h*u(istart+i)*dsdth(istart+i)/pi
	    end do
	    istart = istart + nd(kbod)
	 end do
c
c  discrete integral operator
         istart = 0
         do kbod = 1, k
            hself = 2.d0*pi/nd(kbod)
	    do i = 1, nd(kbod)
	       self = 0.25d0*hself*rkappa(istart+i)*dsdth(istart+i)/pi
	       sum_pv = self*u(istart+i)
	       zcauchy = 2.d0*pi*eye*self*u(istart+i)
	       jstart = 0
	       do nbod = 1, k
	          h = 2.d0*pi/nd(nbod)
	          do j = 1, nd(nbod)
	             if (istart+i.ne.jstart+j) then 
	                delz = z(jstart+j) - z(istart+i)
                        zn = -eye*dz(jstart+j)/dsdth(jstart+j)
	                akern = dreal(delz*dconjg(zn))/(cdabs(delz))**2
	                akern = akern*dsdth(jstart+j)
	                sum_pv  = sum_pv + 0.5d0*h*akern*u(jstart+j)/pi
	                zcauchy = zcauchy + u(jstart+j)
     1                                            *dz(jstart+j)*h/delz
	             end if
                  end do
                  jstart = jstart + nd(nbod)
               end do
               zcauchy = zcauchy/(2.d0*pi*eye)
               w(istart+i) = u(istart+i) - 2.d0*sum 
     1                                   + 2.d0*dreal(zcauchy)
               do nbod = 1, k
                  w(istart+i) = w(istart+i) 
     1              - 2.d0*ak(nbod)*dlog(cdabs(z(istart+i)-zk(nbod)))
               end do
            end do
            istart = istart + nd(kbod)
         end do
c  
c constraints
         w(nbk+k) = 0
         do kbod = 1, k
            w(nbk+k) = w(nbk+k) + ak(kbod)
         end do
         istart = 0
         do kbod = 1, k-1
            h = 2.d0*pi/nd(kbod)
            do i = 1, nd(kbod)
	       w(nbk+kbod) = w(nbk+kbod) 
     1                        + h*dsdth(istart+i)*u(istart+i)
	    end do
	    istart = istart + nd(kbod)
         end do
c
      return
      end
c
c
c********1*********2*********3*********4*********5*********6*********7**
c
      subroutine FASMVP_LAP (k, nd, nbk, zk, z, dz, rkappa, dsdth, 
     1                          ak, u, w)
c---------------
c
      implicit real*8 (a-h,o-z)
      dimension rkappa(nbk), nd(k), dsdth(nbk), u(nbk), ak(k), w(nbk)
      complex*16 zk(k), z(nbk), dz(nbk), eye, delz, zn, zcauchy
c
c  FMM parameters
      parameter (nmax=100000,nsp=500000)
      INTEGER *4 IOUT(2),INFORM(10),IERR(10)
      COMPLEX*16 UCMPLX(NMAX),PHICMPLX(NMAX), ZERO, z2pii,
     *           QA(NMAX), CFIELD(NMAX), WKSP(NSP)
      DIMENSION X(NMAX), Y(NMAX), POTEN(NMAX)
c
         pi = 4.d0*datan(1.d0)
	 eye = dcmplx(0.d0,1.d0)
	 z2pii = 1.d0/(2.d0*pi*eye)
c
c  extract off strength of logarithmic singularities
ccc	 do kbod = 1, k
ccc	    ak(kbod) = u(nbk+kbod)
ccc	 end do
ccc	 sum = 0.d0
ccc	 istart = 0
ccc	 do kbod = 1, k
ccc	    h = 2.d0*pi/nd(kbod)
ccc	    do i = 1, nd(kbod)
ccc	       sum = sum + 0.5d0*h*u(istart+i)*dsdth(istart+i)/pi
ccc	    end do
ccc	    istart = istart + nd(kbod)
ccc	 end do
c 
c set density for fmm call
         istart = 0
         do kbod = 1, k
            h = 2.d0*pi/nd(kbod)
            do i = 1, nd(kbod)
               x(istart+i) = dreal(z(istart+i))
               y(istart+i) = dimag(z(istart+i))
               qa(istart+i) = u(istart+i)*dz(istart+i)*h
            end do
            istart = istart + nd(kbod)
         end do
C
C     set parameters for FMM routine DAPIF2
         IOUT(1) = 0
         IOUT(2) = 0
c      
         IFLAG7 = 3
         NAPB = 30
         NINIRE = 2
         TOL = 1.0d-11
         MEX = 300
         EPS7 = 1.0d-16
	 call DAPIF2 (IOUT,IFLAG7,NBK,NAPB,NINIRE,MEX,IERR,INFORM,
     *                    TOL,EPS7,X,Y,QA,POTEN,CFIELD,WKSP,NSP,CLOSE)
         call PRINI(6,13)
         if (ierr(1).ne.0) then
            write (6,*) '  IERR FROM DAPIF = ',IERR(1)
            write (6,*) '  IERR = ',(ierr(ii),ii=1,4)
            write (6,*) '  INFORM = ',(inform(ii),ii=1,4)
            stop 
         end if
c
c  discrete integral operator
         istart = 0
         do kbod = 1, k
            hself = 2.d0*pi/nd(kbod)
	    do i = 1, nd(kbod)
	       self = 0.25d0*hself*rkappa(istart+i)*dsdth(istart+i)/pi
               zcauchy = self*u(istart+i) -
     1                          dreal(z2pii*cfield(istart+i))
ccc               w(istart+i) = u(istart+i) - 2.d0*sum 
ccc     1                                   + 2.d0*dreal(zcauchy)
               w(istart+i) = u(istart+i) + 2.d0*dreal(zcauchy)
               du = hself*u(istart+i)*dsdth(istart+i)
               w(istart+i) = w(istart+i)-2.d0*du
ccc               do nbod = 1, k
ccc                  w(istart+i) = w(istart+i) 
ccc     1              - 2.d0*ak(nbod)*dlog(cdabs(z(istart+i)-zk(nbod)))
ccc               end do
            end do
            istart = istart + nd(kbod)
         end do

c  
c constraints
ccc         w(nbk+k) = 0
ccc         do kbod = 1, k
ccc            w(nbk+k) = w(nbk+k) + ak(kbod)
ccc         end do
         istart = 0
         do kbod = 1, k
            h = 2.d0*pi/nd(kbod)
            do i = 1, nd(kbod)
	       w(nbk+kbod) = w(nbk+kbod) 
     1                        + h*dsdth(istart+i)*u(istart+i)
	    end do
	    istart = istart + nd(kbod)
         end do
c
      return
      end
c
c
c---------------
      subroutine DIRICH_NEU_MAP (k, nd, nbk, z, u, dz, zk, ak, phi, 
     1                           psi, dphi_dn, dsdth)
c---------------
c
      implicit real*8 (a-h,o-z)
      dimension nd(k), u(nbk), ak(k), psi(nbk), phi(nbk), dphi_dn(nbk),
     1          dsdth(nbk)
      complex*16 zk(k), z(nbk), dz(nbk), zdel, zn
c
c  FMM parameters
      parameter (nmax=100000,nsp=500000)
      COMPLEX*16 UCMPLX(NMAX),PHICMPLX(NMAX), z2pii, eye, zphi,
     *           QA(NMAX), CFIELD(NMAX), WKSP(NSP), u2(nmax)
      DIMENSION X(NMAX), Y(NMAX), x2(nmax), y2(nmax), POTEN(NMAX), 
     *          h(100)
c
c  fft arrays
      complex*16 zf(nmax), zfi(nmax)
c
         pi = 4.d0*datan(1.d0)
	 eye = dcmplx(0.d0,1.d0)
	 z2pii = 1.d0/(2.d0*pi*eye)
c
         istart = 0
         do kbod = 1, k
            h(kbod) = 2.d0*pi/nd(kbod)
            do i = 1, nd(kbod)
               x(istart+i) = dreal(z(istart+i))
               y(istart+i) = dimag(z(istart+i))
               ucmplx(istart+i) = u(istart+i)*dz(istart+i)
            end do
            istart = istart + nd(kbod)
         end do
c
         call PVINTEV(1,K,ND,nbk,X,Y,X2,Y2,h,UCMPLX,
     *                U2,QA,CFIELD,NSP,WKSP,PHICMPLX)
         call PRINI (6, 13)
c
         istart = 0
         do kbod = 1, k
            do i = 1, nd(kbod)
               zphi = -0.5d0*u(istart+i) - phicmplx(istart+i)
     1                                  - z(istart+i)
ccc               zphi = -0.5d0*u(istart+i) - phicmplx(istart+i)
               phi(istart+i) = dreal(zphi)
               psi(istart+i) = dimag(zphi)
               zf(i) = psi(istart+i)
            end do
            call DCFFTI (nd(kbod), wksp)
            call FDIFFF (zf, zfi, nd(kbod), wksp)
            do i = 1, nd(kbod)
               dphi_dn(istart+i) = -dreal(zfi(i))/dsdth(istart+i)
               zn = -eye*dz(istart+i)/dsdth(istart+i)
               alog = 0.d0
ccc               do nbod = 1, k
ccc                  zdel = z(istart+i)-zk(nbod)
ccc                  alog = alog + ak(nbod)*dreal(zdel*dconjg(zn))
ccc     1                            /(cdabs(zdel))**2
ccc               end do
               dphi_dn(istart+i) = dphi_dn(istart+i) + alog
            end do
            istart = istart + nd(kbod)
         end do
c
      return
      end
c
c
c---------------
      subroutine SOLVE_LAPL (nd, k, nbk, rhs, soln, u, ak, 
     1                       rwork, lrwork, iwork, liwork, dsdth)
c---------------
c
      implicit real*8 (a-h,o-z)
      external MATVEC_LAPL, MSOLVE
c
c  System
c
      dimension soln(*), rhs(*), u(*), ak(k), nd(k), dsdth(*)
c
c  DGMRES work arrays
c
      dimension rwork(lrwork),iwork(liwork)
c
c  Timings
c
      real*4 timep(2), etime
c
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
         norder = nbk 
         do i=1,norder
            soln(i) = rhs(i)
ccc            soln(i) = 0.d0
         enddo
c
         t0 = etime(timep)
         call prin2 (' in laplace solve, tol = *', tol, 1)
         call prinF (' in laplace solve, itmax = *', itmax, 1)
         call DGMRES (norder, rhs, soln, nelt, ia, ja, a, isym,
     1               MATVEC_LAPL, MSOLVE, itol, tol, itmax, iter, err,  
     1               ierr, 6, sb, sx, rwork, lrwork, iwork, 
     1               liwork, rw, iw)
         call Prin2 (' after laplace solve, err = *', err, 1)
         call PrinF (' after laplace solve, ierr = *', ierr, 1)
         call PRINI (6,13)
         call PRINF ('  # GMRES ITERATIONS = *',iter,1)
         if (ierr.gt.2) then
            call PRINF ('  SOMETHING WRONG IN GMRES, IERR = *',ierr,1)
            call PRINF ('  iwork = *',iwork,10)
            stop
           elseif (ierr.ge.0) then
            t1 = etime(timep)
            tsec = t1 - t0
c
c  unpack RHS into U
c

            istart = 0
            do kbod = 1, k
            farfield = 0.d0
               h = 2.d0*pi/nd(kbod)
               do i = 1, nd(kbod)
                  u(istart+i) = soln(istart+i)
                  du = u(istart+i)*dsdth(istart+i)*h
                  farfield = farfield + du
               end do
               call PRIn2 (' sum on k = *', farfield, 1)
               istart = istart + nd(kbod)
            end do   
ccc            do kbod = 1, k
ccc               ak(kbod) = soln(nbk+kbod)
ccc            end do
ccc            call PRIN2 (' ak = *', ak, k)
         end if
c
      return
      end
c
c
c---------------
      subroutine MATVEC_LAPL (N, XX, YY, NELT, IA, JA, A, ISYM)
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
      dimension dsdth(nmax), rkappa(nmax), ak(kmax), nd(kmax)
      complex*16 zk(kmax), dz(nmax), z(nmax)
c
c  local work arrays
c
      real*4 timep(2), etime
c
         t0 = etime(timep)
c         
         call FASMVP_LAP (k, nd, nbk, zk, z, dz, rkappa, dsdth, 
     1                       ak, xx, yy)
ccc         call FASMVP (nd, nbk, k, x, y, dz, dsdth, rkappa, 
ccc     *                zk, xx, yy, h, za, pk, stress, dens)
ccc         call MSOLVE(N, yy, xx, NELT, IA, JA, A, ISYM, RWORK, IWORK)
ccc         call PRIn2 ('  xx = *', xx, n)
         t1 = etime(timep)
         tsec = t1 - t0
ccc         WRITE(13,*) 'TIME IN SECONDS FOR MATVEC = ',TSEC
ccc         WRITE(6,*) 'TIME IN SECONDS FOR MATVEC = ',TSEC
c
      RETURN
      END
					 