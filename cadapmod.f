C
C**********************************************************************C
      SUBROUTINE CADAPMOD(IFLAG1,IFLAG2,LEVMAX,XAT,YAT,DZDT,K0,KBOD,ND,
     1           DENS,H,IWORK,LENIW,CWORK,LENW,XTARG,YTARG,NTARG,
     2           GRADW,PHIP,PHIPP,TOL,INFORM,IER,XLEFT,SLENGTH)
C**********************************************************************C
C
c  MODIFICATION TO DENSITY FOR INTERFACE PROBLEMS
c
C     INPUT:
C
C     XAT(i) = x-coordinate of ith source. 
C     YAT(i) = y-coordinate of ith source. 
C     DENS(i) = complex density at ith source point.
C     DZDT(i) = dz/dt on boundary, where curve is equispaced
C               with respect to parameter t. 
C               (i.e. dz = dzdt * dt)
C     K0 = 0 interior problem
C     K0 = 1 for exterior (or wall-bounded flow) problem
C     KBOD = number of inclusions for K0 = 1
C     KBOD = number of interior inclusions for K0 = 0
C     ND(K) = number of points on kth curve (Gamma_k).
C     H   =  spacing in parameter t
C            (i.e. h(k) = length of Gamma_k measured in t, divided by N_k)
C     NATOMS = number of sources
C     XTARG(i) = x-coordinate of ith source. 
C     XLEFT, SLENGTH used for wall-bounded flows to
C                    define a square containing all sources
C                    and targets.
C     YTARG(i) = y-coordinate of ith source. 
C                 _
C                   |
C                 S |
C                 L |
C                 E |
C                 N |
C                 G |
C                 T |
C                 H |
C      -------------|-------------------|--------------------
C              XLEFT <---  SLENGTH  --->
C
C     NTARG = number of targets
C
C     flag to determine what is evaluated
C
C     IFLAG1 = 0 => evaluate GRADW = PHI+z CONJG(PHI')+PSI
C                   at source points only
C     IFLAG1 = 1 => evaluate GRADW and PHIP=PHI' and PHIPP=phi''
C                   at target points only
C
C     flag to determine boundary conditions
C
C     IFLAG2 = 0 => free space
C     IFLAG2 = 1 => wall bounded flow.
C
C     LEVMAX = max number of levels in multipole hierarchy.
C
C              ***    NOTE ***
C              ...LEVMAX must be greater than or equal to one...
C
C*************************************************************
C
C     TOL  is the specified precision.
C          The precision which can be achieved by a double precision
C          version of this program can not be set higher than about
C          1.0D-14. If too high a precision is requested, the
C          subroutine returns with an error code described below.
C
C     CWORK is a workspace array dimensioned as COMPLEX *8.
C          If insufficient workspace is provided, the subroutine
C          returns with an error code described below.
C
C     IWORK is a workspace array dimensioned as INTEGER *4.
C          If insufficient workspace is provided, the subroutine
C          returns with an error code described below.
C
C     INFORM   = information returned to user
C
C              INFORM(1) is the number of levels used in the
C                 calculations (less than or equal to LEVMAX).
C              INFORM(2) is the number of terms used in the
C                 multipole expansions to achieve desired accuracy.
C
C     IER        = error code array
C
C     IER(1) = 0  => no errors encountered
C
C     IER(1) = 4  => TOL set too high. Too many terms in Hermite
C                    expansions needed.
C                    IER(2) returns the number of terms in expansion
C                    required to satisfy TOL, which should not exceed
C                    20.
C
C     IER(1) = 8  => Insufficient workspace allotted.
C                    IER(2) = amount of complex workspace needed
C                      to complete current level.
C                    IER(3) = amount of complex workspace provided.
C                    IER(4) = amount of integer workspace needed
C                      to complete current level.
C                    IER(5) = amount of integer workspace provided.
C
C     (local variables)
C
C     NTERMS     = number of terms used in expansions
C
C**********************************************************************C
      IMPLICIT REAL *8  (A-H,O-Z)
      INTEGER *4 NATOMS,NTARG,LENIW,LENW,IOUT1,IOUT2,IER(1)
      INTEGER *4 IFLAG1,IFLAG2,INFORM(1),IWORK(1),ND(K0:*)
      REAL *8  XAT(1),YAT(1),XTARG(1),YTARG(1),TOL,H(K0:*),PI
      REAL *8  XLEFT,SLENGTH
      COMPLEX *16  DZDT(1),DENS(1)
      COMPLEX *16  CWORK(1),Z2PI
      COMPLEX *16  GRADW(1),PHIP(1),PHIPP(1)
C
C----------------------------------------------
C      determine number of multipole terms needed for
C      specified precision
C
ccc      time1 = second()
      Z2PI = 1.0D0/(2*4*DATAN(1.0D0)*DCMPLX(0.0D0,1.0D0))
      IER(1) = 0
      NTERMS = - NINT(dlog(TOL)/dlog(2.0D0))
      INFORM(2) = NTERMS
      call prinf(' NTERMS IS =*',NTERMS,1)
      IF (NTERMS .GT. 60) THEN
	 IER(1) = 4
	 IER(2) = NTERMS
	 RETURN
      ENDIF
C
C-----  scale sources to 64 x 64 box  and compute scaling constants.
C
      NATOMS = 0
      BOXSIZ = 64.0D0
      DO 20 K = K0,KBOD
         NATOMS = NATOMS + ND(K)
20    CONTINUE
      IF (IFLAG2.EQ.1) THEN
	 RSCAL = BOXSIZ/SLENGTH
         RMINX = XLEFT
         RMINY = 0.0D0
      ELSE
         RMINX = XAT(1)
         RMAXX = XAT(1)
         RMINY = YAT(1)
         RMAXY = YAT(2)
C
C
         DO 100 I=1,NATOMS
	    IF (XAT(I) .GT. RMAXX) RMAXX = XAT(I)
	    IF (YAT(I) .GT. RMAXY) RMAXY = YAT(I)
	    IF (XAT(I) .LT. RMINX) RMINX = XAT(I)
	    IF (YAT(I) .LT. RMINY) RMINY = YAT(I)
 100     CONTINUE
         DO 200 I=1,NTARG
	    IF (XTARG(I) .GT. RMAXX) RMAXX = XTARG(I)
	    IF (YTARG(I) .GT. RMAXY) RMAXY = YTARG(I)
	    IF (XTARG(I) .LT. RMINX) RMINX = XTARG(I)
	    IF (YTARG(I) .LT. RMINY) RMINY = YTARG(I)
 200     CONTINUE
         RSX = RMAXX-RMINX
         RSY = RMAXY-RMINY
         XONEW = (RMAXX+RMINX)/2.0
         YONEW = (RMAXY+RMINY)/2.0
ccc      RSCAL = AMAX1(RSX,RSY)*1.001
         RSCAL = MAX(RSX,RSY)
         RSCAL = BOXSIZ/RSCAL
      call prin2(' RSCAL = *',RSCAL,1)
      ENDIF
C
C
C----- allocate arrays used by algorithm from workspace.
C
      IZSRC = 1
      IQA = IZSRC + NATOMS
      IQB = IQA + NATOMS
      IQC = IQB + NATOMS
      IZSRC2 = IQC + NATOMS
      IQA2 = IZSRC2 + NATOMS
      IQB2 = IQA2 + NATOMS
      IQC2 = IQB2 + NATOMS
      IZTARG = IQC2 + NATOMS
      ITOTP = IZTARG + NTARG
ccc      call prinf(' IZTARG = *',IZTARG,1)
ccc      call prinf(' ITOTP = *',ITOTP,1)
ccc      call prinf(' NATOMS = *',NATOMS,1)
      IF ( ITOTP .GT. LENW) THEN
	 IER(1) = 8
	 IER(2) = ITOTP
	 IER(3) = LENW
	 write(6,*)' insufficient workspace'
	 write(13,*)' insufficient workspace'
	 RETURN
      ELSE
	 J = 0
	 DO 300 K = K0,KBOD
            DO 250 NK=1,ND(K)
	       J = J + 1
	       CWORK(IZSRC+J-1) = 
ccc     1              DCMPLX(XAT(J)-XONEW,YAT(J)-YONEW)*RSCAL
     1              DCMPLX(XAT(J)-RMINX,YAT(J)-RMINY)*RSCAL
c********
c  Leslie's original signs
c
ccc	       CWORK(IQA+J-1) = Z2PI*DENS(J)*DZDT(J)*H(K)
ccc	       CWORK(IQB+J-1) = (DCONJG(DENS(J))*DZDT(J) +
ccc     1              DENS(J)*DCONJG(DZDT(J)))*Z2PI*H(K)
ccc	       CWORK(IQC+J-1) = 
ccc     1          -Z2PI*DCONJG(CWORK(IZSRC+J-1))*DENS(J)*DZDT(J)*H(K) 
c
c*******
c  My fudging
c
	       CWORK(IQA+J-1) = Z2PI*DENS(J)*DZDT(J)*H(K)
	       CWORK(IQB+J-1) = +(-DCONJG(DENS(J))*DZDT(J) +
     1              DENS(J)*DCONJG(DZDT(J)))*Z2PI*H(K)
	       CWORK(IQC+J-1) = 
     1          -Z2PI*DCONJG(CWORK(IZSRC+J-1))*DENS(J)*DZDT(J)*H(K) 
c
c*******
 250        CONTINUE
	    ISTART = ISTART + ND(K)
 300     CONTINUE
         DO 400 J=1,NTARG
	    CWORK(IZTARG+J-1) = 
ccc     1                DCMPLX(XTARG(J)-XONEW,YTARG(J)-YONEW)*RSCAL
     1                DCMPLX(XTARG(J)-RMINX,YTARG(J)-RMINY)*RSCAL
 400     CONTINUE
      ENDIF
      IF (IFLAG1.EQ.0) THEN
         DO 500 J=1,NATOMS
	    GRADW(J) = 0.0D0
 500     CONTINUE
      ELSE
         DO 600 J=1,NTARG
	    GRADW(J) = 0.0D0
	    PHIP(J) = 0.0D0
	    PHIPP(J) = 0.0D0
600     CONTINUE
      ENDIF
C
C     process boxes at each level until LEVMAX is reached or the 
C     number of particles per box is sufficiently small 
C     
C
ccc      time2 = second()
      call prin2(' CADAP initial. time = *',time2-time1,1)
      DO 1000 ILEV = 1,LEVMAX
	 INFORM(1) = ILEV
         call prinf(' LEVMAX = *',LEVMAX,1)
         call prinf(' ILEV = *',ILEV,1)
ccc         time1 = second()
ccc         call prin2(' GRADW = *',GRADW,2*NATOMS)
C
C        allocate memory to integer workspace
C
         NDIM = 2**ILEV
         NALLBX = NDIM*NDIM
C
         LICNT  = NALLBX
         LICNT2 = NALLBX
         LLOC   = NALLBX
         LTOFST = NALLBX+1
         LOFFST = NALLBX+1
         LIBOX  = MAX(NATOMS,NTARG)
         LATADR = NATOMS
         LTRADR = NTARG
C
         ICNT = 1
         ICNT2 = ICNT + LICNT
         IOFFST = ICNT2 + LICNT2
         ITOFST = IOFFST + LOFFST
         IBOX = ITOFST + LTOFST
         IATADR = IBOX + LIBOX
         ITRADR = IATADR + LATADR
         ILOC = ITRADR + LTRADR
         ITOT = ILOC + LLOC
c
         IF ( ITOT .GT. LENIW) THEN
	    IER(1) = 8
	    IER(4) = ITOT
	    IER(5) = LENIW
	    write(6,*)' insufficient workspace'
	    write(13,*)' insufficient workspace'
	    RETURN
	 ENDIF
c
C
C     assign particles to boxes 
C
ccc	 call prin2(' in cadap, before ASSIGN1, ZSRC = *',
ccc     1                CWORK(IZSRC),2*NATOMS)
ccc	 time1 = second()
         CALL ASSIGN1(IFLAG1,ILEV,CWORK(IZSRC),IWORK(IOFFST), 
     1         CWORK(IZTARG),NTARG,IWORK(IBOX),NATOMS,
     2         IWORK(ICNT),IWORK(ICNT2),IWORK(IATADR),
     3         IWORK(ITRADR),IWORK(ITOFST),IWORK(ILOC),NTBOX,NSMAX)
ccc	 time2 = second()
         call prin2(' time for ASSIGN1 = *',time2-time1,1)
	 call prinf(' NTBOX = *',NTBOX,1)
	 call prinf(' NSMAX = *',NSMAX,1)
         LLOCXP = (NTERMS+1)*NTBOX
         ILOCXP1 = ITOTP
         ILOCXP2 = ILOCXP1 + LLOCXP
         ITOTC = ILOCXP2 + LLOCXP
C
ccc         call prin2(' after ASSIGN1 GRADW = *',GRADW,2*NATOMS)
         IF ( ITOTC .GT. LENW) THEN
	    IER(1) = 8
	    IER(2) = ITOTC
	    IER(3) = LENW
	    write(6,*)' insufficient workspace'
	    write(13,*)' insufficient workspace'
	    RETURN
         ENDIF
C
ccc      call prinf('IBOX = *',IWORK(IBOX),NATOMS)
ccc      call prinf('IOFFST = *',IWORK(IOFFST),NALLBX+1)
C
	 IF (IFLAG1.EQ.0) THEN
            CALL DOINT0(IFLAG2,ILEV,CWORK(IZSRC),CWORK(IQA),
     1        CWORK(IQB),CWORK(IQC),CWORK(IZSRC2),CWORK(IQA2),
     2        CWORK(IQB2),CWORK(IQC2),NATOMS,CWORK(IZTARG),NTARG,
     3        CWORK(ILOCXP1),CWORK(ILOCXP2),NTERMS,IWORK(IOFFST),
     4        IWORK(IBOX),IWORK(IATADR),IWORK(ITRADR),IWORK(ITOFST), 
     5        IWORK(ICNT),IWORK(ICNT2),IWORK(ILOC),
     6        NTBOX,GRADW(1))
	 ELSE
            CALL DOINT1(IFLAG2,ILEV,CWORK(IZSRC),CWORK(IQA),
     1        CWORK(IQB),CWORK(IQC),CWORK(IZSRC2),CWORK(IQA2),
     2        CWORK(IQB2),CWORK(IQC2),NATOMS,CWORK(IZTARG),NTARG,
     3        CWORK(ILOCXP1),CWORK(ILOCXP2),NTERMS,IWORK(IOFFST),
     4        IWORK(IBOX),IWORK(IATADR),IWORK(ITRADR),IWORK(ITOFST),
     5        IWORK(ICNT),IWORK(ICNT2),IWORK(ILOC),
     6        NTBOX,GRADW(1),PHIP(1),PHIPP(1))
         ENDIF
ccc         call prin2(' after DOINT0 GRADW = *',GRADW,2*NATOMS)
ccc	 IF ((ILEV.EQ.LEVMAX).OR.(NSMAX.LE.NTERMS)) GOTO 1001
ccc         time2 = second()
         call prin2(' time for current level = *',time2-time1,1)
	 IF (ILEV.EQ.LEVMAX) GOTO 1001
1000  CONTINUE
1001  CONTINUE
      call prin2(' after all far fields GRADW = *',GRADW,2*NATOMS)
      call prinf(' doing near neighbors *',ILEV,0)
      call prinf(' ILEV =  *',ILEV,1)
ccc      timen = second()
      IF (IFLAG1.EQ.0) THEN
         CALL DONN0(IFLAG2,ILEV,CWORK(IZSRC),CWORK(IQA),CWORK(IQB),
     1      CWORK(IQC),NATOMS,CWORK(IZSRC2),CWORK(IQA2),CWORK(IQB2),
     2      CWORK(IQC2),IWORK(IOFFST),IWORK(IBOX),IWORK(IATADR),
     3      GRADW(1))
      ELSE
         CALL DONN1(IFLAG2,ILEV,CWORK(IZSRC),CWORK(IQA),CWORK(IQB),
     1      CWORK(IQC),NATOMS,CWORK(IZSRC2),CWORK(IQA2),CWORK(IQB2),
     2      CWORK(IQC2),CWORK(IZTARG),NTARG,IWORK(IOFFST), 
     3      IWORK(IATADR),IWORK(ITOFST),IWORK(IBOX),IWORK(ITRADR),
     4      GRADW(1),PHIP(1),PHIPP(1)) 
      ENDIF
ccc      timen2 = second()
ccc      call prin2(' after DONN GRADW = *',GRADW,2*NATOMS)
      call prin2(' time for DONN0 is = *',timen2-timen,1)
C
C     rescale computed fields
C
      IF (IFLAG1 .EQ. 0) THEN
	DO I = 1,NATOMS
	   GRADW(I) = GRADW(I)*RSCAL
        ENDDO
      ELSE       
	RSCAL2 = RSCAL**2
	RSCAL3 = RSCAL**3
	DO I = 1,NTARG
	   GRADW(I) = GRADW(I)*RSCAL
c
c   phi = phi(i)*rscal
c
	   PHIP(I) = PHIP(I)*RSCAL2
	   PHIPP(I) = PHIPP(I)*RSCAL3
        ENDDO
      ENDIF
      RETURN
      END
C
