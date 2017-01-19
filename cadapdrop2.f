C********1*********2*********3*********4*********5*********6*********7*C
      SUBROUTINE DOINT0DR(IFLAG2,ILEV,ZSRC,QA,QB,QC,ZSRC2,QA2,QB2,QC2,
     1      NATOMS,ZTARG,NTARG,LOCEXP1,LOCEXP2,NTERMS,IOFFST,IBOX,
     2      IATADR,ITRADR,ITOFST,ICNT,ICNT2,ILOC,NTBOX,GRADW)
C
C     The main subroutine of FMM for IFLAG1 = 0 in CADAP.
C     That is, we evaluate fields at source locations themselves,
C     ignoring self.
C
C     STEP 1)  Sources are assigned to boxes.
C
C     Boxes containing sources are then dealt with sequentially.
C
C     STEP 2)  Far field expansions are formed for each box
C              containing sources.
C
C     The interaction list is then searched for source pts and 
C     processed by means of the following decision analysis,
C     where NT is the number of source pts in interaction box.
C
C         IF (NT > NTERMS/2) convert far field expansion to
C                               local expansion.
C
C         IF (NT .LE. NTERMS/2)  evaluate far field expansion
C                                  directly at target positions.
C
C     IFLAG2 = 0 => free space
C     IFLAG2 = 1 => wall bounded flow.
C
C*****************************************************************
      IMPLICIT REAL *8 (A-H,O-Z)
      INTEGER *4 NTERMS,NATOMS,IOFFST(1),NTBOX,IATADR(1)
      INTEGER *4 IBOX(1),ICNT(1),ICNT2(1),NSIDE,NTARG,ILOC(1)
      INTEGER *4 NLMAX,ITOFST(1),ITRADR(1)
      INTEGER *4 INTACT(30)
      REAL *8 C(120,120),H
      COMPLEX *16  ZSRC(1),QA(1),QB(1),QC(1)
      COMPLEX *16  ZSRC2(1),QA2(1),QB2(1),QC2(1)
      COMPLEX *16  ZTARG(1),GRADW(1),ZCENT,ZCENTR,ZCENT2
      COMPLEX *16  MPOLE1(60),MPOLE2(60),B1(0:60),B2(0:60)
      COMPLEX *16  MPR1(60),MPR2(60)
      COMPLEX *16  LOCEXP1(0:NTERMS,1),LOCEXP2(0:NTERMS,1)
      COMPLEX *16  PHIPT,PHIPPT
      DATA ZERO/0.0d0/
C
C     compute binomial coefficients
C
ccc      call prinf(' In DOINT0, ILEV =  *',ILEV,1)
ccc      call prin2(' In DOINT0, ZSRC =  *',ZSRC,2*NATOMS)
ccc      call prinf(' computing binomial coeffs *',LDC,0)
      NLMAX = NTERMS/2
ccc      NLMAX = 0
      LDC = 120
      C(1,1) = 1.0D0
      DO I=2,120
	 C(1,I) = 0.0D0
	 C(I,1) = 1.0D0
      ENDDO
      DO I=2,120
	 DO J=2,120
	    C(I,J)=C(I-1,J-1)+C(I-1,J)
         ENDDO
      ENDDO
ccc      call prin2(' GRADW = *',GRADW,2*NATOMS)
C
C-----initialize local expansions to zero.
C
ccc      call prinf(' initializing local expansions *',LDC,0)
      IFLAG1 = 0
      NSIDE = 2**(ILEV)
      NBOXES = NSIDE*NSIDE
      H = 64.0D0/NSIDE
      DO 70 I = 1,NTBOX
	 DO 60 J = 0,NTERMS
		LOCEXP1(J,I) = ZERO
		LOCEXP2(J,I) = ZERO
 60      CONTINUE
 70   CONTINUE
c
c---- process all boxes
c
      NSHIFT = 0
      NSLMVAL = 0
      tdir = 0
      tlota1 = 0
      tmpole = 0
      tlota = 0
      DO 1000 I = 1,NBOXES
ccc         call prinf(' processing box I =*',I,1)
	 IOFF =   IOFFST(I)
	 NINBOX = IOFFST(I+1) - IOFF
ccc         call prinf(' NINBOX =*',NINBOX,1)
         DO K = 1,NINBOX
            ZSRC2(K) = ZSRC(IATADR(IOFF+K))
            QA2(K) = QA(IATADR(IOFF+K))
            QB2(K) = QB(IATADR(IOFF+K))
            QC2(K) = QC(IATADR(IOFF+K))
         ENDDO
ccc	 call prin2(' zsrc2 = *',zsrc2,2*NINBOX)
	 call prin2(' QA2 = *',QA2,2*NINBOX)
	 call prin2(' QB2 = *',QB2,2*NINBOX)
	 call prin2(' QC2 = *',QC2,2*NINBOX)
	 IF ( NINBOX .LE. 0 ) GOTO 1000
c
c---- create far field expansion
c
   	 ICOL = 1 + MOD(I-1,NSIDE)
 	 IROW = 1 + (I-1)/NSIDE
ccc	 XC = -32.0D0+ICOL*H-H/2
ccc	 YC = -32.0D0+IROW*H-H/2
	 XC = ICOL*H-H/2
	 YC = IROW*H-H/2
	 ZCENT = DCMPLX(XC,YC)
	 CALL SLMKXP(ZCENT,ZSRC2,QA2,QB2,QC2,
     1              NINBOX,MPOLE1,MPOLE2,NTERMS)
ccc	 call prinf(' ninbox is greater than nfmax *',NLMAX,0)
	 CALL MKINTL(ILEV,I,INTACT,NINTAC)
ccc	 call prinf(' INTACT list is *',INTACT,NINTAC)
	 DO 800 J = 1,NINTAC
	    IADR = INTACT(J)
	    INOFF = IOFFST(IADR)
	    NINNBR = IOFFST(IADR+1) - INOFF
	    IF (NINNBR .LE. NLMAX) THEN
c
c---- evaluate far field expansion
c
ccc	       t1 = second()
	       DO 700 K = 1,NINNBR
		  JT = IATADR(INOFF+K)
		  CALL SLMVALDR(IFLAG1,MPOLE1,MPOLE2,ZCENT,NTERMS,
     1                       ZSRC(JT),GRADW(JT),PHIPT,PHIPPT)
		  NSLMVAL = NSLMVAL+1
700            CONTINUE
ccc	       t2 = second()
ccc	       tmpole = tmpole + (t2-t1)
	    ELSE
c
C     convert far field expansion to local expansion in neighbor box
c
ccc	       t1 = second()
   	       ICOL = 1 + MOD(IADR-1,NSIDE)
 	       IROW = 1 + (IADR-1)/NSIDE
ccc	       XC = -32.0D0+ICOL*H-H/2
ccc	       YC = -32.0D0+IROW*H-H/2
	       XC = ICOL*H-H/2
	       YC = IROW*H-H/2
	       ZCENT2 = DCMPLX(XC,YC)
	       IADLOC = ILOC(IADR)
	       CALL SLSHIFT(MPOLE1,MPOLE2,ZCENT,NTERMS,ZCENT2,
     1               B1,B2,C,LDC)
	       NSHIFT = NSHIFT+1
	       CALL ADDEXP(B1,B2,LOCEXP1(0,IADLOC),
     1                LOCEXP2(0,IADLOC),NTERMS)
ccc	       t2 = second()
ccc	       tlota = tlota + (t2-t1)
	    ENDIF
800      CONTINUE
c
c    include effect of image boxes for wall bounded problems.
c
ccc	 call prinf(' processing image boxes for I = *',I,1)
         IF ((IFLAG2.EQ.1) .AND. (I .LE. 2*NSIDE)) THEN
	    CALL MPREFLECT(ZCENT,MPOLE1,MPOLE2,NTERMS,
     1              MPR1,MPR2)
	    CALL MKINTR(ILEV,I,INTACT,NINTAC)
ccc	 call prinf(' reflected INTACT list is *',INTACT,NINTAC)
	    ZCENTR = DCONJG(ZCENT)
	    DO 900 J = 1,NINTAC
	       IADR = INTACT(J)
	       INOFF = IOFFST(IADR)
	       NINNBR = IOFFST(IADR+1) - INOFF
	       IF (NINNBR .LE. NLMAX) THEN
c
c---- evaluate far field expansion
c
ccc	          t1 = second()
	          DO 850 K = 1,NINNBR
		     JT = IATADR(INOFF+K)
		     CALL SLMVALDR(IFLAG1,MPR1,MPR2,ZCENTR,NTERMS,
     1                       ZSRC(JT),GRADW(JT),PHIPT,PHIPPT)
		     NSLMVAL = NSLMVAL+1
850               CONTINUE
ccc	          t2 = second()
ccc	          tmpole = tmpole + (t2-t1)
	       ELSE
c
C     convert far field expansion to local expansion in neighbor box
c
ccc	          t1 = second()
   	          ICOL = 1 + MOD(IADR-1,NSIDE)
 	          IROW = 1 + (IADR-1)/NSIDE
ccc	          XC = -32.0D0+ICOL*H-H/2
ccc	          YC = -32.0D0+IROW*H-H/2
	          XC = ICOL*H-H/2
	          YC = IROW*H-H/2
	          ZCENT2 = DCMPLX(XC,YC)
	          IADLOC = ILOC(IADR)
	          CALL SLSHIFT(MPR1,MPR2,ZCENTR,NTERMS,ZCENT2,
     1               B1,B2,C,LDC)
	          NSHIFT = NSHIFT+1
	          CALL ADDEXP(B1,B2,LOCEXP1(0,IADLOC),
     1                   LOCEXP2(0,IADLOC),NTERMS)
ccc	          t2 = second()
ccc	          tlota = tlota + (t2-t1)
	       ENDIF
900         CONTINUE
         ENDIF
1000  CONTINUE
ccc      call prinf(' finished doing interactions *',NBOXES,0)
ccc      call prinf(' LOTAs =*',NSHIFT,1)
ccc      call prinf(' number of SLMVAL  =*',NSLMVAL,1)
      call prin2(' GRADW = *',GRADW,2*NATOMS)
ccc      call prinf(' evaluating all local expansions *',NBOXES,0)
ccc      tloc0 = second()
      DO 1200 I = 1,NBOXES
	 INOFF = IOFFST(I)
	 NINBOX = IOFFST(I+1) - INOFF
	 IF (NINBOX .LE. NLMAX) THEN
ccc	 IF (NINBOX .EQ. 0) THEN
	    GOTO 1200
	 ELSE
   	    ICOL = 1 + MOD(I-1,NSIDE)
 	    IROW = 1 + (I-1)/NSIDE
ccc	    XC = -32.0D0+ICOL*H-H/2
ccc	    YC = -32.0D0+IROW*H-H/2
	    XC = ICOL*H-H/2
	    YC = IROW*H-H/2
	    ZCENT = DCMPLX(XC,YC)
	    JADR = ILOC(I)
ccc	    call prinf( 'JADR = *',JADR,1)
ccc	    call prin2( 'locexp1 = *',LOCEXP1(0,JADR),2*(NTERMS+1))
ccc	    call prin2( 'locexp2 = *',LOCEXP2(0,JADR),2*(NTERMS+1))
	    DO 1100 K = 1,NINBOX
	       JT = IATADR(INOFF+K)
	       CALL SLLVALDR(IFLAG1,LOCEXP1(0,JADR),LOCEXP2(0,JADR),
     1               NTERMS,ZCENT,ZSRC(JT),GRADW(JT),PHIPT,PHIPPT)
1100        CONTINUE
	 ENDIF
1200  CONTINUE
ccc      tloc1 = second()
ccc      call prin2(' time for eval local expansions *',tloc1-tloc0,1)
ccc      call prin2(' time for eval mpole expansions *',tmpole,1)
ccc      call prin2(' time for LOTA1 *',tlota1,1)
ccc      call prin2(' time for LOTAs *',tlota,1)
ccc      call prin2(' time for direct int list *',tdir,1)
ccc      call prinf(' finished evaluating local expansions *',NBOXES,0)
      call prin2(' after local expansions GRADW = *',GRADW,2*NATOMS)
      RETURN
      END
C
      SUBROUTINE SLMVALDR(IFLAG1,MPOLE1,MPOLE2,ZCENT,NTERMS,ZTARG,
     1                  GW,PHIP,PHIPP)
C
C     This subroutine evaluates the far field expansions for PHI
C     and PSI about ZCENT to obtain GW (= W_x + i W_y)  at ZTARG.

C     INPUT
C
C     ZCENT = center opf expansions
C     MPOLE1(i) = (i)th coefficient of far field PHI expansion
C     MPOLE2(i) = (i)th coefficient of far field PSI expansion
C     NTERMS = number of terms in expansion
C
C     OUTPUT:
C
C     GW = PHI + Z*DCONJG(PHI') + DCONJG(PSI)
C
      IMPLICIT REAL *8 (A-H,O-Z)
      INTEGER *4 NTERMS
      COMPLEX *16 ZCENT,ZTARG,GW
      COMPLEX *16 MPOLE1(NTERMS),MPOLE2(NTERMS)
      COMPLEX *16 ZSHIFT,ZSINV
      COMPLEX *16 PHI,PHIP,PSI,PHIPP,PHIPLOC,PHIPPLOC
C
c---- evaluate expansion
C
      ZSHIFT = ZTARG - ZCENT
      ZSINV = 1.0D0/ZSHIFT
      PHI = MPOLE1(NTERMS)*ZSINV
      PHIPLOC = NTERMS*MPOLE1(NTERMS)*ZSINV
      PHIPPLOC = NTERMS*(NTERMS+1)*MPOLE1(NTERMS)*ZSINV
      PSI = MPOLE2(NTERMS)*ZSINV
      DO 300 J = NTERMS-1,1,-1
c
c  rename this philoc
c
	 PHI = (PHI + MPOLE1(J))*ZSINV
	 PSI = (PSI + MPOLE2(J))*ZSINV
         PHIPLOC = (PHIPLOC + J*MPOLE1(J))*ZSINV
         PHIPPLOC = (PHIPPLOC + J*(J+1)*MPOLE1(J))*ZSINV
300   CONTINUE
c
c  phi = phi + philoc, philoc doesn't need rescaling
c
      PHIPLOC = -PHIPLOC*ZSINV
      PHIPPLOC = PHIPPLOC*ZSINV*ZSINV
      PHIP = PHIP + PHIPLOC
      PHIPP = PHIPP + PHIPPLOC
      GW = GW + PHI - ZTARG*DCONJG(PHIPLOC) - DCONJG(PSI)
      RETURN
      END
C
C**********************************************************************
      SUBROUTINE SLLVALDR(IFLAG1,LOCAL1,LOCAL2,NTERMS,ZCENT,ZTARG,
     1                  GRADW,PHIP,PHIPP)
C
C     This subroutine evaluates the local expansions LOCAL1 and
C     LOCAL2 for PHI and PSI at the location ZTARG.
C
C     INPUT
C
C     ZCENT = center of the expansion
C     LOCAL1(i) = (i)th coefficient of local PHI expansion
C     LOCAL2(i) = (i)th coefficient of local PSI expansion
C     NTERMS = number of terms in expansion
C
C     OUTPUT:
C
C     GRADW = evaluated velocity
C
      IMPLICIT REAL *8 (A-H,O-Z)
      INTEGER *4 NTERMS
      COMPLEX *16  ZCENT,GRADW,ZTARG,ZSHIFT
      COMPLEX *16 LOCAL1(0:NTERMS),LOCAL2(0:NTERMS)
      COMPLEX *16 PHI,PHIP,PSI,PHIPP,PHIPLOC,PHIPPLOC
C
c---- evaluate expansion
C
ccc      call prin2(' LOCAL1 = *',LOCAL1(0),2*NTERMS+2)
ccc      call prin2(' LOCAL2 = *',LOCAL2(0),2*NTERMS+2)
      ZSHIFT = ZTARG - ZCENT
      PHI = 0.0D0
      PHIPLOC = 0.0D0
      PHIPPLOC = 0.0D0
      PSI = 0.0D0
      DO 300 J = NTERMS,0,-1
	 PHI = LOCAL1(J) + ZSHIFT*PHI
	 PSI = LOCAL2(J) + ZSHIFT*PSI
300   CONTINUE
      DO 500 J = NTERMS,1,-1
         PHIPLOC = J*LOCAL1(J) + ZSHIFT*PHIPLOC
500   CONTINUE
      DO 700 J = NTERMS,2,-1
         PHIPPLOC = J*(J-1)*LOCAL1(J) + ZSHIFT*PHIPPLOC
700   CONTINUE
      GRADW = GRADW+PHI-ZTARG*DCONJG(PHIPLOC)-DCONJG(PSI)
      PHIP = PHIP+PHIPLOC
      PHIPP = PHIPP+PHIPPLOC
      RETURN
      END
C
C****************************************************************
      SUBROUTINE DOINT1DR(IFLAG2,ILEV,ZSRC,QA,QB,QC,ZSRC2,QA2,QB2,QC2,
     1      NATOMS,ZTARG,NTARG,LOCEXP1,LOCEXP2,NTERMS,IOFFST,IBOX,
     2      IATADR,ITRADR,ITOFST,ICNT,ICNT2,ILOC,NTBOX,
     3      GRADW,PHIP,PHIPP)
C
C     The main subroutine of FMM for IFLAG1 = 1 in CADAP.
C     That is, we evaluate fields at target locations only.
C
C     STEP 1)  Sources are assigned to boxes.
C
C     Boxes containing sources are then dealt with sequentially.
C
C     STEP 2)  Far field expansions are formed for each box
C              containing sources.
C
C     The interaction list is then searched for targets and 
C     processed by means of the following decision analysis,
C     where NT is the number of targets in target box.
C
C         IF (NT > NLMAX) convert far field expansion to
C                               local expansion.
C
C         IF (NT .LE. NLMAX)  evaluate far field expansion
C                                  directly at target positions.
C
C*****************************************************************
      IMPLICIT REAL *8 (A-H,O-Z)
      INTEGER *4 NTERMS,NATOMS,IOFFST(1),NTBOX,IATADR(1)
      INTEGER *4 IBOX(1),ICNT(1),ICNT2(1),NSIDE,NTARG,ILOC(1)
      INTEGER *4 NLMAX,ITOFST(1),ITRADR(1)
      INTEGER *4 INTACT(30)
      REAL *8 C(120,120)
      COMPLEX *16  ZSRC(1),QA(1),QB(1),QC(1)
      COMPLEX *16  ZSRC2(1),QA2(1),QB2(1),QC2(1)
      COMPLEX *16  ZTARG(1),GRADW(1),ZCENT,ZCENT2,ZCENTR
      COMPLEX *16  PHIP(1),PHIPP(1)
      COMPLEX *16  MPOLE1(60),MPOLE2(60),B1(60),B2(60)
      COMPLEX *16  MPR1(60),MPR2(60)
      COMPLEX *16  LOCEXP1(0:NTERMS,1),LOCEXP2(0:NTERMS,1)
      DATA ZERO/0.0d0/
C
C     compute binomial coefficients
C
ccc      call prin2(' In DOINT1, ZSRC =  *',ZSRC,2*NATOMS)
ccc      call prin2(' In DOINT1, ZTARG =  *',ZTARG,2*NTARG)
      LDC = 120
ccc      NLMAX = NTERMS/2
      NLMAX = 0
      C(1,1) = 1.0D0
      DO I=2,120
	 C(1,I) = 0.0D0
	 C(I,1) = 1.0D0
      ENDDO
      DO I=2,120
	 DO J=2,120
	    C(I,J)=C(I-1,J-1)+C(I-1,J)
         ENDDO
      ENDDO
C
C-----initialize local expansions to zero.
C
      IFLAG1 = 1
      NSIDE = 2**(ILEV)
      NBOXES = NSIDE*NSIDE
      H = 64.0D0/NSIDE
      DO 70 I = 1,NTBOX
	 DO 60 J = 0,NTERMS
		LOCEXP1(J,I) = ZERO
		LOCEXP2(J,I) = ZERO
 60      CONTINUE
 70   CONTINUE
c
c---- process all boxes
c
      NSHIFT = 0
      NSLMVAL = 0
      tdir = 0
      tlota1 = 0
      tmpole = 0
      tlota = 0
      DO 1000 I = 1,NBOXES
ccc         call prinf(' processing box I =*',I,1)
	 IOFF =   IOFFST(I)
	 NINBOX = IOFFST(I+1) - IOFF
ccc         call prinf(' NINBOX =*',NINBOX,1)
	 DO K = 1,NINBOX
	    ZSRC2(K) = ZSRC(IATADR(IOFF+K))
	    QA2(K) = QA(IATADR(IOFF+K))
	    QB2(K) = QB(IATADR(IOFF+K))
	    QC2(K) = QC(IATADR(IOFF+K))
	 ENDDO
	 IF ( NINBOX .LE. 0 ) GOTO 1000
c
c---- create far field expansion
c
   	 ICOL = 1 + MOD(I-1,NSIDE)
 	 IROW = 1 + (I-1)/NSIDE
ccc	 XC = -32.0D0+ICOL*H-H/2
ccc	 YC = -32.0D0+IROW*H-H/2
	 XC = ICOL*H-H/2
	 YC = IROW*H-H/2
	 ZCENT = DCMPLX(XC,YC)
	 CALL SLMKXP(ZCENT,ZSRC2,QA2,QB2,QC2,
     1           NINBOX,MPOLE1,MPOLE2,NTERMS)
	 CALL MKINTL(ILEV,I,INTACT,NINTAC)
	 DO 800 J = 1,NINTAC
	    IADR = INTACT(J)
	    INOFF = ITOFST(IADR)
	    NINNBR = ITOFST(IADR+1) - INOFF
	    IF (NINNBR .LE. NLMAX) THEN
c
c---- evaluate far field expansion
c
	       DO 700 K = 1,NINNBR
		  JT = ITRADR(INOFF+K)
		  CALL SLMVALDR(IFLAG1,MPOLE1,MPOLE2,ZCENT,NTERMS,
     1                 ZTARG(JT),GRADW(JT),PHIP(JT),PHIPP(JT))
700            CONTINUE
	    ELSE
c
C     convert far field expansion to local expansion in neighbor box
c
ccc	       t1 = second()
   	       ICOL = 1 + MOD(IADR-1,NSIDE)
 	       IROW = 1 + (IADR-1)/NSIDE
ccc	       XC = -32.0D0+ICOL*H-H/2
ccc	       YC = -32.0D0+IROW*H-H/2
	       XC = ICOL*H-H/2
	       YC = IROW*H-H/2
	       ZCENT2 = DCMPLX(XC,YC)
	       IADLOC = ILOC(IADR)
	       CALL SLSHIFT(MPOLE1,MPOLE2,ZCENT,NTERMS,ZCENT2,
     1               B1,B2,C,LDC)
	       NSHIFT = NSHIFT+1
	       CALL ADDEXP(B1,B2,LOCEXP1(0,IADLOC),
     1                 LOCEXP2(0,IADLOC),NTERMS)
ccc	       t2 = second()
ccc	       tlota = tlota + (t2-t1)
	    ENDIF
800      CONTINUE
         IF ((IFLAG2.EQ.1) .AND. (I .LE. 2*NSIDE)) THEN
ccc	 call prinf(' processing image boxes for I = *',I,1)
	    CALL MPREFLECT(ZCENT,MPOLE1,MPOLE2,NTERMS,
     1              MPR1,MPR2)
	    CALL MKINTR(ILEV,I,INTACT,NINTAC)
ccc	 call prinf(' reflected INTACT list is *',INTACT,NINTAC)
	    ZCENTR = DCONJG(ZCENT)
	    DO 900 J = 1,NINTAC
	       IADR = INTACT(J)
	       INOFF = ITOFST(IADR)
	       NINNBR = ITOFST(IADR+1) - INOFF
	       IF (NINNBR .LE. NLMAX) THEN
c
c---- evaluate far field expansion
c
		  DO 850 K = 1,NINNBR
		     JT = ITRADR(INOFF+K)
		     CALL SLMVALDR(IFLAG1,MPR1,MPR2,ZCENTR,NTERMS,
     1                    ZTARG(JT),GRADW(JT),PHIP(JT),PHIPP(JT))
850               CONTINUE
	       ELSE
c
C     convert far field expansion to local expansion in neighbor box
c
ccc		  t1 = second()
   	          ICOL = 1 + MOD(IADR-1,NSIDE)
 	          IROW = 1 + (IADR-1)/NSIDE
ccc		  XC = -32.0D0+ICOL*H-H/2
ccc		  YC = -32.0D0+IROW*H-H/2
		  XC = ICOL*H-H/2
		  YC = IROW*H-H/2
	          ZCENT2 = DCMPLX(XC,YC)
		  IADLOC = ILOC(IADR)
		  CALL SLSHIFT(MPR1,MPR2,ZCENTR,NTERMS,ZCENT2,
     1                  B1,B2,C,LDC)
		  NSHIFT = NSHIFT+1
		  CALL ADDEXP(B1,B2,LOCEXP1(0,IADLOC),
     1                    LOCEXP2(0,IADLOC),NTERMS)
ccc		  t2 = second()
ccc		  tlota = tlota + (t2-t1)
	       ENDIF
900         CONTINUE
         ENDIF
1000  CONTINUE
      DO 1200 I = 1,NBOXES
	 INOFF = ITOFST(I)
	 NINBOX = ITOFST(I+1) - INOFF
	 IF (NINBOX .LE. NLMAX) THEN
	    GOTO 1200
	 ELSE
   	    ICOL = 1 + MOD(I-1,NSIDE)
 	    IROW = 1 + (I-1)/NSIDE
ccc	    XC = -32.0D0+ICOL*H-H/2
ccc	    YC = -32.0D0+IROW*H-H/2
	    XC = ICOL*H-H/2
	    YC = IROW*H-H/2
	    ZCENT = DCMPLX(XC,YC)
	    JADR = ILOC(I)
	    DO 1100 K = 1,NINBOX
	       JT = ITRADR(INOFF+K)
	       CALL SLLVALDR(IFLAG1,LOCEXP1(0,JADR),LOCEXP2(0,JADR),
     1               NTERMS,ZCENT,ZTARG(JT),
     2               GRADW(JT),PHIP(JT),PHIPP(JT))
1100        CONTINUE
	 ENDIF
1200  CONTINUE
      RETURN
      END
C
C
C*********************************************************************C
      SUBROUTINE SLDIRDR(I1,ZTARG,ZSRC,NATOMS,QA,QB,QC,GW,PHIP,PHIPP)
c
c---- direct calculation
c
c
c     qa/(z_i - z) + z * dconjg(qa/(z_i - z)^2) +
c     dconjg(qb/(z_i - z)) + dconjg(qc/(z_i - z)^2)
c
c
      IMPLICIT REAL *8  (A-H,O-Z)
      COMPLEX *16 ZTARG,ZSRC(1),QA(1),QB(1),QC(1),GW,PHIP,PHIPP
      COMPLEX *16 PHILOC,PHIPLOC,PSILOC,ZDIS,ZDIS2
      COMPLEX *16 PHIPPLOC,ZINV,ZINV2C
C
ccc      call prin2(' in SLDIR, ZSRC = *',ZSRC,2*NATOMS)
ccc      call prin2(' in SLDIR, QA = *',QA,2*NATOMS)
ccc      call prin2(' in SLDIR, QB = *',QB,2*NATOMS)
ccc      call prin2(' in SLDIR, QC = *',QC,2*NATOMS)
      PHILOC = 0.0d0
      PHIPLOC = 0.0D0
      PHIPPLOC = 0.0D0
      PSILOC = 0.0d0
      IF (I1.EQ.0) THEN
         DO 30 J = 1,NATOMS
ccc	     ZDIS = ZSRC(J) - ZTARG
ccc	     ZDIS2 = ZDIS*ZDIS
             ZINV = 1.0D0/(ZSRC(J) - ZTARG)
             ZINV2C = DCONJG(ZINV*ZINV)
ccc             call prin2('ZDIS = *',ZDIS,1)
ccc	     PHILOC = PHILOC + QA(J)/ZDIS
ccc	     PHIPLOC = PHIPLOC + QA(J)/ZDIS2
ccc	     PSILOC = PSILOC + QB(J)/ZDIS
ccc	     PSILOC = PSILOC + QC(J)/ZDIS2
	     GW = GW + QA(J)*ZINV - DCONJG(QB(J)*ZINV) -
     1            (ZTARG*DCONJG(QA(J)) + DCONJG(QC(J)))*ZINV2C
30       CONTINUE
ccc         call prin2('PHILOC = *',PHILOC,2)
ccc         call prin2('PHIPLOC = *',PHIPLOC,2)
ccc         call prin2('PSILOC = *',PSILOC,2)
ccc         GW = GW + PHILOC + ZTARG*DCONJG(PHIPLOC) + DCONJG(PSILOC)
ccc         call prin2('GW = *',GW,2)
      ELSE
         DO 50 J = 1,NATOMS
	     ZDIS = ZSRC(J) - ZTARG
	     ZDIS2 = ZDIS*ZDIS
	     PHILOC = PHILOC + QA(J)/ZDIS
	     PHIPLOC = PHIPLOC + QA(J)/ZDIS2
	     PHIPPLOC = PHIPPLOC + 2*QA(J)/(ZDIS2*ZDIS)
	     PSILOC = PSILOC + QB(J)/ZDIS
	     PSILOC = PSILOC + QC(J)/ZDIS2
50       CONTINUE
         GW = GW + PHILOC - ZTARG*DCONJG(PHIPLOC) - DCONJG(PSILOC)
         PHIP = PHIP + PHIPLOC
         PHIPP = PHIPP + PHIPPLOC
      ENDIF
      RETURN
      END
C
C*****************************************************************
      SUBROUTINE DONN0DR(IFLAG2,ILEV,ZSRC,QA,QB,QC,NATOMS,
     2      ZSRC2,QA2,QB2,QC2,IOFFST,IBOX,IATADR,GRADW)
C
C     The nearest neighbor subroutine of FMM for IFLAG1 = 0 in CADAP.
C     That is, we evaluate fields at source locations themselves,
C     ignoring self.
C
C     Sources have already been assigned to boxes at this level
C     in DOINT0.
C
C     Direct calculations are done now for near neighbors.
C*****************************************************************
      IMPLICIT REAL *8 (A-H,O-Z)
      INTEGER *4 NATOMS,IOFFST(1),IATADR(1)
      INTEGER *4 IBOX(1),NSIDE
      INTEGER *4 NBORS(30),NNBORS
      COMPLEX *16  ZSRC(1),QA(1),QB(1),QC(1)
      COMPLEX *16  ZSRC2(1),QA2(1),QB2(1),QC2(1)
      COMPLEX *16  GRADW(1),ZCENT,ZCENT2
      COMPLEX *16  PHIPT,PHIPPT
      DATA ZERO/0.0d0/
c
c---- process all boxes
c
      NSIDE = 2**ILEV
      NBOXES = NSIDE*NSIDE
      IFLAG1 = 0
ccc      call prinf(' in DONN0, NBOXES = *',NBOXES,1)
      DO 800 I = 1,NBOXES
	 IOFF =   IOFFST(I)
	 NINBOX = IOFFST(I+1) - IOFF
ccc	 call prinf(' I box = *',I,1)
ccc	 call prinf(' NINBOX  = *',NINBOX,1)
	 IF ( NINBOX .LE. 0 ) GOTO 800
         DO K = 1,NINBOX
            ZSRC2(K) = ZSRC(IATADR(IOFF+K))
            QA2(K) = QA(IATADR(IOFF+K))
            QB2(K) = QB(IATADR(IOFF+K))
            QC2(K) = QC(IATADR(IOFF+K))
         ENDDO
ccc	 call prin2(' calling SLSELF with ZSRC2 = *',ZSRC2,2*NINBOX)
ccc	 call prin2(' calling SLSELF with QA2 = *',QA2,2*NINBOX)
ccc	 call Prin2(' calling SLSELF with QB2 = *',QB2,2*NINBOX)
ccc	 call prin2(' calling SLSELF with QC2 = *',QC2,2*NINBOX)
         CALL SLSELFDR(ZSRC2,IOFF,IATADR,NINBOX,QA2,QB2,QC2,
     1               GRADW,PHIPT,PHIPPT)
ccc	 call prin2(' after SLSELF GRADW = *',GRADW,2*NATOMS)
         CALL  MKNBOR(ILEV,I,NBORS,NNBORS)
ccc         call prinf(' nbors are *',nbors,nnbors)
	 DO J = 1,NNBORS
	    IADR = NBORS(J)
	    INOFF = IOFFST(IADR)
	    NINNBR = IOFFST(IADR+1) - INOFF
ccc            call prinf(' ninnbr is *',ninnbr,1)
	    IF (IADR.NE.I) THEN
ccc               call prinf(' different boxes .. interact *',J,0)
c
c        interactions with adjacent boxes
c
	       DO K = 1,NINNBR
	          JT = IATADR(INOFF+K)
		  CALL SLDIRDR(IFLAG1,ZSRC(JT),ZSRC2,NINBOX,QA2,
     1                 QB2,QC2,GRADW(JT),PHIPT,PHIPPT)
	       ENDDO
	    ENDIF
	 ENDDO
800   CONTINUE
      IF (IFLAG2.EQ.1) THEN
         DO 1000 I = 1,NSIDE
	    IOFF =   IOFFST(I)
	    NINBOX = IOFFST(I+1) - IOFF
ccc	    call prinf(' I box = *',I,1)
ccc	    call prinf(' NINBOX  = *',NINBOX,1)
	    IF ( NINBOX .LE. 0 ) GOTO 1000
            DO K = 1,NINBOX
               ZSRC2(K) = ZSRC(IATADR(IOFF+K))
               QA2(K) = QA(IATADR(IOFF+K))
               QB2(K) = QB(IATADR(IOFF+K))
               QC2(K) = QC(IATADR(IOFF+K))
            ENDDO
            CALL  MKNBRRFL(ILEV,I,NBORS,NNBORS)
ccc         call prinf(' nbors are *',nbors,nnbors)
	    DO J = 1,NNBORS
	       IADR = NBORS(J)
	       INOFF = IOFFST(IADR)
	       NINNBR = IOFFST(IADR+1) - INOFF
ccc               call prinf(' ninnbr is *',ninnbr,1)
c
c        interactions with adjacent boxes
c
	       DO K = 1,NINNBR
	          JT = IATADR(INOFF+K)
		  CALL SLDREFDR(IFLAG1,ZSRC(JT),ZSRC2,NINBOX,QA2,
     1                 QB2,QC2,GRADW(JT),PHIPT,PHIPPT)
	       ENDDO
	    ENDDO
1000     CONTINUE
      ENDIF
      RETURN
      END
c
      SUBROUTINE SLSELFDR(ZSRC,IOFF,IATADR,NATOMS,QA,QB,QC,
     1                  GW,PHIP,PHIPP)
c
c---- direct calculation
c
      IMPLICIT REAL *8  (A-H,O-Z)
      INTEGER *4 NATOMS,IATADR(*),IOFF
      COMPLEX *16 ZSRC(1)
      COMPLEX *16 QA(1),QB(1),QC(1),GW(1),PHIP,PHIPP
      COMPLEX *16 PHILOC,PHIPLOC,PSILOC,ZDIS,ZDIS2
      COMPLEX *16 PHIPPLOC
C
ccc      call prin2(' in SLSELF with ZSRC = *',ZSRC,2*NATOMS)
ccc      call prin2(' in SLSELF with QA = *',QA,2*NATOMS)
ccc      call prin2(' in SLSELF with QB = *',QB,2*NATOMS)
ccc      call prin2(' in SLSELF with QC = *',QC,2*NATOMS)
ccc      call prin2(' GW = *',GW,2*NATOMS)
      DO 40 K = 1,NATOMS
ccc	 call prinf(' in SLSELF, K = *',K,1)
         PHILOC = 0.0d0
         PHIP = 0.0d0
         PHIPP = 0.0d0
         PSILOC = 0.0d0
	 KT = IATADR(IOFF+K)
ccc	 call prinf(' in SLSELF, KT = *',KT,1)
         DO 30 J = 1,NATOMS
ccc	     call prinf(' trying atom J = *',J,1)
	     IF (J.EQ.K) GOTO 30
	     ZDIS = ZSRC(J) - ZSRC(K)
ccc	     call prin2(' ZDIS  = *',ZDIS,2)
	     ZDIS2 = ZDIS*ZDIS
	     PHILOC = PHILOC + QA(J)/ZDIS
	     PHIP = PHIP + QA(J)/ZDIS2
	     PHIPP = PHIPP + 2*QA(J)/(ZDIS2*ZDIS)
	     PSILOC = PSILOC + QB(J)/ZDIS
	     PSILOC = PSILOC + QC(J)/ZDIS2
30       CONTINUE
ccc	 call prin2(' PHILOC  = *',PHILOC,2)
ccc	 call prin2(' PHIP  = *',PHIP,2)
ccc	 call prin2(' PSILOC  = *',PSILOC,2)
         GW(KT) = GW(KT) + PHILOC - 
     1         ZSRC(K)*DCONJG(PHIP)-DCONJG(PSILOC)
ccc         call prin2(' in SLSELF GW = *',GW,2*NATOMS)
40    CONTINUE
      RETURN
      END
C
C*****************************************************************
      SUBROUTINE DONN1DR(IFLAG2,ILEV,ZSRC,QA,QB,QC,NATOMS,
     2      ZSRC2,QA2,QB2,QC2,ZTARG,NTARG,IOFFST,IATADR,
     3      ITOFST,IBOX,ITRADR,GRADW,PHIP,PHIPP)
C
C     The nearest neighbor subroutine of FMM for IFLAG1 = 1 in CADAP.
C     That is, we evaluate fields at target locations only.
C
C     Sources have already been assigned to boxes at this level
C     in DOINT1.
C
C     Direct calculations are done now for near neighbors.
C*****************************************************************
      IMPLICIT REAL *8 (A-H,O-Z)
      INTEGER *4 NATOMS,IOFFST(1),IATADR(1)
      INTEGER *4 ITOFST(1),ITRADR(1)
      INTEGER *4 IBOX(1),NSIDE
      INTEGER *4 NBORS(30),NNBORS
      COMPLEX *16  ZSRC(1),QA(1),QB(1),QC(1)
      COMPLEX *16  ZSRC2(1),QA2(1),QB2(1),QC2(1)
      COMPLEX *16  ZTARG(1)
      COMPLEX *16  GRADW(1),PHIP(1),PHIPP(1),ZCENT,ZCENT2
      DATA ZERO/0.0d0/
c
c---- process all boxes
c
      NSIDE = 2**ILEV
      NBOXES = NSIDE*NSIDE
      IFLAG1 = 1
      DO 800 I = 1,NBOXES
	 IOFF =   IOFFST(I)
	 NINBOX = IOFFST(I+1) - IOFF
	 IF ( NINBOX .LE. 0 ) GOTO 800
         DO K = 1,NINBOX
            ZSRC2(K) = ZSRC(IATADR(IOFF+K))
            QA2(K) = QA(IATADR(IOFF+K))
            QB2(K) = QB(IATADR(IOFF+K))
            QC2(K) = QC(IATADR(IOFF+K))
         ENDDO
         CALL  MKNBOR(ILEV,I,NBORS,NNBORS)
	 DO J = 1,NNBORS
	    IADR = NBORS(J)
	    INOFF = ITOFST(IADR)
	    NINNBR = ITOFST(IADR+1) - INOFF
	    DO K = 1,NINNBR
	       JT = ITRADR(INOFF+K)
	       CALL SLDIRDR(IFLAG1,ZTARG(JT),ZSRC2,NINBOX,QA2,
     1                 QB2,QC2,GRADW(JT),PHIP(JT),PHIPP(JT))
	    ENDDO
	 ENDDO
800   CONTINUE
      IF (IFLAG2.EQ.1) THEN
         DO 1000 I = 1,NSIDE
	    IOFF =   IOFFST(I)
	    NINBOX = IOFFST(I+1) - IOFF
ccc	    call prinf(' NSIDE = *',NSIDE,1)
ccc	    call prinf(' I box = *',I,1)
ccc	    call prinf(' NINBOX  = *',NINBOX,1)
	    IF ( NINBOX .LE. 0 ) GOTO 1000
            DO K = 1,NINBOX
               ZSRC2(K) = ZSRC(IATADR(IOFF+K))
               QA2(K) = QA(IATADR(IOFF+K))
               QB2(K) = QB(IATADR(IOFF+K))
               QC2(K) = QC(IATADR(IOFF+K))
            ENDDO
            CALL  MKNBRRFL(ILEV,I,NBORS,NNBORS)
ccc            call prinf(' nbors are *',nbors,nnbors)
	    DO J = 1,NNBORS
	       IADR = NBORS(J)
	       INOFF = ITOFST(IADR)
	       NINNBR = ITOFST(IADR+1) - INOFF
ccc                  call prinf(' ninnbr is *',ninnbr,1)
c
c        interactions with adjacent boxes
c
	       DO K = 1,NINNBR
	          JT = ITRADR(INOFF+K)
		  CALL SLDREF(IFLAG1,ZTARG(JT),ZSRC2,NINBOX,QA2,
     1                 QB2,QC2,GRADW(JT),PHIP(JT),PHIPP(JT))
	       ENDDO
	    ENDDO
1000     CONTINUE
      ENDIF
      RETURN
      END
C*********************************************************************C
      SUBROUTINE SLDREFDR(I1,ZTARG,ZSRC,NATOMS,QA,QB,QC,GW,PHIP,PHIPP)
c
c---- direct calculation of reflected sources.
c
c
c     qa/(z_i - z) + z * dconjg(qa/(z_i - z)^2) +
c     dconjg(qb/(z_i - z)) + dconjg(qc/(z_i - z)^2)
c
c
      IMPLICIT REAL *8  (A-H,O-Z)
      COMPLEX *16 ZTARG,ZSRC(1),QA(1),QB(1),QC(1),GW,PHIP,PHIPP
      COMPLEX *16 PHILOC,PHIPLOC,PSILOC,ZDIS,ZDIS2,ZDIS3
      COMPLEX *16 PHIPPLOC,ZINV,ZINV2C
      COMPLEX *16 DELTA1,DELTA2,DELTA3,GAMMA1,GAMMA2
C
ccc      call prin2(' in SLDREF, ZSRC = *',ZSRC,2*NATOMS)
ccc      call prin2(' in SLDREF, QA = *',QA,2*NATOMS)
ccc      call prin2(' in SLDREF, QB = *',QB,2*NATOMS)
ccc      call prin2(' in SLDREF, QC = *',QC,2*NATOMS)
      PHILOC = 0.0d0
      PHIPLOC = 0.0D0
      PHIPPLOC = 0.0D0
      PSILOC = 0.0d0
      IF (I1.EQ.0) THEN
         DO 30 J = 1,NATOMS
	     ZDIS = DCONJG(ZSRC(J)) - ZTARG
	     ZDIS2 = ZDIS*ZDIS
	     ZDIS3 = ZDIS2*ZDIS
ccc             call prin2('ZDIS = *',ZDIS,1)
	     GAMMA1 = -DCONJG(QB(J)) + DCONJG(QA(J))
	     DELTA1 = -DCONJG(QB(J))
	     GAMMA2 = -DCONJG(QC(J)) - DCONJG(QA(J)*ZSRC(J))
	     DELTA2 = 2*GAMMA2 - DCONJG(ZSRC(J))*GAMMA1
ccc	     DELTA3 = -2*DCONJG(ZSRC(J)*ZSRC(J)*QA(J)) +
ccc     1                2*DCONJG(ZSRC(J)*QC(J))
	     DELTA3 = -2*DCONJG(ZSRC(J))*GAMMA2
	     PHILOC = PHILOC + GAMMA1/ZDIS
	     PHILOC = PHILOC + GAMMA2/ZDIS2
	     PHIPLOC = PHIPLOC + GAMMA1/ZDIS2
	     PHIPLOC = PHIPLOC + 2*GAMMA2/ZDIS3
	     PSILOC = PSILOC + DELTA1/ZDIS
	     PSILOC = PSILOC + DELTA2/ZDIS2
	     PSILOC = PSILOC + DELTA3/ZDIS3
30       CONTINUE
ccc         call prin2('PHILOC = *',PHILOC,2)
ccc         call prin2('PHIPLOC = *',PHIPLOC,2)
ccc         call prin2('PSILOC = *',PSILOC,2)
         GW = GW + PHILOC - ZTARG*DCONJG(PHIPLOC) - DCONJG(PSILOC)
ccc         call prin2('GW = *',GW,2)
      ELSE
         DO 50 J = 1,NATOMS
	     ZDIS = DCONJG(ZSRC(J)) - ZTARG
	     ZDIS2 = ZDIS*ZDIS
	     ZDIS3 = ZDIS2*ZDIS
ccc             call prin2('ZDIS = *',ZDIS,1)
	     GAMMA1 = -DCONJG(QB(J)) + DCONJG(QA(J))
	     DELTA1 = -DCONJG(QB(J))
	     GAMMA2 = -DCONJG(QC(J)) - DCONJG(QA(J)*ZSRC(J))
	     DELTA2 = 2*GAMMA2 - DCONJG(ZSRC(J))*GAMMA1
ccc	     DELTA3 = -2*DCONJG(ZSRC(J)*ZSRC(J)*QA(J)) +
ccc     1                2*DCONJG(ZSRC(J)*QC(J))
	     DELTA3 = -2*DCONJG(ZSRC(J))*GAMMA2
	     PHILOC = PHILOC + GAMMA1/ZDIS
	     PHILOC = PHILOC + GAMMA2/ZDIS2
	     PHIPLOC = PHIPLOC + GAMMA1/ZDIS2
	     PHIPLOC = PHIPLOC + 2*GAMMA2/ZDIS3
	     PHIPPLOC = PHIPPLOC + 2*GAMMA1/ZDIS3
	     PHIPPLOC = PHIPPLOC + 6*GAMMA2/(ZDIS3*ZDIS)
	     PSILOC = PSILOC + DELTA1/ZDIS
	     PSILOC = PSILOC + DELTA2/ZDIS2
	     PSILOC = PSILOC + DELTA3/ZDIS3
50       CONTINUE
         GW = GW + PHILOC - ZTARG*DCONJG(PHIPLOC) - DCONJG(PSILOC)
         PHIP = PHIP + PHIPLOC
         PHIPP = PHIPP + PHIPPLOC
      ENDIF
      RETURN
      END
