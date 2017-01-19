ccc      IMPLICIT REAL *8 (A-H,O-Z)
ccc      PARAMETER (NTERMS = 40)  
ccc      COMPLEX *16 ZC(10),ZB(10),ZK(10)
ccc      COMPLEX *16 ZTARG,PHI,PHIP,PSI
ccc      real *8 drand
ccc      external drand
cccC
cccC
ccc      NATOMS = 10
ccc      DO J = 1,NATOMS
ccc         ZK(J) = DCMPLX(drand(0.0),2+drand(0.0))
ccc         ZC(J) = DCMPLX(drand(0.0),drand(0.0))
ccc         ZB(J) = DCMPLX(drand(0.0),drand(0.0))
ccc      ENDDO
ccc      ZTARG = DCMPLX(-0.4D0,0.0D0)
ccc      CALL SLDIRSING(ZTARG,NATOMS,ZK,ZC,ZB,PHI,PHIP,PSI)
cccc
ccc      write(6,*)' GW = ',PHI + ZTARG*DCONJG(PHIP) + DCONJG(PSI)
ccc      write(13,*)' GW = ',PHI + ZTARG*DCONJG(PHIP) + DCONJG(PSI)
cccC
ccc      CALL REFSING(ZTARG,NATOMS,ZK,ZC,ZB,PHI,PHIP,PSI)
ccc      write(6,*)' GW = ',PHI + ZTARG*DCONJG(PHIP) + DCONJG(PSI)
ccc      write(13,*)' GW = ',PHI + ZTARG*DCONJG(PHIP) + DCONJG(PSI)
ccc      STOP
ccc      END
C
C*********************************************************************C
      SUBROUTINE REFSING(ZTARG,NATOMS,ZK,ZC,ZB,PHI,PHIP,PSI)
c
c---- direct calculation of reflected sources.
c
c     original sources are given by:
c
c     PHI = \sum ZC(j) \log( z - ZK(j))
c
c     PSI = \sum dconjg(ZC(j)) \log( z - ZK(j)) +
c           \sum Zb(j)) /(z - ZK(j))
c
c
c     qa/(z_i - z) + z * dconjg(qa/(z_i - z)^2) +
c     dconjg(qb/(z_i - z)) + dconjg(qc/(z_i - z)^2)
c
c
      IMPLICIT REAL *8  (A-H,O-Z)
      COMPLEX *16 ZTARG,ZK(1),ZC(1),ZB(1),PHI,PHIP,PSI
      COMPLEX *16 ZDIS,ZDIS2
      COMPLEX *16 ZINV,ZINV2C,ZTEMP,ZTEMP2
      COMPLEX *16 GAMMA1,DELTA1,DELTA2
      COMPLEX *16 ALPHA1,BETA1,BETA2
C
ccc      PHI = 0.0D0
ccc      PHIP = 0.0D0
ccc      PSI = 0.0D0
      DO 30 J = 1,NATOMS
	 ZDIS = ZTARG - DCONJG(ZK(J))
	 X = DREAL(ZK(J))
	 Y = DIMAG(ZK(J))
	 ZDIS2 = ZDIS*ZDIS
ccc             call prin2('ZDIS = *',ZDIS,1)
c
c        reflect ZB pole terms
c
	 GAMMA1 = -DCONJG(ZB(J))
	 DELTA1 = GAMMA1
	 DELTA2 = DCONJG(ZK(J))*GAMMA1
	 PHI = PHI + GAMMA1/ZDIS
	 PHIP = PHIP - GAMMA1/ZDIS2
	 PSI = PSI + DELTA1/ZDIS
	 PSI = PSI + DELTA2/ZDIS2
c
c        reflect log sources
c
	 PHI = PHI - ZC(J)*CDLOG(ZDIS)
	 PHIP = PHIP - ZC(J)/ZDIS
	 PSI = PSI - DCONJG(ZC(J))*CDLOG(ZDIS)
ccc	 ALPHA1 = DCONJG(ZC(J))*DCMPLX(0.0D0,2*Y)
	 ALPHA1 = DCONJG(ZC(J))*(ZK(J) - DCONJG(ZK(J)))
	 BETA1 = ZC(J)*DCONJG(ZK(J)) + ALPHA1
	 BETA2 = ALPHA1*DCONJG(ZK(J))
	 PHI = PHI + ALPHA1/ZDIS
	 PHIP = PHIP - ALPHA1/ZDIS2
	 PSI = PSI + BETA1/ZDIS
	 PSI = PSI + BETA2/ZDIS2
30    CONTINUE
      RETURN
      END
C*********************************************************************C
      SUBROUTINE SLDIRSING(ZTARG,NATOMS,ZK,ZC,ZB,PHI,PHIP,PSI)
c
c---- direct calculation
c
c     PHI = \sum ZC(j) \log( z - ZK(j))
c
c     PSI = \sum dconjg(ZC(j)) \log( z - ZK(j)) +
c           \sum Zb(j)) \log( z - ZK(j)) +
c
c
      IMPLICIT REAL *8  (A-H,O-Z)
      COMPLEX *16 ZTARG,ZK(1),ZC(1),ZB(1),PHI,PHIP,PSI
      COMPLEX *16 PHILOC,PHIPLOC,PSILOC,ZDIS,ZDIS2
      COMPLEX *16 ZINV,ZINV2C
C
ccc      call prin2(' in SLDIR, ZSRC = *',ZSRC,2*NATOMS)
ccc      call prin2(' in SLDIR, QA = *',QA,2*NATOMS)
ccc      call prin2(' in SLDIR, QB = *',QB,2*NATOMS)
ccc      call prin2(' in SLDIR, QC = *',QC,2*NATOMS)
      PHI = 0.0d0
      PHIP = 0.0D0
      PSI = 0.0d0
      DO 30 J = 1,NATOMS
	 ZDIS = ZTARG - ZK(J)
	 PHI = PHI + ZC(J)*CDLOG(ZDIS)
	 PHIP = PHIP + ZC(J)/ZDIS
	 PSI = PSI + DCONJG(ZC(J))*CDLOG(ZDIS)
	 PSI = PSI - ZC(J)*DCONJG(ZK(J))/ZDIS
	 PSI = PSI + ZB(J)/ZDIS
30    CONTINUE
      RETURN
      END
C

