C==================================================================
C     Subroutines to calculate indirect constraints
C==================================================================

      SUBROUTINE CALCINDIR
C Master subroutine to call the others for computing indirect constraints.
      IMPLICIT NONE
      CALL INITINDIR
      CALL CALCBSMM
      CALL CALCSPAR
      CALL CALCBSGAM
      RETURN
      END

      SUBROUTINE INITINDIR
C Initialization of experimental inputs for indirect constraints
C (B_s -> mu mu, S parameter)
C OUTPUT: common block INDIREXP
      IMPLICIT NONE
      DOUBLE PRECISION BMMEXP, DBMMEXP, SEXP, DSEXP, 
     .     TEXP, DTEXP, RHOST, MHREF
      COMMON/INDIREXP/BMMEXP, DBMMEXP, SEXP, DSEXP, 
     .     TEXP, DTEXP, RHOST, MHREF

C Experimental value of (averaged time-integrated) Bs -> mu mu and its error,
C from CMS + LHCb, CMS-PAS-BPH-13-007
      BMMEXP = 2.9D-9
      DBMMEXP = 0.7D-9

C Experimental values for the S and T parameters, from PDG EW model rev 2018,
C M. Tanabashi et al. (Particle Data Group), Phys. Rev. D 98, 030001 (2018).
C Computed setting U = 0.
      SEXP = 0.02D0
      DSEXP = 0.07D0
      TEXP = 0.06D0
      DTEXP = 0.06D0
      RHOST = 0.92D0
      MHREF = 125.D0

      RETURN
      END

C------------------------------------------------------------------
C B_s -> mu+ mu-:

      SUBROUTINE CALCBSMM
C Computes the time-averaged BR(B_s -> mu mu) divided by its SM value.
C INPUTS: common blocks SM, PHYSPARAMS, INDIREXP
C OUTPUTS: variable RBSMM in common block INDIR
      IMPLICIT NONE
C Common blocks:
      DOUBLE PRECISION MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
      COMMON/PHYSPARAMS/MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
      DOUBLE PRECISION V, MZ, MW, MTPOLE, MBMB, MCMC, MS, MTAU, MMU,
     .     ALPHAEM,ALSMZ, VCB, VUS, VUB, GAMZ, GAMW
      COMMON/SM/V,MZ,MW,MTPOLE,MBMB,MCMC,MS,MTAU,MMU,
     .     ALPHAEM,ALSMZ,VCB,VUS,VUB,GAMZ,GAMW
      DOUBLE PRECISION BMMEXP, DBMMEXP, SEXP, DSEXP, 
     .     TEXP, DTEXP, RHOST, MHREF
      COMMON/INDIREXP/BMMEXP, DBMMEXP, SEXP, DSEXP, 
     .     TEXP, DTEXP, RHOST, MHREF
      DOUBLE PRECISION RBSMM, SPARAM
      INTEGER BSMMOK, SPAROK, BSGAMLOOSEOK, BSGAMTIGHTOK
      COMMON/INDIR/RBSMM, SPARAM, BSMMOK, SPAROK, 
     .     BSGAMLOOSEOK, BSGAMTIGHTOK
C Local variables:
      DOUBLE PRECISION C10SMREF, MTREF, ALSMZREF
      DOUBLE PRECISION TANTH, XT3, FUNCX
      DOUBLE PRECISION C10SM, C10GM
      DOUBLE PRECISION BMMSMREF, DBMMSMREF, BMMSM, DBMMSM, 
     .     BMMGM, DBMMGM, BMMDIFF, DBMM
      DOUBLE PRECISION MT
C Function to be called:
      DOUBLE PRECISION RUNMT

      MT = RUNMT(MTPOLE)

C Reference values for SM calculation of C10, from Li, Lu & Pich, 1404.5865
      C10SMREF = -0.9380D0
      MTREF = 173.1D0
      ALSMZREF = 0.1184D0
      C10SM = C10SMREF*(MTPOLE/MTREF)**1.53 * (ALSMZ/ALSMZREF)**(-0.09)

      TANTH = 2.D0*DSQRT(2.D0)*VCHI/VPHI
      XT3 = MT**2/MH3**2
      IF (DABS(XT3-1.D0).GT.1D-4) THEN
         FUNCX = (XT3/(1.D0-XT3) + XT3*DLOG(XT3)/(1.D0-XT3)**2)
      ELSE
         FUNCX = -0.5D0 - (XT3-1.D0)/6.D0
      ENDIF
      C10GM = C10SM + TANTH**2/8.D0 * MT**2/MW**2 * FUNCX

      RBSMM = (C10GM/C10SM)**2

C Reference values for SM calculation of time-averaged BR(B_s -> mu mu)
C and its theory uncertainty, from Li, Lu & Pich, 1404.5865
      BMMSMREF = 3.67D-9
      DBMMSMREF = 0.25D-9
C Scale reference values according to m_t and alpha_s(M_Z):
      BMMSM = BMMSMREF * (C10SM/C10SMREF)**2
      DBMMSM = DBMMSMREF * (C10SM/C10SMREF)**2
C BR(B_s -> mu mu) in GM model and its theory uncertainty:
      BMMGM = BMMSM * RBSMM
      DBMMGM = DBMMSM * RBSMM

C Expt measurement minus GM prediction:
      BMMDIFF = BMMEXP - BMMGM
C Combined theory + expt uncertainty:
      DBMM = DSQRT(DBMMGM**2 + DBMMEXP**2)

C Check consistency with experimental measurement and set flag (2 sigma):
      IF (DABS(BMMDIFF).LE.2.D0*DBMM) THEN
         BSMMOK = 1
      ELSE
         BSMMOK = 0
      ENDIF

      RETURN
      END

C------------------------------------------------------------------
C S parameter:

      SUBROUTINE CALCSPAR
C Computes the S parameter (relative to the SM with Higgs mass = MH).
C INPUTS: common blocks SM, PHYSPARAMS, INDIREXP
C OUTPUTS: variable SPARAM in common block INDIR
      IMPLICIT NONE
C Common blocks:
      DOUBLE PRECISION MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
      COMMON/PHYSPARAMS/MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
      DOUBLE PRECISION V, MZ, MW, MTPOLE, MBMB, MCMC, MS, MTAU, MMU,
     .     ALPHAEM,ALSMZ, VCB, VUS, VUB, GAMZ, GAMW
      COMMON/SM/V,MZ,MW,MTPOLE,MBMB,MCMC,MS,MTAU,MMU,
     .     ALPHAEM,ALSMZ,VCB,VUS,VUB,GAMZ,GAMW
      DOUBLE PRECISION BMMEXP, DBMMEXP, SEXP, DSEXP, 
     .     TEXP, DTEXP, RHOST, MHREF
      COMMON/INDIREXP/BMMEXP, DBMMEXP, SEXP, DSEXP, 
     .     TEXP, DTEXP, RHOST, MHREF
      DOUBLE PRECISION RBSMM, SPARAM
      INTEGER BSMMOK, SPAROK, BSGAMLOOSEOK, BSGAMTIGHTOK
      COMMON/INDIR/RBSMM, SPARAM, BSMMOK, SPAROK, 
     .     BSGAMLOOSEOK, BSGAMTIGHTOK
C Local variables:
      DOUBLE PRECISION EE, SW, CW, PI
      DOUBLE PRECISION GZHLH3, GZHHH3, GZH5H3, GZH5PH3M
      DOUBLE PRECISION GZZHL, GZZHSM, GZZHH, GZZH5, GZWH5
      DOUBLE PRECISION TMIN, CHISQ 
C Functions:
      DOUBLE PRECISION SF1, SF3

      PI = 4.D0*DATAN(1.D0)
      EE = DSQRT(4.D0*PI*ALPHAEM)
      CW = MW/MZ
      SW = DSQRT(1.D0 - CW**2)

C Couplings (all will be mod-squared, so drop the i and overall sign)
      GZHLH3 = DSQRT(2.D0/3.D0)*EE/(SW*CW)
     .     * (DSQRT(3.D0)*VCHI/V*DCOS(ALPHA) + VPHI/V*DSIN(ALPHA))
      GZHHH3 = DSQRT(2.D0/3.D0)*EE/(SW*CW)
     .     * (VPHI/V*DCOS(ALPHA) - DSQRT(3.D0)*VCHI/V*DSIN(ALPHA))
      GZH5H3 = DSQRT(1.D0/3.D0)*EE/(SW*CW) * VPHI/V
      GZH5PH3M = EE/(2.D0*SW*CW) * VPHI/V

      GZZHL = -EE**2/(6.D0*SW**2*CW**2)
     .     * (8.D0*DSQRT(3.D0)*VCHI*DSIN(ALPHA) 
     .     - 3.D0*VPHI*DCOS(ALPHA))
      GZZHH = EE**2/(6.D0*SW**2*CW**2) 
     .     * (8.D0*DSQRT(3.D0)*VCHI*DCOS(ALPHA)
     .     + 3.D0*VPHI*DSIN(ALPHA))
      GZZH5 = -DSQRT(8.D0/3.D0)*EE**2/(SW**2*CW**2) * VCHI
      GZWH5 = -DSQRT(2.D0)*EE**2/(SW**2*CW) * VCHI

      GZZHSM = EE**2*V/(2.D0*SW**2*CW**2)

C output for debugging
c      PRINT * 
c      PRINT *, "--- debugging of S parameter calculation ---"
c      PRINT *, "GZHLH3 = ", GZHLH3
c      PRINT *, "GZHHH3 = ", GZHHH3
c      PRINT *, "GZH5H3 = ", GZH5H3
c      PRINT *, "GZH5PH3M = ", GZH5PH3M
c      PRINT *, "GZZHL = ", GZZHL
c      PRINT *, "GZZHH = ", GZZHH
c      PRINT *, "ALPHAEM = ", ALPHAEM
c      PRINT *, "MW    = ", MW
c      PRINT *, "GZZH5 = ", GZZH5
c      PRINT *, "EE    = ", EE
c      PRINT *, "SW    = ", SW
c      PRINT *, "CW    = ", CW
c      PRINT *, "VCHI  = ", VCHI
c      PRINT *, "GZWH5 = ", GZWH5
c      PRINT *, "MZ    = ", MZ
c      PRINT *, "MHREF = ", MHREF
c      PRINT *, "f1(MZ,MHREF) = ", SF1(MZ,MHREF)
c      PRINT *, "f3(MZ,MHREF) = ", SF3(MZ,MHREF)
c      PRINT *, "GZZHSM = ", GZZHSM
c      PRINT *
C end output for debugging

      SPARAM = SW**2*CW**2/EE**2/PI *
     .     ( -EE**2/(12.D0*SW**2*CW**2)*(DLOG(MH3**2)+5.D0*DLOG(MH5**2))
     .     + 2.D0*GZHLH3**2*SF1(MHL,MH3)
     .     + 2.D0*GZHHH3**2*SF1(MHH,MH3)
     .     + 2.D0*(GZH5H3**2 + 2.D0*GZH5PH3M**2)*SF1(MH5,MH3)
     .     + GZZHL**2*(SF1(MZ,MHL)/2.D0/MZ**2 - SF3(MZ,MHL))
     .     - GZZHSM**2*(SF1(MZ,MHREF)/2.D0/MZ**2 - SF3(MZ,MHREF))
     .     + GZZHH**2*(SF1(MZ,MHH)/2.D0/MZ**2 - SF3(MZ,MHH))
     .     + GZZH5**2*(SF1(MZ,MH5)/2.D0/MZ**2 - SF3(MZ,MH5))
     .     + 2.D0*GZWH5**2*(SF1(MW,MH5)/2.D0/MW**2 - SF3(MW,MH5)) )

C Marginalize over T parameter:
      TMIN = TEXP + RHOST*(SPARAM - SEXP)*DTEXP/DSEXP

      CHISQ = ((SPARAM-SEXP)**2/DSEXP**2 + (TMIN-TEXP)**2/DTEXP**2
     .     - 2.D0*RHOST*(SPARAM-SEXP)*(TMIN-TEXP)/DSEXP/DTEXP)
     .     /(1.D0-RHOST**2)

C Check consistency with experimental measurements and set flag (2 sigma):
      IF (CHISQ.LE.4.D0) THEN
         SPAROK = 1
      ELSE
         SPAROK = 0
      ENDIF

      RETURN
      END

      DOUBLE PRECISION FUNCTION SF1(M1,M2)
      IMPLICIT NONE
      DOUBLE PRECISION M1, M2
      DOUBLE PRECISION EPS
      EPS = M2**2/M1**2 - 1.D0
      IF (DABS(EPS).GT.1.D-4) THEN
         SF1 = (5.D0*(M2**6-M1**6) 
     .        + 27.D0*(M1**4*M2**2-M1**2*M2**4)
     .        + 12.D0*(M1**6-3.D0*M1**4*M2**2)*DLOG(M1)
     .        + 12.D0*(3.D0*M1**2*M2**4-M2**6)*DLOG(M2))
     .        /36.D0/(M1**2-M2**2)**3
      ELSE
         SF1 = 1.D0/6.D0*DLOG(M1**2) + EPS/12.D0
      ENDIF
      END

      DOUBLE PRECISION FUNCTION SF3(M1,M2)
      IMPLICIT NONE
      DOUBLE PRECISION M1, M2
      DOUBLE PRECISION EPS
      EPS = M2**2/M1**2 - 1.D0
      IF (DABS(EPS).GT.1.D-4) THEN
         SF3 = (M1**4 - M2**4 
     .        + 2.D0*M1**2*M2**2*(DLOG(M2**2) - DLOG(M1**2)))
     .        /2.D0/(M1**2 - M2**2)**3
      ELSE
         SF3 = 1.D0/6.D0/M1**2 - EPS/12.D0/M1**2
      ENDIF
      END

C--------------------------------------------------------------------
C BR(b -> s gamma):

      SUBROUTINE CALCBSGAM
C This is a **first implementation** of the b -> s gamma constraint,
C to be improved in later versions.
C For now, this subroutine reads in a numerical file containing the
C desired limit curve from the b -> s gamma constraint in the MH3 - VCHI
C plane.
C The limit curve file was generated by computing BR(b -> s gamma) using 
C SuperIso v3.3 for the Type-1 2HDM and taking MH3 --> MH+ and 
C cot(theta_H) --> tan(beta).
C The maximum allowed value of VCHI for a given value of MH3 is determined
C using linear interpolation from the data file.  The flag BSGAMOK is then
C set based on whether VCHI is above (=0) or below (=1) the limit.
      IMPLICIT NONE
C Common blocks:
      DOUBLE PRECISION MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
      COMMON/PHYSPARAMS/MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
      DOUBLE PRECISION RBSMM, SPARAM
      INTEGER BSMMOK, SPAROK, BSGAMLOOSEOK, BSGAMTIGHTOK
      COMMON/INDIR/RBSMM, SPARAM, BSMMOK, SPAROK, 
     .     BSGAMLOOSEOK, BSGAMTIGHTOK
C Local variables:
      DOUBLE PRECISION MASS(200), VEV(200)
      DOUBLE PRECISION VCHILIMIT
      INTEGER I, NLOOSE, NTIGHT, NPOINTS

      NLOOSE = 138
      NTIGHT = 138

C Extract the "loose" b -> s gamma constraint:
      NPOINTS = NLOOSE
      OPEN (UNIT = 99, FILE = 'src/bsgloose.data', STATUS = 'OLD')
      READ (99,*)
      READ (99,*)
      DO 10 I = 1,NPOINTS
         READ (99,*) MASS(I), VEV(I)
 10   CONTINUE
      CLOSE (99)

      IF (MH3.LE.MASS(1)) THEN
         VCHILIMIT = VEV(1)
      ELSE IF (MH3.GE.MASS(NPOINTS)) THEN
         VCHILIMIT = VEV(NPOINTS)
      ELSE
         I = 1
 20      IF (MASS(I+1).GE.MH3) THEN
C do the linear interpolation
            CALL LININTERP(MASS(I),VEV(I),MASS(I+1),VEV(I+1),
     .           MH3,VCHILIMIT)
         ELSE
            I = I+1
            GOTO 20
         ENDIF
      ENDIF

      IF (DABS(VCHI).LE.VCHILIMIT) THEN
         BSGAMLOOSEOK = 1
      ELSE
         BSGAMLOOSEOK = 0
      ENDIF

C Extract the "tight" b -> s gamma constraint:
      NPOINTS = NTIGHT
      OPEN (UNIT = 99, FILE = 'src/bsgtight.data', STATUS = 'OLD')
      READ (99,*)
      READ (99,*)
      DO 11 I = 1,NPOINTS
         READ (99,*) MASS(I), VEV(I)
 11   CONTINUE
      CLOSE (99)

      IF (MH3.LE.MASS(1)) THEN
         VCHILIMIT = VEV(1)
      ELSE IF (MH3.GE.MASS(NPOINTS)) THEN
         VCHILIMIT = VEV(NPOINTS)
      ELSE
         I = 1
 21      IF (MASS(I+1).GE.MH3) THEN
C do the linear interpolation
            CALL LININTERP(MASS(I),VEV(I),MASS(I+1),VEV(I+1),
     .           MH3,VCHILIMIT)
         ELSE
            I = I+1
            GOTO 21
         ENDIF
      ENDIF

      IF (DABS(VCHI).LE.VCHILIMIT) THEN
         BSGAMTIGHTOK = 1
      ELSE
         BSGAMTIGHTOK = 0
      ENDIF

      END

