C========================================================================
C       Subroutines to calculate direct collider constraints
C========================================================================
C CALCWWJJ: constraint on H5++ from like-sign WWjj, ATLAS Run 1 1405.6241, recast 
C           for VBF H5++ -> W+W+ in 1407.5053 (numbers courtesy of Cheng-Wei Chiang).
C           Sets the flag WWJJOK.
C CALCLSDM: constraint on H5++ from ATLAS Run 1 anomalous like-sign dimuon production
C           1412.0237, recast for Drell-Yan H++H-- in Higgs Triplet model in 
C           1412.7603, reinterpreted in Georgi-Machacek model in 1502.01275.
C           Sets the flag LSDMOK.
C CALCH5PP: direct experimental constraint on VBF H5++ -> W+W+ from CMS Run 2 
C           (137 fb-1) 2104.04762
C           Sets the flag H5PPOK.
C CALCATLAS8TEVGAGA: constraint from ATLAS pp -> H50 -> gam gam from
C           ATLAS Run 1 (8 TeV, 20.3 fb-1) 1407.6583, using Drell-Yan production
C           of H50 and H5+ or H5-.
C           Sets the flag ATLAS8TEVGAGAOK.
C CALCATLAS13TEVDYGAGA: constraint from ATLAS pp -> H50 -> gam gam from ATLAS
C           Run 2 (37 fb-1) 1707.04147, using VBF production of H50.
C           Sets the flag ATLAS13TEVDYGAGAOK.
C CALCATLAS13TEVVBFGAGA: constraint from ATLAS pp -> H50 -> gam gam from ATLAS
C           Run 2 (37 fb-1) 1707.04147, using VBF production of H50.
C           Sets the flag ATLAS13TEVVBFGAGAOK.
C CALCDYHPP: constraint from ATLAS Drell-Yan pp -> H5++H5-- -> 2,3,4 leptons
C           from ATLAS Run 2 (139 fb-1) 2101.11961
C           Sets the flag DYHPPOK.
C===========================================================================

      SUBROUTINE CALCWWJJ
C Implements the 2-sigma constraint on the m5-vchi plane from the ATLAS
C Run-1 like-sign WW cross section in VBF, as computed in Chiang, Kanemura
C and Yagyu, 1407.5053.  Data file courtesy of Cheng-Wei Chiang.
C When decays of H5 --> V H3 are open, the limit is on sigma*BR(H5->WW),
C so for these points we actually calculate BR(H5++ --> W+W+).
C INPUTS: common block PHYSPARAM
C OUTPUT: common block CONSTR, flag WWJJOK
      IMPLICIT NONE
C Common blocks:
      DOUBLE PRECISION MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
      COMMON/PHYSPARAMS/MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
      DOUBLE PRECISION H5PPBRWW, H5PPBRWH3, H5PPBRH3P, H5PPWDTH
      COMMON/H5PPBRS/H5PPBRWW, H5PPBRWH3, H5PPBRH3P, H5PPWDTH
      INTEGER WWJJOK, LSDMOK, H5PPOK, ATLAS8TEVGAGAOK
      INTEGER ATLAS13TEVDYGAGAOK, ATLAS13TEVVBFGAGAOK
      INTEGER DYHPPOK
      COMMON/CONSTR/WWJJOK, LSDMOK, H5PPOK, ATLAS8TEVGAGAOK,
     .     ATLAS13TEVDYGAGAOK, ATLAS13TEVVBFGAGAOK,
     .     DYHPPOK
C Local variables:
      DOUBLE PRECISION MASS(10), VEV(10)
      DOUBLE PRECISION VCHILIMIT
      INTEGER I, NPOINTS

      NPOINTS = 10
      OPEN (UNIT = 99, FILE = 'src/wwjj.data', STATUS = 'OLD')
      DO 9 I = 1,4
         READ (99, *)
 9    CONTINUE
      DO 10 I = 1, NPOINTS
         READ (99, *) MASS(I), VEV(I)
 10   CONTINUE
      CLOSE (99)

      WWJJOK = 1
C Check if we are within the range of the constraint:
      IF (MH5.GE.MASS(1).AND.MH5.LE.MASS(NPOINTS)) THEN
         I = 1
 20      IF (MASS(I+1).GE.MH5) THEN
C Do the linear interpolation
            CALL LININTERP(MASS(I),VEV(I),MASS(I+1),VEV(I+1),
     .           MH5,VCHILIMIT)
            IF (DABS(VCHI).GT.VCHILIMIT) THEN
C We are in the danger zone, though BR(H5++ -> WW) < 1 might save us.
C Check whether H5 -> H3 V decays are open; if so, calculate BR(H5++ -> WW).
               IF (MH3.GE.MH5) THEN
                  WWJJOK = 0
               ELSE
                  CALL H5PPDECAYS
                  IF (VCHI**2*H5PPBRWW.GT.VCHILIMIT**2) THEN
                     WWJJOK = 0
                  ENDIF
               ENDIF
            ENDIF
         ELSE
            I = I+1
            GOTO 20
         ENDIF
      ENDIF
      RETURN
      END

C===========================================================================

      SUBROUTINE CALCH5PP
C Implements the 95%CL upper bound on VBF H5++/-- -> W+/- W+/- -> likesign 
C dileptons from CMS Run 2 (13 TeV) 137 fb-1 results, 2104.04762.
C Upper limit on sH as a function of m5 is given assuming BR(H5++ -> W+W+) = 100%.
C The cross section is proportional to sH^2, so we scale the exclusion limit
C appropriately when BR(H5++ -> W+W+) < 100%.
C INPUTS: common blocks PHYSPARAMS, H5PPBRS
C OUTPUT: common block CONSTR, flag H5PPOK
      IMPLICIT NONE
C Common blocks:
      DOUBLE PRECISION MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
      COMMON/PHYSPARAMS/MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
      DOUBLE PRECISION H5PPBRWW, H5PPBRWH3, H5PPBRH3P, H5PPWDTH
      COMMON/H5PPBRS/H5PPBRWW, H5PPBRWH3, H5PPBRH3P, H5PPWDTH
      INTEGER WWJJOK, LSDMOK, H5PPOK, ATLAS8TEVGAGAOK
      INTEGER ATLAS13TEVDYGAGAOK, ATLAS13TEVVBFGAGAOK
      INTEGER DYHPPOK
      COMMON/CONSTR/WWJJOK, LSDMOK, H5PPOK, ATLAS8TEVGAGAOK,
     .     ATLAS13TEVDYGAGAOK, ATLAS13TEVVBFGAGAOK,
     .     DYHPPOK
C Local variables:
      DOUBLE PRECISION MASS(12), EXPSH(12)
      DOUBLE PRECISION SHLIMIT, SH, V
      INTEGER I, NPOINTS

      NPOINTS = 12
      OPEN (UNIT = 99, FILE = 'src/cms-VBFHpp-run2-137fb.data', 
     .     STATUS = 'OLD')
      DO 9 I = 1,7
         READ (99, *)
 9    CONTINUE
      DO 10 I = 1, NPOINTS
         READ (99, *) MASS(I), EXPSH(I)
 10   CONTINUE
      CLOSE (99)

      H5PPOK = 1
C Check if we are within the range of the constraint:
      IF (MH5.GE.MASS(1).AND.MH5.LE.MASS(NPOINTS)) THEN
         I = 1
 20      IF (MASS(I+1).GE.MH5) THEN
C Do the linear interpolation
            CALL LININTERP(MASS(I),EXPSH(I),MASS(I+1),EXPSH(I+1),
     .           MH5,SHLIMIT)
            V = DSQRT(VPHI**2 + 8.D0*VCHI**2)
            SH = DSQRT(8.D0)*VCHI/V
            IF (DABS(SH).GT.SHLIMIT) THEN
C We are in the danger zone, though BR(H5++ -> WW) < 1 might save us.
C Check whether H5 -> H3 V decays are open; if so, calculate BR(H5++ -> WW).
               IF (MH3.GE.MH5) THEN
                  H5PPOK = 0
               ELSE
                  CALL H5PPDECAYS
                  IF (SH**2*H5PPBRWW.GT.SHLIMIT**2) THEN
                     H5PPOK = 0
                  ENDIF
               ENDIF
            ENDIF
         ELSE
            I = I+1
            GOTO 20
         ENDIF
      ENDIF
      RETURN
      END

C===========================================================================

      SUBROUTINE CALCLSDM
C Implements the 95%CL upper bound on the fiducial cross section for anomalous
C like-sign dimuon production from ATLAS Run 1 1412.0237.  GM model cross section
C is the lower curve in Fig. 4 of 1502.01275 (-5% thy uncert).
C We scale the exclusion limit appropriately when BR(H5++ -> W+W+) < 100%.
C This requires implementing the full assembly of cross sections, efficiencies,
C and BRs from Kanemura, Kikuchi, Yagyu, & Yokoya 1412.7603.
C Cross sections are NLO QCD, 8 TeV pp.
C INPUTS: common block PHYSPARAM
C OUTPUT: common block CONSTR, flag LSDMOK
      IMPLICIT NONE
C Common blocks:
      DOUBLE PRECISION MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
      COMMON/PHYSPARAMS/MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
      DOUBLE PRECISION H5PPBRWW, H5PPBRWH3, H5PPBRH3P, H5PPWDTH
      COMMON/H5PPBRS/H5PPBRWW, H5PPBRWH3, H5PPBRH3P, H5PPWDTH
      INTEGER WWJJOK, LSDMOK, H5PPOK, ATLAS8TEVGAGAOK
      INTEGER ATLAS13TEVDYGAGAOK, ATLAS13TEVVBFGAGAOK
      INTEGER DYHPPOK
      COMMON/CONSTR/WWJJOK, LSDMOK, H5PPOK, ATLAS8TEVGAGAOK,
     .     ATLAS13TEVDYGAGAOK, ATLAS13TEVVBFGAGAOK,
     .     DYHPPOK
C Local variables:
      DOUBLE PRECISION ATLASLIMIT
      DOUBLE PRECISION M5(6)
      DOUBLE PRECISION SIGHTM22(6), SIGHTM21(6), SIGHTM12(6)
      DOUBLE PRECISION SIGGM22(6), SIGGM21(6), SIGGM12(6)
      DOUBLE PRECISION BRMM(6)
      DOUBLE PRECISION EPS22(6), EPS21(6), EPS12(6)
      DOUBLE PRECISION SIGGM22PT, SIGGM21PT, SIGGM12PT, BRMMPT,
     .     EPS22PT, EPS21PT, EPS12PT
      DOUBLE PRECISION FIDXSPT
      INTEGER I, NPOINTS

C ATALS limit is in fb:
      ATLASLIMIT = 16.D0
      NPOINTS = 6
      LSDMOK = 1

C Fill arrays from Table 1 of arXiv:1412.7603 (note: cross sections are in pb)
      M5 = (/ 50., 60., 70., 80., 90., 100. /)
      SIGHTM22 = (/ 8.52, 3.57, 1.93, 1.16, 0.744, 0.501 /)
      SIGHTM21 = (/ 10.6, 4.47, 2.36, 1.40, 0.891, 0.598 /)
      SIGHTM12 = (/ 6.71, 2.73, 1.40, 0.803, 0.498, 0.326 /)
      BRMM = (/ 2.22, 2.21, 2.19, 2.16, 1.98, 1.61 /)
      EPS22 = (/ 5.1, 9.9, 16., 21., 23., 23. /)
      EPS21 = (/ 4.9, 9.9, 15., 21., 22., 23. /)
      EPS12 = (/ 4.7, 9.7, 15., 21., 23., 22. /)

C Convert BRMM and efficiencies from percents to actual numbers
      DO 10 I=1,6
         BRMM(I) = BRMM(I)*0.01D0
         EPS22(I) = EPS22(I)*0.01D0
         EPS21(I) = EPS21(I)*0.01D0
         EPS12(I) = EPS12(I)*0.01D0
 10   CONTINUE

C Apply scalings for GM subprocess cross sections, and convert to fb
      DO 11 I=1,6
         SIGGM22(I) = SIGHTM22(I)*1D3
         SIGGM21(I) = 0.5D0*SIGHTM21(I)*1D3
         SIGGM12(I) = 0.5D0*SIGHTM12(I)*1D3
 11   CONTINUE

C Check if we are within the range of the constraint:
      IF (MH5.GE.M5(1).AND.MH5.LE.M5(NPOINTS)) THEN
         I = 1
 20      IF (M5(I+1).GE.MH5) THEN
C Use linear interpolation to extract all the relevant numbers at MH5
            CALL LININTERP(M5(I),SIGGM22(I),M5(I+1),SIGGM22(I+1),
     .           MH5,SIGGM22PT)
            CALL LININTERP(M5(I),SIGGM21(I),M5(I+1),SIGGM21(I+1),
     .           MH5,SIGGM21PT)
            CALL LININTERP(M5(I),SIGGM12(I),M5(I+1),SIGGM12(I+1),
     .           MH5,SIGGM12PT)
            CALL LININTERP(M5(I),BRMM(I),M5(I+1),BRMM(I+1),
     .           MH5,BRMMPT)
            CALL LININTERP(M5(I),EPS22(I),M5(I+1),EPS22(I+1),
     .           MH5,EPS22PT)
            CALL LININTERP(M5(I),EPS21(I),M5(I+1),EPS21(I+1),
     .           MH5,EPS21PT)
            CALL LININTERP(M5(I),EPS12(I),M5(I+1),EPS12(I+1),
     .           MH5,EPS12PT)
C To avoid calling H5PPDECAYS unnecessarily, first check the total fiducial 
C cross section assuming BR(H5++ -> WW) = 100%.  
C Also, apply the conservative -5% theory uncertainty.
            FIDXSPT = 0.95D0*(SIGGM22PT*(2.D0*BRMMPT*EPS22PT
     .           - BRMMPT**2*EPS22PT**2)
     .           + SIGGM21PT*BRMMPT*EPS21PT
     .           + SIGGM12PT*BRMMPT*EPS12PT)
            IF (FIDXSPT.GT.ATLASLIMIT) THEN
C We are in the danger zone, though BR(H5++ -> WW) < 1 might save us.
C Check whether H5 -> H3 V decays are open; if so, calculate BR(H5++ -> WW).
               IF (MH3.GE.MH5) THEN
                  LSDMOK = 0
               ELSE
                  CALL H5PPDECAYS
C Now recompute the fiducial cross section taking into account BR(H5++ -> WW).
C Again apply the conservative -5% theory uncertainty.
                  FIDXSPT = 0.95D0
     .                 *(SIGGM22PT*(2.D0*H5PPBRWW*BRMMPT*EPS22PT
     .                 - H5PPBRWW**2*BRMMPT**2*EPS22PT**2)
     .                 + SIGGM21PT*H5PPBRWW*BRMMPT*EPS21PT
     .                 + SIGGM12PT*H5PPBRWW*BRMMPT*EPS12PT)
                  IF (FIDXSPT.GT.ATLASLIMIT) THEN
                     LSDMOK = 0
                  ENDIF
               ENDIF
            ENDIF
         ELSE
            I = I+1
            GOTO 20
         ENDIF
      ENDIF
      RETURN
      END

C===========================================================================

      SUBROUTINE CALCATLAS8TEVGAGA
C By Yongcheng Wu 2018
C Implements the 95%CL exclusion on CSxBR from ATLAS pp --> H50+X --> gam gam + X.
C Analysis and limits from 1407.6583
C INPUTS: common blocks PHYSPARAM, SM
C OUTPUT: common block CONSTR, flag ATLAS8TEVGAGAOK
      IMPLICIT NONE
C Common blocks:
      DOUBLE PRECISION MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
      COMMON/PHYSPARAMS/MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
      DOUBLE PRECISION V, MZ, MW, MTPOLE, MBMB, MCMC, MS, MTAU, MMU,
     .     ALPHAEM, ALSMZ, VCB, VUS, VUB, GAMZ, GAMW
      COMMON/SM/V,MZ,MW,MTPOLE,MBMB,MCMC,MS,MTAU,MMU,
     .     ALPHAEM,ALSMZ,VCB,VUS,VUB,GAMZ,GAMW
      DOUBLE PRECISION H5NBRGA, H5NBRZGA, H5NBRW, H5NBRZ,
     .     H5NBRZH3N, H5NBRWH3P,
     .     H5NBRH3N, H5NBRH3P,
     .     H5NWDTH
      COMMON/H5NBRS/H5NBRGA, H5NBRZGA, H5NBRW, H5NBRZ,
     .     H5NBRZH3N, H5NBRWH3P,
     .     H5NBRH3N, H5NBRH3P,
     .     H5NWDTH
      INTEGER WWJJOK, LSDMOK, H5PPOK, ATLAS8TEVGAGAOK
      INTEGER ATLAS13TEVDYGAGAOK, ATLAS13TEVVBFGAGAOK
      INTEGER DYHPPOK
      COMMON/CONSTR/WWJJOK, LSDMOK, H5PPOK, ATLAS8TEVGAGAOK, 
     .     ATLAS13TEVDYGAGAOK, ATLAS13TEVVBFGAGAOK,
     .     DYHPPOK
C Local variables:
      DOUBLE PRECISION MASSExclusion(117), SIGBR(117)
      DOUBLE PRECISION MASSCSEFFP(109), CSP(109), EFFP(109)
      DOUBLE PRECISION ERRPUP(109),ERRPDO(109)
      DOUBLE PRECISION MASSCSEFFM(109), CSM(109), EFFM(109)
      DOUBLE PRECISION ERRMUP(109),ERRMDO(109)
      DOUBLE PRECISION CSinterP,EFFinterP,SIGBRLIMIT
      DOUBLE PRECISION CSinterM,EFFinterM
      INTEGER I, NPOINTS
      INTEGER ID(109)

C First Read in the Limits from ATLAS
      NPOINTS = 117
      OPEN (UNIT = 99, FILE = 'src/CSxBRUpperAA.dat', 
     .     STATUS = 'OLD')
      DO 10 I = 1,3
         READ (99,*)
 10   CONTINUE
      DO 11 I = 1,NPOINTS
         READ (99,*) MASSExclusion(I), SIGBR(I)
 11   CONTINUE
      CLOSE (99)
C Then Read in the CS 
      NPOINTS = 109
      OPEN (UNIT = 99, FILE = 'src/Eff_H50H5p_NLO.dat',
     .     STATUS = 'OLD')
      DO 12 I = 1,3
        READ (99,*)
 12   CONTINUE
      DO 13 I = 1,NPOINTS
        READ (99,*) ID(I),MASSCSEFFP(I),CSP(I),
     .   ERRPUP(I),ERRPDO(I),EFFP(I)
 13   CONTINUE
      CLOSE (99)

      OPEN (UNIT = 99, FILE = 'src/Eff_H50H5m_NLO.dat',
     .     STATUS = 'OLD')
      DO 14 I = 1,3
        READ (99,*)
 14   CONTINUE
      DO 15 I = 1,NPOINTS
        READ (99,*) ID(I),MASSCSEFFM(I),CSM(I),
     .   ERRMUP(I),ERRMDO(I),EFFM(I)
 15   CONTINUE
      CLOSE (99)


      ATLAS8TEVGAGAOK = 1
C Check if we are within the range of the constraint:
      IF (MH5.GE.MASSExclusion(1).AND.MH5.LE.MASSExclusion(117)) THEN
         I = 1
 20      IF (MASSExclusion(I+1).GE.MH5) THEN
C Do the linear interpolation
            CALL LININTERP(MASSExclusion(I),SIGBR(I),
     .                    MASSExclusion(I+1),SIGBR(I+1),
     .           MH5,SIGBRLIMIT)
C Calculate BR(H50 -> gamma gamma).  This is the time-consuming part.
            CALL H5NDECAYS
         ELSE
            I = I+1
            GOTO 20
         ENDIF
         I = 1
 40      IF (MASSCSEFFP(I+1).GE.MH5) THEN
             CALL LININTERP(MASSCSEFFP(I),CSP(I),
     .              MASSCSEFFP(I+1),CSP(I+1),
     .       MH5,CSinterP)
             CALL LININTERP(MASSCSEFFP(I),EFFP(I),
     .              MASSCSEFFP(I+1),EFFP(I+1),
     .       MH5,EFFinterP)
             CALL LININTERP(MASSCSEFFM(I),CSM(I),
     .              MASSCSEFFM(I+1),CSM(I+1),
     .       MH5,CSinterM)
             CALL LININTERP(MASSCSEFFM(I),EFFM(I),
     .              MASSCSEFFM(I+1),EFFM(I+1),
     .       MH5,EFFinterM)
         ELSE
             I = I+1
             GOTO 40
         ENDIF
         IF(0.95D0*(CSinterP*EFFinterP+CSinterM*EFFinterM)*H5NBRGA
     .     .GT.SIGBRLIMIT) THEN
            ATLAS8TEVGAGAOK = 0
         ENDIF
      ENDIF
      RETURN
      END

C===========================================================================

      SUBROUTINE CALCATLAS13TEVDYGAGA
C Implements the 95%CL exclusion on CSxBR from ATLAS pp --> H50 + X --> gam gam + X
C Analysis and limits from 1707.04147
C INPUTS: common blocks PHYSPARAM, SM
C OUTPUT: common block CONSTR, flag ATLAS13TEVDYGAGAOK
      IMPLICIT NONE
C Common blocks:
      DOUBLE PRECISION MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
      COMMON/PHYSPARAMS/MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
      DOUBLE PRECISION V, MZ, MW, MTPOLE, MBMB, MCMC, MS, MTAU, MMU,
     .     ALPHAEM, ALSMZ, VCB, VUS, VUB, GAMZ, GAMW
      COMMON/SM/V,MZ,MW,MTPOLE,MBMB,MCMC,MS,MTAU,MMU,
     .     ALPHAEM,ALSMZ,VCB,VUS,VUB,GAMZ,GAMW
      DOUBLE PRECISION H5NBRGA, H5NBRZGA, H5NBRW, H5NBRZ,
     .     H5NBRZH3N, H5NBRWH3P,
     .     H5NBRH3N, H5NBRH3P,
     .     H5NWDTH
      COMMON/H5NBRS/H5NBRGA, H5NBRZGA, H5NBRW, H5NBRZ,
     .     H5NBRZH3N, H5NBRWH3P,
     .     H5NBRH3N, H5NBRH3P,
     .     H5NWDTH
      INTEGER WWJJOK, LSDMOK, H5PPOK, ATLAS8TEVGAGAOK
      INTEGER ATLAS13TEVDYGAGAOK, ATLAS13TEVVBFGAGAOK
      INTEGER DYHPPOK
      COMMON/CONSTR/WWJJOK, LSDMOK, H5PPOK, ATLAS8TEVGAGAOK, 
     .     ATLAS13TEVDYGAGAOK, ATLAS13TEVVBFGAGAOK,
     .     DYHPPOK
C Local variables:
      DOUBLE PRECISION MASSExclusion(135), SIGBR(135)
      DOUBLE PRECISION MASSCSEFFP(109), CSP(109), EFFP(109)
      DOUBLE PRECISION ERRPUP(109),ERRPDO(109)
      DOUBLE PRECISION MASSCSEFFM(109), CSM(109), EFFM(109)
      DOUBLE PRECISION ERRMUP(109),ERRMDO(109)
      DOUBLE PRECISION CSinterP,EFFinterP,SIGBRLIMIT
      DOUBLE PRECISION CSinterM,EFFinterM
      INTEGER I, NPOINTS
      INTEGER ID(109)

C First Read in the Limits from ATLAS
      NPOINTS = 135
      OPEN (UNIT = 99, FILE = 'src/CSxBRUpperAA13TeV.dat', 
     .     STATUS = 'OLD')
      DO 10 I = 1,3
         READ (99,*)
 10   CONTINUE
      DO 11 I = 1,NPOINTS
         READ (99,*) MASSExclusion(I), SIGBR(I)
 11   CONTINUE
      CLOSE (99)
C Then Read in the CS 
      NPOINTS = 93
      OPEN (UNIT = 99, FILE = 'src/Eff_H50H5p_NLO_13TeV.dat',
     .     STATUS = 'OLD')
      DO 12 I = 1,3
        READ (99,*)
 12   CONTINUE
      DO 13 I = 1,NPOINTS
        READ (99,*) ID(I),MASSCSEFFP(I),CSP(I),
     .   ERRPUP(I),ERRPDO(I),ERRPUP(I),ERRPDO(I),EFFP(I)
 13   CONTINUE
      CLOSE (99)

      OPEN (UNIT = 99, FILE = 'src/Eff_H50H5m_NLO_13TeV.dat',
     .     STATUS = 'OLD')
      DO 14 I = 1,3
        READ (99,*)
 14   CONTINUE
      DO 15 I = 1,NPOINTS
        READ (99,*) ID(I),MASSCSEFFM(I),CSM(I),
     .   ERRMUP(I),ERRMDO(I),ERRMUP(I),ERRMDO(I),EFFM(I)
 15   CONTINUE
      CLOSE (99)


      ATLAS13TEVDYGAGAOK = 1
C Check if we are within the range of the Efficiency calculation:
      IF (MH5.GE.MASSCSEFFP(1).AND.MH5.LE.MASSCSEFFP(93)) THEN
         I = 1
 20      IF (MASSExclusion(I+1).GE.MH5) THEN
C Do the linear interpolation
            CALL LININTERP(MASSExclusion(I),SIGBR(I),
     .                    MASSExclusion(I+1),SIGBR(I+1),
     .           MH5,SIGBRLIMIT)
C Calculate BR(H50 -> gamma gamma).  This is the time-consuming part.
            CALL H5NDECAYS
         ELSE
            I = I+1
            GOTO 20
         ENDIF
         I = 1
 40      IF (MASSCSEFFP(I+1).GE.MH5) THEN
             CALL LININTERP(MASSCSEFFP(I),CSP(I),
     .              MASSCSEFFP(I+1),CSP(I+1),
     .       MH5,CSinterP)
             CALL LININTERP(MASSCSEFFP(I),EFFP(I),
     .              MASSCSEFFP(I+1),EFFP(I+1),
     .       MH5,EFFinterP)
             CALL LININTERP(MASSCSEFFM(I),CSM(I),
     .              MASSCSEFFM(I+1),CSM(I+1),
     .       MH5,CSinterM)
             CALL LININTERP(MASSCSEFFM(I),EFFM(I),
     .              MASSCSEFFM(I+1),EFFM(I+1),
     .       MH5,EFFinterM)
         ELSE
             I = I+1
             GOTO 40
         ENDIF
         IF(0.95D0*(CSinterP*EFFinterP+CSinterM*EFFinterM)*H5NBRGA
     .     .GT.SIGBRLIMIT) THEN
            ATLAS13TEVDYGAGAOK = 0
         ENDIF
      ENDIF
      RETURN
      END

C===========================================================================

      SUBROUTINE CALCATLAS13TEVVBFGAGA
C Implements the 95%CL exclusion on CSxBR from ATLAS pp --> H50--> gam gam
C Analysis and limits from 1707.04147
C INPUTS: common blocks PHYSPARAM, SM
C OUTPUT: common block CONSTR, flag ATLAS13TEVDYGAGAOK
      IMPLICIT NONE
C Common blocks:
      DOUBLE PRECISION MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
      COMMON/PHYSPARAMS/MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
      DOUBLE PRECISION V, MZ, MW, MTPOLE, MBMB, MCMC, MS, MTAU, MMU,
     .     ALPHAEM, ALSMZ, VCB, VUS, VUB, GAMZ, GAMW
      COMMON/SM/V,MZ,MW,MTPOLE,MBMB,MCMC,MS,MTAU,MMU,
     .     ALPHAEM,ALSMZ,VCB,VUS,VUB,GAMZ,GAMW
      DOUBLE PRECISION H5NBRGA, H5NBRZGA, H5NBRW, H5NBRZ,
     .     H5NBRZH3N, H5NBRWH3P,
     .     H5NBRH3N, H5NBRH3P,
     .     H5NWDTH
      COMMON/H5NBRS/H5NBRGA, H5NBRZGA, H5NBRW, H5NBRZ,
     .     H5NBRZH3N, H5NBRWH3P,
     .     H5NBRH3N, H5NBRH3P,
     .     H5NWDTH
      INTEGER WWJJOK, LSDMOK, H5PPOK, ATLAS8TEVGAGAOK
      INTEGER ATLAS13TEVDYGAGAOK, ATLAS13TEVVBFGAGAOK
      INTEGER DYHPPOK
      COMMON/CONSTR/WWJJOK, LSDMOK, H5PPOK, ATLAS8TEVGAGAOK, 
     .     ATLAS13TEVDYGAGAOK, ATLAS13TEVVBFGAGAOK,
     .     DYHPPOK
C Local variables:
      DOUBLE PRECISION MASSExclusion(135), SIGBR(135)
      DOUBLE PRECISION MASSEFFP(109), MASSCSP(109), CSP(109), EFFP(109)
      DOUBLE PRECISION CSinterP,EFFinterP,SIGBRLIMIT
      DOUBLE PRECISION SHSQ
      INTEGER I, NPOINTS
      INTEGER ID(109)

C First Read in the Limits from ATLAS
      NPOINTS = 135
      OPEN (UNIT = 99, FILE = 'src/CSxBRUpperAA13TeV.dat', 
     .     STATUS = 'OLD')
      DO 10 I = 1,3
         READ (99,*)
 10   CONTINUE
      DO 11 I = 1,NPOINTS
         READ (99,*) MASSExclusion(I), SIGBR(I)
 11   CONTINUE
      CLOSE (99)
C Then Read in the Eff
      NPOINTS = 93
      OPEN (UNIT = 99, FILE = 'src/Eff_H50VBF_LO_13TeV.dat',
     .     STATUS = 'OLD')
      DO 12 I = 1,3
        READ (99,*)
 12   CONTINUE
      DO 13 I = 1,NPOINTS
        READ (99,*) ID(I),MASSEFFP(I),EFFP(I)
 13   CONTINUE
      CLOSE (99)

C Then Read in the CS
      NPOINTS = 49
      OPEN (UNIT = 99, FILE = 'src/VBF_H50_NNLO_13TeV.dat',
     .     STATUS = 'OLD')
      DO 14 I = 1,3
        READ (99,*)
 14   CONTINUE
      DO 15 I = 1,NPOINTS
        READ (99,*) MASSCSP(I),CSP(I)
 15   CONTINUE
      CLOSE (99)


      ATLAS13TEVVBFGAGAOK = 1
C Check if we are within the range of the constraint:
      IF (MH5.GE.MASSCSP(1).AND.MH5.LE.MASSCSP(49)) THEN
         I = 1
 20      IF (MASSExclusion(I+1).GE.MH5) THEN
C Do the linear interpolation
            CALL LININTERP(MASSExclusion(I),SIGBR(I),
     .                    MASSExclusion(I+1),SIGBR(I+1),
     .           MH5,SIGBRLIMIT)
C Calculate BR(H50 -> gamma gamma).  This is the time-consuming part.
            CALL H5NDECAYS
            SHSQ=8.D0*VCHI**2/(VPHI**2+8.D0*VCHI**2)
         ELSE
            I = I+1
            GOTO 20
         ENDIF
         I = 1
 40      IF (MASSCSP(I+1).GE.MH5) THEN
             CALL LININTERP(MASSCSP(I),CSP(I),
     .              MASSCSP(I+1),CSP(I+1),
     .       MH5,CSinterP)
         ELSE
             I = I+1
             GOTO 40
         ENDIF
         I=1
 41      IF (MASSEFFP(I+1).GE.MH5) THEN
             CALL LININTERP(MASSEFFP(I),EFFP(I),
     .              MASSEFFP(I+1),EFFP(I+1),
     .       MH5,EFFinterP)
         ELSE
             I = I+1
             GOTO 41
         ENDIF

         IF((CSinterP*EFFinterP)*H5NBRGA*SHSQ
     .     .GT.SIGBRLIMIT) THEN
            ATLAS13TEVVBFGAGAOK = 0
         ENDIF
      ENDIF
      RETURN
      END

C===========================================================================

      SUBROUTINE CALCDYHPP
C Implements the 95%CL upper bound on Drell-Yan H5++ H5-- --> 4W production
C from ATLAS Run 2 (13 TeV) 139 fb-1 results, 2101.11961.
C Upper limit is on sigma * BR assuming BR(H++ -> W+W+) = 100%.
C Since the ATLAS analysis combines 2ssl, 3l, and 4l final states, we scale
C the exclusion limit by [BR(H5++ -> W+W+)]^2.
C We also incorporate a factor of 0.95 on the cross section to account for
C the +/- 5% theory uncertainty in the NLO QCD Drell-Yan cross section calculation.
C INPUTS: common blocks PHYSPARAMS, H5PPBRS
C OUTPUTS: common block CONSTR, flag DYHPPOK
      IMPLICIT NONE
C Common blocks:
      DOUBLE PRECISION MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
      COMMON/PHYSPARAMS/MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
      DOUBLE PRECISION H5PPBRWW, H5PPBRWH3, H5PPBRH3P, H5PPWDTH
      COMMON/H5PPBRS/H5PPBRWW, H5PPBRWH3, H5PPBRH3P, H5PPWDTH
      INTEGER WWJJOK, LSDMOK, H5PPOK, ATLAS8TEVGAGAOK
      INTEGER ATLAS13TEVDYGAGAOK, ATLAS13TEVVBFGAGAOK
      INTEGER DYHPPOK
      COMMON/CONSTR/WWJJOK, LSDMOK, H5PPOK, ATLAS8TEVGAGAOK,
     .     ATLAS13TEVDYGAGAOK, ATLAS13TEVVBFGAGAOK,
     .     DYHPPOK
C Local variables:
      DOUBLE PRECISION MASS(6), EXPXS(6)
      DOUBLE PRECISION MASSTHY(41), XSTHY(41)
      DOUBLE PRECISION SIGBRLIM, XS
      INTEGER I, NPOINTS
      INTEGER ITHY, NPOINTSTHY

C Read in the experimental limits:
      NPOINTS = 6
      OPEN (UNIT = 99, FILE = 'src/atlas-HppHmm-139fb.data',
     .     STATUS = 'OLD')
      DO 9 I = 1,5
         READ (99, *)
 9    CONTINUE
      DO 10 I = 1, NPOINTS
         READ (99, *) MASS(I), EXPXS(I)
 10   CONTINUE
      CLOSE (99)

C Read in the theoretical cross section curve:
      NPOINTSTHY = 41
      OPEN (UNIT = 98,
     .     FILE = 'src/atlas-HppHmm-139fb-theory.data',
     .     STATUS = 'OLD')
      DO 11 ITHY = 1,6
         READ (98, *)
 11   CONTINUE
      DO 12 ITHY = 1, NPOINTSTHY
         READ (98, *) MASSTHY(ITHY), XSTHY(ITHY)
 12   CONTINUE
      CLOSE (98)

      DYHPPOK = 1
C Check if we are within the range of the constraint:
      IF (MH5.GE.MASS(1).AND.MH5.LE.MASS(NPOINTS)) THEN
         ITHY = 1
 19      IF (MASSTHY(ITHY+1).GE.MH5) THEN
C Do the linear interpolation for the theory cross section
            CALL LININTERP(MASSTHY(ITHY),XSTHY(ITHY),
     .           MASSTHY(ITHY+1),XSTHY(ITHY+1),
     .           MH5,XS)
C Apply the factor of 0.95 to account for theory uncert.
            XS = 0.95D0*XS
         ELSE
            ITHY = ITHY+1
            GOTO 19
         ENDIF
         I = 1
 20      IF (MASS(I+1).GE.MH5) THEN
C Do the linear interpolation for the expt limit
            CALL LININTERP(MASS(I),EXPXS(I),MASS(I+1),EXPXS(I+1),
     .           MH5,SIGBRLIM)
         ELSE
            I = I+1
            GOTO 20
         ENDIF
         IF (XS.GT.SIGBRLIM) THEN
C We are in the danger zone, though BR(H5++ -> WW) < 1 might save us.
C Check whether H5 -> H3 V decays are open; if so, calculate BR(H5++ -> WW).
            IF (MH3.GE.MH5) THEN
               DYHPPOK = 0
            ELSE
               CALL H5PPDECAYS
               IF (XS*H5PPBRWW**2.GT.SIGBRLIM) THEN
                  DYHPPOK = 0
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      RETURN
      END

C=================================================================
