C==================================================================
C     Subroutines to read in & process the input parameters
C==================================================================

      SUBROUTINE LOAD_INPUTS
C Handles the input "case" and calls its own subroutines to read in
C the inputs.  Those subroutines do the necessary calculations
C to compute the entries in LPARAMS.
      IMPLICIT NONE
C Common blocks:
      DOUBLE PRECISION MU2SQ,MU3SQ,LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,
     .     LAMBDA5,M1,M2
      COMMON/LPARAMS/MU2SQ,MU3SQ,LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,
     .     LAMBDA5,M1,M2
      INTEGER INPUTMODE, SILENT
      COMMON/CONTROLS/INPUTMODE,SILENT
      DOUBLE PRECISION MH
      INTEGER INPUTSET
      COMMON/INPUT/MH,INPUTSET
      DOUBLE PRECISION IMHL, IMHH, IMH3, IMH5, ISH, ISA
      COMMON/INPUT3/IMHL,IMHH,IMH3,IMH5,ISH,ISA
C Local variables:

      IF (INPUTMODE.EQ.1) THEN
C Read inputs from the user (interactive mode)
c         PRINT *
c         CALL PRINT_BANNER
c         PRINT *
         PRINT *, " Available input sets: "
         PRINT *, " 1)  mu3^2, lambda1, lambda2, lambda3, lambda4, ",
     .        "lambda5, M1, M2"
         PRINT *, " 2)  mu3^2, mh, lambda2, lambda3, lambda4, ",
     .        "lambda5, M1, M2"
         PRINT *, " 3)  mh, mH, m3, m5, sin(thetaH), sin(alpha), ",
     .        "M1, M2"
         PRINT *, " 4)  mh, m5, sin(thetaH), lambda2, lambda3, ",
     .        "lambda4, M1, M2"
         PRINT *, " 5)  mh, mH, sin(thetaH), sin(alpha), ",
     .        "lambda2, lambda3, lambda4, lambda5"
         PRINT *, " 6)  mh, m5, sin(thetaH), lambda2, lambda3, ",
     .        "lambda4, lambda5, M2"
         PRINT *, " All masses in GeV and masses-squared in GeV^2. "
         PRINT *, " Total v ~= 246 GeV from G_Fermi is hard-coded. "
         PRINT *
         PRINT *, " Enter your input set choice: "
         READ *, INPUTSET
      ELSE IF (INPUTMODE.EQ.0.AND.SILENT.EQ.0) THEN
C Take inputs from the main program
c         PRINT *
c         CALL PRINT_BANNER
c         PRINT *
         PRINT *, " Taking inputs from the main program."
         PRINT *, " INPUTSET = ", INPUTSET
         PRINT *, " Total vev v is computed from G_Fermi. "
         IF (INPUTSET.EQ.1) THEN
            PRINT *, " mu3^2 = ", MU3SQ, " GeV^2"
            PRINT *, " lambda1 = ", LAMBDA1
            PRINT *, " lambda2 = ", LAMBDA2
            PRINT *, " lambda3 = ", LAMBDA3
            PRINT *, " lambda4 = ", LAMBDA4
            PRINT *, " lambda5 = ", LAMBDA5
            PRINT *, " M1 = ", M1, " GeV"
            PRINT *, " M2 =      ", M2, " GeV"
         ELSE IF (INPUTSET.EQ.2) THEN
            PRINT *, " mu3^2 =   ", MU3SQ, " GeV^2 "
            PRINT *, " mh =      ", MH, " GeV"
            PRINT *, " lambda2 = ", LAMBDA2
            PRINT *, " lambda3 = ", LAMBDA3
            PRINT *, " lambda4 = ", LAMBDA4
            PRINT *, " lambda5 = ", LAMBDA5
            PRINT *, " M1 =      ", M1, " GeV"
            PRINT *, " M2 =      ", M2, " GeV"
         ELSE IF (INPUTSET.EQ.3) THEN
            PRINT *, " mh =          ", IMHL, " GeV"
            PRINT *, " mH =          ", IMHH, " GeV"
            PRINT *, " m3 =          ", IMH3, " GeV"
            PRINT *, " m5 =          ", IMH5, " GeV"
            PRINT *, " sin(thetaH) = ", ISH
            PRINT *, " sin(alpha) =  ", ISA
            PRINT *, " M1 =          ", M1, " GeV"
            PRINT *, " M2 =          ", M2, " GeV"
         ELSE IF (INPUTSET.EQ.4) THEN
            PRINT *, " mh =          ", IMHL, " GeV"
            PRINT *, " m5 =          ", IMH5, " GeV"
            PRINT *, " sin(thetaH) = ", ISH
            PRINT *, " lambda2 =     ", LAMBDA2
            PRINT *, " lambda3 =     ", LAMBDA3
            PRINT *, " lambda4 =     ", LAMBDA4
            PRINT *, " M1 =          ", M1, " GeV"
            PRINT *, " M2 =          ", M2, " GeV"
         ELSE IF (INPUTSET.EQ.5) THEN
            PRINT *, " mh =          ", IMHL, " GeV"
            PRINT *, " mH =          ", IMHH, " GeV"
            PRINT *, " sin(thetaH) = ", ISH
            PRINT *, " sin(alpha) =  ", ISA
            PRINT *, " lambda2 =     ", LAMBDA2
            PRINT *, " lambda3 =     ", LAMBDA3
            PRINT *, " lambda4 =     ", LAMBDA4
            PRINT *, " lambda5 =     ", LAMBDA5
         ELSE IF (INPUTSET.EQ.6) THEN
            PRINT *, " mh =          ", IMHL, " GeV"
            PRINT *, " m5 =          ", IMH5, " GeV"
            PRINT *, " sin(thetaH) = ", ISH
            PRINT *, " lambda2 =     ", LAMBDA2
            PRINT *, " lambda3 =     ", LAMBDA3
            PRINT *, " lambda4 =     ", LAMBDA4
            PRINT *, " lambda5 =     ", LAMBDA5
            PRINT *, " M2 =          ", M2, " GeV"
         ENDIF
      ENDIF

      IF (INPUTSET.EQ.1) THEN
         CALL READIN1
      ELSE IF (INPUTSET.EQ.2) THEN
         CALL READIN2
      ELSE IF (INPUTSET.EQ.3) THEN
         CALL READIN3
      ELSE IF (INPUTSET.EQ.4) THEN
         CALL READIN4
      ELSE IF (INPUTSET.EQ.5) THEN
         CALL READIN5
      ELSE IF (INPUTSET.EQ.6) THEN
         CALL READIN6
      ELSE
         PRINT *, " Invalid INPUTSET = ", INPUTSET
      ENDIF

      RETURN
      END

C==================================================================
      SUBROUTINE READIN1
C Reads in the inputs for Input Set 1 and computes the rest of the
C Lagrangian parameters in LPARAMS.
C Input Set 1 is:
C      mu3^2, lambda1, lambda2, lambda3, lambda4, lambda5, M1, M2
C      with v = (sqrt(2) G_F)^{-1/2} ~ 246 GeV set by INITIALIZE_SM.
C
      IMPLICIT NONE
C Common blocks:
      DOUBLE PRECISION MU2SQ,MU3SQ,LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,
     .     LAMBDA5,M1,M2
      COMMON/LPARAMS/MU2SQ,MU3SQ,LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,
     .     LAMBDA5,M1,M2
      INTEGER UNIOK, BFBOK, POSMSQOK, MINOK, INPUTOK
      COMMON/FLAGS/UNIOK, BFBOK, POSMSQOK, MINOK, INPUTOK
      DOUBLE PRECISION V, MZ, MW, MTPOLE, MBMB, MCMC, MS, MTAU, MMU,
     .     ALPHAEM,ALSMZ, VCB, VUS, VUB, GAMZ, GAMW
      COMMON/SM/V,MZ,MW,MTPOLE,MBMB,MCMC,MS,MTAU,MMU,
     .     ALPHAEM,ALSMZ,VCB,VUS,VUB,GAMZ,GAMW
      INTEGER INPUTMODE, SILENT
      COMMON/CONTROLS/INPUTMODE,SILENT
      DOUBLE PRECISION MH
      INTEGER INPUTSET
      COMMON/INPUT/MH,INPUTSET
      DOUBLE PRECISION ZETA, OMEGA, SIGMA, RHO
      COMMON/HYPER/ZETA,OMEGA,SIGMA,RHO
      DOUBLE PRECISION AA, BB, VV
      DIMENSION AA(6), BB(6), VV(6)
      COMMON/MINIMA/AA,BB,VV
C Local variables:
      DOUBLE PRECISION VPHI, VCHI, MU2SQTEST
      DIMENSION VPHI(3), VCHI(3), MU2SQTEST(3)
      INTEGER SOLOK
      DIMENSION SOLOK(3)
      DOUBLE PRECISION GOODVPHI, GOODVCHI, GOODV
      DOUBLE PRECISION ALPHA, BETA, GAMMA, DELTA, RDISC
      DOUBLE COMPLEX DISC, DISC0, DISC1
      DOUBLE COMPLEX II, U1, U2, U3, CROOT, CVCHI1, CVCHI2, CVCHI3
      INTEGER J

      IF (INPUTMODE.EQ.1) THEN
C Read inputs from the user (interactive mode)
      PRINT *
      PRINT *, " Enter the 8 parameters on one line (no commas): "
      READ *, MU3SQ, LAMBDA1, LAMBDA2, LAMBDA3, LAMBDA4, LAMBDA5, M1, M2
      ENDIF

      ZETA = 1.D0/3.D0
      OMEGA = 1.D0/2.D0
      SIGMA = DSQRT(3.D0)/4.D0
      RHO = 2.D0/DSQRT(3.D0)
C Solve the cubic equation for v_chi.
      II = (0.D0,1.D0)
      U1 = (1.D0,0.D0)
      U2 = -0.5D0 + II*DSQRT(3.D0)/2.D0
      U3 = -0.5D0 - II*DSQRT(3.D0)/2.D0
      ALPHA = -48.D0*LAMBDA2 + 24.D0*LAMBDA5 + 12.D0*LAMBDA3 
     .     + 36.D0*LAMBDA4
      BETA = 6.D0*M1 - 18.D0*M2
      GAMMA = 3.D0*MU3SQ + (6.D0*LAMBDA2 - 3.D0*LAMBDA5)*V**2
      DELTA = -3.D0/4.D0*M1*V**2
      DISC = 18.D0 * ALPHA*BETA*GAMMA*DELTA - 4.D0 * BETA**3 *DELTA
     .     + BETA**2 * GAMMA**2 - 4.D0 * ALPHA * GAMMA**3
     .     - 27.D0 * ALPHA**2 * DELTA**2
      DISC0 = BETA**2 - 3.D0 * ALPHA * GAMMA
      DISC1 = 2.D0 * BETA**3 - 9.D0 * ALPHA*BETA*GAMMA 
     .     + 27.D0 * ALPHA**2 * DELTA
      CROOT = (0.5D0 * DISC1 + 0.5D0*CDSQRT(-27.D0*ALPHA**2*DISC)
     .     )**(1.D0/3.D0)
      CVCHI1 = -1.D0/3.D0/ALPHA * (BETA + U1*CROOT + DISC0/U1/CROOT)
      CVCHI2 = -1.D0/3.D0/ALPHA * (BETA + U2*CROOT + DISC0/U2/CROOT)
      CVCHI3 = -1.D0/3.D0/ALPHA * (BETA + U3*CROOT + DISC0/U3/CROOT)
      RDISC = DREAL(DISC)
      IF (RDISC.GE.0.D0) THEN
C The 3 roots are all real (if RDISC = 0 then two are the same)
         VCHI(1) = DREAL(CVCHI1)
         VCHI(2) = DREAL(CVCHI2)
         VCHI(3) = DREAL(CVCHI3)
         DO 10 J=1,3
C Check that v_chi is not too big for each solution.
            IF (8.D0*VCHI(J)**2.LT.V**2) THEN
               SOLOK(J) = 1
               VPHI(J) = DSQRT(V**2 - 8.D0*VCHI(J)**2)
               MU2SQTEST(J) = 2.D0*(-2.D0*LAMBDA1*VPHI(J)**2 
     .              - (3.D0*LAMBDA2 - 3.D0/2.D0*LAMBDA5)*VCHI(J)**2
     .              + 3.D0/4.D0*M1*VCHI(J))
            ELSE
               SOLOK(J) = 0
               MU2SQTEST(J) = 0.D0
            ENDIF
 10      CONTINUE
C Choose among the 3 solutions: check whether each v_chi is actually the 
C minimum for this set of Lagrangian parameters, using GETMINIMA.
         MU2SQ = 0.D0
         INPUTOK = 0
         DO 20 J=1,3
            IF (SOLOK(J).EQ.1) THEN
               MU2SQ = MU2SQTEST(J)
               CALL GETMINIMA
C Identify the deepest of the good extrema:
               GOODVPHI = AA(1)
               GOODVCHI = BB(1)/DSQRT(3.D0)
               GOODV = VV(1)
               IF (VV(2).LT.GOODV) THEN
                  GOODVPHI = AA(2)
                  GOODVCHI = BB(2)/DSQRT(3.D0)
                  GOODV = VV(2)
               ELSE IF (VV(3).LT.GOODV) THEN
                  GOODVPHI = AA(3)
                  GOODVCHI = BB(3)/DSQRT(3.D0)
                  GOODV = VV(3)
               ENDIF
               IF (DABS(VCHI(J)-GOODVCHI).LT.1D-6) THEN
                  INPUTOK = 1
               ELSE
                  SOLOK(J) = 0
               ENDIF
            ENDIF
 20      CONTINUE
         IF ((SOLOK(1) + SOLOK(2) + SOLOK(3)).GT.1.AND.SILENT.EQ.0) THEN
            PRINT *, "*** READIN1: there is more than one",
     .           " valid solution for v_chi! ***"
            INPUTOK = 0
         ELSE
C Now set the value of mu_2^2.
            DO 30 J=1,3
               IF (SOLOK(J).EQ.1) THEN
                  MU2SQ = MU2SQTEST(J)
               ENDIF
 30         CONTINUE
         ENDIF
      ELSE
C Only one root is real.  First determine which one it is.
         IF (DABS(DIMAG(CVCHI1)).LT.1D-6*CDABS(CVCHI1)) THEN
            VCHI(1) = DREAL(CVCHI1)
         ELSE IF (DABS(DIMAG(CVCHI2)).LT.1D-6*CDABS(CVCHI2)) THEN
            VCHI(1) = DREAL(CVCHI2)
         ELSE
            VCHI(1) = DREAL(CVCHI3)
         ENDIF
C Check that the one solution v_chi is not too big.
         IF (8.D0*VCHI(1)**2.LT.V**2) THEN
C Now compute mu2^2.
            VPHI(1) = DSQRT(V**2 - 8.D0*VCHI(1)**2)
            MU2SQ = 2.D0*(-2.D0*LAMBDA1*VPHI(1)**2 
     .           - (3.D0*LAMBDA2 - 3.D0/2.D0*LAMBDA5)*VCHI(1)**2
     .           + 3.D0/4.D0*M1*VCHI(1))
C Check whether v_chi is actually the minimum for this set of Lagrangian 
C parameters, using GETMINIMA.
            CALL GETMINIMA
C Identify the deepest of the good extrema:
            GOODVPHI = AA(1)
            GOODVCHI = BB(1)/DSQRT(3.D0)
            GOODV = VV(1)
            IF (VV(2).LT.GOODV) THEN
               GOODVPHI = AA(2)
               GOODVCHI = BB(2)/DSQRT(3.D0)
               GOODV = VV(2)
            ELSE IF (VV(3).LT.GOODV) THEN
               GOODVPHI = AA(3)
               GOODVCHI = BB(3)/DSQRT(3.D0)
               GOODV = VV(3)
            ENDIF
            IF (DABS(VCHI(1)-GOODVCHI).LT.1D-6) THEN
               INPUTOK = 1
            ELSE
               INPUTOK = 0
            ENDIF
         ELSE
            INPUTOK = 0
            MU2SQ = 0.D0
         ENDIF
      ENDIF

      RETURN
      END

C==================================================================
      SUBROUTINE READIN2
C Reads in the inputs for Input Set 2 and computes the rest of the
C Lagrangian parameters in LPARAMS.
C Input Set 2 is:
C      mu3^2, mh, lambda2, lambda3, lambda4, lambda5, M1, M2
C      with v = 246 GeV hard-coded (it is set by INITIALIZE_SM).
C
      IMPLICIT NONE
C Common blocks:
      DOUBLE PRECISION MU2SQ,MU3SQ,LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,
     .     LAMBDA5,M1,M2
      COMMON/LPARAMS/MU2SQ,MU3SQ,LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,
     .     LAMBDA5,M1,M2
      INTEGER UNIOK, BFBOK, POSMSQOK, MINOK, INPUTOK
      COMMON/FLAGS/UNIOK, BFBOK, POSMSQOK, MINOK, INPUTOK
      DOUBLE PRECISION V, MZ, MW, MTPOLE, MBMB, MCMC, MS, MTAU, MMU,
     .     ALPHAEM,ALSMZ, VCB, VUS, VUB, GAMZ, GAMW
      COMMON/SM/V,MZ,MW,MTPOLE,MBMB,MCMC,MS,MTAU,MMU,
     .     ALPHAEM,ALSMZ,VCB,VUS,VUB,GAMZ,GAMW
      INTEGER INPUTMODE, SILENT
      COMMON/CONTROLS/INPUTMODE,SILENT
      DOUBLE PRECISION MH
      INTEGER INPUTSET
      COMMON/INPUT/MH,INPUTSET
      DOUBLE PRECISION ZETA, OMEGA, SIGMA, RHO
      COMMON/HYPER/ZETA,OMEGA,SIGMA,RHO
      DOUBLE PRECISION AA, BB, VV
      DIMENSION AA(6), BB(6), VV(6)
      COMMON/MINIMA/AA,BB,VV
C Local variables:
      DOUBLE PRECISION VPHI, VCHI, MU2SQTEST, L1TEST
      DIMENSION VPHI(3), VCHI(3), MU2SQTEST(3), L1TEST(3)
      DOUBLE PRECISION MSQ12, MSQ22
      INTEGER SOLOK
      DIMENSION SOLOK(3)
      DOUBLE PRECISION GOODVPHI, GOODVCHI, GOODV
      DOUBLE PRECISION ALPHA, BETA, GAMMA, DELTA, RDISC
      DOUBLE COMPLEX DISC, DISC0, DISC1
      DOUBLE COMPLEX II, U1, U2, U3, CROOT, CVCHI1, CVCHI2, CVCHI3
      INTEGER J
      
      IF (INPUTMODE.EQ.1) THEN
C Read inputs from the user (interactive mode)
      PRINT *
      PRINT *, " Enter the 8 parameters on one line (no commas): "
      READ *, MU3SQ, MH, LAMBDA2, LAMBDA3, LAMBDA4, LAMBDA5, M1, M2
      ENDIF

      ZETA = 1.D0/3.D0
      OMEGA = 1.D0/2.D0
      SIGMA = DSQRT(3.D0)/4.D0
      RHO = 2.D0/DSQRT(3.D0)
C Solve the cubic equation for v_chi.
      II = (0.D0,1.D0)
      U1 = (1.D0,0.D0)
      U2 = -0.5D0 + II*DSQRT(3.D0)/2.D0
      U3 = -0.5D0 - II*DSQRT(3.D0)/2.D0
      ALPHA = -48.D0*LAMBDA2 + 24.D0*LAMBDA5 + 12.D0*LAMBDA3 
     .     + 36.D0*LAMBDA4
      BETA = 6.D0*M1 - 18.D0*M2
      GAMMA = 3.D0*MU3SQ + (6.D0*LAMBDA2 - 3.D0*LAMBDA5)*V**2
      DELTA = -3.D0/4.D0*M1*V**2
      DISC = 18.D0 * ALPHA*BETA*GAMMA*DELTA - 4.D0 * BETA**3 *DELTA
     .     + BETA**2 * GAMMA**2 - 4.D0 * ALPHA * GAMMA**3
     .     - 27.D0 * ALPHA**2 * DELTA**2
      DISC0 = BETA**2 - 3.D0 * ALPHA * GAMMA
      DISC1 = 2.D0 * BETA**3 - 9.D0 * ALPHA*BETA*GAMMA 
     .     + 27.D0 * ALPHA**2 * DELTA
      CROOT = (0.5D0 * DISC1 + 0.5D0*CDSQRT(-27.D0*ALPHA**2*DISC)
     .     )**(1.D0/3.D0)
      CVCHI1 = -1.D0/3.D0/ALPHA * (BETA + U1*CROOT + DISC0/U1/CROOT)
      CVCHI2 = -1.D0/3.D0/ALPHA * (BETA + U2*CROOT + DISC0/U2/CROOT)
      CVCHI3 = -1.D0/3.D0/ALPHA * (BETA + U3*CROOT + DISC0/U3/CROOT)
      RDISC = DREAL(DISC)
      IF (RDISC.GE.0.D0) THEN
C The 3 roots are all real (if RDISC = 0 then two are the same)
         VCHI(1) = DREAL(CVCHI1)
         VCHI(2) = DREAL(CVCHI2)
         VCHI(3) = DREAL(CVCHI3)
         DO 10 J=1,3
C Check that v_chi is not too big for each solution.
            IF (8.D0*VCHI(J)**2.LT.V**2) THEN
               SOLOK(J) = 1
               VPHI(J) = DSQRT(V**2 - 8.D0*VCHI(J)**2)
C Compute lambda1 and mu2sq for each solution:
               MSQ12 = DSQRT(3.D0)*VPHI(J) * (-M1/2.D0 
     .              + (4.D0*LAMBDA2 - 2.D0*LAMBDA5)*VCHI(J))
C               MSQ22 = M1*VPHI(J)**2/4.D0/VCHI(J) - 6.D0*M2*VCHI(J)
C     .              + (8.D0*LAMBDA3 + 24.D0*LAMBDA4)*VCHI(J)**2
               MSQ22 = MU3SQ - 12.D0*M2*VCHI(J) 
     .          + (2.D0*LAMBDA2-LAMBDA5)*VPHI(J)**2
     .          + 12.D0*(LAMBDA3+3.D0*LAMBDA4)*VCHI(J)**2
               L1TEST(J) = 1.D0/8.D0/VPHI(J)**2 
     .              * (MH**2 + MSQ12**2/(MSQ22 - MH**2))
               MU2SQTEST(J) = 2.D0*(-2.D0*L1TEST(J)*VPHI(J)**2 
     .              - (3.D0*LAMBDA2 - 3.D0/2.D0*LAMBDA5)*VCHI(J)**2
     .              + 3.D0/4.D0*M1*VCHI(J))
            ELSE
               SOLOK(J) = 0
               L1TEST(J) = 0.D0
               MU2SQTEST(J) = 0.D0
            ENDIF
 10      CONTINUE
C Choose among the 3 solutions: check whether each v_chi is actually the 
C minimum for this set of Lagrangian parameters, using GETMINIMA.
         LAMBDA1 = 0.D0
         MU2SQ = 0.D0
         INPUTOK = 0
         DO 20 J=1,3
            IF (SOLOK(J).EQ.1) THEN
               LAMBDA1 = L1TEST(J)
               MU2SQ = MU2SQTEST(J)
               CALL GETMINIMA
C Identify the deepest of the good extrema:
               GOODVPHI = AA(1)
               GOODVCHI = BB(1)/DSQRT(3.D0)
               GOODV = VV(1)
               IF (VV(2).LT.GOODV) THEN
                  GOODVPHI = AA(2)
                  GOODVCHI = BB(2)/DSQRT(3.D0)
                  GOODV = VV(2)
               ELSE IF (VV(3).LT.GOODV) THEN
                  GOODVPHI = AA(3)
                  GOODVCHI = BB(3)/DSQRT(3.D0)
                  GOODV = VV(3)
               ENDIF
               IF (DABS(VCHI(J)-GOODVCHI).LT.1D-6) THEN
                  INPUTOK = 1
               ELSE
                  SOLOK(J) = 0
               ENDIF
            ENDIF
 20      CONTINUE
         IF ((SOLOK(1) + SOLOK(2) + SOLOK(3)).GT.1.AND.SILENT.EQ.0) THEN
            PRINT *, "*** READIN2: there is more than one",
     .           " valid solution for v_chi! ***"
            INPUTOK = 0
         ELSE
C Now set the values of lambda1 and mu_2^2.
            DO 30 J=1,3
               IF (SOLOK(J).EQ.1) THEN
                  LAMBDA1 = L1TEST(J)
                  MU2SQ = MU2SQTEST(J)
               ENDIF
 30         CONTINUE
         ENDIF
      ELSE
C Only one root is real.  First determine which one it is.
         IF (DABS(DIMAG(CVCHI1)).LT.1D-6*CDABS(CVCHI1)) THEN
            VCHI(1) = DREAL(CVCHI1)
         ELSE IF (DABS(DIMAG(CVCHI2)).LT.1D-6*CDABS(CVCHI2)) THEN
            VCHI(1) = DREAL(CVCHI2)
         ELSE
            VCHI(1) = DREAL(CVCHI3)
         ENDIF
C Check that the one solution v_chi is not too big.
         IF (8.D0*VCHI(1)**2.LT.V**2) THEN
C Now compute lambda1 and mu2^2.
            VPHI(1) = DSQRT(V**2 - 8.D0*VCHI(1)**2)
            MSQ12 = DSQRT(3.D0)*VPHI(1) * (-M1/2.D0 
     .           + (4.D0*LAMBDA2 - 2.D0*LAMBDA5)*VCHI(1))
C             MSQ22 = M1*VPHI(1)**2/4.D0/VCHI(1) - 6.D0*M2*VCHI(1)
C      .           + (8.D0*LAMBDA3 + 24.D0*LAMBDA4)*VCHI(1)**2
            MSQ22 = MU3SQ - 12.D0*M2*VCHI(1) 
     .          + (2.D0*LAMBDA2-LAMBDA5)*VPHI(1)**2
     .          + 12.D0*(LAMBDA3+3.D0*LAMBDA4)*VCHI(1)**2
            LAMBDA1 = 1.D0/8.D0/VPHI(1)**2 
     .              * (MH**2 + MSQ12**2/(MSQ22 - MH**2))
            MU2SQ = 2.D0*(-2.D0*LAMBDA1*VPHI(1)**2 
     .           - (3.D0*LAMBDA2 - 3.D0/2.D0*LAMBDA5)*VCHI(1)**2
     .           + 3.D0/4.D0*M1*VCHI(1))
C Check whether v_chi is actually the minimum for this set of Lagrangian 
C parameters, using GETMINIMA.
            CALL GETMINIMA
C Identify the deepest of the good extrema:
            GOODVPHI = AA(1)
            GOODVCHI = BB(1)/DSQRT(3.D0)
            GOODV = VV(1)
            IF (VV(2).LT.GOODV) THEN
               GOODVPHI = AA(2)
               GOODVCHI = BB(2)/DSQRT(3.D0)
               GOODV = VV(2)
            ELSE IF (VV(3).LT.GOODV) THEN
               GOODVPHI = AA(3)
               GOODVCHI = BB(3)/DSQRT(3.D0)
               GOODV = VV(3)
            ENDIF
            IF (DABS(VCHI(1)-GOODVCHI).LT.1D-6) THEN
               INPUTOK = 1
            ELSE
               INPUTOK = 0
            ENDIF
         ELSE
            INPUTOK = 0
            LAMBDA1 = 0.D0
            MU2SQ = 0.D0
         ENDIF
      ENDIF

      RETURN
      END


C==================================================================
      SUBROUTINE READIN3
C Reads in the inputs for Input Set 3 and computes the rest of the 
C Lagrangian parameters in LPARAMS.
C Input Set 3 is:
C     mh, mH, m3, m5, sin(thetaH), sin(alpha), M1, M2
C     with v = (sqrt(2) G_F)^{-1/2} ~ 246 GeV set by INITIALIZE_SM.
C
      IMPLICIT NONE
C Common blocks:
      DOUBLE PRECISION MU2SQ,MU3SQ,LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,
     .     LAMBDA5,M1,M2
      COMMON/LPARAMS/MU2SQ,MU3SQ,LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,
     .     LAMBDA5,M1,M2
      INTEGER UNIOK, BFBOK, POSMSQOK, MINOK, INPUTOK
      COMMON/FLAGS/UNIOK, BFBOK, POSMSQOK, MINOK, INPUTOK
      DOUBLE PRECISION V, MZ, MW, MTPOLE, MBMB, MCMC, MS, MTAU, MMU,
     .     ALPHAEM,ALSMZ, VCB, VUS, VUB, GAMZ, GAMW
      COMMON/SM/V,MZ,MW,MTPOLE,MBMB,MCMC,MS,MTAU,MMU,
     .     ALPHAEM,ALSMZ,VCB,VUS,VUB,GAMZ,GAMW
      INTEGER INPUTMODE, SILENT
      COMMON/CONTROLS/INPUTMODE,SILENT
      DOUBLE PRECISION MH
      INTEGER INPUTSET
      COMMON/INPUT/MH,INPUTSET
      DOUBLE PRECISION IMHL, IMHH, IMH3, IMH5, ISH, ISA
      COMMON/INPUT3/IMHL,IMHH,IMH3,IMH5,ISH,ISA
      DOUBLE PRECISION ZETA, OMEGA, SIGMA, RHO
      COMMON/HYPER/ZETA,OMEGA,SIGMA,RHO
      DOUBLE PRECISION AA, BB, VV
      DIMENSION AA(6), BB(6), VV(6)
      COMMON/MINIMA/AA,BB,VV
C Local variables:
      DOUBLE PRECISION VCHI, VPHI, CA
      DOUBLE PRECISION M11SQ, M12SQ, M22SQ
      DOUBLE PRECISION GOODVPHI, GOODVCHI, GOODV
      DOUBLE PRECISION MHL, MHH, MH3, MH5, SH, SA

      IF (INPUTMODE.EQ.1) THEN
C Read inputs from the user (interactive mode)
         PRINT *
         PRINT *, " Enter the 8 parameters on one line (no commas): "
         READ *, IMHL, IMHH, IMH3, IMH5, ISH, ISA, M1, M2
      ENDIF

      MHL = IMHL
      MHH = IMHH
      MH3 = IMH3
      MH5 = IMH5
      SH = ISH
      SA = ISA
      
      INPUTOK = 1

      IF (SH.GE.1.D0.OR.SH.LE.0.D0) THEN
         INPUTOK = 0
      ENDIF
      IF (SA.GT.1.D0.OR.SA.LT.-1.D0) THEN
         INPUTOK = 0
      ENDIF

      VCHI = V * SH / DSQRT(8.D0)
      VPHI = DSQRT(V**2 - 8.D0*VCHI**2)
      LAMBDA5 = 2.D0*MH3**2/V**2 - M1/2.D0/VCHI
      LAMBDA3 = (MH5**2 - M1*VPHI**2/4.D0/VCHI - 12.D0*M2*VCHI
     .     - 3.D0/2.D0*LAMBDA5*VPHI**2) / 8.D0 / VCHI**2
      CA = DSQRT(1.D0 - SA**2)
      M11SQ = CA**2*MHL**2 + SA**2*MHH**2
      M12SQ = SA*CA*(MHH**2 - MHL**2)
      M22SQ = SA**2*MHL**2 + CA**2*MHH**2
      LAMBDA1 = M11SQ/8.D0/VPHI**2
      LAMBDA2 = (2.D0/DSQRT(3.D0)/VPHI*M12SQ + M1 
     .     + 4.D0*LAMBDA5*VCHI) / 8.D0 / VCHI
      LAMBDA4 = (M22SQ - M1*VPHI**2/4.D0/VCHI + 6.D0*M2*VCHI
     .     - 8.D0*LAMBDA3*VCHI**2) / 24.D0 / VCHI**2
      MU3SQ = M1*VPHI**2/4.D0/VCHI + 6.D0*M2*VCHI 
     .     - (2.D0*LAMBDA2 - LAMBDA5)*VPHI**2 
     .     - 4.D0*(LAMBDA3 + 3.D0*LAMBDA4)*VCHI**2
      MU2SQ = -4.D0*LAMBDA1*VPHI**2 
     .     - 3.D0*(2.D0*LAMBDA2 - LAMBDA5)*VCHI**2
     .     + 3.D0/2.D0*M1*VCHI

C Check that the inputted parameters give the right minimum
      ZETA = 1.D0/3.D0
      OMEGA = 1.D0/2.D0
      SIGMA = DSQRT(3.D0)/4.D0
      RHO = 2.D0/DSQRT(3.D0)
      CALL GETMINIMA
C Identify the deepest of the good extrema:
      GOODVPHI = AA(1)
      GOODVCHI = BB(1)/DSQRT(3.D0)
      GOODV = VV(1)
      IF (VV(2).LT.GOODV) THEN
         GOODVPHI = AA(2)
         GOODVCHI = BB(2)/DSQRT(3.D0)
         GOODV = VV(2)
      ELSE IF (VV(3).LT.GOODV) THEN
         GOODVPHI = AA(3)
         GOODVCHI = BB(3)/DSQRT(3.D0)
         GOODV = VV(3)
      ENDIF
      IF (DABS(VCHI-GOODVCHI).GT.1D-6) THEN
         INPUTOK = 0
      ENDIF

      RETURN
      END


C==================================================================
      SUBROUTINE READIN4
C Reads in the inputs for Input Set 4 and computes the rest of the 
C Lagrangian parameters in LPARAMS.
C Input Set 4 is:
C     mh, m5, sin(thetaH), lambda2, lambda3, lambda4, M1, M2
C     with v = (sqrt(2) G_F)^{-1/2} ~ 246 GeV set by INITIALIZE_SM.
C
      IMPLICIT NONE
C Common blocks:
      DOUBLE PRECISION MU2SQ,MU3SQ,LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,
     .     LAMBDA5,M1,M2
      COMMON/LPARAMS/MU2SQ,MU3SQ,LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,
     .     LAMBDA5,M1,M2
      INTEGER UNIOK, BFBOK, POSMSQOK, MINOK, INPUTOK
      COMMON/FLAGS/UNIOK, BFBOK, POSMSQOK, MINOK, INPUTOK
      DOUBLE PRECISION V, MZ, MW, MTPOLE, MBMB, MCMC, MS, MTAU, MMU,
     .     ALPHAEM,ALSMZ, VCB, VUS, VUB, GAMZ, GAMW
      COMMON/SM/V,MZ,MW,MTPOLE,MBMB,MCMC,MS,MTAU,MMU,
     .     ALPHAEM,ALSMZ,VCB,VUS,VUB,GAMZ,GAMW
      INTEGER INPUTMODE, SILENT
      COMMON/CONTROLS/INPUTMODE,SILENT
      DOUBLE PRECISION MH
      INTEGER INPUTSET
      COMMON/INPUT/MH,INPUTSET
      DOUBLE PRECISION IMHL, IMHH, IMH3, IMH5, ISH, ISA
      COMMON/INPUT3/IMHL,IMHH,IMH3,IMH5,ISH,ISA
      DOUBLE PRECISION ZETA, OMEGA, SIGMA, RHO
      COMMON/HYPER/ZETA,OMEGA,SIGMA,RHO
      DOUBLE PRECISION AA, BB, VV
      DIMENSION AA(6), BB(6), VV(6)
      COMMON/MINIMA/AA,BB,VV
C Local variables:
      DOUBLE PRECISION VCHI, VPHI
      DOUBLE PRECISION M12SQ, M22SQ
      DOUBLE PRECISION GOODVPHI, GOODVCHI, GOODV
      DOUBLE PRECISION MHL, MH5, SH

      IF (INPUTMODE.EQ.1) THEN
C Read inputs from the user (interactive mode)
         PRINT *
         PRINT *, " Enter the 8 parameters on one line (no commas): "
         READ *, IMHL, IMH5, ISH, LAMBDA2, LAMBDA3, LAMBDA4, M1, M2
      ENDIF
      
      MHL = IMHL
      MH5 = IMH5
      SH = ISH

      INPUTOK = 1

      IF (SH.GE.1.D0.OR.SH.LE.0.D0) THEN
         INPUTOK = 0
      ENDIF

      VCHI = V * SH / DSQRT(8.D0)
      VPHI = DSQRT(V**2 - 8.D0*VCHI**2)
      LAMBDA5 = 2.D0/3.D0/VPHI**2 * (MH5**2 - M1*VPHI**2/4.D0/VCHI
     .     - 12.D0*M2*VCHI - 8.D0*LAMBDA3*VCHI**2)
      MU3SQ = 2.D0/3.D0*MH5**2 + M1*VPHI**2/12.D0/VCHI 
     .     - 2.D0*LAMBDA2*VPHI**2 - 28.D0/3.D0*LAMBDA3*VCHI**2
     .     - 12.D0*LAMBDA4*VCHI**2 - 2.D0*M2*VCHI

      M12SQ = DSQRT(3.D0)/2.D0*VPHI 
     .     * (-M1 + 4.D0*(2.D0*LAMBDA2 - LAMBDA5)*VCHI)
      M22SQ = M1*VPHI**2/4.D0/VCHI - 6.D0*M2*VCHI
     .     + 8.D0*(LAMBDA3 + 3.D0*LAMBDA4)*VCHI**2
      LAMBDA1 = 1.D0/8.D0/VPHI**2 * (MHL**2 
     .     + M12SQ**2/(M22SQ - MHL**2))
      MU2SQ = -4.D0*LAMBDA1*VPHI**2 
     .     - 3.D0*(2.D0*LAMBDA2 - LAMBDA5)*VCHI**2
     .     + 3.D0/2.D0*M1*VCHI

C Check that the inputted parameters give the right minimum
      ZETA = 1.D0/3.D0
      OMEGA = 1.D0/2.D0
      SIGMA = DSQRT(3.D0)/4.D0
      RHO = 2.D0/DSQRT(3.D0)
      CALL GETMINIMA
C Identify the deepest of the good extrema:
      GOODVPHI = AA(1)
      GOODVCHI = BB(1)/DSQRT(3.D0)
      GOODV = VV(1)
      IF (VV(2).LT.GOODV) THEN
         GOODVPHI = AA(2)
         GOODVCHI = BB(2)/DSQRT(3.D0)
         GOODV = VV(2)
      ELSE IF (VV(3).LT.GOODV) THEN
         GOODVPHI = AA(3)
         GOODVCHI = BB(3)/DSQRT(3.D0)
         GOODV = VV(3)
      ENDIF
      IF (DABS(VCHI-GOODVCHI).GT.1D-6) THEN
         INPUTOK = 0
      ENDIF

      RETURN
      END


C==================================================================
      SUBROUTINE READIN5
C Reads in the inputs for Input Set 5 and computes the rest of the 
C Lagrangian parameters in LPARAMS.
C Input Set 5 is:
C     mh, mH, sin(thetaH), sin(alpha), lambda2, lambda3, lambda4, lambda5
C     with v = (sqrt(2) G_F)^{-1/2} ~ 246 GeV set by INITIALIZE_SM.
C
      IMPLICIT NONE
C Common blocks:
      DOUBLE PRECISION MU2SQ,MU3SQ,LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,
     .     LAMBDA5,M1,M2
      COMMON/LPARAMS/MU2SQ,MU3SQ,LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,
     .     LAMBDA5,M1,M2
      INTEGER UNIOK, BFBOK, POSMSQOK, MINOK, INPUTOK
      COMMON/FLAGS/UNIOK, BFBOK, POSMSQOK, MINOK, INPUTOK
      DOUBLE PRECISION V, MZ, MW, MTPOLE, MBMB, MCMC, MS, MTAU, MMU,
     .     ALPHAEM,ALSMZ, VCB, VUS, VUB, GAMZ, GAMW
      COMMON/SM/V,MZ,MW,MTPOLE,MBMB,MCMC,MS,MTAU,MMU,
     .     ALPHAEM,ALSMZ,VCB,VUS,VUB,GAMZ,GAMW
      INTEGER INPUTMODE, SILENT
      COMMON/CONTROLS/INPUTMODE,SILENT
      DOUBLE PRECISION MH
      INTEGER INPUTSET
      COMMON/INPUT/MH,INPUTSET
      DOUBLE PRECISION IMHL, IMHH, IMH3, IMH5, ISH, ISA
      COMMON/INPUT3/IMHL,IMHH,IMH3,IMH5,ISH,ISA
      DOUBLE PRECISION ZETA, OMEGA, SIGMA, RHO
      COMMON/HYPER/ZETA,OMEGA,SIGMA,RHO
      DOUBLE PRECISION AA, BB, VV
      DIMENSION AA(6), BB(6), VV(6)
      COMMON/MINIMA/AA,BB,VV
C Local variables:
      DOUBLE PRECISION VCHI, VPHI
      DOUBLE PRECISION M11SQ, M12SQ, M22SQ
      DOUBLE PRECISION GOODVPHI, GOODVCHI, GOODV
      DOUBLE PRECISION MHL, MHH, SH, SA, CA

      IF (INPUTMODE.EQ.1) THEN
C Read inputs from the user (interactive mode)
         PRINT *
         PRINT *, " Enter the 8 parameters on one line (no commas): "
         READ *, IMHL, IMHH, ISH, ISA, LAMBDA2, LAMBDA3, 
     .        LAMBDA4, LAMBDA5
      ENDIF
      
      MHL = IMHL
      MHH = IMHH
      SH = ISH
      SA = ISA

      INPUTOK = 1

      IF (SH.GE.1.D0.OR.SH.LE.0.D0) THEN
         INPUTOK = 0
      ENDIF
      IF (SA.GE.1.D0.OR.SA.LE.-1.D0) THEN
         INPUTOK = 0
      ENDIF

      VCHI = V * SH / DSQRT(8.D0)
      VPHI = DSQRT(V**2 - 8.D0*VCHI**2)
      CA = DSQRT(1.D0-SA**2)

      M11SQ = MHL**2*CA**2 + MHH**2*SA**2
      M22SQ = MHL**2*SA**2 + MHH**2*CA**2
      M12SQ = CA*SA*(MHH**2 - MHL**2)

      LAMBDA1 = M11SQ/8.D0/VPHI**2
      M1 = -2.D0*M12SQ/DSQRT(3.D0)/VPHI 
     .     + 4.D0*(2.D0*LAMBDA2-LAMBDA5)*VCHI
      M2 = (-M22SQ*VCHI + 8.D0*(LAMBDA3+3.D0*LAMBDA4)*VCHI**3
     .     + M1*VPHI**2/4.D0)/6.D0/VCHI**2
      MU3SQ = M22SQ - (2.D0*LAMBDA2-LAMBDA5)*VPHI**2 
     .     - 12.D0*(LAMBDA3+3.D0*LAMBDA4)*VCHI**2 + 12.D0*M2*VCHI
      MU2SQ = -4.D0*LAMBDA1*VPHI**2 
     .     - 3.D0*(2.D0*LAMBDA2-LAMBDA5)*VCHI**2 
     .     + 3.D0/2.D0*M1*VCHI

C Check that the inputted parameters give the right minimum
      ZETA = 1.D0/3.D0
      OMEGA = 1.D0/2.D0
      SIGMA = DSQRT(3.D0)/4.D0
      RHO = 2.D0/DSQRT(3.D0)
      CALL GETMINIMA
C Identify the deepest of the good extrema:
      GOODVPHI = AA(1)
      GOODVCHI = BB(1)/DSQRT(3.D0)
      GOODV = VV(1)
      IF (VV(2).LT.GOODV) THEN
         GOODVPHI = AA(2)
         GOODVCHI = BB(2)/DSQRT(3.D0)
         GOODV = VV(2)
      ELSE IF (VV(3).LT.GOODV) THEN
         GOODVPHI = AA(3)
         GOODVCHI = BB(3)/DSQRT(3.D0)
         GOODV = VV(3)
      ENDIF
      IF (DABS(VCHI-GOODVCHI).GT.1D-6) THEN
         INPUTOK = 0
      ENDIF

      RETURN
      END


C==================================================================
      SUBROUTINE READIN6
C Reads in the inputs for Input Set 6 and computes the rest of the 
C Lagrangian parameters in LPARAMS.
C Input Set 6 is:
C     mh, m5, sin(thetaH), lambda2, lambda3, lambda4, lambda5, M2
C     with v = (sqrt(2) G_F)^{-1/2} ~ 246 GeV set by INITIALIZE_SM.
C
      IMPLICIT NONE
C Common blocks:
      DOUBLE PRECISION MU2SQ,MU3SQ,LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,
     .     LAMBDA5,M1,M2
      COMMON/LPARAMS/MU2SQ,MU3SQ,LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,
     .     LAMBDA5,M1,M2
      INTEGER UNIOK, BFBOK, POSMSQOK, MINOK, INPUTOK
      COMMON/FLAGS/UNIOK, BFBOK, POSMSQOK, MINOK, INPUTOK
      DOUBLE PRECISION V, MZ, MW, MTPOLE, MBMB, MCMC, MS, MTAU, MMU,
     .     ALPHAEM,ALSMZ, VCB, VUS, VUB, GAMZ, GAMW
      COMMON/SM/V,MZ,MW,MTPOLE,MBMB,MCMC,MS,MTAU,MMU,
     .     ALPHAEM,ALSMZ,VCB,VUS,VUB,GAMZ,GAMW
      INTEGER INPUTMODE, SILENT
      COMMON/CONTROLS/INPUTMODE,SILENT
      DOUBLE PRECISION MH
      INTEGER INPUTSET
      COMMON/INPUT/MH,INPUTSET
      DOUBLE PRECISION IMHL, IMHH, IMH3, IMH5, ISH, ISA
      COMMON/INPUT3/IMHL,IMHH,IMH3,IMH5,ISH,ISA
      DOUBLE PRECISION ZETA, OMEGA, SIGMA, RHO
      COMMON/HYPER/ZETA,OMEGA,SIGMA,RHO
      DOUBLE PRECISION AA, BB, VV
      DIMENSION AA(6), BB(6), VV(6)
      COMMON/MINIMA/AA,BB,VV
C Local variables:
      DOUBLE PRECISION VCHI, VPHI
      DOUBLE PRECISION M12SQ, M22SQ
      DOUBLE PRECISION GOODVPHI, GOODVCHI, GOODV
      DOUBLE PRECISION MHL, MH5, SH, CH

      IF (INPUTMODE.EQ.1) THEN
C Read inputs from the user (interactive mode)
         PRINT *
         PRINT *, " Enter the 8 parameters on one line (no commas): "
         READ *, IMHL, IMH5, ISH, LAMBDA2, LAMBDA3, LAMBDA4, 
     .        LAMBDA5, M2
      ENDIF
      
      MHL = IMHL
      MH5 = IMH5
      SH = ISH
      CH = DSQRT(1.D0-SH**2)

      INPUTOK = 1

      IF (SH.GE.1.D0.OR.SH.LE.0.D0) THEN
         INPUTOK = 0
      ENDIF

      VCHI = V * SH / DSQRT(8.D0)
      VPHI = DSQRT(V**2 - 8.D0*VCHI**2)
      M1 = DSQRT(2.D0)*MH5**2/V*SH/CH**2
     .     - 3.D0/DSQRT(2.D0)*LAMBDA5*V*SH
     .     - DSQRT(2.D0)*LAMBDA3*V*SH**3/CH**2
     .     - 6.D0*M2*SH**2/CH**2
      MU3SQ = 2.D0/3.D0*MH5**2 + M1*VPHI**2/12.D0/VCHI 
     .     - 2.D0*LAMBDA2*VPHI**2 - 28.D0/3.D0*LAMBDA3*VCHI**2
     .     - 12.D0*LAMBDA4*VCHI**2 - 2.D0*M2*VCHI

      M12SQ = DSQRT(3.D0)/2.D0*VPHI 
     .     * (-M1 + 4.D0*(2.D0*LAMBDA2 - LAMBDA5)*VCHI)
      M22SQ = M1*VPHI**2/4.D0/VCHI - 6.D0*M2*VCHI
     .     + 8.D0*(LAMBDA3 + 3.D0*LAMBDA4)*VCHI**2
      LAMBDA1 = 1.D0/8.D0/VPHI**2 * (MHL**2 
     .     + M12SQ**2/(M22SQ - MHL**2))
      MU2SQ = -4.D0*LAMBDA1*VPHI**2 
     .     - 3.D0*(2.D0*LAMBDA2 - LAMBDA5)*VCHI**2
     .     + 3.D0/2.D0*M1*VCHI

C Check that the inputted parameters give the right minimum
      ZETA = 1.D0/3.D0
      OMEGA = 1.D0/2.D0
      SIGMA = DSQRT(3.D0)/4.D0
      RHO = 2.D0/DSQRT(3.D0)
      CALL GETMINIMA
C Identify the deepest of the good extrema:
      GOODVPHI = AA(1)
      GOODVCHI = BB(1)/DSQRT(3.D0)
      GOODV = VV(1)
      IF (VV(2).LT.GOODV) THEN
         GOODVPHI = AA(2)
         GOODVCHI = BB(2)/DSQRT(3.D0)
         GOODV = VV(2)
      ELSE IF (VV(3).LT.GOODV) THEN
         GOODVPHI = AA(3)
         GOODVCHI = BB(3)/DSQRT(3.D0)
         GOODV = VV(3)
      ENDIF
      IF (DABS(VCHI-GOODVCHI).GT.1D-6) THEN
         INPUTOK = 0
      ENDIF

      RETURN
      END
