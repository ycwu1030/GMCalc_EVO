C=======================================================================
C GMCALC: a calculator for the Georgi-Machacek model 
C (with the most general custodial-symmetry-invariant scalar potential)
C http://people.physics.carleton.ca/~logan/gmcalc/
C========================================================================

      PROGRAM GMPOINT
C Program to compute the spectrum, check theoretical and indirect
C constraints, compute the h and H coupling modification factors, 
C and compute the decay BRs and total width of all the scalars.
      IMPLICIT NONE
C Common blocks:
      DOUBLE PRECISION MU2SQ,MU3SQ,LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,
     .     LAMBDA5,M1,M2
      COMMON/LPARAMS/MU2SQ,MU3SQ,LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,
     .     LAMBDA5,M1,M2
      DOUBLE PRECISION V, MZ, MW, MTPOLE, MBMB, MCMC, MS, MTAU, MMU,
     .     ALPHAEM, ALSMZ, VCB, VUS, VUB, GAMZ, GAMW
      COMMON/SM/V,MZ,MW,MTPOLE,MBMB,MCMC,MS,MTAU,MMU,
     .     ALPHAEM,ALSMZ,VCB,VUS,VUB,GAMZ,GAMW
      INTEGER UNIOK, BFBOK, POSMSQOK, MINOK, INPUTOK
      COMMON/FLAGS/UNIOK, BFBOK, POSMSQOK, MINOK, INPUTOK
      INTEGER INPUTMODE, SILENT
      COMMON/CONTROLS/INPUTMODE,SILENT
      DOUBLE PRECISION MH
      INTEGER INPUTSET
      COMMON/INPUT/MH,INPUTSET
      DOUBLE PRECISION IMHL, IMHH, IMH3, IMH5, ISH, ISA
      COMMON/INPUT3/IMHL,IMHH,IMH3,IMH5,ISH,ISA
      INTEGER OFFSHELL, QCDCORRS
      COMMON/DECAYFLAGS/OFFSHELL, QCDCORRS
      CALL PRINT_BANNER
C==================================================================
C This initialization call must always come before anything else!
      CALL INITIALIZE_SM
C==================================================================

C================== user-defined options ==========================
C INPUTMODE = 0: take inputs from this program - enter values below
C INPUTMODE = 1: read inputs from the command line in interactive mode
      INPUTMODE = 0
C==================================================================
C Options for input sets:
C INPUTSET = 1: mu3^2, lambda1, lambda2, lambda3, lambda4, lambda5, M1, M2
C INPUTSET = 2: mu3^2, mh, lambda2, lambda3, lambda4, lambda5, M1, M2
C INPUTSET = 3: mh, mH, m3, m5, sin(thetaH), sin(alpha), M1, M2
C INPUTSET = 4: mh, m5, sin(thetaH), lambda2, lambda3, lambda4, M1, M2
C INPUTSET = 5: mh, mH, sin(thetaH), sin(alpha), lambda2, lambda3, lambda4, lambda5
C INPUTSET = 6: mh, m5, sin(thetaH), lambda2, lambda3, lambda4, lambda5, M2
      INPUTSET = 4
C==================================================================
C SILENT = 0: echo the inputs to the screen
C SILENT = 1: don't echo the inputs to the screen
      SILENT = 0
C==================================================================
C Decay flags:  Set these to zero to turn off offshell decays
C and/or QCD corrections.  To keep everything on, set both flags to 1.
      OFFSHELL = 1
      QCDCORRS = 1
C==================================================================
C Modify these entries to set the input parameters (used if INPUTMODE = 0). 
      IF (INPUTMODE.EQ.0) THEN
         IF (INPUTSET.EQ.1) THEN
            MU3SQ = 90000.D0
            LAMBDA1 = 4.68D-2
            LAMBDA2 = 0.1D0
            LAMBDA3 = 0.1D0
            LAMBDA4 = 0.1D0
            LAMBDA5 = 0.1D0
            M1 = 100.D0
            M2 = 100.D0
         ELSE IF (INPUTSET.EQ.2) THEN
            MH = 125.D0
            MU3SQ = 250000.D0
            LAMBDA2 = 0.1D0
            LAMBDA3 = 0.1D0
            LAMBDA4 = 0.1D0
            LAMBDA5 = 0.1D0
            M1 = 100.D0
            M2 = 100.D0
         ELSE IF (INPUTSET.EQ.3) THEN
            IMHL = 125.D0
            IMHH = 288.268237D0
            IMH3 = 304.221605D0
            IMH5 = 339.748616D0
            ISH = 0.194487374D0
            ISA = -0.303281383D0
            M1 = 100.D0
            M2 = 100.D0
         ELSE IF (INPUTSET.EQ.4) THEN
C Here the user can set all 8 inputs independently, or vary only IMH5 =
C m_5 in [200,3000] GeV and ISH = s_H in (0,1) for the H5plane
C benchmark.
            IMHL = 125.D0
            IMH5 = 339.748616D0
            ISH = 0.194487374D0
            LAMBDA2 = 0.4D0*IMH5/1000.D0
            LAMBDA3 = -0.1D0
            LAMBDA4 = 0.2D0
            M1 = DSQRT(2.D0)*ISH*(IMH5**2+V**2)/V
            M2 = M1/6.D0
         ELSE IF (INPUTSET.EQ.5) THEN
            IMHL = 125.D0
            IMHH = 288.26779305953880D0
            ISH = 0.19448737400000182D0
            ISA = -0.30328279630518934D0
            LAMBDA2 = 0.1D0
            LAMBDA3 = 0.1D0
            LAMBDA4 = 0.1D0
            LAMBDA5 = 0.1D0
         ELSE IF (INPUTSET.EQ.6) THEN
            IMHL = 125.D0
            IMH5 = 339.748616D0
            ISH = 0.194487374D0
            LAMBDA2 = 0.1D0
            LAMBDA3 = 0.1D0
            LAMBDA4 = 0.1D0
            LAMBDA5 = 0.1D0
            M2 = 100.D0
         ELSE
            PRINT *, "INPUTSET = ", INPUTSET, "is not a valid option."
            PRINT *
            GOTO 10
         ENDIF
      ENDIF

C==================================================================
      CALL LTSTARTER

      CALL LOAD_INPUTS
      IF (INPUTOK.EQ.1) THEN
         CALL THYCHECK
         IF (UNIOK*BFBOK*MINOK.EQ.1) THEN
            CALL CALCPHYS
            CALL CALCINDIR
            CALL HLCOUPS
            CALL HHCOUPS
            CALL H3COUPS
            CALL H5COUPS
            CALL CALCDECAYS
            CALL PRINT_RESULTS
            CALL PRINT_HCOUPS
            CALL PRINT_DECAYS
         ELSE
            PRINT *, "Theory constraint(s) failed!"
         ENDIF
      ELSE
         PRINT *, "Bad input set!"
      ENDIF

      CALL LTENDER

 10   STOP
      END



