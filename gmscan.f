C=======================================================================
C GMCALC: a calculator for the Georgi-Machacek model 
C (with the most general custodial-symmetry-invariant scalar potential)
C http://people.physics.carleton.ca/~logan/gmcalc/
C========================================================================

      PROGRAM GMSCAN
C Program to perform a scan over the allowed GM model parameter ranges,
C subject to theoretical and indirect experimental constraints, with 
C mh set to 125 GeV.  The user sets the maximum value of MU3SQ.
C The output consists of a representative selection
C of masses, couplings, and decay branching ratios.  The full set of
C calculated observables is available to be written out.
      IMPLICIT NONE
C Common blocks:
      DOUBLE PRECISION MU2SQ,MU3SQ,LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,
     .     LAMBDA5,M1,M2
      COMMON/LPARAMS/MU2SQ,MU3SQ,LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,
     .     LAMBDA5,M1,M2
      DOUBLE PRECISION MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
      COMMON/PHYSPARAMS/MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
      INTEGER UNIOK, BFBOK, POSMSQOK, MINOK, INPUTOK
      COMMON/FLAGS/UNIOK, BFBOK, POSMSQOK, MINOK, INPUTOK
      DOUBLE PRECISION KVL, KFL, KGAML, KZGAML, DKGAML, DKZGAML
      COMMON/KAPPASL/KVL,KFL,KGAML,KZGAML,DKGAML,DKZGAML
      DOUBLE PRECISION KVH, KFH, KGAMH, KZGAMH, DKGAMH, DKZGAMH
      COMMON/KAPPASH/KVH,KFH,KGAMH,KZGAMH,DKGAMH,DKZGAMH
      DOUBLE PRECISION V, MZ, MW, MTPOLE, MBMB, MCMC, MS, MTAU, MMU,
     .     ALPHAEM, ALSMZ, VCB, VUS, VUB, GAMZ, GAMW
      COMMON/SM/V,MZ,MW,MTPOLE,MBMB,MCMC,MS,MTAU,MMU,
     .     ALPHAEM,ALSMZ,VCB,VUS,VUB,GAMZ,GAMW
      INTEGER INPUTMODE, SILENT
      COMMON/CONTROLS/INPUTMODE,SILENT
      DOUBLE PRECISION MH
      INTEGER INPUTSET
      COMMON/INPUT/MH,INPUTSET
      DOUBLE PRECISION IMHL, IMHH, IMH3, IMH5, ISH, ISA
      COMMON/INPUT3/IMHL,IMHH,IMH3,IMH5,ISH,ISA
C Common block for indirect constraints
      DOUBLE PRECISION RBSMM, SPARAM
      INTEGER BSMMOK, SPAROK, BSGAMLOOSEOK, BSGAMTIGHTOK
      COMMON/INDIR/RBSMM, SPARAM, BSMMOK, SPAROK, 
     .     BSGAMLOOSEOK, BSGAMTIGHTOK
C Common blocks for decay BRs and total widths
      DOUBLE PRECISION HLBRB, HLBRTA, HLBRMU, HLBRS, HLBRC, HLBRT,
     .     HLBRG, HLBRGA, HLBRZGA, HLBRW, HLBRZ, 
     .     HLBRWH3P, HLBRZH3N,
     .     HLBRH3N, HLBRH3P, HLBRH5N, HLBRH5P, HLBRH5PP, HLWDTH
      COMMON/HLBRS/HLBRB, HLBRTA, HLBRMU, HLBRS, HLBRC, HLBRT,
     .     HLBRG, HLBRGA, HLBRZGA, HLBRW, HLBRZ, 
     .     HLBRWH3P, HLBRZH3N,
     .     HLBRH3N, HLBRH3P, HLBRH5N, HLBRH5P, HLBRH5PP, HLWDTH
      DOUBLE PRECISION HHBRB, HHBRTA, HHBRMU, HHBRS, HHBRC, HHBRT,
     .     HHBRG, HHBRGA, HHBRZGA, HHBRW, HHBRZ, 
     .     HHBRWH3P, HHBRZH3N,
     .     HHBRHL, HHBRH3N, HHBRH3P, HHBRH5N, HHBRH5P, HHBRH5PP, 
     .     HHWDTH
      COMMON/HHBRS/HHBRB, HHBRTA, HHBRMU, HHBRS, HHBRC, HHBRT,
     .     HHBRG, HHBRGA, HHBRZGA, HHBRW, HHBRZ, 
     .     HHBRWH3P, HHBRZH3N,
     .     HHBRHL, HHBRH3N, HHBRH3P, HHBRH5N, HHBRH5P, HHBRH5PP, 
     .     HHWDTH
      DOUBLE PRECISION H3NBRB, H3NBRTA, H3NBRMU, H3NBRS, H3NBRC, H3NBRT,
     .     H3NBRZHL, H3NBRZHH, H3NBRZH5N, H3NBRWH5P,
     .     H3NBRG, H3NBRGA, H3NBRZGA, 
     .     H3NWDTH
      COMMON/H3NBRS/H3NBRB, H3NBRTA, H3NBRMU, H3NBRS, H3NBRC, H3NBRT,
     .     H3NBRZHL, H3NBRZHH, H3NBRZH5N, H3NBRWH5P,
     .     H3NBRG, H3NBRGA, H3NBRZGA,
     .     H3NWDTH
      DOUBLE PRECISION H3PBRBC, H3PBRTA, H3PBRMU, H3PBRSU,
     .     H3PBRCS, H3PBRTB, H3PBRBU,
     .     H3PBRWHL, H3PBRWHH, H3PBRZH5P, H3PBRWH5N, H3PBRWH5PP,
     .     H3PBRWGA,
     .     H3PWDTH
      COMMON/H3PBRS/H3PBRBC, H3PBRTA, H3PBRMU, H3PBRSU,
     .     H3PBRCS, H3PBRTB, H3PBRBU,
     .     H3PBRWHL, H3PBRWHH, H3PBRZH5P, H3PBRWH5N, H3PBRWH5PP,
     .     H3PBRWGA,
     .     H3PWDTH
      DOUBLE PRECISION H5NBRGA, H5NBRZGA, H5NBRW, H5NBRZ,
     .     H5NBRZH3N, H5NBRWH3P, 
     .     H5NBRH3N, H5NBRH3P,
     .     H5NWDTH
      COMMON/H5NBRS/H5NBRGA, H5NBRZGA, H5NBRW, H5NBRZ,
     .     H5NBRZH3N, H5NBRWH3P, 
     .     H5NBRH3N, H5NBRH3P,
     .     H5NWDTH
      DOUBLE PRECISION H5PBRWZ, H5PBRZH3P, H5PBRWH3N, H5PBRH3PN, 
     .     H5PBRWGA, H5PWDTH
      COMMON/H5PBRS/H5PBRWZ, H5PBRZH3P, H5PBRWH3N, H5PBRH3PN, 
     .     H5PBRWGA, H5PWDTH
      DOUBLE PRECISION H5PPBRWW, H5PPBRWH3, H5PPBRH3P, H5PPWDTH
      COMMON/H5PPBRS/H5PPBRWW, H5PPBRWH3, H5PPBRH3P, H5PPWDTH
      DOUBLE PRECISION TOPBRW, TOPBRH3P, TOPWDTH
      COMMON/TOPBRS/TOPBRW, TOPBRH3P, TOPWDTH
      INTEGER OFFSHELL, QCDCORRS
      COMMON/DECAYFLAGS/OFFSHELL, QCDCORRS
C Local variables:
      INTEGER POINTS, TRIALS, NPOINTS
      DOUBLE PRECISION X
      DIMENSION X(7)
      DOUBLE PRECISION LAMBDA2MAX, LAMBDA5MIN, LAMBDA5MAX, M1MAX, M2MAX
      DOUBLE PRECISION MU3MAX
      DOUBLE PRECISION MINMASS
      INTEGER SMLIKE
      DOUBLE PRECISION KV, KF, KGAM, KZGAM, DKGAM, DKZGAM
      DOUBLE PRECISION PI
      PI = 4.D0*DATAN(1.D0)
      CALL PRINT_BANNER
C==================================================================
C This initialization call must always come before anything else!
      CALL INITIALIZE_SM
C==================================================================

C================== user-defined options ==========================
C Decay flags:  Set these to zero to turn off offshell decays
C and/or QCD corrections.  To keep everything on, set both flags to 1.
      OFFSHELL = 1
      QCDCORRS = 1
C==================================================================
C Parameter scan will run until NPOINTS good points are found:
      NPOINTS = 100
C==================================================================
C Set the light Higgs mass and the maximum value of sqrt(MU3SQ) 
C to be scanned over.  
      MH = 125.D0
      MU3MAX = 1200.D0
C==================================================================

C============= flags that should not be changed lightly ===========
C The following flags must be set this way for the scans to 
C function properly.  INPUTMODE = 0 reads the parameters from the 
C main program.  SILENT = 1 avoids dumping details of each of the 
C scan points to the screen.
      INPUTMODE = 0
      SILENT = 1
C==================================================================
C Since in our scan we want to fix mh = 125 GeV, we must use an 
C INPUTSET in which this is taken as an input: i.e., 2, 3, 4, 5, or 6.
      INPUTSET = 2
C==================================================================
      
      CALL LTSTARTER
      OPEN (UNIT = 90, FILE = 'scan-output.data')
      WRITE (90,*) "# kappa_V   kappa_f   kappa_gamma   m_5", 
     .     "    BR(H5+ -> WZ)    H5+ width"

C Start the loop over scan points
      TRIALS = 0
      POINTS = 0
 400  CALL RANDOM_NUMBER(X)
      TRIALS = TRIALS + 1
C============ setting up the scan ranges ==========================
      X(1) = -200.D0 + X(1)*(200.D0+MU3MAX)
      MU3SQ = X(1)*DABS(X(1))
      LAMBDA3 = -0.5D0*PI + X(3)*11.D0/10.D0*PI
      IF (LAMBDA3.LT.0) THEN
         LAMBDA4 = -LAMBDA3 + X(4)*(4.D0/11.D0*LAMBDA3 
     .        + 2.D0/11.D0*PI)
      ELSE
         LAMBDA4 = -LAMBDA3/3.D0 + X(4)*(-10.D0/33.D0*LAMBDA3 
     .        + 2.D0/11.D0*PI)
      ENDIF
      LAMBDA2MAX = DSQRT(4.D0*PI**2 
     .     - 2.D0*PI*(7.D0*LAMBDA3 + 11.D0*LAMBDA4))/3.D0
      LAMBDA2 = -LAMBDA2MAX + X(2)*2.D0*LAMBDA2MAX
      LAMBDA5MIN = -2.D0*PI + LAMBDA2
      LAMBDA5MAX = 2.D0*PI + LAMBDA2
      LAMBDA5 = LAMBDA5MIN + X(5)*(LAMBDA5MAX-LAMBDA5MIN)
      M1MAX = 3.5D0*DSQRT(DABS(MU3SQ))
      IF (M1MAX.LT.3500.D0) THEN
         M1MAX = 3500.D0
      ENDIF
      M1 = X(6)*M1MAX
      M2MAX = 1.3D0*DSQRT(DABS(MU3SQ))
      IF (M2MAX.LT.250.D0) THEN
         M2MAX = 250.D0
      ENDIF
      M2 = -M2MAX + X(7)*2.D0*M2MAX
C==================================================================
 500  CALL LOAD_INPUTS
      IF (INPUTOK.EQ.1) THEN
         CALL THYCHECK
         IF (UNIOK*BFBOK*MINOK.EQ.1) THEN
            CALL CALCPHYS
            IF (VCHI.LT.0.D0) THEN
               M1 = -M1
               M2 = -M2
               GOTO 500
            ENDIF
            CALL CALCINDIR
C Impose the indirect constraints: S parameter and "loose" b -> s gamma
            IF (SPAROK.EQ.1.AND.BSGAMLOOSEOK.EQ.1) THEN
               POINTS = POINTS + 1
               CALL HLCOUPS
               CALL HHCOUPS
C Check whether HL or HH is the 125 GeV Higgs.  If they are both within
C +- 0.5 GeV, take the SM-like one to be the one with larger |KV|.
               IF (MHH.GE.MH+0.5D0) THEN
                  SMLIKE = 1
               ELSE IF (MHL.LE.MH-0.5D0) THEN
                  SMLIKE = 2
               ELSE IF (DABS(KVL).GE.DABS(KVH)) THEN
                  SMLIKE = 1
               ELSE
                  SMLIKE = 2
               ENDIF
               IF (SMLIKE.EQ.1) THEN
                  KV = KVL
                  KF = KFL
                  KGAM = KGAML
                  KZGAM = KZGAML
                  DKGAM = DKGAML
                  DKZGAM = DKZGAML
               ELSE
                  KV = KVH
                  KF = KFH
                  KGAM = KGAMH
                  KZGAM = KZGAMH
                  DKGAM = DKGAMH
                  DKZGAM = DKZGAMH
               ENDIF
               CALL CALCDECAYS
C MINMASS is the mass of the lightest *new* scalar:
               CALL GET_MINMASS(MINMASS)
               WRITE (90, *) KV, KF, KGAM, MH5, H5PBRWZ, H5PWDTH
            ENDIF
         ENDIF
      ENDIF
      IF (POINTS.LT.NPOINTS) THEN
         GOTO 400
      ENDIF
C End of the loop over scan points

      PRINT *, "Number of trials = ", TRIALS
      PRINT *, "Number of good points = ", POINTS
      CLOSE (UNIT = 90)
      CALL LTENDER

      STOP
      END



