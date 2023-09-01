! modified by Ameen on 2018-08-05
C=======================================================================
C GMCALC: a calculator for the Georgi-Machacek model 
C (with the most general custodial-symmetry-invariant scalar potential)
C http://people.physics.carleton.ca/~logan/gmcalc/
C========================================================================

      PROGRAM GMHB5 
C Program to evaluate the GM model against direct constraints from
C HiggsBounds 5. See gmhbhs.f for running GMCALC with HiggsBounds 4
C and HiggsSignals.
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
C Common blocks for couplings:
      DOUBLE PRECISION MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
      COMMON/PHYSPARAMS/MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
      DOUBLE PRECISION KVL, KFL, KGAML, KZGAML, DKGAML, DKZGAML
      COMMON/KAPPASL/KVL,KFL,KGAML,KZGAML,DKGAML,DKZGAML
      DOUBLE PRECISION KVH, KFH, KGAMH, KZGAMH, DKGAMH, DKZGAMH
      COMMON/KAPPASH/KVH,KFH,KGAMH,KZGAMH,DKGAMH,DKZGAMH
      DOUBLE PRECISION KW3, KZ3, KF3
      COMMON/KAPPAS3/KW3,KZ3,KF3
      DOUBLE PRECISION KW5, KZ5, KF5
      COMMON/KAPPAS5/KW5,KZ5,KF5
C Common blocks for decay BRs and total widths:
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
C Variables and common block for HiggsBounds5-calling subroutine:
      integer :: nHzero, nHplus
      integer :: HBresult(7), chan(7), ncombined(7)
      double precision :: obsratio(7)
      COMMON/HB5VARS/obsratio,nHzero,nHplus,HBresult,chan,ncombined
C Variables and common block for HiggsSignals2:
      integer :: nobs
      double precision :: Chisq_mu, Chisq_mh, Chisq, Pvalue
      COMMON/HS2VARS/Chisq_mu,Chisq_mh,Chisq,Pvalue,nobs
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
c      INPUTSET = 6
      INPUTSET = 4
C==================================================================
C SILENT = 0: echo the inputs to the screen
C SILENT = 1: don't echo the inputs to the screen
      SILENT = 1
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
            IMH5 = 1000.D0
            ISH = 0.44D0
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
            IMH5 = 60.D0
            ISH = 0.1D0
            LAMBDA2 = 0.08D0*(IMH5/100.D0)
            LAMBDA3 = -1.5D0
            LAMBDA4 = -LAMBDA3
            LAMBDA5 = -4.D0*LAMBDA2
            M2 = 10.D0
         ELSE
            PRINT *, "INPUTSET = ", INPUTSET, "is not a valid option."
            PRINT *
            GOTO 10
         ENDIF
      ENDIF

C==================================================================
C Initialize LoopTools
      CALL LTSTARTER
C Initialize HiggsBounds
      nHzero = 4
      nHplus = 2
      CALL initialize_HiggsBounds(nHzero, nHplus, "LandH")
C Initialize HiggsSignals (nparam = number of scanned model params for p-value)
      CALL initialize_HiggsSignals_LHC13(nHzero, nHplus)
      CALL setup_output_level(0)
      CALL setup_pdf(2)
      CALL setup_nparam(0)

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
            GOTO 10
         ENDIF
      ELSE
         PRINT *, "Bad input set!"
         GOTO 10
      ENDIF

C Call HiggsBounds and HiggsSignals via the subroutine in src/gmhbhs.f
      CALL CALLHBHS

      print*,' '      
      print*,'*************    HiggsBounds Results  **************'
      print *
      print*,'Allowed combining all Higgses? ', HBresult(1)
      print*,'Most sensitive channel: ', chan(1), ' (see Key.dat)'
      print*,'Ratio of predicted rate to limit: ', obsratio(1)
      print*,'Number of Higgses contributing: ', ncombined(1)
      print *
      print *, 'h, mass = ', MHL, ' GeV'
      print *, 'Allowed by HiggsBounds? ', HBresult(2)
      print *, 'Most sensitive channel: ', chan(2)
      print *, 'Ratio of predicted rate to limit: ', obsratio(2)
      print *
      print *, 'H, mass = ', MHH, ' GeV'
      print *, 'Allowed by HiggsBounds? ', HBresult(3)
      print *, 'Most sensitive channel: ', chan(3)
      print *, 'Ratio of predicted rate to limit: ', obsratio(3)
      print *
      print *, 'H30, mass = ', MH3, ' GeV'
      print *, 'Allowed by HiggsBounds? ', HBresult(4)
      print *, 'Most sensitive channel: ', chan(4)
      print *, 'Ratio of predicted rate to limit: ', obsratio(4)
      print *
      print *, 'H50, mass = ', MH5, ' GeV'
      print *, 'Allowed by HiggsBounds? ', HBresult(5)
      print *, 'Most sensitive channel: ', chan(5)
      print *, 'Ratio of predicted rate to limit: ', obsratio(5)
      print *
      print *, 'H3+, mass = ', MH3, ' GeV'
      print *, 'Allowed by HiggsBounds? ', HBresult(6)
      print *, 'Most sensitive channel: ', chan(6)
      print *, 'Ratio of predicted rate to limit: ', obsratio(6)
      print *
      print *, 'H5+, mass = ', MH5, ' GeV'
      print *, 'Allowed by Higgsbounds? ', HBresult(7)
      print *, 'Most sensitive channel: ', chan(7)
      print *, 'Ratio of predicted rate to limit: ', obsratio(7)
      print *
      print*,'*************    HiggsSignals Results  **************'
      print *
      print *, 'Chisq_mu = ', Chisq_mu
      print *, 'Chisq_mh = ', Chisq_mh
      print *, 'Chisq = ', Chisq
      print *, 'nobs = ', nobs
      print *, 'Pvalue = ', Pvalue
      print *
      print*,'****************************************************'

C Finish up HiggsBounds and HiggsSignals
      CALL finish_HiggsBounds()
      CALL finish_HiggsSignals()
C Finish up LoopTools
      CALL LTENDER

 10   STOP
      END


