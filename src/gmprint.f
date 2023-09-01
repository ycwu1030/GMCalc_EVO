C==================================================================
C     Subroutine to print the banner
C==================================================================

      SUBROUTINE PRINT_BANNER
      IMPLICIT NONE
      PRINT *
      WRITE(*,1)
      WRITE(*,3) "                GMCALC v1.5.3:                   "
      WRITE(*,3) "   A calculator for the Georgi-Machacek model    "
      WRITE(*,3) "http://people.physics.carleton.ca/~logan/gmcalc/ "
      WRITE(*,1)
      PRINT *
 1    FORMAT(70('#'))
 3    FORMAT(2('#'),8X,A50,8X,2('#'))
 4    FORMAT(2('#'),7X,A52,7X,2('#'))
      RETURN
      END


C==================================================================
C     Subroutine to write out the decay BRs and total widths
C==================================================================

      SUBROUTINE PRINT_DECAYS
      IMPLICIT NONE
C Common blocks:
      DOUBLE PRECISION MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
      COMMON/PHYSPARAMS/MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
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

      PRINT *
      PRINT *, "   ===== GMCALC DECAYS ======   "
      PRINT *, "   format: mode  BR # partial width (GeV) "
      PRINT *, "           total width in GeV       "
      PRINT *, "   +hc means hermitian conjugate decay is included "
      PRINT *

      PRINT *, "Top quark decays"
      PRINT *, "W+ b ", TOPBRW, " # ", TOPBRW*TOPWDTH
      PRINT *, "H3+ b", TOPBRH3P, " # ", TOPBRH3P*TOPWDTH
      PRINT *, "WIDTH", TOPWDTH

      PRINT *
      PRINT *, "Decays of h, mass = ", MHL
      PRINT *, "BB   ", HLBRB, " # ", HLBRB*HLWDTH
      PRINT *, "TATA ", HLBRTA, " # ", HLBRTA*HLWDTH
      PRINT *, "MUMU ", HLBRMU, " # ", HLBRMU*HLWDTH
      PRINT *, "SS   ", HLBRS, " # ", HLBRS*HLWDTH
      PRINT *, "CC   ", HLBRC, " # ", HLBRC*HLWDTH
      PRINT *, "TT   ", HLBRT, " # ", HLBRT*HLWDTH
      PRINT *, "GG   ", HLBRG, " # ", HLBRG*HLWDTH
      PRINT *, "GAGA ", HLBRGA, " # ", HLBRGA*HLWDTH
      PRINT *, "ZGA  ", HLBRZGA, " # ", HLBRZGA*HLWDTH
      PRINT *, "WW   ", HLBRW, " # ", HLBRW*HLWDTH
      PRINT *, "ZZ   ", HLBRZ, " # ", HLBRZ*HLWDTH
      PRINT *, "WH3+ +hc", HLBRWH3P, " # ", HLBRWH3P*HLWDTH
      PRINT *, "ZH30 ", HLBRZH3N, " # ", HLBRZH3N*HLWDTH
      PRINT *, "H30  ", HLBRH3N, " # ", HLBRH3N*HLWDTH
      PRINT *, "H3+- ", HLBRH3P, " # ", HLBRH3P*HLWDTH
      PRINT *, "H50  ", HLBRH5N, " # ", HLBRH5N*HLWDTH
      PRINT *, "H5+- ", HLBRH5P, " # ", HLBRH5P*HLWDTH
      PRINT *, "H5++--", HLBRH5PP, " # ", HLBRH5PP*HLWDTH
      PRINT *, "WIDTH", HLWDTH

      PRINT *
      PRINT *, "Decays of H, mass = ", MHH
      PRINT *, "BB   ", HHBRB, " # ", HHBRB*HHWDTH
      PRINT *, "TATA ", HHBRTA, " # ", HHBRTA*HHWDTH
      PRINT *, "MUMU ", HHBRMU, " # ", HHBRMU*HHWDTH
      PRINT *, "SS   ", HHBRS, " # ", HHBRS*HHWDTH
      PRINT *, "CC   ", HHBRC, " # ", HHBRC*HHWDTH
      PRINT *, "TT   ", HHBRT, " # ", HHBRT*HHWDTH
      PRINT *, "GG   ", HHBRG, " # ", HHBRG*HHWDTH
      PRINT *, "GAGA ", HHBRGA, " # ", HHBRGA*HHWDTH
      PRINT *, "ZGA  ", HHBRZGA, " # ", HHBRZGA*HHWDTH
      PRINT *, "WW   ", HHBRW, " # ", HHBRW*HHWDTH
      PRINT *, "ZZ   ", HHBRZ, " # ", HHBRZ*HHWDTH
      PRINT *, "WH3+ +hc", HHBRWH3P, " # ", HHBRWH3P*HHWDTH
      PRINT *, "ZH30 ", HHBRZH3N, " # ", HHBRZH3N*HHWDTH
      PRINT *, "hh   ", HHBRHL, " # ", HHBRHL*HHWDTH
      PRINT *, "H30  ", HHBRH3N, " # ", HHBRH3N*HHWDTH
      PRINT *, "H3+- ", HHBRH3P, " # ", HHBRH3P*HHWDTH
      PRINT *, "H50  ", HHBRH5N, " # ", HHBRH5N*HHWDTH
      PRINT *, "H5+- ", HHBRH5P, " # ", HHBRH5P*HHWDTH
      PRINT *, "H5++--", HHBRH5PP, " # ", HHBRH5PP*HHWDTH
      PRINT *, "WIDTH", HHWDTH

      PRINT *
      PRINT *, "Decays of H3^0, mass = ", MH3
      PRINT *, "BB   ", H3NBRB, " # ", H3NBRB*H3NWDTH
      PRINT *, "TATA ", H3NBRTA, " # ", H3NBRTA*H3NWDTH
      PRINT *, "MUMU ", H3NBRMU, " # ", H3NBRMU*H3NWDTH
      PRINT *, "SS   ", H3NBRS, " # ", H3NBRS*H3NWDTH
      PRINT *, "CC   ", H3NBRC, " # ", H3NBRC*H3NWDTH
      PRINT *, "TT   ", H3NBRT, " # ", H3NBRT*H3NWDTH
      PRINT *, "GG   ", H3NBRG, " # ", H3NBRG*H3NWDTH
      PRINT *, "GAGA ", H3NBRGA, " # ", H3NBRGA*H3NWDTH
      PRINT *, "ZGA  ", H3NBRZGA, " # ", H3NBRZGA*H3NWDTH
      PRINT *, "Zh   ", H3NBRZHL, " # ", H3NBRZHL*H3NWDTH
      PRINT *, "ZH   ", H3NBRZHH, " # ", H3NBRZHH*H3NWDTH
      PRINT *, "ZH50 ", H3NBRZH5N, " # ", H3NBRZH5N*H3NWDTH
      PRINT *, "WH5+ +hc", H3NBRWH5P, " # ", H3NBRWH5P*H3NWDTH
      PRINT *, "WIDTH", H3NWDTH

      PRINT *
      PRINT *, "Decays of H3+, mass = ", MH3
      PRINT *, "BC   ", H3PBRBC, " # ", H3PBRBC*H3PWDTH
      PRINT *, "TAUNU", H3PBRTA, " # ", H3PBRTA*H3PWDTH
      PRINT *, "MU NU", H3PBRMU, " # ", H3PBRMU*H3PWDTH
      PRINT *, "SU   ", H3PBRSU, " # ", H3PBRSU*H3PWDTH
      PRINT *, "CS   ", H3PBRCS, " # ", H3PBRCS*H3PWDTH
      PRINT *, "TB   ", H3PBRTB, " # ", H3PBRTB*H3PWDTH
      PRINT *, "BU   ", H3PBRBU, " # ", H3PBRBU*H3PWDTH
      PRINT *, "W+h  ", H3PBRWHL, " # ", H3PBRWHL*H3PWDTH
      PRINT *, "W+H  ", H3PBRWHH, " # ", H3PBRWHH*H3PWDTH
      PRINT *, "ZH5+ ", H3PBRZH5P, " # ", H3PBRZH5P*H3PWDTH
      PRINT *, "W+H5^0", H3PBRWH5N, " # ", H3PBRWH5N*H3PWDTH
      PRINT *, "W-H5++", H3PBRWH5PP, " # ", H3PBRWH5PP*H3PWDTH
      PRINT *, "W+GA ", H3PBRWGA, " # ", H3PBRWGA*H3PWDTH
      PRINT *, "WIDTH", H3PWDTH

      PRINT *
      PRINT *, "Decays of H5^0, mass = ", MH5
      PRINT *, "GAGA ", H5NBRGA, " # ", H5NBRGA*H5NWDTH
      PRINT *, "ZGA  ", H5NBRZGA, " # ", H5NBRZGA*H5NWDTH
      PRINT *, "WW   ", H5NBRW, " # ", H5NBRW*H5NWDTH
      PRINT *, "ZZ   ", H5NBRZ, " # ", H5NBRZ*H5NWDTH
      PRINT *, "ZH3^0", H5NBRZH3N, " # ", H5NBRZH3N*H5NWDTH
      PRINT *, "WH3+ +hc", H5NBRWH3P, " # ", H5NBRWH3P*H5NWDTH
      PRINT *, "H3^0 ", H5NBRH3N, " # ", H5NBRH3N*H5NWDTH
      PRINT *, "H3+- ", H5NBRH3P, " # ", H5NBRH3P*H5NWDTH
      PRINT *, "WIDTH", H5NWDTH

      PRINT *
      PRINT *, "Decays of H5+, mass = ", MH5
      PRINT *, "WZ   ", H5PBRWZ, " # ", H5PBRWZ*H5PWDTH
      PRINT *, "ZH3+ ", H5PBRZH3P, " # ", H5PBRZH3P*H5PWDTH
      PRINT *, "WH3^0", H5PBRWH3N, " # ", H5PBRWH3N*H5PWDTH
      PRINT *, "H3+H3^0", H5PBRH3PN, " # ", H5PBRH3PN*H5PWDTH
      PRINT *, "W+GA ", H5PBRWGA, " # ", H5PBRWGA*H5PWDTH
      PRINT *, "WIDTH", H5PWDTH

      PRINT *
      PRINT *, "Decays of H5++, mass = ", MH5
      PRINT *, "W+W+ ", H5PPBRWW, " # ", H5PPBRWW*H5PPWDTH
      PRINT *, "W+H3+", H5PPBRWH3, " # ", H5PPBRWH3*H5PPWDTH
      PRINT *, "H3+H3+", H5PPBRH3P, " # ", H5PPBRH3P*H5PPWDTH
      PRINT *, "WIDTH", H5PPWDTH

      RETURN
      END

C==================================================================
C     Subroutine to write out the results in a tidy form
C==================================================================

      SUBROUTINE PRINT_RESULTS
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
C Common block for indirect constraints
      DOUBLE PRECISION RBSMM, SPARAM
      INTEGER BSMMOK, SPAROK, BSGAMLOOSEOK, BSGAMTIGHTOK
      COMMON/INDIR/RBSMM, SPARAM, BSMMOK, SPAROK, 
     .     BSGAMLOOSEOK, BSGAMTIGHTOK

      PRINT *
      PRINT *, "   ===== GMCALC RESULTS =====   "
      PRINT *

      IF (INPUTOK.EQ.1) THEN

      PRINT *, "Scalar potential parameters:"
      PRINT *, "mu_2^2 [GeV^2] = ", MU2SQ
      PRINT *, "mu_3^2 [GeV^2] = ", MU3SQ
      PRINT *, "lambda_1 =       ", LAMBDA1
      PRINT *, "lambda_2 =       ", LAMBDA2
      PRINT *, "lambda_3 =       ", LAMBDA3
      PRINT *, "lambda_4 =       ", LAMBDA4
      PRINT *, "lambda_5 =       ", LAMBDA5
      PRINT *, "M1 [GeV] =       ", M1
      PRINT *, "M2 [GeV] =       ", M2

      PRINT *
      PRINT *, "Theory constraints: 1 = all right, 0 = violated:"
      PRINT *, "Unitarity, |Re a_0| < 1/2:    ", UNIOK
      PRINT *, "Potential bounded from below: ", BFBOK
      PRINT *, "No deeper alternative minima: ", MINOK
      PRINT *, "All masses real and positive: ", POSMSQOK
      PRINT *
      PRINT *, "Indirect constraints:"
      PRINT *, "Bs->mumu/SM = ", RBSMM, " allowed: ", BSMMOK
      PRINT *, "S parameter = ", SPARAM, " allowed: ", SPAROK
      PRINT *, "BR(b -> s gamma) (from SuperIso v3.3) allowed",
     .     " (loose): ",  BSGAMLOOSEOK
      PRINT *, "BR(b -> s gamma) (from SuperIso v3.3) allowed",
     .     " (tight): ", BSGAMTIGHTOK
      PRINT *

      PRINT *, "Physical parameters:"
      PRINT *, "v_phi [GeV] = ", VPHI
      PRINT *, "v_chi [GeV] = ", VCHI
      PRINT *, "M_h [GeV] =   ", MHL
      PRINT *, "M_H [GeV] =   ", MHH
      PRINT *, "M_3 [GeV] =   ", MH3
      PRINT *, "M_5 [GeV] =   ", MH5
      PRINT *, "alpha =       ", ALPHA
      PRINT *, "sin(alpha) =  ", DSIN(ALPHA)
      PRINT *, "cos(alpha) =  ", DCOS(ALPHA)
      PRINT *, "sin(thetaH) = ", DSQRT(8.D0)*VCHI
     .     / DSQRT(VPHI**2+8.D0*VCHI**2)
      PRINT *, "cot(thetaH) = ", VPHI/2.D0/DSQRT(2.D0)/VCHI

      ELSE
C The input parameters were not ok!        
         PRINT *, "!!!! Inputs do not give an acceptable potential !!!!"
      ENDIF

      PRINT *

      RETURN
      END

C==================================================================
C     Subroutine to write out the kappa factors 
C==================================================================

      SUBROUTINE PRINT_HCOUPS
      IMPLICIT NONE
C Common blocks:
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
      DOUBLE PRECISION V, MZ, MW, MTPOLE, MBMB, MCMC, MS, MTAU, MMU,
     .     ALPHAEM,ALSMZ, VCB, VUS, VUB, GAMZ, GAMW
      COMMON/SM/V,MZ,MW,MTPOLE,MBMB,MCMC,MS,MTAU,MMU,
     .     ALPHAEM,ALSMZ,VCB,VUS,VUB,GAMZ,GAMW

      PRINT *
      PRINT *, "Coupling modification factors for h, mass = ", 
     .     MHL, " GeV:"
      PRINT *, "kappa_V =            ", KVL
      PRINT *, "kappa_f =            ", KFL
      PRINT *, "kappa_gamma =        ", KGAML
      PRINT *, "Delta kappa_gamma =  ", DKGAML
      IF (MHL.LT.MZ) THEN
         PRINT *, "(mh < mZ: no h -> Z gamma decay)"
      ELSE
         PRINT *, "kappa_Zgamma =       ", KZGAML
         PRINT *, "Delta kappa_Zgamma = ", DKZGAML
      ENDIF

      PRINT *
      PRINT *, "Coupling modification factors for H, mass = ", 
     .     MHH, " GeV:"
      PRINT *, "kappa_V =            ", KVH
      PRINT *, "kappa_f =            ", KFH
      PRINT *, "kappa_gamma =        ", KGAMH
      PRINT *, "Delta kappa_gamma =  ", DKGAMH
      IF (MHH.LT.MZ) THEN
         PRINT *, "(mH < mZ: no H -> Z gamma decay)"
      ELSE
         PRINT *, "kappa_Zgamma =       ", KZGAMH
         PRINT *, "Delta kappa_Zgamma = ", DKZGAMH
      ENDIF

! modified by Ameen on 2018-08-02
      PRINT *
      PRINT *, "Coupling modification factors for H^3_0, mass = ",
     .     MH3, " GeV:"
      PRINT *, "kappa_W =            ", KW3
      PRINT *, "kappa_Z =            ", KZ3
      PRINT *, "kappa_f =            ", KF3

      PRINT *
      PRINT *, "Coupling modification factors for H^5_0, mass = ",
     .     MH5, " GeV:"
      PRINT *, "kappa_W =            ", KW5
      PRINT *, "kappa_Z =            ", KZ5
      PRINT *, "kappa_f =            ", KF5
! end modified by Ameen

      RETURN
      END

C=====================================================================
C Subroutine to generate the MG5 param_card for use with the FeynRules
C UFO model: LO EFT version
C=====================================================================

      SUBROUTINE WRITE_PARAM_CARD_EFT_LO
      IMPLICIT NONE
C Common blocks:
      DOUBLE PRECISION MU2SQ,MU3SQ,LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,
     .     LAMBDA5,M1,M2
      COMMON/LPARAMS/MU2SQ,MU3SQ,LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,
     .     LAMBDA5,M1,M2
      DOUBLE PRECISION MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
      COMMON/PHYSPARAMS/MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
      DOUBLE PRECISION RxH5PWGA,RxH3PWGA,RxH3PWGATILDE,RxH5NGAGA,
     .     RxH5NZGA,RxH3NGAGA,RxH3NZGA,RxH3NGG,RxHHGAGA,RxHHZGA,RxHHGG,
     .     RxHLGAGA,RxHLZGA,RxHLGG,
     .     IxH5PWGA,IxH3PWGA,IxH3PWGATILDE
      COMMON/LOOPPARAMS/RxH5PWGA,RxH3PWGA,RxH3PWGATILDE,RxH5NGAGA,
     .     RxH5NZGA,RxH3NGAGA,RxH3NZGA,RxH3NGG,RxHHGAGA,RxHHZGA,RxHHGG,
     .     RxHLGAGA,RxHLZGA,RxHLGG,
     .     IxH5PWGA,IxH3PWGA,IxH3PWGATILDE
      DOUBLE PRECISION V, MZ, MW, MTPOLE, MBMB, MCMC, MS, MTAU, MMU,
     .     ALPHAEM,ALSMZ, VCB, VUS, VUB, GAMZ, GAMW
      COMMON/SM/V,MZ,MW,MTPOLE,MBMB,MCMC,MS,MTAU,MMU,
     .     ALPHAEM,ALSMZ,VCB,VUS,VUB,GAMZ,GAMW
      DOUBLE PRECISION GF, CABIBBO, MDPOLE, MUPOLE, MELEC
      COMMON/MG5SM/GF,CABIBBO,MDPOLE,MUPOLE,MELEC
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
C Local variables:
      DOUBLE PRECISION MCPOLE, MBPOLE, TANTHETAH
C Functions to be called:
      DOUBLE PRECISION POLEMQ


      MCPOLE = POLEMQ(MCMC,4)
      MBPOLE = POLEMQ(MBMB,5)
      TANTHETAH = DSQRT(8.D0)*VCHI/VPHI

      OPEN (UNIT = 90, FILE = 'param_card-EFTLO.dat')

      WRITE(90,1)
      WRITE(90,3) "         PARAM_CARD GENERATED BY GMCALC          "
      WRITE(90,1)
      WRITE(90,2)
      WRITE(90,3) "Width set on Auto will be computed following the "
      WRITE(90,3) "information present in the decay.py files of the "
      WRITE(90,3) "model. By default, this is only 1->2 decay modes."
      WRITE(90,2)
      WRITE(90,1)
      WRITE(90,*)

      WRITE(90,4)
      WRITE(90,5) "INFORMATION FOR CKMBLOCK        "
      WRITE(90,4)
      WRITE(90,6) "Block ckmblock                  "
      WRITE(90,7) 1, CABIBBO, "cabi          "
      WRITE(90,*)

c--ADD By Y. Wu @ 26 Feb 2018 

      WRITE(90,4)
      WRITE(90,5) "INFORMATION FOR NEUTRALFORMFACTORS   "
      WRITE(90,4)
      WRITE(90,6) "Block NEUTRALFORMFACTORS               "
      WRITE(90,77) 1, RxH5NGAGA, "H5NGAGA             "
      WRITE(90,77) 2, RxH5NZGA, "H5NZGA              "
      WRITE(90,77) 3, RxH3NGAGA, "H3NGAGA             "
      WRITE(90,77) 4, RxH3NZGA, "H3NZGA            "
      WRITE(90,77) 5, RxH3NGG, "H3NGG               "
      WRITE(90,77) 6, RxHHGAGA, "HHGAGA              "
      WRITE(90,77) 7, RxHHZGA, "HHZGA               "
      WRITE(90,77) 8, RxHHGG, "HHGG                "
      WRITE(90,77) 9, RxHLGAGA, "HLGAGA              "
      WRITE(90,77) 10, RxHLZGA, "HLZGA               "
      WRITE(90,77) 11, RxHLGG, "HLGG                "
      WRITE(90,*)

      WRITE(90,4)
      WRITE(90,5) "INFORMATION FOR FORMFACTORS     "
      WRITE(90,4)
      WRITE(90,6) "Block REXCHARGEDFORMFACTORS     "
      WRITE(90,77) 1, RxH5PWGA, "RxH5PWGA            "
      WRITE(90,77) 2, 0.D0, "RxH5PWGATILDE(DEBUG)"
      WRITE(90,77) 3, RxH3PWGA, "RxH3PWGA            "
      WRITE(90,77) 4, RxH3PWGATILDE, "RxH3PWGATILDE       "
      WRITE(90,*)

      WRITE(90,4)
      WRITE(90,5) "INFORMATION FOR FORMFACTORS     "
      WRITE(90,4)
      WRITE(90,6) "Block IMXCHARGEDFORMFACTORS     "
      WRITE(90,77) 1, IxH5PWGA, "IxH5PWGA            "
      WRITE(90,77) 2, 0.D0, "IxH5PWGATILDE(DEBUG)"
      WRITE(90,77) 3, IxH3PWGA, "IxH3PWGA            "
      WRITE(90,77) 4, IxH3PWGATILDE, "IxH3PWGATILDE       "
      WRITE(90,*)

c--End ADD

      WRITE(90,4)
      WRITE(90,5) "INFORMATION FOR MASS            "
      WRITE(90,4)
      WRITE(90,6) "Block mass                      "
      WRITE(90,7) 1, MDPOLE, "MD            "
      WRITE(90,7) 2, MUPOLE, "MU            "
      WRITE(90,7) 3, MS,     "MS            "
      WRITE(90,7) 4, MCPOLE, "MC            "
      WRITE(90,7) 5, MBPOLE, "MB            "
      WRITE(90,7) 6, MTPOLE, "MT            "
      WRITE(90,7) 11, MELEC, "Me            "
      WRITE(90,7) 13, MMU,   "MM            "
      WRITE(90,7) 15, MTAU,  "MTA           "
      WRITE(90,7) 23, MZ,    "MZ            "
      WRITE(90,7) 25, MHL,   "Mh            "
      WRITE(90,3) "Dependent parameters, given by model restrictions."
      WRITE(90,3) "MG5 ignores these values but they are important   "
      WRITE(90,3) "for interfacing the output of MG5 to external     "
      WRITE(90,3) "programs such as Pythia.                          "
      WRITE(90,7) 12, 0.D0, "ve            "
      WRITE(90,7) 14, 0.D0, "vm            "
      WRITE(90,7) 16, 0.D0, "vt            "
      WRITE(90,7) 21, 0.D0, "g             "
      WRITE(90,7) 22, 0.D0, "a             "
      WRITE(90,7) 24, MW,   "w+   (derived)"
      WRITE(90,7) 252, MHH, "H    (derived)"
      WRITE(90,7) 253, MH3, "H3p  (derived)"
      WRITE(90,7) 255, MH5, "H5pp (derived)"
      WRITE(90,7) 254, MH3, "H3z  (derived)"
      WRITE(90,7) 256, MH5, "H5p  (derived)"
      WRITE(90,7) 257, MH5, "H5z  (derived)"
      WRITE(90,*)

      WRITE(90,4)
      WRITE(90,5) "INFORMATION FOR POTENTIALPARAM  "
      WRITE(90,4)
      WRITE(90,6) "Block potentialparam            "
      WRITE(90,7) 1, LAMBDA2, "lam2          "
      WRITE(90,7) 2, LAMBDA3, "lam3          "
      WRITE(90,7) 3, LAMBDA4, "lam4          "
      WRITE(90,7) 4, LAMBDA5, "lam5          "
      WRITE(90,7) 5, M1,      "M1coeff       "
      WRITE(90,7) 6, M2,      "M2coeff       "
      WRITE(90,*)

      WRITE(90,4)
      WRITE(90,5) "INFORMATION FOR SMINPUTS        "
      WRITE(90,4)
      WRITE(90,6) "Block sminputs                  "
      WRITE(90,7) 1, 1.D0/ALPHAEM, "aEWM1         "
      WRITE(90,7) 2, GF,           "Gf            "
      WRITE(90,7) 3, ALSMZ,        "aS            "
      WRITE(90,*)

      WRITE(90,4)
      WRITE(90,5) "INFORMATION FOR VEV             "
      WRITE(90,4)
      WRITE(90,6) "Block vev                       "
      WRITE(90,7) 1, TANTHETAH, "tanth         "
      WRITE(90,*)

      WRITE(90,4)
      WRITE(90,5) "INFORMATION FOR YUKAWA          "
      WRITE(90,4)
      WRITE(90,6) "Block yukawa                    "
      WRITE(90,7) 1, MDPOLE, "ymdo          "
      WRITE(90,7) 2, MUPOLE, "ymup          "
      WRITE(90,7) 3, MS,     "yms           "
      WRITE(90,7) 4, MCPOLE, "ymc           "
      WRITE(90,7) 5, MBPOLE, "ymb           "
      WRITE(90,7) 6, MTPOLE, "ymt           "
      WRITE(90,7) 11, MELEC, "yme           "
      WRITE(90,7) 13, MMU,   "ymm           "
      WRITE(90,7) 15, MTAU,  "ymtau         "
      WRITE(90,*)

CC Adding Whole decay table for scalars By Yongcheng Wu 2018 Feb

      WRITE(90,4)
      WRITE(90,5) "INFORMATION FOR DECAY           "
      WRITE(90,4)
C       WRITE(90,15) 6, 1.D0,   'WT   '
      WRITE(90,150)
      WRITE(90,15) 23, 2.4952D0,  'WZ   '
      WRITE(90,15) 24, 2.085D0,  'WW   '
C       WRITE(90,16) 25,  'Wh   '
C       WRITE(90,16) 252, 'Wh__2'
C       WRITE(90,16) 253, 'WH3p '
C       WRITE(90,16) 254, 'WH3z '
C       WRITE(90,16) 255, 'WH5pp'
C       WRITE(90,16) 256, 'WH5p '
C       WRITE(90,16) 257, 'WH5z '
      WRITE(90,150)
      WRITE(90,15) 6, TOPWDTH, 'WT   '
      WRITE(90,151)
      WRITE(90,152) TOPBRW,2,24,5,'W+ b      '
      WRITE(90,152) TOPBRH3P,2,253,5,'H3+ b     '
      WRITE(90,153)
      WRITE(90,150)
      WRITE(90,15) 25, HLWDTH, 'Wh   '
      WRITE(90,151)
      WRITE(90,152) HLBRB,2,-5,5,'b~ b      '
      WRITE(90,152) HLBRTA,2,-15,15,'ta~ ta    '
      WRITE(90,152) HLBRMU,2,-13,13,'mu~ mu    '
      WRITE(90,152) HLBRS,2,-3,3,'s~ s      '
      WRITE(90,152) HLBRC,2,-4,4,'c~ c      '
      WRITE(90,152) HLBRT,2,-6,6,'t~ t      '
      WRITE(90,152) HLBRG,2,21,21,'g g       '
      WRITE(90,152) HLBRGA,2,22,22,'a a       '
      WRITE(90,152) HLBRZGA,2,23,22,'Z a      '
      WRITE(90,152) HLBRW,2,-24,24,'W+ W-     '
      WRITE(90,152) HLBRZ,2,23,23,'Z Z       '
      WRITE(90,152) HLBRWH3P,2,-24,253,'W- H3+ +hc'
      WRITE(90,152) HLBRZH3N,2,23,254,'Z H30     '
      WRITE(90,152) HLBRH3N,2,254,254,'H30 H30   '
      WRITE(90,152) HLBRH3P,2,-253,253,'H3m H3p   '
      WRITE(90,152) HLBRH5N,2,257,257,'H50 H50   '
      WRITE(90,152) HLBRH5P,2,-256,256,'H5m H5p   '
      WRITE(90,152) HLBRH5PP,2,-255,255,'H5mm H5pp '
      WRITE(90,153)

      WRITE(90,150)
      WRITE(90,15) 252, HHWDTH, 'WH   '
      WRITE(90,151)
      WRITE(90,152) HHBRB,2,-5,5,'b~ b      '
      WRITE(90,152) HHBRTA,2,-15,15,'ta~ ta    '
      WRITE(90,152) HHBRMU,2,-13,13,'mu~ mu    '
      WRITE(90,152) HHBRS,2,-3,3,'s~ s      '
      WRITE(90,152) HHBRC,2,-4,4,'c~ c      '
      WRITE(90,152) HHBRT,2,-6,6,'t~ t      '
      WRITE(90,152) HHBRG,2,21,21,'g g       '
      WRITE(90,152) HHBRGA,2,22,22,'a a       '
      WRITE(90,152) HHBRZGA,2,23,22,'Z a      '
      WRITE(90,152) HHBRW,2,-24,24,'W+ W-     '
      WRITE(90,152) HHBRZ,2,23,23,'Z Z       '
      WRITE(90,152) HHBRWH3P,2,-24,253,'W- H3+ +hc'
      WRITE(90,152) HHBRZH3N,2,23,254,'Z H30     '
      WRITE(90,152) HHBRHL,2,25,25,'h h       '
      WRITE(90,152) HHBRH3N,2,254,254,'H30 H30   '
      WRITE(90,152) HHBRH3P,2,-253,253,'H3m H3p   '
      WRITE(90,152) HHBRH5N,2,257,257,'H50 H50   '
      WRITE(90,152) HHBRH5P,2,-256,256,'H5m H5p   '
      WRITE(90,152) HHBRH5PP,2,-255,255,'H5mm H5pp '
      WRITE(90,153)

      WRITE(90,150)
      WRITE(90,15) 254, H3NWDTH, 'WH3N '
      WRITE(90,151)
      WRITE(90,152) H3NBRB,2,-5,5,'b~ b      '
      WRITE(90,152) H3NBRTA,2,-15,15,'ta~ ta    '
      WRITE(90,152) H3NBRMU,2,-13,13,'mu~ mu    '
      WRITE(90,152) H3NBRS,2,-3,3,'s~ s      '
      WRITE(90,152) H3NBRC,2,-4,4,'c~ c      '
      WRITE(90,152) H3NBRT,2,-6,6,'t~ t      '
      WRITE(90,152) H3NBRG,2,21,21,'g g       '
      WRITE(90,152) H3NBRGA,2,22,22,'a a       '
      WRITE(90,152) H3NBRZGA,2,23,22,'Z a      '
      WRITE(90,152) H3NBRZHL,2,23,25,'Z h      '
      WRITE(90,152) H3NBRZHH,2,23,252,'Z H      '
      WRITE(90,152) H3NBRZH5N,2,23,257,'Z H50    '
      WRITE(90,152) H3NBRWH5P,2,-24,256,'W- H5p +hc'
      WRITE(90,153)

      WRITE(90,150)
      WRITE(90,15) 253,H3PWDTH, 'WH3P '
      WRITE(90,151)
      WRITE(90,152) H3PBRBC,2,-5,4,'b~ c      '
      WRITE(90,152) H3PBRTA,2,-15,16,'ta~ nu    '
      WRITE(90,152) H3PBRMU,2,-13,14,'mu~ nu    '
      WRITE(90,152) H3PBRSU,2,-3,2,'s~ u      '
      WRITE(90,152) H3PBRCS,2,-3,4,'s~ c      '
      WRITE(90,152) H3PBRTB,2,-5,6,'b~ t      '
      WRITE(90,152) H3PBRBU,2,-5,2,'b~ u      '
      WRITE(90,152) H3PBRWHL,2,24,25,'W+ h      '
      WRITE(90,152) H3PBRWHH,2,24,252,'W+ H      '
      WRITE(90,152) H3PBRZH5P,2,23,256,'Z H5p     '
      WRITE(90,152) H3PBRWH5N,2,24,257,'W+ H50    '
      WRITE(90,152) H3PBRWH5PP,2,-24,255,'W- H5pp   '
      WRITE(90,152) H3PBRWGA,2,24,22,'W+ a      '
      WRITE(90,153)

      WRITE(90,150)
      WRITE(90,15) 257,H5NWDTH, 'WH50 '
      WRITE(90,151)
      WRITE(90,152) H5NBRGA,2,22,22,'a a       '
      WRITE(90,152) H5NBRZGA,2,23,22,'Z a       '
      WRITE(90,152) H5NBRW,2,-24,24,'W+ W-     '
      WRITE(90,152) H5NBRZ,2,23,23,'Z Z       '
      WRITE(90,152) H5NBRZH3N,2,23,254,'Z H30     '
      WRITE(90,152) H5NBRWH3P,2,-24,253,'W- H3p +hc'
      WRITE(90,152) H5NBRH3N,2,254,254,'H30 H30   '
      WRITE(90,152) H5NBRH3P,2,-253,253,'H3m H3p   '
      WRITE(90,153)

      WRITE(90,150)
      WRITE(90,15) 256,H5PWDTH, 'WH5p '
      WRITE(90,152) H5PBRWZ,2,24,23,'W+ Z      '
      WRITE(90,152) H5PBRZH3P,2,23,253,'Z H3p     '
      WRITE(90,152) H5PBRWH3N,2,24,254,'W+ H30    '
      WRITE(90,152) H5PBRH3PN,2,253,254,'H3p H30   '
      WRITE(90,152) H5PBRWGA,2,24,22,'W+ a      '
      WRITE(90,153)

      WRITE(90,150)
      WRITE(90,15) 255,H5PPWDTH, 'WH5pp'
      WRITE(90,151)
      WRITE(90,152) H5PPBRWW,2,24,24,'W+ W+     '
      WRITE(90,152) H5PPBRWH3,2,24,253,'W+ H3p    '
      WRITE(90,152) H5PPBRH3P,2,253,253,'H3p H3p   '
      WRITE(90,153)


      WRITE(90,3) "Dependent parameters, given by model restrictions."
      WRITE(90,3) "MG5 ignores these values but they are important   "
      WRITE(90,3) "for interfacing the output of MG5 to external     "
      WRITE(90,3) "programs such as Pythia.                          "
      WRITE(90,15) 1, 0.D0,  'd    '
      WRITE(90,15) 2, 0.D0,  'u    '
      WRITE(90,15) 3, 0.D0,  's    '
      WRITE(90,15) 4, 0.D0,  'c    '
      WRITE(90,15) 5, 0.D0,  'b    '
      WRITE(90,15) 11, 0.D0, 'e-   '
      WRITE(90,15) 12, 0.D0, 've   '
      WRITE(90,15) 13, 0.D0, 'mu-  '
      WRITE(90,15) 14, 0.D0, 'vm   '
      WRITE(90,15) 15, 0.D0, 'ta-  '
      WRITE(90,15) 16, 0.D0, 'vt   '
      WRITE(90,15) 21, 0.D0, 'g    '
      WRITE(90,15) 22, 0.D0, 'a    '

      WRITE(90,8)
      WRITE(90,9) "QUANTUM NUMBERS OF NEW STATE(S) (NON SM PDG CODE) "
      WRITE(90,8)
      WRITE(90,*)
      WRITE(90,10) 252, "H    "
      WRITE(90,11) 0
      WRITE(90,12) 1
      WRITE(90,13) 1
      WRITE(90,14) 0
      WRITE(90,10) 253, "H3p  "
      WRITE(90,11) 3
      WRITE(90,12) 1
      WRITE(90,13) 1
      WRITE(90,14) 1
      WRITE(90,10) 254, "H3z  "
      WRITE(90,11) 0
      WRITE(90,12) 1
      WRITE(90,13) 1
      WRITE(90,14) 0
      WRITE(90,10) 255, "H5++ "
      WRITE(90,11) 6
      WRITE(90,12) 1
      WRITE(90,13) 1
      WRITE(90,14) 1
      WRITE(90,10) 256, "H5p  "
      WRITE(90,11) 3
      WRITE(90,12) 1
      WRITE(90,13) 1
      WRITE(90,14) 1
      WRITE(90,10) 257, "H5z  "
      WRITE(90,11) 0
      WRITE(90,12) 1
      WRITE(90,13) 1
      WRITE(90,14) 0

 1    FORMAT(70('#'))
 2    FORMAT(2('#'),66X,2('#'))
 3    FORMAT(2('#'),8X,A50,8X,2('#'))
 4    FORMAT(35('#'))
 5    FORMAT(2('#'),1X,A32)
 6    FORMAT(A32)
 7    FORMAT(I5,1X,e12.6,1X,'#',1X,A14)
 77   FORMAT(I5,1X,e12.6,1X,'#',1X,A20)
 8    FORMAT('#',59('='))
 9    FORMAT('#',1X,A50)
 10   FORMAT('Block QNUMBERS ',I3,2X,'#',1X,A5) 
 11   FORMAT(8X,'1',1X,I1,2X,'# 3 times electric charge')
 12   FORMAT(8X,'2',1X,I1,2X,'# number of spin states (2S+1)')
 13   FORMAT(8X,'3',1X,I1,2X,'# colour rep (1: singlet, 3: triplet,',
     .     ' 8: octet)')
 14   FORMAT(8X,'4',1X,I1,2X,'# Particle/Antiparticle distinction',
     .     ' (0=own anti)')
 150  FORMAT('#     PDG',1X,'    Width    ')
 15   FORMAT('DECAY ',I3,1X,e12.6,' # ',A5)
 151  FORMAT('#     ',1X,'  BR  ',1X,'  NDA  ',1X,'  ID1  ',1X,'  ID2')
 152  FORMAT(1X,e12.6,1X,I3,1X,I6,1X,I6,1x,' # ',A10)
 153  FORMAT('#')
 16   FORMAT('DECAY ',I3,1X,'Auto ',' # ',A5)
      
      CLOSE (UNIT = 90)

      RETURN
      END


C=====================================================================
C Subroutine to generate the MG5 param_card for use with the FeynRules
C UFO model: LO version
C=====================================================================

      SUBROUTINE WRITE_PARAM_CARD_LO
      IMPLICIT NONE
C Common blocks:
      DOUBLE PRECISION MU2SQ,MU3SQ,LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,
     .     LAMBDA5,M1,M2
      COMMON/LPARAMS/MU2SQ,MU3SQ,LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,
     .     LAMBDA5,M1,M2
      DOUBLE PRECISION MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
      COMMON/PHYSPARAMS/MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
      DOUBLE PRECISION V, MZ, MW, MTPOLE, MBMB, MCMC, MS, MTAU, MMU,
     .     ALPHAEM,ALSMZ, VCB, VUS, VUB, GAMZ, GAMW
      COMMON/SM/V,MZ,MW,MTPOLE,MBMB,MCMC,MS,MTAU,MMU,
     .     ALPHAEM,ALSMZ,VCB,VUS,VUB,GAMZ,GAMW
      DOUBLE PRECISION GF, CABIBBO, MDPOLE, MUPOLE, MELEC
      COMMON/MG5SM/GF,CABIBBO,MDPOLE,MUPOLE,MELEC
C Local variables:
      DOUBLE PRECISION MCPOLE, MBPOLE, TANTHETAH
C Functions to be called:
      DOUBLE PRECISION POLEMQ

      MCPOLE = POLEMQ(MCMC,4)
      MBPOLE = POLEMQ(MBMB,5)
      TANTHETAH = DSQRT(8.D0)*VCHI/VPHI

      OPEN (UNIT = 90, FILE = 'param_card-LO.dat')

      WRITE(90,1)
      WRITE(90,3) "         PARAM_CARD GENERATED BY GMCALC          "
      WRITE(90,1)
      WRITE(90,2)
      WRITE(90,3) "Width set on Auto will be computed following the "
      WRITE(90,3) "information present in the decay.py files of the "
      WRITE(90,3) "model. By default, this is only 1->2 decay modes."
      WRITE(90,2)
      WRITE(90,1)
      WRITE(90,*)

      WRITE(90,4)
      WRITE(90,5) "INFORMATION FOR CKMBLOCK        "
      WRITE(90,4)
      WRITE(90,6) "Block ckmblock                  "
      WRITE(90,7) 1, CABIBBO, "cabi          "
      WRITE(90,*)

      WRITE(90,4)
      WRITE(90,5) "INFORMATION FOR MASS            "
      WRITE(90,4)
      WRITE(90,6) "Block mass                      "
      WRITE(90,7) 1, MDPOLE, "MD            "
      WRITE(90,7) 2, MUPOLE, "MU            "
      WRITE(90,7) 3, MS,     "MS            "
      WRITE(90,7) 4, MCPOLE, "MC            "
      WRITE(90,7) 5, MBPOLE, "MB            "
      WRITE(90,7) 6, MTPOLE, "MT            "
      WRITE(90,7) 11, MELEC, "Me            "
      WRITE(90,7) 13, MMU,   "MM            "
      WRITE(90,7) 15, MTAU,  "MTA           "
      WRITE(90,7) 23, MZ,    "MZ            "
      WRITE(90,7) 25, MHL,   "Mh            "
      WRITE(90,3) "Dependent parameters, given by model restrictions."
      WRITE(90,3) "MG5 ignores these values but they are important   "
      WRITE(90,3) "for interfacing the output of MG5 to external     "
      WRITE(90,3) "programs such as Pythia.                          "
      WRITE(90,7) 12, 0.D0, "ve            "
      WRITE(90,7) 14, 0.D0, "vm            "
      WRITE(90,7) 16, 0.D0, "vt            "
      WRITE(90,7) 21, 0.D0, "g             "
      WRITE(90,7) 22, 0.D0, "a             "
      WRITE(90,7) 24, MW,   "w+   (derived)"
      WRITE(90,7) 252, MHH, "H    (derived)"
      WRITE(90,7) 253, MH3, "H3p  (derived)"
      WRITE(90,7) 255, MH5, "H5pp (derived)"
      WRITE(90,7) 254, MH3, "H3z  (derived)"
      WRITE(90,7) 256, MH5, "H5p  (derived)"
      WRITE(90,7) 257, MH5, "H5z  (derived)"
      WRITE(90,*)

      WRITE(90,4)
      WRITE(90,5) "INFORMATION FOR POTENTIALPARAM  "
      WRITE(90,4)
      WRITE(90,6) "Block potentialparam            "
      WRITE(90,7) 1, LAMBDA2, "lam2          "
      WRITE(90,7) 2, LAMBDA3, "lam3          "
      WRITE(90,7) 3, LAMBDA4, "lam4          "
      WRITE(90,7) 4, LAMBDA5, "lam5          "
      WRITE(90,7) 5, M1,      "M1coeff       "
      WRITE(90,7) 6, M2,      "M2coeff       "
      WRITE(90,*)

      WRITE(90,4)
      WRITE(90,5) "INFORMATION FOR SMINPUTS        "
      WRITE(90,4)
      WRITE(90,6) "Block sminputs                  "
      WRITE(90,7) 1, 1.D0/ALPHAEM, "aEWM1         "
      WRITE(90,7) 2, GF,           "Gf            "
      WRITE(90,7) 3, ALSMZ,        "aS            "
      WRITE(90,*)

      WRITE(90,4)
      WRITE(90,5) "INFORMATION FOR VEV             "
      WRITE(90,4)
      WRITE(90,6) "Block vev                       "
      WRITE(90,7) 1, TANTHETAH, "tanth         "
      WRITE(90,*)

      WRITE(90,4)
      WRITE(90,5) "INFORMATION FOR YUKAWA          "
      WRITE(90,4)
      WRITE(90,6) "Block yukawa                    "
      WRITE(90,7) 1, MDPOLE, "ymdo          "
      WRITE(90,7) 2, MUPOLE, "ymup          "
      WRITE(90,7) 3, MS,     "yms           "
      WRITE(90,7) 4, MCPOLE, "ymc           "
      WRITE(90,7) 5, MBPOLE, "ymb           "
      WRITE(90,7) 6, MTPOLE, "ymt           "
      WRITE(90,7) 11, MELEC, "yme           "
      WRITE(90,7) 13, MMU,   "ymm           "
      WRITE(90,7) 15, MTAU,  "ymtau         "
      WRITE(90,*)

      WRITE(90,4)
      WRITE(90,5) "INFORMATION FOR DECAY           "
      WRITE(90,4)
      WRITE(90,15) 6, 1.D0,   'WT   '
      WRITE(90,15) 23, 1.D0,  'WZ   '
      WRITE(90,15) 24, 1.D0,  'WW   '
      WRITE(90,16) 25,  'Wh   '
      WRITE(90,16) 252, 'Wh__2'
      WRITE(90,16) 253, 'WH3p '
      WRITE(90,16) 254, 'WH3z '
      WRITE(90,16) 255, 'WH5pp'
      WRITE(90,16) 256, 'WH5p '
      WRITE(90,16) 257, 'WH5z '
      WRITE(90,3) "Dependent parameters, given by model restrictions."
      WRITE(90,3) "MG5 ignores these values but they are important   "
      WRITE(90,3) "for interfacing the output of MG5 to external     "
      WRITE(90,3) "programs such as Pythia.                          "
      WRITE(90,15) 1, 0.D0,  'd    '
      WRITE(90,15) 2, 0.D0,  'u    '
      WRITE(90,15) 3, 0.D0,  's    '
      WRITE(90,15) 4, 0.D0,  'c    '
      WRITE(90,15) 5, 0.D0,  'b    '
      WRITE(90,15) 11, 0.D0, 'e-   '
      WRITE(90,15) 12, 0.D0, 've   '
      WRITE(90,15) 13, 0.D0, 'mu-  '
      WRITE(90,15) 14, 0.D0, 'vm   '
      WRITE(90,15) 15, 0.D0, 'ta-  '
      WRITE(90,15) 16, 0.D0, 'vt   '
      WRITE(90,15) 21, 0.D0, 'g    '
      WRITE(90,15) 22, 0.D0, 'a    '

      WRITE(90,8)
      WRITE(90,9) "QUANTUM NUMBERS OF NEW STATE(S) (NON SM PDG CODE) "
      WRITE(90,8)
      WRITE(90,*)
      WRITE(90,10) 252, "H    "
      WRITE(90,11) 0
      WRITE(90,12) 1
      WRITE(90,13) 1
      WRITE(90,14) 0
      WRITE(90,10) 253, "H3p  "
      WRITE(90,11) 3
      WRITE(90,12) 1
      WRITE(90,13) 1
      WRITE(90,14) 1
      WRITE(90,10) 254, "H3z  "
      WRITE(90,11) 0
      WRITE(90,12) 1
      WRITE(90,13) 1
      WRITE(90,14) 0
      WRITE(90,10) 255, "H5++ "
      WRITE(90,11) 6
      WRITE(90,12) 1
      WRITE(90,13) 1
      WRITE(90,14) 1
      WRITE(90,10) 256, "H5p  "
      WRITE(90,11) 3
      WRITE(90,12) 1
      WRITE(90,13) 1
      WRITE(90,14) 1
      WRITE(90,10) 257, "H5z  "
      WRITE(90,11) 0
      WRITE(90,12) 1
      WRITE(90,13) 1
      WRITE(90,14) 0

 1    FORMAT(70('#'))
 2    FORMAT(2('#'),66X,2('#'))
 3    FORMAT(2('#'),8X,A50,8X,2('#'))
 4    FORMAT(35('#'))
 5    FORMAT(2('#'),1X,A32)
 6    FORMAT(A32)
 7    FORMAT(I5,1X,e12.6,1X,'#',1X,A14)
 8    FORMAT('#',59('='))
 9    FORMAT('#',1X,A50)
 10   FORMAT('Block QNUMBERS ',I3,2X,'#',1X,A5) 
 11   FORMAT(8X,'1',1X,I1,2X,'# 3 times electric charge')
 12   FORMAT(8X,'2',1X,I1,2X,'# number of spin states (2S+1)')
 13   FORMAT(8X,'3',1X,I1,2X,'# colour rep (1: singlet, 3: triplet,',
     .     ' 8: octet)')
 14   FORMAT(8X,'4',1X,I1,2X,'# Particle/Antiparticle distinction',
     .     ' (0=own anti)')
 15   FORMAT('DECAY ',I3,1X,e12.6,' # ',A5)
 16   FORMAT('DECAY ',I3,1X,'Auto ',' # ',A5)
      
      CLOSE (UNIT = 90)

      RETURN
      END


C=====================================================================
C Subroutine to generate the MG5 param_card for use with the FeynRules
C UFO model: NLO version --  A Peterson
C=====================================================================

      SUBROUTINE WRITE_PARAM_CARD_NLO
      IMPLICIT NONE
C Common blocks:
      DOUBLE PRECISION MU2SQ,MU3SQ,LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,
     .     LAMBDA5,M1,M2
      COMMON/LPARAMS/MU2SQ,MU3SQ,LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,
     .     LAMBDA5,M1,M2
      DOUBLE PRECISION MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
      COMMON/PHYSPARAMS/MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
      DOUBLE PRECISION V, MZ, MW, MTPOLE, MBMB, MCMC, MS, MTAU, MMU,
     .     ALPHAEM,ALSMZ, VCB, VUS, VUB, GAMZ, GAMW
      COMMON/SM/V,MZ,MW,MTPOLE,MBMB,MCMC,MS,MTAU,MMU,
     .     ALPHAEM,ALSMZ,VCB,VUS,VUB,GAMZ,GAMW
      DOUBLE PRECISION GF, CABIBBO, MDPOLE, MUPOLE, MELEC
      COMMON/MG5SM/GF,CABIBBO,MDPOLE,MUPOLE,MELEC
C Local variables:
      DOUBLE PRECISION MCPOLE, MBPOLE, TANTHETAH
C Functions to be called:
      DOUBLE PRECISION POLEMQ

      MCPOLE = POLEMQ(MCMC,4)
      MBPOLE = POLEMQ(MBMB,5)
      TANTHETAH = DSQRT(8.D0)*VCHI/VPHI

      OPEN (UNIT = 90, FILE = 'param_card-NLO.dat')

      WRITE(90,1)
      WRITE(90,3) "         PARAM_CARD GENERATED BY GMCALC          "
      WRITE(90,1)
      WRITE(90,2)
      WRITE(90,3) "Width set on Auto will be computed following the "
      WRITE(90,3) "information present in the decay.py files of the "
      WRITE(90,3) "model. By default, this is only 1->2 decay modes."
      WRITE(90,2)
      WRITE(90,1)
      WRITE(90,*)

      WRITE(90,4)
      WRITE(90,5) "INFORMATION FOR LOOP            "
      WRITE(90,4)
      WRITE(90,6) "Block loop                      "
      WRITE(90,7) 1, MW, "MU_R          "
      WRITE(90,*)

      WRITE(90,4)
      WRITE(90,5) "INFORMATION FOR MASS            "
      WRITE(90,4)
      WRITE(90,6) "Block mass                      "
      WRITE(90,7) 6, MTPOLE, "MT            "
      WRITE(90,7) 13, MMU,   "MM            "
      WRITE(90,7) 15, MTAU,  "MTA           "
      WRITE(90,7) 23, MZ,    "MZ            "
      WRITE(90,7) 24, MW,    "MW            "
      WRITE(90,7) 25, MHL,   "Mh            "
      WRITE(90,7) 252, MHH, "H              "
      WRITE(90,7) 253, MH3, "H3p            "
      WRITE(90,7) 255, MH5, "H5pp           "
      WRITE(90,7) 254, MH3, "H3z            "
      WRITE(90,7) 256, MH5, "H5p            "
      WRITE(90,7) 257, MH5, "H5z            "
      WRITE(90,3) "Dependent parameters, given by model restrictions."
      WRITE(90,3) "MG5 ignores these values but they are important   "
      WRITE(90,3) "for interfacing the output of MG5 to external     "
      WRITE(90,3) "programs such as Pythia.                          "
      WRITE(90,7) 1, 0.D0, "MD            "
      WRITE(90,7) 2, 0.D0, "MU            "
      WRITE(90,7) 3, 0.D0, "MS            "
      WRITE(90,7) 4, 0.D0, "MC            "
      WRITE(90,7) 5, 0.D0, "MB            "
      WRITE(90,7) 11, 0.D0, "Me            "
      WRITE(90,7) 12, 0.D0, "ve            "
      WRITE(90,7) 14, 0.D0, "vm            "
      WRITE(90,7) 16, 0.D0, "vt            "
      WRITE(90,7) 21, 0.D0, "g             "
      WRITE(90,7) 22, 0.D0, "a             "
      WRITE(90,*)

      WRITE(90,4)
      WRITE(90,5) "INFORMATION FOR POTENTIALPARAM  "
      WRITE(90,4)
      WRITE(90,6) "Block potentialparam            "
      WRITE(90,7) 1, LAMBDA2, "lam2          "
      WRITE(90,7) 2, LAMBDA3, "lam3          "
      WRITE(90,7) 3, LAMBDA4, "lam4          "
      WRITE(90,7) 4, LAMBDA5, "lam5          "
      WRITE(90,7) 5, M1,      "M1coeff       "
      WRITE(90,7) 6, M2,      "M2coeff       "
      WRITE(90,*)

      WRITE(90,4)
      WRITE(90,5) "INFORMATION FOR SMINPUTS        "
      WRITE(90,4)
      WRITE(90,6) "Block sminputs                  "
      WRITE(90,7) 1, 1.D0/ALPHAEM, "aEWM1         "
      WRITE(90,7) 2, GF,           "Gf            "
      WRITE(90,7) 3, ALSMZ,        "aS            "
      WRITE(90,*)

      WRITE(90,4)
      WRITE(90,5) "INFORMATION FOR VEV             "
      WRITE(90,4)
      WRITE(90,6) "Block vev                       "
      WRITE(90,7) 1, TANTHETAH, "tanth         "
      WRITE(90,*)

      WRITE(90,4)
      WRITE(90,5) "INFORMATION FOR YUKAWA          "
      WRITE(90,4)
      WRITE(90,6) "Block yukawa                    "
      WRITE(90,7) 1, MDPOLE, "ymdo          "
      WRITE(90,7) 2, MUPOLE, "ymup          "
      WRITE(90,7) 3, MS,     "yms           "
      WRITE(90,7) 4, MCPOLE, "ymc           "
      WRITE(90,7) 5, MBPOLE, "ymb           "
      WRITE(90,7) 6, MTPOLE, "ymt           "
      WRITE(90,7) 11, MELEC, "yme           "
      WRITE(90,7) 13, MMU,   "ymm           "
      WRITE(90,7) 15, MTAU,  "ymtau         "
      WRITE(90,*)

      WRITE(90,4)
      WRITE(90,5) "INFORMATION FOR DECAY           "
      WRITE(90,4)
      WRITE(90,15) 6, 1.508336D0,   'WT   '
      WRITE(90,15) 23, 2.495200D0,  'WZ   '
      WRITE(90,15) 24, 2.085000D0,  'WW   '
      WRITE(90,16) 25,  'Wh   '
      WRITE(90,16) 252, 'Wh__2'
      WRITE(90,16) 253, 'WH3p '
      WRITE(90,16) 254, 'WH3z '
      WRITE(90,16) 255, 'WH5pp'
      WRITE(90,16) 256, 'WH5p '
      WRITE(90,16) 257, 'WH5z '
      WRITE(90,3) "Dependent parameters, given by model restrictions."
      WRITE(90,3) "MG5 ignores these values but they are important   "
      WRITE(90,3) "for interfacing the output of MG5 to external     "
      WRITE(90,3) "programs such as Pythia.                          "
      WRITE(90,15) 1, 0.D0,  'd    '
      WRITE(90,15) 2, 0.D0,  'u    '
      WRITE(90,15) 3, 0.D0,  's    '
      WRITE(90,15) 4, 0.D0,  'c    '
      WRITE(90,15) 5, 0.D0,  'b    '
      WRITE(90,15) 11, 0.D0, 'e-   '
      WRITE(90,15) 12, 0.D0, 've   '
      WRITE(90,15) 13, 0.D0, 'mu-  '
      WRITE(90,15) 14, 0.D0, 'vm   '
      WRITE(90,15) 15, 0.D0, 'ta-  '
      WRITE(90,15) 16, 0.D0, 'vt   '
      WRITE(90,15) 21, 0.D0, 'g    '
      WRITE(90,15) 22, 0.D0, 'a    '

      WRITE(90,8)
      WRITE(90,9) "QUANTUM NUMBERS OF NEW STATE(S) (NON SM PDG CODE) "
      WRITE(90,8)
      WRITE(90,*)
      WRITE(90,10) 252, "H    "
      WRITE(90,11) 0
      WRITE(90,12) 1
      WRITE(90,13) 1
      WRITE(90,14) 0
      WRITE(90,10) 253, "H3p  "
      WRITE(90,11) 3
      WRITE(90,12) 1
      WRITE(90,13) 1
      WRITE(90,14) 1
      WRITE(90,10) 254, "H3z  "
      WRITE(90,11) 0
      WRITE(90,12) 1
      WRITE(90,13) 1
      WRITE(90,14) 0
      WRITE(90,10) 255, "H5++ "
      WRITE(90,11) 6
      WRITE(90,12) 1
      WRITE(90,13) 1
      WRITE(90,14) 1
      WRITE(90,10) 256, "H5p  "
      WRITE(90,11) 3
      WRITE(90,12) 1
      WRITE(90,13) 1
      WRITE(90,14) 1
      WRITE(90,10) 257, "H5z  "
      WRITE(90,11) 0
      WRITE(90,12) 1
      WRITE(90,13) 1
      WRITE(90,14) 0

 1    FORMAT(70('#'))
 2    FORMAT(2('#'),66X,2('#'))
 3    FORMAT(2('#'),8X,A50,8X,2('#'))
 4    FORMAT(35('#'))
 5    FORMAT(2('#'),1X,A32)
 6    FORMAT(A32)
 7    FORMAT(I5,1X,e12.6,1X,'#',1X,A14)
 8    FORMAT('#',59('='))
 9    FORMAT('#',1X,A50)
 10   FORMAT('Block QNUMBERS ',I3,2X,'#',1X,A5) 
 11   FORMAT(8X,'1',1X,I1,2X,'# 3 times electric charge')
 12   FORMAT(8X,'2',1X,I1,2X,'# number of spin states (2S+1)')
 13   FORMAT(8X,'3',1X,I1,2X,'# colour rep (1: singlet, 3: triplet,',
     .     ' 8: octet)')
 14   FORMAT(8X,'4',1X,I1,2X,'# Particle/Antiparticle distinction',
     .     ' (0=own anti)')
 15   FORMAT('DECAY ',I3,1X,e12.6,' # ',A5)
 16   FORMAT('DECAY ',I3,1X,'Auto ',' # ',A5)
      
      CLOSE (UNIT = 90)

      RETURN
      END


