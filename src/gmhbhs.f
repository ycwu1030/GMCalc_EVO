C============================================================================
C     Subroutines to call HiggsBounds 5 and HiggsSignals 2
C============================================================================
C Code from Ameen Ismail, 2018 / updated by HEL 2022

      SUBROUTINE CALLHBHS
C Calls HiggsBounds 5 and HiggsSignals 2.
C If you get an error message complaining about the USE statement,
C make sure your compiler supports IEEE intrinsic modules.
      USE, INTRINSIC :: ieee_arithmetic
      USE usefulbits_hs, only : HSres
      IMPLICIT NONE
C Common blocks for masses and couplings:
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
C Common blocks new to this subroutine:
      integer :: nHzero, nHplus
      integer :: HBresult(7), chan(7), ncombined(7)
      double precision :: obsratio(7)
      COMMON/HB5VARS/obsratio,nHzero,nHplus,HBresult,chan,ncombined
      integer :: nobs
      double precision :: Chisq_mu, Chisq_mh, Chisq, Pvalue
      COMMON/HS2VARS/Chisq_mu,Chisq_mh,Chisq,Pvalue,nobs
      double precision :: Chisq_mu_LHCRun1, Chisq_mh_LHCRun1, 
     .     Chisq_LHCRun1, Pvalue_LHCRun1
      double precision :: Chisq_STXS_rates, Chisq_STXS_mh, Chisq_STXS, 
     .     Pvalue_STXS
      integer ::  nobs_LHCRun1,  nobs_STXS
C Local HiggsBounds variables:
      double precision :: MHNeut(4),GammaTotal(4),ghjss_s(4),
     .     ghjss_p(4),ghjcc_s(4),ghjcc_p(4),ghjbb_s(4),ghjbb_p(4),
     .     ghjtt_s(4),ghjtt_p(4),ghjmumu_s(4),ghjmumu_p(4),
     .     ghjtautau_s(4),ghjtautau_p(4),ghjWW(4),ghjZZ(4),
     .     ghjZga(4),ghjgaga(4),ghjgg(4),ghjggZ(4),ghjhiZ(4,4),
     .     BR_hkhjhi(4,4,4),BR_hjinvisible(4),BR_hjhiZ(4,4),BR_hjemu(4),
     .     BR_hjHpiW(4,2)
      double precision :: SMGammaTotal(4), GRat(4), CP_value(4)
      double precision :: MHPlus(2), GammaTotalPlus(2),
     .     CS_lep_HpjHmj_ratio(2),BR_tWpb, BR_tHpjb(2),BR_Hpjcs(2),
     .     BR_Hpjcb(2),BR_Hpjtaunu(2),BR_Hpjtb(2),BR_HpjWZ(2),
     .     BR_HpjhiW(2,4)
C HiggsBounds functions for computing branching ratios:
      double precision :: SMBR_Htoptop,SMBR_Hss, SMBR_Hcc, SMBR_Hbb,
     .     SMBR_Hmumu, SMBR_Htautau,SMBR_HWW, SMBR_HZZ, SMBR_HZgam,
     .     SMBR_Hgamgam, SMBR_Hgg,SMGamma_h
      double precision :: kappa_V, kappa_f, kappa_gamma, kappa_Zgamma
      double precision :: csq_sep
C Local helper functions for couplings:
      DOUBLE PRECISION GET_KHLH3NZ, GET_KHHH3NZ, GET_KH5NH3NZ, invbsq,mt
C Loop counter:
      INTEGER k

C Input variables for the neutral Higgses. See HiggsBounds manual.              
      CP_value = (/ 1.D0, 1.D0, -1.D0, 1.D0 /)
      MHNeut = (/ MHL, MHH, MH3, MH5 /)
      GammaTotal = (/ HLWDTH, HHWDTH, H3NWDTH, H5NWDTH /)
      SMGammaTotal = (/ SMGamma_h(MHL), SMGamma_h(MHH), SMGamma_h(MH3),
     .   SMGamma_h(MH5) /)

      ghjss_s = (/ KFL, KFH, 0.D0, KF5 /)
      ghjss_p = (/ 0.D0, 0.D0, KF3, 0.D0 /)
      ghjcc_s = ghjss_s
      ghjcc_p = -ghjss_p
      ghjbb_s = ghjss_s
      ghjbb_p = ghjss_p
      ghjtt_s = ghjss_s
      ghjtt_p = -ghjss_p
      ghjmumu_s = ghjss_s
      ghjmumu_p = ghjss_p
      ghjtautau_s = ghjss_s
      ghjtautau_p = ghjss_p
      ghjWW = (/ KVL, KVH, KW3, KW5 /)
      ghjZZ = (/ KVL, KVH, KZ3, KZ5 /)
      ghjZga = DSQRT((/
     .     (HLBRZGA*HLWDTH)/(SMBR_HZgam(MHL)*SMGammaTotal(1)),
     .     (HHBRZGA*HHWDTH)/(SMBR_HZgam(MHH)*SMGammaTotal(2)),
     .     (H3NBRZGA*H3NWDTH)/(SMBR_HZgam(MH3)*SMGammaTotal(3)),
     .     (H5NBRZGA*H5NWDTH)/(SMBR_HZgam(MH5)*SMGammaTotal(4)) /))
      ghjgaga = DSQRT((/
     .     (HLBRGA*HLWDTH)/(SMBR_Hgamgam(MHL)*SMGammaTotal(1)),
     .     (HHBRGA*HHWDTH)/(SMBR_Hgamgam(MHH)*SMGammaTotal(2)),
     .     (H3NBRGA*H3NWDTH)/(SMBR_Hgamgam(MH3)*SMGammaTotal(3)),
     .     (H5NBRGA*H5NWDTH)/(SMBR_Hgamgam(MH5)*SMGammaTotal(4)) /))
      ghjgg = DSQRT((/
     .     (HLBRG*HLWDTH)/(SMBR_Hgg(MHL)*SMGammaTotal(1)),
     .     (HHBRG*HHWDTH)/(SMBR_Hgg(MHH)*SMGammaTotal(2)),
     .     (H3NBRG*H3NWDTH)/(SMBR_Hgg(MH3)*SMGammaTotal(3)),
     .     0.D0 /))
      ghjhiZ = RESHAPE( (/ 0.D0, 0.D0, GET_KHLH3NZ(), 0.D0,
     .     0.D0, 0.D0, GET_KHHH3NZ(), 0.D0,
     .     GET_KHLH3NZ(), GET_KHHH3NZ(), 0.D0, GET_KH5NH3NZ(),
     .     0.D0, 0.D0, GET_KH5NH3NZ(), 0.D0 /), (/ 4, 4 /) )

C If the SM BR and the GM model BR are both zero for a loop-induced
C decay, the relevant coupling modiification factor will be NaN at
C this point. If the SM BR is zero but the GM model BR is not, then
C the coupling modification factor will be Inf. Here, we fix this by
C replacing NaNs and Infs with zeros.
      DO k = 1,4
         IF (.NOT.ieee_is_finite(ghjZga(k))) THEN
            ghjZga(k) = 0.D0
         ENDIF
         IF (.NOT.ieee_is_finite(ghjgaga(k))) THEN
            ghjgaga(k) = 0.D0
         ENDIF
         IF (.NOT.ieee_is_finite(ghjgg(k))) THEN
            ghjgg(k) = 0.D0
         ENDIF
      ENDDO

C We must adjust the total widths. If we were to use the GMCALC widths,
C HiggsBounds would warn that the sum of BRs is greater than 1.
C Set kinematic top mass to same value as in HiggsBounds, for pseudoscalar
C ttbar width calculation only.
      mt = 172.6D0
      
      GammaTotal(1) = SMGammaTotal(1)*( ghjWW(1)**2*SMBR_HWW(MHL)
     .   +ghjZZ(1)**2*SMBR_HZZ(MHL) + ghjss_s(1)**2*(SMBR_Hss(MHL)
     .   +SMBR_Hcc(MHL)+SMBR_Hbb(MHL)+SMBR_Htoptop(MHL)
     .   +SMBR_Hmumu(MHL)+SMBR_Htautau(MHL))
     .   + ghjZga(1)**2*SMBR_HZgam(MHL)+ghjgg(1)**2*SMBR_Hgg(MHL)
     .   +ghjgaga(1)**2*SMBR_Hgamgam(MHL) )
     .   +HLWDTH*(HLBRH3N + HLBRH5N + HLBRZH3N + HLBRWH3P + HLBRH3P
     .   + HLBRH5P + HLBRH5PP)

      GammaTotal(2) = SMGammaTotal(2)*( ghjWW(2)**2*SMBR_HWW(MHH)
     .   +ghjZZ(2)**2*SMBR_HZZ(MHH) + ghjss_s(2)**2*(SMBR_Hss(MHH)
     .   +SMBR_Hcc(MHH)+SMBR_Hbb(MHH)+SMBR_Htoptop(MHH)
     .   +SMBR_Hmumu(MHH)+SMBR_Htautau(MHH))
     .   + ghjZga(2)**2*SMBR_HZgam(MHH)+ghjgg(2)**2*SMBR_Hgg(MHH)
     .   +ghjgaga(2)**2*SMBR_Hgamgam(MHH) )
     .   +HHWDTH*(HHBRHL + HHBRH3N + HHBRH5N + HHBRZH3N + HHBRWH3P
     .   + HHBRH3P + HHBRH5P + HHBRH5PP)

      GammaTotal(3) = SMGammaTotal(3)*( ghjWW(3)**2*SMBR_HWW(MH3)
     .   +ghjZZ(3)**2*SMBR_HZZ(MH3) + ghjss_p(3)**2*(SMBR_Hss(MH3)
     .   +SMBR_Hcc(MH3)+SMBR_Hbb(MH3)+invbsq(mt,MH3)*SMBR_Htoptop(MH3)
     .   +SMBR_Hmumu(MH3)+SMBR_Htautau(MH3))
     .   + ghjZga(3)**2*SMBR_HZgam(MH3)+ghjgg(3)**2*SMBR_Hgg(MH3)
     .   +ghjgaga(3)**2*SMBR_Hgamgam(MH3) )
     .   +H3NWDTH*(H3NBRZHL + H3NBRZHH + H3NBRZH5N + H3NBRWH5P)

      GammaTotal(4) = SMGammaTotal(4)*( ghjWW(4)**2*SMBR_HWW(MH5)
     .   +ghjZZ(4)**2*SMBR_HZZ(MH5) + ghjZga(4)**2*SMBR_HZgam(MH5)
     .   +ghjgg(4)**2*SMBR_Hgg(MH5)+ghjgaga(4)**2*SMBR_Hgamgam(MH5) )
     .   +H5NWDTH*(H5NBRH3N + H5NBRZH3N + H5NBRWH3P + H5NBRH3P)

C Calculating non-SM branching ratios.
      GRat = (/ HLWDTH, HHWDTH, H3NWDTH, H5NWDTH /) / GammaTotal
      BR_hjinvisible = 0.D0
      BR_hkhjhi = 0.D0
      BR_hkhjhi(2,1,1) = HHBRHL*GRat(2)
      BR_hkhjhi(1,3,3) = HLBRH3N*GRat(1)
      BR_hkhjhi(2,3,3) = HHBRH3N*GRat(2)
      BR_hkhjhi(4,3,3) = H5NBRH3N*GRat(4)
      BR_hkhjhi(1,4,4) = HLBRH5N*GRat(1)
      BR_hkhjhi(2,4,4) = HHBRH5N*Grat(2)
      BR_hjhiZ = RESHAPE( (/ 0.D0, 0.D0, H3NBRZHL*GRat(3), 0.D0,
     .     0.D0, 0.D0, H3NBRZHH*GRat(3), 0.D0,
     .     HLBRZH3N*GRat(1), HHBRZH3N*GRat(2), 0.D0, H5NBRZH3N*GRat(4),
     .     0.D0, 0.D0, H3NBRZH5N*GRat(3), 0.D0 /), (/ 4, 4 /) )
      BR_hjemu = (/ 0.D0, 0.D0, 0.D0, 0.D0 /) ! will also use this for BR_hjetau and BR_hjmutau
      BR_hjHpiW = RESHAPE( (/ HLBRWH3P*GRat(1), HHBRWH3P*GRat(2), 0.D0,
     .     H5NBRWH3P*GRat(4), 0.D0, 0.D0, H3NBRWH5P*GRat(3), 0.D0 /),
     .     (/ 4, 2 /) )

C Input parameters for charged Higgses. See HiggsBounds manual.
      MHPlus = (/ MH3, MH5 /)
      GammaTotalPlus = (/ H3PWDTH, H5PWDTH /)
      CS_lep_HpjHmj_ratio = (/ 1.D0, 1.D0 /)
      BR_tWpb = TOPBRW
      BR_tHpjb = (/ TOPBRH3P, 0.D0 /)
      BR_Hpjcs = (/ H3PBRCS, 0.D0 /)
      BR_Hpjcb = (/ H3PBRBC, 0.D0 /)
      BR_Hpjtaunu = (/ H3PBRTA, 0.D0 /)
      BR_Hpjtb = (/ H3PBRTB, 0.D0 /)
      BR_HpjWZ = (/ 0.D0, H5PBRWZ /)
      BR_HpjhiW = RESHAPE( (/ H3PBRWHL, 0.D0, H3PBRWHH, 0.D0, 0.D0,
     .     H5PBRWH3N, H3PBRWH5N, 0.D0 /), (/ 2, 4 /) )

C Inputting all the parameters that were calculated
      CALL HiggsBounds_neutral_input_properties(MHNeut, GammaTotal,
     .     CP_value)
      CALL HiggsBounds_neutral_input_nonSMBR(BR_hjinvisible, BR_hkhjhi,
     .     BR_hjhiZ, BR_hjemu, BR_hjemu, BR_hjemu, BR_hjHpiW)
      CALL HiggsBounds_neutral_input_effC(ghjss_s, ghjss_p, ghjcc_s,
     .     ghjcc_p,ghjbb_s,ghjbb_p,ghjtt_s,ghjtt_p,ghjmumu_s,ghjmumu_p,
     .     ghjtautau_s,ghjtautau_p,ghjWW,ghjZZ,ghjZga,ghjgaga,ghjgg,
     .     ghjhiZ)
      CALL HiggsBounds_charged_input(MHPlus, GammaTotalPlus,
     .     CS_lep_HpjHmj_ratio, BR_tWpb, BR_tHpjb, BR_Hpjcs,
     .     BR_Hpjcb, BR_Hpjtaunu, BR_Hpjtb, BR_HpjWZ, BR_HpjhiW)

C Running HiggsBounds
      CALL run_HiggsBounds_full(HBresult, chan, obsratio, ncombined)

C Running HiggsSignals (the first argument is the mode, 1 = peak centred)
      CALL run_HiggsSignals(1, Chisq_mu, Chisq_mh, Chisq, nobs, Pvalue)
      CALL run_HiggsSignals_LHC_Run1_combination(Chisq_mu_LHCRun1, 
     .     Chisq_mh_LHCRun1, Chisq_LHCRun1, nobs_LHCRun1, 
     .     Pvalue_LHCRun1) 
      CALL run_HiggsSignals_STXS(Chisq_STXS_rates, Chisq_STXS_mh, 
     .     Chisq_STXS, nobs_STXS, Pvalue_STXS)
      CALL complete_HS_results()
      Pvalue = HSres(1)%Pvalue

      RETURN
      END

C====================================================================
C Local helper functions for subroutine CALLHB5:
C====================================================================
      FUNCTION GET_KHLH3NZ()
         DOUBLE PRECISION GET_KHLH3NZ, VSM
         VSM = DSQRT(VPHI**2 + 8.D0*VCHI**2)
         GET_KHLH3NZ=-DSQRT(1.D0/6.D0)*(DSQRT(3.D0)*COS(ALPHA)*VCHI/VSM
     .           +SIN(ALPHA)*VPHI/VSM)
         RETURN
      END FUNCTION

      FUNCTION GET_KHHH3NZ()
         DOUBLE PRECISION GET_KHHH3NZ, VSM
         VSM = DSQRT(VPHI**2 + 8.D0*VCHI**2)
         GET_KHHH3NZ=DSQRT(1.D0/6.D0)*(-DSQRT(3.D0)*SIN(ALPHA)*VCHI/VSM
     .           +COS(ALPHA)*VPHI/VSM)
         RETURN
      END FUNCTION

      FUNCTION GET_KH5NH3NZ()
         DOUBLE PRECISION GET_KH5NH3NZ, VSM
         VSM = DSQRT(VPHI**2 + 8.D0*VCHI**2)
         GET_KH5NH3NZ = (1.D0/2.D0)*DSQRT(1.D0/3.D0)*VPHI/VSM
         RETURN
      END FUNCTION
C====================================================================
      function invbsq(mf,mh)
         double precision,intent(in) :: mf,mh
         double precision :: invbsq
         if(mh>2.0D0*mf)then
            invbsq=1.0D0/(1.0D0-4.0D0*(mf/mh)**2.0D0)
         else
            invbsq=0.0D0
         endif
      end function invbsq
