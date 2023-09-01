C==================================================================
C     Subroutines to calculate the scalar branching ratios & widths
C==================================================================

      SUBROUTINE CALCDECAYS
C A quick way to call the 7 subroutines that calculate the decay
C partial widths of the physical scalars.
      CALL HLDECAYS
      CALL HHDECAYS
      CALL H3NDECAYS
      CALL H3PDECAYS
      CALL H5NDECAYS
      CALL H5PDECAYS
      CALL H5PPDECAYS
      CALL TOPDECAYS
      RETURN
      END

C------------------------------------------------------------------
C Subroutine to calculate the decays of h = HL

      SUBROUTINE HLDECAYS
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
      DOUBLE PRECISION RxH5PWGA,RxH3PWGA,RxH3PWGATILDE,RxH5NGAGA,
     .     RxH5NZGA,RxH3NGAGA,RxH3NZGA,RxH3NGG,RxHHGAGA,RxHHZGA,RxHHGG,
     .     RxHLGAGA,RxHLZGA,RxHLGG,
     .     IxH5PWGA,IxH3PWGA,IxH3PWGATILDE
      COMMON/LOOPPARAMS/RxH5PWGA,RxH3PWGA,RxH3PWGATILDE,RxH5NGAGA,
     .     RxH5NZGA,RxH3NGAGA,RxH3NZGA,RxH3NGG,RxHHGAGA,RxHHZGA,RxHHGG,
     .     RxHLGAGA,RxHLZGA,RxHLGG,
     .     IxH5PWGA,IxH3PWGA,IxH3PWGATILDE
      DOUBLE PRECISION HLBRB, HLBRTA, HLBRMU, HLBRS, HLBRC, HLBRT,
     .     HLBRG, HLBRGA, HLBRZGA, HLBRW, HLBRZ, 
     .     HLBRWH3P, HLBRZH3N,
     .     HLBRH3N, HLBRH3P, HLBRH5N, HLBRH5P, HLBRH5PP, HLWDTH
      COMMON/HLBRS/HLBRB, HLBRTA, HLBRMU, HLBRS, HLBRC, HLBRT,
     .     HLBRG, HLBRGA, HLBRZGA, HLBRW, HLBRZ, 
     .     HLBRWH3P, HLBRZH3N,
     .     HLBRH3N, HLBRH3P, HLBRH5N, HLBRH5P, HLBRH5PP, HLWDTH
      INTEGER OFFSHELL, QCDCORRS
      COMMON/DECAYFLAGS/OFFSHELL, QCDCORRS
C Local variables:
      DOUBLE PRECISION HLGAMB, HLGAMTA, HLGAMMU, HLGAMS, HLGAMC, 
     .     HLGAMT, 
     .     HLGAMG, HLGAMGA, HLGAMZGA, HLGAMW, HLGAMZ, 
     .     HLGAMWH3P, HLGAMZH3N,
     .     HLGAMH3N, HLGAMH3P, HLGAMH5N, HLGAMH5P, HLGAMH5PP
      DOUBLE PRECISION KF, KV, CHHZ, CHHW, GHL33, GHL55, CZ11, CZ22
      DOUBLE PRECISION TAUT, TAUW, LAMBDAT, LAMBDAW
      DOUBLE COMPLEX AT, AW, AS, A3P, A5P, A5PP, AMP
      DOUBLE PRECISION EE5
      DOUBLE PRECISION SW, CW, PI
      DOUBLE PRECISION MTRUN, MBRUN, MCRUN, MT, MB, MC
C Functions to be called:
      DOUBLE PRECISION GAMFF, GAMVVOF, GAMVV, GAMWH, GAMZH, GAMHH
      DOUBLE COMPLEX F0, F12, F1, I1, I2
      DOUBLE PRECISION RUNMB, RUNMC, DELTAQCD, DELTATB, DELTATC, RUNALS,
     .     POLEMQ, RUNMT
      PI = 4.D0*DATAN(1.D0)
      CW = MW/MZ
      SW = DSQRT(1.D0 - CW**2)

      IF (QCDCORRS.EQ.1) THEN
C Running quark masses, for Yukawa couplings
         MTRUN = RUNMT(MHL)
         MBRUN = RUNMB(MHL)
         MCRUN = RUNMC(MHL)
C Quark pole masses, for kinematics
         MT = MTPOLE
         MB = POLEMQ(MBMB,5)
         MC = POLEMQ(MCMC,4)
      ELSE
         MTRUN = MTPOLE
         MBRUN = MBMB
         MCRUN = MCMC
         MT = MTPOLE
         MB = MBMB
         MC = MCMC
      ENDIF

C HL decays to fermions:  on-shell for now.
      KF = DCOS(ALPHA)*V/VPHI
      IF (QCDCORRS.EQ.1) THEN
         HLGAMB = GAMFF(3.D0,MHL,MB,MB,MBRUN/V*KF,0.D0)
     .        * (DELTAQCD(MHL) + DELTATB(MHL))
         HLGAMC = GAMFF(3.D0,MHL,MC,MC,MCRUN/V*KF,0.D0)
     .        * (DELTAQCD(MHL) + DELTATC(MHL))
      ELSE
         HLGAMB = GAMFF(3.D0,MHL,MB,MB,MBRUN/V*KF,0.D0)
         HLGAMC = GAMFF(3.D0,MHL,MC,MC,MCRUN/V*KF,0.D0)
      ENDIF
      HLGAMTA = GAMFF(1.D0,MHL,MTAU,MTAU,MTAU/V*KF,0.D0)
      HLGAMMU = GAMFF(1.D0,MHL,MMU,MMU,MMU/V*KF,0.D0)
      HLGAMS = GAMFF(3.D0,MHL,MS,MS,MS/V*KF,0.D0)
      HLGAMT = GAMFF(3.D0,MHL,MT,MT,MTRUN/V*KF,0.D0)

C HL decays to massive gauge bosons: doubly offshell unless switched off
      KV = DCOS(ALPHA)*VPHI/V - 8.D0/DSQRT(3.D0)*DSIN(ALPHA)*VCHI/V
      IF (OFFSHELL.EQ.1) THEN
         HLGAMW = GAMVVOF(1.D0,MHL,MW,MW,GAMW,GAMW,2.D0*MW**2/V*KV)
         HLGAMZ = GAMVVOF(0.5D0,MHL,MZ,MZ,GAMZ,GAMZ,2.D0*MZ**2/V*KV)
      ELSE
C onshell option (OFFSHELL = 0)
         HLGAMW = GAMVV(1.D0,MHL,MW,MW,2.D0*MW**2/V*KV)
         HLGAMZ = GAMVV(0.5D0,MHL,MZ,MZ,2.D0*MZ**2/V*KV)
      ENDIF

C HL decays to vector + scalar (W+H- and W-H+ are summed):
C On-shell above MH2+MV, singly off-shell below (the V is off-shell).
C CHHZ is pure imaginary; we remove the i (since it is mod-squared in GAMVH).
C HLGAMWH3P is the sum of HL -> W+ H3- and HL -> W- H3+
      CHHW = -2.D0*MW/V * (DSQRT(2.D0)*VCHI/V*DCOS(ALPHA)
     .     + DSQRT(2.D0/3.D0)*VPHI/V*DSIN(ALPHA))
      HLGAMWH3P = 2.D0 * GAMWH(MHL,MH3,MW,CHHW,V)
      CHHZ = -2.D0*MZ/V * (DSQRT(2.D0)*VCHI/V*DCOS(ALPHA) 
     .     + DSQRT(2.D0/3.D0)*VPHI/V*DSIN(ALPHA))
      HLGAMZH3N = GAMZH(MHL,MH3,MZ,CHHZ,V,SW)

C HL decays to two scalars: on-shell only.
      GHL33 = 64.D0*LAMBDA1*DCOS(ALPHA)*VCHI**2*VPHI/V**2
     .     - 8.D0/DSQRT(3.D0)*VPHI**2*VCHI/V**2
     .     *DSIN(ALPHA)*(LAMBDA3+3.D0*LAMBDA4)
     .     - 4.D0/DSQRT(3.D0)*VCHI*M1/V**2
     .     *(DSIN(ALPHA)*VCHI - DSQRT(3.D0)*DCOS(ALPHA)*VPHI)
     .     - 16.D0/DSQRT(3.D0)*VCHI**3/V**2*DSIN(ALPHA)
     .     *(6.D0*LAMBDA2+LAMBDA5)
     .     - DCOS(ALPHA)*VPHI**3/V**2*(LAMBDA5-4.D0*LAMBDA2)
     .     + 2.D0*DSQRT(3.D0)*M2*VPHI**2/V**2*DSIN(ALPHA)
     .     - 8.D0/DSQRT(3.D0)*LAMBDA5*VCHI*VPHI/V**2
     .     *(DSIN(ALPHA)*VPHI - DSQRT(3.D0)*DCOS(ALPHA)*VCHI)
      GHL55 = -8.D0*DSQRT(3.D0)*(LAMBDA3+LAMBDA4)*DSIN(ALPHA)*VCHI
     .     + (4.D0*LAMBDA2+LAMBDA5)*DCOS(ALPHA)*VPHI
     .     - 2.D0*DSQRT(3.D0)*M2*DSIN(ALPHA)
      HLGAMH3N = GAMHH(0.5D0,MHL,MH3,MH3,GHL33)
      HLGAMH3P = GAMHH(1.D0,MHL,MH3,MH3,GHL33)
      HLGAMH5N = GAMHH(0.5D0,MHL,MH5,MH5,GHL55)
      HLGAMH5P = GAMHH(1.D0,MHL,MH5,MH5,GHL55)
      HLGAMH5PP = GAMHH(1.D0,MHL,MH5,MH5,GHL55)

C HL decays to massless gauge bosons (loop-induced)
      TAUT = 4.D0*MT**2/MHL**2
      TAUW = 4.D0*MW**2/MHL**2
      LAMBDAT = 4.D0*MT**2/MZ**2
      LAMBDAW = 4.D0*MW**2/MZ**2

C Loop-induced h -> gg: include QCD corrections, Nf = 5 light flavors, 
C renorm scale = MHL.
C SM amplitude for the top loop:
      AT = -F12(TAUT)
      AMP = KF*AT
      IF (QCDCORRS.EQ.1) THEN
         EE5 = 95.D0/4.D0 - 7.D0/6.D0*5.D0
         HLGAMG = RUNALS(MHL,5)**2 * MHL**3 / 128.D0/PI**3 / V**2
     .        * CDABS(AMP)**2
     .        * (1.D0 + EE5*RUNALS(MHL,5)/PI)
         RxHLGG = RUNALS(MHL,5) / 4.D0/PI / V
     .        * CDABS(AMP)
     .        * DSQRT(1.D0 + EE5*RUNALS(MHL,5)/PI)
      ELSE
         HLGAMG = ALSMZ**2 * MHL**3 / 128.D0/PI**3 / V**2
     .        * CDABS(AMP)**2
         RxHLGG = ALSMZ / 4.D0/PI / V
     .        * CDABS(AMP)
      ENDIF

C Loop-induced h -> gamma gamma:
C SM amplitude for the top loop:
      AT = 3.D0 * (2.D0/3.D0)**2 * F12(TAUT)
C SM amplitude for the W loop:
      AW = F1(TAUW)
C Amplitudes for the scalar loops:
      A3P = GHL33*V/2.D0/MH3**2 * F0(4.D0*MH3**2/MHL**2)
      A5P = GHL55*V/2.D0/MH5**2 * F0(4.D0*MH5**2/MHL**2)
      A5PP = GHL55*V/2.D0/MH5**2 * 4.D0 * F0(4.D0*MH5**2/MHL**2)
      AS = A3P + A5P + A5PP
C Combine the amplitudes (all are complex in general):
      AMP = KF*AT + KV*AW + AS
      HLGAMGA = ALPHAEM**2 * MHL**3 / 256.D0/PI**3 / V**2
     .     * CDABS(AMP)**2
      RxHLGAGA = ALPHAEM / 2.D0/PI / V
     .     * CDABS(AMP)
C Loop-induced h -> Z gamma coupling:
      IF (MHL.LE.MZ) THEN
         HLGAMZGA = 0.D0
         RxHLZGA = 0.D0
      ELSE
C SM amplitude for the top loop:
         AT = 3.D0 * (-4.D0/3.D0) * (0.5D0 - 4.D0/3.D0*SW**2)/SW/CW
     .        * (I1(TAUT,LAMBDAT) - I2(TAUT,LAMBDAT))
C SM amplitude for the W loop:
         AW = -CW/SW * (4.D0*(3.D0-SW**2/CW**2) * I2(TAUW,LAMBDAW)
     .        + ((1.D0+2.D0/TAUW)*SW**2/CW**2 - (5.D0+2.D0/TAUW)) 
     .        * I1(TAUW,LAMBDAW))
C Amplitudes for the scalar loops:
         CZ11 = 1.D0/2.D0/SW/CW * (1.D0 - 2.D0*SW**2)
         CZ22 = 1.D0/SW/CW * (1.D0 - 2.D0*SW**2)
         A3P = 2.D0*GHL33*CZ11/MH3**2 
     .        * I1(4.D0*MH3**2/MHL**2,4.D0*MH3**2/MZ**2)
         A5P = 2.D0*GHL55*CZ11/MH5**2
     .        * I1(4.D0*MH5**2/MHL**2,4.D0*MH5**2/MZ**2)
         A5PP = 2.D0*GHL55*CZ22*2.D0/MH5**2
     .        * I1(4.D0*MH5**2/MHL**2,4.D0*MH5**2/MZ**2)
         AS = A3P + A5P + A5PP
C Combine the amplitudes (all are complex in general):
         AMP = KF*AT + KV*AW + V/2.D0*AS
         HLGAMZGA = ALPHAEM**2 * MHL**3 / 128.D0/PI**3 / V**2
     .        * CDABS(AMP)**2 * (1.D0 - MZ**2/MHL**2)**3
         RxHLZGA = ALPHAEM / 2.D0/PI / V
     .        * CDABS(AMP)
      ENDIF

C Compute the total width
      HLWDTH = HLGAMB + HLGAMTA + HLGAMMU + HLGAMS + HLGAMC + HLGAMT
     .     + HLGAMG + HLGAMGA + HLGAMZGA + HLGAMW + HLGAMZ
     .     + HLGAMWH3P + HLGAMZH3N
     .     + HLGAMH3N + HLGAMH3P + HLGAMH5N + HLGAMH5P + HLGAMH5PP

C Compute the BRs
      HLBRB = HLGAMB/HLWDTH
      HLBRTA = HLGAMTA/HLWDTH
      HLBRMU = HLGAMMU/HLWDTH
      HLBRS = HLGAMS/HLWDTH
      HLBRC = HLGAMC/HLWDTH
      HLBRT = HLGAMT/HLWDTH
      HLBRG = HLGAMG/HLWDTH
      HLBRGA = HLGAMGA/HLWDTH
      HLBRZGA = HLGAMZGA/HLWDTH
      HLBRW = HLGAMW/HLWDTH
      HLBRZ = HLGAMZ/HLWDTH
      HLBRWH3P = HLGAMWH3P/HLWDTH
      HLBRZH3N = HLGAMZH3N/HLWDTH
      HLBRH3N = HLGAMH3N/HLWDTH
      HLBRH3P = HLGAMH3P/HLWDTH
      HLBRH5N = HLGAMH5N/HLWDTH
      HLBRH5P = HLGAMH5P/HLWDTH
      HLBRH5PP = HLGAMH5PP/HLWDTH

      RETURN
      END

C------------------------------------------------------------------
C Subroutine to calculate the decays of H = HH

      SUBROUTINE HHDECAYS
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
      INTEGER OFFSHELL, QCDCORRS
      COMMON/DECAYFLAGS/OFFSHELL, QCDCORRS
C Local variables:
      DOUBLE PRECISION HHGAMB, HHGAMTA, HHGAMMU, HHGAMS, HHGAMC, 
     .     HHGAMT, 
     .     HHGAMG, HHGAMGA, HHGAMZGA, HHGAMW, HHGAMZ, 
     .     HHGAMWH3P, HHGAMZH3N,
     .     HHGAMHL, HHGAMH3N, HHGAMH3P, HHGAMH5N, HHGAMH5P, HHGAMH5PP
      DOUBLE PRECISION KF, KV, CHHZ, CHHW, GHLL, GHH33, GHH55, 
     .     CZ11, CZ22
      DOUBLE PRECISION TAUT, TAUW, LAMBDAT, LAMBDAW
      DOUBLE COMPLEX AT, AW, AS, A3P, A5P, A5PP, AMP
      DOUBLE PRECISION EE5
      DOUBLE PRECISION SW, CW, PI
      DOUBLE PRECISION MTRUN, MBRUN, MCRUN, MT, MB, MC
C Functions to be called:
      DOUBLE PRECISION GAMFF, GAMVVOF, GAMVV, GAMWH, GAMZH, GAMHH
      DOUBLE COMPLEX F0, F12, F1, I1, I2
      DOUBLE PRECISION RUNMB, RUNMC, DELTAQCD, DELTATB, DELTATC, RUNALS,
     .     POLEMQ, RUNMT
      PI = 4.D0*DATAN(1.D0)
      CW = MW/MZ
      SW = DSQRT(1.D0 - CW**2)

      IF (QCDCORRS.EQ.1) THEN
C Running quark masses, for Yukawa couplings
         MTRUN = RUNMT(MHH)
         MBRUN = RUNMB(MHH)
         MCRUN = RUNMC(MHH)
C     Quark pole masses, for kinematics
         MT = MTPOLE
         MB = POLEMQ(MBMB,5)
         MC = POLEMQ(MCMC,4)
      ELSE
         MTRUN = MTPOLE
         MBRUN = MBMB
         MCRUN = MCMC
         MT = MTPOLE
         MB = MBMB
         MC = MCMC
      ENDIF

C HH decays to fermions: on-shell for now.
      KF = DSIN(ALPHA)*V/VPHI
      IF (QCDCORRS.EQ.1) THEN
         HHGAMB = GAMFF(3.D0,MHH,MB,MB,MBRUN/V*KF,0.D0)
     .        * (DELTAQCD(MHH) + DELTATB(MHH))
         HHGAMC = GAMFF(3.D0,MHH,MC,MC,MCRUN/V*KF,0.D0)
     .        * (DELTAQCD(MHH) + DELTATC(MHH))
      ELSE
         HHGAMB = GAMFF(3.D0,MHH,MB,MB,MBRUN/V*KF,0.D0)
         HHGAMC = GAMFF(3.D0,MHH,MC,MC,MCRUN/V*KF,0.D0)
      ENDIF
      HHGAMTA = GAMFF(1.D0,MHH,MTAU,MTAU,MTAU/V*KF,0.D0)
      HHGAMMU = GAMFF(1.D0,MHH,MMU,MMU,MMU/V*KF,0.D0)
      HHGAMS = GAMFF(3.D0,MHH,MS,MS,MS/V*KF,0.D0)
      HHGAMT = GAMFF(3.D0,MHH,MT,MT,MTRUN/V*KF,0.D0)

C HH decays to massive gauge bosons: doubly offshell unless switched off
      KV = DSIN(ALPHA)*VPHI/V + 8.D0/DSQRT(3.D0)*DCOS(ALPHA)*VCHI/V
      IF (OFFSHELL.EQ.1) THEN
         HHGAMW = GAMVVOF(1.D0,MHH,MW,MW,GAMW,GAMW,2.D0*MW**2/V*KV)
         HHGAMZ = GAMVVOF(0.5D0,MHH,MZ,MZ,GAMZ,GAMZ,2.D0*MZ**2/V*KV)
      ELSE
C onshell option (OFFSHELL = 0)
         HHGAMW = GAMVV(1.D0,MHH,MW,MW,2.D0*MW**2/V*KV)
         HHGAMZ = GAMVV(0.5D0,MHH,MZ,MZ,2.D0*MZ**2/V*KV)
      ENDIF

C HH decays to vector + scalar (W+H- and W-H+ are summed):
C On-shell above threshold, V off-shell below threshold.
C CHHZ is pure imaginary; we remove the i (since it is mod-squared in GAMVH).
C HHGAMWH3P is the sum of HH -> W+ H3- and HH -> W- H3+
      CHHW = -2.D0*MW/V * (DSQRT(2.D0)*VCHI/V*DSIN(ALPHA)
     .     - DSQRT(2.D0/3.D0)*VPHI/V*DCOS(ALPHA))
      HHGAMWH3P = 2.D0 * GAMWH(MHH,MH3,MW,CHHW,V)
      CHHZ = -2.D0*MZ/V * (DSQRT(2.D0)*VCHI/V*DSIN(ALPHA) 
     .     - DSQRT(2.D0/3.D0)*VPHI/V*DCOS(ALPHA))
      HHGAMZH3N = GAMZH(MHH,MH3,MZ,CHHZ,V,SW)

C HH decays to two scalars: on-shell only.
      GHLL = 24.D0*LAMBDA1*DCOS(ALPHA)**2*DSIN(ALPHA)*VPHI
     .     + 2.D0*(DSQRT(3.D0)*DCOS(ALPHA)*VCHI
     .     *(3.D0*DCOS(ALPHA)**2 - 2.D0)
     .     + DSIN(ALPHA)*VPHI*(1.D0 - 3.D0*DCOS(ALPHA)**2))
     .     * (2.D0*LAMBDA2 - LAMBDA5)
     .     + 8.D0*DSQRT(3.D0)*DCOS(ALPHA)*DSIN(ALPHA)**2
     .     *VCHI*(LAMBDA3 + 3.D0*LAMBDA4)
     .     - DSQRT(3.D0)/2.D0*M1*DCOS(ALPHA)
     .     *(3.D0*DCOS(ALPHA)**2 - 2.D0)
     .     - 4.D0*DSQRT(3.D0)*M2*DCOS(ALPHA)*DSIN(ALPHA)**2
      GHH33 = 64.D0*LAMBDA1*DSIN(ALPHA)*VCHI**2*VPHI/V**2
     .     + 8.D0/DSQRT(3.D0)*VPHI**2*VCHI/V**2 
     .     *DCOS(ALPHA)*(LAMBDA3 + 3.D0*LAMBDA4)
     .     + 4.D0/DSQRT(3.D0)*VCHI*M1/V**2
     .     *(DCOS(ALPHA)*VCHI + DSQRT(3.D0)*DSIN(ALPHA)*VPHI)
     .     + 16.D0/DSQRT(3.D0)*VCHI**3/V**2*DCOS(ALPHA)
     .     *(6.D0*LAMBDA2 + LAMBDA5)
     .     + DSIN(ALPHA)*VPHI**3/V**2*(4.D0*LAMBDA2 - LAMBDA5)
     .     - 2.D0*DSQRT(3.D0)*M2*VPHI**2/V**2*DCOS(ALPHA)
     .     + 8.D0/DSQRT(3.D0)*LAMBDA5*VCHI*VPHI/V**2
     .     * (DCOS(ALPHA)*VPHI + DSQRT(3.D0)*DSIN(ALPHA)*VCHI)
      GHH55 = 8.D0*DSQRT(3.D0)*(LAMBDA3+LAMBDA4)*DCOS(ALPHA)*VCHI
     .     + (4.D0*LAMBDA2+LAMBDA5)*DSIN(ALPHA)*VPHI
     .     + 2.D0*DSQRT(3.D0)*M2*DCOS(ALPHA)
      HHGAMHL = GAMHH(0.5D0,MHH,MHL,MHL,GHLL)
      HHGAMH3N = GAMHH(0.5D0,MHH,MH3,MH3,GHH33)
      HHGAMH3P = GAMHH(1.D0,MHH,MH3,MH3,GHH33)
      HHGAMH5N = GAMHH(0.5D0,MHH,MH5,MH5,GHH55)
      HHGAMH5P = GAMHH(1.D0,MHH,MH5,MH5,GHH55)
      HHGAMH5PP = GAMHH(1.D0,MHH,MH5,MH5,GHH55)

C HH decays to massless gauge bosons (loop-induced)
      TAUT = 4.D0*MT**2/MHH**2
      TAUW = 4.D0*MW**2/MHH**2
      LAMBDAT = 4.D0*MT**2/MZ**2
      LAMBDAW = 4.D0*MW**2/MZ**2

C Loop-induced H -> gg:
C SM amplitude for the top loop:
      AT = -F12(TAUT)
      AMP = KF*AT
      IF (QCDCORRS.EQ.1) THEN
         EE5 = 95.D0/4.D0 - 7.D0/6.D0*5.D0
         HHGAMG = RUNALS(MHH,5)**2 * MHH**3 / 128.D0/PI**3 / V**2
     .        * CDABS(AMP)**2
     .        * (1.D0 + EE5*RUNALS(MHH,5)/PI)
         RxHHGG = RUNALS(MHH,5) / 4.D0/PI / V
     .        * CDABS(AMP)
     .        * DSQRT(1.D0 + EE5*RUNALS(MHH,5)/PI)
      ELSE
         HHGAMG = ALSMZ**2 * MHH**3 / 128.D0/PI**3 / V**2
     .        * CDABS(AMP)**2
         RxHHGG = ALSMZ / 4.D0/PI / V
     .        * CDABS(AMP)
      ENDIF

C Loop-induced H -> gamma gamma:
C SM amplitude for the top loop:
      AT = 3.D0 * (2.D0/3.D0)**2 * F12(TAUT)
C SM amplitude for the W loop:
      AW = F1(TAUW)
C Amplitudes for the scalar loops:
      A3P = GHH33*V/2.D0/MH3**2 * F0(4.D0*MH3**2/MHH**2)
      A5P = GHH55*V/2.D0/MH5**2 * F0(4.D0*MH5**2/MHH**2)
      A5PP = GHH55*V/2.D0/MH5**2 * 4.D0 * F0(4.D0*MH5**2/MHH**2)
      AS = A3P + A5P + A5PP
C Combine the amplitudes (all are complex in general):
      AMP = KF*AT + KV*AW + AS
      HHGAMGA = ALPHAEM**2 * MHH**3 / 256.D0/PI**3 / V**2
     .     * CDABS(AMP)**2
      RxHHGAGA = ALPHAEM / 2.D0/PI / V
     .     * CDABS(AMP)
C Loop-induced H -> Z gamma coupling:
      IF (MHH.LE.MZ) THEN
         HHGAMZGA = 0.D0
         RxHHZGA = 0.D0
      ELSE
C SM amplitude for the top loop:
         AT = 3.D0 * (-4.D0/3.D0) * (0.5D0 - 4.D0/3.D0*SW**2)/SW/CW
     .        * (I1(TAUT,LAMBDAT) - I2(TAUT,LAMBDAT))
C SM amplitude for the W loop:
         AW = -CW/SW * (4.D0*(3.D0-SW**2/CW**2) * I2(TAUW,LAMBDAW)
     .        + ((1.D0+2.D0/TAUW)*SW**2/CW**2 - (5.D0+2.D0/TAUW)) 
     .        * I1(TAUW,LAMBDAW))
C Amplitudes for the scalar loops:
         CZ11 = 1.D0/2.D0/SW/CW * (1.D0 - 2.D0*SW**2)
         CZ22 = 1.D0/SW/CW * (1.D0 - 2.D0*SW**2)
         A3P = 2.D0*GHH33*CZ11/MH3**2 
     .        * I1(4.D0*MH3**2/MHH**2,4.D0*MH3**2/MZ**2)
         A5P = 2.D0*GHH55*CZ11/MH5**2
     .        * I1(4.D0*MH5**2/MHH**2,4.D0*MH5**2/MZ**2)
         A5PP = 2.D0*GHH55*CZ22*2.D0/MH5**2
     .        * I1(4.D0*MH5**2/MHH**2,4.D0*MH5**2/MZ**2)
         AS = A3P + A5P + A5PP
C Combine the amplitudes (all are complex in general):
         AMP = KF*AT + KV*AW + V/2.D0*AS
         HHGAMZGA = ALPHAEM**2 * MHH**3 / 128.D0/PI**3 / V**2
     .        * CDABS(AMP)**2 * (1.D0 - MZ**2/MHH**2)**3
         RxHHZGA = ALPHAEM / 2.D0/PI / V
     .        * CDABS(AMP)
      ENDIF

C Compute the total width
      HHWDTH = HHGAMB + HHGAMTA + HHGAMMU + HHGAMS + HHGAMC + HHGAMT
     .     + HHGAMG + HHGAMGA + HHGAMZGA + HHGAMW + HHGAMZ
     .     + HHGAMWH3P + HHGAMZH3N
     .     + HHGAMHL + HHGAMH3N + HHGAMH3P + HHGAMH5N + HHGAMH5P 
     .     + HHGAMH5PP

C Compute the BRs
      HHBRB = HHGAMB/HHWDTH
      HHBRTA = HHGAMTA/HHWDTH
      HHBRMU = HHGAMMU/HHWDTH
      HHBRS = HHGAMS/HHWDTH
      HHBRC = HHGAMC/HHWDTH
      HHBRT = HHGAMT/HHWDTH
      HHBRG = HHGAMG/HHWDTH
      HHBRGA = HHGAMGA/HHWDTH
      HHBRZGA = HHGAMZGA/HHWDTH
      HHBRW = HHGAMW/HHWDTH
      HHBRZ = HHGAMZ/HHWDTH
      HHBRWH3P = HHGAMWH3P/HHWDTH
      HHBRZH3N = HHGAMZH3N/HHWDTH
      HHBRHL = HHGAMHL/HHWDTH
      HHBRH3N = HHGAMH3N/HHWDTH
      HHBRH3P = HHGAMH3P/HHWDTH
      HHBRH5N = HHGAMH5N/HHWDTH
      HHBRH5P = HHGAMH5P/HHWDTH
      HHBRH5PP = HHGAMH5PP/HHWDTH

      RETURN
      END

C------------------------------------------------------------------
C Subroutine to calculate the decays of H_3^0 = H3N (pseudoscalar)

      SUBROUTINE H3NDECAYS
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
      DOUBLE PRECISION H3NBRB, H3NBRTA, H3NBRMU, H3NBRS, H3NBRC, H3NBRT,
     .     H3NBRZHL, H3NBRZHH, H3NBRZH5N, H3NBRWH5P,
     .     H3NBRG, H3NBRGA, H3NBRZGA, 
     .     H3NWDTH
      COMMON/H3NBRS/H3NBRB, H3NBRTA, H3NBRMU, H3NBRS, H3NBRC, H3NBRT,
     .     H3NBRZHL, H3NBRZHH, H3NBRZH5N, H3NBRWH5P,
     .     H3NBRG, H3NBRGA, H3NBRZGA,
     .     H3NWDTH
      INTEGER OFFSHELL, QCDCORRS
      COMMON/DECAYFLAGS/OFFSHELL, QCDCORRS
C Local variables:
      DOUBLE PRECISION H3NGAMB, H3NGAMTA, H3NGAMMU, H3NGAMS, H3NGAMC,
     .     H3NGAMT,
     .     H3NGAMZHL, H3NGAMZHH, H3NGAMZH5N, H3NGAMWH5P,
     .     H3NGAMG, H3NGAMGA, H3NGAMZGA
      DOUBLE PRECISION TANH
      DOUBLE PRECISION CZ3L, CZ3H, CZ35, CW35
      DOUBLE PRECISION TAUT, LAMBDAT
      DOUBLE COMPLEX AT, AMP
      DOUBLE PRECISION SW, CW, PI
      DOUBLE PRECISION MTRUN, MBRUN, MCRUN, MT, MB, MC
C Functions to be called:
      DOUBLE PRECISION GAMFF, GAMWH, GAMZH
      DOUBLE COMPLEX FA12, I2
      DOUBLE PRECISION RUNMB, RUNMC, RUNALS, POLEMQ, RUNMT
      PI = 4.D0*DATAN(1.D0)
      CW = MW/MZ
      SW = DSQRT(1.D0 - CW**2)

      IF (QCDCORRS.EQ.1) THEN
C Running quark masses, for Yukawa couplings
         MTRUN = RUNMT(MH3)
         MBRUN = RUNMB(MH3)
         MCRUN = RUNMC(MH3)
C Quark pole masses, for kinematics
         MT = MTPOLE
         MB = POLEMQ(MBMB,5)
         MC = POLEMQ(MCMC,4)
      ELSE
         MTRUN = MTPOLE
         MBRUN = MBMB
         MCRUN = MCMC
         MT = MTPOLE
         MB = MBMB
         MC = MCMC
      ENDIF

C H3N decays to fermions: computing at tree level for now, on-shell for now.
C Dropping an overall i on the pseudoscalar couplings, because they are 
C mod-squared anyway in GAMFF.
      TANH = DSQRT(8.D0)*VCHI/VPHI
      H3NGAMB = GAMFF(3.D0,MH3,MB,MB,0.D0,MBRUN/V*TANH)
      H3NGAMTA = GAMFF(1.D0,MH3,MTAU,MTAU,0.D0,MTAU/V*TANH)
      H3NGAMMU = GAMFF(1.D0,MH3,MMU,MMU,0.D0,MMU/V*TANH)
      H3NGAMS = GAMFF(3.D0,MH3,MS,MS,0.D0,MS/V*TANH)
      H3NGAMC = GAMFF(3.D0,MH3,MC,MC,0.D0,-MCRUN/V*TANH)
      H3NGAMT = GAMFF(3.D0,MH3,MT,MT,0.D0,-MTRUN/V*TANH)

C H3N decays to vector + scalar (W+H- and W-H+ are summed): 
C On-shell above threshold, V off-shell below threshold.
C CHHV's are pure imaginary: we remove the i since it's mod-squared.
C H3NGAMWH5P is the sum of H3^0 -> W+ H5- and H3^0 -> W- H5+.
      CZ3L = 2.D0*MZ/V * (DSQRT(2.D0)*DCOS(ALPHA)*VCHI/V
     .     + DSQRT(2.D0/3.D0)*DSIN(ALPHA)*VPHI/V)
      H3NGAMZHL = GAMZH(MH3,MHL,MZ,CZ3L,V,SW)
      CZ3H = -2.D0*MZ/V * (DSQRT(2.D0/3.D0)*DCOS(ALPHA)*VPHI/V
     .     - DSQRT(2.D0)*DSIN(ALPHA)*VCHI/V)
      H3NGAMZHH = GAMZH(MH3,MHH,MZ,CZ3H,V,SW)
      CZ35 = 2.D0*MZ*VPHI/DSQRT(3.D0)/V**2
      H3NGAMZH5N = GAMZH(MH3,MH5,MZ,CZ35,V,SW)
      CW35 = -MW*VPHI/V**2
      H3NGAMWH5P = 2.D0 * GAMWH(MH3,MH5,MW,CW35,V)

C H3N decays to massless gauge bosons (loop-induced)
      TAUT = 4.D0*MT**2/MH3**2
      LAMBDAT = 4.D0*MT**2/MZ**2

C Loop-induced H3^0 -> gg:
C Top quark loop (note: be careful with signs when adding other fermions!)
      AT = -FA12(TAUT)
      AMP = TANH*AT
      IF (QCDCORRS.EQ.1) THEN
         H3NGAMG = RUNALS(MH3,5)**2 * MH3**3 / 128.D0/PI**3 / V**2
     .        * CDABS(AMP)**2
         RxH3NGG = RUNALS(MH3,5) / 4.D0/PI / V
     .        * CDABS(AMP)
      ELSE
         H3NGAMG = ALSMZ**2 * MH3**3 / 128.D0/PI**3 / V**2
     .        * CDABS(AMP)**2
         RxH3NGG = ALSMZ / 4.D0/PI / V
     .        * CDABS(AMP)
      ENDIF

C Loop-induced H3^0 -> gamma gamma:
C Top quark loop:
      AT = 3.D0 * (2.D0/3.D0)**2 * FA12(TAUT)
      AMP = TANH*AT
      H3NGAMGA = ALPHAEM**2 * MH3**3 / 256.D0/PI**3 / V**2
     .     * CDABS(AMP)**2
      RxH3NGAGA = ALPHAEM / 2.D0/PI / V
     .     * CDABS(AMP)
C Loop-induced H3^0 -> Z gamma:
      IF (MH3.LE.MZ) THEN
         H3NGAMZGA = 0.D0
         RxH3NZGA = 0.D0
      ELSE
C Top quark loop:
         AT = 3.D0 * (-4.D0/3.D0) * (0.5D0 - 4.D0/3.D0*SW**2)/SW/CW
     .        * (-I2(TAUT,LAMBDAT))
         AMP = TANH*AT
         H3NGAMZGA = ALPHAEM**2 * MH3**3 / 128.D0/PI**3 / V**2
     .        * CDABS(AMP)**2 * (1.D0 - MZ**2/MH3**2)**3
         RxH3NZGA = ALPHAEM / 2.D0/PI / V
     .        * CDABS(AMP)
      ENDIF

C Compute the total width
      H3NWDTH = H3NGAMB + H3NGAMTA + H3NGAMMU + H3NGAMS + H3NGAMC
     .     + H3NGAMT
     .     + H3NGAMZHL + H3NGAMZHH + H3NGAMZH5N + H3NGAMWH5P
     .     + H3NGAMG + H3NGAMGA + H3NGAMZGA

C Compute the BRs
      H3NBRB = H3NGAMB/H3NWDTH
      H3NBRTA = H3NGAMTA/H3NWDTH
      H3NBRMU = H3NGAMMU/H3NWDTH
      H3NBRS = H3NGAMS/H3NWDTH
      H3NBRC = H3NGAMC/H3NWDTH
      H3NBRT = H3NGAMT/H3NWDTH
      H3NBRZHL = H3NGAMZHL/H3NWDTH
      H3NBRZHH = H3NGAMZHH/H3NWDTH
      H3NBRZH5N = H3NGAMZH5N/H3NWDTH
      H3NBRWH5P = H3NGAMWH5P/H3NWDTH
      H3NBRG = H3NGAMG/H3NWDTH
      H3NBRGA = H3NGAMGA/H3NWDTH
      H3NBRZGA = H3NGAMZGA/H3NWDTH

      RETURN
      END

C------------------------------------------------------------------
C Subroutine to calculate the decays of H_3^+ = H3P

      SUBROUTINE H3PDECAYS
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
      INTEGER OFFSHELL, QCDCORRS
      COMMON/DECAYFLAGS/OFFSHELL, QCDCORRS
C Local variables:
      DOUBLE PRECISION H3PGAMBC, H3PGAMTA, H3PGAMMU, H3PGAMSU,
     .     H3PGAMCS, H3PGAMTB, H3PGAMBU,
     .     H3PGAMWHL, H3PGAMWHH, H3PGAMZH5P, H3PGAMWH5N, H3PGAMWH5PP,
     .     H3PGAMWGA
      DOUBLE PRECISION TANH
      DOUBLE PRECISION CW3L, CW3H, CZ35, CW35N, CW35PP
      DOUBLE PRECISION CW, SW
      DOUBLE PRECISION MTRUN, MBRUN, MCRUN, MT, MB, MC
C Functions to be called:
      DOUBLE PRECISION GAMFF, GAMWH, GAMZH
      DOUBLE PRECISION RUNMB, RUNMC, DELTAQCD, DELTATB, DELTATC, RUNALS,
     .     POLEMQ, RUNMT
      DOUBLE PRECISION HETLOOPH3WGA
      CW = MW/MZ
      SW = DSQRT(1.D0 - CW**2)

      IF (QCDCORRS.EQ.1) THEN
C Running quark masses, for Yukawa couplings
         MTRUN = RUNMT(MH3)
         MBRUN = RUNMB(MH3)
         MCRUN = RUNMC(MH3)
C Quark pole masses, for kinematics
         MT = MTPOLE
         MB = POLEMQ(MBMB,5)
         MC = POLEMQ(MCMC,4)
      ELSE
         MTRUN = MTPOLE
         MBRUN = MBMB
         MCRUN = MCMC
         MT = MTPOLE
         MB = MBMB
         MC = MCMC
      ENDIF

C H3P decays to fermions: computing at tree-level for now, on-shell for now.
C Treating V_cs = 1, V_tb = 1.
      TANH = DSQRT(8.D0)*VCHI/VPHI
      H3PGAMBC = VCB**2 * GAMFF(3.D0,MH3,MB,MC,
     .     -(MCRUN-MBRUN)/DSQRT(2.D0)/V*TANH,
     .     (MCRUN+MBRUN)/DSQRT(2.D0)/V*TANH)
      H3PGAMTA = GAMFF(1.D0,MH3,MTAU,0.D0,
     .     MTAU/DSQRT(2.D0)/V*TANH,MTAU/DSQRT(2.D0)/V*TANH)
      H3PGAMMU = GAMFF(1.D0,MH3,MMU,0.D0,
     .     MMU/DSQRT(2.D0)/V*TANH,MMU/DSQRT(2.D0)/V*TANH)
      H3PGAMSU = VUS**2 * GAMFF(3.D0,MH3,MS,0.D0,
     .     MS/DSQRT(2.D0)/V*TANH,MS/DSQRT(2.D0)/V*TANH)
      H3PGAMCS = GAMFF(3.D0,MH3,MS,MC,
     .     -(MCRUN-MS)/DSQRT(2.D0)/V*TANH,
     .     (MCRUN+MS)/DSQRT(2.D0)/V*TANH)
      H3PGAMTB = GAMFF(3.D0,MH3,MB,MT,
     .     -(MTRUN-MBRUN)/DSQRT(2.D0)/V*TANH,
     .     (MTRUN+MBRUN)/DSQRT(2.D0)/V*TANH)
      H3PGAMBU = VUB**2 * GAMFF(3.D0,MH3,MB,0.D0,
     .     MBRUN/DSQRT(2.D0)/V*TANH,MBRUN/DSQRT(2.D0)/V*TANH)

C H3P decays to vector + scalar: 
C On-shell above threshold, V off-shell below threshold.
      CW3L = 2.D0*MW/V**2 * (DSQRT(2.D0)*DCOS(ALPHA)*VCHI 
     .     + DSQRT(2.D0/3.D0)*DSIN(ALPHA)*VPHI)
      H3PGAMWHL = GAMWH(MH3,MHL,MW,CW3L,V)
      CW3H = 2.D0*MW/V**2 * (DSQRT(2.D0)*DSIN(ALPHA)*VCHI
     .     - DSQRT(2.D0/3.D0)*DCOS(ALPHA)*VPHI)
      H3PGAMWHH = GAMWH(MH3,MHH,MW,CW3H,V)
      CZ35 = -MZ*VPHI/V**2
      H3PGAMZH5P = GAMZH(MH3,MH5,MZ,CZ35,V,SW)
      CW35N = -MW*VPHI/DSQRT(3.D0)/V**2
      H3PGAMWH5N = GAMWH(MH3,MH5,MW,CW35N,V)
      CW35PP = -DSQRT(2.D0)*MW*VPHI/V**2
      H3PGAMWH5PP = GAMWH(MH3,MH5,MW,CW35PP,V)

C H3P decay to W+ gamma (loop induced)
      H3PGAMWGA = HETLOOPH3WGA()

C Compute the total width
      H3PWDTH = H3PGAMBC + H3PGAMTA + H3PGAMMU + H3PGAMSU
     .     + H3PGAMCS + H3PGAMTB + H3PGAMBU
     .     + H3PGAMWHL + H3PGAMWHH + H3PGAMZH5P 
     .     + H3PGAMWH5N + H3PGAMWH5PP
     .     + H3PGAMWGA

C Compute the BRs
      H3PBRBC = H3PGAMBC/H3PWDTH
      H3PBRTA = H3PGAMTA/H3PWDTH
      H3PBRMU = H3PGAMMU/H3PWDTH
      H3PBRSU = H3PGAMSU/H3PWDTH
      H3PBRCS = H3PGAMCS/H3PWDTH
      H3PBRTB = H3PGAMTB/H3PWDTH
      H3PBRBU = H3PGAMBU/H3PWDTH
      H3PBRWHL = H3PGAMWHL/H3PWDTH
      H3PBRWHH = H3PGAMWHH/H3PWDTH
      H3PBRZH5P = H3PGAMZH5P/H3PWDTH
      H3PBRWH5N = H3PGAMWH5N/H3PWDTH
      H3PBRWH5PP = H3PGAMWH5PP/H3PWDTH
      H3PBRWGA = H3PGAMWGA/H3PWDTH

      RETURN
      END

C------------------------------------------------------------------
C Subroutine to calculate the decays of H_5^0 = H5N (scalar)

      SUBROUTINE H5NDECAYS
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
      DOUBLE PRECISION H5NBRGA, H5NBRZGA, H5NBRW, H5NBRZ,
     .     H5NBRZH3N, H5NBRWH3P, 
     .     H5NBRH3N, H5NBRH3P,
     .     H5NWDTH
      COMMON/H5NBRS/H5NBRGA, H5NBRZGA, H5NBRW, H5NBRZ,
     .     H5NBRZH3N, H5NBRWH3P, 
     .     H5NBRH3N, H5NBRH3P,
     .     H5NWDTH
      INTEGER OFFSHELL, QCDCORRS
      COMMON/DECAYFLAGS/OFFSHELL, QCDCORRS
C Local variables:
      DOUBLE PRECISION H5NGAMGA, H5NGAMZGA, H5NGAMW, H5NGAMZ,
     .     H5NGAMZH3N, H5NGAMWH3P,
     .     H5NGAMH3N, H5NGAMH3P
      DOUBLE PRECISION KW, KZ
      DOUBLE PRECISION TAUW
      DOUBLE COMPLEX AW, AS, A3P, A5P, A5PP, AMP
      DOUBLE PRECISION SW, CW, PI
      DOUBLE PRECISION CZ53, CW53, G533N, G533P
      DOUBLE PRECISION G555P, G555PP, CZ11, CZ22
C Functions to be called:
      DOUBLE PRECISION GAMVVOF, GAMVV, GAMWH, GAMZH, GAMHH
      DOUBLE COMPLEX F0, F1, I1, I2
      DOUBLE PRECISION HETLOOPH5ZGA
      PI = 4.D0*DATAN(1.D0)
      CW = MW/MZ
      SW = DSQRT(1.D0 - CW**2)

C H5N decays to massive gauge bosons: doubly offshell unless switched off
      KW = 2.D0*DSQRT(2.D0/3.D0)*VCHI/V
      KZ = -4.D0*DSQRT(2.D0/3.D0)*VCHI/V
      IF (OFFSHELL.EQ.1) THEN
         H5NGAMW = GAMVVOF(1.D0,MH5,MW,MW,GAMW,GAMW,2.D0*MW**2/V*KW)
         H5NGAMZ = GAMVVOF(0.5D0,MH5,MZ,MZ,GAMZ,GAMZ,2.D0*MZ**2/V*KZ)
      ELSE
C onshell option (OFFSHELL = 0)
         H5NGAMW = GAMVV(1.D0,MH5,MW,MW,2.D0*MW**2/V*KW)
         H5NGAMZ = GAMVV(0.5D0,MH5,MZ,MZ,2.D0*MZ**2/V*KZ)
      ENDIF

C H5N decays to vector + scalar (W+H3- and W-H3+ are summed):
C On-shell above threshold, V off-shell below threshold.
C CZ53 is pure imaginary; remove the i (since it gets mod-squared).
C H5NGAMWH3P is the sum of H5^0 -> W+ H3- and H5^0 -> W- H3+.
      CZ53 = -MZ/V**2*VPHI*2.D0/DSQRT(3.D0)
      H5NGAMZH3N = GAMZH(MH5,MH3,MZ,CZ53,V,SW)
      CW53 = -MW/V**2*VPHI/DSQRT(3.D0)
      H5NGAMWH3P = 2.D0 * GAMWH(MH5,MH3,MW,CW53,V)

C H5N decays to two scalars: on-shell only.
      G533N = -2.D0*DSQRT(2.D0/3.D0)/V**2 
     .     * (-8.D0*LAMBDA5*VCHI**3 + 4.D0*M1*VCHI**2 
     .     + (-4.D0*LAMBDA5 + 2.D0*LAMBDA3)*VPHI**2*VCHI
     .     + 3.D0*M2*VPHI**2)
      G533P = DSQRT(2.D0/3.D0)/V**2 
     .     * (-8.D0*LAMBDA5*VCHI**3 + 4.D0*M1*VCHI**2 
     .     + (-4.D0*LAMBDA5 + 2.D0*LAMBDA3)*VPHI**2*VCHI
     .     + 3.D0*M2*VPHI**2)
      H5NGAMH3N = GAMHH(0.5D0,MH5,MH3,MH3,G533N)
      H5NGAMH3P = GAMHH(1.D0,MH5,MH3,MH3,G533P)

C H5N decays to massless gauge bosons (loop-induced)
      TAUW = 4.D0*MW**2/MH5**2
      G555P = DSQRT(6.D0)*(2.D0*LAMBDA3*VCHI - M2)
      G555PP = -2.D0*DSQRT(6.D0)*(2.D0*LAMBDA3*VCHI - M2)

C Loop-induced H5^0 -> gamma gamma:
C SM amplitude for the W loop:
      AW = F1(TAUW)
C Amplitudes for the scalar loops:
      A3P = G533P*V/2.D0/MH3**2 * F0(4.D0*MH3**2/MH5**2)
      A5P = G555P*V/2.D0/MH5**2 * F0(4.D0*MH5**2/MH5**2)
      A5PP = G555PP*V/2.D0/MH5**2 * 4.D0 * F0(4.D0*MH5**2/MH5**2)
      AS = A3P + A5P + A5PP
C Combine the amplitudes (all are complex in general):
      AMP = KW*AW + AS
      H5NGAMGA = ALPHAEM**2 * MH5**3 / 256.D0/PI**3 / V**2
     .     * CDABS(AMP)**2
      RxH5NGAGA = ALPHAEM / 2.D0/PI / V
     .     * CDABS(AMP)
C Loop-induced H5^0 -> Z gamma:
      H5NGAMZGA = HETLOOPH5ZGA()

C Compute the total width
      H5NWDTH = H5NGAMGA + H5NGAMZGA + H5NGAMW + H5NGAMZ 
     .     + H5NGAMZH3N + H5NGAMWH3P
     .     + H5NGAMH3N + H5NGAMH3P
C Compute the BRs
      H5NBRGA = H5NGAMGA/H5NWDTH
      H5NBRZGA = H5NGAMZGA/H5NWDTH
      H5NBRW = H5NGAMW/H5NWDTH
      H5NBRZ = H5NGAMZ/H5NWDTH
      H5NBRZH3N = H5NGAMZH3N/H5NWDTH
      H5NBRWH3P = H5NGAMWH3P/H5NWDTH
      H5NBRH3N = H5NGAMH3N/H5NWDTH
      H5NBRH3P = H5NGAMH3P/H5NWDTH

      RETURN
      END

C------------------------------------------------------------------
C Subroutine to calculate the decays of H_5^+ = H5P

      SUBROUTINE H5PDECAYS
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
      DOUBLE PRECISION H5PBRWZ, H5PBRZH3P, H5PBRWH3N, H5PBRH3PN, 
     .     H5PBRWGA, H5PWDTH
      COMMON/H5PBRS/H5PBRWZ, H5PBRZH3P, H5PBRWH3N, H5PBRH3PN, 
     .     H5PBRWGA, H5PWDTH
      INTEGER OFFSHELL, QCDCORRS
      COMMON/DECAYFLAGS/OFFSHELL, QCDCORRS
C Local variables:
      DOUBLE PRECISION H5PGAMWZ, H5PGAMZH3P, H5PGAMWH3N, H5PGAMH3PN,
     .     H5PGAMWGA
      DOUBLE PRECISION KV, SV
      DOUBLE PRECISION CW, SW
      DOUBLE PRECISION CZ35, CW35, G5P3P3N
C Functions to be called:
      DOUBLE PRECISION GAMVVOF, GAMVV, GAMWH, GAMZH, GAMHH
      DOUBLE PRECISION HETLOOPH5WGA
      CW = MW/MZ
      SW = DSQRT(1.D0 - CW**2)

C H5P decays to massive gauge bosons: doubly offshell unless switched off
C WZ non-identical bosons: symmetry factor SV = 1.
      KV = -2.D0*DSQRT(2.D0)*VCHI/V
      IF (OFFSHELL.EQ.1) THEN
         H5PGAMWZ = GAMVVOF(1.D0,MH5,MW,MZ,GAMW,GAMZ,2.D0*MW*MZ/V*KV)
      ELSE
         H5PGAMWZ = GAMVV(1.D0,MH5,MW,MZ,2.D0*MW*MZ/V*KV)
      ENDIF

C H5P decays to vector + scalar:
C On-shell above threshold, V off-shell below threshold.
C CW35 is pure imaginary; remove the i (since it gets mod-squared).
      CZ35 = -MZ/V**2*VPHI
      H5PGAMZH3P = GAMZH(MH5,MH3,MZ,CZ35,V,SW)
      CW35 = -MW/V**2*VPHI
      H5PGAMWH3N = GAMWH(MH5,MH3,MW,CW35,V)

C H5P decays to two scalars: on-shell only.
C G5P3P3N is pure imaginary; remove the i (since it gets mod-squared).
      G5P3P3N = -DSQRT(2.D0)/V**2 
     .     * (-8.D0*LAMBDA5*VCHI**3 + 4.D0*M1*VCHI**2
     .     + (-4.D0*LAMBDA5 + 2.D0*LAMBDA3)*VPHI**2*VCHI
     .     + 3.D0*M2*VPHI**2)
      H5PGAMH3PN = GAMHH(1.D0,MH5,MH3,MH3,G5P3P3N)

C H5P decay to W gamma:
      H5PGAMWGA = HETLOOPH5WGA()

C Compute the total width
      H5PWDTH = H5PGAMWZ + H5PGAMZH3P + H5PGAMWH3N + H5PGAMH3PN
     .     + H5PGAMWGA

C Compute the BRs
      H5PBRWZ = H5PGAMWZ/H5PWDTH
      H5PBRZH3P = H5PGAMZH3P/H5PWDTH
      H5PBRWH3N = H5PGAMWH3N/H5PWDTH
      H5PBRH3PN = H5PGAMH3PN/H5PWDTH
      H5PBRWGA = H5PGAMWGA/H5PWDTH

      RETURN
      END

C------------------------------------------------------------------
C Subroutine to calculate the decays of H_5^++ = H5PP

      SUBROUTINE H5PPDECAYS
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
      DOUBLE PRECISION H5PPBRWW, H5PPBRWH3, H5PPBRH3P, H5PPWDTH
      COMMON/H5PPBRS/H5PPBRWW, H5PPBRWH3, H5PPBRH3P, H5PPWDTH
      INTEGER OFFSHELL, QCDCORRS
      COMMON/DECAYFLAGS/OFFSHELL, QCDCORRS
C Local variables:
      DOUBLE PRECISION H5PPGAMWW, H5PPGAMWH3, H5PPGAMH3P
      DOUBLE PRECISION KW
      DOUBLE PRECISION CW53, G5PP3P3P
C Functions to be called:
      DOUBLE PRECISION GAMVVOF, GAMVV, GAMWH, GAMHH

C H5PP decays to massive gauge bosons: doubly offshell unless switched off
C Symmetry factor SV = 1/2 for identical W+W+
      KW = 4.D0*VCHI/V
      IF (OFFSHELL.EQ.1) THEN
         H5PPGAMWW = GAMVVOF(0.5D0,MH5,MW,MW,GAMW,GAMW,2.D0*MW**2/V*KW)
      ELSE
         H5PPGAMWW = GAMVV(0.5D0,MH5,MW,MW,2.D0*MW**2/V*KW)
      ENDIF

C H5PP decays to vector + scalar: 
C On-shell above threshold, V off-shell below threshold.
      CW53 = DSQRT(2.D0)*MW*VPHI/V**2
      H5PPGAMWH3 = GAMWH(MH5,MH3,MW,CW53,V)

C H5PP decays to two scalars: on-shell only.
      G5PP3P3P = -2.D0/V**2 
     .     * (-8.D0*LAMBDA5*VCHI**3 + 4.D0*M1*VCHI**2 
     .     + (-4.D0*LAMBDA5 + 2.D0*LAMBDA3)*VPHI**2*VCHI
     .     + 3.D0*M2*VPHI**2)
      H5PPGAMH3P = GAMHH(0.5D0,MH5,MH3,MH3,G5PP3P3P)

C Compute the total width
      H5PPWDTH = H5PPGAMWW + H5PPGAMWH3 + H5PPGAMH3P

C Compute the BRs
      H5PPBRWW = H5PPGAMWW/H5PPWDTH
      H5PPBRWH3 = H5PPGAMWH3/H5PPWDTH
      H5PPBRH3P = H5PPGAMH3P/H5PPWDTH

      RETURN
      END

C------------------------------------------------------------------
C Subroutine to calculate the decays of the top quark (can have t -> b H3+)

      SUBROUTINE TOPDECAYS
      IMPLICIT NONE
C Common blocks:
      DOUBLE PRECISION MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
      COMMON/PHYSPARAMS/MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
      DOUBLE PRECISION V, MZ, MW, MTPOLE, MBMB, MCMC, MS, MTAU, MMU,
     .     ALPHAEM,ALSMZ, VCB, VUS, VUB, GAMZ, GAMW
      COMMON/SM/V,MZ,MW,MTPOLE,MBMB,MCMC,MS,MTAU,MMU,
     .     ALPHAEM,ALSMZ,VCB,VUS,VUB,GAMZ,GAMW
      DOUBLE PRECISION TOPBRW, TOPBRH3P, TOPWDTH
      COMMON/TOPBRS/TOPBRW, TOPBRH3P, TOPWDTH
C External functions to be called:
      DOUBLE PRECISION RUNMB
C Local variables:
      DOUBLE PRECISION TOPGAMW, TOPGAMH3P
      DOUBLE PRECISION TANH, PI
      DOUBLE PRECISION MT, MB
      PI = 4.D0*DATAN(1.D0)
      TANH = 2.D0*DSQRT(2.D0)*VCHI/VPHI

C Should probably use the running top mass here too... need to improve 
C QCD corrections.
      MB = RUNMB(MH3)
      MT = MTPOLE

C These formulae neglect mb in the kinematics but keep it in the H3+ coupling.
C We're setting Vtb = 1.
      TOPGAMW = MT/16.D0/PI/V**2 * (MT**2 + 2.D0*MW**2)
     .     * (1 - MW**2/MT**2)**2
      IF (MT.GT.MH3) THEN
         TOPGAMH3P = MT/16.D0/PI/V**2 * TANH**2 * (MT**2 + MB**2)
     .        * (1 - MH3**2/MT**2)**2
      ELSE
         TOPGAMH3P = 0.D0
      ENDIF

      TOPWDTH = TOPGAMW + TOPGAMH3P

      TOPBRW = TOPGAMW/TOPWDTH
      TOPBRH3P = TOPGAMH3P/TOPWDTH

      RETURN
      END

C======================================================================
C Functions used by the decay table subroutines

      DOUBLE PRECISION FUNCTION GAMFF(NC,MH,MF1,MF2,CS,CP)
C Function for Gamma(H -> f1 f2)
C X1 = mf1/mH, X2 = mf2/mH
C NC = number of colors of fermions f1 and f2
C CS = scalar and CP = pseudoscalar parts of H coupling, 
C Feynman rule = -i (CS + CP gamma5)
      IMPLICIT NONE
      DOUBLE PRECISION NC, MH, MF1, MF2, CS, CP
      DOUBLE PRECISION X1, X2, LAMBDA, PI
      PI = 4.D0*DATAN(1.D0)
      X1 = MF1/MH
      X2 = MF2/MH
      IF ((X1+X2).GE.1.D0) THEN
         GAMFF = 0.D0
      ELSE
         GAMFF = NC*MH/8.D0/PI * DSQRT(LAMBDA(X1**2,X2**2))
     .        * ( CS**2 * (1.D0-(X1+X2)**2) 
     .        + CP**2 * (1.D0-(X1-X2)**2) )
      ENDIF
      END

      DOUBLE PRECISION FUNCTION LAMBDA(X,Y)
C Kinematic function lambda
      IMPLICIT NONE
      DOUBLE PRECISION X,Y
      LAMBDA = (1.D0 - X - Y)**2 - 4.D0*X*Y
      END


      DOUBLE PRECISION FUNCTION GAMVV(SV,MH,MV1,MV2,C)
C Function for Gamma(H -> V1 V2) above threshold (V1 & V2 on shell).
C Symmetry factor: SV = 1 for V1 != V2, SV = 1/2 for V1 = V2.
C C is the HVV coupling, Feynman rule = i C g_{munu}.
      IMPLICIT NONE
      DOUBLE PRECISION SV, MH, MV1, MV2, C
      DOUBLE PRECISION K1, K2, PI
      DOUBLE PRECISION LAMBDA
      PI = 4.D0*DATAN(1.D0)
      K1 = MV1**2/MH**2
      K2 = MV2**2/MH**2
      IF ((MV1+MV2).GE.MH) THEN
         GAMVV = 0.D0
      ELSE
         GAMVV = SV * C**2 * MH**3 / 64.D0/PI / MV1**2 / MV2**2 
     .        * (1.D0 - 2.D0*K1 - 2.D0*K2 + 10.D0*K1*K2 + K1**2 + K2**2) 
     .        * DSQRT(LAMBDA(K1,K2))
      ENDIF
      END

C Function GAMVVOF is in the file gamvvof_final.f.

C Function rho for change of variables in the Breit-Wigner integration
      DOUBLE PRECISION FUNCTION RHO(Q,M,GA)
      IMPLICIT NONE
      DOUBLE PRECISION Q, M, GA
      DOUBLE PRECISION PI
      PI = 4.D0*DATAN(1.D0)
      RHO = 1.D0/PI * DATAN((Q**2-M**2)/M/GA)
      END

C Inverse function for finding Q in terms of rho
      DOUBLE PRECISION FUNCTION QFROMRHO(RHO,M,GA)
      IMPLICIT NONE
      DOUBLE PRECISION RHO, M, GA
      DOUBLE PRECISION PI
      PI = 4.D0*DATAN(1.D0)
      QFROMRHO = DSQRT(M*GA*DTAN(PI*RHO) + M**2)
      END
C====================================================================

      DOUBLE PRECISION FUNCTION GAMWH(MH1,MH2,MW,C,V)
C Function for Gamma(H1 -> H2 W), for one choice of the sign of the W.
C Calls GAMVH for on-shell decays when mH1 > mH2 + MW.
C Calls GAMVSTARH for singly off-shell decays when mH1 < mH2 + MW.
      IMPLICIT NONE
      DOUBLE PRECISION MH1, MH2, MW, C, V
      DOUBLE PRECISION DV
      DOUBLE PRECISION GAMVH, GAMVSTARH
      INTEGER OFFSHELL, QCDCORRS
      COMMON/DECAYFLAGS/OFFSHELL, QCDCORRS
      IF (MH1.GT.MH2+MW) THEN
C Compute the H1 -> H2 W on-shell width.
         GAMWH = GAMVH(MH1,MH2,MW,C)
      ELSE IF (MH1.GT.MH2.AND.OFFSHELL.EQ.1) THEN
C Compute the H1 -> H2 W* singly off-shell width.  
C DV = 3/2 counts only one sign choice for the W*.
         DV = 3.D0/2.D0
         GAMWH = GAMVSTARH(DV,MH1,MH2,MW,C,V)
      ELSE
         GAMWH = 0.D0
      ENDIF
      END

      DOUBLE PRECISION FUNCTION GAMZH(MH1,MH2,MZ,C,V,SW)
C Function for Gamma(H1 -> H2 Z).
C Calls GAMVH for on-shell decays when mH1 > mH2 + MZ.
C Calls GAMVSTARH for singly off-shell decays when mH1 < mH2 + MZ.
      IMPLICIT NONE
      DOUBLE PRECISION MH1, MH2, MZ, C, V, SW
      DOUBLE PRECISION DV
      DOUBLE PRECISION GAMVH, GAMVSTARH
      INTEGER OFFSHELL, QCDCORRS
      COMMON/DECAYFLAGS/OFFSHELL, QCDCORRS
      IF (MH1.GT.MH2+MZ) THEN
C Compute the H1 -> H2 Z on-shell width.
         GAMZH = GAMVH(MH1,MH2,MZ,C)
      ELSE IF (MH1.GT.MH2.AND.OFFSHELL.EQ.1) THEN
C Compute the H1 -> H2 Z* singly off-shell width.
         DV = 3.D0*(7.D0/12.D0 - 10.D0/9.D0*SW**2 + 40.D0/27.D0*SW**4) 
         GAMZH = GAMVSTARH(DV,MH1,MH2,MZ,C,V)
      ELSE
         GAMZH = 0.D0
      ENDIF
      END

      DOUBLE PRECISION FUNCTION GAMVH(MH1,MH2,MV,C)
C Function for Gamma(H1 -> V H2) above threshold (V & H2 on shell).
C V denotes Z, W+, or W-: the decays H -> W+H- and H -> W-H+ are distinct.
C C is the V* H1 H2* coupling, Feynman rule = i C (p1-p2)_{mu}.
      IMPLICIT NONE
      DOUBLE PRECISION MH1, MH2, MV, C
      DOUBLE PRECISION X1, X2, Y1, Y2, LAMBDA, PI
      PI = 4.D0*DATAN(1.D0)
      X1 = MH1**2/MV**2
      X2 = MH2**2/MV**2
      Y1 = MV**2/MH1**2
      Y2 = MH2**2/MH1**2
      IF ((MH2 + MV).GE.MH1) THEN
         GAMVH = 0.D0
      ELSE
         GAMVH = C**2 * MV**2 / 16.D0/PI / MH1
     .        * LAMBDA(X1,X2) * DSQRT(LAMBDA(Y1,Y2))
      ENDIF
      END

      DOUBLE PRECISION FUNCTION GAMVSTARH(DV,MH1,MH2,MV,C,V)
C Function for Gamma(H1 -> H2 V*) below threshold (singly off-shell).
C C is the V* H1 H2* coupling, Feynman rule = i C (p1-p2)_{mu}.
C Formula is taken from Djouadi et al, hep-ph/9511342, Eq. (43-46).
C We correct the typing error in the kinematic function G_ij pointed
C out in Akeroyd hep-ph/9806337 (last term should be +2 lambdaij/kj 
C instead of -2 lambdaij/kj).
C Note that this function does not work when the 2-body decay V -> H1 H2 
C is possible: in that case we set the H1 -> V* H2 width to zero.
      IMPLICIT NONE
      DOUBLE PRECISION DV, MH1, MH2, MV, C, V
      DOUBLE PRECISION COEFF, PI
      DOUBLE PRECISION KINEMATIC, KI, KJ, LAMIJ, ATANARG
      IF (MV.GT.(MH1+MH2)) THEN
         GAMVSTARH = 0.D0
      ELSE
         PI = 4.D0*DATAN(1.D0)
         COEFF = 3.D0/8.D0/PI**3 * C**2 * MV**2/2.D0/V**2 * MH1 * DV
         KI = MH2**2/MH1**2
         KJ = MV**2/MH1**2
         LAMIJ = -1.D0 + 2.D0*KI + 2.D0*KJ - (KI - KJ)**2
         ATANARG = (KJ*(1.D0-KJ+KI) - LAMIJ)/(1.D0-KI)/DSQRT(LAMIJ)
         KINEMATIC = 0.25D0 * ( 2.D0*(-1.D0+KJ-KI) * DSQRT(LAMIJ)
     .        * (PI/2.D0 + DATAN(ATANARG))
     .        + (LAMIJ - 2.D0*KI)*DLOG(KI) 
     .        + (1.D0-KI)/3.D0*(5.D0*(1.D0+KI) 
     .        - 4.D0*KJ + 2.D0/KJ*LAMIJ))
         GAMVSTARH = COEFF * KINEMATIC
      ENDIF
      END

      DOUBLE PRECISION FUNCTION GAMHH(SH,MH1,MH2,MH3,G123)
C Function for Gamma(H1 -> H2 H3)
C SH is a symmetry factor: SH = 1/2 for H2 = H3; 
C    SH = 1 for H2 != H3, including when H2/3 are h.c. of each other.
C G123 is the H1 H2* H3* coupling, Feynman rule = -i G123.
      IMPLICIT NONE
      DOUBLE PRECISION SH, MH1, MH2, MH3, G123
      DOUBLE PRECISION X2, X3, LAMBDA, PI
      PI = 4.D0*DATAN(1.D0)
      X2 = MH2**2/MH1**2
      X3 = MH3**2/MH1**2
      IF ((MH2 + MH3).GE.MH1) THEN
         GAMHH = 0.D0
      ELSE
         GAMHH = SH * G123**2 / 16.D0/PI / MH1
     .        * DSQRT(LAMBDA(X2,X3))
      ENDIF
      END


C==================================================================
C     Subroutines to calculate the couplings of h and H (kappas)
C==================================================================

      SUBROUTINE HLCOUPS
C Computes kappa_V, kappa_f, kappa_gamma, kappa_Zgamma,
C Delta kappa_gamma, and Delta kappa_Zgamma for the lighter 
C Higgs mass eigenstate h (of mass MHL).
C INPUTS: common blocks LPARAMS, PHYSPARAMS
C OUTPUTS: common block KAPPASL
      IMPLICIT NONE
C Common blocks:
      DOUBLE PRECISION MU2SQ,MU3SQ,LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,
     .     LAMBDA5,M1,M2
      COMMON/LPARAMS/MU2SQ,MU3SQ,LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,
     .     LAMBDA5,M1,M2
      DOUBLE PRECISION MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
      COMMON/PHYSPARAMS/MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
      DOUBLE PRECISION KVL, KFL, KGAML, KZGAML, DKGAML, DKZGAML
      COMMON/KAPPASL/KVL,KFL,KGAML,KZGAML,DKGAML,DKZGAML
      DOUBLE PRECISION V, MZ, MW, MTPOLE, MBMB, MCMC, MS, MTAU, MMU,
     .     ALPHAEM,ALSMZ, VCB, VUS, VUB, GAMZ, GAMW
      COMMON/SM/V,MZ,MW,MTPOLE,MBMB,MCMC,MS,MTAU,MMU,
     .     ALPHAEM,ALSMZ,VCB,VUS,VUB,GAMZ,GAMW
C Local variables:
      DOUBLE PRECISION VSM
      DOUBLE COMPLEX F0, F12, F1
      DOUBLE COMPLEX AT, AW, AS, A3P, A5P, A5PP
      DOUBLE PRECISION TAUT, TAUW
      DOUBLE PRECISION GHL33, GHL55, CZ11, CZ22
      DOUBLE COMPLEX I1, I2
      DOUBLE COMPLEX ATZ, AWZ, ASZ, A3PZ, A5PZ, A5PPZ
      DOUBLE PRECISION LAMBDAT, LAMBDAW
      DOUBLE PRECISION SW, CW
      DOUBLE PRECISION MT
      CW = MW/MZ
      SW = DSQRT(1.D0 - CW**2)

      VSM = DSQRT(VPHI**2 + 8.D0*VCHI**2)
      KVL = DCOS(ALPHA)*VPHI/VSM 
     .     - 8.D0/DSQRT(3.D0)*DSIN(ALPHA)*VCHI/VSM
      KFL = DCOS(ALPHA)*VSM/VPHI

      MT = MTPOLE

C Loop-induced h -> gamma gamma coupling:
C SM amplitude for the top loop:
      TAUT = 4.D0*MT**2/MHL**2
      AT = 3.D0 * (2.D0/3.D0)**2 * F12(TAUT)
C SM amplitude for the W loop:
      TAUW = 4.D0*MW**2/MHL**2
      AW = F1(TAUW)
C Amplitudes for the scalar loops:
      GHL33 = 64.D0*LAMBDA1*DCOS(ALPHA)*VCHI**2*VPHI/VSM**2
     .     - 8.D0/DSQRT(3.D0)*VPHI**2*VCHI/VSM**2
     .     *DSIN(ALPHA)*(LAMBDA3+3.D0*LAMBDA4)
     .     - 4.D0/DSQRT(3.D0)*VCHI*M1/VSM**2
     .     *(DSIN(ALPHA)*VCHI - DSQRT(3.D0)*DCOS(ALPHA)*VPHI)
     .     - 16.D0/DSQRT(3.D0)*VCHI**3/VSM**2*DSIN(ALPHA)
     .     *(6.D0*LAMBDA2+LAMBDA5)
     .     - DCOS(ALPHA)*VPHI**3/VSM**2*(LAMBDA5-4.D0*LAMBDA2)
     .     + 2.D0*DSQRT(3.D0)*M2*VPHI**2/VSM**2*DSIN(ALPHA)
     .     - 8.D0/DSQRT(3.D0)*LAMBDA5*VCHI*VPHI/VSM**2
     .     *(DSIN(ALPHA)*VPHI - DSQRT(3.D0)*DCOS(ALPHA)*VCHI)
      GHL55 = -8.D0*DSQRT(3.D0)*(LAMBDA3+LAMBDA4)*DSIN(ALPHA)*VCHI
     .     + (4.D0*LAMBDA2+LAMBDA5)*DCOS(ALPHA)*VPHI
     .     - 2.D0*DSQRT(3.D0)*M2*DSIN(ALPHA)
      A3P = GHL33*VSM/2.D0/MH3**2 * F0(4.D0*MH3**2/MHL**2)
      A5P = GHL55*VSM/2.D0/MH5**2 * F0(4.D0*MH5**2/MHL**2)
      A5PP = GHL55*VSM/2.D0/MH5**2 * 4.D0 * F0(4.D0*MH5**2/MHL**2)
      AS = A3P + A5P + A5PP
C Combine the amplitudes:
      DKGAML = DREAL(AS/(AT + AW))
      KGAML = DREAL((KFL*AT + KVL*AW + AS)/(AT + AW))

C Loop-induced h -> Z gamma coupling:
      IF (MHL.LT.MZ) THEN
         KZGAML = 1.D0
         DKZGAML = 0.D0
      ELSE
C SM amplitude for the top loop:
      LAMBDAT = 4.D0*MT**2/MZ**2
      ATZ = 3.D0 * (-4.D0/3.D0) * (0.5D0 - 4.D0/3.D0*SW**2)/SW/CW
     .     * (I1(TAUT,LAMBDAT) - I2(TAUT,LAMBDAT))
C SM amplitude for the W loop:
      LAMBDAW = 4.D0*MW**2/MZ**2
      AWZ = -CW/SW * (4.D0*(3.D0-SW**2/CW**2) * I2(TAUW,LAMBDAW)
     .     + ((1.D0+2.D0/TAUW)*SW**2/CW**2 - (5.D0+2.D0/TAUW)) 
     .     * I1(TAUW,LAMBDAW))
C Amplitudes for the scalar loops:
      CZ11 = 1.D0/2.D0/SW/CW * (1.D0 - 2.D0*SW**2)
      CZ22 = 1.D0/SW/CW * (1.D0 - 2.D0*SW**2)
      A3PZ = 2.D0*GHL33*CZ11/MH3**2 
     .     * I1(4.D0*MH3**2/MHL**2,4.D0*MH3**2/MZ**2)
      A5PZ = 2.D0*GHL55*CZ11/MH5**2
     .     * I1(4.D0*MH5**2/MHL**2,4.D0*MH5**2/MZ**2)
      A5PPZ = 2.D0*GHL55*CZ22*2.D0/MH5**2
     .     * I1(4.D0*MH5**2/MHL**2,4.D0*MH5**2/MZ**2)
      ASZ = A3PZ + A5PZ + A5PPZ
C Combine the amplitudes:
      DKZGAML = DREAL(VSM*ASZ/2.D0/(ATZ + AWZ))
      KZGAML = DREAL((2.D0*KFL*ATZ + 2.D0*KVL*AWZ + VSM*ASZ)
     .     /2.D0/(ATZ + AWZ))
      ENDIF

      RETURN
      END

C------------------------------------------------------------------
      SUBROUTINE HHCOUPS
C Computes kappa_V, kappa_f, kappa_gamma, kappa_Zgamma,
C Delta kappa_gamma, and Delta kappa_Zgamma for the heavier
C Higgs mass eigenstate H (of mass MHH).
C INPUTS: common blocks LPARAMS, PHYSPARAMS
C OUTPUTS: common block KAPPASH
      IMPLICIT NONE
C Common blocks:
      DOUBLE PRECISION MU2SQ,MU3SQ,LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,
     .     LAMBDA5,M1,M2
      COMMON/LPARAMS/MU2SQ,MU3SQ,LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,
     .     LAMBDA5,M1,M2
      DOUBLE PRECISION MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
      COMMON/PHYSPARAMS/MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
      DOUBLE PRECISION KVH, KFH, KGAMH, KZGAMH, DKGAMH, DKZGAMH
      COMMON/KAPPASH/KVH,KFH,KGAMH,KZGAMH,DKGAMH,DKZGAMH
      DOUBLE PRECISION V, MZ, MW, MTPOLE, MBMB, MCMC, MS, MTAU, MMU,
     .     ALPHAEM,ALSMZ, VCB, VUS, VUB, GAMZ, GAMW
      COMMON/SM/V,MZ,MW,MTPOLE,MBMB,MCMC,MS,MTAU,MMU,
     .     ALPHAEM,ALSMZ,VCB,VUS,VUB,GAMZ,GAMW
C Local variables:
      DOUBLE PRECISION VSM
      DOUBLE COMPLEX F0, F12, F1
      DOUBLE COMPLEX AT, AW, AS, A3P, A5P, A5PP
      DOUBLE PRECISION TAUT, TAUW
      DOUBLE PRECISION GHH33, GHH55, CZ11, CZ22
      DOUBLE COMPLEX I1, I2
      DOUBLE COMPLEX ATZ, AWZ, ASZ, A3PZ, A5PZ, A5PPZ
      DOUBLE PRECISION LAMBDAT, LAMBDAW
      DOUBLE PRECISION SW, CW
      DOUBLE PRECISION MT
      CW = MW/MZ
      SW = DSQRT(1.D0 - CW**2)

      VSM = DSQRT(VPHI**2 + 8.D0*VCHI**2)
      KVH = DSIN(ALPHA)*VPHI/VSM 
     .     + 8.D0/DSQRT(3.D0)*DCOS(ALPHA)*VCHI/VSM
      KFH = DSIN(ALPHA)*VSM/VPHI

      MT = MTPOLE
      
C Loop-induced H -> gamma gamma coupling:
C SM amplitude for the top loop:
      TAUT = 4.D0*MT**2/MHH**2
      AT = 3.D0 * (2.D0/3.D0)**2 * F12(TAUT)
C SM amplitude for the W loop:
      TAUW = 4.D0*MW**2/MHH**2
      AW = F1(TAUW)
C Amplitudes for the scalar loops:
      GHH33 = 64.D0*LAMBDA1*DSIN(ALPHA)*VCHI**2*VPHI/VSM**2
     .     + 8.D0/DSQRT(3.D0)*VPHI**2*VCHI/VSM**2 
     .     *DCOS(ALPHA)*(LAMBDA3 + 3.D0*LAMBDA4)
     .     + 4.D0/DSQRT(3.D0)*VCHI*M1/VSM**2
     .     *(DCOS(ALPHA)*VCHI + DSQRT(3.D0)*DSIN(ALPHA)*VPHI)
     .     + 16.D0/DSQRT(3.D0)*VCHI**3/VSM**2*DCOS(ALPHA)
     .     *(6.D0*LAMBDA2 + LAMBDA5)
     .     + DSIN(ALPHA)*VPHI**3/VSM**2*(4.D0*LAMBDA2 - LAMBDA5)
     .     - 2.D0*DSQRT(3.D0)*M2*VPHI**2/VSM**2*DCOS(ALPHA)
     .     + 8.D0/DSQRT(3.D0)*LAMBDA5*VCHI*VPHI/VSM**2
     .     * (DCOS(ALPHA)*VPHI + DSQRT(3.D0)*DSIN(ALPHA)*VCHI)
      GHH55 = 8.D0*DSQRT(3.D0)*(LAMBDA3+LAMBDA4)*DCOS(ALPHA)*VCHI
     .     + (4.D0*LAMBDA2+LAMBDA5)*DSIN(ALPHA)*VPHI
     .     + 2.D0*DSQRT(3.D0)*M2*DCOS(ALPHA)
      A3P = GHH33*VSM/2.D0/MH3**2 * F0(4.D0*MH3**2/MHH**2)
      A5P = GHH55*VSM/2.D0/MH5**2 * F0(4.D0*MH5**2/MHH**2)
      A5PP = GHH55*VSM/2.D0/MH5**2 * 4.D0 * F0(4.D0*MH5**2/MHH**2)
      AS = A3P + A5P + A5PP
C Combine the amplitudes:
      DKGAMH = DREAL(AS/(AT + AW))
      KGAMH = DREAL((KFH*AT + KVH*AW + AS)/(AT + AW))

C Loop-induced H -> Z gamma coupling:
      IF (MHH.LT.MZ) THEN
         KZGAMH = 1.D0
         DKZGAMH = 0.D0
      ELSE
C SM amplitude for the top loop:
      LAMBDAT = 4.D0*MT**2/MZ**2
      ATZ = 3.D0 * (-4.D0/3.D0) * (0.5D0 - 4.D0/3.D0*SW**2)/SW/CW
     .     * (I1(TAUT,LAMBDAT) - I2(TAUT,LAMBDAT))
C SM amplitude for the W loop:
      LAMBDAW = 4.D0*MW**2/MZ**2
      AWZ = -CW/SW * (4.D0*(3.D0-SW**2/CW**2) * I2(TAUW,LAMBDAW)
     .     + ((1.D0+2.D0/TAUW)*SW**2/CW**2 - (5.D0+2.D0/TAUW)) 
     .     * I1(TAUW,LAMBDAW))
C Amplitudes for the scalar loops:
      CZ11 = 1.D0/2.D0/SW/CW * (1.D0 - 2.D0*SW**2)
      CZ22 = 1.D0/SW/CW * (1.D0 - 2.D0*SW**2)
      A3PZ = 2.D0*GHH33*CZ11/MH3**2 
     .     * I1(4.D0*MH3**2/MHH**2,4.D0*MH3**2/MZ**2)
      A5PZ = 2.D0*GHH55*CZ11/MH5**2
     .     * I1(4.D0*MH5**2/MHH**2,4.D0*MH5**2/MZ**2)
      A5PPZ = 2.D0*GHH55*CZ22*2.D0/MH5**2
     .     * I1(4.D0*MH5**2/MHH**2,4.D0*MH5**2/MZ**2)
      ASZ = A3PZ + A5PZ + A5PPZ
C Combine the amplitudes:
      DKZGAMH = DREAL(VSM*ASZ/2.D0/(ATZ + AWZ))
      KZGAMH = DREAL((2.D0*KFH*ATZ + 2.D0*KVH*AWZ + VSM*ASZ)
     .     /2.D0/(ATZ + AWZ))
      ENDIF

      RETURN
      END


! modified by Ameen on 2018-08-02
C------------------------------------------------------------------
      SUBROUTINE H3COUPS
C Computes kappa_W, kappa_Z, and kappa_f
C for the triplet Higgs mass eigenstate H^3_0 (of mass MH3).
C Note that kappa_f as calculated here has the correct sign for
C down-type quarks. Multiply by -1 for up-type quarks.
C INPUTS: common blocks LPARAMS, PHYSPARAMS
C OUTPUTS: common block KAPPAS3
      IMPLICIT NONE
C Common blocks:
      DOUBLE PRECISION MU2SQ,MU3SQ,LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,
     .     LAMBDA5,M1,M2
      COMMON/LPARAMS/MU2SQ,MU3SQ,LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,
     .     LAMBDA5,M1,M2
      DOUBLE PRECISION MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
      COMMON/PHYSPARAMS/MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
      DOUBLE PRECISION KW3,KZ3,KF3
      COMMON/KAPPAS3/KW3,KZ3,KF3
C Local variables:
      DOUBLE PRECISION VSM

      VSM = DSQRT(VPHI**2 + 8.D0*VCHI**2)
      KW3 = 0.D0
      KZ3 = 0.D0
      KF3 = 2.D0*DSQRT(2.D0)*VCHI/VPHI 

      RETURN
      END


C------------------------------------------------------------------
      SUBROUTINE H5COUPS
C Computes kappa_W, kappa_Z, and kappa_f
C for the fiveplet Higgs mass eigenstate H^5_0 (of mass MH5).
C INPUTS: common blocks LPARAMS, PHYSPARAMS
C OUTPUTS: common block KAPPAS5
      IMPLICIT NONE
C Common blocks:
      DOUBLE PRECISION MU2SQ,MU3SQ,LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,
     .     LAMBDA5,M1,M2
      COMMON/LPARAMS/MU2SQ,MU3SQ,LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,
     .     LAMBDA5,M1,M2
      DOUBLE PRECISION MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
      COMMON/PHYSPARAMS/MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
      DOUBLE PRECISION KW5,KZ5,KF5
      COMMON/KAPPAS5/KW5,KZ5,KF5
C Local variables:
      DOUBLE PRECISION VSM

      VSM = DSQRT(VPHI**2 + 8.D0*VCHI**2)
      KW5 = 2.D0*DSQRT(2.D0)/DSQRT(3.D0)*VCHI/VSM
      KZ5 = -2.D0*KW5
      KF5 = 0.D0

      RETURN
      END
! end modified by Ameen


C Functions for loop integrals used by HLCOUPS and HHCOUPS:

      DOUBLE COMPLEX FUNCTION F0(TAU)
C This is the scalar loop function F_0 from the Higgs Hunter's Guide
C discussion of H -> gamma gamma.  Its value asymptotes to -1/3 
C in the limit of tau >> 1.  TAU should be 4*ms^2/mH^2.
      IMPLICIT NONE
      DOUBLE PRECISION TAU
      DOUBLE COMPLEX FF
      F0 = TAU*(1.D0 - TAU*FF(TAU))
      END

      DOUBLE COMPLEX FUNCTION F12(TAU)
C This is the fermion loop function F_{1/2} from the Higgs Hunter's
C Guide discussion of H -> gamma gamma.  Its vaule asymptotes to -4/3
C in the limit of tau >> 1.  TAU should be 4*mf^2/mH^2.
      IMPLICIT NONE
      DOUBLE PRECISION TAU
      DOUBLE COMPLEX FF
      F12 = -2.D0*TAU*(1.D0 + (1.D0-TAU)*FF(TAU))
      END

      DOUBLE COMPLEX FUNCTION FA12(TAU)
C This is the fermion loop function F^A_{1/2} for the pseudoscalar
C from the Higgs Hunter's Guide Appendix C.  
C TAU should be 4*mf^2/mA^2.
      IMPLICIT NONE
      DOUBLE PRECISION TAU
      DOUBLE COMPLEX FF
      FA12 = -2.D0*TAU*FF(TAU)
      END

      DOUBLE COMPLEX FUNCTION F1(TAU)
C This is the spin-1 boson loop function F_1 from the Higgs Hunter's 
C Guide discussion of H -> gamma gamma.  Its value asymptotes to 7 when
C tau >> 1.  TAU should be 4*MW^2/MH^2.
      IMPLICIT NONE
      DOUBLE PRECISION TAU
      DOUBLE COMPLEX FF
      F1 = 2.D0 + 3.D0*TAU + 3.D0*TAU*(2.D0-TAU)*FF(TAU)
      END

      DOUBLE COMPLEX FUNCTION FF(TAU)
C This is the function f(tau) from the Higgs Hunter's Guide discussion
C of H -> gamma gamma.
      IMPLICIT NONE
      DOUBLE PRECISION TAU
      DOUBLE PRECISION PI
      DOUBLE PRECISION ETAP, ETAM
      DOUBLE COMPLEX I
      IF (TAU.GE.1.D0) THEN
         FF = (DASIN(1.D0/DSQRT(TAU)))**2
      ELSE
         PI = 4.D0*DATAN(1.D0)
         ETAP = 1.D0 + DSQRT(1.D0 - TAU)
         ETAM = 1.D0 - DSQRT(1.D0 - TAU)
         I = (0.D0,1.D0)
         FF = -0.25D0*(DLOG(ETAP/ETAM) - I*PI)**2
      ENDIF
      END

      DOUBLE COMPLEX FUNCTION GG(TAU)
C This is the function g(tau) used in the calculation of h -> Z gamma.
      IMPLICIT NONE
      DOUBLE PRECISION TAU
      DOUBLE PRECISION PI
      DOUBLE PRECISION ETAP, ETAM
      DOUBLE COMPLEX I
      IF (TAU.GE.1.D0) THEN
         GG = DSQRT(TAU - 1.D0) * DASIN(1.D0/DSQRT(TAU))
      ELSE
         PI = 4.D0*DATAN(1.D0)
         ETAP = 1.D0 + DSQRT(1.D0 - TAU)
         ETAM = 1.D0 - DSQRT(1.D0 - TAU)
         I = (0.D0,1.D0)
         GG = 0.5D0 * DSQRT(1.D0 - TAU) * (DLOG(ETAP/ETAM) - I*PI)
      ENDIF
      END

      DOUBLE COMPLEX FUNCTION I1(A,B)
C This is the function I_2(a,b) used in the calculation of scalar, top quark,
C and W boson contributions to h -> Z gamma.
      IMPLICIT NONE
      DOUBLE PRECISION A, B
      DOUBLE COMPLEX FF, GG
      I1 = A*B/2.D0/(A-B) + A**2*B**2/2.D0/(A-B)**2 * (FF(A) - FF(B))
     .     + A**2*B/(A-B)**2 * (GG(A) - GG(B))
      END

      DOUBLE COMPLEX FUNCTION I2(A,B)
C This is the function I_2(a,b) used in the calculation of top and W boson
C contributions to h -> Z gamma.
      IMPLICIT NONE
      DOUBLE PRECISION A, B
      DOUBLE COMPLEX FF
      I2 = -A*B/2.D0/(A-B) * (FF(A) - FF(B))
      END

