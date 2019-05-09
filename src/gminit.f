C==================================================================
C     Subroutine to initialize the SM parameters common block
C==================================================================

      SUBROUTINE INITIALIZE_SM
      IMPLICIT NONE
      DOUBLE PRECISION V, MZ, MW, MTPOLE, MBMB, MCMC, MS, MTAU, MMU, 
     .     ALPHAEM,ALSMZ, VCB, VUS, VUB, GAMZ, GAMW
      COMMON/SM/V,MZ,MW,MTPOLE,MBMB,MCMC,MS,MTAU,MMU,
     .     ALPHAEM,ALSMZ,VCB,VUS,VUB,GAMZ,GAMW
      DOUBLE PRECISION PI, CW, SW, GG, EE
      DOUBLE PRECISION GF, CABIBBO, MDPOLE, MUPOLE, MELEC
      COMMON/MG5SM/GF,CABIBBO,MDPOLE,MUPOLE,MELEC
      PI = 4.D0*DATAN(1.D0)
      GF = 1.1663787D-5
      MZ = 91.1876D0
      MW = 80.385D0
      GAMZ = 2.4952D0
      GAMW = 2.085D0
c      MT = 172.D0
c      MB = 4.49D0
c      MC = 1.42D0
c      MS = 0.100D0
c      MTAU = 1.77682D0
c      MMU = 0.105658367D0
c      ALPHAEM = 1.D0/137.036D0
c      ALPHAS = 0.119D0
c      VCB = 0.0410D0
c      VUS = 0.2253D0
c      VUB = 0.0846D0*VCB
      MTPOLE = 1.725D2
      MBMB = 4.18D0
      MCMC = 1.42D0
      MS = 0.D0
      MTAU = 1.77682D0
      MMU = 1.056583715D-1
c      ALPHAEM = 1.D0/1.279D2
      ALSMZ = 1.18D-1
      CABIBBO = 2.277360D-1
      VCB = 0.D0
      VUB = 0.D0
C Derived parameters:
      V = DSQRT(1.D0/DSQRT(2.D0)/GF)
      VUS = DSIN(CABIBBO)
      GG = 2.D0*MW/V
      CW = MW/MZ
      SW = DSQRT(1.D0 - CW**2)
      EE = GG*SW
      ALPHAEM = EE**2/4.D0/PI
c      MW = DSQRT(MZ**2/2.D0 + DSQRT(MZ**4/4.D0 
c     .     - PI*ALPHAEM*MZ**2/DSQRT(2.D0)/GF))
C More pole masses, for MG5
      MDPOLE = 5.04D-3
      MUPOLE = 2.55D-3
      MELEC = 5.11D-4

      CALL INIT_QCD
      RETURN
      END

