C Subroutines to initialize and end LoopTools
C (These are defined to facilitate switching off the LoopTools use by instead
C linking heteroloops-dummy.F.)

      SUBROUTINE LTSTARTER
      IMPLICIT NONE
#include "looptools.h"
      CALL LTINI
      PRINT *
      RETURN
      END

      SUBROUTINE LTENDER
      IMPLICIT NONE
#include "looptools.h"      
      CALL LTEXI
      RETURN
      END

C Calculations of the one-loop decay widths:
C=====================================================================
C H5+ --> W+ gamma
C=====================================================================
      DOUBLE PRECISION FUNCTION HETLOOPH5WGA()
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
C Functions to be called:
      DOUBLE COMPLEX LOOPSSS, LOOPVSS, LOOPSVV, LOOPVVV,
     .     LOOPFFF, TLOOPFFF
C Local variables:
      DOUBLE COMPLEX S
      DOUBLE PRECISION GH3NH3PH5M, GH5PH5MH5N, GH5PH5PH5MM
      DOUBLE PRECISION L1, L2, L3, L4, L5
      DOUBLE PRECISION SW, CW
      DOUBLE PRECISION PI
      PI = 4.D0*DATAN(1.D0)

      CW = MW/MZ
      SW = DSQRT(1.D0 - CW**2)

      L1 = LAMBDA1
      L2 = LAMBDA2
      L3 = LAMBDA3
      L4 = LAMBDA4
      L5 = LAMBDA5

      IF (MH5.LE.MW) THEN
         HETLOOPH5WGA = 0.D0
         RxH5PWGA = 0.D0
         IxH5PWGA = 0.D0
      ELSE

         S = (0.D0,0.D0)

C This first coupling has an i removed:
      GH3NH3PH5M = -DSQRT(2.D0)/V**2*(-8.D0*L5*VCHI**3 + 4.D0*M1*VCHI**2
     .     + (-4.D0*L5 + 2.D0*L3)*VPHI**2*VCHI + 3.D0*M2*VPHI**2)
      GH5PH5MH5N = DSQRT(6.D0)*(2.D0*L3*VCHI - M2)
      GH5PH5PH5MM = -6.D0*(2.D0*L3*VCHI - M2)

C H30 H3- H3- loop:
      S = S - ALPHAEM*(-1.D0)/PI*(1.D0/2.D0/SW)*GH3NH3PH5M
     .     * LOOPSSS(MW**2,MH5**2,MH3**2,MH3**2)
C H50 H5- H5- loop:
      S = S - ALPHAEM*(-1.D0)/PI*(DSQRT(3.D0)/2.D0/SW)*GH5PH5MH5N
     .     * LOOPSSS(MW**2,MH5**2,MH5**2,MH5**2)
C H5- H5-- H5-- loop:
      S = S - ALPHAEM*(-2.D0)/PI*(-1.D0/DSQRT(2.D0)/SW)*GH5PH5PH5MM
     .     * LOOPSSS(MW**2,MH5**2,MH5**2,MH5**2)
C H5++ H5+ H5+ loop:
      S = S - ALPHAEM*(1.D0)/PI*(1.D0/DSQRT(2.D0)/SW)*GH5PH5PH5MM
     .     * LOOPSSS(MW**2,MH5**2,MH5**2,MH5**2)

C Z H5- H5- loop:
      S = S - ALPHAEM**2*(-1.D0)*(1.D0-2.D0*SW**2)/CW/SW
     .     * (-DSQRT(2.D0)*VCHI/2.D0/CW/SW**2) / MZ**2
     .     * LOOPVSS(MW**2,MH5**2,MZ**2,MH5**2)
C W- H5-- H5-- loop:
      S = S - ALPHAEM**2*(-2.D0)*(1.D0/DSQRT(2.D0)/SW)
     .     * 2.D0*VCHI/SW**2 / MW**2
     .     * LOOPVSS(MW**2,MH5**2,MW**2,MH5**2)

C H50 W- W- loop:
      S = S + ALPHAEM**2*(-1.D0)*(-DSQRT(3.D0)/2.D0/SW)
     .     * DSQRT(2.D0/3.D0)*VCHI/SW**2 / 2.D0/MW**2
     .     * LOOPSVV(MW**2,MH5**2,MH5**2,MW**2)
C H5++ W+ W+ loop:
      S = S + ALPHAEM**2*(1.D0)*(1.D0/DSQRT(2.D0)/SW)
     .     * 2.D0*VCHI/SW**2 / 2.D0/MW**2
     .     * LOOPSVV(MW**2,MH5**2,MH5**2,MW**2)

C Z W- W- loop:
      S = S + ALPHAEM**2/2.D0*(-1.D0)*(CW/SW)
     .     *(-DSQRT(2.D0)*VCHI/CW/SW**2)
     .     * LOOPVVV(MW**2,MH5**2,MZ**2,MW**2)

      RxH5PWGA = REAL(S)
      IxH5PWGA = AIMAG(S)
      HETLOOPH5WGA = MH5**3/32.D0/PI * (1.D0 - MW**2/MH5**2)**3
     .     * CDABS(S)**2

      ENDIF

      RETURN
      END


C=====================================================================
C H3+ --> W+ gamma
C=====================================================================
      DOUBLE PRECISION FUNCTION HETLOOPH3WGA()
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
C Functions to be called:
      DOUBLE COMPLEX LOOPSSS, LOOPVSS, LOOPSVV, LOOPVVV,
     .     LOOPFFF, TLOOPFFF
C Local variables:
      DOUBLE COMPLEX S, TS
      DOUBLE COMPLEX S3S, SVSS, SSVV, STB, TSTB
      DOUBLE PRECISION GH3NH3PH5M
      DOUBLE PRECISION GH3PH3MH5N, GHLH3PH3M, GHHH3PH3M, GH3PH3PH5MM
      DOUBLE PRECISION L1, L2, L3, L4, L5
      DOUBLE PRECISION SA, CA, SH, CH, SW, CW, MB, MT
      DOUBLE PRECISION PI
      PI = 4.D0*DATAN(1.D0)

      MT = MTPOLE
      MB = MBMB

      SA = DSIN(ALPHA)
      CA = DCOS(ALPHA)
      CH = VPHI/V
      SH = DSQRT(1.D0 - CH**2)
      CW = MW/MZ
      SW = DSQRT(1.D0 - CW**2)

      L1 = LAMBDA1
      L2 = LAMBDA2
      L3 = LAMBDA3
      L4 = LAMBDA4
      L5 = LAMBDA5

      IF (MH3.LE.MW) THEN
         HETLOOPH3WGA = 0.D0
         RxH3PWGA = 0.D0
         RxH3PWGATILDE = 0.D0
         IxH3PWGA = 0.D0
         IxH3PWGATILDE = 0.D0
      ELSE

C This first coupling has an i removed:
      GH3NH3PH5M = -DSQRT(2.D0)/V**2*(-8.D0*L5*VCHI**3 + 4.D0*M1*VCHI**2
     .     + (-4.D0*L5 + 2.D0*L3)*VPHI**2*VCHI + 3.D0*M2*VPHI**2)
      GH3PH3MH5N = DSQRT(2.D0/3.D0)/V**2 * (-8.D0*L5*VCHI**3
     .     + 4.D0*M1*VCHI**2 + (-4.D0*L5 + 2.D0*L3)*VPHI**2*VCHI
     .     + 3.D0*M2*VPHI**2)
      GHLH3PH3M = 64.D0*L1*CA*VCHI**2*VPHI/V**2 
     .     - 8.D0/DSQRT(3.D0)*VPHI**2*VCHI/V**2*SA*(L3+3.D0*L4)
     .     - 4.D0/DSQRT(3.D0)*VCHI*M1/V**2
     .     *(SA*VCHI - DSQRT(3.D0)*CA*VPHI)
     .     - 16.D0/DSQRT(3.D0)*VCHI**3/V**2*SA*(6.D0*L2+L5)
     .     - CA*VPHI**3/V**2*(L5-4.D0*L2)
     .     + 2.D0*DSQRT(3.D0)*M2*VPHI**2/V**2*SA
     .     - 8.D0/DSQRT(3.D0)*L5*VCHI*VPHI/V**2
     .     *(SA*VPHI-DSQRT(3.D0)*CA*VCHI)
      GHHH3PH3M = 64.D0*L1*SA*VCHI**2*VPHI/V**2
     .     + 8.D0/DSQRT(3.D0)*VPHI**2*VCHI/V**2*CA*(L3+3.D0*L4)
     .     + 4.D0/DSQRT(3.D0)*VCHI*M1/V**2
     .     *(CA*VCHI + DSQRT(3.D0)*SA*VPHI)
     .     + 16.D0/DSQRT(3.D0)*VCHI**3/V**2*CA*(6.D0*L2+L5)
     .     + SA*VPHI**3/V**2*(4.D0*L2 - L5)
     .     - 2.D0*DSQRT(3.D0)*M2*VPHI**2/V**2*CA
     .     + 8.D0/DSQRT(3.D0)*L5*VCHI*VPHI/V**2
     .     *(CA*VPHI + DSQRT(3.D0)*SA*VCHI)
      GH3PH3PH5MM = -2.D0/V**2 * (-8.D0*L5*VCHI**3 + 4.D0*M1*VCHI**2
     .     + (-4.D0*L5 + 2.D0*L3)*VPHI**2*VCHI + 3.D0*M2*VPHI**2)

      S = (0.D0,0.D0)

      S3S = (0.D0,0.D0)
      SVSS = (0.D0,0.D0)
      SSVV = (0.D0,0.D0)

C H50 H3- H3- loop:
      S3S = -ALPHAEM*(-1.D0)/PI*(-CH/2.D0/DSQRT(3.D0)/SW)
     .     * GH3PH3MH5N * LOOPSSS(MW**2,MH3**2,MH5**2,MH3**2)
C HL H3- H3- loop:
      S3S = S3S - ALPHAEM*(-1.D0)/PI
     .     *(DSQRT(2.D0/3.D0)*SA*CH + CA*SH/2.D0)/SW * GHLH3PH3M
     .     * LOOPSSS(MW**2,MH3**2,MHL**2,MH3**2)
C HH H3- H3- loop:
      S3S = S3S - ALPHAEM*(-1.D0)/PI
     .     *(-DSQRT(2.D0/3.D0)*CA*CH + SA*SH/2.D0)/SW * GHHH3PH3M
     .     * LOOPSSS(MW**2,MH3**2,MHH**2,MH3**2)
C H30 H5- H5- loop:
      S3S = S3S - ALPHAEM*(-1.D0)/PI * (-CH/2.D0/SW) * GH3NH3PH5M
     .     * LOOPSSS(MW**2,MH3**2,MH3**2,MH5**2)
C H3- H5-- H5-- loop:
      S3S = S3S - ALPHAEM*(-2.D0)/PI * (-CH/DSQRT(2.D0)/SW) 
     .     * GH3PH3PH5MM
     .     * LOOPSSS(MW**2,MH3**2,MH3**2,MH5**2)
C H5++ H3+ H3+ loop:
      S3S = S3S - ALPHAEM*(1.D0)/PI * (CH/DSQRT(2.D0)/SW) 
     .     * GH3PH3PH5MM
     .     * LOOPSSS(MW**2,MH3**2,MH5**2,MH3**2)

C Z H5- H5- loop:
      SVSS = -ALPHAEM**2*(-1.D0)*(-VPHI/2.D0/V/CW/SW) 
     .     * (-DSQRT(2.D0)*VCHI/CW/SW**2) / MZ**2
     .     * LOOPVSS(MW**2,MH3**2,MZ**2,MH5**2)
C W- H5-- H5-- loop:
      SVSS = SVSS - ALPHAEM**2*(-2.D0)*(VPHI/DSQRT(2.D0)/V/SW)
     .     * (2.D0*VCHI/SW**2) / MW**2
     .     * LOOPVSS(MW**2,MH3**2,MW**2,MH5**2)

C HL W- W- loop:
      SSVV = ALPHAEM**2*(-1.D0)
     .     *(-DSQRT(2.D0/3.D0)*(DSQRT(3.D0)*CA*VCHI/V + SA*VPHI/V)/SW)
     .     *(-8.D0*DSQRT(3.D0)*SA*VCHI + 3.D0*CA*VPHI)/6.D0/SW**2 
     .     / 2.D0/MW**2
     .     * LOOPSVV(MW**2,MH3**2,MHL**2,MW**2)
C HH W- W- loop:
      SSVV = SSVV + ALPHAEM**2*(-1.D0)
     .     *(-DSQRT(2.D0/3.D0)*(DSQRT(3.D0)*SA*VCHI/V - CA*VPHI/V)/SW)
     .     *(8.D0*DSQRT(3.D0)*CA*VCHI + 3.D0*SA*VPHI)/6.D0/SW**2
     .     / 2.D0/MW**2
     .     * LOOPSVV(MW**2,MH3**2,MHH**2,MW**2)
C H50 W- W- loop:
      SSVV = SSVV + ALPHAEM**2*(-1.D0)*(VPHI/2.D0/DSQRT(3.D0)/V/SW)
     .     *(DSQRT(2.D0/3.D0)*VCHI/SW**2) / 2.D0/MW**2
     .     * LOOPSVV(MW**2,MH3**2,MH5**2,MW**2)
C H5++ W+ W+ loop:
      SSVV = SSVV + ALPHAEM**2*(1.D0)*(VPHI/DSQRT(2.D0)/V/SW)
     .     *(2.D0*VCHI/SW**2) / 2.D0/MW**2
     .     * LOOPSVV(MW**2,MH3**2,MH5**2,MW**2)

C top-bottom loop: CP even part
      STB = ALPHAEM*3.D0*(SH/CH)/2.D0/PI/V/SW
     .     * LOOPFFF(MW**2,MH3**2,MT**2,MB**2)

C top-bottom loop: CP odd part
      TSTB = ALPHAEM*3.D0*(SH/CH)/2.D0/PI/V/SW
     .     * TLOOPFFF(MW**2,MH3**2,MT**2,MB**2)

      S = S3S + SVSS + SSVV + STB
      TS = TSTB

      RxH3PWGA = REAL(S)
      RxH3PWGATILDE = REAL(TS)
      IxH3PWGA = AIMAG(S)
      IxH3PWGATILDE = AIMAG(TS)
      HETLOOPH3WGA = MH3**3/32.D0/PI * (1.D0 - MW**2/MH3**2)**3
     .     * (CDABS(S)**2 + CDABS(TS)**2)

      ENDIF

      RETURN
      END


C=====================================================================
C H50 --> Z gamma
C=====================================================================
      DOUBLE PRECISION FUNCTION HETLOOPH5ZGA()
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
C Functions to be called:
C (The functions I1 and I2 are defined in gmdecays.f.)
      DOUBLE COMPLEX LOOPSSS, LOOPVSS, LOOPSVV, LOOPVVV,
     .     LOOPFFF, TLOOPFFF
      DOUBLE COMPLEX I1, I2
C Local variables:
      DOUBLE COMPLEX S
      DOUBLE PRECISION TAUW, LAMBDAW, G533P, G555P, G555PP, 
     .     CZ11, CZ22, KW
      DOUBLE COMPLEX AW, A3P, A5P, A5PP, AS, AMP
      DOUBLE PRECISION L1, L2, L3, L4, L5
      DOUBLE PRECISION CW, SW
      DOUBLE PRECISION PI
      PI = 4.D0*DATAN(1.D0)

      CW = MW/MZ
      SW = DSQRT(1.D0 - CW**2)

      L1 = LAMBDA1
      L2 = LAMBDA2
      L3 = LAMBDA3
      L4 = LAMBDA4
      L5 = LAMBDA5

      IF (MH5.LE.MZ) THEN
         HETLOOPH5ZGA = 0.D0
         RxH5NZGA = 0.D0
      ELSE

      TAUW = 4.D0*MW**2/MH5**2
      LAMBDAW = 4.D0*MW**2/MZ**2
      G533P = DSQRT(2.D0/3.D0)/V**2
     .     * (-8.D0*L5*VCHI**3 + 4.D0*M1*VCHI**2
     .     + (-4.D0*L5 + 2.D0*L3)*VPHI**2*VCHI
     .     + 3.D0*M2*VPHI**2)
      G555P = DSQRT(6.D0)*(2.D0*L3*VCHI - M2)
      G555PP = -2.D0*DSQRT(6.D0)*(2.D0*L3*VCHI - M2)
      KW = 2.D0*DSQRT(2.D0/3.D0)*VCHI/V

C SM amplitude for the W loop: 
      AW = -CW/SW * (4.D0*(3.D0-SW**2/CW**2) * I2(TAUW,LAMBDAW)
     .     + ((1.D0+2.D0/TAUW)*SW**2/CW**2 - (5.D0+2.D0/TAUW))
     .     * I1(TAUW,LAMBDAW))
C Amplitudes for the scalar loops:
      CZ11 = 1.D0/2.D0/SW/CW * (1.D0 - 2.D0*SW**2)
      CZ22 = 1.D0/SW/CW * (1.D0 - 2.D0*SW**2)
      A3P = 2.D0*G533P*CZ11/MH3**2
     .     * I1(4.D0*MH3**2/MH5**2,4.D0*MH3**2/MZ**2)
      A5P = 2.D0*G555P*CZ11/MH5**2
     .     * I1(4.D0*MH5**2/MH5**2,4.D0*MH5**2/MZ**2)
      A5PP = 2.D0*G555PP*CZ22*2.D0/MH5**2
     .     * I1(4.D0*MH5**2/MH5**2,4.D0*MH5**2/MZ**2)
      AS = A3P + A5P + A5PP
C Combine the amplitudes (all are complex in general): 
      AMP = KW*AW + V/2.D0*AS
      S = -ALPHAEM/2.D0/PI/V * AMP

C H5+ W+ W+ and H5- W- W- loops:
      S = S + 2.D0 * ALPHAEM**2*(1.D0)*(-DSQRT(3.D0)/2.D0/SW)
     .     * (-DSQRT(2.D0)*VCHI/CW/SW**2) / 2.D0/MW**2
     .     * LOOPSVV(MZ**2,MH5**2,MH5**2,MW**2)

C W+ H5+ H5+ and W- H5- H5- loops:
      S = S - 2.D0 * ALPHAEM**2*(1.D0)*(DSQRT(3.D0)/2.D0/SW)
     .     * (-DSQRT(2.D0)*VCHI/CW/SW**2) / MW**2
     .     * LOOPVSS(MZ**2,MH5**2,MW**2,MH5**2)

      RxH5NZGA = CDABS(S)
      HETLOOPH5ZGA = MH5**3 / 32.D0/PI * (1.D0 - MZ**2/MH5**2)**3
     .     * CDABS(S)**2

      ENDIF

      RETURN
      END


C############################################################################
C Subroutines for the new "heterogeneous" loop functions for H5+ -> W+ gamma, 
C H50 -> Z gamma, and H3+ -> W+ gamma, computed using LoopTools.
C K. Hartling and H.E. Logan, December 2015

C Three scalars in the loop:
      DOUBLE COMPLEX FUNCTION LOOPSSS(KSQ,MHSQ,M1SQ,M2SQ)
      IMPLICIT NONE
      DOUBLE PRECISION KSQ, MHSQ, M1SQ, M2SQ
      DOUBLE PRECISION QSQ
#include "looptools.h"
      QSQ = 0.D0
      LOOPSSS = C0I(cc12,KSQ,QSQ,MHSQ,M1SQ,M2SQ,M2SQ)
     .     + C0I(cc22,KSQ,QSQ,MHSQ,M1SQ,M2SQ,M2SQ)
     .     + C0I(cc2,KSQ,QSQ,MHSQ,M1SQ,M2SQ,M2SQ)
      END

C One vector and two scalars in the loop:
      DOUBLE COMPLEX FUNCTION LOOPVSS(KSQ,MHSQ,MVSQ,MSSQ)
      IMPLICIT NONE
      DOUBLE PRECISION KSQ, MHSQ, MSSQ, MVSQ
      DOUBLE PRECISION QSQ
#include "looptools.h"
      QSQ = 0.D0
      LOOPVSS = MVSQ * ((-2.D0) *( 
     .     C0I(cc12,KSQ,QSQ,MHSQ,MVSQ,MSSQ,MSSQ)
     .     + C0I(cc22,KSQ,QSQ,MHSQ,MVSQ,MSSQ,MSSQ)
     .     + 2.D0*C0I(cc1,KSQ,QSQ,MHSQ,MVSQ,MSSQ,MSSQ)
     .     + 3.D0*C0I(cc2,KSQ,QSQ,MHSQ,MVSQ,MSSQ,MSSQ)
     .     + 2.D0*C0I(cc0,KSQ,QSQ,MHSQ,MVSQ,MSSQ,MSSQ))
     .     - 2.D0 * (MHSQ - MSSQ)/MVSQ *(
     .     C0I(cc12,KSQ,QSQ,MHSQ,MVSQ,MSSQ,MSSQ)
     .     + C0I(cc22,KSQ,QSQ,MHSQ,MVSQ,MSSQ,MSSQ)
     .     + C0I(cc2,KSQ,QSQ,MHSQ,MVSQ,MSSQ,MSSQ) ))
      END

C One scalar and two vectors in the loop:
      DOUBLE COMPLEX FUNCTION LOOPSVV(KSQ,MHSQ,MSSQ,MVSQ)
      IMPLICIT NONE
      DOUBLE PRECISION KSQ, MHSQ, MSSQ, MVSQ
      DOUBLE PRECISION QSQ
#include "looptools.h"
      QSQ = 0.D0
      LOOPSVV = 2.D0*MVSQ
     .     * ( -2.D0*C0I(cc12,KSQ,QSQ,MHSQ,MSSQ,MVSQ,MVSQ)
     .     - 2.D0*C0I(cc22,KSQ,QSQ,MHSQ,MSSQ,MVSQ,MVSQ)
     .     + 4.D0*C0I(cc1,KSQ,QSQ,MHSQ,MSSQ,MVSQ,MVSQ)
     .     + 2.D0*C0I(cc2,KSQ,QSQ,MHSQ,MSSQ,MVSQ,MVSQ)
     .     - (MHSQ - MSSQ)/MVSQ * 2.D0
     .     * (C0I(cc12,KSQ,QSQ,MHSQ,MSSQ,MVSQ,MVSQ)
     .     + C0I(cc22,KSQ,QSQ,MHSQ,MSSQ,MVSQ,MVSQ)
     .     + C0I(cc2,KSQ,QSQ,MHSQ,MSSQ,MVSQ,MVSQ)))
      END

C Three vectors in the loop:
      DOUBLE COMPLEX FUNCTION LOOPVVV(KSQ,MHSQ,MZSQ,MWSQ)
      IMPLICIT NONE
      DOUBLE PRECISION KSQ, MHSQ, MZSQ, MWSQ
      DOUBLE PRECISION QSQ
      DOUBLE PRECISION MZ, MW
#include "looptools.h"
      QSQ = 0.D0
      MZ = DSQRT(MZSQ)
      MW = DSQRT(MWSQ)
      LOOPVVV = -2.D0 * (
     .     (10.D0*C0I(cc12,KSQ,QSQ,MHSQ,MZSQ,MWSQ,MWSQ)
     .     + 10.D0*C0I(cc22,KSQ,QSQ,MHSQ,MZSQ,MWSQ,MWSQ)
     .     + C0I(cc1,KSQ,QSQ,MHSQ,MZSQ,MWSQ,MWSQ)
     .     + 10.D0*C0I(cc2,KSQ,QSQ,MHSQ,MZSQ,MWSQ,MWSQ)
     .     + 5.D0*C0I(cc0,KSQ,QSQ,MHSQ,MZSQ,MWSQ,MWSQ))
     .     + 0.5D0*(C0I(cc12,KSQ,QSQ,MHSQ,MZSQ,MWSQ,MWSQ)
     .     + C0I(cc22,KSQ,QSQ,MHSQ,MZSQ,MWSQ,MWSQ)
     .     - 2.D0*C0I(cc1,KSQ,QSQ,MHSQ,MZSQ,MWSQ,MWSQ)
     .     + 3.D0*C0I(cc2,KSQ,QSQ,MHSQ,MZSQ,MWSQ,MWSQ)
     .     + 2.D0*C0I(cc0,KSQ,QSQ,MHSQ,MZSQ,MWSQ,MWSQ))
     .     + 0.5D0*(C0I(cc12,KSQ,QSQ,MHSQ,MZSQ,MWSQ,MWSQ)
     .     + C0I(cc22,KSQ,QSQ,MHSQ,MZSQ,MWSQ,MWSQ)
     .     - C0I(cc2,KSQ,QSQ,MHSQ,MZSQ,MWSQ,MWSQ))
     .     + MZ**2/MW**2*(1.D0-MWSQ/MZSQ) 
     .     * (C0I(cc12,KSQ,QSQ,MHSQ,MZSQ,MWSQ,MWSQ)
     .     + C0I(cc22,KSQ,QSQ,MHSQ,MZSQ,MWSQ,MWSQ)
     .     + 2.D0*C0I(cc1,KSQ,QSQ,MHSQ,MZSQ,MWSQ,MWSQ)
     .     + 3.D0*C0I(cc2,KSQ,QSQ,MHSQ,MZSQ,MWSQ,MWSQ)
     .     + 2.D0*C0I(cc0,KSQ,QSQ,MHSQ,MZSQ,MWSQ,MWSQ))
     .     + MHSQ/MW**2 
     .     * (C0I(cc12,KSQ,QSQ,MHSQ,MZSQ,MWSQ,MWSQ)
     .     + C0I(cc22,KSQ,QSQ,MHSQ,MZSQ,MWSQ,MWSQ)
     .     + C0I(cc2,KSQ,QSQ,MHSQ,MZSQ,MWSQ,MWSQ))
     .     + (C0I(cc12,KSQ,QSQ,MHSQ,MZSQ,MWSQ,MWSQ)
     .     + C0I(cc22,KSQ,QSQ,MHSQ,MZSQ,MWSQ,MWSQ)
     .     + C0I(cc2,KSQ,QSQ,MHSQ,MZSQ,MWSQ,MWSQ)))
      END

C Three fermions in the loop: CP-even part
      DOUBLE COMPLEX FUNCTION LOOPFFF(KSQ,MHSQ,MTSQ,MBSQ)
      IMPLICIT NONE
      DOUBLE PRECISION KSQ, MHSQ, MTSQ, MBSQ
      DOUBLE PRECISION QSQ, QT, QB
      DOUBLE COMPLEX LOOPTBB, LOOPBTT
#include "looptools.h"
      QSQ = 0.D0
      QT = 2.D0/3.D0
      QB = -1.D0/3.D0
      LOOPTBB = MTSQ * (2.D0*C0I(cc12,KSQ,QSQ,MHSQ,MTSQ,MBSQ,MBSQ)
     .     + 2.D0*C0I(cc22,KSQ,QSQ,MHSQ,MTSQ,MBSQ,MBSQ)
     .     + 3.D0*C0I(cc2,KSQ,QSQ,MHSQ,MTSQ,MBSQ,MBSQ)
     .     + C0I(cc1,KSQ,QSQ,MHSQ,MTSQ,MBSQ,MBSQ)
     .     + C0I(cc0,KSQ,QSQ,MHSQ,MTSQ,MBSQ,MBSQ))
     .     - MBSQ * (2.D0*C0I(cc12,KSQ,QSQ,MHSQ,MTSQ,MBSQ,MBSQ)
     .     + 2.D0*C0I(cc22,KSQ,QSQ,MHSQ,MTSQ,MBSQ,MBSQ)
     .     + C0I(cc2,KSQ,QSQ,MHSQ,MTSQ,MBSQ,MBSQ)
     .     - C0I(cc1,KSQ,QSQ,MHSQ,MTSQ,MBSQ,MBSQ))
      LOOPBTT = -MBSQ * (2.D0*C0I(cc12,KSQ,QSQ,MHSQ,MBSQ,MTSQ,MTSQ)
     .     + 2.D0*C0I(cc22,KSQ,QSQ,MHSQ,MBSQ,MTSQ,MTSQ)
     .     + 3.D0*C0I(cc2,KSQ,QSQ,MHSQ,MBSQ,MTSQ,MTSQ)
     .     + C0I(cc1,KSQ,QSQ,MHSQ,MBSQ,MTSQ,MTSQ)
     .     + C0I(cc0,KSQ,QSQ,MHSQ,MBSQ,MTSQ,MTSQ))
     .     + MTSQ * (2.D0*C0I(cc12,KSQ,QSQ,MHSQ,MBSQ,MTSQ,MTSQ)
     .     + 2.D0*C0I(cc22,KSQ,QSQ,MHSQ,MBSQ,MTSQ,MTSQ)
     .     + C0I(cc2,KSQ,QSQ,MHSQ,MBSQ,MTSQ,MTSQ)
     .     - C0I(cc1,KSQ,QSQ,MHSQ,MBSQ,MTSQ,MTSQ))
      LOOPFFF = QB*LOOPTBB + QT*LOOPBTT
      END

C Three fermions in the loop: CP-odd part (tilde'd amplitude)
      DOUBLE COMPLEX FUNCTION TLOOPFFF(KSQ,MHSQ,MTSQ,MBSQ)
      IMPLICIT NONE
      DOUBLE PRECISION KSQ, MHSQ, MTSQ, MBSQ
      DOUBLE PRECISION QSQ, QT, QB
      DOUBLE COMPLEX TLOOPTBB, TLOOPBTT
#include "looptools.h"
      QSQ = 0.D0
      QT = 2.D0/3.D0
      QB = -1.D0/3.D0
      TLOOPTBB = -MTSQ * (C0I(cc1,KSQ,QSQ,MHSQ,MTSQ,MBSQ,MBSQ)
     .     + C0I(cc2,KSQ,QSQ,MHSQ,MTSQ,MBSQ,MBSQ)
     .     + C0I(cc0,KSQ,QSQ,MHSQ,MTSQ,MBSQ,MBSQ))
     .     + MBSQ * (C0I(cc1,KSQ,QSQ,MHSQ,MTSQ,MBSQ,MBSQ)
     .     + C0I(cc2,KSQ,QSQ,MHSQ,MTSQ,MBSQ,MBSQ))
      TLOOPBTT = -MBSQ * (C0I(cc1,KSQ,QSQ,MHSQ,MBSQ,MTSQ,MTSQ)
     .     + C0I(cc2,KSQ,QSQ,MHSQ,MBSQ,MTSQ,MTSQ)
     .     + C0I(cc0,KSQ,QSQ,MHSQ,MBSQ,MTSQ,MTSQ))
     .     + MTSQ * (C0I(cc1,KSQ,QSQ,MHSQ,MBSQ,MTSQ,MTSQ)
     .     + C0I(cc2,KSQ,QSQ,MHSQ,MBSQ,MTSQ,MTSQ))
      TLOOPFFF = QB*TLOOPTBB + QT*TLOOPBTT
      END
