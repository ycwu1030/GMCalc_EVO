      SUBROUTINE INIT_QCD
C Initializes the QCD strong-coupling scales LAMBDA4 and LAMBDA5.
C Algorithm idea based on 2HDMC, Eriksson, Rathsman & Stal, 0902.0851      
      IMPLICIT NONE
C Common blocks:
      DOUBLE PRECISION LAMBDA4, LAMBDA5
      COMMON/LAMQCD/LAMBDA4,LAMBDA5
      DOUBLE PRECISION V, MZ, MW, MTPOLE, MBMB, MCMC, MS, MTAU, MMU, 
     .     ALPHAEM,ALSMZ, VCB, VUS, VUB, GAMZ, GAMW
      COMMON/SM/V,MZ,MW,MTPOLE,MBMB,MCMC,MS,MTAU,MMU,
     .     ALPHAEM,ALSMZ,VCB,VUS,VUB,GAMZ,GAMW
C External functions to be called:
      DOUBLE PRECISION RUNALS
C Local variables:
      DOUBLE PRECISION TOLERANCE
      DOUBLE PRECISION LAMBDAMIN, LAMBDAMAX, LAMBDAMID
      DOUBLE PRECISION ALSMIN, ALSMID
      DOUBLE PRECISION ALSMATCH

      TOLERANCE = 1.D-6

C First initialize LAMBDA5.  
      LAMBDAMIN = 1.D-3
      LAMBDAMAX = 1.D0

 60   LAMBDA5 = LAMBDAMIN
      ALSMIN = RUNALS(MZ,5)
      LAMBDAMID = (LAMBDAMIN+LAMBDAMAX)/2.D0
      LAMBDA5 = LAMBDAMID
      ALSMID = RUNALS(MZ,5)
      IF ((ALSMZ.GE.ALSMIN).AND.(ALSMZ.LE.ALSMID)) THEN
         LAMBDAMAX = LAMBDAMID
      ELSE
         LAMBDAMIN = LAMBDAMID
      ENDIF
      IF (DABS(1.D0-ALSMIN/ALSMZ).GE.TOLERANCE) THEN
         GOTO 60
      ENDIF
      LAMBDA5 = LAMBDAMIN

C Now initialize LAMBDA4.  Match by requiring continuity of alpha_s at
C the b quark mass threshold, as done in SuperIso v3.3 manual appendix A.
      LAMBDAMIN = 1.D-3
      LAMBDAMAX = 1.D0
      ALSMATCH = RUNALS(MBMB,5)

 61   LAMBDA4 = LAMBDAMIN
      ALSMIN = RUNALS(MBMB,4)
      LAMBDAMID = (LAMBDAMIN+LAMBDAMAX)/2.D0
      LAMBDA4 = LAMBDAMID
      ALSMID = RUNALS(MBMB,4)
      IF ((ALSMATCH.GE.ALSMIN).AND.(ALSMATCH.LE.ALSMID)) THEN
         LAMBDAMAX = LAMBDAMID
      ELSE
         LAMBDAMIN = LAMBDAMID
      ENDIF
      IF (DABS(1.D0-ALSMIN/ALSMATCH).GE.TOLERANCE) THEN
         GOTO 61
      ENDIF
      LAMBDA4 = LAMBDAMIN

      RETURN
      END


      DOUBLE PRECISION FUNCTION RUNALS(MU,NF)
C Running strong coupling constant alpha_s at scale MU, evaluated in the 
C NF flavor scheme.  NF = 4 and 5 are implemented.  Uses LAMBDA4 and 
C LAMBDA5, which are the (4fs and 5fs) scales at which alpha_s diverges;
C these parameters are set by subroutine INIT_QCD.
C From Eq.(10) of Djouadi, Spira & Zerwas, hep-ph/9511344
      IMPLICIT NONE
      DOUBLE PRECISION MU
      INTEGER NF
C Common blocks:
      DOUBLE PRECISION LAMBDA4, LAMBDA5
      COMMON/LAMQCD/LAMBDA4,LAMBDA5
C Local variables:
      DOUBLE PRECISION PI, LAMBDA
      PI = 4.D0*DATAN(1.D0)

      IF (NF.EQ.4) THEN 
         LAMBDA = LAMBDA4
      ELSE IF (NF.EQ.5) THEN
         LAMBDA = LAMBDA5
      ENDIF

      RUNALS = 12.D0*PI/(33.D0-2.D0*NF)/DLOG(MU**2/LAMBDA**2)
     .     * (1.D0 - 6.D0*(153.D0-19.D0*NF)/(33.D0-2.D0*NF)**2
     .       * DLOG(DLOG(MU**2/LAMBDA**2))/DLOG(MU**2/LAMBDA**2))
      
      RETURN
      END

      DOUBLE PRECISION FUNCTION RUNMB(MU)
C Running bottom quark mass evaluated at scale MU.
C Only scales MU at or above the charm mass mc are implemented. 
C From Eq.(4) of Djouadi, Spira & Zerwas, hep-ph/9511344
      IMPLICIT NONE
      DOUBLE PRECISION MU
C Common blocks:
      DOUBLE PRECISION V, MZ, MW, MTPOLE, MBMB, MCMC, MS, MTAU, MMU, 
     .     ALPHAEM, ALSMZ, VCB, VUS, VUB, GAMZ, GAMW
      COMMON/SM/V,MZ,MW,MTPOLE,MBMB,MCMC,MS,MTAU,MMU,
     .     ALPHAEM,ALSMZ,VCB,VUS,VUB,GAMZ,GAMW
C External functions to be called:
      DOUBLE PRECISION C4QCD, C5QCD
      DOUBLE PRECISION RUNALS
C Local variables:
      DOUBLE PRECISION XMU, XMB
      DOUBLE PRECISION PI
      PI = 4.D0*DATAN(1.D0)

      IF (MU.GE.MBMB) THEN
         XMU = RUNALS(MU,5)/PI
         XMB = RUNALS(MBMB,5)/PI
         RUNMB = MBMB * C5QCD(XMU)/C5QCD(XMB)
      ELSE 
         XMU = RUNALS(MU,4)/PI
         XMB = RUNALS(MBMB,4)/PI
         RUNMB = MBMB * C4QCD(XMU)/C4QCD(XMB)
      ENDIF

      RETURN
      END

      DOUBLE PRECISION FUNCTION RUNMC(MU)
C Running charm quark mass evaluated at scale MU.
C Only scales MU at or above the charm mass mc are implemented.
C From Eq.(4) of Djouadi, Spira & Zerwas, hep-ph/9511344
      IMPLICIT NONE
      DOUBLE PRECISION MU
C Common blocks:
      DOUBLE PRECISION V, MZ, MW, MTPOLE, MBMB, MCMC, MS, MTAU, MMU, 
     .     ALPHAEM, ALSMZ, VCB, VUS, VUB, GAMZ, GAMW
      COMMON/SM/V,MZ,MW,MTPOLE,MBMB,MCMC,MS,MTAU,MMU,
     .     ALPHAEM,ALSMZ,VCB,VUS,VUB,GAMZ,GAMW
C External functions to be called:
      DOUBLE PRECISION C4QCD, C5QCD
      DOUBLE PRECISION RUNALS
C Local variables:
      DOUBLE PRECISION XMU, XMC
      DOUBLE PRECISION MCMB, XMB
      DOUBLE PRECISION PI
      PI = 4.D0*DATAN(1.D0)

      IF (MU.LE.MBMB) THEN
         XMU = RUNALS(MU,4)/PI
         XMC = RUNALS(MCMC,4)/PI
         RUNMC = MCMC * C4QCD(XMU)/C4QCD(XMC)
      ELSE
         XMB = RUNALS(MBMB,4)/PI
         XMC = RUNALS(MCMC,4)/PI
         MCMB = MCMC * C4QCD(XMB)/C4QCD(XMC)
         XMU = RUNALS(MU,5)/PI
         XMB = RUNALS(MBMB,5)/PI
         RUNMC = MCMB * C5QCD(XMU)/C5QCD(XMB)
      ENDIF

      RETURN
      END

      DOUBLE PRECISION FUNCTION RUNMT(MU)
C Running top quark mass evaluated at scale MU.
C Only scales MU at or above the charm mass mc are implemented. 
C From Eq.(4) of Djouadi, Spira & Zerwas, hep-ph/9511344
      IMPLICIT NONE
      DOUBLE PRECISION MU
C Common blocks:
      DOUBLE PRECISION V, MZ, MW, MTPOLE, MBMB, MCMC, MS, MTAU, MMU, 
     .     ALPHAEM, ALSMZ, VCB, VUS, VUB, GAMZ, GAMW
      COMMON/SM/V,MZ,MW,MTPOLE,MBMB,MCMC,MS,MTAU,MMU,
     .     ALPHAEM,ALSMZ,VCB,VUS,VUB,GAMZ,GAMW
C External functions to be called:
      DOUBLE PRECISION C4QCD, C5QCD
      DOUBLE PRECISION RUNALS, MSBARMQ
C Local variables:
      DOUBLE PRECISION XMU, XMT, MTMT
      DOUBLE PRECISION XMB, MTMB
      DOUBLE PRECISION PI
      PI = 4.D0*DATAN(1.D0)

      MTMT = MSBARMQ(MTPOLE,5)

      IF (MU.GE.MBMB) THEN
         XMU = RUNALS(MU,5)/PI
         XMT = RUNALS(MTMT,5)/PI
         RUNMT = MTMT * C5QCD(XMU)/C5QCD(XMT)
      ELSE 
         XMT = RUNALS(MTMT,5)/PI
         XMB = RUNALS(MBMB,5)/PI
         MTMB = MTMT * C5QCD(XMB)/C5QCD(XMT)
         XMB = RUNALS(MBMB,4)/PI
         XMU = RUNALS(MU,4)/PI
         RUNMT = MTMB * C4QCD(XMU)/C4QCD(XMB)
      ENDIF

      RETURN
      END

      
      DOUBLE PRECISION FUNCTION C4QCD(X)
      IMPLICIT NONE
      DOUBLE PRECISION X
      C4QCD = (25.D0/6.D0*X)**(12.D0/25.D0)
     .     * (1.D0 + 1.014D0*X + 1.389D0*X**2)
      RETURN
      END

      DOUBLE PRECISION FUNCTION C5QCD(X)
      IMPLICIT NONE
      DOUBLE PRECISION X
      C5QCD = (23.D0/6.D0*X)**(12.D0/23.D0)
     .     * (1.D0 + 1.175D0*X + 1.501D0*X**2)
      RETURN
      END


      DOUBLE PRECISION FUNCTION DELTAQCD(MH)
C QCD correction to h -> qqbar decays
C From Eq.(3) of Djouadi, Spira & Zerwas, hep-ph/9511344
      IMPLICIT NONE
      DOUBLE PRECISION MH
C Common blocks:
      DOUBLE PRECISION V, MZ, MW, MTPOLE, MBMB, MCMC, MS, MTAU, MMU, 
     .     ALPHAEM, ALSMZ, VCB, VUS, VUB, GAMZ, GAMW
      COMMON/SM/V,MZ,MW,MTPOLE,MBMB,MCMC,MS,MTAU,MMU,
     .     ALPHAEM,ALSMZ,VCB,VUS,VUB,GAMZ,GAMW
C External functions to be called:
      DOUBLE PRECISION RUNALS
C Local variables:
      DOUBLE PRECISION PI
      INTEGER NF

      PI = 4.D0*DATAN(1.D0)

      IF (MH.GE.MBMB) THEN 
         NF = 5
      ELSE
         NF = 4
      ENDIF

      DELTAQCD = 1.D0 + 5.67D0*RUNALS(MH,NF)/PI
     .     + (35.94 - 1.36*NF)*(RUNALS(MH,NF)/PI)**2

      RETURN
      END

      DOUBLE PRECISION FUNCTION DELTATB(MH)
C Top-mass-dependent two-loop QCD correction to h -> bbbar
C From Eq.(3) of Djouadi, Spira & Zerwas, hep-ph/9511344
      IMPLICIT NONE
      DOUBLE PRECISION MH
C Common blocks:
      DOUBLE PRECISION V, MZ, MW, MTPOLE, MBMB, MCMC, MS, MTAU, MMU, 
     .     ALPHAEM, ALSMZ, VCB, VUS, VUB, GAMZ, GAMW
      COMMON/SM/V,MZ,MW,MTPOLE,MBMB,MCMC,MS,MTAU,MMU,
     .     ALPHAEM,ALSMZ,VCB,VUS,VUB,GAMZ,GAMW
C External functions to be called:
      DOUBLE PRECISION RUNALS, RUNMB
C Local variables:
      DOUBLE PRECISION PI
      INTEGER NF

      PI = 4.D0*DATAN(1.D0)

      IF (MH.GE.MBMB) THEN 
         NF = 5
      ELSE
         NF = 4
      ENDIF

      DELTATB = (RUNALS(MH,NF)/PI)**2 
     .     * (1.57D0 - 2.D0/3.D0*DLOG(MH**2/MTPOLE**2)
     .       + 1.D0/9.D0 * (DLOG(RUNMB(MH)**2/MH**2))**2)

      RETURN
      END

      DOUBLE PRECISION FUNCTION DELTATC(MH)
C Top-mass-dependent two-loop QCD correction to h -> ccbar
C From Eq.(3) of Djouadi, Spira & Zerwas, hep-ph/9511344
      IMPLICIT NONE
      DOUBLE PRECISION MH
C Common blocks:
      DOUBLE PRECISION V, MZ, MW, MTPOLE, MBMB, MCMC, MS, MTAU, MMU, 
     .     ALPHAEM, ALSMZ, VCB, VUS, VUB, GAMZ, GAMW
      COMMON/SM/V,MZ,MW,MTPOLE,MBMB,MCMC,MS,MTAU,MMU,
     .     ALPHAEM,ALSMZ,VCB,VUS,VUB,GAMZ,GAMW
C External functions to be called:
      DOUBLE PRECISION RUNALS, RUNMC
C Local variables:
      DOUBLE PRECISION PI
      INTEGER NF

      PI = 4.D0*DATAN(1.D0)

      IF (MH.GE.MBMB) THEN 
         NF = 5
      ELSE
         NF = 4
      ENDIF

      DELTATC = (RUNALS(MH,NF)/PI)**2 
     .     * (1.57D0 - 2.D0/3.D0*DLOG(MH**2/MTPOLE**2)
     .       + 1.D0/9.D0 * (DLOG(RUNMC(MH)**2/MH**2))**2)

      RETURN
      END

      DOUBLE PRECISION FUNCTION POLEMQ(MQMQ,NF)
C Computes the pole mass of quark Q, given the running mass at scale MQ.
C Based on the caption of Table 1 in Djouadi, Spira & Zerwas, hep-ph/9511344
      IMPLICIT NONE
      DOUBLE PRECISION MQMQ
      INTEGER NF
      DOUBLE PRECISION RUNALS
      DOUBLE PRECISION PI
      PI = 4.D0*DATAN(1.D0)

      POLEMQ = MQMQ * (1.D0 + 4.D0*RUNALS(MQMQ,NF)/3.D0/PI)

      RETURN
      END

      DOUBLE PRECISION FUNCTION MSBARMQ(MQPOLE,NF)
C Reverses the function of POLEMQ by calculating the MSbar running
C mass at scale MQ from the pole mass.
      IMPLICIT NONE
      DOUBLE PRECISION MQPOLE
      INTEGER NF
      DOUBLE PRECISION RUNALS
      DOUBLE PRECISION PI
      PI = 4.D0*DATAN(1.D0)
     
      MSBARMQ = MQPOLE / (1.D0 + 4.D0*RUNALS(MQPOLE,NF)/3.D0/PI)

      RETURN
      END

