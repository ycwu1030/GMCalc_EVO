C Code in this file written by Mark B. Reimer, winter 2018.
C Note: In order for this new GAMVVOF function to work, it uses the
C following functions from the original code:
C DOUBLE PRECISION FUNCTION LAMBDA(X,Y)
C DOUBLE PRECISION FUNCTION GAMVV(SV,MH,MV1,MV2,C)
C DOUBLE PRECISION FUNCTION RHO(Q,M,GA)
C DOUBLE PRECISION FUNCTION QFROMRHO(RHO,M,GA)

      DOUBLE PRECISION FUNCTION GAMVVOF(SV,MH,MV1,MV2,GA1,GA2,C)
C Doubly-offshell Gamma(H -> V1* V2*) numerical integration routine.
C Symmetry factor: SV = 1 for V1 != V2, SV = 1/2 for V1 = V2.
C C is the HVV coupling, Feynman rule = i C g_{munu}.
C This function uses the VEGAS algorithm or 16-point
C Gauss-Legendre quadrature, depending on the value of MH
C relative to MV1+MV2.
      IMPLICIT NONE
      DOUBLE PRECISION SV, MH, MV1, MV2, GA1, GA2, C
C Functions to be called
      DOUBLE PRECISION GAMVVOFVEGAS
      DOUBLE PRECISION GAMVVOFGL16
      
      IF (MH.GT.(MV1+MV2)) THEN
C We are above threshold.
      GAMVVOF = GAMVVOFGL16(SV,MH,MV1,MV2,GA1,GA2,C)
      ELSE
C We are below threshold; need more integration points.
      GAMVVOF = GAMVVOFVEGAS(SV,MH,MV1,MV2,GA1,GA2,C)
      ENDIF
      END FUNCTION



      DOUBLE PRECISION FUNCTION GAMVVOFVEGAS(SV,MH,MV1,MV2,GA1,GA2,C)
C Doubly-offshell Gamma(H -> V1* V2*) numerical integration routine.
C Symmetry factor: SV = 1 for V1 != V2, SV = 1/2 for V1 = V2.
C C is the HVV coupling, Feynman rule = i C g_{munu}.
C This function uses the VEGAS algorithm for integration.
      IMPLICIT NONE
C variables to be passed to gvvofintfnq
C Common block is initialized here before one calls gvvofintfnq
      DOUBLE PRECISION SV,MH,MV1,MV2,GA1,GA2,C
      DOUBLE PRECISION SVP,MHP,MV1P,MV2P,GA1P,GA2P,CP
      COMMON/GVVOFPARAMS/SVP,MHP,MV1P,MV2P,GA1P,GA2P,CP
      
C variables used by the INTEG subroutine from vegas
      INTEGER IDIM
      INTEGER IPOINT, ITER, IPOINT1, ITER1
      DOUBLE PRECISION ACC, RES
      DOUBLE PRECISION VAR(2)
      
      DOUBLE PRECISION DUMMYVAR

C functions called
      DOUBLE PRECISION GVVOFINTFNQ

      SVP = SV
      MHP = MH
      MV1P = MV1
      MV2P = MV2
      GA1P = GA1
      GA2P = GA2
      CP = C
      
      IDIM = 2
      IPOINT = 10000
      ITER = 3
      IPOINT1 = 10000
      ITER1 = 3
      ACC = 1.D-03
      
      DATA VAR /0.4,0.4/
      
C Make a call to the function to be integrated, so that INTEG
C recognizes it as a function.
      DUMMYVAR = GVVOFINTFNQ(VAR)
      
C Initialize vegas
      CALL RSTART(12,34,56,78)
      
      CALL INTEG(GVVOFINTFNQ,IDIM,IPOINT,ITER,IPOINT1,ITER1,ACC,RES)
      
      GAMVVOFVEGAS = RES
      END FUNCTION


      DOUBLE PRECISION FUNCTION GVVOFINTFNQ(VAR)
C VAR(1) and VAR(2) are between 0 and 1.
      IMPLICIT NONE
      DOUBLE PRECISION VAR(2)
      DOUBLE PRECISION PI
C Denormalized variables
      DOUBLE PRECISION Q1,Q2,Q1MIN,Q2MIN,Q1MAX,Q2MAX
C Parts of overall function to be integrated
      DOUBLE PRECISION DR1DQ1, DR2DQ2
C Ratios between VARs and Qs
      DOUBLE PRECISION DQ1DVAR1, DQ2DVAR2
C Denormalized function value
      DOUBLE PRECISION GAMVVOF_FN
C functions to be called
      DOUBLE PRECISION RHO,QFROMRHO,GAMVV
C Common block must be initialized by calling program
      DOUBLE PRECISION SV,MH,MV1,MV2,GA1,GA2,C
      COMMON/GVVOFPARAMS/SV,MH,MV1,MV2,GA1,GA2,C
      
      PI = 4.D0*DATAN(1.D0)
      
      Q1MIN = 0.D0
      Q2MIN = 0.D0
      Q1MAX = MH
      Q2MAX = MH
      
      DQ1DVAR1 = Q1MAX - Q1MIN
      DQ2DVAR2 = Q2MAX - Q2MIN
      
      Q1 = Q1MIN + DQ1DVAR1 * VAR(1)
      Q2 = Q2MIN + DQ2DVAR2 * VAR(2)
      
      DR1DQ1 = (2/PI)*(MV1*GA1*Q1) /
     .            (MV1**2 * GA1**2 + (Q1**2-MV1**2)**2)
      DR2DQ2 = (2/PI)*(MV2*GA2*Q2) /
     .            (MV2**2 * GA2**2 + (Q2**2-MV2**2)**2)
      
      IF ((Q1 .GT. 0.D0) .AND.
     .      (Q2 .GT. 0.D0) .AND.
     .         (MH - Q1 - Q2 .GT. 0.D0)) THEN
         GAMVVOF_FN = DR1DQ1 * DR2DQ2
     .                 * Q1**2/MV1**2 * Q2**2/MV2**2
     .                  * GAMVV(SV,MH,Q1,Q2,C)
         GVVOFINTFNQ = DQ1DVAR1 * DQ2DVAR2 * GAMVVOF_FN
      ELSE
         GVVOFINTFNQ = 0.D0
      ENDIF 
      END FUNCTION




      DOUBLE PRECISION FUNCTION GAMVVOFGL16(SV,MH,MV1,MV2,GA1,GA2,C)
C Doubly-offshell Gamma(H -> V1* V2*) numerical integration routine.
C Symmetry factor: SV = 1 for V1 != V2, SV = 1/2 for V1 = V2.
C C is the HVV coupling, Feynman rule = i C g_{munu}.
C This function uses 16-point Gauss-Legendre quadrature.
      IMPLICIT NONE
      DOUBLE PRECISION SV, MH, MV1, MV2, GA1, GA2, C
      DOUBLE PRECISION Q1, Q2
      DOUBLE PRECISION R1, R2
      DOUBLE PRECISION R1I, R1F, R2I, R2F
      DOUBLE PRECISION SLICE
      INTEGER I, J, NSTEP
      DOUBLE PRECISION X1M, X1R, DX1, X2M, X2R, DX2
C Half the number of points used in the quadrature
      INTEGER GLPTSDIV2
      PARAMETER (GLPTSDIV2=8)
      DOUBLE PRECISION X(GLPTSDIV2),W(GLPTSDIV2)
C Sample points chosen
      DATA X/.0950125098376374D0, .2816035507792589D0,
     .   .4580167776572274D0, .6178762444026438D0,
     .   .7554044083550030D0, .8656312023878318D0,
     .   .9445750230732326D0, .9894009349916499D0/
C Weights of sample points
      DATA W/.1894506104550685D0, .1826034150449236D0,
     .   .1691565193950025D0, .1495959888165767D0,
     .   .1246289712555339D0, .0951585116824928D0,
     .   .0622535239386479D0, .0271524594117541D0/
C Functions to be called
      DOUBLE PRECISION GAMVV
      DOUBLE PRECISION RHO, QFROMRHO
      
C Set the limits of integration:
         R1I = RHO(0.D0,MV1,GA1)
         R1F = RHO(MH,MV1,GA1)
         
C The integral ranges from X1M-X1R to X1M+X1R
         X1M=0.5D0*(R1F+R1I)
         X1R=0.5D0*(R1F-R1I)
         
         R2I = RHO(0.D0,MV2,GA2)
         GAMVVOFGL16 = 0.D0
         
         DO 10 I = 1, GLPTSDIV2
            DX1 = X1R*X(I)
C Sample point above the midpoint
            R1 = X1M+DX1
            Q1 = QFROMRHO(R1,MV1,GA1)
            R2F = RHO(MH-Q1,MV2,GA2)
            X2M=0.5D0*(R2F+R2I)
            X2R=0.5D0*(R2F-R2I)
            SLICE = 0.D0
            DO 20 J = 1, GLPTSDIV2
               DX2 = X2R*X(J)
               
               R2 = X2M+DX2
               Q2 = QFROMRHO(R2,MV2,GA2)
               SLICE = SLICE + W(J)
     .              * Q1**2/MV1**2 * Q2**2/MV2**2
     .              * GAMVV(SV,MH,Q1,Q2,C)
               
               R2 = X2M-DX2
               Q2 = QFROMRHO(R2,MV2,GA2)
               SLICE = SLICE + W(J)
     .              * Q1**2/MV1**2 * Q2**2/MV2**2
     .              * GAMVV(SV,MH,Q1,Q2,C)
 20         CONTINUE
            SLICE = X2R*SLICE
            GAMVVOFGL16 = GAMVVOFGL16 + W(I) * SLICE
C Sample point below the midpoint
            R1 = X1M-DX1
            Q1 = QFROMRHO(R1,MV1,GA1)
            R2F = RHO(MH-Q1,MV2,GA2)
            X2M=0.5D0*(R2F+R2I)
            X2R=0.5D0*(R2F-R2I)
            SLICE = 0.D0
            DO 30 J = 1, GLPTSDIV2
               DX2 = X2R*X(J)
               
               R2 = X2M+DX2
               Q2 = QFROMRHO(R2,MV2,GA2)
               SLICE = SLICE + W(J)
     .              * Q1**2/MV1**2 * Q2**2/MV2**2
     .              * GAMVV(SV,MH,Q1,Q2,C)
               
               R2 = X2M-DX2
               Q2 = QFROMRHO(R2,MV2,GA2)
               SLICE = SLICE + W(J)
     .              * Q1**2/MV1**2 * Q2**2/MV2**2
     .              * GAMVV(SV,MH,Q1,Q2,C)
 30         CONTINUE
            SLICE = X2R*SLICE
            GAMVVOFGL16 = GAMVVOFGL16 + W(I) * SLICE
 10      CONTINUE
C Normalize final result
         GAMVVOFGL16 = X1R*GAMVVOFGL16
      END
