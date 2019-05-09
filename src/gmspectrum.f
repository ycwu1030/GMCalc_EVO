
C==================================================================
C     Subroutine to calculate the physical masses & parameters
C==================================================================

      SUBROUTINE CALCPHYS
C Computes the physical masses, vevs, and the mixing angle between the two 
C custodial singlets, starting from the Lagrangian parameters.
C INPUTS: common block LPARAMS
C OUTPUTS: common block PHYSPARAMS
C          POSMSQOK flag
C
C MHL is the lighter custodial singlet mass and  MHH is the heavier 
C custodial singlet mass.
C The mixing angle alpha is defined according to
C    h = cos(alpha) H1 - sin(alpha) H1'
C    H = sin(alpha) H1 + cos(alpha) H1',
C where H1 is the doublet and H1' is the custodial singlet of the triplet.
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
      DOUBLE PRECISION ZETA, OMEGA, SIGMA, RHO
      COMMON/HYPER/ZETA,OMEGA,SIGMA,RHO
      DOUBLE PRECISION AA, BB, VV
      DIMENSION AA(6), BB(6), VV(6)
      COMMON/MINIMA/AA,BB,VV
C Local variables:
      DOUBLE PRECISION MH3SQ, MH5SQ, MHLSQ, MHHSQ
      DOUBLE PRECISION M1OVERVCHI
      DOUBLE PRECISION GOODV
      DOUBLE PRECISION VSM
      DOUBLE PRECISION MM11, MM12, MM22
      DOUBLE PRECISION SIN2A, COS2A
      INTEGER I
      DOUBLE PRECISION PI
      PI = 4.D0*DATAN(1.D0)

      POSMSQOK = 1

C First minimize the potential and find the vevs.  Use GETMINIMA.
      ZETA = 1.D0/3.D0
      OMEGA = 1.D0/2.D0
      SIGMA = DSQRT(3.D0)/4.D0
      RHO = 2.D0/DSQRT(3.D0)
      CALL GETMINIMA
C Debug!
c      DO 99 I = 1,3
c         PRINT *, "subroutine CALCPHYS: 3 extrema:"
c         PRINT *, "v_chi = ", BB(I)/DSQRT(3.D0)
c         PRINT *, "v_phi = ", AA(I)
c         PRINT *, "V = ", VV(I)
c 99   CONTINUE
C End debug block!
C Identify the deepest of the good extrema:
      VPHI = AA(1)
      VCHI = BB(1)/DSQRT(3.D0)
      GOODV = VV(1)
      IF (VV(2).LT.GOODV) THEN
         VPHI = AA(2)
         VCHI = BB(2)/DSQRT(3.D0)
         GOODV = VV(2)
      ELSE IF (VV(3).LT.GOODV) THEN
         VPHI = AA(3)
         VCHI = BB(3)/DSQRT(3.D0)
         GOODV = VV(3)
      ENDIF
      VSM = DSQRT(VPHI**2 + 8.D0*VCHI**2)

C Compute the masses and mixing angle alpha.
C Alpha is in the range (-pi/2, pi/2].
C For numerical safety when M1 and VCHI are both very small, use 
C a minimization-condition identity for M1/VCHI:
      M1OVERVCHI = 4.D0/VPHI**2 * (MU3SQ 
     .     + (2.D0*LAMBDA2 - LAMBDA5)*VPHI**2 
     .     + 4.D0*(LAMBDA3 + 3.D0*LAMBDA4)*VCHI**2
     .     - 6.D0*M2*VCHI)
      MH3SQ = (M1OVERVCHI/4.D0 + LAMBDA5/2.D0)*VSM**2
      MH5SQ = M1OVERVCHI/4.D0*VPHI**2 + 12.D0*M2*VCHI
     .     + 3.D0/2.D0*LAMBDA5*VPHI**2 + 8.D0*LAMBDA3*VCHI**2
      IF (MH3SQ.LT.0.D0.OR.MH5SQ.LT.0.D0) THEN
         POSMSQOK = 0
      ELSE
         MH3 = DSQRT(MH3SQ)
         MH5 = DSQRT(MH5SQ)
      ENDIF
C Mass-squared matrix elements for two custodial singlet states:
      MM11 = 8.D0 * LAMBDA1 * VPHI**2
      MM12 = DSQRT(3.D0)*VPHI 
     .     * (-M1/2.D0 + (4.D0*LAMBDA2 - 2.D0*LAMBDA5)*VCHI)
      MM22 = M1OVERVCHI*VPHI**2/4.D0 - 6.D0*M2*VCHI
     .     + 8.D0*VCHI**2*(LAMBDA3 + 3.D0*LAMBDA4)
      MHLSQ = 0.5D0*(MM11 + MM22 - DSQRT((MM11-MM22)**2 + 4.D0*MM12**2))
      MHHSQ = 0.5D0*(MM11 + MM22 + DSQRT((MM11-MM22)**2 + 4.D0*MM12**2))
      IF (MHLSQ.LT.0.D0.OR.MHHSQ.LT.0.D0) THEN
         POSMSQOK = 0
      ELSE
         MHL = DSQRT(MHLSQ)
         MHH = DSQRT(MHHSQ)
      ENDIF
C Mixing angle for two custodial singlet states:
C Arcsin will give 2alpha in the range [-pi/2,pi/2]:
      SIN2A = -2.D0*MM12/(MHLSQ-MHHSQ)
      ALPHA = 0.5D0*DASIN(SIN2A)
C To determine whether we ought to be in (-pi,-pi/2] or [pi/2,pi],
C check the sign of cos(2alpha):
      COS2A = (MM11-MM22)/(MHLSQ-MHHSQ)
      IF (COS2A.LT.0.D0) THEN
         IF (SIN2A.GE.0.D0) THEN
            ALPHA = PI/2.D0 - ALPHA
         ELSE
            ALPHA = -PI/2.D0 - ALPHA
         ENDIF
      ENDIF

      RETURN
      END

C==================================================================
C     Subroutines to check theoretical constraints on the
C     scalar potential parameters
C==================================================================

      SUBROUTINE THYCHECK
C Parent subroutine for performing the theory checks:
C perturbative unitarity of the lambda_i, 
C bounded-from-belowness of the scalar potential,
C and the absence of deeper alternative minima.
      IMPLICIT NONE
C Common blocks:
      INTEGER UNIOK, BFBOK, POSMSQOK, MINOK, INPUTOK
      COMMON/FLAGS/UNIOK, BFBOK, POSMSQOK, MINOK, INPUTOK

      CALL UNICHECK
      CALL BFBCHECK
C If the potential is not bounded from below, then we are definitely
C not going to be living in the proper EWSB minimum!
      IF (BFBOK.EQ.0) THEN
         MINOK = 0
      ELSE
         CALL MINCHECK
      ENDIF

      RETURN
      END


C======================================================================
      SUBROUTINE UNICHECK
C Checks the partial-wave unitarity conditions on the five lambda
C parameters of the scalar potential, using |Re(a_0)| \leq 1/2.
C Up to date as of Jan 3, 2014 (after comparison with Aoki & Kanemura, 
C arXiv:0712.4053).
C INPUTS: the 5 lambda parameters
C OUTPUTS: the flag UNIOK 
C          = 1 if the unitarity bound is satisfied 
C          = 0 if unitarity bound is violated
      IMPLICIT NONE
C Common blocks:
      DOUBLE PRECISION MU2SQ,MU3SQ,LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,
     .     LAMBDA5,M1,M2
      COMMON/LPARAMS/MU2SQ,MU3SQ,LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,
     .     LAMBDA5,M1,M2
      INTEGER UNIOK, BFBOK, POSMSQOK, MINOK, INPUTOK
      COMMON/FLAGS/UNIOK, BFBOK, POSMSQOK, MINOK, INPUTOK
C Local variables: 
      DOUBLE PRECISION X1P, X1M, X2P, X2M, Y1, Y2, Y3, Y4, Y5
      DOUBLE PRECISION PI
      DOUBLE PRECISION A, B, C, D

      PI = 4.D0*DATAN(1.D0)

      UNIOK = 1

C We use the same notation for the unitarity conditions as 
C Aoki & Kanemura, 0712.4053.
C The zeroth-partial-wave unitarity condition is given by 
C |Re a_0| > 1/2, or each of |x_i|, |y_i| < 8 pi.

      A = 12.D0*LAMBDA1 + 14.D0*LAMBDA3 + 22.D0*LAMBDA4
      B = 12.D0*LAMBDA1 - 14.D0*LAMBDA3 - 22.D0*LAMBDA4
      X1P = A + DSQRT(B**2 + 144.D0*LAMBDA2**2)
      X1M = A - DSQRT(B**2 + 144.D0*LAMBDA2**2)
      C = 4.D0*LAMBDA1 - 2.D0*LAMBDA3 + 4.D0*LAMBDA4
      D = 4.D0*LAMBDA1 + 2.D0*LAMBDA3 - 4.D0*LAMBDA4
      X2P = C + DSQRT(D**2 + 4.D0*LAMBDA5**2)
      X2M = C - DSQRT(D**2 + 4.D0*LAMBDA5**2)
      Y1 = 16.D0*LAMBDA3 + 8.D0*LAMBDA4
c      Y2 = 4.D0*LAMBDA3 + 8.D0*LAMBDA4
c      Y3 = 4.D0*LAMBDA2 - LAMBDA5
c      Y4 = 4.D0*LAMBDA2 + 2.D0*LAMBDA5
      Y5 = 4.D0*LAMBDA2 - 4.D0*LAMBDA5

C Note that the constraints on y_2, y_3, and y_4 do not provide any 
C additional information.  Therefore we do not compute them, for speed.

      IF (DABS(X1P).GT.8.D0*PI) THEN
         UNIOK = 0
      ELSE IF (DABS(X1M).GT.8.D0*PI) THEN
         UNIOK = 0
      ELSE IF (DABS(X2P).GT.8.D0*PI) THEN
         UNIOK = 0
      ELSE IF (DABS(X2M).GT.8.D0*PI) THEN
         UNIOK = 0
      ELSE IF (DABS(Y1).GT.8.D0*PI) THEN
         UNIOK = 0
      ELSE IF (DABS(Y5).GT.8.D0*PI) THEN
         UNIOK = 0
      ENDIF

      RETURN
      END


C=======================================================================
      SUBROUTINE BFBCHECK
C Checks the bounded-from-below conditions on the five lambda
C parameters of the scalar potential.
C Up to date as of Jan 3, 2014.
C INPUTS: the 5 lambda parameters
C OUTPUTS: the flag BFBOK
C          = 1 if the bounded-from-below conditions are satisfied
C          = 0 if the bounded-from-below conditions are violated
      IMPLICIT NONE
C Common blocks:
      DOUBLE PRECISION MU2SQ,MU3SQ,LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,
     .     LAMBDA5,M1,M2
      COMMON/LPARAMS/MU2SQ,MU3SQ,LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,
     .     LAMBDA5,M1,M2
      INTEGER UNIOK, BFBOK, POSMSQOK, MINOK, INPUTOK
      COMMON/FLAGS/UNIOK, BFBOK, POSMSQOK, MINOK, INPUTOK
C Local variables:
      DOUBLE PRECISION RHS
      DOUBLE PRECISION ZETA, OMEGAP, OMEGAM
      INTEGER I

      BFBOK = 1

C Condition 1: lambda_1 > 0
      IF (LAMBDA1.LT.0.D0) THEN
         BFBOK = 0
         GOTO 99
      ENDIF

C Condition 2: lambda_4 > -(1/3) lambda_3 (for lambda_3 >= 0)
C              lambda_4 > -lambda_3 (for lambda_3 < 0)
      IF (LAMBDA3.GE.0.D0) THEN
         IF (LAMBDA4.LT.(-LAMBDA3/3.D0)) THEN
            BFBOK = 0
            GOTO 99
         ENDIF
      ELSE
         IF (LAMBDA4.LT.(-LAMBDA3)) THEN
            BFBOK = 0
            GOTO 99
         ENDIF
      ENDIF

C Condition 3: 
C for lambda_5 >= 0 and lambda_3 >= 0:
C     lambda_2 > (1/2)lambda_5 - 2 sqrt(lambda_1(lambda_3/3 + lambda4))
C for lambda_5 >= 0 and lambda_3 < 0: have to scan over zeta
C     lambda_2 > omega+ lambda_5 - 2 sqrt(lambda_1(zeta*lambda_3 + lambda_4))
C for lambda_5 < 0: have to scan over zeta
C     lambda_2 > omega- lambda_5 - 2 sqrt(lambda_1(zeta*lambda_3 + lambda_4))
      IF (LAMBDA5.GE.0.D0) THEN
         IF (LAMBDA3.GE.0.D0) THEN
            RHS = LAMBDA5/2.D0 
     .           - 2.D0*DSQRT(LAMBDA1*(LAMBDA3/3.D0 + LAMBDA4))
            IF (LAMBDA2.LT.RHS) THEN
               BFBOK = 0
            ENDIF
         ELSE
C Case 2: lambda_5 >= 0 but lambda_3 < 0.  Scan over zeta in (1/3, 1).
            DO 10 I=0,1000
               ZETA = 1.D0/3.D0 + I*2.D0/3.D0/1000.D0
               RHS = OMEGAP(ZETA)*LAMBDA5 
     .              - 2.D0*DSQRT(LAMBDA1*(ZETA*LAMBDA3 + LAMBDA4))
               IF (LAMBDA2.LT.RHS) THEN
                  BFBOK = 0
                  GOTO 99
               ENDIF
 10         CONTINUE
         ENDIF
      ELSE
C Case 3: lambda_5 < 0.  Scan over zeta in (1/3, 1).
         DO 20 I=0,1000
            ZETA = 1.D0/3.D0 + I*2.D0/3.D0/1000.D0
            RHS = OMEGAM(ZETA)*LAMBDA5
     .           - 2.D0*DSQRT(LAMBDA1*(ZETA*LAMBDA3 + LAMBDA4))
            IF (LAMBDA2.LT.RHS) THEN
               BFBOK = 0
               GOTO 99
            ENDIF
 20      CONTINUE
      ENDIF

 99   RETURN
      END

C Functions omega+ and omega- of zeta used in the bounded-from-below
C conditions on lambda_2: 
C These two functions are called by the subroutine BFBCHECK.
      DOUBLE PRECISION FUNCTION OMEGAP(ZETA)
      IMPLICIT NONE
      DOUBLE PRECISION ZETA
      DOUBLE PRECISION B
      B = DSQRT(3.D0/2.D0*(ZETA - 1.D0/3.D0))
      OMEGAP = (1.D0-B)/6.D0 + DSQRT(2.D0)/3.D0 
     .     * DSQRT((1.D0-B)*(0.5D0+B))
      END

      DOUBLE PRECISION FUNCTION OMEGAM(ZETA)
      IMPLICIT NONE
      DOUBLE PRECISION ZETA
      DOUBLE PRECISION B
      B = DSQRT(3.D0/2.D0*(ZETA - 1.D0/3.D0))
      OMEGAM = (1.D0-B)/6.D0 - DSQRT(2.D0)/3.D0 
     .     * DSQRT((1.D0-B)*(0.5D0+B))
      END

C======================================================================
      SUBROUTINE MINCHECK
C Checks for the presence of deeper alternative minima.  This subroutine
C only gets called if the potential is bounded from below.
C INPUTS: the 5 lambda parameters, mu2sq, mu3sq, M1, and M2.
C OUTPUTS: the flag MINOK
C          = 1 if there are no deeper alternative minima
C          = 0 if there are deeper alternative minima (bad!)
      IMPLICIT NONE
C Common blocks:
      DOUBLE PRECISION MU2SQ,MU3SQ,LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,
     .     LAMBDA5,M1,M2
      COMMON/LPARAMS/MU2SQ,MU3SQ,LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,
     .     LAMBDA5,M1,M2
      INTEGER UNIOK, BFBOK, POSMSQOK, MINOK, INPUTOK
      COMMON/FLAGS/UNIOK, BFBOK, POSMSQOK, MINOK, INPUTOK
      DOUBLE PRECISION ZETA, OMEGA, SIGMA, RHO
      COMMON/HYPER/ZETA,OMEGA,SIGMA,RHO
      DOUBLE PRECISION AA, BB, VV
      DIMENSION AA(6), BB(6), VV(6)
      COMMON/MINIMA/AA,BB,VV
C Local variables:
      DOUBLE PRECISION THETA
      DOUBLE PRECISION FZETA, FOMEGA, FSIGMA, FRHO
      DOUBLE PRECISION GOODV
      INTEGER I, J
      DOUBLE PRECISION PI

      PI = 4.D0*DATAN(1.D0)

      MINOK = 1

C Compute the 6 (potential) extrema for theta = a and identify the deepest one:
      THETA = DACOS(1.D0/DSQRT(3.D0))
      ZETA = FZETA(THETA)
      OMEGA = FOMEGA(THETA)
      SIGMA = FSIGMA(THETA)
      RHO = FRHO(THETA)
      CALL GETMINIMA

C Identify the good vacuum and record its depth:
      GOODV = VV(1)
      IF (VV(2).LT.GOODV) THEN
         GOODV = VV(2)
      ELSE IF (VV(3).LT.GOODV) THEN
         GOODV = VV(3)
      ENDIF
c      PRINT *, "Good vacuum: V = ", GOODV
c      PRINT *

C Check the depths of the (bad) a=0 vacua:
      DO 10 J=4,6
         IF (VV(J).LT.GOODV) THEN
            MINOK = 0
c            PRINT *, "Bad vacuum found!", AA(J), BB(J), VV(J)
            GOTO 99
         ENDIF
 10   CONTINUE

C If lambda3, lambda5, M1, and M2 are all positive, then the desired 
C vacuum at theta = a = acos(1/sqrt(3)) is the global minimum.
C If lambda3 and lambda5 are positive and M1 and M2 are negative, then
C the other desired vacuum at theta = pi + a is the global minimum.
C We have already checked them both because we include negative values
C of the triplet vev parameter b.
      IF (LAMBDA3.GT.0.D0.AND.LAMBDA5.GT.0.D0.AND.M1*M2.GT.0.D0) THEN
         GOTO 99
      ENDIF

C Scan over theta for the other cases and check all the minima.
C Avoid stepping on theta = a or theta = pi + a, which we have already
C checked.
C This do loop covers a --> pi + a
      DO 20 I=1,999
         THETA = DACOS(1.D0/DSQRT(3.D0)) + PI*I/1000.D0
         ZETA = FZETA(THETA)
         OMEGA = FOMEGA(THETA)
         SIGMA = FSIGMA(THETA)
         RHO = FRHO(THETA)
         CALL GETMINIMA
         DO 21 J=1,6
            IF(VV(J).LT.GOODV) THEN
               MINOK = 0
c               PRINT *, "Bad vacuum found!", THETA/PI, VV(J)
               GOTO 99
            ENDIF
 21      CONTINUE
 20   CONTINUE
C This do loop covers pi + a --> 2pi + a
      DO 30 I=1,999
         THETA = PI + DACOS(1.D0/DSQRT(3.D0)) + PI*I/1000.D0
         ZETA = FZETA(THETA)
         OMEGA = FOMEGA(THETA)
         SIGMA = FSIGMA(THETA)
         RHO = FRHO(THETA)
         CALL GETMINIMA
         DO 31 J=1,6
            IF(VV(J).LT.GOODV) THEN
               MINOK = 0
c               PRINT *, "Bad vacuum found!", THETA/PI, VV(J)
               GOTO 99
            ENDIF
 31      CONTINUE
 30   CONTINUE

 99   RETURN
      END

C These four functions are used by the subroutine MINCHECK.
      DOUBLE PRECISION FUNCTION FZETA(THETA)
      IMPLICIT NONE
      DOUBLE PRECISION THETA
      FZETA = 0.5D0*DSIN(THETA)**4 + DCOS(THETA)**4
      END

      DOUBLE PRECISION FUNCTION FOMEGA(THETA)
      IMPLICIT NONE
      DOUBLE PRECISION THETA
      FOMEGA = 0.25D0*DSIN(THETA)**2 
     .     + 1.D0/DSQRT(2.D0)*DSIN(THETA)*DCOS(THETA)
      END

      DOUBLE PRECISION FUNCTION FSIGMA(THETA)
      IMPLICIT NONE
      DOUBLE PRECISION THETA
      FSIGMA = 1.D0/2.D0/DSQRT(2.D0)*DSIN(THETA) 
     .     + 0.25D0*DCOS(THETA)
      END
      
      DOUBLE PRECISION FUNCTION FRHO(THETA)
      IMPLICIT NONE
      DOUBLE PRECISION THETA
      FRHO = 3.D0*DSIN(THETA)**2*DCOS(THETA)
      END

C===================================================================
      SUBROUTINE GETMINIMA
C Computes the field values and potential at each of the 6 possible
C minima, for given fixed values of zeta, omega, sigma, & rho.
C The field values are defined as:
C     a^2 = Tr(Phi^dagger Phi) = v_phi^2
C     b^2 = Tr(X^dagger X) = 2 v_chi^2 + v_xi^2 
C                          = 3 v_chi^2 in the desired vacuum.
C INPUTS: the 5 lambda parameters, mu2sq, mu3sq, M1, and M2;
C         and the values of zeta, omega, sigma, and rho.
C OUTPUT: three arrays with 6 entries each, containing the solutions
C for a, b, and V.  Only the first three entries can be good vacua.
      IMPLICIT NONE
C Common blocks:
      DOUBLE PRECISION MU2SQ,MU3SQ,LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,
     .     LAMBDA5,M1,M2
      COMMON/LPARAMS/MU2SQ,MU3SQ,LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,
     .     LAMBDA5,M1,M2
      DOUBLE PRECISION ZETA, OMEGA, SIGMA, RHO
      COMMON/HYPER/ZETA,OMEGA,SIGMA,RHO
      DOUBLE PRECISION AA, BB, VV
      DIMENSION AA(6), BB(6), VV(6)
      COMMON/MINIMA/AA,BB,VV
C Local variables:
      DOUBLE PRECISION CALCV
      DOUBLE PRECISION ALPHA, BETA, GAMMA, DELTA, RDISC, ASQ
      DOUBLE COMPLEX DISC, DISC0, DISC1
      DOUBLE COMPLEX II, U1, U2, U3, CROOT, BC1, BC2, BC3
      INTEGER J

      II = (0.D0,1.D0)
      U1 = (1.D0,0.D0)
      U2 = -0.5D0 + II*DSQRT(3.D0)/2.D0
      U3 = -0.5D0 - II*DSQRT(3.D0)/2.D0

C Solve for the values of a and b at the extrema: 

C 1) a != 0 case; cubic equation for b.
      ALPHA = 4.D0*(ZETA*LAMBDA3 + LAMBDA4)
     .     - 1.D0/LAMBDA1 * (LAMBDA2 - OMEGA*LAMBDA5)**2
      BETA = -3.D0*M2*RHO + 3.D0/2.D0/LAMBDA1 * M1*SIGMA 
     .     * (LAMBDA2 - OMEGA*LAMBDA5)
      GAMMA = MU3SQ - 1.D0/2.D0/LAMBDA1 * MU2SQ 
     .     * (LAMBDA2 - OMEGA*LAMBDA5)
     .     - 1.D0/2.D0/LAMBDA1 * M1**2 * SIGMA**2
      DELTA = 1.D0/4.D0/LAMBDA1 * MU2SQ * M1 * SIGMA
      DISC = 18.D0 * ALPHA*BETA*GAMMA*DELTA - 4.D0 * BETA**3 *DELTA
     .     + BETA**2 * GAMMA**2 - 4.D0 * ALPHA * GAMMA**3
     .     - 27.D0 * ALPHA**2 * DELTA**2
      DISC0 = BETA**2 - 3.D0 * ALPHA * GAMMA
      DISC1 = 2.D0 * BETA**3 - 9.D0 * ALPHA*BETA*GAMMA 
     .     + 27.D0 * ALPHA**2 * DELTA
      CROOT = (0.5D0 * DISC1 + 0.5D0*CDSQRT(-27.D0*ALPHA**2*DISC)
     .     )**(1.D0/3.D0)
C Handling the special case DISC0 = 0:
C Thanks to R. Ruiz for pointing out this bug. Fixed 2015/01/20.
      IF (CDABS(4.D0*DISC0**3/DISC1**2).LT.1D-12) THEN
         CROOT = (0.5D0 * DISC1 + 0.5D0 * DISC1)**(1.D0/3.D0)
      ENDIF
      BC1 = -1.D0/3.D0/ALPHA * (BETA + U1*CROOT + DISC0/U1/CROOT)
      BC2 = -1.D0/3.D0/ALPHA * (BETA + U2*CROOT + DISC0/U2/CROOT)
      BC3 = -1.D0/3.D0/ALPHA * (BETA + U3*CROOT + DISC0/U3/CROOT)
      
      RDISC = DREAL(DISC)
      IF (RDISC.GE.0.D0) THEN
C The 3 roots are all real (if RDISC = 0 then two are the same)
         IF (M1.EQ.0.D0.AND.M2.EQ.0.D0) THEN
C Z2 is preserved; use the simpler quadratic formula.
            BB(1) = 0.D0
            BB(2) = DSQRT(-GAMMA/ALPHA)
            BB(3) = -DSQRT(-GAMMA/ALPHA)
         ELSE
            BB(1) = DREAL(BC1)
            BB(2) = DREAL(BC2)
            BB(3) = DREAL(BC3)
         ENDIF
         DO 10 J=1,3
            ASQ = 1.D0/4.D0/LAMBDA1
     .           * (-MU2SQ + 2.D0*M1*SIGMA*BB(J)
     .           - 2.D0*(LAMBDA2 - OMEGA*LAMBDA5) * BB(J)**2)
C "a" must also be real, or this is not a solution.
            IF (ASQ.GE.0.D0) THEN
               AA(J) = DSQRT(ASQ)
            ELSE
               AA(J) = 0.D0
               BB(J) = 0.D0
            ENDIF
            VV(J) = CALCV(AA(J),BB(J))
 10      CONTINUE
      ELSE
C Discriminant < 0, so only one root is real.  First determine which one it is.
         IF (DABS(DIMAG(BC1)).LT.1D-6*CDABS(BC1)) THEN
            BB(1) = DREAL(BC1)
         ELSE IF (DABS(DIMAG(BC2)).LT.1D-6*CDABS(BC2)) THEN
            BB(1) = DREAL(BC2)
         ELSE
            BB(1) = DREAL(BC3)
         ENDIF
         AA(1) = DSQRT(1.D0/4.D0/LAMBDA1
     .           * (-MU2SQ + 2.D0*M1*SIGMA*BB(1)
     .           - 2.D0*(LAMBDA2 - OMEGA*LAMBDA5) * BB(1)**2))
         VV(1) = CALCV(AA(1),BB(1))
         DO 20 J=2,3
            BB(J) = 0.D0
            AA(J) = 0.D0
            VV(J) = CALCV(AA(J),BB(J))
 20      CONTINUE
      ENDIF

C 2) a = 0 and b = 0 case (bad vacuum)
      AA(4) = 0.D0
      BB(4) = 0.D0
      VV(4) = CALCV(AA(4),BB(4))

C 3) a = 0 and b != 0 case: quadratic equation for b (bad vacua)
      ALPHA = 4.D0*(ZETA*LAMBDA3 + LAMBDA4)
      BETA = -3.D0*M2*RHO
      GAMMA = MU3SQ
      RDISC = BETA**2 - 4.D0*ALPHA*GAMMA
      IF (RDISC.GE.0.D0) THEN
C The quadratic equation has two real roots (equal if RDISC = 0).
         BB(5) = 1.D0/2.D0/ALPHA * (-BETA + DSQRT(RDISC))
         BB(6) = 1.D0/2.D0/ALPHA * (-BETA - DSQRT(RDISC))
      ELSE
C The quadratic equation does not have any real roots.
         BB(5) = 0.D0
         BB(6) = 0.D0
      ENDIF
      DO 30 J=5,6
         AA(J) = 0.D0
         VV(J) = CALCV(AA(J),BB(J))
 30   CONTINUE

      RETURN
      END

C This function is used by the subroutine GETMINIMA.
      DOUBLE PRECISION FUNCTION CALCV(A,B)
C Computes the value of the potential given a, b, the Lagrangian
C parameters, zeta, omega, sigma, and rho.
C INPUTS: the 5 lambda parameters, mu2sq, mu3sq, M1, and M2;
C         the values of zeta, omega, sigma, and rho;
C         and the values of a and b at the extremum of interest.
C OUTPUT: the value of the potential at this point.
      IMPLICIT NONE
C Inputs:
      DOUBLE PRECISION A, B
C Common blocks:
      DOUBLE PRECISION MU2SQ,MU3SQ,LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,
     .     LAMBDA5,M1,M2
      COMMON/LPARAMS/MU2SQ,MU3SQ,LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,
     .     LAMBDA5,M1,M2
      DOUBLE PRECISION ZETA, OMEGA, SIGMA, RHO
      COMMON/HYPER/ZETA,OMEGA,SIGMA,RHO
      CALCV = MU2SQ/2.D0 * A**2 + MU3SQ/2.D0 * B**2
     .     + LAMBDA1 * A**4 + LAMBDA2 * A**2 * B**2
     .     + LAMBDA3 * ZETA * B**4 + LAMBDA4 * B**4
     .     - LAMBDA5 * OMEGA * A**2 * B**2
     .     - M1 * SIGMA * A**2 * B - M2 * RHO * B**3
      END

