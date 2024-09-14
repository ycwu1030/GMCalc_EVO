      SUBROUTINE GET_MINMASS(MINMASS)
C Returns the mass of the lightest *new* particle.  If HH turns out to
C be the custodial singlet whose mass gets set to the measured MH, then
C HL is considered a *new* particle.
      IMPLICIT NONE
      DOUBLE PRECISION MINMASS
C Common blocks:
      DOUBLE PRECISION MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
      COMMON/PHYSPARAMS/MHL,MHH,MH3,MH5,ALPHA,VPHI,VCHI
      DOUBLE PRECISION MH
      INTEGER INPUTSET,Z2
      COMMON/INPUT/MH,INPUTSET,Z2
C Local variables:
      DOUBLE PRECISION EPS

      EPS = 1.D-2

      IF (MHH.LT.MH3.AND.MHH.LT.MH5) THEN
         MINMASS = MHH
      ELSE IF (MH3.LT.MHH.AND.MH3.LT.MH5) THEN
         MINMASS = MH3
      ELSE
         MINMASS = MH5
      ENDIF
      IF (MHL.LT.MH-EPS) THEN
         IF (MHL.LT.MH3.AND.MHL.LT.MH5) THEN
            MINMASS = MHL
         ELSE IF (MH3.LT.MHL.AND.MH3.LT.MH5) THEN
            MINMASS = MH3
         ELSE
            MINMASS = MH5
         ENDIF
      ENDIF

      RETURN
      END


      FUNCTION RANDOM_RANGE (LBOUND,UBOUND)

      IMPLICIT NONE
      DOUBLE PRECISION LBOUND,UBOUND
      DOUBLE PRECISION LEN
      DOUBLE PRECISION RANDOM_RANGE
      DOUBLE PRECISION T


      LEN=UBOUND-LBOUND
      CALL RANDOM_NUMBER(T)
      RANDOM_RANGE=LBOUND+LEN*T

      return
      end

      SUBROUTINE INIT_RANDOM_SEED()
          INTEGER :: i, n, clock
          INTEGER, DIMENSION(:), ALLOCATABLE :: seed
          CALL RANDOM_SEED(size=n)
          ALLOCATE(seed(n))
          CALL SYSTEM_CLOCK(COUNT=clock)
          seed = clock + 37 * (/ (i-1, i=1, n) /)
          CALL RANDOM_SEED(PUT=seed)
          DEALLOCATE(seed)
      END SUBROUTINE
