C Dummy routines for use without LoopTools

C Subroutines to initialize and end LoopTools
C (These are defined to facilitate switching off the LoopTools use by instead
C linking heteroloops-dummy.F.)

      SUBROUTINE LTSTARTER
      PRINT *
      PRINT *, "*************************************************",
     .     "******************"
      PRINT *, "* Note: LoopTools is not linked.                 ",
     .     "                 *"
      PRINT *, "* Decays of H3+ -> W+GA, H50 -> ZGA, H5+ -> W+GA ",
     .     "are set to zero. *"
      PRINT *, "*************************************************",
     .     "******************"
      PRINT *
      RETURN
      END

      SUBROUTINE LTENDER
      RETURN
      END

C Dummy routines for non-standard loop-induced decays

      DOUBLE PRECISION FUNCTION HETLOOPH5WGA()
      HETLOOPH5WGA = 0.D0
      RETURN
      END

      DOUBLE PRECISION FUNCTION HETLOOPH3WGA()
      HETLOOPH3WGA = 0.D0
      RETURN
      END

      DOUBLE PRECISION FUNCTION HETLOOPH5ZGA()
      HETLOOPH5ZGA = 0.D0
      RETURN
      END

