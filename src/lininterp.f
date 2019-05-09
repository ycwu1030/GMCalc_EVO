* H.Logan 3/25/2003
**********************************************************
*                                                        *
      SUBROUTINE LININTERP(X1,Y1,X2,Y2,X,Y)
*                                                        *
* PERFORMS LINEAR INTERPOLATION BETWEEN THE INPUT POINTS *
* (X1,Y1) AND (X2,Y2) TO FIND THE Y VALUE AT POINT X.    *
*                                                        *
* INPUTS:                                                *
*   REFERENCE POINTS X1, Y1, X2, Y2                      *
*   DESIRED X VALUE -- MUST HAVE X1 < X < X2             *
*                                                        *
* OUTPUTS:                                               *
*   INTERPOLATED Y VALUE                                 *
*                                                        *
**********************************************************
      IMPLICIT NONE
      DOUBLE PRECISION X1,Y1,X2,Y2,X,Y
      IF (X1.LE.X.AND.X2.GE.X.AND.X1.NE.X2) THEN
         Y = Y1 + (X-X1)*(Y2-Y1)/(X2-X1)
      ELSE IF (X1.EQ.X2) THEN
         PRINT *, "** LININTERP: ERROR: X1 = X2"
         Y = 0.D0
      ELSE
         PRINT *, "** LININTERP: ERROR: X OUTSIDE INTERVAL (X1,X2)"
         Y = 0.D0
      ENDIF
      END
