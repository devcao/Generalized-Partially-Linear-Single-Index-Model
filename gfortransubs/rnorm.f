      SUBROUTINE RNORM(U1, U2)
C
C     ALGORITHM AS 53.1  APPL. STATIST. (1972) VOL.21, NO.3
C
C     Sets U1 and U2 to two independent standardized random normal
C     deviates.   This is a Fortran version of the method given in
C     Knuth(1969).
C
C     Function RAN must give a result randomly and rectangularly
C     distributed between the limits 0 and 1 exclusive.
C
      REAL U1, U2
      REAL RAND
C
C     Local variables
C
      REAL X, Y, S, ONE, TWO
      DATA ONE /1.0/, TWO /2.0/
C
    1 X = RAND(0)
      Y = RAND(0)
      X = TWO * X - ONE
      Y = TWO * Y - ONE
      S = X * X + Y * Y
      IF (S .GT. ONE) GO TO 1
      S = SQRT(- TWO * LOG(S) / S)
      U1 = X * S
      U2 = Y * S
      RETURN
      END
