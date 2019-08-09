C ** Correction report.
C ** Correction 1. 15/06/93: 5 lines corrected **
C ** Correction 2. 13/01/94: 2 lines added
C ** End of Correction report.
C  THIS VERSION: 13/01/1994 AT 04:41:59 PM.
C ** Correction 1. 15/06/93: 5 lines corrected **
CS    SUBROUTINE SBNDFA( N, NSEMIB, DIAG, OFFDIA, LOFFDI )
CD    SUBROUTINE DBNDFA( N, NSEMIB, DIAG, OFFDIA, LOFFDI )
      INTEGER N, NSEMIB, LOFFDI
CS    REAL             DIAG( N ), OFFDIA( LOFFDI, N )
CD    DOUBLE PRECISION DIAG( N ), OFFDIA( LOFFDI, N )
C ** Correction 1. 15/06/93:  end of correction **
C
C  ***************************************************************
C
C  COMPUTE THE L D L(TRANSPOSE) FACTORIZATION OF A BANDED MATRIX.
C
C  NICK GOULD, 17TH OCTOBER, 1990.
C  FOR CGT PRODUCTIONS.
C
C  ***************************************************************
C
      INTEGER          I,      IPJM1,  J,      K,      M
CS    REAL             ZERO,   OFFD,   TAU1,   TAU2,   GAMMA,  OFFSUM,
CD    DOUBLE PRECISION ZERO,   OFFD,   TAU1,   TAU2,   GAMMA,  OFFSUM,
     *                 ONE
      LOGICAL          PHASE1
      INTRINSIC        MIN,    MAX,    ABS
C
C  COMMON VARIABLES.
C
CS    REAL             EPSMCH, EPSNEG, TINY, BIG
CD    DOUBLE PRECISION EPSMCH, EPSNEG, TINY, BIG
CS    COMMON / SMACHN / EPSMCH, EPSNEG, TINY, BIG
CD    COMMON / DMACHN / EPSMCH, EPSNEG, TINY, BIG
C ** Correction 2. 13/01/94: 2 lines added
CS    SAVE   / SMACHN /
CD    SAVE   / DMACHN /
C ** Correction 2. 13/01/94: end of correction **
C
C  SET CONSTANT VALUES.
C
CS    PARAMETER ( ZERO = 0.0E+0, ONE = 1.0E+0 )
CD    PARAMETER ( ZERO = 0.0D+0, ONE = 1.0E+0 )
      PHASE1 = .TRUE.
      GAMMA  = ONE
C
C  CHECK THAT THE INITIAL DIAGONALS ARE POSITIVE.
C
      DO 10 I = 1, N
C
C  THE MATRIX IS INDEFINITE. ENTER PHASE-2.
C
         IF ( DIAG( I ) .LE. ZERO ) PHASE1 = .FALSE.
C
C   FIND THE LARGEST DIAGONAL ENTRY.
C
         GAMMA = MAX( GAMMA, ABS( DIAG( I ) ) ) 
   10 CONTINUE
C
C  SET PIVOT TOLERANCES.
C
      TAU1 = GAMMA * EPSMCH ** 0.333
      TAU2 = TAU1
C
C  LOOP OVER THE COLUMNS OF THE MATRIX.
C
      DO 500 I = 1, N
         M = MIN( NSEMIB, N - I )
C
C  PERFORM THE I-TH ELIMINATION. FIRST, CHECK THAT THE RESULTING
C  PIVOTS ARE POSITIVE.
C
         IF ( PHASE1 ) THEN
            DO 100 J = 1, M
               IF ( DIAG( I + J ) - ( OFFDIA( J, I ) / 
     *              DIAG( I ) ) * OFFDIA( J, I )  .LE. TAU1 ) THEN
C
C  THE MATRIX IS INDEFINITE. ENTER PHASE-2.
C
                  PHASE1 = .FALSE.
                  GO TO 110
               END IF
  100       CONTINUE   
         END IF
C
C  IF THE MATRIX IS INDEFINITE, MODIFY THE DIAGONALS.
C
  110    CONTINUE
C
C  COMPUTE THE GERSHGORIN RADIUS FOR THE PRINCIPAL DIAGONAL OF THE
C  LAST N - I BY N - I SUBMATRIX.
C
         IF ( .NOT. PHASE1 ) THEN
            OFFSUM = ZERO
            DO 120 J = 1, M
               OFFSUM = OFFSUM + ABS( OFFDIA( J, I ) )
  120       CONTINUE
C
C  PERTURB THE DIAGONAL SO THAT THE GERSHGORIN DISK LIES IN THE
C  POSITIVE HALF OF THE COMPLEX PLANE.
C
            OFFSUM = MAX( ZERO, - DIAG( I ) + MAX( OFFSUM, TAU2 ) ) 
C           WRITE(6,*) ' DIAGONAL ', I, ' MODIFIED BY ', OFFSUM
            DIAG( I ) = DIAG( I ) + OFFSUM
         END IF
C
C  PERFORM THE I-TH STEP OF THE FACTORIZATION.
C
         DO 230 J = 1, M
            OFFD  = OFFDIA( J, I ) 
C
C  UPDATE THE SCHUR COMPLEMENT. (1) OFF DIAGONAL TERMS.
C
            IPJM1    = J
            DO 220 K = 1, J - 1
               OFFDIA( IPJM1 - K, I + K ) = 
     *         OFFDIA( IPJM1 - K, I + K ) - OFFD * OFFDIA( K, I ) 
  220       CONTINUE
C
C  (2) DIAGONAL TERMS.
C
            OFFD          = OFFD / DIAG ( I )
            DIAG( I + J ) = DIAG( I + J ) - OFFD * OFFDIA( J, I ) 
C
C  FIND THE SUBDIAGONAL OF THE I-TH COLUMN OF THE FACTOR L.
C
            OFFDIA( J, I ) = OFFD 
  230    CONTINUE   
  500 CONTINUE   
      RETURN
C
C  END OF SUBROUTINE BNDFA.
C
      END
C ** Correction report.
C ** Correction 1. 15/06/93: 5 lines corrected **
C ** End of Correction report.
C  THIS VERSION: 13/01/1994 AT 04:41:59 PM.
C ** Correction 1. 15/06/93: 5 lines corrected **
CS    SUBROUTINE SBNDSL( N, NSEMIB, DIAG, OFFDIA, LOFFDI, RHS )
CD    SUBROUTINE DBNDSL( N, NSEMIB, DIAG, OFFDIA, LOFFDI, RHS )
      INTEGER N, NSEMIB, LOFFDI
CS    REAL             DIAG( N ), OFFDIA( LOFFDI, N ), RHS( N )
CD    DOUBLE PRECISION DIAG( N ), OFFDIA( LOFFDI, N ), RHS( N )
C ** Correction 1. 15/06/93:  end of correction **
C
C  ***************************************************************
C
C  SOLVE THE SYSTEM OF LINEAR EQUATIONS 
C
C  L D L(TRANSPOSE) X = RHS,
C
C  PUTTING THE SOLUTION IN RHS.
C
C  NICK GOULD, 17TH OCTOBER, 1990.
C  FOR CGT PRODUCTIONS.
C
C  ***************************************************************
C
      INTEGER          I,      J,     M
      INTRINSIC        MIN
C
C  FORWARD SOLVE TO OBTAIN THE SOLUTION TO L Y = RHS, PUTTING THE
C  SOLUTION IN RHS.
C
      DO 100 I    = 1, N
         M        = MIN( NSEMIB, N - I )
         DO 20 J  = 1, M
            RHS( I + J ) = RHS( I + J ) - OFFDIA( J, I ) * RHS( I )
   20    CONTINUE   
C
C  OBTAIN THE SOLUTION TO THE DIAGONAL SYSTEM D Y = RHS, PUTTING THE
C  SOLUTION IN RHS.
C
         RHS( I ) = RHS( I ) / DIAG( I )
  100 CONTINUE   
C
C  BACK SOLVE TO OBTAIN THE SOLUTION TO L(TRANSPOSE) Y = RHS, 
C  PUTTING THE SOLUTION IN RHS.
C
      DO 300 I    = N, 1, - 1
         M        = MIN( NSEMIB, N - I )
         DO 220 J = 1, M
            RHS( I ) = RHS( I  ) - OFFDIA( J, I ) * RHS( I + J )
  220    CONTINUE       
  300 CONTINUE   
      RETURN
C
C  END OF SUBROUTINE BNDSL.
C
      END
