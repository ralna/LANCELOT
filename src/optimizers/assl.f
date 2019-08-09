C  THIS VERSION: 14/08/1991 AT 12:00:44 PM.
CS    SUBROUTINE SASSLB( N, M, ICLASS, LBD, BD, IRNBD, IPBD, LCD,
CD    SUBROUTINE DASSLB( N, M, ICLASS, LBD, BD, IRNBD, IPBD, LCD,
     *                   CD, JCNCD, IPCD, R, LR, Q, LQ, RHS, X,
     *                   VECTOR, INFORM )
C
C  SOLVE THE SYSTEM OF EQUATIONS
C
C     / A  B \ / X1 \ _ / RHS1 \
C     \ C  D / \ X2 / - \ RHS2 /
C
C  NICK GOULD, AUGUST 3RD 1988.
C  FOR CGT PRODUCTIONS.
C
      INTEGER          INFORM, N, M, ICLASS, LBD, LCD, LR, LQ
      INTEGER          IRNBD( LBD ), JCNCD( LCD ), IPBD( * ), IPCD( * )
CS    REAL             BD( LBD ), CD( LCD ), RHS( * ), X( * ),
CD    DOUBLE PRECISION BD( LBD ), CD( LCD ), RHS( * ), X( * ),
     *                 R( LR ), Q( LQ, * ), VECTOR( * )
C
C  LOCAL VARIABLES.
C
      INTEGER          I, ICOL, IROW, J, K, K1, K2
CS    REAL             SCALAR, ZERO, ONE
CD    DOUBLE PRECISION SCALAR, ZERO, ONE
CS    EXTERNAL         SAXPY, SCOPY, SASSLF
CD    EXTERNAL         DAXPY, DCOPY, DASSLF
      SAVE
C
C  DATA STATEMENTS.
C
CS    DATA ZERO, ONE / 0.0E+0, 1.0E+0 /
CD    DATA ZERO, ONE / 0.0D+0, 1.0D+0 /
      IF ( INFORM .LE. 0 ) THEN
         INFORM = - 5
         RETURN
      END IF
      GO TO ( 100, 200, 300 ), INFORM
C
C  COMPUTE THE PART OF THE SOLUTION X2.
C
  100 CONTINUE
      DO 110 I      = 1, M
         X( N + I ) = ZERO
  110 CONTINUE
      IF ( N .GT. 0 ) THEN
CS       CALL SCOPY( N, RHS, 1, VECTOR, 1 )
CD       CALL DCOPY( N, RHS, 1, VECTOR, 1 )
         INFORM = 2
C
C  RETURN TO GET THE PRODUCT A(INVERSE) * VECTOR.
C
         RETURN
      END IF
  200 CONTINUE
C
C  COPY VECTOR ONTO X1.
C
CS    CALL SCOPY( N, VECTOR, 1, X, 1 )
CD    CALL DCOPY( N, VECTOR, 1, X, 1 )
      IF ( M .LE. 0 ) GO TO 400
C
C  TRANSFORM THE RIGHT HAND SIDE.
C
      IF ( ICLASS .EQ. 4 ) THEN
C
C  S IS NEGATIVE DEFINITE AND THE FACTORIZATION OF - S IS USED.
C
         DO 220 IROW = 1, M
            SCALAR   = - RHS( N + IROW )
            IF ( N .GT. 0 ) THEN
               K1 = IPBD( IROW )
               K2 = IPBD( IROW + 1 ) - 1
               IF ( K1 .LE. K2 ) THEN
                  DO 210 K = K1, K2
                     J     = IRNBD( K )
                     IF ( J .LE. N )
     *                  SCALAR = SCALAR + BD( K ) * VECTOR( J )
  210             CONTINUE
               END IF
            END IF
            X( N + IROW ) = SCALAR
  220    CONTINUE
      ELSE
C
C  THE FACTORIZATION OF S IS USED.
C
         DO 250 IROW = 1, M
            SCALAR   = RHS( N + IROW )
            IF ( N .GT. 0 ) THEN
C
C  S IS SYMMETRIC.
C
               IF ( ICLASS .GT. 1 ) THEN
                  K1 = IPBD( IROW )
                  K2 = IPBD( IROW + 1 ) - 1
                  IF ( K1 .LE. K2 ) THEN
                     DO 230 K = K1, K2
                        J     = IRNBD( K )
                        IF ( J .LE. N )
     *                     SCALAR = SCALAR - BD( K ) * VECTOR( J )
  230                CONTINUE
                  END IF
               ELSE
C
C  S IS UNSYMMETRIC.
C
                  K1 = IPCD( IROW )
                  K2 = IPCD( IROW + 1 ) - 1
                  IF ( K1 .LE. K2 ) THEN
                     DO 240 K = K1, K2
                        J     = JCNCD( K )
                        IF ( J .LE. N )
     *                     SCALAR = SCALAR - CD( K ) * VECTOR( J )
  240                CONTINUE
                  END IF
               END IF
            END IF
C
C  IF A QR FACTORIZATION IS USED, TRANSFORM BY Q(TRANSPOSE).
C
            IF ( ICLASS .LE. 2 ) THEN
CS             CALL SAXPY( M, SCALAR, Q( IROW, 1 ), LQ, X( N + 1 ), 1 )
CD             CALL DAXPY( M, SCALAR, Q( IROW, 1 ), LQ, X( N + 1 ), 1 )
            ELSE
               X( N + IROW ) = SCALAR
            END IF
  250    CONTINUE
      END IF
C
C  IF AN R(TRANSPOSE) R FACTORIZATION IS USED, TRANSFORM
C  BY R(INVERSE_TRANSPOSE).
C
      IF ( ICLASS .GE. 3 ) THEN
CS       CALL SASSLF( M, R, LR, X( N + 1 ), .TRUE. )
CD       CALL DASSLF( M, R, LR, X( N + 1 ), .TRUE. )
      END IF
C
C  TRANSFORM BY R(INVERSE).
C
CS    CALL SASSLF( M, R, LR, X( N + 1 ), .FALSE. )
CD    CALL DASSLF( M, R, LR, X( N + 1 ), .FALSE. )
C
C  COMPUTE THE PART OF THE SOLUTION X1.
C
      IF ( N .LE. 0 ) GO TO 400
      DO 260 I       = 1, N
         VECTOR( I ) = ZERO
  260 CONTINUE
      DO 280 ICOL = 1, M
         SCALAR   = X( N + ICOL )
         K1       = IPBD( ICOL )
         K2       = IPBD( ICOL + 1 ) - 1
         IF ( K1 .LE. K2 ) THEN
            DO 270 K = K1, K2
               J     = IRNBD( K )
               IF ( J .LE. N )
     *             VECTOR( J ) = VECTOR( J ) + BD( K ) * SCALAR
  270       CONTINUE
         END IF
  280 CONTINUE
      INFORM = 3
C
C  RETURN TO GET THE PRODUCT A(INVERSE) * VECTOR.
C
      RETURN
  300 CONTINUE
CS    CALL SAXPY( N, - ONE, VECTOR, 1, X, 1 )
CD    CALL DAXPY( N, - ONE, VECTOR, 1, X, 1 )
  400 CONTINUE
      INFORM = 0
      RETURN
C
C  END OF ASSLB.
C
      END
C  THIS VERSION: 14/08/1991 AT 12:00:44 PM.
CS    SUBROUTINE SASSLC( N, M, ICLASS, LBD, BD, IRNBD, IPBD,
CD    SUBROUTINE DASSLC( N, M, ICLASS, LBD, BD, IRNBD, IPBD,
     *                   LCD, CD, JCNCD, IPCD, R, LR, Q, LQ,
     *                   SPIKE, VECTOR, INFORM )
C
C  FORM AND FACTORIZE THE SCHUR COMPLEMENT
C
C      S = D - C * A(INVERSE) * B
C
C  OF THE MATRIX A IN THE SYMMETRIC OR UNSYMMETRIC BLOCK MATRIX
C
C     / A  B \
C     \ C  D /
C
C  WHEN AN EXTRA ROW AND COLUMN ARE APPENDED.
C
C  NICK GOULD, AUGUST 3RD 1988.
C  FOR CGT PRODUCTIONS.
C
      INTEGER          INFORM, N, M, ICLASS, LBD, LCD, LR, LQ
      INTEGER          IRNBD( LBD ), JCNCD( LCD ), IPBD( * ), IPCD( * )
CS    REAL             BD( LBD ), CD( LCD ), R( LR ), Q( LQ, * ),
CD    DOUBLE PRECISION BD( LBD ), CD( LCD ), R( LR ), Q( LQ, * ),
     *                 VECTOR( * ), SPIKE( * )
C
C  LOCAL VARIABLES.
C
      INTEGER          I, IROW, J, K, K1, K2, NEWCLR, NEWDIR
      INTEGER          KIRN, NEWDIA, MNEW
CS    REAL             DIANEW, SCALAR, ZERO
CD    DOUBLE PRECISION DIANEW, SCALAR, ZERO
C
C  MACHINE FUNCTIONS.
C
      INTRINSIC        SQRT
CS    REAL             SDOT
CD    DOUBLE PRECISION DDOT
CS    EXTERNAL         SDOT, SASSLF
CD    EXTERNAL         DDOT, DASSLF
      COMMON / PRECNN / PRNTER
      LOGICAL PRNTER
      SAVE
C
C  DATA STATEMENTS.
C
CS    DATA ZERO / 0.0E+0 /
CD    DATA ZERO / 0.0D+0 /
      IF ( INFORM .LE. 0 ) THEN
         INFORM = - 5
         RETURN
      END IF
      GO TO ( 100, 200, 200 ), INFORM
  100 CONTINUE
C
C  THE NEW ROW AND COLUMN WILL BE THE MNEW-TH.
C
      MNEW   = M + 1
      NEWCLR = M * MNEW / 2
      NEWDIR = MNEW * ( MNEW + 1 ) / 2
C
C  ENSURE THAT THERE IS SUFFICIENT SPACE.
C
      IF ( M .LT. 0 .OR. N .LT. 0 .OR. ICLASS .LT. 1 .OR. ICLASS
     *     .GT. 4 .OR. LQ .LT. MNEW .OR. LR .LT. NEWDIR ) THEN
         INFORM = - 1
         RETURN
      END IF
C
C  REMOVE SUBDIAGONAL ELEMENTS FROM THE DATA STRUCTURES ASSOCIATED
C  WITH THE NEW ROW OF B AND THE UPPER TRIANGULAR PART OF D.
C
      K1       = IPBD( MNEW )
      K2       = IPBD( MNEW + 1 ) - 1
      KIRN     = K1
      NEWDIA   = N + MNEW
      DO 110 K = K1, K2
         I     = IRNBD( K )
         IF ( I .LE. NEWDIA ) THEN
            IRNBD( KIRN ) = I
            BD( KIRN )    = BD( K )
            KIRN          = KIRN + 1
         END IF
  110 CONTINUE
      IPBD( MNEW + 1 ) = KIRN
C
C  RETURN TO OBTAIN THE INVERSE OF A TIMES THE FIRST N COMPONENTS
C  OF COLNEW. FIRST CALCULATE COLNEW.
C
      IF ( M .GT. 0 ) THEN
         DO 130 I           = 1, M
            R( NEWCLR + I ) = ZERO
            SPIKE( I )      = ZERO
  130    CONTINUE
      END IF
      SPIKE( MNEW ) = ZERO
      IF ( N .GT. 0 ) THEN
         DO 140 I       = 1, N
            VECTOR( I ) = ZERO
  140    CONTINUE
      END IF
C
C  S IS SYMMETRIC.
C
      K1 = IPBD( MNEW )
      K2 = IPBD( MNEW + 1 ) - 1
      IF ( K1 .LE. K2 ) THEN
         DO 150 K = K1, K2
            J     = IRNBD( K )
            IF ( J .LE. N ) THEN
               VECTOR( J ) = BD( K )
            ELSE
               SPIKE( J - N ) = BD( K )
            END IF
  150    CONTINUE
      END IF
      DIANEW = SPIKE( MNEW )
      IF ( N .GT. 0 ) THEN
         INFORM = 2
         RETURN
      END IF
  200 CONTINUE
C
C  FORM THE NEW LAST COLUMN OF S.
C
      IF ( M .GT. 0 ) THEN
C
C  S IS NEGATIVE DEFINITE AND THE FACTORS OF - S ARE USED.
C
            DO 220 IROW = 1, M
               SCALAR   = - SPIKE( IROW )
               IF ( N .GT. 0 ) THEN
                  K1 = IPBD( IROW )
                  K2 = IPBD( IROW + 1 ) - 1
                  IF ( K1 .LE. K2 ) THEN
                     DO 210 K = K1, K2
                        J     = IRNBD( K )
                        IF ( J .LE. N )
     *                     SCALAR = SCALAR + BD( K ) * VECTOR( J )
  210                CONTINUE
                  END IF
               END IF
C
C  TRANSFORM THE NEW COLUMN OF R.
C
               R( NEWCLR + IROW ) = SCALAR
  220       CONTINUE
      END IF
C
C  FORM THE NEW DIAGONAL OF S PRIOR TO THE R(TRANS) R FACTORIZATION.
C
      IF ( N .GT. 0 ) THEN
         K1 = IPBD( MNEW )
         K2 = IPBD( MNEW + 1 ) - 1
         IF ( K1 .LE. K2 ) THEN
            DO 520 K = K1, K2
               J     = IRNBD( K )
               IF ( J .LE. N ) DIANEW = DIANEW - BD( K ) * VECTOR( J )
  520       CONTINUE
         END IF
      END IF
            SCALAR = - DIANEW
C
C  FIND THE NEW COLUMN OF R.
C
         IF ( M .GT. 0 ) THEN
            NEWCLR = NEWCLR + 1
CS          CALL SASSLF( M, R, LR, R( NEWCLR ), .TRUE. )
CD          CALL DASSLF( M, R, LR, R( NEWCLR ), .TRUE. )
CS          SCALAR = SCALAR - SDOT( M, R( NEWCLR ), 1, R( NEWCLR ), 1 )
CD          SCALAR = SCALAR - DDOT( M, R( NEWCLR ), 1, R( NEWCLR ), 1 )
         END IF
C
C  CHECK THAT THE MATRIX IS INDEED POSITIVE DEFINITE.
C
         IF ( SCALAR .LE. ZERO ) THEN
            INFORM = - ICLASS
            RETURN
         END IF
C
C  FIND THE NEW DIAGONAL OF R.
C
         R( NEWDIR ) = SQRT( SCALAR )
      INFORM = 0
      RETURN
C
C  END OF ASSLC.
C
      END
C  THIS VERSION: 14/08/1991 AT 12:00:44 PM.
CS    SUBROUTINE SASSLE( MSOFAR, JBEGIN, R, LR, Q, LQ, SPIKE,
CD    SUBROUTINE DASSLE( MSOFAR, JBEGIN, R, LR, Q, LQ, SPIKE,
     *                   QUSED, INFORM )
C
C  USE PLANE-ROTATION MATRICES TO REDUCE AN UPPER TRIANGULAR
C  MATRIX WITH A NEW HORIZONTAL SPIKE ROW TO UPPER TRIANGULAR FORM.
C  THE SPIKE IS INPUT IN THE ARRAY SPIKE AND ITS
C  FIRST NONZERO OCCURS IN POSITION IBEGIN.
C  IF THE LOGICAL QUSED IS TRUE, APPLY THE PLANE-ROTATIONS TO Q.
C
C  NICK GOULD, AUGUST 5TH 1988.
C  FOR CGT PRODUCTIONS.
C
      INTEGER          MSOFAR, JBEGIN, LR, LQ, INFORM
      LOGICAL          QUSED
CS    REAL             R( LR ), Q( LQ, * ), SPIKE( * )
CD    DOUBLE PRECISION R( LR ), Q( LQ, * ), SPIKE( * )
C
C  LOCAL VARIABLES.
C
      INTEGER          J, K, MNEW, NEXTR, NEXTW
CS    REAL             C, S, X, Y, ZERO
CD    DOUBLE PRECISION C, S, X, Y, ZERO
C
C  EXTERNAL FUNCTIONS.
C
CS    EXTERNAL         SROT, SROTG
CD    EXTERNAL         DROT, DROTG
C
C  SET DATA.
C
CS    DATA ZERO / 0.0E+0 /
CD    DATA ZERO / 0.0D+0 /
      MNEW = MSOFAR + 1
C
C  REDUCE THE NEW ROW TO ZERO BY APPLYING PLANE-ROTATION MATRICES.
C
      IF ( JBEGIN .LE. MSOFAR ) THEN
         DO 20 J  = JBEGIN, MSOFAR
            NEXTW = J + 1
            NEXTR = J * NEXTW / 2
C
C  USE A PLANE-ROTATION IN THE PLANE ( J, MNEW ).
C
CS          CALL SROTG( R( NEXTR ), SPIKE( J ), C, S )
CD          CALL DROTG( R( NEXTR ), SPIKE( J ), C, S )
C
C  APPLY THE PLANE-ROTATIONS TO THE REMAINING ELEMENTS IN ROWS
C  J AND MNEW OF R.
C
            NEXTR         = NEXTR + J
            DO 10 K       = J + 1, MNEW
               X          = R( NEXTR )
               Y          = SPIKE( K )
               R( NEXTR ) = C * X + S * Y
               SPIKE( K ) = C * Y - S * X
               NEXTR      = NEXTR + K
   10       CONTINUE
C
C  APPLY THE PLANE-ROTATIONS TO THE REMAINING ELEMENTS IN COLUMNS
C  J AND MNEW OF Q.
C
CS          IF ( QUSED ) CALL SROT( MNEW, Q( 1, J ), 1,
CD          IF ( QUSED ) CALL DROT( MNEW, Q( 1, J ), 1,
     *                              Q( 1, MNEW ), 1, C, S )
   20    CONTINUE
      END IF
C
C  CHECK THAT THE NEW DIAGONAL ENTRY OF R IS NON-ZERO.
C
      IF ( SPIKE( MNEW ) .EQ. ZERO ) THEN
         INFORM = - 2
      ELSE
         R( MNEW * ( MNEW + 1 ) / 2  ) = SPIKE( MNEW )
      END IF
      RETURN
      END
C  THIS VERSION: 14/08/1991 AT 12:00:44 PM.
CS    SUBROUTINE SASSLF( NR, R, LR, X, TRANS )
CD    SUBROUTINE DASSLF( NR, R, LR, X, TRANS )
C
C  COMPUTE THE SOLUTION TO THE TRIANGULAR SYSTEMS
C
C  R * X(OUTPUT) = X(INPUT) (TRANS = .FALSE.) OR
C
C  R(TRANSPOSE) * X(OUTPUT) = X(INPUT) (TRANS = .TRUE.),
C
C  WHERE R IS AN UPPER TRIANGULAR MATRIX STORED BY COLUMNS.
C
C  NICK GOULD, AUGUST 15TH 1988.
C  FOR CGT PRODUCTIONS.
C
      INTEGER          NR, LR
      LOGICAL          TRANS
CS    REAL             R( LR ), X( * )
CD    DOUBLE PRECISION R( LR ), X( * )
C
C  LOCAL VARIABLES.
C
      INTEGER          I, II, NEXTR
CS    REAL             SCALAR
CD    DOUBLE PRECISION SCALAR
C
C  EXTERNAL FUNCTIONS.
C
CS    EXTERNAL         SAXPY, SDOT
CD    EXTERNAL         DAXPY, DDOT
CS    REAL             SDOT
CD    DOUBLE PRECISION DDOT
      IF ( NR .LE. 0 ) RETURN
      IF ( TRANS ) THEN
C
C  SOLVE R(TRANSPOSE) * X(OUTPUT) = X(INPUT).
C
         X( 1 ) = X( 1 ) / R( 1 )
         IF ( NR .GT. 1 ) THEN
            NEXTR     = 2
            DO 100 I  = 2, NR
CS             SCALAR = X( I ) - SDOT( I - 1, R( NEXTR ), 1, X, 1 )
CD             SCALAR = X( I ) - DDOT( I - 1, R( NEXTR ), 1, X, 1 )
               NEXTR  = NEXTR + I
               X( I ) = SCALAR / R( NEXTR - 1 )
  100       CONTINUE
         END IF
      ELSE
C
C  SOLVE R * X(OUTPUT) = X(INPUT).
C
         NEXTR     = NR * ( NR + 1 ) / 2
         DO 200 II = 1, NR
            I      = NR + 1 - II
            SCALAR = X( I ) / R( NEXTR )
            X( I ) = SCALAR
            NEXTR  = NEXTR - I
CS          CALL SAXPY( I - 1, - SCALAR, R( NEXTR + 1 ), 1, X, 1 )
CD          CALL DAXPY( I - 1, - SCALAR, R( NEXTR + 1 ), 1, X, 1 )
  200    CONTINUE
      END IF
      RETURN
C
C  END OF ASSLF.
C
      END
