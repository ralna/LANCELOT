C     ******************************************************************
C
CS    SUBROUTINE SDSCDR( N, IW, A, RHS, NULROW, CONSIS,
CD    SUBROUTINE DDSCDR( N, IW, A, RHS, NULROW, CONSIS,
     *                   IW1, LIW1, WORK, EPSMCH )
C
C     ******************************************************************
C
      INTEGER N, LIW1
      LOGICAL CONSIS
      INTEGER IW( * ), NULROW( N ), IW1( LIW1 )
CS    REAL             A( * ), RHS( N ), WORK( N ), EPSMCH
CD    DOUBLE PRECISION A( * ), RHS( N ), WORK( N ), EPSMCH
C
C     ******************************************************************
C
C     PROGRAMMING :    M. LESCRENIER ( AUG 1987 ).
C     =============
C
C     DESCRIPTION :
C     =============
C
C     GIVEN THE GAUSSIAN FACTORIZED FORM OF A SYMMETRIC MATRIX
C     KNOWN TO BE SINGULAR
C                                        T        T
C                                   P A P  = L D L
C
C     WHERE D IS A BLOCK DIAGONAL MATRIX WITH 1X1 OR 2X2 BLOCKS
C     COMPUTE THE SOLUTION S OF THE SYSTEM
C                                           A * S = RHS     (1)
C     IF THIS SYSTEM IS CONSISTENT,
C     OTHERWISE,
C
C     COMPUTE A DESCENT DIRECTION S FOR THE QUADRATIC Q
C                           T
C             Q(X) = 0.5 * X * A * X - RHS * X
C                                                 T
C     I.E. A DIRECTION S SUCH THAT A * S = 0 AND S * RHS > 0
C
C     METHOD:
C     =======
C     SEE "ON THE LOCATION OF DIRECTIONS OF INFINITE DESCENT FOR
C          NONLINEAR PROGRAMMING ALGORITHMS" BY
C         A.R. CONN AND N.I.M. GOULD,
C         SIAM J. NUM. ANAL. VOL 21, NO 6, DEC 1984.
C
C     MORE SPECIFICALLY, LEMMA 2.1, P 1167.
C
C     PARAMETERS :
C     ============
C
C     N
C            I: ORDER OF THE MATRIX A.
C            O: UNMODIFIED.
C
C     IW
C            I: ABS(IW(1)) = NUMBER OF BLOCK PIVOT ROWS.
C               IW(1) < 0 MEANS THAT AT LEAST ONE 2X2 PIVOT HAS BEEN
C                         USED.
C               INTEGER INFORMATION ON EACH BLOCK PIVOT ROW FOLLOWS.
C               FOR EACH BLOCK PIVOT ROW, WE STORE
C                   1) THE NUMBER OF COLUMNS
C                   2) THE NUMBER OF ROWS
C                   3) THE COLUMN INDICES
C               IF THE BLOCK PIVOT ROW IS REDUCED TO 1 ROW, WE STORE
C                   1) -(THE NUMBER OF COLUMNS)
C                   3) THE COLUMN INDICES
C               THE FIRST COLUMN INDEX OF A 2X2 PIVOT IS FLAGGED
C               NEGATIVE.
C            O: UNMODIFIED.
C
C     A
C            I: THE FACTORS STORED ROW BY ROW,
C               THE INVERSES OF THE PIVOTS ARE STORED INSTEAD OF THE
C                   PIVOTS THEMSELVES,
C               THE OTHER ENTRIES ARE NEGATED.
C            O: UNMODIFIED.
C
C     RHS
C            I: SEE RHS IN DESCRIPTION ABOVE.
C            O: SEE S   IN DESCRIPTION ABOVE.
C
C     NULROW
C            I: IF THE I-TH ROW OF THE MATRIX IS ENTIRELY NUL
C               NULROW(I)=1;
C               NULROW(I)=0 OTHERWISE.
C            O: UNMODIFIED.
C
C     CONSIS
C            I: UNDEFINED.
C            O: .TRUE. IF THE SYSTEM (1) IS CONSISTENT
C               .FALSE. OTHERWISE.
C
C     IW1
C            I: UNDEFINED WORK VECTOR OF LENGTH LIW1.
C            O: IW1(I) = POINTER TO THE BEGINNING OF EACH BLOCK PIVOT
C                        ROW IN IW (I=1,...,ABS(IW(1)))
C     LIW1
C            I: LENGTH OF THE VECTOR IW1 (N IS SAFE)
C            O: UNMODIFIED
C
C     WORK
C            I: UNDEFINED WORK VECTOR OF LENGTH N.
C            O: MEANINGLESS.
C
C     EPSMCH
C            I: RELATIVE PRECISION OF THE MACHINE.
C            O: UNMODIFIED.
C
C     ******************************************************************
C
      INTEGER I, IVAR1, IVAR2, LATOP, NBLK, NCOLS, NROWS, IS, IAPOS2
      INTEGER IAPOS, IPOS, J, J1, J2, IBLK, IROW, K, KK1, KK2, I1S, I2S 
      INTEGER IWORK, LG, K1, K2, K3, IVAR, LOOP, IIS, JPOS, JJ1, JJ2
CS    REAL             V1, V2, V3, W1, W2, TOLER, FACTOR, ZERO
CD    DOUBLE PRECISION V1, V2, V3, W1, W2, TOLER, FACTOR, ZERO
      INTRINSIC ABS
CS    DATA FACTOR, ZERO / 1.0E+2, 0.0E+0 /
CD    DATA FACTOR, ZERO / 1.0D+4, 0.0D+0 /
      CONSIS = .TRUE.
      TOLER  = FACTOR * EPSMCH
C
C     FIND THE VECTOR WORK SUCH THAT L * WORK = P * RHS
C     CHECK CONSISTENCY IN THE CASE OF ENTIRELY NUL ROWS.
C
      DO 5 I = 1, N
         WORK( I ) = RHS( I )
         IF ( NULROW( I ) .EQ. 1 ) THEN
            IF ( ABS( RHS( I ) ) .GT. TOLER ) CONSIS = .FALSE.
         END IF
5     CONTINUE
      NBLK = ABS( IW( 1 ) + 0 )
      IF ( NBLK .GT. LIW1 ) THEN
         WRITE( 6, 1000 ) NBLK
         STOP
      END IF
C
C IAPOS. RUNNING POINTER TO CURRENT PIVOT POSITION IN ARRAY A.
C IPOS.  RUNNING POINTER TO BEGINNING OF BLOCK PIVOT ROW IN IW.
C
      IAPOS = 1
      IPOS  = 2
      J2    = 1
      IBLK  = 0
      NROWS = 0
      DO 70 IROW = 1, N
        IF ( NROWS .GT. 0 ) GO TO 20
        IBLK = IBLK + 1
        IF ( IBLK .GT. NBLK ) GO TO 80
        IPOS = J2 + 1
C ABS(NCOLS) IS NUMBER OF VARIABLES (COLUMNS) IN BLOCK PIVOT ROW.
        NCOLS = - IW( IPOS )
C NROWS IS NUMBER OF ROWS IN BLOCK PIVOT.
        NROWS = 1
        IF ( NCOLS .GT. 0 ) GO TO 10
        NCOLS = - NCOLS
        IPOS = IPOS + 1
        NROWS = IW( IPOS )
   10   CONTINUE
        J1 = IPOS + 1
        J2 = IPOS + NCOLS
        IF ( IBLK .EQ. NBLK ) THEN
           KK1 = IPOS + NROWS + 1
           KK2 = IPOS + NCOLS
        END IF
C JUMP IF WE HAVE A 2 BY 2 PIVOT.
   20   CONTINUE
        IF ( IW( J1 ) .LT. 0 ) GO TO 40
C PERFORM FORWARD SUBSTITUTION USING 1 BY 1 PIVOT.
        NROWS = NROWS - 1
        IAPOS = IAPOS + 1
        J1 = J1 + 1
        IF (J1.GT.J2) GO TO 70
        IWORK = IW(J1-1)
        W1 = WORK(IWORK)
C
C       CHECK CONSISTENCY OF THE LINEAR SYSTEM (1)
C
        IF (A(IAPOS-1).EQ.ZERO.AND.ABS(W1).GT.TOLER) CONSIS = .FALSE.
        IF (ABS(W1).GT.EPSMCH) THEN
           K = IAPOS
           DO 30 J=J1,J2
             IWORK = ABS(IW(J)+0)
             WORK(IWORK) = WORK(IWORK) + A(K)*W1
             K = K + 1
  30       CONTINUE
        END IF
        IAPOS = IAPOS + J2 - J1 + 1
        GO TO 70
C PERFORM OPERATIONS WITH 2 BY 2 PIVOT
  40    NROWS = NROWS - 2
        J1 = J1 + 2
        IAPOS = IAPOS + 2
        IF (J1.GT.J2) GO TO 60
        IWORK = -IW(J1-2)
        W1 = WORK(IWORK)
        IWORK = IW(J1-1)
        W2 = WORK(IWORK)
        K1 = IAPOS
        K3 = IAPOS + J2 - J1 + 2
        DO 50 J=J1,J2
          IWORK = ABS(IW(J)+0)
          WORK(IWORK) = WORK(IWORK) + W1*A(K1) + W2*A(K3)
          K1 = K1 + 1
          K3 = K3 + 1
  50    CONTINUE
  60    IAPOS = IAPOS + 2*(J2-J1+1) + 1
  70  CONTINUE
  80  LATOP = IAPOS - 1
C
C     CHECK CONSISTENCY OF THE LINEAR SYSTEM (1)
C
      DO 85 K=KK1,KK2
         IF (ABS(WORK(IW(K))).GT.TOLER) CONSIS = .FALSE.
85    CONTINUE
      IF ( CONSIS ) THEN
C
C       REDEFINE RHS AS THE SOLUTION X OF D * X = WORK
C
        IPOS = 2
        IAPOS = 1
        DO 100 IBLK=1,NBLK
C
C          CONSIDER THE IBLK-TH BLOCK PIVOT.
C
           IW1(IBLK) = IPOS
           IF (IW(IPOS).GT.0) THEN
              NCOLS = IW(IPOS)
              NROWS = IW(IPOS+1)
              J1 = IPOS + 2
           ELSE
              NCOLS = -IW(IPOS)
              NROWS = 1
              J1 = IPOS + 1
           END IF
           J2 = J1 + NROWS - 1
           IPOS = J1 + NCOLS
           LG = NCOLS
           J = J1
C
C          CONSIDER THE NEXT PIVOT (1X1 OR 2X2) OF THIS BLOCK PIVOT.
C
90         IF (IW(J).GT.0) THEN
C
C             IT IS A 1X1 PIVOT.
C
              RHS(IW(J)) = A(IAPOS)*WORK(IW(J))
              IAPOS = IAPOS + LG
              LG = LG - 1
              J = J + 1
           ELSE
C
C             IT IS A 2X2 PIVOT.
C
              V1 = A(IAPOS)
              V2 = A(IAPOS+1)
              IAPOS = IAPOS + LG
              LG = LG - 1
              V3 = A(IAPOS)
              IAPOS = IAPOS + LG
              LG = LG - 1
              IVAR1 = -IW(J)
              IVAR2 = IW(J+1)
              RHS(IVAR1) = V1*WORK(IVAR1) + V2*WORK(IVAR2)
              RHS(IVAR2) = V2*WORK(IVAR1) + V3*WORK(IVAR2)
              J = J + 2
           END IF
           IF (J.LE.J2) GO TO 90
100     CONTINUE
        LATOP = IAPOS - 1
      ELSE
C
C       REDEFINE RHS AS A VECTOR X <> 0 SUCH THAT D * X = 0
C
        IPOS = 2
        IAPOS = 1
        DO 120 IBLK=1,NBLK
C
C          CONSIDER THE IBLK-TH BLOCK PIVOT.
C
           IW1(IBLK) = IPOS
           IF (IW(IPOS).GT.0) THEN
              NCOLS = IW(IPOS)
              NROWS = IW(IPOS+1)
              J1 = IPOS + 2
           ELSE
              NCOLS = -IW(IPOS)
              NROWS = 1
              J1 = IPOS + 1
           END IF
           J2 = J1 + NROWS - 1
           IF (IBLK.EQ.NBLK) THEN
              K1 = J2 + 1
              K2 = J1 + NCOLS - 1
              DO 105 K=K1,K2
                 IWORK = IW(K)
                 RHS(IWORK) = WORK(IWORK)
105           CONTINUE
           END IF
           IPOS = J1 + NCOLS
           LG = NCOLS
           J = J1
C
C          CONSIDER THE NEXT PIVOT (1X1 OR 2X2) OF THIS BLOCK PIVOT.
C
110        IF (IW(J).GT.0) THEN
C
C             IT IS A 1X1 PIVOT.
C
              IF ( A( IAPOS ) .EQ. ZERO ) THEN
                 RHS( IW( J ) ) = WORK( IW( J ) )
              ELSE
                 RHS(IW(J)) = ZERO
              END IF
              IAPOS = IAPOS + LG
              LG = LG - 1
              J = J + 1
           ELSE
C
C             IT IS A 2X2 PIVOT.
C
              IAPOS = IAPOS + LG
              LG = LG - 1
              IAPOS = IAPOS + LG
              LG = LG - 1
              IVAR = -IW(J)
              RHS(IVAR) = ZERO
              IVAR = IW(J+1)
              RHS(IVAR) = ZERO
              J = J + 2
           END IF
           IF (J.LE.J2) GO TO 110
120     CONTINUE
        LATOP = IAPOS - 1
      END IF
C
C     LATOP IS A POINTER TO THE LAST ENTRY OF THE FACTORS OF A.
C
C     REDEFINE RHS AS THE SOLUTION X OF THE UPPER TRIANGULAR SYSTEM
C                   T
C                 (L *P) * X = RHS
C
C     IAPOS : RUNNING POINTER TO CURRENT PIVOT POSITION IN ARRAY A.
C     IPOS  : RUNNING POINTER TO BEGINNING OF CURRENT BLOCK PIVOT ROW.
C
      IAPOS = LATOP + 1
      NROWS = 0
      IBLK = NBLK + 1
C
C     RUN THROUGH BLOCK PIVOT ROWS IN THE REVERSE ORDER.
C
      DO 190 LOOP=1,N
        IF (NROWS.GT.0) GO TO 140
C
C       CONSIDER NEXT BLOCK PIVOT ROW.
C
        IBLK = IBLK - 1
C
C       STOPPING CRITERIA.
C
        IF (IBLK.LT.1) GO TO 200
        IPOS = IW1(IBLK)
C
C       ABS(NCOLS) IS NUMBER OF VARIABLES (COLUMNS) IN BLOCK PIVOT
C       ROW.
C       NROWS IS NUMBER OF ROWS IN BLOCK PIVOT ROW.
C
        NCOLS = -IW(IPOS)
        NROWS = 1
        IF (NCOLS.GT.0) GO TO 130
        NCOLS = -NCOLS
        IPOS = IPOS + 1
        NROWS = IW(IPOS)
  130   JPOS = IPOS + NROWS
        J2 = IPOS + NCOLS
C
C       CONSIDER NEXT PIVOT IN THE IBLK-TH BLOCK PIVOT ROW.
C
  140   IF (NROWS.EQ.1) GO TO 150
C
C       JUMP IF WE HAVE A 2 BY 2 PIVOT.
C
        IF (IW(JPOS-1).LT.0) GO TO 170
C
C       PERFORM BACK-SUBSTITUTION USING 1 BY 1 PIVOT.
C
  150   NROWS = NROWS - 1
        IAPOS = IAPOS - (J2-JPOS+1)
        IIS = IW(JPOS)
        W1 = RHS(IIS)
        J1 = JPOS + 1
        K = IAPOS + 1
        DO 160 J=J1,J2
          IS = ABS(IW(J)+0)
          W1 = W1 + A(K)*RHS(IS)
          K = K + 1
  160   CONTINUE
        RHS(IIS) = W1
        JPOS = JPOS - 1
        GO TO 190
C
C       PERFORM OPERATIONS WITH 2 BY 2 PIVOT
C
  170   NROWS = NROWS - 2
        IAPOS2 = IAPOS - (J2-JPOS+1)
        IAPOS = IAPOS2 - (J2-JPOS+2)
        I1S = -IW(JPOS-1)
        I2S = IW(JPOS)
        W1 = RHS(I1S) + RHS(I2S)
        W2 = RHS(I1S) + RHS(I2S)
        J1 = JPOS + 1
        JJ1 = IAPOS + 2
        JJ2 = IAPOS2 + 1
        DO 180 J=J1,J2
          IS = ABS(IW(J)+0)
          W1 = W1 + RHS(IS)*A(JJ1)
          W2 = W2 + RHS(IS)*A(JJ2)
          JJ1 = JJ1 + 1
          JJ2 = JJ2 + 1
  180   CONTINUE
        RHS(I1S) = W1
        RHS(I2S) = W2
        JPOS = JPOS - 2
  190 CONTINUE
  200 CONTINUE
1000  FORMAT(' Message from -DSCDIR-',/,
     *       ' Increase the parameter -LIW1- to ',I5)
      RETURN
      END
C     ******************************************************************
C
CS    SUBROUTINE SNEGCR( N, IW, A, IW1, LIW1, NUMBER, NEXT, S, SAS )
CD    SUBROUTINE DNEGCR( N, IW, A, IW1, LIW1, NUMBER, NEXT, S, SAS )
C
C     ******************************************************************
C
      LOGICAL NEXT
      INTEGER N, LIW1, IW1( LIW1 ), NUMBER
      INTEGER IW( * )
CS    REAL             A( * ), S( * ), SAS
CD    DOUBLE PRECISION A( * ), S( * ), SAS
C
C     ******************************************************************
C
C     PROGRAMMING :    M. LESCRENIER ( AUG 1987 ).
C     =============
C
C     DESCRIPTION :
C     =============
C
C     GIVEN THE GAUSSIAN FACTORIZED FORM OF A SYMMETRIC MATRIX
C     KNOWN TO BE NON POSITIVE DEFINITE,
C                                        T        T
C                                   P A P  = L D L
C
C     WHERE D IS A BLOCK DIAGONAL MATRIX WITH 1X1 OR 2X2 BLOCKS
C     COMPUTE A DIRECTION OF NEGATIVE CURVATURE OF A.
C
C     METHOD:
C     -------
C     1)CHOOSE A NEGATIVE EIGENVALUE OF D, L, AND
C       ITS CORRESPONDING EIGENVECTOR V
C     2)S IS CHOSEN AS THE SOLUTION OF THE UPPER TRIANGULAR SYSTEM
C
C                     T
C                 (D*L *P)*S = D*V
C
C     WE KNOW THEN THAT
C                   T       T  T     T         T           2
C                  S A S = S  P L D L  P S  = V D V = L ]V]  < 0
C
C     PARAMETERS :
C     ============
C
C     N
C            I: ORDER OF THE MATRIX A.
C            O: UNMODIFIED.
C
C     IW
C            I: ABS(IW(1)) = NUMBER OF BLOCK PIVOT ROWS.
C               IW(1) < 0 MEANS THAT AT LEAST ONE 2X2 PIVOT HAS BEEN
C                         USED.
C               INTEGER INFORMATION ON EACH BLOCK PIVOT ROW FOLLOWS.
C               FOR EACH BLOCK PIVOT ROW, WE STORE
C                   1) THE NUMBER OF COLUMNS
C                   2) THE NUMBER OF ROWS
C                   3) THE COLUMN INDICES
C               IF THE BLOCK PIVOT ROW IS REDUCED TO 1 ROW, WE STORE
C                   1) -(THE NUMBER OF COLUMNS)
C                   3) THE COLUMN INDICES
C               THE FIRST COLUMN INDEX OF A 2X2 PIVOT IS FLAGGED
C               NEGATIVE.
C            O: UNMODIFIED.
C
C     A
C            I: THE FACTORS STORED ROW BY ROW,
C               THE INVERSES OF THE PIVOTS ARE STORED INSTEAD OF THE
C                   PIVOTS THEMSELVES,
C               THE OTHER ENTRIES ARE NEGATED.
C            O: UNMODIFIED.
C
C     IW1
C            I: UNDEFINED WORK VECTOR OF LENGTH LIW1.
C            O: IW1(I) = POINTER TO THE BEGINNING OF EACH BLOCK PIVOT
C                        ROW IN IW (I=1,...,ABS(IW(1)))
C     LIW1
C            I: LENGTH OF THE VECTOR IW1 (N IS SAFE)
C            O: UNMODIFIED
C
C     NUMBER
C            I: NUMBER OF THE PREVIOUS CHOSEN NEGATIVE EIGENVALUE.
C            I: NUMBER OF THE NEW CHOSEN NEGATIVE EIGENVALUE.
C
C     NEXT
C            I: IF NEXT=.TRUE. CHOOSE THE NEXT NEGATIVE EIGENVALUE
C               OTHERWISE, CHOOSE THE SMALLEST NEGATIVE EIGENVALUE.
C            O: UNMODIFIED
C
C     S
C            I: UNDEFINED.
C            O: DIRECTION OF NEGATIVE CURVATURE OF A.
C
C     SAS
C            I: UNDEFINED.
C                                       T
C            O: THE NEGATIVE CURVATURE S A S
C
C     ******************************************************************
C
      LOGICAL TWO, FIRST
      INTEGER IAPOS, IPOS, J, J1, J2, IBLK,K, INUM
      INTEGER IS, I1S, I2S, LG, LOOP, IIS, JPOS, JJ1, JJ2, IAPOS2
      INTEGER IVAR1, IVAR2, IVAR3, IVAR4, LATOP, NBLK, NCOLS, NROWS
      INTEGER NEWNUM, I
CS    REAL             LMIN, L, DET, RHO, AUX
CS    REAL             W1, W2, V1, V2, V3, V4, V5, V6, ZERO
CD    DOUBLE PRECISION LMIN, L, DET, RHO, AUX
CD    DOUBLE PRECISION W1, W2, V1, V2, V3, V4, V5, V6, ZERO
      INTRINSIC ABS, SQRT
CS    DATA ZERO / 0.0E+0 /
CD    DATA ZERO / 0.0D+0 /
C
C     INUM COUNTS THE ENCOUNTERED NEGATIVE EIGENVALUES
C
      INUM = 0
      IF ( .NOT. NEXT ) THEN
C
C        CHOOSE THE SMALLEST EIGENVALUE OF D
C        -----------------------------------
C        TWO IS SET TO FALSE IF THE SMALLEST EIGENVALUE OF D IS A
C        1X1 BLOCK, TWO IS SET TO TRUE OTHERWISE.
C
C        V1, ...,V6 ARE USED TO STORE TWO 2X2 BLOCKS OF D.
C
C        IF TWO=FALSE,
C            IVAR1 = THE VARIABLE INDEX CORRESPONDING TO THE
C                    SMALLEST EIGENVALUE OF D
C        OTHERWISE,
C            IVAR1 AND IVAR2 ARE THE VARIABLES INDICES
C            OF THE 2X2 BLOCK OF D HOLDING THE SMALLEST EIGENVALUE.
C            IVAR3 AND IVAR4 MAY BE THE VARIABLES INDICES CORRESPONDING
C            TO THE SECOND 2X2 BLOCK HELD IN V4,V5,V6.
C         LATOP IS A POINTER TO THE LAST ENTRY OF THE FACTORS OF A.
C
         NBLK = ABS(IW(1)+0)
         IF (NBLK.GT.LIW1) THEN
            WRITE(6,1000) NBLK
            STOP
         END IF
         LMIN = ZERO
         IPOS = 2
         IAPOS = 1
         DO 20 IBLK=1,NBLK
C
C           CONSIDER THE IBLK-TH BLOCK PIVOT.
C
            IW1(IBLK) = IPOS
            IF (IW(IPOS).GT.0) THEN
               NCOLS = IW(IPOS)
               NROWS = IW(IPOS+1)
               J1 = IPOS + 2
            ELSE
               NCOLS = -IW(IPOS)
               NROWS = 1
               J1 = IPOS + 1
            END IF
            J2 = J1 + NROWS - 1
            IPOS = J1 + NCOLS
            LG = NCOLS
            J = J1
C
C           CONSIDER THE NEXT PIVOT (1X1 OR 2X2) OF THIS BLOCK PIVOT.
C
10          IF (IW(J).GT.0) THEN
C
C              IT IS A 1X1 PIVOT.
C
               L = 1.0D0/A(IAPOS)
               IAPOS = IAPOS + LG
               LG = LG - 1
               IF (L.LT.ZERO) INUM = INUM + 1
               IF (L.LT.LMIN) THEN
                  LMIN = L
                  IVAR1 = IW(J)
                  TWO = .FALSE.
                  NEWNUM = INUM
               END IF
               J = J + 1
            ELSE
C
C              IT IS A 2X2 PIVOT.
C
               V4 = A(IAPOS)
               V5 = A(IAPOS+1)
               IAPOS = IAPOS + LG
               LG = LG - 1
               V6 = A(IAPOS)
               IAPOS = IAPOS + LG
               LG = LG - 1
               IVAR3 = -IW(J)
               IVAR4 = IW(J+1)
C
C              COMPUTE THE 2X2 PIVOT.
C
               DET  = 1.0D0/(V4*V6 - V5**2)
               AUX  = V6
               V6 = V4*DET
               V4 = AUX*DET
               V5 = -V5*DET
C
C              COMPUTE THE SMALLEST EIGENVALUE OF THIS 2X2 PIVOT.
C
               RHO = (V4+V6)**2 - 4.0D0*(V4*V6-V5**2)
               L = 0.5*(V4+V6-SQRT(RHO))
               IF (L.LT.ZERO) INUM = INUM + 1
               IF (L.LT.LMIN) THEN
                  LMIN = L
                  IVAR1 = IVAR3
                  IVAR2 = IVAR4
                  V1 = V4
                  V2 = V5
                  V3 = V6
                  TWO = .TRUE.
                  NEWNUM = INUM
               END IF
               J = J + 2
            END IF
            IF (J.LE.J2) GO TO 10
20       CONTINUE
         LATOP = IAPOS - 1
         NUMBER = NEWNUM
      ELSE
C
C        CHOOSE THE NEXT NEGATIVE EIGENVALUE
C        -----------------------------------
C        LATOP IS A POINTER TO THE LAST ENTRY OF THE FACTORS OF A.
C
         FIRST = .TRUE.
         NBLK = ABS(IW(1)+0)
         IF (NBLK.GT.LIW1) THEN
            WRITE(6,1000) NBLK
            STOP
         END IF
         IPOS = 2
         IAPOS = 1
         DO 21 IBLK=1,NBLK
C
C           CONSIDER THE IBLK-TH BLOCK PIVOT.
C
            IW1(IBLK) = IPOS
            IF (IW(IPOS).GT.0) THEN
               NCOLS = IW(IPOS)
               NROWS = IW(IPOS+1)
               J1 = IPOS + 2
            ELSE
               NCOLS = -IW(IPOS)
               NROWS = 1
               J1 = IPOS + 1
            END IF
            J2 = J1 + NROWS - 1
            IPOS = J1 + NCOLS
            LG = NCOLS
            J = J1
C
C           CONSIDER THE NEXT PIVOT (1X1 OR 2X2) OF THIS BLOCK PIVOT.
C
11          IF (IW(J).GT.0) THEN
C
C              IT IS A 1X1 PIVOT.
C
               L = 1.0D0/A(IAPOS)
               IAPOS = IAPOS + LG
               LG = LG - 1
               IF (L.LT.ZERO) THEN
                  INUM = INUM + 1
                  IF (FIRST) THEN
                     LMIN = L
                     IVAR1 = IW(J)
                     TWO = .FALSE.
                     FIRST = .FALSE.
                     NEWNUM = 1
                  END IF
                  IF (INUM.EQ.NUMBER+1) THEN
                     NEWNUM = NUMBER + 1
                     LMIN = L
                     IVAR1 = IW(J)
                     TWO = .FALSE.
                  END IF
               END IF
               J = J + 1
            ELSE
C
C              IT IS A 2X2 PIVOT.
C
               V1 = A(IAPOS)
               V2 = A(IAPOS+1)
               IAPOS = IAPOS + LG
               LG = LG - 1
               V3 = A(IAPOS)
               IAPOS = IAPOS + LG
               LG = LG - 1
               IVAR1 = -IW(J)
               IVAR2 = IW(J+1)
C
C              COMPUTE THE 2X2 PIVOT.
C
               DET  = 1.0D0/(V1*V3 - V2**2)
               AUX  = V3
               V3 = V1*DET
               V1 = AUX*DET
               V2 = -V2*DET
C
C              COMPUTE THE SMALLEST EIGENVALUE OF THIS 2X2 PIVOT.
C
               RHO = (V1+V3)**2 - 4.0D0*(V1*V3-V2**2)
               L = 0.5*(V1+V3-SQRT(RHO))
               IF (L.LT.ZERO) THEN
                  INUM = INUM + 1
                  IF (FIRST) THEN
                     NEWNUM = 1
                     LMIN = L
                     TWO = .TRUE.
                     FIRST = .FALSE.
                  END IF
                  IF (INUM.EQ.NUMBER+1) THEN
                     NEWNUM = NUMBER + 1
                     LMIN = L
                     TWO = .TRUE.
                  END IF
               END IF
               J = J + 2
            END IF
            IF (J.LE.J2) GO TO 11
21       CONTINUE
         LATOP = IAPOS - 1
         NUMBER = NEWNUM
      END IF
      DO 30 I = 1, N
         S(I) = ZERO
30    CONTINUE
      IF (TWO) THEN
C
C        COMPUTE THE PRODUCT D*V WHERE D IS THE BLOCK DIAGONAL MATRIX
C        IS V IS THE EIGENVECTOR CORRESPONDING TO THE SMALLEST
C        EIGENVALUE OF D. NOTE THAT D*V = LMIN*V.
C
         IF (ABS(V1-LMIN).GE.ABS(V2)) THEN
            S(IVAR1) = LMIN*V2/(LMIN-V1)
            S(IVAR2) = LMIN
         ELSE
            S(IVAR1) = LMIN
            S(IVAR2) = LMIN*(LMIN-V1)/V2
         END IF
         SAS = (S(IVAR1)**2 + S(IVAR2)**2)/LMIN
      ELSE
         S(IVAR1) = LMIN
         SAS      = LMIN
      END IF
C
C
C     COMPUTE S AS THE SOLUTION OF THE UPPER TRIANGULAR SYSTEM
C                     T
C                 (D*L *P) * S = (D*V)
C
C     IAPOS : RUNNING POINTER TO CURRENT PIVOT POSITION IN ARRAY A.
C     IPOS  : RUNNING POINTER TO BEGINNING OF CURRENT BLOCK PIVOT ROW.
C
      IAPOS = LATOP + 1
      NROWS = 0
      IBLK = NBLK + 1
C
C     RUN THROUGH BLOCK PIVOT ROWS IN THE REVERSE ORDER.
C
      DO 120 LOOP=1,N
        IF (NROWS.GT.0) GO TO 50
C
C       CONSIDER NEXT BLOCK PIVOT ROW.
C
        IBLK = IBLK - 1
C
C       STOPPING CRITERIA.
C
        IF (IBLK.LT.1) GO TO 130
        IPOS = IW1(IBLK)
C
C       ABS(NCOLS) IS NUMBER OF VARIABLES (COLUMNS) IN BLOCK PIVOT
C       ROW.
C       NROWS IS NUMBER OF ROWS IN BLOCK PIVOT ROW.
C
        NCOLS = -IW(IPOS)
        NROWS = 1
        IF (NCOLS.GT.0) GO TO 40
        NCOLS = -NCOLS
        IPOS = IPOS + 1
        NROWS = IW(IPOS)
   40   JPOS = IPOS + NROWS
        J2 = IPOS + NCOLS
C
C       CONSIDER NEXT PIVOT IN THE IBLK-TH BLOCK PIVOT ROW.
C
  50    IF (NROWS.EQ.1) GO TO 60
C
C       JUMP IF WE HAVE A 2 BY 2 PIVOT.
C
        IF (IW(JPOS-1).LT.0) GO TO 90
C
C       PERFORM BACK-SUBSTITUTION USING 1 BY 1 PIVOT.
C
  60    NROWS = NROWS - 1
        IAPOS = IAPOS - (J2-JPOS+1)
        IIS = IW(JPOS)
        W1 = S(IIS)*A(IAPOS)
        J1 = JPOS + 1
        K = IAPOS + 1
        DO 70 J=J1,J2
          IS = ABS(IW(J)+0)
          W1 = W1 + A(K)*S(IS)
          K = K + 1
  70    CONTINUE
        S(IIS) = W1
        JPOS = JPOS - 1
        GO TO 120
C
C       PERFORM OPERATIONS WITH 2 BY 2 PIVOT
C
  90    NROWS = NROWS - 2
        IAPOS2 = IAPOS - (J2-JPOS+1)
        IAPOS = IAPOS2 - (J2-JPOS+2)
        I1S = -IW(JPOS-1)
        I2S = IW(JPOS)
        W1 = S(I1S)*A(IAPOS) + S(I2S)*A(IAPOS+1)
        W2 = S(I1S)*A(IAPOS+1) + S(I2S)*A(IAPOS2)
        J1 = JPOS + 1
        JJ1 = IAPOS + 2
        JJ2 = IAPOS2 + 1
        DO 100 J=J1,J2
          IS = ABS(IW(J)+0)
          W1 = W1 + S(IS)*A(JJ1)
          W2 = W2 + S(IS)*A(JJ2)
          JJ1 = JJ1 + 1
          JJ2 = JJ2 + 1
  100   CONTINUE
        S(I1S) = W1
        S(I2S) = W2
        JPOS = JPOS - 2
  120 CONTINUE
  130 CONTINUE
1000  FORMAT(' Message from -NEGCUR-',/,
     *       ' Increase the parameter -LIW1- to ',I5)
      RETURN
      END
C     ******************************************************************
C
CCS    SUBROUTINE PRINT( IW, A )
CCD    SUBROUTINE PRINT( IW, A )
C
C     ******************************************************************
C
C      INTEGER IW( * )
CCS    REAL             A( * )
CCD    DOUBLE PRECISION A( * )
C
C     ******************************************************************
C
C     PROGRAMMING :    M. LESCRENIER ( AUG 1987 ).
C     =============
C
C     DESCRIPTION :
C     =============
C
C     PRINT THE FACTORS OF A (COMPUTED BY MA27)
C
C     PARAMETERS :
C     ============
C
C     IW
C            I: ABS(IW(1)) = NUMBER OF BLOCK PIVOTS.
C               IW(1) < 0 MEANS THAT AT LEAST ONE 2X2 PIVOT HAS BEEN
C                         USED.
C               INTEGER INFORMATION ON EACH BLOCK PIVOT FOLLOWS.
C               FOR EACH BLOCK PIVOT, WE STORE
C                   1) THE NUMBER OF COLUMNS
C                   2) THE NUMBER OF ROWS
C                   3) THE COLUMN INDICES
C               IF THE BLOCK PIVOT IS REDUCED TO 1 ROW, WE STORE
C                   1) -(THE NUMBER OF COLUMNS)
C                   3) THE COLUMN INDICES
C               THE FIRST COLUMN INDEX OF A 2X2 PIVOT IS FLAGGED
C               NEGATIVE.
C            O: UNMODIFIED.
C
C     A
C            I: THE FACTORS STORED ROW BY ROW,
C               THE INVERSES OF THE PIVOTS ARE STORED INSTEAD OF THE
C                   PIVOTS THEMSELVES,
C               THE OTHER ENTRIES ARE NEGATED.
C            O: UNMODIFIED.
C
C     ******************************************************************
C
C     INTEGER NBLK, NCOLS, NROWS, IPOS, IAPOS, IBLK, J1, J2, J
C     INTEGER LG, IROWS
C     INTRINSIC ABS
C     NBLK  = ABS( IW( 1 ) + 0 )
C     IPOS  = 2
C     IAPOS = 1
C     DO 20 IBLK = 1, NBLK
C        IF ( IW( IPOS ) .GT. 0 ) THEN
C           NCOLS = IW( IPOS )
C           NROWS = IW( IPOS + 1 )
C           J1 = IPOS + 2
C        ELSE
C           NCOLS = - IW( IPOS )
C           NROWS = 1
C           J1 = IPOS + 1
C        END IF
C        WRITE( 6, 1000 ) IBLK, NROWS, NCOLS
C        J2 = J1 + NCOLS - 1
C        IPOS = J2 + 1
C        WRITE( 6, 2000 ) ( IW( J ), J = J1, J2)
C        WRITE(6,3000)
C        LG = NCOLS
C        DO 10 IROWS = 1, NROWS
C           J1 = IAPOS
C           J2 = IAPOS + LG - 1
C           WRITE( 6, 4000 ) ( A( J ), J = J1, J2 )
C           LG = LG - 1
C           IAPOS = J2 + 1
C  10    CONTINUE
C  20 CONTINUE
C1000 FORMAT(/, ' Block pivot = ',I3, ' NROWS = ', I3, ' NCOLS = ', I3 )
C2000 FORMAT(' Column indices = ', 10I6, / , ( 18X, 10I6 ) )
C3000 FORMAT(' Real entries ... each row starts on a new line' )
C4000 FORMAT( 1P, 5D16.8 )
C     RETURN
C     END
C     ******************************************************************
C
CS    SUBROUTINE SSYSCH( N, IW, A, EPSMCH, NULROW, INFORM )
CD    SUBROUTINE DSYSCH( N, IW, A, EPSMCH, NULROW, INFORM )
C
C     ******************************************************************
C
      INTEGER N, INFORM
      INTEGER IW( * ), NULROW( N )
CS    REAL             A( * ), EPSMCH
CD    DOUBLE PRECISION A( * ), EPSMCH
C
C     ******************************************************************
C
C     PROGRAMMING :    M. LESCRENIER ( AUG 1987 ).
C     =============
C
C     DESCRIPTION :
C     =============
C
C     GIVEN THE GAUSSIAN FACTORIZED FORM OF A SYMMETRIC MATRIX A
C
C                                        T        T
C                                   P A P  = L D L
C
C     WHERE D IS A BLOCK DIAGONAL MATRIX WITH 1X1 OR 2X2 BLOCKS
C     SET INFORM TO 1 IF THE MATRIX IS POSITIVE DEFINITE,
C                   2 IF THE MATRIX IS INDEFINITE,
C                   3 IF THE MATRIX IS SINGULAR.
C
C     PARAMETERS :
C     ============
C
C     N
C            I: ORDER OF THE MATRIX A.
C            O: UNMODIFIED.
C
C     IW
C            I: ABS(IW(1)) = NUMBER OF BLOCK PIVOT ROWS.
C               IW(1) < 0 MEANS THAT AT LEAST ONE 2X2 PIVOT HAS BEEN
C                         USED.
C               INTEGER INFORMATION ON EACH BLOCK PIVOT ROW FOLLOWS.
C               FOR EACH BLOCK PIVOT ROW, WE STORE
C                   1) THE NUMBER OF COLUMNS
C                   2) THE NUMBER OF ROWS
C                   3) THE COLUMN INDICES
C               IF THE BLOCK PIVOT ROW IS REDUCED TO 1 ROW, WE STORE
C                   1) -(THE NUMBER OF COLUMNS)
C                   3) THE COLUMN INDICES
C               THE FIRST COLUMN INDEX OF A 2X2 PIVOT IS FLAGGED
C               NEGATIVE.
C            O: UNMODIFIED.
C
C     A
C            I: THE FACTORS STORED ROW BY ROW,
C               THE INVERSES OF THE PIVOTS ARE STORED INSTEAD OF THE
C                   PIVOTS THEMSELVES,
C               THE OTHER ENTRIES ARE NEGATED.
C            O: UNMODIFIED.
C
C     EPSMCH
C            I: RELATIVE PRECISION OF THE MACHINE
C            O: UNMODIFIED.
C
C     INFORM
C            I: UNDEFINED.
C            O: 1 IF THE MATRIX IS POSITIVE DEFINITE,
C               2 IF THE MATRIX IS INDEFINITE,
C               3 IF THE MATRIX IS SINGULAR.
C
C     NULROW
C            I: UNDEFINED.
C            O: IF THE I-TH ROW OF THE MATRIX IS ENTIRELY NUL,
C               NULROW(I) = 1;
C               NULROW(I) = 0 OTHERWISE
C
C     ******************************************************************
C
      LOGICAL SINGUL, POSDEF
      INTEGER I, IPOS, IAPOS, J1, NVAR, NBLK, NCOLS, NROWS
      INTEGER J2, J, LG, IBLK
CS    REAL             L, DET, RHO, AUX, TOLER, V1, V2, V3, ZERO, ONE
CD    DOUBLE PRECISION L, DET, RHO, AUX, TOLER, V1, V2, V3, ZERO, ONE
      INTRINSIC ABS
CS    DATA ZERO, ONE / 0.0E+0, 1.0E+0 /
CD    DATA ZERO, ONE / 0.0D+0, 1.0D+0 /
      DO 5 I = 1, N
         NULROW( I ) = 1

5     CONTINUE
      NVAR   = 0
      SINGUL = .FALSE.
      POSDEF = .TRUE.
      TOLER  = 1.0D+4 * EPSMCH
      NBLK   = ABS( IW( 1 ) + 0 )
      IPOS   = 2
      IAPOS  = 1
      DO 20 IBLK = 1, NBLK
C
C        CONSIDER THE IBLK-TH BLOCK PIVOT.
C
         IF ( IW( IPOS ) .GT. 0 ) THEN
            NCOLS = IW( IPOS )
            NROWS = IW( IPOS + 1 )
            J1    = IPOS + 2
         ELSE
            NCOLS = - IW( IPOS )
            NROWS = 1
            J1    = IPOS + 1
         END IF
         NVAR = NVAR + NROWS
         IF ( IBLK .EQ. NBLK .AND. NROWS .NE. NCOLS ) SINGUL = .TRUE.
         J2   = J1 + NROWS - 1
         IPOS = J1 + NCOLS
         LG   = NCOLS
         J   = J1
C
C  CONSIDER THE NEXT PIVOT (1X1 OR 2X2) OF THIS BLOCK PIVOT.
C
   10    CONTINUE
         IF ( IW( J ) .GT. 0 ) THEN
C
C  IT IS A 1X1 PIVOT.
C
            NULROW( IW( J ) ) = 0
C           WRITE( 6, * ) ' A(IAPOS) = ', A(IAPOS)
            IF ( ABS( A( IAPOS ) ) * TOLER .GT. ONE ) THEN
               SINGUL     = .TRUE.
               A( IAPOS ) = ZERO
            END IF
            IF ( A( IAPOS ) .LE. ZERO ) POSDEF = .FALSE.
            IAPOS = IAPOS + LG
            LG    = LG - 1
            J     = J + 1
         ELSE
C
C  IT IS A 2X2 PIVOT.
C
            NULROW( - IW( J ) )   = 0
            NULROW( IW( J + 1 ) ) = 0
            V1 = A( IAPOS )
            V2 = A( IAPOS + 1 )
            IAPOS = IAPOS + LG
            LG = LG - 1
            V3 = A( IAPOS )
            IAPOS = IAPOS + LG
            LG = LG - 1
C
C  COMPUTE THE 2X2 PIVOT.
C
            DET = ONE / ( V1 * V3 - V2 ** 2 )
            AUX = V3
            V3  = V1 * DET
            V1  = AUX * DET
            V2  = - V2 * DET
C
C           COMPUTE THE SMALLEST EIGENVALUE OF THIS 2X2 PIVOT.
C
            RHO = ( V1 + V3 ) ** 2 - 4.0D0 * ( V1 * V3 - V2 ** 2 )
            L   = 0.5D0 * ( V1 + V3 - SQRT( RHO ) )
C           WRITE( 6, * ) ' L = ', L
            IF ( ABS( L ) .LT. TOLER ) THEN
               SINGUL = .TRUE.
               L      = ZERO
            END IF
            IF ( L .LE. ZERO ) POSDEF = .FALSE.
            J = J + 2
         END IF
         IF ( J .LE. J2 ) GO TO 10
   20 CONTINUE
      IF ( NVAR .NE. N ) THEN
C
C  THE MATRIX HAS ENTIRELY NUL ROWS.
C
         SINGUL = .TRUE.
      END IF
      IF ( SINGUL ) THEN
         INFORM = 3
      ELSE
         IF ( POSDEF ) THEN
            INFORM = 1
         ELSE
            INFORM = 2
         END IF
      END IF
      RETURN
      END
