C  THIS VERSION: 14/08/1991 AT 11:58:46 AM.
CS    SUBROUTINE SDRCHG( NG, X , GVALS , LGVALS, XT    , GTVALS, LGTVAL,
CD    SUBROUTINE DDRCHG( NG, X , GVALS , LGVALS, XT    , GTVALS, LGTVAL,
     *                   ITESTG, NTESTG, RELPR , IPRINT, IOUT  ,
     *                   WARNNG, DEBUG , JUMPTO )
      INTEGER            NG    , LGVALS, LGTVAL, NTESTG, JUMPTO
      INTEGER            IOUT  , IPRINT
CS    REAL               RELPR
CD    DOUBLE PRECISION   RELPR
      LOGICAL            WARNNG, DEBUG
      INTEGER            ITESTG( NG )
CS    REAL               X( NG ), XT( NG ), GVALS( LGVALS, 3 ),
CD    DOUBLE PRECISION   X( NG ), XT( NG ), GVALS( LGVALS, 3 ),
     *                   GTVALS( LGTVAL, 3 )
C
C  GIVEN A VECTOR OF FUNCTIONALS GVALS, CHECK THEIR ANALYTICAL
C  DERIVATIVES AGAINST APPROXIMATIONS BY
C  DIFFERENCES AT THE GIVEN POINT X.
C
C  NICK GOULD, 26 OCTOBER 1988.
C  FOR CGT PRODUCTIONS.
C
      INTEGER          I, J
CS    REAL             EPSQRT, COMP, GTOL, ONE
CD    DOUBLE PRECISION EPSQRT, COMP, GTOL, ONE
C
C  MACHINE FUNCTIONS.
C
      INTRINSIC        SQRT, MAX, ABS
      SAVE
C
C  SET CONSTANT REAL PARAMETERS.
C
CS    PARAMETER ( ONE = 1.0E+0 )
CD    PARAMETER ( ONE = 1.0D+0 )
      IF ( JUMPTO .GT. 0 ) GO TO ( 100, 200, 300 ), JUMPTO
      IF ( IOUT .GT. 0 .AND. IPRINT .GT. 0 ) WRITE( IOUT, 2000 )
      WARNNG = .FALSE.
CS    EPSQRT =  SQRT( RELPR )
CD    EPSQRT = DSQRT( RELPR )
C
C  RETURN TO OBTAIN THE VALUES OF THE FUNCTIONS AND THEIR DERIVATIVES
C  AT THE POINT X.
C
      JUMPTO = 1
      RETURN
  100 CONTINUE
C
C  COPY THE COMPONENTS OF GVALS FOR GROUPS THAT ARE
C  TO BE TESTED INTO GTVALS.
C
      DO 110 J = 1, NTESTG
         I     = ITESTG( J )
         GTVALS( I, 1 ) = GVALS( I, 1 )
         GTVALS( I, 2 ) = GVALS( I, 2 )
         GTVALS( I, 3 ) = GVALS( I, 3 )
C
C  COPY THE COMPONENTS OF X FOR GROUPS THAT ARE TO BE
C  TESTED INTO XT AND PERTURB X BY A SMALL QUANTITY.
C
         XT( I ) = X( I )
         X( I )  = X( I ) + EPSQRT
  110 CONTINUE
C
C  EVALUATE THE REQUIRED GROUPS AT THE PERTURBED POINT.
C
      JUMPTO = 2
      RETURN
  200 CONTINUE
C
C  ESTIMATE THE FIRST DERIVATIVE OF THE ITH GROUP AND
C  TEST IT W.R.T. ITS ANALYTICAL VALUE.
C
      DO 210 J = 1, NTESTG
         I     = ITESTG( J )
         COMP  = ( GVALS( I, 1 ) - GTVALS( I, 1 ) ) / EPSQRT
C
C  TEST AGREEMENT BETWEEN ANALYTICAL AND APPROX. VALUES.
C
CS       GTOL = 1.0E+2 * EPSQRT * MAX( ONE, ABS( GTVALS( I, 2 ) ) )
CD       GTOL = 1.0D+2 * EPSQRT * MAX( ONE, ABS( GTVALS( I, 2 ) ) )
         IF ( DEBUG .OR.
     *        ABS( GTVALS( I, 2 ) - COMP ) .GT. GTOL ) THEN
            IF ( ABS( GTVALS( I, 2 ) - COMP ) .GT. GTOL ) THEN
               WARNNG = .TRUE.
               IF ( IOUT .GT. 0 .AND. IPRINT .GT. 0 )
     *              WRITE( IOUT, 2040 )
            END IF
            IF ( IOUT .GT. 0 .AND. IPRINT .GT. 0 ) THEN
               WRITE( IOUT, 2010 ) I
               WRITE( IOUT, 2020 ) GTVALS( I, 2 ), COMP
            END IF
         END IF
  210 CONTINUE
C
C  EVALUATE THE REQUIRED GROUP DERIVATIVES AT THE
C  PERTURBED POINT.
C
      JUMPTO = 3
      RETURN
  300 CONTINUE
C
C  ESTIMATE THE SECOND DERIVATIVE OF THE I-TH GROUP AND
C  TEST IT W.R.T. ITS ANALYTICAL VALUE.
C
      DO 310 J = 1, NTESTG
         I     = ITESTG( J )
         COMP = ( GVALS( I, 2 ) - GTVALS( I, 2 ) ) / EPSQRT
C
C  TEST AGREEMENT BETWEEN ANALYTICAL AND APPROX. VALUES.
C
CS       GTOL = 1.0E+2 * EPSQRT * MAX( ONE, ABS( GTVALS( I, 3 ) ) )
CD       GTOL = 1.0D+2 * EPSQRT * MAX( ONE, ABS( GTVALS( I, 3 ) ) )
         IF ( DEBUG .OR.
     *        ABS( GTVALS( I, 3 ) - COMP ) .GT. GTOL ) THEN
            IF ( ABS( GTVALS( I, 3 ) - COMP ) .GT. GTOL ) THEN
                WARNNG = .TRUE.
                IF ( IOUT .GT. 0 .AND. IPRINT .GT. 0 )
     *               WRITE( IOUT, 2050 )
            END IF
            IF ( IOUT .GT. 0 .AND. IPRINT .GT. 0 ) THEN
               WRITE( IOUT, 2030 ) I
               WRITE( IOUT, 2020 ) GTVALS( I, 3 ), COMP
            END IF
         END IF
  310 CONTINUE
C
C  RESET THE COMPONENTS OF GVALS AND X FOR GROUPS THAT HAVE BEEN TESTED.
C
      DO 610 J         = 1, NTESTG
         I             = ITESTG( J )
         GVALS( I, 1 ) = GTVALS( I, 1 )
         GVALS( I, 2 ) = GTVALS( I, 2 )
         GVALS( I, 3 ) = GTVALS( I, 3 )
         X( I )        = XT( I )
  610 CONTINUE
      IF ( IOUT .GT. 0 .AND. IPRINT .GT. 0 ) THEN
         IF ( WARNNG ) THEN
            WRITE( IOUT, 2060 )
         ELSE
            WRITE( IOUT, 2070 )
         END IF
      END IF
      JUMPTO = 0
      RETURN
C
C  NON-EXECUTABLE STATEMENTS.
C
 2000 FORMAT( /, ' *********** Checking group derivatives **********',
     *        /, ' *                                               *',/)
 2010 FORMAT(' 1st derivative of group function ', I4, ' :',
     *       /,42('-'))
 2020 FORMAT(' Anal. deriv. = ',1P, D14.6,' Approx. deriv. = ',
     *       1P, D14.6, / )
 2030 FORMAT(' 2nd derivative of group function ', I4, ' :',
     *       /,42('-'))
 2040 FORMAT(' Possible mistake in the computation of the 1st',
     *       ' derivative')
 2050 FORMAT(' Possible mistake in the computation of the 2nd',
     *       ' derivative')
 2060 FORMAT( ' *                                               * ', /,
     *        ' *********** Derivatives checked - warnings ****** ' )
 2070 FORMAT( ' *                                               * ', /,
     *        ' *********** Derivatives checked O.K. ************ ' )
C
C  END OF DRCHG.
C
      END
