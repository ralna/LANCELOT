C ** Correction report.
C ** Correction 1. 13/01/94: 2 lines added
C ** End of Correction report.
C  THIS VERSION: 13/01/1994 AT 04:37:34 PM.
CS    SUBROUTINE SDRCHE( N, NEL, ISTAEV, LSTAEV, ISTADH, LSTADH, IELVAR, 
CD    SUBROUTINE DDRCHE( N, NEL, ISTAEV, LSTAEV, ISTADH, LSTADH, IELVAR,
     *                   LELVAR, IELVRT, LELVRT, INTVAR, LNTVAR,
     *                   IWORK , LIWORK, INTREP, LINTRE, ICALCF, LCALCF, 
     *                   NCALCF, X , LX, XT    , LXT   , FUVALS,
     *                   LFUVAL, FTUVAL, LFTUVA, XEL   , LXEL  , XINT  ,
     *                   LXINT , WORK  , LWORK , RELPR , SECOND, ITESTL,
     *                   LTESTL, IPRINT, IOUT  , RANGES, WARNNG, 
     *                   DEBUG , JUMPTO )
      INTEGER            N, NEL, NCALCF, LFUVAL, JUMPTO, IPRINT, IOUT
      INTEGER            LSTAEV, LSTADH, LELVAR, LELVRT, LNTVAR, LINTRE
      INTEGER            LCALCF, LIWORK, LWORK , LX    , LXT   , LTESTL
      INTEGER            LFTUVA, LXEL  , LXINT
CS    REAL               RELPR
CD    DOUBLE PRECISION   RELPR
      LOGICAL            SECOND, WARNNG, DEBUG
      INTEGER            ISTADH( LSTADH ), INTVAR( LNTVAR )
      INTEGER            ISTAEV( LSTAEV ), IELVAR( LELVAR )
      INTEGER            ICALCF( LCALCF ), ITESTL( LTESTL )
      INTEGER            IWORK ( LIWORK ), IELVRT( LELVRT )
      LOGICAL            INTREP( LINTRE )
CS    REAL               X     ( LX     ), XT    ( LXT    ),
CD    DOUBLE PRECISION   X     ( LX     ), XT    ( LXT    ),
     *                   FTUVAL( LFTUVA ), XEL   ( LXEL   ),
     *                   XINT  ( LXINT  ), WORK  ( LWORK  ),
     *                   FUVALS( LFUVAL )
      EXTERNAL           RANGES
C
C     ******************************************************************
C
C     PROGRAMMING :    N. GOULD AND M. LESCRENIER ( OCTOBER 1987 ).
C     =============
C
C     DESCRIPTION :
C     =============
C
C     GIVEN A PARTIALLY SEPARABLE FUNCTION, CHECK THE ANALYTICAL
C     GRADIENTS (AND HESSIANS IF REQUIRED) AGAINST APPROXIMATIONS BY
C     DIFFERENCES AT THE GIVEN POINT X.
C
C     PARAMETERS :
C     ============
C
C     N
C            I: NUMBER OF VARIABLES OF THE PARTIALLY SEPAR. FUNCTION.
C            O: UNMODIFIED.
C
C     NEL
C            I: NUMBER OF ELEMENTS OF THE PARTIALLY SEPAR. FUNCTION.
C            O: UNMODIFIED.
C
C     ISTAEV
C            I: POINTERS TO THE LISTS OF VARIABLES OF THE ELEMENTS
C               (STORED IN INVAR).
C            O: UNMODIFIED.
C
C     ISTADH
C            I: ISTADH(I)-ISTADH(1)+1 POINTS TO THE I-TH ELEMENT
C               HESSIAN IN HX
C            O: UNMODIFIED.
C
C     IELVAR I: LISTS OF VARIABLES OF THE ELEMENTS
C            O: UNMODIFIED.
C
C     INTVAR
C            I: INTVAR(I) GIVES THE NUMBER OF INTERNAL VARIABLES
C               IN ELEMENT I.
C            O: UNMODIFIED.
C
C     ICALCF
C            I: INTEGER VECTOR WHOSE COMPONENTS GIVE THE ELEMENTS TESTED
C            O: UNMODIFIED.
C
C     X
C            I: POINT AT WHICH THE GRADIENTS MUST BE CHECKED.
C            O: UNMODIFIED.
C
C
C     RELPR
C            I: RELATIVE PRECISION OF THE MACHINE.
C            O: UNMODIFIED.
C
C     SECOND
C            I: .FALSE. IF ONLY THE ANALYTICAL GRADIENTS MUST BE
C                          CHECKED.
C               .TRUE.  IF THE ANALYTICAL GRADIENTS AND HESSIANS MUST
C                          BE CHECKED.
C            O: UNMODIFIED.
C
C     ******************************************************************
C
      INTEGER          I , II, IELF  , INFORM, IP1, L, J , K , KK    , 
     *                 LFXI  , LG, LH, LGXI  , LHXI  , NELVAR, NEL1  ,
     *                 NINVAR, NSIZEH, NTESTL, LEND  , ITEST
CS    REAL             EPSQRT, COMP  , GTOL  , ZERO  , ONE
CD    DOUBLE PRECISION EPSQRT, COMP  , GTOL  , ZERO  , ONE
      LOGICAL          INTRE
C
C  MACHINE FUNCTIONS.
C
      INTRINSIC        SQRT  , MAX1  , ABS
CS    REAL             SDOT
CD    DOUBLE PRECISION DDOT
CS    EXTERNAL         SCOPY , SSETVL, SDOT  , SGNINV
CD    EXTERNAL         DCOPY , DSETVL, DDOT  , DGNINV
      SAVE
C
C  SET CONSTANT REAL PARAMETERS.
C
CS    PARAMETER      ( ZERO = 0.0E+0, ONE = 1.0E+0 )
CD    PARAMETER      ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      IF ( JUMPTO .GT. 0 ) GO TO ( 100, 240, 250 ), JUMPTO
      IF ( IOUT .GT. 0 .AND. IPRINT .GT. 0 ) WRITE( IOUT, 2000 )
      WARNNG = .FALSE.
      EPSQRT = SQRT( RELPR )
C
C  SET THE STARTING ADDRESSES FOR THE PARTITIONS WITHIN FUVALS.
C
      LFXI = 0
      LGXI = LFXI + NEL
C
C  SET UP THE STARTING ADDRESSES FOR THE ELEMENT GRADIENTS
C  WITH RESPECT TO THEIR INTERNAL VARIABLES.
C
      NEL1             = NEL + 1
      INTVAR( NEL1 )   = 0
      K                = INTVAR( 1 )
      INTVAR( 1 )      = LGXI + 1
      DO 10 I          = 1, NEL
         IP1           = I + 1
         L             = INTVAR( IP1 )
         INTVAR( IP1 ) = INTVAR( I ) + K
         K             = L
   10 CONTINUE
C
C  ENSURE THAT ALL THE ELEMENT FUNCTIONS ARE EVALUATED
C  AT THE INITIAL POINT.
C
      NTESTL         = NCALCF
      DO 15 I        = 1, NCALCF
         ICALCF( I ) = ITESTL( I )
   15 CONTINUE
      LHXI   = INTVAR( NEL1 ) - 1
      NINVAR = INTVAR( NEL1 ) - INTVAR( 1 )
C
C  SET UP THE STARTING ADDRESSES FOR THE ELEMENT HESSIANS
C  WITH RESPECT TO THEIR INTERNAL VARIABLES.
C
      K = LHXI + 1
      IF ( SECOND ) THEN
         DO 20 I        = 1, NEL
            ISTADH( I ) = K
            NSIZEH      = INTVAR( I + 1 ) - INTVAR( I )
            K           = K + NSIZEH * ( NSIZEH + 1 ) / 2
   20    CONTINUE
      END IF
      LEND = K - 1
C
C  INITIALIZE FUVALS TO ZERO.
C
CS    CALL SSETVL( LEND, FUVALS, 1, ZERO )
CD    CALL DSETVL( LEND, FUVALS, 1, ZERO )
      JUMPTO = 1
      RETURN
  100 CONTINUE
C
C  COPY FUVALS INTO FTUVAL.
C
CS    CALL SCOPY( LEND, FUVALS, 1, FTUVAL, 1 )
CD    CALL DCOPY( LEND, FUVALS, 1, FTUVAL, 1 )
C
C  CHECK THE ANALYTICAL GRADIENT (IN INTERNAL REPRESENTATION)
C
      JUMPTO = 2
C
C  MOCK DO LOOP TO ALLOW REVERSE COMMUNICATION.
C
      ITEST  = 0
      NCALCF = 1
  200 CONTINUE
         ITEST = ITEST + 1
         IF ( ITEST .GT. NTESTL ) GO TO 290
         IELF = ITESTL( ITEST )
C        IF ( IOUT .GT. 0 .AND. IPRINT .GT. 0 ) WRITE( IOUT, 2010 ) IELF
         NINVAR         = INTVAR( IELF + 1 ) - INTVAR( IELF )
         NELVAR         = ISTAEV( IELF + 1 ) - ISTAEV( IELF )
         K              = 0
         KK             = INTVAR( IELF )
         INTRE          = INTREP( IELF )
         ICALCF( 1 )    = IELF
C
C  ASSIGN TEMPORARY INDICES TO THE VARIABLES.
C  STORE THE ELEMENTAL VARIABLES IN XEL.
C
         II              = ISTAEV( IELF )
         DO 210 I        = 1, NELVAR
            IELVRT( II ) = I
            XEL( I )     = X( IELVAR( II ) )
            II           = II + 1
  210    CONTINUE
         IF ( INTRE ) THEN
C
C  COMPUTE THE VALUES OF THE INTERNAL VARIABLES.
C
            CALL RANGES( IELF, .FALSE., XEL, XINT, NELVAR, NINVAR )
CS          CALL SGNINV( IELF, WORK, NELVAR, NINVAR, IWORK,
CD          CALL DGNINV( IELF, WORK, NELVAR, NINVAR, IWORK,
     *                   WORK( NELVAR * NINVAR + 1 ), INFORM,
     *                   IOUT, RANGES )
            IF ( INFORM .EQ. 1 ) THEN
               IF ( IOUT .GT. 0 .AND. IPRINT .GT. 0 )
     *              WRITE( IOUT, 2080 ) IELF
               GO TO 280
            END IF
         END IF
         IF ( SECOND ) LH = ISTADH( IELF )
C
C  MOCK DO LOOP TO ALLOW REVERSE COMMUNICATION.
C
         J = 0
  220    CONTINUE
            J  = J + 1
            IF ( J .GT. NINVAR ) GO TO 280
C
C  CHECK THE K-TH COMPONENT OF THE IELF-TH GRADIENT.
C
            IF ( INTRE ) THEN
               XINT( J )   = XINT( J ) + EPSQRT
C
C  PUT THE ELEMENTAL VARIABLES IN XEL.
C
               L           = 1
               DO 230 I    = 1, NELVAR
                  XT( I )  = XEL( I )
CS                XEL( I ) = SDOT( NINVAR, WORK( L ), 1, XINT, 1 )
CD                XEL( I ) = DDOT( NINVAR, WORK( L ), 1, XINT, 1 )
                  L        = L + NINVAR
  230          CONTINUE
            ELSE
               XT( 1 )  = XEL( J )
               XEL( J ) = XEL( J ) + EPSQRT
C              WRITE(6,*) ' VARIABLE ', J, ' WAS ', XT( 1 ),
C    *                    ' BECOMES ', XEL( J )
            END IF
C
C  EVALUATE THE IELF-TH ELEMENT FUNCTION AT THE PERTURBED POINT.
C
            JUMPTO = 2
            RETURN
  240       CONTINUE
C
C  ESTIMATE THE K-TH COMPONENT OF THE GRADIENT AND
C  TEST IT W.R.T. ITS ANALYTICAL VALUE.
C
            COMP = ( FUVALS( IELF ) - FTUVAL( IELF ) ) / EPSQRT
C
C  PRINT THE COMPONENTS FOR COMPARISON.
C
C           IF ( IOUT .GT. 0 .AND. IPRINT .GT. 0 )
C    *           WRITE( IOUT, 2020 ) J, FTUVAL( KK ), J, COMP
C
C  TEST AGREEMENT BETWEEN ANALYTICAL AND APPROX. VALUES.
C
CS          GTOL = 1.0E+2 * EPSQRT * MAX( ONE, ABS( FTUVAL( KK ) ) )
CD          GTOL = 1.0D+2 * EPSQRT * MAX( ONE, ABS( FTUVAL( KK ) ) )
            IF ( DEBUG .OR. ABS( FTUVAL( KK ) - COMP ) .GT. GTOL ) THEN
               IF ( ABS( FTUVAL( KK ) - COMP ) .GT. GTOL ) THEN
                  WARNNG = .TRUE.
                  IF ( IOUT .GT. 0 .AND. IPRINT .GT. 0 )
     *                 WRITE( IOUT, 2050 )
               END IF
               IF ( IOUT .GT. 0 .AND. IPRINT .GT. 0 ) THEN
                  WRITE( IOUT, 2010 ) IELF
                  WRITE( IOUT, 2020 ) J, FTUVAL( KK ), J, COMP
               END IF
            END IF
            KK = KK + 1
C
C  CHECK THE ANALYTICAL HESSIAN (IN INTERNAL REPRESENTATION)
C
            IF ( SECOND ) THEN
C
C  COMPUTE THE IELF-TH PERTURBED GRADIENT IN INTERNAL REPRESENTATION
C
               JUMPTO = 3
               RETURN
            END IF
  250       CONTINUE
            IF ( SECOND ) THEN
C
C  ESTIMATE EACH COMPONENT OF THE HESSIAN'S COLUMN AND
C  TEST IT W.R.T. ITS ANALYTICAL VALUE
C
               LG       = INTVAR( IELF )
               DO 260 K = 1, J
                  COMP  = ( FUVALS( LG ) - FTUVAL( LG ) ) / EPSQRT
C
C  PRINT THE COMPONENTS FOR COMPARISON
C
C                 IF ( IOUT .GT. 0 .AND. IPRINT .GT. 0 )
C     *              WRITE( IOUT, 2040 ) K, J,
C     *                                  FTUVAL( LH ), K, J, COMP
C
C  TEST AGREEMENT BETWEEN ANALYTICAL AND APPROX. VALUES
C
CS                GTOL = 1.0E+2 * EPSQRT * 
CD                GTOL = 1.0D+2 * EPSQRT *
     *                   MAX( ONE, ABS( FTUVAL( LH ) ) )
                  IF ( ABS( FTUVAL( LH ) - COMP ) .GT. GTOL
     *                 .OR. DEBUG ) THEN
                     IF ( ABS( FTUVAL( LH ) - COMP ) .GT. GTOL ) THEN
                        WARNNG = .TRUE.
                        IF ( IOUT .GT. 0 .AND. IPRINT .GT. 0 )
     *                       WRITE( IOUT, 2060 )
                     END IF
                     IF ( IOUT .GT. 0 .AND. IPRINT .GT. 0 ) THEN
                        WRITE( IOUT, 2030 ) IELF
                        WRITE( IOUT, 2040 ) K, J, FTUVAL( LH ),
     *                                      K, J, COMP
                     END IF
                  END IF
                  LG = LG + 1
                  LH = LH + 1
  260          CONTINUE
C
C  RESET THE POINT X.
C
            END IF
            IF ( INTRE ) THEN
               XINT( J ) = XINT( J ) - EPSQRT
               DO 270 I  = 1, NELVAR
                  XEL( I ) = XT( I )
  270           CONTINUE
            ELSE
               XEL( J ) = XT( 1 )
            END IF
            GO TO 220
  280    CONTINUE
         GO TO 200
  290 CONTINUE
C
C  PREPARE TO EXIT. RESET THE POINTERS TO NUMBER OF INTERNAL VARIABLES.
C
      DO 510 I       = 1, NEL
         INTVAR( I ) = INTVAR( I + 1 ) - INTVAR( I )
  510 CONTINUE
C
C  RESET FUVALS TO FTUVAL.
C
CS    CALL SCOPY( LEND, FTUVAL, 1, FUVALS, 1 )
CD    CALL DCOPY( LEND, FTUVAL, 1, FUVALS, 1 )
      IF ( WARNNG ) THEN
         IF ( IOUT .GT. 0 .AND. IPRINT .GT. 0 ) WRITE( IOUT, 2090 )
      ELSE
         IF ( IOUT .GT. 0 .AND. IPRINT .GT. 0 ) WRITE( IOUT, 2070 )
      END IF
      JUMPTO = 0
      RETURN
C
C  NON-EXECUTABLE STATEMENTS.
C
 2000 FORMAT( /, ' ********* Checking element derivatives **********',
     *        /, ' *                                               *',/)
 2010 FORMAT( ' Gradient of element function ', I4, ' :', /, 42('-') )
 2020 FORMAT( ' Anal. grad.(', I4, ')      = ', 1P, D14.6,
     *        ' Approx.(', I4, '    )      = ', 1P, D14.6, / )
 2030 FORMAT( ' Hessian of element function ', I4, ' :', /, 41('-') )
 2040 FORMAT( ' Anal. hess.(', I4, ',', I4, ') = ', 1P, D14.6,
     *        ' Appr. hess.(', I4, ',', I4, ') = ', 1P, D14.6, / )
 2050 FORMAT( ' Possible mistake in the computation of the gradient ')
 2060 FORMAT( ' Possible mistake in the computation of the Hessian  ')
 2070 FORMAT( ' *                                               * ', /,
     *        ' *********** Derivatives checked O.K. ************ ' )
 2080 FORMAT( ' The range transformation for element ', I5, ' is',
     *        ' singular ' )
 2090 FORMAT( ' *                                               * ', /,
     *        ' ******** Derivatives checked - warnings ********* ' )
C
C  END OF DRCHE.
C
      END
C  THIS VERSION: 13/01/1994 AT 04:37:34 PM.
CS    SUBROUTINE SGNINV( IELF  , A     , NELVAR, NINVAR, IWORK , WORK  ,
CD    SUBROUTINE DGNINV( IELF  , A     , NELVAR, NINVAR, IWORK , WORK  ,
     *                   INFORM, IOUT  , RANGES )
      INTEGER            IELF  , NELVAR, NINVAR, INFORM, IOUT
      INTEGER            IWORK ( *    )
CS    REAL               A     ( NINVAR, NELVAR )     , WORK( * )
CD    DOUBLE PRECISION   A     ( NINVAR, NELVAR )     , WORK( * )
      EXTERNAL           RANGES
C
C  TO FIND THE INVERSE TRANSFORMATION TO THE "GATHER" USED BY RANGES.
C
C  NICK GOULD, 15TH MARCH 1990.
C  FOR CGT PRODUCTIONS.
C
      INTEGER           I , IR, J , K , KP1, L
CS    REAL              SIZEA , SNRM2 , ZERO  , HALF  , ONE
CD    DOUBLE PRECISION  SIZEA , DNRM2 , ZERO  , HALF  , ONE
CS    REAL              BSQ   , RMAX  , SIGMA , SUM
CD    DOUBLE PRECISION  BSQ   , RMAX  , SIGMA , SUM
CS    EXTERNAL          SCOPY , SSETVL, SNRM2
CD    EXTERNAL          DCOPY , DSETVL, DNRM2
CS    REAL              EPSMCH, EPSNEG, TINY  , BIG
CD    DOUBLE PRECISION  EPSMCH, EPSNEG, TINY  , BIG
CS    COMMON / SMACHN / EPSMCH, EPSNEG, TINY  , BIG
CD    COMMON / DMACHN / EPSMCH, EPSNEG, TINY  , BIG
C ** Correction 1. 13/01/94: 2 lines added
CS    SAVE   / SMACHN /
CD    SAVE   / DMACHN /
C ** Correction 1. 13/01/94: end of correction **
CS    PARAMETER       ( ZERO = 0.0E+0, HALF = 5.0E-1, ONE = 1.0E+0 )
CD    PARAMETER       ( ZERO = 0.0D+0, HALF = 5.0D-1, ONE = 1.0D+0 )
C
C  ASSEMBLE THE "GATHER" MATRIX INTO A.
C
      INFORM = 0
CS    CALL SSETVL( NINVAR, WORK, 1, ZERO )
CD    CALL DSETVL( NINVAR, WORK, 1, ZERO )
      DO 10 I      = 1, NINVAR
         WORK( I ) = ONE
         CALL RANGES( IELF, .TRUE., WORK, WORK( NINVAR + 1 ),
     *                NELVAR, NINVAR )
CS       CALL SCOPY( NELVAR, WORK( NINVAR + 1 ), 1, A( I, 1 ), NINVAR )
CD       CALL DCOPY( NELVAR, WORK( NINVAR + 1 ), 1, A( I, 1 ), NINVAR )
CS       SIZEA = SNRM2( NELVAR, WORK( NINVAR + 1 ), 1 )
CD       SIZEA = DNRM2( NELVAR, WORK( NINVAR + 1 ), 1 )
         IF ( SIZEA .LT. EPSMCH ) THEN
            INFORM = 1
            RETURN
         END IF
         WORK( I ) = ZERO
   10 CONTINUE
C
C  FORM THE (MOORE-PENROSE) GENERALIZED INVERSE OF A USING THE
C  STORAGE MISERLY METHOD OF POWELL (AERE REPORT R-6072).
C
C  TRANSFORM A TO LOWER TRIANGULAR FORM BY A SEQUENCE OF ELEMENTARY
C  TRANSFORMATIONS WITH ROW AND COLUMN INTERCHANGES FOR STABILITY.
C  RECORD ALL ROW AND COLUMN INTERCHANGES IN IWORK.
C
      DO 20 I       = 1, NINVAR
         IWORK( I ) = I
   20 CONTINUE
      DO 30 I                = 1, NELVAR
         IWORK( NINVAR + I ) = I
   30 CONTINUE   
C
C  PERFORM THE K-TH ELEMENTARY TRANSFORMATION.
C
      DO 200 K = 1, NINVAR
         KP1   = K + 1
C
C  FIND THE LARGEST ROW.
C
         RMAX     = ZERO
         DO 120 I = K, NINVAR
            SUM   = ZERO
            DO 110 J = K, NELVAR
               SUM   = SUM + A( I, J ) ** 2
  110       CONTINUE
            IF ( SUM .GT. RMAX ) THEN
               IR   = I
               RMAX = SUM
            END IF
  120    CONTINUE
C
C  IF THE MATRIX IS NOT OF FULL RANK, STOP.
C
         IF ( RMAX .EQ. ZERO ) THEN
            IF ( IOUT .GT. 0 ) WRITE( IOUT, 2000 ) NINVAR - K
            STOP
         END IF
C
C  IF THE CURRENT ROW IS NOT THE LARGEST, SWOP IT WITH THE LARGEST.
C
         IF ( IR .GT. K ) THEN
            L           = IWORK( K )
            IWORK( K )  = IWORK( IR )
            IWORK( IR ) = L
            DO 130 J      = 1, NELVAR
               SUM        = A( K, J )
               A( K, J )  = A( IR, J )
               A( IR, J ) = SUM
  130       CONTINUE
         END IF
C
C  FIND LARGEST ELEMENT IN THE PIVOTAL ROW.
C
         RMAX     = ZERO
         SUM      = ZERO
         DO 140 J = K, NELVAR
            SUM   = SUM + A( K, J ) ** 2
            IF ( RMAX .LT. ABS( A( K, J ) ) ) THEN
               IR   = J
               RMAX = ABS( A( K, J ) )
            END IF
  140    CONTINUE
C
C  IF THE CURRENT COLUMN IS NOT THE LARGEST, SWOP IT WITH THE LARGEST.
C
         IF ( IR .GT. K ) THEN
            I          = NINVAR + K
            J          = NINVAR + IR
            L          = IWORK( I )
            IWORK( I ) = IWORK( J )
            IWORK( J ) = L
            DO 150 I      = 1, NINVAR
               RMAX       = A( I, K )
               A( I, K )  = A( I, IR )
               A( I, IR ) = RMAX
  150      CONTINUE
         END IF
C
C  REPLACE THE PIVOTAL ROW BY THE HOUSHOLDER TRANSFORMATION VECTOR.
C
         SIGMA     = SQRT( SUM )
         BSQ       = SQRT( SUM + SIGMA * ABS( A( K, K ) ) )
         WORK( K ) = SIGN( SIGMA + ABS( A( K, K ) ),
     *                       A( K, K ) ) / BSQ
         A( K, K ) = - SIGN( SIGMA, A( K, K ) )
         IF ( KP1 .LE. NELVAR ) THEN
            DO 160 J     = KP1, NELVAR
               A( K, J ) = A( K, J ) / BSQ
  160       CONTINUE
C
C  APPLY THE TRANSFORMATION TO THE REMAINING ROWS OF A.
C
            DO 190 I    = KP1, NINVAR
               SUM      = WORK( K ) * A( I, K )
               DO 170 J = KP1, NELVAR
                  SUM   = SUM + A( K, J ) * A( I, J )
  170          CONTINUE
               A( I, K )    = A( I, K ) - SUM * WORK( K )
               DO 180 J     = KP1, NELVAR
                  A( I, J ) = A( I, J ) - SUM * A( K, J )
  180          CONTINUE
  190       CONTINUE
         END IF
  200 CONTINUE
C
C  THE REDUCTION OF A IS COMPLETE. BUILD THE GENERALIZED INVERSE.
C  FIRSTLY, APPLY THE FIRST ELEMENTARY TRANSFORMATION.
C
      SUM      = WORK( NINVAR ) / A( NINVAR, NINVAR )
      DO 210 J = NINVAR + 1, NELVAR
         A( NINVAR, J ) = - SUM * A( NINVAR, J )
  210 CONTINUE
      A( NINVAR, NINVAR ) = ONE / A( NINVAR, NINVAR )
     *                      - SUM * WORK( NINVAR )
C
C  NOW APPLY THE REMAINING NINVAR - 1 TRANSFORMATIONS.
C
      DO 300 K = NINVAR - 1, 1, - 1
         KP1   = K + 1
C
C  FIRST TRANSFORM THE LAST NINVAR - K ROWS.
C
         DO 240 I = KP1, NINVAR
            SUM   = ZERO
            DO 220 J = KP1, NELVAR
               SUM  = SUM + A( K, J ) * A( I, J )
  220       CONTINUE
            DO 230 J      = KP1, NELVAR
               A( I, J ) = A( I, J ) - SUM * A( K, J )
  230       CONTINUE
            WORK( I ) = - SUM * WORK( K )
  240    CONTINUE
C
C  THEN CALCULATE THE NEW K-TH ROW.
C
         DO 260 J    = KP1, NELVAR
            SUM      = - WORK( K ) * A( K, J )
            DO 250 I = KP1, NINVAR
               SUM   = SUM - A( I, K ) * A( I, J )
  250       CONTINUE
            A( K, J ) = SUM / A( K, K )
  260    CONTINUE
C
C  UPDATE THE K-TH COLUMN.
C
         SUM          = ONE - WORK( K ) ** 2
         DO 270 I     = KP1, NINVAR
            SUM       = SUM - A( I, K ) * WORK( I )
            A( I, K ) = WORK( I )
  270    CONTINUE
         A( K, K ) = SUM / A( K, K )
  300 CONTINUE   
C
C  UNDO THE ROW INTERCHANGES.
C
      DO 430 I = 1, NINVAR
  410    CONTINUE
         IR = IWORK( I )
         IF ( I .LT. IR ) THEN
            IWORK( I )  = IWORK( IR )
            IWORK( IR ) = IR
C
C  SWAP ROWS I AND IR.
C
            DO 420 J      = 1, NELVAR
               SUM        = A( I, J )
               A( I, J )  = A( IR, J )
               A( IR, J ) = SUM
  420       CONTINUE
            GO TO 410
         END IF
  430 CONTINUE
C
C  UNDO THE COLUMN INTERCHANGES.
C
      DO 460 J = 1, NELVAR
  440    CONTINUE
         I  = NINVAR + J
         IR = IWORK( I )
         IF ( J .LT. IR ) THEN
            K       = NINVAR + IR
            IWORK( I ) = IWORK( K )
            IWORK( K ) = IR
C
C  SWAP COLUMNS J AND IR.
C
            DO 450 I      = 1, NINVAR
               SUM        = A( I, J )
               A( I, J )  = A( I, IR )
               A( I, IR ) = SUM
  450       CONTINUE
            GO TO 440
         END IF
  460 CONTINUE
      RETURN
C
C  NON-EXECUTABLE STATEMENTS.
C
 2000 FORMAT( /, ' *** Error message from GENINV *** ', I3, 
     *       ' reduced rows found to be zero ' )
C
C  END OF GNINV.
C
      END
