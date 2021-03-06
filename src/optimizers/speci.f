C ** Correction report.
C ** Correction 1. 19/08/96: 2 lines replaced **
C ** Correction 2. 19/08/96: 4 lines replaced by 3 **
C ** Correction 3. 19/08/96: 1 line replaced **
C ** Correction 4. 19/08/96: 1 line replaced **
C ** Correction 5. 20/08/96: 1 line replaced **
C ** End of Correction report.
C  THIS VERSION: 20/08/1996 AT 17:27:00 AM.
CS    SUBROUTINE SSPECI( IISPEC, IPRNT,  IPSTRT, IPSTOP, IPGAP,  MAXIT, 
CD    SUBROUTINE DSPECI( IISPEC, IPRNT,  IPSTRT, IPSTOP, IPGAP,  MAXIT, 
     *                   NSEMIB, RMU,    RMUTOL, GETSCA, DECHKE, DECHKG, 
     *                   FATALE, FATALG, ICHOSE, IPRNTS, SCALEG, SCALEV,
     *                   FINDMX, STOPGA, STOPCA, FIRSTC, FIRSTG, TESTAL,
     *                   WARMST, RADIUS, ISTORE, ITEST , IOUT )
      INTEGER            IISPEC, IPRNT,  IPSTRT, IPSTOP, IPGAP,  MAXIT 
      INTEGER            IPRNTS, ITEST,  ISTORE, IOUT  , ICHOSE( 6 ) 
CS    REAL               RMU,    RMUTOL, STOPGA, STOPCA, FIRSTC, FIRSTG
CD    DOUBLE PRECISION   RMU,    RMUTOL, STOPGA, STOPCA, FIRSTC, FIRSTG
CS    REAL               RADIUS, FINDMX
CD    DOUBLE PRECISION   RADIUS, FINDMX
      LOGICAL            GETSCA, DECHKE, DECHKG, FATALE, FATALG
      LOGICAL            SCALEG, SCALEV, TESTAL, WARMST
C
C  READ THE ALGORITHM SPECIFICATION FILE AND SET UP APPROPRIATE
C  INPUT PARAMETERS FOR LANCELOT.
C
C  NICK GOULD 25/05/90
C  FOR CGT PRODUCTIONS.
C
      INTEGER            LDICT,  ISTART, ISTOP,  LSTART, LSTOP,  IVALUE
      INTEGER            ITYPE,  IFIELD, NSEMIB
CS    REAL               VALUE
CD    DOUBLE PRECISION   VALUE
      PARAMETER        ( LDICT = 54,     LSTART = 1,     LSTOP = 80 )
      CHARACTER * 20     FIELD
      CHARACTER * 80     LISPEC
      CHARACTER * 8      FIRST8, DICT( LDICT )
CS    EXTERNAL           SGETVL, SGETIN
CD    EXTERNAL           DGETVL, DGETIN
C
C  PROBLEM DATA.
C 
      DATA DICT(  1 ) / 'BEGIN   ' /  
      DATA DICT(  2 ) / 'BEGIN-SP' / 
      DATA DICT(  3 ) / 'BFGS-APP' /
      DATA DICT(  4 ) / 'GILL-MUR' /
      DATA DICT(  5 ) / 'CG-METHO' /
      DATA DICT(  6 ) / 'SCHNABEL' /
      DATA DICT(  7 ) / 'CHECK-EL' /
      DATA DICT(  8 ) / 'CHECK-GR' /
      DATA DICT(  9 ) / 'DFP-APPR' /
      DATA DICT( 10 ) / 'MAXIMIZE' /
      DATA DICT( 11 ) / 'TWO-NORM' /
      DATA DICT( 12 ) / 'DIAGONAL' /
      DATA DICT( 13 ) / 'END     ' /
      DATA DICT( 14 ) / 'END-SPEC' /
      DATA DICT( 15 ) / 'EXACT-CA' /
      DATA DICT( 16 ) / 'EXACT-SE' /
      DATA DICT( 17 ) / 'EXPANDIN' /
      DATA DICT( 18 ) / 'FULL-MAT' /
      DATA DICT( 19 ) / 'IGNORE-E' /
      DATA DICT( 20 ) / 'INITIAL-' /
      DATA DICT( 21 ) / 'DIRECT-M' /
      DATA DICT( 22 ) / 'MAXIMUM-' /
      DATA DICT( 23 ) / 'MUNKSGAA' /
      DATA DICT( 24 ) / 'PRINT-LE' /
      DATA DICT( 25 ) / 'START-PR' /
      DATA DICT( 26 ) / 'STOP-PRI' /
      DATA DICT( 27 ) / 'PSB-APPR' /
      DATA DICT( 28 ) / 'BANDSOLV' /  
      DATA DICT( 29 ) / 'SR1-APPR' /
      DATA DICT( 30 ) / 'CHECK-AL' /
      DATA DICT( 31 ) / 'MULTIFRO' /
      DATA DICT( 32 ) / 'USERS-PR' /
      DATA DICT( 33 ) / 'IGNORE-G' /
      DATA DICT( 34 ) / 'INFINITY' /
      DATA DICT( 35 ) / 'CHECK-DE' /
      DATA DICT( 36 ) / 'IGNORE-D' /
      DATA DICT( 37 ) / 'INEXACT-' /
      DATA DICT( 38 ) / 'GET-SCAL' /
      DATA DICT( 39 ) / 'USE-CONS' /
      DATA DICT( 40 ) / 'USE-SCAL' /
      DATA DICT( 41 ) / 'USE-VARI' /
      DATA DICT( 42 ) / 'PRINT-SC' /
      DATA DICT( 43 ) / 'GRADIENT' /
      DATA DICT( 44 ) / 'CONSTRAI' /
      DATA DICT( 45 ) / 'DECREASE' /
      DATA DICT( 46 ) / 'FIRST-CO' /
      DATA DICT( 47 ) / 'FIRST-GR' /
      DATA DICT( 48 ) / 'TRUST-RE' /
      DATA DICT( 49 ) / 'ITERATIO' /
      DATA DICT( 50 ) / 'MODIFIED' /
      DATA DICT( 51 ) / 'SOLVE-BQ' /
      DATA DICT( 52 ) / 'RESTART-' /
      DATA DICT( 53 ) / 'SAVE-DAT' /
      DATA DICT( 54 ) / 'FINITE-D' /
C
C  SET DEFAULT VALUES FOR LANCELOT PARAMETERS.
C
CS    FINDMX = 1.0E+0
CD    FINDMX = 1.0D+0
      IPRNT  = 1
      IPGAP  = 1
      IPRNTS = 0
      IPSTRT = 0
      IPSTOP = 100
      MAXIT  = 100
      ISTORE = 0
      NSEMIB = 5
      DECHKE = .FALSE.
      DECHKG = .FALSE.
      FATALE = .TRUE.
      FATALG = .TRUE.
      GETSCA = .FALSE.
      SCALEV = .FALSE.
      SCALEG = .FALSE.
      TESTAL = .FALSE.
CS    STOPCA = 1.0E-3
CD    STOPCA = 1.0D-5
CS    STOPGA = 1.0E-3
CD    STOPGA = 1.0D-5
CS    RMU    = 1.0E-1
CD    RMU    = 1.0D-1
      RMUTOL = RMU
      FIRSTC = RMU
      FIRSTG = RMU
      WARMST = .FALSE.
CS    RADIUS = - 1.0E+0
CD    RADIUS = - 1.0D+0
C
C  ICHOSE( 1 ) = 1 SPECIFIES 2-NORM TRUST-REGION, ANYTHING ELSE IS
C                FOR THE INFINITY NORM.
C  ICHOSE( 2 ) = IMETH GIVES THE METHOD TO BE USED FOR SOLVING THE
C                LINEAR SYSTEM. 1=CG, 2=DIAGONAL PRECONDITIONED CG,
C                3=USER-PROVIDED PRECONDITIONED CG, 4=EXPANDING BAND
C                PRECONDITIONED CG, 5=MUNKSGAARD'S PRECONDITIONED CG,
C                6=SCHNABEL-ESKOW MODIFIED CHOLESKY PRECONDITIONED CG, 
C                7=GILL-MURRAY-PONCELEON-SAUNDERS MODIFIED CHOLESKY
C                PRECONDITIONED CG, 8=BAND MATRIX
c                PRECONDITIONED CG, 11=MULTIFRONTAL DIRECT
C                METHOD, 12=DIRECT MODIFIED MULTIFRONTAL METHOD.
C  ICHOSE( 3 ) = 0 IF EXACT FIRST DERIVATIVES ARE GIVEN AND = 1 IF
C                FINITE DIFFERENCE APPROXIMATIONS ARE TO BE CALCULATED.
C  ICHOSE( 4 ) = IUPDAT, THE APPROXIMATION TO THE SECOND DERIVATIVES
C                USED. 0=EXACT, 1=BFGS, 2=DFP, 3=PSB, 4=SR1.
C  ICHOSE( 5 ) = 1 IF THE EXACT CAUCHY POINT IS REQUIRED, =2 IF AN
C                APPROXIMATION SUFFICES.
C  ICHOSE( 6 ) = 1 IF THE THE MINIMIZER OF THE QUADRATIC MODEL WITHIN
C                THE INTERSECTION OF THE TRUST-REGION AND FEASIBLE BOX
C                IS TO BE SOUGHT (TO A PRESCRIBED ACCURACY), =2 IF AN
C                APPROXIMATION SUFFICES.
C
      ICHOSE( 1 ) = 2
      ICHOSE( 2 ) = 8
      ICHOSE( 3 ) = 0
      ICHOSE( 4 ) = 4
      ICHOSE( 5 ) = 1
      ICHOSE( 6 ) = 2
C
C  READ THE NEXT LINE.
C  CHECK THAT THE END OF THE FILE HAS NOT BEEN REACHED. IF IT HAS,
C  NO CALL TO LANCELOT WILL BE MADE. OTHERWISE, LOOK FOR THE POSITION
C  OF THE FIRST NON-BLANK CHARACTER, ISTART, AND THE POSITION OF THE
C  END OF THE FOLLOWING STRING, ISTOP. IF THE FIRST CHARACTER  IS A
C  '*', THE CARD IS A COMMENT AND IS IGNORED.
C
    1 CONTINUE
      READ( UNIT = IISPEC, FMT = 1000, END = 600, ERR = 600 ) 
     *      LISPEC( LSTART: LSTOP )
      IF ( LISPEC( 1: 1 ) .EQ. '*' ) GO TO 1 
      DO 2 ISTART = LSTART, LSTOP
         IF ( LISPEC( ISTART: ISTART ) .NE. ' ' ) GO TO 3
    2 CONTINUE
      GO TO 1      
    3 CONTINUE
      DO 4 ISTOP = ISTART, LSTOP - 1
         IF ( LISPEC( ISTOP + 1: ISTOP + 1 ) .EQ. ' ' ) GO TO 5
    4 CONTINUE
      ISTOP = LSTOP
    5 CONTINUE
C
C  OBTAIN THE FIRST 8 CHARACTERS FROM THE STRING, FIRST8.
C
      FIRST8 =  '       '
      FIRST8( 1: MIN( 8, ISTOP - ISTART + 1 ) ) = 
     *   LISPEC( ISTART: ISTART + MIN( 7, ISTOP - ISTART ) )
C
C  CONVERT LOWER CASE CHARACTERS IN FIRTS8 TO UPPER CASE.
C
C ** Correction 1. 19/08/96: 2 lines replaced **
CS    CALL SCNVRT( FIRST8 )
CD    CALL DCNVRT( FIRST8 )
C ** Correction 1. 19/08/96: end of correction **
C
C  LOOK UP THE FIRST 8 CHARACTERS FROM THE STRING IN THE DICTIONARY.
C
      DO 6 ITYPE = 1, LDICT
         IF ( FIRST8 .EQ. DICT( ITYPE ) ) GO TO 7
    6 CONTINUE
      WRITE( IOUT, 2000 ) LISPEC( LSTART: LSTOP )
      GO TO 1
    7 CONTINUE
C
C  READ ADDITIONAL NUMERICAL VALUES FROM CERTAIN CARDS. FIRSTLY, 
C  DETERMINE THE POSITION OF THE START OF THE NUMERICAL VALUE.
C  THEN TRANSFER THE VALUE INTO THE CHARACTER STRING FIELD AND 
C  FINALLY CONVERT FIELD INTO A NUMERICAL VALUE. 
C 
      IF ( ITYPE .EQ. 20 .OR. ITYPE. EQ. 22 .OR. ITYPE .EQ. 24 .OR.
     *     ITYPE .EQ. 25 .OR. ITYPE. EQ. 26 .OR. ITYPE .EQ. 28 .OR.
     *     ITYPE .EQ. 43 .OR. ITYPE. EQ. 44 .OR. ITYPE .EQ. 45 .OR.
     *     ITYPE .EQ. 46 .OR. ITYPE .EQ. 47 .OR. ITYPE .EQ. 48 .OR.
     *     ITYPE .EQ. 49 .OR. ITYPE .EQ. 53 ) THEN
         DO 8 IFIELD = ISTOP + 1, LSTOP
            IF ( LISPEC( IFIELD: IFIELD ) .NE. ' ' ) GO TO 9
    8    CONTINUE
    9    CONTINUE
         FIELD =  '                    '
         FIELD( 1: MIN( 20, LSTOP - IFIELD + 1 ) ) = 
     *      LISPEC( IFIELD: IFIELD + MIN( 19, LSTOP - IFIELD ) )
         IF ( ITYPE .EQ. 20 .OR. ITYPE. EQ. 43 .OR. ITYPE .EQ. 44 .OR.
     *        ITYPE .EQ. 45 .OR. ITYPE. EQ. 46 .OR. 
     *        ITYPE .EQ. 47 .OR. ITYPE .EQ. 48 ) THEN
CS          CALL SGETVL( FIELD, VALUE )
CD          CALL DGETVL( FIELD, VALUE )
         ELSE
CS          CALL SGETIN( FIELD, IVALUE )
CD          CALL DGETIN( FIELD, IVALUE )
         END IF
      END IF
C
C  BRANCH ACCORDING TO THE STRING FOUND. IGNORE ENTRIES FOR WHICH
C  NO MATCH IS FOUND. 
C
      GO TO ( 10, 10, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130,
     *        130, 150, 160, 170, 40, 190, 200, 210, 220, 230, 240,
     *        250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350,
     *        360, 370, 380, 390, 400, 410, 420, 430, 440, 450, 460,
     *        470, 480, 490,  60, 510, 520, 530, 540 ), ITYPE
C
C   LINE CONTAINS STRING 'BEGIN   ' OR 'BEGIN-SP'.
C
   10 CONTINUE
      GO TO 1 
C
C   LINE CONTAINS STRING 'BFGS-APP'.
C
   30 CONTINUE
      ICHOSE( 4 ) = 1
      GO TO 1 
C
C   LINE CONTAINS STRING 'FULL-MAT' OR 'GILL-MUR'.
C
   40 CONTINUE
      ICHOSE( 2 ) = 7
      GO TO 1 
C
C   LINE CONTAINS STRING 'CG-METHO'.
C
   50 CONTINUE
      ICHOSE( 2 ) = 1
      GO TO 1 
C
C   LINE CONTAINS STRING 'MODIFIED' OR 'SCHNABEL'.
C
   60 CONTINUE
      ICHOSE( 2 ) = 6
      GO TO 1
C
C   LINE CONTAINS STRING 'CHECK-EL'.
C
   70 CONTINUE
      DECHKE = .TRUE.
      GO TO 1 
C
C   LINE CONTAINS STRING 'CHECK-GR'.
C
   80 CONTINUE
      DECHKG = .TRUE.
      GO TO 1 
C
C   LINE CONTAINS STRING 'DFP-APPR'.
C
   90 CONTINUE
      ICHOSE( 4 ) = 2
      GO TO 1 
C
C   LINE CONTAINS STRING 'MAXIMIZE'.
C
  100 CONTINUE
CS    FINDMX = - 1.0E+0
CD    FINDMX = - 1.0D+0
      GO TO 1 
C
C   LINE CONTAINS STRING 'TWO-NORM'.
C
  110 CONTINUE
      ICHOSE( 1 ) = 1
      GO TO 1 
C
C   LINE CONTAINS STRING 'DIAGONAL'.
C
  120 CONTINUE
      ICHOSE( 2 ) = 2
      GO TO 1 
C
C   LINE CONTAINS STRING 'END     ' OR 'END-SPEC'.
C
  130 CONTINUE
      ITEST = 1
      RETURN
C
C   LINE CONTAINS STRING 'EXACT-CA'.
C
  150 CONTINUE
      ICHOSE( 5 ) = 1
      GO TO 1 
C
C   LINE CONTAINS STRING 'EXACT-SE'.
C
  160 CONTINUE
      IF ( ICHOSE( 3 ) .EQ. 0 ) THEN
         ICHOSE( 4 ) = 0
      ELSE
         WRITE( IOUT, 2010 )
      END IF
      GO TO 1 
C
C   LINE CONTAINS STRING 'EXPANDIN'.
C
  170 CONTINUE
      ICHOSE( 2 ) = 4
      GO TO 1 
C
C   LINE CONTAINS STRING 'IGNORE-E'.
C
  190 CONTINUE
      FATALE = .FALSE.
      GO TO 1 
C
C   LINE CONTAINS STRING 'INITIAL-'.
C
  200 CONTINUE
      RMU = VALUE
      GO TO 1 
C
C   LINE CONTAINS STRING 'DIRECT-M'.
C
  210 CONTINUE
      ICHOSE( 2 ) = 12
      GO TO 1 
C
C   LINE CONTAINS STRING 'MAXIMUM-'.
C
  220 CONTINUE
      MAXIT = IVALUE
      GO TO 1 
C
C   LINE CONTAINS STRING 'MUNKSGAA'.
C
  230 CONTINUE
      ICHOSE( 2 ) = 5
      GO TO 1 
C
C   LINE CONTAINS STRING 'PRINT-LE'.
C
  240 CONTINUE
      IPRNT = IVALUE
      GO TO 1 
C
C   LINE CONTAINS STRING 'START-PR'.
C
  250 CONTINUE
      IPSTRT = IVALUE
      GO TO 1 
C
C   LINE CONTAINS STRING 'STOP-PRI'.
C
  260 CONTINUE
      IPSTOP = IVALUE
      GO TO 1 
C
C   LINE CONTAINS STRING 'PSB-APPR'.
C
  270 CONTINUE
      ICHOSE( 4 ) = 3
      GO TO 1 
C
C   LINE CONTAINS STRING 'BANDSOLV'.
C
  280 CONTINUE
      ICHOSE( 2 ) = 8
      NSEMIB = IVALUE
      GO TO 1 
C
C   LINE CONTAINS STRING 'SR1-APPR'.
C
  290 CONTINUE
      ICHOSE( 4 ) = 4
      GO TO 1 
C
C   LINE CONTAINS STRING 'CHECK-AL'.
C
  300 CONTINUE
      DECHKE = .TRUE.
      DECHKG = .TRUE.
      TESTAL = .TRUE.
      GO TO 1 
C
C   LINE CONTAINS STRING 'MULTIFRO'.
C
  310 CONTINUE
      ICHOSE( 2 ) = 11
      GO TO 1 
C
C   LINE CONTAINS STRING 'USERS-PR'.
C
  320 CONTINUE
      ICHOSE( 2 ) = 3
      GO TO 1 
C
C   LINE CONTAINS STRING 'IGNORE-G'.
C
  330 CONTINUE
      FATALG = .FALSE.
      GO TO 1 
C
C   LINE CONTAINS STRING 'INFINITY'.
C
  340 CONTINUE
      ICHOSE( 1 ) = 2
      GO TO 1 
C
C   LINE CONTAINS STRING 'CHECK-DE'.
C
  350 CONTINUE
      DECHKE = .TRUE.
      DECHKG = .TRUE.
      GO TO 1 
C
C   LINE CONTAINS STRING 'IGNORE-D'.
C
  360 CONTINUE
      FATALE = .FALSE.
      FATALG = .FALSE.
      GO TO 1 
C
C   LINE CONTAINS STRING 'INEXACT-'.
C
  370 CONTINUE
      ICHOSE( 5 ) = 2
      GO TO 1 
C
C   LINE CONTAINS STRING 'GET-SCAL'.
C
  380 CONTINUE
      GETSCA = .TRUE.
      IPRNTS = 1
      GO TO 1 
C
C   LINE CONTAINS STRING 'USE-CONS'.
C
  390 CONTINUE
      GETSCA = .TRUE.
      SCALEG = .TRUE.
      GO TO 1 
C
C   LINE CONTAINS STRING 'USE-SCAL'.
C
  400 CONTINUE
      GETSCA = .TRUE.
      SCALEV = .TRUE.
      SCALEG = .TRUE.
      GO TO 1 
C
C   LINE CONTAINS STRING 'USE-VARI'.
C
  410 CONTINUE
      GETSCA = .TRUE.
      SCALEV = .TRUE.
      GO TO 1 
C
C   LINE CONTAINS STRING 'PRINT-SC'.
C
  420 CONTINUE
      IPRNTS = 1
      GO TO 1 
C
C   LINE CONTAINS STRING 'GRADIENT'.
C
  430 CONTINUE
      STOPGA = VALUE
      GO TO 1 
C
C   LINE CONTAINS STRING 'CONSTRAI'.
C
  440 CONTINUE
      STOPCA = VALUE
      GO TO 1 
C
C   LINE CONTAINS STRING 'DECREASE'.
C
  450 CONTINUE
      RMUTOL = VALUE
      GO TO 1 
C
C   LINE CONTAINS STRING 'FIRST-CO'.
C
  460 CONTINUE
      FIRSTC = VALUE
      GO TO 1 
C
C   LINE CONTAINS STRING 'FIRST-GR'.
C
  470 CONTINUE
      FIRSTG = VALUE
      GO TO 1 
C
C   LINE CONTAINS STRING 'TRUST-RE'.
C
  480 CONTINUE
      RADIUS = VALUE
      GO TO 1 
C
C   LINE CONTAINS STRING 'ITERATIO'.
C
  490 CONTINUE
      IPGAP = IVALUE 
      GO TO 1 
C
C   LINE CONTAINS STRING 'SOLVE-BQ'.
C
  510 CONTINUE
      ICHOSE( 6 ) = 1
      GO TO 1 
C
C   LINE CONTAINS STRING 'RESTART-'.
C
  520 CONTINUE
      WARMST = .TRUE.
      GO TO 1 
C
C   LINE CONTAINS STRING 'SAVE-DAT'.
C
  530 CONTINUE
      ISTORE = IVALUE
      GO TO 1 
C
C   LINE CONTAINS STRING 'FINITE-D'.
C
  540 CONTINUE
      ICHOSE( 3 ) = 1
      IF ( ICHOSE( 4 ) .EQ. 0 ) THEN
         WRITE( IOUT, 2010 )
         ICHOSE( 4 ) = 4
      END IF
C
C  BRANCH BACK TO READ THE NEXT INPUT LINE
C
      GO TO 1 
  600 CONTINUE
      ITEST = 0
      RETURN
C
C  NON-EXECUTABLE STATEMENTS.
C
 1000 FORMAT( A80 )
 2000 FORMAT( ' Line not recognised in specification file. Line reads:',
     *        /, A80 )
 2010 FORMAT( ' ** Warning. As finite-difference gradients have been',
     *        ' specified, exact second ', /,
     *        '    derivatives may not be used. A secant update will',
     *        ' be substituted.' )
C
C  END OF SUBROUTINE SPECI.
C
      END
C  THIS VERSION: 15/11/1991 AT 09:47:14 AM.
CS    SUBROUTINE SGETIN( FIELD, IVALUE )
CD    SUBROUTINE DGETIN( FIELD, IVALUE )
      INTEGER        IVALUE
      CHARACTER * 20 FIELD
C
C  READ THE INTEGER NUMBER IVALUE STORED IN THE CHARACTER FIELD.
C
C  NICK GOULD 02/08/89
C  FOR CGT PRODUCTIONS.
C
      INTEGER        I, J
      CHARACTER * 20 FIELD2
C
C  RIGHT-SHIFT THE FIELD, ELIMINATING BLANKS.
C
      FIELD2  = '            '
            J = 20
      DO 10 I = 20, 1, - 1
         IF ( FIELD( I : I ) .EQ. ' ' ) GO TO 10
         FIELD2( J : J ) = FIELD( I : I )
         J = J - 1
   10 CONTINUE
      READ( UNIT = FIELD2, FMT = 2000 ) IVALUE
      RETURN
C
C  NON-EXECUTABLE STATEMENTS.
C
 2000 FORMAT( I20 )
C
C  END OF SUBROUTINE GETIN.
C
      END
C  THIS VERSION: 15/11/1991 AT 09:47:14 AM.
CS    SUBROUTINE SGETVL( FIELD, VALUE )
CD    SUBROUTINE DGETVL( FIELD, VALUE )
CS    REAL             VALUE
CD    DOUBLE PRECISION VALUE
      CHARACTER * 20 FIELD
C
C  READ THE REAL NUMBER VALUE STORED IN THE CHARACTER FIELD.
C
C  NICK GOULD 02/08/89
C  FOR CGT PRODUCTIONS.
C
      READ( UNIT = FIELD, FMT = 2000 ) VALUE
      RETURN
C
C  NON-EXECUTABLE STATEMENTS.
C
 2000 FORMAT( BN, F20.0 )
      END
C  THIS VERSION: 15/11/1991 AT 09:47:14 AM.
C ** Correction 2. 19/08/96: 4 lines replaced by 3 **
CS    SUBROUTINE SCNVRT( STRING )
CD    SUBROUTINE DCNVRT( STRING )
      CHARACTER * 8 STRING
C ** Correction 2. 19/08/96: end of correction **
C
C  CONVERT ANY LOWER CASE CHARACTERS IN THE CHARACTER ARRAY STRING TO
C  UPPER CASE. THIS IS NOT VERY EFFICIENT AND MAYBE SHOULD BE REPLACED
C  BY A HASHING ROUTINE.
C
C  NICK GOULD 25/05/90
C  FOR CGT PRODUCTIONS.
C
      INTEGER       I, LETTER
      CHARACTER * 1 LOWER( 26 ), UPPER( 26 )
      DATA LOWER / 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 
     *             'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 
     *             'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z' /
      DATA UPPER / 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 
     *             'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 
     *             'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z' /
C
C  LOOP OVER EACH CHARACTER IN THE STRING.
C
C ** Correction 5. 20/08/96: 1 line replaced **
      DO 100 I = 1, 8
C ** Correction 5. 20/08/96: end of correction **
C
C  SEE IF THE CURRENT LETTER IS LOWER CASE. IF SO REPLACE IT BY ITS 
C  UPPER CASE COUNTERPART.
C
         DO 10 LETTER = 1, 26
C ** Correction 3. 19/08/96: 1 line replaced **
            IF ( STRING( I : I ) .EQ. LOWER( LETTER ) ) GO TO 20
C ** Correction 3. 19/08/96: end of correction **
   10    CONTINUE
         GO TO 100
   20    CONTINUE
C ** Correction 4. 19/08/96: 1 line replaced **
         STRING( I : I ) = UPPER( LETTER )
C ** Correction 4. 19/08/96: end of correction **
  100 CONTINUE
      RETURN
C
C  END OF SUBROUTINE CNVRT.
C
      END
