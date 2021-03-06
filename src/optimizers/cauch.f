C ** Correction report.
C ** Correction 1. 13/01/94: 2 lines corrected **
C ** Correction 2. 18/06/96: 1 line replaced by 20 **
C ** End of Correction report.
C  THIS VERSION: 13/01/1994 AT 04:37:01 PM.
C  VERSION OF 20/4/89 INCLUDING HEAP SORT.
C
C  ** FOR THE CRAY 2, LINES STARTING CDIR$ IVDEP TELL THE COMPILER TO
C     IGNORE POTENTIAL VECTOR DEPENDENCIES AS THEY ARE KNOWN TO BE O.K.
C
CS    SUBROUTINE SCAUCH( N, X0, XT, G, BND, INDEX, F, EPSTOL, BOUNDX,
CD    SUBROUTINE DCAUCH( N, X0, XT, G, BND, INDEX, F, EPSTOL, BOUNDX,
     *                   DXSQR, FXT, P, Q, IVAR, NFREE,
     *                   NVAR1, NVAR2, NNONNZ, INONNZ, BREAKP,
     *                   IOUT, JUMPTO, IDEBUG )
C
C  *******************************************************************
C
C  FIND THE GENERALIZED CAUCHY POINT IN THE DIRECTION P FROM X0
C  FOR A GIVEN QUADRATIC FUNCTION WITHIN A BOX SHAPED REGION.
C
C  IF WE DEFINE THE 'CAUCHY ARC' X(T) = PROJECTION OF X0 + T*P INTO
C  THE BOX REGION BND(*,1) <= X(*) <= BND(*,2), THE GENERALIZED CAUCHY
C  POINT IS THE FIRST LOCAL MINIMIZER OF THE QUADRATIC FUNCTION
C
C    0.5 (X - X0) (TRANSPOSE ) B (X - X0) + G (TRANSPOSE) (X - X0) + F
C
C  FOR POINTS LYING ON X(T), WITH T >= 0. OPTIONALLY, THE SEARCH
C  FOR THE GENERALIZED CAUCHY POINT MAY BE TERMINATED AT THE FIRST
C  POINT ENCOUNTERED ON THE BOUNDARY OF THE SPHERICAL REGION
C  ||X - X0|| <= R, WHERE ||.|| DENOTES THE 2-NORM. CONTROL IS PASSED
C  FROM THE ROUTINE WHENEVER A PRODUCT OF THE VECTOR P WITH B IS
C  REQUIRED, AND THE USER IS RESPONSIBLE FOR FORMING THE PRODUCT
C  IN THE VECTOR Q.
C
C  PROGRESS THROUGH THE ROUTINE IS CONTROLLED BY THE PARAMETER
C  JUMPTO.
C
C  IF JUMPTO = 0, THE GENERALIZED CAUCHY POINT HAS BEEN FOUND.
C  IF JUMPTO = 1, AN INITIAL ENTRY HAS BEEN MADE.
C  IF JUMPTO = 2, 3, 4 THE VECTOR Q = B * P IS REQUIRED.
C
C  THE STATUS OF THE ARRAY INDEX GIVES THE STATUS OF THE VARIABLES.
C
C  IF INDEX( I ) = 0, THE I-TH VARIABLE IS FREE.
C  IF INDEX( I ) = 1, THE I-TH VARIABLE IS ON ITS LOWER BOUND.
C  IF INDEX( I ) = 2, THE I-TH VARIABLE IS ON ITS UPPER BOUND.
C  IF INDEX( I ) = 3, THE I-TH VARIABLE IS FIXED.
C  IF INDEX( I ) = 4, THE I-TH VARIABLE IS TEMPORARILY FIXED.
C
C  THE ADDRESSES OF THE FREE VARIABLES ARE GIVEN IN THE FIRST
C  NVAR2 ENTRIES OF THE ARRAY IVAR.
C  IF THE PRODUCT OF B * P IS REQUIRED (JUMPTO = 2, 3, 4), THE
C  NONZERO COMPONENTS OF P OCCUR IN POSITIONS IVAR(I) FOR
C  I LYING BETWEEN NVAR1 AND NVAR2.
C
C  AT THE INITIAL POINT, VARIABLES WITHIN EPSTOL OF THEIR BOUNDS AND
C  FOR WHICH THE SEARCH DIRECTION POINTS OUT OF THE BOX WILL BE FIXED.
C
C  NICK GOULD, 20TH APRIL 1989.
C  FOR CGT PRODUCTIONS.
C
C  VERSIONS AVAILABLE :
C  --------------------
C
C  THERE IS A SINGLE AND A DOUBLE PRECISION VERSION. TO USE THE
C  SINGLE PRECISION VERSION, REPLACE THE CHARACTERS 'CS' WITH '  ' IN
C  THE FIRST TWO COLUMNS OF ALL LINES. TO USE THE DOUBLE PRECISION
C  VERSION, REPLACE THE CHARACTERS 'CD' WITH '  ' IN THE FIRST TWO
C  COLUMNS OF EVERY LINE.
C
C  PARAMETERS : (REAL MEANS DOUBLE PRECISION IN THE DOUBLE VERSION)
C  ------------
C
C  N      (INTEGER) THE NUMBER OF INDEPENDENT VARIABLES.
C          ** THIS VARIABLE IS NOT ALTERED BY THE SUBROUTINE.
C  X0     (REAL ARRAY OF LENGTH AT LEAST N) THE POINT X0 FROM
C          WHICH THE 'CAUCHY ARC' COMMENCES.
C          ** THIS VARIABLE IS NOT ALTERED BY THE SUBROUTINE.
C  XT     (REAL ARRAY OF LENGTH AT LEAST N) THE CURRENT ESTIMATE
C          OF THE GENERALIZED CAUCHY POINT.
C  G      (REAL ARRAY OF LENGTH AT LEAST N) THE COEFFICIENTS OF
C          THE LINEAR TERM IN THE QUADRATIC FUNCTION.
C          ** THIS VARIABLE IS NOT ALTERED BY THE SUBROUTINE.
C  BND    (TWO DIMENSIONAL REAL ARRAY WITH LEADING DIMENSION N AND
C          SECOND DIMENSION 2) THE LOWER (BND(*,1)) AND
C          UPPER (BND(*,2)) BOUNDS ON THE VARIABLES.
C          ** THIS VARIABLE IS NOT ALTERED BY THE SUBROUTINE.
C  INDEX  (INTEGER ARRAY OF LENGTH AT LEAST N) GIVES THE STATUS OF
C          THE VARIABLES AT THE CURRENT ESTIMATE OF THE GENERALIZED
C          CAUCHY POINT AS FOLLOWS:
C          IF INDEX( I ) = 0, THE I-TH VARIABLE IS FREE.
C          IF INDEX( I ) = 1, THE I-TH VARIABLE IS ON ITS LOWER BOUND.
C          IF INDEX( I ) = 2, THE I-TH VARIABLE IS ON ITS UPPER BOUND.
C          IF INDEX( I ) = 3, THE I-TH VARIABLE IS FIXED.
C          IF INDEX( I ) = 4, THE I-TH VARIABLE IS TEMPORARILY FIXED.
C  F      (REAL) THE VALUE OF THE QUADRATIC AT X0, SEE ABOVE.
C          ** THIS VARIABLE IS NOT ALTERED BY THE SUBROUTINE.
C  EPSTOL (REAL) A TOLERANCE ON FEASIBILITY OF X0, SEE ABOVE.
C          ** THIS VARIABLE IS NOT ALTERED BY THE SUBROUTINE.
C  BOUNDX (LOGICAL) THE SEARCH FOR THE GENERALIZED CAUCHY POINT
C          WILL BE TERMINATED ON THE BOUNDARY OF THE SPHERICAL
C          REGION ||X-X0|| <= R IF AND ONLY IF BOUNDX IS SET TO
C          .TRUE. ON INITIAL (JUMPTO=1) ENTRY.
C          ** THIS VARIABLE IS NOT ALTERED BY THE SUBROUTINE.
C  DXSQR  (REAL) THE SQUARE OF THE TWO NORM OF THE DISTANCE BETWEEN
C          THE CURRENT ESTIMATE OF THE GENERALIZED CAUCHY POINT AND
C          X0. DXSQR WILL ONLY BE SET IF BOUNDX IS .TRUE.
C  FXT    (REAL) THE VALUE OF THE PIECEWISE QUADRATIC FUNCTION
C          AT THE CURRENT ESTIMATE OF THE GENERALIZED CAUCHY POINT.
C  P      (REAL ARRAY OF LENGTH AT LEAST N) CONTAINS THE VALUES OF THE
C          COMPONENTS OF THE VECTOR P. ON INITIAL (JUMPTO=1) ENTRY,
C          P MUST CONTAIN THE INITIAL DIRECTION
C          OF THE 'CAUCHY ARC'. ON A NON OPTIMAL EXIT,
C          (JUMPTO=2,3,4), P IS THE VECTOR FOR WHICH THE PRODUCT B*P
C          IS REQUIRED BEFORE THE NEXT RE-ENTRY. ON A TERMINAL EXIT
C          (JUMPTO=0), P CONTAINS THE STEP XT - X0. THE COMPONENTS
C          IVAR(I) = NVAR1, ... , NVAR2 OF P CONTAIN THE VALUES OF THE
C          NONZERO COMPONENTS OF P (SEE, IVAR, NVAR1, NVAR2).
C  Q      (REAL ARRAY OF LENGTH AT LEAST N) ON A NON INITIAL ENTRY
C         (JUMPTO=2,3,4), Q MUST CONTAIN THE VECTOR B * P. ONLY THE
C          COMPONENTS IVAR(I), I=1,...,NFREE, OF Q NEED BE SET (THE
C          OTHER COMPONENTS ARE NOT USED).
C  IVAR   (INTEGER ARRAY OF LENGTH AT LEAST N) ON ALL NORMAL EXITS
C         (JUMPTO=0,2,3,4), IVAR(I), I=NVAR1,...,NVAR2, GIVES THE INDICES
C          OF THE NONZERO COMPONENTS OF P.
C  NFREE  (INTEGER) THE NUMBER OF FREE VARIABLES AT THE INITIAL POINT.
C  NVAR1  (INTEGER) SEE IVAR, ABOVE.
C  NVAR2  (INTEGER) SEE IVAR, ABOVE.
C  NNONNZ (INTEGER) THE NUMBER OF NONZERO COMPONENTS OF Q ON A
C         JUMPTO=3 ENTRY. NNONNZ NEED NOT BE SET ON OTHER ENTRIES.
C  INONNZ (INTEGER ARRAY OF LENGTH AT LEAST NNONNZ) ON JUMPTO = 3
C         ENTRIES, INONN(I), I = 1,....,NNONNZ, MUST GIVE THE
C         INDICES OF THE NONZERO COMPONENTS OF Q. ON OTHER ENTRIES,
C         INNONNZ NEED NOT BE SET.
C  BREAKP (REAL ARRAY OF LENGTH AT LEAST N) USED AS WORKSPACE.
C  IOUT   (INTEGER) THE FORTRAN OUTPUT CHANNEL NUMBER TO BE USED.
C  JUMPTO (INTEGER) CONTROLS FLOW THROUGH THE SUBROUTINE.
C          IF JUMPTO = 0, THE GENERALIZED CAUCHY POINT HAS BEEN FOUND.
C          IF JUMPTO = 1, AN INITIAL ENTRY HAS BEEN MADE.
C          IF JUMPTO = 2, 3, 4, THE VECTOR Q = B * P IS REQUIRED.
C  IDEBUG (INTEGER) ALLOWS DETAILED PRINTING. IF IDEBUG IS LARGER
C          THAN 4, DETAILED OUTPUT FROM THE ROUTINE WILL BE GIVEN.
C          OTHERWISE, NO OUTPUT OCCURS.
C
C  COMMON BLOCKS.
C  --------------
C
C  THE COMMON BLOCK SCOMSB/DCOMSB IS USED TO PASS FURTHER INFORMATION.
C  THE BLOCK CONTAINS THE VARIABLES
C
C    ACCCG, RATIO, RADIUS, RADMAX, FINDMX, CMA31, ITERCG, ITCGMX, 
C    NGEVAL, ISKIP, IFIXED, NBANDW
C
C  IN THAT ORDER.  VARIABLES ITERCG, ITCGMX, NGEVAL, ISKIP, IFIXED AND
C  NBANDW ARE INTEGER AND ACCCG, RATIO, RADIUS, RADMAX, FINDMX
C  AND CMA31 ARE REAL. ONLY RADIUS AND FINDMX NEED BE SET AS FOLLOWS:
C
C  RADIUS (REAL) THE RADIUS, R, OF THE SPHERICAL REGION. RADIUS
C          NEED NOT BE SET IF BOUNDX IS .FALSE. ON INITIAL ENTRY.
C          ** THIS VARIABLE IS NOT ALTERED BY THE SUBROUTINE.
C
C  FINDMX (REAL) WHEN PRINTING THE VALUE OF THE OBJECTIVE FUNCTION,
C          THE VALUE CALCULATED IN FMODEL WILL BE MULTIPLIED BY THE
C          SCALE FACTOR FINDMX. THIS ALLOWS THE USER, FOR INSTANCE,
C          TO FIND THE MAXIMUM OF A QUADRATIC FUNCTION F, BY MINIMIZING
C          THE FUNCTION - F, WHILE MONITORING THE PROGRESS AS
C          IF A MAXIMIZATION WERE ACTUALLY TAKING PLACE, BY SETTING
C          FINDMX TO - 1.0. NORMALLY FINDMX SHOULD BE SET TO 1.0.
C
C  THE COMMON BLOCK SMACHN/DMACHN IS USED TO PASS MACHINE DEPENDENT 
C  CONSTANTS. THE BLOCK CONTAINS THE REAL VARIABLES
C
C  EPSMCH, EPSNEG, TINY, BIG
C
C  IN THAT ORDER. EPSMCH AND EPSNEG ARE, RESPECTIVELY, THE SMALLEST
C  POSITIVE NUMBERS SUCH THAT 1.0 + EPSMCH > 1.0 AND 1.0 - EPSNEG < 1.0.
C  TINY IS THE SMALLEST POSITIVE FULL-PRECISION NUMBER AND BIG IS THE
C  LARGEST REPRESENTABLE NUMBER. ** THESE VALUES MUST BE SET BY THE USER.
C                                   ------------------------------------      
C
C  INITIALIZATION.
C  ---------------
C
C  ON THE INITIAL CALL TO THE SUBROUTINE THE FOLLOWING VARIABLES MUST
C  BE SET BY THE USER: -
C
C      N, X0, G, P, BND, F, EPSTOL, BOUNDX, IOUT, JUMPTO, IDEBUG.
C
C  JUMPTO MUST HAVE THE VALUE 1.
C  IN ADDITION, IF THE I-TH VARIABLE IS REQUIRED TO BE FIXED AT
C  ITS INITIAL VALUE, X0(I), INDEX(I) MUST BE SET TO 3.
C  RADIUS MUST BE SPECIFIED IF BOUNDX IS .TRUE. ON INITIAL ENTRY.
C
C  RE-ENTRY.
C  ---------
C
C  IF THE VARIABLE JUMPTO HAS THE VALUE 2, 3 OR 4 ON EXIT, THE
C  SUBROUTINE MUST BE RE-ENTERED WITH THE VECTOR Q CONTAINING
C  THE PRODUCT OF THE SECOND DERIVATIVE MATRIX B AND THE OUTPUT
C  VECTOR P. ALL OTHER PARAMETERS MUST NOT BE ALTERED.
C                                 -------------------
C
C  ****************************************************************
C
      INTEGER          JUMPTO, N, NFREE, NVAR1, NVAR2
      INTEGER          NNONNZ, IOUT, IDEBUG
CS    REAL             DXSQR, EPSTOL, F, FXT
CD    DOUBLE PRECISION DXSQR, EPSTOL, F, FXT
      LOGICAL          BOUNDX
      INTEGER          INDEX( N ), IVAR( N ), INONNZ( N )
CS    REAL             X0( N ), XT( N ), BND( N, 2 ), G( N ), P( N ),
CD    DOUBLE PRECISION X0( N ), XT( N ), BND( N, 2 ), G( N ), P( N ),
     *                 Q( N ), BREAKP( N )
C
C  LOCAL VARIABLES.
C
      INTEGER          I, IBREAK, ITER, J, INSORT, NBREAK, NFREED, NZERO
CS    REAL             T, TK, TSTAR, GXT, HXT,
CD    DOUBLE PRECISION T, TK, TSTAR, GXT, HXT,
     *                 TWO, HALF, ZERO, GZERO, HZERO, EPSTL2,
     *                 TPTTP, PTP, TBNDSQ, FEASEP, TEN, TCAUCH, QIPI,
     *                 TBREAK, DELTAT, GP, PBP, EPSQRT, GXTOLD
      LOGICAL          XLOWER, XUPPER, PRNTER, PRONEL, RECOMP
C
C  COMMON VARIABLES.
C
      INTEGER           ITERCG, ITCGMX, NGEVAL, ISKIP , IFIXED, NBANDW
CS    REAL              ACCCG , RATIO , RADIUS, RADMAX, FINDMX, CMA31
CD    DOUBLE PRECISION  ACCCG , RATIO , RADIUS, RADMAX, FINDMX, CMA31
CS    COMMON / SCOMSB / ACCCG , RATIO , RADIUS, RADMAX, FINDMX, CMA31,
CD    COMMON / DCOMSB / ACCCG , RATIO , RADIUS, RADMAX, FINDMX, CMA31,
     *                  ITERCG, ITCGMX, NGEVAL, ISKIP , IFIXED, NBANDW
CS    REAL              EPSMCH, EPSNEG, TINY  , BIG
CD    DOUBLE PRECISION  EPSMCH, EPSNEG, TINY  , BIG
CS    COMMON / SMACHN / EPSMCH, EPSNEG, TINY  , BIG
CD    COMMON / DMACHN / EPSMCH, EPSNEG, TINY  , BIG
C
C  MACHINE FUNCTIONS.
C
      INTRINSIC         ABS, MIN, SQRT
C
C  EXTERNAL FUNCTIONS AND SUBROUTINES.
C
CS    EXTERNAL          SSORTA, SSORTB
CD    EXTERNAL          DSORTA, DSORTB
C
C  THE FOLLOWING VARIABLES ARE SAVED BETWEEN CALLS OF THE ROUTINE.
C
      SAVE              ITER,  PRNTER, PRONEL, GXT,    HXT,    TK
      SAVE              PTP,   TPTTP,  NBREAK, DELTAT, TCAUCH
      SAVE              NZERO, EPSTL2, EPSQRT, GXTOLD, TBREAK
C ** Correction 1. 13/01/94: 2 lines corrected **
CS    SAVE            / SCOMSB /, / SMACHN /
CD    SAVE            / DCOMSB /, / DMACHN /
C ** Correction 1. 13/01/94:  end of correction **
C
C  SET CONSTANT REAL PARAMETERS.
C
CS    PARAMETER ( ZERO   = 0.0E+0, HALF   = 5.0E-1, TWO    = 2.0E+0 )
CD    PARAMETER ( ZERO   = 0.0D+0, HALF   = 5.0D-1, TWO    = 2.0E+0 )
CS    PARAMETER ( TEN    = 1.0E+1                                   )
CD    PARAMETER ( TEN    = 1.0D+1                                   )
C
C  GZERO IS THE SMALLEST VALUE THAT THE FIRST DERIVATIVE OF
C  THE PIECEWISE QUADRATIC MAY HAVE AND STILL BE CONSIDERED NONZERO.
C  - HZERO IS THE SMALLEST VALUE THAT THE SECOND DERIVATIVE OF
C  THE PIECEWISE QUADRATIC MAY HAVE AND STILL BE CONSIDERED CONVEX.
C
CS    PARAMETER ( GZERO  = 1.0E-5,  HZERO  = 1.0E-5                 )
CD    PARAMETER ( GZERO  = 1.0D-10, HZERO  = 1.0D-10                )
C
C  JUMP TO DIFFERENT PARTS OF THE SUBROUTINE DEPENDING ON
C  THE VALUE OF JUMPTO.
C
      GO TO ( 100, 200, 400, 450 ), JUMPTO
C
C  ON INITIAL ENTRY, SET CONSTANTS.
C
  100 CONTINUE
      PRNTER = IDEBUG .GE. 4 .AND. IOUT .GT. 0
      PRONEL = IDEBUG .EQ. 2 .AND. IOUT .GT. 0
      IF ( PRNTER ) WRITE( IOUT, 2010 ) 
      IF ( PRONEL ) WRITE( IOUT, 2120 ) 
C     ITA    = 0
      NBREAK = 0
      NFREED = 0
      NZERO  = N + 1
      EPSTL2 = TEN * EPSMCH
      EPSQRT = SQRT( EPSMCH )
      TBREAK = ZERO
C
C  IF NECESSARY, INITIALIZE THE DISTANCES TO THE SPHERICAL BOUNDARY.
C  PTP IS THE SUM OF THE SQUARES OF THE FREE COMPONENTS OF THE CAUCHY
C  DIRECTION. TPTTP IS THE SUM OF THE SQUARES OF THE DISTANCES TO THE 
C  BOUNDARY OF THE BOX FOR THE VARIABLES WHICH ARE FIXED AS THE
C  COMPUTATION PROCEEDS. 
C
      IF ( BOUNDX ) THEN
         DXSQR = ZERO
         PTP   = ZERO
         TPTTP = ZERO
      END IF
      IF ( IDEBUG .GE. 100 ) THEN
         DO 105 I = 1, N
            WRITE( IOUT, 2150 ) I, BND( I, 1 ), X0( I ), 
     *                             BND( I, 2 ), P( I )
  105    CONTINUE
      END IF
C
C  FIND THE STATUS OF THE VARIABLES.
C
CDIR$ IVDEP
      DO 120 I = 1, N
C
C  CHECK TO SEE WHETHER THE VARIABLE IS FIXED.
C
         IF ( INDEX( I ) .LE. 2 ) THEN
            INDEX( I ) = 0
            XUPPER     = BND( I, 2 ) - X0( I ) .LE. EPSTOL
            XLOWER     = X0( I ) - BND( I, 1 ) .LE. EPSTOL
            IF ( .NOT. ( XUPPER .OR. XLOWER ) ) THEN
C
C  THE VARIABLE LIES BETWEEN ITS BOUNDS. CHECK TO SEE IF THE SEARCH
C  DIRECTION IS ZERO.
C
               IF ( ABS( P( I ) ) .GT.   EPSMCH ) GO TO 110
               NZERO         = NZERO - 1
               IVAR( NZERO ) = I
            ELSE
               IF ( XLOWER ) THEN
C
C  THE VARIABLE LIES CLOSE TO ITS LOWER BOUND.
C
                  IF (    P( I )   .GT.   EPSMCH ) THEN
                     NFREED = NFREED + 1
                     GO TO 110
                  END IF
                  INDEX( I ) = 1
               ELSE
C
C  THE VARIABLE LIES CLOSE TO ITS UPPER BOUND.
C
                  IF (    P( I )   .LT. - EPSMCH ) THEN
                     NFREED = NFREED + 1
                     GO TO 110
                  END IF
                  INDEX( I ) = 2
               END IF
            END IF
         END IF
C
C  SET THE SEARCH DIRECTION TO ZERO.
C
         XT( I ) = X0( I )
         P( I )  = ZERO
C        IF ( PRNTER ) WRITE( IOUT, 2020 ) I, ZERO
         GO TO 120
  110    CONTINUE
C
C  IF THE VARIABLE IS FREE, SET UP THE POINTERS TO THE NONZEROS IN
C  THE VECTOR P READY FOR CALCULATING Q = B * P.
C
         NBREAK            = NBREAK + 1
         IVAR( NBREAK )    = I
         IF ( BOUNDX ) PTP = PTP + P( I ) * P( I )
  120 CONTINUE
C
C  RECORD THE NUMBER OF FREE VARIABLES AT THE STARTING POINT.
C
      NFREE = NBREAK
      NVAR2 = NFREE
      FXT   = F
C
C  IF ALL OF THE VARIABLES ARE FIXED, EXIT.
C
      IF ( PRONEL ) WRITE( IOUT, 2070 ) NFREED, N - NBREAK
      IF ( PRNTER ) WRITE( IOUT, 2110 ) NFREED, N - NBREAK
      IF ( NBREAK .EQ. 0 ) GO TO 600
      ITER = 0
C
C  FIND THE BREAKPOINTS FOR THE PIECEWISE LINEAR ARC (THE DISTANCES
C  TO THE BOUNDARY).
C
      DO 130 J = 1, NBREAK
         I     = IVAR( J )
         IF ( P( I ) .GT. EPSMCH ) THEN
            T = ( BND( I, 2 ) - X0( I ) ) / P( I )
         ELSE
            T = ( BND( I, 1 ) - X0( I ) ) / P( I )
         END IF
         BREAKP( J ) = T
  130 CONTINUE
C
C  ORDER THE BREAKPOINTS IN INCREASING SIZE USING A HEAPSORT.
C  BUILD THE HEAP.
C
CS    CALL SSORTA( NBREAK, BREAKP, IVAR, NFREE, .TRUE., INSORT )
CD    CALL DSORTA( NBREAK, BREAKP, IVAR, NFREE, .TRUE., INSORT )
C
C  RETURN TO THE MAIN ROUTINE TO EVALUATE Q = B * P.
C
      JUMPTO = 2
      NVAR1  = 1
      NVAR2  = NFREE
      RETURN
C
C  CALCULATE The FIRST DERIVATIVE (GXT) AND
C  SECOND DERIVATIVE (HXT) OF THE UNIVARIATE PIECEWISE QUADRATIC
C  FUNCTION AT THE START OF THE PIECEWISE LINEAR ARC.
C
  200 CONTINUE
      GXT      = ZERO
      HXT      = ZERO
      DO 210 J = 1, NFREE
         I     = IVAR( J )
         GXT   = GXT + G( I ) * P( I )
         HXT   = HXT + Q( I ) * P( I )
  210 CONTINUE
C
C  START THE MAIN LOOP TO FIND THE FIRST LOCAL MINIMIZER OF
C  THE PIECEWISE QUADRATIC FUNCTION. CONSIDER THE PROBLEM
C  OVER SUCCESSIVE PIECES.
C
  300 CONTINUE
C     ITA = ITA + 1
C
C  PRINT DETAILS OF THE PIECEWISE QUADRATIC IN THE NEXT INTERVAL.
C
      ITER = ITER + 1
      IF ( PRNTER ) WRITE( IOUT, 2030 )
     *                 ITER, FXT * FINDMX, GXT * FINDMX, HXT * FINDMX
      IF ( PRONEL ) WRITE( IOUT, 2080 )
     *                 ITER, FXT * FINDMX, GXT * FINDMX, HXT * FINDMX
C
C  IF THE GRADIENT OF THE UNIVARIATE FUNCTION INCREASES, EXIT.
C
      IF ( GXT .GT. GZERO ) GO TO 600
C
C  RECORD THE VALUE OF THE LAST BREAKPOINT.
C
      TK = TBREAK
C
C  FIND THE NEXT BREAKPOINT (END OF THE PIECE).
C
      TBREAK = BREAKP( 1 )
CS    CALL SSORTB( NBREAK, BREAKP, IVAR, NFREE, .TRUE., INSORT )
CD    CALL DSORTB( NBREAK, BREAKP, IVAR, NFREE, .TRUE., INSORT )
C
C  COMPUTE THE LENGTH OF THE CURRENT PIECE.
C
      DELTAT = TBREAK - TK
C
C  IF NECESSARY, COMPUTE THE DISTANCE TO THE SPHERICAL BOUNDARY.
C
      IF ( BOUNDX ) THEN
         TBNDSQ = SQRT( ( RADIUS ** 2 - TPTTP ) / PTP )
         IF ( TBNDSQ .LT. TBREAK ) DELTAT = TBNDSQ - TK
      END IF
C
C  PRINT DETAILS OF THE BREAKPOINT.
C
      IF ( PRNTER ) THEN
         IF ( BOUNDX ) THEN
            WRITE( IOUT, 2140 ) TBREAK, TBNDSQ
         ELSE
            WRITE( IOUT, 2040 ) TBREAK
         END IF
      END IF
C
C  IF THE GRADIENT OF THE UNIVARIATE FUNCTION IS SMALL AND
C  ITS CURVATURE IS POSITIVE, EXIT.
C
      IF ( ABS( GXT ) .LE. GZERO ) THEN
         IF ( HXT .GT. - HZERO .OR. DELTAT .GE. BIG ) THEN
            TCAUCH = TK
            GO TO 600
         ELSE
            TCAUCH = TBREAK
         END IF
      ELSE
C
C  IF THE GRADIENT OF THE UNIVARIATE FUNCTION IS NONZERO AND
C  ITS CURVATURE IS POSITIVE, COMPUTE THE LINE MINIMUM.
C
         IF ( HXT .GT. ZERO ) THEN
            TSTAR = - GXT / HXT
            IF ( PRNTER ) WRITE( IOUT, 2050 ) TSTAR
C
C  IF THE LINE MINIMUM OCCURS BEFORE THE BREAKPOINT, THE
C  LINE MINIMUM GIVES THE GENERALIZED CAUCHY POINT. EXIT.
C
            TCAUCH = MIN( TK + TSTAR, TBREAK )
            IF ( TSTAR .LT. DELTAT ) THEN
               DELTAT = TSTAR
               GO TO 500
            END IF
         ELSE
            TCAUCH = TBREAK
         END IF
      END IF
C
C  IF THE BREAKPOINT OCCURS ON THE SPHERICAL BOUNDARY, EXIT.
C
      IF ( BOUNDX ) THEN
         IF ( TBNDSQ .LE. TCAUCH ) THEN
            TCAUCH = TBNDSQ
            IF ( PRNTER ) WRITE( IOUT, 2130 ) 
            IF ( PRONEL ) WRITE( IOUT, 2130 ) 
            GO TO 500
         END IF
      END IF
C
C  UPDATE THE UNIVARIATE FUNCTION AND GRADIENT VALUES.
C
      FXT    = FXT + DELTAT * ( GXT + HALF * DELTAT * HXT )
      GXTOLD = GXT
      GXT    = GXT + DELTAT * HXT
C
C  RECORD THE NEW BREAKPOINT AND THE AMOUNT BY WHICH OTHER BREAKPOINTS
C  ARE ALLOWED TO VARY FROM THIS ONE AND STILL BE CONSIDERED TO BE
C  WITHIN THE SAME CLUSTER.
C
      FEASEP = TBREAK + EPSTL2
C
C  MOVE THE APPROPRIATE VARIABLE(S) TO THEIR BOUND(S).
C
  320 CONTINUE
      IBREAK = IVAR( NBREAK )
      NBREAK = NBREAK - 1
      IF ( PRNTER ) WRITE( IOUT, 2020 ) IBREAK, TBREAK
C
C  INDICATE THE STATUS OF THE NEWLY FIXED VARIABLE - THE VALUE
C  IS NEGATED TO INDICATE THAT THE VARIABLE HAS JUST BEEN FIXED.
C
      IF ( P( IBREAK ) .LT. ZERO ) THEN
         INDEX( IBREAK ) = - 1
      ELSE
         INDEX( IBREAK ) = - 2
      END IF
C
C  IF ALL OF THE REMAINING SEARCH DIRECTION IS ZERO, RETURN.
C
C ** Correction 2. 18/06/96: 1 line replaced by 20 **
      IF( NBREAK .EQ. 0 ) THEN
CDIR$ IVDEP
         DO 330 J = 1, NVAR2
            I     = IVAR( J )
C
C  RESTORE INDEX TO ITS CORRECT SIGN.
C
            INDEX( I ) = - INDEX( I )
C
C  MOVE THE VARIABLE ONTO ITS BOUND.
C
            XT( I ) = BND( I, INDEX( I ) )
C
C  STORE THE STEP FROM THE INITIAL POINT TO THE CAUCHY POINT IN P. 
C
            P( I ) = XT( I ) - X0( I )
  330    CONTINUE   
         NVAR2 = 0
         GO TO 600
      END IF
C ** Correction 2. 18/06/96:  end of correction **
C
C  DETERMINE IF OTHER VARIABLES HIT THEIR BOUNDS AT THE BREAKPOINT.
C
      IF ( BREAKP( 1 ) .LT. FEASEP ) THEN
CS       CALL SSORTB( NBREAK, BREAKP, IVAR, NFREE, .TRUE., INSORT )
CD       CALL DSORTB( NBREAK, BREAKP, IVAR, NFREE, .TRUE., INSORT )
         GO TO 320
      END IF
C
C  RETURN TO THE MAIN ROUTINE TO EVALUATE Q = B * P.
C
      JUMPTO = 3
      NVAR1  = NBREAK + 1
      RETURN
C
C  UPDATE THE FIRST AND SECOND DERIVATIVES OF THE UNIVARIATE FUNCTION.
C
  400 CONTINUE
C
C  START WITH THE SECOND-ORDER TERMS. ONLY PROCESS NONZERO
C  COMPONENTS OF Q.
C
CDIR$ IVDEP
      DO 410 J = 1, NNONNZ
         I     = INONNZ( J )
         QIPI  = Q( I ) * P( I )
         IF ( INDEX( I ) .EQ. 0 ) THEN
C
C  INCLUDE CONTRIBUTIONS FROM THE FREE COMPONENTS OF Q.
C
            GXT = GXT - QIPI * TBREAK
            HXT = HXT - QIPI * TWO
         ELSE
            IF ( INDEX( I ) .LT. 0 ) THEN
C
C  INCLUDE CONTRIBUTIONS FROM THE COMPONENTS OF Q WHICH WERE JUST FIXED.
C
               GXT  = GXT - QIPI * TBREAK
               HXT  = HXT - QIPI
            ELSE
C
C  INCLUDE CONTRIBUTIONS FROM THE COMPONENTS OF Q WERE PREVIOUSLY FIXED.
C
               GXT = GXT - QIPI
            END IF
         END IF
  410 CONTINUE
C
C  NOW INCLUDE THE CONTRIBUTIONS FROM THE VARIABLES WHICH
C  HAVE JUST BEEN FIXED.
C
CDIR$ IVDEP
      DO 420 J = NVAR1, NVAR2
         I     = IVAR( J )
         GXT   = GXT - P( I ) * G( I )
C
C  CONTINUE UPDATING THE DISTANCES TO THE SPHERICAL BOUNDARY.
C
         IF ( BOUNDX ) THEN
            TPTTP = TPTTP + ( TBREAK * P( I ) ) ** 2
            PTP   = PTP   - (          P( I ) ) ** 2
         END IF
C
C  RESTORE INDEX TO ITS CORRECT SIGN.
C
         INDEX( I ) = - INDEX( I )
C
C  MOVE THE VARIABLE ONTO ITS BOUND.
C
         XT( I ) = BND( I, INDEX( I ) )
C
C  STORE THE STEP FROM THE INITIAL POINT TO THE CAUCHY POINT IN P. 
C
         P( I ) = XT( I ) - X0( I )
  420 CONTINUE
C
C  COMPUTE THE SQUARE OF THE DISTANCE TO THE CURRENT POINT.
C
      IF ( BOUNDX ) DXSQR = TPTTP + PTP * TCAUCH ** 2
C
C  RESET THE NUMBER OF FREE VARIABLES.
C
      NVAR2 = NBREAK
C
C  CHECK THAT THE SIZE OF THE LINE GRADIENT HAS NOT SHRUNK 
C  SIGNIFICANTLY IN THE CURRENT SEGMENT OF THE  PIECEWISE ARC.
C  IF IT HAS, THERE MAY BE A LOSS OF ACCURACY, SO THE LINE
C  DERIVATIVES WILL BE RECOMPUTED.    
C
      RECOMP = ABS( GXT ) .LT. - EPSQRT * GXTOLD
C
C  IF REQUIRED, COMPUTE THE TRUE LINE GRADIENT AND CURVATURE.
C  FIRSTLY, COMPUTE THE MATRIX-VECTOR PRODUCT B * P.
C
      IF ( RECOMP .OR. PRNTER ) THEN
         JUMPTO = 4
         NVAR1  = 1
         RETURN
      END IF
  450 CONTINUE
C
C  CALCULATE THE LINE GRADIENT AND CURVATURE.
C
      IF ( RECOMP .OR. PRNTER ) THEN
         PBP      = ZERO
         GP       = ZERO
         IF ( IDEBUG .GT. 100 .AND. IOUT .GT. 0 ) WRITE( IOUT, 2100 )
     *      ( IVAR( I ), P(IVAR( I ) ), I = 1, NVAR2 ) 
         DO 460 J = 1, NVAR2
            I     = IVAR( J )
            QIPI  = P( I ) * Q( I )
            PBP   = PBP + QIPI
            GP    = GP +  P( I ) * G( I ) + TBREAK * QIPI
  460    CONTINUE
         DO 470 J = NVAR2 + 1, NFREE
            I     = IVAR( J )
            GP    = GP +  P( I ) * Q( I )
  470    CONTINUE
         IF ( PRNTER ) WRITE( IOUT, 2090 )
     *        GP * FINDMX, PBP * FINDMX, GXT * FINDMX, HXT * FINDMX
         IF ( RECOMP ) THEN
            GXT = GP
            HXT = PBP
         END IF
      END IF
C
C  JUMP BACK TO CALCULATE THE NEXT BREAKPOINT.
C
      GO TO 300
C
C  STEP TO THE GENERALIZED CAUCHY POINT.
C
  500 CONTINUE
C
C  CALCULATE THE FUNCTION VALUE FOR THE PIECEWISE QUADRATIC.
C
      FXT = FXT + DELTAT * ( GXT + HALF * DELTAT * HXT )
C
C  IF NECESSARY, UPDATE THE DISTANCES TO THE SPHERICAL BOUNDARY.
C
      IF ( BOUNDX ) DXSQR = TPTTP + PTP * DELTAT ** 2
      IF ( PRNTER ) WRITE( IOUT, 2060 ) FXT * FINDMX
C
C  THE GENERALIZED CAUCHY POINT HAS BEEN FOUND. SET THE ARRAY P
C  TO THE STEP FROM THE INITIAL POINT TO THE CAUCHY POINT.
C
  600 CONTINUE
CDIR$ IVDEP
      DO 610 J   = 1, NVAR2
         I       = IVAR( J )
         P( I )  = P( I ) * TCAUCH
         XT( I ) = X0( I ) + P( I )
  610 CONTINUE
C
C  RECORD THAT VARIABLES WHOSE GRADIENTS WERE ZERO AT THE INITIAL
C  POINT ARE FREE VARIABLES.
C
      DO 620 J = NZERO, N
         NFREE = NFREE + 1
         IVAR( NFREE ) = IVAR( J )
  620 CONTINUE
C
C  SET RETURN CONDITIONS.
C
      JUMPTO = 0
      NVAR1  = 1
      NVAR2  = NFREE
      RETURN
C
C  NON-EXECUTABLE STATEMENTS.
C
 2010 FORMAT(/ ' ----------------- CAUCH entered -------------------')
 2020 FORMAT(    ' Variable ', I4, ' is fixed, step = ', 1P, D12.4 )
 2030 FORMAT( /, ' Piece', I5, ' - F, G and H at start point ',
     *        1P, 3D12.4 )
 2040 FORMAT( /, ' Next break point = ', 1P, D12.4 )
 2050 FORMAT( /, ' Stationary point = ', 1P, D12.4 )
 2060 FORMAT( /, ' Function value at the Cauchy point ', 1P, D12.4 )
 2070 FORMAT(  /, '  ', I7, ' vars. freed ', I7, ' vars. remain fixed' )
 2080 FORMAT( 25X, I7, 1P, 3D12.4 )
 2090 FORMAT( /, ' Calculated GXT and HXT = ', 1P, 2D12.4,
     *        /, ' Recurred   GXT and HXT = ',     2D12.4 )
 2100 FORMAT( ' Current search direction ', /, ( 4( I6, 1P, D12.4 ) ) )
 2110 FORMAT( /, I8, ' variables freed from their bounds ', 
     *        /, I8, ' variables remain fixed ', / )
 2120 FORMAT( /, '    ** CAUCH entered ** Segment    Model   ',
     *           '   Gradient   Curvature ')
 2130 FORMAT( /, ' Spherical trust region encountered ', / )
 2140 FORMAT( /, ' Next break point = ', 1P, D12.4, 
     *           ' Spherical boundary = ', 1P, D12.4 )
 2150 FORMAT( ' VAR L X U P ', I6, 1P, 4D12.4 )
C
C  END OF SUBROUTINE CAUCH.
C
      END
