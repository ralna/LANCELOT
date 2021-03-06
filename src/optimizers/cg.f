C ** Correction report.
C ** Correction 1. 13/01/94: 2 lines corrected **
C ** Correction 2. 27/08/96: 1 line corrected **
C ** End of Correction report.
C  THIS VERSION: 27/08/1996 AT 08:00:00 AM.
C
C  ** FOR THE CRAY 2, LINES STARTING'CDIR$ IVDEP' TELL THE COMPILER TO
C     IGNORE POTENTIAL VECTOR DEPENDENCIES AS THEY ARE KNOWN TO BE O.K.
C
CS    SUBROUTINE SCG   ( N , X0, XT, G , BND   , INDEX , EPSGRD, FMODEL, 
CD    SUBROUTINE DCG   ( N , X0, XT, G , BND   , INDEX , EPSGRD, FMODEL,
     *                   XSCALE, GRAD  , BREAKP, INFORM, P , Q , IVAR  ,
     *                   NVAR  , NRESRT, BOUNDX, GNRMSQ, DXSQR , IOUT  ,
     *                   JUMPTO, IDEBUG )
C
C  *******************************************************************
C
C  FIND AN APPROXIMATION TO THE MINIMIZER OF THE QUADRATIC FUNCTION
C
C     0.5 (X-X0) (TRANSPOSE ) B (X-X0) + G (TRANSPOSE) (X-X0) + F
C
C  WITHIN THE BOX REGION BND(*,1) <= X(*) <= BND(*,2) SUBJECT TO THE
C  FURTHER RESTRICTION THAT SPECIFIED VARIABLES REMAIN FIXED AT
C  THEIR INITIAL VALUES.
C
C  ** VERSION TO ALLOW RESTARTS WHEN BOUNDS ARE ENCOUNTERED.
C
C  OPTIONALLY, THE SEARCH MAY BE TERMINATED AT THE FIRST POINT
C  ENCOUNTERED ON THE BOUNDARY OF THE SPHERICAL REGION
C  ||X - X0|| <= R, WHERE ||.|| DENOTES THE 2-NORM. FURTHERMORE,
C  THE MINIMIZATION MAY CONTINUE ALONG THE EDGES OF THE BOX IF
C  SO DESIRED.
C
C  CONTROL IS PASSED
C  FROM THE ROUTINE WHENEVER A PRODUCT OF THE VECTOR P WITH B IS
C  REQUIRED, AND THE USER IS RESPONSIBLE FOR FORMING THE PRODUCT
C  IN THE VECTOR Q. CONTROL IS ALSO PASSED TO THE USER WHEN THE
C  'PRECONDITIONED' GRADIENT W ** (-1) GRAD IS TO BE FORMED IN Q.
C  HERE, W IS ANY POSITIVE DEFINITE AND SYMMETRIC APPROXIMATION TO
C  THE HESSIAN B.
C
C  PROGRESS THROUGH THE ROUTINE IS CONTROLLED BY THE PARAMETER
C  JUMPTO.
C
C  IF JUMPTO = 0, NO FURTHER ENTRY IS CALLED FOR.
C  IF JUMPTO = 1, AN INITIAL ENTRY HAS BEEN MADE.
C  IF JUMPTO = 2, THE VECTOR Q = W ** (-1) * GRAD IS REQUIRED.
C  IF JUMPTO = 3, THE VECTOR Q = B * P IS REQUIRED.
C  IF JUMPTO = 4, THE NORM OF THE GRADIENT OF THE CURRENT
C                 ITERATE IS SMALLER THAN EPSGRD. THE USER MAY
C                 WISH TO PERFORM ADDITIONAL CONVERGENCE TESTS
C                 IN THE CALLING PROGRAM AND RE-ENTER CG WITH A
C                 SMALLER VALUE OF EPSGRD. IF SUCH A RE-ENTRY IS
C                 REQUIRED, JUMPTO SHOULD NOT BE ALTERED.
C                 OTHERWISE, JUMPTO SHOULD BE RESET TO 0.
C  IF JUMPTO = 5, AN EDGE OF THE BOX HAS BEEN ENCOUNTERED.
C                 IF THE USER WISHES TO CONTINUE THE MINIMIZ-
C                 -ATION ALONG THE EDGE, JUMPTO SHOULD BE RESET
C                 TO 2 AND THE THE VECTOR Q CALCULATED AS ABOVE
C                 BEFORE RE-ENTRY. OTHERWISE, JUMPTO SHOULD BE
C                 RESET TO 0.
C
C  THE STATUS OF THE ARRAY INDEX GIVES THE STATUS OF THE VARIABLES.
C
C  IF INDEX( I ) = 0, THE I-TH VARIABLE IS FREE.
C  IF INDEX( I ) = 1, THE I-TH VARIABLE IS FIXED ON ITS LOWER BOUND.
C  IF INDEX( I ) = 2, THE I-TH VARIABLE IS FIXED ON ITS UPPER BOUND.
C  IF INDEX( I ) = 3, THE I-TH VARIABLE IS PERMANENTLY FIXED.
C  IF INDEX( I ) = 4, THE I-TH VARIABLE IS FIXED AT SOME OTHER VALUE.
C
C  THE ADDRESSES OF THE FREE VARIABLES ARE GIVEN IN THE FIRST
C  NVAR ENTRIES OF THE ARRAY IVAR.
C
C  IF THE PRODUCT W ** (-1) * GRAD IS REQUIRED (JUMPTO = 2), THE
C  NONZERO COMPONENTS OF GRAD OCCUR IN POSITIONS IVAR(I) FOR
C  I LYING BETWEEN ONE AND NVAR AND HAVE THE VALUES GRAD(I).
C
C  IF THE PRODUCT B * P IS REQUIRED (JUMPTO = 3), THE
C  NONZERO COMPONENTS OF P OCCUR IN POSITIONS IVAR(I) FOR
C  I LYING BETWEEN ONE AND NVAR.
C
C  NICK GOULD, 11TH AUGUST 1987.
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
C  THE PRECONDITIONER :
C  --------------------
C
C  WHEN JUMPTO = 2, THE USER IS ASKED TO PROVIDE THE PRODUCT OF
C  THE INVERSE OF A 'PRECONDITIONING' MATRIX W AND THE VECTOR GRAD.
C  THE MATRIX W SHOULD HAVE THE FOLLOWING PROPERTIES:
C
C  A) W MUST BE SYMMETRIC AND POSITIVE DEFINITE.
C  B) W SHOULD BE AN APPROXIMATION TO THE SECOND DERIVATIVE MATRIX
C     B IN THE SENSE THAT THE EIGENVALUES OF W ** (-1) B SHOULD
C     BE CLOSE TO UNITY WHENEVER B IS POSITIVE DEFINITE. IDEALLY, THESE
C     EIGENVALUES SHOULD BE CLUSTERED ABOUT A SMALL NUMBER OF DISTINCT
C     VALUES.
C  C) W SHOULD HAVE THE PROPERTY THAT THE PRODUCT W ** (-1) * GRAD
C     IS CHEAP TO OBTAIN. BY THIS WE MEAN THAT THE AMOUNT OF
C     ARITHMETIC WORK INVOLVED IN FORMING THE PRODUCT SHOULD BE A
C     SMALL MULTIPLE OF N.
C
C  POPULAR PRECONDITIONERS INCLUDE THE DIAGONAL PRECONDITIONER,
C  IN WHICH W IS A DIAGONAL MATRIX WITH DIAGONAL VALUES EQUAL TO
C  THOSE OF B (MODIFIED TO BE POSITIVE IF NECESSARY) AND THE
C  INCOMPLETE FACTORIZATION PRECONDITIONER, IN WHICH ELEMENTS
C  OF B ARE CHANGED TO ZEROS IN W IN ORDER TO REDUCE THE NUMBER OF
C  NONZEROS IN THE CHOLESKY FACTORIZATION OF W. UNFORTUNATELY,
C  THE CHOICE OF A GOOD PRECONDITIONER IS PROBLEM DEPENDENT.
C
C  PARAMETERS : (REAL MEANS DOUBLE PRECISION IN THE DOUBLE VERSION)
C  ------------
C
C  N      (INTEGER) THE NUMBER OF INDEPENDENT VARIABLES.
C          ** THIS VARIABLE IS NOT ALTERED BY THE SUBROUTINE.
C  X0     (REAL ARRAY OF LENGTH AT LEAST N) THE POINT X0.
C          ** THIS VARIABLE IS NOT ALTERED BY THE SUBROUTINE.
C  XT     (REAL ARRAY OF LENGTH AT LEAST N) THE CURRENT ESTIMATE
C          OF THE MINIMIZER FOR THE PROBLEM. XT SHOULD BE INITIALIZED
C          AT THE USER'S ESTIMATE OF THE MINIMIZER AND MUST SATISFY
C          THE RESTRICTION BND(1,I) <= XT(I) <= BND(2,I), I=1,...,N
C  G      (REAL ARRAY OF LENGTH AT LEAST N) THE COEFFICIENTS OF
C          THE LINEAR TERM IN THE QUADRATIC FUNCTION.
C          ** THIS VARIABLE IS NOT ALTERED BY THE SUBROUTINE.
C  BND    (TWO DIMENSIONAL REAL ARRAY WITH LEADING DIMENSION N AND
C          SECOND DIMENSION 2) THE LOWER (BND(*,1)) AND
C          UPPER (BND(*,2)) BOUNDS ON THE VARIABLES.
C          ** THIS VARIABLE IS NOT ALTERED BY THE SUBROUTINE.
C  INDEX  (INTEGER ARRAY OF LENGTH AT LEAST N) SPECIFIES WHICH
C          OF THE VARIABLES ARE TO BE FIXED AT THE START OF THE
C          MINIMIZATION. INDEX SHOULD BE SET AS FOLLOWS:
C          IF INDEX( I ) = 0, THE I-TH VARIABLE IS FREE.
C          IF INDEX( I ) = 1, THE I-TH VARIABLE IS ON ITS LOWER BOUND.
C          IF INDEX( I ) = 2, THE I-TH VARIABLE IS ON ITS UPPER BOUND.
C          IF INDEX( I ) = 3, 4, THE I-TH VARIABLE IS FIXED AT XT(I).
C  EPSGRD (REAL) THE MINIMIZATION WILL CONTINUE UNTIL THE NORM OF
C          THE SCALED GRADIENT IS LESS THAN SQRT(EPSGRD).
C          ** THIS VARIABLE IS NOT ALTERED BY THE SUBROUTINE.
C  FMODEL (REAL) THE VALUE OF THE QUADRATIC FUNCTION AT THE CURRENT
C          ESTIMATE OF THE MINIMIZER. FMODEL SHOULD BE INITIALIZED AT
C          THE VALUE OF THE QUADRATIC AT THE INITIAL POINT.
C  XSCALE (REAL ARRAY OF LENGTH N) THE SCALE FACTORS USED IN
C         COMPUTING THE NORM OF THE SCALED GRADIENT, GNORM =
C         SQRT( SUM(FROM I=1)(TO NVAR) (GRAD(I) * XSCALE(IVAR(I)))**2 ).
C          ** THIS VARIABLE IS NOT ALTERED BY THE SUBROUTINE.
C  GRAD   (REAL ARRAY OF LENGTH N) THE VALUE OF NONZERO COMPONENTS
C          OF THE GRADIENT OF THE QUADRATIC FUNCTION AT THE CURRENT
C          ESTIMATE OF THE MINIMIZER. THE COMPONENT GRAD(I), I = 1,..,
C          NVAR, GIVES THE VALUE OF THE IVAR(I)-TH COMPONENT OF THE
C          GRADIENT OF THE QUADRATIC FUNCTION (SEE IVAR, NVAR).
C  BREAKP (REAL ARRAY OF LENGTH N) IS USED AS WORKSPACE.
C  INFORM (INTEGER) THE VALUE OF INFORM ON TERMINAL (JUMPTO=0) EXIT
C          DESCRIBES THE PROGRESS OF THE ALGORITHM AS FOLLOWS:
C          IF INFORM = 10, THE OUTPUT VALUE XT CONTAINS THE REQUIRED
C                          MINIMIZER (SEE EPSGRD).
C          IF INFORM = 11, TOO MANY ITERATIONS HAVE BEEN PERFORMED
C                          BY THE ROUTINE WITHOUT SOLVING THE PROBLEM.
C                          THE BEST ESTIMATE OF THE SOLUTION IS
C                          GIVEN IN XT.
C          IF INFORM = 12, A FREE VARIABLE HAS ENCOUNTERED ONE OF ITS
C                          BOUNDS. THE BEST ESTIMATE OF THE SOLUTION
C                          FOUND TO DATE IS GIVEN IN XT.
C          IF INFORM = 13, A DIRECTION OF NEGATIVE CURVATURE HAS BEEN
C                          DETERMINED AND A STEP TO ONE OF THE BOUNDS
C                          HAS BEEN TAKEN. THE BEST ESTIMATE OF THE
C                          SOLUTION OBTAINED IS GIVEN IN XT.
C          IF INFORM = 14, THE LAST STEP TAKEN IS SO SMALL THAT FURTHER
C                          PROGRESS IS UNLIKELY. THE BEST ESTIMATE OF
C                          THE SOLUTION IS GIVEN IN XT.
C  P      (REAL ARRAY OF LENGTH AT LEAST N) CONTAINS THE VALUES OF THE
C          COMPONENTS OF THE VECTOR P. ON AN INITIAL ENTRY (JUMPTO=1),
C          THE COMPONENTS OF PMUST BE SET TO ZERO. ON A NON OPTIMAL EXIT
C          (JUMPTO=3), P IS THE VECTOR FOR WHICH THE PRODUCT B*P
C          IS REQUIRED BEFORE THE NEXT RE-ENTRY. THE COMPONENTS
C          IVAR(I), I = 1, ... , NVAR OF P CONTAIN THE VALUES OF THE
C          NONZERO COMPONENTS OF P (SEE, IVAR, NVAR).
C  Q      (REAL ARRAY OF LENGTH AT LEAST N) ON THE INITIAL ENTRY
C          (JUMPTO=1), Q MUST CONTAIN THE PRODUCT B * ( XT - X0 ).
C          Q MUST CONTAIN THE VECTOR W ** (-1) *  P (JUMPTO=2) OR
C          B * P (JUMPTO=3), ON A NON INITIAL ENTRY. FOR NON INITIAL
C          ENTRIES, ONLY THE COMPONENTS WITH INDICES IVAR(I) I=1,..,
C          NVAR NEED BE SET (THE OTHER COMPONENTS ARE NOT USED).
C  IVAR   (INTEGER ARRAY OF LENGTH AT LEAST N) ON ALL NORMAL EXITS
C         (JUMPTO=0,2,3), IVAR(I), I=1,...,NVAR, GIVES THE INDICES
C          OF THE NONZERO COMPONENTS OF P AND GRAD.
C  NVAR   (INTEGER) SEE IVAR, ABOVE.
C  NRESRT (INTEGER) THE SEARCH FOR THE MINIMIZER WILL BE RESTARTED IN 
C          THE PRECONDITIONED STEEPEST-DESCENT DIRECTION EVERY NRESRT
C          ITERATIONS.
C          ** THIS VARIABLE IS NOT ALTERED BY THE SUBROUTINE.
C  BOUNDX (LOGICAL) THE SEARCH FOR THE MINIMIZER WILL
C          BE TERMINATED ON THE BOUNDARY OF THE SPHERICAL
C          REGION ||X-X0|| <= R IF AND ONLY IF BOUNDX IS SET TO
C          .TRUE. ON INITIAL (JUMPTO=1) ENTRY.
C          ** THIS VARIABLE IS NOT ALTERED BY THE SUBROUTINE.
C  GNRMSQ (REAL) THE NORM OF THE PRECONDITIONED GRADIENT AT THE CURRENT
C          ESTIMATE OF THE MINIMIZER.
C  DXSQR  (REAL) THE SQUARE OF THE TWO NORM OF THE DISTANCE BETWEEN
C          THE CURRENT ESTIMATE OF THE MINIMIZER AND X0. DXSQR WILL
C          ONLY BE SET IF BOUNDX IS .TRUE.
C  IOUT   (INTEGER) THE STANDARD FORTRAN OUTPUT CHANNEL TO BE USED.
C  JUMPTO (INTEGER) CONTROLS FLOW THROUGH THE SUBROUTINE.
C          IF JUMPTO = 0, NO FURTHER ENTRY IS CALLED FOR.
C          IF JUMPTO = 1, AN INITIAL ENTRY HAS BEEN MADE.
C          IF JUMPTO = 2, THE VECTOR Q = W ** (-1) GRAD IS REQUIRED.
C          IF JUMPTO = 3, THE VECTOR Q = B * P IS REQUIRED.
C          IF JUMPTO = 4, THE NORM OF THE GRADIENT OF THE CURRENT
C                         ITERATE IS SMALLER THAN EPSGRD. THE USER MAY
C                         WISH TO PERFORM ADDITIONAL CONVERGENCE TESTS
C                         IN THE CALLING PROGRAM AND RE-ENTER CG WITH A
C                         SMALLER VALUE OF EPSGRD. IF SUCH A RE-ENTRY IS
C                         REQUIRED, JUMPTO SHOULD NOT BE ALTERED.
C                         OTHERWISE, JUMPTO SHOULD BE RESET TO 0.
C          IF JUMPTO = 5, AN EDGE OF THE BOX HAS BEEN ENCOUNTERED.
C                         IF THE USER WISHES TO CONTINUE THE MINIMIZ-
C                         -ATION ALONG THE EDGE, JUMPTO SHOULD BE RESET
C                         TO 2 AND THE THE VECTOR Q CALCULATED AS ABOVE
C                         BEFORE RE-ENTRY. OTHERWISE, JUMPTO SHOULD BE
C                         RESET TO 0.
C  IDEBUG (INTEGER) ALLOWS DETAILED PRINTING. IF IDEBUG IS LARGER
C          THAN 4, DETAILED OUTPUT FROM THE ROUTINE WILL BE GIVEN.
C          OTHERWISE, NO OUTPUT OCCURS.
C
C  COMMON VARIABLES.
C  -----------------
C
C  THE COMMON BLOCK SCOMSB/DCOMSB IS USED TO PASS FURTHER INFORMATION.
C  THE BLOCK CONTAINS THE VARIABLES
C
C    ACCCG, RATIO, RADIUS, RADMAX, FINDMX, CMA31, ITERCG, NGEVAL, 
C    ISKIP, IFIXED, NBANDW
C
C  IN THAT ORDER.  VARIABLES ITERCG, NGEVAL, ISKIP, IFIXED AND
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
C  ITERCG GIVES THE NUMBER OF CONJUGATE GRADIENT ITERATIONS TAKEN
C  WHILE IFIXED GIVES THE VARIABLE WHICH MOST RECENTLY HITS A BOUND.
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
C  --------------
C
C  ON THE INITIAL CALL TO THE SUBROUTINE THE FOLLOWING VARIABLES MUST
C  BE SET BY THE USER: -
C
C      N, X0, XT, G, BND, FMODEL, EPSGRD, BOUNDX, P, Q, IOUT, JUMPTO,
C      XSCALE.
C
C  JUMPTO MUST HAVE THE VALUE 1.
C  IN ADDITION, IF THE I-TH VARIABLE IS REQUIRED TO BE FIXED AT
C  ITS INITIAL VALUE, X0(I), INDEX(I) MUST BE SET TO 3 OR 4.
C  RADIUS MUST BE SPECIFIED IF BOUNDX IS .TRUE. ON INITIAL ENTRY.
C
C  RE-ENTRY.
C  ---------
C
C  IF THE VARIABLE JUMPTO HAS THE VALUE 2 ON EXIT, THE
C  SUBROUTINE MUST BE RE-ENTERED WITH THE VECTOR Q CONTAINING
C  THE PRODUCT OF THE INVERSE OF THE PRECONDITIONING MATRIX W AND
C  THE OUTPUT VECTOR GRAD. ALL OTHER PARAMETERS MUST NOT BE ALTERED.
C                                               -------------------
C  IF THE VARIABLE JUMPTO HAS THE VALUE 3 ON EXIT, THE
C  SUBROUTINE MUST BE RE-ENTERED WITH THE VECTOR Q CONTAINING
C  THE PRODUCT OF THE SECOND DERIVATIVE MATRIX B AND THE OUTPUT
C  VECTOR P. ALL OTHER PARAMETERS MUST NOT BE ALTERED.
C                                 -------------------
C
C  ****************************************************************
C
      INTEGER           N     , INFORM, NVAR  , NRESRT, IOUT  , JUMPTO
      INTEGER           IDEBUG
CS    REAL              EPSGRD, FMODEL, GNRMSQ, DXSQR
CD    DOUBLE PRECISION  EPSGRD, FMODEL, GNRMSQ, DXSQR
      LOGICAL           BOUNDX
      INTEGER           INDEX ( N ), IVAR  ( N )
CS    REAL              G     ( N ), GRAD  ( N ), X0   ( N ), 
CD    DOUBLE PRECISION  G     ( N ), GRAD  ( N ), X0   ( N ),
     *                  XT    ( N ), XSCALE( N ), BND  ( N, 2 ),
     *                  P     ( N ), Q     ( N ), BREAKP( N )
C
C  LOCAL VARIABLES.
C
      INTEGER           I , J , ITER  , ITSLE , NEWVAR
CS    REAL              STEPMX, STEPIB, STEPBD, STEPTR, ALPHA , BETA,
CD    DOUBLE PRECISION  STEPMX, STEPIB, STEPBD, STEPTR, ALPHA , BETA,
     *                  DXTP  , OLDGNS, HALF  , PI    , PTP   , CURVAT,
     *                  TWO   , ZERO  , ONE   , ONEPEP, TEN   , REALGR
      LOGICAL           PRNTER, PRONEL
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
      INTRINSIC         MOD   , MIN   , SQRT
C
C  THE FOLLOWING VARIABLES ARE SAVED BETWEEN CALLS OF THE ROUTINE.
C
      SAVE              ITER  , ITSLE , PRNTER, PRONEL 
      SAVE              OLDGNS, ONEPEP, ALPHA 
C ** Correction 1. 13/01/94: 2 lines corrected **
CS    SAVE            / SCOMSB /, / SMACHN /
CD    SAVE            / DCOMSB /, / DMACHN /
C ** Correction 1. 13/01/94:  end of correction **
C
C  SET CONSTANT REAL PARAMETERS.
C
CS    PARAMETER ( ZERO   = 0.0E+0, HALF   = 5.0E-1, ONE    = 1.0E+0 )
C ** Correction 2. 27/08/96: 1 line corrected **
CD    PARAMETER ( ZERO   = 0.0D+0, HALF   = 5.0D-1, ONE    = 1.0D+0 )
C ** Correction 2. 27/08/96: end of correction **
CS    PARAMETER ( TWO    = 2.0E+0, TEN    = 1.0E+1                  )
CD    PARAMETER ( TWO    = 2.0D+0, TEN    = 1.0D+1                  )
C
C  PASS CONTROL TO DIFFERENT PARTS OF THE SUBROUTINE DEPENDING
C  ON THE VALUE OF JUMPTO.
C
      IF ( JUMPTO .LE. 0 ) RETURN
      GO TO ( 100, 300, 400, 330 ), JUMPTO
C
C  ON INITIAL ENTRY, SET CONSTANTS.
C
  100 CONTINUE
      ITER   = 0
      ITSLE  = 0
      NVAR   = 0
      ONEPEP = ONE + TEN * EPSMCH
      ALPHA  = ZERO
      GNRMSQ = ZERO
      PRNTER = IDEBUG .GE. 4 .AND. IOUT .GT. 0
      PRONEL = IDEBUG .EQ. 2 .AND. IOUT .GT. 0
C
C  IF NECESSARY, CHECK THAT XT LIES WITHIN THE SPHERICAL BOUNDARY.
C
      IF ( BOUNDX ) THEN
         IF ( DXSQR .GT. RADIUS ** 2 ) THEN
            IF ( PRNTER ) WRITE( IOUT, 2080 )
            IF ( PRONEL ) WRITE( IOUT, 2190 )
            INFORM = 12
            JUMPTO = 0
            GO TO 500
         END IF
      END IF
      DO 110 I  = 1, N
C
C  PLACE THE ADDRESSES OF THE FREE VARIABLES IN IVAR AND THE GRADIENT
C  OF THE QUADRATIC WITH RESPECT TO THE FREE VARIABLES IN GRAD.
C
         IF ( INDEX( I ) .EQ. 0 ) THEN
            NVAR         = NVAR   + 1
            IVAR( NVAR ) = I
            GRAD( NVAR ) = G( I ) + Q( I )
         END IF
  110 CONTINUE
C
C  IF THERE ARE NO FREE VARIABLES, EXIT.
C
      IF ( NVAR .EQ. 0 ) THEN
C
C  PRINT DETAILS OF THE CURRENT ITERATION.
C
         IF ( PRNTER ) THEN
            IF ( ITSLE .EQ. 0 ) WRITE( IOUT, 2000 )
            WRITE( IOUT, 2160 ) ITER, FMODEL * FINDMX, ZERO, ZERO,
     *                          EPSGRD
         END IF
         IF ( PRONEL ) THEN
            IF ( ITSLE .EQ. 0 ) WRITE( IOUT, 2110 ) EPSGRD
            WRITE( IOUT, 2170 ) ITER, FMODEL * FINDMX, ZERO, ZERO, ZERO
         END IF
         IF ( PRNTER ) WRITE( IOUT, 2050 ) ITER
         INFORM = 10
         JUMPTO = 4
         GO TO 500
      END IF
C
C  SET AN UPPER BOUND ON THE NUMBER OF C.G. ITERATIONS ALLOWED.
C
      ITCGMX = 3 * NVAR
C     ITCGMX = 50
C
C  START THE MAIN ITERATION.
C
  200 CONTINUE
C
C  RETURN TO THE MAIN ROUTINE TO OBTAIN THE PRECONDITIONED
C  GRADIENT Q = M ** (-1) GRAD.
C
      JUMPTO = 2
      RETURN
C
C  COMPUTE THE NORM OF THE PRECONDITIONED GRADIENT.
C
  300 CONTINUE
      GNRMSQ    = ZERO
      DO 310 I  = 1, NVAR
         GNRMSQ = GNRMSQ + Q( IVAR( I ) ) * GRAD( I )
  310 CONTINUE
      IF ( IDEBUG .GE. 100 ) THEN
         DO 315 I = 1, NVAR
            J     = IVAR( I )
            WRITE( IOUT, 2210 ) J, BND( J, 1 ), XT( J ), 
     *                             BND( J, 2 ), GRAD( I )
  315    CONTINUE
      END IF
C
C  IF THE PRECONDITIONED GRADIENT IS SUFFICIENTLY SMALL, EVALUATE
C  THE NORM OF THE TRUE GRADIENT OF THE QUADRATIC MODEL.
C
      IF ( GNRMSQ .LE. EPSGRD ) THEN
         REALGR    = ZERO
         DO 320 I  = 1, NVAR
            REALGR = MAX( REALGR,
     *                  ( GRAD( I ) * XSCALE( IVAR( I ) ) ) ** 2 )
  320    CONTINUE
C
C  PRINT DETAILS OF THE CURRENT ITERATION.
C
         IF ( PRNTER ) THEN
            IF ( ITSLE .EQ. 0 ) WRITE( IOUT, 2000 )
            WRITE( IOUT, 2160 )
     *         ITER, FMODEL * FINDMX, GNRMSQ, REALGR, EPSGRD
         END IF
         IF ( PRONEL ) THEN
            IF ( ITSLE .EQ. 0 ) WRITE( IOUT, 2110 ) EPSGRD
            WRITE( IOUT, 2170 ) ITER, FMODEL * FINDMX, GNRMSQ, REALGR,
     *                          ALPHA
         END IF
C
C  IF THE GRADIENT OF THE MODEL IS SUFFICIENTLY SMALL, EXIT.
C
         IF ( REALGR .LE. EPSGRD ) THEN
            IF ( PRNTER ) WRITE( IOUT, 2050 ) ITER
            INFORM = 10
            JUMPTO = 4
            GO TO 500
         END IF
      ELSE
C
C  PRINT DETAILS OF THE CURRENT ITERATION.
C
         IF ( PRNTER ) THEN
            IF ( ITSLE .EQ. 0 ) WRITE( IOUT, 2000 )
            WRITE( IOUT, 2020 ) ITER, FMODEL * FINDMX, GNRMSQ, EPSGRD
         END IF
         IF ( PRONEL ) THEN
            IF ( ITSLE .EQ. 0 ) WRITE( IOUT, 2110 ) EPSGRD
            WRITE( IOUT, 2120 ) ITER, FMODEL * FINDMX, GNRMSQ, ALPHA 
         END IF
      END IF
C
C  START OF THE NEW ITERATION.
C
  330 CONTINUE
      ITER  = ITER + 1
      ITSLE = ITSLE + 1
C
C  ---------------------------------------------------------------------
C                COMPUTE THE SEARCH DIRECTION
C  ---------------------------------------------------------------------
C
      IF ( MOD( ITER, NRESRT + 1 ) .EQ. 1 ) THEN
         DO 340 J  = 1, NVAR
            I      = IVAR( J )
            P( I ) = - Q( I )
  340    CONTINUE
      ELSE
C
C  CALCULATE THE SCALE FACTOR, BETA, WHICH MAKES THE SEARCH DIRECTION
C  B-CONJUGATE TO ITS PREDECESSORS.
C
         BETA      = GNRMSQ / OLDGNS
CDIR$ IVDEP
         DO 350 J  = 1, NVAR
            I      = IVAR( J )
            P( I ) = - Q( I ) + BETA * P( I )
  350    CONTINUE
      END IF
C
C  SAVE THE NORM OF THE PRECONDITIONED GRADIENT.
C
      OLDGNS = GNRMSQ
C
C  RETURN TO THE MAIN ROUTINE TO EVALUATE Q = B * P.
C
      JUMPTO = 3
      RETURN
C
C  COMPUTE THE MAXIMUM POSSIBLE STEP, STEPMX.
C
  400 CONTINUE
      CURVAT = ZERO
      STEPBD = BIG
C
C  ---------------------------------------------------------------------
C                COMPUTE THE CURVATURE
C  ---------------------------------------------------------------------
C
CDIR$ IVDEP
      DO 410 J  = 1, NVAR
C
C  COMPUTE THE CURVATURE ALONG THE CURRENT SEARCH DIRECTION.
C
         I      = IVAR( J )
         PI     = P( I )
         CURVAT = CURVAT + Q( I ) * PI
C
C  FIND THE DISTANCE TO THE I-TH BOUND, STEPIB, AND SAVE IT FOR LATER.
C  TAKE PRECAUTIONS IF THE CURRENT ITERATE IS VERY CLOSE TO A BOUND.
C
         IF ( PI .GT.   EPSMCH ) THEN
            STEPIB = ( EPSMCH + ( BND( I, 2 ) - XT( I ) ) ) / PI
         ELSE
            IF ( PI .LT. - EPSMCH ) THEN
               STEPIB = ( - EPSMCH + ( BND( I, 1 ) - XT( I ) ) ) / PI
            ELSE
               STEPIB = BIG
            END IF
         END IF
         IF ( STEPIB .LE. ZERO ) STEPIB = EPSMCH / ABS( PI )
         BREAKP( I )                    = STEPIB
C        IF ( PRNTER ) WRITE( 6, 2180 ) I, XT( I ), P( I ), BREAKP( I )
C
C  FIND THE DISTANCE TO NEAREST BOUND, STEPBD.
C
         STEPBD = MIN( STEPBD, STEPIB )
  410 CONTINUE
C
C  IF NECESSARY, COMPUTE THE TERMS NEEDED TO FIND THE DISTANCE
C  TO THE SPHERICAL BOUNDARY.
C
      STEPMX = STEPBD
      IF ( BOUNDX ) THEN
         DXTP = ZERO
         PTP  = ZERO
CDIR$ IVDEP
         DO 420 J  = 1, NVAR
            I      = IVAR( J )
            PI     = P( I )
            DXTP = DXTP + PI * ( XT( I ) - X0( I ) )
            PTP  = PTP  + PI * PI
  420    CONTINUE
C
C  NOW COMPUTE THE DISTANCE TO THE SPHERICAL BOUNDARY, STEPTR, AND
C  FIND THE SMALLER OF THIS AND THE DISTANCE TO THE BOUNDARY OF THE
C  FEASIBLE BOX, STEPMX.
C
         STEPTR = ( SQRT( DXTP * DXTP - PTP * ( DXSQR - RADIUS ** 2 ) )
     *              - DXTP ) / PTP
         STEPMX = MIN( STEPTR, STEPMX )
      END IF
C
C  IF THE CURVATURE IS POSITIVE, COMPUTE THE STEP TO THE MINIMIZER OF
C   THE QUADRATIC ALONG THE SEARCH DIRECTION, ALPHA.
C
      IF ( CURVAT .GT. ZERO ) THEN
         ALPHA = OLDGNS / CURVAT
         IF ( PRNTER ) WRITE( IOUT, 2030 ) ALPHA, STEPMX
C
C  ---------------------------------------------------------------------
C                LINE MINIMIZER ENCOUNTERED
C  ---------------------------------------------------------------------
C
C  IF THE MINIMIZER LIES BEFORE THE BOUNDARY, UPDATE THE SOLUTION 
C  AND PREPARE FOR ANOTHER CG STEP.
C
         IF ( ALPHA .LE. STEPMX ) THEN
C
C  UPDATE THE FUNCTION VALUE.
C
            FMODEL = FMODEL - HALF * ALPHA * OLDGNS
C
C  THE STEP LIES WITHIN THE BOUNDS. UPDATE THE SOLUTION.
C
CDIR$ IVDEP
            DO 430 J     = 1, NVAR
               I         = IVAR( J )
               XT( I )   = XT( I )   + ALPHA * P( I )
               GRAD( J ) = GRAD( J ) + ALPHA * Q( I )
  430       CONTINUE
C
C  IF NECESSARY, UPDATE THE SQUARE OF THE DISTANCE FROM X0 TO XT.
C
            IF (BOUNDX) DXSQR = DXSQR + ALPHA *
     *                          ( TWO * DXTP + ALPHA * PTP)
C
C  PERFORM ANOTHER ITERATION.
C
            IF ( ITER .LT. ITCGMX ) THEN
               IF ( ALPHA .GE. EPSMCH ) THEN
                  GO TO 200
               ELSE
C
C  THE STEP TAKEN IS VERY SMALL. FURTHER PROGRESS IS UNLIKELY. EXIT.
C
                  IF ( PRNTER ) WRITE( IOUT, 2010 ) ALPHA
                  IF ( PRONEL ) WRITE( IOUT, 2220 ) ITER, ALPHA
                  INFORM = 14
                  JUMPTO = 0
                  GO TO 500
               END IF
C
C  MORE THAN MAXIT ITERATIONS HAVE BEEN PERFORMED. EXIT.
C
            ELSE
               IF ( PRNTER ) WRITE( IOUT, 2040 )
               INFORM = 11
               JUMPTO = 0
               GO TO 500
            END IF
         END IF
C
C  IF NEGATIVE CURVATURE IS ENCOUNTERED, PREPARE TO EXIT.
C
      ELSE 
         INFORM = 13   
         IF ( PRNTER ) WRITE( IOUT, 2060 )
         IF ( PRONEL ) WRITE( IOUT, 2150 ) ITER
      END IF
C
C  ---------------------------------------------------------------------
C                BOUNDARY ENCOUNTERED
C  ---------------------------------------------------------------------
C
C  UPDATE THE FUNCTION VALUE.
C
      FMODEL = FMODEL + STEPMX * ( - OLDGNS + HALF * STEPMX * CURVAT )
C
C  IF A STEP IS MADE OUTSIDE THE SPHERICAL BOUNDARY, EXIT.
C
      IF ( BOUNDX .AND. STEPTR .LE. STEPBD ) THEN
C
C  UPDATE THE SOLUTION VALUES.
C
CDIR$ IVDEP
         DO 440 J   = 1, NVAR
            I       = IVAR( J )
            XT( I ) = MIN( BND( I, 2 ), MAX( BND( I, 1 ),
     *                     XT( I ) + STEPMX * P( I ) ) )
  440    CONTINUE
C
C  UPDATE THE SQUARE OF THE DISTANCE FROM X0 TO XT.
C
         DXSQR = DXSQR + STEPMX * ( TWO * DXTP + STEPMX * PTP )

C  SET THE EXIT CONDITIONS.
C
         INFORM = 12
         JUMPTO = 0
         IF ( PRNTER ) WRITE( IOUT, 2200 ) ITER
         GO TO 500
      END IF
C
C  A FREE VARIABLE ENCOUNTERS A BOUND, PREPARE FOR A RESTART.
C  COMPUTE WHICH VARIABLES ENCOUNTER THEIR BOUND.
C
      NEWVAR   = 0
CDIR$ IVDEP
      DO 450 J = 1, NVAR
         I     = IVAR( J )
C     
C  VARIABLE I LIES OFF ITS BOUND. UPDATE THE POINT.
C
         IF ( BREAKP( I ) .GE. STEPBD * ONEPEP ) THEN
            XT( I )   = XT( I ) + STEPBD * P( I )
C
C  SHIFT THE COMPONENTS OF IVAR AND GRAD AND UPDATE THE GRADIENT
C  TO ACCOUNT FOR THE REDUCED SET OF FREE VARIABLES.
C
            NEWVAR         = NEWVAR + 1
            IVAR( NEWVAR ) = I
            GRAD( NEWVAR ) = GRAD( J ) + STEPBD * Q( I )
         ELSE
C
C  VARIABLE I ENCOUNTERS ITS BOUND. FLAG THE VARIABLE.
C
            IVAR( J ) = - I
C
C  AN UPPER BOUND IS ENCOUNTERED.
C
            IF ( P( I ) .GT. EPSMCH ) THEN
C
C  MOVE TO THE POINT AT WHICH THE BOUNDARY IS ENCOUNTERED.
C
               XT( I ) = BND( I, 2 )
C
C  THE MODULUS OF IFIXED GIVES THE INDEX OF ONE OF THE VARIABLES
C  WHICH ENCOUNTERS A BOUND. IFIXED IS SET POSITIVE, TO INDICATE
C  THE VARIABLE HITS ITS UPPER BOUND.
C
               IFIXED = I
               IF ( PRNTER ) WRITE( IOUT, 2090 ) BND( I, 2 ), I
               IF ( PRONEL ) WRITE( IOUT, 2130 ) ITER, FMODEL, I
C
C  A LOWER BOUND IS ENCOUNTERED.
C
            ELSE
               XT( I ) = BND( I, 1 )
C
C  THE MODULUS OF IFIXED GIVES THE INDEX OF ONE OF THE VARIABLES
C  WHICH ENCOUNTERS A BOUND. IFIXED IS SET NEGATIVE, TO INDICATE
C  THE VARIABLE HITS ITS LOWER BOUND.
C
               IFIXED  = - I
               IF ( PRNTER ) WRITE( IOUT, 2100 ) BND( I, 1 ), I
               IF ( PRONEL ) WRITE( IOUT, 2140 ) ITER, FMODEL, I
            END IF
            P( I ) = ZERO
         END IF
  450 CONTINUE
      NVAR = NEWVAR
C
C  THERE ARE STILL SOME FREE VARIABLES.
C
      IF ( NVAR .GT. 0 ) THEN
         IF ( PRNTER ) WRITE( IOUT, 2070 ) ITER
         ITER   = 0
         IF ( INFORM .NE. 13 ) INFORM = 12
         JUMPTO = 5
      ELSE
C
C  IF THERE ARE NO FREE VARIABLES, EXIT.
C
         IF ( PRNTER ) WRITE( IOUT, 2050 ) ITER
         IF ( INFORM .NE. 13 ) INFORM = 12 
         JUMPTO = 0
      END IF
C
C  PREPARE TO EXIT.
C
  500 CONTINUE
      ITERCG = ITERCG + ITSLE
      ITSLE  = 0
      RETURN
C
C  NON-EXECUTABLE STATEMENTS.
C
 2000 FORMAT(/ ' --------------------  CG entered -------------------')
 2010 FORMAT(/ ' Step ', 1P, D12.4, ' too small in C.G.' )
 2020 FORMAT(/ ' CG iteration ',I3 ,
     *       / ' Function/gradient norm**2 are ',1P, 2D12.4,
     *         ' EPSGRD = ', 1P, D12.4 )
 2030 FORMAT(  ' Minimizing step and step to bound = ', 1P, 2D12.4 )
 2040 FORMAT(/ ' Maximum number of iterations exceeded in C.G.' )
 2050 FORMAT(/ ' Convergence criterion satisfied after ',
     *          I3,' C.G. iterations ' )
 2060 FORMAT(/ ' Direction of negative curvature ' )
 2070 FORMAT(/ ' Bound encountered in C.G. after ', I4, ' iterations ' )
 2080 FORMAT(/ ' Boundary of spherical region reached ', / )
 2090 FORMAT(/ ' Upper bound of ', 1P, D12.4, ' on variable ',
     *        I3, ' encountered ' )
 2100 FORMAT(/ ' Lower bound of ', 1P, D12.4, ' on variable ',
     *        I3, ' encountered ' )
 2110 FORMAT(/ '                                           ',
     *         '   Grad tol  =', 1P, D8.1,
     *       / '    **  CG entered **    Iteration  Model  ',
     *         '   Proj.Grad.  True Grad.  Step ' )
 2120 FORMAT( 25X, I7, 1P, 2D12.4, 12X, D12.4 )
 2130 FORMAT( 25X, I7, 1P, D12.4, '  -- Upper bound ', 
     *        I5, ' encountered' )
 2140 FORMAT( 25X, I7, 1P, D12.4, '  -- Lower bound ', 
     *        I5, ' encountered' )
 2150 FORMAT( 25X, I7, 13X, ' -- Negative curvature encountered ' )
 2160 FORMAT(/ ' CG iteration ',I3 ,
     *       /, ' Function/pr. gradient norm**2 are ', 1P, 2D12.4,
     *       /, ' Model gradient norm**2/EPSGRD are ', 1P, 2D12.4 )
 2170 FORMAT( 25X, I7, 1P, 4D12.4 )
C2180 FORMAT( ' I, XT, P, BREAKP = ', I5, 1P, 3D12.4 )
 2190 FORMAT(/ '    **  CG entered **    Iteration  Model  ',
     *       / '   Boundary of the spherical region reached ' )
 2200 FORMAT(/ ' Spherical bound encountered in C.G. after ', 
     *         I4, ' iterations ' )
 2210 FORMAT( ' VAR L X U G ', I6, 1P, 4D12.4 )
 2220 FORMAT( 25X, I7, 13X, ' -- Step ', 1P, D12.4, ' too small ' )
C
C  END OF SUBROUTINE CG.
C
      END

