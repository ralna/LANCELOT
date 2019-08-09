C ** Correction report.
C ** Correction 1. 13/01/94: 4 lines moved below **
C ** Correction 2. 13/01/94: 4 lines moved from above **
C ** Correction 3. 14/01/94: 1 line moved above **
C ** Correction 4. 14/01/94: 1 line moved from below **
C ** Correction 5. 09/08/96: 1 line replaced by 5 **
C ** Correction 6. 09/08/96: 1 line replaced by 5 **
C ** Correction 7. 30/05/97: 2 lines added **
C ** Correction 8. 30/05/97: 7 lines replaced by 11 **
C ** Correction 9. 30/05/97: 2 lines added **
C ** Correction 10. 30/05/97: 3 lines replaced by 8 **
C ** Correction 11. 30/05/97: 5 lines added **
C ** End of Correction report.
C  THIS VERSION: 30/05/1997 AT 09:15:00 AM.
CS    SUBROUTINE SAUGLG( N,  NG, NEL,    IELING, LELING, ISTADG, LSTADG,
CD    SUBROUTINE DAUGLG( N,  NG, NEL,    IELING, LELING, ISTADG, LSTADG,
     *                   IELVAR, LELVAR, ISTAEV, LSTAEV, INTVAR, LNTVAR,
     *                   ISTADH, LSTADH, ICNA  , LICNA , ISTADA, LSTADA,
     *                   A , LA, B , LB, BL    , LBL   , BU    , LBU   ,
     *                   GSCALE, LGSCAL, ESCALE, LESCAL, VSCALE, LVSCAL,
     *                   GXEQX,  LGXEQX, INTREP, LINTRE, KNDOFC, LKNDOF,
     *                   RANGES, INFORM, FOBJ  , X , LX, U , LU, GVALS , 
     *                   LGVALS, FT    , LFT   , FUVALS, LFUVAL, XT    , 
     *                   LXT   , ICALCF, LCALCF, NCALCF, ICALCG, LCALCG,
     *                   NCALCG, IVAR  , LIVAR , NVAR  , Q , LQ, DGRAD ,
     *                   LDGRAD, ICHOSE, ITER  , MAXIT , QUADRT, VNAMES,
     *                   LVNAME, GNAMES, LGNAME, STOPG , STOPC , IWK   ,
     *                   LIWK  , WK    , LWK   , IPRINT, IOUT  )
C
C  *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
C
C  AUGLG, A METHOD FOR FINDING A LOCAL MINIMIZER OF A FUNCTION
C  SUBJECT TO GENERAL CONSTRAINTS AND SIMPLE BOUNDS ON THE
C  SIZES OF THE VARIABLES. THE METHOD IS DESCRIBED IN THE PAPER
C  'A GLOBALLY CONVERGENT AUGMENTED LAGRANGIAN ALGORITHM FOR 
C  OPTIMIZATION WITH GENERAL CONSTRAINTS AND SIMPLE BOUNDS' BY
C  A.R.CONN, N.I.M.GOULD and PH.L.TOINT, 
C  SIAM J. NUM. ANAL. 28 (1991) PP.545-572.
C
C  THE OBJECTIVE FUNCTION IS ASSUMED TO BE OF THE FORM
C
C                             ISTADG(OBJ+1)-1
C      F( X ) = SUM GS   * G(   SUM    ES * F( X ) + A (TRANS) X - B )
C               OBJ   OBJ   OBJ  J=      J   J        M+1           M+1
C            IN OBJSET        ISTADG(OBJ)
C
C  AND THE CONSTRAINTS (I = 1, .... , NG, I .NE. OBJSET) OF THE FORM
C
C                        ISTADG(I+1)-1
C      CI( X ) = GS * G(   SUM     ES * F ( X ) + A (TRANS) X - B ) = 0
C                  I   I J=ISTADG(I) J   J         I             I
C
C  WHERE THE FJ( X ) ARE KNOWN AS NONLINEAR ELEMENT FUNCTIONS,
C  THE AI(TRANS) X + BI ARE THE LINEAR ELEMENT FUNCTIONS,
C  THE GSI ARE GROUP WEIGHTS AND THE ESI ARE ELEMENT WEIGHTS.
C  EACH FJ IS EXPECTED ONLY TO INVOLVE A FEW 'INTERNAL'
C  VARIABLES, THAT IS THERE IS A LINEAR TRANSFORMATION FROM THE
C  PROBLEM VARIABLES TO THE VARIABLES ACTUALLY NEEDED TO DEFINE
C  THE ELEMENT (AND ITS DERIVATIVES) WHOSE RANGE SPACE IS VERY SMALL.
C
C  CONTENTS OF THE ARRAY FUVALS:
C  -----------------------------
C
C     <- NEL -> <-- NINVAR --> <-- NHEL --> <---- N ---> <-
C
C    ---------------------------------------------------------
C    |  FJ(X)  |  GRAD FJ(X)  | HESS FJ(X) | GRAD F(X)) | ....
C    ---------------------------------------------------------
C   |         |              |            |            |
C  LFXI     LGXI           LHXI         LGGFX         LDX
C  (=0)
C        -> <------ N -----> <-- NVARGP ->
C
C       -----------------------------------
C      ... | DIAG HESS F(X) |  GRAD GI(X) |
C       -----------------------------------
C         |                |             |
C        LDX             LGFX           LEND
C
C  ONLY THE UPPER TRIANGULAR PART OF EACH ELEMENT HESSIAN IS STORED;
C  THE STORAGE IS BY COLUMNS.
C
C  CONTENTS OF THE ARRAY ISTADG:
C  -----------------------------
C
C         -> <------------------------ NEL -----------------------> <-
C
C           --------------------------------------------------------
C  PART OF  | ELEMENTS  | ELEMENTS  | ............... | ELEMENTS  | .
C  FUVALS   | GROUP 1   | GROUP 2   |                 | GROUP NG  |
C           ---------------------------------------------------------
C          | |           |                             |         | |
C       LFXI | |--- > ---|                             |       LGXI|
C            | |   |-------------------- > ------------|           |
C            | |   | |------------------------- > -----------------|
C            ---------
C  ISTADG:   | ..... |    POINTER TO THE POSITION OF THE 1ST ELEMENT
C            ---------    OF EACH CONSTRAINT AND THE OBJECTIVE
C                         FUNCTION WITHIN THE ARRAY.
C            <-NG+1->
C
C  CONTENTS OF THE ARRAYS IELVAR AND ISTAEV:
C  ----------------------------------------
C
C          <--------------------- NELVAR -------------------------->
C
C          ---------------------------------------------------------
C          | VARIABLES | VARIABLES | ............... |  VARIABLES  |
C  IELVAR: | ELEMENT 1 | ELEMENT 2 |                 | ELEMENT NEL |
C          ---------------------------------------------------------
C           |           |                             |             |
C           | |--- > ---|                             |             |
C           | |    |------------------- > ------------|             |
C           | |    | |----------------- > --------------------------|
C           ----------
C  ISTAEV:  | ...... |    POINTER TO THE POSITION OF THE 1ST VARIABLE
C           ----------    IN EACH ELEMENT (INCLUDING ONE TO THE END).
C
C          <- NEL+1 ->
C
C  CONTENTS OF THE ARRAY INTVAR:
C  -----------------------------
C
C  ON INITIAL ENTRY, INTVAR( I ), I = 1, ... , NEL, GIVES THE NUMBER OF
C  INTERNAL VARIABLES FOR ELEMENT I. THEREAFTER, INTVAR PROVIDES
C  POINTERS TO THE START OF EACH ELEMENT GRADIENT WITH RESPECT TO ITS
C  INTERNAL VARIABLES AS FOLLOWS:
C
C         -> <---------------------- NINVAR -----------------------> <-
C
C         -------------------------------------------------------------
C  PART OF  | GRADIENT  | GRADIENT  | ............... |  GRADIENT   | .
C  FUVALS   | ELEMENT 1 | ELEMENT 2 |                 | ELEMENT NEL |
C         -------------------------------------------------------------
C          | |           |                             |           | |
C       LGXI | |--- > ---|                             |         LHXI|
C            | |   |-------------------- > ------------|             |
C            | |   | |------------------------- > -------------------|
C            ---------
C  INTVAR:   | ..... |    POINTER TO THE POSITION OF THE 1ST ENTRY OF
C            ---------    THE GRADIENT FOR EACH ELEMENT.
C
C            <-NEL+1->
C
C  CONTENTS OF THE ARRAY ISTADH:
C  -----------------------------
C
C         -> <---------------------- NHEL -------------------------> <-
C
C         -------------------------------------------------------------
C  PART OF  | HESSIAN   | HESSIAN   | ............... |  HESSIAN    | .
C  FUVALS   | ELEMENT 1 | ELEMENT 2 |                 | ELEMENT NEL |
C         -------------------------------------------------------------
C          | |           |                             |           |
C       LHXI | |--- > ---|                             |         LGGFX
C            | |    |------------------- > ------------|
C            | |    |
C            ---------
C  ISTADH:   | ..... |    POINTER TO THE POSITION OF THE 1ST ENTRY OF
C            ---------    THE HESSIAN FOR EACH ELEMENT, WITH RESPECT
C                         TO ITS INTERNAL VARIABLES.
C            <- NEL ->
C
C  CONTENTS OF THE A, ARRAYS ICNA AND ISTADA:
C  ----------------------------------------
C
C          <--------------------- NA ----------------------------->
C
C          ---------------------------------------------------------
C          |   VALUES  |   VALUES  | ............... |    VALUES   |
C  A:      |    A(1)   |    A(2)   |                 |     A(NG)   |
C          ---------------------------------------------------------
C          | VARIABLES | VARIABLES | ............... |  VARIABLES  |
C  ICNA:   |    A(1)   |    A(2)   |                 |     A(NG)   |
C          ---------------------------------------------------------
C           |           |                             |             |
C           | |--- > ---|                             |             |
C           | |    |------------------- > ------------|             |
C           | |    | |----------------- > --------------------------|
C           ----------
C  ISTADA:  | ...... |    POINTER TO THE POSITION OF THE 1ST VARIABLE IN
C           ----------    THE LINEAR ELEMENT FOR EACH GROUP (INCLUDING
C                         ONE TO THE END).
C
C           <- NG+1 ->
C
C  IF THE ROUTINE TERMINATES WITH A NEGATIVE VALUE OF INFORM,
C  THE USER IS ASKED TO RE-ENTER THE SUBROUTINE WITH FURTHER
C  INFORMATION.
C
C  IF INFORM  = - 1, THE USER MUST SUPPLY THE FUNCTION
C                    AND DERIVATIVE VALUES OF EACH FJ AT THE POINT XT.
C  IF INFORM  = - 2, THE USER MUST SUPPLY THE FUNCTION AND DERIVATIVE
C                    OF EACH GI FOR THE ARGUMENT FT(I) .
C  IF INFORM  = - 3, THE USER MUST SUPPLY THE VALUE, ALONE, OF EACH
C                    FUNCTION FJ EVALUATED AT THE POINT XT.
C  IF INFORM  = - 4, THE USER MUST SUPPLY THE VALUE OF EACH FUNCTION
C                    GI, ALONE, FOR THE ARGUMENT FT(I).
C  IF INFORM  = - 5, THE USER MUST SUPPLY THE DERIVATIVES, ALONE, OF THE
C                    FUNCTIONS FJ AND GI AT THE POINT XT AND ARGUMENT
C                    FT(I) RESPECTIVELY.
C  IF INFORM  = - 6, THE USER MUST SUPPLY THE DERIVATIVES, ALONE, OF THE
C                    FUNCTIONS FJ AT THE POINT XT.
C  IF INFORM  = - 7, THE USER MUST SUPPLY THE VALUE, ALONE, OF EACH
C                    FUNCTION FJ EVALUATED AT THE POINT XT.
C  IF INFORM <= - 8, THE USER MUST PROVIDE THE PRODUCT OF THE
C                    INVERSE OF THE PRECONDITIONING MATRIX AND
C                    THE VECTOR D. THE NONZERO COMPONENTS OF
C                    GRAD OCCUR IN POSITIONS IVAR(I), I = 1,..,NVAR2
C                    AND HAVE THE VALUES DGRAD(I). THE PRODUCT MUST
C                    BE RETURNED IN THE VECTOR Q. THIS RETURN IS ONLY
C                    POSSIBLE IF THE LOGICAL ICHOSE( 2 ) IS 3
C
C  NICK GOULD, 19TH OF JANUARY 1990.
C  FOR CGT PRODUCTIONS.
C
C  *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
C
      INTEGER          N , NG, MAXIT , INFORM, ITER  , NVAR  , LWK   
      INTEGER          NEL   , LELVAR, LSTAEV, LSTADH, LNTVAR, LCALCF
      INTEGER          LELING, LINTRE, LFT   , LGXEQX, LSTADG, LGVALS
      INTEGER          LGSCAL, LESCAL, LVSCAL, LKNDOF, LCALCG, NCALCG
      INTEGER          LA, LB, LICNA , LSTADA, NOBJGR, NCALCF, IOUT
      INTEGER          IPRINT, LFUVAL, LIWK  , LU    
      INTEGER          LIVAR , LX    , LBL   , LBU   , LQ    , LDGRAD
      INTEGER          LGNAME, LXT   , LVNAME
CS    REAL             STOPG , STOPC , FOBJ
CD    DOUBLE PRECISION STOPG , STOPC , FOBJ
      INTEGER          IELVAR( LELVAR ), ISTAEV( LSTAEV )
      INTEGER          ISTADH( LSTADH ), IWK( LIWK ), IELING( LELING )
      INTEGER          INTVAR( LNTVAR ), ISTADG( LSTADG )
      INTEGER          ICNA  ( LICNA  ), ISTADA( LSTADA ), ICHOSE( 6 )
      INTEGER          ICALCF( LCALCF ), KNDOFC( LKNDOF )
      INTEGER          ICALCG( LCALCG ), IVAR  ( LIVAR  )
      LOGICAL          QUADRT, GXEQX( LGXEQX ), INTREP( LINTRE )
CS    REAL             X( LX ), FUVALS( LFUVAL ), BL( LBL ), BU( LBU ),
CD    DOUBLE PRECISION X( LX ), FUVALS( LFUVAL ), BL( LBL ), BU( LBU ),
     *                 A( LA ), B( LB ), U( LU ),
     *                 Q( LQ ), DGRAD( LDGRAD ),
     *                 XT( LXT ), FT( LFT ), GVALS( LGVALS, 3 ),
     *                 GSCALE( LGSCAL ), ESCALE( LESCAL ),
     *                 VSCALE( LVSCAL ), WK( LWK )
      CHARACTER * 10   VNAMES ( LVNAME ), GNAMES( LGNAME )
      EXTERNAL         RANGES
C
C  COMMON VARIABLES.
C
      INTEGER           ITERCG, ITCGMX, NGEVAL, ISKIP , IFIXED, NSEMIB
CS    REAL              ACCCG , RATIO , RADIUS, RADMAX, FINDMX, CMA31
CD    DOUBLE PRECISION  ACCCG , RATIO , RADIUS, RADMAX, FINDMX, CMA31
CS    COMMON / SCOMSB / ACCCG , RATIO , RADIUS, RADMAX, FINDMX, CMA31,
CD    COMMON / DCOMSB / ACCCG , RATIO , RADIUS, RADMAX, FINDMX, CMA31,
     *                  ITERCG, ITCGMX, NGEVAL, ISKIP , IFIXED, NSEMIB
CS    REAL              RMU   , RMUTOL, FIRSTC, FIRSTG
CD    DOUBLE PRECISION  RMU   , RMUTOL, FIRSTC, FIRSTG
      LOGICAL           NEWSOL
CS    COMMON / SCOMAL / RMU   , RMUTOL, FIRSTC, FIRSTG, NEWSOL
CD    COMMON / DCOMAL / RMU   , RMUTOL, FIRSTC, FIRSTG, NEWSOL
      INTEGER           IGETFD, LSPTRS, LNPTRS, LSELTS, LNELTS
      INTEGER           LSTAJC, LNSTJC, LGCOLJ, LNGCLJ
      LOGICAL           UNSUCC
      COMMON / INWKSP / IGETFD, LSPTRS, LNPTRS, LSELTS, LNELTS,
     *                  LSTAJC, LNSTJC, LGCOLJ, LNGCLJ, UNSUCC
C
C  LOCAL VARIABLES.
C
      INTEGER          I , IG, J , M , LGFX  , ICRIT , NCRIT
      INTEGER          IFIXD , IPDGEN, IDDGEN, ISTATE, IC
CS    REAL             ONE, F, ZERO  , HALF  , POINT1, POIN01, SMACHR,
CD    DOUBLE PRECISION ONE, F, ZERO  , HALF  , POINT1, POIN01, DMACHR,
     *                 CNORM , ETAK  , ETA0  , OMEGAK, OMEGA0, TAU   ,
     *                 GAMMA1, ALPHAE, BETAE , ALPHAK, RMUINV, YIUI  ,
     *                 RADSAV, WMIN  , ALPHAO, BETAO , OCNORM, EPSMCH,
     *                 OMEMIN, ETAMIN, HDASH , EPSTOL, EPSGRD, CTT   ,
     *                 SCALEG, EPSLAM, TEN   , FLOWER, PJGNRM, THETA
      LOGICAL          REEVAL, ITZERO 
      CHARACTER * 5    STATE( 5 )
C
C  FUNCTIONS.
C
      INTRINSIC        MIN, MAX, ABS, ANINT, LOG10
CS    EXTERNAL         SSBMIN, SMACHR
CD    EXTERNAL         DSBMIN, DMACHR
C ** Correction 1. 13/01/94: 4 lines moved below**
C ** Correction 1. 13/01/94:  end of correction **
C
C  SET CONSTANT PARAMETERS.
C
CS    PARAMETER ( ZERO   = 0.0E+0, HALF   = 5.0E-1, ONE    = 1.0E+0 )
CD    PARAMETER ( ZERO   = 0.0D+0, HALF   = 5.0D-1, ONE    = 1.0E+0 )
CS    PARAMETER ( POINT1 = 1.0E-1, POIN01 = 1.0E-2, TEN    = 1.0E+1 )
CD    PARAMETER ( POINT1 = 1.0D-1, POIN01 = 1.0D-2, TEN    = 1.0D+1 )
CS    PARAMETER ( WMIN   = 1.0E-1, THETA  = 1.0E-1                  )
CD    PARAMETER ( WMIN   = 1.0D-1, THETA  = 1.0D-1                  )
C ** Correction 4. 14/01/94: 1 line moved from below **
      SAVE
C ** Correction 4. 14/01/94: end of correction **
C ** Correction 2. 13/01/94: 4 lines moved from above **
C
C  SET DATA VALUES.
C
      DATA STATE  / ' FREE', 'LOWER', 'UPPER', 'FIXED', 'DEGEN' /
C ** Correction 2. 13/01/94: end of correction **
C
C  FACTORS WHICH INFLUENCE THE INTERMEDIATE CONVERGENCE TOLERANCES.
C
CS    DATA GAMMA1, TAU           / 1.0E-1, 1.0E-1 /
CD    DATA GAMMA1, TAU           / 1.0D-1, 1.0D-1 /
CS    DATA ALPHAE, BETAE         / 1.0E-1, 9.0E-1 /
CD    DATA ALPHAE, BETAE         / 1.0D-1, 9.0D-1 /
CS    DATA ALPHAO, BETAO         / 1.0E+0, 1.0E+0 /
CD    DATA ALPHAO, BETAO         / 1.0D+0, 1.0D+0 /
C ** Correction 3. 14/01/94: 1 line moved above **
C ** Correction 3. 14/01/94: end of correction **
C
C  BRANCH TO THE INTERIOR OF THE CODE IF A RE-ENTRY IS BEING MADE.
C
      IF ( INFORM .LT. 0 ) GO TO 310
C
C  SET INITIAL INTEGER VALUES.
C
      NOBJGR   = 0
      M        = 0
      DO 10 IG = 1, NG
         IF ( KNDOFC( IG ) .GE. 2 ) THEN
            IF ( KNDOFC( IG ) .GT. 4 ) THEN
               INFORM = 7
               RETURN
            ELSE
               M = M + 1
            END IF
         ELSE
            NOBJGR = NOBJGR + 1
         END IF
   10 CONTINUE
      ICRIT  = 0
      NCRIT  = 9
      ITZERO = ITER .EQ. 0
C
C  SET INITIAL REAL VALUES.
C
CS    EPSMCH = SMACHR( 1 )
CD    EPSMCH = DMACHR( 1 )
      EPSTOL = EPSMCH ** 0.75
      OMEMIN = STOPG 
      ETAMIN = STOPC 
      EPSGRD = STOPG
      EPSLAM = STOPC
      OMEGA0 = FIRSTG / ( MIN( RMU, GAMMA1 ) ** ALPHAO )
      ETA0   = FIRSTC / ( MIN( RMU, GAMMA1 ) ** ALPHAE )
CS    FLOWER = - SMACHR( 5 )
CD    FLOWER = - DMACHR( 5 )
C
C  SET INITIAL VALUES FOR THE LAGRANGE MULTIPLIER ESTIMATES, U, INTERNAL 
C  GROUP SCALINGS, WK,  AND THE ARRAY, GXEQX,  WHICH TELLS IF EACH 
C  GROUP IS TRIVIAL.
C
      IF ( NG .GT. 0 ) THEN
         DO 20 IG            = 1, NG
C           U( IG )          = ZERO
            IF ( KNDOFC( IG ) .GT. 1 ) THEN
               WK( IG )         = ONE
               GXEQX( NG + IG ) = .FALSE.
            ELSE 
               WK( IG )         = GSCALE( IG )
               GXEQX( NG + IG ) = GXEQX( IG )
            END IF
   20    CONTINUE
      END IF
      IF ( IOUT .GT. 0 .AND. IPRINT .GT. 0 ) WRITE( IOUT, 2000 )
C
C  SET THE CONVERGENCE TOLERENCES.
C
      ALPHAK = MIN( RMU, GAMMA1 )
      ETAK   = MAX( ETAMIN, ETA0   * ( ALPHAK ** ALPHAE ) )
      OMEGAK = MAX( OMEMIN, OMEGA0 * ( ALPHAK ** ALPHAO ) )
      ACCCG  = POIN01
      RADSAV = RADIUS
      IF ( IOUT .GT. 0 .AND. IPRINT .GT. 0 )
     *     WRITE( IOUT, 2050 ) RMU, OMEGAK, ETAK
      PJGNRM = OMEGAK
C
C  START OF THE OUTER ITERATION.
C
  100 CONTINUE
C
C  START INNER LOOP TO MINIMIZE PENALTY FUNCTION FOR THE
C  CURRENT VALUES OF RMU AND U.
C
  200 CONTINUE
CS    CALL SSBMIN( N,  NG, NEL   , IELING, LELING, ISTADG, LSTADG,
CD    CALL DSBMIN( N,  NG, NEL   , IELING, LELING, ISTADG, LSTADG,
     *             IELVAR, LELVAR, ISTAEV, LSTAEV, INTVAR, LNTVAR,
     *             ISTADH, LSTADH, ICNA  , LICNA , ISTADA, LSTADA,
     *             A     , LA    , B     , LB    , BL    , LBL   , 
     *             BU    , LBU   , WK( 1), NG    , ESCALE, LESCAL, 
     *             VSCALE, LVSCAL, GXEQX( NG + 1),
     *             NG    , INTREP, LINTRE, 
     *             RANGES, INFORM, F , X , LX    , GVALS , LGVALS, 
     *             FT    , LFT   , FUVALS, LFUVAL, XT    , LXT   , 
     *             ICALCF, LCALCF, NCALCF, ICALCG, LCALCG, NCALCG, 
     *             IVAR  , LIVAR , NVAR  , Q , LQ, DGRAD , LDGRAD, 
     *             ICHOSE, ITER  , MAXIT , QUADRT, PJGNRM, FLOWER,
     *             IWK   , LIWK  , 
     *             WK( 2 * NG + 1 ), LWK - 2 * NG , IPRINT, IOUT )
C
C  RETURN TO THE USER TO OBTAIN PROBLEM DEPENDENT INFORMATION.
C
  300 CONTINUE
      IF ( INFORM .LE. - 1 .AND. INFORM .GE. - 10 ) RETURN
C
C  RE-ENTRY POINT.
C
  310 CONTINUE
      IF ( INFORM .LT. 0 ) THEN
C
C  CALCULATE PROBLEM RELATED INFORMATION.
C
         IF ( INFORM .EQ. - 1 ) THEN
C
C  IF THERE ARE SLACK VARIABLES, INITIALIZE  THEM TO MINIMIZE THE
C  INFEASIBILITY OF THEIR ASSOCIATED CONSTRAINTS.
C
            IF ( ITZERO ) THEN
               DO 340 IG = 1, NG
                  IF ( KNDOFC( IG ) .GE. 3 ) THEN
C
C  CALCULATE THE CONSTRAINT VALUE FOR THE INEQUALITY CONSTRAINTS.
C  IT IS ASSUMED THAT THE SLACK VARIABLE OCCURS LAST IN THE LIST 
C  OF VARIABLES IN THE LINEAR ELEMENT.
C
                     CTT   = - B( IG )
C
C  INCLUDE THE CONTRIBUTION FROM THE LINEAR ELEMENT.
C
                     DO 320 J = ISTADA( IG ), ISTADA( IG + 1 ) - 2
                        CTT  = CTT + A( J ) * X( ICNA( J ) )
  320                CONTINUE
C
C  INCLUDE THE CONTRIBUTIONS FROM THE NONLINEAR ELEMENTS.
C
                     DO 330 J = ISTADG( IG ), ISTADG( IG + 1 ) - 1
                        CTT  = CTT + ESCALE( J ) *
     *                               FUVALS( IELING( J ) )
  330                CONTINUE
                     J  = ISTADA( IG + 1 ) - 1
                     IC = ICNA( J )
C
C  THE SLACK VARIABLE CORRESPONDS TO A LESS-THAN-OR-EQUAL-TO
C  CONSTRAINT. SET ITS VALUE AS CLOSE AS POSSIBLE TO THE 
C  CONSTRAINT VALUE.
C
                     IF ( KNDOFC( IG ) .EQ. 3 ) THEN
                        X( IC ) = MIN( MAX( BL( IC ), - CTT ),
     *                                 BU( IC ) )

C
C  THE SLACK VARIABLE CORRESPONDS TO A GREATER-THAN-OR-EQUAL-TO
C  CONSTRAINT. SET ITS VALUE AS CLOSE AS POSSIBLE TO THE 
C  CONSTRAINT VALUE.
C
                     ELSE
                        X( IC ) = MIN( MAX( BL( IC ),   CTT ),
     *                                 BU( IC ) )
                     END IF
C 
C  COMPUTE A SUITABLE SCALE FACTOR FOR THE SLACK VARIABLE.
C
                     IF ( X( IC ) .GT. ONE ) THEN
                        VSCALE( IC ) = TEN ** ANINT( LOG10( X( IC ) ) )
                     ELSE
                        VSCALE( IC ) = ONE
                     END IF
                  END IF
  340          CONTINUE
               ITZERO = .FALSE.
            END IF
         END IF
         IF ( INFORM .EQ. - 2 .OR. INFORM .EQ. - 4 ) THEN
C
C  RECORD THE UNSCALED CONSTRAINT VALUES IN WK( NG + IG ).
C
            DO 350 IG = 1, NG
               IF ( GXEQX( IG ) ) THEN
                  WK( NG + IG ) = FT( IG )
               ELSE
                  WK( NG + IG ) = GVALS( IG, 1 )
               END IF
  350       CONTINUE
            IF ( M .GT. 0 ) THEN
C
C  PRINT THE CONSTRAINT VALUES ON THE FIRST ITERATION.
C
               IF ( ITER .EQ. 0 ) THEN
                  IF ( IOUT .GT. 0 .AND. IPRINT .GE. 1 ) THEN
                     J     = 1
                     FOBJ  = ZERO
                     CNORM = ZERO
                     IF ( IPRINT .GE. 3 ) WRITE( IOUT, 2120 )
                     DO 351 I = 1, NG
                        IF ( KNDOFC( I ) .EQ. 1 ) THEN
                           IF ( I - 1 .GE. J .AND. IPRINT .GE. 3 )
     *                       WRITE( IOUT, 2090 ) ( GNAMES( IG ), IG, 
     *                          WK( NG + IG ) * GSCALE( IG ),
     *                              IG = J, I - 1 )
                           J    = I + 1
                           FOBJ = FOBJ + WK( NG + I ) * GSCALE( I )
                        END IF
 351                 CONTINUE   
                     IF ( NG .GE. J .AND. IPRINT .GE. 3 )
     *                 WRITE( IOUT, 2090 ) ( GNAMES( IG ), IG, 
     *                       WK( NG + IG ) * GSCALE( IG ), IG = J, NG )
C
C  PRINT THE OBJECTIVE FUNCTION VALUE ON THE FIRST ITERATION.
C
                     IF ( NOBJGR .GT. 0 ) THEN
                        WRITE( IOUT, 2010 ) FOBJ * FINDMX
                     ELSE
                        WRITE( IOUT, 2020 ) 
                     END IF   
                  END IF
C
C  CALCULATE THE CONSTRAINT NORM.
C
                  CNORM     = ZERO
                  DO 360 IG = 1, NG
                     IF ( KNDOFC( IG ) .NE. 1 ) CNORM = MAX( CNORM,
     *                    ABS( GSCALE( IG ) * WK( NG + IG ) ) )
  360             CONTINUE
                  IF ( IOUT .GT. 0 .AND. IPRINT .EQ. 1 ) THEN
                     IF ( ITER .EQ. 0 ) THEN
                        WRITE( IOUT, 2200 ) CNORM
                     ELSE
                        WRITE( IOUT, 2180 ) RMU, OMEGAK, CNORM, ETAK
                     END IF
                  END IF
               END IF
C
C  CALCULATE THE TERMS INVOLVING THE CONSTRAINTS FOR THE
C  AUGMENTED LAGRANGIAN FUNCTION.
C
               RMUINV    = HALF / RMU
               DO 370 IG = 1, NG
                  IF ( KNDOFC( IG ) .NE. 1 ) THEN
                     YIUI           = GSCALE( IG ) * WK( NG + IG )
     *                                + RMU * U( IG )
                     GVALS( IG, 1 ) = ( RMUINV * YIUI ) * YIUI
                  END IF
  370          CONTINUE
            END IF
         END IF
         IF ( INFORM .EQ. - 2 .OR. INFORM .EQ. - 5 ) THEN
            IF ( M .GT. 0 ) THEN
C
C  CALCULATE THE DERIVATIVES OF THE TERMS INVOLVING THE CONSTRAINTS
C  FOR THE AUGMENTED LAGRANGIAN FUNCTION.
C
               DO 380 IG = 1, NG
                  IF ( KNDOFC( IG ) .NE. 1 ) THEN
                     SCALEG = GSCALE( IG )
                     IF ( GXEQX( IG ) ) THEN
                        GVALS( IG, 3 ) = SCALEG * ( SCALEG / RMU )
                        GVALS( IG, 2 ) = SCALEG * ( FT( IG ) *
     *                                 ( SCALEG / RMU ) + U( IG ) )
                     ELSE
                        HDASH          = SCALEG * ( U( IG ) +
     *                                   WK( NG + IG ) * ( SCALEG
     *                                   / RMU ) )
                        GVALS( IG, 3 ) = HDASH * GVALS( IG, 3 ) + 
     *                                   ( ( SCALEG * GVALS( IG, 2 )
     *                                     ) ** 2 ) / RMU 
                        GVALS( IG, 2 ) = HDASH * GVALS( IG, 2 )
                     END IF
                  END IF
  380          CONTINUE
            END IF
         END IF
         GO TO 200
      ELSE
C
C  TEST WHETHER THE MAXIMUM ALLOWED NUMBER OF ITERATIONS HAS
C  BEEN REACHED.
C
         IF ( INFORM .GT. 3 ) RETURN
C
C  PRINT THE VALUES OF THE CONSTRAINT FUNCTIONS.
C
         IF ( IOUT .GT. 0 ) THEN
            J    = 1
            FOBJ = ZERO
            IF ( IPRINT .GE. 3 ) WRITE( IOUT, 2120 )
            DO 381 I = 1, NG
               IF ( KNDOFC( I ) .EQ. 1 ) THEN
                  IF ( I - 1 .GE. J .AND. IPRINT .GE. 3 )
     *               WRITE( IOUT, 2090 ) ( GNAMES( IG ), IG, 
     *                   WK( NG + IG ) * GSCALE( IG ),
     *                       IG = J, I - 1 )
                  J    = I + 1
C ** Correction 6. 09/08/96: 1 line replaced by 5 **
                  IF ( GXEQX( I ) ) THEN
                     FOBJ = FOBJ + FT( I ) * GSCALE( I )
                  ELSE
                     FOBJ = FOBJ + GVALS( I, 1 ) * GSCALE( I )
                  END IF                     
C ** Correction 6. 09/08/96: end of correction 5 **
               END IF
 381        CONTINUE   
            IF ( NG .GE. J .AND. IPRINT .GE. 3 ) 
     *           WRITE( IOUT, 2090 ) ( GNAMES( IG ), IG, 
     *                 WK( NG + IG ) * GSCALE( IG ), IG = J, NG )
C
C  PRINT THE OBJECTIVE FUNCTION VALUE.
C
            IF ( IPRINT .GE. 1 ) THEN
               IF ( NOBJGR .GT. 0 ) THEN
                  WRITE( IOUT, 2010 ) FOBJ * FINDMX
               ELSE
                  WRITE( IOUT, 2020 ) 
               END IF   
            END IF   
         END IF
C
C  CALCULATE THE CONSTRAINT NORM.
C
         OCNORM    = CNORM
         CNORM     = ZERO
         DO 390 IG = 1, NG
            IF ( KNDOFC( IG ) .NE. 1 ) CNORM = MAX( CNORM, 
     *           ABS( WK( NG + IG ) * GSCALE(IG ) ) )
  390    CONTINUE
         IF ( INFORM .EQ. 1 ) GO TO 500
         IF ( IOUT .GT. 0 .AND. IPRINT .GT. 0 )
     *        WRITE( IOUT, 2060 ) RMU, PJGNRM, OMEGAK, CNORM, ETAK
C
C  TEST FOR CONVERGENCE OF THE OUTER ITERATION.
C
         IF ( ( OMEGAK .LE. STOPG .OR. PJGNRM .LE. STOPG ) .AND.
     *          CNORM  .LE. STOPC ) GO TO 500
C
C  COMPUTE THE RATIO OF SUCCESSIVE NORMS OF CONSTRAINT VIOLATIONS.
C  IF THIS RATIO IS NOT SUBSTANTIALLY DECREASED OVER NCRIT ITERATIONS,
C  EXIT WITH THE WARNING THAT NO FEASIBLE POINT CAN BE FOUND.
C
         IF ( CNORM .GT. 9.9D-1 * OCNORM ) THEN
            ICRIT = ICRIT + 1
            IF ( ICRIT .GE. NCRIT ) THEN
               INFORM = 8
               IF ( IOUT .GT. 0 .AND. IPRINT .GT. 0 )
     *              WRITE( IOUT, 2160 ) NCRIT
              GO TO 500
            END IF
         ELSE
            ICRIT = 0
         END IF
C
C  RECORD THAT AN APPROXIMATE MINIMIZER OF THE AUGMENTED LAGRANGIAN 
C  FUNCTION HAS BEEN FOUND.
C
         NEWSOL = .TRUE.
         IF ( IOUT .GT. 0 .AND. IPRINT .GE. 3 )
     *       WRITE( IOUT, 2070 ) ( X( I ), I = 1, N )
C
C  ANOTHER ITERATION WILL BE PERFORMED.
C
         INFORM = - 1
C
C  CHECK TO SEE IF THE CONSTRAINT HAS BEEN SUFFICIENTLY REDUCED.
C   
         IF ( CNORM .LT. ETAK .AND. RMU .LE. RMUTOL ) THEN
            IF ( OCNORM .GT. 1.0D-10 .AND. IOUT .GT. 0 .AND.
     *           IPRINT .GE. 3 )
     *         WRITE( IOUT, 2080 ) CNORM / OCNORM, ALPHAK ** BETAE
C
C  THE CONSTRAINT NORM HAS BEEN REDUCED SUFFICIENTLY.
C  UPDATE THE LAGRANGE MULTIPLIER ESTIMATES, U.
C
            IF ( M .GT. 0 ) THEN
               DO 420 IG = 1, NG
                  IF ( KNDOFC( IG ) .NE. 1 ) U( IG ) = U( IG ) + 
     *                 WK( NG + IG ) * ( GSCALE( IG ) / RMU )
  420          CONTINUE
            END IF
            IF ( IOUT .GT. 0 .AND. IPRINT .GE. 1 ) WRITE( IOUT, 2040)
            IF ( IOUT .GT. 0 .AND. IPRINT .GE. 3 ) THEN
               J = 1
               DO 430 I = 1, NG
                  IF ( KNDOFC( I ) .EQ. 1 ) THEN
                     IF ( I - 1 .GE. J ) WRITE( IOUT, 2030 ) 
     *                  ( U( IG ), IG = J, I - 1 )
                     J = I + 1
                  END IF
 430          CONTINUE   
              IF ( NG .GE. J ) WRITE( IOUT, 2030 ) 
     *                      ( U( IG ), IG = J, NG )
            END IF
C
C  DECREASE THE CONVERGENCE TOLERANCES.
C
            ALPHAK = MIN( RMU, GAMMA1 )
            ETAK   = MAX( ETAMIN, ETAK   * ( ALPHAK ** BETAE ) )
            OMEGAK = MAX( OMEMIN, OMEGAK * ( ALPHAK ** BETAO ) )
C
C  PREPARE FOR THE NEXT OUTER ITERATION.
C
            IF ( IOUT .GT. 0 .AND. IPRINT .GT. 0 )
     *           WRITE( IOUT, 2050 ) RMU, OMEGAK, ETAK
            PJGNRM = OMEGAK
            RADIUS = RADSAV
C
C  MOVE VARIABLES WHICH ARE CLOSE TO THEIR BOUNDS ONTO THE BOUND.
C
            LGFX   = ISTADH( NEL + 1 ) - 1
            REEVAL = .FALSE.
            DO 440 I = 1, N
               XT( I ) = X( I )
               IF ( X( I ) .NE. BL( I ) .AND. X( I ) - BL( I ) .LE. 
     *              THETA * FUVALS( LGFX + I ) ) THEN
                  REEVAL  = .TRUE.
                  XT( I ) = BL( I )
               END IF
               IF ( X( I ) .NE. BU( I ) .AND. X( I ) - BU( I ) .GE. 
     *              THETA * FUVALS( LGFX + I ) ) THEN
                  REEVAL  = .TRUE.
                  XT( I ) = BU( I )
               END IF
 440        CONTINUE
         ELSE
C
C  REDUCE THE PENALTY PARAMETER AND RESET THE CONVERGENCE TOLERANCES.
C
            RMU    = TAU * RMU
            ALPHAK = MIN( RMU, GAMMA1 )
            ETAK   = MAX( ETAMIN, ETA0   * ( ALPHAK ** ALPHAE ) )
            OMEGAK = MAX( OMEMIN, OMEGA0 * ( ALPHAK ** ALPHAO ) )
            ACCCG  = POIN01
C
C  PREPARE FOR THE NEXT OUTER ITERATION.
C
            IF ( IOUT .GT. 0 .AND. IPRINT .GT. 0 ) WRITE( IOUT, 2150)
            IF ( IOUT .GT. 0 .AND. IPRINT .GT. 0 )
     *           WRITE( IOUT, 2050 ) RMU, OMEGAK, ETAK
            PJGNRM = OMEGAK
            RADIUS = RADSAV
C
C  MOVE VARIABLES WHICH ARE CLOSE TO THEIR BOUNDS ONTO THE BOUND.
C
            REEVAL = .FALSE.
            LGFX   = ISTADH( NEL + 1 ) - 1
            DO 450 I = 1, N
               XT( I ) = X( I )
               IF ( X( I ) .NE. BL( I ) .AND. X( I ) - BL( I ) .LE. 
     *              THETA * FUVALS( LGFX + I ) ) THEN
                  REEVAL  = .TRUE.
                  XT( I ) = BL( I )
               END IF
               IF ( X( I ) .NE. BU( I ) .AND. X( I ) - BU( I ) .GE. 
     *              THETA * FUVALS( LGFX + I ) ) THEN
                  REEVAL  = .TRUE.
                  XT( I ) = BU( I )
               END IF
 450        CONTINUE
C
C  IF FINITE-DIFFERENCE GRADIENTS ARE USED, USE CENTRAL DIFFERENCES
C  WHENEVER THE PENALTY PARAMETER IS SMALL.
C
            IF ( ICHOSE( 3 ) .GE. 1 .AND. RMU .LT. EPSMCH ** 0.25 )
     *         ICHOSE( 3 ) = 2 
         END IF
C
C  SEE IF WE NEED TO RE-EVALUATE THE PROBLEM FUNCTIONS.
C
         IF ( REEVAL ) THEN
            IF ( IOUT .GT. 0 .AND. IPRINT .GT. 0 ) WRITE( IOUT, 2190)
            NGEVAL = NGEVAL + 1
            IF ( ICHOSE( 3 ) .EQ. 1 ) IGETFD = 0
CS          CALL SCLCFG( UNSUCC, N, NCALCF, NCALCG, ISTAEV, LSTAEV,
CD          CALL DCLCFG( UNSUCC, N, NCALCF, NCALCG, ISTAEV, LSTAEV,
     *                   ISTADG, LSTADG, IELING, LELING,
     *                   ICALCF, LCALCF, ICALCG, LCALCG,
     *                   IWK( LSPTRS + 1 ), LNPTRS,
     *                   IWK( LSELTS + 1 ), LNELTS,
     *                   IWK( LSTAJC + 1 ), LNSTJC,
     *                   IWK( LGCOLJ + 1 ), LNGCLJ, X, XT )
            DO 460 I  = 1, N
               X( I ) = XT( I )
  460       CONTINUE   
            GO TO 300
         END IF
         GO TO 100
      END IF
C
C  END OF THE MAIN ITERATION.
C
  500 CONTINUE
C
C  RECORD THE NORMS OF THE CONSTRAINTS AND PROJECTED GRADIENT.
C
      STOPG = PJGNRM
      STOPC = CNORM
C
C  COMPUTE THE FINAL, UNWEIGHTED, MULTIPLIER ESTIMATES.
C
      IPDGEN = 0
      IF ( INFORM .EQ. 0 .AND. M .GT. 0 ) THEN
         DO 510 IG  = 1, NG
            IF ( KNDOFC( IG ) .NE. 1 ) THEN
               SCALEG  = GSCALE( IG )
               U( IG ) = SCALEG * ( U( IG ) +
     *                   WK( NG + IG ) * ( SCALEG / RMU ) )
               IF ( ABS( U( IG ) ) .LE. EPSLAM ) IPDGEN = IPDGEN + 1
            END IF
  510    CONTINUE
      END IF
      IF ( IOUT .GT. 0 .AND. IPRINT .GT. 0 ) THEN
C     IF ( IOUT .GT. 0 .AND. IPRINT .GT. 1 ) THEN
         LGFX   = ISTADH( NEL + 1 ) - 1
         IFIXD  = 0
         IDDGEN = 0
         WRITE( IOUT, 2100 )
C ** Correction 7. 30/05/97: 2 lines added **
         J = N
C        J = 2
C ** Correction 7. 30/05/97: end of correction 7 **
         DO 520 I  = 1, N
            ISTATE = 1
            IF ( X( I ) .GE. BU( I ) - EPSTOL *
     *           MAX( ONE, ABS( BU( I ) ) ) ) ISTATE = 3
            IF ( X( I ) .LE. BL( I ) + EPSTOL *
     *           MAX( ONE, ABS( BL( I ) ) ) ) ISTATE = 2
            IF ( BU( I ) - BL( I ) .LE. 2.0D+0 * EPSMCH ) ISTATE = 4
            IF ( MAXIT .GE. 0 ) THEN
              IF ( ( ISTATE .EQ. 2 .OR. ISTATE .EQ. 3 ) .AND.
     *               ABS( FUVALS( LGFX + I ) ) .LE. EPSGRD ) ISTATE = 5
               IF ( ISTATE .EQ. 5 ) IDDGEN = IDDGEN + 1
            END IF
            IF ( ISTATE .GT. 1 ) IFIXD  = IFIXD + 1
C ** Correction 8. 30/05/97: 7 lines replaced by 11 **
            IF ( I .LE. J .OR. I .GE. N + 1 - J ) THEN
               IF ( MAXIT  .GE. 0 ) THEN
                  WRITE( IOUT, 2110 ) VNAMES( I ), I, STATE( ISTATE ), 
     *                     X( I ), BL( I ), BU( I ), FUVALS( LGFX + I )
               ELSE
                  WRITE( IOUT, 2210 ) VNAMES( I ), I, STATE( ISTATE ), 
     *                     X( I ), BL( I ), BU( I )
               END IF
            ELSE
               IF ( I .EQ. J + 1 ) WRITE( IOUT, 2220 ) 
            END IF
C ** Correction 8. 30/05/97: end of correction 8 **
  520    CONTINUE
      END IF
C
C  COMPUTE THE OBJECTIVE FUNCTION VALUE.
C 
      FOBJ = ZERO
      IF ( IOUT .GT. 0 .AND. IPRINT .GT. 0 ) WRITE( IOUT, 2130 )
C ** Correction 9. 30/05/97: 2 lines added **
      J = N
C     J = 2
C ** Correction 9. 30/05/97: end of correction 9 **
      DO 530 IG = 1, NG
         IF ( KNDOFC( IG ) .NE. 1 ) THEN 
C ** Correction 10. 30/05/97: 3 lines replaced by 8 **
            IF ( IOUT .GT. 0 .AND. IPRINT .GT. 0 ) THEN
               IF ( IG .LE. J + 1 .OR. IG .GE. NG + 1 - J ) THEN
                    WRITE( IOUT, 2140 ) GNAMES( IG ), IG,
     *                  WK( NG + IG ), GSCALE( IG ), U( IG )
               ELSE
                  IF ( IG .EQ. J + 2 ) WRITE( IOUT, 2230 ) 
               END IF
            END IF
C ** Correction 10. 30/05/97: end of correction 10 **
         ELSE
C ** Correction 5. 09/08/96: 1 line replaced by 5 **
            IF ( GXEQX( IG ) ) THEN
               FOBJ = FOBJ + FT( IG ) * GSCALE( IG )
            ELSE
               FOBJ = FOBJ + GVALS( IG, 1 ) * GSCALE( IG )
            END IF
C ** Correction 5. 09/08/96: end of correction 5 **
         END IF
  530 CONTINUE
      IF ( IOUT .GT. 0 .AND. IPRINT .GT. 0 ) THEN
         IF ( NOBJGR .GT. 0 ) THEN
            WRITE( IOUT, 2010 ) FOBJ * FINDMX
         ELSE
            WRITE( IOUT, 2020 ) 
         END IF   
         WRITE( IOUT, 2170 ) N, M, IPDGEN, IFIXD, IDDGEN
      END IF
      RETURN
C
C  NON EXECUTABLE STATEMENTS.
C
 2000 FORMAT( /, ' *********  Starting optimization  ************** ' )
 2010 FORMAT( /, ' Objective function value  ', 1P, D22.14 )
 2020 FORMAT( /, ' There is no objective function ' )
 2030 FORMAT( /, ' Multiplier values ', /, ( 1P, 5D12.4 ) )
 2040 FORMAT( /, ' ******** Updating multiplier estimates ********** ' )
 2050 FORMAT( /, ' Penalty parameter ', 1P, D12.4,
     *           ' Required projected gradient norm = ', 1P, D12.4, /,
     *           '                   ', 12X,
     *           ' Required constraint         norm = ', 1P, D12.4 )
 2060 FORMAT( /, ' Penalty parameter       = ', 1P, D12.4, /,
     *           ' Projected gradient norm = ', 1P, D12.4,
     *           ' Required gradient   norm = ', 1P, D12.4, /,
     *           ' Constraint         norm = ', 1P, D12.4, 
     *           ' Required constraint norm = ', 1P, D12.4 )
 2070 FORMAT( /, ' Solution   values ', /, ( 1P, 5D12.4 ) )
 2080 FORMAT( /, ' ||c|| / ||c(old)|| = ', 1P, D12.4,
     *           ' vs ALPHA ** BETAE = ', D12.4 )
 2090 FORMAT( ( 4X, A10, I7, 6X, 1P, D22.14 ) )
 2100 FORMAT( /, ' Variable name Number Status     Value',
     *           '   Lower bound Upper bound  |  Dual value ',
     *        /, ' ------------- ------ ------     -----',
     *           '   ----------- -----------  |  ----------' )
 2110 FORMAT( 2X, A10, I7, 4X, A5, 1P, 3D12.4, '  |', D12.4 )
 2120 FORMAT( /, ' Constraint name Number        Value ' )
 2130 FORMAT( /, ' Constraint name Number    Value    Scale factor ',
     *           '| Lagrange multiplier',
     *        /, ' --------------- ------    -----    ----- ------ ',
     *           '| -------------------')
 2140 FORMAT( 4X, A10, I7, 2X, 1P, 2D12.4, '  |   ', D12.4 )
 2150 FORMAT( /, ' ***********    Reducing mu    *************** ' )
 2160 FORMAT( /, ' Constraint violations have not decreased',
     *           ' substantially over ', I4, ' major iterations. ',
     *        /, ' Problem possibly infeasible, terminating run. ' )
 2170 FORMAT( /, ' There are ', I7, ' variables in total. ',
     *        /, ' There are ', I7, ' equality constraints. ',
     *        /, ' Of these  ', I7, ' are primal degenerate. ',
     *        /, ' There are ', I7, ' variables on their bounds. ',
     *        /, ' Of these  ', I7, ' are dual degenerate. ' )
 2180 FORMAT( /, ' Penalty parameter       = ', 1P, D12.4, /,
     *           '                           ', 12X, 
     *           ' Required gradient norm   = ', 1P, D12.4, /,
     *           ' Constraint norm         = ', 1P, D12.4, 
     *           ' Required constraint norm = ', 1P, D12.4 )
 2190 FORMAT( /, ' Using the shifted starting point. ' )
 2200 FORMAT(    ' Constraint norm           ', 1P, D22.14 )
 2210 FORMAT( 2X, A10, I7, 4X, A5, 1P, 3D12.4, '  |      - ' )
C ** Correction 11. 30/05/97: 5 lines added **
 2220 FORMAT( '  .               .    .....',
     *        ' ........... ........... ...........',
     *        '  | ...........' )
 2230 FORMAT( '    .               .   ........... ...........',
     *        '  |    ........... ' )
C ** Correction 11. 30/05/97: end of correction 11 **
C
C  END OF AUGLG.
C
      END
