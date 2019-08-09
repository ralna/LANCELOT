C ** Correction report.
C ** Correction 1. 22/07/92: 5 lines added **
C ** Correction 2. 22/07/92: 3 lines added **
C ** Correction 3. 29/01/93: 2 lines corrected **
C ** Correction 4. 05/02/93: 1 line corrected **
C ** Correction 5. 05/02/93: 1 line corrected **
C ** Correction 6. 07/07/94: 2 lines added for VAXs **
C ** Correction 7. 06/12/94: 1 line corrected **
C ** Correction 8. 19/01/95: 1 line corrected **
C ** Correction 9. 03/12/98: 2 lines corrected **
C ** Correction 10. 03/12/98: 4 lines corrected **
C ** Correction 11. 03/12/98: 7 lines corrected **
C ** Correction 12. 03/12/98: 2 lines changed to 3 **
C ** End of Correction report.
C  THIS VERSION: 03/12/98 08:10:00 GMT 1998
C  VERSION OF 03/03/88 JACOBIAN STORED BY COLUMNS.
C  VERSION OF 10/03/89 DIRECT METHOD ALLOWED FOR LINEAR SYSTEMS.
C  VERSION OF 15/03/89 DIRECT METHOD PRECONDITIONER ALLOWED.
C  VERSION OF 28/03/89 MUNKSGAARD'S PRECONDITIONER ALLOWED.
C  VERSION OF 05/05/89 JACOBIAN STORED BY BOTH ROWS AND COLUMNS.
C  VERSION OF 22/11/89 LIST RATHER THAN FLAG CALCULATION OF ELEMENTS.
C  VERSION OF 19/01/90 ALLOW SCALING OF GROUPS, VARIABLES AND ELEMENTS.
C  VERSION OF 22/08/90 ALLOW ACCURATE SOLUTION OF MODEL PROBLEM.
C  VERSION OF 16/10/91 ALLOW FINITE-DIFFERENCE GRADIENTS.
C
C  ** FOR THE CRAY 2, LINES STARTING 'CDIR$ IVDEP' TELL THE COMPILER TO
C     IGNORE POTENTIAL VECTOR DEPENDENCIES AS THEY ARE KNOWN TO BE O.K.
C
CS    SUBROUTINE SSBMIN( N,  NG, NEL   , IELING, LELING, ISTADG, LSTADG,
CD    SUBROUTINE DSBMIN( N,  NG, NEL   , IELING, LELING, ISTADG, LSTADG,
     *                   IELVAR, LELVAR, ISTAEV, LSTAEV, INTVAR, LNTVAR,
     *                   ISTADH, LSTADH, ICNA  , LICNA,  ISTADA, LSTADA,
     *                   A , LA, B , LB, BL    , LBL   , BU    , LBU   , 
     *                   GSCALE, LGSCAL, ESCALE, LESCAL, VSCALE, LVSCAL, 
     *                   GXEQX , LGXEQX, INTREP, LINTRE, RANGES, INFORM, 
     *                   F , X , LX    , GVALS , LGVALS, FT    , LFT   , 
     *                   FUVALS, LFUVAL, XT    , LXT   , ICALCF, LCALCF, 
     *                   NCALCF, ICALCG, LCALCG, NCALCG, IVAR  , LIVAR ,
     *                   NVAR  , Q , LQ, DGRAD , LDGRAD, ICHOSE, ITER  , 
     *                   MAXIT , QUADRT, STOPG , FLOWER, IWK   , LIWK  ,
     *                   WK    , LWK   , IPRINT, IOUT  )
C
C  *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
C
C  SBMIN, A METHOD FOR FINDING A LOCAL MINIMIZER OF A FUNCTION
C  SUBJECT TO RESTRICTIONS ON THE SIZES OF THE VARIABLES.
C  THE ALGORITHM IS DESCRIBED IN THE PAPER 'TESTING A CLASS OF
C  METHODS FOR SOLVING MINIMIZATION PROBLEMS WITH SIMPLE BOUNDS
C  ON THEIR VARIABLES' BY A.R.CONN, N.I.M.GOULD AND PH.L.TOINT,
C  MATHEMATICS OF COMPUTATION, 50 (1988) PP.399-430.
C
C  NICK GOULD, 16TH OF OCTOBER 1991.
C  FOR CGT PRODUCTIONS.
C
C  ------ THE GROUP PARTIALLY SEPARABLE VERSION -------
C
C  THE OBJECTIVE FUNCTION IS ASSUMED TO BE OF THE FORM
C
C               NG          
C     F( X ) = SUM  GS * G(   SUM     ES  * F ( X ) + A (TRANS) X - B )
C              I=1    I   I J IN J(I)   IJ   J         I             I
C                                 
C
C  WHERE THE F(J)( X ) ARE KNOWN AS NONLINEAR ELEMENT FUNCTIONS,
C  THE A(I)(TRANS) X - BI ARE THE LINEAR ELEMENT FUNCTIONS, THE GS(I)
C  ARE GROUP WEIGHTS, THE ES(IJ) ARE ELEMENT WEIGHTS, THE G(I)
C  ARE CALLED GROUP FUNCTIONS AND EACH SET J(I) IS A SUBSET OF THE SET
C  OF THE FIRST NEL INTEGERS. EACH F(J) IS EXPECTED ONLY TO INVOLVE A
C  FEW 'INTERNAL' VARIABLES, THAT IS THERE IS A LINEAR TRANSFORMATION
C  FROM THE PROBLEM VARIABLES TO THE VARIABLES ACTUALLY NEEDED TO
C  DEFINE THE ELEMENT (AND ITS DERIVATIVES) WHOSE RANGE SPACE IS VERY
C  SMALL. EACH GROUP FUNCTION IS A DIFFERENTIABLE UNIVARIATE FUNCTION.
C
C  CONTENTS OF THE ARRAY FUVALS:
C  -----------------------------
C
C     <-NEL-><-- NINVAR --> <-- NHEL --> <---- N ---> <-
C
C    ---------------------------------------------------------
C    |  FJ(X)  |  GRAD FJ(X)  | HESS FJ(X) | GRAD F(X)) | ....
C    ---------------------------------------------------------
C   |         |              |            |            |
C  LFXI     LGXI           LHXI         LGGFX         LDX
C  (=0)
C        -> <------ N -----> <-- NVARGP ->
C
C       --------------------------------------
C      ... | DIAG SCALING F(X) |  GRAD GI(X) |
C       --------------------------------------
C         |                   |             |
C        LDX                LGRJAC         LEND
C
C  ONLY THE UPPER TRIANGULAR PART OF EACH ELEMENT HESSIAN IS STORED;
C  THE STORAGE IS BY COLUMNS.
C
C  CONTENTS OF THE ARRAYS ISTADG, ESCALE AND IELING:
C  ------------------------------------------------
C
C           <---------------------- NGEL -------------------------->
C
C           --------------------------------------------------------
C  ESCALE:  | EL.SCALES | EL.SCALES | ............... | EL.SCALES  |
C           | GROUP 1   | GROUP 2   |                 | GROUP NG   |
C           --------------------------------------------------------
C  IELING:  | ELEMENTS  | ELEMENTS  | ............... | ELEMENTS   |
C           | GROUP 1   | GROUP 2   |                 | GROUP NG   |
C           --------------------------------------------------------
C            |           |                             |            |
C            | |--- > ---|                             |            |
C            | |   |-------------------- > ------------|            |
C            | |   | |------------------------- > ------------------|
C            ---------
C  ISTADG:   | ..... |    POINTER TO THE POSITION OF THE 1ST ELEMENT
C            ---------    OF EACH GROUP WITHIN THE ARRAY.
C
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
C  *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
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
C  IF INFORM  = - 8, 9, 10, THE USER MUST PROVIDE THE PRODUCT OF THE
C                    INVERSE OF THE PRECONDITIONING MATRIX AND
C                    THE VECTOR D. THE NONZERO COMPONENTS OF
C                    GRAD OCCUR IN POSITIONS IVAR(I), I = 1,..,NVAR2
C                    AND HAVE THE VALUES DGRAD(I). THE PRODUCT MUST
C                    BE RETURNED IN THE VECTOR Q. THIS RETURN IS ONLY
C                    POSSIBLE IF ICHOSE( 2 ) IS 3.
C
C  IF THE USER DOES NOT WISH TO COMPUTE AN ELEMENT OR GROUP FUNCTION AT
C  A PARTICULAR ARGUMENT RETURNED FROM SBMIN, SHE MAY RESET INFORM TO
C  -11 AND RE-ENTER. SBMIN WILL TREAT SUCH A RE-ENTRY AS IF THE CURRENT
C  ITERATION HAD BEEN UNSUCCESSFUL AND REDUCE THE TRUST REGION. THIS 
C  FACILITY IS USEFUL WHEN, FOR INSTANCE, THE USER IS ASKED TO EVALUATE
C  A FUNCTION AT A POINT OUTSIDE ITS DOMAIN OF DEFINITION.
C
C  *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
C
      INTEGER           N,  NG, NEL   , MAXIT , INFORM, ITER  , NVAR
      INTEGER           LFUVAL, LELVAR, LSTAEV, LSTADH, LELING, LNTVAR
      INTEGER           LCALCF, NCALCF, LINTRE, LFT   , LGXEQX, LSTADG
      INTEGER           LGSCAL, LESCAL, LVSCAL, LCALCG, NCALCG, IPRINT
      INTEGER           LIVAR , LX    , LBL   , LBU   , LQ    , LDGRAD
      INTEGER           LXT   , LGVALS, LA, LB, LICNA , LSTADA, IOUT  
      INTEGER           LWK   , LIWK
CS    REAL              STOPG , FLOWER, F
CD    DOUBLE PRECISION  STOPG , FLOWER, F
      LOGICAL           QUADRT
      INTEGER           IVAR  ( LIVAR  ), IELVAR( LELVAR ), ICHOSE( 6 )
      INTEGER           ISTAEV( LSTAEV ), ICALCF( LCALCF )
      INTEGER           ISTADH( LSTADH ), IWK   ( LIWK   )
      INTEGER           INTVAR( LNTVAR ), ISTADG( LSTADG ) 
      INTEGER           ICNA  ( LICNA  ), ISTADA( LSTADA )
      INTEGER           ICALCG( LCALCG ), IELING( LELING )
      LOGICAL           GXEQX ( LGXEQX ), INTREP( LINTRE )
CS    REAL              X     ( LX     ), FUVALS( LFUVAL ), 
CD    DOUBLE PRECISION  X     ( LX     ), FUVALS( LFUVAL ), 
     *                  BL    ( LBL    ), BU    ( LBU    ), 
     *                  XT    ( LXT    ), A     ( LA     ), Q    ( LQ ),
     *                  B     ( LB     ), DGRAD ( LDGRAD ), 
     *                  FT    ( LFT    ), GVALS ( LGVALS  , 3 ),
     *                  GSCALE( LGSCAL ), ESCALE( LESCAL ),
     *                  VSCALE( LVSCAL ), WK    ( LWK    )
      EXTERNAL          RANGES
C
C  COMMON VARIABLES.
C
      INTEGER           ITERCG, ITCGMX, NGEVAL, ISKIP , IFIXED, NSEMIB
CS    REAL              ACCCG , RATIO , RADIUS, RADMAX, FINDMX, CMA31
CD    DOUBLE PRECISION  ACCCG , RATIO , RADIUS, RADMAX, FINDMX, CMA31
CS    COMMON / SCOMSB / ACCCG , RATIO , RADIUS, RADMAX, FINDMX, CMA31,
CD    COMMON / DCOMSB / ACCCG , RATIO , RADIUS, RADMAX, FINDMX, CMA31,
     *                  ITERCG, ITCGMX, NGEVAL, ISKIP , IFIXED, NSEMIB
CS    REAL              EPSMCH, EPSNEG, TINY  , BIG
CD    DOUBLE PRECISION  EPSMCH, EPSNEG, TINY  , BIG
CS    COMMON / SMACHN / EPSMCH, EPSNEG, TINY  , BIG
CD    COMMON / DMACHN / EPSMCH, EPSNEG, TINY  , BIG
      INTEGER           IGETFD, LSPTRS, LNPTRS, LSELTS, LNELTS
      INTEGER           LSTAJC, LNSTJC, LGCOLJ, LNGCLJ
      LOGICAL           UNSUCC
      COMMON / INWKSP / IGETFD, LSPTRS, LNPTRS, LSELTS, LNELTS,
     *                  LSTAJC, LNSTJC, LGCOLJ, LNGCLJ, UNSUCC
C
C  LOCAL VARIABLES.
C
      INTEGER           I,  J , LX0   , IEL   , IS, LL, LWKSTR, IFACTR
      INTEGER           JUMPTO, K,  KK, LDX   , LFXI  , LGRJAC, LGXI 
      INTEGER           NBPROD, NVAR1 , K1, K2, NVAR2 , LNWK  , LNWKB
      INTEGER           MAXSEL, NINVAR, NEL1  , NFREEC, LNWKC , IBQPST
      INTEGER           INFOR , NUMBER, LIWK2 , LGGFX , NG1   , NFIXED
      INTEGER           LINDEX, LSVGRP, LSWKSP, LWK2  , LP, IG, LBND
      INTEGER           NFREEF, NN    , NFREE , LNGUVL, LNGRJC, LDELTX
      INTEGER           LHXI  , LVALJR, LSYMMD, LSYMMH, MAXSIN, NGEL
      INTEGER           LNNDEX, LNWKSP, LNSTGV, LSTAGV, LSEND , LSLGRP
      INTEGER           LNIUSE, LNNNON, LNNNO2, LNSYMD, LNSYMH, LNLGRP
      INTEGER           LNVGRP, LNVLJR, LWKB  , LWKC  , LXCP  , LGX0
      INTEGER           LNQGRD, LNBRAK, LNP   , LNBND , LNFREC, LFREEC
      INTEGER           LNFXI , LNGXI , LNHXI , LNHUVL, LNGGFX, LNDX
      INTEGER           NNONNZ, LIUSED, LNNONZ, LNONZ2, LBREAK, LQGRAD
      INTEGER           LSTYPE, LSSWTR, LSSIWT, LSIWTR, LSWTRA
      INTEGER           LNTYPE, LNSWTR, LNSIWT, LNIWTR, LNWTRA, NTYPE
      INTEGER           LNISET, LNSVSE, LSISET, LSSVSE, NSETS 
C     INTEGER           ICORAS
      INTEGER           ISYS( 5 )
CS    REAL              PJGNRM, EPSTLP, GMODEL, EPSGRD, FTT   , XI, GI, 
CD    DOUBLE PRECISION  PJGNRM, EPSTLP, GMODEL, EPSGRD, FTT   , XI, GI,
     *                  FNEW  , RADMIN, CGSTOP, ONE   , DIAMIN, DIAMAX,
     *                  ARED  , PRERED, RHO   , FMODEL, GAMMA0, CURV  ,
     *                  BLI   , BUI   , DXSQR , HALF  , GAMMA1, STEPMX, 
     *                  GAMMA2, ETA   , RMU   , POINT1, SMALLH, RESMIN,
     *                  GNORM , QGNORM, OLDRAD, DISTAN, EPSCNS, SLOPE ,
     *                  RADTOL, FILL  , STEP  , STPTOL, ZERO  , TENEPS,
C    *                  ASTEP ,
     *                  STPMIN, EPSTLN, VSCMAX, DLTNRM, F0    , TENM2 ,
     *                  TEN   , HUNDRD, TENTEN, TENM4
      LOGICAL           PRCOND, TWONRM, DIRECT, MYPREC, FIRSUP, FIRSTC
      LOGICAL           ALLLIN, DIAG  , ALTRIV, NEXT  , SECOND, PRINT1
      LOGICAL           MODCHL, DPRCND, IPRCND, MUNKS , SEPREC, DENSEP
      LOGICAL           CALCDI, XACTCP, REUSEC, GMPSPR, SLVBQP
      LOGICAL           REFACT, BAND  , FDGRAD, CENTRL
C     LOGICAL           SAMEAS, INTERA
C     CHARACTER * 1     INCREA
      CHARACTER * 2     LS    , LSYS( 5 )
      CHARACTER * 6     CGEND , LISEND, CGENDS( 5 ), LSENDS( 5 )
      CHARACTER * 7     ATIME , CHRTIM
C
C  EXTERNAL SUBROUTINES AND FUNCTIONS USED.
C
CS    EXTERNAL          SCG   , SCAUCH, SELGRD, SHSPRD, SFRNTL, SGTPGR  
CS    EXTERNAL          SPRECN, SPRGRA, SSECNT, SCLCFG                  
CS    EXTERNAL          SINXAC, SNRM2 , SINITW, SDOT  , SMACHR, CPUTIM  
CS    EXTERNAL          SSCALH, SFDGRD, CHRTIM
CD    EXTERNAL          DCG   , DCAUCH, DELGRD, DHSPRD, DFRNTL, DGTPGR
CD    EXTERNAL          DPRECN, DPRGRA, DSECNT, DCLCFG
CD    EXTERNAL          DINXAC, DNRM2 , DINITW, DDOT  , DMACHR, CPUTIM
CD    EXTERNAL          DSCALH, DFDGRD, CHRTIM
C
C  ** FOR TESTING PURPOSES ONLY, REMOVE ON GENERAL RELEASE **
C
C     INTEGER           ISTATE( 10000 )
C
C  INTRINSIC FUNCTIONS.
C
      INTRINSIC         ABS,    MIN,    MAX,    SQRT
C
C  BASIC LINEAR ALGEBRA AND MACHINE FUNCTIONS.
C
CS    REAL              SNRM2,  SDOT,   SMACHR, SDNRM
CD    DOUBLE PRECISION  DNRM2,  DDOT,   DMACHR, DDNRM
CS    EXTERNAL          SCOPY,  SSETVL, SSETVI
CD    EXTERNAL          DCOPY,  DSETVL, DSETVI
C ** Correction 6. 07/07/94: 2 lines added for VAXs **
CSVMS EXTERNAL         SDATSB
CDVMS EXTERNAL         DDATSB
C ** Correction 6. 07/07/94: end of correction **
C
C  TIMING PARAMETERS. TCA, TLS, TMV AND TUP ARE, RESPECTIVELY, THE 
C  TIMES SPENT IN FINDING THE CAUCHY STEP, FINDING THE APPROXIMATE
C  MINIMIZER OF THE MODEL, FORMING THE PRODUCT OF THE HESSIAN WITH 
C  A SPECIFIED VECTOR AND IN UPDATING THE SECOND DERIVATIVE APPROX-
C  -IMATIONS. TIME GIVES THE CLOCK TIME ON INITIAL ENTRY TO SBMIN. 
C  T GIVES THE INSTANTANEOUS CLOCK TIME AND DUM IS A DUMMY ARGUMENT.
C
      REAL              DUM,    T, TUP, TIME,   TCA,    TLS,    TMV
      SAVE
C
C  TIMING FUNCTION.
C
      REAL             CPUTIM
C
C  SET CONSTANT PARAMETERS.
C
CS    PARAMETER ( ZERO   = 0.0E+0, HALF   = 5.0E-1, ONE    = 1.0E+0  )
C ** Correction 7. 06/12/94: 1 line corrected **
CD    PARAMETER ( ZERO   = 0.0D+0, HALF   = 5.0D-1, ONE    = 1.0D+0  )
C ** Correction 7. 06/12/94: end of correction **
CS    PARAMETER ( POINT1 = 1.0E-1, GAMMA1 = 5.0E-1, GAMMA2 = 2.0E+0  )
CD    PARAMETER ( POINT1 = 1.0D-1, GAMMA1 = 5.0D-1, GAMMA2 = 2.0D+0  )
CS    PARAMETER ( RMU    = 2.5E-1, ETA    = 7.5E-1, GAMMA0 = 6.25E-2 )
CD    PARAMETER ( RMU    = 2.5D-1, ETA    = 7.5D-1, GAMMA0 = 6.25D-2 )
C ** Correction 8. 19/01/95: 1 line corrected **
CS    PARAMETER ( STPTOL = 1.0E-1, TENM2  = 1.0E-2, TENTEN = 1.0E+10 )
C ** Correction 8. 19/01/95: end of correction **
CD    PARAMETER ( STPTOL = 1.0D-1, TENM2  = 1.0D-2, TENTEN = 1.0D+10 )
CS    PARAMETER ( TEN    = 1.0E+1, HUNDRD = 1.0E+2, TENM4  = 1.0E-4  )
CD    PARAMETER ( TEN    = 1.0D+1, HUNDRD = 1.0D+2, TENM4  = 1.0D-4  )
C
C  INITIALIZE DATA.
C
      DATA        ISYS   / 0, 0, 0, 0, 0 /
      DATA        CGENDS / ' CONVR', ' MAXIT', ' BOUND', ' -CURV',
     *                     ' S<EPS' /
      DATA        LSENDS / ' PSDEF', ' INDEF', ' SINGC', ' SINGI',
     *                     ' PRTRB' /
      DATA        LSYS   / 'PD', 'IN', 'SC', 'SI', 'PP' /
      DATA        TCA, TLS, TMV, TUP / 0.0E+0, 0.0E+0, 0.0E+0, 0.0E+0 /
C
C  BRANCH TO DIFFERENT PARTS OF THE CODE DEPENDING ON THE
C  INPUT VALUE OF INFORM.
C
      I = - INFORM
      IF ( INFORM .LT. 0 )
     *   GO TO ( 40, 70, 430, 470, 540, 540, 540, 90, 320, 570, 
     *           590 ), I
C     INTERA = .FALSE.
C
C  TEST WHETHER THE MAXIMUM ALLOWED NUMBER OF ITERATIONS HAS
C  BEEN REACHED.
C
      IF ( MAXIT .LT. 0 ) THEN
         IF ( IPRINT .GE. 0 ) WRITE( IOUT, 2090 )
         INFORM = 1
         GO TO 700
      END IF
C
C  RECORD THE NUMBER OF NONLINEAR ELEMENTS, NGEL. 
C
      NG1  = NG            + 1
      NEL1 = NEL           + 1
      NGEL = ISTADG( NG1 ) - 1
      IF ( IPRINT .GE. 1 ) WRITE( IOUT, 2020 ) N, NG, NEL
C
C  CHECK THE LENGTHS OF INPUT ARRAYS. EXIT IF THERE IS INSUFFICIENT
C  ROOM IN ANY ARRAY.
C
C  INTEGER ARRAYS.
C
      INFORM = 4
      IF ( LNTVAR .LT. NEL1 ) THEN
         WRITE( IOUT, 2400 ) NEL1, LNTVAR
         RETURN
      END IF
      IF ( LSTADH .LT. NEL1 ) THEN
         WRITE( IOUT, 2410 ) NEL1, LSTADH
         RETURN
      END IF
      IF ( LSTADG .LT. NG1 ) THEN
         WRITE( IOUT, 2450 ) NG1, LSTADG
         RETURN
      END IF
      IF ( LSTAEV .LT. NEL1 ) THEN
         WRITE( IOUT, 2460 ) NEL1, LSTAEV
         RETURN
      END IF
      IF ( LSTADA .LT. NG1 ) THEN
         WRITE( IOUT, 2470 ) NG1, LSTADA
         RETURN
      END IF
      IF ( LICNA .LT. ISTADA( NG1 ) - 1 ) THEN
         WRITE( IOUT, 2290 ) ISTADA( NG1 ) - 1, LICNA
         RETURN
      END IF
      IF ( LCALCF .LT. NEL ) THEN
         WRITE( IOUT, 2520 ) NEL, LCALCF
         RETURN
      END IF
      IF ( LCALCG .LT. NG ) THEN
         WRITE( IOUT, 2420 ) NG, LCALCG
         RETURN
      END IF
C
C  REAL ARRAYS.
C
      INFORM = 5
      IF ( LGVALS .LT. NG ) THEN
         WRITE( IOUT, 2500 ) NG, LGVALS
         RETURN
      END IF
      IF ( LFT .LT. NG ) THEN
         WRITE( IOUT, 2510 ) NG, LFT
         RETURN
      END IF
      IF ( LA .LT. ISTADA( NG1 ) - 1 ) THEN
         WRITE( IOUT, 2140 ) ISTADA( NG1 ) - 1, LA
         RETURN
      END IF
      IF ( LB .LT. NG ) THEN
         WRITE( IOUT, 2690 ) NG, LB
         RETURN
      END IF
C
C  LOGICAL ARRAYS.
C
      INFORM = 6
      IF ( LGXEQX .LT. NG ) THEN
         WRITE( IOUT, 2480 ) NG, LGXEQX
         RETURN
      END IF
      IF ( LINTRE .LT. NEL ) THEN
         WRITE( IOUT, 2490 ) NEL, LINTRE
         RETURN
      END IF
C
C  INITIALIZE INTEGER PARAMETERS.
C
C  ISKIP GIVES THE TOTAL NUMBER OF SECANT UPDATES THAT ARE SKIPPED.
C  ITERCG IS THE NUMBER OF CG ITERATIONS PERFORMED. NGEVAL IS THE 
C  NUMBER OF DERIVATIVE EVALUATIONS MADE. NUMBER IS USED TO CONTROL
C  WHICH NEGATIVE EIGENVALUE IS PICKED WHEN THE NEGATIVE CURVATURE
C  EXPLOITING MULTIFRONTAL SCHEME IS USED.
C
      ITER   = 0
      ISKIP  = 0
      ITERCG = 0
      NGEVAL = 0
      NUMBER = 0
C
C  INITIALIZE FLOATING-POINT PARAMETERS.
C
C  EPSMCH IS THE (POSITIVE) MACHINE PRECISION, EPSNEG IS THE (NEGATIVE)
C  MACHINE PRECISION. TINY AND BIG ARE THE SMALLEST AND LARGEST 
C  POSITIVE FLOATING POINT NUMBERS. EPSTLP AND EPSTLN ARE
C  TOLERANCES ON HOW FAR A VARIABLE MAY LIE AWAY FROM ITS
C  BOUND AND STILL BE CONSIDERED ACTIVE. RADTOL IS THE SMALLEST VALUE
C  THAT THE TRUST REGION RADIUS IS ALLOWED AND STPMIN IS THE SMALLEST
C  ALLOWABLE STEP BETWEEN CONSECUTIVE ITERATES. TENEPS IS 10 TIMES THE
C  MACHINE PRECISION(!) AND EPSGRD IS THE SMALLEST THAT THE CG 
C  RESIDUALS WILL BE REQUIRED TO BE. VSCMAX IS THE LARGEST SPECIFIED
C  VARIABLE SCALING AND FILL IS THE AMOUNT OF FILL-IN WHEN A DIRECT
C  METHOD IS USED TO FIND AN APPROXIMATE SOLUTION TO THE MODEL PROBLEM.
C
CS    EPSMCH = SMACHR( 1 )
CD    EPSMCH = DMACHR( 1 )
C ** Correction 9. 03/12/98: 2 lines corrected **
CS    EPSNEG = SMACHR( 1 )
CD    EPSNEG = DMACHR( 1 )
C ** Correction 9. End of correction **
C ** Correction 3. 29/01/93: 2 lines corrected **
CS    TINY   = SMACHR( 4 )
CD    TINY   = DMACHR( 4 )
C ** Correction 3. 29/01/93: end of correction **
CS    BIG    = SMACHR( 5 )
CD    BIG    = DMACHR( 5 )
      EPSTLP = EPSMCH
      EPSTLN = EPSNEG
      EPSGRD = HUNDRD  * EPSMCH ** 2
      TENEPS = TEN     * EPSMCH
      RADTOL = POINT1  * EPSMCH
      SMALLH = EPSMCH ** 0.3333
C     STPMIN = EPSMCH
      STPMIN = EPSMCH ** 0.75
      VSCMAX = ZERO
      FILL   = ZERO
      TIME   = CPUTIM( DUM )
C
C  INITIALIZE LOGICAL PARAMETERS.
C
C  FIRSUP IS .TRUE. IF INITIAL SECOND DERIVATIVE APPROXIMATIONS ARE
C  TO BE RESCALED USING THE SHANNO-PHUA SCALING.
C
      FIRSUP = .FALSE.
      NEXT   = .FALSE.
C
C  ALLLIN IS .TRUE. IF THERE ARE NO NONLINEAR ELEMENTS AND .FALSE.
C  OTHERWISE.
C
      ALLLIN = NEL .EQ. 0
C
C  TWONRM IS .TRUE. IF THE TWO-NORM TRUST REGION IS TO BE USED, AND
C  IS .FALSE. IF THE INFINITY-NORM TRUST REGION IS REQUIRED.
C
      TWONRM = ICHOSE( 1 ) .EQ. 1
C
C  DIRECT IS .TRUE. IF THE LINEAR SYSTEM IS TO BE SOLVED USING A
C  DIRECT METHOD (MA27). OTHERWISE, THE LINEAR SYSTEM WILL BE SOLVED
C  USING CONJUGATE GRADIENTS.
C
      DIRECT = ICHOSE( 2 ) .GE. 11
C
C  MODCHL IS .TRUE. IF THE HESSIAN IS TO BE FORCED TO BE POSITIVE DEF.
C  PRCOND IS .TRUE. IF PRECONDITIONING IS TO BE USED IN THE CONJUGATE
C  GRADIENT ITERATION.
C
      MODCHL = ICHOSE( 2 ) .EQ. 12
      PRCOND = .NOT. DIRECT .AND. ICHOSE( 2 ) .GE. 2
C
C  DPRCND IS .TRUE. IF THE USER WISHES TO USE A DIAGONAL PRECONDITIONER.
C
      DPRCND = ICHOSE( 2 ) .EQ. 2
      CALCDI = ICHOSE( 2 ) .EQ. 2
C
C  MYPREC IS .TRUE. IF THE USER IS TO TAKE RESPONSIBILITY FOR
C  PROVIDING THE PRECONDITIONER.
C
      MYPREC = ICHOSE( 2 ) .EQ. 3
C
C  IPRCND IS .TRUE. IF THE USER WISHES TO USE A POSITIVE DEFINITE
C  PERTURBATION OF THE INNER BAND OF THE TRUE MATRIX AS A PRECONDITIONER
C
      IPRCND = ICHOSE( 2 ) .EQ. 4
C
C  MUNKS IS .TRUE. IF THE MUNKSGAARDS PRECONDITIONER IS TO BE USED.
C
      MUNKS = ICHOSE( 2 ) .EQ. 5
C
C  SEPREC IS .TRUE. IF THE USER WISHES TO USE THE SCHNABEL-ESKOW
C  POSITIVE DEFINITE PERTURBATION OF THE COMPLETE MATRIX AS 
C  A PRECONDITIONER.
C
      SEPREC = ICHOSE( 2 ) .EQ. 6
C
C  GMPSPR IS .TRUE. IF THE USER WISHES TO USE THE GILL-MURRAY-PONCELEON-
C  SAUNDERS POSITIVE DEFINITE PERTURBATION OF THE COMPLETE MATRIX AS 
C  A PRECONDITIONER.
C
      GMPSPR = ICHOSE( 2 ) .EQ. 7
C
C  BAND IS .TRUE. IF THE USER WISHES TO USE A BANDSOLVER AS
C  A PRECONDITIONER.
C
      BAND   = ICHOSE( 2 ) .EQ. 8
      DIAG   = PRCOND .AND. .NOT. MYPREC
C
C  FDGRAD IS .FALSE. IF THE USER PROVIDES EXACT FIRST DERIVATIVES
C  OF THE NONLINEAR ELEMENT FUNCTIONS AND .TRUE. OTHERWISE.
C
      FDGRAD = ICHOSE( 3 ) .GE. 1
C
C  SECOND IS .TRUE. IF THE USER PROVIDES EXACT SECOND DERIVATIVES
C  OF THE NONLINEAR ELEMENT FUNCTIONS AND .FALSE. OTHERWISE.
C
      SECOND = ICHOSE( 4 ) .EQ. 0
C
C  XACTCP IS .TRUE, IF THE USER WISHES TO CALCULATE THE EXACT CAUCHY
C  POINT IN THE FASHION OF CONN, GOULD AND TOINT (1988). IF AN
C  APPROXIMATION SUFFICES, XACTCP WILL BE .FALSE.
C
      XACTCP = ICHOSE( 5 ) .EQ. 1
C
C  SLVBQP IS .TRUE. IF A GOOD APPROXIMATION TO THE MINIMUM OF THE 
C  QUADRATIC MODEL IS TO BE SOUGHT AT EACH ITERATION, WHILE SLVBQP
C  IS .FALSE. IF A LESS ACCURATE SOLUTION IS DESIRED.
C
      SLVBQP = ICHOSE( 6 ) .EQ. 1
C
C  UNSUCC IS .TRUE. IF THE LAST ATTEMPTED STEP PROVED UNSUCCESSFUL.
C
      UNSUCC = .FALSE.
C
C  PRINT DETAILS OF THE OBJECTIVE FUNCTION CHARACTERISTICS.
C
      IF ( IPRINT .GE. 3 ) THEN
         WRITE( IOUT, 2310 )
         IF ( IPRINT .GE. 10 .OR. NG .LE. 100 ) THEN
            DO 10 IG = 1, NG
               K1    = ISTADG( IG )
               K2    = ISTADG( IG + 1 ) - 1
C
C  PRINT DETAILS OF THE GROUPS.
C
               IF ( K1 .LE. K2 ) THEN
                  IF ( K1 .EQ. K2 ) THEN
                     WRITE( IOUT, 2380 ) IG, 1, IELING( K1 )
                  ELSE
                     WRITE( IOUT, 2320 ) IG, K2 - K1 + 1,
     *                                 ( IELING( K ), K = K1, K2 )
                  END IF
               ELSE
                  WRITE( IOUT, 2330 ) IG
               END IF
               IF ( .NOT. GXEQX( IG ) ) WRITE( IOUT, 2390 )
C
C  PRINT DETAILS OF THE NONLINEAR ELEMENTS.
C
               K1 = ISTADA( IG )
               K2 = ISTADA( IG + 1 ) - 1
               IF ( K1 .LE. K2 )  WRITE( IOUT, 2340 )
     *                          ( ICNA( J ), J = K1, K2 )
   10       CONTINUE
         END IF
         IF ( .NOT. ALLLIN ) THEN
            WRITE( IOUT, 2350 )
            IF ( IPRINT .GE. 10 .OR. NEL .LE. 100 ) THEN
               DO 20 IEL = 1, NEL
                  K1     = ISTAEV( IEL )
                  K2     = ISTAEV( IEL + 1 ) - 1
C
C  PRINT DETAILS OF THE NONLINEAR ELEMENTS.
C
                  IF ( K1 .LE. K2 ) THEN
                     WRITE( IOUT, 2360 ) IEL, INTVAR( IEL ), 
     *                  K2 - K1 + 1, ( IELVAR( K ), K = K1, K2 )
                  ELSE
                     WRITE( IOUT, 2370 ) IEL
                  END IF
   20          CONTINUE
            END IF
         END IF
      END IF
C
C  PARTITION THE WORKSPACE ARRAYS FUVALS, IWK AND WK. INITIALIZE 
C  CERTAIN PORTIONS OF IWK AND WK.
C
CS    CALL SINITW( N , NG, NEL   , IELING, LELING, ISTADG, LSTADG,
CD    CALL DINITW( N , NG, NEL   , IELING, LELING, ISTADG, LSTADG,
     *             IELVAR, LELVAR, ISTAEV, LSTAEV, INTVAR, LNTVAR,
     *             ISTADH, LSTADH, ICNA  , LICNA , ISTADA, LSTADA,
     *             GXEQX , LGXEQX, INTREP, LINTRE, LFUVAL, ALTRIV, 
     *             DIRECT, FDGRAD, LFXI  , LGXI  , LHXI  , LGGFX , 
     *             LDX   , LGRJAC, LQGRAD, LBREAK, LP    , LXCP  , 
     *             LX0   , LGX0  , LDELTX, LBND  , LWKSTR, LSPTRS, 
     *             LSELTS, LINDEX, LSWKSP, LSTAGV, LSTAJC, LIUSED, 
     *             LFREEC, LNNONZ, LNONZ2, LSYMMD, LSYMMH, LSLGRP, 
     *             LSVGRP, LGCOLJ, LVALJR, LSEND , LNPTRS, LNELTS, 
     *             LNNDEX, LNWKSP, LNSTGV, LNSTJC, LNIUSE, LNFREC, 
     *             LNNNON, LNNNO2, LNSYMD, LNSYMH, LNLGRP, LNVGRP, 
     *             LNGCLJ, LNVLJR, LNQGRD, LNBRAK, LNP   , LNBND , 
     *             LNFXI , LNGXI , LNGUVL, LNHXI , LNHUVL, LNGGFX, 
     *             LNDX  , LNGRJC, LIWK2 , LWK2  , MAXSIN, NINVAR, 
     *             NTYPE , NSETS , MAXSEL, LSTYPE, LSSWTR, LSSIWT, 
     *             LSIWTR, LSWTRA, LNTYPE, LNSWTR, LNSIWT, LNIWTR, 
     *             LNWTRA, LSISET, LSSVSE, LNISET, LNSVSE, RANGES,
     *             IWK   , LIWK  , WK    , LWK   ,
     *             IPRINT, IOUT  , INFORM )
      IF ( INFORM .NE. 0 ) RETURN 
      NN    = NINVAR + N
      LNWK  = MAX( NG, MAXSEL )
      LNWKB = MAXSIN
      LNWKC = MAXSIN
      LWKB  = LNWK + 1
      LWKC  = LWKB + LNWKB
C
C  ---------------------------------------------------------------------
C  STEP 0 OF THE ALGORITHM (SEE PAPER).
C  ---------------------------------------------------------------------
C
      DO 30 I = 1, N
C
C  PROJECT USER SUPPLIED STARTING POINT INTO THE FEASIBLE BOX.
C
         BLI = BL( I )
         BUI = BU( I )
C ** Correction 1. 22/07/92: 5 lines added **
         IF ( BLI .GT. BUI ) THEN
            IF ( IPRINT .GE. 0 ) WRITE( IOUT, 2960 ) BLI, I, BUI
            INFORM = 8
            RETURN
         END IF
C ** Correction 1. 22/07/92: end of correction **
         XI  = MAX( BLI, MIN( BUI, X( I ) ) )
C
C  FIND THE MAXIMUM VARIABLE SCALE FACTOR.
C
         VSCMAX = MAX( VSCMAX, VSCALE( I ) )
C
C  FIND INITIAL ACTIVE SET.
C
         IS = 0
C ** Correction 10. 03/12/98: 4 lines corrected **
         IF ( XI .LE. BLI * ( ONE + SIGN( EPSTLP, BLI ) ) ) IS = 1
         IF ( XI .GE. BUI * ( ONE - SIGN( EPSTLN, BUI ) ) ) IS = 2
         IF ( BUI * ( ONE - SIGN( EPSTLN, BUI ) ) .LE. 
     *        BLI * ( ONE + SIGN( EPSTLP, BLI ) ) ) IS = 3
C ** Correction 10. 03/12/98: End of correction **
         IWK( LINDEX + I ) = IS
C        IF ( IPRINT .GE. 0 ) ISTATE( I ) = IS
         IF ( IS .EQ. 3 ) XI = HALF * ( BLI + BUI )
C
C  COPY THE INITIAL POINT INTO XT PRIOR TO CALCULATING FUNCTION VALUES.
C
         X( I )  = XI
         XT( I ) = XI
   30 CONTINUE
C
C  THE ACTIVE SET (THE SET OF VARIABLES ON THEIR BOUNDS) REMAINS FIXED
C  AFTER ICORAS ITERATIONS.
C
C     ICORAS = 0
C
C  ENSURE THAT ALL THE ELEMENT FUNCTIONS ARE EVALUATED
C  AT THE INITIAL POINT.
C
      NCALCF         = NEL
      DO 35 I        = 1, NCALCF
         ICALCF( I ) = I
   35 CONTINUE
C
C  RETURN TO THE CALLING PROGRAM TO OBTAIN THE ELEMENT FUNCTION
C  AND, IF POSSIBLE, DERIVATIVE VALUES AT THE INITIAL POINT.
C
      IF ( FDGRAD ) IGETFD = 0
      INFORM = - 1
      NGEVAL = NGEVAL + 1
      RETURN
   40 CONTINUE
C
C  IF FINITE-DIFFERENCE GRADIENTS ARE USED, COMPUTE THEIR VALUES.
C
      IF ( FDGRAD .AND. .NOT. ALLLIN ) THEN
C
C  STORE THE VALUES OF THE NONLINEAR ELEMENTS FOR FUTURE USE.
C
         IF ( IGETFD .EQ. 0 ) THEN
CS          CALL SCOPY( NEL, FUVALS, 1, WK( LWKSTR + 1 ), 1 )
CD          CALL DCOPY( NEL, FUVALS, 1, WK( LWKSTR + 1 ), 1 )
            CENTRL = ICHOSE( 3 ) .EQ. 2
         END IF 
C
C  OBTAIN A FURTHER SET OF DIFFERENCES.
C
CS       CALL SFDGRD( N, NEL, IELVAR, LELVAR, ISTAEV, LSTAEV,
CD       CALL DFDGRD( N, NEL, IELVAR, LELVAR, ISTAEV, LSTAEV,
     *                IWK   ( LSELTS + 1   ), LNELTS,
     *                IWK   ( LSPTRS + 1   ), LNPTRS,
     *                IELING, LELING,
     *                IWK   ( LSSVSE + 1   ), LNSVSE,
     *                IWK   ( LSISET + 1   ), LNISET, NSETS , 
     *                IWK   ( LSEND  + 1   ), NEL   , 
     *                ICALCF, LCALCF, NCALCF, INTVAR, LNTVAR,
     *                IWK   ( LSSWTR + 1   ), LNSWTR,
     *                IWK   ( LSSIWT + 1   ), LNSIWT, 
     *                IWK   ( LSTYPE + 1   ), LNTYPE, NTYPE ,
     *                IWK   ( LSIWTR + 1   ), LNIWTR,
     *                WK    ( LSWTRA + 1   ), LNWTRA, X , XT, 
     *                FUVALS, LFUVAL, CENTRL, IGETFD )
         IF ( IGETFD .GT. 0 ) RETURN
C
C  RESTORE THE VALUES OF THE NONLINEAR ELEMENTS AT X..
C
         IGETFD = NSETS + 1
CS       CALL SCOPY( NEL, WK( LWKSTR + 1 ), 1, FUVALS, 1 )
CD       CALL DCOPY( NEL, WK( LWKSTR + 1 ), 1, FUVALS, 1 )
      END IF
C
C  THE CONVERGENCE TOLERANCE IS MODIFIED TO REFLECT THE SCALING.
C
      EPSCNS = STOPG * VSCMAX
C
C  COMPUTE THE MINIMUM NORM OF THE RESIDUAL THAT IS TO BE REQUIRED
C  WHEN OBTAINING THE APPROXIMATE MINIMIZER OF THE MODEL PROBLEM.
C
CS    RESMIN = MIN( TENM4, MAX( EPSGRD, EPSCNS ** 2.02 ) )
CD    RESMIN = MIN( TENM4, MAX( EPSGRD, EPSCNS ** 2.02 ) )
C
C  COMPUTE THE GROUP ARGUMENT VALUES FT.
C
      DO 55 IG = 1, NG
         FTT   = - B( IG )
C
C  INCLUDE THE CONTRIBUTION FROM THE LINEAR ELEMENT.
C
         DO 45 J = ISTADA( IG ), ISTADA( IG + 1 ) - 1
            FTT  = FTT + A( J ) * X( ICNA( J ) )
   45    CONTINUE
C
C  INCLUDE THE CONTRIBUTIONS FROM THE NONLINEAR ELEMENTS.
C
         DO 50 J = ISTADG( IG ), ISTADG( IG + 1 ) - 1
            FTT  = FTT + ESCALE( J ) * FUVALS( IELING( J ) )
   50    CONTINUE
         FT( IG ) = FTT
   55 CONTINUE
C
C  COMPUTE THE GROUP FUNCTION VALUES.
C
      IF ( ALTRIV ) THEN
CS       F = SDOT( NG, GSCALE, 1, FT, 1 )
CD       F = DDOT( NG, GSCALE, 1, FT, 1 )
CS       CALL SCOPY( NG, FT, 1, GVALS( 1, 1 ), 1 )
CD       CALL DCOPY( NG, FT, 1, GVALS( 1, 1 ), 1 )
CS       CALL SSETVL( NG, GVALS( 1, 2 ), 1, ONE )
CD       CALL DSETVL( NG, GVALS( 1, 2 ), 1, ONE )
      ELSE
C
C  IF NECESSARY, RETURN TO THE CALLING PROGRAM TO OBTAIN THE GROUP
C  FUNCTION AND DERIVATIVE VALUES AT THE INITIAL POINT.
C  ENSURE THAT ALL THE GROUP FUNCTIONS ARE EVALUATED
C  AT THE INITIAL POINT.
C
         NCALCG          = NG
         DO 65 IG        = 1, NG
            ICALCG( IG ) = IG
            IF ( GXEQX( IG ) ) THEN
               GVALS( IG, 1 ) = FT( IG )
               GVALS( IG, 2 ) = ONE
            END IF
   65    CONTINUE
         INFORM = - 2
         RETURN
      END IF
   70 CONTINUE
      REUSEC = .FALSE.
      IF ( .NOT. ALTRIV ) THEN
         F        = ZERO
CDIR$ IVDEP
         DO 75 IG = 1, NG
            F     = F + GSCALE( IG ) * GVALS( IG, 1 )
   75    CONTINUE
      END IF
C
C  IF A SECANT METHOD IS TO BE USED, INITIALIZE THE SECOND
C  DERIVATIVES OF EACH ELEMENT AS A SCALED IDENTITY MATRIX.
C
      T  = CPUTIM( DUM )
      IF ( .NOT. SECOND .AND. .NOT. ALLLIN ) THEN
CS       CALL SSCALH( .TRUE., N, NEL, NCALCF, LHXI, ISTAEV, LSTAEV,
CD       CALL DSCALH( .TRUE., N, NEL, NCALCF, LHXI, ISTAEV, LSTAEV,
     *                ISTADH, LSTADH, ICALCF, LCALCF, 
     *                INTVAR, LNTVAR, IELVAR, LELVAR, IWK( LSYMMD + 1 ),
     *                MAXSIN, INTREP, LINTRE, FUVALS, LFUVAL, 
     *                WK( NINVAR + 1 ), LNWK, WK( NN + 1 ), LNWKB,
     *                WK( LP + 1 ), WK, NINVAR, RANGES )
      END IF
C
C  IF A TWO-NORM TRUST REGION IS TO BE USED, INITIALIZE THE VECTOR P.
C
CS    IF ( TWONRM ) CALL SSETVL( N, WK( LP + 1 ), 1, ZERO )
CD    IF ( TWONRM ) CALL DSETVL( N, WK( LP + 1 ), 1, ZERO )
      TUP = TUP + CPUTIM( DUM ) - T
C
C  COMPUTE THE GRADIENT VALUE.
C
CS    CALL SELGRD( N, NG, .TRUE., ICNA, LICNA, ISTADA, LSTADA,
CD    CALL DELGRD( N, NG, .TRUE., ICNA, LICNA, ISTADA, LSTADA,
     *             IELING, LELING, ISTADG, LSTADG, ISTAEV, LSTAEV,
     *             IELVAR, LELVAR, INTVAR, LNTVAR, IWK( LSVGRP + 1 ),
     *             LNVGRP, IWK( LSTAJC + 1 ), LNSTJC, 
     *             IWK( LSTAGV + 1 ), LNSTGV, A, LA, GVALS( 1, 2 ), 
C ** Correction 4. 05/02/93: 1 line corrected **
     *             LGVALS, FUVALS, LHXI, FUVALS( LGGFX + 1 ),
C ** Correction 4. 05/02/93: end of correction **
     *             GSCALE, LGSCAL, 
     *             ESCALE, LESCAL, FUVALS( LGRJAC + 1 ),
     *             LNGRJC, WK, WK( N + 1 ), MAXSEL,
     *             GXEQX, LGXEQX, INTREP, LINTRE, RANGES )
C
C  FIND THE INITIAL PROJECTED GRADIENT AND ITS NORM.
C
CS    CALL SPRGRA( N, X, FUVALS( LGGFX + 1 ), VSCALE, BL, BU, DGRAD,
CD    CALL DPRGRA( N, X, FUVALS( LGGFX + 1 ), VSCALE, BL, BU, DGRAD,
     *             IVAR, NVAR, PJGNRM )
      NFREE = NVAR
C
C  FIND THE NORM OF THE PROJECTED GRADIENT.
C
      IF ( PRCOND .AND. NVAR .GT. 0 .AND. MYPREC ) THEN
C
C  USE THE USERS PRECONDITIONER.
C
         INFORM = - 8
         RETURN
      END IF
   90 CONTINUE
C
C  FIND THE NORM OF THE 'PRECONDITIONED' PROJECTED GRADIENT. ALSO,
C  FIND THE DIAGONAL ELEMENTS OF THE ASSEMBLED HESSIAN AS SCALINGS,
C  IF REQUIRED.
C
CS    CALL SGTPGR( N, NG, NGEL, NN, NVAR, MAXSEL, QGNORM, SMALLH,
CD    CALL DGTPGR( N, NG, NGEL, NN, NVAR, MAXSEL, QGNORM, SMALLH,
     *             PJGNRM, CALCDI, DPRCND, MYPREC, IVAR, ISTADH, LSTADH,
     *             ISTAEV, LSTAEV, IELVAR, LELVAR, INTVAR, LNTVAR,
     *             IELING, LELING, IWK( LSLGRP + 1 ), LNLGRP,
     *             IWK( LGCOLJ + 1 ), LNGCLJ, IWK( LSTAGV + 1 ),
     *             LNSTGV, IWK( LSVGRP + 1 ), LNVGRP, IWK( LVALJR + 1 ),
     *             LNVLJR, IWK( LSYMMD + 1 ), IWK( LSYMMH + 1 ), MAXSIN,
     *             DGRAD, Q, GVALS( 1, 2 ), GVALS( 1, 3 ), 
     *             FUVALS( LDX + 1 ), GSCALE, ESCALE, LESCAL,
     *             FUVALS( LGRJAC + 1 ), LNGRJC, FUVALS, LNHUVL,
     *             WK, LNWK, WK( LWKB ), LNWKB, WK( LWKC ), LNWKC,
     *             GXEQX, LGXEQX, INTREP, LINTRE, RANGES ) 
C
C  SET INITIAL TRUST REGION RADIUS.
C
      IF ( RADIUS .LT. ZERO ) THEN
C
C  ENSURE THAT THE INITIAL CAUCHY STEP IS OF ORDER UNITY.
C
CS       GNORM  = SNRM2( N, FUVALS( LGGFX + 1 ), 1 )
CD       GNORM  = DNRM2( N, FUVALS( LGGFX + 1 ), 1 )
         RADIUS = MIN( RADMAX, POINT1 * GNORM ) 
C
C  ENSURE THAT A LINEARIZED MODEL OF THE OBJECTIVE FUNCTION GIVES A
C  VALUE NO SMALLER THAN THE KNOWN LOWER BOUND, FLOWER.
C
C        IF ( ( F .GT. FLOWER ) .AND. ( FLOWER .GT. - BIG ) )
C    *      RADIUS = MIN( RADIUS, ( F - FLOWER ) / GNORM ** 2 )
      END IF
      OLDRAD = RADIUS 
      PRINT1 = .TRUE.
C
C  ---------------------------------------------------------------------
C  MAIN ITERATION LOOP OF THE ALGORITHM (SEE PAPER).
C  ---------------------------------------------------------------------
C
  100 CONTINUE
C
C  IF REQUIRED, PRINT ONE LINE OF DETAILS OF THE CURRENT ITERATION.
C
      IF ( IPRINT .EQ. 1 .OR. IPRINT .EQ. 2 ) THEN
C
C  IF NEEDED, PRINT THE ITERATION HEADER.
C
         IF ( PRINT1 .OR. IPRINT .EQ. 2 ) THEN
            IF ( DIRECT ) THEN
               WRITE( IOUT, 2600 )
            ELSE
               WRITE( IOUT, 2000 )
            END IF
         END IF
C
C  PRINT THE ITERATION DETAILS.
C
         ATIME = CHRTIM( CPUTIM( DUM ) - TIME ) 
         IF ( DIRECT ) THEN
            IF ( PRINT1 ) THEN
               WRITE( IOUT, 2810 ) ITER, NGEVAL,         F * FINDMX,
     *                           PJGNRM,
     *                           NFREE, ATIME
            ELSE IF ( INFORM .EQ. - 11 ) THEN
               WRITE( IOUT, 2870 ) ITER, NGEVAL, FILL,   F * FINDMX, 
     *                           PJGNRM,      OLDRAD, STEP, LISEND, 
     *                           NFREE, ATIME
            ELSE
               WRITE( IOUT, 2610 ) ITER, NGEVAL, FILL,   F * FINDMX, 
     *                           PJGNRM, RHO, OLDRAD, STEP, LISEND, 
     *                           NFREE, ATIME
            END IF
         ELSE
            IF ( PRINT1 ) THEN
               WRITE( IOUT, 2800 ) ITER, NGEVAL, ITERCG, F * FINDMX,
     *                           PJGNRM, NFREE, ATIME
            ELSE IF ( INFORM .EQ. - 11 ) THEN
               WRITE( IOUT, 2860 ) ITER, NGEVAL, ITERCG, F * FINDMX,
     *                           PJGNRM, OLDRAD, STEP, CGEND, NFREE,
     *                           ATIME
            ELSE
               WRITE( IOUT, 2010 ) ITER, NGEVAL, ITERCG, F * FINDMX, 
     *                           PJGNRM, RHO, OLDRAD, STEP, CGEND, 
     *                           NFREE, ATIME
            END IF
         END IF
         PRINT1 = .FALSE.
         IF ( IPRINT .EQ. 2 ) THEN
            WRITE( IOUT, 2760 ) EPSCNS
         END IF
      END IF
C
C  ---------------------------------------------------------------------
C  STEP 1 OF THE ALGORITHM (SEE PAPER).
C  ---------------------------------------------------------------------
C
C  IF REQUIRED, PRINT MORE THOROUGH DETAILS OF THE CURRENT ITERATION.
C
      IF ( IPRINT .GE. 3 ) THEN
         IF ( ITER .EQ. 0 ) THEN
            WRITE( IOUT, 2830 ) ITER, F * FINDMX, NGEVAL, PJGNRM, 
     *                 ITERCG,          ISKIP
C    *               , ICORAS
         ELSE
            WRITE( IOUT, 2170 ) ITER, F * FINDMX, NGEVAL, PJGNRM, 
     *                  ITERCG, OLDRAD, ISKIP
C    *                , ICORAS
         END IF
         WRITE( IOUT, 2050 ) ( X( KK ), KK = 1, N )
         IF ( IPRINT .GE. 4 ) THEN
            WRITE( IOUT, 2060 ) ( FUVALS( LGGFX + KK ) * FINDMX,
     *                            KK = 1, N )
            IF ( IPRINT .GE. 5 ) THEN
               WRITE( IOUT, 2270 ) ( FUVALS( KK ),   KK = 1, NEL )
               WRITE( IOUT, 2280 ) ( GVALS( KK, 1 ), KK = 1, NG )
               IF ( CALCDI ) WRITE( IOUT, 2080 )
     *                   ( FUVALS( LDX + KK ) * FINDMX, KK = 1, N )
               IF ( IPRINT .GE. 6 ) THEN
                  LL = LHXI - LGXI
                  IF ( LL .GT. 0 ) WRITE( IOUT, 2790 )
     *               ( FUVALS( LGXI + KK ), KK = 1, LL )
                  WRITE( IOUT, 2780 ) ( GVALS( KK, 2 ), KK = 1, NG )
                  LL = LGGFX - LHXI
                  IF ( LL .GT. 0 ) WRITE( IOUT, 2180 )
     *               ( FUVALS( LHXI + KK ), KK = 1, LL )
                  IF ( IPRINT .GE. 7 ) THEN
                     WRITE( IOUT, 2901 ) ( IWK( K ), K=LSPTRS+1, LSELTS)
                     WRITE( IOUT, 2902 ) ( IWK( K ), K=LSELTS+1, LINDEX)
                     WRITE( IOUT, 2904 ) ( IWK( K ), K=LINDEX+1, LSWKSP)
                     WRITE( IOUT, 2905 ) ( IWK( K ), K=LSWKSP+1, LIUSED)
                     WRITE( IOUT, 2908 ) ( IWK( K ), K=LIUSED+1, LFREEC)
                     WRITE( IOUT, 2917 ) ( IWK( K ), K=LFREEC+1, LNNONZ)
                     WRITE( IOUT, 2909 ) ( IWK( K ), K=LNNONZ+1, LNONZ2)
                     WRITE( IOUT, 2910 ) ( IWK( K ), K=LNONZ2+1, LSYMMD)
                     WRITE( IOUT, 2911 ) ( IWK( K ), K=LSYMMD+1, LSYMMH)
                     WRITE( IOUT, 2912 ) ( IWK( K ), K=LSYMMH+1, LSLGRP)
                     WRITE( IOUT, 2913 ) ( IWK( K ), K=LSLGRP+1, LSTAJC)
                     WRITE( IOUT, 2907 ) ( IWK( K ), K=LSTAJC+1, LSTAGV)
                     WRITE( IOUT, 2906 ) ( IWK( K ), K=LSTAGV+1, LSVGRP)
                     WRITE( IOUT, 2914 ) ( IWK( K ), K=LSVGRP+1, LGCOLJ)
                     WRITE( IOUT, 2915 ) ( IWK( K ), K=LGCOLJ+1, LVALJR)
                     WRITE( IOUT, 2916 ) ( IWK( K ), K=LVALJR+1, LSTYPE)
                     WRITE( IOUT, 2923 ) ( IWK( K ), K=LSTYPE+1, LSSWTR)
                     WRITE( IOUT, 2918 ) ( IWK( K ), K=LSSWTR+1, LSSIWT)
                     WRITE( IOUT, 2919 ) ( IWK( K ), K=LSSIWT+1, LSIWTR)
                     WRITE( IOUT, 2920 ) ( IWK( K ), K=LSIWTR+1, LSISET)
                     WRITE( IOUT, 2921 ) ( IWK( K ), K=LSISET+1, LSSVSE)
                     WRITE( IOUT, 2922 ) ( IWK( K ), K=LSSVSE+1, LSEND )
                  END IF
               END IF
            END IF
         END IF
      END IF
C
C  TEST FOR CONVERGENCE.
C
      IF ( PJGNRM .LE. EPSCNS ) THEN
         INFORM = 0
         GO TO 600
      END IF
C
C  TEST WHETHER THE MAXIMUM ALLOWED NUMBER OF ITERATIONS HAS
C  BEEN REACHED.
C
      IF ( ITER .GE. MAXIT ) THEN
         IF ( IPRINT .GE. 0 ) WRITE( IOUT, 2090 )
         INFORM = 1
         GO TO 600
      END IF
C
C  TEST WHETHER THE TRUST REGION RADIUS IS TOO SMALL FOR PROGRESS.
C
      IF ( RADIUS .LT. RADTOL ) THEN
         IF ( IPRINT .GE. 0 ) WRITE( IOUT, 2130 )
         INFORM = 2
         GO TO 600
      END IF
      ITER = ITER + 1
C
C  ---------------------------------------------------------------------
C  STEP 2 OF THE ALGORITHM (SEE PAPER).
C  ---------------------------------------------------------------------
C
C  USE IWK TO INDICATE WHICH ELEMENTS ARE NEEDED FOR THE MATRIX-
C  VECTOR PRODUCT B * P. IF IWK(I) = NBPROD, THE I-TH ELEMENT IS USED.
C
      NBPROD = 0
      IF ( .NOT. ALLLIN ) THEN
         DO 210 I             = 1, NGEL
            IWK( LSWKSP + I ) = 0
  210    CONTINUE
      END IF
C
C  ESTIMATE THE NORM OF THE PRECONDITIONING MATRIX BY COMPUTING
C  ITS SMALLEST AND LARGEST (IN MAGNITUDE) DIAGONALS.
C
      DIAMIN = BIG
      DIAMAX = ZERO
      IF ( CALCDI ) THEN
         DO 220 I = 1, N
            IF ( IWK( LINDEX + I ) .EQ. 0 ) THEN
               DIAMIN = MIN( DIAMIN, FUVALS( LDX + I ) )
               DIAMAX = MAX( DIAMAX, FUVALS( LDX + I ) )
            END IF
  220    CONTINUE
      END IF
C
C  IF ALL THE DIAGONALS ARE SMALL, THE NORM WILL BE ESTIMATED AS ONE.
C
      IF ( DIAMAX .LE. EPSMCH ) THEN
         DIAMIN = ONE
         DIAMAX = ONE
      END IF
C
C  INITIALIZE VALUES FOR THE GENERALIZED CAUCHY POINT CALCULATION.
C
      STEPMX = ZERO
      F0     = F
      IBQPST = 1
CDIR$ IVDEP
      DO 230 I = 1, N
C
C  SET THE BOUNDS ON THE VARIABLES FOR THE MODEL PROBLEM. IF A TWO-NORM
C  TRUST REGION IS TO BE USED, THE BOUNDS ARE JUST THE BOX CONSTRAINTS.  
C
         IF ( TWONRM ) THEN
            WK( LBND + I     ) = BL( I )
            WK( LBND + I + N ) = BU( I )
         ELSE
C
C  IF AN INFINITY-NORM TRUST REGION IS TO BE USED, THE BOUNDS ARE THE
C  INTERSECTION OF THE TRUST REGION WITH THE BOX CONSTRAINTS.  
C
            IF ( CALCDI ) THEN
               DISTAN = RADIUS / SQRT( FUVALS( LDX + I ) )
            ELSE
               DISTAN = RADIUS * VSCALE( I )
            END IF
            WK( LBND + I     ) = MAX( X( I ) - DISTAN, BL( I ) )
            WK( LBND + I + N ) = MIN( X( I ) + DISTAN, BU( I ) )
         END IF
C
C  COMPUTE THE CAUCHY DIRECTION, DGRAD, AS A SCALED STEEPEST-DESCENT
C  DIRECTION. NORMALIZE THE DIAGONAL SCALINGS IF NECESSARY.
C
         WK( LX0    + I ) = X( I ) 
         WK( LDELTX + I ) = ZERO
         WK( LGX0   + I ) = FUVALS( LGGFX + I )
         WK( LP     + I ) = ZERO
         IF ( REUSEC ) GO TO 230
         IF ( CALCDI ) THEN
            J           = LDX + I
            DGRAD( I )  = - FUVALS( LGGFX + I ) / FUVALS( J )
            FUVALS( J ) =   FUVALS( J ) / DIAMAX
         ELSE
            DGRAD( I ) = - FUVALS( LGGFX + I ) *
     *                     ( VSCALE( I ) / VSCMAX ) ** 2
         END IF
C
C  IF AN APPROXIMATION TO THE CAUCHY POINT IS TO BE USED, CALCULATE A
C  SUITABLE INITIAL ESTIMATE OF THE LINE MINIMUM, STEPMX.
C
         IF ( .NOT. XACTCP ) THEN
            IF ( DGRAD( I ) .NE. ZERO ) THEN
               IF ( DGRAD( I ) .GT. ZERO ) THEN
                  STEPMX = MAX( STEPMX, ( BU( I ) - X( I ) )/DGRAD( I ))
               ELSE
                  STEPMX = MAX( STEPMX, ( BL( I ) - X( I ) )/DGRAD( I ))
               END IF
            END IF
         END IF
C
C  RELEASE ANY ARTIFICIALLY FIXED VARIABLES FROM THEIR BOUNDS,
C
         IF ( IWK( LINDEX + I ) .EQ. 4 ) IWK( LINDEX + I ) = 0
  230 CONTINUE
C
C  THE VALUE OF THE INTEGER IFACTR CONTROLS WHETHER A NEW FACTORIZATION
C  OF THE HESSIAN OF THE MODEL IS OBTAINED (IFACTR = 1) OR WHETHER A
C  SCHUR-COMPLEMENT UPDATE TO AN EXISTING FACTORIZATION IS REQUIRED
C  (IFACTR = 2) WHEN FORMING THE PRECONDITIONER.
C
      IFACTR = 1
      REFACT = .FALSE.
C
C  IF A PREVIOUSLY CALCULATED GENERALIZED CAUCHY POINT STILL LIES
C  WITHIN THE TRUST-REGION BOUNDS, IT WILL BE REUSED.
C
      IF ( REUSEC ) THEN
C
C  RETRIEVE THE CAUCHY POINT.
C
CS       CALL SCOPY( N, WK( LXCP + 1 ), 1, XT, 1 )
CD       CALL DCOPY( N, WK( LXCP + 1 ), 1, XT, 1 )
C
C  RETRIEVE THE SET OF FREE VARIABLES.
C
         NVAR         = NFREEC
         DO 235 I     = 1, NVAR
            IVAR( I ) = IWK( LFREEC + I )
            IWK( LINDEX + IVAR( I ) ) = 0
  235    CONTINUE
C
C  SKIP THE REMAINDER OF STEP 2.
C
         REUSEC = .FALSE.
         IF ( IPRINT .GT. 1 ) WRITE( IOUT, 2300 )
         GO TO 290   
      END IF
C
C  EVALUATE THE GENERALISED CAUCHY POINT, XT.
C
      JUMPTO = 1
      FIRSTC = .TRUE.
  240 CONTINUE
      T = CPUTIM( DUM )
      IF ( XACTCP ) THEN
C
C  THE EXACT CAUCHY POINT IS REQUIRED.
C
CS       CALL SCAUCH( N, WK( LX0 + 1 ), XT, WK( LGX0 + 1 ), 
CD       CALL DCAUCH( N, WK( LX0 + 1 ), XT, WK( LGX0 + 1 ), 
     *                WK( LBND + 1 ), IWK( LINDEX + 1 ), F0, EPSTLP, 
     *                TWONRM, DXSQR, FMODEL, DGRAD, Q, IVAR, NVAR, 
     *                NVAR1, NVAR2, NNONNZ, IWK( LNNONZ + 1 ), 
     *                WK( LBREAK + 1 ), IOUT, JUMPTO, IPRINT )
      ELSE
C
C  AN APPROXIMATION TO THE CAUCHY POINT SUFFICES.
C
CS       CALL SINXAC( N, WK( LX0 + 1 ), XT, WK( LGX0 + 1 ), 
CD       CALL DINXAC( N, WK( LX0 + 1 ), XT, WK( LGX0 + 1 ), 
     *                WK( LBND + 1 ), IWK( LINDEX + 1 ), F0, EPSTLP, 
     *                STEPMX, POINT1, TWONRM, DXSQR, FMODEL, DGRAD, Q,
     *                IVAR, NVAR, NVAR1, NVAR2,
     *                NNONNZ, IWK( LNNONZ + 1 ), WK( LBREAK + 1 ),
     *                WK( LQGRAD + 1 ), IOUT, JUMPTO, IPRINT )
      END IF
      TCA = TCA + CPUTIM( DUM ) - T
C
C  SCATTER THE NONZEROS IN DGRAD ONTO WK( LP ).
C
      DO 250 J        = NVAR1, NVAR2
         I            = IVAR( J )
         WK( LP + I ) = DGRAD( I )
  250 CONTINUE
C
C  A FURTHER MATRIX-VECTOR PRODUCT IS REQUIRED.
C
      IF ( JUMPTO .GT. 0 ) THEN
         T      = CPUTIM( DUM )
         NBPROD = NBPROD + 1
C
C  CALCULATE THE PRODUCT OF THE HESSIAN WITH THE VECTOR P.
C
         DENSEP = JUMPTO .EQ. 2 .OR. ( XACTCP .AND. JUMPTO .EQ. 4 )
CS       CALL SHSPRD( N, NN, NG, NGEL, NVAR, NVAR1, NVAR2, NBPROD,
CD       CALL DHSPRD( N, NN, NG, NGEL, NVAR, NVAR1, NVAR2, NBPROD,
     *                ALLLIN, IVAR, ISTAEV, LSTAEV, ISTADH, LSTADH, 
     *                INTVAR, LNTVAR, IELING, LELING, IELVAR, LELVAR,
     *                IWK( LSTAJC + 1 ), LNSTJC, IWK( LSELTS + 1 ),
     *                LNELTS, IWK( LSPTRS + 1 ), LNPTRS, 
     *                IWK( LGCOLJ + 1 ), LNGCLJ, IWK( LSLGRP + 1 ),
     *                LNLGRP, IWK( LSWKSP + 1 ), LNWKSP, 
     *                IWK( LSVGRP + 1 ), LNVGRP, IWK( LSTAGV + 1 ),
     *                LNSTGV, IWK( LVALJR + 1 ), LNVLJR, 
     *                NNONNZ, IWK( LNNONZ + 1 ), LNNNON,
     *                IWK( LIUSED + 1 ), LNIUSE, IWK( LNONZ2 + 1 ),
     *                LNNNO2, IWK( LSYMMH + 1 ), MAXSIN, 
     *                WK( LP + 1 ), Q, GVALS( 1, 2 ),
     *                GVALS( 1, 3 ), FUVALS( LGRJAC + 1 ), LNGRJC,
     *                GSCALE, ESCALE, LESCAL, FUVALS, LNHUVL,
     *                WK, LNWK, WK( LWKB ), LNWKB, WK( LWKC ), LNWKC,
     *                GXEQX, LGXEQX, INTREP, LINTRE,
     *                DENSEP, RANGES )
         IF ( IPRINT .GE. 5 .AND. JUMPTO .EQ. 3 )
     *       WRITE( IOUT, 2710 ) ( IWK( LNNONZ + I ), I = 1, NNONNZ )
         TMV = TMV + CPUTIM( DUM ) - T
C
C  RESET THE COMPONENTS OF WK( LP ) THAT HAVE CHANGED TO ZERO.
C
         DO 280 J                = NVAR1, NVAR2
            WK( LP + IVAR( J ) ) = ZERO
  280    CONTINUE
C
C  IF REQUIRED, PRINT A LIST OF THE NONZEROS OF WK( LP ).
C 
         IF ( JUMPTO .EQ. 3 .AND. IPRINT .GE. 5 .AND. .NOT. ALLLIN )
     *      WRITE( IOUT, 2580 ) NBPROD,
     *           ( IWK( LSWKSP + I ), I = 1, NGEL )
C
C  CONTINUE THE CAUCHY POINT CALCULATION.
C
         GO TO 240
      END IF
C
C  STORE THE CAUCHY POINT FOR FUTURE USE.
C
CS    CALL SCOPY( N, XT, 1, WK( LXCP + 1 ), 1 )
CD    CALL DCOPY( N, XT, 1, WK( LXCP + 1 ), 1 )
C
C  STORE THE SET OF FREE VARIABLES AT CAUCHY POINT FOR FUTURE USE.
C
      NFREEC               = NVAR
      DO 285 I             = 1, NFREEC
         IWK( LFREEC + I ) = IVAR( I )
  285 CONTINUE   
C
C  SEE IF AN ACCURATE APPROXIMATION TO THE MINIMUM OF THE QUADRATIC
C  MODEL IS TO BE SOUGHT.
C
      IF ( SLVBQP ) THEN
C
C  FIX THE VARIABLES WHICH THE CAUCHY POINT PREDICTS ARE ACTIVE AT THE
C  SOLUTION.
C
         IF ( FIRSTC ) THEN
            FIRSTC = .FALSE.
            DO 286 I = 1, N
               IF ( IWK( LINDEX + I ) .EQ. 1 .OR.
     *              IWK( LINDEX + I ) .EQ. 2 ) THEN
                  IWK( LINDEX + I )  = 4
                  WK( LBND + I )     = XT( I )
                  WK( LBND + N + I ) = XT( I )
               END IF
  286       CONTINUE   
         ELSE
C
C  UPDATE THE STEP TAKEN AND THE SET OF VARIABLES WHICH ARE CONSIDERED
C  FREE.
C
            NVAR     = 0
            DO 288 I = 1, N
               WK( LP + I ) = WK( LP + I ) + WK( LDELTX + I )
               IF ( WK( LP + I ) .NE. ZERO .OR.
     *              IWK( LINDEX + I ) .EQ. 0 ) THEN
                  NVAR = NVAR + 1
                  IVAR( NVAR ) = I
               END IF
  288       CONTINUE   
            NVAR2 = NVAR
         END IF
      END IF
C
C  IF REQUIRED, PRINT THE ACTIVE SET AT THE GENERALIZED CAUCHY POINT.
C
  290 CONTINUE   
      IF ( IPRINT .GE. 4 ) THEN
         WRITE( IOUT, 2030 )
         DO 295 I = 1, N
            IF ( IWK( LINDEX + I ) .EQ. 2 .AND. XT( I ) .GE. BU( I )
     *           - ABS( BU( I ) ) * EPSTLN ) WRITE( IOUT, 2100 ) I
            IF ( IWK( LINDEX + I ) .EQ. 1 .AND. XT( I ) .LE. BL( I )
     *           + ABS( BL( I ) ) * EPSTLP ) WRITE( IOUT, 2110 ) I
            IF ( IWK( LINDEX + I ) .EQ. 4 )  WRITE( IOUT, 2770 ) I
  295    CONTINUE
      END IF
C
C  ---------------------------------------------------------------------
C  STEP 3 OF THE ALGORITHM (SEE PAPER).
C  ---------------------------------------------------------------------
C
      JUMPTO = 1
C
C  IF AN ITERATIVE METHOD IS TO BE USED, SET UP CONVERGENCE TOLERANCES.
C
      CGSTOP = MAX( RESMIN, MIN( ACCCG, QGNORM ) * QGNORM * QGNORM ) *
     *                                             DIAMIN / DIAMAX
CS    IF ( TWONRM .AND. .NOT. DIRECT ) DXSQR = SDOT( N, WK( LP + 1 ), 1,
CD    IF ( TWONRM .AND. .NOT. DIRECT ) DXSQR = DDOT( N, WK( LP + 1 ), 1,
     *                                               WK( LP + 1 ), 1 )
      IF ( IPRINT .GE. 5 .AND. IOUT .GT. 0 .AND.
     *     TWONRM .AND. .NOT. DIRECT ) WRITE( IOUT, 2900 ) SQRT( DXSQR )
      STEP = RADIUS
C
C  IF AN INCOMPLETE FACTORIZATION PRECONDITIONER IS TO BE USED, DECIDE
C  ON THE SEMI-BANDWIDTH, NSEMIB, OF THE PRECONDITIONER. 
C  FOR THE EXPANDING BAND METHOD, THE ALLOWABLE BANDWIDTH INCREASES AS
C  THE SOLUTION IS APPROACHED.
C
      IF ( IPRCND ) THEN
         NSEMIB = N / 5
         IF ( PJGNRM .LE. POINT1 ) NSEMIB = N / 2
         IF ( PJGNRM .LE. TENM2  ) NSEMIB = N
      ELSE
         IF ( .NOT. BAND ) NSEMIB = N
      END IF
C
C  IF MUNKSGAARDS PRECONDITIONER IS TO BE USED, SET THE STABILITY FACTOR
C  REQUIRED BY MC31, TO BE MORE STRINGENT AS THE SOLUTION IS APPROACHED.
C 
      IF ( MUNKS ) THEN
         CMA31 = POINT1
         IF ( PJGNRM .LE. POINT1 ) CMA31 = TENM2
         IF ( PJGNRM .LE. TENM2  ) CMA31 = ZERO
      ELSE
         CMA31 = ZERO
      END IF
C
C  SET A LIMIT ON THE NUMBER OF CG ITERATIONS THAT ARE TO BE ALLOWED.
C
      ITCGMX = N
      IF ( IPRCND .OR. BAND   ) ITCGMX = MAX( 5, N  / ( NSEMIB + 1 ) )
      IF ( SEPREC .OR. GMPSPR ) ITCGMX = 5
C
C  CALCULATE AN APPROXIMATE MINIMIZER OF THE MODEL WITHIN THE SPECIFIED
C  BOUNDS.
C
  300 CONTINUE
      IF ( JUMPTO .EQ. 4 ) GO TO 320
C
C  THE PRODUCT OF THE HESSIAN WITH THE VECTOR WK( LP ) IS REQUIRED.
C
      IF ( JUMPTO .NE. 2 ) THEN
         NBPROD = NBPROD + 1
         NVAR1  = 1
         T      = CPUTIM( DUM )
C
C  SET THE REQUIRED COMPONENTS OF Q TO ZERO.
C
         IF ( JUMPTO .EQ. 1 ) THEN
CS          CALL SSETVL( N, Q, 1, ZERO )
CD          CALL DSETVL( N, Q, 1, ZERO )
         ELSE
CS          CALL SSETVI( NVAR2, Q, IVAR, ZERO )
CD          CALL DSETVI( NVAR2, Q, IVAR, ZERO )
         END IF
C
C  COMPUTE THE MATRIX-VECTOR PRODUCT WITH THE DENSE VECTOR WK( LP ).
C
CS       CALL SHSPRD( N, NN, NG, NGEL, NVAR, NVAR1, NVAR2, NBPROD,
CD       CALL DHSPRD( N, NN, NG, NGEL, NVAR, NVAR1, NVAR2, NBPROD,
     *                ALLLIN, IVAR, ISTAEV, LSTAEV, ISTADH, LSTADH, 
     *                INTVAR, LNTVAR, IELING, LELING, IELVAR, LELVAR,
     *                IWK( LSTAJC + 1 ), LNSTJC, IWK( LSELTS + 1 ),
     *                LNELTS, IWK( LSPTRS + 1 ), LNPTRS, 
     *                IWK( LGCOLJ + 1 ), LNGCLJ, IWK( LSLGRP + 1 ),
     *                LNLGRP, IWK( LSWKSP + 1 ), LNWKSP, 
     *                IWK( LSVGRP + 1 ), LNVGRP, IWK( LSTAGV + 1 ),
     *                LNSTGV, IWK( LVALJR + 1 ), LNVLJR, 
     *                NNONNZ, IWK( LNNONZ + 1 ), LNNNON,
     *                IWK( LIUSED + 1 ), LNIUSE, IWK( LNONZ2 + 1 ),
     *                LNNNO2, IWK( LSYMMH + 1 ), MAXSIN, 
     *                WK( LP + 1 ), Q, GVALS( 1, 2 ),
     *                GVALS( 1, 3 ), FUVALS( LGRJAC + 1 ), LNGRJC,
     *                GSCALE, ESCALE, LESCAL, FUVALS, LNHUVL,
     *                WK, LNWK, WK( LWKB ), LNWKB, WK( LWKC ), LNWKC,
     *                GXEQX, LGXEQX, INTREP, LINTRE,
     *                .TRUE., RANGES )
         TMV = TMV + CPUTIM( DUM ) - T
C
C  IF REQUIRED, PRINT A LIST OF THE NONZEROS OF WK( LP ).
C 
         IF ( IPRINT .GE. 5 .AND. .NOT. ALLLIN ) WRITE( IOUT, 2580 )
     *        NBPROD, ( IWK( LSWKSP + I ), I = 1, NGEL )
         IF ( JUMPTO .EQ. 1 ) THEN
C
C  IF REQUIRED, PRINT THE STEP TAKEN.
C
            IF ( IPRINT .GE. 20 ) WRITE( IOUT, 2880 )
     *                ( WK( LP + I ), I = 1, N )
C
C  COMPUTE THE VALUE OF THE MODEL AT THE GENERALIZED CAUCHY POINT AND
C  THEN RESET WK( LP ) TO ZERO.
C
            FNEW   = FMODEL
            FMODEL = F
CDIR$ IVDEP
            DO 310 J = 1, NVAR2
               I     = IVAR( J )
               FMODEL = FMODEL +
     *            ( FUVALS( LGGFX + I ) + HALF * Q( I ) ) * WK( LP + I )
               WK( LP + I ) = ZERO
  310       CONTINUE
C
C  IF REQUIRED, COMPARE THE RECURRED AND CALCULATED MODEL VALUES.
C
            IF ( IPRINT .GE. 5 ) WRITE( IOUT, 2190 )
     *                FMODEL * FINDMX, FNEW * FINDMX
         END IF
      ELSE
C
C  EVALUATE THE 'PRECONDITIONED' GRADIENT. IF THE USER HAS
C  SUPPLIED A PRECONDITIONER, RETURN TO THE CALLING PROGRAM.
C
         IF ( MYPREC ) THEN
            INFORM = - 9
            RETURN
         ELSE
            IF ( IPRCND .OR. MUNKS .OR. SEPREC .OR. GMPSPR .OR.
     *           BAND ) THEN
C
C  IF REQUIRED, OBTAIN A FULL INVERSE PRECONDITIONER.
C
               T = CPUTIM( DUM )
CS             CALL SPRECN( IFACTR, MUNKS, BAND, SEPREC, N, NG, 
CD             CALL DPRECN( IFACTR, MUNKS, BAND, SEPREC, N, NG, 
     *                      MAXSEL, NFREEF, NFIXED, REFACT, NVAR2, IVAR,
     *                      ISTADH, LSTADH, ICNA, LICNA, 
     *                      ISTADA, LSTADA, INTVAR, LNTVAR, IELVAR,
     *                      LELVAR, IELING, LELING, ISTADG, LSTADG,
     *                      ISTAEV, LSTAEV, IWK( LSTAGV + 1 ),
     *                      LNSTGV, IWK( LSVGRP + 1 ), LNVGRP,
     *                      IWK( LSEND + 1 ), LIWK2, A, LA, FUVALS,
     *                      LNGUVL, FUVALS, LNHUVL, GVALS( 1, 2 ),
     *                      GVALS( 1, 3 ), DGRAD, Q, GSCALE, ESCALE, 
     *                      LESCAL, WK( LWKSTR + 1 ), LWK2, 
     *                      GXEQX, LGXEQX, INTREP, LINTRE, RANGES,
     *                      IPRINT, IOUT, INFOR )
               TLS    = TLS + CPUTIM( DUM ) - T
               IFACTR = 0
C
C  CHECK FOR ERROR RETURNS.
C
               IF ( INFOR .EQ. 10 ) THEN
                  INFORM = 4
                  RETURN
               END IF
               IF ( INFOR .EQ. 11 ) THEN
                  INFORM = 5
                  RETURN
               END IF
            ELSE
C
C  IF REQUIRED, USE A DIAGONAL PRECONDITIONER.
C
CDIR$ IVDEP
               DO 315 J = 1, NVAR2
                  I     = IVAR( J )
                  IF ( DPRCND ) THEN
                     Q( I ) = DGRAD( J ) / FUVALS( LDX + I )
                  ELSE
C
C  NO PRECONDITIONER IS REQUIRED.
C
                     Q( I ) = DGRAD( J ) * VSCALE( I )
                  END IF
  315          CONTINUE
            END IF
         END IF
      END IF
  320 CONTINUE
C
C  THE MINIMIZATION WILL TAKE PLACE OVER ALL VARIABLES WHICH ARE
C  NOT ON THE TRUST REGION BOUNDARY WITH NEGATIVE GRADIENTS PUSHING 
C  OVER THE BOUNDARY.
C
      IF ( DIRECT ) THEN
C
C  - - - - - - - - - - - - DIRECT METHOD - - - - - - - - - - - - - - - -
C
C  MINIMIZE THE QUADRATIC USING A DIRECT METHOD. THE METHOD USED IS
C  A MULTIFRONTAL SYMMETRIC INDEFINITE FACTORIZATION SCHEME.
C  EVALUATE THE GRADIENT OF THE QUADRATIC AT XT.
C
         NVAR     = 0
         GMODEL   = ZERO
         DO 325 I = 1, N
            IF ( IWK( LINDEX + I ) .EQ. 0 ) THEN
               NVAR          = NVAR  + 1
               IVAR( NVAR )  = I
               GI            = FUVALS( LGGFX + I ) + Q( I )
               DGRAD( NVAR ) = GI
               GMODEL        = MAX( GMODEL, ABS( GI ) )
            ELSE
               GI = ZERO
            END IF
            WK( LP + I )     = ZERO
            WK( LQGRAD + I ) = GI
  325    CONTINUE
         NVAR2 = NVAR
C
C  CHECK IF THE GRADIENT OF THE MODEL AT THE GENERALIZED CAUCHY POINT
C  IS ALREADY SMALL ENOUGH. COMPUTE THE (SCALED) STEP
C  MOVED FROM THE PREVIOUS TO THE CURRENT ITERATE.
C
CS       STEP = SDNRM( N, XT, X, TWONRM, VSCALE, .TRUE. )
CD       STEP = DDNRM( N, XT, X, TWONRM, VSCALE, .TRUE. )
C
C  IF THE STEP TAKEN IS SMALL RELATIVE TO THE TRUST REGION RADIUS,
C  ENSURE THAT AN ACCURATE APPROXIMATION TO THE MINIMIZER OF THE
C  MODEL IS FOUND.
C
         IF ( STEP .LE. STPTOL * RADIUS ) THEN
            IF ( MAX( RESMIN, STEP * CGSTOP / ( RADIUS * STPTOL ) )
     *           .GE. GMODEL ) GO TO 400
         ELSE
            IF ( GMODEL * GMODEL .LT. CGSTOP ) GO TO 400
         END IF
C
C  FACTORIZE THE MATRIX AND OBTAIN THE SOLUTION TO THE LINEAR SYSTEM,
C  A DIRECTION OF NEGATIVE CURVATURE OR A DESCENT DIRECTION FOR THE
C  QUADRATIC MODEL.
C
         T = CPUTIM( DUM )
CS       CALL SFRNTL( N , NG, MAXSEL, INTVAR, LNTVAR, IELVAR, LELVAR,
CD       CALL DFRNTL( N , NG, MAXSEL, INTVAR, LNTVAR, IELVAR, LELVAR,
     *                INTREP, LINTRE, IELING, LELING, ISTADG, LSTADG,
     *                ISTAEV, LSTAEV, IWK   ( LSTAGV + 1 )  , LNSTGV,
     *                A , LA, ICNA  , LICNA , ISTADA, LSTADA,
     *                FUVALS, LNGUVL, IWK   ( LSVGRP + 1 )  , LNVGRP,
     *                FUVALS, LNHUVL, ISTADH, LSTADH, GXEQX , LGXEQX,
     *                GVALS( 1, 2 ) , GVALS( 1, 3 ) , IVAR  , NVAR2 ,
     *                WK   ( LQGRAD + 1 )   , WK( LP + 1 )  , XT    ,
     *                WK( LBND + 1 ), FMODEL, GSCALE, ESCALE,
     *                LESCAL, CGSTOP, NUMBER, NEXT  , MODCHL,
     *                IWK( LSEND + 1 )      , LIWK2 , WK( LWKSTR + 1 ),
     *                LWK2  , RANGES, IPRINT, IOUT  , INFOR  )
         TLS  = TLS + CPUTIM( DUM ) - T
         NVAR = NVAR2
C
C  CHECK FOR ERROR RETURNS.
C
         IF ( INFOR .EQ. 10 ) THEN
            INFORM = 4
            RETURN
         END IF
         IF ( INFOR .EQ. 11 ) THEN
            INFORM = 5
            RETURN
         END IF
C
C  SAVE DETAILS OF THE SYSTEM SOLVED.
C
         FILL          = MAX( FILL, RATIO )
         ISYS( INFOR ) = ISYS( INFOR ) + 1
         LS            = LSYS( INFOR )
         LISEND        = LSENDS( INFOR )
C
C  COMPUTE THE (SCALED) STEP FROM THE PREVIOUS TO THE CURRENT ITERATE.
C  IN THE APPROPRIATE NORM.
C
CS       STEP = SDNRM( N, XT, X, TWONRM, VSCALE, .TRUE. )
CD       STEP = DDNRM( N, XT, X, TWONRM, VSCALE, .TRUE. )
C
C  FOR DEBUGGING, COMPUTE THE DIRECTIONAL DERIVATIVE AND CURVATURE 
C  ALONG THE DIRECTION WK( LP ).
C
         IF ( IPRINT .GE. 3 ) THEN
            IF ( .NOT. ALLLIN ) THEN
               DO 350 J             = 1, NGEL
                  IWK( LSWKSP + J ) = 0
  350          CONTINUE
            END IF
            NVAR1    = 0
            DO 355 I = 1, NVAR2
               IF ( IVAR( I ) .GT. 0 ) THEN
                  NVAR1         = NVAR1 + 1
                  IVAR( NVAR1 ) = IVAR ( I )
               END IF
  355       CONTINUE
            NVAR2  = NVAR1
            NVAR   = NVAR2
            NVAR1  = 1
            NBPROD = 1
C
C  EVALUATE THE PRODUCT OF THE HESSIAN WITH THE DENSE VECTOR WK( LP ). 
C
            T = CPUTIM( DUM )
CS          CALL SSETVL( N, Q, 1, ZERO )
CD          CALL DSETVL( N, Q, 1, ZERO )
CS          CALL SHSPRD( N, NN, NG, NGEL, NVAR, NVAR1, NVAR2, NBPROD,
CD          CALL DHSPRD( N, NN, NG, NGEL, NVAR, NVAR1, NVAR2, NBPROD,
     *                   ALLLIN, IVAR, ISTAEV, LSTAEV, ISTADH, LSTADH, 
     *                   INTVAR, LNTVAR, IELING, LELING, IELVAR, LELVAR,
     *                   IWK( LSTAJC + 1 ), LNSTJC, IWK( LSELTS + 1 ),
     *                   LNELTS, IWK( LSPTRS + 1 ), LNPTRS, 
     *                   IWK( LGCOLJ + 1 ), LNGCLJ, IWK( LSLGRP + 1 ),
     *                   LNLGRP, IWK( LSWKSP + 1 ), LNWKSP, 
     *                   IWK( LSVGRP + 1 ), LNVGRP, IWK( LSTAGV + 1 ),
     *                   LNSTGV, IWK( LVALJR + 1 ), LNVLJR, 
     *                   NNONNZ, IWK( LNNONZ + 1 ), LNNNON,
     *                   IWK( LIUSED + 1 ), LNIUSE, IWK( LNONZ2 + 1 ),
     *                   LNNNO2, IWK( LSYMMH + 1 ), MAXSIN, 
     *                   WK( LP + 1 ), Q, GVALS( 1, 2 ),
     *                   GVALS( 1, 3 ), FUVALS( LGRJAC + 1 ),  LNGRJC,
     *                   GSCALE, ESCALE, LESCAL, FUVALS, LNHUVL, 
     *                   WK, LNWK, WK( LWKB ), LNWKB, WK( LWKC ), LNWKC,
     *                   GXEQX, LGXEQX, INTREP, LINTRE, .TRUE., RANGES )
            TMV = TMV + CPUTIM( DUM ) - T
C
C  COMPUTE THE CURVATURE.
C
            CURV     = ZERO
            DO 360 J = 1, NVAR2
               I     = IVAR( J )
               CURV  = CURV  + Q( I ) * WK( LP + I )
  360       CONTINUE
C
C  COMPARE THE CALCULATED AND RECURRED CURVATURE.
C
            WRITE( IOUT, 2550 ) CURV
            WRITE( IOUT, 2530 ) INFOR
            IF ( INFOR .EQ. 1 .OR. INFOR .EQ. 3 .OR. INFOR .EQ. 5 ) THEN
               DO 365 J = 1, NVAR2
                  I     = IVAR( J )
                  WRITE( IOUT, 2540 ) I, I, WK( LP + I ), Q( I ),
     *                                WK( LQGRAD + I )
  365          CONTINUE
            END IF
         END IF
      ELSE
C
C  - - - - - - - - - - - - ITERATIVE METHOD - - - - - - - - - - - - - -
C
C   MINIMIZE THE QUADRATIC USING AN ITERATIVE METHOD. THE METHOD USED
C   IS A SAFEGUARDED PRECONDITIONED CONJUGATE GRADIENT SCHEME.
C
         T = CPUTIM( DUM )
CS       CALL SCG( N, X, XT, FUVALS( LGGFX + 1 ), WK( LBND + 1 ),
CD       CALL DCG( N, X, XT, FUVALS( LGGFX + 1 ), WK( LBND + 1 ),
     *             IWK( LINDEX + 1 ), CGSTOP, FMODEL, VSCALE, DGRAD,
     *             WK( LBREAK + 1 ), INFORM, WK( LP + 1 ), Q, IVAR,
     *             NVAR, NVAR2, TWONRM, GMODEL, DXSQR,
     *             IOUT, JUMPTO, IPRINT )
         TLS = TLS + CPUTIM( DUM ) - T
C
C  COMPUTE THE STEP MOVED FROM THE PREVIOUS TO THE CURRENT ITERATE.
C
         IF ( JUMPTO .EQ. 0 .OR. JUMPTO .EQ. 4 .OR.
     *        JUMPTO .EQ. 5 ) THEN
CS          STEP = SDNRM( N, XT, X, TWONRM, VSCALE, .TRUE. )
CD          STEP = DDNRM( N, XT, X, TWONRM, VSCALE, .TRUE. )
         END IF
C
C  THE NORM OF THE GRADIENT OF THE QUADRATIC MODEL IS SMALLER THAN
C  CGSTOP. PERFORM ADDITIONAL TESTS TO SEE IF THE CURRENT ITERATE
C  IS ACCEPTABLE.
C
         NVAR2 = NVAR
         IF ( JUMPTO .EQ. 4 ) THEN
C
C  IF THE (SCALED) STEP TAKEN IS SMALL RELATIVE TO THE TRUST REGION
C  RADIUS, ENSURE THAT AN ACCURATE APPROXIMATION TO THE MINIMIZER OF THE
C  MODEL IS FOUND.
C
            IF ( STEP .LE. STPTOL * RADIUS .AND. .NOT. QUADRT
     *                                     .AND. .NOT. SLVBQP ) THEN
               IF ( MAX( RESMIN, STEP * CGSTOP / ( RADIUS * STPTOL ) )
     *              .GE. GMODEL ) THEN
                  IF ( IPRINT .GE. 5 ) WRITE( IOUT, 2820 ) STEP
                  JUMPTO = 0
               ELSE
                  GI = STEP * CGSTOP / ( RADIUS * STPTOL )
                  IF ( IPRINT .GE. 5 ) WRITE( IOUT, 2260)
     *                                    GI, STEP, RADIUS
                  JUMPTO = 4
               END IF
            ELSE
               JUMPTO = 0
            END IF
         END IF
C
C  A BOUND HAS BEEN ENCOUNTERED IN CG. IF THE BOUND IS A TRUST
C  REGION BOUND, STOP THE MINIMIZATION.
C
         IF ( JUMPTO .EQ. 5 ) THEN
            IFACTR = 2
            IF ( TWONRM ) THEN
               JUMPTO = 2
            ELSE
              IF ( SLVBQP ) THEN
                  JUMPTO = 2
               ELSE
                  JUMPTO = 0
C
C  THE BOUND ENCOUNTERED IS AN UPPER BOUND.
C
                  IF ( IFIXED .GT. 0 ) THEN
                     IF ( CALCDI  ) THEN
                        IF ( BU( IFIXED ) .LT. X( IFIXED ) + RADIUS
     *                     / SQRT( FUVALS( LDX + IFIXED ) ) ) JUMPTO = 2
                     ELSE
                        IF ( BU( IFIXED ) .LT. X( IFIXED )
     *                       + RADIUS * VSCALE( IFIXED ) ) JUMPTO = 2
                     END IF
                  ELSE
C
C  THE BOUND ENCOUNTERED IS A LOWER BOUND.
C
                     IF ( CALCDI ) THEN
                        IF ( BL( - IFIXED ) .GT. X( - IFIXED ) - RADIUS
     *                     / SQRT( FUVALS( LDX - IFIXED ) ) ) JUMPTO = 2
                     ELSE
                        IF ( BL( - IFIXED ) .GT. X( - IFIXED )
     *                       - RADIUS * VSCALE( - IFIXED ) ) JUMPTO = 2
                     END IF
                  END IF
               END IF
            END IF
            IF ( IPRINT .GE. 5 .AND. JUMPTO .EQ. 2 ) WRITE( IOUT, 2200 )
         END IF
C
C  IF THE BOUND ENCOUNTERED WAS A PROBLEM BOUND, CONTINUE MINIMIZING 
C  THE MODEL
C
         IF ( JUMPTO .GT. 0 ) GO TO 300
         CGEND = CGENDS( INFORM - 9 )
      END IF
C
C  IF REQUIRED, COMPUTE THE VALUE OF THE MODEL FROM FIRST PRINCIPLES.
C
      IF ( IPRINT .GE. 10 ) THEN
         NBPROD = NBPROD + 1
         NVAR1  = 1
         NVAR   = N
         NVAR2  = NVAR
C
C  COMPUTE THE STEP TAKEN, WK( LP ).
C
         DO 370 I        = 1, N
            IVAR( I )    = I
            WK( LP + I ) = XT( I ) - X( I )
  370    CONTINUE
C
C  EVALUATE THE PRODUCT OF THE HESSIAN WITH THE DENSE VECTOR WK( LP ).
C
         T = CPUTIM( DUM )
CS       CALL SSETVL( N, Q, 1, ZERO )
CD       CALL DSETVL( N, Q, 1, ZERO )
CS       CALL SHSPRD( N, NN, NG, NGEL, NVAR, NVAR1, NVAR2, NBPROD,
CD       CALL DHSPRD( N, NN, NG, NGEL, NVAR, NVAR1, NVAR2, NBPROD,
     *                ALLLIN, IVAR, ISTAEV, LSTAEV, ISTADH, LSTADH, 
     *                INTVAR, LNTVAR, IELING, LELING, IELVAR, LELVAR,
     *                IWK( LSTAJC + 1 ), LNSTJC, IWK( LSELTS + 1 ),
     *                LNELTS, IWK( LSPTRS + 1 ), LNPTRS, 
     *                IWK( LGCOLJ + 1 ), LNGCLJ, IWK( LSLGRP + 1 ),
     *                LNLGRP, IWK( LSWKSP + 1 ), LNWKSP, 
     *                IWK( LSVGRP + 1 ), LNVGRP, IWK( LSTAGV + 1 ),
     *                LNSTGV, IWK( LVALJR + 1 ), LNVLJR, 
     *                NNONNZ, IWK( LNNONZ + 1 ), LNNNON,
     *                IWK( LIUSED + 1 ), LNIUSE, IWK( LNONZ2 + 1 ),
     *                LNNNO2, IWK( LSYMMH + 1 ), MAXSIN, 
     *                WK( LP + 1 ), Q, GVALS( 1, 2 ),
     *                GVALS( 1, 3 ), FUVALS( LGRJAC + 1 ), LNGRJC,
     *                GSCALE, ESCALE, LESCAL, FUVALS, LNHUVL,  
     *                WK, LNWK, WK( LWKB ), LNWKB, WK( LWKC ), LNWKC,
     *                GXEQX, LGXEQX, INTREP, LINTRE,
     *                .TRUE., RANGES )
         TMV  = TMV + CPUTIM( DUM ) - T
C
C  IF REQUIRED, PRINT THE STEP TAKEN.
C
         IF ( IPRINT .GE. 20 ) WRITE( IOUT, 2880 )
     *             ( WK( LP + I ), I = 1, N )
C
C  COMPUTE THE MODEL VALUE, FNEW, AND RESET WK( LP ) TO ZERO.
C
         FNEW  = F
CDIR$ IVDEP
         DO 375 J = 1, NVAR2
            I     = IVAR( J )
            FNEW  = FNEW +
     *         ( FUVALS( LGGFX + I ) + HALF * Q( I ) ) * WK( LP + I )
            WK( LP + I ) = ZERO
  375    CONTINUE
         WRITE( IOUT, 2620 ) FNEW * FINDMX, FMODEL * FINDMX
      END IF
  380 CONTINUE
C
C  ---------------------------------------------------------------------
C  STEP 3.5 OF THE ALGORITHM (SEE PAPER).
C  ---------------------------------------------------------------------
C
C  AN ACCURATE APPROXIMATION TO THE MINIMUM OF THE QUADRATIC
C  MODEL IS TO BE SOUGHT.
C
      IF ( SLVBQP .AND. IBQPST .LE. 2 ) THEN
C        IBQPST = IBQPST + 1
C
C  COMPUTE THE GRADIENT VALUE.
C
         NVAR  = N
         NVAR2 = NVAR
C
C  COMPUTE THE STEP TAKEN.
C
         DO 390 I            = 1, N
            IVAR( I )        = I
            WK( LDELTX + I ) = XT( I ) - X( I )
  390    CONTINUE   
C
C  EVALUATE THE PRODUCT OF THE HESSIAN WITH THE DENSE STEP VECTOR 
C
         T = CPUTIM( DUM )
CS       CALL SSETVL( N, Q, 1, ZERO )
CD       CALL DSETVL( N, Q, 1, ZERO )
CS       CALL SHSPRD( N, NN, NG, NGEL, NVAR, NVAR1, NVAR2, NBPROD,
CD       CALL DHSPRD( N, NN, NG, NGEL, NVAR, NVAR1, NVAR2, NBPROD,
     *                ALLLIN, IVAR, ISTAEV, LSTAEV, ISTADH, LSTADH, 
     *                INTVAR, LNTVAR, IELING, LELING, IELVAR, LELVAR,
     *                IWK( LSTAJC + 1 ), LNSTJC, IWK( LSELTS + 1 ),
     *                LNELTS, IWK( LSPTRS + 1 ), LNPTRS, 
     *                IWK( LGCOLJ + 1 ), LNGCLJ, IWK( LSLGRP + 1 ),
     *                LNLGRP, IWK( LSWKSP + 1 ), LNWKSP, 
     *                IWK( LSVGRP + 1 ), LNVGRP, IWK( LSTAGV + 1 ),
     *                LNSTGV, IWK( LVALJR + 1 ), LNVLJR, 
     *                NNONNZ, IWK( LNNONZ + 1 ), LNNNON,
     *                IWK( LIUSED + 1 ), LNIUSE, IWK( LNONZ2 + 1 ),
     *                LNNNO2, IWK( LSYMMH + 1 ), MAXSIN, 
     *                WK( LDELTX + 1 ), Q, GVALS( 1, 2 ),
     *                GVALS( 1, 3 ), FUVALS( LGRJAC + 1 ), LNGRJC,
     *                GSCALE, ESCALE, LESCAL, FUVALS, LNHUVL, 
     *                WK, LNWK, WK( LWKB ), LNWKB, WK( LWKC ), LNWKC,
     *                GXEQX, LGXEQX, INTREP, LINTRE,
     *                .TRUE., RANGES )
         TMV = TMV + CPUTIM( DUM ) - T
C
C  COMPUTE THE MODEL GRADIENT AT XT.
C
         DLTNRM   = ZERO
CDIR$ IVDEP
         DO 395 I = 1, N
C           IF ( IWK( LINDEX + I ) .NE. 0 ) 
C     *           WK( LGX0 + I ) = FUVALS( LGGFX + I ) + Q( I )
            WK( LGX0 + I ) = FUVALS( LGGFX + I ) + Q( I )
            DLTNRM = MAX( DLTNRM, ABS( Q( I ) ) )
  395    CONTINUE
C
C  SAVE THE VALUES OF THE NONZERO COMPONENTS OF THE GRADIENT.
C
      K        = 0
      DO 396 J = 1, NFREEF
         I     = IWK( LSEND + J ) 
         IF ( I .GT. 0 ) THEN
            K = K + 1
            WK( LGX0 + I ) = DGRAD( K )
         END IF
  396 CONTINUE
      IF ( IPRINT .GE. 1000 ) THEN
         DO 397 I = 1, N
            WRITE( 6, * ) WK( LGX0 + I ), FUVALS( LGGFX + I ) + Q( I )
  397    CONTINUE
      END IF   
C
C  FIND THE PROJECTED GRADIENT OF THE MODEL AND ITS NORM.
C
CS       CALL SPRGRA( N, XT, WK( LGX0 + 1 ), VSCALE, WK( LBND + 1 ),
CD       CALL DPRGRA( N, XT, WK( LGX0 + 1 ), VSCALE, WK( LBND + 1 ),
     *                WK( LBND + N + 1 ), DGRAD, IVAR, NVAR, GMODEL )
C
C  CHECK FOR CONVERGENCE OF THE INNER ITERATION.
C
         IF ( IPRINT .GT. 1 ) WRITE( IOUT, 2750 ) GMODEL, SQRT( CGSTOP )
         IF ( GMODEL * GMODEL .GT. CGSTOP .AND.
     *        DLTNRM .GT. EPSMCH ) THEN
C
C  THE APPROXIMATION TO THE MINIMIZER OF THE QUADRATIC MODEL IS NOT YET
C  GOOD ENOUGH. PERFORM ANOTHER ITERATION.
C
C  STORE THE FUNCTION VALUE AT THE STARTING POINT FOR THE CAUCHY SEARCH.
C
            F0       = FMODEL
            DO 398 I = 1, N
C
C  SET THE STARING POINT FOR THE CAUCHY STEP.
C
               WK( LX0 + I ) = XT( I )
C
C  SET THE CAUCHY DIRECTION.
C
               DGRAD( I )   = - WK( LGX0 + I ) * 
     *                        ( VSCALE( I ) / VSCMAX ) ** 2
               WK( LP + I ) = ZERO
  398       CONTINUE
C
C  IF POSSIBLE, USE THE EXISTING PRECONDITIONER.
C
            IF ( REFACT ) THEN
               IFACTR = 1
            ELSE
C
C  ENSURE THAT A NEW SCHUR COMPLEMENT IS CALCULATED. RESTORE THE COMPLETE
C  LIST OF VARIABLES THAT WERE FREE WHEN THE FACTORIZATION WAS CALCULATED.
C
               IFACTR = 2
               NFIXED = 0
               DO 399 I            = 1, NFREEF
                  IWK( LSEND + I ) = ABS( IWK( LSEND + I ) )
  399          CONTINUE
            END IF   
            JUMPTO = 1
            GO TO 240
         END IF
      END IF
C
C  ---------------------------------------------------------------------
C  STEP 4 OF THE ALGORITHM (SEE PAPER).
C  ---------------------------------------------------------------------
C
C  TEST FOR ACCEPTANCE OF NEW POINT AND TRUST REGION
C  MANAGEMENT.
C
  400 CONTINUE
C
C  DETERMINE WHICH NONLINEAR ELEMENTS AND NON-TRIVIAL GROUPS NEED TO
C  BE RE-EVALUATED BY CONSIDERING WHICH OF THE VARIABLES HAVE CHANGED.
C
CS    CALL SCLCFG( UNSUCC, N, NCALCF, NCALCG,
CD    CALL DCLCFG( UNSUCC, N, NCALCF, NCALCG,
     *             ISTAEV, LSTAEV, ISTADG, LSTADG, IELING,
     *             LELING, ICALCF, LCALCF, ICALCG, LCALCG,
     *             IWK( LSPTRS + 1 ), LNPTRS, IWK( LSELTS + 1 ),
     *             LNELTS, IWK( LSTAJC + 1 ), LNSTJC,
     *             IWK( LGCOLJ + 1 ), LNGCLJ, X, XT )
C
C  IF REQUIRED, PRINT A LIST OF THE NONLINEAR ELEMENTS AND GROUPS
C  WHICH HAVE CHANGED.
C
      IF ( IPRINT .GE. 5 .AND. .NOT. ALLLIN ) WRITE( IOUT, 2590 )
     *      ( ICALCF( I ), I = 1, NCALCF )
      IF ( IPRINT .GE. 5 .AND. .NOT. ALTRIV ) WRITE( IOUT, 2720 )
     *      ( ICALCG( I ), I = 1, NCALCG )
C
C  IF THE STEP TAKEN IS RIDICULOUSLY SMALL, EXIT.
C
      IF ( STEP .LE. STPMIN ) THEN
         INFORM = 3
         GO TO 600
      END IF
C
C  RETURN TO THE CALLING PROGRAM TO OBTAIN THE FUNCTION
C  VALUE AT THE NEW POINT.
C
      INFORM = - 3
      RETURN
  430 CONTINUE
C
C  COMPUTE THE GROUP ARGUMENT VALUES FT.
C
      DO 460 IG = 1, NG
         FTT    = - B( IG )
C
C  INCLUDE THE CONTRIBUTION FROM THE LINEAR ELEMENT.
C
         DO 440 J = ISTADA( IG ), ISTADA( IG + 1 ) - 1
            FTT   = FTT + A( J ) * XT( ICNA( J ) )
  440    CONTINUE
C
C  INCLUDE THE CONTRIBUTIONS FROM THE NONLINEAR ELEMENTS.
C
         DO 450 J = ISTADG( IG ), ISTADG( IG + 1 ) - 1
            FTT   = FTT + ESCALE( J ) * FUVALS( IELING( J ) )
  450    CONTINUE
         FT( IG ) = FTT
  460 CONTINUE
C
C  COMPUTE THE GROUP FUNCTION VALUES.
C
      IF ( ALTRIV ) THEN
CS       FNEW = SDOT( NG, GSCALE, 1, FT, 1 )
CD       FNEW = DDOT( NG, GSCALE, 1, FT, 1 )
      ELSE
C
C  IF NECESSARY, RETURN TO THE CALLING PROGRAM TO OBTAIN THE GROUP
C  FUNCTION AND DERIVATIVE VALUES AT THE INITIAL POINT.
C
         INFORM = - 4
         RETURN
      END IF
  470 CONTINUE
      IF ( .NOT. ALTRIV ) THEN
         FNEW      = ZERO
CDIR$ IVDEP
         DO 480 IG = 1, NG
            IF ( GXEQX( IG ) ) THEN
               FNEW = FNEW + GSCALE( IG ) * FT( IG )
            ELSE
               FNEW = FNEW + GSCALE( IG ) * GVALS( IG, 1 )
            END IF
  480    CONTINUE
      END IF
C
C  COMPUTE THE ACTUAL AND PREDICTED REDUCTIONS IN THE FUNCTION VALUE.
C  ENSURE THAT ROUNDING ERRORS DO NOT DOMINATE.
C
      ARED   = ( F - FNEW   ) + MAX( ONE, ABS( F ) ) * TENEPS
      PRERED = ( F - FMODEL ) + MAX( ONE, ABS( F ) ) * TENEPS
      IF ( ABS( ARED ) .LT. TENEPS .AND. ABS( F ) .GT. TENEPS )
     *   ARED = PRERED
      IF ( QUADRT ) THEN
         RHO = ONE
      ELSE
         RHO = ARED / PRERED
      END IF
      IF ( IPRINT .GE. 100 ) WRITE( IOUT, 2940 ) F, FNEW, F, FMODEL
      IF ( IPRINT .GE.   3 ) WRITE( IOUT, 2070 ) ARED, PRERED, RHO
C
C  IF IN INTERACTIVE MODE, RETURN TO THE USER TO DECIDE WHETHER TO 
C  ACCEPT THE CURRENT STEP.
C
C     IF ( INTERA ) THEN
C        WRITE( IOUT, 4446 ) RHO
C        READ( 5, 4447 ) INCREA         
C        IF ( INCREA .EQ. 'H' .OR. INCREA .EQ. 'D' ) THEN
C           IF ( INCREA .EQ. 'D' ) ASTEP = 2.0D+0
C           IF ( INCREA .EQ. 'H' ) ASTEP = 5.0D-1 
C           DO 490 I   = 1, N
C              XT( I ) = X( I ) + ASTEP * ( XT( I ) - X( I ) )
C 490       CONTINUE
C           GO TO 380
C        END IF
C     END IF
C
C  -------------------------------------------------------------------
C  STEP 5 OF THE ALGORITHM (SEE PAPER).
C  -------------------------------------------------------------------
C
      OLDRAD = RADIUS
C
C  - - - - STEP MANAGEMENT WHEN THE ITERATION HAS PROVED UNSUCCESSFUL -
C
      IF ( RHO .LE. RMU .OR. PRERED .LE. ZERO ) THEN
C
C  UNSUCCESSFUL STEP. CALCULATE THE RADIUS WHICH WOULD JUST INCLUDE
C  THE NEWLY FOUND POINT, XT.
C
         UNSUCC = .TRUE.
         IF ( RHO .GE. ZERO .AND. PRERED .GT. ZERO ) THEN
            RADMIN = STEP
         ELSE
C
C  VERY UNSUCCESSFUL STEP. OBTAIN AN ESTIMATE OF THE RADIUS
C  REQUIRED TO OBTAIN A SUCCESSFUL STEP ALONG THE STEP TAKEN,
C  RADMIN, IF SUCH A STEP WERE TAKEN AT THE NEXT ITERATION.
C
            SLOPE    = ZERO
            DO 510 I = 1, N
               SLOPE = SLOPE + FUVALS( LGGFX + I ) * ( XT( I ) - X( I ))
  510       CONTINUE
            CURV   = FMODEL - F - SLOPE 
            RADMIN = STEP * ( ETA - ONE ) * SLOPE /
     *               ( FNEW - F - SLOPE - ETA * CURV )
         END IF
C
C  COMPUTE AN UPPER BOUND ON THE NEW TRUST REGION RADIUS. RADMIN, 
C  THE ACTUAL RADIUS WILL BE THE CURRENT RADIUS MULTIPLIED BY THE 
C  LARGEST POWER OF GAMMA1 FOR WHICH THE PRODUCT IS SMALLER THAN
C  RADMIN.
C
         RADMIN = MAX( RADIUS * GAMMA0, RADMIN )
C
C  IF THE TRUST REGION RADIUS HAS SHRUNK TOO MUCH, EXIT. THIS 
C  MAY INDICATE A DERIVATIVE BUG OR THAT THE USER IS ASKING FOR 
C  TOO MUCH ACCURACY IN THE FINAL GRADIENT.
C
         IF ( RADMIN .LT. RADTOL ) THEN
            IF ( IPRINT .GE. 0 ) WRITE( IOUT, 2130 )
            INFORM = 2
            GO TO 600
         END IF
C
C  CONTINUE REDUCING THE RADIUS BY THE FACTOR GAMMA1 UNTIL IT IS
C  SMALLER THAN RADMIN.
C
  520    CONTINUE
         RADIUS = GAMMA1 * RADIUS
         IF ( RADIUS .GE. RADMIN ) GO TO 520
C
C  COMPUTE THE DISTANCE OF THE GENERALIZED CAUCHY POINT FROM THE
C  INITIAL POINT.
C
         IF ( CALCDI ) THEN
            DO 530 I = 1, N
               WK( I ) = ONE / SQRT( FUVALS( LDX + I ) )
  530       CONTINUE   
CS          STEP = SDNRM( N, XT, X, TWONRM, WK, .TRUE. )
CD          STEP = DDNRM( N, XT, X, TWONRM, WK, .TRUE. )
         ELSE
CS          STEP = SDNRM( N, XT, X, TWONRM, VSCALE, .TRUE. )
CD          STEP = DDNRM( N, XT, X, TWONRM, VSCALE, .TRUE. )
         END IF
C
C  IF THE GENERALIZED CAUCHY POINT LIES WITHIN THE NEW TRUST REGION,
C  IT MAY BE REUSED.
C
C        REUSEC = STEP .LT. RADIUS
C
C  START A FURTHER ITERATION USING THE NEWLY REDUCED TRUST REGION.
C
         IF ( DIRECT ) NEXT = .FALSE.
         GO TO 100
C
C  - - - - - - - - - - - SUCCESSFUL STEP - - - - - - - - - - - - - - - -
C
      ELSE
         UNSUCC = .FALSE.
C
C  - - STEP MANAGEMENT WHEN THE ITERATION HAS PROVED VERY SUCCESSFUL - -
C
         IF ( RHO .GE. ETA ) THEN
C
C  IF IN INTERACTIVE MODE, RETURN TO THE USER TO DECIDE WHETHER TO 
C  INCREASE THE TRUST REGION RADIUS.
C
C           IF ( INTERA ) THEN
C              WRITE( IOUT, 4444 )
C              IF ( INCREA .EQ. '>' )
C    *         RADIUS = MIN( MAX( RADIUS, GAMMA2 * STEP ), RADMAX )
C           ELSE
C
C  INCREASE THE TRUST REGION RADIUS.
C  NOTE THAT WE REQUIRE THE STEP TAKEN TO BE AT LEAST A CERTAIN
C  MULTIPLE OF THE DISTANCE TO THE TRUST REGION BOUNDARY.
C
               RADIUS = MIN( MAX( RADIUS, GAMMA2 * STEP ), RADMAX )
C           END IF         
         END IF
C
C  - - DERIVATIVE EVALUATIONS WHEN THE ITERATION HAS PROVED SUCCESSFUL -
C
C  EVALUATE THE GRADIENT AND APPROXIMATE HESSIAN. FIRSTLY,
C  SAVE THE OLD ELEMENT GRADIENTS IF APPROXIMATE HESSIANS
C  ARE TO BE USED.
C
         T = CPUTIM( DUM )
CS       IF ( .NOT. SECOND .AND. .NOT. ALLLIN ) CALL SCOPY( NINVAR, 
CD       IF ( .NOT. SECOND .AND. .NOT. ALLLIN ) CALL DCOPY( NINVAR,
     *                                 FUVALS( LGXI + 1 ), 1, WK, 1 )
         TUP = TUP + CPUTIM( DUM ) - T
C
C  IF THEY ARE USED, UPDATE THE SECOND DERIVATIVE APPROXIMATIONS.
C
         IF ( .NOT. SECOND .AND. .NOT. ALLLIN ) THEN
C
C  FORM THE DIFFERENCES IN THE ITERATES, WK( LP ).
C
CDIR$ IVDEP
            DO 550 J        = 1, N
               WK( LP + J ) = XT( J ) - X( J )
  550       CONTINUE
         END IF
C
C  ACCEPT THE COMPUTED POINT AND FUNCTION VALUE.
C
         F = FNEW
CS       CALL SCOPY( N, XT, 1, X, 1 )
CD       CALL DCOPY( N, XT, 1, X, 1 )
C
C  RETURN TO THE CALLING PROGRAM TO OBTAIN THE DERIVATIVE
C  VALUES AT THE NEW POINT.
C
         IF ( FDGRAD ) IGETFD = 0
         IF ( .NOT. ( ALTRIV .AND. ALLLIN ) ) THEN
            IF ( ALTRIV ) THEN
               INFORM = - 6
            ELSE
               INFORM = - 5
            END IF
            NGEVAL = NGEVAL + 1
            RETURN
         END IF
      END IF
  540 CONTINUE
C
C  IF FINITE-DIFFERENCE GRADIENTS ARE USED, COMPUTE THEIR VALUES.
C
      IF ( FDGRAD .AND. .NOT. ALLLIN ) THEN
C
C  STORE THE VALUES OF THE NONLINEAR ELEMENTS FOR FUTURE USE.
C
         IF ( IGETFD .EQ. 0 ) THEN
CS          CALL SCOPY( NEL, FUVALS, 1, WK( LWKSTR + 1 ), 1 )
CD          CALL DCOPY( NEL, FUVALS, 1, WK( LWKSTR + 1 ), 1 )
            CENTRL = ICHOSE( 3 ) .EQ. 2 .OR. PJGNRM .LT. EPSMCH ** 0.25
         END IF
C
C  OBTAIN A FURTHER SET OF DIFFERENCES.
C
CS       CALL SFDGRD( N, NEL, IELVAR, LELVAR, ISTAEV, LSTAEV,
CD       CALL DFDGRD( N, NEL, IELVAR, LELVAR, ISTAEV, LSTAEV,
     *                IWK   ( LSELTS + 1   ), LNELTS,
     *                IWK   ( LSPTRS + 1   ), LNPTRS,
     *                IELING, LELING,
     *                IWK   ( LSSVSE + 1   ), LNSVSE,
     *                IWK   ( LSISET + 1   ), LNISET, NSETS , 
     *                IWK   ( LSEND  + 1   ), NEL   , 
     *                ICALCF, LCALCF, NCALCF, INTVAR, LNTVAR,
     *                IWK   ( LSSWTR + 1   ), LNSWTR,
     *                IWK   ( LSSIWT + 1   ), LNSIWT, 
     *                IWK   ( LSTYPE + 1   ), LNTYPE, NTYPE ,
     *                IWK   ( LSIWTR + 1   ), LNIWTR,
     *                WK    ( LSWTRA + 1   ), LNWTRA, X , XT, 
     *                FUVALS, LFUVAL, CENTRL, IGETFD )
         IF ( IGETFD .GT. 0 ) THEN
            INFORM = - 7
            RETURN
         END IF
C
C  RESTORE THE VALUES OF THE NONLINEAR ELEMENTS AT X..
C
         IGETFD = NSETS + 1
CS       CALL SCOPY( NEL, WK( LWKSTR + 1 ), 1, FUVALS, 1 )
CD       CALL DCOPY( NEL, WK( LWKSTR + 1 ), 1, FUVALS, 1 )
      END IF
C
C  COMPUTE THE GRADIENT VALUE.
C
      T = CPUTIM( DUM )
CS    CALL SELGRD( N, NG, .FALSE., ICNA, LICNA, ISTADA, LSTADA,
CD    CALL DELGRD( N, NG, .FALSE., ICNA, LICNA, ISTADA, LSTADA,
     *             IELING, LELING, ISTADG, LSTADG, ISTAEV, LSTAEV,
     *             IELVAR, LELVAR, INTVAR, LNTVAR, IWK( LSVGRP + 1 ),
     *             LNVGRP, IWK( LSTAJC + 1 ), LNSTJC, 
     *             IWK( LSTAGV + 1 ), LNSTGV, A, LA, GVALS( 1, 2 ), 
C ** Correction 5. 05/02/93: 1 line corrected **
     *             LGVALS, FUVALS, LHXI, FUVALS( LGGFX + 1 ),
C ** Correction 5. 05/02/93: end of correction **
     *             GSCALE, LGSCAL, 
     *             ESCALE, LESCAL, FUVALS( LGRJAC + 1 ),
     *             LNGRJC, WK( NINVAR + 1 ), WK( NN + 1 ), MAXSEL,
     *             GXEQX, LGXEQX, INTREP, LINTRE, RANGES )
C
C  IF THEY ARE USED, UPDATE THE SECOND DERIVATIVE APPROXIMATIONS.
C
      IF ( .NOT. SECOND .AND. .NOT. ALLLIN ) THEN
C
C  FORM THE DIFFERENCES IN THE GRADIENTS, WK( ).
C
CDIR$ IVDEP
         DO 560 J   = 1, NINVAR
            WK( J ) = FUVALS( LGXI + J ) - WK( J )
  560    CONTINUE
         IF ( FIRSUP ) THEN
C
C  IF A SECANT METHOD IS TO BE USED, SCALE THE INITIAL SECOND
C  DERIVATIVE MATRIX FOR EACH ELEMENT SO AS TO SATISFY
C  THE WEAK SECANT CONDITION.
C
CS       CALL SSCALH( .FALSE., N, NEL, NCALCF, LHXI, ISTAEV, LSTAEV,
CD       CALL DSCALH( .FALSE., N, NEL, NCALCF, LHXI, ISTAEV, LSTAEV,
     *                ISTADH, LSTADH, ICALCF, LCALCF, 
     *                INTVAR, LNTVAR, IELVAR, LELVAR, IWK( LSYMMD + 1 ),
     *                MAXSIN, INTREP, LINTRE, FUVALS, LFUVAL, 
     *                WK( NINVAR + 1 ), LNWK, WK( NN + 1 ), LNWKB,
     *                WK( LP + 1 ), WK, NINVAR, RANGES )
            FIRSUP = .FALSE.
         END IF
C
C  UPDATE THE SECOND DERIVATIVE APPROXIMATIONS USING ONE OF FOUR
C  POSSIBLE SECANT UPDATING FORMULAE, BFGS, DFP, PSB AND SR1.
C
CS       CALL SSECNT( N, IELVAR, LELVAR, ISTAEV, LSTAEV, INTVAR,
CD       CALL DSECNT( N, IELVAR, LELVAR, ISTAEV, LSTAEV, INTVAR,
     *                LNTVAR, INTREP, LINTRE, ISTADH, LSTADH, NEL,
     *                NINVAR, FUVALS, LFUVAL, ICALCF, LCALCF, NCALCF,
     *                WK( LP + 1 ), WK, WK( NINVAR + 1 ), 2 * MAXSEL,
     *                ICHOSE( 4 ), IPRINT, IOUT, RANGES )
      END IF
      TUP = TUP + CPUTIM( DUM ) - T
C
C  COMPUTE THE PROJECTED GRADIENT AND ITS NORM.
C
CS    CALL SPRGRA( N, X, FUVALS( LGGFX + 1 ), VSCALE, BL, BU, DGRAD,
CD    CALL DPRGRA( N, X, FUVALS( LGGFX + 1 ), VSCALE, BL, BU, DGRAD,
     *             IVAR, NVAR, PJGNRM )
      NFREE = NVAR
C
C  IF REQUIRED, USE THE USERS PRECONDITIONER.
C
      IF ( PRCOND .AND. NVAR .GT. 0 .AND. MYPREC ) THEN
         INFORM = - 10
         RETURN
      END IF
  570 CONTINUE
C
C  FIND THE NORM OF THE 'PRECONDITIONED' PROJECTED GRADIENT. ALSO,
C  IF REQUIRED, FIND THE DIAGONAL ELEMENTS OF THE ASSEMBLED HESSIAN.
C
CS    CALL SGTPGR( N, NG, NGEL, NN, NVAR, MAXSEL, QGNORM, SMALLH,
CD    CALL DGTPGR( N, NG, NGEL, NN, NVAR, MAXSEL, QGNORM, SMALLH,
     *             PJGNRM, CALCDI, DPRCND, MYPREC, IVAR, ISTADH, LSTADH,
     *             ISTAEV, LSTAEV, IELVAR, LELVAR, INTVAR, LNTVAR,
     *             IELING, LELING, IWK( LSLGRP + 1 ), LNLGRP,
     *             IWK( LGCOLJ + 1 ), LNGCLJ, IWK( LSTAGV + 1 ),
     *             LNSTGV, IWK( LSVGRP + 1 ), LNVGRP, IWK( LVALJR + 1 ),
     *             LNVLJR, IWK( LSYMMD + 1 ), IWK( LSYMMH + 1 ), MAXSIN,
     *             DGRAD, Q, GVALS( 1, 2 ), GVALS( 1, 3 ), 
     *             FUVALS( LDX + 1 ), GSCALE, ESCALE, LESCAL,
     *             FUVALS( LGRJAC + 1 ), LNGRJC, FUVALS, LNHUVL,
     *             WK, LNWK, WK( LWKB ), LNWKB, WK( LWKC ), LNWKC,
     *             GXEQX, LGXEQX, INTREP, LINTRE, RANGES ) 
C
C  FOR TESTING PURPOSES, CHECK IF THE ACTIVE SET HAS CHANGED.
C
C     IF ( IPRINT .GE. 0 ) THEN
C        SAMEAS   = .TRUE.
C        DO 580 I = 1, N
C           IS    = 0
C ** Correction 11. 03/12/98: 7 lines corrected **
C           XI = X( I )
C           BLI = BL( I )
C           BUI = BU( I )
C           IF ( XI .LE. BLI * ( ONE + SIGN( EPSTLP, BLI ) ) ) IS = 1
C           IF ( XI .GE. BUI * ( ONE - SIGN( EPSTLN, BUI ) ) ) IS = 2
C           IF ( BUI * ( ONE - SIGN( EPSTLN, BUI ) ) .LE. 
C    *           BLI * ( ONE + SIGN( EPSTLP, BLI ) ) ) IS = 3
C ** Correction 11. 03/12/98: End of correction **
C           IF ( ISTATE( I ) .NE. IS ) SAMEAS = .FALSE.
C           ISTATE( I ) = IS
C 580    CONTINUE
C        IF ( .NOT. SAMEAS ) ICORAS = ITER
C     END IF
      IF ( DIRECT ) NEXT = STEP .LT. TENTEN * EPSMCH .AND. INFOR .EQ. 2
      GO TO 100
C
C  IF THE USER'S COMPUTED GROUP FUNCTION VALUES ARE INADEQUATE, REDUCE
C  THE TRUST REGION RADIUS.
C
  590 CONTINUE
      UNSUCC = .TRUE.
      OLDRAD = RADIUS
      RADIUS = GAMMA1 * RADIUS
C
C  IF THE TRUST REGION RADIUS HAS SHRUNK TOO MUCH, EXIT. THIS 
C  MAY INDICATE A DERIVATIVE BUG OR THAT THE USER IS ASKING FOR 
C  TOO MUCH ACCURACY IN THE FINAL GRADIENT.
C
      IF ( RADIUS .LT. RADTOL ) THEN
         IF ( IPRINT .GE. 0 ) WRITE( IOUT, 2130 )
         INFORM = 2
         GO TO 600
      END IF
      GO TO 100
C
C ---------------------------------------------------------------------
C   END THE MAIN LOOP.
C ---------------------------------------------------------------------
C
  600 CONTINUE
C
C  RECORD THE VALUE OF THE PROJECTED GRADIENT ON EXIT.
C
      STOPG = PJGNRM
C
C  PRINT DETAILS OF THE SOLUTION.
C
      IF ( IPRINT .GT. 0 ) THEN
         IF ( ITER .GT. MAXIT ) ITER = MAXIT
         IF ( ITER .EQ. 0 ) THEN
            WRITE( IOUT, 2830 ) ITER, F * FINDMX, NGEVAL, PJGNRM, 
     *                 ITERCG,          ISKIP
C    *               , ICORAS
         ELSE
            WRITE( IOUT, 2170 ) ITER, F * FINDMX, NGEVAL, PJGNRM, 
     *                  ITERCG, OLDRAD, ISKIP
C    *                , ICORAS
         END IF
         K        = 0
         DO 610 I = 1, N
C ** Correction 12. 03/12/98: 2 lines changed to 5 **
            XI = X( I )
            BLI = BL( I )
            BUI = BU( I )
            IF ( XI .LE. BLI * ( ONE + SIGN( EPSTLP, BLI ) ) .OR.
     *           XI .GE. BUI * ( ONE - SIGN( EPSTLN, BUI ) ) ) K = K + 1
C ** Correction 12. 03/12/98: end of correction **
  610    CONTINUE
         WRITE( IOUT, 2040 ) N, K
         IF ( IPRINT .GE. 3 ) THEN
            WRITE( IOUT, 2050 ) ( X( KK ), KK = 1, N )
            WRITE( IOUT, 2060 ) ( FUVALS( LGGFX + KK ) * FINDMX,
     *                            KK = 1, N )
         END IF
         WRITE( IOUT, 2120 ) TCA, TLS, TMV, TUP
         IF ( XACTCP ) THEN
            WRITE( IOUT, 2840 )
         ELSE
            WRITE( IOUT, 2850 )
         END IF
         IF ( SLVBQP ) WRITE( IOUT, 2950 )
         IF ( ICHOSE( 2 ) .EQ.  1 ) WRITE( IOUT, 2630 )
         IF ( ICHOSE( 2 ) .EQ.  2 ) WRITE( IOUT, 2640 )
         IF ( ICHOSE( 2 ) .EQ.  3 ) WRITE( IOUT, 2650 )
         IF ( ICHOSE( 2 ) .EQ.  4 ) WRITE( IOUT, 2660 )
         IF ( ICHOSE( 2 ) .EQ.  5 ) WRITE( IOUT, 2700 )
         IF ( ICHOSE( 2 ) .EQ.  6 ) WRITE( IOUT, 2730 )
         IF ( ICHOSE( 2 ) .EQ.  7 ) WRITE( IOUT, 2740 )
         IF ( ICHOSE( 2 ) .EQ.  8 ) WRITE( IOUT, 2890 ) NSEMIB
         IF ( ICHOSE( 2 ) .EQ. 11 ) WRITE( IOUT, 2670 )
         IF ( ICHOSE( 2 ) .EQ. 12 ) WRITE( IOUT, 2680 )
         IF ( DIRECT ) THEN
            IF ( MODCHL ) THEN
               WRITE( IOUT, 2560 ) ISYS( 1 ), ISYS( 5 ), FILL
            ELSE
               WRITE( IOUT, 2570 ) ( ISYS( I ), I = 1, 4 ), FILL
            END IF
         END IF
         IF ( TWONRM ) THEN
            WRITE( IOUT, 2150 )
         ELSE
            WRITE( IOUT, 2160 )
         END IF
         IF ( ICHOSE( 3 ) .GE. 1 ) WRITE( IOUT, 2930 )
         IF ( ICHOSE( 4 ) .EQ. 0 ) WRITE( IOUT, 2210 )
         IF ( ICHOSE( 4 ) .EQ. 1 ) WRITE( IOUT, 2220 )
         IF ( ICHOSE( 4 ) .EQ. 2 ) WRITE( IOUT, 2230 )
         IF ( ICHOSE( 4 ) .EQ. 3 ) WRITE( IOUT, 2240 )
         IF ( ICHOSE( 4 ) .EQ. 4 ) WRITE( IOUT, 2250 )
      END IF
  700 CONTINUE
      RETURN
C
C  NON-EXECUTABLE STATEMENTS.
C
 2000 FORMAT( /, '  Iter #g.ev c.g.it      f    proj.g    rho   ',
     *           ' radius   step  cgend #free   time' )
 2010 FORMAT( 3I6, 1P, D10.2, D8.1, D9.1, 2D8.1, A6, I6, 0P, A7 )
 2020 FORMAT( /, ' There are ', I8, ' variables',
     *        /, ' There are ', I8, ' groups',
     *        /, ' There are ', I8, ' nonlinear elements ' )
 2030 FORMAT( / )
 2040 FORMAT( /, ' There are ', I6, ' variables and ', I6,
     *           ' active bounds ' )
 2050 FORMAT( /, ' X = ', / ( 1P, 6D12.4 ) )
 2060 FORMAT( /, ' G = ', / ( 1P, 6D12.4 ) )
 2070 FORMAT( /, ' Actual change    = ', 1P, D20.12, /
     *           ' Predicted change = ', 1P, D20.12, /
     *           ' Ratio (rho)      = ', 1P, D20.12 )
 2080 FORMAT( /, ' Diagonals of second derivatives ', / ( 1P, 6D12.4 ) )
 2090 FORMAT( /, ' SBMIN : maximum number of iterations reached ' )
 2100 FORMAT(    ' The variable number ',I3, ' is at its upper bound ' )
 2110 FORMAT(    ' The variable number ',I3, ' is at its lower bound ' )
 2120 FORMAT( /, ' Times for Cauchy, systems, products and',
     *           ' updates ', 0P, 4F8.2 )
 2130 FORMAT( /, ' SBMIN : trust region radius too small ' )
 2140 FORMAT( /, ' A must have dimension at least ', I8, /,
     *           ' ------------- At present LA is ', I8 )
 2150 FORMAT( /, ' Two-norm trust region used ' )
 2160 FORMAT( /, ' Infinity-norm trust region used ' )
 2170 FORMAT( /, ' Iteration number           ', I5,
     *           '  Merit function value    = ', 1P, D19.11,
     *        /, ' No. derivative evaluations ', I5,
     *           '  Projected gradient norm = ', 1P, D19.11,
     *        /, ' C.G. iterations       ', I10,
     *           '  Trust region radius     = ', 1P, D19.11,
     *        /, ' Number of updates skipped  ', I5
     *        )
C    *         ,  '  Correct active set after', I5, ' iteration(s) ' )
 2180 FORMAT( /, ' Element hessians ', / ( 1P, 6D12.4 ) )
 2190 FORMAT( ' *** Calculated quadratic at CP ', 1P, D22.14, /,
     *        ' *** Recurred   quadratic at CP ', 1P, D22.14 )
 2200 FORMAT( /, ' Restarting the conjugate gradient iteration ' )
 2210 FORMAT( /, ' Exact second derivatives used ' )
 2220 FORMAT( /, ' B.F.G.S. approximation to second derivatives used ' )
 2230 FORMAT( /, ' D.F.P. approximation to second derivatives used ' )
 2240 FORMAT( /, ' P.S.B. approximation to second derivatives used ' )
 2250 FORMAT( /, ' S.R.1 Approximation to second derivatives used ' )
 2260 FORMAT( /, ' C.G. tolerance of ', 1P, D12.4, ' has not been',
     *           ' achieved. ', /, ' Actual step length = ', 1P, D12.4,
     *           ' Radius = ', 1P, D12.4 )
 2270 FORMAT( /, ' Element values ', / ( 1P, 6D12.4 ) )
 2280 FORMAT( /, ' Group values ', / ( 1P, 6D12.4 ) )
 2290 FORMAT( /, ' ICNA must have dimension at least ', I8, /,
     *           ' ------------- At present LICNA is ', I8 )
 2300 FORMAT( /, ' Reusing previous generalized Cauchy point ' )
 2310 FORMAT( /, ' ------- Group information ------ ' )
 2320 FORMAT( /, ' Group ', I5, ' contains ', I5, ' nonlinear',
     *           ' element(s). These are element(s)', 2I5,
     *        /, ( 16I5 ) )
 2330 FORMAT( /, ' Group ', I5, ' contains     no nonlinear',
     *           ' elements. ' )
 2340 FORMAT( '  * The group has a linear element with variable(s)',
     *           ' X( I ), I =', 3I5, /, ( 3X, 19I5 ) )
 2350 FORMAT( /, ' ------ Element information ----- ' )
 2360 FORMAT( /, ' Nonlinear element', I5, ' has ', I4, ' internal',
     *           ' and ', I4, ' elemental variable(s),',
     *        /, ' X( I ), I =   ', 13I5, /, ( 16I5 ) )
 2370 FORMAT( /, ' Nonlinear element', I5, ' has   no internal',
     *           ' or       elemental variables.' )
 2380 FORMAT( /, ' Group ', I5, ' contains ', I5, ' nonlinear',
     *           ' element.  This  is  element ', I5 )
 2390 FORMAT( '  * The group function is non-trivial ' )
 2400 FORMAT( /, ' INTVAR must have dimension at least ', I8, /,
     *           ' -------------- At present LNTVAR is ', I8 )
 2410 FORMAT( /, ' ISTADH must have dimension at least ', I8, /,
     *           ' -------------- At present LSTADH is ', I8 )
 2420 FORMAT( /, ' ICALCG must have dimension at least ', I8, /,
     *           ' -------------- At present LCALCG is ', I8 )
 2450 FORMAT( /, ' ISTADG must have dimension at least ', I8, /,
     *           ' -------------- At present LSTADG is ', I8 )
 2460 FORMAT( /, ' ISTAEV must have dimension at least ', I8, /,
     *           ' -------------- At present LSTAEV is ', I8 )
 2470 FORMAT( /, ' ISTADA must have dimension at least ', I8, /,
     *           ' -------------- At present LSTADA is ', I8 )
 2480 FORMAT( /, ' GXEQX must have dimension at least ', I8, /,
     *           ' ------------- At present LGXEQX is ', I8 )
 2490 FORMAT( /, ' INTREP must have dimension at least ', I8, /,
     *           ' -------------- At present LINTRE is ', I8 )
 2500 FORMAT( /, ' GVALS must have leading dimension at least ', I8, /,
     *           ' -------------- At present LGVALS is ', I8 )
 2510 FORMAT( /, ' FT must have dimension at least ', I8, /,
     *           ' ------------- At present LFT is ', I8 )
 2520 FORMAT( /, ' ICALCF must have dimension at least ', I8, /,
     *           ' -------------- At present LCALCF is ', I8 )
 2530 FORMAT( ' FRONTL - INFOR = ', I1 )
 2540 FORMAT( ' P, H * P(',I6,'), RHS(',I6,') = ', 1P, 3D15.7 )
 2550 FORMAT( ' CURV  = ', 1P, D12.4 )
 2560 FORMAT( ' No. pos. def. systems = ', I4, ' No. indef. systems = ',
     *        I4, / ' Ratio ( fill-in ) = ', 1P, D11.2 )
 2570 FORMAT( ' PD = ', I4, ' UD = ', I4, ' SC = ', I4, ' SI = ', I4,
     *        / ' Ratio ( fill-in ) = ', 1P, D11.2 )
 2580 FORMAT( /, ' The matrix-vector product used elements',
     *         ' marked ', I5, ' in the following list ', /, ( 20I4 ) )
 2590 FORMAT( /, ' Functions for the following elements need to be',
     *           ' re-evaluated ', /, ( 12I6 ) )
 2600 FORMAT( /, '  Iter #g.ev fill-in     f    proj.g    rho   ',
     *           ' radius   step  lsend #free   time' )
 2610 FORMAT( 2I6, 0P, F6.1, 1P, D10.2, D8.1, D9.1,
     *        2D8.1, A6, I6, 0P, A7 )
 2620 FORMAT( ' *** Calculated quadratic at end CG ', 1P, D22.14, /,
     *        ' *** Recurred   quadratic at end CG ', 1P, D22.14 )
 2630 FORMAT( /, ' Conjugate gradients without preconditioner used ' )
 2640 FORMAT( /, ' Conjugate gradients with diagonal preconditioner',
     *           ' used ' )
 2650 FORMAT( /, ' Conjugate gradients with user-supplied',
     *           ' preconditioner used ' )
 2660 FORMAT( /, ' Conjugate gradients with band inverse',
     *           ' preconditioner used ' )
 2670 FORMAT( /, ' Exact matrix factorization used ' )
 2680 FORMAT( /, ' Modified matrix factorization used ' )
 2690 FORMAT( /, ' B must have dimension at least ', I8, /,
     *           ' ------------- At present LB IS ', I8 )
 2700 FORMAT( /, ' Conjugate gradients with Munksgaards',
     *           ' preconditioner used ' )
 2710 FORMAT( ' Nonzeros of Hessian * P are in positions ', /, ( 24I3 ))
 2720 FORMAT( /, ' Functions for the following groups need to be',
     *           ' re-evaluated ', /, ( 12I6 ) )
 2730 FORMAT( /, ' Conjugate gradients with Schnabel-Eskow ',
     *           ' modified Cholesky preconditioner used ' )
 2740 FORMAT( /, ' Conjugate gradients with GMPS modified Cholesky',
     *           ' preconditioner used ' )
 2750 FORMAT( /, '    ** Model gradient is ', 1P, D12.4,
     *        ' Required accuracy is ', 1P, D12.4 )
 2760 FORMAT( /, ' Required gradient accuracy ', 1P, D8.1 )
 2770 FORMAT(    ' The variable number ',I3, ' is temporarily fixed ' )
 2780 FORMAT( /, ' Group derivatives ', / ( 1P, 6D12.4 ) )
 2790 FORMAT( /, ' Element gradients ', / ( 1P, 6D12.4 ) )
 2800 FORMAT( 3I6, 1P, D10.2, D8.1,  '     -       -       -   ',
     *        '   -  ', I6, 0P, A7 )
 2810 FORMAT( 2I6, 6X, 1P, D10.2, D8.1,  '     -       -       -   ',
     *        '   -  ', I6, 0P, A7 )
 2820 FORMAT( ' Norm of trial step ', 1P, D12.4 )
 2830 FORMAT( /, ' Iteration number           ', I5,
     *           '  Merit function value    = ', 1P, D19.11,
     *        /, ' No. derivative evaluations ', I5,
     *           '  Projected gradient norm = ', 1P, D19.11,
     *        /, ' C.G. iterations       ', I10,
     *        /, ' Number of updates skipped  ', I5
     *        )
C    *         , '  Correct active set after', I5, ' iteration(s) ' )
 2840 FORMAT( /, ' Exact Cauchy step computed ' )
 2850 FORMAT( /, ' Approximate Cauchy step computed ' )
 2860 FORMAT( 3I6, 1P, D10.2, D8.1, '     -   ',
     *        2D8.1, A6, I6, 0P, A7 )
 2870 FORMAT( 2I6, 0P, F6.1, 1P, D10.2, D8.1, '     -   ',
     *        2D8.1, A6, I6, 0P, A7 )
 2880 FORMAT( /, ' Change in X = ', / ( 1P, 6D12.4 ) )
 2890 FORMAT( /, ' Bandsolver preconditioned C.G. used', 
     *           ' (semi-bandwidth = ', I6, ') ' )
 2900 FORMAT( /, ' Two-norm of step to Cauchy point = ', 1P, D12.4 )
 2901 FORMAT( ' IWK( LSPTRS ) = ', /, ( 16I5 ) )
 2902 FORMAT( ' IWK( LSELTS ) = ', /, ( 16I5 ) )
 2904 FORMAT( ' IWK( LINDEX ) = ', /, ( 16I5 ) )
 2905 FORMAT( ' IWK( LSWKSP ) = ', /, ( 16I5 ) )
 2906 FORMAT( ' IWK( LSTAGV ) = ', /, ( 16I5 ) )
 2907 FORMAT( ' IWK( LSTAJC ) = ', /, ( 16I5 ) )
 2908 FORMAT( ' IWK( LIUSED ) = ', /, ( 16I5 ) )
 2909 FORMAT( ' IWK( LNNONZ ) = ', /, ( 16I5 ) )
 2910 FORMAT( ' IWK( LNONZ2 ) = ', /, ( 16I5 ) )
 2911 FORMAT( ' IWK( LSYMMD ) = ', /, ( 16I5 ) )
 2912 FORMAT( ' IWK( LSYMMH ) = ', /, ( 16I5 ) )
 2913 FORMAT( ' IWK( LSLGRP ) = ', /, ( 16I5 ) )
 2914 FORMAT( ' IWK( LSVGRP ) = ', /, ( 16I5 ) )
 2915 FORMAT( ' IWK( LGCOLJ ) = ', /, ( 16I5 ) )
 2916 FORMAT( ' IWK( LVALJR ) = ', /, ( 16I5 ) )
 2917 FORMAT( ' IWK( LFREEC ) = ', /, ( 16I5 ) )
 2918 FORMAT( ' IWK( LSSWTR ) = ', /, ( 16I5 ) )
 2919 FORMAT( ' IWK( LSSIWT ) = ', /, ( 16I5 ) )
 2920 FORMAT( ' IWK( LSIWTR ) = ', /, ( 16I5 ) )
 2921 FORMAT( ' IWK( LSISET ) = ', /, ( 16I5 ) )
 2922 FORMAT( ' IWK( LSSVSE ) = ', /, ( 16I5 ) )
 2923 FORMAT( ' IWK( LSTYPE ) = ', /, ( 16I5 ) )
 2930 FORMAT( /, ' Finite-difference approximations to',  
     *           ' nonlinear-element gradients used' )
 2940 FORMAT( /, ' Old F = ', 1P, D20.12, ' New   F = ', 1P, D20.12,
     *        /, ' Old F = ', 1P, D20.12, ' Model F = ', 1P, D20.12 )
 2950 FORMAT( /, ' Accuarate solution of BQP computed ' )
C ** Correction 2. 22/07/92: 3 lines added **
 2960 FORMAT( /, ' Lower bound ', 1P, D12.4, ' on variable ', I8,
     *           ' larger than upper bound ', 1P, D12.4, //,  
     *           ' Execution terminating ' )
C ** Correction 2. 22/07/92: end of correction **
C4444 FORMAT( ' V successful step: increase(>) or leave(=) radius?')
C              READ(5,4445) INCREA         
C4445 FORMAT( A1 )
C4446 FORMAT( ' Rho = ', 1P, D10.2, ' Half(H), Double(D) or',
C    *        ' Accept(A) step' )
C4447 FORMAT( A1 )
C
C  END OF SUBROUTINE SBMIN.
C
      END
CS    BLOCK DATA SDATSB
CD    BLOCK DATA DDATSB
C
C  DEFINE INITIAL VALUES FOR COMMON VARIABLES
C  FOR THE SUBROUTINE SBMIN VIA BLOCK DATA.
C
C  NICK GOULD, 15TH MARCH 1990.
C  FOR CGT PRODUCTIONS.
C
      INTEGER           ITERCG, ITCGMX, NGEVAL, ISKIP , IFIXED, NSEMIB
CS    REAL              ACCCG , RATIO , RADIUS, RADMAX, FINDMX, CMA31
CD    DOUBLE PRECISION  ACCCG , RATIO , RADIUS, RADMAX, FINDMX, CMA31
CS    COMMON / SCOMSB / ACCCG , RATIO , RADIUS, RADMAX, FINDMX, CMA31 ,
CD    COMMON / DCOMSB / ACCCG , RATIO , RADIUS, RADMAX, FINDMX, CMA31 ,
     *                  ITERCG, ITCGMX, NGEVAL, ISKIP , IFIXED, NSEMIB
CS    SAVE / SCOMSB /
CD    SAVE / DCOMSB /
CS    DATA FINDMX / 1.0E+0 /
CD    DATA FINDMX / 1.0D+0 /
CS    DATA ACCCG / 1.0E-2 /, RADIUS / - 1.0E+0 /, RADMAX / 1.0E+20 /
CD    DATA ACCCG / 1.0D-2 /, RADIUS / - 1.0D+0 /, RADMAX / 1.0D+20 /
      END
