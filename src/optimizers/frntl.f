C ** Correction report.
C ** Correction 1. 13/01/94: 2 lines corrected **
C ** Correction 2. 16/07/94: 1 line corrected **
C ** Correction 3. 16/07/94: 1 line corrected **
C ** Correction 4. 05/11/97: 1 line corrected **
C ** Correction 5. 19/08/98: added final argument to call to ASMBL.
C ** End of Correction report.
C  THIS VERSION: 19/08/1998 AT 18:00:51 PM.
CS    SUBROUTINE SFRNTL( N , NG, MAXSEL, INTVAR, LNTVAR, IELVAR, LELVAR, 
CD    SUBROUTINE DFRNTL( N , NG, MAXSEL, INTVAR, LNTVAR, IELVAR, LELVAR,
     *                   INTREP, LINTRE, IELING, LELING,
     *                   ISTADG, LSTADG, ISTAEV, LSTAEV, ISTAGV, LNSTGV,
     *                   A , LA, ICNA  , LICNA , ISTADA, LSTADA,
     *                   GUVALS, LNGUVL, ISVGRP, LNVGRP, HUVALS, LNHUVL,
     *                   ISTADH, LSTADH, GXEQX , LGXEQX,
     *                   GVALS2, GVALS3, IFREE , NFREE , GMODEL, P     ,
     *                   XT    , BND   , FMODEL, GSCALE, ESCALE, LESCAL, 
     *                   GSTOP , NUMBER, NEXT  , MODCHL, IWK   , LIWK  ,
     *                   WK, LWK, RANGES, IPRINT, MP, INFORM )
      INTEGER           N,      NG,     MAXSEL, NFREE,   MP   , NUMBER
      INTEGER           LSTADH, LICNA,  LSTADA, LNTVAR, LELVAR, LELING
      INTEGER           LSTADG, LSTAEV, LNSTGV, LNVGRP, LIWK
      INTEGER           LNGUVL, LNHUVL, LESCAL, LWK,    LGXEQX, LINTRE
      INTEGER           INFORM, IPRINT, LA
CS    REAL              FMODEL, GSTOP
CD    DOUBLE PRECISION  FMODEL, GSTOP
      LOGICAL           MODCHL, NEXT
      INTEGER           IFREE( N ), ISTADH( LSTADH ), ICNA( LICNA )  
      INTEGER           ISTADA( LSTADA ), INTVAR( LNTVAR )
      INTEGER           IELVAR( LELVAR ), IELING( LELING )
      INTEGER           ISTADG( LSTADG ), ISTAEV( LSTAEV )
      INTEGER           ISTAGV( LNSTGV ), ISVGRP( LNVGRP ), IWK( LIWK )
      LOGICAL           GXEQX( LGXEQX ), INTREP( LINTRE )
CS    REAL              A( LA ), GVALS2( NG ), GVALS3( NG ),
CD    DOUBLE PRECISION  A( LA ), GVALS2( NG ), GVALS3( NG ),
     *                  GUVALS( LNGUVL ), HUVALS( LNHUVL ), GSCALE( NG), 
     *                  ESCALE( LESCAL ), XT( N ), WK( LWK ),
     *                  GMODEL( N ), P( N ), BND( N, 2 )
      EXTERNAL RANGES
C
C     ******************************************************************
C
C     PROGRAMMING :    M. LESCRENIER ( SEP 1987 ).
C     =============    N. GOULD      ( FEB 1989 ).
C
C     DESCRIPTION :
C     =============
C
C     LET US CALL H, THE ASSEMBLED HESSIAN OF THE GROUP PARTIALLY
C     SEPARABLE FUNCTION, RESTRICTED TO THE ONLY IFREE VARIABLES,
C
C     IF H IS POSITIVE DEFINITE,
C        SOLVE H*P = GMODEL,
C
C     IF H IS INDEFINITE,
C        FIND A DIRECTION P SUCH THAT
C                    T                       T
C                 P * H * P < 0  AND  P * GMODEL > 0
C
C     IF H IS SINGULAR BUT H*P = GMODEL IS CONSISTENT
C        SOLVE H*P = GMODEL,
C
C     IF H IS SINGULAR BUT H*P = GMODEL IS INCONSISTENT
C        FIND A DESCENT DIRECTION SOL SUCH THAT
C                                          T
C                             H*P = 0 AND P * GMODEL > 0
C
C     PARAMETERS :
C     ============
C
C     N
C            I: NUMBER OF VARIABLES.
C            O: UNMODIFIED.
C
C     NG
C            I: NUMBER OF GROUPS.
C            O: UNMODIFIED.
C
C     NFREE
C            I: NUMBER OF IFREE VARIABLES (RECORDED IN IFREE).
C            O: UNMODIFIED.
C
C     MAXSEL
C            I: MAXIMUM NUMBER OF ELEMENTAL VARIABLES OF AN ELEMENT.
C            O: UNMODIFIED.
C
C     IFREE
C            I: THE FREE VARIABLES.
C            O: UNMODIFIED.
C
C     ISTAEV
C            I: ISTAEV(I) POINTS TO THE BEGIN OF THE LIST OF VARIABLES
C               OF THE I-TH ELEMENT FUNCTION (STORED IN IFREE)
C            O: UNMODIFIED.
C
C     IFREE
C            I: LISTS OF THE ELEMENTAL VARIABLES OF THE ELEMENTS.
C            O: UNMODIFIED.
C
C     HUVALS
C            I: ELEMENT HESSIANS STORED IN INTERNAL REPRESENTATION
C               (THE ONLY LOWER PARTS ARE STORED BY ROWS).
C            O: UNMODIFIED.
C
C     ISTADH
C            I: ISTADH(I)-ISTADH(1)+1 POINTS TO THE I-TH ELEMENT
C               HESSIAN IN HUVALS.
C            O: UNMODIFIED.
C
C     GMODEL
C            I: RIGHT HAND SIDE OF THE LINEAR SYSTEM
C            O: UNMODIFIED.
C
C     P
C            I: UNDEFINED.
C            O: ONLY THE IFREE VARIABLES ARE SET,
C               THE OTHER REMAINS UNMODIFIED.
C
C               IF H IS POSITIVE DEFINITE,
C                  SOLVE H*P = GMODEL,
C
C               IF H IS INDEFINITE,
C                  FIND A DIRECTION P SUCH THAT
C                    T                  T
C                  P * H * P < 0  AND  P * GMODEL > 0
C
C               IF H IS SINGULAR BUT H*P = GMODEL IS CONSISTENT
C                  SOLVE H*P = GMODEL,
C
C               IF H IS SINGULAR BUT H*P = GMODEL IS INCONSISTENT
C                  FIND A DESCENT DIRECTION P SUCH THAT
C                               T
C                  H*P = 0 AND P * GMODEL > 0
C
C     STHS
C            I: UNDEFINED
C                              T
C            O: THE CURVATURE P * H * P
C
C     INFORM
C            I: UNDEFINED.
C            O:  1   IF H, THE HESSIAN RESTRICTED TO THE IFREE
C                    VARIABLES, IS POSITIVE DEFINITE
C                2   IF H IF INDEFINITE,
C                3   IF H IS SINGULAR BUT H*P = GMODEL CONSISTENT
C                4   IF H IS SINGULAR BUT H*P = GMODEL INCONSISTENT
C
C     WK
C            I: WORK VECTOR.
C            O: MEANINGLESS.
C
C     LWK
C            I: LENGTH OF WK (GREATER THAN THE NUMBER OF NON-ZEROS
C               IN THE FACTOR OF H + N + 3 * MAXSEL).
C            O: UNMODIFIED.
C
C     IWK
C            I: WORK VECTOR.
C            O: MEANINGLESS.
C
C     LIWK
C            I: LENGTH OF IWK.
C            O: UNMODIFIED.
C
C     NUMBER, NEXT : SEE THE ROUTINE -NEGCUR-
C
C  PARTITION OF INTEGER WORKSPACE ARRAY IWK
C  ----------------------------------------
C
C   <-                        LIWK                                 ->
C  -------------------------------------------------------------------
C  |<-NNZH->|<-NNZH->]<-3*NFREE->|<-2*NFREE->| <-NIRBDU-> |          |
C  -------------------------------------------------------------------
C   |        |        |           |           |            |
C  LIRNH    LJCNH    LIKEEP      LIW1        LIMA27
C
C  PARTITION OF REAL WORKSPACE ARRAY WK
C  ------------------------------------
C
C   <-                          LWK                                ->
C  -------------------------------------------------------------------
C  |   <- NRLBDU  - >  |               |<-NFREE->|<-NFREE->|<-NFREE->|
C  -------------------------------------------------------------------
C   |                   |               |         |         |
C  LWMA27                              LRHS      LPERT     LIW
C
C     ******************************************************************
C
      INTEGER I, IFLAG, LH, LIKEEP, LIMA27, LIRNH, LIWKH, LIWKUS
      INTEGER LIWORK, LIW1, LJCNH, LNXTRW, LRHS, LW, LWKH, LWKUSE
      INTEGER J, LPERT, MLP, LWMA27, LWORK, MAXFRT, NNZH, NSTEPS
      INTEGER NFIXED, IPCD( 1 ), JCNCD( 1 ), LBD, LIPBD, LIRNBD, LP2
      INTEGER LQ, LQMAX, LR, LRHS2, LRMAX, NFMAX, INFO39
C ** Correction 2. 16/07/94: 1 line corrected **
      INTEGER IDUM1, IDUM2, NRLBDU, NIRBDU, LIH, NSEMIW, LINXTR, MAXSBW
C ** Correction 2. 16/07/94:  end of correction **
      LOGICAL PRONEL, PRNTER, CONSIS, NEGATE
CS    REAL DOTPRD, OPS, ZERO, HALF, PERTUR, ALPHA, SGRAD
CS    REAL STHS, ONE, ONEPEP, TEN, HP
CS    REAL AMODEL, ATEMP, GNRMSQ, CD( 1 )
CD    DOUBLE PRECISION DOTPRD, OPS, ZERO, HALF, PERTUR, ALPHA, SGRAD
CD    DOUBLE PRECISION STHS, ONE, ONEPEP, TEN, HP
CD    DOUBLE PRECISION AMODEL, ATEMP, GNRMSQ, CD( 1 )
      INTRINSIC ABS, DBLE, FLOAT, INT, MAX, MIN, SQRT
CS    EXTERNAL SASMBL, SDSCDR, MA27A , MA27B , MA27C , SNEGCR, SSYSCH
CD    EXTERNAL DASMBL, DDSCDR, MA27AD, MA27BD, MA27CD, DNEGCR, DSYSCH
CS    EXTERNAL SMCFA, SASSLB, SASSLC
CD    EXTERNAL DMCFA, DASSLB, DASSLC
CS    COMMON / MA27E  / OPS, IDUM1( 7 ), NRLBDU, NIRBDU, IDUM2( 5 )
CD    COMMON / MA27ED / OPS, IDUM1( 7 ), NRLBDU, NIRBDU, IDUM2( 5 )
      INTEGER           ITERCG, ITCGMX, NGEVAL, ISKIP , IFIXED, NSEMIB
CS    REAL              ACCCG , RATIO , RADIUS, RADMAX, FINDMX, CMA31
CD    DOUBLE PRECISION  ACCCG , RATIO , RADIUS, RADMAX, FINDMX, CMA31
CS    COMMON / SCOMSB / ACCCG , RATIO , RADIUS, RADMAX, FINDMX, CMA31,
CD    COMMON / DCOMSB / ACCCG , RATIO , RADIUS, RADMAX, FINDMX, CMA31,
     *                  ITERCG, ITCGMX, NGEVAL, ISKIP , IFIXED, NSEMIB
CS    REAL             EPSMCH, EPSNEG, TINY, BIG
CD    DOUBLE PRECISION EPSMCH, EPSNEG, TINY, BIG
CS    COMMON / SMACHN / EPSMCH, EPSNEG, TINY, BIG
CD    COMMON / DMACHN / EPSMCH, EPSNEG, TINY, BIG
C ** Correction 1. 13/01/94: 2 lines corrected **
CS    SAVE            / SCOMSB /, / SMACHN /
CD    SAVE            / DCOMSB /, / DMACHN /
C ** Correction 1. 13/01/94:  end of correction **
CS    DATA ZERO, HALF, ONE, TEN / 0.0E+0, 5.0E-1, 1.0E+0, 1.0E+1 /
CD    DATA ZERO, HALF, ONE, TEN / 0.0D+0, 5.0D-1, 1.0D+0, 1.0D+1 /
      PRONEL = IPRINT .EQ. 2
      PRNTER = IPRINT .GE. 5
      IF ( PRONEL ) WRITE( MP, 2120 ) GSTOP
      ONEPEP = ONE + TEN * EPSMCH
C
C  DEFINE THE WORK SPACE FOR ASEMBL.
C
      LWKUSE = N + 3 * MAXSEL
      LWKH   = LWK - LWKUSE
      IF ( LWKH .LE. 0 ) THEN
         INFORM = 10
         WRITE( MP, 2000 )
         RETURN
      END IF
      LIWKUS = N
      LIWKH  = LIWK - LIWKUS
      LH = MIN( LWKH, ( LIWKH - 3 * NFREE ) / 4 )
      IF ( LH .LE. 0 ) THEN
         INFORM = 11
         WRITE( MP, 2010 )
         RETURN
      END IF
      LIRNH  = 1
      LJCNH  = LIRNH + LH
      LNXTRW = LJCNH + LH
C
C  ASSEMBLE THE HESSIAN RESTRICTED TO THE IFREE VARIABLES.
C
      LIH    = LH
      LINXTR = LH + NFREE     
      NSEMIW = NSEMIB 
C ** Correction 3. 16/07/94: 1 line corrected **
C ** Correction 4. 05/11/97: 1 line corrected **
C ** Correction 5. 19/08/98: added another final argument to call to ASMBL.
CS    CALL SASMBL( N , NG, MAXSEL, NSEMIW, LH    , LH    , NNZH  , 
CD    CALL DASMBL( N , NG, MAXSEL, NSEMIW, LH    , LH    , NNZH  , 
     *             NFREE , IFREE , ISTADH, LSTADH, ICNA  , LICNA , 
     *             ISTADA, LSTADA, INTVAR, LNTVAR, IELVAR, LELVAR,
     *             IELING, LELING, ISTADG, LSTADG, ISTAEV, LSTAEV,
     *             ISTAGV, LNSTGV, ISVGRP, LNVGRP, IWK  ( LIRNH ),
     *             IWK ( LJCNH ) , IWK ( LNXTRW ), LINXTR,
     *             IWK( LIWKH + 1 )      , LIWK - LIWKH  ,
     *             A , LA, GUVALS, LNGUVL, HUVALS, 
     *             LNHUVL, GVALS2, GVALS3, GSCALE, ESCALE, LESCAL, 
     *             WK    , WK (LWKH + 1 ), LWK - LWKH    , GXEQX , 
     *             LGXEQX, INTREP, LINTRE, RANGES, IPRINT, MP    ,
     *             .FALSE.       , MAXSBW, INFORM, .TRUE., .FALSE. )
C ** Correction 5. 19/08/98: end of correction **
C ** Correction 4. 05/11/97: end of correction **
C ** Correction 3. 16/07/94: end of correction **
      IF ( INFORM .GT. 0 ) THEN
         INFORM = 11
         WRITE( MP, 2010 )
         RETURN
      END IF
C
C  COMPRESS THE VECTOR IWK.
C
      DO 10 I = 1, NNZH
         IWK( NNZH + I ) = IWK( LH + I )
   10 CONTINUE
      LJCNH = LIRNH + NNZH
C
C  PRINT THE ASSEMBLED HESSIAN.
C
C      DO 25 I = 1, NNZH
C         J = NNZH + I
C         WRITE( MP ,2500 ) IWK( I ), IWK( J ), WK( I )
C25    CONTINUE
C
C  DEFINE THE WORK SPACE FOR MA27AD.
C
      LIWORK = LIWK - 2 * NNZH - 5 * NFREE
      IF ( LIWORK .LE. 0 ) THEN
         INFORM = 11
         WRITE( MP, 2010 )
         RETURN
      END IF
      LIKEEP = LJCNH + NNZH
      LIW1   = LIKEEP + 3 * NFREE
      LIMA27 = LIW1   + 2 * NFREE
C
C  CHOOSE PIVOTS FOR GAUSSIAN ELIMINATION.
C
      IFLAG = 0
CS    CALL MA27A ( NFREE, NNZH, IWK( LIRNH ), IWK( LJCNH ),
CD    CALL MA27AD( NFREE, NNZH, IWK( LIRNH ), IWK( LJCNH ),
     *             IWK( LIMA27 ), LIWORK, IWK( LIKEEP ),
     *             IWK( LIW1 ), NSTEPS, IFLAG )
      IF ( IFLAG .LT. 0 ) THEN
         WRITE( MP, 2020 ) IFLAG
         IF ( IFLAG .EQ. - 3 ) THEN
            INFORM = 11
            WRITE( MP, 2010 )
         ELSE
            INFORM = 10
            WRITE( MP, 2000 )
         END IF
         RETURN
      END IF
C
C  DEFINE THE WORK SPACE FOR MA27BD OR MCFAD, MA27CD, NEGCUR OR DSCDIR.
C
      LWORK = LWK - 3 * NFREE
      IF ( LWORK .LE. 0 ) THEN
         INFORM = 10
         WRITE( MP, 2000 )
         RETURN
      END IF
      LWMA27 = 1
      LRHS   = LWMA27 + LWORK
      LPERT  = LRHS   + NFREE
      LW     = LPERT  + NFREE
C
C  FACTORIZE THE RESTRICTED HESSIAN.
C
      PERTUR = ZERO
      IF ( MODCHL ) THEN
         MLP = 0
         IF ( IPRINT .GE. 20 ) MLP = MP
CS       CALL SMCFA( NFREE, NNZH, IWK( LIRNH ), IWK( LJCNH ),
CD       CALL DMCFA( NFREE, NNZH, IWK( LIRNH ), IWK( LJCNH ),
     *                WK( LWMA27 ), LWORK, IWK( LIMA27 ), LIWORK,
     *                IWK( LIKEEP ), NSTEPS, MAXFRT, IWK( LIW1 ), IFLAG,
     *                WK( LRHS ), WK( LPERT ), MLP )
         DO 20 J = 1, NFREE
            PERTUR = MAX( PERTUR, ABS( WK( LPERT + J - 1 ) ) )
   20    CONTINUE
         IF ( PRNTER ) THEN
            DO 21 J = 1, NFREE
               IF ( WK( LPERT + J - 1 ) .NE. ZERO ) WRITE( MP, 2110 )
     *                         WK( LPERT + J - 1 ), IFREE( J )
   21       CONTINUE
         END IF
      ELSE
CS       CALL MA27B ( NFREE, NNZH, IWK( LIRNH ), IWK( LJCNH ),
CD       CALL MA27BD( NFREE, NNZH, IWK( LIRNH ), IWK( LJCNH ),
     *                WK( LWMA27 ), LWORK, IWK( LIMA27 ), LIWORK,
     *                IWK( LIKEEP ), NSTEPS, MAXFRT, IWK( LIW1 ), IFLAG)
      END IF
C
C  RECORD THE RELATIVE FILL-IN.
C
      RATIO = DBLE( FLOAT( NRLBDU ) ) / DBLE( FLOAT( NNZH ) )
      IF ( IFLAG .NE. 0 .AND. IFLAG .NE. 3 ) THEN
         WRITE( MP, 2030 ) IFLAG
         IF ( IFLAG .EQ. - 3 ) THEN
            INFORM = 11
            WRITE( MP, 2010 )
         ELSE
            INFORM = 10
            WRITE( MP, 2000 )
         END IF
         RETURN
      END IF
C
C     CALL PRINT( IWK( LIMA27 ), WK )
C
      IF ( MODCHL ) THEN
         INFORM = 1
      ELSE
CS       CALL SSYSCH( NFREE, IWK( LIMA27 ), WK( LWMA27 ), EPSMCH,
CD       CALL DSYSCH( NFREE, IWK( LIMA27 ), WK( LWMA27 ), EPSMCH,
     *                IWK( LIKEEP ), INFORM )
      END IF
C
C  DEFINE THE MODEL OF THE GRADIENT RESTRICTED TO THE
C  FREE VARIABLES IN WK(LRHS).
C
      DO 30 J = 1, NFREE
         WK( LRHS + J - 1 ) = - GMODEL( IFREE( J ) )
   30 CONTINUE
      IF ( INFORM .EQ. 1 ) THEN
C
C  THE RESTRICTED HESSIAN IS POSITIVE DEFINITE -
C  SOLVE THE LINEAR SYSTEM H * S = - GMODEL. S IS OUTPUT IN GMODEL.
C
CS       CALL MA27C ( NFREE, WK( LWMA27 ), LWORK, IWK( LIMA27 ), LIWORK,
CD       CALL MA27CD( NFREE, WK( LWMA27 ), LWORK, IWK( LIMA27 ), LIWORK,
     *                WK( LW ), MAXFRT, WK( LRHS ), IWK( LIW1 ), NSTEPS)
         NEGATE = .FALSE.
      END IF
      IF ( INFORM .EQ. 2 ) THEN
C
C  THE RESTRICTED HESSIAN IS NOT POSITIVE DEFINITE
C  COMPUTE A DIRECTION OF NEGATIVE CURVATURE.
C
CS       CALL SNEGCR( NFREE, IWK( LIMA27 ), WK( 1 ), IWK( LIW1 ),
CD       CALL DNEGCR( NFREE, IWK( LIMA27 ), WK( 1 ), IWK( LIW1 ),
     *                NFREE, NUMBER, NEXT, WK( LRHS ), STHS )
C
C  CHECK IF THE DOT PRODUCT OF THIS DIRECTION TIMES GMODEL IS > 0.
C
         DOTPRD = ZERO
         DO 40 J = 1, NFREE
            DOTPRD = DOTPRD - GMODEL( IFREE( J ) ) * WK( LRHS + J - 1 )
   40    CONTINUE
         NEGATE = DOTPRD .LE. ZERO
      END IF
      IF ( INFORM .EQ. 3 ) THEN
CS       CALL SDSCDR( NFREE, IWK( LIMA27 ), WK( 1 ), WK( LRHS ),
CD       CALL DDSCDR( NFREE, IWK( LIMA27 ), WK( 1 ), WK( LRHS ),
     *                IWK( LIKEEP ), CONSIS, IWK( LIW1 ), NFREE,
     *                WK( LW ), EPSMCH )
         NEGATE = .FALSE.
         IF ( .NOT. CONSIS ) THEN
            INFORM = 4
            STHS   = ZERO
         END IF
      END IF
C
C  CONSTRUCT THE SOLUTION P FROM WK(LRHS).
C
      IF ( NEGATE ) THEN
         DO 50 J = 1, NFREE
            P( IFREE( J ) ) = - WK( LRHS + J - 1 )
   50    CONTINUE
      ELSE
         DO 60 J = 1, NFREE
            P( IFREE( J ) ) = WK( LRHS + J - 1 )
   60    CONTINUE
      END IF
C
C  CALCULATE THE SLOPE OF THE QUADRATIC MODEL FROM THE CAUCHY POINT.
C
      SGRAD    = ZERO
      DO 70 J  = 1, NFREE
         SGRAD = SGRAD + P( IFREE( J ) ) * GMODEL( IFREE( J ) )
   70 CONTINUE
C
C  OBTAIN THE CURVATURE.
C
      IF ( INFORM .EQ. 1 .OR. INFORM .EQ. 3 ) STHS = - SGRAD
C
C  MODIFY THE CURVATURE IF THE HESSIAN HAS BEEN PERTURBED.
C
      IF ( PERTUR .GT. ZERO ) THEN
         INFORM = 5
         DO 80 J = 1, NFREE
            STHS = STHS - WK( LPERT + J - 1 ) * WK( LRHS + J - 1 ) ** 2
   80    CONTINUE
      END IF
      NFIXED = 0
  100 CONTINUE
C
C  COMPUTE THE STEPLENGTH TO THE MINIMIZER OF THE MODEL ALONG THE
C  CURRENT SEARCH DIRECTION.
C
      IF ( INFORM .EQ. 1 .OR. INFORM .EQ. 3 .OR. INFORM .EQ. 5 ) THEN
         IF ( INFORM .EQ. 5 ) THEN
            IF ( STHS .GT. ZERO ) THEN
               ALPHA = - SGRAD / STHS
            ELSE
               ALPHA = BIG
            END IF
         ELSE
            ALPHA = ONE
         END IF
      ELSE
         ALPHA = BIG
      END IF
C
C  IF REQUIRED, PRINT DETAILS OF THE CURRENT STEP.
C
      IF ( IPRINT .GE. 20 ) THEN
         WRITE( MP, 2060 )
         DO 110 J = 1, NFREE
            I     = IFREE( J )
            IF ( I .GT. 0 ) WRITE( MP, 2070 ) I, XT( I ), GMODEL( I ), 
     *         P( I ), BND( I, 1 ), BND( I, 2 )
  110    CONTINUE
      END IF
C
C  FIND THE LARGEST FEASIBLE STEP IN THE DIRECTION P FROM XT.
C
      AMODEL = ALPHA
      DO 120 J = 1, NFREE
         I     = IFREE( J )
         IF ( I .GT. 0 ) THEN
            IF ( ABS( P( I ) ) .GT. EPSMCH ) THEN
               IF ( P( I ) .GT. ZERO ) THEN
                  ALPHA = MIN(ALPHA, ( BND( I, 2 ) - XT( I ) ) / P( I ))
               ELSE
                  ALPHA = MIN(ALPHA, ( BND( I, 1 ) - XT( I ) ) / P( I ))
               END IF
            END IF
         END IF
  120 CONTINUE
C
C  UPDATE THE MODEL FUNCTION VALUE.
C
      FMODEL = FMODEL + ALPHA * ( SGRAD + ALPHA * HALF * STHS )
C
C  IF REQUIRED, PRINT THE MODEL VALUE, SLOPE AND THE STEP TAKEN.
C
      IF ( PRNTER ) WRITE( MP, 2040 ) AMODEL, ALPHA, FMODEL, SGRAD, STHS
C
C  CHECK TO SEE IF THE BOUNDARY  IS ENCOUNTERED.
C
      GNRMSQ = ZERO
      IF ( ALPHA .LT. AMODEL ) THEN
C
C  A FREE VARIABLE HAS ENCOUNTERED A BOUND. MAKE A SECOND PASS TO
C  DETERMINE WHICH VARIABLE AND COMPUTE THE SIZE OF THE MODEL GRADIENT
C  AT THE NEW POINT.
C
         DO 130 J = 1, NFREE
            I     = IFREE( J )
            IF ( I .GT. 0 ) THEN
               IF ( ABS ( P( I ) ) .GT. EPSMCH ) THEN
                  ATEMP = BIG
                  IF ( P( I ) .GT. ZERO ) THEN
                     ATEMP = ( BND( I, 2 ) - XT( I ) ) / P( I )
                  ELSE
                     ATEMP = ( BND( I, 1 ) - XT( I ) ) / P( I )
                  END IF
C
C  VARIABLE I ENCOUNTERS ITS BOUND.
C
                  IF ( ATEMP .LE. ALPHA * ONEPEP ) THEN
                     IF ( PRNTER ) WRITE( MP, 2090 ) I
                     IFREE( J ) = - I
                  ELSE
C
C  UPDATE THE GRADIENT.
C
                     IF ( INFORM .EQ. 1 .OR. INFORM .GE. 4 )
     *                  HP = -GMODEL( I )
                     IF ( INFORM .EQ. 4 ) HP = ZERO
                     IF ( INFORM .EQ. 5 )
     *                  HP = HP - WK( LPERT + J - 1 )*P( I )
                     GMODEL( I ) = GMODEL( I ) + ALPHA * HP
                     GNRMSQ      = GNRMSQ + GMODEL( I ) ** 2
                  END IF
                  XT( I ) = XT( I ) + ALPHA * P( I )
               END IF
            END IF
  130    CONTINUE
      ELSE
C
C  STEP TO THE NEW POINT.
C
         DO 140 J = 1, NFREE
            I     = IFREE( J )
            IF ( I .GT. 0 ) THEN
               IF ( ABS( P( I ) ) .GT. EPSMCH ) THEN
                  XT( I ) = XT( I ) + ALPHA * P( I )
               END IF
            END IF
  140    CONTINUE
      END IF
C     DO 140 J = 1, NFREE
C        I     = IFREE( J )
C        IF ( I .GT. 0 ) THEN
C           IF ( ABS( P( I ) ) .GT. EPSMCH ) THEN
C              XT( I ) = XT( I ) + ALPHA * P( I )
C
C  UPDATE THE GRADIENT.
C
C              IF ( INFORM .EQ. 1 .OR. INFORM .GE. 4 ) HP = -GMODEL( I )
C              IF ( INFORM .EQ. 4 ) HP = ZERO
C              IF ( INFORM .EQ. 5 ) HP = HP - WK( LPERT + J - 1 )*P( I )
C              GMODEL( I ) = GMODEL( I ) + ALPHA * HP
C              GNRMSQ      = GNRMSQ + GMODEL( I ) ** 2
C           END IF
C        END IF
C 140 CONTINUE
C
C  IF THE MODEL GRADIENT IS SUFFICIENTLY SMALL, EXIT.
C
      IF ( PRNTER ) WRITE( MP, 2050 ) GNRMSQ, GSTOP
      IF ( PRONEL ) WRITE( MP, 2130 ) NFIXED, FMODEL, GNRMSQ
      IF ( GNRMSQ .LE. GSTOP ) RETURN
C
C  CONTINUE THE MINIMIZATION IN THE RESTRICTED SUBSPACE.
C
      IF ( .NOT. MODCHL ) RETURN
C
C  SET UP FURTHER WORKSPACE ADDRESSES TO PARTITION THE UNUSED PORTIONS
C  OF WK AND IWK. CALCULATE THE LARGEST NUMBER OF VARIABLES WHICH CAN
C  BE FIXED.
C
      IF ( NFIXED .EQ. 0 ) THEN
         NFMAX = MIN( NFREE, ( LIWORK - NIRBDU - 1 ) / 2 , INT( SQRT(
     *             FLOAT( 2 * ( LWORK - NRLBDU - 2 * NFREE ) ) ) ) - 5 )
         IF ( PRNTER ) WRITE( MP, 2080 ) NFMAX
         IF ( NFMAX .LE. 0 ) RETURN
C
C  SET INTEGER WORKSPACE ADDRESSES.
C
         LIRNBD = LIMA27 + NIRBDU
         LIPBD  = LIRNBD + NFMAX
         IWK( LIPBD ) = 1
C
C  SET REAL WORKSPACE ADDRESSES.
C
         LQMAX  = NFMAX
         LRMAX  = NFMAX * ( NFMAX + 1 ) / 2
         LBD    = LWMA27 + NRLBDU
         LQ     = LBD    + NFMAX
         LR     = LQ     + LQMAX
         LRHS2  = LR     + LRMAX
         LP2    = LRHS2  + NFREE + NFMAX
      END IF
C
C  DETERMINE WHICH VARIABLES ARE TO BE FIXED.
C
      DO 220 J = 1, NFREE
         I     = IFREE( J )
         IF ( I .LT. 0 ) THEN
C
C  IF MORE THAN NFMAX VARIABLES HAVE BEEN FIXED, RETURN.
C
            IF ( NFIXED .GE. NFMAX ) RETURN
C
C  UPDATE THE FACTORIZATION OF THE SCHUR COMPLEMENT TO ALLOW FOR
C  THE REMOVAL OF THE J-TH ROW AND COLUMN OF THE ORIGINAL HESSIAN -
C  THIS REMOVAL IS EFFECTED BY APPENDING THE J-TH ROW AND COLUMN
C  OF THE IDENTITY MATRIX TO THE HESSIAN.
C
            WK( LBD + NFIXED ) = ONE
            IWK( LIRNBD + NFIXED ) = J
            IWK( LIPBD + NFIXED + 1 ) = NFIXED + 2
            INFO39 = 1
  210       CONTINUE
CS          CALL SASSLC( NFREE, NFIXED, 4, NFMAX, WK( LBD ),
CD          CALL DASSLC( NFREE, NFIXED, 4, NFMAX, WK( LBD ),
     *                   IWK( LIRNBD ), IWK( LIPBD ), 1, CD, JCNCD,
     *                   IPCD, WK( LR ), LRMAX, WK( LQ ), LQMAX,
     *                   WK( LRHS2 ), WK( LRHS ), INFO39 )
            IF ( INFO39 .GT. 0 ) THEN
C
C  ASSL REQUIRES ADDITIONAL INFORMATION. COMPUTE THE SOLUTION TO
C  THE EQUATION H * S = RHS, RETURNING THE SOLUTION S IN RHS.
C
CS             CALL MA27C ( NFREE, WK( LWMA27 ), LWORK, IWK( LIMA27 ),
CD             CALL MA27CD( NFREE, WK( LWMA27 ), LWORK, IWK( LIMA27 ),
     *                      LIWORK, WK( LW ), MAXFRT, WK( LRHS ),
     *                      IWK( LIW1 ), NSTEPS )
               GO TO 210
            END IF
            IF ( INFO39 .LT. 0 ) THEN
               IF ( PRNTER ) WRITE( MP, 2100 ) INFO39
               RETURN
            END IF
            NFIXED = NFIXED + 1
         END IF
  220 CONTINUE
C
C  DEFINE THE MODEL OF THE GRADIENT RESTRICTED TO THE
C  FREE VARIABLES IN RHS.
C
      DO 230 J = 1, NFREE
         I     = IFREE( J )
         IF ( I .GT. 0 ) THEN
            WK( LRHS2 + J - 1 ) = - GMODEL( I )
         ELSE
            WK( LRHS2 + J - 1 ) = ZERO
         END IF
  230 CONTINUE
      DO 240 J = 1, NFIXED
         WK( LRHS2 + NFREE + J - 1 ) = ZERO
  240 CONTINUE
C
C  SOLVE THE NEW LINEAR SYSTEM H * P2 = RHS2.
C
      INFO39 = 1
  310 CONTINUE
CS    CALL SASSLB( NFREE, NFIXED, 4, NFMAX, WK( LBD ),
CD    CALL DASSLB( NFREE, NFIXED, 4, NFMAX, WK( LBD ),
     *             IWK( LIRNBD ), IWK( LIPBD ), 1, CD, JCNCD,
     *             IPCD, WK( LR ), LRMAX, WK( LQ ), LQMAX,
     *             WK( LRHS2 ), WK( LP2 ), WK( LRHS ), INFO39 )
      IF ( INFO39 .GT. 0 ) THEN
C
C  ASSL REQUIRES ADDITIONAL INFORMATION. COMPUTE THE SOLUTION TO
C  THE EQUATION H * S = RHS, RETURNING THE SOLUTION S IN RHS.
C
CS       CALL MA27C ( NFREE, WK( LWMA27 ), LWORK, IWK( LIMA27 ),
CD       CALL MA27CD( NFREE, WK( LWMA27 ), LWORK, IWK( LIMA27 ),
     *                LIWORK, WK( LW ), MAXFRT, WK( LRHS ),
     *                IWK( LIW1 ), NSTEPS )
         GO TO 310
      END IF
C
C  STORE THE CORRECTIONS TO THE FREE VARIABLES.
C
      DO 320 J = 1, NFREE
         I     = IFREE( J )
         IF ( I .LT. 0 ) IFREE( J ) = 0
         IF ( I .GT. 0 ) P( I )     = WK( LP2 + J - 1 )
  320 CONTINUE
C
C  CALCULATE THE SLOPE OF THE QUADRATIC MODEL FROM THE CAUCHY POINT.
C
      SGRAD    = ZERO
      DO 330 J = 1, NFREE
         I     = IFREE( J )
         IF ( I .GT. 0 ) SGRAD = SGRAD + P( I ) * GMODEL( I )
  330 CONTINUE
C
C  OBTAIN THE CURVATURE.
C
      IF ( INFORM .EQ. 1 .OR. INFORM .EQ. 3 .OR. INFORM .EQ. 5 )
     *   STHS = - SGRAD
C
C  MODIFY THE CURVATURE IF THE HESSIAN HAS BEEN PERTURBED.
C
      IF ( PERTUR .GT. ZERO ) THEN
         DO 340 J = 1, NFREE
            IF ( IFREE( J ) .GT. 0 ) STHS = STHS -
     *         WK( LPERT + J - 1 ) * WK( LP2 + J - 1 ) ** 2
  340    CONTINUE
      END IF
      GO TO 100
C
C NON-EXECUTABLE STATEMENTS.
C
 2000 FORMAT( ' ** MESSAGE FROM -FRNTL-', /,
     *        '    INCREASE THE PARAMETER -LWK-' )
 2010 FORMAT( ' ** MESSAGE FROM -FRNTL-', /,
     *        '    INCREASE THE PARAMETER -LIWK-' )
 2020 FORMAT( ' ** MESSAGE FROM -FRNTL-', /,
     *        '    VALUE OF IFLAG AFTER MA27A = ', I3 )
 2030 FORMAT( ' ** MESSAGE FROM -FRNTL-', /,
     *        '    VALUE OF IFLAG AFTER MA27B = ', I3 )
 2040 FORMAT( /, ' MODEL STEP AND ACTUAL STEP = ', 1P, 2D12.4, /,
     *           ' FMODEL, SGRAD, STHS = ', 3D12.4 )
 2050 FORMAT( /, ' MODEL GRADIENT ** 2 ', 1P, D12.4, ' GSTOP ', D12.4 )
 2060 FORMAT( /, '    I      XT          G           P           BL',
     *           '          BU' )
 2070 FORMAT( I5, 1P, 5D12.4 )
 2080 FORMAT( /, ' AT MOST ', I6, ' VARIABLES WILL BE FIXED ' )
 2090 FORMAT( ' VARIABLE ', I6, ' HAS ENCOUNTERED A BOUND IN FRONTL ' )
 2100 FORMAT( ' ** MESSAGE FROM -FRONTL-', /,
     *        '    VALUE OF IFLAG AFTER ASSLC = ', I3 )
 2110 FORMAT( ' PERTUBATION ', 1P, D12.4, ' FOR DIAGONAL ', I6 )
 2120 FORMAT( /, '    ** FRONTL ENTERED ** NFIXED     MODEL   ',
     *           '   GRADIENT GSTOP = ', 1P, D12.4 )
 2130 FORMAT( 25X, I7, 1P, 3D12.4 )
C2500 FORMAT( ' H(', I2, ',', I2, ') = ', 1P, D17.9 )
C
C  END OF FRNTL.
C
      END
