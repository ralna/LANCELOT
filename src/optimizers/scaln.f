C  THIS VERSION: Sat Dec 13 11:20:42 EST 1997

C ** Correction report.
C ** Corrections 1-11: Caused by replacement of MC19 by MC29
C **                   in the Harwell Subroutine Library **
C ** Correction 1. 10/08/93: 3 lines replaced **
C ** Correction 2. 10/08/93: 1 line corrected **
C ** Correction 3. 10/08/93: 4 lines replaced by 3 **
C ** Correction 4. 10/08/93: 1 line corrected **
C ** Correction 5. 10/08/93: 3 lines replaced by 2 **
C ** Correction 6. 10/08/93: 3 lines replaced **
C ** Correction 7. 10/08/93: 3 lines replaced **
C ** Correction 8. 10/08/93: 2 lines replaced by 3 **
C ** Correction 9. 10/08/93: 2 lines replaced by 1 **
C ** Correction 10. 10/08/93: 2 lines replaced by 1 **
C ** Correction 11. 10/08/93: 3 lines replaced **
C ** Correction 12. 19/08/96: 1 line corrected **
C ** Correction 13. 13/12/97: 1 line corrected **
C ** Correction 14. 13/12/97: 1 line corrected **
C ** Correction 15. 13/12/97: 1 line corrected **
C ** Correction 16. 13/12/97: 1 line corrected **
C ** End of Correction report.
CS    SUBROUTINE SSCALN( N , NG, NSLACK, NEL   , X , LX, NCALCF, ICALCF, 
CD    SUBROUTINE DSCALN( N , NG, NSLACK, NEL   , X , LX, NCALCF, ICALCF,
     *                   LCALCF, ICNA  , LICNA , ISTADA, LSTADA, IELING,
     *                   LELING, IRNGRJ, ICNGRJ, INTVAR, LNTVAR,
     *                   ISTADG, LSTADG, ISTAEV, LSTAEV, IELVAR, LELVAR,
     *                   ISTINV, LSTINV, ISVGRP, LSVGRP, ISTADH, LSTADH,
     *                   ISTAGV, LNSTAG, A , LA, B , LB, FT    , LFT   ,
     *                   GVALS2, LGVALS, FUVALS, LFUVAL, STOPGA, STOPCA, 
     *                   GSCALE, LGSCAL, ESCALE, LESCAL, VSCALE, LVSCAL, 
     *                   SCALEG, SCALEV, KNDOFC, LKNDOF, GRJAC , LGRJAC, 
C ** Correction 13. 13/12/97: 1 line corrected **
     *                   VARSCA, LVARSC, GRPSCA, LGRPSC, GETSCA, 
     *                   WK    , LWK   , WRK   , LWRK  , IWRK  , LIWRK ,
     *                   GXEQX , LGXEQX, INTREP, LINTRE, RANGES,
     *                   IOUT  , IPRINT, FDGRAD, INFORM )
      INTEGER            N , NG, NEL   , LWK   , LGRJAC, LSVGRP, LSTADH
      INTEGER            NSLACK, IOUT  , IPRINT, INFORM, NCALCF, LIWRK
      INTEGER            LICNA , LSTADA, LSTAEV, LNSTAG, LSTINV, LSTADG 
      INTEGER            LNTVAR, LCALCF, LKNDOF, LELVAR, LELING
      INTEGER            LGSCAL, LESCAL, LVSCAL, LGXEQX, LINTRE, LWRK
      INTEGER            LFT   , LA    , LB    , LGVALS, LFUVAL
      INTEGER            LVARSC, LGRPSC, LX    
CS    REAL               STOPGA, STOPCA   
CD    DOUBLE PRECISION   STOPGA, STOPCA   
C ** Correction 14. 13/12/97: 1 line corrected **
      LOGICAL            GETSCA, SCALEG, SCALEV, FDGRAD
      INTEGER            ICNA  ( LICNA  ), ISTADA( LSTADA ) 
      INTEGER            ISTAEV( LSTAEV ), ISTAGV( LNSTAG )
      INTEGER            INTVAR( LNTVAR ), ICALCF( LCALCF )
      INTEGER            KNDOFC( LKNDOF ), IELVAR( LELVAR )
      INTEGER            ISTINV( LSTINV ), ISTADG( LSTADG ) 
      INTEGER            IELING( LELING ), ISVGRP( LSVGRP ) 
      INTEGER            IRNGRJ( LGRJAC ), ICNGRJ( LGRJAC )
      INTEGER            IWRK  ( LIWRK  ), ISTADH( LSTADH )
C ** Correction 1. 10/08/93: 3 lines replaced **
CS    REAL               X     ( LX     ), GRJAC ( LGRJAC ), 
CD    DOUBLE PRECISION   X     ( LX     ), GRJAC ( LGRJAC ),
     *                   VARSCA( LVARSC ), GRPSCA( LGRPSC ),
C ** Correction 1. 10/08/93: end of correction **
     *                   FT    ( LFT    ), A     ( LA     ),    
     *                   B     ( LB     ), GVALS2( LGVALS ),
     *                   FUVALS( LFUVAL ), WK    ( LWK    ), 
     *                   GSCALE( LGSCAL ), ESCALE( LESCAL ),
     *                   VSCALE( LVSCAL ), WRK   ( LWRK   )
      LOGICAL            GXEQX ( LGXEQX ), INTREP( LINTRE )
      EXTERNAL           RANGES
C
C  TO CALCULATE SUITABLE VARIABLE AND CONSTRAINT SCALINGS
C  FOR A NONLINEAR PROGRAMMING PROBLEM INPUT TO LANCELOT.
C
C  NICK GOULD, 8TH JUNE 1990.
C  FOR CGT PRODUCTIONS.
C
C  LOCAL VARIABLES.
C
C ** Correction 2. 10/08/93: 1 line corrected **
      INTEGER            I , IG, NNZGRJ, LSTAGV, LENDGV, IESTEV
C ** Correction 2. 10/08/93: end of correction **
      INTEGER            NVAREL, IEL, L, IG1, J, II, K , IEINTV, NVAR
      INTEGER            IGETFD, L1    , L2    , L3    , L4    
      INTEGER            L5    , L6    , L7    , L8    , LXTEMP, L9     
      INTEGER            L10   , L11   , L12   , L13   , LWKSTR, LSPTRS 
      INTEGER            LSELTS, L14   , L15   , L16   , L17   , L18    
      INTEGER            L19   , L20   , L21   , L22   , L23   , L24    
      INTEGER            L25   , L26   , L27   , LSEND , LNPTRS, LNELTS 
      INTEGER            L28   , L29   , L30   , L31   , L32   , L33   
      INTEGER            L34   , L35   , L36   , L37   , L38   , L39    
      INTEGER            L40   , L41   , L42   , L43   , LNXTEM, L44    
      INTEGER            L45   , L46   , L47   , L48   , L49   , L50    
      INTEGER            L51   , L52   , L53   , L54   , L55   , L56    
      INTEGER            NTYPE , NSETS , L57   , LSTYPE, LSSWTR, LSSIWT 
      INTEGER            LSIWTR, LSWTRA, LNTYPE, LNSWTR, LNSIWT, LNIWTR 
      INTEGER            LNWTRA, LSISET, LSSVSE, LNISET, LNSVSE
C ** Correction 3. 10/08/93: 4 lines replaced by 3 **
CS    REAL               SMACHR, ZERO  , GI    , SCALE , SCALEE, FTT   , 
CD    DOUBLE PRECISION   DMACHR, ZERO  , GI    , SCALE , SCALEE, FTT   , 
     *                   RTEN  , GRPMAX, VARMIN, ONE   , TOL
C ** Correction 3. 10/08/93: end of correction **
      LOGICAL            ALTRIV, CENTRL
      INTRINSIC          ABS   , ANINT , LOG10 , MIN   , MAX   , EXP   
      INTRINSIC          SQRT
C ** Correction 4. 10/08/93: 1 line corrected **
      SAVE               ALTRIV, NVAR
C ** Correction 4. 10/08/93: end of correction **
CS    EXTERNAL           SSETVI, SSETVL, MC29A , SCOPY , SMACHR, SINITW
CD    EXTERNAL           DSETVI, DSETVL, MC29AD, DCOPY , SMACHR, DINITW
CS    EXTERNAL           SFDGRD
CD    EXTERNAL           DFDGRD, DMACHR
C
C  SET CONSTANT REAL PARAMETERS.
C
C ** Correction 5. 10/08/93: 3 lines replaced by 2 **
CS    PARAMETER        ( ZERO = 0.0E+0, ONE = 1.0E+0, RTEN  = 1.0E+1 )
CD    PARAMETER        ( ZERO = 0.0D+0, ONE = 1.0D+0, RTEN  = 1.0D+1 )
C ** Correction 5. 10/08/93: end of correction **
      IF ( INFORM .LT. 0 ) GO TO ( 110, 160 ), - INFORM
C
C  0. INITIALIZE DATA.
C
      NVAR = N - NSLACK
      IF ( LSTINV .LT. MAX( NVAR, NEL ) ) THEN
         WRITE( IOUT, 2000 ) MAX( NVAR, NEL ) - LSTINV
         INFORM = 4
         RETURN
      END IF
C ** Correction 6. 10/08/93: 3 lines replaced **
      IF ( NVAR .LE. 0 .OR. NG .LE. 0 ) RETURN
      IF ( LWK .LT. 2 * NG + 3 * NVAR ) THEN
         WRITE( IOUT, 2020 ) 2 * NG + 3 * NVAR - LWK
C ** Correction 6. 10/08/93: end of correction **
         INFORM = 5
         RETURN
      END IF
C
C  STORE THE INDICES OF VARIABLES WHICH APPEARS IN EACH GROUP IN ISVGRP.
C  THOSE FROM GROUP IG OCCUR FROM ISTAGV(IG) TO ISTAGV(IG+1) - 1.
C
      ALTRIV      = .TRUE.
      LSTAGV      = 1
      ISTAGV( 1 ) = LSTAGV
CDIR$ IVDEP
      DO 20 J        = 1, NVAR
         ISTINV( J ) = 0
   20 CONTINUE
      DO 70 IG  = 1, NG
C
C  CHECK TO SEE IF ALL OF THE GROUPS ARE TRIVIAL.
C
         IF ( .NOT. GXEQX( IG ) ) ALTRIV = .FALSE.
C
C  RUN THROUGH ALL THE ELEMENTAL VARIABLES CHANGING THE I-TH ENTRY OF
C  ISTINV( LSWKSP ) FROM ZERO TO ONE IF VARIABLE I APPEARS IN AN ELEMENT.
C  CONSIDER VARIABLES WHICH ARISE FROM NONLINEAR ELEMENTS.
C
         DO 40 K    = ISTADG( IG ), ISTADG( IG + 1 ) - 1 
            IEL     = IELING( K )
            DO 30 J = ISTAEV( IEL ), ISTAEV( IEL + 1 ) - 1
               I    = IELVAR( J )
               IF ( ISTINV( I ) .EQ. 0 ) THEN
                  IF ( LSTAGV .GT. LSVGRP ) THEN
                     WRITE( IOUT, 2010 ) LSVGRP - LSTAGV
                     INFORM = 4
                     RETURN
                  END IF
                  ISTINV( I ) = 1
C
C  RECORD THE NONLINEAR VARIABLES FROM THE IG-TH GROUP.
C
                  ISVGRP( LSTAGV ) = I
                  LSTAGV           = LSTAGV + 1
               END IF
   30       CONTINUE
   40    CONTINUE
C
C  CONSIDER VARIABLES WHICH ARISE FROM THE LINEAR ELEMENT.
C
         DO 50 J = ISTADA( IG ), ISTADA( IG + 1 ) - 1
            I    = ICNA( J )
            IF ( I .LE. NVAR ) THEN
               IF ( ISTINV( I ) .EQ. 0 ) THEN
                  IF ( LSTAGV .GT. LSVGRP ) THEN
                     WRITE( IOUT, 2010 ) LSVGRP - LSTAGV
                     INFORM = 4
                     RETURN
                  END IF
                  ISTINV( I ) = 1
C
C  RECORD THE LINEAR VARIABLES FROM THE IG-TH GROUP.
C
                  ISVGRP( LSTAGV ) = I
                  LSTAGV           = LSTAGV + 1
               END IF
            END IF
   50    CONTINUE
C
C  RESET ISTINV TO ZERO.
C
         DO 60 J                  = ISTAGV( IG ), LSTAGV - 1
            ISTINV( ISVGRP( J ) ) = 0
   60    CONTINUE
         ISTAGV( IG + 1 ) = LSTAGV
   70 CONTINUE
C
C  SET UP THE STARTING ADDRESSES FOR THE ELEMENT GRADIENTS
C  WITH RESPECT TO THEIR INTERNAL VARIABLES.
C
      ISTINV( 1 ) = NEL + 1
      IF ( NEL .GT. 0 ) THEN
         DO 80 I            = 1, NEL
            ISTINV( I + 1 ) = ISTINV( I ) + INTVAR( I )
            ICALCF( I ) = I
   80    CONTINUE
      END IF
C
C  IF FINITE-DIFFERENCE GRADIENTS ARE REQUIRED, SET UP THE NECESSARY
C  DATA STRUCTURES (THIS IS A BIT WASTEFUL AS NOT ALL THE PARTITIONS 
C  ARE USED).
C
      IF ( FDGRAD ) THEN
         DO 90 I         = 1, NEL
            ISTINV( I  ) = INTVAR( I )
   90    CONTINUE
CS       CALL SINITW( N , NG, NEL   , IELING, LELING, ISTADG, LSTADG,
CD       CALL DINITW( N , NG, NEL   , IELING, LELING, ISTADG, LSTADG,
     *                IELVAR, LELVAR, ISTAEV, LSTAEV, ISTINV, LNTVAR,
     *                ISTADH, LSTADH, ICNA  , LICNA , ISTADA, LSTADA,
     *                GXEQX , LGXEQX, INTREP, LINTRE, LFUVAL, ALTRIV, 
     *               .FALSE., FDGRAD, L1    , L2    , L3    , L4    , 
     *                L5    , L6    , L7    , L8    , LXTEMP, L9    , 
     *                L10   , L11   , L12   , L13   , LWKSTR, LSPTRS, 
     *                LSELTS, L14   , L15   , L16   , L17   , L18   , 
     *                L19   , L20   , L21   , L22   , L23   , L24   , 
     *                L25   , L26   , L27   , LSEND , LNPTRS, LNELTS, 
     *                L28   , L29   , L30   , L31   , L32   , L33   ,
     *                L34   , L35   , L36   , L37   , L38   , L39   , 
     *                L40   , L41   , L42   , L43   , LNXTEM, L44   , 
     *                L45   , L46   , L47   , L48   , L49   , L50   , 
     *                L51   , L52   , L53   , L54   , L55   , L56   , 
     *                NTYPE , NSETS , L57   , LSTYPE, LSSWTR, LSSIWT, 
     *                LSIWTR, LSWTRA, LNTYPE, LNSWTR, LNSIWT, LNIWTR, 
     *                LNWTRA, LSISET, LSSVSE, LNISET, LNSVSE, RANGES,
     *                IWRK  , LIWRK , WRK   , LWRK  ,
     *                IPRINT, IOUT  , INFORM )
         IF ( INFORM .GT. 0 ) THEN
            WRITE( IOUT, 2100 ) 
            RETURN
         END IF
         IGETFD = 0
      END IF
C
C  1. CALCULATE THE REQUIRED ELEMENT AND GROUP DERIVATIVES.
C
C  RETURN TO THE CALLING PROGRAM TO OBTAIN THE ELEMENT FUNCTION
C  AND DERIVATIVE VALUES AT THE INITIAL POINT.
C
      NCALCF = NEL
      INFORM = - 1
      RETURN
  110 CONTINUE
C
C  IF FINITE-DIFFERENCE GRADIENTS ARE USED, COMPUTE THEIR VALUES.
C
      IF ( FDGRAD .AND. NEL .GT. 0 ) THEN
C
C  STORE THE VALUES OF THE NONLINEAR ELEMENTS FOR FUTURE USE.
C
         IF ( IGETFD .EQ. 0 ) THEN
CS          CALL SCOPY( NEL, FUVALS, 1, WRK( LWKSTR + 1 ), 1 )
CD          CALL DCOPY( NEL, FUVALS, 1, WRK( LWKSTR + 1 ), 1 )
CS          CALL SCOPY( N, X, 1, WRK( LXTEMP + 1 ), 1 )
CD          CALL DCOPY( N, X, 1, WRK( LXTEMP + 1 ), 1 )
            CENTRL = .TRUE.
         END IF 
C
C  OBTAIN A FURTHER SET OF DIFFERENCES.
C
CS       CALL SFDGRD( N, NEL, IELVAR, LELVAR, ISTAEV, LSTAEV,
CD       CALL DFDGRD( N, NEL, IELVAR, LELVAR, ISTAEV, LSTAEV,
     *                IWRK  ( LSELTS + 1   ), LNELTS,
     *                IWRK  ( LSPTRS + 1   ), LNPTRS,
     *                IELING, LELING,
     *                IWRK  ( LSSVSE + 1   ), LNSVSE,
     *                IWRK  ( LSISET + 1   ), LNISET, NSETS , 
     *                IWRK  ( LSEND  + 1   ), NEL   , 
     *                ICALCF, LCALCF, NCALCF, INTVAR, LNTVAR,
     *                IWRK  ( LSSWTR + 1   ), LNSWTR,
     *                IWRK  ( LSSIWT + 1   ), LNSIWT, 
     *                IWRK  ( LSTYPE + 1   ), LNTYPE, NTYPE ,
     *                IWRK  ( LSIWTR + 1   ), LNIWTR,
     *                WRK   ( LSWTRA + 1   ), LNWTRA,
     *                WRK   ( LXTEMP + 1   ), X, 
     *                FUVALS, LFUVAL, CENTRL, IGETFD )
         IF ( IGETFD .GT. 0 ) RETURN
C
C  RESTORE THE VALUES OF THE NONLINEAR ELEMENTS AT X.
C
         IGETFD = NSETS + 1
CS       CALL SCOPY( NEL, WRK( LWKSTR + 1 ), 1, FUVALS, 1 )
CD       CALL DCOPY( NEL, WRK( LWKSTR + 1 ), 1, FUVALS, 1 )
      END IF
C
C  COMPUTE THE GROUP ARGUMENT VALUES FT.
C
      DO 140 IG = 1, NG
         FTT   = - B( IG )
C
C  INCLUDE THE CONTRIBUTION FROM THE LINEAR ELEMENT.
C
         DO 120 J = ISTADA( IG ), ISTADA( IG + 1 ) - 1
            IF ( ICNA( J ) .LE. NVAR ) FTT  = FTT + 
     *                                        A( J ) * X( ICNA( J ) )
  120    CONTINUE
C
C  INCLUDE THE CONTRIBUTIONS FROM THE NONLINEAR ELEMENTS.
C
         DO 130 J = ISTADG( IG ), ISTADG( IG + 1 ) - 1
            FTT  = FTT + ESCALE( J ) * FUVALS( IELING( J ) )
  130    CONTINUE
         FT( IG )     = FTT
         ICALCF( IG ) = IG
  140 CONTINUE
C
C  COMPUTE THE GROUP FUNCTION VALUES.
C
      IF ( ALTRIV ) THEN
         DO 150 IG       = 1, NG
            GVALS2( IG ) = FT( IG )
  150    CONTINUE
      ELSE
C
C  IF NECESSARY, RETURN TO THE CALLING PROGRAM TO OBTAIN THE GROUP
C  FUNCTION AND DERIVATIVE VALUES AT THE INITIAL POINT.
C
         NCALCF = NG
         INFORM = - 2
         RETURN
      END IF
  160 CONTINUE
      IF ( .NOT. ALTRIV ) THEN
CDIR$ IVDEP
         DO 170 IG = 1, NG
            IF ( GXEQX( IG ) ) GVALS2( IG ) = FT( IG )
  170    CONTINUE   
      END IF
C
C  2. FORM THE JACOBIAN MATRIX OF THE NONLINEAR FUNCTION
C     ( GROUP(1), .... , GROUP(NG) )(TRANSPOSE).
C
      NNZGRJ = 0
C
C  CONSIDER THE IG-TH GROUP.
C
      DO 270 IG = 1, NG
         IG1    = IG + 1
         LSTAGV = ISTAGV( IG )
         LENDGV = ISTAGV( IG1 ) - 1
C
C  INITIALIZE THE GROUP DERIVATIVE TO ZERO.
C
CS       CALL SSETVI( LENDGV - LSTAGV + 1, WK( 1 ), ISVGRP( LSTAGV ),
CD       CALL DSETVI( LENDGV - LSTAGV + 1, WK( 1 ), ISVGRP( LSTAGV ),
     *                ZERO )
C
C  SEE IF THE GROUP HAS ANY NONLINEAR ELEMENTS.
C
         DO 230 II = ISTADG( IG ), ISTADG( IG1 ) - 1
            IEL    = IELING( II )
            IEINTV = ISTINV( IEL )
            IESTEV = ISTAEV( IEL )
            NVAREL = ISTAEV( IEL + 1 ) - IESTEV
            SCALEE = ESCALE( II )
            IF ( INTREP( IEL ) ) THEN
C
C  THE IEL-TH ELEMENT HAS AN INTERNAL REPRESENTATION.
C
               CALL RANGES( IEL, .TRUE., FUVALS( IEINTV ),
     *                      WK( NVAR + 1 ), NVAREL,
     *                      ISTINV( IEL + 1 ) - IEINTV )
CDIR$ IVDEP
               DO 210 I   = 1, NVAREL
                  J       = IELVAR( IESTEV )
                  WK( J ) = WK( J ) + SCALEE * WK( NVAR + I )
                  IESTEV  = IESTEV + 1
  210          CONTINUE
            ELSE
C
C  THE IEL-TH ELEMENT HAS NO INTERNAL REPRESENTATION.
C
CDIR$ IVDEP
               DO 220 I   = 1, NVAREL
                  J       = IELVAR( IESTEV )
                  WK( J ) = WK( J ) + SCALEE * FUVALS( IEINTV )
                  IEINTV  = IEINTV + 1
                  IESTEV  = IESTEV + 1
  220          CONTINUE
            END IF
  230    CONTINUE
C
C  INCLUDE THE CONTRIBUTION FROM THE LINEAR ELEMENT.
C
CDIR$ IVDEP
         DO 240 K = ISTADA( IG ), ISTADA( IG1 ) - 1
            IF ( ICNA( K ) .LE. NVAR ) THEN
               J       = ICNA( K )
               WK( J ) = WK( J ) + A( K )
            END IF
  240   CONTINUE
C
C  STORE THE VALUES OF THE NONZERO ENTRIES OF THE GRADIENT OF
C  THE IG-TH GROUP IN GRJAC, ALONG WITH THEIR ROW AND COLUMN
C  INDICES IN THE OVERALL JACOBIAN.
C
         IF ( .NOT. GXEQX( IG ) ) THEN
            GI = GSCALE( IG ) * GVALS2( IG )
         ELSE
            GI = GSCALE( IG )
         END IF
CDIR$ IVDEP
         DO 250 J   = LSTAGV, LENDGV
            I       = ISVGRP( J )
            NNZGRJ  = NNZGRJ + 1
            IF ( NNZGRJ .GT. LGRJAC ) THEN
               WRITE( IOUT, 2030 ) NNZGRJ - LGRJAC
               INFORM = 5
               RETURN
            END IF
            IRNGRJ( NNZGRJ ) = IG
            ICNGRJ( NNZGRJ ) = I   
            GRJAC ( NNZGRJ ) = GI * WK( I )
            IF ( IPRINT .GE. 100 ) WRITE( IOUT, 2070 )
     *           IG, I, GRJAC( NNZGRJ )
  250    CONTINUE
  270 CONTINUE
C
C  3. REMOVE ALL ENTRIES WHICH ARE SMALLER THAN TOL TIMES THE
C     LARGEST ENTRIES IN BOTH THEIR ROWS AND COLUMNS.
C
C  FIND THE LARGEST ENTRIES, IN ABSOLUTE VALUE, IN EACH ROW AND COLUMN.
C
CS    TOL =       SMACHR ( 1 )
CD    TOL = SQRT( DMACHR( 1 ) )
CS    CALL SSETVL( NG,   WK( 1 ), 1, ZERO )
CD    CALL DSETVL( NG,   WK( 1 ), 1, ZERO )
C ** Correction 12. 19/08/96: 1 line corrected **
CS    CALL SSETVL( NVAR, WK( NG + 1 ), 1, ZERO )
CD    CALL DSETVL( NVAR, WK( NG + 1 ), 1, ZERO )
C ** Correction 12. 19/08/96: end of correction **
      DO 310 K = 1, NNZGRJ
         I     = IRNGRJ( K )
         J     = ICNGRJ( K )
         GI    = ABS( GRJAC( K ) )
         WK( I ) = MAX( WK( I ), GI )
         WK( NG + J ) = MAX( WK( NG + J ), GI )
  310 CONTINUE
C
C  REMOVE SMALL ENTRIES.
C
      L        = 0
      DO 320 K = 1, NNZGRJ
         I     = IRNGRJ( K )
         J     = ICNGRJ( K )
         GI    = ABS( GRJAC( K ) ) / TOL
         IF ( GI .GT. WK( I ) .OR. 
     *        GI .GT. WK( NG + J ) ) THEN
            L           = L + 1
            IRNGRJ( L ) = I
            ICNGRJ( L ) = J
            GRJAC( L )  = GRJAC( K )
         ELSE
            IF ( IPRINT .GE. 100 ) WRITE( IOUT, 2080 )
     *           I, J, GRJAC( K )
         END IF
  320 CONTINUE
      NNZGRJ = L 
C
C  4. APPLY THE ROW AND COLUMN EQUILIBRATING SCHEME OF CURTIS
C     AND REID TO THE JACOBIAN.
C
C ** Correction 7. 10/08/93: 3 lines replaced **
CS    CALL MC29A ( NG, NVAR, NNZGRJ, GRJAC, IRNGRJ, ICNGRJ, GRPSCA,
CD    CALL MC29AD( NG, NVAR, NNZGRJ, GRJAC, IRNGRJ, ICNGRJ, GRPSCA,
     *             VARSCA, WK, 0, L1 )
C ** Correction 7. 10/08/93: end of correction **
C
C  5. CALCULATE THE SCALE FACTORS.
C
C  OBTAIN THE SMALLEST VARIABLE AND GROUP SCALE FACTOR.
C
C ** Correction 8. 10/08/93: 2 lines replaced by 3 **
      GRPMAX = ZERO
CS    VARMIN = SMACHR( 5 )
CD    VARMIN = DMACHR( 5 )
C ** Correction 8. 10/08/93: end of correction **
CDIR$ IVDEP
      DO 510 I       = 1, NVAR
         VARSCA( I ) = EXP( VARSCA( I ) )
         VARMIN      = MIN( VARSCA( I ), VARMIN )
  510 CONTINUE   
CDIR$ IVDEP
      DO 520 IG = 1, NG
         IF ( KNDOFC( IG ) .NE. 1 ) THEN
            GRPSCA( IG ) = EXP( GRPSCA( IG ) )
            GRPMAX       = MAX( GRPSCA( IG ), GRPMAX )
         END IF
  520 CONTINUE   
C
C  SCALE THE FACTORS RELATIVE TO THEIR LARGEST MEMBERS.
C
CDIR$ IVDEP
      DO 530 I       = 1, NVAR
         VARSCA( I ) = RTEN ** ANINT( LOG10( VARSCA( I ) / VARMIN ) )
  530 CONTINUE   
CDIR$ IVDEP
      DO 540 IG = 1, NG
         IF ( KNDOFC( IG ) .NE. 1 ) GRPSCA( IG ) =
     *        RTEN ** ANINT( LOG10( GRPSCA( IG ) / GRPMAX ) )
  540 CONTINUE   
C
C  USE THE VARIABLE SCALINGS WITHIN LANCELOT.
C
      IF ( SCALEV ) THEN
         SCALE = ZERO
CDIR$ IVDEP
         DO 550 I       = 1, NVAR
C ** Correction 9. 10/08/93: 2 lines replaced by 1 **
            VSCALE( I ) =  VARSCA( I )  
C ** Correction 9. 10/08/93: end of correction **
            SCALE = MAX( SCALE, VSCALE( I ) )
  550    CONTINUE   
         STOPGA = STOPGA / SCALE
      END IF
C
C  USE THE GROUP SCALINGS WITHIN LANCELOT.
C
      IF ( SCALEG ) THEN
         SCALE = ONE
CDIR$ IVDEP
         DO 560 IG = 1, NG
            IF ( KNDOFC( IG ) .NE. 1 ) THEN
C ** Correction 10. 10/08/93: 2 lines replaced by 1 **
               GSCALE( IG ) = GSCALE( IG ) * GRPSCA( IG ) 
C ** Correction 10. 10/08/93: end of correction **
               SCALE = MIN( SCALE, GSCALE( IG ) )
            END IF
  560    CONTINUE   
         STOPCA = STOPCA * SCALE
      END IF
C
C  SUCCESSFUL CONCLUSION TO THE CALCULATION.
C
      IF ( IPRINT .GT. 0 ) THEN
C ** Correction 15. 13/12/97: 1 line corrected **
        IF ( SCALEG .OR. GETSCA ) THEN
            WRITE( IOUT, 2040 )
            J = 1
            DO 580 I = 1, NG
               IF ( KNDOFC( I ) .EQ. 1 ) THEN
                  IF ( I - 1 .GE. J ) WRITE( IOUT, 2060 )
     *                   ( IG, GRPSCA( IG ), IG = J, I - 1 )
                  J = I + 1
               END IF
 580        CONTINUE   
            IF ( NG .GE. J ) WRITE( IOUT, 2060 )
     *                   ( IG, GRPSCA( IG ), IG = J, NG )
         END IF
C ** Correction 16. 13/12/97: 1 line corrected **
         IF ( SCALEV .OR. GETSCA  ) THEN
            WRITE( IOUT, 2050 ) ( I, VARSCA( I ), I = 1, NVAR )
         END IF
         WRITE( IOUT, 2090 ) STOPGA, STOPCA
      END IF
      INFORM = 0
      RETURN
C
C  NON-EXECUTABLE STATEMENTS.
C
 2000 FORMAT( ' ** Return from SCALN. Size of ISTINV too small. ',
     *        '  Increase by at least ', I8 )
 2010 FORMAT( ' ** Return from SCALN. Size of ISVGRP too small.', 
     *        ' increase by at least ', I8 )
 2020 FORMAT( ' ** Return from SCALN. Size of WK too small. Increase',
     *        ' by at least ', I8 )
 2030 FORMAT( ' ** Return from SCALN. Size of GRJAC too small.', 
     *        ' Increase by at least ', I8 )
 2040 FORMAT( /, ' Multiply groups by the factors:' )
C ** Correction 11. 10/08/93: 3 lines replaced **
 2050 FORMAT( /, ' Divide variables by the factors:', 
     *        / 4( I6, 1P, D12.4 ) )
 2060 FORMAT( 4( I6, 1P, D12.4 ) )
C ** Correction 11. 10/08/93: end of correction **
 2070 FORMAT( ' Row ', I6, ' column ', I6, ' value ', 1P, D12.4 ) 
 2080 FORMAT( ' Row ', I6, ' column ', I6, ' value ', 1P, D12.4, 
     *        ' removed ' ) 
 2090 FORMAT( /, ' Scaled projected gradient tolerance = ', 1P, D12.4,
     *        /, ' Scaled constraint tolerance         = ', 1P, D12.4 )
 2100 FORMAT( ' ** Return from SCALN. Insufficient workspace ' )
C
C  END OF SCALN.
C
      END 
