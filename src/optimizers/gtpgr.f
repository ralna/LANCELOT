C  THIS VERSION: 14/08/1991 AT 11:59:37 AM.
C
C  ** FOR THE CRAY 2, LINES STARTING CDIR$ IVDEP TELL THE COMPILER TO
C     IGNORE POTENTIAL VECTOR DEPENDENCIES AS THEY ARE KNOWN TO BE O.K.
C
CS    SUBROUTINE SGTPGR( N, NG, NGEL, NN, NVAR, MAXSEL, QGNORM, SMALLH,
CD    SUBROUTINE DGTPGR( N, NG, NGEL, NN, NVAR, MAXSEL, QGNORM, SMALLH,
     *                   PJGNRM, CALCDI, DPRCND, MYPREC, IVAR, ISTADH,
     *                   LSTADH, ISTAEV, LSTAEV, IELVAR, LELVAR,
     *                   INTVAR, LNTVAR, IELING, LELING, ISLGRP, LNLGRP,
     *                   IGCOLJ, LNGCLJ, ISTAGV, LNSTGV, ISVGRP, LNVGRP,
     *                   IVALJR, LNVLJR, ISYMMD, ISYMMH, MAXSZH,
     *                   DGRAD, Q, GVALS2, GVALS3, DIAG, GSCALE,
     *                   ESCALE, LESCAL, GRJAC, LNGRJC, HUVALS, LNHUVL,
     *                   WK, LNWK, WKB, LNWKB, WKC, LNWKC,
     *                   GXEQX, LGXEQX, INTREP, LINTRE, RANGES )
C
C  FIND THE NORM OF THE PROJECTED GRADIENT, SCALED IF DESIRED.
C  IF REQUIRED, ALSO FIND THE DIAGONAL ELEMENTS OF THE HESSIAN MATRIX.
C
C  NICK GOULD, 15TH MARCH 1990.
C  FOR CGT PRODUCTIONS.
C
      INTEGER          N, NG, NN, NVAR, MAXSEL, MAXSZH
      INTEGER          LSTADH, LSTAEV, LELVAR, LNTVAR, LELING, LNLGRP
      INTEGER          LNGCLJ, LNSTGV, LNVGRP, LNVLJR, LESCAL, LNGRJC
      INTEGER          LNHUVL, LNWK, LNWKB, LNWKC, LGXEQX, LINTRE
      LOGICAL          CALCDI, DPRCND, MYPREC
CS    REAL             QGNORM, SMALLH, PJGNRM
CD    DOUBLE PRECISION QGNORM, SMALLH, PJGNRM
      INTEGER          IVAR( N ), ISTADH( LSTADH ), ISTAEV( LSTAEV )
      INTEGER          IELVAR( LELVAR ), INTVAR( LNTVAR )
      INTEGER          IELING( LELING ), ISLGRP( LNLGRP )
      INTEGER          IGCOLJ( LNGCLJ ), ISTAGV( LNSTGV )
      INTEGER          ISVGRP( LNVGRP ), IVALJR( LNVLJR )
      INTEGER          ISYMMH( MAXSZH, MAXSZH ), ISYMMD( MAXSZH )
CS    REAL             DGRAD( N ), Q( N ), GVALS2( NG ), GVALS3( NG ),
CD    DOUBLE PRECISION DGRAD( N ), Q( N ), GVALS2( NG ), GVALS3( NG ),
     *                 DIAG( N ), GSCALE( NG ), ESCALE( LESCAL ), 
     *                 GRJAC( LNGRJC ), HUVALS( LNHUVL ), 
     *                 WK( LNWK ), WKB( LNWKB ), WKC( LNWKC )
      LOGICAL          GXEQX( LGXEQX ), INTREP( LINTRE )
      EXTERNAL         RANGES
C
C  LOCAL VARIABLES.
C
      INTEGER          I, IEL, IG, IROW, J, IJHESS, K, KK, L, LL
      INTEGER          IELL, NIN, NVAREL, JCOL, IELHST, NGEL
CS    REAL             GDASH, G2DASH, PI, ZERO, ONE
CD    DOUBLE PRECISION GDASH, G2DASH, PI, ZERO, ONE
C
C  MACHINE FUNCTIONS.
C
      INTRINSIC        ABS, SQRT, MAX
CS    REAL             SDOT
CD    DOUBLE PRECISION DDOT
CS    EXTERNAL         SDOT, SSETVL
CD    EXTERNAL         DDOT, DSETVL
C
C  SET CONSTANT REAL PARAMETERS.
C
CS    PARAMETER ( ZERO   = 0.0E+0, ONE    = 1.0E+0 )
CD    PARAMETER ( ZERO   = 0.0D+0, ONE    = 1.0D+0 )
      QGNORM = ZERO
      IF ( MYPREC ) THEN
         DO 10 J   = 1, NVAR
            QGNORM = QGNORM + DGRAD( J ) * Q( IVAR( J ) )
   10    CONTINUE
         QGNORM = SQRT( QGNORM )
      ELSE
         IF ( CALCDI ) THEN
C
C  OBTAIN THE DIAGONAL ELEMENTS OF THE HESSIAN.
C  INITIALIZE THE DIAGONALS AS ZERO.
C
CS          CALL SSETVL( N, DIAG, 1, ZERO )
CD          CALL DSETVL( N, DIAG, 1, ZERO )
CS          CALL SSETVL( MAXSEL, WK, 1, ZERO )
CD          CALL DSETVL( MAXSEL, WK, 1, ZERO )
C
C  OBTAIN THE CONTRIBUTIONS FROM THE SECOND DERIVATIVES OF THE ELEMENTS.
C
            DO 80 IELL = 1, NGEL
               IEL     = IELING( IELL )
               IG      = ISLGRP( IELL )
               NVAREL  = ISTAEV( IEL + 1 ) - ISTAEV( IEL )
               LL      = ISTAEV( IEL )
               IF ( GXEQX( IG ) ) THEN
                  GDASH = ESCALE( IELL ) * GSCALE( IG )
               ELSE
                  GDASH = ESCALE( IELL ) * GSCALE( IG ) * GVALS2( IG)
               END IF
               IF ( INTREP( IEL ) ) THEN
                  NIN      = INTVAR( IEL + 1 ) - INTVAR( IEL )
                  DO 50 KK = 1, NVAREL
C
C  THE IEL-TH ELEMENT HESSIAN HAS AN INTERNAL REPRESENTATION.
C  SET WK AS THE K-TH COLUMN OF THE IDENTITY MATRIX.
C
                     WK( KK ) = ONE
C
C  GATHER WK INTO ITS INTERNAL VARIABLES, WKB.
C
                     CALL RANGES( IEL, .FALSE., WK, WKB, NVAREL, NIN )
                     WK( KK ) = ZERO
C
C  MULTIPLY THE INTERNAL VARIABLES BY THE ELEMENT HESSIAN.
C  CONSIDER THE FIRST COLUMN OF THE HESSIAN.
C
                     IELHST = ISTADH( IEL )
                     PI     = WKB( 1 )
CDIR$ IVDEP
                     DO 20 IROW     = 1, NIN
                        IJHESS      = ISYMMH( 1, IROW ) + IELHST
                        WKC( IROW ) = PI * HUVALS( IJHESS )
  20                 CONTINUE
C
C  NOW CONSIDER THE REMAINING COLUMNS OF THE HESSIAN.
C
                     DO 40 JCOL = 2, NIN
                        PI      = WKB( JCOL )
CDIR$ IVDEP
                        DO 30 IROW     = 1, NIN
                           IJHESS      = ISYMMH( JCOL, IROW ) + IELHST
                           WKC( IROW ) = WKC( IROW ) +
     *                                   PI * HUVALS( IJHESS )
  30                    CONTINUE
  40                 CONTINUE
C
C  ADD THE KK-TH DIAGONAL OF THE IEL-TH ELEMENT HESSIAN.
C
                     J         = IELVAR( LL )
                     LL        = LL + 1
CS                   DIAG( J ) = DIAG( J ) + GDASH * SDOT( NIN, WKB, 1,
CD                   DIAG( J ) = DIAG( J ) + GDASH * DDOT( NIN, WKB, 1,
     *                                                          WKC, 1 )
   50             CONTINUE
               ELSE
C
C  THE IEL-TH ELEMENT HESSIAN HAS NO INTERNAL REPRESENTATION.
C
                  IELHST = ISTADH( IEL )
CDIR$ IVDEP
                  DO 70 IROW   = 1, NVAREL
                     IJHESS    = ISYMMD( IROW ) + IELHST
                     J         = IELVAR( LL )
                     LL        = LL + 1
                     DIAG( J ) = DIAG( J ) + GDASH * HUVALS( IJHESS )
   70             CONTINUE
               END IF
   80       CONTINUE
C
C  IF THE GROUP IS NON-TRIVIAL, ADD ON RANK-ONE FIRST ORDER TERMS.
C
            DO 120 IG = 1, NG
               IF ( .NOT. GXEQX( IG ) ) THEN
                  G2DASH = GSCALE( IG ) * GVALS3( IG )
CDIR$ IVDEP
                  DO 110 K     = ISTAGV( IG ), ISTAGV( IG + 1 ) - 1
                     L         = ISVGRP( K )
                     DIAG( L ) = DIAG( L ) + G2DASH *
     *                           GRJAC( IVALJR( K ) ) ** 2
  110             CONTINUE
               END IF
  120       CONTINUE
C
C  TAKE THE ABSOLUTE VALUES OF ALL THE DIAGONAL ENTRIES, SETTING
C  ENSURING THAT ALL ENTRIES ARE LARGER THAN THE TOLERANCE SMALLH.
C
            DO 210 I     = 1, N
               DIAG( I ) = MAX( SMALLH, ABS( DIAG( I ) ) )
  210       CONTINUE
         END IF
C
C  USE THE DIAGONALS TO CALCULATE A SCALED NORM OF THE GRADIENT.
C
         IF ( DPRCND ) THEN
CDIR$ IVDEP
            DO 310 J  = 1, NVAR
               QGNORM = QGNORM + DGRAD( J ) *
     *                  ( DGRAD( J ) / DIAG( IVAR( J ) ) )
  310       CONTINUE
            QGNORM = SQRT( QGNORM )
         ELSE
            QGNORM = PJGNRM
         END IF
      END IF
      RETURN
C
C  END OF GTPGR.
C
      END
