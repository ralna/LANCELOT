C ** Correction report.
C ** Correction 1. 06/08/93: Routines *DNRM, *SETVL, *SETVI appended
C                            from linpac.f 
C ** Correction 2. 13/01/94: 2 lines added
C ** Correction 3. 13/01/94: 2 lines added
C ** Correction 4. 04/01/2000: 1 line altered
C ** Correction 5. 04/01/2000: 1 line removed
C ** End of Correction report.
C  THIS VERSION: 04/01/2000 AT 16:00:00 PM.
C  MODIFIED BY PH TOINT 30/04/1992.
CS    SUBROUTINE SSECNT( N     , INVAR , LINVAR, ISTAEV, LSTAEV, INTVAR,
CD    SUBROUTINE DSECNT( N     , INVAR , LINVAR, ISTAEV, LSTAEV, INTVAR,
     *                   LNTVAR, INTREP, LINTRE, ISTADH, LSTADH, NEL   ,
     *                   NINVAR, FUVALS, LFUVAL, ICALCF, LCALCF, NCALCF,
     *                   S , Y , WORK  , LWORK , IUPDAT, IDEBUG, IOUT  ,
     *                   RANGES )
      INTEGER            N     , LINVAR, LSTAEV, LSTADH, NEL   , NINVAR
      INTEGER            LCALCF, LINTRE, LNTVAR, LWORK , IUPDAT, IDEBUG
      INTEGER            NCALCF, IOUT  , LFUVAL
      INTEGER            INVAR ( LINVAR       ), ISTAEV( LSTAEV       )
      INTEGER            INTVAR( LNTVAR       ), ISTADH( LSTADH       )
      INTEGER            ICALCF( LCALCF       )
CS    REAL               FUVALS( LFUVAL       ), S     ( N            ),
CD    DOUBLE PRECISION   FUVALS( LFUVAL       ), S     ( N            ),
     *                   Y     ( NINVAR       ), WORK  ( LWORK        )
      LOGICAL            INTREP( LINTRE       )
      EXTERNAL           RANGES
C
C  COMPUTES THE SECANT UPDATE TO THE SECOND DERIVATIVE MATRIX.
C  FOR EACH ELEMENT FUNCTION.
C
C  IF IUPDAT = 1, THE B.F.G.S. UPDATE IS USED.
C  IF IUPDAT = 2, THE D.F.P. UPDATE IS USED.
C  IF IUPDAT = 3, THE P.S.B. UPDATE IS USED.
C  IF IUPDAT = 4, THE S.R.1 UPDATE IS USED.
C
C  NICK GOULD, 15TH OF MARCH 1990.
C  FOR CGT PRODUCTIONS.
C
      INTEGER           I, IEL, II, LL, IPOS  , IPOS1 , J , JJ, K , KK,
     *                  NIN   , NVAREL
CS    REAL              STS   , WTS   , WTW   , YTS   , YTY   , YJ, SJ,
CD    DOUBLE PRECISION  STS   , WTS   , WTW   , YTS   , YTY   , YJ, SJ,
     *                  SKIPR1, SKIPBD, ZERO  , WJ
      LOGICAL           INTRNL, PRNTER
      INTRINSIC         ABS
CS    EXTERNAL          SSETVL
CD    EXTERNAL          DSETVL
C
C  SET COMMON BLOCK.
C
      INTEGER           ITERCG, ITCGMX, NGEVAL, ISKIP , IFIXED, NBANDW
CS    REAL              ACCCG , RATIO , RADIUS, RADMAX, FINDMX, CMA31
CD    DOUBLE PRECISION  ACCCG , RATIO , RADIUS, RADMAX, FINDMX, CMA31
CS    COMMON / SCOMSB / ACCCG , RATIO , RADIUS, RADMAX, FINDMX, CMA31 ,
CD    COMMON / DCOMSB / ACCCG , RATIO , RADIUS, RADMAX, FINDMX, CMA31 ,
     *                  ITERCG, ITCGMX, NGEVAL, ISKIP , IFIXED, NBANDW
CS    SAVE / SCOMSB /
CD    SAVE / DCOMSB /
C
C  SET CONSTANT REAL PARAMETERS.
C
CS    PARAMETER ( ZERO   = 0.0E+0, SKIPBD = 1.0E-8, SKIPR1 = 1.0E-8 )
CD    PARAMETER ( ZERO   = 0.0D+0, SKIPBD = 1.0E-8, SKIPR1 = 1.0D-8 )
      PRNTER = IDEBUG .GE. 4 .AND. IOUT .GT. 0
C
C  CONSIDER THE IEL-TH ELEMENT.
C
      DO 500 I   = 1, NCALCF
         IEL     = ICALCF( I )
         LL      = ISTAEV( IEL ) - 1
         KK      = ISTADH( IEL ) - 1
         JJ      = INTVAR( IEL ) - 1 - NEL
         NVAREL  = ISTAEV( IEL + 1 ) - ISTAEV( IEL )
         NIN     = INTVAR( IEL + 1 ) - INTVAR( IEL )
         INTRNL  = INTREP( IEL )
C
C  IF THE ELEMENT HAS AN INTERNAL REPRESENTATION, TRANSFORM S.
C
         IF ( INTRNL ) THEN
            DO 10 J      = 1, NVAREL
               WORK( J ) = S( INVAR( LL + J ) )
   10       CONTINUE
            CALL RANGES( IEL, .FALSE., WORK( 1 ),
     *                   WORK( NVAREL + 1 ), NVAREL, NIN )
         END IF
C
C  COMPUTE SCALARS FOR THE BROYDEN-FLETCHER-GOLDFARB-SHANNO
C  AND DAVIDON-FLETCHER-POWELL UPDATES.
C
         IF ( IUPDAT .LE. 2 ) THEN
            YTS     = ZERO
            YTY     = ZERO
            IF ( INTRNL ) THEN
               DO 20 J = 1, NIN
                  YJ   = Y( JJ + J )
                  YTS  = YTS + YJ * WORK( NVAREL + J )
                  YTY  = YTY + YJ * YJ
   20          CONTINUE
            ELSE
               DO 30 J = 1, NVAREL
                  YJ   = Y( JJ + J )
                  YTS  = YTS + YJ * S( INVAR( LL + J ) )
                  YTY  = YTY + YJ * YJ
   30          CONTINUE
            END IF
            IF ( YTS .LE. SKIPBD * YTY ) THEN
               IF ( IUPDAT .EQ. 1 ) THEN
                  IF ( PRNTER ) WRITE( IOUT, 2010 ) IEL
               ELSE
                  IF ( PRNTER ) WRITE( IOUT, 2020 ) IEL
               END IF
               ISKIP  = ISKIP + 1
               GO TO 500
            END IF
         END IF
C
C  CALCULATE THE PRODUCT OF THE ELEMENT HESSIAN WITH THE VECTOR S AND
C  FOR THE DFP, PSB AND SR1 UPDATES, ADD THE ANSWER TO Y.
C  FIRST INITIALIZE WORK.
C
         IF ( IUPDAT .EQ. 1 ) THEN
CS          CALL SSETVL( NIN, WORK, 1, ZERO )
CD          CALL DSETVL( NIN, WORK, 1, ZERO )
         ELSE
            DO 40 IPOS      = 1, NIN
               WORK( IPOS ) = - Y( JJ + IPOS )
   40       CONTINUE
         END IF
         DO 70 IPOS = 1, NIN
            IF ( INTRNL ) THEN
               SJ = WORK( NVAREL + IPOS )
            ELSE
               SJ = S( INVAR( LL + IPOS ) )
            END IF
            K = KK + IPOS * ( IPOS - 1 ) / 2
C
C  FORM THE PRODUCT OF THE IPOS-TH COMPONENT WITH THE ELEMENT HESSIAN.
C
            DO 50 II      = 1, IPOS
               K          = K + 1
               WORK( II ) = WORK( II ) + SJ * FUVALS( K )
   50       CONTINUE
            IPOS1 = IPOS + 1
            IF ( IPOS1 .LE. NIN ) THEN
               DO 60 II      = IPOS1, NIN
                  K          = K + ( II - 1 )
                  WORK( II ) = WORK( II ) + SJ * FUVALS( K )
   60          CONTINUE
            END IF
   70    CONTINUE
C
C  COMPUTE THE INNER PRODUCT OF THIS VECTOR WITH S.
C
         WTS = ZERO
         IF ( INTRNL ) THEN
            DO 80 J = 1, NIN
               WTS  = WTS + WORK( J ) * WORK( NVAREL + J )
   80       CONTINUE
         ELSE
            DO 90 J = 1, NIN
               WTS  = WTS + WORK( J ) * S( INVAR( LL + J ) )
   90       CONTINUE
         END IF
C
C  COMPUTE S(TRANS) S FOR ALL UPDATES, BUT SYMMETRIC RANK ONE.
C
         IF ( IUPDAT .NE. 4 ) THEN
            STS = ZERO
            IF ( INTRNL ) THEN
               DO 100 J = 1, NIN
                  STS  = STS + WORK( NVAREL + J )**2
  100          CONTINUE
            ELSE
               DO 105 J = 1, NIN
                  STS  = STS + S( INVAR( LL + J ) )**2
  105          CONTINUE
            END IF
C
C  SKIP THE POSITIVE DEFINITE UPDATES IF THE ELEMENT HESSIAN IS NOT
C  POSITIVE DEFINITE ENOUGH, DUE TO ROUNDING ERRORS.
C
            IF ( WTS .LE. SKIPBD * STS .AND. IUPDAT .NE. 3 ) THEN
               IF ( IUPDAT .EQ. 1 ) THEN
                  IF ( PRNTER ) WRITE( IOUT, 2050 ) IEL
               ELSE
                  IF ( PRNTER ) WRITE( IOUT, 2060 ) IEL
               END IF
               ISKIP  = ISKIP + 1
               GO TO 500
            END IF
         END IF
C
C  BROYDEN FLETCHER GOLDFARB SHANNO UPDATE.
C
         IF ( IUPDAT .EQ. 1 ) THEN
            DO 120 J           = 1, NIN
               YJ              = Y( JJ + J ) / YTS
               WJ              = WORK( J ) / WTS
               DO 110 K        = 1, J
                  KK           = KK + 1
                  FUVALS( KK ) = FUVALS( KK ) - WJ * WORK( K ) +
     *                                          YJ * Y( JJ + K )
  110          CONTINUE
  120       CONTINUE
         END IF
C
C  DAVIDON FLETCHER POWELL UPDATE.
C
         IF ( IUPDAT .EQ. 2 ) THEN
            DO 220 J           = 1, NIN
               WJ              =  Y( JJ + J ) / YTS
               YJ              = ( WORK( J ) - WJ * WTS ) / YTS
               DO 210 K        = 1, J
                  KK           = KK + 1
                  FUVALS( KK ) = FUVALS( KK ) - WJ * WORK( K ) -
     *                                          YJ * Y( JJ + K )
  210          CONTINUE
  220       CONTINUE
         END IF
C
C  POWELL SYMMETRIC BROYDEN UPDATE.
C
         IF ( IUPDAT .EQ. 3 ) THEN
C
C  UPDATE FOR AN INTERNAL ELEMENTAL REPRESENTATION.
C
            IF ( INTRNL ) THEN
               IF ( STS .NE. ZERO ) THEN
                  DO 330 J           = 1, NIN
                     WJ              = - WORK( NVAREL + J ) / STS
                     SJ              = ( - WORK( J ) - WJ * WTS ) / STS
                     DO 320 K        = 1, J
                        KK           = KK + 1
                        FUVALS( KK ) = FUVALS( KK ) + WJ * WORK( K )
     *                                 + SJ * WORK( NVAREL + K )
  320                CONTINUE
  330             CONTINUE
               ELSE
                  IF ( PRNTER ) WRITE( IOUT, 2030 ) IEL
                  ISKIP = ISKIP + 1
                  GO TO 500
               END IF
            ELSE
C
C  UPDATE FOR AN ELEMENTAL REPRESENTATION.
C
               IF ( STS .NE. ZERO ) THEN
                  DO 360 J           = 1, NIN
                     WJ              = - S( INVAR( LL + J ) ) / STS
                     SJ              = ( - WORK( J ) - WJ * WTS ) / STS
                     DO 350 K        = 1, J
                        KK           = KK + 1
                        FUVALS( KK ) = FUVALS( KK ) + WJ * WORK( K )
     *                                 + SJ * S( INVAR( LL + K ) )
  350                CONTINUE
  360             CONTINUE
               ELSE
                  IF ( PRNTER ) WRITE( IOUT, 2030 ) IEL
                  ISKIP = ISKIP + 1
                  GO TO 500
               END IF
            END IF
         END IF
C
C  SYMMETRIC RANK 1 UPDATE.
C
         IF ( IUPDAT .EQ. 4 ) THEN
            WTW      = ZERO
            DO 410 J = 1, NIN
               WJ    = WORK( J )
               WTW   = WTW + WJ * WJ
  410       CONTINUE
            IF ( ABS( WTS ) .LE. SKIPR1 * WTW ) THEN
               IF ( PRNTER ) WRITE( IOUT, 2040 ) IEL
               ISKIP = ISKIP + 1
               GO TO 500
            ELSE
               DO 430 J           = 1, NIN
                  WJ              = WORK( J ) / WTS
                  DO 420 K        = 1, J
                     KK           = KK + 1
                     FUVALS( KK ) = FUVALS( KK ) - WJ * WORK( K )
  420             CONTINUE
  430          CONTINUE
            END IF
         END IF
  500 CONTINUE
      RETURN
C
C  NON-EXECUTABLE STATEMENTS.
C
 2010 FORMAT( /, ' BFGS update skipped in element ', I5,
     *           ' Y(TRANS) S is negative' )
 2020 FORMAT( /, ' DFP update skipped in element ', I5,
     *           ' Y(TRANS) S is negative' )
 2030 FORMAT( /, ' PSB update skipped in element ', I5,
     *           ' S(TRANS) S is zero ' )
 2040 FORMAT( /, ' SR1 update skipped in element ', I5,
     *           ' W(TRANS) S is too small ' )
 2050 FORMAT( /, ' BFGS update skipped in element ', I5,
     *           ' S(TRANS) B S is negative' )
 2060 FORMAT( /, ' DFP update skipped in element ', I5,
     *           ' S(TRANS) B S is negative' )
C
C  END OF SUBROUTINE SECNT.
C
      END
C  THIS VERSION: 13/01/1994 AT 04:41:23 PM.
CS    SUBROUTINE SSCALH( INITH, N, NEL, NCALCF, LHXI, ISTAEV, LSTAEV,
CD    SUBROUTINE DSCALH( INITH, N, NEL, NCALCF, LHXI, ISTAEV, LSTAEV,
     *                   ISTADH, LSTADH, ICALCF, LCALCF, INTVAR,
     *                   LNTVAR, IELVAR, LELVAR, ISYMMD, MAXSZH,
     *                   INTREP, LINTRE, FUVALS, LFUVAL,
     *                   WK, LNWK, WKB, LNWKB, S, Y, LY, RANGES )
      INTEGER          N, NEL, NCALCF, LHXI, LFUVAL, LSTAEV, LSTADH
      INTEGER          LCALCF, MAXSZH, LELVAR, LNTVAR, LY, LINTRE
      INTEGER          LNWK, LNWKB
      LOGICAL          INITH
      INTEGER          ISTAEV( LSTAEV ), ISTADH( LSTADH )
      INTEGER          ICALCF( LCALCF ), INTVAR( LNTVAR )
      INTEGER          IELVAR( LELVAR ), ISYMMD( MAXSZH )
      LOGICAL          INTREP( LINTRE )
CS    REAL             FUVALS( LFUVAL ), WK( LNWK ), WKB( LNWKB ),
CD    DOUBLE PRECISION FUVALS( LFUVAL ), WK( LNWK ), WKB( LNWKB ),
     *                 S( N ), Y( LY )
      EXTERNAL         RANGES
C
C  INITIALIZE THE APPROXIMATE SECOND DERIVATIVE MATRIX FOR EACH
C  NONLINEAR ELEMENT AS A (SCALED) IDENTITY MATRIX, IF INITH IS
C  .TRUE., AND SCALE THESE MATRICES TO SATISFY THE WEAK SECANT
C  EQUATION, IF INITH IS .FALSE.
C
C  NICK GOULD, 9TH JULY 1990.
C  FOR CGT PRODUCTIONS.
C
C  LOCAL VARIABLES.
C
      INTEGER          I, J, IEL, NIN, LHUVAL, NEL1, NVAREL
CS    REAL             ZERO, ONE, SCALE, YTS, STHS, SI
CD    DOUBLE PRECISION ZERO, ONE, SCALE, YTS, STHS, SI
C
C  EXTERNAL SUBROUTINES AND FUNCTIONS USED.
C
CS    EXTERNAL         SINISH, SSETVL
CD    EXTERNAL         DINISH, DSETVL
C
C  SET CONSTANT REAL PARAMETERS.
C
CS    PARAMETER ( ZERO   = 0.0E+0, ONE    = 1.0E+0 )
CD    PARAMETER ( ZERO   = 0.0D+0, ONE    = 1.0D+0 )
C
C  IF A SECANT METHOD IS TO BE USED, INITIALIZE THE SECOND
C  DERIVATIVES OF EACH ELEMENT AS A SCALED IDENTITY MATRIX.
C
      NEL1   = NEL + 1
      LHUVAL = ISTADH( NEL + 1 ) - ISTADH( 1 )
      IF ( INITH ) THEN
C
C  SET ALL VALUES TO ZERO.
C
CS       CALL SSETVL( LHUVAL, FUVALS( LHXI + 1 ), 1, ZERO )
CD       CALL DSETVL( LHUVAL, FUVALS( LHXI + 1 ), 1, ZERO )
C
C  RESET THE DIAGONALS TO THE VALUE SCALE.
C
         SCALE     = ONE
         DO 10 IEL = 1, NEL
            NIN    = INTVAR( IEL + 1 ) - INTVAR( IEL )
CS          CALL SINISH( IEL, NIN, SCALE, ISTADH,
CD          CALL DINISH( IEL, NIN, SCALE, ISTADH,
     *                   FUVALS, ISYMMD, MAXSZH )
   10    CONTINUE
      ELSE
C
C  AT THE END OF THE FIRST SUCCESSFUL ITERATION, SCALE THE INITIAL
C  SECOND DERIVATIVE MATRIX FOR EACH ELEMENT SO AS TO SATISFY
C  THE WEAK SECANT CONDITION OF SHANNO AND PHUA.
C
         DO 200 I  = 1, NCALCF
            IEL    = ICALCF( I )
            NIN    = INTVAR( IEL + 1 ) - INTVAR( IEL )
            NVAREL = ISTAEV( IEL + 1 ) - ISTAEV( IEL )
            YTS    = ZERO
            STHS   = ZERO
C
C  IF THE ELEMENT HAS AN INTERNAL REPRESENTATION, TRANSFORM S
C  INTO ITS INTERNAL VARIABLES, WKB.
C
            IF ( INTREP( IEL ) ) THEN
CDIR$ IVDEP
               DO 110 J   = 1, NVAREL
                  WK( J ) = S( IELVAR( ISTAEV( IEL ) + J - 1 ) )
  110          CONTINUE
               CALL RANGES( IEL, .FALSE., WK, WKB, NVAREL, NIN )
C
C  COMPUTE THE SCALARS YTS = Y(TRANS) S AND STHS =
C  S(TRANS) H S, REMEMBERING THAT H = I HERE.
C
               DO 120 J = 1, NIN
                  SI    = WKB( J )
                  YTS   = YTS + Y( INTVAR( IEL ) - NEL1 + J ) * SI
                  STHS  = STHS + SI * SI
  120          CONTINUE
            ELSE
               DO 130 J = 1, NVAREL
                  SI    = S( IELVAR( ISTAEV( IEL ) + J - 1 ) )
                  YTS   = YTS + Y( INTVAR( IEL ) - NEL1 + J ) * SI
                  STHS  = STHS + SI * SI
  130          CONTINUE
            END IF
C
C  SCALE THE ELEMENT HESSIANS BY THE QUANTITY YTS / STHS AS SUGGESTED
C  BY SHANNO AND PHUA.
C
            SCALE = YTS / STHS
CS          CALL SINISH( IEL, NIN, SCALE, ISTADH, FUVALS,
CD          CALL DINISH( IEL, NIN, SCALE, ISTADH, FUVALS,
     *                   ISYMMD, MAXSZH )
  200    CONTINUE
      END IF
      RETURN
C
C  END OF SUBROUTINE SCALH.
C
      END
C  THIS VERSION: 13/01/1994 AT 04:41:23 PM.
CS    SUBROUTINE SFDGRD( N, NEL, IELVAR, LELVAR, ISTAEV, LSTAEV,
CD    SUBROUTINE DFDGRD( N, NEL, IELVAR, LELVAR, ISTAEV, LSTAEV,
     *                   ISELTS, LNELTS, ISPTRS, LNPTRS, IELING,
     *                   LELING, ISVSET, LSVSET, ISET  , LSET  ,
     *                   NSETS , INVSET, LNVSET, ICALCF, LCALCF,
     *                   NCALCF, INTVAR, LNTVAR, ISWTRA,
     *                   LSWTRA, ISIWTR, LSIWTR, ITYPEE, LTYPEE,
     *                   NTYPE , IWTRAN, LIWTRA, WTRANS, LWTRAN,
     *                   X , XT, FUVALS, LFUVAL, CENTRL, IGETFD )
      INTEGER            N, NEL, LELVAR, LSTAEV, LELING, LTYPEE, NTYPE
      INTEGER            LSWTRA, LSIWTR, LIWTRA, LWTRAN, NSETS
      INTEGER            LNELTS, LNPTRS, LSVSET, LSET  , LNVSET
      INTEGER            IGETFD, LCALCF, NCALCF, LNTVAR, LFUVAL
      LOGICAL            CENTRL
      INTEGER            ISELTS( LNELTS       ), ISPTRS( LNPTRS        )
      INTEGER            IELVAR( LELVAR       ), ISTAEV( LSTAEV        )
      INTEGER            IELING( LELING       ), ISVSET( LSVSET        )
      INTEGER            ISET  ( LSET         ), INVSET( LNVSET        )
      INTEGER            ICALCF( LCALCF       ), INTVAR( LNTVAR        )
      INTEGER            ISWTRA( LSWTRA       ), IWTRAN( LIWTRA        )
      INTEGER            ISIWTR( LSIWTR       ), ITYPEE( LTYPEE        )
CS    REAL               X     ( N            ), XT    ( N             )
CS    REAL               FUVALS( LFUVAL       ), WTRANS( LWTRAN        )
CD    DOUBLE PRECISION   X     ( N            ), XT    ( N             )
CD    DOUBLE PRECISION   FUVALS( LFUVAL       ), WTRANS( LWTRAN        )
C
C  --------------------------------------------------------------------
C
C  Obtain finite-difference estimates of the first derivatives of the
C  nonlinear element functions.
C
C  Nick Gould, 23rd September 1991.
C  For CGT productions.
C
C  --------------------------------------------------------------------
C
      INTEGER          I , J , K , L , IEL   , IPT   , ITYPE
      INTEGER          NINVAR, LWFREE, LIWFRE, IELL  , IVAR  , INTV
      LOGICAL          BACKWD
CS    REAL             DIFF  , ZERO  , ONE   , TWODIF, TWO
CD    DOUBLE PRECISION DIFF  , ZERO  , ONE   , TWODIF, TWO
CS    PARAMETER      ( ZERO = 0.0E+0 , ONE = 1.0E+0  , TWO = 2.0E+0 )
CD    PARAMETER      ( ZERO = 0.0D+0 , ONE = 1.0D+0  , TWO = 2.0D+0 )
C
C  COMMON VARIABLES.
C
CS    REAL              EPSMCH, EPSNEG, TINY, BIG
CD    DOUBLE PRECISION  EPSMCH, EPSNEG, TINY, BIG
CS    COMMON / SMACHN / EPSMCH, EPSNEG, TINY, BIG
CD    COMMON / DMACHN / EPSMCH, EPSNEG, TINY, BIG
C
C  INTRINSIC FUNCTIONS.
C
C ** Correction 5. 04/01/2000: 1 line removed
      SAVE
      IF ( IGETFD .GT. NSETS ) THEN
         IGETFD = - 1
         RETURN
      END IF
C     WRITE( 6, * ) ' XT ', ( XT( I ), I = 1, N )
C     WRITE( 6, * ) ' FUVALS ', ( FUVALS( I ), I = 1, NEL )
C
C  Calculate the difference intervals.
C
      IF ( CENTRL ) THEN
         DIFF   = EPSMCH ** 0.33333
         TWODIF = TWO * DIFF
      ELSE
         DIFF   = EPSMCH ** 0.5
      END IF
C
C  ---------------------------------------------------------------------
C  Compute the finite differences corresponding to the variables
C  from the IGETFD-th set.
C  ---------------------------------------------------------------------
C
      IF ( IGETFD .GT. 0 ) THEN
         DO 490 I = ISVSET( IGETFD ), ISVSET( IGETFD + 1 ) - 1
            IVAR  = ISET( I )
C
C  Loop over the elements which use variable IVAR.
C  The elements are obtained from a linked-list.
C
            IPT = ISPTRS( IVAR )
            IF ( IPT .GE. 0 ) THEN
               IELL = ISELTS( IVAR )
  410          CONTINUE
               IEL  = IELING( IELL )
C
C  If the element has an internal representation, check that the
C  variable IVAR belongs in the "independence" set.
C
               ITYPE = ITYPEE( IEL )
               IF ( ITYPE .GT. 0 ) THEN
                  LIWFRE = ISIWTR( ITYPE )
                  NINVAR = IWTRAN( LIWFRE     )
                  DO 420 J = 1, NINVAR
                     K = J - 1
                     L = IWTRAN( LIWFRE + NINVAR + 1 + J ) - 1
                     IF ( IVAR .EQ. IELVAR( ISTAEV( IEL ) + L ) )
     *                  GO TO 440
  420             CONTINUE
                  GO TO 470
               ELSE
C
C  Find which internal variable is used.
C
                  K = 0
                  DO 430 J = ISTAEV( IEL ), ISTAEV( IEL + 1 ) - 1
                     IF ( IELVAR( J ) .EQ. IVAR ) GO TO 440
                     K = K + 1
  430             CONTINUE
               END IF
  440          CONTINUE
C
C  Form a central difference.
C
               INTV = INTVAR( IEL )
               L    = INTV + K
               IF ( CENTRL ) THEN
                  IF ( BACKWD ) THEN
                     IF ( INVSET( IEL ) .EQ. 2 ) THEN
                        FUVALS( L ) = ( FUVALS( L ) -
     *                                  FUVALS( IEL ) ) / TWODIF
C                       write( 6, * ) ' el ', IEL, ' gradient is ',
C    *                     FUVALS( L )

C  Ensure that derivatives for repeated variables are set to zero.
C
                        IF ( ITYPE .EQ. 0 ) THEN
                           K = K + 1
                           L = K
                           DO 450 J = ISTAEV( IEL     ) + L,
     *                                ISTAEV( IEL + 1 ) - 1
                              IF ( IELVAR( J ) .EQ. IVAR )
     *                           FUVALS( INTV + K ) = ZERO
                              K = K + 1
  450                      CONTINUE
                        END IF
                        INVSET( IEL ) = 0
                     END IF
                  ELSE
                     IF ( INVSET( IEL ) .EQ. 1 ) THEN
C                       write( 6, * ) ' el ', IEL,
C    *                     ' gradient would be ',
C    *                     ( FUVALS( IEL ) - FUVALS( L ) ) / DIFF
                        FUVALS( L )   = FUVALS( IEL )
                        INVSET( IEL ) = 2
                     END IF
                  END IF
C
C  Form a forward difference.
C
               ELSE
                  IF ( INVSET( IEL ) .EQ. 1 ) THEN
C                    write(6,*) ' ELEMENT ', IEL, ' F = ', FUVALS( L )
                     FUVALS( L ) = ( FUVALS( IEL ) -
     *                               FUVALS( L ) ) / DIFF
C                       write(6,*) ' NEW F = ', FUVALS( IEL ),
C    *                  ' grad ', FUVALS( L )
C
C  Ensure that derivatives for repeated variables are set to zero.
C
                     IF ( ITYPE .EQ. 0 ) THEN
                        K = K + 1
                        L = K
                        DO 460 J = ISTAEV( IEL     ) + L,
     *                             ISTAEV( IEL + 1 ) - 1
                           IF ( IELVAR( J ) .EQ. IVAR )
     *                          FUVALS( INTV + K ) = ZERO
                           K = K + 1
  460                   CONTINUE
                     END IF
                     INVSET( IEL ) = 0
                  END IF
               END IF
  470          CONTINUE
C
C  See if further variables use variable IVAR.
C
               IF ( IPT .GT. 0 ) THEN
                  IELL = ISELTS( IPT )
                  IPT  = ISPTRS( IPT )
                  GO TO 410
              END IF
            END IF
C
C  Reset the variables from the IGETFD-th set to their initial values.
C
            IF ( CENTRL ) THEN
               IF ( BACKWD ) THEN
                  XT( IVAR ) = X( IVAR )
               ELSE
                  XT( IVAR ) = X( IVAR ) - DIFF
               END IF
            ELSE
               XT( IVAR ) = X( IVAR )
            END IF
  490    CONTINUE
      ELSE
C
C  ---------------------------------------------------------------------
C  Initialise the point, XT,  at which the elements are to be computed
C  as the current point X.
C  ---------------------------------------------------------------------
C
         IF ( CENTRL ) BACKWD = .TRUE.
         DO 510 IVAR   = 1, N
            XT( IVAR ) = X( IVAR )
  510    CONTINUE
C
C  Empty the list, INVSET, of elements which need to be re-evaluated.
C  Those marked - 1 need not be re-evaluated.
C
         DO 520 J       = 1, NEL
            INVSET( J ) = - 1
  520    CONTINUE
         DO 540 K = 1, NCALCF
            J     = ICALCF( K )
C           IF ( .NOT. CENTRL ) THEN
               DO 530 I = INTVAR( J ), INTVAR( J + 1 ) - 1
                  FUVALS( I ) = FUVALS( J )
  530          CONTINUE
C           END IF
            INVSET( J ) = 0
  540    CONTINUE
      END IF
C
C  Compute the next set of finite difference intervals.
C
      IF ( CENTRL ) THEN
         IF ( BACKWD ) THEN
            BACKWD = .FALSE.
            IGETFD = IGETFD + 1
         ELSE
            BACKWD = .TRUE.
         END IF
      ELSE
         IGETFD = IGETFD + 1
      END IF
C
C  If all the difference have been computed, prepare to return.
C
      IF ( IGETFD .GT. NSETS ) THEN
C
C  Run through the list of elements with internal variables,
C  transforming the differences by the matrix W(-T).
C

         IF ( NTYPE .GT. 0 ) THEN
            NCALCF = 0
            DO 620 IEL = 1, NEL
               ITYPE   = ITYPEE( IEL )
               IF  ( ITYPE .GT. 0 ) THEN
                  LIWFRE = ISIWTR( ITYPE  )
                  LWFREE = ISWTRA( ITYPE  )
                  NINVAR = IWTRAN( LIWFRE )
C
C  Transform the differences.
C
CS                CALL SGESLV( NINVAR, IWTRAN( LIWFRE + 2 ),
CD                CALL DGESLV( NINVAR, IWTRAN( LIWFRE + 2 ),
     *                   WTRANS( LWFREE ), FUVALS( INTVAR( IEL ) ) )
               END IF
C
C  Reset ICALCF to its original value.
C
               IF ( INVSET( IEL ) .NE. - 1 ) THEN
                  NCALCF = NCALCF + 1
                  ICALCF( NCALCF ) = IEL
               END IF
  620       CONTINUE
         END IF
         IGETFD = - 1
         RETURN
      END IF
C
C  Prepare to return to obtain additional element function values.
C  Compute the difference intervals for the IGETFD-th set.
C
      IF ( .NOT. ( CENTRL .AND. BACKWD ) ) THEN
         NCALCF = 0
         DO 680 I = ISVSET( IGETFD ), ISVSET( IGETFD + 1 ) - 1
            IVAR  = ISET( I )
            XT( IVAR ) = X( IVAR ) + DIFF
C
C  Loop over the elements which use variable IVAR.
C  The elements are obtained from a linked-list.
C
            IPT = ISPTRS( IVAR )
            IF ( IPT .GE. 0 ) THEN
               IELL = ISELTS( IVAR )
  630          CONTINUE
               IEL  = IELING( IELL )
               ITYPE = ITYPEE( IEL )
C
C  Check that the variable belongs to the "indendence" set of an
C  element with an internal representation.
C
               IF ( ITYPE .GT. 0 ) THEN
                  LIWFRE = ISIWTR( ITYPE )
                  NINVAR = IWTRAN( LIWFRE     )
                  DO 640 J = 1, NINVAR
                     K = IWTRAN( LIWFRE + NINVAR + 1 + J ) - 1
                     IF ( IVAR .EQ. IELVAR( ISTAEV( IEL ) + K ) )
     *                  GO TO 660
  640             CONTINUE
                  GO TO 670
               END IF
C
C  Flag the nonlinear elements which will need to be recalculated.
C
  660          CONTINUE
               IF ( INVSET( IEL ) .EQ. 0 ) THEN
                  INVSET( IEL ) = 1
                  NCALCF = NCALCF + 1
                  ICALCF( NCALCF ) = IEL
               END IF
  670          CONTINUE
               IF ( IPT .GT. 0 ) THEN
                  IELL = ISELTS( IPT )
                  IPT  = ISPTRS( IPT )
                  GO TO 630
               END IF
            END IF
  680    CONTINUE
C        write(6, * ) ' EVALUATE ', ( ICALCF( I ), I = 1, NCALCF )
      END IF
      RETURN
C
C  End of subroutine FDGRD.
C
      END
C  THIS VERSION: 13/01/1994 AT 04:41:23 PM.
CS    SUBROUTINE SGELIM( M , N , IPVT  , JCOL  , A      )
CD    SUBROUTINE DGELIM( M , N , IPVT  , JCOL  , A      )
      INTEGER            M , N
      INTEGER            IPVT  ( M    ), JCOL  ( N     )
CS    REAL               A     ( M     , N     )
CD    DOUBLE PRECISION   A     ( M     , N     )
C
C  Perform the first M steps of Gaussian Elimination with
C  complete pivoting on the M by N ( M <= N) matrix A.
C
C  Nick Gould, 23rd September 1991.
C  For CGT productions.
C
      INTEGER            I , J , K     , IPIVOT, JPIVOT
CS    REAL               APIVOT, ONE   , ATEMP
CD    DOUBLE PRECISION   APIVOT, ONE   , ATEMP
CS    PARAMETER        ( ONE = 1.0E+0 )
CD    PARAMETER        ( ONE = 1.0D+0 )
C
C  Initialize the column indices.
C
      DO 10 J = 1, N
         JCOL( J ) = J
   10 CONTINUE
C
C  Main loop.
C
      DO 100 K = 1, M
C
C  Compute the K-th pivot.
C
         APIVOT = - ONE
         DO 30 J = K, N
            DO 20 I = K, M
               IF ( ABS( A( I, J ) ) .GT. APIVOT ) THEN
                  APIVOT = ABS( A( I, J ) )
                  IPIVOT = I
                  JPIVOT = J
               END IF
   20       CONTINUE
   30    CONTINUE
C
C  Interchange rows I and IPIVOT.
C
         IPVT( K ) = IPIVOT
         IF ( IPIVOT .GT. K ) THEN
            DO 40 J   = K, N
               ATEMP          = A( IPIVOT, J )
               A( IPIVOT, J ) = A( K     , J )
               A( K     , J ) = ATEMP
   40       CONTINUE
         END IF
C
C  Interchange columns J and JPIVOT.
C
         IF ( JPIVOT .GT. K ) THEN
            J              = JCOL( JPIVOT )
            JCOL( JPIVOT ) = JCOL( K      )
            JCOL( K      ) = J
            DO 50 I = 1, M
               ATEMP          = A( I, JPIVOT )
               A( I, JPIVOT ) = A( I, K )
               A( I, K      ) = ATEMP
   50       CONTINUE
         END IF
C
C  Perform the elimination.
C
         APIVOT       = A( K, K )
         DO 70 I      = K + 1, M
            ATEMP     = A( I, K ) / APIVOT
            A( I, K ) = ATEMP
            DO 60 J      = K + 1, N
               A( I, J ) = A( I, J ) - ATEMP * A( K, J )
   60       CONTINUE
   70    CONTINUE
  100 CONTINUE
      RETURN
C
C  End of subroutine GELIM.
C
      END
C  THIS VERSION: 13/01/1994 AT 04:41:23 PM.
CS    SUBROUTINE SGESLV( M     , IPVT  , A , X  )
CD    SUBROUTINE DGESLV( M     , IPVT  , A , X  )
      INTEGER            M
      INTEGER            IPVT  ( M    )
CS    REAL               A     ( M     , M     ), X    ( M       )
CD    DOUBLE PRECISION   A     ( M     , M     ), X    ( M       )
C
C  Solve the equations A(T)x = b. The vector b is input in X.
C  The LU factors of P A are input in A; The permutation P is stored
C  in IPVT. The solution x is output in X.
C
C  Nick Gould, 23rd September 1991.
C  For CGT productions.
C
      INTEGER            I , K
CS    REAL               XTEMP, ZERO
CD    DOUBLE PRECISION   XTEMP, ZERO
CS    PARAMETER        ( ZERO = 0.0E+0 )
CD    PARAMETER        ( ZERO = 0.0D+0 )
C
C  Solve U(T)y = b. The vector b is input in X; y is output in X.
C
      DO 20 K  = 1, M
         XTEMP = ZERO
         DO 10 I   = 1, K - 1
            XTEMP  = XTEMP + A( I, K ) * X( I )
   10    CONTINUE
         X( K ) = ( X( K ) - XTEMP ) / A( K, K )
   20 CONTINUE
C
C  Solve L(T) x = y. The vector y is input in X; x is output in X.
C
      DO 40 K  = M - 1, 1, - 1
         XTEMP = ZERO
         DO 30 I   = K + 1, M
            XTEMP  = XTEMP + A( I, K ) * X( I )
   30    CONTINUE
         X( K ) = X( K ) - XTEMP
         I      = IPVT( K )
         IF ( I .NE. K ) THEN
            XTEMP  = X( I )
            X( I ) = X( K )
            X( K ) = XTEMP
         END IF
   40 CONTINUE
      RETURN
C
C  End of subroutine GESLV.
C
      END

C  THIS VERSION: 13/01/1994 AT 04:41:23 PM.
C
C  ** FOR THE CRAY 2, LINES STARTING 'CDIR$ IVDEP' TELL THE COMPILER TO
C     IGNORE POTENTIAL VECTOR DEPENDENCIES AS THEY ARE KNOWN TO BE O.K.
C
CS    SUBROUTINE SINISH( IEL, NIN, SCALE, ISTADH, HUVALS,
CD    SUBROUTINE DINISH( IEL, NIN, SCALE, ISTADH, HUVALS,
     *                   ISYMMD, MAXSZH )
C
C  SET THE DIAGONAL ENTRIES OF THE IEL-TH ELEMENT HESSIAN TO
C  THE VALUE SCALE.
C
C  NICK GOULD, 19TH OF MAY 1989.
C  FOR CGT PRODUCTIONS.
C
      INTEGER          IEL, NIN, MAXSZH
CS    REAL             SCALE
CD    DOUBLE PRECISION SCALE
      INTEGER          ISTADH( * ), ISYMMD( MAXSZH )
CS    REAL             HUVALS( * )
CD    DOUBLE PRECISION HUVALS( * )
C
C  LOCAL VARIABLES.
C
      INTEGER          IELHST, IROW
C
C  INITIALIZE ALL ENTRIES TO ZERO.
C
      IELHST = ISTADH( IEL )
CDIR$ IVDEP
      DO 10 IROW = 1, NIN
         HUVALS( ISYMMD( IROW ) + IELHST ) = SCALE
   10 CONTINUE
      RETURN
C
C  END OF SUBROUTINE INISH.
C
      END
C  THIS VERSION: 13/01/1994 AT 04:41:23 PM.
C
C  ** FOR THE CRAY 2, LINES STARTING 'CDIR$ IVDEP' TELL THE COMPILER TO
C     IGNORE POTENTIAL VECTOR DEPENDENCIES AS THEY ARE KNOWN TO BE O.K.
C
CS    SUBROUTINE SSYMMH( MAXSZH, ISYMMH, ISYMMD )
CD    SUBROUTINE DSYMMH( MAXSZH, ISYMMH, ISYMMD )
      INTEGER MAXSZH
      INTEGER ISYMMH( MAXSZH, MAXSZH ), ISYMMD( MAXSZH )
C
C  GIVEN A COLUMNWISE STORAGE SCHEME OF THE UPPER TRIANGLE OF A
C  SYMMETRIC MATRIX OF ORDER MAXSZH, COMPUTE THE POSITION OF THE
C  I,J-TH ENTRY OF THE SYMMETRIC MATRIX IN THIS SCHEME.
C
C  THE VALUE ISYMMH( I, J ) + 1 GIVES THE POSITION OF THE I,J-TH
C  ENTRY OF THE MATRIX IN THE UPPER TRIANGULAR SCHEME.
C
C  NICK GOULD, 10TH OF MAY 1989.
C  FOR CGT PRODUCTIONS.
C
      INTEGER I, J, K
      K       = 0
      DO 20 J = 1, MAXSZH
CDIR$ IVDEP
         DO 10 I           = 1, J - 1
            ISYMMH( I, J ) = K
            ISYMMH( J, I ) = K
            K              = K + 1
   10    CONTINUE
         ISYMMD( J )    = K
         ISYMMH( J, J ) = K
         K              = K + 1
   20 CONTINUE
      RETURN
C
C  END OF SYMMH.
C
      END
C  THIS VERSION: 13/01/1994 AT 04:41:23 PM.
CS    SUBROUTINE SPRGRA( N, X, G, XSCALE, BL, BU, GRAD,
CD    SUBROUTINE DPRGRA( N, X, G, XSCALE, BL, BU, GRAD,
     *                   IVAR, NVAR, PJGNRM )
C
C  COMPUTE THE PROJECTION OF THE GRADIENT INTO THE FEASIBLE BOX.
C
C  NICK GOULD, 10TH OF MAY 1989.
C  FOR CGT PRODUCTIONS.
C
      INTEGER          N, NVAR
CS    REAL             PJGNRM
CD    DOUBLE PRECISION PJGNRM
      INTEGER          IVAR( N )
CS    REAL             X( N ), G( N ), BL( N ), BU( N ), GRAD( N ),
CD    DOUBLE PRECISION X( N ), G( N ), BL( N ), BU( N ), GRAD( N ),
     *                 XSCALE( N )
C
C  LOCAL VARIABLES.
C
      INTEGER          I
CS    REAL             GI, ZERO
CD    DOUBLE PRECISION GI, ZERO
C
C  COMMON VARIABLES.
C
CS    REAL             EPSMCH, EPSNEG, TINY, BIG
CD    DOUBLE PRECISION EPSMCH, EPSNEG, TINY, BIG
CS    COMMON / SMACHN / EPSMCH, EPSNEG, TINY, BIG
CD    COMMON / DMACHN / EPSMCH, EPSNEG, TINY, BIG
C ** Correction 2. 13/01/94: 2 lines added
CS    SAVE   / SMACHN /
CD    SAVE   / DMACHN /
C ** Correction 2. 13/01/94: end of correction **
C
C  MACHINE FUNCTIONS.
C
C ** Correction 4. 04/01/2000: 1 line altered
      INTRINSIC        ABS, MIN
C
C  SET CONSTANT REAL PARAMETERS.
C
CS    PARAMETER ( ZERO   = 0.0E+0 )
CD    PARAMETER ( ZERO   = 0.0D+0 )
      NVAR    = 0
      PJGNRM  = ZERO
      DO 20 I = 1, N
         GI   = G( I ) * XSCALE( I )
         IF ( GI .EQ. ZERO ) GO TO 20
C
C  COMPUTE THE PROJECTION OF THE GRADIENT WITHIN THE BOX.
C
         IF ( GI .LT. ZERO ) THEN
            GI = - MIN( ABS( BU( I ) - X( I ) ), - GI )
         ELSE
            GI =   MIN( ABS( BL( I ) - X( I ) ),   GI )
         END IF
C
C  RECORD THE NONZERO COMPONENTS OF THE CAUCHY DIRECTION IN GRAD.
C
         IF ( ABS( GI ) .GT. EPSMCH ) THEN
            NVAR         = NVAR + 1
            PJGNRM       = MAX( PJGNRM, ABS( GI ) )
            IVAR( NVAR ) = I
            GRAD( NVAR ) = GI
         END IF
   20 CONTINUE
      RETURN
      END
C  THIS VERSION: 13/01/1994 AT 04:41:23 PM.
CS    SUBROUTINE SCLCFG( UNSUCC, N     , NCALCF, NCALCG,
CD    SUBROUTINE DCLCFG( UNSUCC, N     , NCALCF, NCALCG,
     *                   ISTAEV, LSTAEV, ISTADG, LSTADG, IELING,
     *                   LELING, ICALCF, LCALCF, ICALCG, LCALCG,
     *                   ISPTRS, LNPTRS, ISELTS, LNELTS,
     *                   ISTAJC, LNSTJC, IGCOLJ, LNGCLJ, X, XNEW )
C
C  *********************************************************************
C
C  THIS ROUTINE SELECTS THE ELEMENTS AND GROUPS WHOSE VALUE MUST BE
C  RECOMPUTED, DUE TO THE CHANGE IN THE VARIABLE VECTOR FROM X TO XNEW.
C  IT IS ASSUMED THAT THE VALUES OF THE ELEMENTS AND GROUPS ARE
C  AVAILABLE AT X.  THESE ELEMENTS (GROUPS) ARE STORED IN THE VECTOR
C  ICALCF (ICALCG) FROM POSITION 1 TO NCALCF (NCALCG).
C
C  PH TOINT (WITH A FEW MODS BY NICK GOULD), SEPTEMBER 1990.
C  FOR CGT PRODUCTIONS.
C
C  ARGUMENTS
C
      INTEGER          N     , NCALCF, NCALCG
      INTEGER          LNPTRS, LNELTS, LSTAEV, LSTADG
      INTEGER          LNSTJC, LNGCLJ, LCALCF, LCALCG, LELING
      LOGICAL          UNSUCC
      INTEGER          ISTAEV( LSTAEV ), ISTADG( LSTADG )
      INTEGER          IELING( LELING )
      INTEGER          ICALCF( LCALCF ), ICALCG( LCALCG )
      INTEGER          ISPTRS( LNPTRS ), ISELTS( LNELTS )
      INTEGER          ISTAJC( LNSTJC ), IGCOLJ( LNGCLJ )
CS    REAL             X( N ), XNEW( N )
CD    DOUBLE PRECISION X( N ), XNEW( N )
C
C  INTERNAL VARIABLES
C
      INTEGER          I, K, IEL, IELL, IPT, IG
CS    REAL             DIFF, ONE, ZERO, XI
CD    DOUBLE PRECISION DIFF, ONE, ZERO, XI
      INTRINSIC        ABS
C
C  COMMON VARIABLES.
C
CS    REAL             EPSMCH, EPSNEG, TINY, BIG
CD    DOUBLE PRECISION EPSMCH, EPSNEG, TINY, BIG
CS    COMMON / SMACHN / EPSMCH, EPSNEG, TINY, BIG
CD    COMMON / DMACHN / EPSMCH, EPSNEG, TINY, BIG
C ** Correction 3. 13/01/94: 2 lines added
CS    SAVE   / SMACHN /
CD    SAVE   / DMACHN /
C ** Correction 3. 13/01/94: end of correction **
C
C  SET CONSTANT REAL PARAMETERS.
C
CS    PARAMETER ( ZERO   = 0.0E+0, ONE    = 1.0E+0 )
CD    PARAMETER ( ZERO   = 0.0D+0, ONE    = 1.0D+0 )
C
C  INITIALIZE THE NUMBER OF ELEMENTS AND GROUPS TO BE RECOMPUTED
C
      IF ( .NOT. UNSUCC ) THEN
         NCALCF = 0
         NCALCG = 0
      ELSE
C
C  RESET THE ELEMENT POINTERS TO THEIR CORRECT SIGNS
C
         DO 10 I = 1, NCALCF
            IEL = ICALCF( I )
            ISTAEV( IEL ) = - ISTAEV( IEL )
  10     CONTINUE
C
C  RESET THE GROUP POINTERS TO THEIR CORRECT SIGNS
C
         DO 20 I = 1, NCALCG
            IG = ICALCG( I )
            ISTADG( IG ) = - ISTADG( IG )
  20     CONTINUE
      END IF
C
C  DETECT THE VARIABLES THAT HAVE CHANGED SIGNIFICANTLY FROM X
C
      DO 100 I = 1, N
         XI   = XNEW( I )
         DIFF = ABS( XI - X( I ) )
         IF ( XI .NE. ZERO ) THEN
            IF ( DIFF .LT. EPSMCH * ABS( XI ) ) GO TO 100
         ELSE
            IF ( DIFF .LT. TINY ) GO TO 100
         END IF
C
C  THE I-TH VARIABLE HAS BEEN MODIFIED: FLAG ALL ELEMENTS WHERE IT
C  APPEARS
C
         IPT = ISPTRS( I )
         IF ( IPT .GE. 0) THEN
            IELL = ISELTS( I )
 200        CONTINUE
            IEL = IELING( IELL )
            IF ( ISTAEV( IEL ) .GT. 0 ) THEN
               NCALCF = NCALCF + 1
               ICALCF( NCALCF ) = IEL
               ISTAEV( IEL ) = - ISTAEV( IEL )
            END IF
            IF ( IPT .GT. 0 ) THEN
               IELL = ISELTS( IPT )
               IPT = ISPTRS( IPT )
               GO TO 200
            END IF
         END IF
C
C  FLAG ALL GROUPS THAT CONTAIN THE I-TH VARIABLE
C
         DO 300 K = ISTAJC( I ), ISTAJC( I + 1 ) - 1
            IG = IGCOLJ( K )
            IF ( ISTADG( IG ) .GT. 0 ) THEN
               NCALCG = NCALCG + 1
               ICALCG( NCALCG ) = IG
               ISTADG( IG ) = - ISTADG( IG )
            END IF
 300     CONTINUE
C
C  END OF THE LOOP ON THE VARIABLES
C
 100  CONTINUE
C
C  RESET THE ELEMENT POINTERS TO THEIR CORRECT SIGNS
C
      DO 400 I = 1, NCALCF
         IEL = ICALCF( I )
         ISTAEV( IEL ) = - ISTAEV( IEL )
 400  CONTINUE
C
C  RESET THE GROUP POINTERS TO THEIR CORRECT SIGNS
C
      DO 500 I = 1, NCALCG
         IG = ICALCG( I )
         ISTADG( IG ) = - ISTADG( IG )
 500  CONTINUE
      RETURN
C
C  END OF SUBROUTINE CLCFG.
C
      END
C  THIS VERSION: 13/01/1994 AT 04:41:23 PM.
      REAL FUNCTION RANDOM( I )
C
C  GIVES A RANDOM REAL IN (0, 1) (IF I > 0) OR (-1, 1) (IF I < 0).
C
      INTEGER I
      DOUBLE PRECISION GL, GR, R, S, BIGINT, ONE, TWOBIG, RLARGE
      INTRINSIC AINT, MOD
      PARAMETER( BIGINT = 3.2768D+4  , ONE    = 1.0D+0 )
      PARAMETER( RLARGE = 9.228907D+6, TWOBIG = 2.0D+0 * BIGINT )
      SAVE GL, GR
      DATA GL / 2.1845D+4 /
      DATA GR / 2.1845D+4 /
      R  = GR * RLARGE / TWOBIG
      S  = AINT( R )
      GL = MOD( S + GL * RLARGE, TWOBIG )
      GR = R - S
      IF ( I .GE. 0 ) THEN
         RANDOM = ( GL + GR ) / TWOBIG
      ELSE
         RANDOM = ( GL + GR ) / BIGINT - ONE
      END IF
      GR = GR * TWOBIG
      RETURN
C
C  END OF RANDOM.
C
      END
C  THIS VERSION: 13/01/1994 AT 04:41:23 PM.
      CHARACTER * 7 FUNCTION CHRTIM( TIME )
      REAL TIME
C
C  Obtain a 7 character representation of the time TIME.
C
C  Nick Gould, 13th September 1991.
C  For CGT Productions.
C
      INTEGER ITIM
      REAL TIM, TIMM, TIMH, TIMD
      CHARACTER * 7 CTIM
      CHRTIM( 1: 7 ) = '       '
      TIM = TIME
      TIMM = TIME / 6.0D+1
      TIMH = TIMM / 6.0D+1
      TIMD = TIMH / 2.4D+1
      IF ( TIM .LE. 9999.9 ) THEN
         TIM = TIME
         WRITE( UNIT = CTIM, FMT = 2000 ) TIM
         CHRTIM = CTIM
      ELSE IF ( TIM .LE. 99999.9 ) THEN
         TIM = TIME
         WRITE( UNIT = CTIM, FMT = 2000 ) TIM
         CHRTIM( 1: 1 ) = ' '
         CHRTIM( 2: 7 ) = CTIM( 1: 6 )
      ELSE IF ( TIM .LE. 9.99999E+5 ) THEN
         ITIM = INT( TIME )
         WRITE( UNIT = CTIM, FMT = 2010 ) ITIM
         CHRTIM = CTIM
      ELSE IF ( TIMM .LE. 9.99999E+4 ) THEN
         ITIM = INT( TIMM )
         WRITE( UNIT = CTIM( 1: 6 ), FMT = 2020 ) ITIM
         CHRTIM = CTIM( 1: 6 ) // 'm'
      ELSE IF ( TIMH .LE. 9.99999E+4 ) THEN
         ITIM = INT( TIMH )
         WRITE( UNIT = CTIM( 1: 6 ), FMT = 2020 ) ITIM
         CHRTIM = CTIM( 1: 6 ) // 'h'
      ELSE IF ( TIMD .LE. 9.99999E+4 ) THEN
         ITIM = INT( TIMD )
         WRITE( UNIT = CTIM( 1: 6 ), FMT = 2020 ) ITIM
         CHRTIM = CTIM( 1: 6 ) // 'd'
      ELSE
         CHRTIM = ' ******'
      END IF
 2000 FORMAT( 0P, F7.1 )
 2010 FORMAT( I7 )
 2020 FORMAT( I6 )
      RETURN
      END
C ** Correction 1. 06/08/93: Routines *DNRM, *SETVL, *SETVI appended
C                            from linpac.f 
CS    REAL             FUNCTION SDNRM( N, X, Y, TWONRM, SCALE, SCALED )
CD    DOUBLE PRECISION FUNCTION DDNRM( N, X, Y, TWONRM, SCALE, SCALED )
      INTEGER N
      LOGICAL TWONRM, SCALED
CS    REAL             X( N ), Y( N ), SCALE( * )
CD    DOUBLE PRECISION X( N ), Y( N ), SCALE( * )
C
C  COMPUTE THE SCALED (OR UNSCALED) TWO (OR INFINITY) NORM DISTANCE
C  BETWEEN THE VECTORS X AND Y.
C
C  NICK GOULD, JULY 10TH 1990.
C  FOR CGT PRODUCTIONS.
C
C  LOCAL VARIABLES.
C
      INTEGER          I
CS    REAL             ZERO
CD    DOUBLE PRECISION ZERO
      INTRINSIC        MAX, SQRT, ABS
C
C  SET CONSTANT REAL PARAMETERS.
C
CS    PARAMETER ( ZERO   = 0.0E+0 )
CD    PARAMETER ( ZERO   = 0.0D+0 )
CS    SDNRM = ZERO
CD    DDNRM = ZERO
      IF ( SCALED ) THEN
C
C  COMPUTE THE SCALED TWO-NORM DISTANCE BETWEEN X AND Y.
C
         IF ( TWONRM ) THEN
            DO 10 I = 1, N
CS             SDNRM  = SDNRM + ( ( X( I ) - Y( I ) ) / SCALE( I ) ) **2
CD             DDNRM  = DDNRM + ( ( X( I ) - Y( I ) ) / SCALE( I ) ) **2
   10       CONTINUE
CS          SDNRM = SQRT( SDNRM )
CD          DDNRM = SQRT( DDNRM )
C
C  COMPUTE THE SCALED INFINITY-NORM DISTANCE BETWEEN X AND Y.
C
         ELSE
            DO 20 I = 1, N
CS             SDNRM  = MAX( SDNRM,
CD             DDNRM  = MAX( DDNRM,
     *                   ABS( ( X( I ) - Y( I ) ) / SCALE( I ) ) )
   20       CONTINUE
         END IF
      ELSE 
C
C  COMPUTE THE TWO-NORM DISTANCE BETWEEN X AND Y.
C
         IF ( TWONRM ) THEN
            DO 30 I = 1, N
CS             SDNRM = SDNRM + ( X( I ) - Y( I ) ) ** 2
CD             DDNRM = DDNRM + ( X( I ) - Y( I ) ) ** 2
   30       CONTINUE
CS          SDNRM = SQRT( SDNRM )
CD          DDNRM = SQRT( DDNRM )
C
C  COMPUTE THE INFINITY-NORM DISTANCE BETWEEN X AND Y.
C
         ELSE
            DO 40 I = 1, N
CS             SDNRM = MAX( SDNRM, ABS( X( I ) - Y( I ) ) )
CD             DDNRM = MAX( DDNRM, ABS( X( I ) - Y( I ) ) )
   40       CONTINUE
         END IF
      END IF
      RETURN
C
C  END OF FUNCTION DNRM.
C
      END
CS    SUBROUTINE SSETVL( N, X, INCX, VL )
CD    SUBROUTINE DSETVL( N, X, INCX, VL )
C
C     ******************************************************************
C
      INTEGER          INCX, N
CS    REAL             VL, X( * )
CD    DOUBLE PRECISION VL, X( * )
      INTEGER          IX, M, I, J
      IF ( N .LE. 0 ) RETURN
      IF ( INCX .NE. 1 ) THEN
         IF( INCX .LT. 0 ) THEN
             IX = ( - N + 1 ) * INCX + 1
         ELSE
             IX = 1
         END IF
         J         = IX + ( N - 1 ) * INCX
         DO 100 I  = IX, J, INCX
            X( I ) = VL
  100    CONTINUE
      ELSE
         M = MOD( N, 5 )
         IF ( M .NE. 0 ) THEN
            DO 200 I  = 1, M
               X( I ) = VL
  200       CONTINUE
         END IF
         IF ( N .GE. 5 ) THEN
            IX            = M + 1
            DO 300 I      = IX, N, 5
               X( I )     = VL
               X( I + 1 ) = VL
               X( I + 2 ) = VL
               X( I + 3 ) = VL
               X( I + 4 ) = VL
  300       CONTINUE
         END IF
      END IF
      RETURN
      END
CS    SUBROUTINE SSETVI( NVAR, X, IVAR, VL )
CD    SUBROUTINE DSETVI( NVAR, X, IVAR, VL )
C
C     ******************************************************************
C
      INTEGER          NVAR, IVAR( * ), I
CS    REAL             VL, X( * )
CD    DOUBLE PRECISION VL, X( * )
      IF ( NVAR .LE. 0 ) RETURN
      DO 10 I           = 1, NVAR
         X( IVAR( I ) ) = VL
   10 CONTINUE
      RETURN
      END
C ** Correction 1. 06/08/93: end of correction **
