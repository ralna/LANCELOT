C  THIS VERSION: 18/07/1994 AT 10:51:57 AM.
C ** Correction report.
C ** Correction 1. 18/07/94: whole routine replaced **
C ** End of Correction report.
CS    SUBROUTINE SBARIA( N,  NG, NEL,    IELING, LELING, ISTADG, LSTADG,
CD    SUBROUTINE DBARIA( N,  NG, NEL,    IELING, LELING, ISTADG, LSTADG,
     *                   IELVAR, LELVAR, ISTAEV, LSTAEV, INTVAR, LNTVAR,
     *                   ISTADH, LSTADH, ICNA  , LICNA , ISTADA, LSTADA,
     *                   A , LA, B , LB, BL    , LBL   , BU    , LBU   ,
     *                   GSCALE, LGSCAL, ESCALE, LESCAL, VSCALE, LVSCAL,
     *                   GXEQX,  LGXEQX, INTREP, LINTRE, KNDOFC, LKNDOF,
     *                   RANGES, INFORM, FOBJ  , X , LX, U , LU, GVALS , 
     *                   LGVALS, FT    , LFT,    FUVALS, LFUVAL, XT    , 
     *                   LXT   , ICALCF, LCALCF, NCALCF, ICALCG, LCALCG,
     *                   NCALCG, IVAR  , LIVAR , NVAR  , Q , LQ, DGRAD , 
     *                   LDGRAD, ICHOSE, ITER  , MAXIT , QUADRT, VNAMES,
     *                   LVNAME, GNAMES, LGNAME, STOPG , STOPC , IWK   , 
     *                   LIWK  , WK    , LWK   , IPRINT, IOUT  )
C
      INTEGER          N, NEL, MAXIT , INFORM, ITER  , NVAR  , LU
      INTEGER          LELVAR, NG    , LSTAEV, LSTADH, LNTVAR, LCALCF
      INTEGER          LELING, LINTRE, LFT   , LGXEQX, LSTADG, LGVALS
      INTEGER          LGSCAL, LESCAL, LVSCAL, LCALCG, NCALCG, IPRINT
      INTEGER          LIVAR , LX    , LBL   , LBU   , LQ    , LDGRAD
      INTEGER          LGNAME, LXT   , LVNAME, LFUVAL, LWK   , LIWK
      INTEGER          LA, LB, LICNA , LSTADA, LKNDOF, NCALCF, IOUT
CS    REAL             STOPG , STOPC , FOBJ
CD    DOUBLE PRECISION STOPG , STOPC , FOBJ
      LOGICAL          QUADRT
      INTEGER          IELVAR( LELVAR       ), ISTAEV( LSTAEV         )
      INTEGER          ISTADH( LSTADH       ), IWK   ( LIWK           )
      INTEGER          INTVAR( LNTVAR       ), ISTADG( LSTADG         )
      INTEGER          ICNA  ( LICNA        ), ISTADA( LSTADA         )
      INTEGER          ICALCG( LCALCG       ), ICHOSE( 6              )
      INTEGER          KNDOFC( LKNDOF       ), ICALCF( LCALCF         )
      INTEGER          IVAR  ( LIVAR        ), IELING( LELING         )
      LOGICAL          INTREP( LINTRE       ), GXEQX ( LGXEQX         )
CS    REAL             X     ( LX           ), FUVALS( LFUVAL         ), 
CD    DOUBLE PRECISION X     ( LX           ), FUVALS( LFUVAL         ), 
     *                 BL    ( LBL          ), BU    ( LBU            ),
     *                 A     ( LA           ), B     ( LB             ), 
     *                 Q     ( LQ           ), DGRAD ( LDGRAD         ),
     *                 XT    ( LXT          ), FT    ( LFT            ),
     *                 GVALS ( LGVALS , 3   ), U     ( LU             ),
     *                 GSCALE( LGSCAL       ), ESCALE( LESCAL         ),
     *                 VSCALE( LVSCAL       ), WK    ( LWK            )
      CHARACTER * 10   VNAMES( LVNAME       ), GNAMES( LGNAME         )      
      EXTERNAL         RANGES
C
C  DUMMY SUBROUTINE AVAILABLE WITH LANCELOT.
C
C  NICK GOULD, 6TH DECEMBER 1990.
C  FOR CGT PRODUCTIONS.
C
      WRITE( 6, 2000 )
      STOP
C
C  NON-EXECUTABLE STATEMENTS.
C
 2000 FORMAT( /, ' We regret that the barrier function option',
     *           ' that you have ',
     *        /, ' chosen is not available with LANCELOT at',
     *           ' this time. ', 
     *        //, ' *** EXECUTION TERMINATING *** ', / )
      END
  
