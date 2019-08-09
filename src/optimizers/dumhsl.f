C ** Correction report.
C ** Correction 1. 10/08/93: 30 MC19 lines replaced by 29 MC29 lines **
C ** End of Correction report.
C  THIS VERSION: 11/08/1993 AT 08:58:48 AM.
C ** Correction 1. 10/08/93: 30 MC19 lines replaced by 29 MC29 lines **
CS    SUBROUTINE MC29A ( M, N, NE, A, IRN, ICN, R, C, W, LP, IFAIL )
CD    SUBROUTINE MC29AD( M, N, NE, A, IRN, ICN, R, C, W, LP, IFAIL )
      INTEGER          M, N, NE, LP, IFAIL
      INTEGER          IRN( * ), ICN( * )
CS    REAL             A( * ), R( * ), C( * ), W( * )
CD    DOUBLE PRECISION A( * ), R( * ), C( * ), W( * )
C
C  DUMMY SUBROUTINE AVAILABLE WITH LANCELOT.
C
C  NICK GOULD, 10TH AUGUST 1993.
C  FOR CGT PRODUCTIONS.
C
      WRITE( 6, 2000 )
      STOP
C
C  NON-EXECUTABLE STATEMENTS.
C
 2000 FORMAT( /, ' We regret that the solution options that you have ',
     *        /, ' chosen are not all freely available with LANCELOT. ', 
     *        //, ' If you have the Harwell Subroutine Library, this ',
     *        /, ' option may be enabled by replacing the dummy ', /,
     *        ' subroutine MC29A (single precision) or MC29AD', /,
     *        ' (double precision) with its H.S.L. namesake ', /,
     *        ' and dependencies.', //,
     *        ' *** EXECUTION TERMINATING *** ', / )
C
C  END OF DUMMY SUBROUTINE.
C
      END
C ** Correction 1. 10/08/93: end of correction **
CS    SUBROUTINE MA27A ( N, NZ, IRN, ICN, IW, LIW, IKEEP, IW1, NSTEPS,
CD    SUBROUTINE MA27AD( N, NZ, IRN, ICN, IW, LIW, IKEEP, IW1, NSTEPS,
     *                   IFLAG )
      INTEGER N, NZ, LIW, NSTEPS, IFLAG
      INTEGER IW1( N, 2 )
      INTEGER IRN( * ), ICN( * ), IW( LIW ), IKEEP( N, 3 )
C
C  DUMMY SUBROUTINE AVAILABLE WITH LANCELOT.
C
C  NICK GOULD, 7TH NOVEMBER 1990.
C  FOR CGT PRODUCTIONS.
C
      WRITE( 6, 2000 )
      STOP
C
C  NON-EXECUTABLE STATEMENTS.
C
 2000 FORMAT( /, ' We regret that the solution options that you have ',
     *        /, ' chosen are not all freely available with LANCELOT. ', 
     *        //, ' If you have the Harwell Subroutine Library, this ',
     *        /, ' option may be enabled by replacing the dummy ', /,
     *        ' subroutine MA27A (single precision) or MA27AD', /,
     *        ' (double precision) with its H.S.L. namesake ', /,
     *        ' and dependencies.', //,
     *        ' *** EXECUTION TERMINATING *** ', / )
C
C  END OF DUMMY SUBROUTINE.
C
      END
CS    SUBROUTINE MA27B ( N, NZ, IRN, ICN, A, LA, IW, LIW, IKEEP, NSTEPS,
CD    SUBROUTINE MA27BD( N, NZ, IRN, ICN, A, LA, IW, LIW, IKEEP, NSTEPS,
     *                   MAXFRT, IW1, IFLAG )
      INTEGER            N, NZ, LA, NSTEPS, MAXFRT, IFLAG, LIW
      INTEGER            IRN( * ), ICN( * ), IW( LIW )
      INTEGER            IKEEP( N, 3 ), IW1( N )
CS    REAL               A( LA )
CD    DOUBLE PRECISION   A( LA )
C
C  DUMMY SUBROUTINE AVAILABLE WITH LANCELOT.
C
C  NICK GOULD, 7TH NOVEMBER 1990.
C  FOR CGT PRODUCTIONS.
C
      WRITE( 6, 2000 )
      STOP
C
C  NON-EXECUTABLE STATEMENTS.
C
 2000 FORMAT( /, ' We regret that the solution options that you have ',
     *        /, ' chosen are not all freely available with LANCELOT. ', 
     *        //, ' If you have the Harwell Subroutine Library, this ',
     *        /, ' option may be enabled by replacing the dummy ', /,
     *        ' subroutine MA27B (single precision) or MA27BD', /,
     *        ' (double precision) with its H.S.L. namesake ', /,
     *        ' and dependencies.', //,
     *        ' *** EXECUTION TERMINATING *** ', / )
C
C  END OF DUMMY SUBROUTINE.
C
      END
CS    SUBROUTINE MA27C ( N, A, LA, IW, LIW, W, MAXFRT, RHS, IW1, 
CD    SUBROUTINE MA27CD( N, A, LA, IW, LIW, W, MAXFRT, RHS, IW1, 
     *                   NSTEPS )
      INTEGER            N, LA, NSTEPS, MAXFRT, LIW
      INTEGER            IW1( NSTEPS ), IW( LIW )
CS    REAL               A( LA ), W( MAXFRT ), RHS( N )
CD    DOUBLE PRECISION   A( LA ), W( MAXFRT ), RHS( N )
C
C  DUMMY SUBROUTINE AVAILABLE WITH LANCELOT.
C
C  NICK GOULD, 7TH NOVEMBER 1990.
C  FOR CGT PRODUCTIONS.
C
      WRITE( 6, 2000 )
      STOP
C
C  NON-EXECUTABLE STATEMENTS.
C
 2000 FORMAT( /, ' We regret that the solution options that you have ',
     *        /, ' chosen are not all freely available with LANCELOT. ', 
     *        //, ' If you have the Harwell Subroutine Library, this ',
     *        /, ' option may be enabled by replacing the dummy ', /,
     *        ' subroutine MA27C (single precision) or MA27CD', /,
     *        ' (double precision) with its H.S.L. namesake ', /,
     *        ' and dependencies.', //,
     *        ' *** EXECUTION TERMINATING *** ', / )
C
C  END OF DUMMY SUBROUTINE.
C
      END
CS    SUBROUTINE MA31D ( A, IRN, IA, N, IK, IP, ROW )
CD    SUBROUTINE MA31DD( A, IRN, IA, N, IK, IP, ROW )
      INTEGER           IA, N
      LOGICAL           ROW
      INTEGER           IK( N ), IP( N ), IRN( IA )
CS    REAL              A( IA )
CD    DOUBLE PRECISION  A( IA )
C
C  DUMMY SUBROUTINE AVAILABLE WITH LANCELOT.
C
C  NICK GOULD, 7TH NOVEMBER 1990.
C  FOR CGT PRODUCTIONS.
C
      WRITE( 6, 2000 )
      STOP
C
C  NON-EXECUTABLE STATEMENTS.
C
 2000 FORMAT( /, ' We regret that the solution options that you have ',
     *        /, ' chosen are not all freely available with LANCELOT. ', 
     *        //, ' If you have the Harwell Subroutine Library, this ',
     *        /, ' option may be enabled by replacing the dummy ', /,
     *        ' subroutine MA31D (single precision) or MA31DD', /,
     *        ' (double precision) with its H.S.L. namesake ', /,
     *        ' and dependencies.', //,
     *        ' *** EXECUTION TERMINATING *** ', / )
C
C  END OF DUMMY SUBROUTINE.
C
      END
