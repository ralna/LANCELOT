C  THIS VERSION: 11/08/1993 AT 08:58:48 AM.
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
