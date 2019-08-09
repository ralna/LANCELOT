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
