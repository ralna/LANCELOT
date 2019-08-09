C  THIS VERSION: 20/02/2001 AT 15:00:00 PM.
C ** Correction report.
C ** Correction -1. 03/03/00: Integer formats increased
C ** Correction 1. 29/01/93: 1 line corrected **
C ** Correction 2. 29/01/93: 4 lines replaced by 7 **
C ** Correction 3. 29/01/93: 1 line deleted **
C ** Correction 4. 28/07/93: 20 lines replaced by 49 **
C ** Correction 5. 28/07/93: 3 lines added **
C ** Correction 6. 28/07/93: 1 line corrected **
C ** Correction 7. 28/07/93: 3 lines added **
C ** Correction 8. 10/08/93: 1 line corrected **
C ** Correction 9. 13/01/94: 1 line added **
C ** Correction 10. 06/08/97: 4 lines replaced by 8 **
C ** Correction 11. 13/12/97: 2 lines corrected **
C ** Correction 12. 16/01/98: 9 lines added **
C ** Correction 13. 22/12/99: 1 line added **
C ** Correction 14. 01/05/2000: 1 line corrected **
C ** Correction 15. 01/05/2000: 2 lines added **
C ** Correction 16. 01/05/2000: 1 line corrected **
C ** Correction 17. 20/02/2000: 2 lines updated **
C ** End of Correction report.
CS    SUBROUTINE SWRKSP( INPUT , IOUT  , IINRUN, IOUTSL, IOUT1L, IINPB ,
CD    SUBROUTINE DWRKSP( INPUT , IOUT  , IINRUN, IOUTSL, IOUT1L, IINPB ,
     *                   IOUTPB, IALIVE, DEBUG , DEBUGF, DEBUGG )
      INTEGER            INPUT , IOUT  , IINRUN, IOUTSL, IOUT1L
      INTEGER            IINPB , IOUTPB, IALIVE
      LOGICAL            DEBUG , DEBUGF, DEBUGG
C
C  PARTITION THE ARRAYS FOR THE MAIN PROGRAM FOR LANCELOT
C  SUITE OF PROGRAMS AND THEN CALL THE MAIN PROGRAM.
C
C  NICK GOULD, FOR CGT PRODUCTIONS.
C  DECEMBER 4TH 1990.
C
      INTEGER          LIWK, LWK , LFUVAL, LLOGIC, LCHARA
C
C----------------------------------------------------------------------------
C
C          PARAMETERS WHOSE VALUE MIGHT BE CHANGED BY THE USERS
C          ----------------------------------------------------
C
C  THE FOLLOWING PARAMETERS DEFINE THE SIZES OF VARIOUS PROBLEM
C  DEPENDENT ARRAYS.
C
C  THESE MAY BE CHANGED BY THE USE TO SUIT A PARTICULAR PROBLEM OR
C  SYSTEM CONFIGURATION.
C
C  HOWEVER, THEY SHOULD BE CHANGED WITH CARE.  THE USERS ARE URGED TO
C  READ THE LANCELOT MANUAL AND UNDERSTAND WHAT THEY ARE DOING!
C
C  THE PACKAGE WILL ISSUE ERROR MESSAGES IF ANY OF THESE SIZES IS TOO SMALL,
C  TELLING WHICH PARAMETER TO INCREASE.
C
C----------------------------------------------------------------------------
C
C  INTEGER WORKSPACE
C
CHUG  PARAMETER      ( LIWK   =  30000000 )
CBIG  PARAMETER      ( LIWK   =  3000000  )
CMED  PARAMETER      ( LIWK   =  300000   )
CTOY  PARAMETER      ( LIWK   =  10000    )
C
C  REAL/DOUBLE PRECISION WORKSPACE
C
CHUG  PARAMETER      ( LWK    =  30000000 )
CBIG  PARAMETER      ( LWK    =  3000000  )
CMED  PARAMETER      ( LWK    =  300000   )
CTOY  PARAMETER      ( LWK    =  10000    )
C
C  LOGICAL WORKSPACE
C
CHUG  PARAMETER      ( LLOGIC =  1500000 )
CBIG  PARAMETER      ( LLOGIC =  150000  )
CMED  PARAMETER      ( LLOGIC =  15000   )
CTOY  PARAMETER      ( LLOGIC =  1500    )
C
C  CHARACTER WORKSPACE
C
CHUG  PARAMETER      ( LCHARA =  3000000 )
CBIG  PARAMETER      ( LCHARA =  300000  )
CMED  PARAMETER      ( LCHARA =  30000   )
CTOY  PARAMETER      ( LCHARA =  3000    )
C
C  WORKSPACE TO STORE THE PROBLEM'S FUNCTION AND DERIVATIVES VALUES
C  (PROBLEM DEPENDENT)
C
CHUG  PARAMETER      ( LFUVAL =  5000000 )
CBIG  PARAMETER      ( LFUVAL =  500000  )
CMED  PARAMETER      ( LFUVAL =  50000   )
CTOY  PARAMETER      ( LFUVAL =  5000    )
C
C-----------------------------------------------------------------------------
C
C  END OF PARAMETERS WHICH MIGHT BE CHANGED BY THE USERS
C
C-----------------------------------------------------------------------------
C
      INTEGER          IWK( LIWK )
CS    REAL             WK ( LWK  ), FUVALS( LFUVAL )
CD    DOUBLE PRECISION WK ( LWK  ), FUVALS( LFUVAL )
      LOGICAL          LOGI ( LLOGIC )
      CHARACTER * 10   CHA  ( LCHARA )
      INTEGER          N,  NG, NELNUM, NGEL  , NVARS , NNZA  , NGPVLU,
     *                 NEPVLU, NG1   , NEL1  , ISTADG, ISTGP , ISTADA,
     *                 ISTAEV, ISTEP , ITYPEG, KNDOFC, ITYPEE,
     *                 IELING, IELVAR, ICNA  , ISTADH, INTVAR, IVAR  ,
     *                 ICALCF, ICALCG, IWRK  , A     , B     , BL    ,
     *                 BU    , X     , U     , GPVALU, EPVALU,
     *                 ESCALE, GSCALE, VSCALE, GVALS , XT    , DGRAD ,
     *                 Q     , WRK   , INTREP, GXEQX , GNAMES, VNAMES,
     *                 LO    , CH    , LIWORK, LWORK , NGNG  , FT    ,
     *                 LSTADG, LSTGP , LSTADA, LSTAEV, LSTEP , LTYPEG,
     *                 LKNDOF, LTYPEE, LELING, LELVAR, LICNA , LSTADH,
     *                 LNTVAR, LIVAR , LCALCF, LCALCG, ETYPES, GTYPES,
     *                 LA, LB, LBL   , LBU   , LX, LU, LGPVLU, LETYPE,
     *                 LEPVLU, LESCAL, LGSCAL, LVSCAL, LGTYPE, LGVALS,
     *                 LXT   , LDGRAD, LQ    , LFT   , LINTRE, LGXEQX,
     *                 LGNAME, LVNAME, NELTYP, NGRTYP, IALGOR
      CHARACTER * 8    PNAME
C
C  INPUT THE PROBLEM DIMENSIONS.
C
      READ( INPUT, 1000 ) N,  NG, NELNUM, NGEL, NVARS, NNZA, NGPVLU,
     *                    NEPVLU, NELTYP, NGRTYP
      READ( INPUT, 1010 ) IALGOR, PNAME
      NG1    = NG     + 1
      NGNG   = NG     + NG
      NEL1   = NELNUM + 1
C
C  SET DIMENSIONS FOR THE PARTITIONED INTEGER WORKSPACE.
C
      LSTADG = MAX( 1, NG1 )
      LSTGP  = MAX( 1, NG1 )
      LSTADA = MAX( 1, NG1 )
      LSTAEV = MAX( 1, NEL1 )
      LSTEP  = MAX( 1, NEL1 )
      LTYPEG = MAX( 1, NG )
      LKNDOF = MAX( 1, NG )
      LTYPEE = MAX( 1, NELNUM )
      LELING = MAX( 1, NGEL )
      LELVAR = MAX( 1, NVARS )
      LICNA  = MAX( 1, NNZA )
      LSTADH = MAX( 1, NEL1 )
      LNTVAR = MAX( 1, NEL1 )
      LIVAR  = MAX( 1, N )
      LCALCF = MAX( 1, NELNUM, NG )
      LCALCG = MAX( 1, NG )
C
C  PARTITION THE INTEGER WORKSPACE.
C
      ISTADG = 1
      ISTGP  = ISTADG + NG1
      ISTADA = ISTGP  + NG1
      ISTAEV = ISTADA + NG1
      ISTEP  = ISTAEV + NEL1
      ITYPEG = ISTEP  + NEL1
      KNDOFC = ITYPEG + NG
      ITYPEE = KNDOFC + NG
      IELING = ITYPEE + NELNUM
      IELVAR = IELING + NGEL
      ICNA   = IELVAR + NVARS
      ISTADH = ICNA   + NNZA
      INTVAR = ISTADH + NEL1
      IVAR   = INTVAR + NEL1
      ICALCF = IVAR   + N
      ICALCG = ICALCF + LCALCF
      IWRK   = ICALCG + NG
      LIWORK = LIWK   - IWRK
C
C  ENSURE THERE IS SUFFICIENT ROOM.
C
      IF ( LIWORK .LT. 0 ) THEN
         WRITE( IOUT, 2000 ) 'IWK   ', 'LIWK  ', - LIWORK
         RETURN
      END IF
C
C  SET DIMENSIONS FOR THE PARTITIONED REAL WORKSPACE.
C
      LA     = MAX( 1, NNZA )
      LB     = MAX( 1, NG )
      IF ( IALGOR .LE. 2 ) THEN
         LBL = MAX( 1, N )
         LBU = MAX( 1, N )
      ELSE
         LBL = MAX( 1, N + NG )
         LBU = MAX( 1, N + NG )
      END IF
      LX     = MAX( 1, N )
      LU     = MAX( 1, NG )
      LGPVLU = MAX( 1, NGPVLU )
      LEPVLU = MAX( 1, NEPVLU )
      LESCAL = MAX( 1, NGEL )
      LGSCAL = MAX( 1, NG )
      LVSCAL = MAX( 1, N )
      LGVALS = MAX( 1, NG )
      LXT    = MAX( 1, N )
      LDGRAD = MAX( 1, N )
      LQ     = MAX( 1, N )
      LFT    = MAX( 1, NG )
C
C  PARTITION THE REAL WORKSPACE.
C
      A      = 1
      B      = A      + NNZA
      BL     = B      + NG
      IF ( IALGOR .LE. 2 ) THEN
         BU  = BL     + N
         X   = BU     + N
      ELSE
         BU  = BL     + N + NG
         X   = BU     + N + NG
      END IF
      U      = X      + N
      GPVALU = U      + NG
      EPVALU = GPVALU + NGPVLU
      ESCALE = EPVALU + NEPVLU
      GSCALE = ESCALE + NGEL
      VSCALE = GSCALE + NG
      GVALS  = VSCALE + N
      XT     = GVALS  + 3 * NG
      DGRAD  = XT     + N
      Q      = DGRAD  + N
      FT     = Q      + N
      WRK    = FT     + NG
      LWORK  = LWK    - WRK
C
C  ENSURE THERE IS SUFFICIENT ROOM.
C
      IF ( LWORK .LT. 0 ) THEN
         WRITE( IOUT, 2000 ) 'WK   ', 'LWK   ', - LWORK
         RETURN
      END IF
C
C  SET DIMENSIONS FOR THE PARTITIONED LOGICAL WORKSPACE.
C
      LINTRE = MAX( 1, NELNUM )
      LGXEQX = MAX( 1, NGNG )
C
C  PARTITION THE LOGICAL WORKSPACE.
C
      INTREP = 1
      GXEQX  = INTREP + NELNUM
      LO     = GXEQX  + NGNG
C
C  ENSURE THERE IS SUFFICIENT ROOM.
C
      IF ( LLOGIC .LT. LO ) THEN
         WRITE( IOUT, 2000 ) 'LOGI  ', 'LLOGIC', LO - LLOGIC
         RETURN
      END IF
C
C  SET DIMENSIONS FOR THE PARTITIONED CHARACTER WORKSPACE.
C
      LGNAME = MAX( 1, NG )
      LVNAME = MAX( 1, N )
      LETYPE = MAX( 1, NELTYP )
      LGTYPE = MAX( 1, NGRTYP )
C
C  PARTITION THE CHARACTER WORKSPACE.
C
      GNAMES = 1
      VNAMES = GNAMES + NG
      ETYPES = VNAMES + N
      GTYPES = ETYPES + NELTYP
      CH     = GTYPES + NGRTYP
C
C  ENSURE THERE IS SUFFICIENT ROOM.
C
      IF ( LCHARA .LT. CH ) THEN
         WRITE( IOUT, 2000 ) 'CHA   ', 'LCHARA', CH - LCHARA
         RETURN
      END IF
C
C  CALL THE MAIN PROGRAM.
C
CS    CALL SLANCE( N , NG, NELNUM, NG1, NEL1, NGNG, NGEL, NVARS, NNZA,
CD    CALL DLANCE( N , NG, NELNUM, NG1, NEL1, NGNG, NGEL, NVARS, NNZA,
     *             NGPVLU,NEPVLU,NELTYP,NGRTYP,LSTADG,LSTGP ,LSTADA,
     *             LSTAEV,LSTEP ,LTYPEG,LKNDOF,LTYPEE,LELING,LELVAR,
     *             LICNA ,LSTADH,LNTVAR,LIVAR ,LCALCF,LCALCG,LIWORK, LA, 
     *             LB, LBL, LBU, LX, LU, LGPVLU, LEPVLU, LESCAL, LGSCAL, 
     *             LVSCAL, LGVALS, LXT, LDGRAD, LQ, LFT, LWORK , LFUVAL,
     *             LINTRE, LGXEQX, LGNAME, LVNAME, LETYPE, LGTYPE,
     *             IWK( ISTADG ), IWK( ISTGP  ), IWK( ISTADA ),
     *             IWK( ISTAEV ), IWK( ISTEP  ), IWK( ITYPEG ),
     *             IWK( KNDOFC ), IWK( ITYPEE ), IWK( IELING ),
     *             IWK( IELVAR ), IWK( ICNA   ), IWK( ISTADH ),
     *             IWK( INTVAR ), IWK( IVAR   ), IWK( ICALCF ),
     *             IWK( ICALCG ), IWK( IWRK   ), WK ( A ), WK ( B ),
     *             WK( BL ), WK( BU ), WK ( X ), WK( U ), WK( GPVALU ),
     *             WK ( EPVALU ), WK( ESCALE ), WK( GSCALE ),
     *             WK ( VSCALE ), WK( GVALS  ), WK( XT ), WK ( DGRAD ),
     *             WK( Q ), WK( FT ), WK( WRK ), FUVALS, LOGI( INTREP ), 
     *             LOGI(GXEQX), CHA(GNAMES), CHA(VNAMES), CHA( ETYPES ), 
     *             CHA( GTYPES ),INPUT,IOUT,IINRUN,IOUTSL,IOUT1L,IINPB,
     *             IOUTPB, IALIVE, IALGOR,PNAME, DEBUG, DEBUGF, DEBUGG )
      RETURN
C
C  NON-EXECUTABLE STATEMENTS.
C
C ** Correction -1a. 03/03/00: Integer formats increased
 1000 FORMAT( 10I8 )
 1010 FORMAT( I2, A8 )
 2000 FORMAT( /, ' ** PROGRAM WRKSP: ARRAY LENGTH ', A6, ' TOO SMALL.',
     *        /, ' -- MINIMIZATION ABANDONED.',
     *        /, ' -- INCREASE THE PARAMETER ', A6, ' BY AT LEAST ', I8,
     *           ' AND RESTART.'  )
C
C  END OF WRKSP.
C
      END
C  THIS VERSION: 13/01/1994 AT 04:37:53 PM.
CS    SUBROUTINE SLANCE( N , NG, NELNUM, NG1   , NEL1  , NGNG  , NGEL  ,
CD    SUBROUTINE DLANCE( N , NG, NELNUM, NG1   , NEL1  , NGNG  , NGEL  ,
     *                   NVARS , NNZA  , NGPVLU, NEPVLU, NELTYP, NGRTYP,
     *                   LSTADG, LSTGP , LSTADA, LSTAEV, LSTEP , LTYPEG, 
     *                   LKNDOF, LTYPEE, LELING, LELVAR, LICNA , LSTADH, 
     *                   LNTVAR, LIVAR , LCALCF, LCALCG, LIWK  ,
     *                   LA, LB, LBL   , LBU   , LX, LU, LGPVLU,
     *                   LEPVLU, LESCAL, LGSCAL, LVSCAL, LGVALS,
     *                   LXT   , LDGRAD, LQ    , LFT   , LWK   , LFUVAL, 
     *                   LINTRE, LGXEQX, LGNAME, LVNAME, LETYPE, LGTYPE,
     *                   ISTADG, ISTGP , ISTADA, ISTAEV, ISTEP , ITYPEG,
     *                   KNDOFC, ITYPEE, IELING, IELVAR, ICNA  ,
     *                   ISTADH, INTVAR, IVAR  , ICALCF, ICALCG, IWK   ,
     *                   A , B , BL, BU, X , U , GPVALU, EPVALU,
     *                   ESCALE, GSCALE, VSCALE, GVALS , XT    , DGRAD ,
     *                   Q     , FT    , WK    , FUVALS, INTREP, GXEQX ,
     *                   GNAMES, VNAMES, ETYPES, GTYPES, INPUT , IOUT  ,
     *                   IINRUN, IOUTSL, IOUT1L, IINPB , IOUTPB, IALIVE,
     *                   IALGOR, PNAME , DEBUG , DEBUGF, DEBUGG )
C
C  MAIN PROGRAM FOR LANCELOT SUITE OF PROGRAMS
C
C  NICK GOULD, FOR CGT PRODUCTIONS.
C  DECEMBER 4TH 1990.
C
      INTEGER            N,  NG, NELNUM, NGEL  , NVARS , NNZA  , NGPVLU,
     *                   NEPVLU, NG1   , NEL1  , NGNG  , INPUT , IOUT  , 
     *                   IINRUN, IOUTSL, IOUT1L, IINPB , IOUTPB, IALIVE, 
     *                   LSTADG, LSTGP , LSTADA, LSTAEV, LSTEP , LTYPEG, 
     *                   LKNDOF, LTYPEE, LELING, LELVAR, LICNA , LSTADH, 
     *                   LNTVAR, LIVAR , LCALCF, LCALCG, LIWK  , NELTYP,
     *                   LA, LB, LBL   , LBU   , LX, LU, LGPVLU, NGRTYP,
     *                   LEPVLU, LESCAL, LGSCAL, LVSCAL, LGVALS, LXT   ,
     *                   LDGRAD, LQ    , LFT   , LWK   , LFUVAL, LINTRE,
     *                   LGXEQX, LGNAME, LVNAME, LETYPE, LGTYPE, IALGOR
      LOGICAL            DEBUG , DEBUGF, DEBUGG
      CHARACTER * 8      PNAME
      INTEGER            ISTADG( LSTADG ), ISTGP ( LSTGP  ),
     *                   ISTADA( LSTADA ), ISTAEV( LSTAEV ),
     *                   ISTEP ( LSTAEV ), ITYPEG( LTYPEG ),
     *                   KNDOFC( LKNDOF ), ITYPEE( LTYPEE ),
     *                   IELING( LELING ),
     *                   IELVAR( LELVAR ), ICNA  ( LICNA  ),
     *                   ISTADH( LSTADH ), INTVAR( LNTVAR ),
     *                   IVAR  ( LIVAR  ), ICALCF( LCALCF ),
     *                   ICALCG( LCALCG ), IWK   ( LIWK   )
CS    REAL               A     ( LA     ), B     ( LB     ),
CD    DOUBLE PRECISION   A     ( LA     ), B     ( LB     ),
     *                   BL    ( LBL    ), BU    ( LBU    ),
     *                   X     ( LX     ), U     ( LU     ),
     *                   GPVALU( LGPVLU ), EPVALU( LEPVLU ),
     *                   ESCALE( LESCAL ), GSCALE( LGSCAL ),
     *                   VSCALE( LVSCAL ), GVALS ( LGVALS  , 3 ),
     *                   XT    ( LXT    ), DGRAD ( LDGRAD ),
     *                   Q     ( LQ     ), FT    ( LFT    ),
     *                   WK    ( LWK    ), FUVALS( LFUVAL )
      LOGICAL            INTREP( LINTRE ), GXEQX ( LGXEQX )
      CHARACTER * 10     GNAMES( LGNAME ), VNAMES( LVNAME ),
     *                   ETYPES( LETYPE ), GTYPES( LGTYPE )
C
C  LOCAL VARIABLES.
C
      INTEGER            NOBJGR, ITEST , MAXIT , NELNM1, NCALCG
      INTEGER            ITER  , JUMPTO, LFTUV , IPRNTS, IPRINT
      INTEGER            LWRKD , LXEL  , LXINT , LXTT  , NCALCF
C ** Correction 6. 28/07/93: 1 line corrected **
      INTEGER            NIN   , NINMAX, NLMAX , LW, L , IR, IC
C ** Correction 6. 28/07/93: end of correction **
      INTEGER            NVAR  , LGTSTR, LGTVAL, LFTT  , IELTYP
      INTEGER            I , J , IPRNT , IPSTOP, IPSTRT, IFFLAG
      INTEGER            LSSTIV, LSTINV, LSTAGV, LSSVGR, LSVGRP, LGRJAC
      INTEGER            LSGRJA, NUMVAR, LGFX  , IGRTYP, IPGAP , ISTORE
      INTEGER            LRNGRJ, LCNGRJ, LVRSCA, NORDER, LGPSCA, LWKST
      INTEGER            NFREE , NFIXED, NLOWER, NUPPER, NBOTH , NSLACK
      INTEGER            NLINOB, NNLNOB, NLINEQ, NNLNEQ, NLININ, NNLNIN
      INTEGER            INFORM, ICHOSE( 6 )
      REAL               DUM   , TIME  , TIMM  , TTOTAL
CS    REAL               RMU   , STOPGA, STOPCA, RMUTOL, FOBJ
CD    DOUBLE PRECISION   RMU   , STOPGA, STOPCA, RMUTOL, FOBJ
CS    REAL               FIRSTC, FIRSTG, OBJFBN, OBFBND( 2 )
CD    DOUBLE PRECISION   FIRSTC, FIRSTG, OBJFBN, OBFBND( 2 )
C ** Correction 1. 29/01/93: 1 lines corrected **
      LOGICAL            TESTAL, DECHKE, DECHKG, WARNNG, SCALEG, DSAVE
C ** Correction 1. 29/01/93: end of correction **
      LOGICAL            QUADRT, FATALE, FATALG, SECOND, GETSCA, FDGRAD
      LOGICAL            SCALEV, NEWSOL, SQUARE, SCALED, WARMST, ALIVE
      CHARACTER * 3      MINMAX
      CHARACTER * 5      ST    , OPTIMI
      CHARACTER * 8      PNAME2
      CHARACTER * 10     TEMPNA
C
C  INTRINSICS, FUNCTIONS AND SUBROUTINES.
C
      INTRINSIC          DBLE, MIN, MAX, MOD
      REAL               RANDOM, CPUTIM
CS    REAL               SMACHR
CD    DOUBLE PRECISION   DMACHR
      EXTERNAL           RANDOM, CPUTIM, SETTYP, RANGES, ELFUNS, GROUPS
CS    EXTERNAL           SMACHR, SSPECI, SBARIA, SSBMIN, SAUGLG, SDRCHE
CD    EXTERNAL           DMACHR, DSPECI, DBARIA, DSBMIN, DAUGLG, DDRCHE
CS    EXTERNAL           SDRCHG, SSCALN
CD    EXTERNAL           DDRCHG, DSCALN
C
C  COMMON VARIABLES.
C
      COMMON / NUMEL  / NELNM1
CS    COMMON / SCOMAL / RMU   , RMUTOL, FIRSTC, FIRSTG, NEWSOL
CD    COMMON / DCOMAL / RMU   , RMUTOL, FIRSTC, FIRSTG, NEWSOL
      INTEGER           ITERCG, ITCGMX, NGEVAL, ISKIP , IFIXED, NBANDW
CS    REAL              ACCCG , RATIO , RADIUS, RADMAX, FINDMX, CMA31
CD    DOUBLE PRECISION  ACCCG , RATIO , RADIUS, RADMAX, FINDMX, CMA31
CS    COMMON / SCOMSB / ACCCG , RATIO , RADIUS, RADMAX, FINDMX, CMA31,
CD    COMMON / DCOMSB / ACCCG , RATIO , RADIUS, RADMAX, FINDMX, CMA31,
     *                  ITERCG, ITCGMX, NGEVAL, ISKIP , IFIXED, NBANDW
CS    REAL              EPSMCH, EPSNEG, TINY, BIG
CD    DOUBLE PRECISION  EPSMCH, EPSNEG, TINY, BIG
CS    COMMON / SMACHN / EPSMCH, EPSNEG, TINY, BIG
CD    COMMON / DMACHN / EPSMCH, EPSNEG, TINY, BIG
C ** Correction 9. 13/01/94: 1 lines added
      SAVE
C ** Correction 9. 13/01/94: end of correction **
      NELNM1 = NELNUM
C
C  SET UP INITIAL DATA.
C
CS    EPSMCH = SMACHR( 1 )
CD    EPSMCH = DMACHR( 1 )
      QUADRT = .FALSE.
      NEWSOL = .FALSE.
C
C  SET UP DATA FOR NEXT PROBLEM.
C
   10 CONTINUE
CS    CALL SSPECI( IINRUN, IPRNT , IPSTRT, IPSTOP, IPGAP , MAXIT ,
CD    CALL DSPECI( IINRUN, IPRNT , IPSTRT, IPSTOP, IPGAP , MAXIT ,
     *             NBANDW, RMU   , RMUTOL, GETSCA, DECHKE, DECHKG,
     *             FATALE, FATALG, ICHOSE, IPRNTS, SCALEG, SCALEV,
     *             FINDMX, STOPGA, STOPCA, FIRSTC, FIRSTG, TESTAL,
     *             WARMST, RADIUS, ISTORE, ITEST , IOUT )
      FDGRAD = ICHOSE( 3 ) .GE. 1
      IF ( ITEST .EQ. 0 ) GO TO 900
      ITER   = 0
      IPRINT = 0
      SECOND = ICHOSE( 4 ) .EQ. 0
C
C  PRINT OUT PROBLEM DATA. INPUT THE NUMBER OF VARIABLES, GROUPS,
C  ELEMENTS AND THE IDENTITY OF THE OBJECTIVE FUNCTION GROUP.
C
      IF ( IALGOR .EQ. 2 ) THEN
         READ( INPUT, 1002 ) NSLACK
      ELSE
         NSLACK = 0
      END IF
      NUMVAR = N - NSLACK
      IF ( IALGOR .EQ. 1 ) SCALEG = .FALSE.
      IF ( DEBUG ) WRITE( 6, 1100 ) PNAME, N, NG, NELNUM
C
C  INPUT THE STARTING ADDRESSES OF THE ELEMENTS IN EACH GROUP,
C  OF THE PARAMETERS USED FOR EACH GROUP AND
C  OF THE NONZEROS OF THE LINEAR ELEMENT IN EACH GROUP.
C
      READ( INPUT, 1010 ) ( ISTADG( I ), I = 1, NG1 )
      IF ( DEBUG ) WRITE( 6, 1110 ) 'ISTADG',
     *   ( ISTADG( I ), I = 1, NG1 )
      READ( INPUT, 1010 ) ( ISTGP ( I ), I = 1, NG1 )
      IF ( DEBUG ) WRITE( 6, 1110 ) 'ISTGP ',
     *   ( ISTGP ( I ), I = 1, NG1 )
      READ( INPUT, 1010 ) ( ISTADA( I ), I = 1, NG1 )
      IF ( DEBUG ) WRITE( 6, 1110 ) 'ISTADA',
     *   ( ISTADA( I ), I = 1, NG1 )
C
C  INPUT THE STARTING ADDRESSES OF THE VARIABLES AND PARAMETERS
C  IN EACH ELEMENT.
C
      READ( INPUT, 1010 ) ( ISTAEV( I ), I = 1, NEL1 )
      IF ( DEBUG ) WRITE( 6, 1110 ) 'ISTAEV',
     *   ( ISTAEV( I ), I = 1, NEL1 )
      READ( INPUT, 1010 ) ( ISTEP( I ), I = 1, NEL1 )
      IF ( DEBUG ) WRITE( 6, 1110 ) 'ISTEP ',
     *   ( ISTEP( I ), I = 1, NEL1 )
C
C  INPUT THE GROUP TYPE OF EACH GROUP
C
      READ( INPUT, 1010 ) ( ITYPEG( I ), I = 1, NG )
      IF ( DEBUG ) WRITE( 6, 1110 ) 'ITYPEG',
     *   ( ITYPEG( I ), I = 1, NG )
      IF ( IALGOR .GE. 2 ) THEN
         READ( INPUT, 1010 ) ( KNDOFC( I ), I = 1, NG )
         IF ( DEBUG ) WRITE( 6, 1110 ) 'KNDOFC',
     *      ( KNDOFC( I ), I = 1, NG )
      END IF
C
C  INPUT THE ELEMENT TYPE OF EACH ELEMENT
C
      READ( INPUT, 1010 ) ( ITYPEE( I ), I = 1, NELNUM )
      CALL SETTYP( NELNUM, ITYPEE )
      IF ( DEBUG ) WRITE( 6, 1110 ) 'ITYPEE',
     *   ( ITYPEE( I ), I = 1, NELNUM )
C
C  INPUT THE NUMBER OF INTERNAL VARIABLES FOR EACH ELEMENT.
C
      READ( INPUT, 1010 ) ( INTVAR( I ), I = 1, NELNUM )
      IF ( DEBUG ) WRITE( 6, 1110 ) 'INTVAR',
     *   ( INTVAR( I ), I = 1, NELNUM )
C
C  INPUT THE IDENTITY OF EACH INDIVIDUAL ELEMENT.
C
      READ( INPUT, 1010 ) ( IELING( I ), I = 1, NGEL )
      IF ( DEBUG ) WRITE( 6, 1110 ) 'IELING',
     *   ( IELING( I ), I = 1, NGEL )
C
C  INPUT THE VARIABLES IN EACH GROUP'S ELEMENTS.
C
      NVARS = ISTAEV( NEL1 ) - 1
      READ( INPUT, 1010 ) ( IELVAR( I ), I = 1, NVARS )
      IF ( DEBUG ) WRITE( 6, 1110 ) 'IELVAR',
     *   ( IELVAR( I ), I = 1, NVARS )
C
C  INPUT THE COLUMN ADDRESSES OF THE NONZEROS IN EACH LINEAR ELEMENT.
C
      READ( INPUT, 1010 ) ( ICNA( I ), I = 1, NNZA )
      IF ( DEBUG ) WRITE( 6, 1110 ) 'ICNA  ',
     *   ( ICNA( I ), I = 1, NNZA )
C
C  INPUT THE VALUES OF THE NONZEROS IN EACH LINEAR ELEMENT, THE
C  CONSTANT TERM IN EACH GROUP, THE LOWER AND UPPER BOUNDS ON
C  THE VARIABLES AND THE STARTING POINT FOR THE MINIMIZATION.
C
      READ( INPUT, 1020 ) ( A( I ), I = 1, NNZA )
      IF ( DEBUG ) WRITE( 6, 1120 ) 'A     ',
     *   ( A( I ), I = 1, NNZA )
      READ( INPUT, 1020 ) ( B( I ), I = 1, NG )
      IF ( DEBUG ) WRITE( 6, 1120 ) 'B     ',
     *   ( B( I ), I = 1, NG )
      IF ( IALGOR .LE. 2 ) THEN
         READ( INPUT, 1020 ) ( BL( I ), I = 1, N )
         IF ( DEBUG ) WRITE( 6, 1120 ) 'BL    ',
     *      ( BL( I ), I = 1, N )
         READ( INPUT, 1020 ) ( BU( I ), I = 1, N )
         IF ( DEBUG ) WRITE( 6, 1120 ) 'BU    ',
     *      ( BU( I ), I = 1, N )
      ELSE
         READ( INPUT, 1020 ) ( BL( I ), I = 1, N + NG )
         IF ( DEBUG ) WRITE( 6, 1120 ) 'BL    ',
     *      ( BL( I ), I = 1, N + NG )
         READ( INPUT, 1020 ) ( BU( I ), I = 1, N + NG )
         IF ( DEBUG ) WRITE( 6, 1120 ) 'BU    ',
     *      ( BU( I ), I = 1, N + NG )
      END IF
      READ( INPUT, 1020 ) ( X( I ), I = 1, N )
      IF ( DEBUG ) WRITE( 6, 1120 ) 'X     ',
     *   ( X( I ), I = 1, N )
      IF ( IALGOR .GE. 2 ) READ( INPUT, 1020 )( U( I ), I = 1, NG )
      IF ( DEBUG ) WRITE( 6, 1120 ) 'U     ',
     *   ( U( I ), I = 1, NG )
C
C  INPUT THE PARAMETERS IN EACH GROUP.
C
      READ( INPUT, 1020 ) ( GPVALU( I ), I = 1, NGPVLU )
      IF ( DEBUG ) WRITE( 6, 1120 ) 'GPVALU',
     *   ( GPVALU( I ), I = 1, NGPVLU )
C
C  INPUT THE PARAMETERS IN EACH INDIVIDUAL ELEMENT.
C
      READ( INPUT, 1020 ) ( EPVALU( I ), I = 1, NEPVLU )
      IF ( DEBUG ) WRITE( 6, 1120 ) 'EPVALU',
     *   ( EPVALU( I ), I = 1, NEPVLU )
C
C  INPUT THE SCALE FACTORS FOR THE NONLINEAR ELEMENTS.
C
      READ( INPUT, 1020 ) ( ESCALE( I ), I = 1, NGEL )
      IF ( DEBUG ) WRITE( 6, 1120 ) 'ESCALE',
     *   ( ESCALE( I ), I = 1, NGEL )
C
C  INPUT THE SCALE FACTORS FOR THE GROUPS.
C
      READ( INPUT, 1020 ) ( GSCALE( I ), I = 1, NG )
      IF ( DEBUG ) WRITE( 6, 1120 ) 'GSCALE',
     *   ( GSCALE( I ), I = 1, NG )
C
C  IF A MAXIMUM IS REQUIRED, CHANGE THE SIGN OF THE WEIGHT(S)
C  ASSOCIATED WITH THE OBJECTIVE FUNCTION GROUP(S).
C
      IF ( FINDMX .LT. 0.0 ) THEN
         IF ( IALGOR .EQ. 1 ) THEN
            DO 101 I = 1, NG
               GSCALE( I ) = - GSCALE( I )
  101       CONTINUE
         ELSE
            DO 102 I = 1, NG
               IF ( KNDOFC( I ) .EQ. 1 ) GSCALE( I ) = - GSCALE( I )
  102       CONTINUE
         END IF
      END IF
C
C  INPUT THE SCALE FACTORS FOR THE VARIABLES.
C
      READ( INPUT, 1020 ) ( VSCALE( I ), I = 1, N )
      IF ( DEBUG ) WRITE( 6, 1120 ) 'VSCALE',
     *   ( VSCALE( I ), I = 1, N )
C
C  INPUT THE LOWER AND UPPER BOUNDS ON THE OBJECTIVE FUNCTION.
C
      READ( INPUT, 1080 ) OBFBND( 1 ), OBFBND( 2 )
      IF ( DEBUG ) WRITE( 6, 1180 ) 'OBFBND',
     *    OBFBND( 1 ), OBFBND( 2 )
C
C  SET THE LOWER BOUND FOR THE MINIMIZATION.
C
      IF ( FINDMX .GT. 0.0 ) THEN
         OBJFBN = OBFBND( 1 )
      ELSE
         OBJFBN = - OBFBND( 2 )
      END IF
C
C  INPUT A LOGICAL ARRAY WHICH SAYS WHETHER AN ELEMENT HAS INTERNAL
C  VARIABLES.
C
      READ( INPUT, 1030 ) ( INTREP( I ), I = 1, NELNUM )
      IF ( DEBUG ) WRITE( 6, 1130 ) 'INTREP',
     *   ( INTREP( I ), I = 1, NELNUM )
C
C  INPUT A LOGICAL ARRAY WHICH SAYS WHETHER A GROUP IS TRIVIAL.
C
      READ( INPUT, 1030 ) ( GXEQX( I ), I = 1, NG )
      IF ( DEBUG ) WRITE( 6, 1130 ) 'GXEQX ',
     *   ( GXEQX( I ), I = 1, NG )
C
C  INPUT THE NAMES GIVEN TO THE GROUPS AND TO THE VARIABLES.
C
      READ( INPUT, 1040 ) ( GNAMES( I ), I = 1, NG )
      IF ( DEBUG ) WRITE( 6, 1140 ) 'GNAMES',
     *   ( GNAMES( I ), I = 1, NG )
      READ( INPUT, 1040 ) ( VNAMES( I ), I = 1, N )
      IF ( DEBUG ) WRITE( 6, 1140 ) 'VNAMES',
     *   ( VNAMES( I ), I = 1, N )
C
C  INPUT THE NAMES GIVEN TO THE ELEMENT AND GROUP TYPES.
C
      READ( INPUT, 1040 ) ( ETYPES( I ), I = 1, NELTYP )
      IF ( DEBUG ) WRITE( 6, 1140 ) 'ETYPES',
     *   ( ETYPES( I ), I = 1, NELTYP )
      READ( INPUT, 1040 ) ( GTYPES( I ), I = 1, NGRTYP )
      IF ( DEBUG ) WRITE( 6, 1140 ) 'GTYPES',
     *   ( GTYPES( I ), I = 1, NGRTYP )
C ** Correction 13. 22/12/99: 1 line added **
      READ( INPUT, 1010 ) ( L, I = 1, N )
C
C  IF THE PROBLEM HAS NO OBJECTIVE FUNCTION, IMPLICITLY RESTATE IT AS
C  A LEAST-SQUARES PROBLEM. ENSURE THAT ALL GROUPS ARE NON-TRIVIAL.
C  REMEMBER TO SQUARE THE SCALING FACTORS.
C
      SQUARE = .FALSE.
      IF ( IALGOR .GE. 2 ) THEN
         DO 111 I = 1, NG
            IF ( KNDOFC( I ) .EQ. 1 ) GO TO 112
 111     CONTINUE
         IALGOR = 1
         SQUARE = .TRUE.
         FINDMX = 1.0
         STOPGA = MIN( STOPGA, STOPCA )
         DO 114 I      = 1, NG
            GXEQX( I ) = .FALSE.
            GSCALE( I ) = GSCALE( I ) ** 2
 114     CONTINUE
 112     CONTINUE
      END IF
      IF ( IPRNT .GT. 0 ) THEN
         IF ( IALGOR .EQ. 1 ) WRITE( IOUT, 2100 )
         IF ( IALGOR .EQ. 2 ) WRITE( IOUT, 2110 )
         IF ( IALGOR .EQ. 3 ) WRITE( IOUT, 2120 )
         WRITE( IOUT, 2280 )
         WRITE( IOUT, 2010 ) PNAME
      END IF
C
C  READ A PREVIOUS SOLUTION FILE FOR A RE-ENTRY.
C
      IF ( WARMST .AND. IINPB .GT. 0 ) THEN
         REWIND( IINPB )
         READ( IINPB, 2510 ) PNAME2
         IF ( PNAME2 .NE. PNAME ) THEN
            WRITE( IOUT, 2500 )
            WRITE( IOUT, 2550 ) PNAME2, PNAME
            RETURN
         END IF
         READ( IINPB, 2520 ) I
         IF ( I .NE. N ) THEN
            WRITE( IOUT, 2500 )
            WRITE( IOUT, 2560 ) N, I
            RETURN
         END IF
         READ( IINPB, 2520 ) I
         IF ( I .NE. NG ) THEN
            WRITE( IOUT, 2500 )
            WRITE( IOUT, 2570 ) NG, I
            RETURN
         END IF
         READ( IINPB, 2530 ) RMU
         READ( IINPB, 2590 )
         DO 22 I = 1, N
            READ( IINPB, 2540 ) X( I ), TEMPNA
            IF ( TEMPNA .NE. VNAMES( I ) ) THEN
               WRITE( IOUT, 2500 )
               WRITE( IOUT, 2580 ) TEMPNA
               RETURN
            END IF
   22    CONTINUE
         IF ( IALGOR .GE. 2 ) THEN
            READ( IINPB, 2590 )
            DO 23 I = 1, NG
               IF ( KNDOFC( I ) .GT. 1 ) THEN
                  READ( IINPB, 2540 ) U( I ), TEMPNA
                  IF ( TEMPNA .NE. GNAMES( I ) ) THEN
                     WRITE( IOUT, 2500 )
                     WRITE( IOUT, 2600 ) TEMPNA
                     RETURN
                  END IF
               END IF
   23       CONTINUE
         END IF
      END IF
C
C  IF REQUIRED, TEST THE DERIVATIVES OF THE ELEMENT FUNCTIONS.
C
      IF ( NELNUM .GT. 0 .AND. DECHKE .AND. .NOT. FDGRAD ) THEN
         JUMPTO     = 0
C
C  CHECK THE DERIVATIVES OF THE ELEMENT FUNCTIONS AT THE POINT XT.
C
         DO 30 J = 1, N
CS          IF ( 1.0E+0 .LE. BL( J ) ) THEN
CD          IF ( 1.0D+0 .LE. BL( J ) ) THEN
CS             XT( J ) = BL( J ) + 1.0E-1 *       RANDOM( J )   *
CS   *                   MIN( 1.0E+0, BU( J ) - BL( J ) )
CD             XT( J ) = BL( J ) + 1.0D-1 * DBLE( RANDOM( J ) ) *
CD   *                   MIN( 1.0D+0, BU( J ) - BL( J ) )
            ELSE
CS             IF ( 1.0E+0 .GE. BU( J ) ) THEN
CD             IF ( 1.0D+0 .GE. BU( J ) ) THEN
CS                XT( J ) = BU( J ) - 1.0E-1 *       RANDOM( J )   *
CS   *                      MIN( 1.0E+0, BU( J ) - BL( J ) )
CD                XT( J ) = BU( J ) - 1.0D-1 * DBLE( RANDOM( J ) ) *
CD   *                      MIN( 1.0D+0, BU( J ) - BL( J ) )
               ELSE
CS                XT( J ) = 1.0E+0 + 1.0E-1 *       RANDOM( J )   *
CS   *                      MIN( 1.0E+0, BU( J ) - 1.0E+0 )
CD                XT( J ) = 1.0D+0 + 1.0D-1 * DBLE( RANDOM( J ) ) *
CD   *                      MIN( 1.0D+0, BU( J ) - 1.0D+0 )
               END IF
            END IF
   30    CONTINUE
         IF ( TESTAL ) THEN
C
C  TEST ALL THE NONLINEAR ELEMENT FUNCTIONS.
C
            NCALCF      = NELNUM
            DO 40 J     = 1, NELNUM
               IWK( J ) = J
   40       CONTINUE
         ELSE
C
C  TEST ONE NONLINEAR ELEMENT FUNCTION OF EACH TYPE.
C
            NCALCF               = 0
            DO 50 J              = 1, NELNUM
               IWK( NELNUM + J ) = 0
   50       CONTINUE
            IF ( IPRNT .GT. 0 .AND. IOUT .GT. 0 ) WRITE( IOUT, 2060 )
            DO 60 J   = 1, NELNUM
               IELTYP = ITYPEE( J )
               IF ( IWK( NELNUM + IELTYP ) .EQ. 0 ) THEN
                  IF ( IPRNT .GT. 0 .AND. IOUT .GT. 0 )
     *               WRITE( IOUT, 2040 ) J, ETYPES( IELTYP )
                  IWK( NELNUM + IELTYP ) = 1
                  NCALCF                 = NCALCF + 1
                  IWK( NCALCF )          = J
               END IF
   60       CONTINUE
         END IF
C
C  PARTITION THE WORKSPACE ARRAY.
C
         NINMAX    = 0
         NLMAX     = 0
         DO 70 J   = 1, NELNUM
            NIN    = INTVAR( J )
            NINMAX = MAX( NINMAX, NIN )
            NLMAX  = MAX( NLMAX, ISTAEV( J + 1 ) - ISTAEV( J ) )
   70    CONTINUE
         LXTT  = 1
         LXEL  = LXTT  + NLMAX
         LXINT = LXEL  + NLMAX
         LWRKD = LXINT + NINMAX
         LFTUV = LWRKD + NINMAX * NLMAX + 2 * NINMAX + NLMAX
C
C  CHECK THE DERIVATIVES OF THE NONLINEAR ELEMENT FUNCTIONS.
C  ---------------------------------------------------------
C
  100    CONTINUE
CS       CALL SDRCHE( N, NELNUM, ISTAEV, LSTAEV, ISTADH, LSTADH, IELVAR,
CD       CALL DDRCHE( N, NELNUM, ISTAEV, LSTAEV, ISTADH, LSTADH, IELVAR,
     *                LELVAR, IWK( NELNUM + 1 ), NVARS , INTVAR, LNTVAR,
     *                IWK( NELNUM + NVARS + 1 ), LIWK - NELNUM - NVARS ,
     *                INTREP, LINTRE, ICALCF, LCALCF,
     *                NCALCF, XT    , LXT   , WK( LXTT ), LXEL - LXTT,
     *                FUVALS, LFUVAL, WK( LFTUV    ), LWK - LFTUV,
     *                WK( LXEL ) , LXINT - LXEL, WK( LXINT ),
     *                LWRKD - LXINT, WK( LWRKD ), LFTUV - LWRKD, EPSMCH, 
     *                SECOND, IWK( 1 ), NELNUM, IPRNT , IOUT  ,
     *                RANGES, WARNNG, DEBUGF, JUMPTO )
         IF ( JUMPTO .EQ. 1 ) THEN
            CALL ELFUNS( FUVALS, XT, EPVALU, NCALCF, ITYPEE,
     *                   ISTAEV, IELVAR, INTVAR, ISTADH, ISTEP,
     *                   ICALCF, 1 )
            J = 2
            IF ( SECOND ) J = 3
            CALL ELFUNS( FUVALS, XT, EPVALU, NCALCF, ITYPEE,
     *                   ISTAEV, IELVAR, INTVAR, ISTADH, ISTEP,
     *                   ICALCF, J )
         END IF
         IF ( JUMPTO .EQ. 2 ) THEN
            CALL ELFUNS( FUVALS, WK( LXEL ), EPVALU, NCALCF, ITYPEE,
     *                   ISTAEV, IWK( NELNUM + 1 ), INTVAR, ISTADH,
     *                   ISTEP, ICALCF, 1 )
         END IF
         IF ( JUMPTO .EQ. 3 ) THEN
            CALL ELFUNS( FUVALS, WK( LXEL ), EPVALU, NCALCF, ITYPEE,
     *                   ISTAEV, IWK( NELNUM + 1 ), INTVAR, ISTADH,
     *                   ISTEP, ICALCF, 2 )
         END IF
         IF ( JUMPTO .GT. 0 ) GO TO 100
C
C  STOP IF THERE WERE ANY WARNING MESSAGES.
C
         IF ( WARNNG .AND. FATALE ) THEN
            IF ( IOUT .GT. 0 .AND. IPRNT .EQ. 0 ) WRITE( IOUT, 2370 )
            RETURN
         END IF
      END IF
      IF ( DECHKG .AND. NG .GT. 0 ) THEN
         IF ( TESTAL ) THEN
C
C  CHECK THE DERIVATIVES OF THE GROUP FUNCTIONS AT THE POINTS FT.
C  TEST ALL THE NONTRIVIAL GROUP FUNCTIONS.
C
C ** Correction 14. 01/05/2000: 1 line corrected **
            NCALCG           = 0
            DO 110 J         = 1, NG
C ** Correction 15. 01/05/2000: 2 lines added **
               IF (  ITYPEG( J ) .LE. 0 ) GO TO 110
               NCALCG        = NCALCG + 1
               IWK( NG + J ) = 1
C ** Correction 16. 01/05/2000: 1 line corrected **
               IWK( NCALCG ) = J
CS             FT( J )       =       RANDOM( J )   + 1.0E-1
CD             FT( J )       = DBLE( RANDOM( J ) ) + 1.0D-1
  110       CONTINUE
         ELSE
C
C  TEST ONE NONTRIVIAL GROUP FUNCTION OF EACH TYPE.
C
            NCALCG           = 0
            DO 120 J         = 1, NG
               IWK( NG + J ) = 0
  120       CONTINUE
            IF ( IPRNT .GT. 0 ) WRITE( IOUT, 2060 )
            DO 130 J  = 1, NG
               IGRTYP = ITYPEG( J )
               IF ( IGRTYP .LE. 0 ) GO TO 130
               IF ( IWK( NG + IGRTYP ) .EQ. 0 ) THEN
                  IWK( NG + IGRTYP ) = 1
                  NCALCG        = NCALCG + 1
                  IWK( NCALCG ) = J
                  IF ( IPRNT .GT. 0 ) WRITE( IOUT, 2050 )
     *                                J, GTYPES( IGRTYP )
CS                FT( J ) =       RANDOM( J )   + 1.0E-1
CD                FT( J ) = DBLE( RANDOM( J ) ) + 1.0D-1
               END IF
  130       CONTINUE
         END IF
         LFTT   = 1
         LGTSTR = LFTT + NG
         LGTVAL = 3 * NG
C
C  CHECK THE DERIVATIVES OF THE GROUP FUNCTIONS.
C  ---------------------------------------------
C
  200    CONTINUE
CS       CALL SDRCHG( NG, FT( 1 ), GVALS, LGVALS, WK( LFTT ),
CD       CALL DDRCHG( NG, FT( 1 ), GVALS, LGVALS, WK( LFTT ),
     *                WK( LGTSTR ), LGTVAL, IWK, NCALCG, EPSMCH,
     *                IPRNT, IOUT, WARNNG, DEBUGG, JUMPTO )
         IF ( JUMPTO .EQ. 1 .OR. JUMPTO .EQ. 2 )
     *      CALL GROUPS( GVALS, LGVALS, FT, GPVALU,
     *                   NCALCG, ITYPEG, ISTGP, IWK, .FALSE. )
         IF ( JUMPTO .EQ. 1 .OR. JUMPTO .EQ. 3 )
     *      CALL GROUPS( GVALS, LGVALS, FT, GPVALU,
     *                   NCALCG, ITYPEG, ISTGP, IWK, .TRUE. )
         IF ( JUMPTO .GT. 0 ) GO TO 200
C
C  STOP IF THERE WERE ANY WARNING MESSAGES.
C
         IF ( WARNNG .AND. FATALG ) THEN
            IF ( IOUT .GT. 0 .AND. IPRNT .EQ. 0 ) WRITE( IOUT, 2380 )
            RETURN
         END IF
      END IF
C
C  OBTAIN APPROPRIATE VARIABLE AND GROUP SCALINGS, IF REQUIRED.
C
      IF ( GETSCA ) THEN
C
C  SET UP REAL WORKSPACE ADDRESSES.
C
         LVRSCA = 1
         NVAR   = NUMVAR
         NORDER = MAX( NVAR, NG )
         LGPSCA = LVRSCA + NORDER
         LWKST  = LGPSCA + NORDER
C ** Correction 8. 10/08/93: 1 line corrected **
         LW     = 2 * NG + 3 * NVAR
C ** Correction 8. 10/08/93: end of correction **
         LSGRJA = LWKST + LW
C
C  SET UP INTEGER WORKSPACE ADDRESSES.
C
         LSSTIV = 1
         LSTINV = MAX( NVAR, NELNUM + 1 )
         LSTAGV = LSSTIV + LSTINV
         LSSVGR = LSTAGV + NG + 1
         LSVGRP = ( LIWK - LSSVGR ) / 3
         LGRJAC = MIN( LSVGRP, LWK - LSGRJA )
         LRNGRJ = LSSVGR + LSVGRP
         LCNGRJ = LRNGRJ + LGRJAC
         SCALED = .TRUE.
C
C  CHECK THE DERIVATIVES OF THE ELEMENT FUNCTIONS AT THE POINT XT.
C
         DO 210 J   = 1, N
            XT( J ) = MAX( BL( J ), MIN( BU( J ), X( J ) ) )
  210    CONTINUE
         INFORM = 0
C
C  CALCULATE THE SCALINGS.
C
  220    CONTINUE
         IPRNTS = IPRNT
CS       CALL SSCALN( N , NG, NSLACK, NELNUM, XT    , LXT   , NCALCF,
CD       CALL DSCALN( N , NG, NSLACK, NELNUM, XT    , LXT   , NCALCF,
     *                ICALCF, LCALCF, ICNA  , LICNA , ISTADA, LSTADA,
     *                IELING, LELING, IWK( LRNGRJ ) , IWK( LCNGRJ ) ,
     *                INTVAR, LNTVAR, ISTADG, LSTADG, ISTAEV, LSTAEV,
     *                IELVAR, LELVAR, IWK( LSSTIV ) , LSTINV,
     *                IWK( LSSVGR ) , LSVGRP, ISTADH, LSTADH,
     *                IWK( LSTAGV ) , NG + 1, A , LA, B , LB, FT,
     *                LFT   , GVALS( 1, 2 ) , LGVALS, FUVALS, LFUVAL,
     *                STOPGA, STOPCA, GSCALE, LGSCAL, ESCALE, LESCAL,
     *                VSCALE, LVSCAL, SCALEG, SCALEV, KNDOFC, LKNDOF,
     *                WK ( LSGRJA ) , LGRJAC, WK ( LVRSCA ) , NORDER,
C ** Correction 11. 13/12/97: 2 lines corrected **
     *                WK ( LGPSCA ) , NORDER, GETSCA,  WK ( LWKST  ) ,
     *                LW    , WK ( LWKST  ) , LWK  - LWKST  ,
     *                IWK( LRNGRJ ) , LIWK - LRNGRJ ,
     *                GXEQX , LGXEQX, INTREP, LINTRE,  RANGES,
     *                IOUT  , IPRNTS, FDGRAD, INFORM )
         IF ( INFORM .LT. 0 ) THEN
C
C  FURTHER PROBLEM INFORMATION IS REQUIRED.
C
            IF ( INFORM .EQ. - 1 ) THEN
C
C  EVALUATE THE ELEMENT FUNCTION AND DERIVATIVE VALUE.
C
               CALL ELFUNS( FUVALS, XT, EPVALU, NCALCF, ITYPEE, ISTAEV,
     *                      IELVAR, IWK( LSSTIV ), ISTADH, ISTEP,
     *                      ICALCF, 1 )
               IF ( .NOT. FDGRAD ) CALL
     *              ELFUNS( FUVALS, XT, EPVALU, NCALCF, ITYPEE, ISTAEV,
     *                      IELVAR, IWK( LSSTIV ), ISTADH, ISTEP,
     *                      ICALCF, 2 )
            END IF
            IF ( INFORM .EQ. - 2 ) THEN
C
C  EVALUATE THE GROUP FUNCTION DERIVATIVES.
C
               IF ( SQUARE ) THEN
                  DO 230 J         = 1, NCALCF
                     I             = ICALCF( J )
                     GVALS( I, 1 ) = FT( I )
 230              CONTINUE
                  CALL GROUPS( GVALS, LGVALS, FT, GPVALU, NCALCF,
     *                         ITYPEG, ISTGP, ICALCF, .FALSE. )
                  DO 240 J         = 1, NCALCF
                     I             = ICALCF( J )
                     U( I )        = GVALS( I, 1 )
                     GVALS( I, 1 ) = GVALS( I, 1 ) ** 2
                     GVALS( I, 2 ) = 1.0
                     GVALS( I, 3 ) = 0.0
 240              CONTINUE
               END IF
               CALL GROUPS( GVALS, LGVALS, FT, GPVALU,
     *                      NCALCF, ITYPEG, ISTGP, ICALCF, .TRUE. )
               IF ( SQUARE ) THEN
                  DO 250 J         = 1, NCALCF
                     I             = ICALCF( J )
                     GVALS( I, 3 ) = 2.0 * U( I ) * GVALS( I, 3 )
     *                               + 2.0 * GVALS( I, 2 ) ** 2
                     GVALS( I, 2 ) = 2.0 * U( I ) * GVALS( I, 2 )
 250              CONTINUE
               END IF
            END IF
            GO TO 220
         END IF
         IF ( INFORM .GT. 0 ) RETURN
C ** Correction 12. 16/01/98: 9 lines added **
C
C  IF REQUIRED, SCALE THE INITIAL LAGRANGE MULTIPLIER ESTIMATES.
C
         IF ( SCALEG ) THEN
            DO 260 IG = 1, NG
               IF ( KNDOFC( IG ) .NE. 1 )
     *              U( IG ) = U( IG ) / WK ( LGPSCA + IG - 1 )
  260       CONTINUE
         END IF
C ** Correction 12. 16/01/98: 1 line added **
      ELSE
         SCALED = .FALSE.
      END IF
C
C  PREPARE FOR THE MINIMIZATION.
C
      INFORM = 0
      TIME   = 0.0
      TTOTAL = CPUTIM( DUM )
      IF ( FINDMX .GT. 0.0 ) THEN
         MINMAX = 'MIN'
         IF ( IPRNT .GT. 0 ) WRITE( IOUT, 2340 )
      ELSE
         MINMAX = 'MAX'
         IF ( IPRNT .GT. 0 ) WRITE( IOUT, 2350 )
      END IF
C
C  CALL THE MINIMIZER.
C  -------------------
C
  300 CONTINUE
      TIMM = CPUTIM( DUM )
      IF ( ITER .GE. IPSTRT ) IPRINT = IPRNT
      IF ( MOD( ITER - IPSTRT, IPGAP ) .NE. 0 .OR.
     *     ITER .GE. IPSTOP ) IPRINT = 0
      IF ( IALGOR .EQ. 1 ) THEN
CS       CALL SSBMIN( N , NG, NELNUM, IELING, LELING, ISTADG, LSTADG,
CD       CALL DSBMIN( N , NG, NELNUM, IELING, LELING, ISTADG, LSTADG,
     *                IELVAR, LELVAR, ISTAEV, LSTAEV, INTVAR, LNTVAR,
     *                ISTADH, LSTADH, ICNA  , LICNA , ISTADA, LSTADA,
     *                A , LA, B , LB, BL    , LBL   , BU    , LBU   ,
     *                GSCALE, LGSCAL, ESCALE, LESCAL, VSCALE, LVSCAL,
     *                GXEQX , LGXEQX, INTREP, LINTRE, RANGES, INFORM,
     *                FOBJ  , X , LX, GVALS , LGVALS, FT    , LFT   ,
     *                FUVALS, LFUVAL, XT    , LXT   , ICALCF, LCALCF,
     *                NCALCF, ICALCG, LCALCG, NCALCG, IVAR  , LIVAR ,
     *                NVAR  , Q , LQ, DGRAD , LDGRAD, ICHOSE, ITER  ,
     *                MAXIT , QUADRT, STOPGA, OBJFBN, IWK   , LIWK  ,
     *                WK    , LWK   , IPRINT, IOUT  )
      END IF
      IF ( IALGOR .EQ. 2 ) THEN
CS       CALL SAUGLG( N, NG , NELNUM, IELING, LELING, ISTADG, LSTADG,
CD       CALL DAUGLG( N, NG , NELNUM, IELING, LELING, ISTADG, LSTADG,
     *                IELVAR, LELVAR, ISTAEV, LSTAEV, INTVAR, LNTVAR,
     *                ISTADH, LSTADH, ICNA  , LICNA , ISTADA, LSTADA,
     *                A , LA, B , LB, BL    , LBL   , BU    , LBU   ,
     *                GSCALE, LGSCAL, ESCALE, LESCAL, VSCALE, LVSCAL,
     *                GXEQX , LGXEQX, INTREP, LINTRE, KNDOFC, LKNDOF,
     *                RANGES, INFORM, FOBJ  , X , LX, U , LU, GVALS ,
     *                LGVALS, FT    , LFT   , FUVALS, LFUVAL, XT    ,
     *                LXT   , ICALCF, LCALCF, NCALCF, ICALCG, LCALCG,
     *                NCALCG, IVAR  , LIVAR , NVAR  , Q , LQ, DGRAD ,
     *                LDGRAD, ICHOSE, ITER  , MAXIT , QUADRT, VNAMES,
     *                LVNAME, GNAMES, LGNAME, STOPGA, STOPCA, IWK   ,
     *                LIWK  , WK    , LWK   , IPRINT, IOUT  )
      END IF
      IF ( IALGOR .EQ. 3 ) THEN
CS       CALL SBARIA( N, NG , NELNUM, IELING, LELING, ISTADG, LSTADG,
CD       CALL DBARIA( N, NG , NELNUM, IELING, LELING, ISTADG, LSTADG,
     *                IELVAR, LELVAR, ISTAEV, LSTAEV, INTVAR, LNTVAR,
     *                ISTADH, LSTADH, ICNA  , LICNA , ISTADA, LSTADA,
     *                A , LA, B , LB, BL    , LBL   , BU    , LBU   ,
     *                GSCALE, LGSCAL, ESCALE, LESCAL, VSCALE, LVSCAL,
     *                GXEQX , LGXEQX, INTREP, LINTRE, KNDOFC, LKNDOF,
     *                RANGES, INFORM, FOBJ  , X , LX, U , LU, GVALS ,
     *                LGVALS, FT    , LFT   , FUVALS, LFUVAL, XT    ,
     *                LXT   , ICALCF, LCALCF, NCALCF, ICALCG, LCALCG,
     *                NCALCG, IVAR  , LIVAR , NVAR  , Q , LQ, DGRAD ,
     *                LDGRAD, ICHOSE, ITER  , MAXIT , QUADRT, VNAMES,
     *                LVNAME, GNAMES, LGNAME, STOPGA, STOPCA, IWK   ,
     *                LIWK  , WK    , LWK   , IPRINT, IOUT  )
      END IF
      TIMM = CPUTIM( DUM ) - TIMM
      TIME = TIME + TIMM
C
C  CHECK TO SEE IF THE USER WISHES TO CONTINUE THE MINIMIZATION.
C
CIBM  INQUIRE( UNIT =  IALIVE    , EXIST = ALIVE )
CUNIX INQUIRE( FILE = 'ALIVE.d'  , EXIST = ALIVE )
CVMS  INQUIRE( FILE = 'ALIVE.DAT', EXIST = ALIVE )
      IF ( INFORM .LT. 0 .AND. .NOT. ALIVE ) INFORM = 10
C
C  WRITE THE SOLUTION FILE FOR POSSIBLE RE-ENTRY.
C
C ** Correction 2. 29/01/93: 4 lines replaced by 7 **
      IF ( ISTORE .GT. 0 ) THEN
         DSAVE = IOUTPB .GT. 0 .AND. ( MOD( ITER, ISTORE ) .EQ. 0
     *           .AND. ( INFORM .EQ. - 2 .OR. INFORM .EQ. - 4 ) )
      ELSE
         DSAVE = IOUTPB .GT. 0 .AND. INFORM .GE. 0 .AND. ISTORE .EQ. 0
      END IF
      IF ( DSAVE ) THEN
C ** Correction 2. 29/01/93: end of correction **
            REWIND( IOUTPB )
            WRITE( IOUTPB, 2450 ) PNAME, N, NG, RMU
            WRITE( IOUTPB, 2300 )
            DO 310 I = 1, N
               WRITE( IOUTPB, 2460 ) X( I ), VNAMES( I )
  310       CONTINUE
            IF ( IALGOR .GE. 2 ) THEN
               WRITE( IOUTPB, 2290 )
               DO 320 I = 1, NG
                  IF ( KNDOFC( I ) .GT. 1 )
     *               WRITE( IOUTPB, 2460 ) U( I ), GNAMES( I )
  320          CONTINUE
            END IF
C ** Correction 3. 29/01/93: 1 line deleted **
C ** Correction 3. 29/01/93: end of correction **
      END IF
C
C  THE USER WISHES TO TERMINATE THE MINIMIZATION.
C
      IF ( .NOT. ALIVE ) THEN
         IF ( IOUT .GT. 0 ) WRITE( IOUT, 2070 )
         GO TO 500
      END IF
C
C  IF THERE IS INSUFFICIENT WORKSPACE, TERMINATE EXECUTION.
C
      IF ( INFORM .EQ. 4 .OR. INFORM .EQ. 5 .OR.
     *     INFORM .EQ. 6 .OR. INFORM .EQ. 7 ) GO TO 500
C
C  IF THE APPROXIMATION TO THE SOLUTION HAS CHANGED APPRECIABLY,
C  PRINT OUT THE SOLUTION IN SIF FORMAT.
C
      IF ( NEWSOL .OR. INFORM .EQ. 0 .OR. INFORM .EQ. 3 ) THEN
         NEWSOL = .FALSE.
         IF ( IOUTSL .GT. 0 ) THEN
            REWIND( IOUTSL )
            WRITE( IOUTSL, 2250 ) PNAME, RMU
            WRITE( IOUTSL, 2300 )
            DO 330 I = 1, N
               WRITE( IOUTSL, 2260 ) VNAMES( I ), X( I )
  330       CONTINUE
            IF ( IALGOR .GE. 2 ) THEN
               WRITE( IOUTSL, 2290 )
               NOBJGR   = 0
               DO 340 I = 1, NG
                  IF ( KNDOFC( I ) .GT. 1 ) THEN
                     WRITE( IOUTSL, 2260 ) GNAMES( I ), U( I )
                  ELSE
                     NOBJGR = NOBJGR + 1
                  END IF
  340          CONTINUE
               IF ( NOBJGR .GT. 0 ) THEN
                  IF ( FINDMX .GT. 0.0 ) THEN
                     WRITE( IOUTSL, 2270 ) FOBJ
                  ELSE
                     WRITE( IOUTSL, 2360 ) - FOBJ
                  END IF
               END IF
            ELSE
               IF ( FINDMX .GT. 0.0 ) THEN
                  WRITE( IOUTSL, 2270 ) FOBJ
               ELSE
                  WRITE( IOUTSL, 2360 ) - FOBJ
               END IF
            END IF
         END IF
      END IF
      IF ( INFORM .LT. 0 ) THEN
C
C  FURTHER PROBLEM INFORMATION IS REQUIRED.
C
         IF ( INFORM .EQ. - 1 .OR. INFORM .EQ. - 3 .OR.
     *        INFORM .EQ. - 7 ) THEN
            IF ( IPRINT .GE. 10 ) WRITE( IOUT, 2200 )
C
C  EVALUATE THE ELEMENT FUNCTION VALUES.
C
            CALL ELFUNS( FUVALS, XT, EPVALU, NCALCF, ITYPEE, ISTAEV,
     *                   IELVAR, INTVAR, ISTADH, ISTEP,
     *                   ICALCF, 1 )
         END IF
         IF ( ( INFORM .EQ. - 1 .OR. INFORM .EQ. - 5 .OR.
     *        INFORM .EQ. - 6 ) .AND. .NOT. FDGRAD ) THEN
            IFFLAG = 2
            IF ( SECOND ) IFFLAG = 3
C
C  EVALUATE THE ELEMENT FUNCTION DERIVATIVES.
C
            IF ( IPRINT .GE. 10 ) WRITE( IOUT, 2210 )
            CALL ELFUNS( FUVALS, XT, EPVALU, NCALCF, ITYPEE, ISTAEV,
     *                   IELVAR, INTVAR, ISTADH, ISTEP,
     *                   ICALCF, IFFLAG )
         END IF
         IF ( INFORM .EQ. - 2 .OR. INFORM .EQ. - 4 ) THEN
C
C  EVALUATE THE GROUP FUNCTION VALUES.
C
            IF ( IPRINT .GE. 10 ) WRITE( IOUT, 2220 )
            IF ( IPRINT .GE. 100 ) WRITE( IOUT, 2390 )
     *           ( FT( I ) , I = 1, NG )
            IF ( SQUARE ) THEN
               DO 350 J         = 1, NCALCG
                  I             = ICALCG( J )
                  GVALS( I, 1 ) = FT( I )
 350           CONTINUE
            END IF
            CALL GROUPS( GVALS, LGVALS, FT, GPVALU,
     *                   NCALCG, ITYPEG, ISTGP, ICALCG, .FALSE. )
            IF ( SQUARE ) THEN
               DO 360 J         = 1, NCALCG
                  I             = ICALCG( J )
                  U( I )        = GVALS( I, 1 )
                  GVALS( I, 1 ) = GVALS( I, 1 ) ** 2
 360           CONTINUE
            END IF
            IF ( IPRINT .GE. 100 ) WRITE( IOUT, 2400 )
     *           ( GVALS( I, 1 ) , I = 1, NG )
         END IF
         IF ( INFORM .EQ. - 2 .OR. INFORM .EQ. - 5 ) THEN
C
C  EVALUATE THE GROUP FUNCTION DERIVATIVES.
C
            IF ( IPRINT .GE. 10 ) WRITE( IOUT, 2230 )
            IF ( SQUARE ) THEN
               DO 370 J         = 1, NCALCG
                  I             = ICALCG( J )
                  GVALS( I, 2 ) = 1.0
                  GVALS( I, 3 ) = 0.0
 370           CONTINUE
            END IF
            CALL GROUPS( GVALS, LGVALS, FT, GPVALU,
     *                   NCALCG, ITYPEG, ISTGP, ICALCG, .TRUE. )
            IF ( SQUARE ) THEN
               DO 380 J         = 1, NCALCG
                  I             = ICALCG( J )
                  GVALS( I, 3 ) = 2.0 * U( I ) * GVALS( I, 3 )
     *                            + 2.0 * GVALS( I, 2 ) ** 2
                  GVALS( I, 2 ) = 2.0 * U( I ) * GVALS( I, 2 )
 380           CONTINUE
            END IF
         END IF
         IF ( INFORM .LE. - 8 ) THEN
C
C  EVALUATE THE PRECONDITIONED GRADIENT- USE A DIAGONAL PRECONDITIONER.
C
            IF ( IPRINT .GE. 10 ) WRITE( IOUT, 2240 )
            DO 390 J = 1, NVAR
               Q( IVAR( J ) ) = DGRAD( J )
  390       CONTINUE
         END IF
         GO TO 300
      END IF
C
C  THE MINIMIZATION HAS BEEN TERMINATED.
C  -------------------------------------
C
      IF ( SQUARE ) THEN
         IF ( FOBJ .GT. STOPCA ) THEN
            INFORM = 8
            IF ( IPRNT .GT. 0 .AND. IOUT .GT. 0 ) WRITE( IOUT, 2160 )
         END IF
      END IF
      IF ( IPRNT .GT. 0 .AND. IOUT .GT. 0 .AND. IALGOR .EQ. 1 ) THEN
         LGFX = ISTADH( NELNUM + 1 ) - 1
         WRITE( IOUT, 2020 )
C ** Correction 4. 28/07/93: 20 lines replaced by 49 **
         L = N
C        L = 2
         DO 420 J = 1, 2
            IF ( J .EQ. 1 ) THEN
               IR = 1
               IC = MIN( L, N )
            ELSE
               IF ( IC. LT. N - L ) WRITE( IOUT, 2080 )
               IR = MAX( IC + 1, N - IC + 1 )
               IC = N
            END IF
            DO 410 I = IR, IC
               ST    = ' FREE'
               IF ( X( I ) .GE. BU( I ) - EPSMCH ) ST = 'UPPER'
               IF ( X( I ) .LE. BL( I ) + EPSMCH ) ST = 'LOWER'
CS             IF ( BU( I ) - BL( I ) .LE. 2.0E+0 * EPSMCH) ST = 'FIXED'
CD             IF ( BU( I ) - BL( I ) .LE. 2.0D+0 * EPSMCH) ST = 'FIXED'
               IF ( MAXIT .GE. 0 ) THEN
                  WRITE( IOUT, 2030 ) VNAMES( I ), ST, X( I ),
     *                           BL( I ), BU( I ), FUVALS( LGFX + I)
               ELSE
                  WRITE( IOUT, 2150 ) VNAMES( I ), ST, X( I ),
     *                           BL( I ), BU( I )
               END IF
  410       CONTINUE
  420    CONTINUE
         IF ( SQUARE ) THEN
            WRITE( IOUT, 2130 )
            L = NG
C           L = 2
            DO 440 J = 1, 2
               IF ( J .EQ. 1 ) THEN
                  IR = 1
                  IC = MIN( L, NG )
               ELSE
                  IF ( IC. LT. NG - L ) WRITE( IOUT, 2440 )
                  IR = MAX( IC + 1, NG - IC + 1 )
                  IC = NG
               END IF
               DO 430 I = IR, IC
                  IF ( MAXIT .GE. 0 ) THEN
                     WRITE( IOUT, 2140 ) GNAMES( I ), I,
     *                    U( I ), GSCALE( I ), 2.0 * U( I )
                  ELSE
                     WRITE( IOUT, 2430 ) GNAMES( I ), I,
     *                    U( I ), GSCALE( I )
                  END IF
  430          CONTINUE
  440       CONTINUE
C ** Correction 4. 28/07/93: end of correction **
         ELSE
            WRITE( IOUT, 2180 ) FOBJ * FINDMX
         END IF
      END IF
  500 CONTINUE
      TTOTAL = CPUTIM( DUM ) - TTOTAL
      IF ( IOUT .GT. 0 .AND. IPRNT .EQ. 0 ) THEN
         WRITE( IOUT, 2180 ) FOBJ * FINDMX
         DO 510 I = 1, N
            WRITE( IOUT, 2170 ) VNAMES( I ), X( I )
  510    CONTINUE
         IF ( INFORM .NE. 0 ) WRITE( IOUT, 2190 ) INFORM
      END IF
      IF ( IOUT .GT. 0 .AND. IPRNT .GT. 0 ) THEN
         WRITE( IOUT, 2000 ) INFORM, ITER, TIME, TTOTAL - TIME
         IF ( FINDMX .GT. 0.0 ) THEN
            WRITE( IOUT, 2340 )
         ELSE
            WRITE( IOUT, 2350 )
         END IF
         WRITE( IOUT, 2010 ) PNAME
      END IF
C
C  WRITE THE SOLUTION SUMMARY FILE.
C
      IF ( IOUT1L .GT. 0 ) THEN
         READ( INPUT, 1090 ) PNAME,
     *          NFREE , NFIXED, NLOWER, NUPPER, NBOTH , NSLACK,
     *          NLINOB, NNLNOB, NLINEQ, NNLNEQ, NLININ, NNLNIN
         WRITE( IOUT1L, 2090 ) PNAME,
     *          NFREE , NFIXED, NLOWER, NUPPER, NBOTH , NSLACK,
     *          NLINOB, NNLNOB, NLINEQ, NNLNEQ, NLININ, NNLNIN
         IF ( IALGOR .EQ. 1 ) THEN
            OPTIMI = 'SBMIN'
         END IF
         IF ( IALGOR .EQ. 2 ) THEN
            OPTIMI = 'AUGLG'
         END IF
         IF ( IALGOR .EQ. 3 ) THEN
            OPTIMI = 'BARIA'
         END IF
         WRITE( IOUT1L, 2310 ) PNAME, NUMVAR,
     *      MINMAX, OPTIMI, ICHOSE( 1 ), ICHOSE( 2 ), ICHOSE( 3 ),
     *      ICHOSE( 4 ), ICHOSE( 5 ), ICHOSE( 6 ), ITER,
     *      NGEVAL, ITERCG, INFORM, PNAME, TTOTAL, FOBJ * FINDMX
      END IF
      IF ( INFORM .LT. 10 ) GO TO 10
C
C  END OF EXECUTION.
C
  900 CONTINUE
      RETURN
C
C  NON-EXECUTABLE STATEMENTS.
C
C ** Correction -1b. 03/03/00: Integer formats increased
 1002 FORMAT( I8 )
 1010 FORMAT( ( 10I8 ) )
 1020 FORMAT( ( 1P, 4D16.8 ) )
 1030 FORMAT( ( 72L1 ) )
 1040 FORMAT( ( 8A10 ) )
 1080 FORMAT( 1P, 2D16.8 )
C ** Correction -1c. 03/03/00: Integer formats increased
 1090 FORMAT( A8, 12I8 )
 1100 FORMAT( A8, 3I8 )
 1110 FORMAT( 1X, A6, /, ( 1X, 10I8 ) )
 1120 FORMAT( 1X, A6, /, ( 1X, 1P, 4D16.8 ) )
 1130 FORMAT( 1X, A6, /, ( 1X, 72L1 ) )
 1140 FORMAT( 1X, A6, /, ( 1X, 8A10 ) )
 1180 FORMAT( 1X, A6, /, 1P, 2D16.6 )
 2000 FORMAT( /, ' INFORM         = ', I16,
     *           ' Number of iterations = ', I16,
     *        /, ' Time(LANCELOT) = ', 0P, F16.2,
     *           ' Time(other)          = ', 0P, F16.2 )
 2010 FORMAT( /, ' ************* Problem ', A8, ' *****************' )
 2020 FORMAT( /, ' Variable name Status    Value',
     *           '    Lower bound Upper bound |  Dual value ',
     *        /, ' ------------- ------    -----',
     *           '    ----------- ----------- |  ----------' )
 2030 FORMAT( 2X, A10, 4X, A5, 1P, 3D12.4, ' |', D12.4 )
 2040 FORMAT( ' Element number ', I6, ' chosen as representative of',
     *        ' element type ', A10 )
 2050 FORMAT( ' Group number ', I6, ' chosen as representative of',
     *        ' group type ', A10 )
 2060 FORMAT( / )
 2070 FORMAT( /, ' ** User wishes to terminate execution. Stopping ' )
C ** Correction 5. 28/07/93: 3 lines added **
 2080 FORMAT( '  .             .....',
     *        ' ........... ........... ...........',
     *        ' | ...........' )
C ** Correction 5. 28/07/93: end of correction **
C ** Correction -1d. 03/03/00: Integer formats increased
 2090 FORMAT( A8, 12I8 )
 2100 FORMAT( /, ' *-*-*-*-* LANCELOT A -*- SBMIN Minimizer *-*-*-*-*' )
 2110 FORMAT( /, ' *-*-*-*-* LANCELOT A -*- AUGLG Minimizer *-*-*-*-*' )
 2120 FORMAT( /, ' *-*-*-*-* LANCELOT A -*- BARIA Minimizer *-*-*-*-*' )
 2130 FORMAT( /, ' Constraint name Number    Value    Scale factor ',
     *           '| Lagrange multiplier',
     *        /, ' --------------- ------    -----    ----- ------ ',
     *           '| -------------------')
 2140 FORMAT( 4X, A10, I7, 2X, 1P, 2D12.4, '  |   ', D12.4 )
 2150 FORMAT( 2X, A10, 4X, A5, 1P, 3D12.4, ' |      - ' )
 2160 FORMAT( /, ' Constraint violations are large. Problem',
     *           ' possibly infeasible. ' )
 2170 FORMAT( 12X, A10, 6X, 1P, D22.14 )
 2180 FORMAT( /, ' objective function value = ', 1P, D22.14, / )
 2190 FORMAT( /, ' ** Warning. Exit from lancelot with INFORM = ', I3 )
 2200 FORMAT( /, ' Evaluating element functions ' )
 2210 FORMAT( /, ' Evaluating derivatives of element functions ' )
 2220 FORMAT( /, ' Evaluating group functions ' )
 2230 FORMAT( /, ' Evaluating derivatives of group functions ' )
 2240 FORMAT( /, ' Evaluating user supplied  preconditioner ' )
 2250 FORMAT( '*   Lancelot solution for problem name: ', A8, /,
     *        '*   penalty parameter value is ', 1P, D12.4 )
 2260 FORMAT( '    SOLUTION  ', A10, 1P, D12.5 )
 2270 FORMAT( /, ' XL SOLUTION  ', 10X, 1P, D12.5 )
C ** Correction 10. 06/08/97: 4 lines replaced by 8 **
 2280 FORMAT( /, ' Copyright CGT productions, 1991',
     *        //, ' - Use of this code is restricted to those who',
     *        /,  ' - agree to abide by the conditions-of-use',
     *        /,  ' - set out in the CONDIT.USE file distributed',
     *        /,  ' - with the source to the LANCELOT codes or from',
     *            ' the WWW at',
     *        /,  ' - http://www.numerical.rl.ac.uk',
     *            '/lancelot/blurb.html' )
C ** Correction 17. 20/02/2000: previous 2 lines updated **
C ** Correction 10. 06/08/97: end of correction **
 2290 FORMAT( /, '*   Lagrange multipliers ', / )
 2300 FORMAT( /, '*   variables ', / )
 2310 FORMAT( A8, I6, 1X, A3, 1X, A5, '( ', I1, ' ', I2, ' ', I1,
     *        ' ', I1, ' ', I1, ' ', I1, ' ) ', 3I8, I4,
     *        /, A8, ' TIME ', 0P, F12.2,  ' FINAL FUNCTION VALUE ',
     *           1P, D12.4 )
 2340 FORMAT( /, ' *-*-*-*-*-*-* Minimizer sought *-*-*-*-*-*-*-*-*' )
 2350 FORMAT( /, ' *-*-*-*-*-*-* Maximizer sought *-*-*-*-*-*-*-*-*' )
 2360 FORMAT( /, ' XU SOLUTION  ', 10X, 1P, D12.5 )
 2370 FORMAT( /, ' Possible error in element derivative. Stopping ' )
 2380 FORMAT( /, ' Possible error in group derivative. Stopping ' )
 2390 FORMAT( /, ' Group values ', /, ( 1P, 6D12.4 ) )
 2400 FORMAT( /, ' Group function values ', /, ( 1P, 6D12.4 ) )
C ** Correction 7. 28/07/93: 3 lines added **
 2430 FORMAT( 4X, A10, I7, 2X, 1P, 2D12.4, '  |        - ' )
 2440 FORMAT( '    .               .   ........... ...........',
     *        '  |    ........... ' )
C ** Correction 7. 28/07/93: end of correction **
 2450 FORMAT( 16X, A8, '  problem name ', /,
     *        16X, I8, '  number of variables ', /,
     *        16X, I8, '  number of groups ', /,
     *        1P, D24.16, '  penalty parameter value ' )
 2460 FORMAT( 1P, D24.16, 2X, A10 )
 2500 FORMAT( /, ' Re-entry requested ' )
 2510 FORMAT( 16X, A8 )
 2520 FORMAT( 16X, I8 )
 2530 FORMAT( 1P, D24.16 )
 2540 FORMAT( 1P, D24.16, 2X, A10 )
 2550 FORMAT( /, ' *** Exit from LANCE: re-entry requested with data ',
     *        ' for problem ', A8, /,
     *        '     while the most recently decoded problem is ', A8 )
 2560 FORMAT( /, ' *** Exit from LANCE: number of variables changed',
     *        ' from ', I8, ' to ', I8, ' on re-entry ' )
 2570 FORMAT( /, ' *** Exit from LANCE: number of groups changed ',
     *        ' from ', I8, ' to ', I8, ' on re-entry ' )
 2580 FORMAT( /, ' *** Exit from LANCE: variable named ', A10,
     *        ' out of order on re-entry ' )
 2590 FORMAT( /, / )
 2600 FORMAT( /, ' *** Exit from LANCE: Lagrange multiplier for ', A10,
     *        ' out of order on re-entry ' )
      END
