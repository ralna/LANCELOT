C ** Correction report.
C ** Correction 1. 23/03/93: 2 lines corrected **
C ** Correction 2. 15/06/93: 1 lines corrected **
C ** Correction 3. 15/06/93: 1 line added **
C ** Correction 4. 15/06/93: 1 line added **
C ** Correction 5. 15/06/93: 1 line corrected **
C ** Correction 6. 15/06/93: 1 line corrected **
C ** Correction 7. 15/06/93: 1 line corrected **
C ** Correction 8. 15/06/93: 1 line corrected **
C ** Correction 9. 13/01/94: 2 lines corrected **
C ** Correction 10. 19/08/96: 5 lines added **
C ** Correction 11. 05/11/97: 1 line corrected **
C ** Correction 12. 19/08/98: added final argument to call to ASMBL.
C ** End of Correction report.
C  THIS VERSION: 19/08/1998 AT 18:00:51 PM.
CS    SUBROUTINE SPRECN( IFACTR, MUNKS , BAND  , SEPREC, N     , NG    ,     
CD    SUBROUTINE DPRECN( IFACTR, MUNKS , BAND  , SEPREC, N     , NG    ,     
     *                   MAXSEL, NFREEF, NFIXED, REFACT, NVAR  , IVAR  ,   
     *                   ISTADH, LSTADH, ICNA  , LICNA , ISTADA, LSTADA, 
     *                   INTVAR, LNTVAR, IELVAR, LELVAR, IELING, LELING, 
     *                   ISTADG, LSTADG, ISTAEV, LSTAEV, ISTAGV, LNSTGV, 
     *                   ISVGRP, LNVGRP, IWK   , LIWK  , A     , LA    , 
     *                   GUVALS, LNGUVL, HUVALS, LNHUVL, GVALS2, GVALS3, 
     *                   GRAD  , Q     , GSCALE, ESCALE, LESCAL, WK    , 
     *                   LWK   , GXEQX , LGXEQX, INTREP, LINTRE, RANGES,
     *                   IPRINT, IOUT  , INFORM )
C
C  ********************************************************************
C
C  FORM THE PRODUCT Q = M ** -1 GRAD, WHERE M IS A PRECONDITIONER
C  FOR THE LINEAR SYSTEM H * X = B AND H IS THE HESSIAN OF A
C  GROUPS PARTIALLY SEPARABLE FUNCTION.
C
C  THIS IS ACHIEVED IN THREE STAGES.
C  1) ASSEMBLE THE MATRIX H.
C  2) FIND AND FACTORIZE THE POSITIVE DEFINITE MATRIX M.
C  3) SOLVE M * Q = GRAD.
C
C  NICK GOULD, MARCH 28TH 1989.
C  FOR CGT PRODUCTIONS.
C
C  ********************************************************************
C
      INTEGER           IFACTR, N , NG, MAXSEL, NVAR  , IOUT
      INTEGER           LSTADH, LICNA , LSTADA, LNTVAR, LELVAR, LELING
      INTEGER           LSTADG, LSTAEV, LNSTGV, LNVGRP, LIWK  , NFREEF
      INTEGER           LNGUVL, LNHUVL, LESCAL, LWK   , LGXEQX, LINTRE
      INTEGER           NFIXED, INFORM, IPRINT, LA   
      LOGICAL           MUNKS,  SEPREC, REFACT, BAND
      INTEGER           IVAR( N ), ISTADH( LSTADH ), ICNA( LICNA )  
      INTEGER           ISTADA( LSTADA ), INTVAR( LNTVAR )
      INTEGER           IELVAR( LELVAR ), IELING( LELING )
      INTEGER           ISTADG( LSTADG ), ISTAEV( LSTAEV )
      INTEGER           ISTAGV( LNSTGV ), ISVGRP( LNVGRP ), IWK( LIWK )
      LOGICAL           GXEQX( LGXEQX ), INTREP( LINTRE )
CS    REAL              A( LA ), GVALS2( NG ), GVALS3( NG ),
CD    DOUBLE PRECISION  A( LA ), GVALS2( NG ), GVALS3( NG ),
     *                  GUVALS( LNGUVL ), HUVALS( LNHUVL ), GSCALE( NG), 
     *                  ESCALE( LESCAL ), GRAD( N ), Q( N ), WK( LWK )
      EXTERNAL RANGES
C
C  ...................................................................
C
C  PARTITION OF INTEGER WORKSPACE ARRAY IWK
C  ----------------------------------------
C
C  FOR THE MUNKSGAARD OPTION:
C
C   <-                            LIWK                       ->
C  ------------------------------------------------------------
C  |<-NFREEF->|<-2*NNZH->|<-3*NNZH->|<-4*NFREEF->|<-4*NFREEF->|
C  ------------------------------------------------------------
C   |          |          |          |            | 
C  LIFREE     LIRNH      LJCNH      LIICCG       LIK
C
C  FOR THE MA27 OPTION:
C
C   <-                        LIWK                                 ->
C  -----------------------------------------------------------------------
C  |<-NFREEF->|<-NNZH->|<-NNZH->|<-3*NFREEF->|<-2*NFREEF->| <-NIRBDU-> | |
C  -----------------------------------------------------------------------
C   |          |        |        |            |            |            |
C  LIFREE     LIRNH    LJCNH    LIKEEP       LIW1         LIMA27
C
C  PARTITION OF REAL WORKSPACE ARRAY WK
C  ------------------------------------
C
C   <-                          LWK                                ->
C  -------------------------------------------------------------------
C  |   <- NRLBDU  - >  |               |<-NFREEF->|<-NFREEF->|<-NFREEF->|
C  -------------------------------------------------------------------
C   |                   |               |          |          |
C  LWMA27                              LW         LPERT      LRHS
C
C  ...................................................................
C
C  LOCAL VARIABLES.
C
      INTEGER           I,  J,  IFLAG,  LH, LW, LIKEEP, LIMA27, LIRNH
      INTEGER           LIWKUS, IFICCG, NSTEPS, LWICCG, LIRNBD, LP2
      INTEGER           LIWORK, LIW1,   LJCNH,  LNXTRW, LRHS,   LWKH
      INTEGER           LWKUSE, LIWKH,  NEG2,   LQ, LG, NSEMIW, MAXSBW
      INTEGER           LPERT,  MLP,    LWMA27, LWORK,  MAXFRT, NNZH
      INTEGER           LR,     LRHS2,  LRMAX,  NFMAX,  LIFREE, LIH
      INTEGER           LDIAG,  LIK,    LIICCG, LIWFAC, LICCGG, LWKFAC
      INTEGER           IAI,    IAJ,    NZ0,    NUPDAT, NUPMAX, NEXTRA
      INTEGER           IDUM1,  IDUM2,  IDUM3,  LHSCAL, LBAND,  NEG1
      INTEGER           LIVUSE, INFOAS, LQMAX,  LBD,    LIPBD,  LOFFDI
C ** Correction 2. 15/06/93: 1 lines corrected **
      INTEGER           LINXTR, NVNNZH , NVLH,  LENOFF
C ** Correction 2. 15/06/93:  end of correction **
      INTEGER           IPCD(   1 ),    JCNCD(  1 )
      LOGICAL           PRNTER
      REAL              TFACTR, T1STSL, TUPDAT, TSOLVE, CPUTIM, DUM
      REAL              TT, T
CS    REAL              OPS,    ZERO,   PERTUR, ONE,
CD    DOUBLE PRECISION  OPS,    ZERO,   PERTUR, ONE,
     *                  DIAMIN, DIAMAX, DIATOL, CD( 1 )
      INTRINSIC         ABS,    DBLE,   FLOAT,  INT,    MAX,    MIN
      INTRINSIC         SQRT
CS    EXTERNAL          SASMBL, MA27A , MA27C,  SMCFA,  SASSLB, SASSLC
CD    EXTERNAL          DASMBL, MA27AD, MA27CD, DMCFA,  DASSLB, DASSLC
CS    EXTERNAL          SICCGA, SICCGB, SSYPRC, MA27B,  CPUTIM
CD    EXTERNAL          DICCGA, DICCGB, DSYPRC, MA27BD, CPUTIM
CS    EXTERNAL          SBNDSL, SBNDFA
CD    EXTERNAL          DBNDSL, DBNDFA
C
C  COMMON VARIABLES.
C
      INTEGER           LMA271, LMA272, LMA273
CS    REAL              U
CD    DOUBLE PRECISION  U
CS    COMMON / MA27D  / U, LMA271, LMA272, LMA273
CD    COMMON / MA27DD / U, LMA271, LMA272, LMA273
      INTEGER           NRLNEC, NIRNEC, NRLBDU, NIRBDU
CS    COMMON / MA27E  / OPS,    IDUM1(  3 ),    NRLNEC, NIRNEC, 
CD    COMMON / MA27ED / OPS,    IDUM1(  3 ),    NRLNEC, NIRNEC, 
     *                  IDUM2(  2 ),    NRLBDU, NIRBDU, IDUM3(  5 )
      INTEGER           LICCG1, LICCG2
      REAL              DD
CS    COMMON / MA31I  / DD, LICCG1, LICCG2
CD    COMMON / MA31ID / DD, LICCG1, LICCG2
      INTEGER           IAJMAX, IAIMAX, IDUM4,  ND,     IDUM5
CS    COMMON / MA31J  / IAJMAX, IAIMAX, IDUM4,  ND,     IDUM5
CD    COMMON / MA31JD / IAJMAX, IAIMAX, IDUM4,  ND,     IDUM5
      INTEGER           ITERCG, ITCGMX, NGEVAL, ISKIP , IFIXED, NSEMIB
CS    REAL              ACCCG , RATIO , RADIUS, RADMAX, FINDMX, CICCG
CD    DOUBLE PRECISION  ACCCG , RATIO , RADIUS, RADMAX, FINDMX, CICCG
CS    COMMON / SCOMSB / ACCCG , RATIO , RADIUS, RADMAX, FINDMX, CICCG,
CD    COMMON / DCOMSB / ACCCG , RATIO , RADIUS, RADMAX, FINDMX, CICCG,
     *                  ITERCG, ITCGMX, NGEVAL, ISKIP , IFIXED, NSEMIB
CS    REAL              EPSMCH, EPSNEG, TINY,   BIG
CD    DOUBLE PRECISION  EPSMCH, EPSNEG, TINY,   BIG
CS    COMMON / SMACHN / EPSMCH, EPSNEG, TINY,   BIG
CD    COMMON / DMACHN / EPSMCH, EPSNEG, TINY,   BIG
      COMMON / PRECNN / PRNTER
C
C  SAVE INTERNAL VALUES BETWEEN CALLS.
C
      SAVE              LWKUSE, LWKH  , LIWKUS, LIWKH , LH    , LIRNH
      SAVE              LJCNH , LWORK , LWKFAC, LG, LW, LNXTRW, LOFFDI
      SAVE              IAI   , IAJ   , LIFREE, LIK   , LIICCG
      SAVE              LHSCAL, LDIAG , LRHS  , NZ0   , LICCGG, LIWFAC
      SAVE              LIWORK, LIKEEP, LIW1  , LIMA27, LWMA27, LP2
      SAVE              LPERT , NFMAX , LIVUSE, LIRNBD, LIPBD , LWICCG
      SAVE              LQMAX , LRMAX , LBD   , LQ, LR, LRHS2 , NSEMIW
      SAVE              TFACTR, T1STSL, TUPDAT, TSOLVE, LBAND , MAXFRT
C ** Correction 6. 15/06/93: 1 lines corrected **
      SAVE              NNZH  , LIH   , NUPMAX, NUPDAT, LENOFF
C ** Correction 6. 15/06/93:  end of correction **
C ** Correction 9. 13/01/94: 2 lines corrected **
CS    SAVE            / SCOMSB /, / SMACHN /
CD    SAVE            / DCOMSB /, / DMACHN /
C ** Correction 9. 13/01/94:  end of correction **
C
C  SET CONSTANT REAL PARAMETERS.
C
CS    PARAMETER       ( ZERO  = 0.0E+0, ONE   = 1.0E+0 )
CD    PARAMETER       ( ZERO  = 0.0D+0, ONE   = 1.0D+0 )
      PRNTER = IPRINT .GE. 5 .AND. IOUT .GT. 0
      IF ( IPRINT .GE. 1000 ) THEN
         LMA271 = 6
         LMA272 = 6
         LMA273 = 1
         LICCG1 = 6
         LICCG2 = 6
       END IF
C
C  ---------------------------------------------------------------------
C
C  STAGE 1A - FORM H. THIS STAGE NEEDS ONLY
C  BE PERFORMED WHEN THE INTEGER IFACTR IS 1.
C
C  ---------------------------------------------------------------------
C
  100 CONTINUE
      IF ( IFACTR .EQ. 1 ) THEN 
        IF ( IPRINT .GE. 200 .AND. IOUT .GT. 0 ) WRITE( IOUT, 2210 ) 
        T = CPUTIM( DUM )
        NUPDAT = 0
C
C  DEFINE THE REAL WORK SPACE NEEDED FOR ASEMBL. 
C  ENSURE THAT THERE IS SUFFICIENT SPACE.
C
         LWKUSE = N + 3 * MAXSEL
         LWKH   = LWK - LWKUSE
         IF ( LWKH .LE. 0 ) THEN
            INFORM = 10
            WRITE( IOUT, 2000 )
            RETURN
         END IF
C
C  DEFINE THE INTEGER WORK SPACE NEEDED FOR ASEMBL. 
C  ENSURE THAT THERE IS SUFFICIENT SPACE.
C
         LIWKUS = N
         LIWKH  = LIWK - LIWKUS
C ** Correction 10. 19/08/96: 5 lines added **
         IF ( LIWKH .LE. 0 ) THEN
            INFORM = 11
            WRITE( IOUT, 2010 )
            RETURN
         END IF
C ** Correction 10. 19/08/96: end of correction **
         IF ( BAND ) THEN
            LH     = LWKH
            LIH    = 1
            LINXTR = 1
         ELSE
            LH     = MIN( LWKH, ( LIWKH - 3 * NVAR ) / 4 )
            LIH    = LH
            LINXTR = LH + NVAR
         END IF
         IF ( LH .LE. 0 ) THEN
            INFORM = 11
            WRITE( IOUT, 2010 )
            RETURN
         END IF
C
C  SET STARTING ADDRESSES FOR PARTITIONS OF THE INTEGER WORKSPACE.
C
         LIFREE = 1
         LIRNH  = LIFREE + NVAR
         LJCNH  = LIRNH  + LIH
         LNXTRW = LJCNH  + LIH
C
C  ASSEMBLE THE HESSIAN RESTRICTED TO THE VARIABLES IVAR( I ), I = 1,..,
C  NVAR. REMOVE THE NONZEROS WHICH LIE OUTSIDE A BAND WITH 
C  SEMI-BANDWIDTH NSEMIB.
C
         TT = CPUTIM( DUM )
         IF ( BAND ) THEN
            NSEMIW = MIN( NSEMIB, MAX( NVAR - 1, 0 ) )
C ** Correction 3. 15/06/93: 1 line added **
            LENOFF = MAX( 1, NSEMIW )
C ** Correction 3. 15/06/93:  end of correction **
         ELSE
            NSEMIW = NSEMIB
         END IF
C ** Correction 11. 05/11/97: 1 line corrected **
C ** Correction 12. 19/08/98: added final argument to call to ASMBL.
CS       CALL SASMBL( N , NG, MAXSEL, NSEMIW, LH    , LIH   , NNZH  , 
CD       CALL DASMBL( N , NG, MAXSEL, NSEMIW, LH    , LIH   , NNZH  , 
     *                NVAR  , IVAR  , ISTADH, LSTADH, ICNA  , LICNA , 
     *                ISTADA, LSTADA, INTVAR, LNTVAR, IELVAR, LELVAR,
     *                IELING, LELING, ISTADG, LSTADG, ISTAEV, LSTAEV,
     *                ISTAGV, LNSTGV, ISVGRP, LNVGRP, IWK  ( LIRNH ),
     *                IWK ( LJCNH ) , IWK ( LNXTRW ), LINXTR,
     *                IWK( LIWKH + 1 )      , LIWK - LIWKH  ,
     *                A , LA, GUVALS, LNGUVL, HUVALS, 
     *                LNHUVL, GVALS2, GVALS3, GSCALE, ESCALE, LESCAL, 
     *                WK    , WK (LWKH + 1 ), LWK - LWKH    , GXEQX , 
     *                LGXEQX, INTREP, LINTRE, RANGES, IPRINT, IOUT    ,
     *                BAND  , MAXSBW, INFORM, .TRUE., .FALSE. )
C ** Correction 12. 19/08/98: end of correction **
C ** Correction 11. 05/11/97: end of correction **
         IF ( IPRINT .GE. 200 .AND. IOUT .GT. 0 ) 
     *       WRITE( IOUT, 2220 ) CPUTIM( DUM ) - TT
C
C  CHECK THAT THERE IS SUFFICIENT INTEGER WORKSPACE.
C
         IF ( INFORM .GT. 0 ) THEN
            INFORM = 11
            WRITE( IOUT, 2010 )
            RETURN
         END IF
C
C  ------------------------------------------------------
C
C  STAGE 2 - FORM AND FACTORIZE M. THIS STAGE NEEDS ONLY
C  BE PERFORMED WHEN THE INTEGER IFACTR IS 1.
C
C  - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C
C  MUNKSGAARD'S PRECONDITIONER, ICCG, IS TO BE USED.
C
C  - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
         IF ( MUNKS ) THEN
C
C  DECIDE HOW MUCH ROOM IS AVAILABLE FOR THE INCOMPLETE FACTORIZATION.
C  NEXTRA GIVES THE AMOUNT OF WORKSPACE ABOVE THE MINIMUM REQUIRED
C  FOR THE FACTORIZATION WHICH IS TO BE ALLOWED FOR FILL-IN.
C
            NEXTRA = NNZH
            IF ( LIWK - 9 * NVAR - 3 * NNZH - 2 * NEXTRA .LT. 0 ) THEN
               IF ( LIWK - 9 * NVAR - 3 * NNZH .LT. 0 ) THEN
                  INFORM = 11
                  WRITE( IOUT, 2010 )
                  RETURN
               ELSE
                  IAI = NNZH
                  IAJ = 2 * NNZH
               END IF
            ELSE
               IAI = NNZH + NEXTRA
               IAJ = 2 * NNZH + NEXTRA
            END IF
            NVNNZH = NVAR + IAI
            NVLH   = NVAR + LH
C
C  COMPRESS THE VECTOR IWK TO REMOVE ZEROS CREATED IN THE ASSEMBLY.
C
            IF ( LH .GE. IAI ) THEN
               DO 110 I = 1, NNZH
                  IWK( NVNNZH + I ) = IWK( NVLH + I )
  110          CONTINUE
            ELSE
               DO 120 I = NNZH, 1, -1
                  IWK( NVNNZH + I ) = IWK( NVLH + I )
  120          CONTINUE
            END IF
            LJCNH = LIRNH + IAI
C
C  DEFINE THE INTEGER WORK SPACE FOR ICCGA.
C
            LIK    = LIWK - 4 * NVAR + 1
            LIICCG = LIK  - 4 * NVAR
C
C  DEFINE THE REAL WORK SPACE FOR ICCGA.
C
            LWORK = LWK - 3 * NVAR
            IF ( LWORK .LE. 3 * NNZH ) THEN
               INFORM = 10
               WRITE( IOUT, 2000 )
               RETURN
            END IF
            LWICCG = 1
            LHSCAL = LWICCG + LWORK
            LDIAG  = LHSCAL + NVAR
            LRHS   = LDIAG  + NVAR
C
C  PRINT DETAILS OF THE PARTITIONS.
C
C           WRITE( 6, 1100 ) NVAR, NNZH, IAI, IAJ, CICCG
C           DO 99 I = 1, NNZH
C              WRITE( 6, 1120 ) I, IWK( LIRNH+I-1), IWK( LJCNH+I-1), 
C    *                          WK( I )
C  99       CONTINUE
C           DO 97 I = 1, NVAR
C              WRITE( 6, 1130 ) I, GRAD( I )
C  97       CONTINUE
C
C  FORM AND FACTORIZE MUNKSGAARD'S PRECONDITIONER.
C
CS          CALL SICCGA( NVAR, NNZH, WK( LWICCG ), IWK( LIRNH ),
CD          CALL DICCGA( NVAR, NNZH, WK( LWICCG ), IWK( LIRNH ),
     *                   IWK( LJCNH ), IAI, IAJ, IWK( LIK ),
     *                   IWK( LIICCG ), WK( LHSCAL ), CICCG, IFICCG )
            IF ( PRNTER .OR. ( IOUT .GT. 0 .AND. IPRINT .EQ. 2 ) ) THEN
               IF ( PRNTER ) THEN
                  WRITE( IOUT, 2140 ) NVAR, NNZH, IAJMAX
               ELSE
                  WRITE( IOUT, 2150 ) NVAR, NNZH, IAJMAX
               END IF
            END IF
            IF ( ( IFICCG .LT. 0 .AND. IOUT .GT. 0 ) .OR.
     *           ( IFICCG .GT. 0 .AND. PRNTER ) )
     *           WRITE( IOUT, 2160 ) IFICCG
C
C  COMPRESS THE VECTOR IWK TO REMOVE UNUSED LOCATIONS.
C
            NZ0    = NNZH   - ND
            LICCGG = IAJ    - NZ0
            LIWFAC = LJCNH  + NZ0 - 1
            DO 140 I           = 1, LICCGG
               IWK( NVAR + I ) = IWK( LIWFAC + I )
  140       CONTINUE
C
C  STORE STARTING ADDRESSES FOR THE FACTORS OF THE PRECONDITIONER.
C
            LIWFAC = NVAR   + 1
            LWKFAC = LWICCG + NZ0
            LG     = LWICCG + IAJ
            LIWORK = LIWK - LICCGG - 4 * NVAR
C
C  RECORD THE RELATIVE FILL-IN.
C
            IF ( NNZH .GT. 0 ) THEN  
CS             RATIO =       FLOAT( LICCGG )   /       FLOAT( NNZH )  
CD             RATIO = DBLE( FLOAT( LICCGG ) ) / DBLE( FLOAT( NNZH ) )
            ELSE
               RATIO = ONE
            END IF 
C 
C  - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  
C  A BAND PRECONDITIONER IS TO BE USED.
C
C  - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
         ELSE
            IF ( BAND ) THEN
C
C  CHECK THAT THERE WAS SUFFICIENT REAL WORKSPACE FOR THE FACTORIZATION.
C
               LWORK = LWK - NVAR * ( NSEMIW + 2 ) 
               IF ( LWORK .LE. 0 ) THEN
                  INFORM = 10
                  WRITE( IOUT, 2000 )
                  RETURN
               END IF
C
C  DEFINE THE REAL WORK SPACE FOR THE BAND SOLVER.
C
               LDIAG  = 1
               LOFFDI = LDIAG  + NVAR
               LG     = LOFFDI + NSEMIW * NVAR
               LRHS   = LG     + LWORK
C
C  DEFINE THE INTEGER WORK SPACE FOR THE BAND SOLVER.
C
               LIWORK = LIWK - NVAR
               LBAND  = NVAR
C
C  FACTORIZE THE BAND MATRIX.
C
C ** Correction 4. 15/06/93: 1 line added **
CS             CALL SBNDFA( NVAR, NSEMIW, WK( LDIAG ), WK( LOFFDI ),
CD             CALL DBNDFA( NVAR, NSEMIW, WK( LDIAG ), WK( LOFFDI ),
     *                      LENOFF )
C ** Correction 4. 15/06/93:  end of correction **
               IF ( PRNTER .OR. IPRINT .EQ. 2 ) THEN
                  IF ( PRNTER ) THEN
                     WRITE( IOUT, 2120 ) NVAR, NSEMIW, MAXSBW
                  ELSE
                     WRITE( IOUT, 2130 ) NVAR, NSEMIW, MAXSBW
                  END IF
               END IF
            ELSE
C 
C  - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  
C  A MULTI-FRONTAL PRECONDITIONER, MA27, IS TO BE USED.
C
C  - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
C  COMPRESS THE VECTOR IWK TO REMOVE ZEROS CREATED IN THE ASSEMBLY.
C
               NVNNZH = NVAR + NNZH
               NVLH   = NVAR + LH
               DO 150 I = 1, NNZH
                  IWK( NVNNZH + I ) = IWK( NVLH + I )
  150          CONTINUE
               LJCNH = LIRNH + NNZH
C
C  DEFINE THE STARTING ADDRESSES FOR INTEGER WORK SPACE FOR MA27A.
C
               LIWORK = LIWK - 2 * NNZH - 6 * NVAR
               IF ( LIWORK .LE. 0 ) THEN
                  INFORM = 11
                  WRITE( IOUT, 2010 )
                  RETURN
               END IF
               LIKEEP = LJCNH  + NNZH
               LIW1   = LIKEEP + 3 * NVAR
               LIMA27 = LIW1   + 2 * NVAR
C
C  CHOOSE THE PIVOT SEQUENCE FOR THE FACTORIZATION BY ANALYZING THE
C  SPARSITY PATTERN OF M.
C
               IFLAG = 0
               TT = CPUTIM( DUM )
CS             CALL MA27A ( NVAR, NNZH, IWK( LIRNH ), IWK( LJCNH ),
CD             CALL MA27AD( NVAR, NNZH, IWK( LIRNH ), IWK( LJCNH ),
     *                      IWK( LIMA27 ), LIWORK, IWK( LIKEEP ),
     *                      IWK( LIW1 ), NSTEPS, IFLAG )
               IF ( IPRINT .GE. 200 .AND. IOUT .GT. 0 ) 
     *              WRITE( IOUT, 2230 ) CPUTIM(DUM) - TT
C
C  CHECK THAT THERE WAS SUFFICIENT INTEGER WORKSPACE FOR THE ANALYSIS.
C
               IF ( IFLAG .LT. 0 ) THEN
                  WRITE( IOUT, 2020 ) IFLAG
                  IF ( IFLAG .EQ. - 3 ) THEN
                     INFORM = 11
                     WRITE( IOUT, 2010 )
                  ELSE
C
C  CHECK THAT THERE WAS SUFFICIENT REAL WORKSPACE FOR THE ANALYSIS.
C
                     INFORM = 10
                     WRITE( IOUT, 2000 )
                  END IF
                  RETURN
               END IF
C
C  CHECK THAT THERE IS SUFFICIENT REAL WORKSPACE FOR THE FACTORIZATION.
C
               LWORK = LWK - 3 * NVAR
               IF ( LWORK .LE. NRLNEC ) THEN
                  INFORM = 10
                  WRITE( IOUT, 2000 )
                  RETURN
               END IF
C
C  DEFINE THE WORK SPACE FOR THE FACTORIZATION, MCFA, 
C  AND THE SOLVER, MA27C, SUBPROGRAMS.
C
               LWMA27 = 1
               LW     = LWMA27 + LWORK
               LPERT  = LW     + NVAR
               LRHS   = LPERT  + NVAR
C
C  FACTORIZE THE MATRIX M.
C
               MLP = 0
               IF ( IPRINT .GE. 100 ) MLP = IOUT
               IF ( SEPREC ) THEN
C
C  USE THE SCHNABEL-ESKOW MODIFIED CHOLESKY FACTORIZATION.
C
CS                CALL SMCFA( NVAR, NNZH, IWK( LIRNH ), IWK( LJCNH ),
CD                CALL DMCFA( NVAR, NNZH, IWK( LIRNH ), IWK( LJCNH ),
     *                      WK( LWMA27 ), LWORK, IWK( LIMA27 ), LIWORK,
     *                      IWK( LIKEEP ), NSTEPS, MAXFRT, IWK( LIW1 ),
     *                      IFLAG, WK( LRHS ), WK( LPERT ), MLP )
C
C  CALCULATE THE MAXIMUM AND MINIMUM DIAGONALS IN THE FACTORIZATION.
C
                  DIAMAX    = ZERO
                  DIAMIN    = BIG
                  DO 160 J  = 1, NVAR
                     DIAMIN = MIN( DIAMIN, WK( LRHS + J - 1 ) )
                     DIAMAX = MAX( DIAMAX, WK( LRHS + J - 1 ) )
  160             CONTINUE
C
C  ENSURE THAT NO DIAGONAL IS SMALLER THAN 10**-8 TIMES THE LARGEST.
C
                  DIATOL                = 1.0D-8 * DIAMAX
                  DO 170 J              = 1, NVAR
                     WK( LRHS + J - 1 ) = MAX( DIATOL,
     *                                         WK( LRHS + J - 1 ) )
  170             CONTINUE
C
C  CALCULATE THE MAXIMUM PERTURBATION MADE TO THE DIAGONALS OF H.
C
                  IF ( PRNTER .OR. IPRINT .EQ. 2 ) THEN
                     PERTUR    = ZERO
                     DO 180 J  = 1, NVAR
                        PERTUR = MAX( PERTUR,
     *                                ABS( WK( LPERT + J - 1 ) ) )
                        IF ( PRNTER ) THEN
                           IF ( WK( LPERT + J - 1 ) .NE. ZERO ) 
     *                        WRITE( IOUT, 2040)
     *                               WK( LPERT + J - 1 ), IVAR( J )
                        END IF
  180                CONTINUE
                     IF ( PRNTER ) THEN
                        WRITE( IOUT, 2060 ) PERTUR, NVAR, NNZH, 
     *                         IDUM2( 1 ), NRLBDU, DIAMIN, DIAMAX
                     ELSE
                        WRITE( IOUT, 2070 ) NVAR, NNZH, PERTUR, 
     *                         IDUM2( 1 ), NRLBDU, DIAMIN, DIAMAX
                     END IF
                  END IF
               ELSE
C
C  USE THE GILL-MURRAY-PONCELEON-SAUNDERS MODIFICATION TO THE
C  SYMMETRIC INDEFINITE FACTORIZATION.
C
                  TT = CPUTIM( DUM )
CS                CALL MA27B ( NVAR, NNZH, IWK( LIRNH ), IWK( LJCNH ),
CD                CALL MA27BD( NVAR, NNZH, IWK( LIRNH ), IWK( LJCNH ),
     *                      WK( LWMA27 ), LWORK, IWK( LIMA27 ), LIWORK,
     *                      IWK( LIKEEP ), NSTEPS, MAXFRT, IWK( LIW1 ),
     *                      IFLAG )
                  IF ( IPRINT .GE. 200 .AND. IOUT .GT. 0 ) 
     *            WRITE( IOUT, 2240 ) CPUTIM( DUM ) - TT
                  TT = CPUTIM( DUM )
CS                CALL SSYPRC( NVAR, LWORK, LIWORK, WK( LWMA27 ),
CD                CALL DSYPRC( NVAR, LWORK, LIWORK, WK( LWMA27 ),
     *                         IWK( LIMA27 ), NEG1, NEG2 )
                  IF ( IPRINT .GE. 200 .AND. IOUT .GT. 0 ) 
     *            WRITE( IOUT, 2250 ) CPUTIM( DUM ) - TT
                  IF ( PRNTER .OR. IPRINT .EQ. 2 ) THEN
                     IF ( PRNTER ) THEN
                        WRITE( IOUT, 2090 ) NVAR, NNZH, IDUM2( 1 ),
     *                                    NRLBDU, NEG1, NEG2
                     ELSE
                        WRITE( IOUT, 2100 ) NVAR, NNZH, IDUM2( 1 ),
     *                                    NRLBDU, NEG1, NEG2
                     END IF
                  END IF
               END IF   
C
C  RECORD THE RELATIVE FILL-IN.
C
               IF ( NNZH .GT. 0 ) THEN  
C ** Correction 1. 23/03/93: 2 lines corrected **
CS                RATIO =       FLOAT( NRLBDU )   /       FLOAT( NNZH )  
CD                RATIO = DBLE( FLOAT( NRLBDU ) ) / DBLE( FLOAT( NNZH ))
C ** Correction 1. 23/03/93: end of correction **
               ELSE
                  RATIO = ONE
               END IF 
C
C  CHECK THAT THERE WAS SUFFICIENT WORKSPACE FOR THE FACTORIZATION.
C
               IF ( IFLAG .NE. 0 .AND. IFLAG .NE. 3 ) THEN
                  WRITE( IOUT, 2030 ) IFLAG
                  IF ( IFLAG .EQ. - 3 ) THEN
                     INFORM = 11
                     WRITE( IOUT, 2010 )
                  ELSE
                     INFORM = 10
                     WRITE( IOUT, 2000 )
                  END IF
                  RETURN
               END IF
C
C  CHECK THAT THERE WAS SUFFICIENT REAL WORKSPACE FOR THE SOLVES.
C
               IF ( LWORK - N .LE. 0 ) THEN
                  INFORM = 12
                  WRITE( IOUT, 2020 )
                  RETURN
               END IF
C
C  SET THE STARTING ADDRESSES FOR FURTHER PARTITIONS OF THE 
C  REAL WORKSPACE.
C
               LG = LWMA27 + NRLBDU
            END IF
         END IF
C
C  - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  
C  THE FACTORIZATION IS COMPLETE.
C  
C  - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
C  STORE THE NUMBER OF FREE VARIABLES, NFREEF, AND LIST OF FREE
C  VARIABLES WHEN THE FACTORIZATION IS PERFORMED IN CASE THEY
C  THEY ARE NEEDED ON A SUBSEQUENT ENTRY.
C
         NFIXED                   = 0
         NFREEF                   = NVAR
         DO 190 J                 = 1, NVAR
            IWK( LIFREE + J - 1 ) = IVAR( J )
  190    CONTINUE
C
C  RECORD THE TIME TAKEN TO ASSEMBLE AND FACTORIZE THE PRECONDITIONER.
C
         TFACTR = CPUTIM( DUM ) - T
         TUPDAT = 0.0
         TSOLVE = 0.0
         NUPMAX = 100
      END IF
C
C  ---------------------------------------------------------------------
C
C  STAGE 2 (ALTERNATIVE) - UPDATE THE SCHUR COMPLEMENT OF M. THIS STAGE
C  NEEDS ONLY BE PERFORMED WHEN THE INTEGER IFACTR IS 2.
C
C  ---------------------------------------------------------------------
C
      IF ( IFACTR .EQ. 2 ) THEN
        IF ( IPRINT .GE. 200 .AND. IOUT .GT. 0 ) WRITE( IOUT, 2260 ) 
C
C  REFACTORIZE THE REMAINING COEFFICIENT MATRIX IF THE LAST 
C  UPDATE AND SOLVE TOOK LONGER THAN THE PREVIOUS REFACTORIZATION AND
C  SUBSEQUENT SOLVE.
C
         NUPDAT = NUPDAT + 1
         IF ( NUPDAT .GT. NUPMAX ) THEN
C        IF ( NUPDAT .GT. NUPMAX .OR. TUPDAT + TSOLVE .GT.
C    *        TFACTR + T1STSL ) THEN
            IFACTR = 1
            IF ( IPRINT .GE. 2 .AND. IOUT .GT. 0 ) WRITE( IOUT, 2170 )
     *        NUPDAT, NUPMAX, TUPDAT + TSOLVE, TFACTR + T1STSL
            REFACT = .TRUE.
            GO TO 100
         END IF
C
C  SET UP FURTHER WORKSPACE ADDRESSES TO PARTITION THE UNUSED PORTIONS
C  OF WK AND IWK. CALCULATE THE LARGEST NUMBER OF VARIABLES, NFMAX,
C   WHICH CAN BE FIXED. 
C
         T = CPUTIM( DUM )
         IF ( NFIXED .EQ. 0 ) THEN
            IF ( MUNKS ) THEN
               NFMAX = MIN( NFREEF, ( LIWORK - N - 2 ) / 2,
     *                      INT( SQRT( FLOAT( 2 * ( LWORK - IAJ ) -
     *                      3 * N - NFREEF ) ) ) - 5 )
               LIVUSE = LICCGG + 1
            ELSE
               IF ( BAND ) THEN
                  NFMAX = MIN( NFREEF, ( LIWORK - N - 2 ) / 2,
     *                         INT( SQRT( FLOAT( 2 * ( LWORK - N )
     *                         - 3 * ( N + NFREEF ) ) ) ) - 5 )
                  LIVUSE = LBAND + 1
               ELSE
                  NFMAX = MIN( NFREEF, ( LIWORK - NIRBDU - N - 2 ) / 2,
     *                         INT( SQRT( FLOAT( 2 * ( LWORK - NRLBDU )
     *                         - 3 * ( N + NFREEF ) ) ) ) - 5 )
                  LIVUSE = LIMA27 + NIRBDU
               END IF
            END IF
C
C  IF THERE IS INSUFFICIENT ROOM, THE MATRIX WILL HAVE TO BE REFACTORIZED
C  EVERY TIME A VARIABLE IS FIXED.
C
            IF ( NFMAX .LE. 0 ) THEN
               IFACTR = 1
               REFACT = .TRUE.
               GO TO 100
            END IF
C
C  SET INTEGER WORKSPACE ADDRESSES FOR THE SOLVE.
C
            LIRNBD = LIVUSE + N
            LIPBD  = LIRNBD + NFMAX
            IWK( LIPBD ) = 1
C
C  SET REAL WORKSPACE ADDRESSES FOR THE SOLVE.
C
            LQMAX  = NFMAX
            LRMAX  = NFMAX * ( NFMAX + 1 ) / 2
            LBD    = LG     + N
            LQ     = LBD    + NFMAX
            LR     = LQ     + LQMAX
            LRHS2  = LR     + LRMAX
            LP2    = LRHS2  + NFREEF + NFMAX
         END IF
C
C  RECORD THE VARIABLES WHICH ARE STILL FREE BY FIRST SETTING
C  THE APPROPRIATE PARTITION OF IWK TO ZERO.
C
         DO 210 I = 1, N
            IWK( LIVUSE + I - 1 ) = 0
  210    CONTINUE
C
C  NOW SET THE ENTRIES IN IWK CORRESPONDING TO THE FREE VARIABLES
C  TO ONE.
C
         DO 220 J = 1, NVAR
            I     = IVAR( J )
            IWK( LIVUSE + I - 1 ) = 1
  220    CONTINUE
C
C  COMPARE THIS WITH THE LIST OF THOSE WHICH WERE FREE AT THE LAST
C  FACTORIZATION.  FIX ANY VARIABLES WHICH WAS FREE BUT NO LONGER
C  APPEARS IN THE LIST.
C
         DO 240 J = 1, NFREEF
            I     = IWK( LIFREE + J - 1 )
            IF ( I .GT. 0 ) THEN
               IF ( IWK( LIVUSE + I - 1 ) .EQ. 0 ) THEN
C
C  IF MORE THAN NFMAX VARIABLES HAVE BEEN FIXED, REFACTORIZE THE
C  MATRIX
C
                  IF ( NFIXED .GE. NFMAX ) THEN
                     IFACTR = 1
                     REFACT = .TRUE.
                     GO TO 100
                  END IF
C
C  UPDATE THE FACTORIZATION OF THE SCHUR COMPLEMENT TO ALLOW FOR
C  THE REMOVAL OF THE J-TH ROW AND COLUMN OF THE ORIGINAL HESSIAN -
C  THIS REMOVAL IS EFFECTED BY APPENDING THE J-TH ROW AND COLUMN
C  OF THE IDENTITY MATRIX TO THE HESSIAN.
C
                  WK( LBD + NFIXED )        = ONE
                  IWK( LIRNBD + NFIXED )    = J
                  IWK( LIPBD + NFIXED + 1 ) = NFIXED + 2
                  INFOAS = 1
C                 DUM = FLOAT( NVAR - IWK( LIKEEP + NFIXED - 1 ) + 1 )
C     *                     / FLOAT( NVAR )
C                 IF ( IPRINT .GE. 200 .AND. IOUT .GT. 0 ) 
C     *              WRITE( IOUT, 2270 )
C     *              ' SPARSE RHS, PERCENTAGE SKIPPED ', DUM
  230             CONTINUE
C
C  CALL ASSLC TO UPDATE THE SCHUR-COMPLEMENT.
C
                  TT = CPUTIM( DUM )
CS                CALL SASSLC( NFREEF, NFIXED, 4, NFMAX, WK( LBD ),
CD                CALL DASSLC( NFREEF, NFIXED, 4, NFMAX, WK( LBD ),
     *                         IWK( LIRNBD ), IWK( LIPBD ), 1, CD,
     *                         JCNCD, IPCD, WK( LR ), LRMAX, WK( LQ ),
     *                         LQMAX, WK( LRHS2 ), WK( LRHS ), INFOAS )
                  IF ( IPRINT .GE. 200 .AND. IOUT .GT. 0 ) 
     *               WRITE( IOUT, 2270 ) CPUTIM(DUM) - TT
                  IF ( INFOAS .GT. 0 ) THEN
C
C  ASSL REQUIRES ADDITIONAL INFORMATION. COMPUTE THE SOLUTION TO
C  THE EQUATION H * S = RHS, RETURNING THE SOLUTION S IN RHS.
C
C  FOR MUNSKGAARD'S FACTORIZATION.
C
                     IF ( MUNKS ) THEN
CS                      CALL SICCGB( NFREEF, WK( LWKFAC ), IWK( LIWFAC),
CD                      CALL DICCGB( NFREEF, WK( LWKFAC ), IWK( LIWFAC),
     *                               LICCGG, WK( LDIAG ), WK( LHSCAL ),
     *                               IWK( LIK ), WK( LRHS ) )
                     ELSE
C
C  FOR THE BAND FACTORIZATION.
C
                        IF ( BAND ) THEN
C ** Correction 5. 15/06/93: 1 line corrected **
CS                         CALL SBNDSL( NFREEF, NSEMIW, WK( LDIAG ),
CD                         CALL DBNDSL( NFREEF, NSEMIW, WK( LDIAG ),
     *                              WK( LOFFDI ), LENOFF, WK( LRHS ) )
C ** Correction 5. 15/06/93:  end of correction **
                        ELSE   
C
C  FOR THE MULTIFRONTAL FACTORIZATION.
C
                           TT = CPUTIM( DUM )
CS                         CALL MA27C ( NFREEF, WK( LWMA27 ), NRLBDU,
CD                         CALL MA27CD( NFREEF, WK( LWMA27 ), NRLBDU,
     *                         IWK( LIMA27 ), NIRBDU, WK( LW ), MAXFRT,
     *                               WK( LRHS ), IWK( LIW1 ), NSTEPS )
                           IF ( IPRINT .GE. 200 .AND. IOUT .GT. 0 ) 
     *                     WRITE( IOUT, 2280 ) CPUTIM( DUM ) - TT
                        END IF
                     END IF
                     GO TO 230
                  END IF
C
C  IF THE SCHUR-COMPLEMENT IS NUMERICALLY INDEFINITE, REFACTORIZE
C  THE PRECONDITIONING MATRIX TO ALLEVIATE THE EFFECT OF ROUNDING.
C
                  IF ( INFOAS .LT. 0 ) THEN
                     WRITE( IOUT, 2050 ) INFOAS
C                    IF ( PRNTER ) WRITE( IOUT, 2050 ) INFOAS
                     IFACTR = 1
                     REFACT = .TRUE.
                     GO TO 100
                  END IF
C                 WRITE( 6, 2080 ) TFACTR, TMA27C, TASSLB, TASSLC
C
C  RECORD THAT THE RELEVANT VARIABLE IS NOW FIXED.
C
                  IWK( LIFREE + J - 1 ) = - I
                  NFIXED = NFIXED + 1
               END IF
            END IF
  240    CONTINUE
         TUPDAT = CPUTIM( DUM ) - T
         TSOLVE = 0.0
      END IF
C
C  -------------------------------------------------
C
C  STAGE 3 - SOLVE FOR THE PRECONDITIONED GRADIENT.
C
C  - - - - - - - - - - - - - - - - - - - - - - - - -
C
C  INITIAL SOLVE USING THE ORIGINAL FACTORIZATION.
C
C  - - - - - - - - - - - - - - - - - - - - - - - - -
C
C  PUT THE COMPONENTS OF GRAD INTO WK(LRHS).
C
      IF ( IPRINT .GE. 200 .AND. IOUT .GT. 0 ) WRITE( IOUT, 2200 ) 
      T = CPUTIM( DUM )
      IF ( NFIXED .EQ. 0 ) THEN
         DO 310 J = 1, NVAR
            WK( LRHS + J - 1 ) = GRAD( J )
  310    CONTINUE
C
C  COMPUTE THE SOLUTION TO THE EQUATION M * S = RHS,
C  RETURNING THE SOLUTION S IN RHS.
C
C  USING MUNSKGAARD'S FACTORIZATION.
C
         IF ( MUNKS ) THEN
CS          CALL SICCGB( NFREEF, WK( LWKFAC ), IWK( LIWFAC ),
CD          CALL DICCGB( NFREEF, WK( LWKFAC ), IWK( LIWFAC ),
     *                   LICCGG, WK( LDIAG ), WK( LHSCAL ),
     *                   IWK( LIK ), WK( LRHS ) )
         ELSE
C
C  FOR THE BAND FACTORIZATION.
C
            IF ( BAND ) THEN
C ** Correction 7. 15/06/93: 1 line corrected **
CS             CALL SBNDSL( NFREEF, NSEMIW, WK( LDIAG ),
CD             CALL DBNDSL( NFREEF, NSEMIW, WK( LDIAG ),
     *                      WK( LOFFDI ), LENOFF, WK( LRHS ) )
C ** Correction 7. 15/06/93:  end of correction **
            ELSE   
C
C  USING THE MULTIFRONTAL FACTORIZATION.
C
            TT = CPUTIM( DUM )
CS             CALL MA27C ( NFREEF, WK( LWMA27 ), NRLBDU, IWK( LIMA27 ),
CD             CALL MA27CD( NFREEF, WK( LWMA27 ), NRLBDU, IWK( LIMA27 ),
     *                      NIRBDU, WK( LW ), MAXFRT, WK( LRHS ),
     *                      IWK( LIW1 ), NSTEPS )
            IF ( IPRINT .GE. 200 .AND. IOUT .GT. 0 ) 
     *          WRITE( IOUT, 2280 ) CPUTIM( DUM ) - TT
            END IF
         END IF
C
C  SCATTER THE FREE COMPONENTS OF THE SOLUTION INTO Q.
C
         DO 320 J = 1, NVAR
            Q( IVAR( J ) ) = WK( LRHS + J - 1 )
  320    CONTINUE
         T1STSL = CPUTIM( DUM ) - T
      ELSE
C
C  - - - - - - - - - - - - - - - - - - - - - - - - -
C
C  SUBSEQUENT SOLVES USING THE ORIGINAL FACTORIZATION
C  AND THE FACTORIZATION OF THE SCHUR-COMPLEMENT.
C
C  - - - - - - - - - - - - - - - - - - - - - - - - -
C
C  SOLVE FOR THE PRECONDITIONED GRADIENT USING THE SCHUR COMPLEMENT
C  UPDATE. PUT THE COMPONENTS OF GRAD INTO WK(LRHS2).
C
         DO 330 J = 1, NVAR
            WK( LG + IVAR( J ) - 1 ) = GRAD( J )
  330    CONTINUE
         DO 340 J = 1, NFREEF
            I     = IWK( LIFREE + J - 1 )
            IF ( I .GT. 0 ) THEN
               WK( LRHS2 + J - 1 ) = WK( LG + I - 1 )
            ELSE
               WK( LRHS2 + J - 1 ) = ZERO
            END IF
  340    CONTINUE
         DO 350 J = 1, NFIXED
            WK( LRHS2 + NFREEF + J - 1 ) = ZERO
  350    CONTINUE
C
C  SOLVE THE LINEAR SYSTEM H * P2 = RHS2.
C
         INFOAS = 1
  360    CONTINUE
C
C  CALL ASSLB TO SOLVE THE SYSTEM.
C
         TT = CPUTIM( DUM )
CS       CALL SASSLB( NFREEF, NFIXED, 4, NFMAX, WK( LBD ),
CD       CALL DASSLB( NFREEF, NFIXED, 4, NFMAX, WK( LBD ),
     *                IWK( LIRNBD ), IWK( LIPBD ), 1, CD, JCNCD,
     *                IPCD, WK( LR ), LRMAX, WK( LQ ), LQMAX,
     *                WK( LRHS2 ), WK( LP2 ), WK( LRHS ), INFOAS )
         IF ( IPRINT .GE. 200 .AND. IOUT .GT. 0 ) 
     *        WRITE( IOUT, 2290 ) CPUTIM( DUM ) - TT
         IF ( INFOAS .GT. 0 ) THEN
C
C  ASSL REQUIRES ADDITIONAL INFORMATION. COMPUTE THE SOLUTION TO
C  THE EQUATION H * S = RHS, RETURNING THE SOLUTION S IN RHS.
C
C  USING MUNSKGAARD'S FACTORIZATION.
C
            IF ( MUNKS ) THEN
CS             CALL SICCGB( NFREEF, WK( LWKFAC ), IWK( LIWFAC ),
CD             CALL DICCGB( NFREEF, WK( LWKFAC ), IWK( LIWFAC ),
     *                      LICCGG, WK( LDIAG ), WK( LHSCAL ),
     *                      IWK( LIK ), WK( LRHS ) )
            ELSE
C
C  FOR THE BAND FACTORIZATION.
C
               IF ( BAND ) THEN
C ** Correction 8. 15/06/93: 1 line corrected **
CS                CALL SBNDSL( NFREEF, NSEMIW, WK( LDIAG ),
CD                CALL DBNDSL( NFREEF, NSEMIW, WK( LDIAG ),
     *                      WK( LOFFDI ), LENOFF, WK( LRHS ) )
C ** Correction 8. 15/06/93:  end of correction **
               ELSE   
C
C  USING THE MULTIFRONTAL FACTORIZATION.
C
                  TT = CPUTIM( DUM )
CS                CALL MA27C ( NFREEF, WK( LWMA27 ), NRLBDU,
CD                CALL MA27CD( NFREEF, WK( LWMA27 ), NRLBDU, 
     *                         IWK( LIMA27 ), NIRBDU, WK( LW ), 
     *                         MAXFRT, WK( LRHS ),IWK( LIW1 ), NSTEPS )
                  IF ( IPRINT .GE. 200 .AND. IOUT .GT. 0 ) 
     *            WRITE( IOUT, 2280 ) CPUTIM( DUM ) - TT
               END IF
            END IF
            GO TO 360
         END IF
C
C  SCATTER THE FREE COMPONENTS OF THE SOLUTION INTO Q.
C
         DO 370 J = 1, NFREEF
            I     = IWK( LIFREE + J - 1 )
            IF ( I .GT. 0 ) Q( I ) = WK( LP2 + J - 1 )
  370    CONTINUE
         TSOLVE = CPUTIM( DUM ) - T
         IF( IPRINT .GE. 10 .AND. IOUT. GT. 0 ) WRITE( IOUT, 2110 )
     *       TUPDAT + TSOLVE, TFACTR + T1STSL
      END IF
      INFORM = 0
      RETURN
C
C NON-EXECUTABLE STATEMENTS.
C
 2000 FORMAT( ' ** Message from -PRECN-', /,
     *        '    Increase the parameter -LWK-' )
 2010 FORMAT( ' ** Message from -PRECN-', /,
     *        '    Increase the parameter -LIWK-' )
 2020 FORMAT( ' ** Message from -PRECN-', /,
     *        '    Value of IFLAG after MA27A = ', I3 )
 2030 FORMAT( ' ** Message from -PRECN-', /,
     *        '    Value of IFLAG after MA27B = ', I3 )
 2040 FORMAT( ' Perturbation ', 1P, D12.4, ' for diagonal ', I6 )
 2050 FORMAT( ' ** Message from -PRECN-', /,
     *        '    Value of INFOAS after ASSLC = ', I3 )
 2060 FORMAT( /,
     *        ' ** Preconditioner: diagonals are perturbed by at most ',
     *        1P, D12.4, /,
     *        '    Order of preconditioner         = ', I8, /,
     *        '    # nonzeros in preconditioner    = ', I8, /,
     *        '    Predicted # nonzeros in factors = ', I8, /,
     *        '    Actual    # nonzeros in factors = ', I8, /,
     *        '    Smallest diagonal               = ', 1P, D12.4, /,
     *        '    Largest diagonal                = ', 1P, D12.4 )
 2070 FORMAT( /, '    -- Preconditioner formed. N = ', I8, /,
     *        '    -- NNZ  = ', I8, ' Max pert  = ', 1P, D8.1, /,
     *        '    -- Pred fill = ', I8, ' fill = ', I8, /,
     *        '    -- Diag min = ', 1P, D8.1, ' Dia max = ', 1P, D8.1 )
 2090 FORMAT( /,
     *        ' ** Preconditioner: ', /,
     *        '    order of preconditioner         = ', I8, /,
     *        '    # nonzeros in preconditioner    = ', I8, /,
     *        '    Predicted # nonzeros in factors = ', I8, /,
     *        '    actual    # nonzeros in factors = ', I8, /,
     *        '    # negative 1 x 1 block pivots   = ', I8, /,
     *        '    #          2 x 2 block pivots   = ', I8 )
 2100 FORMAT( /, '    -- Preconditioner formed.', /,
     *        '    -- N = ', I8, ' NNZ  = ', I8, 
     *        '    Pred fill = ', I8, ' Fill = ', I8, /,
     *        '    -- # negative 1x1 pivots = ', I8, 
     *              ' # 2x2 pivots = ', I8 )
 2110 FORMAT( /, ' T(updated) = ', F7.3, ' vs T(factored) = ', F7.3 )
 2120 FORMAT( /,
     *        ' ** Preconditioner: ', /,
     *        '    Order of preconditioner         = ', I8, /,
     *        '    Semi-bandwidth used             = ', I8, /,
     *        '    True semi-bandwidth             = ', I8 )
 2130 FORMAT( /, '    -- Preconditioner formed.', /,
     *        '    -- N = ', I8, ' Semi-bandwidth = ', I8, 
     *        '    true semi-bandwith = ', I8 )
 2140 FORMAT( /,
     *        ' ** Preconditioner: ', /,
     *        '    order of preconditioner         = ', I8, /,
     *        '    # nonzeros in preconditioner    = ', I8, /,
     *        '    # nonzeros in factors           = ', I8 )
 2150 FORMAT( /, '    -- Preconditioner formed.', /,
     *        '    -- N = ', I8, ' NNZ  = ', I8, ' fill = ', I8 )
 2160 FORMAT( /, '  Warning message from ICCG. IFLAG = ', I2 )
 2170 FORMAT( /, ' Refactorizing: update ', I6,
     *           ' out of an allowed total of ', I6, 
     *        /, '                Time to update = ', F7.1, 
     *           ' v.s. time to refactorize = ',F7.1 )
 2200 FORMAT( /, ' Solve ' )
 2210 FORMAT( /, ' Factorization ' )
 2220 FORMAT( /, ' time(ASMBL) = ', F9.1 )
 2230 FORMAT( /, ' time(MA27A) = ', F9.1 )
 2240 FORMAT( /, ' time(MA27B) = ', F9.1 )
 2250 FORMAT( /, ' time(SYPRC) = ', F9.1 )
 2260 FORMAT( /, ' Update ' )
 2270 FORMAT( /, ' time(ASSLC) = ', F9.1 )
 2280 FORMAT( /, ' time(MA27C) = ', F9.1 )
 2290 FORMAT( /, ' time(ASSLB) = ', F9.1 )
C1100 FORMAT( 4I6, 1P, D12.4 )
C1120 FORMAT( 3I6, 1P, D12.4 )
C1130 FORMAT( I6, 1P, D12.4 )
C
C  END OF PRECN.
C
      END
