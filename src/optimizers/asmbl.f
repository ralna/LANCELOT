C ** Correction report.
C ** Correction 1. 12/10/92: 3 lines corrected **
C ** Correction 2. 25/05/94: 1 line added **
C ** Correction 3. 25/05/94: 2 lines added **
C ** Correction 4. 25/05/94: 2 lines added **
C ** Correction 5. 20/03/96: 2 lines interchanged **
C ** Correction 6. 05/11/97: 1 line corrected **
C ** Correction 7. 05/11/97: 1 line corrected **
C ** Correction 8. 05/11/97: 3 lines changed to 4 **
C ** Correction 9. 05/11/97: 1 line corrected **
C ** Correction 10. 05/11/97: 8 lines changed to 9 **
C ** Correction 11. 05/11/97: 1 line corrected **
C ** Correction 12. 19/08/98: 1 line corrected **
C ** Correction 13. 19/08/98: 1 line corrected **
C ** Correction 14. 19/08/98: 1 line corrected **
C ** Correction 15. 19/08/98: 1 line corrected **
C  THIS VERSION: 19/08/1998 AT 18:00:51 PM.
CS    SUBROUTINE SASMBL( N , NG, MAXSEL, NSEMIB, LH    , LIH   , NNZH  ,
CD    SUBROUTINE DASMBL( N , NG, MAXSEL, NSEMIB, LH    , LIH   , NNZH  ,
     *                   NFREE , IFREE , ISTADH, LSTADH, ICNA  , LICNA ,  
     *                   ISTADA, LSTADA, INTVAR, LNTVAR, IELVAR, LELVAR, 
     *                   IELING, LELING, ISTADG, LSTADG, ISTAEV, LSTAEV, 
     *                   ISTAGV, LNSTGV, ISVGRP, LNVGRP, IRNH  , JCNH  ,   
     *                   NXTROW, LNXTRW, IWK   , LIWK  , A , LA, GUVALS, 
     *                   LNGUVL, HUVALS, LNHUVL, GVALS2, GVALS3, GSCALE, 
     *                   ESCALE, LESCAL, H , WK, LWK   , GXEQX , LGXEQX, 
     *                   INTREP, LINTRE, RANGES, IPRINT, IOUT  , BAND  ,   
C ** Correction 6. 05/11/97: 1 line corrected **
C ** Correction 12. 19/08/98: 1 line corrected **
     *                   MAXSBW, INFORM, NOZERO, FIXSTR )
C ** Correction 12. 19/08/98: end of correction **
C ** Correction 6. 05/11/97: end of correction **
C
C  ********************************************************************
C
C  ASSEMBLE THE SECOND DERIVATIVE MATRIX OF A GROUPS PARTIALLY
C  SEPARABLE FUNCTION IN EITHER CO-ORDINATE OR BAND FORMAT.
C
C  NICK GOULD, FEBRUARY 20TH 1989.
C  FOR CGT PRODUCTIONS.
C
C  ********************************************************************
C
      INTEGER          N , NG, MAXSEL, NFREE , LNXTRW, NNZH
      INTEGER          LA, LH, IOUT  , INFORM, NSEMIB, IPRINT, LIH
      INTEGER          LSTADH, LICNA , LSTADA, LNTVAR, LELVAR, LELING
      INTEGER          LSTADG, LSTAEV, LNSTGV, LNVGRP, LIWK,   MAXSBW
      INTEGER          LNGUVL, LNHUVL, LESCAL, LWK,    LGXEQX, LINTRE
C ** Correction 7. 05/11/97: 1 line corrected **
C ** Correction 13. 19/08/98: 1 line corrected **
      LOGICAL          BAND  , NOZERO, FIXSTR
C ** Correction 13. 19/08/98: end of correction **
C ** Correction 7. 05/11/97: end of correction **
      INTEGER          IFREE(  N ),      ISTADH( LSTADH ), ICNA( LICNA )  
      INTEGER          ISTADA( LSTADA ), INTVAR( LNTVAR ), JCNH( LIH )
      INTEGER          IELVAR( LELVAR ), IELING( LELING ), IRNH( LIH )
      INTEGER          ISTADG( LSTADG ), ISTAEV( LSTAEV )
      INTEGER          ISTAGV( LNSTGV ), ISVGRP( LNVGRP ), IWK( LIWK )
      INTEGER          NXTROW( 2,        LNXTRW )
      LOGICAL          GXEQX( LGXEQX ),  INTREP( LINTRE )
CS    REAL             A(     LA ),      GVALS2( NG ),     GVALS3( NG ),
CD    DOUBLE PRECISION A(     LA ),      GVALS2( NG ),     GVALS3( NG ),
     *                 GUVALS( LNGUVL ), HUVALS( LNHUVL ), GSCALE( NG ), 
     *                 ESCALE( LESCAL ), H(      LH ),     WK(     LWK )
      EXTERNAL         RANGES
C
C  LOCAL VARIABLES.
C
      INTEGER          I, II,  IG, J,  JJ, K,  L, IP,  NN,     NNN
      INTEGER          NEWPT,  NIN,    IELL,   IEL,    IHNEXT
      INTEGER          NVAREL, IG1,    LISTVS, LISTVE, IELH
      INTEGER          INEXT,  IJHESS, IROW,   JCOL,   JCOLST, ISTART
CS    REAL             ONE,    ZERO,   WKI,    HESNEW, GDASH,  G2DASH,
CD    DOUBLE PRECISION ONE,    ZERO,   WKI,    HESNEW, GDASH,  G2DASH,
     *                 SCALEE
      CHARACTER * 2    MATRIX( 36, 36 )
CS    EXTERNAL         SSETVL, SSETVI
CD    EXTERNAL         DSETVL, DSETVI
      INTRINSIC        ABS,    MAX,    MIN
C
C  SET CONSTANT REAL PARAMETERS.
C
CS    PARAMETER ( ZERO   = 0.0E+0, ONE    = 1.0E+0 )
CD    PARAMETER ( ZERO   = 0.0D+0, ONE    = 1.0D+0 )
      IF ( .NOT. BAND .AND. NFREE .GT. LNXTRW ) GO TO 620
C
C  RENUMBER THE FREE VARIABLES SO THAT THEY ARE VARIABLES 1 TO NFREE.
C
      DO 10 I     = 1, N
         IWK( I ) = 0
   10 CONTINUE
      DO 20 I = 1, NFREE
         IWK( IFREE( I ) ) = I
C
C  INITIALIZE THE LINK LIST WHICH POINTS TO THE ROW NUMBERS WHICH
C  ARE USED IN THE COLUMNS OF THE ASSEMBLED HESSIAN.
C
C  NXTROW( 1, . ) GIVES THE LINK LIST. THE LIST FOR COLUMN J STARTS
C                 IN NXTROW( 1, J ) AND ENDS WHEN NXTROW( 1, K ) = - 1.
C  NXTROW( 2, . ) GIVES THE POSITION IN H OF THE CURRENT LINK.
C
         IF ( .NOT. BAND ) NXTROW( 1, I ) = - 1
   20 CONTINUE
      IF ( IPRINT .GE. 10 ) THEN
         WRITE( IOUT, 2060 ) NFREE, ( IFREE( I ), I = 1, NFREE )
      END IF
C
C  IF A BAND STORAGE SCHEME IS TO BE USED, INITIALIZE THE ENTRIES
C  WITHIN THE BAND AS ZERO.
C
      IF ( BAND ) THEN
         MAXSBW  = 0
         IF ( NFREE * ( NSEMIB + 1 ) .GT. LH ) GO TO 610
         DO 30 I = 1, NFREE * ( NSEMIB + 1 )
            H( I ) = ZERO
   30    CONTINUE   
      ELSE
         NNZH  = 0
         NEWPT = NFREE + 1
         IF ( NEWPT .GT. LNXTRW ) GO TO 620
      END IF
C
C  ------------------------------------------------------
C  FORM THE RANK-ONE SECOND ORDER TERM FOR THE ITH GROUP.
C  ------------------------------------------------------
C
      DO 200 IG = 1, NG
         IF ( IPRINT .GE. 100 ) WRITE( IOUT, 2070 ) IG
         IF ( GXEQX( IG ) ) GO TO 200
C ** Correction 8. 05/11/97: 3 lines changed to 4 **
C ** Correction 14. 19/08/98: 1 line corrected **
         IF ( .NOT. FIXSTR .AND. GSCALE( IG ) .EQ. ZERO ) GO TO 200
C ** Correction 14. 19/08/98: end of correction **
         G2DASH = GSCALE( IG ) * GVALS3( IG )
C        IF ( IPRINT .GE. 100 )WRITE(6,*) ' GVALS3(IG) ', GVALS3(IG)
         IF ( NOZERO .AND. G2DASH .EQ. ZERO ) GO TO 200
C ** Correction 8. 05/11/97: end of correction **
         IG1    = IG + 1
         LISTVS = ISTAGV( IG )
         LISTVE = ISTAGV( IG1 ) - 1
C
C  FORM THE GRADIENT OF THE IG-TH GROUP.
C
CS       CALL SSETVI( LISTVE - LISTVS + 1, WK, ISVGRP( LISTVS ), ZERO )
CD       CALL DSETVI( LISTVE - LISTVS + 1, WK, ISVGRP( LISTVS ), ZERO )
C
C  CONSIDER ANY NONLINEAR ELEMENTS FOR THE GROUP.
C
         DO 130 IELL = ISTADG( IG ), ISTADG( IG1 ) - 1
            IEL      = IELING( IELL )
            K        = INTVAR( IEL )
            L        = ISTAEV( IEL )
            NVAREL   = ISTAEV( IEL + 1 ) - L
            SCALEE   = ESCALE( IELL )
            IF ( INTREP( IEL ) ) THEN
C
C  THE IEL-TH ELEMENT HAS AN INTERNAL REPRESENTATION.
C
               NIN = INTVAR( IEL + 1 ) - K
               CALL RANGES( IEL, .TRUE., GUVALS( K ),
     *                      WK( N + 1 ), NVAREL, NIN )
               DO 110 I   = 1, NVAREL
                  J       = IELVAR( L )
                  WK( J ) = WK( J ) + SCALEE * WK( N + I )
                  L       = L + 1
  110          CONTINUE
            ELSE
C
C  THE IEL-TH ELEMENT HAS NO INTERNAL REPRESENTATION.
C
               DO 120 I   = 1, NVAREL
                  J       = IELVAR( L )
                  WK( J ) = WK( J ) + SCALEE * GUVALS( K )
                  K       = K + 1
                  L       = L + 1
  120          CONTINUE
            END IF
  130    CONTINUE
C
C  INCLUDE THE CONTRIBUTION FROM THE LINEAR ELEMENT.
C
         DO 140 K   = ISTADA( IG ), ISTADA( IG1 ) - 1
            J       = ICNA( K )
            WK( J ) = WK( J ) + A( K )
  140    CONTINUE
C
C  THE GRADIENT IS COMPLETE. FORM THE J-TH COLUMN OF THE RANK-ONE MATRIX
C
         DO 190 L = LISTVS, LISTVE
            JJ    = ISVGRP( L )
            J     = IWK( JJ )
            IF ( J .EQ. 0 ) GO TO 190
C
C  FIND THE ENTRY IN ROW I OF THIS COLUMN.
C
            DO 180 K = LISTVS, LISTVE
               II    = ISVGRP( K )
               I     = IWK( II )
               IF ( I .EQ. 0 .OR. I .GT. J ) GO TO 180
C
C  SKIP ALL ELEMENTS WHICH LIE OUTSIDE A BAND OF WIDTH NSEMIB.
C
               IF ( BAND ) MAXSBW = MAX( MAXSBW, J - I )
               IF ( J - I .GT. NSEMIB ) GO TO 180
               HESNEW = WK( II ) * WK( JJ ) * G2DASH
               IF ( IPRINT .GE. 100 ) WRITE( IOUT, 2090 ) I, J, HESNEW
C ** Correction 9. 05/11/97: 1 line corrected **
               IF ( NOZERO .AND. HESNEW .EQ. ZERO ) GO TO 180
C ** Correction 9. 05/11/97: end of correction **
C
C  OBTAIN THE APPROPRIATE STORAGE LOCATION IN H FOR THE NEW ENTRY.
C
C
C  CASE 1: BAND MATRIX STORAGE SCHEME.
C
               IF ( BAND ) THEN
C
C  THE ENTRY BELONGS ON THE DIAGONAL.
C
                  IF ( I .EQ. J ) THEN
                     H( I ) = H( I ) + HESNEW
C
C  THE ENTRY BELONGS OFF THE DIAGONAL.
C
                  ELSE
                     H( NFREE + J - I + NSEMIB * ( I - 1 ) ) = 
     *         H( NFREE + J - I + NSEMIB * ( I - 1 ) ) + HESNEW
                  END IF   
C
C  CASE 2: COORDINATE STORAGE SCHEME.
C
               ELSE
                  ISTART = J
  150             CONTINUE
                  INEXT = NXTROW( 1, ISTART )
                  IF ( INEXT .EQ. - 1 ) THEN
                     IF ( NEWPT .GT. LNXTRW ) GO TO 620
C
C  THE (I,J)-TH LOCATION IS EMPTY. PLACE THE NEW ENTRY IN THIS LOCATION
C  AND ADD ANOTHER LINK TO THE LIST.
C
                     NNZH = NNZH + 1
                     IF ( NNZH .GT. LH .OR. NNZH .GT. LIH ) GO TO 610
                     IRNH( NNZH )        = I
                     JCNH( NNZH )        = J
                     H( NNZH )           = HESNEW
                     NXTROW( 1, ISTART ) = NEWPT
                     NXTROW( 2, ISTART ) = NNZH
                     NXTROW( 1, NEWPT )  = - 1
                     NEWPT               = NEWPT + 1
                  ELSE
C
C  CONTINUE SEARCHING THE LINKED LIST FOR AN ENTRY IN ROW I, COLUMN J.
C
                     IF ( IRNH( NXTROW( 2, ISTART ) ) .EQ. I ) THEN
                        IP = NXTROW( 2, ISTART )
C                       WRITE( 6, 2110 ) IP, H( IP )
                        H( IP ) = H( IP ) + HESNEW
C                       WRITE( 6, 2120 ) IP, H( IP )
                     ELSE
                        ISTART = INEXT
                        GO TO 150
                     END IF
                  END IF
               END IF
  180       CONTINUE
  190    CONTINUE
  200 CONTINUE
C
C  RESET THE WORKSPACE ARRAY TO ZERO.
C
CS    CALL SSETVL( MAXSEL, WK, 1, ZERO )
CD    CALL DSETVL( MAXSEL, WK, 1, ZERO )
C
C  --------------------------------------------------------
C  ADD ON THE LOW RANK FIRST ORDER TERMS FOR THE I-TH GROUP.
C  --------------------------------------------------------
C
      DO 300 IG = 1, NG
         IF ( IPRINT .GE. 100 ) WRITE( IOUT, 2100 ) IG
C ** Correction 10. 05/11/97: 8 lines changed to 9 **
C ** Correction 15. 19/08/98: 1 line corrected **
         IF ( .NOT. FIXSTR .AND. GSCALE( IG ) .EQ. ZERO ) GO TO 300
C ** Correction 15. 19/08/98: end of correction **
         IF ( GXEQX( IG ) ) THEN
            GDASH = GSCALE( IG )
         ELSE
            GDASH = GSCALE( IG ) * GVALS2( IG )
C           IF ( IPRINT .GE. 100 )WRITE(6,*) ' GVALS2(IG) ', GVALS2(IG)
C ** Correction 5. 20/03/96: 2 lines interchanged **
            IF ( NOZERO .AND. GDASH .EQ. ZERO ) GO TO 300
         END IF
C ** Correction 5. 20/03/96: end of correction **
C ** Correction 10. 05/11/97: end of correction **
         IG1 = IG + 1
C
C  SEE IF THE GROUP HAS ANY NONLINEAR ELEMENTS.
C
         DO 290 IELL = ISTADG( IG ), ISTADG( IG1 ) - 1
            IEL      = IELING( IELL )
            LISTVS   = ISTAEV( IEL )
            LISTVE   = ISTAEV( IEL + 1 ) - 1
            NVAREL   = LISTVE - LISTVS + 1
            IELH     = ISTADH( IEL )
            IHNEXT   = IELH
            SCALEE   = ESCALE( IELL )
            DO 250 L = LISTVS, LISTVE
               J     = IWK( IELVAR( L ) )
               IF ( J .NE. 0 ) THEN
C
C  THE IEL-TH ELEMENT HAS AN INTERNAL REPRESENTATION.
C  COMPUTE THE J-TH COLUMN OF THE ELEMENT HESSIAN MATRIX.
C
                  IF ( INTREP( IEL ) ) THEN
C
C  COMPUTE THE J-TH COLUMN OF THE HESSIAN.
C
                     WK( L - LISTVS + 1 ) = ONE
C
C  FIND THE INTERNAL VARIABLES.
C
                     NIN = INTVAR( IEL + 1 ) - INTVAR( IEL )
                     CALL RANGES( IEL, .FALSE., WK( 1 ),
     *                            WK( MAXSEL + 1 ), NVAREL, NIN )
C
C  MULTIPLY THE INTERNAL VARIABLES BY THE ELEMENT HESSIAN.
C
                     NN = MAXSEL + NIN
CS                   CALL SSETVL( NIN, WK( NN + 1 ), 1, ZERO )
CD                   CALL DSETVL( NIN, WK( NN + 1 ), 1, ZERO )
C
C  ONLY THE UPPER TRIANGLE OF THE ELEMENT HESSIAN IS STORED.
C
                     JCOLST = IELH - 1
                     DO 220 JCOL = 1, NIN
                        IJHESS   = JCOLST
                        JCOLST   = JCOLST + JCOL
                        WKI      = WK( MAXSEL + JCOL ) * GDASH
                        DO 210 IROW = 1, NIN
                           IF ( IROW .LE. JCOL ) THEN
                              IJHESS = IJHESS + 1
                           ELSE
                              IJHESS = IJHESS + IROW - 1
                           END IF
                           WK( NN + IROW ) = WK( NN + IROW ) +
     *                                       WKI * HUVALS( IJHESS )
  210                   CONTINUE
  220                CONTINUE
C
C  SCATTER THE PRODUCT BACK ONTO THE ELEMENTAL VARIABLES.
C
                     NNN = NN + NIN
                     CALL RANGES( IEL, .TRUE., WK( NN + 1 ),
     *                            WK( NNN + 1 ), NVAREL, NIN )
                     WK( L - LISTVS + 1 ) = ZERO
                  END IF
C
C  FIND THE ENTRY IN ROW I OF THIS COLUMN.
C
                  DO 240 K = LISTVS, L
                     I     = IWK( IELVAR( K ) )
C
C  SKIP ALL ELEMENTS WHICH LIE OUTSIDE A BAND OF WIDTH NSEMIB.
C
C ** Correction 1. 12/10/92: 3 lines corrected **
                     IF ( BAND .AND. I .NE. 0 ) MAXSBW = 
     *                    MAX( MAXSBW, ABS( J - I ) )
                     IF ( ABS( I - J ) .LE. NSEMIB .AND. I .NE. 0 ) THEN
C ** Correction 1. 12/10/92: end of correction **
C
C  ONLY THE UPPER TRIANGLE OF THE MATRIX IS STORED.
C
                        IF ( I .LE. J ) THEN
                           II = I
                           JJ = J
                        ELSE
                           II = J
                           JJ = I
                        END IF
C
C  OBTAIN THE APPROPRIATE STORAGE LOCATION IN H FOR THE NEW ENTRY.
C
                        IF ( INTREP( IEL ) ) THEN
                           HESNEW = SCALEE * WK( NNN + K - LISTVS + 1 )
                        ELSE
                           HESNEW = SCALEE * HUVALS( IHNEXT ) * GDASH
                        END IF
                        IF ( IPRINT .GE. 100 )
     *                     WRITE( IOUT, 2080 ) II, JJ, IEL, HESNEW
C
C  CASE 1: BAND MATRIX STORAGE SCHEME.
C
                        IF ( BAND ) THEN
C
C  THE ENTRY BELONGS ON THE DIAGONAL.
C
                           IF ( II .EQ. JJ ) THEN
                              H( II ) = H( II ) + HESNEW
C ** Correction 2. 25/05/94: 1 line added **
                              IF ( K .NE. L ) H( II ) = H( II ) + HESNEW
C ** Correction 2. 25/05/94: end of correction **
C
C  THE ENTRY BELONGS OFF THE DIAGONAL.
C
                           ELSE
                              H( NFREE + JJ - II + NSEMIB * ( II - 1 ) )
     *                     =  H( NFREE + JJ - II + NSEMIB * ( II - 1 ) )
     *                       + HESNEW
                           END IF   
C
C  CASE 2: COORDINATE STORAGE SCHEME.
C
                        ELSE
C ** Correction 11. 05/11/97: 1 line corrected **
                           IF ( .NOT. NOZERO .OR. HESNEW .NE. ZERO) THEN
C ** Correction 11. 05/11/97: end of correction **
                              ISTART = JJ
  230                         CONTINUE
                              INEXT = NXTROW( 1, ISTART )
                              IF ( INEXT .EQ. - 1 ) THEN
                                 IF ( NEWPT .GT. LNXTRW ) GO TO 620
C
C  THE (I,J)-TH LOCATION IS EMPTY. PLACE THE NEW ENTRY IN THIS LOCATION
C  AND ADD ANOTHER LINK TO THE LIST.
C
                                 NNZH = NNZH + 1
                                 IF ( NNZH .GT. LH .OR. 
     *                                NNZH .GT. LIH ) GO TO 610
                                 IRNH( NNZH )        = II
                                 JCNH( NNZH )        = JJ
                                 H( NNZH )           = HESNEW
C ** Correction 3. 25/05/94: 2 lines added **
                                 IF ( K .NE. L .AND. II .EQ. JJ )  
     *                              H( NNZH ) = HESNEW + HESNEW
C ** Correction 3. 25/05/94: end of correction **
                                 NXTROW( 1, ISTART ) = NEWPT
                                 NXTROW( 2, ISTART ) = NNZH
                                 NXTROW( 1, NEWPT )  = - 1
                                 NEWPT               = NEWPT + 1
                              ELSE
C
C  CONTINUE SEARCHING THE LINKED LIST FOR AN ENTRY IN ROW I, COLUMN J.
C
                              IF ( IRNH( NXTROW( 2, ISTART ) )
     *                             .EQ. II ) THEN
                                    IP      = NXTROW( 2, ISTART )
C                                   WRITE( 6, 2110 ) IP, H( IP )
                                    H( IP ) = H( IP ) + HESNEW
C ** Correction 4. 25/05/94: 2 lines added **
                                    IF ( K .NE. L .AND. II .EQ. JJ )  
     *                                 H( IP ) = H( IP ) + HESNEW
C ** Correction 4. 25/05/94: end of correction **
C                                   WRITE( 6, 2120 ) IP, H( IP )
                                 ELSE
                                    ISTART = INEXT
                                    GO TO 230
                                 END IF
                              END IF
                           END IF
                        END IF
                     END IF
                     IHNEXT = IHNEXT + 1
  240             CONTINUE
               ELSE
                  IHNEXT = IHNEXT + L - LISTVS + 1
               END IF
  250       CONTINUE
  290    CONTINUE
  300 CONTINUE
C
C  ----------------------------------------
C  FOR DEBUGGING, PRINT THE NONZERO VALUES.
C  ----------------------------------------
C
      IF ( IPRINT .GE. 10 ) THEN
         IF ( .NOT. BAND ) WRITE( IOUT, 2000 )
     *      ( IRNH( I ), JCNH( I ), H( I ), I = 1, NNZH )
C        if(iprint.ge.1000) then
C           write(62,*) 'ELEMENT USES'
C           DO 320 I = 1, NNZH
C              if ( irnh(i).eq.jcnh(i))then
C                 write(62,7000)irnh(i),jcnh(i),
C    *                          irnh(i),jcnh(i),irnh(i)
C              else
C                 write(62,7010)irnh(i),jcnh(i),
C    *                          irnh(i),jcnh(i),irnh(i),
C    *                          irnh(i),jcnh(i),jcnh(i)
C              end if   
C7000       format( 1x, 'T', 2x, 'D', 2i4, 1x, 'DIAG',/,
C    *              1x, 'V', 2x, 'D', 2i4, 1x, 'X', 24x, i8 )
C7010       format( 1x, 'T', 2x, 'O', 2i4, 1x, 'OFFDIAG',/,
C    *              1x, 'V', 2x, 'O', 2i4, 1x, 'X', 24x, i8, /,
C    *              1x, 'V', 2x, 'O', 2i4, 1x, 'Y', 24x, i8 )
C 320       CONTINUE
C           write(62,*) 'GROUP USES'
C           DO 330 I = 1, NNZH
C              if ( irnh(i).eq.jcnh(i))then
C                 write(62,8000)irnh(i),jcnh(i),h(i)
C              else
C                 write(62,8010)irnh(i),jcnh(i),h(i)
C              end if   
C8000       format( 1x, 'E', 2x, 'LINGROUP', 2x, 'D', 2i4, 1x,
C    *                1p, d12.4 )
C8010       format( 1x, 'E', 2x, 'LINGROUP', 2x, 'O', 2i4, 1x,
C    *                1p, d12.4 )
C 330       CONTINUE
C        end if   
C
C  FOR DEBUGGING, PRINT THE NONZERO PATTERN OF THE MATRIX.
C
         IF ( NFREE .LE. 36 ) THEN
            DO 420 I = 1, NFREE
               DO 410 J = 1, NFREE
                  MATRIX( I, J ) = '  '
  410          CONTINUE
  420       CONTINUE
            IF ( BAND ) THEN
               DO 440 I = 1, NFREE
                  IF ( H( I ) .NE. ZERO ) MATRIX( I, I ) = ' *'
                  DO 430 J = 1, MIN( NSEMIB, NFREE - I )
                     K = NFREE + J + NSEMIB * ( I - 1 )
                     IF ( H( K ) .NE. ZERO ) THEN
                        MATRIX( I + J, I ) = ' *'
                        MATRIX( I, I + J ) = ' *'
                     END IF
  430             CONTINUE
  440          CONTINUE
            ELSE
               DO 450 I = 1, NNZH
                  IF ( IRNH( I ) .GT. NFREE ) THEN
                     WRITE( IOUT, * ) ' ENTRY OUT OF BOUNDS IN ASEMBL ',
     *                                ' ROW NUMBER = ', IRNH( I )
                     STOP
                  END IF
                  IF ( JCNH( I ) .GT. NFREE ) THEN
                     WRITE( IOUT, * ) ' ENTRY OUT OF BOUNDS IN ASEMBL ',
     *                                ' COL NUMBER = ', JCNH( I )
                     STOP
                  END IF
                  MATRIX( IRNH( I ), JCNH( I ) ) = ' *'
                  MATRIX( JCNH( I ), IRNH( I ) ) = ' *'
  450          CONTINUE
            END IF
            WRITE( IOUT, 2040 ) ( I, I = 1, NFREE )
            DO 460 I = 1, NFREE
               WRITE( IOUT, 2050 ) I, ( MATRIX( I, J ), J = 1, NFREE )
  460       CONTINUE
         END IF
      END IF
      INFORM = 0
      RETURN
C
C  UNSUCCESSFUL RETURNS.
C
  610 CONTINUE
      INFORM = 1
      IF ( LH .LE. LIH ) THEN
         WRITE( IOUT, 2010 ) LH
      ELSE
         WRITE( IOUT, 2030 ) LIH
      END IF
      RETURN
  620 CONTINUE
      INFORM = 2
      WRITE( IOUT, 2020 ) LNXTRW
      RETURN
C
C  NON-EXECUTABLE STATEMENTS
C
 2000 FORMAT( '    Row  Column    Value        Row  Column    Value ', /
     *        '    ---  ------    -----        ---  ------    ----- ', /
     *        ( 2I6, 1P, D24.16, 2I6, 1P, D24.16 ) )
 2010 FORMAT( ' ** Array dimension LH =', I6, ' too small in ASMBL. ' )
 2020 FORMAT( ' ** Array dimension LNXTRW =',I8, ' too small in ASMBL.')
 2030 FORMAT( ' ** Array dimension LIH =', I6, ' too small in ASMBL.' )
 2040 FORMAT( /, 5X, 36I2 )
 2050 FORMAT( I3, 2X, 36A2 )
 2060 FORMAT( /, I6, ' free variables. They are ', 8I5, /, ( 14I5 ) )
 2070 FORMAT( ' Group ', I5, ' rank-one terms ' )
 2080 FORMAT( ' Row ', I6, ' Column ', I6, ' used from element ', I6,
     *        ' value = ', 1P, D24.16 )
 2090 FORMAT( ' Row ', I6, ' column ', I6, ' used. Value = ',1P,D24.16)
 2100 FORMAT( ' Group ', I5, ' second-order terms ' )
C2110 FORMAT( ' Modifying existing H. IP,    H = ', I3, 1P, D12.4 )
C2120 FORMAT( ' Modifying existing H. IP, NEWH = ', I3, 1P, D12.4 )
C
C  END OF ASMBL.
C
      END
