* used to be tobig.f

      PROGRAM TOBIG

*  SPECIALIZES SOURCE FILES TO LARGE SIZE LANCELOT

*  PH.L. TOINT AND A. CONN, JUNE 1991.

      INTEGER      INDEV, OUTDEV, EOL, K
      CHARACTER*1  LINE(80)
      CHARACTER*5  HEAD
      PARAMETER    ( INDEV = 5, OUTDEV = 6 )

*  LOOP ON THE LINES

 10   CONTINUE
      READ( INDEV, '(80A1)', END = 99 ) ( LINE(K), K = 1, 80 )

*  DETECT END OF LINE

      DO 20 EOL = 80, 1, -1
         IF( LINE(EOL) .NE. ' ' ) GO TO 30
 20   CONTINUE

*  ACTIVATE THE STATEMENTS SUITABLE FOR  THE DESIRED ENVIRONMENT

 30   CONTINUE
      IF( LINE(1) .EQ. 'C' ) THEN
        DO 25 K = 1, 5
          HEAD(K:K) = LINE(K)
 25     CONTINUE
        IF( HEAD .EQ. 'CBIG ' ) GO TO 50
      END IF
      GO TO 40

*  REPLACE HEAD BY 5 BLANKS

 50   CONTINUE
      DO 60 K = 1, 5
        LINE(K) = ' '
 60   CONTINUE

*  PRINT LINE

 40   CONTINUE
      WRITE( OUTDEV, '(80A1)' ) ( LINE(K), K = 1, EOL )
      GO TO 10

*  END OF FILE

 99   CONTINUE
      STOP
      END