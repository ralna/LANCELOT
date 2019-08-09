C  THIS VERSION: 11/08/1993 AT 08:58:48 AM.
CS    SUBROUTINE MA31D ( A, IRN, IA, N, IK, IP, ROW )
CD    SUBROUTINE MA31DD( A, IRN, IA, N, IK, IP, ROW )
      INTEGER           IA, N
      LOGICAL           ROW
      INTEGER           IK( N ), IP( N ), IRN( IA )
CS    REAL              A( IA )
CD    DOUBLE PRECISION  A( IA )
C
C  DUMMY SUBROUTINE AVAILABLE WITH LANCELOT.
C  NB. MA31 IS NO LONGER AVAILABLE FROM HSL OR ITS ARCHIVE
C
C  NICK GOULD, 7TH NOVEMBER 1990.
C  FOR CGT PRODUCTIONS.
C
C     WRITE( 6, 2000 )
C     STOP
C
C  NON-EXECUTABLE STATEMENTS.
C
C2000 FORMAT( /, ' We regret that the solution options that you have ',
C    *        /, ' chosen are not all freely available with LANCELOT. ',
C    *        //, ' If you have the Harwell Subroutine Library, this ',
C    *        /, ' option may be enabled by replacing the dummy ', /,
C    *        ' subroutine MA31D (single precision) or MA31DD', /,
C    *        ' (double precision) with its H.S.L. namesake ', /,
C    *        ' and dependencies.', //,
C    *        ' *** EXECUTION TERMINATING *** ', / )
C
      COMMON/MA31JD/LROW,LCOL,NCP,ND,IPD
      COMMON/MA31KD/NURL,NUCL,NUAL
      NCP=NCP+1
C
C STORE FIRST ELEMENT OF ENTRY IN IK(J). THEN OVERWRITE IT BY -J
C
      DO 5 J=1,N
       NZ=IK(J)
       NN=NUAL
       IF (.NOT.ROW) NN=NUCL
       IF (NZ.GT.0 .AND. IP(J).GE.NN) THEN
         K=IP(J)
         IK(J)=IRN(K)
         IRN(K)=-J
       END IF
    5 CONTINUE
C
C KN IS POSITION OF NEXT ENTRY IN COMPRESSED FILE.
C
      KN=IA+1
      IPI=IA+1
      KL=IA-NUCL+1
      IF (ROW) KL=IA-NUAL+1
C
C LOOP THROUGH OLD FILE SKIPPING ZERO ELEMENTS AND MOVING GENUINE
C ELEMENTS FORWARD. THE ENTRY NUMBER BECOMES KNOWN ONLY WHEN
C ITS END IS DETECTED BY PRESENCE OF A NEGATIVE INTEGER.
C
      DO 15 KK=1,KL
       K=IA+1-KK
       IF (IRN(K).NE.0) THEN
        KN=KN-1
        IF (ROW) A(KN)=A(K)
C
C END OF ENTRY. RESTORE IRN(K), SET POINTERS TO START OF ENTRY AND
C STORE CURRENT KN IN IPI READY FOR USE WHEN NEXT LAST ENTRY IS
C DETECTED.
C
         IF (IRN(K).LT.0) THEN
          J=-IRN(K)
          IRN(K)=IK(J)
          IP(J)=KN
          IK(J)=IPI-KN
          IPI=KN
          END IF
         IRN(KN)=IRN(K)
       END IF
   15 CONTINUE
      IF (ROW) THEN
       NUAL=KN
      ELSE
       NUCL=KN
      END IF
      RETURN
C
C  END OF DUMMY SUBROUTINE.
C
      END
