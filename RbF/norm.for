	PROGRAM NORM
	INTEGER I,J,N,IFAIL
	PARAMETER (N=1115)
	DOUBLE PRECISION KOFF,PI,Y(N),ER,RESULT,OMEGA(N),GRE(N),GIM(N)
	PI=3.1415926D0
	KOFF=1.0D0
        OPEN(1,FILE='DATA1.TXT')
	DO J=1,N
	  READ(1,1)OMEGA(J),GIM(J),GRE(J)
	ENDDO
	CLOSE(1)
10	DO J=1,N
	  Y(J)=OMEGA(J)*GIM(J)
        ENDDO
	CALL D01GAF(OMEGA,Y,N,RESULT,ER,IFAIL)
        IF (DABS(RESULT-PI/2.0D0).GT.0.01D0) THEN
        DO I=1,N
	  GIM(I)=GIM(I)*1.002D0
	ENDDO
        KOFF=KOFF*1.002
        PRINT*,result
	GOTO 10
        ENDIF
        PRINT*,KOFF
1       FORMAT(3E16.8)
	END
      
      SUBROUTINE D01GAF(X,Y,N,ANS,ER,IFAIL)
C
C     THIS SUBROUTINE INTEGRATES A FUNCTION (Y) SPECIFIED
C     NUMERICALLY AT N POINTS (X), WHERE N IS AT LEAST 4,
C     OVER THE RANGE X(1) TO X(N).  THE POINTS NEED NOT BE
C     EQUALLY SPACED, BUT SHOULD BE DISTINCT AND IN ASCENDING
C     OR DESCENDING ORDER.  AN ERROR ESTIMATE IS RETURNED.
C     THE METHOD IS DUE TO GILL AND MILLER.
C
C     NAG COPYRIGHT 1975
C     MARK 5 RELEASE
C     MARK 7 REVISED IER-154 (DEC 1978)
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D01GAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ANS, ER
      INTEGER           IFAIL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  X(N), Y(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  C, D1, D2, D3, H1, H2, H3, H4, R1, R2, R3, R4, S
      INTEGER           I, NN
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Executable Statements ..
      ANS = 0.0D0
      ER = 0.0D0
      IF (N.GE.4) GO TO 20
      IFAIL = P01ABF(IFAIL,1,SRNAME,0,P01REC)
      RETURN
C
C     CHECK POINTS ARE STRICTLY INCREASING OR DECREASING
C
   20 H2 = X(2) - X(1)
      DO 80 I = 3, N
         H3 = X(I) - X(I-1)
         IF (H2*H3) 40, 60, 80
   40    IFAIL = P01ABF(IFAIL,2,SRNAME,0,P01REC)
         RETURN
   60    IFAIL = P01ABF(IFAIL,3,SRNAME,0,P01REC)
         RETURN
   80 CONTINUE
C
C     INTEGRATE OVER INITIAL INTERVAL
C
      D3 = (Y(2)-Y(1))/H2
      H3 = X(3) - X(2)
      D1 = (Y(3)-Y(2))/H3
      H1 = H2 + H3
      D2 = (D1-D3)/H1
      H4 = X(4) - X(3)
      R1 = (Y(4)-Y(3))/H4
      R2 = (R1-D1)/(H4+H3)
      H1 = H1 + H4
      R3 = (R2-D2)/H1
      ANS = H2*(Y(1)+H2*(D3/2.0D0-H2*(D2/6.0D0-(H2+2.0D0*H3)*R3/12.0D0))
     *      )
      S = -(H2**3)*(H2*(3.0D0*H2+5.0D0*H4)+10.0D0*H3*H1)/60.0D0
      R4 = 0.0D0
C
C     INTEGRATE OVER CENTRAL PORTION OF RANGE
C
      NN = N - 1
      DO 120 I = 3, NN
         ANS = ANS + H3*((Y(I)+Y(I-1))/2.0D0-H3*H3*(D2+R2+(H2-H4)*R3)
     *         /12.0D0)
         C = H3**3*(2.0D0*H3*H3+5.0D0*(H3*(H4+H2)+2.0D0*H4*H2))/120.0D0
         ER = ER + (C+S)*R4
         IF (I.NE.3) S = C
         IF (I.EQ.3) S = S + 2.0D0*C
         IF (I-N+1) 100, 140, 100
  100    H1 = H2
         H2 = H3
         H3 = H4
         D1 = R1
         D2 = R2
         D3 = R3
         H4 = X(I+2) - X(I+1)
         R1 = (Y(I+2)-Y(I+1))/H4
         R4 = H4 + H3
         R2 = (R1-D1)/R4
         R4 = R4 + H2
         R3 = (R2-D2)/R4
         R4 = R4 + H1
         R4 = (R3-D3)/R4
  120 CONTINUE
C
C     INTEGRATE OVER FINAL INTERVAL
C
  140 CONTINUE
      ANS = ANS + H4*(Y(N)-H4*(R1/2.0D0+H4*(R2/6.0D0+(2.0D0*H3+H4)
     *      *R3/12.0D0)))
      ER = ER - H4**3*R4*(H4*(3.0D0*H4+5.0D0*H2)+10.0D0*H3*(H2+H3+H4))
     *     /60.0D0 + S*R4
      ANS = ANS + ER
      IFAIL = 0
      RETURN
      END

        
      INTEGER FUNCTION P01ABF(IFAIL,IERROR,SRNAME,NREC,REC)
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C     MARK 13 REVISED. IER-621 (APR 1988).
C     MARK 13B REVISED. IER-668 (AUG 1988).
C
C     P01ABF is the error-handling routine for the NAG Library.
C
C     P01ABF either returns the value of IERROR through the routine
C     name (soft failure), or terminates execution of the program
C     (hard failure). Diagnostic messages may be output.
C
C     If IERROR = 0 (successful exit from the calling routine),
C     the value 0 is returned through the routine name, and no
C     message is output
C
C     If IERROR is non-zero (abnormal exit from the calling routine),
C     the action taken depends on the value of IFAIL.
C
C     IFAIL =  1: soft failure, silent exit (i.e. no messages are 
C                 output)
C     IFAIL = -1: soft failure, noisy exit (i.e. messages are output)
C     IFAIL =-13: soft failure, noisy exit but standard messages from
C                 P01ABF are suppressed
C     IFAIL =  0: hard failure, noisy exit
C
C     For compatibility with certain routines included before Mark 12
C     P01ABF also allows an alternative specification of IFAIL in which
C     it is regarded as a decimal integer with least significant digits
C     cba. Then
C
C     a = 0: hard failure  a = 1: soft failure
C     b = 0: silent exit   b = 1: noisy exit
C
C     except that hard failure now always implies a noisy exit.
C
C     S.Hammarling, M.P.Hooper and J.J.du Croz, NAG Central Office.
C
C     .. Scalar Arguments ..
      INTEGER                 IERROR, IFAIL, NREC
      CHARACTER*(*)           SRNAME
C     .. Array Arguments ..
      CHARACTER*(*)           REC(*)
C     .. Local Scalars ..
      INTEGER                 I, NERR
      CHARACTER*72            MESS
C     .. External Subroutines ..
      EXTERNAL                P01ABZ, X04AAF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC               ABS, MOD
C     .. Executable Statements ..
      IF (IERROR.NE.0) THEN
C        Abnormal exit from calling routine
         IF (IFAIL.EQ.-1 .OR. IFAIL.EQ.0 .OR. IFAIL.EQ.-13 .OR.
     *       (IFAIL.GT.0 .AND. MOD(IFAIL/10,10).NE.0)) THEN
C           Noisy exit
            CALL X04AAF(0,NERR)
            DO 20 I = 1, NREC
               CALL X04BAF(NERR,REC(I))
   20       CONTINUE
            IF (IFAIL.NE.-13) THEN
               WRITE (MESS,FMT=99999) SRNAME, IERROR
               CALL X04BAF(NERR,MESS)
               IF (ABS(MOD(IFAIL,10)).NE.1) THEN
C                 Hard failure
                  CALL X04BAF(NERR,
     *                     ' ** NAG hard failure - execution terminated'
     *                        )
                  CALL P01ABZ
               ELSE
C                 Soft failure
                  CALL X04BAF(NERR,
     *                        ' ** NAG soft failure - control returned')
               END IF
            END IF
         END IF
      END IF
      P01ABF = IERROR
      RETURN
C
99999 FORMAT (' ** ABNORMAL EXIT from NAG Library routine ',A,': IFAIL',
     *  ' =',I6)
      END

      SUBROUTINE P01ABZ
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C
C     Terminates execution when a hard failure occurs.
C
C     ******************** IMPLEMENTATION NOTE ********************
C     The following STOP statement may be replaced by a call to an
C     implementation-dependent routine to display a message and/or
C     to abort the program.
C     *************************************************************
C     .. Executable Statements ..
      STOP
      END

      SUBROUTINE X04AAF(I,NERR)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978
C     MARK 7C REVISED IER-190 (MAY 1979)
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     IF I = 0, SETS NERR TO CURRENT ERROR MESSAGE UNIT NUMBER
C     (STORED IN NERR1).
C     IF I = 1, CHANGES CURRENT ERROR MESSAGE UNIT NUMBER TO
C     VALUE SPECIFIED BY NERR.
C
C     .. Scalar Arguments ..
      INTEGER           I, NERR
C     .. Local Scalars ..
      INTEGER           NERR1
C     .. Save statement ..
      SAVE              NERR1
C     .. Data statements ..
      DATA              NERR1/0/
C     .. Executable Statements ..
      IF (I.EQ.0) NERR = NERR1
      IF (I.EQ.1) NERR1 = NERR
      RETURN
      END

      SUBROUTINE X04BAF(NOUT,REC)
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C
C     X04BAF writes the contents of REC to the unit defined by NOUT.
C
C     Trailing blanks are not output, except that if REC is entirely
C     blank, a single blank character is output.
C     If NOUT.lt.0, i.e. if NOUT is not a valid Fortran unit identifier,
C     then no output occurs.
C
C     .. Scalar Arguments ..
      INTEGER           NOUT
      CHARACTER*(*)     REC
C     .. Local Scalars ..
      INTEGER           I
C     .. Intrinsic Functions ..
      INTRINSIC         LEN
C     .. Executable Statements ..
      IF (NOUT.GE.0) THEN
C        Remove trailing blanks
         DO 20 I = LEN(REC), 2, -1
            IF (REC(I:I).NE.' ') GO TO 40
   20    CONTINUE
C        Write record to external file
   40    WRITE (NOUT,FMT=99999) REC(1:I)
      END IF
      RETURN
C
99999 FORMAT (A)
      END



