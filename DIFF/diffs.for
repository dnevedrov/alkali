	program diff

      INTEGER          NOUT
      PARAMETER        (NOUT=6)
*     .. Local Scalars ..
      DOUBLE PRECISION HBASE, XVAL
      INTEGER          I, IFAIL, J, K, L, NDER, I1
*     .. Local Arrays ..
      DOUBLE PRECISION DER(14), EREST(14)
*     .. External Functions ..
      DOUBLE PRECISION FUN, F, Y
      EXTERNAL         FUN
*     .. External Subroutines ..
      EXTERNAL         D04AAF
*     .. Intrinsic Functions ..
      INTRINSIC        ABS
*     .. Executable Statements ..
      OPEN(1,FILE='PLOT.DAT')
	DO I1 = 1, 1000
	   Y = 1.5D-10+4.0D-10/1000*I1
	F = FUN(Y)

         WRITE(1,99)Y,F
	ENDDO
	CLOSE(1)
99    FORMAT(2E16.8)
      HBASE = 2.0D-12
      NDER = 8
      XVAL = 2.82D-10
      X = XVAL


C      DO 40 K = 1, 1
         IFAIL = 0
C        
         CALL D04AAF(XVAL,NDER,HBASE,DER,EREST,FUN,IFAIL)
C
         WRITE (NOUT,*)
         WRITE (NOUT,99999) 'with step length', HBASE,
     *     '  the results are'
         WRITE (NOUT,*) 'Order        Derivative       Error estimate'
         DO 20 I = 1, 8
            WRITE (NOUT,99998) I, DER(I), EREST(I)

   20    CONTINUE
         HBASE = HBASE*0.1D0
   40 CONTINUE
    
C	PRINT*,DER(3)**2/(DER(2)*DER(4))
C      DO I=1,20
C C     Y=FUN(XVAL+0.1 D-10*I)
C	ENDDO
      STOP

99999 FORMAT (1X,A,F9.4,A)
99998 FORMAT (1X,I2,2D21.4)
      END

      DOUBLE PRECISION FUNCTION FUN(X)
c     .. Scalar Arguments ..
      DOUBLE PRECISION DAT(5,16),X0,C,D,B,K,X
	DOUBLE PRECISION Z1, Z2, EL, EPS0, DIFF2, DIFF3, DIFF4, PI
	INTEGER M,N,I,J,L
      CHARACTER*4 CR(16), CRYST
      DATA CR /'LIF', 'LICL', 'LIBR', 'LII',
     *         'NAF', 'NACL', 'NABR', 'NAI',
     *         'KF', 'KCL', 'KBR', 'KI',
     *         'RBF', 'RBCL', 'RBBR', 'RBI'/
	DATA ((DAT(J,I),J=1,5),I=1,16)/
     *2.014d-10,7.09D-79,29.7D-99,1.15D-58, 1.12D10,
     *2.57d-10,17.4D-79,75.5D-99,1.78D-58, 1.38D10,
     *2.751d-10,20.7D-79,93.1D-99,1.99D-58, 1.42D10,
     *3.0d-10,25.0D-79,114.0D-99,2.50D-58, 1.13D10,
     *2.317d-10,17.1D-79,45.0D-99,1.52D-58, 1.23D10,
     *2.82d-10,41.8D-79,117.0D-99,2.26D-58, 1.58D10,
     *2.989d-10,49.9D-79,147.0D-99,2.59D-58, 1.47D10,
     *3.237d-10,59.7D-79,183.0D-99,2.88D-58, 1.56D10,
     *2.674d-10,45.1D-79,137.0D-99,2.43D-58, 1.37D10,
     *3.147d-10,110.0D-79,353.0D-99,3.52D-58, 1.44D10,
     *3.298d-10,132.0D-79,442.0D-99,3.77D-58, 1.52D10,
     *3.533d-10,158.0D-79,548.0D-99,4.16D-58, 1.51D10,
     *2.815d-10,62.8D-79,203.0D-99,2.92D-58, 1.39D10,
     *3.291d-10,153.0D-79,520.0D-99,3.99D-58, 1.59D10,
     *3.445d-10,184.0D-79,650.0D-99,4.38D-58, 1.55D10,
     *3.671d-10,221.0D-79,806.0D-99,4.68D-58, 1.67D10/
   
     	CRYST = 'NACL'
	
	DO L=1,14
	  IF (CRYST.EQ.CR(L)) THEN 

	  GOTO 30
	ENDIF
	ENDDO

30	M=1
	N=4
	Z1 = 1
	Z2 = -1
	PI = 3.1415926D0
	EL = 1.602176462D-19
	EPS0 = 8.854187817D-12
	X0 = DAT(1,L)
	C = DAT(2,L)
	D = DAT(3,L)
	B = DAT(4,L)
	K = DAT(5,L)
C	PRINT*,X0,C,D,B,K
C	DIFF2 = -72.0D0*D/X**10 - 42.0D0*C/X**8 +
C     *      8.0D0*B*DEXP(-K*(X-X0))*K/X**5 + 
C     *     EL**2*Z1*Z2/(2.0D0*EPS0*PI*X**3)+
C     *      20.0D0*B*DEXP(-K*(X-X0))/X**6+
C     *    B*DEXP(-K*(X-X0))*K**2/X**4 C
C	 PRINT*, x,DIFF2
C	DIFF3 = -720D0*D/X**11 + 336D0*C/X**9 -
C     *      60.0D0*B*DEXP(K)/X**6 + 3.0D0*EL**2/(2.0D0*EPS0*PI*X**4)-
C     *      20.0D0*B*DEXP(K*(X-X0))/X**6

C	DIFF4 = -7920.0D0*D/X**12 + 3024.0D0*C/X**10 -
C     *      480.0D0*B*DEXP(K)/X**7 - 6.0D0*EL**2/(EPS0*PI*X**5)-
C     *      840.0D0*B*DEXP(K*(X-X0))/X**8
C	PRINT*, DIFF2, DIFF3, DIFF4
c     .. Intrinsic Functions ..
C      INTRINSIC                     EXP
c     .. Executable Statements ..
      FUN = Z1*Z2*EL*EL/(4.0D0*3.1415926D0*EPS0*X)+
     *     B/(X**N)*DEXP(-K*(X**M-X0**M))-C/(X**6)-D/(X**8)    
C       PRINT*,X
      RETURN
      END

	SUBROUTINE D04AAF(XVAL,NDER,HBASE,DER,EREST,FUN,IFAIL)
C
C     ***   PURPOSE   ***
C
C     A SUBROUTINE FOR NUMERICAL DIFFERENTIATION AT A POINT.  IT
C     RETURNS A SET OF APPROXIMATIONS TO THE J-TH ORDER DERIVATIVE
C     (J=1,2,...14) OF FUN(X) EVALUATED AT X = XVAL  AND, FOR EACH
C     DERIVATIVE, AN ERROR ESTIMATE ( WHICH INCLUDES THE EFFECT OF
C     AMPLIFICATION OF ROUND- OFF ERRORS).
C
C
C     ***   INPUT  PARAMETERS   ***
C
C     (1) XVAL   REAL.  THE ABSCISSA AT WHICH THE SET OF
C     DERIVATIVES IS REQUIRED.
C     (2) NDER   INTEGER. THE HIGHEST ORDER DERIVATIVE REQUIRED.
C     IF(NDER.GT.0)  ALL DERIVATIVES UP TO MIN(NDER,14) ARE
C     CALCULATED.
C     IF(NDER.LT.0 AND NDER EVEN)  EVEN ORDER DERIVATIVES
C     UP TO MIN(-NDER,14) ARE CALCULATED.
C     IF(NDER.LT.0 AND NDER ODD )  ODD  ORDER DERIVATIVES
C     UP TO MIN(-NDER,13) ARE CALCULATED.
C     (3) HBASE  REAL.  A STEP LENGTH.
C     (6) FUN    THE NAME OF A REAL FUNCTION SUBPROGRAMME,
C     WHICH IS REQUIRED BY THE ROUTINE AS A SUBPROGRAMME
C     AND WHICH REPRESENTS THE FUNCTION BEING DIFFERENTIATED
C     THE ROUTINE REQUIRES 21 FUNCTION EVALUATIONS FUN(X)
C     LOCATED AT    X = XVAL    AND AT
C     X = XVAL + (2*J-1)*HBASE,  J=-9,-8, .... +9,+10.
C     THE FUNCTION VALUE AT  X = XVAL IS DISREGARDED WHEN
C     ONLY ODD ORDER DERIVATIVES ARE REQUIRED.
C
C     ***   OUTPUT PARAMETERS   ***
C
C     (4) DER(J)   J=1,2,...14.  REAL. A VECTOR OF
C     APPROXIMATIONS TO THE J-TH DERIVATIVE OF FUN(X)
C     EVALUATED AT X = XVAL.
C     (5) EREST(J) J=1,2,...14.  REAL. A VECTOR OF
C     ESTIMATES OF THE DABSOLUTE ACCURACY OF DER(J).  THESE
C     ARE NEGATIVE WHEN EREST(J).GT.ABS(DER(J)), OR WHEN,
C     FOR SOME OTHER REASON THE ROUTINE IS DOUBTFUL ABOUT
C     THE VALIDITY OF THE RESULT.
C
C     ***   IMPORTANT WARNING   ***
C
C     EREST IS AN ESTIMATE OF THE OVERALL ERROR.  IT IS ESSENTIAL
C     FOR
C     PROPER USE THAT THE USER CHECKS EACH VALUE OF DER
C     SUBSEQUENTLY
C     USED TO SEE THAT IT IS ACCURATE ENOUGH FOR HIS PURPOSES.
C     FAILURE TO DO THIS MAY RESULT IN THE CONTAMINATION OF ALL
C     SUBSEQUENT RESULTS.
C     IT IS TO BE EXPECTED THAT IN NEARLY ALL CASES  DER(14) WILL
C     BE
C     UNUSABLE.( 14 REPRESENTS A LIMIT IN WHICH FOR THE EASIEST
C     FUNCTION LIKELY TO BE ENCOUNTERED THIS ROUTINE MIGHT JUST
C     OBTAIN
C     AN APPROXIMATION OF THE CORRECT SIGN.)
C
C     ***   NOTE ON CALCULATION   ***
C
C     THE CALCULATION IS BASED ON THE EXTENDED T-TABLE (T SUB
C     K,P,S)
C     DESCRIBED IN LYNESS AND MOLER, NUM. MATH., 8, 458-464,(1966).
C     REFERENCES IN COMMENT CARDS TO TTAB(NK,NP,NS) REFER TO
C     T-TABLE
C     WITH NK=K+1,  NP=P+1,  NS=S+1.   SINCE ONLY PART OF THE
C     EXTENDED
C     T-TABLE IS NEEDED AT ONE TIME, THAT PART IS STORED IN
C     RTAB(10,7)
C     AND IS SUBSEQUENTLY OVERWRITTEN.  HERE
C     RTAB(NK,NP-NS+1)  =  TTAB(NK,NP,NS).
C
C     NAG COPYRIGHT 1976
C     MARK 5 RELEASE
C     MARK 6B REVISED  IER-116 (MAR 1978)
C     MARK 7 REVISED IER-139 (DEC 1978)
C     MARK 8 REVISED. IER-219 (MAR 1980)
C     MARK 8D REVISED. IER-271 (DEC 1980).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D04AAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  HBASE, XVAL
      INTEGER           IFAIL, NDER
C     .. Array Arguments ..
      DOUBLE PRECISION  DER(14), EREST(14)
C     .. Function Arguments ..
      DOUBLE PRECISION  FUN
      EXTERNAL          FUN
C     .. Local Scalars ..
      DOUBLE PRECISION  A, B, BIG, ERNOW, ERPREV, FACT, FTEST, FZ, H,
     *                  HFACT, HSQ, HTEST, ONE, RANMIN, SMALL, TEMP,
     *                  TERM, THRON2, TWO, XJNS, XMAX, XMIN, XNUM, ZERO
      INTEGER           J, JNS, N, NDERP, NDPAR, NEGER, NK, NKTOP, NLO,
     *                  NMAX, NMIN, NN1, NP, NPAR, NPMIN, NPR, NS,
     *                  NSMAX, NSQ, NTOP, NTOPM2, NUP
C     .. Local Arrays ..
      DOUBLE PRECISION  ACOF(10), RANGE(7), RTAB(10,7), TZEROM(10)
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AKF, X02ALF
      INTEGER           P01ABF
      EXTERNAL          X02AKF, X02ALF, P01ABF
C     .. Data statements ..
      DATA              ZERO, ONE, TWO, THRON2/0.0D0, 1.0D0, 2.0D0,
     *                  1.5D0/
C     .. Executable Statements ..
C
C     TO ALTER THE PRECISION OF THE WHOLE ROUTINE,ALTER THE
C     PRECISION
C     IN THE ABOVE DECLARATION AND DATA STATEMENTS
C
C     QUIT ON ZERO HBASE AND ZERO NDER.  SET CONTROL PARAMETER
C     NDPAR.
C     NDPAR = +1  ODD-ORDER DERIVATIVES ONLY.
C     NDPAR =  0  ALL DERIVATIVES
C     NDPAR = -1  EVEN-ORDER DERIVATIVES ONLY.
      BIG = X02ALF()
      SMALL = 1000.0D0*X02AKF()
      IF (HBASE.EQ.ZERO) GO TO 20
      NDERP = NDER
      NDPAR = 0
      IF (NDER) 40, 20, 60
   20 IFAIL = P01ABF(IFAIL,1,SRNAME,0,P01REC)
      RETURN
   40 NDPAR = 4*(NDER/2) - 2*NDER - 1
      NDERP = -NDER
   60 CONTINUE
      IF (NDERP.GT.14) NDERP = 14
C
C     NEXT, EVALUATE 21 FUNCTION VALUES , SET ACOF(N), AND SET
C     FIRST
C     COLUMN OF RTAB FOR ODD ORDER DERIVATIVES AND STORE IN TZEROM
C     THE
C     FIRST COLUMN OF RTAB FOR EVEN ORDER DERIVATIVES.
      FZ = FUN(XVAL)
      DO 80 N = 1, 10
         NN1 = 2*N - 1
         NSQ = NN1*NN1
         ACOF(N) = NSQ
         XNUM = NN1
         H = HBASE*XNUM
         TEMP = FUN(XVAL+H)
         TERM = FUN(XVAL-H)
         RTAB(N,1) = (TEMP-TERM)/(TWO*H)
         TZEROM(N) = (TEMP-FZ-FZ+TERM)/(TWO*H**2)
   80 CONTINUE
C
C     SET UP FOR ODD-ORDER DERIVATIVES.
      IF (NDPAR.EQ.-1) GO TO 100
      NPAR = 1
      NSMAX = (NDERP+1)/2
      IF (NSMAX.LE.0) GO TO 440
      GO TO 140
C
C     SET UP FOR EVEN-ORDER DERIVATIVES.
  100 CONTINUE
      NPAR = -1
      IF (NDPAR.EQ.+1) GO TO 520
      NSMAX = NDERP/2
      IF (NSMAX.LE.0) GO TO 520
C
C     SET THE FIRST COLUMN OF THE T-TABLE FOR EVEN ORDER
C     DERIVATIVES.
      DO 120 J = 1, 10
         RTAB(J,1) = TZEROM(J)
  120 CONTINUE
  140 CONTINUE
C
C     ODD-ORDER AND EVEN ORDER PATHS JOIN UP HERE.
C
      NEGER = 0
      HSQ = HBASE*HBASE
      FACT = ONE
      HFACT = ONE
      DO 420 NS = 1, NSMAX
C
C        FROM HERE ON TO   STATEMENT 360 (NEAR END)  WE ARE DEALING
C        WITH A
C        SPECIFIED VALUE OF NS.
C
C        FOR EACH VALUE OF NP WE CALCULATE RANGE(NP)  WHICH IS THE
C        DIFFERE
C        BETWEEN THE GREATEST AND THE LEAST OF THE ELEMENTS
C        TTAB(NK,NP,NS)
C        AS WE GO ALONG WE DETERMINE ALSO THE MINIMUM OF RANGE(NP)
C        WHICH I
C        RANMIN = RANGE(NPMIN).
C        WE RETAIN NUP AND NLO WHICH ARE THE VALUES  OF NK GIVING THE
C        EXTREME VALUES OF TTAB(NK,NPMIN,NS).
C        THIS PART OF NS LOOP CONCLUDES AT STATEMENT NUMBER 280.
C
C
C        FIRST CALCULATE ELEMENTS OF T-TABLE  FOR CURRENT NS VALUE.
         NPMIN = NS
         IF (NS.EQ.1) NPMIN = 2
         DO 180 NP = NPMIN, 7
            NKTOP = 10 - NP + 1
            NPR = NP - NS + 1
            DO 160 NK = 1, NKTOP
               J = NP + NK - 1
               A = ACOF(J)
               B = ACOF(NK)
               TERM = ZERO
               TEMP = ZERO
               IF (NP.NE.NS) TEMP = A*RTAB(NK,NPR-1) -
     *                              B*RTAB(NK+1,NPR-1)
               IF (NS.NE.1) TERM = -RTAB(NK,NPR) + RTAB(NK+1,NPR)
               RTAB(NK,NPR) = (TERM+TEMP)/(A-B)
  160       CONTINUE
  180    CONTINUE
C
C        NOW CALCULATE NUP,NLO AND RANMIN.
         DO 280 NP = NS, 7
            NTOP = 11 - NP
            NPR = NP - NS + 1
            XMAX = RTAB(1,NPR)
            XMIN = RTAB(1,NPR)
            NMAX = 1
            NMIN = 1
C
            DO 200 NK = 2, NTOP
               TEMP = RTAB(NK,NPR)
               IF (TEMP.GT.XMAX) NMAX = NK
               IF (TEMP.GT.XMAX) XMAX = TEMP
               IF (TEMP.LT.XMIN) NMIN = NK
               IF (TEMP.LT.XMIN) XMIN = TEMP
  200       CONTINUE
C
            RANGE(NP) = XMAX - XMIN
            IF (NP.NE.NS) GO TO 220
            RANMIN = RANGE(NP)
            GO TO 240
  220       CONTINUE
            IF (RANGE(NP).GE.RANMIN) GO TO 260
            RANMIN = RANGE(NP)
  240       CONTINUE
            NPMIN = NP
            NUP = NMAX
            NLO = NMIN
            IF (NLO.EQ.NUP) NLO = NLO + 1
  260       CONTINUE
  280    CONTINUE
C
C        NEXT WE TAKE THE AVERAGE OF ALL EXCEPT THE EXTREME VALUES OF
C        TTAB(NK,NPMIN,NS) AS THE RESULT AND RANMIN AS THE ERROR
C        ESTIMATE.
         TERM = ZERO
         NTOP = 11 - NPMIN
         NTOPM2 = NTOP - 2
         XNUM = NTOPM2
         J = NPMIN - NS + 1
         DO 320 NK = 1, NTOP
            IF (NK.EQ.NUP) GO TO 300
            IF (NK.EQ.NLO) GO TO 300
            TERM = TERM + RTAB(NK,J)
  300       CONTINUE
  320    CONTINUE
         TERM = TERM/XNUM
C
C        THE ABOVE RESULT AND ERROR ESTIMATE REFER TO TAYLOR
C        COEFFICIENT.
C        NEXT WE DO SCALING TO OBTAIN DERIVATIVE RESULTS INSTEAD.
C        JNS AND XJNS ARE ACTUAL ORDER OF DERIVATIVE BEING TREATED.
C        SAFETY FACTORS 1.5,1.5,2.0,2.0 AND 2.0 ARE MULTIPLIED INTO
C        THE
C        ERROR ESTIMATES FOR J = 10,11,12,13 AND 14 RESPECTIVELY.
C        THESE
C        HAVE BEEN CHOSEN BY INSPECTION OF PERFORMANCE STATISTICS.
C        THERE
C        IS NO ANALYTIC JUSTIFICATION FOR THESE PARTICULAR FACTORS.
         JNS = 2*NS - 1
         IF (NPAR.LT.0) JNS = 2*NS
         XJNS = JNS
         FACT = FACT*XJNS
C        TEST FOR UNDERFLOW OF HFACT
         IF (HFACT.GE.1.0D0) GO TO 340
         HTEST = HFACT*BIG
         FTEST = TERM*FACT
         IF ((TWO*RANMIN*FACT).GT.FTEST) FTEST = TWO*RANMIN*FACT
         IF (HTEST.GT.FTEST) GO TO 340
         DER(JNS) = BIG
         EREST(JNS) = -DER(JNS)
         NEGER = NEGER + 1
         GO TO 360
C        END OF CODE TO HANDLE ZERO HFACT
  340    CONTINUE
         DER(JNS) = TERM*FACT/HFACT
         IF (JNS.EQ.10) RANMIN = THRON2*RANMIN
         IF (JNS.EQ.11) RANMIN = THRON2*RANMIN
         IF (JNS.GE.12) RANMIN = TWO*RANMIN
         EREST(JNS) = RANMIN*FACT/(HFACT)
C
C        SET SIGN OF EREST.  EREST NEGATIVE EITHER IF IT SWAMPS DER,
C        OR IF
C        TWO PREVIOUS CONSECUTIVE ERESTS OF THIS PARITY ARE NEGATIVE.
C        IT MAY ALSO BE SET NEGATIVE AT END (750).
         IF (NEGER.GE.2) EREST(JNS) = -EREST(JNS)
         IF (NEGER.GE.2) GO TO 360
         IF (TERM.LT.RANMIN) EREST(JNS) = -EREST(JNS)
         IF (-TERM.GE.RANMIN) EREST(JNS) = -EREST(JNS)
         IF (EREST(JNS).LT.ZERO) NEGER = NEGER + 1
         IF (EREST(JNS).GE.ZERO) NEGER = 0
  360    CONTINUE
C
C        SECOND TEST FOR UNDERFLOW OF HFACT
         IF (HFACT.GE.1.0D0) GO TO 380
         IF (HFACT.EQ.0.0D0) GO TO 400
         IF (HSQ.GE.SMALL/HFACT) GO TO 380
         HFACT = 0.0D0
         GO TO 400
  380    HFACT = HSQ*HFACT
  400    FACT = FACT*(XJNS+ONE)
  420 CONTINUE
C
  440 CONTINUE
      IF (NPAR.GT.0) GO TO 100
C
C     SET SIGN OF EREST NEGATIVE IF TWO PREVIOUS CONSECUTIVE ERESTS
C     ARE NEGATIVE.
      IF (NDPAR.NE.0) GO TO 520
      IF (NDERP.LE.2) GO TO 520
      NEGER = 0
      ERPREV = EREST(1)
      DO 500 J = 2, NDERP
         ERNOW = EREST(J)
         IF (NEGER.EQ.2) GO TO 460
         IF (ERPREV.GE.ZERO) GO TO 480
         IF (ERNOW.GE.ZERO) GO TO 480
  460    NEGER = 2
         IF (ERNOW.GT.ZERO) EREST(J) = -ERNOW
  480    ERPREV = ERNOW
  500 CONTINUE
C
  520 CONTINUE
      IFAIL = 0
      RETURN
C     END OF D04AAF   ***   ***   ***   ***   ***
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
      DOUBLE PRECISION FUNCTION X02AKF()
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS  B**(EMIN-1)  (THE SMALLEST POSITIVE MODEL NUMBER)
C
      DOUBLE PRECISION X02CON
      DATA X02CON /Z'0010000000000000' /
C     .. Executable Statements ..
      X02AKF = X02CON
      RETURN
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
      DOUBLE PRECISION FUNCTION X02ALF()
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS  (1 - B**(-P)) * B**EMAX  (THE LARGEST POSITIVE MODEL
C     NUMBER)
C
      DOUBLE PRECISION X02CON
      DATA X02CON /Z'7FEFFFFFFFFFFFFF' /
C     .. Executable Statements ..
      X02ALF = X02CON
      RETURN
      END
