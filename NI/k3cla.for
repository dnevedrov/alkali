      PROGRAM KCL

      INTEGER K,J,I,NT,P01ABF,IFAIL,N1,N2,NL,M,NM,NW,NT1
      PARAMETER (NW=1375,NT=1000,NM=9)
      COMPLEX*16 DC,G(36,NW),G1(9,9),G2(9,9),G3(9,9),G4(9,9),INVERSE,
     *   GG(9,9),TMP1(9,9),TMP2(9,9),TMP3(9,9),V(9,9),P1(9,9),P2(9,9),
     *   AL,BE,SPUR,V1(3,3),V2(3,3),V3(3,3),W(9,9),W2(3,3),W1(3,3),
     *   W3(3,3),IMG(9,9),VIG1(9,9,NW),VIG4(9,9,NW)
      DOUBLE PRECISION OMEGA(NW),ER,RESULT,F1(NW),E(NT),HBAR,PI,WL(NT),
     *    X01AAF,T(NT),DT,GAMMA(NT),B,K1,K2,K3,K4,AMP(NT),INC,X,
     *    AU,DK2,MD,MCL,MK,WMAX,R0,GAM,KK3,Y, xtemp, wmax1, wmax2,
     *    k20, k30, k40  
      CHARACTER*1 TR
      CHARACTER*3 DIREC
C     ------------------------------------------------
      DIREC='100'
      NT1=1
C     ------------------------------------------------ 
      AU=1.6605655D-27
      TR='N'
      ZETA=0.0D0
      AL=DCMPLX(1.0D0,0.0D0)
      BE=DCMPLX(0.0D0,0.0D0)
      
	HBAR = 6.626d-34
      PI=3.1415926d0
      R0=3.237D-10
      WMAX=1.4336D13
	wmax1 = 2.1736d13
	wmax2 = 3.19d13
      K1= -0.1652D-8
	K2 = 0.9813D2
	K3 = -0.3113D13
	K4 = 0.9232D23
      DK2=0.0D-4
      MD=23.0D0*AU
      MCL=127.0D0*AU
      MK=23.0D0*AU
		K20 = K1/R0
	K30 = (K2 - K20)/R0
	K40 = K3/R0 - (K2-K20)/(2.0D0*R0*R0)
	IF (DYREC.EQ.'100') THEN
	K3 = K3
	K4 = K4
	ELSEIF (DYREC.EQ.'110') THEN
		K3 = K3 + K30
	K4 = K4 + K40
	ELSEIF (DYREC.EQ.'111') THEN
	K3 = K3 + 2.0D0*K30
	K4 = K4 + 2.0D0*K40
	endif
      KK3=(K3/2.0D0)/MK
C     ------------------------------------------------  
      DO J=1,9
         DO I=1,9
            V(J,I)=DCMPLX(0.0D0,0.0D0)
         ENDDO
      ENDDO
      DO J=1,3
         DO I=1,3
            V1(J,I)=DCMPLX(0.0D0,0.0D0)
            V2(J,I)=DCMPLX(0.0D0,0.0D0)
            V3(J,I)=DCMPLX(0.0D0,0.0D0)
         ENDDO
      ENDDO                
      IF (DIREC.EQ.'100') THEN 
			nfac = 4.0d0
	sfac = -4.004d0
         AMP(1)=0.60D-10 
         V1(1,2)=DCMPLX(-DSQRT(MK/MD),0.0D0)
         V1(1,3)=DCMPLX(DSQRT(MK/MD),0.0D0)
         V1(2,1)=DCMPLX(-DSQRT(MK/MD),0.0D0)
         V1(3,1)=DCMPLX(DSQRT(MK/MD),0.0D0)
         V1(2,2)=DCMPLX(1.0D0,0.0D0)
         V1(3,3)=DCMPLX(1.0D0,0.0D0)
      ELSEIF (DIREC.EQ.'110') THEN
	nfac = 8.0d0
	sfac = -3.527d0
         AMP(1)=0.80D-10      
         V2(1,2)=DCMPLX(-DSQRT(MK/MD),0.0D0)
         V2(1,3)=DCMPLX(DSQRT(MK/MD),0.0D0)
         V2(2,1)=DCMPLX(-DSQRT(MK/MD),0.0D0)
         V2(3,1)=DCMPLX(DSQRT(MK/MD),0.0D0)
         V2(2,2)=DCMPLX(1.0D0,0.0D0)
         V2(3,3)=DCMPLX(1.0D0,0.0D0)
      ELSEIF (DIREC.EQ.'111') THEN
	sfac = -3.0507d0
	nfac = 12.0d0
         AMP(1)=1.0D-10
         V3(1,2)=DCMPLX(-DSQRT(MK/MD),0.0D0)
         V3(1,3)=DCMPLX(DSQRT(MK/MD),0.0D0)
         V3(2,1)=DCMPLX(-DSQRT(MK/MD),0.0D0)
         V3(3,1)=DCMPLX(DSQRT(MK/MD),0.0D0)
         V3(2,2)=DCMPLX(1.0D0,0.0D0)
         V3(3,3)=DCMPLX(1.0D0,0.0D0)
      ENDIF
      DO J=1,3
         DO I=1,3
            V(J,I)=V1(J,I)
            V(J+3,I+3)=V2(J,I)
            V(J+6,I+6)=V3(J,I)
         ENDDO
      ENDDO
C     ================================================
      CALL KIRGUDU(OMEGA,G)
      DO J=1,9
         DO I=1,9
            W(J,I)=DCMPLX(0.0D0,0.0D0)
         ENDDO
      ENDDO
C     ================================================
      T(1)=0.0D0 
C     ================================================    
      DO M=1,NT
         X=AMP(M)*1.0D10
	if (x.eq.0.0) then 
	amp(m) = xtemp
	x = xtemp*1.0d10
	endif

         IF (DIREC.EQ.'100') THEN
	  y = 2.5799d0*x**4 - 7.7743d0*x**3 + 8.7701d0*x**2 - 4.4474d0*x + 
     *	  1.6719d0
            WL(M)=Y*WMAX1
	
	
         ELSEIF (DIREC.EQ.'110') THEN
	y = -1.9193d0*x**4 + 4.453d0*x**3 - 2.905d0*x**2 - 0.0709d0*x 
     * 	+ 1.2522d0
            WL(M)=Y*WMAX1
         ELSEIF (DIREC.EQ.'111') THEN
	 y = -8.2993d0*x**4 + 26.137d0*x**3 - 29.915d0*x**2 + 14.351d0*x
     *	  - 1.4313d0

            WL(M)=Y*WMAX1

         ENDIF

         NT1=M
         IF (WL(M)/WMAX1.GT.0.99d0) GOTO 20
         E(M)=(AMP(M)*WL(M))**2*MD/2.0D0
c     *      +  (1.0D0/64.0D0)*K4*AMP(M)**4.
c         N1=(-(WL(M)-WMAX)/WMAX)*(NW/2)
		n1 = (wl(m)-wmax)/wmax*308
	  
c         N2=(NW/2)
		n2 = wmax/wmax*308
c	 print*,'n1 = ', n1,n2
         NL=WL(M)/WMAX*308
	   GAM=DK2/MCL+(1.0D0/nfac)*K4*AMP(M)**2/MCL*sfac
         DO K=N1,N2            
c            GAM=DK2/MCL+(1.0D0/nfac)*K4*AMP(M)**2/MCL
            B=(1.0D0-MD/MCL)*OMEGA(K)**2+2.0D0*GAM   
            W1(1,1)=DCMPLX(B,0.0D0)
            W1(2,2)=DCMPLX(GAM,0.0D0)
            W1(3,3)=DCMPLX(GAM,0.0D0)
            W1(1,2)=DCMPLX(-GAM,0.0D0)
            W1(1,3)=DCMPLX(-GAM,0.0D0)
            W1(2,1)=DCMPLX(-GAM,0.0D0) 
            W1(2,3)=DCMPLX(0.0D0,0.0D0)
            W1(3,1)=DCMPLX(-GAM,0.0D0)
            W1(3,2)=DCMPLX(0.0D0,0.0D0)
            W2(1,1)=DCMPLX(B,0.0D0)
            W2(2,2)=DCMPLX(GAM,0.0D0)
            W2(3,3)=DCMPLX(GAM,0.0D0)
            W2(1,2)=DCMPLX(-GAM,0.0D0)
            W2(1,3)=DCMPLX(-GAM,0.0D0)
            W2(2,1)=DCMPLX(-GAM,0.0D0) 
            W2(2,3)=DCMPLX(0.0D0,0.0D0)
            W2(3,1)=DCMPLX(-GAM,0.0D0)
            W2(3,2)=DCMPLX(0.0D0,0.0D0)
            W3(1,1)=DCMPLX(B,0.0D0)
            W3(2,2)=DCMPLX(GAM,0.0D0)
            W3(3,3)=DCMPLX(GAM,0.0D0)
            W3(1,2)=DCMPLX(-GAM,0.0D0)
            W3(1,3)=DCMPLX(-GAM,0.0D0)
            W3(2,1)=DCMPLX(-GAM,0.0D0) 
            W3(2,3)=DCMPLX(0.0D0,0.0D0)
            W3(3,1)=DCMPLX(-GAM,0.0D0)
            W3(3,2)=DCMPLX(0.0D0,0.0D0)
            DO J=1,3
               DO I=1,3
                  W(J,I)=W1(J,I)
                  W(J+3,I+3)=W2(J,I)
                  W(J+6,I+6)=W3(J,I)
               ENDDO
            ENDDO                              
C     ----------------------------------------------------------------
	
            CALL TIGR(GG,G,K,W)
            DO J=1,NM
               DO I=1,NM
                  VIG1(J,I,K)=GG(J,I)
                  VIG4(J,I,K)=DCMPLX(DREAL(GG(J,I)),-DIMAG(GG(J,I)))
               ENDDO
            ENDDO            
         ENDDO
	
         DO K=N1,N2
            DO J=1,NM
               DO I=1,NM
                  G1(J,I)=VIG1(J,I,K)
                  G4(J,I)=VIG4(J,I,K)
                  G3(J,I)=VIG1(J,I,NL-K)
                  G2(J,I)=VIG4(J,I,NL-K)
               ENDDO
            ENDDO           
C     ...P(w)............................................................. 
            CALL F06ZAF(TR,TR,NM,NM,NM,AL,V,NM,G2,NM,BE,TMP1,NM)
            CALL F06ZAF(TR,TR,NM,NM,NM,AL,TMP1,NM,V,NM,BE,TMP2,NM)
            CALL F06ZAF(TR,TR,NM,NM,NM,AL,TMP2,NM,G1,NM,BE,TMP1,NM)
            DO J=1,NM
               DO I=1,NM
                  TMP2(J,I)=DC(J,I)-
     *                 (AMP(M)*KK3)**2*TMP1(J,I)
                  IMG(J,I)=DCMPLX(DIMAG(G1(J,I)),0.0D0)
               ENDDO
            ENDDO
            CALL INVERSE1(TMP2)
            CALL F06ZAF(TR,TR,NM,NM,NM,AL,V,NM,IMG,NM,BE,TMP1,NM)         
            CALL F06ZAF(TR,TR,NM,NM,NM,AL,TMP1,NM,TMP2,NM,BE,P1,NM)         
C     ----P(w_L - w)-------------------------------------------------------
            CALL F06ZAF(TR,TR,NM,NM,NM,AL,V,NM,G4,NM,BE,TMP1,NM)
            CALL F06ZAF(TR,TR,NM,NM,NM,AL,TMP1,NM,V,NM,BE,TMP2,NM)
            CALL F06ZAF(TR,TR,NM,NM,NM,AL,TMP2,NM,G3,NM,BE,TMP1,NM)
            DO J=1,NM
               DO I=1,NM
                  TMP2(J,I)=DC(J,I)-
     *                 (AMP(M)*KK3)**2*TMP1(J,I)
                  IMG(J,I)=DCMPLX(DIMAG(G3(J,I)),0.0D0)
               ENDDO
            ENDDO
            CALL INVERSE1(TMP2)
            CALL F06ZAF(TR,TR,NM,NM,NM,AL,V,NM,IMG,NM,BE,TMP1,NM)
            CALL F06ZAF(TR,TR,NM,NM,NM,AL,TMP1,NM,TMP2,NM,BE,P2,NM)
C     ----------------------------------------------------------------
            CALL F06ZAF(TR,TR,NM,NM,NM,AL,P1,NM,P2,NM,BE,TMP1,NM)
            F1(K)=DREAL(SPUR(TMP1))
       ENDDO
C     ================================================================
         DO K=1,N1-1
            F1(K)=0.0D0
         ENDDO
         DO K=N2+1,NW
            F1(K)=0.0D0
         ENDDO
	
         CALL D01GAF(OMEGA,F1,NW,RESULT,ER,IFAIL)
         GAMMA(M)=(RESULT)*(HBAR*WL(M)*0.25D0/PI)*(AMP(M)*KK3)**2
C         E(M+1)=E(M)*0.9999D0
         IF (GAMMA(M).EQ.0.0D0) THEN
            DT=0.0D0
         ELSEIF (GAMMA(M).GT.0.0D0) THEN
            DT=(E(M)-E(M+1))/GAMMA(M)
         ENDIF
	if (amp(m).eq.0.0) amp(m) = xtemp
         AMP(M+1)=AMP(M)-1.0d-10/1000.0d0
	if (amp(m+1).gt.0.0d0) xtemp= amp(m+1)

c         AMP(M+1)=DSQRT(2.0D0*E(M+1)/(MD*WL(M)**2))
         IF (M.LE.NT-1) T(M+1)=(T(M)+DT*WMAX/(2.0D0*PI))
         PRINT*,M,AMP(M),gamma(m)/e(m)/wmax2
20    ENDDO
C     ================================================================      
      CALL RYBA(NT1,T,GAMMA,E,AMP,DIREC)
C     ================================================================      
      END
      
      COMPLEX*16 FUNCTION DC(J,I)
      INTEGER J,I
      DC=DCMPLX(0.0D0,0.0D0)
      IF (J.EQ.I) DC=DCMPLX(1.0D0,0.0D0)
      END
      
      COMPLEX*16 FUNCTION SPUR(P)
      INTEGER J
      COMPLEX*16 P(9,9)
      SPUR=DCMPLX(0.0D0,0.0D0)
      DO J=1,9
         SPUR=SPUR+P(J,J)
      ENDDO
      END
      
      SUBROUTINE RYBA(NT,T,G,E,A,DIREC)
      INTEGER K,NT
      DOUBLE PRECISION T(NT),G(NT),A(NT),E(NT)
      CHARACTER*3 DIREC
      OPEN(1,FILE=DIREC)
      DO K=1,NT
         WRITE(1,3)T(K),G(K)/E(K)/3.19d13,A(K)*1.0d10,E(K)
      ENDDO
      CLOSE(1)
 3    FORMAT(4E16.8)
      RETURN
      END
          
      SUBROUTINE KIRGUDU(OMEGA,G)
      INTEGER J,I,NR
      PARAMETER (NR=1375)
      COMPLEX*16 G(36,NR)
      DOUBLE PRECISION GIM(36,NR),GRE(36,NR),OMEGA(NR)
      OPEN(2,FILE='DATA1.TXT')
      DO I=1,NR
         READ(2,2)OMEGA(I),GIM(1,I),GRE(1,I)
      ENDDO
      CLOSE(2)
      OPEN(3,FILE='DATA2.TXT')
      DO I=1,NR
         READ(3,2)OMEGA(I),GIM(2,I),GRE(2,I)
      ENDDO
      CLOSE(3)
      OPEN(4,FILE='DATA3.TXT')
      DO I=1,NR
         READ(4,2)OMEGA(I),GIM(3,I),GRE(3,I)
      ENDDO
      CLOSE(4)
      OPEN(5,FILE='DATA4.TXT')
      DO I=1,NR
         READ(5,2)OMEGA(I),GIM(4,I),GRE(4,I)
      ENDDO
      CLOSE(5)
      OPEN(36,FILE='DATA5.TXT')
      DO I=1,NR
        READ(36,2)OMEGA(I),GIM(5,I),GRE(5,I)
      ENDDO
      CLOSE(36)
      OPEN(7,FILE='DATA6.TXT')
      DO I=1,NR
         READ(7,2)OMEGA(I),GIM(6,I),GRE(6,I)
      ENDDO
      CLOSE(7)
      OPEN(8,FILE='DATA7.TXT')
      DO I=1,NR
         READ(8,2)OMEGA(I),GIM(7,I),GRE(7,I)
      ENDDO
      CLOSE(8)
      OPEN(9,FILE='DATA8.TXT')
      DO I=1,NR
         READ(9,2)OMEGA(I),GIM(8,I),GRE(8,I)
      ENDDO
      CLOSE(9)
      OPEN(10,FILE='DATA9.TXT')
      DO I=1,NR
         READ(10,2)OMEGA(I),GIM(9,I),GRE(9,I)
      ENDDO
      CLOSE(10)
      OPEN(11,FILE='DATA10.TXT')
      DO I=1,NR
         READ(11,2)OMEGA(I),GIM(10,I),GRE(10,I)
      ENDDO
      CLOSE(11) 
      OPEN(12,FILE='DATA11.TXT')
      DO I=1,NR
         READ(12,2)OMEGA(I),GIM(11,I),GRE(11,I)
      ENDDO
      CLOSE(12)
      OPEN(13,FILE='DATA12.TXT')
      DO I=1,NR
         READ(13,2)OMEGA(I),GIM(12,I),GRE(12,I)
      ENDDO
      CLOSE(13)
      OPEN(14,FILE='DATA13.TXT')
      DO I=1,NR
         READ(14,2)OMEGA(I),GIM(13,I),GRE(13,I)
      ENDDO
      CLOSE(14)
      OPEN(15,FILE='DATA14.TXT')
      DO I=1,NR
         READ(15,2)OMEGA(I),GIM(14,I),GRE(14,I)
      ENDDO
      CLOSE(15)
      OPEN(16,FILE='DATA15.TXT')
      DO I=1,NR
         READ(16,2)OMEGA(I),GIM(15,I),GRE(15,I)
      ENDDO
      CLOSE(16)
      OPEN(17,FILE='DATA16.TXT')
      DO I=1,NR
         READ(17,2)OMEGA(I),GIM(16,I),GRE(16,I)
      ENDDO
      CLOSE(17)
      OPEN(18,FILE='DATA17.TXT')
      DO I=1,nr
         READ(18,2)OMEGA(I),GIM(17,I),GRE(17,I)
      ENDDO
      CLOSE(18)
      OPEN(19,FILE='DATA18.TXT')
      DO I=1,nr
         READ(19,2)OMEGA(I),GIM(18,I),GRE(18,I)
      ENDDO
      CLOSE(19)
      OPEN(20,FILE='DATA19.TXT')
      DO I=1,nr
         READ(20,2)OMEGA(I),GIM(19,I),GRE(19,I)
      ENDDO
      CLOSE(20)
      OPEN(21,FILE='DATA20.TXT')
      DO I=1,nr
         READ(21,2)OMEGA(I),GIM(20,I),GRE(20,I)
      ENDDO
      CLOSE(21)
      OPEN(22,FILE='DATA21.TXT')
      DO I=1,nr
         READ(22,2)OMEGA(I),GIM(21,I),GRE(21,I)
      ENDDO
      CLOSE(22)
      OPEN(23,FILE='DATA22.TXT')
      DO I=1,nr
         READ(23,2)OMEGA(I),GIM(22,I),GRE(22,I)
      ENDDO
      CLOSE(23)
      OPEN(24,FILE='DATA23.TXT')
      DO I=1,nr
         READ(24,2)OMEGA(I),GIM(23,I),GRE(23,I)
      ENDDO
      CLOSE(24)
      OPEN(25,FILE='DATA24.TXT')
      DO I=1,nr
         READ(25,2)OMEGA(I),GIM(24,I),GRE(24,I)
      ENDDO
      CLOSE(25)
      OPEN(26,FILE='DATA25.TXT')
      DO I=1,nr
         READ(26,2)OMEGA(I),GIM(25,I),GRE(25,I)
      ENDDO
      CLOSE(26)
      OPEN(27,FILE='DATA26.TXT')
      DO I=1,nr
         READ(27,2)OMEGA(I),GIM(26,I),GRE(26,I)
      ENDDO
      CLOSE(27)
      OPEN(28,FILE='DATA27.TXT')
      DO I=1,nr
         READ(28,2)OMEGA(I),GIM(27,I),GRE(27,I)
      ENDDO
      CLOSE(28)
      OPEN(29,FILE='DATA28.TXT')
      DO I=1,nr
         READ(29,2)OMEGA(I),GIM(28,I),GRE(28,I)
      ENDDO
      CLOSE(29)
      OPEN(30,FILE='DATA29.TXT')
      DO I=1,nr
         READ(30,2)OMEGA(I),GIM(29,I),GRE(29,I)
      ENDDO
      CLOSE(30)
      OPEN(31,FILE='DATA30.TXT')
      DO I=1,nr
         READ(31,2)OMEGA(I),GIM(30,I),GRE(30,I)
      ENDDO
      CLOSE(31)
      OPEN(32,FILE='DATA31.TXT')
      DO I=1,nr
         READ(32,2)OMEGA(I),GIM(31,I),GRE(31,I)
      ENDDO
      CLOSE(32)
      OPEN(33,FILE='DATA32.TXT')
      DO I=1,nr
         READ(33,2)OMEGA(I),GIM(32,I),GRE(32,I)
      ENDDO
      CLOSE(33)
      OPEN(34,FILE='DATA33.TXT')
      DO I=1,nr
         READ(34,2)OMEGA(I),GIM(33,I),GRE(33,I)
      ENDDO
      CLOSE(34)
      OPEN(35,FILE='DATA34.TXT')
      DO I=1,nr
         READ(35,2)OMEGA(I),GIM(34,I),GRE(34,I)
      ENDDO
      CLOSE(35)    
      OPEN(36,FILE='DATA35.TXT')
      DO I=1,nr
         READ(36,2)OMEGA(I),GIM(35,I),GRE(35,I)
      ENDDO
      CLOSE(36)
      OPEN(37,FILE='DATA36.TXT')
      DO I=1,nr
         READ(37,2)OMEGA(I),GIM(36,I),GRE(36,I)
      ENDDO
      CLOSE(37)    
      DO J=1,36
         DO I=1,NR
            G(J,I)=DCMPLX(GRE(J,I),GIM(J,I))
         ENDDO
      ENDDO
 2    FORMAT(3E16.8)
      END

      SUBROUTINE TIGR(GG,G,I,W)
      INTEGER I1,I2,I
      COMPLEX*16 GG(9,9),TMP(9,9),GG0(9,9),G(36,1375),DC,W(9,9),
     *     AL,BE,KOFF,GG1(3,3),GG2(3,3),GG3(3,3)
      CHARACTER*1 TR
      AL=DCMPLX(1.0D0,0.0D0)
      BE=DCMPLX(0.0D0,0.0D0)
      TR='N'
      KOFF=DCMPLX(95.33D0,0.0D0)
      DO I1=1,9
         DO I2=1,9
            GG0(I1,I2)= DCMPLX(0.0D0,0.0D0)
         ENDDO
      ENDDO
      GG1(1,1) = KOFF*G(12,I)
      GG1(1,2) = KOFF*G(10,I)
      GG1(1,3) = KOFF*G(10,I)
      GG1(2,1) = KOFF*G(10,I)
      GG1(2,2) = KOFF*G(1,I)
      GG1(2,3) = DCMPLX(0.0D0,0.0D0)
      GG1(3,1) = KOFF*G(10,I)
      GG1(3,2) = DCMPLX(0.0D0,0.0D0)
      GG1(3,3) = KOFF*G(1,I)
      GG2(1,1) = KOFF*G(12,I)
      GG2(1,2) = KOFF*G(10,I)
      GG2(1,3) = KOFF*G(10,I)
      GG2(2,1) = KOFF*G(10,I)
      GG2(2,2) = KOFF*G(1,I)
      GG2(2,3) = DCMPLX(0.0D0,0.0D0)
      GG2(3,1) = KOFF*G(10,I)
      GG2(3,2) = DCMPLX(0.0D0,0.0D0)
      GG2(3,3) = KOFF*G(1,I)
      GG3(1,1) = KOFF*G(12,I)
      GG3(1,2) = KOFF*G(10,I)
      GG3(1,3) = KOFF*G(10,I)
      GG3(2,1) = KOFF*G(10,I)
      GG3(2,2) = KOFF*G(1,I)
      GG3(2,3) = DCMPLX(0.0D0,0.0D0)
      GG3(3,1) = KOFF*G(10,I)
      GG3(3,2) = DCMPLX(0.0D0,0.0D0)
      GG3(3,3) = KOFF*G(1,I)
      DO I1=1,3
         DO I2=1,3
            GG0(I1,I2)=GG1(I1,I2)
            GG0(I1+3,I2+3)=GG2(I1,I2)
            GG0(I1+6,I2+6)=GG2(I1,I2)
         ENDDO
      ENDDO      
      CALL F06ZAF(TR,TR,9,9,9,AL,GG0,9,W,9,BE,TMP,9)
      DO I1=1,9
         DO I2=1,9
            TMP(I1,I2)=DC(I1,I2)-TMP(I1,I2)
         ENDDO
      ENDDO
      CALL INVERSE1(TMP)
      CALL F06ZAF(TR,TR,9,9,9,AL,GG0,9,TMP,9,BE,GG,9)
 2    FORMAT(3E16.8)
      RETURN
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
 20   H2 = X(2) - X(1)
      DO 80 I = 3, N
         H3 = X(I) - X(I-1)
         IF (H2*H3) 40, 60, 80
 40      IFAIL = P01ABF(IFAIL,2,SRNAME,0,P01REC)
         RETURN
 60      IFAIL = P01ABF(IFAIL,3,SRNAME,0,P01REC)
         RETURN
 80   CONTINUE
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
     *     )
      S = -(H2**3)*(H2*(3.0D0*H2+5.0D0*H4)+10.0D0*H3*H1)/60.0D0
      R4 = 0.0D0
C     
C     INTEGRATE OVER CENTRAL PORTION OF RANGE
C     
      NN = N - 1
      DO 120 I = 3, NN
         ANS = ANS + H3*((Y(I)+Y(I-1))/2.0D0-H3*H3*(D2+R2+(H2-H4)*R3)
     *        /12.0D0)
         C = H3**3*(2.0D0*H3*H3+5.0D0*(H3*(H4+H2)+2.0D0*H4*H2))/120.0D0
         ER = ER + (C+S)*R4
         IF (I.NE.3) S = C
         IF (I.EQ.3) S = S + 2.0D0*C
         IF (I-N+1) 100, 140, 100
 100     H1 = H2
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
 120  CONTINUE
C     
C     INTEGRATE OVER FINAL INTERVAL
C     
 140  CONTINUE
      ANS = ANS + H4*(Y(N)-H4*(R1/2.0D0+H4*(R2/6.0D0+(2.0D0*H3+H4)
     *     *R3/12.0D0)))
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
C     Abnormal exit from calling routine
         IF (IFAIL.EQ.-1 .OR. IFAIL.EQ.0 .OR. IFAIL.EQ.-13 .OR.
     *        (IFAIL.GT.0 .AND. MOD(IFAIL/10,10).NE.0)) THEN
C     Noisy exit
            CALL X04AAF(0,NERR)
            DO 20 I = 1, NREC
               CALL X04BAF(NERR,REC(I))
 20         CONTINUE
            IF (IFAIL.NE.-13) THEN
               WRITE (MESS,FMT=99999) SRNAME, IERROR
               CALL X04BAF(NERR,MESS)
               IF (ABS(MOD(IFAIL,10)).NE.1) THEN
C     Hard failure
                  CALL X04BAF(NERR,
     *                 ' ** NAG hard failure - execution terminated'
     *                 )
                  CALL P01ABZ
               ELSE
C     Soft failure
                  CALL X04BAF(NERR,
     *                 ' ** NAG soft failure - control returned')
               END IF
            END IF
         END IF
      END IF
      P01ABF = IERROR
      RETURN
C     
99999 FORMAT (' ** ABNORMAL EXIT from NAG Library routine ',A,': IFAIL',
     *     ' =',I6)
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
C     Remove trailing blanks
         DO 20 I = LEN(REC), 2, -1
            IF (REC(I:I).NE.' ') GO TO 40
 20      CONTINUE
C     Write record to external file
   40    WRITE (NOUT,FMT=99999) REC(1:I)
      END IF
      RETURN
C     
99999 FORMAT (A)
      END


      DOUBLE PRECISION FUNCTION X01AAF(X)
C     MARK 8 RE-ISSUE. NAG COPYRIGHT 1980.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     RETURNS THE VALUE OF THE MATHEMATICAL CONSTANT PI.
C
C     X IS A DUMMY ARGUMENT
C
C     IT MAY BE NECESSARY TO ROUND THE REAL CONSTANT IN THE
C     ASSIGNMENT STATEMENT TO A SMALLER NUMBER OF SIGNIFICANT
C     DIGITS IN ORDER TO AVOID COMPILATION PROBLEMS.  IF SO, THEN
C     THE NUMBER OF DIGITS RETAINED SHOULD NOT BE LESS THAN
C     .     2 + INT(FLOAT(IT)*ALOG10(IB))
C     WHERE  IB  IS THE BASE FOR THE REPRESENTATION OF FLOATING-
C     .             POINT NUMBERS
C     . AND  IT  IS THE NUMBER OF IB-ARY DIGITS IN THE MANTISSA OF
C     .             A FLOATING-POINT NUMBER.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X
C     .. Executable Statements ..
      X01AAF = 3.14159265358979323846264338328D0
      RETURN
      END

      SUBROUTINE INVERSE1(A)
      INTEGER NP
      PARAMETER (NP=9)
      INTEGER IPVT(NP),JOB
      DOUBLE PRECISION RCOND 
      COMPLEX*16 A(NP,NP),Z(NP),DET(2),WORK(NP)
      JOB=01
      CALL CGECO(A,NP,NP,IPVT,RCOND,Z)
      CALL CGEDI(A,NP,NP,IPVT,DET,WORK,JOB)
      RETURN
      END

    
      SUBROUTINE F06ZAF(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,
     *                  LDC)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C  Purpose
C  =======
C
C  ZGEMM  performs one of the matrix-matrix operations
C
C     C := alpha*op( A )*op( B ) + beta*C,
C
C  where  op( X ) is one of
C
C     op( X ) = X   or   op( X ) = X'   or   op( X ) = conjg( X' ),
C
C  alpha and beta are scalars, and A, B and C are matrices, with op( A )
C  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
C
C  Parameters
C  ==========
C
C  TRANSA - CHARACTER*1.
C           On entry, TRANSA specifies the form of op( A ) to be used in
C           the matrix multiplication as follows:
C
C              TRANSA = 'N' or 'n',  op( A ) = A.
C
C              TRANSA = 'T' or 't',  op( A ) = A'.
C
C              TRANSA = 'C' or 'c',  op( A ) = conjg( A' ).
C
C           Unchanged on exit.
C
C  TRANSB - CHARACTER*1.
C           On entry, TRANSB specifies the form of op( B ) to be used in
C           the matrix multiplication as follows:
C
C              TRANSB = 'N' or 'n',  op( B ) = B.
C
C              TRANSB = 'T' or 't',  op( B ) = B'.
C
C              TRANSB = 'C' or 'c',  op( B ) = conjg( B' ).
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C           On entry,  M  specifies  the number  of rows  of the  matrix
C           op( A )  and of the  matrix  C.  M  must  be at least  zero.
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry,  N  specifies the number  of columns of the matrix
C           op( B ) and the number of columns of the matrix C. N must be
C           at least zero.
C           Unchanged on exit.
C
C  K      - INTEGER.
C           On entry,  K  specifies  the number of columns of the matrix
C           op( A ) and the number of rows of the matrix op( B ). K must
C           be at least  zero.
C           Unchanged on exit.
C
C  ALPHA  - COMPLEX         .
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  A      - COMPLEX          array of DIMENSION ( LDA, ka ), where ka is
C           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
C           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
C           part of the array  A  must contain the matrix  A,  otherwise
C           the leading  k by m  part of the array  A  must contain  the
C           matrix A.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
C           LDA must be at least  max( 1, m ), otherwise  LDA must be at
C           least  max( 1, k ).
C           Unchanged on exit.
C
C  B      - COMPLEX          array of DIMENSION ( LDB, kb ), where kb is
C           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
C           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
C           part of the array  B  must contain the matrix  B,  otherwise
C           the leading  n by k  part of the array  B  must contain  the
C           matrix B.
C           Unchanged on exit.
C
C  LDB    - INTEGER.
C           On entry, LDB specifies the first dimension of B as declared
C           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
C           LDB must be at least  max( 1, k ), otherwise  LDB must be at
C           least  max( 1, n ).
C           Unchanged on exit.
C
C  BETA   - COMPLEX         .
C           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
C           supplied as zero then C need not be set on input.
C           Unchanged on exit.
C
C  C      - COMPLEX          array of DIMENSION ( LDC, n ).
C           Before entry, the leading  m by n  part of the array  C must
C           contain the matrix  C,  except when  beta  is zero, in which
C           case C need not be set on entry.
C           On exit, the array  C  is overwritten by the  m by n  matrix
C           ( alpha*op( A )*op( B ) + beta*C ).
C
C  LDC    - INTEGER.
C           On entry, LDC specifies the first dimension of C as declared
C           in  the  calling  (sub)  program.   LDC  must  be  at  least
C           max( 1, m ).
C           Unchanged on exit.
C
C
C  Level 3 Blas routine.
C
C  -- Written on 8-February-1989.
C     Jack Dongarra, Argonne National Laboratory.
C     Iain Duff, AERE Harwell.
C     Jeremy Du Croz, Numerical Algorithms Group Ltd.
C     Sven Hammarling, Numerical Algorithms Group Ltd.
C
C
C     .. Entry Points ..
      ENTRY             ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,
     *                  BETA,C,LDC)
C     .. Parameters ..
      COMPLEX*16        ONE
      PARAMETER         (ONE=(1.0D+0,0.0D+0))
      COMPLEX*16        ZERO
      PARAMETER         (ZERO=(0.0D+0,0.0D+0))
C     .. Scalar Arguments ..
      COMPLEX*16        ALPHA, BETA
      INTEGER           K, LDA, LDB, LDC, M, N
      CHARACTER*1       TRANSA, TRANSB
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), B(LDB,*), C(LDC,*)
C     .. Local Scalars ..
      COMPLEX*16        TEMP
      INTEGER           I, INFO, J, L, NCOLA, NROWA, NROWB
      LOGICAL           CONJA, CONJB, NOTA, NOTB
C     .. External Subroutines ..
      EXTERNAL          F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC         DCONJG, MAX
C     .. Executable Statements ..
C
C     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
C     conjugated or transposed, set  CONJA and CONJB  as true if  A  and
C     B  respectively are to be  transposed but  not conjugated  and set
C     NROWA, NCOLA and  NROWB  as the number of rows and  columns  of  A
C     and the number of rows of  B  respectively.
C
      NOTA = (TRANSA.EQ.'N' .OR. TRANSA.EQ.'n')
      NOTB = (TRANSB.EQ.'N' .OR. TRANSB.EQ.'n')
      CONJA = (TRANSA.EQ.'C' .OR. TRANSA.EQ.'c')
      CONJB = (TRANSB.EQ.'C' .OR. TRANSB.EQ.'c')
      IF (NOTA) THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF (NOTB) THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
C
C     Test the input parameters.
C
      INFO = 0
      IF (( .NOT. NOTA) .AND. ( .NOT. CONJA)
     *    .AND. ( .NOT. (TRANSA.EQ.'T' .OR. TRANSA.EQ.'t'))) THEN
         INFO = 1
      ELSE IF (( .NOT. NOTB) .AND. ( .NOT. CONJB)
     *         .AND. ( .NOT. (TRANSB.EQ.'T' .OR. TRANSB.EQ.'t'))) THEN
         INFO = 2
      ELSE IF (M.LT.0) THEN
         INFO = 3
      ELSE IF (N.LT.0) THEN
         INFO = 4
      ELSE IF (K.LT.0) THEN
         INFO = 5
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
         INFO = 8
      ELSE IF (LDB.LT.MAX(1,NROWB)) THEN
         INFO = 10
      ELSE IF (LDC.LT.MAX(1,M)) THEN
         INFO = 13
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F06ZAF/ZGEMM ',INFO)
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR. (((ALPHA.EQ.ZERO) .OR. (K.EQ.0))
     *     .AND. (BETA.EQ.ONE))) RETURN
C
C     And when  alpha.eq.zero.
C
      IF (ALPHA.EQ.ZERO) THEN
         IF (BETA.EQ.ZERO) THEN
            DO 40 J = 1, N
               DO 20 I = 1, M
                  C(I,J) = ZERO
   20          CONTINUE
   40       CONTINUE
         ELSE
            DO 80 J = 1, N
               DO 60 I = 1, M
                  C(I,J) = BETA*C(I,J)
   60          CONTINUE
   80       CONTINUE
         END IF
         RETURN
      END IF
C
C     Start the operations.
C
      IF (NOTB) THEN
         IF (NOTA) THEN
C
C           Form  C := alpha*A*B + beta*C.
C
            DO 180 J = 1, N
               IF (BETA.EQ.ZERO) THEN
                  DO 100 I = 1, M
                     C(I,J) = ZERO
  100             CONTINUE
               ELSE IF (BETA.NE.ONE) THEN
                  DO 120 I = 1, M
                     C(I,J) = BETA*C(I,J)
  120             CONTINUE
               END IF
               DO 160 L = 1, K
                  IF (B(L,J).NE.ZERO) THEN
                     TEMP = ALPHA*B(L,J)
                     DO 140 I = 1, M
                        C(I,J) = C(I,J) + TEMP*A(I,L)
  140                CONTINUE
                  END IF
  160          CONTINUE
  180       CONTINUE
         ELSE IF (CONJA) THEN
C
C           Form  C := alpha*conjg( A' )*B + beta*C.
C
            DO 240 J = 1, N
               DO 220 I = 1, M
                  TEMP = ZERO
                  DO 200 L = 1, K
                     TEMP = TEMP + DCONJG(A(L,I))*B(L,J)
  200             CONTINUE
                  IF (BETA.EQ.ZERO) THEN
                     C(I,J) = ALPHA*TEMP
                  ELSE
                     C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                  END IF
  220          CONTINUE
  240       CONTINUE
         ELSE
C
C           Form  C := alpha*A'*B + beta*C
C
            DO 300 J = 1, N
               DO 280 I = 1, M
                  TEMP = ZERO
                  DO 260 L = 1, K
                     TEMP = TEMP + A(L,I)*B(L,J)
  260             CONTINUE
                  IF (BETA.EQ.ZERO) THEN
                     C(I,J) = ALPHA*TEMP
                  ELSE
                     C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                  END IF
  280          CONTINUE
  300       CONTINUE
         END IF
      ELSE IF (NOTA) THEN
         IF (CONJB) THEN
C
C           Form  C := alpha*A*conjg( B' ) + beta*C.
C
            DO 400 J = 1, N
               IF (BETA.EQ.ZERO) THEN
                  DO 320 I = 1, M
                     C(I,J) = ZERO
  320             CONTINUE
               ELSE IF (BETA.NE.ONE) THEN
                  DO 340 I = 1, M
                     C(I,J) = BETA*C(I,J)
  340             CONTINUE
               END IF
               DO 380 L = 1, K
                  IF (B(J,L).NE.ZERO) THEN
                     TEMP = ALPHA*DCONJG(B(J,L))
                     DO 360 I = 1, M
                        C(I,J) = C(I,J) + TEMP*A(I,L)
  360                CONTINUE
                  END IF
  380          CONTINUE
  400       CONTINUE
         ELSE
C
C           Form  C := alpha*A*B'          + beta*C
C
            DO 500 J = 1, N
               IF (BETA.EQ.ZERO) THEN
                  DO 420 I = 1, M
                     C(I,J) = ZERO
  420             CONTINUE
               ELSE IF (BETA.NE.ONE) THEN
                  DO 440 I = 1, M
                     C(I,J) = BETA*C(I,J)
  440             CONTINUE
               END IF
               DO 480 L = 1, K
                  IF (B(J,L).NE.ZERO) THEN
                     TEMP = ALPHA*B(J,L)
                     DO 460 I = 1, M
                        C(I,J) = C(I,J) + TEMP*A(I,L)
  460                CONTINUE
                  END IF
  480          CONTINUE
  500       CONTINUE
         END IF
      ELSE IF (CONJA) THEN
         IF (CONJB) THEN
C
C           Form  C := alpha*conjg( A' )*conjg( B' ) + beta*C.
C
            DO 560 J = 1, N
               DO 540 I = 1, M
                  TEMP = ZERO
                  DO 520 L = 1, K
                     TEMP = TEMP + DCONJG(A(L,I))*DCONJG(B(J,L))
  520             CONTINUE
                  IF (BETA.EQ.ZERO) THEN
                     C(I,J) = ALPHA*TEMP
                  ELSE
                     C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                  END IF
  540          CONTINUE
  560       CONTINUE
         ELSE
C
C           Form  C := alpha*conjg( A' )*B' + beta*C
C
            DO 620 J = 1, N
               DO 600 I = 1, M
                  TEMP = ZERO
                  DO 580 L = 1, K
                     TEMP = TEMP + DCONJG(A(L,I))*B(J,L)
  580             CONTINUE
                  IF (BETA.EQ.ZERO) THEN
                     C(I,J) = ALPHA*TEMP
                  ELSE
                     C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                  END IF
  600          CONTINUE
  620       CONTINUE
         END IF
      ELSE
         IF (CONJB) THEN
C
C           Form  C := alpha*A'*conjg( B' ) + beta*C
C
            DO 680 J = 1, N
               DO 660 I = 1, M
                  TEMP = ZERO
                  DO 640 L = 1, K
                     TEMP = TEMP + A(L,I)*DCONJG(B(J,L))
  640             CONTINUE
                  IF (BETA.EQ.ZERO) THEN
                     C(I,J) = ALPHA*TEMP
                  ELSE
                     C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                  END IF
  660          CONTINUE
  680       CONTINUE
         ELSE
C
C           Form  C := alpha*A'*B' + beta*C
C
            DO 740 J = 1, N
               DO 720 I = 1, M
                  TEMP = ZERO
                  DO 700 L = 1, K
                     TEMP = TEMP + A(L,I)*B(J,L)
  700             CONTINUE
                  IF (BETA.EQ.ZERO) THEN
                     C(I,J) = ALPHA*TEMP
                  ELSE
                     C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                  END IF
  720          CONTINUE
  740       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06ZAF (ZGEMM ).
C
      END


      SUBROUTINE CGECO (A, LDA, N, IPVT, RCOND, Z)
C***  BEGIN PROLOGUE  CGECO
C***  PURPOSE  Factor a matrix using Gaussian elimination and estimate
C     the condition number of the matrix.
C***  LIBRARY   SLATEC (LINPACK)
C***  CATEGORY  D2C1
C***  TYPE      COMPLEX (SGECO-S, DGECO-D, CGECO-C)
C***  KEYWORDS  CONDITION NUMBER, GENERAL MATRIX, LINEAR ALGEBRA, LINPACK,
C     MATRIX FACTORIZATION
C***  AUTHOR  Moler, C. B., (U. of New Mexico)
C***  DESCRIPTION
C     
C     CGECO factors a complex matrix by Gaussian elimination
C     and estimates the condition of the matrix.
C     
C     If  RCOND  is not needed, CGEFA is slightly faster.
C     To solve  A*X = B , follow CGECO By CGESL.
C     To Compute  INVERSE(A)*C , follow CGECO by CGESL.
C     To compute  DETERMINANT(A) , follow CGECO by CGEDI.
C     To compute  INVERSE(A) , follow CGECO by CGEDI.
C     
C     On Entry
C     
C     A       COMPLEX(LDA, N)
C     the matrix to be factored.
C     
C     LDA     INTEGER
C     the leading dimension of the array  A .
C     
C     N       INTEGER
C     the order of the matrix  A .
C     
C     On Return
C     
C     A       an upper triangular matrix and the multipliers
C             which were used to obtain it.
C             The factorization can be written  A = L*U  where
C             L  is a product of permutation and unit lower
C             triangular matrices and  U  is upper triangular.
C     
C     IPVT    INTEGER(N)
C                an integer vector of pivot indices.
C
C        RCOND   REAL
C                an estimate of the reciprocal condition of  A .
C                For the system  A*X = B , relative perturbations
C                in  A  and  B  of size  EPSILON  may cause
C                relative perturbations in  X  of size  EPSILON/RCOND .
C                If  RCOND  is so small that the logical expression
C                           1.0 + RCOND .EQ. 1.0
C                is true, then  A  may be singular to working
C                precision.  In particular,  RCOND  is zero  if
C                exact singularity is detected or the estimate
C                underflows.
C
C        Z       COMPLEX(N)
C                a work vector whose contents are usually unimportant.
C                If  A  is close to a singular matrix, then  Z  is
C                an approximate null vector in the sense that
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  CAXPY, CDOTC, CGEFA, CSSCAL, SCASUM
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  CGECO
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LDA,N,IPVT(*)
      COMPLEX*16 A(LDA,*),Z(*)
      DOUBLE PRECISION RCOND
      COMPLEX*16 CDOTC,EK,T,WK,WKM
      DOUBLE PRECISION ANORM,S,SCASUM,SM,YNORM
      INTEGER INFO,J,K,KB,KP1,L
      COMPLEX*16 ZDUM,ZDUM1,ZDUM2,CSIGN1
      DOUBLE PRECISION CABS1
      CABS1(ZDUM) = DABS(DREAL(ZDUM)) + DABS(DIMAG(ZDUM))
      CSIGN1(ZDUM1,ZDUM2) = CABS1(ZDUM1)*(ZDUM2/CABS1(ZDUM2))
C
C     COMPUTE 1-NORM OF A
C
C***FIRST EXECUTABLE STATEMENT  CGECO
      ANORM = 0.0D0
      DO 10 J = 1, N
         ANORM = MAX(ANORM,SCASUM(N,A(1,J),1))
   10 CONTINUE
C
C     FACTOR
C
      CALL CGEFA(A,LDA,N,IPVT,INFO)
C
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  CTRANS(A)*Y = E .
C     CTRANS(A)  IS THE CONJUGATE TRANSPOSE OF A .
C     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
C     GROWTH IN THE ELEMENTS OF W  WHERE  CTRANS(U)*W = E .
C     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
C
C     SOLVE CTRANS(U)*W = E
C
      EK = DCMPLX(1.0D0,0.0D0)
      DO 20 J = 1, N
         Z(J) = DCMPLX(0.0D0,0.0D0)
   20 CONTINUE
      DO 100 K = 1, N
         IF (CABS1(Z(K)) .NE. 0.0D0) EK = CSIGN1(EK,-Z(K))
         IF (CABS1(EK-Z(K)) .LE. CABS1(A(K,K))) GO TO 30
            S = CABS1(A(K,K))/CABS1(EK-Z(K))
            CALL CSSCAL(N,S,Z,1)
            EK = CMPLX(S,0.0D0)*EK
   30    CONTINUE
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = CABS1(WK)
         SM = CABS1(WKM)
         IF (CABS1(A(K,K)) .EQ. 0.0D0) GO TO 40
            WK = WK/CONJG(A(K,K))
            WKM = WKM/CONJG(A(K,K))
         GO TO 50
   40    CONTINUE
            WK = DCMPLX(1.0D0,0.0D0)
            WKM = DCMPLX(1.0D0,0.0D0)
   50    CONTINUE
         KP1 = K + 1
         IF (KP1 .GT. N) GO TO 90
            DO 60 J = KP1, N
               SM = SM + CABS1(Z(J)+WKM*CONJG(A(K,J)))
               Z(J) = Z(J) + WK*CONJG(A(K,J))
               S = S + CABS1(Z(J))
   60       CONTINUE
            IF (S .GE. SM) GO TO 80
               T = WKM - WK
               WK = WKM
               DO 70 J = KP1, N
                  Z(J) = Z(J) + T*CONJG(A(K,J))
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
         Z(K) = WK
  100 CONTINUE
      S = 1.0D0/SCASUM(N,Z,1)
      CALL CSSCAL(N,S,Z,1)
C
C     SOLVE CTRANS(L)*Y = W
C
      DO 120 KB = 1, N
         K = N + 1 - KB
         IF (K .LT. N) Z(K) = Z(K) + CDOTC(N-K,A(K+1,K),1,Z(K+1),1)
         IF (CABS1(Z(K)) .LE. 1.0E0) GO TO 110
            S = 1.0D0/CABS1(Z(K))
            CALL CSSCAL(N,S,Z,1)
  110    CONTINUE
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
  120 CONTINUE
      S = 1.0D0/SCASUM(N,Z,1)
      CALL CSSCAL(N,S,Z,1)
C
      YNORM = 1.0D0
C
C     SOLVE L*V = Y
C
      DO 140 K = 1, N
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
         IF (K .LT. N) CALL CAXPY(N-K,T,A(K+1,K),1,Z(K+1),1)
         IF (CABS1(Z(K)) .LE. 1.0D0) GO TO 130
            S = 1.0D0/CABS1(Z(K))
            CALL CSSCAL(N,S,Z,1)
            YNORM = S*YNORM
  130    CONTINUE
  140 CONTINUE
      S = 1.0D0/SCASUM(N,Z,1)
      CALL CSSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
C     SOLVE  U*Z = V
C
      DO 160 KB = 1, N
         K = N + 1 - KB
         IF (CABS1(Z(K)) .LE. CABS1(A(K,K))) GO TO 150
            S = CABS1(A(K,K))/CABS1(Z(K))
            CALL CSSCAL(N,S,Z,1)
            YNORM = S*YNORM
  150    CONTINUE
         IF (CABS1(A(K,K)) .NE. 0.0D0) Z(K) = Z(K)/A(K,K)
         IF (CABS1(A(K,K)) .EQ. 0.0D0) Z(K) = DCMPLX(1.0D0,0.0D0)
         T = -Z(K)
         CALL CAXPY(K-1,T,A(1,K),1,Z(1),1)
  160 CONTINUE
C     MAKE ZNORM = 1.0
      S = 1.0D0/SCASUM(N,Z,1)
      CALL CSSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
      IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0
      RETURN
      END
*DECK CAXPY
      SUBROUTINE CAXPY (N, CA, CX, INCX, CY, INCY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C***BEGIN PROLOGUE  CAXPY
C***PURPOSE  Compute a constant times a vector plus a vector.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A7
C***TYPE      COMPLEX (SAXPY-S, DAXPY-D, CAXPY-C)
C***KEYWORDS  BLAS, LINEAR ALGEBRA, TRIAD, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       CA  complex scalar multiplier
C       CX  complex vector with N elements
C     INCX  storage spacing between elements of CX
C       CY  complex vector with N elements
C     INCY  storage spacing between elements of CY
C
C     --Output--
C       CY  complex result (unchanged if N .LE. 0)
C
C     Overwrite complex CY with complex  CA*CX + CY.
C     For I = 0 to N-1, replace  CY(LY+I*INCY) with CA*CX(LX+I*INCX) +
C       CY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
C     defined in a similar way using INCY.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C   920801  Removed variable CANORM.  (RWC, WRB)
C***END PROLOGUE  CAXPY
      COMPLEX*16 CX(*), CY(*), CA
C***FIRST EXECUTABLE STATEMENT  CAXPY
      IF (N.LE.0 .OR. CA.EQ.DCMPLX(0.0D0,0.0D0)) RETURN
      IF (INCX.EQ.INCY .AND. INCX.GT.0) GO TO 20
C
C     Code for unequal or nonpositive increments.
C
      KX = 1
      KY = 1
      IF (INCX .LT. 0) KX = 1+(1-N)*INCX
      IF (INCY .LT. 0) KY = 1+(1-N)*INCY
      DO 10 I = 1,N
        CY(KY) = CY(KY) + CA*CX(KX)
        KX = KX + INCX
        KY = KY + INCY
   10 CONTINUE
      RETURN
C
C     Code for equal, positive, non-unit increments.
C
   20 NS = N*INCX
      DO 30 I = 1,NS,INCX
        CY(I) = CA*CX(I) + CY(I)
   30 CONTINUE
      RETURN
      END
*DECK CDOTC
      COMPLEX*16 FUNCTION CDOTC (N, CX, INCX, CY, INCY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C***BEGIN PROLOGUE  CDOTC
C***PURPOSE  Dot product of two complex vectors using the complex
C            conjugate of the first vector.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A4
C***TYPE      COMPLEX (CDOTC-C)
C***KEYWORDS  BLAS, INNER PRODUCT, LINEAR ALGEBRA, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       CX  complex vector with N elements
C     INCX  storage spacing between elements of CX
C       CY  complex vector with N elements
C     INCY  storage spacing between elements of CY
C
C     --Output--
C    CDOTC  complex result (zero if N .LE. 0)
C
C     Returns the dot product of complex CX and CY, using CONJUGATE(CX)
C     CDOTC = SUM for I = 0 to N-1 of CONJ(CX(LX+I*INCX))*CY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
C     defined in a similar way using INCY.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  CDOTC
      COMPLEX*16 CX(*),CY(*)
C***FIRST EXECUTABLE STATEMENT  CDOTC
      CDOTC = DCMPLX(0.0D0,0.0D0)
      IF (N .LE. 0) RETURN
      IF (INCX.EQ.INCY .AND. INCX.GT.0) GO TO 20
C
C     Code for unequal or nonpositive increments.
C
      KX = 1
      KY = 1
      IF (INCX .LT. 0) KX = 1+(1-N)*INCX
      IF (INCY .LT. 0) KY = 1+(1-N)*INCY
      DO 10 I = 1,N
        CDOTC = CDOTC + CONJG(CX(KX))*CY(KY)
        KX = KX + INCX
        KY = KY + INCY
   10 CONTINUE
      RETURN
C
C     Code for equal, positive increments.
C
   20 NS = N*INCX
      DO 30 I = 1,NS,INCX
      CDOTC = CDOTC + CONJG(CX(I))*CY(I)
   30 CONTINUE
      RETURN
      END
*DECK CGEFA
      SUBROUTINE CGEFA (A, LDA, N, IPVT, INFO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C***BEGIN PROLOGUE  CGEFA
C***PURPOSE  Factor a matrix using Gaussian elimination.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D2C1
C***TYPE      COMPLEX (SGEFA-S, DGEFA-D, CGEFA-C)
C***KEYWORDS  GENERAL MATRIX, LINEAR ALGEBRA, LINPACK,
C             MATRIX FACTORIZATION
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     CGEFA factors a complex matrix by Gaussian elimination.
C
C     CGEFA is usually called by CGECO, but it can be called
C     directly with a saving in time if  RCOND  is not needed.
C     (Time for CGECO) = (1 + 9/N)*(Time for CGEFA) .
C
C     On Entry
C
C        A       COMPLEX(LDA, N)
C                the matrix to be factored.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C     On Return
C
C        A       an upper triangular matrix and the multipliers
C                which were used to obtain it.
C                The factorization can be written  A = L*U  where
C                L  is a product of permutation and unit lower
C                triangular matrices and  U  is upper triangular.
C
C        IPVT    INTEGER(N)
C                an integer vector of pivot indices.
C
C        INFO    INTEGER
C                = 0  normal value.
C                = K  if  U(K,K) .EQ. 0.0 .  This is not an error
C                     condition for this subroutine, but it does
C                     indicate that CGESL or CGEDI will divide by zero
C                     if called.  Use  RCOND  in CGECO for a reliable
C                     indication of singularity.
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  CAXPY, CSCAL, ICAMAX
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  CGEFA
      INTEGER LDA,N,IPVT(*),INFO
      COMPLEX*16 A(LDA,*)
C
      COMPLEX*16 T
      INTEGER ICAMAX,J,K,KP1,L,NM1
      COMPLEX*16 ZDUM
      DOUBLE PRECISION CABS1
      CABS1(ZDUM) = DABS(DREAL(ZDUM)) + DABS(DIMAG(ZDUM))
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
C***FIRST EXECUTABLE STATEMENT  CGEFA
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
C
C        FIND L = PIVOT INDEX
C
         L = ICAMAX(N-K+1,A(K,K),1) + K - 1
         IPVT(K) = L
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         IF (CABS1(A(L,K)) .EQ. 0.0D0) GO TO 40
C
C           INTERCHANGE IF NECESSARY
C
            IF (L .EQ. K) GO TO 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE
C
C           COMPUTE MULTIPLIERS
C
            T = -(1.0E0,0.0E0)/A(K,K)
            CALL CSCAL(N-K,T,A(K+1,K),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .EQ. K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE
               CALL CAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (CABS1(A(N,N)) .EQ. 0.0D0) INFO = N
      RETURN
      END
*DECK CSCAL
      SUBROUTINE CSCAL (N, CA, CX, INCX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C***BEGIN PROLOGUE  CSCAL
C***PURPOSE  Multiply a vector by a constant.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A6
C***TYPE      COMPLEX (SSCAL-S, DSCAL-D, CSCAL-C)
C***KEYWORDS  BLAS, LINEAR ALGEBRA, SCALE, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       CA  complex scale factor
C       CX  complex vector with N elements
C     INCX  storage spacing between elements of CX
C
C     --Output--
C       CX  complex result (unchanged if N .LE. 0)
C
C     Replace complex CX by complex CA*CX.
C     For I = 0 to N-1, replace CX(IX+I*INCX) with CA*CX(IX+I*INCX),
C     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900821  Modified to correct problem with a negative increment.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  CSCAL
      COMPLEX*16 CA, CX(*)
      INTEGER I, INCX, IX, N
C***FIRST EXECUTABLE STATEMENT  CSCAL
      IF (N .LE. 0) RETURN
C
      IF (INCX .EQ. 1) GOTO 20
C
C     Code for increment not equal to 1.
C
      IX = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      DO 10 I = 1,N
        CX(IX) = CA*CX(IX)
        IX = IX + INCX
   10 CONTINUE
      RETURN
C
C     Code for increment equal to 1.
C
   20 DO 30 I = 1,N
        CX(I) = CA*CX(I)
   30 CONTINUE
      RETURN
      END
*DECK CSSCAL
      SUBROUTINE CSSCAL (N, SA, CX, INCX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C***BEGIN PROLOGUE  CSSCAL
C***PURPOSE  Scale a complex vector.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A6
C***TYPE      COMPLEX (CSSCAL-C)
C***KEYWORDS  BLAS, LINEAR ALGEBRA, SCALE, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       SA  single precision scale factor
C       CX  complex vector with N elements
C     INCX  storage spacing between elements of CX
C
C     --Output--
C       CX  scaled result (unchanged if N .LE. 0)
C
C     Replace complex CX by (single precision SA) * (complex CX)
C     For I = 0 to N-1, replace CX(IX+I*INCX) with  SA * CX(IX+I*INCX),
C     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900821  Modified to correct problem with a negative increment.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  CSSCAL
      COMPLEX*16 CX(*)
      DOUBLE PRECISION SA
      INTEGER I, INCX, IX, N
C***FIRST EXECUTABLE STATEMENT  CSSCAL
      IF (N .LE. 0) RETURN
C
      IF (INCX .EQ. 1) GOTO 20
C
C     Code for increment not equal to 1.
C
      IX = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      DO 10 I = 1,N
        CX(IX) = SA*CX(IX)
        IX = IX + INCX
   10 CONTINUE
      RETURN
C
C     Code for increment equal to 1.
C
   20 DO 30 I = 1,N
        CX(I) = SA*CX(I)
   30 CONTINUE
      RETURN
      END
*DECK ICAMAX
      INTEGER FUNCTION ICAMAX (N, CX, INCX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C***BEGIN PROLOGUE  ICAMAX
C***PURPOSE  Find the smallest index of the component of a complex
C            vector having the maximum sum of magnitudes of real
C            and imaginary parts.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A2
C***TYPE      COMPLEX (ISAMAX-S, IDAMAX-D, ICAMAX-C)
C***KEYWORDS  BLAS, LINEAR ALGEBRA, MAXIMUM COMPONENT, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       CX  complex vector with N elements
C     INCX  storage spacing between elements of CX
C
C     --Output--
C   ICAMAX  smallest index (zero if N .LE. 0)
C
C     Returns the smallest index of the component of CX having the
C     largest sum of magnitudes of real and imaginary parts.
C     ICAMAX = first I, I = 1 to N, to maximize
C     ABS(REAL(CX(IX+(I-1)*INCX))) + ABS(IMAG(CX(IX+(I-1)*INCX))),
C     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900821  Modified to correct problem with a negative increment.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  ICAMAX
      COMPLEX*16 CX(*)
      DOUBLE PRECISION SMAX, XMAG
      INTEGER I, INCX, IX, N
      COMPLEX*16 ZDUM
      DOUBLE PRECISION CABS1
      CABS1(ZDUM) = DABS(DREAL(ZDUM)) + DABS(DIMAG(ZDUM))
C***FIRST EXECUTABLE STATEMENT  ICAMAX
      ICAMAX = 0
      IF (N .LE. 0) RETURN
      ICAMAX = 1
      IF (N .EQ. 1) RETURN
C
      IF (INCX .EQ. 1) GOTO 20
C
C     Code for increment not equal to 1.
C
      IX = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      SMAX = CABS1(CX(IX))
      IX = IX + INCX
      DO 10 I = 2,N
        XMAG = CABS1(CX(IX))
        IF (XMAG .GT. SMAX) THEN
          ICAMAX = I
          SMAX = XMAG
        ENDIF
        IX = IX + INCX
   10 CONTINUE
      RETURN
C
C     Code for increment equal to 1.
C
   20 SMAX = CABS1(CX(1))
      DO 30 I = 2,N
        XMAG = CABS1(CX(I))
        IF (XMAG .GT. SMAX) THEN
          ICAMAX = I
          SMAX = XMAG
        ENDIF
   30 CONTINUE
      RETURN
      END
*DECK SCASUM
      DOUBLE PRECISION FUNCTION SCASUM (N, CX, INCX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C***BEGIN PROLOGUE  SCASUM
C***PURPOSE  Compute the sum of the magnitudes of the real and
C            imaginary elements of a complex vector.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A3A
C***TYPE      COMPLEX (SASUM-S, DASUM-D, SCASUM-C)
C***KEYWORDS  BLAS, LINEAR ALGEBRA, SUM OF MAGNITUDES OF A VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       CX  complex vector with N elements
C     INCX  storage spacing between elements of CX
C
C     --Output--
C   SCASUM  single precision result (zero if N .LE. 0)
C
C     Returns sums of magnitudes of real and imaginary parts of
C     components of CX.  Note that this is not the L1 norm of CX.
C     CASUM = sum from 0 to N-1 of ABS(REAL(CX(IX+I*INCX))) +
C             ABS(IMAG(CX(IX+I*INCX))),
C     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900821  Modified to correct problem with a negative increment.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  SCASUM
      COMPLEX*16 CX(*)
      INTEGER I, INCX, IX, N
C***FIRST EXECUTABLE STATEMENT  SCASUM
      SCASUM = 0.0D0
      IF (N .LE. 0) RETURN
C
      IF (INCX .EQ. 1) GOTO 20
C
C     Code for increment not equal to 1.
C
      IX = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      DO 10 I = 1,N
        SCASUM = SCASUM + DABS(DREAL(CX(IX))) + DABS(DIMAG(CX(IX)))
        IX = IX + INCX
   10 CONTINUE
      RETURN
C
C     Code for increment equal to 1.
C
   20 DO 30 I = 1,N
        SCASUM = SCASUM + DABS(DREAL(CX(I))) + DABS(DIMAG(CX(I)))
   30 CONTINUE
      RETURN
      END



*DECK CGEDI
      SUBROUTINE CGEDI (A, LDA, N, IPVT, DET, WORK, JOB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C***BEGIN PROLOGUE  CGEDI
C***PURPOSE  Compute the determinant and inverse of a matrix using the
C            factors computed by CGECO or CGEFA.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D2C1, D3C1
C***TYPE      COMPLEX (SGEDI-S, DGEDI-D, CGEDI-C)
C***KEYWORDS  DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK, MATRIX
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     CGEDI computes the determinant and inverse of a matrix
C     using the factors computed by CGECO or CGEFA.
C
C     On Entry
C
C        A       COMPLEX(LDA, N)
C                the output from CGECO or CGEFA.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C        IPVT    INTEGER(N)
C                the pivot vector from CGECO or CGEFA.
C
C        WORK    COMPLEX(N)
C                work vector.  Contents destroyed.
C
C        JOB     INTEGER
C                = 11   both determinant and inverse.
C                = 01   inverse only.
C                = 10   determinant only.
C
C     On Return
C
C        A       inverse of original matrix if requested.
C                Otherwise unchanged.
C
C        DET     COMPLEX(2)
C                determinant of original matrix if requested.
C                Otherwise not referenced.
C                Determinant = DET(1) * 10.0**DET(2)
C                with  1.0 .LE. CABS1(DET(1)) .LT. 10.0
C                or  DET(1) .EQ. 0.0 .
C
C     Error Condition
C
C        A division by zero will occur if the input factor contains
C        a zero on the diagonal and the inverse is requested.
C        It will not occur if the subroutines are called correctly
C        and if CGECO has set RCOND .GT. 0.0 or CGEFA has set
C        INFO .EQ. 0 .
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  CAXPY, CSCAL, CSWAP
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  CGEDI
      INTEGER LDA,N,IPVT(*),JOB
      COMPLEX*16 A(LDA,*),DET(2),WORK(*),T
      DOUBLE PRECISION TEN
      INTEGER I,J,K,KB,KP1,L,NM1
      COMPLEX*16 ZDUM
      DOUBLE PRECISION CABS1
      CABS1(ZDUM) = DABS(DREAL(ZDUM)) + DABS(DIMAG(ZDUM))
C***FIRST EXECUTABLE STATEMENT  CGEDI
C
C     COMPUTE DETERMINANT
C
      IF (JOB/10 .EQ. 0) GO TO 70
         DET(1) = DCMPLX(1.0D0,0.0D0)
         DET(2) = DCMPLX(0.0D0,0.0D0)
         TEN = 10.0D0
         DO 50 I = 1, N
            IF (IPVT(I) .NE. I) DET(1) = -DET(1)
            DET(1) = A(I,I)*DET(1)
            IF (CABS1(DET(1)) .EQ. 0.0D0) GO TO 60
   10       IF (CABS1(DET(1)) .GE. 1.0D0) GO TO 20
               DET(1) = CMPLX(TEN,0.0D0)*DET(1)
               DET(2) = DET(2) - DCMPLX(1.0D0,0.0D0)
            GO TO 10
   20       CONTINUE
   30       IF (CABS1(DET(1)) .LT. TEN) GO TO 40
               DET(1) = DET(1)/CMPLX(TEN,0.0D0)
               DET(2) = DET(2) + DCMPLX(1.0D0,0.0D0)
            GO TO 30
   40       CONTINUE
   50    CONTINUE
   60    CONTINUE
   70 CONTINUE
C
C     COMPUTE INVERSE(U)
C
      IF (MOD(JOB,10) .EQ. 0) GO TO 150
         DO 100 K = 1, N
            A(K,K) = DCMPLX(1.0D0,0.0D0)/A(K,K)
            T = -A(K,K)
            CALL CSCAL(K-1,T,A(1,K),1)
            KP1 = K + 1
            IF (N .LT. KP1) GO TO 90
            DO 80 J = KP1, N
               T = A(K,J)
               A(K,J) = DCMPLX(0.0D0,0.0D0)
               CALL CAXPY(K,T,A(1,K),1,A(1,J),1)
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
C
C        FORM INVERSE(U)*INVERSE(L)
C
         NM1 = N - 1
         IF (NM1 .LT. 1) GO TO 140
         DO 130 KB = 1, NM1
            K = N - KB
            KP1 = K + 1
            DO 110 I = KP1, N
               WORK(I) = A(I,K)
               A(I,K) = DCMPLX(0.0D0,0.0D0)
  110       CONTINUE
            DO 120 J = KP1, N
               T = WORK(J)
               CALL CAXPY(N,T,A(1,J),1,A(1,K),1)
  120       CONTINUE
            L = IPVT(K)
            IF (L .NE. K) CALL CSWAP(N,A(1,K),1,A(1,L),1)
  130    CONTINUE
  140    CONTINUE
  150 CONTINUE
      RETURN
      END

*DECK CSWAP
      SUBROUTINE CSWAP (N, CX, INCX, CY, INCY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C***BEGIN PROLOGUE  CSWAP
C***PURPOSE  Interchange two vectors.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A5
C***TYPE      COMPLEX (SSWAP-S, DSWAP-D, CSWAP-C, ISWAP-I)
C***KEYWORDS  BLAS, INTERCHANGE, LINEAR ALGEBRA, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       CX  complex vector with N elements
C     INCX  storage spacing between elements of CX
C       CY  complex vector with N elements
C     INCY  storage spacing between elements of CY
C
C     --Output--
C       CX  input vector CY (unchanged if N .LE. 0)
C       CY  input vector CX (unchanged if N .LE. 0)
C
C     Interchange complex CX and complex CY
C     For I = 0 to N-1, interchange  CX(LX+I*INCX) and CY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
C     defined in a similar way using INCY.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  CSWAP
      COMPLEX*16 CX(*),CY(*),CTEMP
C***FIRST EXECUTABLE STATEMENT  CSWAP
      IF (N .LE. 0) RETURN
      IF (INCX.EQ.INCY .AND. INCX.GT.0) GO TO 20
C
C     Code for unequal or nonpositive increments.
C
      KX = 1
      KY = 1
      IF (INCX .LT. 0) KX = 1+(1-N)*INCX
      IF (INCY .LT. 0) KY = 1+(1-N)*INCY
      DO 10 I = 1,N
        CTEMP = CX(KX)
        CX(KX) = CY(KY)
        CY(KY) = CTEMP
        KX = KX + INCX
        KY = KY + INCY
   10 CONTINUE
      RETURN
C
C     Code for equal, positive, non-unit increments.
C
   20 NS = N*INCX
      DO 30 I = 1,NS,INCX
        CTEMP = CX(I)
        CX(I) = CY(I)
        CY(I) = CTEMP
   30 CONTINUE
      RETURN
      END


      SUBROUTINE F06AAZ ( SRNAME, INFO )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      INTEGER            INFO
      CHARACTER*13       SRNAME
C     ..
C
C  Purpose
C  =======
C
C  F06AAZ  is an error handler for the Level 2 BLAS routines.
C
C  It is called by the Level 2 BLAS routines if an input parameter is
C  invalid.
C
C  Parameters
C  ==========
C
C  SRNAME - CHARACTER*13.
C           On entry, SRNAME specifies the name of the routine which
C           called F06AAZ.
C
C  INFO   - INTEGER.
C           On entry, INFO specifies the position of the invalid
C           parameter in the parameter-list of the calling routine.
C
C
C  Auxiliary routine for Level 2 Blas.
C
C  Written on 20-July-1986.
C
C     .. Local Scalars ..
      INTEGER            IFAIL
      CHARACTER*80       REC (1)
C     .. External Functions ..
      INTEGER            P01ABF
      EXTERNAL           P01ABF
C     ..
C     .. Executable Statements ..
      WRITE (REC (1),99999) SRNAME, INFO
      IFAIL = 0
      IFAIL = P01ABF (IFAIL, -1, SRNAME(1:6), 1, REC)
C
      RETURN
C
99999 FORMAT ( ' ** On entry to ', A13, ' parameter number ', I2,
     $         ' had an illegal value' )
C
C     End of F06AAZ.
C
      END


