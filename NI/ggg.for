      PROGRAM GGG
C     ======================================================
C     THE PROGRAM CALCULATES THE DYNAMICAL GREEN'S FUNCTIONS
C     IN ALKALI HALIDES.
C     COPYRIGHT DMITRI NEVEDROV
C     INSTITUTE OF THEORETICAL PHYSICS, UNIVERSITY OF TARTU,
C     TAHE 4, EE2400 TARTU, ESTONIA
C     PHONE:  00372 7 428164
C     FAX:    00372 7 383033
C     E-MAIL: dmitri@eeter.fi.tartu.ee
C     HTTP:/eeter.fi.tartu.ee/~dmitri
C     ======================================================
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      INTEGER NR,NDIV,IMP
      DOUBLE PRECISION STEP,GIM(37,1670),GRE(36,1670),DEN(1670),FM(2)
      CHARACTER*4,CR
      CHARACTER*80,SOURCE
 10   CR='NAI'
      CALL INPCR(CR,SOURCE,FM,IOUT)
c      print*,fm
c      stop
      NDIV=100
      IMP=2
      NEIG=1
      STEP=0.125D0
      CALL GFUN(36,1670,NDIV,IMP,NEIG,STEP,NR,GIM,GRE,DEN,FM(1))
      CALL CONVER(NR,FM(1),GRE,GIM)
      END
      
      SUBROUTINE DIAG(EVAL,EVEC,GR,Q,FFF,GGG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      INTEGER NM
      PARAMETER (NM=6)
      DOUBLE PRECISION A(NM,NM),B(NM,NM),C(NM,NM),D(NM,NM),E(NM,NM),
     *     G(NM,NM),F(NM,NM),BC(NM,NM),CB(NM,NM),EVAL(NM),EVAL1(NM),
     *     EVEC(NM,NM),SM(NM,NM),GR(3,NM),Q(3),C1(4,NM),C2(4,NM),
     *     R1(4,NM),R2(4,NM),R3(4,NM),VEC(NM),FFF(50),GGG(50)
      INTEGER IFAIL,IZ,I,J,L,K
      PARAMETER (IZ=6)
      EXTERNAL F01AAF,F01CKF,F01AJF,F02AMF
      DOUBLE PRECISION HP(12),AK1,AK2,HR(6),AM1,AM2,AM12,
     *     ZZ,Y1Z,Y2Z,YY1,YY2,YY12,Z(IZ),ATOL,EPS,X02AJF
      ATOL=0.0D0
      EPS=X02AJF()
      OPEN(346,FILE='BBB')
      READ(346,222)(HP(I),I=1,12),AK1,AK2,(HR(I),I=1,6),
     *     AM1,AM2,AM12,ZZ,Y1Z,Y2Z,YY1,YY2,YY12
      CLOSE(346)
      CALL COUL(C1,C2,Q,1.6D0,FFF,GGG)
      CALL REPUL(R1,R2,R3,Q)
      DO 180 K=1,4
         VEC(K)=0.0D0
         DO 80 I=1,6
            DO 80 J=I,6
               IF (I-3) 5,5,50
 5             IF (J-3) 10,10,30
 10            IF (I-J) 20,15,20
 15            L=I
               VEC(1)=AK1
               GOTO 25
 20            L=I+J+1
               VEC(1)=0.0D0
 25            D2=Y1Z
               D3=YY1
               V1=R1(K,L)
               GOTO 70
 30            IF (J-I-3) 40,35,40
 35            L=I
               GOTO 45
 40            L=I+J-2
 45            V1=R3(K,L)
               P1=C2(K,L)
               B(I,J)=V1+Y2Z*P1
               B(J,I)=V1-Y1Z*P1
               D1=-ZZ
               D3=YY12
               VEC(1)=0.0D0
               GOTO 75
 50            IF(I-J) 60,55,60
 55            L=I-3
               VEC(1)=AK2
               GOTO 65
 60            L=I+J-5
               VEC(1)=0.0d0
 65            V1=R2(K,L)
               D2=-Y2Z
               D3=YY2
 70            D1=ZZ
               P1=-C1(K,L)
               V2=V1+D2*P1
               B(I,J)=V2
               B(J,I)=V2
 75            V2=V1+D1*P1
               A(I,J)=V2
               A(J,I)=V2
               V2=V1+VEC(K)+D3*P1
               C(I,J)=V2
               C(J,I)=V2
 80         CONTINUE
            DO I=1,6
               DO J=1,6
                  D(I,J)=B(J,I)
               ENDDO
            ENDDO
            IF (K-1) 90,90,100
 90         CALL F01AAF(C,NM,NM,F,NM,EVAL,IFAIL)
            DO I3=1,6
               DO J3=1,6
                  C(I3,J3)=F(I3,J3)
               ENDDO
            ENDDO
            IFAIL=1
            CALL F01CKF(CB,C,D,NM,NM,NM,Z,IZ,1,IFAIL)
            CALL F01CKF(D,B,CB,NM,NM,NM,Z,IZ,1,IFAIL)
            DO I3=1,6
               DO J3=1,6
                  BC(I3,J3)=CB(J3,I3)
               ENDDO
            ENDDO
            GOTO 110
 100        CALL F01CKF(F,B,CB,NM,NM,NM,Z,IZ,1,IFAIL)
            CALL F01CKF(E,BC,D,NM,NM,NM,Z,IZ,1,IFAIL)
            CALL F01CKF(D,C,CB,NM,NM,NM,Z,IZ,1,IFAIL)
            CALL F01CKF(G,BC,D,NM,NM,NM,Z,IZ,1,IFAIL)
            DO I3=1,NM
               DO J3=1,NM
                  D(I3,J3)=F(I3,J3)+E(I3,J3)-G(I3,J3)
               ENDDO
            ENDDO
 110        DO 150 I=1,NM
               DO 150 J=1,I
                  IF (I-3) 120,120,115
 115              IF (J-3) 125,125,130
 120              V1=AM1
                  GOTO 135
 125              V1=AM12
                  GOTO 135
 130              V1=AM2
 135              V2=V1*(A(I,J)-D(I,J))
                  IF (K.GT.1) GOTO 140
                  SM(I,J)=V2
                  SM(J,I)=V2
                  GOTO 150
 140              D(I,J)=V2
                  D(J,I)=V2
 150           CONTINUE
               IF (K.GT.1) GOTO 160
               CALL F01AJF(6,ATOL,SM,6,EVAL,EVAL1,EVEC,6)
               CALL F02AMF(6,EPS,EVAL,EVAL1,EVEC,6,IFAIL)
               GOTO 180
 160           DO 175 I=1,6
                  DO 165 J=1,6
 165                 VEC(J)=EVEC(J,I)
                     D1=0.0d0
                     DO 170 J=1,6
                        DO 170 M=1,6
 170                       D1=D1+VEC(J)*D(J,M)*VEC(M)
 175                       GR(K-1,I)=D1
 180                    CONTINUE
                        RETURN
 221                    FORMAT(1E16.8)
 222                    FORMAT(29D30.20)
 224                    FORMAT(6E12.3)
                        END
      
      SUBROUTINE GFUN(MF,NRES,NDIV,IMP,NEIGH,STEP,NR,GIM,GRE,DEN,FM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      INTEGER J,I,NR,NRES,NDIV,IMP,NEIGH,MF,M1,M2,M3
      DOUBLE PRECISION STEP,FM,GILTR,EVAL(6),EVEC(6,6),GG(10),WW(9),CC,
     *     GIM(MF,NRES),GRE(MF,NRES),DEN(NRES),Q(3),SI(3),CO(3),VI(3),
     *     VN(3),G(3),GRAD(3,6),F(37),DF(37,3),FFF(50),GGG(50),AIM,W,X,
     *     X01AAF
      CALL LATSUM(1.6D0,FFF,GGG)
      B=0.5D0/DBLE(NDIV)
      NR=INT(FM/STEP)+3
      IF (((NR/2)*2).EQ.NR) NR=NR+1
      A=0.018836537D0*STEP
      DO I=1,MF
		DO J=1,NR
			GIM(I,J)=0.0D0
		ENDDO
      ENDDO
      L11=3*(2-IMP)
      L22=3*(IMP-1)
      N1=3*NDIV
      N2=N1+5
      AN=0.25d0/DBLE(NDIV**3)
      DO 80 K1=1,NDIV
      PRINT*,k1
      Q(1)=B*DBLE(2*K1-1)
      DO 75 K2=1,K1
c      IF ((2*(K1+K2)).GE.N2) GOTO 80
      Q(2)=B*DBLE(2*K2-1)
      DO 70 K3=1,K2
      KK=2*(K1+K2+K3)
c      IF (KK.GT.N2) GOTO 75
      Q(3)=B*DBLE(2*K3-1)
C      WEIGHT=VES(K1,K2,K3)
      WEIGHT=1.0d0
      IF (KK.GT.N1) WEIGHT=WEIGHT/2.0d0
      WEIGHT=WEIGHT*AN/24.0d0
      CALL DIAG(EVAL,EVEC,GRAD,Q,FFF,GGG)
      DO 20 I=1,3
      C=X01AAF(X)*(2.0D0*Q(I)-1.0D0)
      SI(I)=DSIN(C)
20    CO(I)=DCOS(C)
      DO 65 K=1,6
      E=(EVAL(K))
      DO 25 I=1,3
      LA=I+L11
      LB=I+L22
      VN(I)=EVEC(LA,K)
25    VI(I)=EVEC(LB,K)
      CALL FUNTRI(NEIGH,F,DF,VN,VI,SI,CO)
      GRM=DSQRT(GRAD(1,K)**2+GRAD(2,K)**2+GRAD(3,K)**2)
      DO 30 I=1,3
30    G(I)=B*GRAD(I,K)/GRM
      CALL AREB(G,M1,M2,M3,GG,WW,CC)
      DO 35 I=1,MF
      VN(1)=DF(I,M1)
      VN(2)=DF(I,M2)
      VN(3)=DF(I,M3)
      DO 35 J=1,3
35    DF(I,J)=VN(J)
      A1=WEIGHT/GRM
      X1=GRM*(DABS(G(1))+DABS(G(2))+DABS(G(3)))
      X2=1.0d0/A
      A2=E-X1
      IF (A2) 45,45,40
40    IF (A2.LT.0.0D0) THEN
      ENDIF
      LOW=X2*DSQRT(A2)
      IF (LOW-1) 45,50,50
45    LOW=1
50    LUP=X2*DSQRT((E+X1))+2.0D0
      IF (LUP.GT.NR) LUP=NR
      DO 65 I=LOW,LUP
      W=(((A*DBLE(I-1))**2)-E)/GRM
      CALL AREA(AIM,VI,W,G,GG,WW,CC)
      DO 65 LA=1,MF
		X1=0.0d0
		DO 60 LB=1,3
			W=B*DF(LA,LB)
60		X1=X1+W*VI(LB)
		GIM(LA,I)=GIM(LA,I)+A1*(AIM*F(LA)+X1)
65    CONTINUE
70    CONTINUE
75    CONTINUE
80    CONTINUE
      DO 100 I=1,MF
      DO 90 J=1,NR
      X1=2.0D0*(J-1)
      X2=GIM(I,J)
      DEN(J)=X1*X2
90    GIM(I,J)=X01AAF(X)*X2
      DO 95 J=1,NR
95    GRE(I,J)=GILTR(DEN,NR,J*2)
100   CONTINUE
      RETURN
      END

      DOUBLE PRECISION FUNCTION GILTR(F,N,L)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      INTEGER N,L
      DOUBLE PRECISION f(n),sum(2)
      IF(L-N) 5,15,15
 5    IF(L-1) 15,15,10
10    P=DBLE(L-1)
      P=0.5D0/P
      A=F(L)
      B=DBLE(N-L)
      C=DBLE(N+L-2)
      R=A*P*DLOG(B/C)
      GOTO 20
15    A=0.0D0
      R=0.0D0
20    N1=N-1
      DO 40 M=1,2
      B=0.0D0
      M1=M+1
      DO 35 K=M1,N1,2
      IF (K-L) 30,25,30
25    B=B+0.5D0*P*(F(K+1)-F(K-1))
      GOTO 35
30    C=DBLE((K-L)*(K+L-2))
      B=B+(F(K)-A)/C
35    continue
40    sum(m)=b
      if(l.eq.1) goto 45
      b1=0.5d0*a/(DBLE(l-1)**2)
      goto 50
45    b1=0.25d0*(f(3)-2.0d0*f(2))
50    if(l.eq.n) goto 55
      b2=0.5d0*(f(n)-a)/DBLE((n-l)*(n+l-2))
      goto 60
55    b2=0.125d0*(3.0d0*f(n)-4.0d0*f(n-1)+f(n-2))/DBLE(N1)
60    giltr=-r-2.0d0*(2.0d0*sum(1)+sum(2)+b1+b2)/3.0d0
      return
      end

      SUBROUTINE INPCR(CRYST,SOURCE,FM,IOUT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      INTEGER IOUT,J,K
      DOUBLE PRECISION FM(2),DAT(14,25),FMM(25),ALG1(10),ALG2(4)
      CHARACTER*80 S(25), SOURCE
      CHARACTER*4 CR(25), CRYST
      DATA CR /'LI7H', 'LI6H', 'LI7D', 'LI6D',
     *         'LIF','NAF', 'NACL', 'NABR', 'NAI',
     *         'KF', 'KCL', 'KBR', 'KI',
     *         'RBF', 'RBCL', 'RBBR', 'RBI',
     *         'CSF','MGO', 'CAO', 'SRO', 4*'XXXX'/
      DATA (S(J),J=1,10)/
     *'4 by J.L.Verbie et al. Phys. Rev. 168, 980 (1980) RT',
     *'4 by J.L.Verbie et al. Phys. Rev. 168, 980 (1980) RT',
     *'4 by J.L.Verbie et al. Phys. Rev. 168, 980 (1980) RT',
     *'4 by J.L.Verbie et al. Phys. Rev. 168, 980 (1980) RT',
     *'2 by G.Dolling et al. Phys. Rev. 168, 970 (1968) T=300K',
     *'6 by W.J.L.Buyers, Phys. Rev. 153, 923 (1967) T=295K',
     *'2 by G.Raunio & S.Rolandson, Phys. Rev. B2, 2098 (1970)
     * T=80K',
     *'2 by J.S.Reid, T.Smith & W.J.L.Buyers, Phys. Rev. B1,
     * 1833 (1970) RT',
     *'6 by R.A.Cowley et al. Phys. Rev. 131, 1030 (1963) 100K',
     *'4 by W.Buhrer, phys. stat. sol. 41, 789 (1970) T=300K'/
      DATA (S(J),J=11,24)/
     *'6 by J.R.D.Copley et al. Phys. Rev. 182, 965 (1969) T=115K',
     *'6 by R.A.Cowley et al. Phys. Rev. 131, 1030 (1963) 100K', 
     *'3 by G.Dolling et al. Phys. Rev. 147, 577 (1966) T=95K',
     *'2 by G.Raunio & S.Rolandson, Phys. Rev. B2, 2098 (1970) T=80K',
     *'2 by G.Raunio & S.Rolandson, Phys. Rev. B2, 2098 (1970) T=80K',
     *'4 by S.Rolandson & G.Raunio, J. Phys. C 4, 958 (1971) T=80K',
     *'4 by G.Raunio & S.Rolandson, phys. stat. sol., 40, 749 (1970)
     * T=80K',
     *'4 by S.Rolandson, phys. stat. sol. 52, 643 (1972) T=80K',
     *'8 by M.J.L.Snagster et al. J. Phys. C3, 1026 (1970) T=293K',
     *'2 by D.H.Saunderson & G.Pecham, J. Phys. C 4, 2009 (1971) T=293K',
     *'4 by K.H.Rieder et al. Phys. Rev. B12, 3374 (1975) T=300K',
     *4*'XXXX'/

      DATA ((DAT(J,K),J=1,14),K=1,9)/
     *2.0418d0,7.0d0,1.0d0,6.146d0,-0.748d0,
     *0.0d0,0.0d0,-0.4d0,0.253d0,0.878d0,0.0d0,0.0d0,1.55d0,0.362d0,
     *2.0418d0,6.0d0,1.0d0,6.146d0,-0.748d0,
     *0.0d0,0.0d0,-0.4d0,0.253d0,0.878d0,0.0d0,0.0d0,1.55d0,0.362d0,
     *2.0346d0,7.0d0,2.0d0,6.146d0,-0.748d0,
     *0.0d0,0.0d0,-0.4d0,0.253d0,0.878d0,0.0d0,0.0d0,1.55d0,0.362d0,
     *2.0346d0,6.d0,2.d0,6.146d0,-0.748d0,
     *0.0d0,0.0d0,-0.4d0,0.253d0,0.878d0,0.0d0,0.0d0,1.55d0,0.362d0,
     *2.009d0,6.939d0,18.9984d0,6.797d0,-0.653d0,
     *0.0d0,0.0d0,1.048d0,-0.049d0,0.902d0,0.0d0,
     *0.0d0,0.725d0,0.1115d0,
     *2.310d0,22.9898d0,18.9984d0,9.26d0,-0.77d0,
     *0.0d0,0.0d0,0.34d0,-0.02d0,0.907d0,0.27d0,0.01d0,0.7d0,0.116d0,
     *2.8000d0,22.9898d0,35.4530d0,9.77d0,-0.864d0,
     *0.56d0,-0.006d0,-0.006d0,0.004d0,0.89d0,
     *0.18d0,-0.031d0,2.35d0,0.131d0,
     *2.9866d0,22.9898d0,79.904d0,11.39d0,-1.04d0,
     *0.0d0,0.0d0,0.0d0,0.0d0,1.0d0,0.91d0,-0.16d0,3.94d0,0.19d0,
     *3.208d0,22.9898d0,126.9044d0,9.94d0,-0.867d0,
     *0.0d0,0.0d0,0.618d0,-0.041d0,0.89d0,1.98d0,-0.112d0,
     *4.34d0,0.136d0/
      DATA ((DAT(J,K),J=1,14),K=10,17)/
     *2.6735d0,39.102d0,18.9984d0,10.76d0,-1.130d0,
     *0.317d0,-0.054d0,-0.247d0,0.085d0,0.961d0,1.233d0,
     *-0.060d0,0.455d0,0.04d0,
     *3.12d0,39.102d0,35.453d0,12.12d0,-1.095d0,
     *0.0d0,0.0d0,-0.1d0,0.058d0,0.918d0,1.701d0,
     *-0.025d0,2.126d0,0.121d0,
     *3.293d0,39.102d0,79.904d0,13.15d0,-1.227d0,
     *0.0d0,0.0d0,-0.399d0,0.054d0,0.965d0,2.12d0,
     *-0.101d0,3.05d0,0.141d0,
     *3.51d0,39.102d0,126.9044d0,13.4d0,-1.0d0,
     *-0.16d0,-0.01d0,-0.29d0,0.05d0,0.92d0,2.28d0,
     *-0.11d0,4.51d0,0.13d0,
     *2.805d0,85.47d0,18.9984d0,12.31d0,-1.416d0,
     *-0.17d0,0.05d0,0.02d0,0.0003d0,1.0d0,1.54d0,
     *-0.008d0,1.15d0,0.099d0,
     *3.27d0,85.47d0,35.453d0,12.01d0,-0.827d0,
     *0.31d0,-0.01d0,0.36d0,-0.04d0,0.79d0,0.66d0,
     *0.066d0,2.31d0,0.103d0,
     *3.42d0,85.47d0,79.904d0,12.0d0,-0.634d0,
     *0.5d0,-0.18d0,0.11d0,0.036d0,0.79d0,3.04d0,
     *0.011d0,0.56d0,0.042d0,
     *3.629d0,85.47d0,126.9044d0,13.82d0,-0.938d0,
     *-0.28d0,0.19d0,0.64d0,-0.21d0,0.87d0,2.29d0,
     *-0.099d0,5.64d0,0.045/
      DATA ((DAT(J,K),J=1,14),K=18,25)/
     *2.984d0,132.905d0,18.9984d0,11.56d0,-1.16d0,
     *0.78d0,-0.29d0,-0.38d0,0.16d0,0.946d0,
     *2.48d0,-0.052d0,0.46d0,0.054d0,
     *2.1065d0,24.312d0,15.9994d0,29.13d0,-3.51d0,
     *2*0.0d0,-1.353d0,0.256d0,1.885d0,2*0.0d0,1.76d0,0.685d0,
     *2.4d0,40.08d0,15.9994d0,35.117d0,-4.482d0,
     *2*0.0d0,-2.143d0,0.693d0,2.0d0,0.32d0,0.0391d0,2.51d0,0.817d0,
     *2.58d0,87.62d0,15.9994d0,31.97d0,-4.224d0,
     *4*0.0d0,1.78d0,0.198d0,0.095d0,1.615d0,0.522d0,
     *56*0.0d0/

      DATA FMM/8*0.0D0,171.48D0,16*0.0D0/
      DO K=1,25
      IF (CRYST.EQ.CR(K)) THEN
      GOTO 30
      ENDIF
      ENDDO
      IOUT=1
      PRINT 1, CRYST, CR
      RETURN
30    IOUT=0
      SOURCE=S(K)
      DO J=1,10
      ALG1(J)=DAT(J,K)
      ENDDO
      DO J=1,4
      ALG2(J)=DAT(J+10,K)
      ENDDO
      CALL ALG(FM,ALG1,ALG2)
      IF (FMM(K).NE.0.0D0) FM(1)=FMM(K)
      RETURN
1     FORMAT(//'DATA FOR CRYSTAL ',A4,' ARE NOT FOUND'/
     *      'THE FOLLOWING CRYSTALS CAN BE USED:'/
     *      /2(10A6/),5A6//'REPEAT INPUT'/)
      END

      SUBROUTINE ALG(F0,ALG1,ALG2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DOUBLE PRECISION R0,AM1,AM2,A12,B12,A11,B11,A22,B22,Z,Y1,Y2,AK1,
     *     AK2,F1,F2,G1,G2,H1,H2,R1,R2,R12,ZZ,YZ1,YZ2,YY1,YY2,YY12,ALP1,
     *     ALP2,D1,D2,F0(2),C1(2),ALG1(10),ALG2(4),C
      DATA C1/8.377580d0, -4.188790d0/
      R0=ALG1(1)
      AM1=ALG1(2)
      AM2=ALG1(3)
      A12=ALG1(4)
      B12=ALG1(5)
      A11=ALG1(6)
      B11=ALG1(7)
      A22=ALG1(8)
      B22=ALG1(9)
      Z=ALG1(10)
      ALP1=ALG2(1)
      D1=ALG2(2)
      ALP2=ALG2(3)
      D2=ALG2(4)
      C=A12+2.0d0*B12
      F1=C+2.0d0*A11+4.0d0*B11
      F2=C+2.0d0*A22+4.0d0*B22
      G1=A11+B11
      G2=A22+B22
      H1=A11-B11
      H2=A22-B22
      P=1.0D0/(R0*R0*R0)
      P1=695.6167d0*P
      R1=P1/AM1
      R2=P1/AM2
      R12=P1/DSQRT(AM1*AM2)
      P=0.5d0*P
      IF (ALP1.NE.0.0d0) GOTO 5
      Y1=0.0d0
      AK1=1.0d8
      GOTO 10
5     Y1=-P*ALP1*C/D1
      AK1=-C*(1.0D0+Y1/D1)
10    IF (ALP2.NE.0.0D0) GOTO 15
      Y2=0.0d0
      AK2=1.0d9
      GOTO 20
15    Y2=-P*ALP2*C/D2
      AK2=-C*(1.0D0+Y2/D2)
20    DEL=AK1*AK2+C*(AK1+AK2)
      CA=C*(1.0d0-C*(AK1+AK2)/DEL)
      ZA=Z+C*(Y2*AK1-Y1*AK2)/DEL
      ALP=(AK2*Y1*Y1+AK1*Y2*Y2+C*(Y1+Y2)**2)/DEL
      DO 25 K=1,2
      S=C1(K)
25    F0(K)=53.0883d0*DSQRT((R1+R2)*(CA+S*ZA*ZA/(1.0d0+S*ALP)))
      ZZ=Z*Z
      YZ1=Y1*Z
      YZ2=Y2*Z
      YY1=Y1*Y1
      YY2=Y2*Y2
      YY12=Y1*Y2
      OPEN(345,FILE='BBB')
      WRITE(345,222)R0,AM1,AM2,A12,B12,A11,B11,A22,B22,Z,
     *      Y1,Y2,AK1,AK2,F1,F2,G1,G2,H1,H2,R1,R2,
     *      R12,ZZ,YZ1,YZ2,YY1,YY2,YY12
      CLOSE(345)
      OPEN(245,FILE='CCC')
      WRITE(245,111)ALP1,D1,ALP2,D2
      CLOSE(245)
111   FORMAT(4D30.20)
222   FORMAT(29D30.20)
      RETURN
      END

      subroutine area(s,es,p,x,G,W,C)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DOUBLE PRECISION ES(3),X(3),G(10),W(9),C,S,P
      P1=DABS(p)
      if(p1.ge.w(4)) goto 55
      p2=p*p
      p3=p1*p2
      if(p1.gt.w(1)) goto 15
      if(c) 10,5,5
5     s=4.0d0/g(1)
      es(1)=p1*s/g(1)
      es(2)=0.0d0
      es(3)=0.0d0
      goto 30
10    s=g(7)*(4.0d0*(g(4)+g(5)+g(6))-p2-w(9))
      r=p1*g(7)/3.0d0
      do 12 k=1,3
12    es(k)=r*(3.0d0*w(8-k)-p2)/g(k)
      goto 30
15    if(p1.gt.w(2)) goto 20
      s=0.5D0*g(7)*(-p2+2.0d0*c*w(1)*p1-w(9)+
     *4.0D0*(2.0D0*g(6)+g(5)+g(4)))
      r =g(7)/6.0d0
      es(1)=r*(-p3-3.0d0*p2*(g(2)+g(3))+
     *3.0d0*p1*(w(7)+4.0d0*g(6))-
     *   w(8)*(w(4)+g(1)))/g(1)
      es(2)=r*(-p3+3.0d0*p2*(g(1)-g(3))-3.0d0*c*w(6)*p1+
     *   w(8)*(w(3)+g(2)))/g(2)
      es(3)=r*(-p3+3.0d0*p2*(g(1)-g(2))-3.0d0*c*w(5)*p1+
     *w(8)*(w(2)+g(3)))/g(3)
      goto 30
20    if(p1.gt.w(3)) goto 25
      s=2.0d0*(g(1)+g(2)-p1)/g(4)
      r=1.0d0/(3.0d0*g(4))
      r1=3.0d0*(g(8)-g(9))-g(10)
      es(1)=r*(-3.0d0*p2+6.0d0*g(2)*p1+r1)/g(1)
      es(2)=r*(-3.0d0*p2+6.0d0*g(1)*p1-r1-2.0d0*g(10))/g(2)
      es(3)= 2.0d0*r*g(3)
      goto 30
25    s=0.5d0*g(7)*(p1-w(4))**2
      r=s/3.0d0
      es(1)=r*(p1+c*w(1)+g(1))/g(1)
      es(2)=r*(p1-w(2)+g(2))/g(2)
      es(3)=r*(p1-w(3)+g(3))/g(3)
30    if(p.ge.0.0d0) goto 40
      do 35 i=1,3
35    es(i)=-es(i)
40    do 50 i=1,3
      if(x(i)) 45,50,50
45    es(i)=-es(i)
50    continue
      goto 65
55    s=0.0d0
      do 60 i=1,3
60    es(i)=0.0d0
65    return
      end

      subroutine areb(x,m1,m2,m3,G,W,C)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      INTEGER M1,M2,M3,K
      DOUBLE PRECISION X(3),G(10),W(9),C
      do 5 k=1,3
5     g(k)=DABS(x(k))
      m1=1
      m3=1
      x1=g(1)
      y=x1
      do 25 k=2,3
      a=g(k)
      if(a-x1) 15,15,10
10    m1=k
      x1=a
      goto 25
15    if(a-y) 20,25,25
20    m3=k
      y=a
25    continue
      do 30 k=1,3
      if(k.eq.m1) goto 30
      if(k.eq.m3) goto 30
      M2=k
      g(2)=g(k)
30    continue
      g(1)=x1
      g(3)=y
      do 35 k=1,3
35    g(k+7)=g(k)**2
      g(4)=g(1)*g(2)
      g(5)=g(1)*g(3)
      g(6)=g(2)*g(3)
      b=g(1)*g(2)*g(3)
      if(b.ne.0.0d0) g(7)=1.0d0/b
      a=g(1)-g(2)-g(3)
      if(a) 40,45,45
40    c=-1.0d0
      goto 50
45    c=1.0d0
50    w(1)=DABS(a)
      w(2)=g(1)-g(2)+g(3)
      w(3)=g(1)+g(2)-g(3)
      w(4)=g(1)+g(2)+g(3)
      w(5)=w(1)*w(2)
      w(6)=w(1)*w(3)
      w(7)=w(2)*w(3)
      w(8)=w(1)**2
      w(9)=w(4)**2
      return
      end

      SUBROUTINE FUNTRI(NEIGH,F,DF,E1,E2,S,C)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      INTEGER NEIGH
      DOUBLE PRECISION X,F(37),DF(36,3),E1(3),E2(3),C(3),S(3),
     * X01AAF,EC(3,3),ES(3,3),TC(3,3),TS(3,3),P,FN(13)
      P=X01AAF(X)
      DO 10 I=1,3
      X=C(I)
      Y=S(I)
      DO 10 K=1,3
      V=E1(K)
      EC(K,I)=V*X
      ES(K,I)=V*Y
      IF (NEIGH.EQ.1) GOTO 10
      Z=E2(K)
      TC(K,I)=Z*(X*X-Y*Y)
      TS(K,I)=2.0d0*Z*X*Y
10    CONTINUE

      DO 170 I=1,1
      GOTO 161
C      GOTO (20,25), NEIGH
C20    GOTO (30,40,70,85,90,95,100,110,120,135,140,155,161), I
C25    GOTO (40,50,60,70,75,80,100,110,125,135,145,160,161), I

30    F(1)=0.50D0
      DO 35 K=1,3
35    DF(1,K)=0.0D0
      GOTO 170

40    X=ES(1,1)+ES(2,2)+ES(3,3)
      F(I)=2.0D0*X*X
      DO 45 K=1,3
45    DF(I,K)=4.0D0*X*EC(K,K)
      GOTO 170

50    Y=TS(1,1)+TS(2,2)+TS(3,3)
      F(I)=2.0D0*Y*Y
      DO 55 K=1,3
55    DF(I,K)=8.0D0*Y*TC(K,K)
      GOTO 170

60    F(3)=2.0D0*X*Y
      DO 65 K=1,3
65    DF(3,K)=2.0D0*EC(K,K)*Y+4.0D0*TC(K,K)*X
      GOTO 170

70    X=ES(1,1)-ES(2,2)
      Y=ES(1,1)-ES(3,3)
      Z=ES(2,2)-ES(3,3)
      F(I)=X*X+Y*Y+Z*Z
      DF(I,1)=2.0D0*EC(1,1)*(X+Y)
      DF(I,2)=2.0D0*EC(2,2)*(Z-X)
      DF(I,3)=2.0D0*EC(3,3)*(-Y-Z)
      GOTO 170

75    X1=TS(1,1)-TS(2,2)
      Y1=TS(1,1)-TS(3,3)
      Z1=TS(2,2)-TS(3,3)
      F(5)=X1*X1+Y1*Y1+Z1*Z1
      DF(5,1)=4.0d0*TC(1,1)*(X1+Y1)
      DF(5,2)=4.0d0*TC(2,2)*(Z1-X1)
      DF(5,3)=4.0d0*TC(3,3)*(-Y1-Z1)
      GOTO 170

80    F(6)=X*X1+Y*Y1+Z*Z1
      DF(6,1)=EC(1,1)*(X1+Y1)+2.0d0*TC(1,1)*(X+Y)
      DF(6,2)=EC(2,2)*(Z1-X1)+2.0d0*TC(2,2)*(Z-X)
      DF(6,3)=-EC(3,3)*(Y1+Z1)-2.0d0*TC(3,3)*(Y+Z)
      GOTO 170

85    X=ES(1,2)+ES(2,1)
      Y=ES(1,3)+ES(3,1)
      Z=ES(2,3)+ES(3,2)
      F(4)=X*X+Y*Y+Z*Z
      DF(4,1)=2.0d0*(EC(2,1)*X+EC(3,1)*Y)
      DF(4,2)=2.0d0*(EC(1,2)*X+EC(3,2)*Z)
      DF(4,3)=2.0d0*(EC(1,3)*Y+EC(2,3)*Z)
      GOTO 170

90    X=ES(1,2)-ES(2,1)
      Y=ES(1,3)-ES(3,1)
      Z=ES(2,3)-ES(3,2)
      F(5)=X*X+Y*Y+Z*Z
      DF(5,1)=-2.0d0*(EC(2,1)*X+EC(3,1)*Y)
      DF(5,2)=2.0d0*(EC(1,2)*X-EC(3,2)*Z)
      DF(5,3)=2.0d0*(EC(1,3)*Y+EC(2,3)*Z)
      GOTO 170

95    X=EC(1,2)-EC(1,3)
      Y=EC(2,1)-EC(2,3)
      Z=EC(3,1)-EC(3,2)
      F(6)=X*X+Y*Y+Z*Z
      DF(6,1)=-2.0d0*(ES(2,1)*Y+ES(3,1)*Z)
      DF(6,2)=2.0d0*(ES(3,2)*Z-ES(1,2)*X)
      DF(6,3)=2.0d0*(ES(1,3)*X+ES(2,3)*Y)
      GOTO 170

100   F(7)=E2(1)**2+E2(2)**2+E2(3)**2
      DO 105  K=1,3
105   DF(7,K)=0.0d0
      GOTO 170

110   F(8)=2.0d0*(EC(1,1)**2+EC(2,2)**2+EC(3,3)**2)
      DO 115 K=1,3
115   DF(8,K)=-4.0d0*EC(K,K)*ES(K,K)
      GOTO 170

120   X=EC(1,2)+EC(1,3)
      Y=EC(2,1)+EC(2,3)
      Z=EC(3,1)+EC(3,2)
      F(9)=X*X+Y*Y+Z*Z
      DF(9,1)=-2.0d0*(ES(2,1)*Y+ES(3,1)*Z)
      DF(9,2)=-2.0d0*(ES(1,2)*X+ES(3,2)*Z)
      DF(9,3)=-2.0d0*(ES(1,3)*X+ES(2,3)*Y)
      GOTO 170

125   F(9)=2.0d0*(TC(1,1)**2+TC(2,2)**2+TC(3,3)**2)
      DO 130 K=1,3
130   DF(9,K)=-8.0d0*TC(K,K)*TS(K,K)
      GOTO 170

135   X1=E1(1)*E2(1)
      Y1=E1(2)*E2(2)
      Z1=E1(3)*E2(3)
      F(10)=P*(X1*C(1)+Y1*C(2)+Z1*C(3))
      DF(10,1)=-P*X1*S(1)
      DF(10,2)=-P*Y1*S(2)
      DF(10,3)=-P*Z1*S(3)
      GOTO 170

140   F(11)=X1*(C(2)+C(3))+Y1*(C(1)+C(3))+Z1*(C(1)+C(2))
      DF(11,1)=-S(1)*(Y1+Z1)
      DF(11,2)=-S(2)*(X1+Z1)
      DF(11,3)=-S(3)*(X1+Y1)
      GOTO 170

145   F(11)=P*(E2(1)*TC(1,1)+E2(2)*TC(2,2)+E2(3)*TC(3,3))
      DO 150 K=1,3
150   DF(11,K)=-2.0d0*P*TS(K,K)*E2(K)
      GOTO 170

155   F(12)=P*(EC(1,1)*X+EC(2,2)*Y+EC(3,3)*Z)
      DF(12,1)=-P*(ES(1,1)*X+ES(2,1)*EC(2,2)+ES(3,1)*EC(3,3))
      DF(12,2)=-P*(ES(2,2)*Y+ES(1,2)*EC(1,1)+ES(3,2)*EC(3,3))
      DF(12,3)=-P*(ES(3,3)*Z+ES(1,3)*EC(1,1)+ES(2,3)*EC(2,2))
      GOTO 170

 160  F(12)=2.0d0*(EC(1,1)*TC(1,1)+EC(2,2)*TC(2,2)+EC(3,3)*TC(3,3))
      DO 165 K=1,3
 165     DF(12,K)=-2.0d0*(ES(K,K)*TC(K,K)+EC(K,K)*TS(K,K))
         
 161  F(1)=E2(1)*E2(1)
      F(2)=E2(1)*E2(1)*(C(1)*C(1)-S(1)*S(1))
      F(3)=E2(1)*E2(2)*(C(1)*C(1)-S(1)*S(1))
      F(4)=E2(1)*E2(2)*C(1)*C(2)
      F(5)=E2(1)*E2(1)*C(1)*C(2)
      F(6)=E2(1)*E2(3)*C(1)*C(2)
      F(7)=E2(1)*E2(1)*(C(2)*C(2)-S(2)*S(2))
      F(8)=E2(1)*E2(3)*(C(2)*C(2)-S(2)*S(2))
      F(9)=E2(1)*E2(1)*C(2)*C(3)
      F(10)=E1(1)*E2(1)*C(1)
      F(11)=E2(1)*E1(2)*C(1)
      F(12)=E1(1)*E1(1)

      F(13)=E2(1)*E1(2)*C(2)
      F(14)=E2(1)*E1(1)*C(2)
      F(15)=E2(1)*E1(3)*C(2)

      F(16)=E2(1)*E1(3)*C(3)
      F(17)=E2(1)*E1(1)*C(3)
      F(18)=E2(1)*E1(2)*C(3)

      F(19)=E2(2)*E1(1)*C(1)
      F(20)=E2(2)*E1(2)*C(1)
      F(21)=E2(2)*E1(3)*C(1)

      F(22)=E2(2)*E1(2)*C(2)
      F(23)=E2(2)*E1(1)*C(2)
      F(24)=E2(2)*E1(3)*C(2)

      F(25)=E2(2)*E1(3)*C(3)
      F(26)=E2(2)*E1(1)*C(3)
      F(27)=E2(2)*E1(2)*C(3)

      F(28)=E2(3)*E1(1)*C(1)
      F(29)=E2(3)*E1(2)*C(1)
      F(30)=E2(3)*E1(3)*C(1)

      F(31)=E2(3)*E1(2)*C(2)
      F(32)=E2(3)*E1(1)*C(2)
      F(33)=E2(3)*E1(3)*C(2)

      F(34)=E2(3)*E1(3)*C(3)
      F(35)=E2(3)*E1(1)*C(3)
      F(36)=E2(3)*E1(2)*C(3)

      F(37)=0.5D0
      
      GOTO 170

170   CONTINUE
      RETURN
      END

      SUBROUTINE LATSUM(E,FFF,GGG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)  
      DOUBLE PRECISION E,S15ADF,R,R1,R2,V,A,B,C,P
      EXTERNAL S15ADF
      INTEGER IFAIL,L
      INTRINSIC DSQRT, DBLE, DEXP
      DOUBLE PRECISION FFF(50),GGG(50)
      DATA C/1.1283791670955D0/
      V=E*E
      DO L=1,50
      R=DSQRT(DBLE(L))
      R1=E*R
      R2=V*DBLE(L)
      A=1.0D0/DBLE(L)
      IFAIL=1
      P=S15ADF(R1,IFAIL)*A/R
      B=C*E*DEXP(-R2)
      FFF(L)=B*A+P
      GGG(L)=2.0d0*B*V+3.0d0*FFF(L)
      ENDDO
      RETURN
      END

      DOUBLE PRECISION FUNCTION S15ADF(X,IFAIL)
C     MARK 5A REVISED - NAG COPYRIGHT 1976
C     MARK 5C REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     COMPLEMENT OF ERROR FUNCTION ERFC(X)
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X
      INTEGER                          IFAIL
C     .. Local Scalars ..
      DOUBLE PRECISION                 T, XHI, XLO, Y
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, EXP
C     .. Data statements ..
C     PRECISION DEPENDENT CONSTANTS
C08   DATA XLO/-4.5D0/
C12   DATA XLO/-5.25D0/
      DATA XLO/-5.75D0/
C     DATA XLO/-6.25D0/
C18   DATA XLO/-6.5D0/
C
C     RANGE DEPENDENT CONSTANTS
      DATA XHI/ 2.66D+1 /
C     XHI = LARGEST X SUCH THAT EXP(-X*X) .GT. MINREAL (ROUNDED DOWN)
CR1   DATA XHI/13.0D0/
CR2   DATA XHI/9.5D0/
CR3   DATA XHI/13.0D0/
CR4   DATA XHI/25.0D0/
CR5   DATA XHI/26.0D0/
C     .. Executable Statements ..
C
C     NO FAILURE EXITS
      IFAIL = 0
C     TEST EXTREME EXITS
      IF (X.GE.XHI) GO TO 20
      IF (X.LE.XLO) GO TO 40
C
C     EXPANSION ARGUMENT
      T = 1.0D0 - 7.5D0/(ABS(X)+3.75D0)
C
C      * EXPANSION (0021) *
C
C     EXPANSION (0021) EVALUATED AS Y(T)  --PRECISION 08E
C08   Y = (((((((((((+3.1475326D-5)*T-1.3874589D-4)*T-6.4127909D-6)
C08  *    *T+1.7866301D-3)*T-8.2316935D-3)*T+2.4151896D-2)
C08  *    *T-5.4799165D-2)*T+1.0260225D-1)*T-1.6357229D-1)
C08  *    *T+2.2600824D-1)*T-2.7342192D-1)*T + 1.4558972D-1
C
C     EXPANSION (0021) EVALUATED AS Y(T)  --PRECISION 12E
C12   Y = ((((((((((((((((-4.21661579602D-8*T-8.63384346353D-8)
C12  *    *T+6.06038693567D-7)*T+5.90655413508D-7)
C12  *    *T-6.12872971594D-6)*T+3.73223486059D-6)
C12  *    *T+4.78645837248D-5)*T-1.52546487034D-4)
C12  *    *T-2.55222360474D-5)*T+1.80299061562D-3)
C12  *    *T-8.22062412199D-3)*T+2.41432185990D-2)
C12  *    *T-5.48023263289D-2)*T+1.02604312548D-1)
C12  *    *T-1.63571895545D-1)*T+2.26008066898D-1)
C12  *    *T-2.73421931495D-1)*T + 1.45589721275D-1
C
C     EXPANSION (0021) EVALUATED AS Y(T)  --PRECISION 14E
      Y = (((((((((((((((-2.2356173494379D-9
     *    *T+4.5302502889845D-9)*T+2.5918103316137D-8)
     *    *T-6.3684846832869D-8)*T-1.7642194353331D-7)
     *    *T+6.4907607131235D-7)*T+7.4296952017617D-7)
     *    *T-6.1758018478516D-6)*T+3.5866167916231D-6)
     *    *T+4.7895180610590D-5)*T-1.5246364229106D-4)
     *    *T-2.5534256252531D-5)*T+1.8029626230333D-3)
     *    *T-8.2206213481002D-3)*T+2.4143223946968D-2)
     *    *T-5.4802326675661D-2)*T+1.0260431203382D-1
      Y = (((Y*T-1.6357189552481D-1)*T+2.2600806691658D-1)
     *    *T-2.7342193149541D-1)*T + 1.4558972127504D-1
C
C     EXPANSION (0021) EVALUATED AS Y(T)  --PRECISION 16E
C     Y = (((((((((((((((+3.328130055126039D-10
C    *    *T-5.718639670776992D-10)*T-4.066088879757269D-9)
C    *    *T+7.532536116142436D-9)*T+3.026547320064576D-8)
C    *    *T-7.043998994397452D-8)*T-1.822565715362025D-7)
C    *    *T+6.575825478226343D-7)*T+7.478317101785790D-7)
C    *    *T-6.182369348098529D-6)*T+3.584014089915968D-6)
C    *    *T+4.789838226695987D-5)*T-1.524627476123466D-4)
C    *    *T-2.553523453642242D-5)*T+1.802962431316418D-3)
C    *    *T-8.220621168415435D-3)*T+2.414322397093253D-2
C     Y = (((((Y*T-5.480232669380236D-2)*T+1.026043120322792D-1)
C    *    *T-1.635718955239687D-1)*T+2.260080669166197D-1)
C    *    *T-2.734219314954260D-1)*T + 1.455897212750385D-1
C
C     EXPANSION (0021) EVALUATED AS Y(T)  --PRECISION 18E
C18   Y = (((((((((((((((-1.58023488119651697D-11
C18  *    *T-4.94972069009392927D-11)*T+1.86424953544623784D-10)
C18  *    *T+6.29796246918239617D-10)*T-1.34751340973493898D-9)
C18  *    *T-4.84566988844706300D-9)*T+9.22474802259858004D-9)
C18  *    *T+3.14410318645430670D-8)*T-7.26754673242913196D-8)
C18  *    *T-1.83380699508554268D-7)*T+6.59488268069175234D-7)
C18  *    *T+7.48541685740064308D-7)*T-6.18344429012694168D-6)
C18  *    *T+3.58371497984145357D-6)*T+4.78987832434182054D-5)
C18  *    *T-1.52462664665855354D-4)*T-2.55353311432760448D-5
C18   Y = ((((((((Y*T+1.80296241673597993D-3)
C18  *    *T-8.22062115413991215D-3)
C18  *    *T+2.41432239724445769D-2)*T-5.48023266949776152D-2)
C18  *    *T+1.02604312032198239D-1)*T-1.63571895523923969D-1)
C18  *    *T+2.26008066916621431D-1)*T-2.73421931495426482D-1)*T +
C18  *     1.45589721275038539D-1
C
      S15ADF = EXP(-X*X)*Y
      IF (X.LT.0.0D0) S15ADF = 2.0D0 - S15ADF
      RETURN
C
   20 S15ADF = 0.0D0
      RETURN
   40 S15ADF = 2.0D0
      RETURN
C
      END

      DOUBLE PRECISION FUNCTION S15AEF(X,IFAIL)
C     MARK 4 RELEASE NAG COPYRIGHT 1974.
C     MARK 4.5 REVISED
C     MARK 8 REVISED. IER-221 (MAR 1980)
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 14 REVISED. IER-755 (DEC 1989).
C     ERF(X)
C
C     **************************************************************
C
C     TO EXTRACT THE CORRECT CODE FOR A PARTICULAR MACHINE-RANGE,
C     ACTIVATE THE STATEMENTS CONTAINED IN COMMENTS BEGINNING  CDD ,
C     WHERE  DD  IS THE APPROXIMATE NUMBER OF SIGNIFICANT DECIMAL
C     DIGITS REPRESENTED BY THE MACHINE
C     DELETE THE ILLEGAL DUMMY STATEMENT OF THE FORM
C     * EXPANSION (DATA) *
C
C     **************************************************************
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X
      INTEGER                          IFAIL
C     .. Local Scalars ..
      DOUBLE PRECISION                 BJ, BJP1, BJP2, HALF, ONE,
     *                                 SQRTPI, THREE, TWENTY, TWO,
     *                                 X2, XUP, XV, ZERO
      INTEGER                          J, NCFC, NCFD
C     .. Local Arrays ..
C07   DOUBLE PRECISION                 C(8), D(8)
C12   DOUBLE PRECISION                 C(11), D(12)
C14   DOUBLE PRECISION                 C(15), D(15)
      DOUBLE PRECISION                 C(18), D(17)
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, EXP, SIGN
C     .. Data statements ..
C      * EXPANSION (DATA) *
C07   DATA NCFC,NCFD/8,8/,XUP/4.0D0/,SQRTPI/1.772454D0/
C07  A,C(1),C(2),C(3),C(4),C(5),C(6),C(7),C(8)
C07  A/1.944907D0,4.2019D-2,-1.8687D-2,5.129D-3,-1.068D-3
C07  A,1.74D-4,-2.1D-5,2.0D-6/
C07  A,D(1),D(2),D(3),D(4),D(5),D(6),D(7),D(8)
C07  A/1.483110D0,-3.01071D-1,6.8995D-2,-1.3916D-2,2.421D-3
C07  A,-3.66D-4,4.9D-5,-6.0D-6/
C
      DATA NCFC,NCFD/11,12/,XUP/5.0D0/,SQRTPI/1.7724538509D0/,C(1),
     * C(2),C(3),C(4),C(5),C(6),C(7),C(8),C(9),C(10),C(11)
     * /1.9449071068D0,4.20186582D-2,-1.86866104D-2,5.1281062D-3,
     * -1.0683107D-3,1.744738D-4,-2.15642D-5,1.7283D-6,-2.D-8,-1.65D-8,
     * 2.D-9/,D(1),D(2),D(3),D(4),D(5),D(6),D(7),D(8),D(9),
     * D(10),D(11),D(12)/1.4831105641D0,-3.010710734D-1,6.89948307D-2,
     * -1.39162713D-2,2.4207995D-3,-3.658640D-4,4.86210D-5,-5.7493D-6,
     * 6.113D-7,-5.90D-8,5.2D-9,-4.D-10/
C
C14   DATA NCFC,NCFD/15,15/,XUP/5.75D0/,SQRTPI/1.7724538509055D0/
C14  A,C(1),C(2),C(3),C(4),C(5),C(6),C(7),C(8),C(9),C(10)
C14  A,C(11),C(12),C(13),C(14),C(15)
C14  A/1.9449071068179D0,4.20186582324D-2,-1.86866103977D-2
C14  A,5.1281061839D-3,-1.0683107462D-3,1.744737872D-4
C14  A,-2.15642056D-5,1.7282658D-6,-2.00479D-8,-1.64782D-8
C14  A,2.0008D-9,2.58D-11,-3.06D-11,1.9D-12,4.0D-13/
C14  A,D(1),D(2),D(3),D(4),D(5),D(6),D(7),D(8),D(9),D(10)
C14  A,D(11),D(12),D(13),D(14),D(15)
C14  A/1.4831105640848D0,-3.010710733866D-1,6.89948306898D-2
C14  A,-1.39162712647D-2,2.4207995224D-3,-3.658639686D-4
C14  A,4.86209844D-5,-5.7492565D-6,6.113243D-7,-5.89910D-8
C14  A,5.2070D-9,-4.233D-10,3.19D-11,-2.2D-12,1.0D-13/
C
      DATA NCFC,NCFD/18,17/,XUP/6.25D0/,SQRTPI/1.7724538509055160D0/
     A,C(1),C(2),C(3),C(4),C(5),C(6),C(7),C(8),C(9),C(10),C(11)
     A,C(12),C(13),C(14),C(15),C(16),C(17),C(18)
     A/1.9449071068178803D0,4.20186582324414D-2,-1.86866103976769D-2
     A,5.1281061839107D-3,-1.0683107461726D-3,1.744737872522D-4
     A,-2.15642065714D-5,1.7282657974D-6,-2.00479241D-8
     A,-1.64782105D-8,2.0008475D-9,2.57716D-11,-3.06343D-11
     A,1.9158D-12,3.703D-13,-5.43D-14,-4.0D-15,1.2D-15/
     A,D(1),D(2),D(3),D(4),D(5),D(6),D(7),D(8),D(9),D(10),D(11)
     A,D(12),D(13),D(14),D(15),D(16),D(17)
     A/1.4831105640848036D0,-3.010710733865950D-1,6.89948306898316D-2
     A,-1.39162712647222D-2,2.4207995224335D-3,-3.658639685849D-4
     A,4.86209844323D-5,-5.7492565580D-6,6.113243578D-7
     A,-5.89910153D-8,5.2070091D-9,-4.232976D-10,3.18811D-11
     A,-2.2361D-12,1.467D-13,-9.0D-15,5.0D-16/
C
      DATA                             ZERO, ONE, TWO, THREE, TWENTY,
     *                                 HALF/0.0D0, 1.0D0, 2.0D0, 3.0D0,
     *                                 20.0D0, 0.5D0/
C     .. Executable Statements ..
C
C     NO FAILURE EXITS
      IFAIL = 0
      XV = ABS(X)
      IF (XV.GE.XUP) GO TO 120
      IF (XV.LE.TWO) GO TO 60
      X2 = TWO - TWENTY/(XV+THREE)
C
C     SUMMATION
      BJP2 = ZERO
      BJP1 = C(NCFC)
      J = NCFC - 1
   20 BJ = X2*BJP1 - BJP2 + C(J)
      IF (J.EQ.1) GO TO 40
      BJP2 = BJP1
      BJP1 = BJ
      J = J - 1
      GO TO 20
   40 X2 = HALF*(BJ-BJP2)/XV*EXP(-X*X)/SQRTPI
      S15AEF = (ONE-X2)*SIGN(ONE,X)
      GO TO 140
C
   60 X2 = X*X - TWO
C     SUMMATION
      BJP2 = ZERO
      BJP1 = D(NCFD)
      J = NCFD - 1
   80 BJ = X2*BJP1 - BJP2 + D(J)
      IF (J.EQ.1) GO TO 100
      BJP2 = BJP1
      BJP1 = BJ
      J = J - 1
      GO TO 80
  100 S15AEF = HALF*(BJ-BJP2)*X
      GO TO 140
 
  120 S15AEF = SIGN(ONE,X)
  140 RETURN
      END

      SUBROUTINE REPUL(R1,R2,R3,Q)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      INTEGER I,J,K,N,M
      DOUBLE PRECISION TEMP(9),SL,PI,RAM1(3),A12,B12,A11,B11,
     * A22,B22,Z,Y1,Y2,AK1,AK2,F1,F2,G1,G2,H1,H2,Q(3),C(3),
     * S(3),DC(3),CP(3),R1(4,6),R2(4,6),R3(4,6),X01AAF,X
      INTRINSIC DCOS, DSIN
      PI=X01AAF(X)
      OPEN(347,FILE='BBB')
      READ(347,222)
     * RAM1(1),RAM1(2),RAM1(3),A12,B12,A11,B11,A22,B22,Z,Y1,
     * Y2,AK1,AK2,F1,F2,G1,G2,H1,H2,TEMP(1),TEMP(2),TEMP(3),
     * TEMP(4),TEMP(5),TEMP(6),TEMP(7),TEMP(8),TEMP(9)
      CLOSE(347)
      DO K=1,3
      SL=PI*Q(K)
      C(K)=DCOS(SL)
      S(K)=DSIN(SL)
      DC(K)=PI*S(K)
      ENDDO
      P1=2.0d0*B11
      P2=2.0d0*B22
      DO 90 M=1,3
      L=1
      DO 20 I=1,3
      IF (M-I) 15,10,15
10    CP(1)=C(M)
      SL=DC(M)
      GOTO 20
15    L=L+1
      CP(L)=C(I)
20    CONTINUE
      D1=CP(2)+CP(3)
      D2=D1*CP(1)
      D3=CP(2)*CP(3)
      V1=SL*D1
      DO 90 N=M,3
      IF (M-N) 60,25,60
25    R1(1,M)=F1-G1*D2-P1*D3
      R2(1,M)=F2-G2*D2-P2*D3
      R3(1,M)=-A12*CP(1)-B12*D1
      DO 55 K=1,3
      I=K+1
      IF (K-M) 35,30,35
30    R1(I,M)=G1*V1
      R2(I,M)=G2*V1
      R3(I,M)=A12*SL
      GOTO 55
35    DO 50 J=1,3
      IF (J-M) 40,50,40
40    IF (J-K) 45,50,45
45    V2=C(J)
50    CONTINUE
      V3=DC(K)
      R1(I,M)=V3*(G1*CP(1)+P1*V2)
      R2(I,M)=V3*(G2*CP(1)+P2*V2)
      R3(I,M)=B12*V3
55    CONTINUE
      GOTO 90
60    L=M+N+1
      V2=S(M)*S(N)
      R1(1,L)=H1*V1
      R2(1,L)=H2*V2
      R3(1,L)=0.0D0
      DO 86 K=1,3
      I=K+1
      IF (K-M) 70,65,70
65    V3=C(M)*DC(N)
      GOTO 85
70    IF (K-N) 80,75,80
75    V3=C(N)*DC(M)
      GOTO 85
80    V3=0.0d0
85    R1(I,L)=H1*V3
      R2(I,L)=H2*V3
86    R3(I,L)=0.0D0
90    CONTINUE
222   FORMAT(29D30.20)
      RETURN
      END

      SUBROUTINE COUL(C1,C2,Q,E,FFF,GGG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      INTEGER K,I,J,MA
      DOUBLE PRECISION E,C1(4,6),C2(4,6)
      DOUBLE PRECISION FFF(50),GGG(50),PI,PI2,X01AAF,X
      DOUBLE PRECISION Q(3),S(3),G1(6),G2(6),H1(6),H2(6),
     *          DG1(3,6),DG2(3,6),DH1(3,6),DH2(3,6),SQ(3)

      PI=X01AAF(X)
      PI2=PI*PI
      V=E*E
      IF (Q(1)) 5,130,5
5     DO 10 K=1,6
      G1(K)=0.0d0
      G2(K)=0.0d0
      H1(K)=0.0d0
      H2(K)=0.0d0
      DO 10 I=1,3
      DG1(I,K)=0.0d0
      DG2(I,K)=0.0d0
      DH1(I,K)=0.0d0
      DH2(I,K)=0.0d0
10    CONTINUE
      DO 50 MA=1,2
      DO 50 I1=MA,6,2
      I=I1-1
11    S(1)=DBLE(I)
      DO 45 J1=MA,6,2
      J=J1-1
12    S(2)=DBLE(J)
      DO 40 K1=MA,6,2
      K=K1-1
15    S(3)=DBLE(K)
      LR=I**2+J**2+K**2
      IF (LR.GT.35) GOTO 40
      B=(S(1)+Q(1))**2+(S(2)+Q(2))**2+(S(3)+Q(3))**2
      R=0.25D0*B*PI2/V
      IF (R.GT.25.0D0) GOTO 35
      A=DCOS(PI*(S(1)+S(2)+S(3)))
      DO M=1,3
      SQ(M)=S(M)+Q(M)
      ENDDO
      R1=DEXP(-R)/B
      DO 30 M=1,3
      DO 30 N=1,M
      IF (M.EQ.N) GOTO 20
      L=M+N+1
      GOTO 25
20    L=M
25    SL=SQ(M)*SQ(N)
      GA=R1*SL
      GB=GA*A
      G1(L)=G1(L)+GA
      G2(L)=G2(L)+GB
      P=2.0d0*SL*(1.0d0+R)/B
      DO 30 KA=1,3
      IF (KA-M) 22,21,22
21    SL=SQ(N)
      GOTO 23
22    SL=0.0d0
23    IF (KA-N) 26,24,26
24    SL1=SQ(M)
      GOTO 27
26    SL1=0.0d0
27    GA=R1*(SL+SL1-P*SQ(KA))
      GB=GA*A
      DG1(KA,L)=DG1(KA,L)+GA
      DG2(KA,L)=DG2(KA,L)+GB
30    CONTINUE
35    K=-K
      IF (K.LT.0) GOTO 15
40    CONTINUE
      J=-J
      IF (J.LT.0) GOTO 12
45    CONTINUE
      I=-I
      IF (I.LT.0) GOTO 11
50    CONTINUE
      DO 105 I1=1,4
      I=I1-1
56    S(1)=DBLE(I)
      DO 100 J1=1,4
      J=J1-1
57    S(2)=DBLE(J)
      DO 95 K1=1,4
      K=K1-1
60    S(3)=DBLE(K)
      LD=I*I+J*J+K*K
      IF (LD.EQ.0) GOTO 95
      IF (LD.GT.10) GOTO 95
      S(1)=DBLE(I)
      S(2)=DBLE(J)
      S(3)=DBLE(K)
      LA=I+J+K
      IF (LA.EQ.(LA/2)*2) GOTO 65
      KP=2
      GOTO 70
65    KP=1
70    A=1.0d0/DBLE(LD)
      B=PI*(Q(1)*S(1)+Q(2)*S(2)+Q(3)*S(3))
      DH=-PI*DSIN(B)
      H=DCOS(B)
      DO 90 M=1,3
      DO 90 N=1,M
      IF (M.EQ.N) GOTO 75
      L=M+N+1
      FA=0.0d0
      GOTO 80
75    L=M
      FA=-FFF(LD)
80    B=FA+GGG(LD)*S(M)*S(N)*A
      SL=H*B
      SL1=DH*B
      IF (KP.EQ.2) GOTO 85
      H1(L)=H1(L)+SL
      DO 82 KA=1,3
82    DH1(KA,L)=DH1(KA,L)+SL1*S(KA)
      GOTO 90
85    H2(L)=H2(L)+SL
      DO 87 KA=1,3
87    DH2(KA,L)=DH2(KA,L)+SL1*S(KA)
90    CONTINUE
92    K=-K
      IF (K.LT.0) GOTO 60
95    CONTINUE
      J=-J
      IF (J.LT.0) GOTO 57
100   CONTINUE
      I=-I
      IF (I.LT.0) GOTO 56
105   CONTINUE
      A=1.504505561435D0*E*V
      SL=4.0D0*PI
      DO 125 K=1,6
      IF (K.GT.3) GOTO 115
      B=A
      GOTO 120
115   B=0.0D0
120   R=B+2.0D0*H1(K)-SL*G1(K)
      R1=SL*G2(K)-2.0D0*H2(K)
      C1(1,K)=R
      C2(1,K)=R1
      DO 125 I1=1,3
      I=I1+1
      C1(I,K)=2.0d0*DH1(I1,K)-SL*DG1(I1,K)
      C2(I,K)=-2.0d0*DH2(I1,K)+SL*DG2(I1,K)
125   CONTINUE
      GOTO 150
130   A=4.188790D0
      DO 145 K=1,6
      IF (K.GT.3) GOTO 140
      IF (K.EQ.1) GOTO 135
      SL=A
      SL1=-A
      GOTO 144
135   SL=-2.0D0*A
      SL1=-SL
      GOTO 144
140   SL=0.0D0
      SL1=0.0D0
144   C1(1,K)=SL
      C2(1,K)=SL1
      DO 145 I=2,4
      C1(I,K)=0.0D0
      C2(I,K)=0.0D0
145   CONTINUE
150   RETURN
      END

      SUBROUTINE SCOUL(C1,C2,Q,E,F,G)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)  
      DOUBLE PRECISION E,C1(6),C2(6)
      DOUBLE PRECISION F(50),G(50),PI,PI2
      DOUBLE PRECISION Q(3),S(3),G1(6),G2(6),H1(6),H2(6)
      PI=3.14159465358D0
      PI2=PI**2
      MP=0
      MP1=0
      V=E*E
      IF (Q(1)) 7,130,7
7     DO 10 K=1,6
      G1(K)=0.0D0
      G2(K)=0.0D0
      H1(K)=0.0D0
10    H2(K)=0.0D0
      DO 50 MA=1,2
      DO 50 I1=MA,6,2
      I=I1-1
11    S(1)=DBLE(I)
      DO 45 J1=MA,6,2
      J=J1-1
12    S(2)=DBLE(J)
      DO 40 K1=MA,6,2
      K=K1-1
15    S(3)=DBLE(K)
      LR=I**2+J**2+K**2
      IF (LR.GT.35) GOTO 40
      B=(S(1)+Q(1))**2+(S(2)+Q(2))**2+(S(3)+Q(3))**2
      R=0.25D0*B*PI2/V
      IF (R.GT.25.0D0) GOTO 35
      MP=MP+1
      R1=DEXP(-R)/B
      DO 30 M=1,3
      DO 30 N=1,M
      IF (M.EQ.N) GOTO 20
      L=M+N+1
      GOTO 25
20    L=M
25    GA=R1*(S(M)+Q(M))*(S(N)+Q(N))
      GB=GA*DCOS(PI*(S(1)+S(2)+S(3)))
      G1(L)=G1(L)+GA
30    G2(L)=G2(L)+GB
35    K=-K
      IF (K.LT.0) GOTO 15
40    CONTINUE
      J=-J
      IF (J.LT.0) GOTO 12
45    CONTINUE
      I=-I
      IF (I.LT.0) GOTO 11
50    CONTINUE
      DO 105 I1=1,4
      I=I1-1
56    S(1)=DBLE(I)
      DO 100 J1=1,4
      J=J1-1
57    S(2)=DBLE(J)
      DO 95 K1=1,4
      K=K1-1
60    S(3)=DBLE(K)
      LD=I**2+J**2+K**2
      IF (LD.EQ.0) GOTO 95
      IF (LD.GT.10) GOTO 95
      MP1=MP1+1
      S(1)=DBLE(I)
      S(2)=DBLE(J)
      S(3)=DBLE(K)
      LA=I+J+K
      LB=LA/2
      LC=2*LB
      IF (LA.EQ.LC) GOTO 65
      KP=2
      GOTO 70
65    KP=1
70    A=1.0D0/DBLE(LD)
      H=DCOS(PI*(Q(1)*S(1)+Q(2)*S(2)+Q(3)*S(3)))
      DO 90 M=1,3
      DO 90 N=1,M
      IF (M.EQ.N) GOTO 75
      L=M+N+1
      FA=0.0D0
      GOTO 80
75    L=M
      FA=-F(LD)
80    SL=H*(FA+G(LD)*S(M)*S(N)*A)
      IF (KP.EQ.2) GOTO 85
      H1(L)=H1(L)+SL
      GOTO 90
85    H2(L)=H2(L)+SL
90    CONTINUE
92    K=-K
      IF (K.LT.0) GOTO 60
95    CONTINUE
      J=-J
      IF (J.LT.0) GOTO 57
100   CONTINUE
      I=-I
      IF (I.LT.0) GOTO 56
105   CONTINUE
      A=1.504505556D0*E*V
      DO 125 K=1,6
      IF (K.GT.3) GOTO 115
      B=A
      GOTO 120
115   B=0.0D0
120   C1(K)=B+2.0D0*H1(K)-4.0D0*PI*G1(K)
125   C2(K)=4.0D0*PI*G2(K)-2.0D0*H2(K)
      GOTO 150
130   A=4.18879020D0
      DO 145 K=1,6
      IF (K.GT.3) GOTO 140
      IF (K.EQ.1) GOTO 135
      C1(K)=A
      C2(K)=-A
      GOTO 145
135   C1(K)=-2.0D0*A
      C2(K)=2.0D0*A
      GOTO 145
140   C1(K)=0.0D0
      C2(K)=0.0D0
145   CONTINUE
150   RETURN
      END



      SUBROUTINE F01CKF(A,B,C,N,P,M,Z,IZ,OPT,IFAIL)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C
C     RETURNS WITH THE RESULT OF THE MATRIX MULTIPLICATION OF B AND
C     C
C     IN THE MATRIX A, WITH THE OPTION TO OVERWRITE B OR C.
C     **************************************************************
C     ****
C     MATRIX MULTIPLICATION
C     A=B*C          IF OPT= 1 A,B,C ASSUMED DISTINCT
C     IF OPT=2 B IS OVERWRITTEN
C     IF OPT=3 C IS OVERWRITTEN
C     **************************************************************
C     ****
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01CKF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, IZ, M, N, OPT, P
C     .. Array Arguments ..
      DOUBLE PRECISION  A(N,P), B(N,M), C(M,P), Z(IZ)
C     .. Local Scalars ..
      INTEGER           I, IB1, J
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          DGEMM, DGEMV
C     .. Executable Statements ..
      IF (M.GE.1 .AND. N.GE.1 .AND. P.GE.1) GO TO 20
C     *****ERROR CHECK 1 FAILS:- NEGATIVE DIMENSIONS******
      IFAIL = P01ABF(IFAIL,1,SRNAME,0,P01REC)
      RETURN
   20 IF (OPT.EQ.1) GO TO 60
      IF (IZ.GE.M) GO TO 40
      IFAIL = P01ABF(IFAIL,4,SRNAME,0,P01REC)
      RETURN
   40 IF (OPT.EQ.2) GO TO 80
      IF (OPT.EQ.3) GO TO 160
   60 IFAIL = 0
      CALL DGEMM('No transpose','No Transpose',N,P,M,1.0D0,B,N,C,M,
     *           0.0D0,A,N)
      RETURN
C     ******OVERWRITE B******
   80 IF (M.EQ.P) GO TO 100
C     ******ERROR CHECK 2 FAILS:- OVERWRITING INCOMPATIBLE
C     ARRAYS******
      IFAIL = P01ABF(IFAIL,2,SRNAME,0,P01REC)
      RETURN
  100 IFAIL = 0
      IB1 = N
      IF (M.EQ.1) IB1 = 1
      DO 140 I = 1, N
         CALL DGEMV('T',M,P,1.0D0,C,M,B(I,1),IB1,0.0D0,Z,1)
         DO 120 J = 1, M
            B(I,J) = Z(J)
  120    CONTINUE
  140 CONTINUE
      RETURN
C     ******OVERWRITE C******
  160 IF (N.EQ.M) GO TO 180
C     ******ERROR CHECK 3 FAILS:- OVERWRITING INCOMPATIBLE
C     ARRAYS******
      IFAIL = P01ABF(IFAIL,3,SRNAME,0,P01REC)
      RETURN
  180 IFAIL = 0
      DO 220 J = 1, P
         CALL DGEMV('N',M,M,1.0D0,B,N,C(1,J),1,0.0D0,Z,1)
         DO 200 I = 1, M
            C(I,J) = Z(I)
  200    CONTINUE
  220 CONTINUE
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

      SUBROUTINE F06PAF( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     AXP4 VERSION FOR VECTOR MACHINES
C     .. Entry Points ..
      ENTRY      DGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N
      CHARACTER*1        TRANS
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
C     ..
C
C  Purpose
C  =======
C
C  DGEMV  performs one of the matrix-vector operations
C
C     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
C
C  where alpha and beta are scalars, x and y are vectors and A is an
C  m by n matrix.
C
C  Parameters
C  ==========
C
C  TRANS  - CHARACTER*1.
C           On entry, TRANS specifies the operation to be performed as
C           follows:
C
C              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
C
C              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
C
C              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C           On entry, M specifies the number of rows of the matrix A.
C           M must be at least zero.
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the number of columns of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - DOUBLE PRECISION.
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
C           Before entry, the leading m by n part of the array A must
C           contain the matrix of coefficients.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, m ).
C           Unchanged on exit.
C
C  X      - DOUBLE PRECISION array of DIMENSION at least
C           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
C           and at least
C           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
C           Before entry, the incremented array X must contain the
C           vector x.
C           Unchanged on exit.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C  BETA   - DOUBLE PRECISION.
C           On entry, BETA specifies the scalar beta. When BETA is
C           supplied as zero then Y need not be set on input.
C           Unchanged on exit.
C
C  Y      - DOUBLE PRECISION array of DIMENSION at least
C           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
C           and at least
C           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
C           Before entry with BETA non-zero, the incremented array Y
C           must contain the vector y. On exit, Y is overwritten by the
C           updated vector y.
C
C  INCY   - INTEGER.
C           On entry, INCY specifies the increment for the elements of
C           Y. INCY must not be zero.
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C  -- DO-loops unrolled on 20-November-1986.
C     Peter Mayes, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP, TEMP1, TEMP2, TEMP3, TEMP4
      INTEGER            I, INFO, IY, J, JX, KX, KY, LENX, LENY, M4, N4
C     .. External Subroutines ..
      EXTERNAL           F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( .NOT.(TRANS.EQ.'N' .OR. TRANS.EQ.'n').AND.
     $         .NOT.(TRANS.EQ.'T' .OR. TRANS.EQ.'t').AND.
     $         .NOT.(TRANS.EQ.'C' .OR. TRANS.EQ.'c')      )THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL F06AAZ( 'F06PAF/DGEMV ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
C
C     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
C     up the start points in  X  and  Y.
C
      IF( (TRANS.EQ.'N' .OR. TRANS.EQ.'n') )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
C
C     Start the operations. In this version the inner loops are all
C     equivalent to AXPY operations.
C
C     First form  y := beta*y.
C
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      JX = KX
      IF( (TRANS.EQ.'N' .OR. TRANS.EQ.'n') )THEN
C
C        Form  y := alpha*A*x + y.
C
         IF( INCY.EQ.1 )THEN
C**** U n r o l l   t o   d e p t h   4 ********************************
            N4 = 4*( N/4 )
            DO 60, J = 1, N4, 4
               TEMP1 = ALPHA*X( JX )
               TEMP2 = ALPHA*X( JX + INCX )
               TEMP3 = ALPHA*X( JX + 2*INCX )
               TEMP4 = ALPHA*X( JX + 3*INCX )
               IF( TEMP1.NE.ZERO.OR.TEMP2.NE.ZERO.OR.TEMP3.NE.ZERO.OR.
     $             TEMP4.NE.ZERO )THEN
                  DO 50, I = 1, M
                     Y( I ) = ( ( ( ( Y( I ) + TEMP1*A( I, J ) )
     $                        + TEMP2*A( I, J + 1 ) )
     $                        + TEMP3*A( I, J + 2 ) )
     $                        + TEMP4*A( I, J + 3 ) )
   50             CONTINUE
               END IF
               JX = JX + 4*INCX
   60       CONTINUE
C**** Clean-up loop ****************************************************
            DO 80, J = N4 + 1, N, 1
               TEMP = ALPHA*X( JX )
               IF( TEMP.NE.ZERO )THEN
                  DO 70, I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         ELSE
C**** U n r o l l   t o   d e p t h   4 ********************************
            N4 = 4*( N/4 )
            DO 100, J = 1, N4, 4
               TEMP1 = ALPHA*X( JX )
               TEMP2 = ALPHA*X( JX + INCX )
               TEMP3 = ALPHA*X( JX + 2*INCX )
               TEMP4 = ALPHA*X( JX + 3*INCX )
               IF( TEMP1.NE.ZERO.OR.TEMP2.NE.ZERO.OR.TEMP3.NE.ZERO.OR.
     $             TEMP4.NE.ZERO )THEN
                  IY = KY
                  DO 90, I = 1, M
                     Y( IY ) = ( ( ( ( Y( IY ) + TEMP1*A( I, J ) )
     $                         + TEMP2*A( I, J + 1 ) )
     $                         + TEMP3*A( I, J + 2 ) )
     $                         + TEMP4*A( I, J + 3 ) )
                     IY = IY + INCY
   90             CONTINUE
               END IF
               JX = JX + 4*INCX
  100       CONTINUE
C**** Clean-up loop ****************************************************
            DO 120, J = N4 + 1, N, 1
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY = KY
                  DO 110, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY = IY + INCY
  110             CONTINUE
               END IF
               JX = JX + INCX
  120       CONTINUE
         END IF
      ELSE
C
C        Form  y := alpha*A'*x + y.
C
         IF( INCY.EQ.1 )THEN
C**** U n r o l l   t o   d e p t h   4 ********************************
            M4 = 4*( M/4 )
            DO 140, J = 1, M4, 4
               TEMP1 = ALPHA*X( JX )
               TEMP2 = ALPHA*X( JX + INCX )
               TEMP3 = ALPHA*X( JX + 2*INCX )
               TEMP4 = ALPHA*X( JX + 3*INCX )
               IF( TEMP1.NE.ZERO.OR.TEMP2.NE.ZERO.OR.TEMP3.NE.ZERO.OR.
     $             TEMP4.NE.ZERO )THEN
                  DO 130, I = 1, N
                     Y( I ) = ( ( ( ( Y( I ) + TEMP1*A( J, I ) )
     $                        + TEMP2*A( J + 1, I ) )
     $                        + TEMP3*A( J + 2, I ) )
     $                        + TEMP4*A( J + 3, I ) )
  130             CONTINUE
               END IF
               JX = JX + 4*INCX
  140       CONTINUE
C**** Clean-up loop ****************************************************
            DO 160, J = M4 + 1, M, 1
               TEMP = ALPHA*X( JX )
               IF( TEMP.NE.ZERO )THEN
                  DO 150, I = 1, N
                     Y( I ) = Y( I ) + TEMP*A( J, I )
  150             CONTINUE
               END IF
               JX = JX + INCX
  160       CONTINUE
         ELSE
C**** U n r o l l   t o   d e p t h   4 ********************************
            M4 = 4*( M/4 )
            DO 180, J = 1, M4, 4
               TEMP1 = ALPHA*X( JX )
               TEMP2 = ALPHA*X( JX + INCX )
               TEMP3 = ALPHA*X( JX + 2*INCX )
               TEMP4 = ALPHA*X( JX + 3*INCX )
               IF( TEMP1.NE.ZERO.OR.TEMP2.NE.ZERO.OR.TEMP3.NE.ZERO.OR.
     $             TEMP4.NE.ZERO )THEN
                  IY = KY
                  DO 170, I = 1, N
                     Y( IY ) = ( ( ( ( Y( IY ) + TEMP1*A( J, I ) )
     $                         + TEMP2*A( J + 1, I ) )
     $                         + TEMP3*A( J + 2, I ) )
     $                         + TEMP4*A( J + 3, I ) )
                     IY = IY + INCY
  170             CONTINUE
               END IF
               JX = JX + 4*INCX
  180       CONTINUE
C**** Clean-up loop ****************************************************
            DO 200, J = M4 + 1, M, 1
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY = KY
                  DO 190, I = 1, N
                     Y( IY ) = Y( IY ) + TEMP*A( J, I )
                     IY = IY + INCY
  190             CONTINUE
               END IF
               JX = JX + INCX
  200       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06PAF (DGEMV ).
C
      END

      SUBROUTINE F06YAF(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,
     *                  LDC)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C  Purpose
C  =======
C
C  DGEMM  performs one of the matrix-matrix operations
C
C     C := alpha*op( A )*op( B ) + beta*C,
C
C  where  op( X ) is one of
C
C     op( X ) = X   or   op( X ) = X',
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
C              TRANSA = 'C' or 'c',  op( A ) = A'.
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
C              TRANSB = 'C' or 'c',  op( B ) = B'.
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
C  ALPHA  - REAL            .
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  A      - REAL             array of DIMENSION ( LDA, ka ), where ka is
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
C  B      - REAL             array of DIMENSION ( LDB, kb ), where kb is
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
C  BETA   - REAL            .
C           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
C           supplied as zero then C need not be set on input.
C           Unchanged on exit.
C
C  C      - REAL             array of DIMENSION ( LDC, n ).
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
      ENTRY             DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,
     *                  BETA,C,LDC)
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALPHA, BETA
      INTEGER           K, LDA, LDB, LDC, M, N
      CHARACTER*1       TRANSA, TRANSB
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*)
C     .. Local Scalars ..
      DOUBLE PRECISION  TEMP
      INTEGER           I, INFO, J, L, NCOLA, NROWA, NROWB
      LOGICAL           NOTA, NOTB
C     .. External Subroutines ..
      EXTERNAL          F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
C     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
C     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
C     and  columns of  A  and the  number of  rows  of  B  respectively.
C
      NOTA = (TRANSA.EQ.'N' .OR. TRANSA.EQ.'n')
      NOTB = (TRANSB.EQ.'N' .OR. TRANSB.EQ.'n')
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
      IF (( .NOT. NOTA) .AND. ( .NOT. (TRANSA.EQ.'C' .OR. TRANSA.EQ.'c')
     *    ) .AND. ( .NOT. (TRANSA.EQ.'T' .OR. TRANSA.EQ.'t'))) THEN
         INFO = 1
      ELSE IF (( .NOT. NOTB) .AND. ( .NOT. (TRANSB.EQ.'C' .OR.
     *         TRANSB.EQ.'c')) .AND. ( .NOT. (TRANSB.EQ.'T' .OR.
     *         TRANSB.EQ.'t'))) THEN
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
         CALL F06AAZ('F06YAF/DGEMM ',INFO)
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR. (((ALPHA.EQ.ZERO) .OR. (K.EQ.0))
     *     .AND. (BETA.EQ.ONE))) RETURN
C
C     And if  alpha.eq.zero.
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
         ELSE
C
C           Form  C := alpha*A'*B + beta*C
C
            DO 240 J = 1, N
               DO 220 I = 1, M
                  TEMP = ZERO
                  DO 200 L = 1, K
                     TEMP = TEMP + A(L,I)*B(L,J)
  200             CONTINUE
                  IF (BETA.EQ.ZERO) THEN
                     C(I,J) = ALPHA*TEMP
                  ELSE
                     C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                  END IF
  220          CONTINUE
  240       CONTINUE
         END IF
      ELSE
         IF (NOTA) THEN
C
C           Form  C := alpha*A*B' + beta*C
C
            DO 340 J = 1, N
               IF (BETA.EQ.ZERO) THEN
                  DO 260 I = 1, M
                     C(I,J) = ZERO
  260             CONTINUE
               ELSE IF (BETA.NE.ONE) THEN
                  DO 280 I = 1, M
                     C(I,J) = BETA*C(I,J)
  280             CONTINUE
               END IF
               DO 320 L = 1, K
                  IF (B(J,L).NE.ZERO) THEN
                     TEMP = ALPHA*B(J,L)
                     DO 300 I = 1, M
                        C(I,J) = C(I,J) + TEMP*A(I,L)
  300                CONTINUE
                  END IF
  320          CONTINUE
  340       CONTINUE
         ELSE
C
C           Form  C := alpha*A'*B' + beta*C
C
            DO 400 J = 1, N
               DO 380 I = 1, M
                  TEMP = ZERO
                  DO 360 L = 1, K
                     TEMP = TEMP + A(L,I)*B(J,L)
  360             CONTINUE
                  IF (BETA.EQ.ZERO) THEN
                     C(I,J) = ALPHA*TEMP
                  ELSE
                     C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                  END IF
  380          CONTINUE
  400       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06YAF (DGEMM ).
C
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

      SUBROUTINE F01AJF(N,ATOL,A,IA,D,E,Z,IZ)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C
C     TRED2
C     THIS SUBROUTINE REDUCES THE GIVEN LOWER TRIANGLE OF A
C     SYMMETRIC MATRIX, A, STORED IN THE ARRAY A(N,N), TO
C     TRIDIAGONAL FORM USING HOUSEHOLDERS REDUCTION. THE DIAGONAL
C     OF THE RESULT IS STORED IN THE ARRAY D(N) AND THE
C     SUB-DIAGONAL IN THE LAST N - 1 STORES OF THE ARRAY E(N)
C     (WITH THE ADDITIONAL ELEMENT E(1) = 0). THE TRANSFORMATION
C     MATRICES ARE ACCUMULATED IN THE ARRAY Z(N,N). THE ARRAY
C     A IS LEFT UNALTERED UNLESS THE ACTUAL PARAMETERS
C     CORRESPONDING TO A AND Z ARE IDENTICAL.
C     1ST AUGUST 1971
C
C     REVISED BY VINCE FERNANDO AT MARK 14 TO INTRODUCE SCALING INTO
C     THE GENERATION OF HOUSEHOLDER MATRICES AS PROPOSED BY
C     G.W. STEWART, INTRODUCTION TO MATRIX COMPUTATIONS, CHAPTER 7.
C     ATOL IS NOW A DUMMY PARAMETER.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ATOL
      INTEGER           IA, IZ, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N), D(N), E(N), Z(IZ,N)
C     .. Local Scalars ..
      DOUBLE PRECISION  F, G, H, HH, SCALE
      INTEGER           I, II, J, K, L
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      INTEGER           IDAMAX
      EXTERNAL          DDOT, IDAMAX
C     .. External Subroutines ..
      EXTERNAL          DGEMV, DGER, DSCAL, DSYMV, DSYR2
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Executable Statements ..
      DO 40 I = 1, N
         DO 20 J = I, N
            Z(J,I) = A(J,I)
   20    CONTINUE
         D(I) = A(N,I)
   40 CONTINUE
      IF (N.EQ.1) GO TO 340
      DO 240 II = 2, N
         I = N - II + 2
         L = I - 1
         IF (L.EQ.1) GO TO 60
C        FIND THE ELEMENT OF LARGEST ABSOLUTE VALUE IN D
         K = IDAMAX(L,D,1)
         SCALE = ABS(D(K))
C        IF D IS A NULL VECTOR THEN SKIP THE TRANSFORMATION
         IF (SCALE.GT.0.0D0) GO TO 120
   60    E(I) = D(L)
         H = 0.0D0
         DO 80 J = 1, L
            Z(J,I) = 0.0D0
   80    CONTINUE
         DO 100 J = 1, L
            Z(I,J) = 0.0D0
            D(J) = Z(I-1,J)
  100    CONTINUE
         GO TO 220
  120    CALL DSCAL(L,1.0D0/SCALE,D,1)
         H = DDOT(L,D,1,D,1)
         F = D(I-1)
         G = SQRT(H)
         IF (F.GE.0.0D0) G = -G
         E(I) = G*SCALE
         H = H - F*G
         D(I-1) = F - G
C        COPY U
         DO 140 J = 1, L
            Z(J,I) = D(J)
  140    CONTINUE
C        FORM A*U
         CALL DSYMV('L',L,1.0D0/H,Z,IZ,D,1,0.0D0,E,1)
C        FORM P
         F = 0.0D0
         DO 160 J = 1, L
            F = F + E(J)*D(J)
  160    CONTINUE
C        FORM K
         HH = F/(H+H)
C        FORM Q
         DO 180 J = 1, L
            E(J) = E(J) - HH*D(J)
  180    CONTINUE
C        FORM REDUCED A
         CALL DSYR2('L',L,-1.0D0,D,1,E,1,Z,IZ)
         DO 200 J = 1, L
            D(J) = Z(L,J)
            Z(I,J) = 0.0D0
  200    CONTINUE
  220    D(I) = H
  240 CONTINUE
C     ACCUMULATION OF TRANSFORMATION MATRICES
      DO 300 I = 2, N
         L = I - 1
         Z(N,L) = Z(L,L)
         Z(L,L) = 1.0D0
         H = D(I)
         IF (H.EQ.0.0D0) GO TO 260
         CALL DGEMV('T',L,L,1.0D0/H,Z,IZ,Z(1,I),1,0.0D0,D,1)
         CALL DGER(L,L,-1.0D0,Z(1,I),1,D,1,Z,IZ)
  260    DO 280 J = 1, L
            Z(J,I) = 0.0D0
  280    CONTINUE
  300 CONTINUE
      DO 320 I = 1, N
         D(I) = Z(N,I)
         Z(N,I) = 0.0D0
  320 CONTINUE
  340 Z(N,N) = 1.0D0
      E(1) = 0.0D0
      RETURN
      END

      SUBROUTINE F02AMF(N,EPS,D,E,Z,IZ,IFAIL)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C
C     TQL2
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS OF A
C     TRIDIAGONAL MATRIX, T, GIVEN WITH ITS DIAGONAL ELEMENTS IN
C     THE ARRAY D(N) AND ITS SUB-DIAGONAL ELEMENTS IN THE LAST N
C     - 1 STORES OF THE ARRAY E(N), USING QL TRANSFORMATIONS. THE
C     EIGENVALUES ARE OVERWRITTEN ON THE DIAGONAL ELEMENTS IN THE
C     ARRAY D IN ASCENDING ORDER. THE EIGENVECTORS ARE FORMED IN
C     THE ARRAY Z(N,N), OVERWRITING THE ACCUMULATED
C     TRANSFORMATIONS AS SUPPLIED BY THE SUBROUTINE F01AJF. THE
C     SUBROUTINE WILL FAIL IF ALL EIGENVALUES TAKE MORE THAN 30*N
C     ITERATIONS.
C     1ST APRIL 1972
C
C     .. Parameters ..
      INTEGER           VLEN
      PARAMETER         (VLEN=128)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02AMF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EPS
      INTEGER           IFAIL, IZ, N
C     .. Array Arguments ..
      DOUBLE PRECISION  D(N), E(N), Z(IZ,N)
C     .. Local Scalars ..
      DOUBLE PRECISION  B, C, F, G, H, P, R, S
      INTEGER           I, I1, IPOS, ISAVE, ISEG, J, K, L, M
C     .. Local Arrays ..
      DOUBLE PRECISION  CC(VLEN), SS(VLEN)
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          F06QXF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, SQRT
C     .. Executable Statements ..
      ISAVE = IFAIL
      IF (N.EQ.1) GO TO 40
      DO 20 I = 2, N
         E(I-1) = E(I)
   20 CONTINUE
   40 E(N) = 0.0D0
      B = 0.0D0
      F = 0.0D0
      J = 30*N
      DO 300 L = 1, N
         H = EPS*(ABS(D(L))+ABS(E(L)))
         IF (B.LT.H) B = H
C        LOOK FOR SMALL SUB-DIAG ELEMENT
         DO 60 M = L, N
            IF (ABS(E(M)).LE.B) GO TO 80
   60    CONTINUE
   80    IF (M.EQ.L) GO TO 280
  100    IF (J.LE.0) GO TO 400
         J = J - 1
C        FORM SHIFT
         G = D(L)
         H = D(L+1) - G
         IF (ABS(H).GE.ABS(E(L))) GO TO 120
         P = H*0.5D0/E(L)
         R = SQRT(P*P+1.0D0)
         H = P + R
         IF (P.LT.0.0D0) H = P - R
         D(L) = E(L)/H
         GO TO 140
  120    P = 2.0D0*E(L)/H
         R = SQRT(P*P+1.0D0)
         D(L) = E(L)*P/(1.0D0+R)
  140    H = G - D(L)
         I1 = L + 1
         IF (I1.GT.N) GO TO 180
         DO 160 I = I1, N
            D(I) = D(I) - H
  160    CONTINUE
  180    F = F + H
C        QL TRANSFORMATION
         P = D(M)
         C = 1.0D0
         S = 0.0D0
         DO 260 K = M - 1, L, -VLEN
            ISEG = MAX(K-VLEN+1,L)
            DO 240 I = K, ISEG, -1
               G = C*E(I)
               H = C*P
               IF (ABS(P).LT.ABS(E(I))) GO TO 200
               C = E(I)/P
               R = SQRT(C*C+1.0D0)
               E(I+1) = S*P*R
               S = C/R
               C = 1.0D0/R
               GO TO 220
  200          C = P/E(I)
               R = SQRT(C*C+1.0D0)
               E(I+1) = S*E(I)*R
               S = 1.0D0/R
               C = C/R
  220          P = C*D(I) - S*G
               D(I+1) = H + S*(C*G+S*D(I))
C           STORE ROTATIONS
               CC(VLEN-K+I) = C
               SS(VLEN-K+I) = -S
  240       CONTINUE
C        UPDATE VECTORS
            IPOS = VLEN - K + ISEG
            CALL F06QXF('Right','Variable','Backward',N,K-ISEG+2,1,
     *                  K-ISEG+2,CC(IPOS),SS(IPOS),Z(1,ISEG),IZ)
  260    CONTINUE
         E(L) = S*P
         D(L) = C*P
         IF (ABS(E(L)).GT.B) GO TO 100
  280    D(L) = D(L) + F
  300 CONTINUE
C     ORDER EIGENVALUES AND EIGENVECTORS
      DO 380 I = 1, N
         K = I
         P = D(I)
         I1 = I + 1
         IF (I1.GT.N) GO TO 340
         DO 320 J = I1, N
            IF (D(J).GE.P) GO TO 320
            K = J
            P = D(J)
  320    CONTINUE
  340    IF (K.EQ.I) GO TO 380
         D(K) = D(I)
         D(I) = P
         DO 360 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  360    CONTINUE
  380 CONTINUE
      IFAIL = 0
      RETURN
  400 IFAIL = P01ABF(ISAVE,1,SRNAME,0,P01REC)
      RETURN
      END

      DOUBLE PRECISION FUNCTION F06EAF( N, X, INCX, Y, INCY )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      DOUBLE PRECISION          DDOT
      ENTRY                     DDOT  ( N, X, INCX, Y, INCY )
C     .. Scalar Arguments ..
      INTEGER                           INCX, INCY, N
C     .. Array Arguments ..
      DOUBLE PRECISION                  X( * ), Y( * )
C     ..
C
C  F06EAF returns the value
C
C     F06EAF = x'y
C
C
C  Nag Fortran 77 version of the Blas routine DDOT.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 21-September-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION      ZERO
      PARAMETER           ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION      SUM
      INTEGER               I, IX, IY
C     ..
C     .. Executable Statements ..
      SUM = ZERO
      IF( N.GT.0 )THEN
         IF( ( INCX.EQ.INCY ).AND.( INCX.GT.0 ) )THEN
            DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
               SUM = SUM + X( IX )*Y( IX )
   10       CONTINUE
         ELSE
            IF( INCY.GE.0 )THEN
               IY = 1
            ELSE
               IY = 1 - ( N - 1 )*INCY
            END IF
            IF( INCX.GT.0 )THEN
               DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
                  SUM = SUM + X( IX )*Y( IY )
                  IY  = IY  + INCY
   20          CONTINUE
            ELSE
               IX = 1 - ( N - 1 )*INCX
               DO 30, I = 1, N
                  SUM = SUM + X( IX )*Y( IY )
                  IX  = IX  + INCX
                  IY  = IY  + INCY
   30          CONTINUE
            END IF
         END IF
      END IF
C
      F06EAF = SUM
      RETURN
C
C     End of F06EAF. ( DDOT )
C
      END
      SUBROUTINE F06EDF( N, ALPHA, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      ENTRY      DSCAL ( N, ALPHA, X, INCX )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
C     ..
C
C  F06EDF performs the operation
C
C     x := alpha*x
C
C
C  Nag Fortran 77 version of the Blas routine DSCAL.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 26-November-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      INTEGER            IX
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ALPHA.EQ.ZERO )THEN
            DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = ZERO
   10       CONTINUE
         ELSE IF( ALPHA.EQ.( -ONE ) )THEN
            DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = -X( IX )
   20       CONTINUE
         ELSE IF( ALPHA.NE.ONE )THEN
            DO 30, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = ALPHA*X( IX )
   30       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06EDF. ( DSCAL )
C
      END
      INTEGER FUNCTION F06JLF( N, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      INTEGER          IDAMAX
      ENTRY            IDAMAX( N, X, INCX )
C     .. Scalar Arguments ..
      INTEGER                  INCX, N
C     .. Array Arguments ..
      DOUBLE PRECISION         X( * )
C     ..
C
C  F06JLF returns the smallest value of i such that
C
C     abs( x( i ) ) = max( abs( x( j ) ) )
C                      j
C
C
C  Nag Fortran 77 version of the Blas routine IDAMAX.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 31-May-1983.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      DOUBLE PRECISION         XMAX
      INTEGER                  I, IMAX, IX
C     .. Intrinsic Functions ..
      INTRINSIC                ABS
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IMAX = 1
         IF( N.GT.1 )THEN
            XMAX = ABS( X( 1 ) )
            IX   = 1
            DO 10, I = 2, N
               IX = IX + INCX
               IF( XMAX.LT.ABS( X( IX ) ) )THEN
                  XMAX = ABS( X( IX ) )
                  IMAX = I
               END IF
   10       CONTINUE
         END IF
      ELSE
         IMAX = 0
      END IF
C
      F06JLF = IMAX
      RETURN
C
C     End of F06JLF. ( IDAMAX )
C
      END


      SUBROUTINE F06PCF( UPLO, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     DOT VERSION FOR VECTOR MACHINES
C     .. Entry Points ..
      ENTRY      DSYMV ( UPLO, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA, BETA
      INTEGER            INCX, INCY, LDA, N
      CHARACTER*1        UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
C     ..
C
C  Purpose
C  =======
C
C  DSYMV  performs the matrix-vector  operation
C
C     y := alpha*A*x + beta*y,
C
C  where alpha and beta are scalars, x and y are n element vectors and
C  A is an n by n symmetric matrix.
C
C  Parameters
C  ==========
C
C  UPLO   - CHARACTER*1.
C           On entry, UPLO specifies whether the upper or lower
C           triangular part of the array A is to be referenced as
C           follows:
C
C              UPLO = 'U' or 'u'   Only the upper triangular part of A
C                                  is to be referenced.
C
C              UPLO = 'L' or 'l'   Only the lower triangular part of A
C                                  is to be referenced.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the order of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - DOUBLE PRECISION.
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
C           Before entry with  UPLO = 'U' or 'u', the leading n by n
C           upper triangular part of the array A must contain the upper
C           triangular part of the symmetric matrix and the strictly
C           lower triangular part of A is not referenced.
C           Before entry with UPLO = 'L' or 'l', the leading n by n
C           lower triangular part of the array A must contain the lower
C           triangular part of the symmetric matrix and the strictly
C           upper triangular part of A is not referenced.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, n ).
C           Unchanged on exit.
C
C  X      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the n
C           element vector x.
C           Unchanged on exit.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C  BETA   - DOUBLE PRECISION.
C           On entry, BETA specifies the scalar beta. When BETA is
C           supplied as zero then Y need not be set on input.
C           Unchanged on exit.
C
C  Y      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCY ) ).
C           Before entry, the incremented array Y must contain the n
C           element vector y. On exit, Y is overwritten by the updated
C           vector y.
C
C  INCY   - INTEGER.
C           On entry, INCY specifies the increment for the elements of
C           Y. INCY must not be zero.
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY
C     .. External Subroutines ..
      EXTERNAL           F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( .NOT.(UPLO.EQ.'U' .OR. UPLO.EQ.'u').AND.
     $         .NOT.(UPLO.EQ.'L' .OR. UPLO.EQ.'l')      )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 5
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 7
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 10
      END IF
      IF( INFO.NE.0 )THEN
         CALL F06AAZ( 'F06PCF/DSYMV ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( ( N.EQ.0 ).OR.( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
C
C     Set up the start points in  X  and  Y.
C
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( N - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( N - 1 )*INCY
      END IF
C
C     Start the operations. In this version the inner loops are all
C     equivalent to dot-product operations.
C
C     First form  y := beta*y.
C
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, N
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, N
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, N
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, N
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      JY = KY
      IF( (UPLO.EQ.'U' .OR. UPLO.EQ.'u') )THEN
C
C        Form  y  when A is stored in upper triangle.
C
         IF( INCX.EQ.1 )THEN
            DO 60, J = 1, N
               TEMP = ZERO
               DO 50, I = 1, J
                  TEMP = TEMP + A( I, J )*X( I )
   50          CONTINUE
               DO 55, I = J + 1, N
                  TEMP = TEMP + A( J, I )*X( I )
   55          CONTINUE
               Y(JY) = Y( JY ) + ALPHA*TEMP
               JY    = JY      + INCY
   60       CONTINUE
         ELSE
            JX = KX
            DO 80, J = 1, N
               TEMP = ZERO
               IX   = KX
               DO 70, I = 1, J
                  TEMP = TEMP + A( I, J )*X( IX )
                  IX   = IX   + INCX
   70          CONTINUE
               IX = JX + INCX
               DO 75, I = J + 1, N
                  TEMP = TEMP + A( J, I )*X( IX )
                  IX   = IX   + INCX
   75          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
               JX      = JX      + INCX
   80       CONTINUE
         END IF
      ELSE
C
C        Form  y  when A is stored in lower triangle.
C
         IF( INCX.EQ.1 )THEN
            DO 100, J = 1, N
               TEMP = ZERO
               DO 85, I = 1, J - 1
                  TEMP = TEMP + A( J, I )*X( I )
   85          CONTINUE
               DO 90, I = J, N
                  TEMP = TEMP + A( I, J )*X( I )
   90          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  100       CONTINUE
         ELSE
            JX = KX
            DO 120, J = 1, N
               TEMP = ZERO
               IX   = KX
               DO 105, I = 1, J - 1
                  TEMP = TEMP + A( J, I )*X( IX )
                  IX   = IX   + INCX
  105          CONTINUE
               IX = JX
               DO 110, I = J, N
                  TEMP = TEMP + A( I, J )*X( IX )
                  IX   = IX   + INCX
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
               JX      = JX      + INCX
  120       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06PCF (DSYMV ).
C
      END
      SUBROUTINE F06PMF( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     .. Entry Points ..
      ENTRY      DGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, INCY, LDA, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
C     ..
C
C  Purpose
C  =======
C
C  DGER   performs the rank 1 operation
C
C     A := alpha*x*y' + A,
C
C  where alpha is a scalar, x is an m element vector, y is an n element
C  vector and A is an m by n matrix.
C
C  Parameters
C  ==========
C
C  M      - INTEGER.
C           On entry, M specifies the number of rows of the matrix A.
C           M must be at least zero.
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the number of columns of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - DOUBLE PRECISION.
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  X      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( m - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the m
C           element vector x.
C           Unchanged on exit.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C  Y      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCY ) ).
C           Before entry, the incremented array Y must contain the n
C           element vector y.
C           Unchanged on exit.
C
C  INCY   - INTEGER.
C           On entry, INCY specifies the increment for the elements of
C           Y. INCY must not be zero.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
C           Before entry, the leading m by n part of the array A must
C           contain the matrix of coefficients. On exit, A is
C           overwritten by the updated matrix.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, m ).
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JY, KX
C     .. External Subroutines ..
      EXTERNAL           F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( M.LT.0 )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL F06AAZ( 'F06PMF/DGER  ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through A.
C
      IF( INCY.GT.0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.EQ.1 )THEN
         DO 20, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
            JY = JY + INCY
   20    CONTINUE
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO 40, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
C
      RETURN
C
C     End of F06PMF (DGER  ).
C
      END
      SUBROUTINE F06PRF( UPLO, N, ALPHA, X, INCX, Y, INCY, A, LDA )
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     .. Entry Points ..
      ENTRY      DSYR2 ( UPLO, N, ALPHA, X, INCX, Y, INCY, A, LDA )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, INCY, LDA, N
      CHARACTER*1        UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
C     ..
C
C  Purpose
C  =======
C
C  DSYR2  performs the symmetric rank 2 operation
C
C     A := alpha*x*y' + alpha*y*x' + A,
C
C  where alpha is a scalar, x and y are n element vectors and A is an n
C  by n symmetric matrix.
C
C  Parameters
C  ==========
C
C  UPLO   - CHARACTER*1.
C           On entry, UPLO specifies whether the upper or lower
C           triangular part of the array A is to be referenced as
C           follows:
C
C              UPLO = 'U' or 'u'   Only the upper triangular part of A
C                                  is to be referenced.
C
C              UPLO = 'L' or 'l'   Only the lower triangular part of A
C                                  is to be referenced.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the order of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - DOUBLE PRECISION.
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  X      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the n
C           element vector x.
C           Unchanged on exit.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C  Y      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCY ) ).
C           Before entry, the incremented array Y must contain the n
C           element vector y.
C           Unchanged on exit.
C
C  INCY   - INTEGER.
C           On entry, INCY specifies the increment for the elements of
C           Y. INCY must not be zero.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
C           Before entry with  UPLO = 'U' or 'u', the leading n by n
C           upper triangular part of the array A must contain the upper
C           triangular part of the symmetric matrix and the strictly
C           lower triangular part of A is not referenced. On exit, the
C           upper triangular part of the array A is overwritten by the
C           upper triangular part of the updated matrix.
C           Before entry with UPLO = 'L' or 'l', the leading n by n
C           lower triangular part of the array A must contain the lower
C           triangular part of the symmetric matrix and the strictly
C           upper triangular part of A is not referenced. On exit, the
C           lower triangular part of the array A is overwritten by the
C           lower triangular part of the updated matrix.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, n ).
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP1, TEMP2
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY
C     .. External Subroutines ..
      EXTERNAL           F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( .NOT.(UPLO.EQ.'U' .OR. UPLO.EQ.'u').AND.
     $         .NOT.(UPLO.EQ.'L' .OR. UPLO.EQ.'l')      )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL F06AAZ( 'F06PRF/DSYR2 ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( ( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
C
C     Set up the start points in X and Y if the increments are not both
C     unity.
C
      IF( ( INCX.NE.1 ).OR.( INCY.NE.1 ) )THEN
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( N - 1 )*INCX
         END IF
         IF( INCY.GT.0 )THEN
            KY = 1
         ELSE
            KY = 1 - ( N - 1 )*INCY
         END IF
         JX = KX
         JY = KY
      END IF
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through the triangular part
C     of A.
C
      IF( (UPLO.EQ.'U' .OR. UPLO.EQ.'u') )THEN
C
C        Form  A  when A is stored in the upper triangle.
C
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 20, J = 1, N
               IF( ( X( J ).NE.ZERO ).OR.( Y( J ).NE.ZERO ) )THEN
                  TEMP1 = ALPHA*Y( J )
                  TEMP2 = ALPHA*X( J )
                  DO 10, I = 1, J
                     A( I, J ) = A( I, J ) + X( I )*TEMP1 + Y( I )*TEMP2
   10             CONTINUE
               END IF
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               IF( ( X( JX ).NE.ZERO ).OR.( Y( JY ).NE.ZERO ) )THEN
                  TEMP1 = ALPHA*Y( JY )
                  TEMP2 = ALPHA*X( JX )
                  IX    = KX
                  IY    = KY
                  DO 30, I = 1, J
                     A( I, J ) = A( I, J ) + X( IX )*TEMP1
     $                                     + Y( IY )*TEMP2
                     IX        = IX        + INCX
                     IY        = IY        + INCY
   30             CONTINUE
               END IF
               JX = JX + INCX
               JY = JY + INCY
   40       CONTINUE
         END IF
      ELSE
C
C        Form  A  when A is stored in the lower triangle.
C
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 60, J = 1, N
               IF( ( X( J ).NE.ZERO ).OR.( Y( J ).NE.ZERO ) )THEN
                  TEMP1 = ALPHA*Y( J )
                  TEMP2 = ALPHA*X( J )
                  DO 50, I = J, N
                     A( I, J ) = A( I, J ) + X( I )*TEMP1 + Y( I )*TEMP2
   50             CONTINUE
               END IF
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( ( X( JX ).NE.ZERO ).OR.( Y( JY ).NE.ZERO ) )THEN
                  TEMP1 = ALPHA*Y( JY )
                  TEMP2 = ALPHA*X( JX )
                  IX    = JX
                  IY    = JY
                  DO 70, I = J, N
                     A( I, J ) = A( I, J ) + X( IX )*TEMP1
     $                                     + Y( IY )*TEMP2
                     IX        = IX        + INCX
                     IY        = IY        + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
               JY = JY + INCY
   80       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06PRF (DSYR2 ).
C
      END
      SUBROUTINE F06QXF( SIDE, PIVOT, DIRECT, M, N, K1, K2, C, S, A,
     $                   LDA )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER            K1, K2, LDA, M, N
      CHARACTER*1        DIRECT, PIVOT, SIDE
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
C     ..
C
C  F06QXF  performs the transformation
C
C     A := P*A,   when   SIDE = 'L' or 'l'  (  Left-hand side )
C
C     A := A*P',  when   SIDE = 'R' or 'r'  ( Right-hand side )
C
C  where A is an m by n matrix and P is an orthogonal matrix, consisting
C  of a  sequence  of  plane  rotations,  applied  in  planes  k1 to k2,
C  determined by the parameters PIVOT and DIRECT as follows:
C
C     When  PIVOT  = 'V' or 'v'  ( Variable pivot )
C     and   DIRECT = 'F' or 'f'  ( Forward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
C
C        where  P( k )  is a plane rotation matrix for the  ( k, k + 1 )
C        plane.
C
C     When  PIVOT  = 'V' or 'v'  ( Variable pivot )
C     and   DIRECT = 'B' or 'b'  ( Backward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C        where  P( k )  is a plane rotation matrix for the  ( k, k + 1 )
C        plane.
C
C     When  PIVOT  = 'T' or 't'  ( Top pivot )
C     and   DIRECT = 'F' or 'f'  ( Forward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k2 - 1 )*P( k2 - 2 )*...*P( k1 ),
C
C        where  P( k )  is a plane rotation matrix for the ( k1, k + 1 )
C        plane.
C
C     When  PIVOT  = 'T' or 't'  ( Top pivot )
C     and   DIRECT = 'B' or 'b'  ( Backward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C        where  P( k )  is a plane rotation matrix for the ( k1, k + 1 )
C        plane.
C
C     When  PIVOT  = 'B' or 'b'  ( Bottom pivot )
C     and   DIRECT = 'F' or 'f'  ( Forward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k2 - 1 )*P( k2 - 2 )*...*P( k1 ),
C
C        where  P( k )  is a  plane rotation  matrix  for the  ( k, k2 )
C        plane.
C
C     When  PIVOT  = 'B' or 'b'  ( Bottom pivot )
C     and   DIRECT = 'B' or 'b'  ( Backward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C        where  P( k )  is a  plane rotation  matrix  for the  ( k, k2 )
C        plane.
C
C  c( k ) and s( k )  must contain the  cosine and sine  that define the
C  matrix  P( k ).  The  two by two  plane rotation  part of the  matrix
C  P( k ), R( k ), is assumed to be of the form
C
C     R( k ) = (  c( k )  s( k ) ).
C              ( -s( k )  c( k ) )
C
C  If m, n or k1 are less than unity,  or k2 is not greater than k1,  or
C  SIDE = 'L' or 'l'  and  k2  is greater than  m, or  SIDE = 'R' or 'r'
C  and  k2  is greater than  n,  then an  immediate return  is effected.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 20-November-1986.
C     Sven Hammarling and Mick Pont, Nag Central Office.
C
C     Alternative code used in this version when SIDE.eq.'L' to
C     enable vectorisation but working with rows of the array A.
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   CTEMP, STEMP, TEMP
      INTEGER            I, J
      LOGICAL            LEFT, RIGHT
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      LEFT = ( SIDE.EQ.'L' ).OR.( SIDE.EQ.'l' )
      RIGHT = ( SIDE.EQ.'R' ).OR.( SIDE.EQ.'r' )
      IF( ( MIN( M, N, K1 ).LT.1 ).OR.( K2.LE.K1 ).OR.
     $    ( ( LEFT ).AND.( K2.GT.M ) ).OR.
     $    ( ( RIGHT ).AND.( K2.GT.N ) ) )RETURN
      IF( LEFT )THEN
         IF( ( PIVOT.EQ.'V' ).OR.( PIVOT.EQ.'v' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 20 J = K1, K2 - 1
                  IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
                     CTEMP = C( J )
                     STEMP = S( J )
                     DO 10 I = 1, N
                        TEMP = A( J + 1, I )
                        A( J + 1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
   10                CONTINUE
                  END IF
   20          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 40 J = K2 - 1, K1, -1
                  IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
                     CTEMP = C( J )
                     STEMP = S( J )
                     DO 30 I = 1, N
                        TEMP = A( J + 1, I )
                        A( J + 1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
   30                CONTINUE
                  END IF
   40          CONTINUE
            END IF
         ELSE IF( ( PIVOT.EQ.'T' ).OR.( PIVOT.EQ.'t' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 60 J = K1 + 1, K2
                  CTEMP = C( J - 1 )
                  STEMP = S( J - 1 )
                  IF( ( CTEMP.NE.ONE ).OR.( STEMP.NE.ZERO ) )THEN
                     DO 50 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( K1, I )
                        A( K1, I ) = STEMP*TEMP + CTEMP*A( K1, I )
   50                CONTINUE
                  END IF
   60          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 80 J = K2, K1 + 1, -1
                  CTEMP = C( J - 1 )
                  STEMP = S( J - 1 )
                  IF( ( CTEMP.NE.ONE ).OR.( STEMP.NE.ZERO ) )THEN
                     DO 70 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( K1, I )
                        A( K1, I ) = STEMP*TEMP + CTEMP*A( K1, I )
   70                CONTINUE
                  END IF
   80          CONTINUE
            END IF
         ELSE IF( ( PIVOT.EQ.'B' ).OR.( PIVOT.EQ.'b' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 100 J = K1, K2 - 1
                  IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
                     CTEMP = C( J )
                     STEMP = S( J )
                     DO 90 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( K2, I ) + CTEMP*TEMP
                        A( K2, I ) = CTEMP*A( K2, I ) - STEMP*TEMP
   90                CONTINUE
                  END IF
  100          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 120 J = K2 - 1, K1, -1
                  IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
                     CTEMP = C( J )
                     STEMP = S( J )
                     DO 110 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( K2, I ) + CTEMP*TEMP
                        A( K2, I ) = CTEMP*A( K2, I ) - STEMP*TEMP
  110                CONTINUE
                  END IF
  120          CONTINUE
            END IF
         END IF
      ELSE IF( RIGHT )THEN
         IF( ( PIVOT.EQ.'V' ).OR.( PIVOT.EQ.'v' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 140 J = K1, K2 - 1
                  IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
                     CTEMP = C( J )
                     STEMP = S( J )
                     DO 130 I = 1, M
                        TEMP = A( I, J + 1 )
                        A( I, J + 1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  130                CONTINUE
                  END IF
  140          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 160 J = K2 - 1, K1, -1
                  IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
                     CTEMP = C( J )
                     STEMP = S( J )
                     DO 150 I = 1, M
                        TEMP = A( I, J + 1 )
                        A( I, J + 1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  150                CONTINUE
                  END IF
  160          CONTINUE
            END IF
         ELSE IF( ( PIVOT.EQ.'T' ).OR.( PIVOT.EQ.'t' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 180 J = K1 + 1, K2
                  CTEMP = C( J - 1 )
                  STEMP = S( J - 1 )
                  IF( ( CTEMP.NE.ONE ).OR.( STEMP.NE.ZERO ) )THEN
                     DO 170 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, K1 )
                        A( I, K1 ) = STEMP*TEMP + CTEMP*A( I, K1 )
  170                CONTINUE
                  END IF
  180          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 200 J = K2, K1 + 1, -1
                  CTEMP = C( J - 1 )
                  STEMP = S( J - 1 )
                  IF( ( CTEMP.NE.ONE ).OR.( STEMP.NE.ZERO ) )THEN
                     DO 190 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, K1 )
                        A( I, K1 ) = STEMP*TEMP + CTEMP*A( I, K1 )
  190                CONTINUE
                  END IF
  200          CONTINUE
            END IF
         ELSE IF( ( PIVOT.EQ.'B' ).OR.( PIVOT.EQ.'b' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 220 J = K1, K2 - 1
                  IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
                     CTEMP = C( J )
                     STEMP = S( J )
                     DO 210 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, K2 ) + CTEMP*TEMP
                        A( I, K2 ) = CTEMP*A( I, K2 ) - STEMP*TEMP
  210                CONTINUE
                  END IF
  220          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 240 J = K2 - 1, K1, -1
                  IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
                     CTEMP = C( J )
                     STEMP = S( J )
                     DO 230 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, K2 ) + CTEMP*TEMP
                        A( I, K2 ) = CTEMP*A( I, K2 ) - STEMP*TEMP
  230                CONTINUE
                  END IF
  240          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06QXF. ( SGESRC )
C
      END

      SUBROUTINE F01AAF(A,LDA,N,X,LDX,P,IFAIL)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C
C     APPROXIMATE INVERSE OF A REAL MATRIX
C     1ST AUGUST 1971
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01AAF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDA, LDX, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,N), P(N), X(LDX,N)
C     .. Local Scalars ..
      DOUBLE PRECISION  D1, EPS
      INTEGER           I, ISAVE, J, JP
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F03AFF, DGEMV, DSCAL, DSWAP, DTRMV
C     .. Executable Statements ..
      ISAVE = IFAIL
      IFAIL = 1
C     Compute LU factorization of A
      EPS = X02AJF()
      CALL F03AFF(N,EPS,A,LDA,D1,I,P,IFAIL)
      IF (IFAIL.EQ.0) GO TO 20
      IFAIL = P01ABF(ISAVE,IFAIL,SRNAME,0,P01REC)
      RETURN
   20 CONTINUE
C     Copy U to array X
      DO 80 J = 1, N
         DO 40 I = 1, J - 1
            X(I,J) = A(I,J)
   40    CONTINUE
         X(J,J) = 1.0D0
         DO 60 I = J + 1, N
            X(I,J) = 0.0D0
   60    CONTINUE
   80 CONTINUE
C     Compute inverse of U in array X, overwriting U
      DO 120 J = 1, N
         CALL DTRMV('U','N','U',J-1,X,LDX,X(1,J),1)
         DO 100 I = 1, J - 1
            X(I,J) = -X(I,J)
  100    CONTINUE
  120 CONTINUE
C     Compute  X * inv(L)
      DO 140 J = N, 1, -1
         IF (J.LT.N) CALL DGEMV('N',N,N-J,-1.0D0,X(1,J+1),LDX,A(J+1,J),
     *                          1,1.0D0,X(1,J),1)
         CALL DSCAL(N,1.0D0/A(J,J),X(1,J),1)
  140 CONTINUE
C     Permute columns of X
      DO 160 J = N, 1, -1
         JP = P(J) + 0.5D0
         IF (JP.NE.J) CALL DSWAP(N,X(1,J),1,X(1,JP),1)
  160 CONTINUE
      RETURN
      END
      SUBROUTINE F03AFF(N,EPS,A,IA,D1,ID,P,IFAIL)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 3 REVISED.
C     MARK 4.5 REVISED
C     MARK 11 REVISED. VECTORISATION (JAN 1984).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 12 REVISED (LEVEL 2 BLAS) (MAR 1986)
C     MARK 12 REVISED. EXTENDED BLAS (JUNE 1986)
C
C     UNSYMDET
C     THE UNSYMMETRIC MATRIX, A, IS STORED IN THE N*N ARRAY A(I,J),
C     I=1,N, J=1,N. THE DECOMPOSITION A=LU, WHERE L IS A
C     LOWER TRIANGULAR MATRIX AND U IS A UNIT UPPER TRIANGULAR
C     MATRIX, IS PERFORMED AND OVERWRITTEN ON A, OMITTING THE UNIT
C     DIAGONAL OF U. A RECORD OF ANY INTERCHANGES MADE TO THE ROWS
C     OF A IS KEPT IN P(I), I=1,N, SUCH THAT THE I-TH ROW AND
C     THE P(I)-TH ROW WERE INTERCHANGED AT THE I-TH STEP. THE
C     DETERMINANT, D1 * 2.0**ID, OF A IS ALSO COMPUTED. THE
C     SUBROUTINE
C     WILL FAIL IF A, MODIFIED BY THE ROUNDING ERRORS, IS SINGULAR
C     OR ALMOST SINGULAR. SETS IFAIL = 0 IF SUCCESSFUL ELSE IFAIL =
C     1.
C     1ST DECEMBER 1971
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F03AFF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  D1, EPS
      INTEGER           IA, ID, IFAIL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N), P(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  X, Y
      INTEGER           I, J, K, L
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          DGEMV, DTRSV
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Executable Statements ..
      DO 20 I = 1, N
         P(I) = 0.0D0
   20 CONTINUE
      DO 60 J = 1, N
         DO 40 I = 1, N
            P(I) = P(I) + A(I,J)**2
   40    CONTINUE
   60 CONTINUE
      DO 80 I = 1, N
         IF (P(I).LE.0.0D0) GO TO 240
         P(I) = 1.0D0/SQRT(P(I))
   80 CONTINUE
      D1 = 1.0D0
      ID = 0
      DO 220 K = 1, N
         L = K
         X = 0.0D0
         DO 100 I = K, N
            Y = ABS(A(I,K)*P(I))
            IF (Y.LE.X) GO TO 100
            X = Y
            L = I
  100    CONTINUE
         IF (L.EQ.K) GO TO 140
         D1 = -D1
         DO 120 J = 1, N
            Y = A(K,J)
            A(K,J) = A(L,J)
            A(L,J) = Y
  120    CONTINUE
         P(L) = P(K)
  140    P(K) = L
         D1 = D1*A(K,K)
         IF (X.LT.8.0D0*EPS) GO TO 240
  160    IF (ABS(D1).LT.1.0D0) GO TO 180
         D1 = D1*0.0625D0
         ID = ID + 4
         GO TO 160
  180    IF (ABS(D1).GE.0.0625D0) GO TO 200
         D1 = D1*16.0D0
         ID = ID - 4
         GO TO 180
  200    IF (K.LT.N) THEN
            CALL DTRSV('L','N','N',K,A,IA,A(1,K+1),1)
            CALL DGEMV('N',N-K,K,-1.0D0,A(K+1,1),IA,A(1,K+1),1,1.0D0,
     *                 A(K+1,K+1),1)
         END IF
  220 CONTINUE
      IFAIL = 0
      RETURN
  240 IFAIL = P01ABF(IFAIL,1,SRNAME,0,P01REC)
      RETURN
      END
      SUBROUTINE F06EGF( N, X, INCX, Y, INCY )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      ENTRY      DSWAP ( N, X, INCX, Y, INCY )
C     .. Scalar Arguments ..
      INTEGER            INCX, INCY, N
C     .. Array Arguments ..
      DOUBLE PRECISION   X( * ), Y( * )
C     ..
C
C  F06EGF performs the operations
C
C     temp := x,   x := y,   y := temp.
C
C
C  Nag Fortran 77 version of the Blas routine DSWAP.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 26-November-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, IX, IY
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ( INCX.EQ.INCY ).AND.( INCY.GT.0 ) )THEN
            DO 10, IY = 1, 1 + ( N - 1 )*INCY, INCY
               TEMP    = X( IY )
               X( IY ) = Y( IY )
               Y( IY ) = TEMP
   10       CONTINUE
         ELSE
            IF( INCX.GE.0 )THEN
               IX = 1
            ELSE
               IX = 1 - ( N - 1 )*INCX
            END IF
            IF( INCY.GT.0 )THEN
               DO 20, IY = 1, 1 + ( N - 1 )*INCY, INCY
                  TEMP    = X( IX )
                  X( IX ) = Y( IY )
                  Y( IY ) = TEMP
                  IX      = IX      + INCX
   20          CONTINUE
            ELSE
               IY = 1 - ( N - 1 )*INCY
               DO 30, I = 1, N
                  TEMP    = X( IX )
                  X( IX ) = Y( IY )
                  Y( IY ) = TEMP
                  IY      = IY      + INCY
                  IX      = IX      + INCX
   30          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06EGF. ( DSWAP )
C
      END


      SUBROUTINE F06PFF( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     AXP4 VERSION FOR VECTOR MACHINES
C     .. Entry Points ..
      ENTRY      DTRMV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
C     .. Scalar Arguments ..
      INTEGER            INCX, LDA, N
      CHARACTER*1        DIAG, TRANS, UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * )
C     ..
C
C  Purpose
C  =======
C
C  DTRMV  performs one of the matrix-vector operations
C
C     x := A*x,   or   x := A'*x,
C
C  where x is n element vector and A is an n by n unit, or non-unit,
C  upper or lower triangular matrix.
C
C  Parameters
C  ==========
C
C  UPLO   - CHARACTER*1.
C           On entry, UPLO specifies whether the matrix is an upper or
C           lower triangular matrix as follows:
C
C              UPLO = 'U' or 'u'   A is an upper triangular matrix.
C
C              UPLO = 'L' or 'l'   A is a lower triangular matrix.
C
C           Unchanged on exit.
C
C  TRANS  - CHARACTER*1.
C           On entry, TRANS specifies the operation to be performed as
C           follows:
C
C              TRANS = 'N' or 'n'   x := A*x.
C
C              TRANS = 'T' or 't'   x := A'*x.
C
C              TRANS = 'C' or 'c'   x := A'*x.
C
C           Unchanged on exit.
C
C  DIAG   - CHARACTER*1.
C           On entry, DIAG specifies whether or not A is unit
C           triangular as follows:
C
C              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
C
C              DIAG = 'N' or 'n'   A is not assumed to be unit
C                                  triangular.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the order of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
C           Before entry with  UPLO = 'U' or 'u', the leading n by n
C           upper triangular part of the array A must contain the upper
C           triangular matrix and the strictly lower triangular part of
C           A is not referenced.
C           Before entry with UPLO = 'L' or 'l', the leading n by n
C           lower triangular part of the array A must contain the lower
C           triangular matrix and the strictly upper triangular part of
C           A is not referenced.
C           Note that when  DIAG = 'U' or 'u', the diagonal elements of
C           A are not referenced either, but are assumed to be unity.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, n ).
C           Unchanged on exit.
C
C  X      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the n
C           element vector x. On exit, X is overwritten with the
C           tranformed vector x.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C  -- DO-loops unrolled on 20-November-1986.
C     Peter Mayes, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP1, TEMP2, TEMP3, TEMP4
      INTEGER            I, IND1, IND2, IND3, INFO, IX, J, JX, KX, N4
      LOGICAL            NOUNIT
C     .. External Subroutines ..
      EXTERNAL           F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC          MAX, MOD
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( .NOT.(UPLO .EQ.'U' .OR. UPLO .EQ.'u').AND.
     $         .NOT.(UPLO .EQ.'L' .OR. UPLO .EQ.'l')      )THEN
         INFO = 1
      ELSE IF( .NOT.(TRANS.EQ.'N' .OR. TRANS.EQ.'n').AND.
     $         .NOT.(TRANS.EQ.'T' .OR. TRANS.EQ.'t').AND.
     $         .NOT.(TRANS.EQ.'C' .OR. TRANS.EQ.'c')      )THEN
         INFO = 2
      ELSE IF( .NOT.(DIAG .EQ.'U' .OR. DIAG .EQ.'u').AND.
     $         .NOT.(DIAG .EQ.'N' .OR. DIAG .EQ.'n')      )THEN
         INFO = 3
      ELSE IF( N.LT.0 )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      END IF
      IF( INFO.NE.0 )THEN
         CALL F06AAZ( 'F06PFF/DTRMV ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 )
     $   RETURN
C
      NOUNIT = (DIAG.EQ.'N' .OR. DIAG.EQ.'n')
C
C     Set up the start point in X if the increment is not unity. This
C     will be  ( N - 1 )*INCX  too small for descending loops.
C
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
C
C     Start the operations. In this version the inner loops are all
C     equivalent to AXPY operations.
C
      IF( (TRANS.EQ.'N' .OR. TRANS.EQ.'n') )THEN
C
C        Form  x := A*x.
C
         IF( (UPLO.EQ.'U' .OR. UPLO.EQ.'u') )THEN
            IF( INCX.EQ.1 )THEN
C**** U n r o l l   t o   D e p t h   4 ********** L O O P  A **********
               N4 = MOD( N, 4 ) + 1
C**** Take care of the beginning triangle ******************************
               IF( N4.GE.4 )THEN
                  IF( NOUNIT )
     $               X( N4 - 3 ) = X( N4 - 3 )*A( N4 - 3, N4 - 3 )
                  X( N4 - 3 ) = X( N4 - 3 ) + X( N4 - 2 )
     $                          *A( N4 - 3, N4 - 2 ) + X( N4 - 1 )
     $                          *A( N4 - 3, N4 - 1 )
               END IF
               IF( N4.GE.3 )THEN
                  IF( NOUNIT )
     $               X( N4 - 2 ) = X( N4 - 2 )*A( N4 - 2, N4 - 2 )
                  X( N4 - 2 ) = X( N4 - 2 ) + X( N4 - 1 )
     $                          *A( N4 - 2, N4 - 1 )
               END IF
               IF( N4.GE.2 )THEN
                  IF( NOUNIT )
     $               X( N4 - 1 ) = X( N4 - 1 )*A( N4 - 1, N4 - 1 )
               END IF
C**** Unrolled main loop ***********************************************
               DO 20, J = N4, N, 4
                  TEMP1 = X( J )
                  TEMP2 = X( J + 1 )
                  TEMP3 = X( J + 2 )
                  TEMP4 = X( J + 3 )
                  IF( TEMP1.NE.ZERO.OR.TEMP2.NE.ZERO.OR.TEMP3.NE.
     $                ZERO.OR.TEMP4.NE.ZERO )THEN
                     DO 10, I = 1, J - 1
                        X( I ) = ( ( ( ( X( I ) + TEMP1*A( I, J ) )
     $                           + TEMP2*A( I, J + 1 ) )
     $                           + TEMP3*A( I, J + 2 ) )
     $                           + TEMP4*A( I, J + 3 ) )
   10                CONTINUE
                  END IF
C**** 4 by 4 triangle left over at end of block ************************
                  IF( NOUNIT )
     $               X( J ) = X( J )*A( J, J )
                  X( J ) = X( J ) + X( J + 1 )*A( J, J + 1 ) +
     $                     X( J + 2 )*A( J, J + 2 ) + X( J + 3 )
     $                     *A( J, J + 3 )
                  IF( NOUNIT )
     $               X( J + 1 ) = X( J + 1 )*A( J + 1, J + 1 )
                  X( J + 1 ) = X( J + 1 ) + X( J + 2 )*A( J + 1, J + 2 )
     $                          + X( J + 3 )*A( J + 1, J + 3 )
                  IF( NOUNIT )
     $               X( J + 2 ) = X( J + 2 )*A( J + 2, J + 2 )
                  X( J + 2 ) = X( J + 2 ) + X( J + 3 )*A( J + 2, J + 3 )
                  IF( NOUNIT )
     $               X( J + 3 ) = X( J + 3 )*A( J + 3, J + 3 )
   20          CONTINUE
            ELSE
C**** U n r o l l   t o   D e p t h   4 ********** L O O P  B **********
               JX = KX
               N4 = MOD( N, 4 ) + 1
C**** Take care of the beginning triangle ******************************
               IF( N4.EQ.4 )THEN
                  IND1 = JX + 2*INCX
                  IND2 = JX + INCX
                  IND3 = JX
                  JX = JX + 3*INCX
               ELSE IF( N4.EQ.3 )THEN
                  IND1 = JX + INCX
                  IND2 = JX
                  JX = JX + 2*INCX
               ELSE IF( N4.EQ.2 )THEN
                  IND1 = JX
                  JX = JX + INCX
               END IF
               IF( N4.GE.4 )THEN
                  IF( NOUNIT )
     $               X( IND3 ) = X( IND3 )*A( N4 - 3, N4 - 3 )
                  X( IND3 ) = X( IND3 ) + X( IND2 )*A( N4 - 3, N4 - 2 )
     $                         + X( IND1 )*A( N4 - 3, N4 - 1 )
               END IF
               IF( N4.GE.3 )THEN
                  IF( NOUNIT )
     $               X( IND2 ) = X( IND2 )*A( N4 - 2, N4 - 2 )
                  X( IND2 ) = X( IND2 ) + X( IND1 )*A( N4 - 2, N4 - 1 )
               END IF
               IF( N4.GE.2 )THEN
                  IF( NOUNIT )
     $               X( IND1 ) = X( IND1 )*A( N4 - 1, N4 - 1 )
               END IF
C**** Unrolled main loop ***********************************************
               DO 40, J = N4, N, 4
                  TEMP1 = X( JX )
                  TEMP2 = X( JX + INCX )
                  TEMP3 = X( JX + 2*INCX )
                  TEMP4 = X( JX + 3*INCX )
                  IF( TEMP1.NE.ZERO.OR.TEMP2.NE.ZERO.OR.TEMP3.NE.
     $                ZERO.OR.TEMP4.NE.ZERO )THEN
                     IX = KX
                     DO 30, I = 1, J - 1
                        X( IX ) = ( ( ( ( X( IX ) + TEMP1*A( I, J ) )
     $                            + TEMP2*A( I, J + 1 ) )
     $                            + TEMP3*A( I, J + 2 ) )
     $                            + TEMP4*A( I, J + 3 ) )
                        IX = IX + INCX
   30                CONTINUE
                  END IF
C**** 4 by 4 triangle left over at end of block ************************
                  IF( NOUNIT )
     $               X( JX ) = X( JX )*A( J, J )
                  X( JX ) = X( JX ) + X( JX + INCX )*A( J, J + 1 ) +
     $                      X( JX + 2*INCX )*A( J, J + 2 ) +
     $                      X( JX + 3*INCX )*A( J, J + 3 )
                  IF( NOUNIT )
     $               X( JX + INCX ) = X( JX + INCX )*A( J + 1, J + 1 )
                  X( JX + INCX ) = X( JX + INCX ) + X( JX + 2*INCX )
     $                             *A( J + 1, J + 2 ) + X( JX + 3*INCX )
     $                             *A( J + 1, J + 3 )
                  IF( NOUNIT )
     $               X( JX + 2*INCX ) = X( JX + 2*INCX )
     $                                  *A( J + 2, J + 2 )
                  X( JX + 2*INCX ) = X( JX + 2*INCX ) + X( JX + 3*INCX )
     $                               *A( J + 2, J + 3 )
                  IF( NOUNIT )
     $               X( JX + 3*INCX ) = X( JX + 3*INCX )
     $                                  *A( J + 3, J + 3 )
                  JX = JX + 4*INCX
   40          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
C**** U n r o l l   t o   D e p t h   4 ********** L O O P  C **********
               N4 = 4*( N/4 )
C**** Take care of the beginning triangle ******************************
               IF( N - N4.GE.3 )THEN
                  IF( NOUNIT )
     $               X( N4 + 3 ) = X( N4 + 3 )*A( N4 + 3, N4 + 3 )
                  X( N4 + 3 ) = X( N4 + 3 ) + X( N4 + 2 )
     $                          *A( N4 + 3, N4 + 2 ) + X( N4 + 1 )
     $                          *A( N4 + 3, N4 + 1 )
               END IF
               IF( N - N4.GE.2 )THEN
                  IF( NOUNIT )
     $               X( N4 + 2 ) = X( N4 + 2 )*A( N4 + 2, N4 + 2 )
                  X( N4 + 2 ) = X( N4 + 2 ) + X( N4 + 1 )
     $                          *A( N4 + 2, N4 + 1 )
               END IF
               IF( N - N4.GE.1 )THEN
                  IF( NOUNIT )
     $               X( N4 + 1 ) = X( N4 + 1 )*A( N4 + 1, N4 + 1 )
               END IF
C**** Unrolled main loop ***********************************************
               DO 60, J = N4, 1, -4
                  TEMP1 = X( J )
                  TEMP2 = X( J - 1 )
                  TEMP3 = X( J - 2 )
                  TEMP4 = X( J - 3 )
                  IF( TEMP1.NE.ZERO.OR.TEMP2.NE.ZERO.OR.TEMP3.NE.
     $                ZERO.OR.TEMP4.NE.ZERO )THEN
                     DO 50, I = N, J + 1, -1
                        X( I ) = ( ( ( ( X( I ) + TEMP1*A( I, J ) )
     $                           + TEMP2*A( I, J - 1 ) )
     $                           + TEMP3*A( I, J - 2 ) )
     $                           + TEMP4*A( I, J - 3 ) )
   50                CONTINUE
C**** 4 by 4 triangle left over at end of block ************************
                     IF( NOUNIT )
     $                  X( J ) = X( J )*A( J, J )
                     X( J ) = X( J ) + X( J - 1 )*A( J, J - 1 ) +
     $                        X( J - 2 )*A( J, J - 2 ) + X( J - 3 )
     $                        *A( J, J - 3 )
                     IF( NOUNIT )
     $                  X( J - 1 ) = X( J - 1 )*A( J - 1, J - 1 )
                     X( J - 1 ) = X( J - 1 ) + X( J - 2 )
     $                            *A( J - 1, J - 2 ) + X( J - 3 )
     $                            *A( J - 1, J - 3 )
                     IF( NOUNIT )
     $                  X( J - 2 ) = X( J - 2 )*A( J - 2, J - 2 )
                     X( J - 2 ) = X( J - 2 ) + X( J - 3 )
     $                            *A( J - 2, J - 3 )
                     IF( NOUNIT )
     $                  X( J - 3 ) = X( J - 3 )*A( J - 3, J - 3 )
                  END IF
   60          CONTINUE
            ELSE
C**** U n r o l l   t o   D e p t h   4 ********** L O O P  D **********
               KX = KX + ( N - 1 )*INCX
               JX = KX
               N4 = 4*( N/4 )
C**** Take care of the beginning triangle ******************************
               IF( N - N4.EQ.3 )THEN
                  IND1 = JX - 2*INCX
                  IND2 = JX - INCX
                  IND3 = JX
                  JX = JX - 3*INCX
               ELSE IF( N - N4.EQ.2 )THEN
                  IND1 = JX - INCX
                  IND2 = JX
                  JX = JX - 2*INCX
               ELSE IF( N - N4.EQ.1 )THEN
                  IND1 = JX
                  JX = JX - INCX
               END IF
               IF( N - N4.GE.3 )THEN
                  IF( NOUNIT )
     $               X( IND3 ) = X( IND3 )*A( N4 + 3, N4 + 3 )
                  X( IND3 ) = X( IND3 ) + X( IND2 )*A( N4 + 3, N4 + 2 )
     $                         + X( IND1 )*A( N4 + 3, N4 + 1 )
               END IF
               IF( N - N4.GE.2 )THEN
                  IF( NOUNIT )
     $               X( IND2 ) = X( IND2 )*A( N4 + 2, N4 + 2 )
                  X( IND2 ) = X( IND2 ) + X( IND1 )*A( N4 + 2, N4 + 1 )
               END IF
               IF( N - N4.GE.1 )THEN
                  IF( NOUNIT )
     $               X( IND1 ) = X( IND1 )*A( N4 + 1, N4 + 1 )
               END IF
C**** Unrolled main loop ***********************************************
               DO 80, J = N4, 1, -4
                  TEMP1 = X( JX )
                  TEMP2 = X( JX - INCX )
                  TEMP3 = X( JX - 2*INCX )
                  TEMP4 = X( JX - 3*INCX )
                  IF( TEMP1.NE.ZERO.OR.TEMP2.NE.ZERO.OR.TEMP3.NE.
     $                ZERO.OR.TEMP4.NE.ZERO )THEN
                     IX = KX
                     DO 70, I = N, J + 1, -1
                        X( IX ) = ( ( ( ( X( IX ) + TEMP1*A( I, J ) )
     $                            + TEMP2*A( I, J - 1 ) )
     $                            + TEMP3*A( I, J - 2 ) )
     $                            + TEMP4*A( I, J - 3 ) )
                        IX = IX - INCX
   70                CONTINUE
C**** 4 by 4 triangle left over at end of block ************************
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )*A( J, J )
                     X( JX ) = X( JX ) + X( JX - INCX )*A( J, J - 1 ) +
     $                         X( JX - 2*INCX )*A( J, J - 2 ) +
     $                         X( JX - 3*INCX )*A( J, J - 3 )
                     IF( NOUNIT )
     $                  X( JX - INCX ) = X( JX - INCX )
     $                                   *A( J - 1, J - 1 )
                     X( JX - INCX ) = X( JX - INCX ) + X( JX - 2*INCX )
     $                                *A( J - 1, J - 2 ) +
     $                                X( JX - 3*INCX )*A( J - 1, J - 3 )
                     IF( NOUNIT )
     $                  X( JX - 2*INCX ) = X( JX - 2*INCX )
     $                                     *A( J - 2, J - 2 )
                     X( JX - 2*INCX ) = X( JX - 2*INCX ) +
     $                                  X( JX - 3*INCX )
     $                                  *A( J - 2, J - 3 )
                     IF( NOUNIT )
     $                  X( JX - 3*INCX ) = X( JX - 3*INCX )
     $                                     *A( J - 3, J - 3 )
                  END IF
                  JX = JX - 4*INCX
   80          CONTINUE
            END IF
         END IF
      ELSE
C
C        Form  x := A'*x.
C
         IF( (UPLO.EQ.'U' .OR. UPLO.EQ.'u') )THEN
            IF( INCX.EQ.1 )THEN
C**** U n r o l l   t o   D e p t h   4 ********** L O O P  E **********
               N4 = 4*( N/4 )
C**** Take care of the beginning triangle ******************************
               IF( N - N4.GE.3 )THEN
                  IF( NOUNIT )
     $               X( N4 + 3 ) = X( N4 + 3 )*A( N4 + 3, N4 + 3 )
                  X( N4 + 3 ) = X( N4 + 3 ) + X( N4 + 2 )
     $                          *A( N4 + 2, N4 + 3 ) + X( N4 + 1 )
     $                          *A( N4 + 1, N4 + 3 )
               END IF
               IF( N - N4.GE.2 )THEN
                  IF( NOUNIT )
     $               X( N4 + 2 ) = X( N4 + 2 )*A( N4 + 2, N4 + 2 )
                  X( N4 + 2 ) = X( N4 + 2 ) + X( N4 + 1 )
     $                          *A( N4 + 1, N4 + 2 )
               END IF
               IF( N - N4.GE.1 )THEN
                  IF( NOUNIT )
     $               X( N4 + 1 ) = X( N4 + 1 )*A( N4 + 1, N4 + 1 )
               END IF
C**** Unrolled main loop ***********************************************
               DO 100, J = N4, 1, -4
                  TEMP1 = X( J )
                  TEMP2 = X( J - 1 )
                  TEMP3 = X( J - 2 )
                  TEMP4 = X( J - 3 )
                  IF( TEMP1.NE.ZERO.OR.TEMP2.NE.ZERO.OR.TEMP3.NE.
     $                ZERO.OR.TEMP4.NE.ZERO )THEN
                     DO 90, I = N, J + 1, -1
                        X( I ) = ( ( ( ( X( I ) + TEMP1*A( J, I ) )
     $                           + TEMP2*A( J - 1, I ) )
     $                           + TEMP3*A( J - 2, I ) )
     $                           + TEMP4*A( J - 3, I ) )
   90                CONTINUE
C**** 4 by 4 triangle left over at end of block ************************
                     IF( NOUNIT )
     $                  X( J ) = X( J )*A( J, J )
                     X( J ) = X( J ) + X( J - 1 )*A( J - 1, J ) +
     $                        X( J - 2 )*A( J - 2, J ) + X( J - 3 )
     $                        *A( J - 3, J )
                     IF( NOUNIT )
     $                  X( J - 1 ) = X( J - 1 )*A( J - 1, J - 1 )
                     X( J - 1 ) = X( J - 1 ) + X( J - 2 )
     $                            *A( J - 2, J - 1 ) + X( J - 3 )
     $                            *A( J - 3, J - 1 )
                     IF( NOUNIT )
     $                  X( J - 2 ) = X( J - 2 )*A( J - 2, J - 2 )
                     X( J - 2 ) = X( J - 2 ) + X( J - 3 )
     $                            *A( J - 3, J - 2 )
                     IF( NOUNIT )
     $                  X( J - 3 ) = X( J - 3 )*A( J - 3, J - 3 )
                  END IF
  100          CONTINUE
            ELSE
C**** U n r o l l   t o   D e p t h   4 ********** L O O P  F **********
               KX = KX + ( N - 1 )*INCX
               JX = KX
               N4 = 4*( N/4 )
C**** Take care of the beginning triangle ******************************
               IF( N - N4.EQ.3 )THEN
                  IND1 = JX - 2*INCX
                  IND2 = JX - INCX
                  IND3 = JX
                  JX = JX - 3*INCX
               ELSE IF( N - N4.EQ.2 )THEN
                  IND1 = JX - INCX
                  IND2 = JX
                  JX = JX - 2*INCX
               ELSE IF( N - N4.EQ.1 )THEN
                  IND1 = JX
                  JX = JX - INCX
               END IF
               IF( N - N4.GE.3 )THEN
                  IF( NOUNIT )
     $               X( IND3 ) = X( IND3 )*A( N4 + 3, N4 + 3 )
                  X( IND3 ) = X( IND3 ) + X( IND2 )*A( N4 + 2, N4 + 3 )
     $                         + X( IND1 )*A( N4 + 1, N4 + 3 )
               END IF
               IF( N - N4.GE.2 )THEN
                  IF( NOUNIT )
     $               X( IND2 ) = X( IND2 )*A( N4 + 2, N4 + 2 )
                  X( IND2 ) = X( IND2 ) + X( IND1 )*A( N4 + 1, N4 + 2 )
               END IF
               IF( N - N4.GE.1 )THEN
                  IF( NOUNIT )
     $               X( IND1 ) = X( IND1 )*A( N4 + 1, N4 + 1 )
               END IF
C**** Unrolled main loop ***********************************************
               DO 120, J = N4, 1, -4
                  TEMP1 = X( JX )
                  TEMP2 = X( JX - INCX )
                  TEMP3 = X( JX - 2*INCX )
                  TEMP4 = X( JX - 3*INCX )
                  IF( TEMP1.NE.ZERO.OR.TEMP2.NE.ZERO.OR.TEMP3.NE.
     $                ZERO.OR.TEMP4.NE.ZERO )THEN
                     IX = KX
                     DO 110, I = N, J + 1, -1
                        X( IX ) = ( ( ( ( X( IX ) + TEMP1*A( J, I ) )
     $                            + TEMP2*A( J - 1, I ) )
     $                            + TEMP3*A( J - 2, I ) )
     $                            + TEMP4*A( J - 3, I ) )
                        IX = IX - INCX
  110                CONTINUE
C**** 4 by 4 triangle left over at end of block ************************
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )*A( J, J )
                     X( JX ) = X( JX ) + X( JX - INCX )*A( J - 1, J ) +
     $                         X( JX - 2*INCX )*A( J - 2, J ) +
     $                         X( JX - 3*INCX )*A( J - 3, J )
                     IF( NOUNIT )
     $                  X( JX - INCX ) = X( JX - INCX )
     $                                   *A( J - 1, J - 1 )
                     X( JX - INCX ) = X( JX - INCX ) + X( JX - 2*INCX )
     $                                *A( J - 2, J - 1 ) +
     $                                X( JX - 3*INCX )*A( J - 3, J - 1 )
                     IF( NOUNIT )
     $                  X( JX - 2*INCX ) = X( JX - 2*INCX )
     $                                     *A( J - 2, J - 2 )
                     X( JX - 2*INCX ) = X( JX - 2*INCX ) +
     $                                  X( JX - 3*INCX )
     $                                  *A( J - 3, J - 2 )
                     IF( NOUNIT )
     $                  X( JX - 3*INCX ) = X( JX - 3*INCX )
     $                                     *A( J - 3, J - 3 )
                  END IF
                  JX = JX - 4*INCX
  120          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
C**** U n r ol l   t o   D e p t h  4 ********** L O O P  G ************
               N4 = MOD( N, 4 ) + 1
C**** Take care of the beginning triangle ******************************
               IF( N4.GE.4 )THEN
                  IF( NOUNIT )
     $               X( N4 - 3 ) = X( N4 - 3 )*A( N4 - 3, N4 - 3 )
                  X( N4 - 3 ) = X( N4 - 3 ) + X( N4 - 2 )
     $                          *A( N4 - 2, N4 - 3 ) + X( N4 - 1 )
     $                          *A( N4 - 1, N4 - 3 )
               END IF
               IF( N4.GE.3 )THEN
                  IF( NOUNIT )
     $               X( N4 - 2 ) = X( N4 - 2 )*A( N4 - 2, N4 - 2 )
                  X( N4 - 2 ) = X( N4 - 2 ) + X( N4 - 1 )
     $                          *A( N4 - 1, N4 - 2 )
               END IF
               IF( N4.GE.2 )THEN
                  IF( NOUNIT )
     $               X( N4 - 1 ) = X( N4 - 1 )*A( N4 - 1, N4 - 1 )
               END IF
C**** Unrolled main loop ***********************************************
               DO 140, J = N4, N, 4
                  TEMP1 = X( J )
                  TEMP2 = X( J + 1 )
                  TEMP3 = X( J + 2 )
                  TEMP4 = X( J + 3 )
                  IF( TEMP1.NE.ZERO.OR.TEMP2.NE.ZERO.OR.TEMP3.NE.
     $                ZERO.OR.TEMP4.NE.ZERO )THEN
                     DO 130, I = 1, J - 1
                        X( I ) = ( ( ( ( X( I ) + TEMP1*A( J, I ) )
     $                           + TEMP2*A( J + 1, I ) )
     $                           + TEMP3*A( J + 2, I ) )
     $                           + TEMP4*A( J + 3, I ) )
  130                CONTINUE
                  END IF
C**** 4 by 4 triangle left over at end of block ************************
                  IF( NOUNIT )
     $               X( J ) = X( J )*A( J, J )
                  X( J ) = X( J ) + X( J + 1 )*A( J + 1, J ) +
     $                     X( J + 2 )*A( J + 2, J ) + X( J + 3 )
     $                     *A( J + 3, J )
                  IF( NOUNIT )
     $               X( J + 1 ) = X( J + 1 )*A( J + 1, J + 1 )
                  X( J + 1 ) = X( J + 1 ) + X( J + 2 )*A( J + 2, J + 1 )
     $                          + X( J + 3 )*A( J + 3, J + 1 )
                  IF( NOUNIT )
     $               X( J + 2 ) = X( J + 2 )*A( J + 2, J + 2 )
                  X( J + 2 ) = X( J + 2 ) + X( J + 3 )*A( J + 3, J + 2 )
                  IF( NOUNIT )
     $               X( J + 3 ) = X( J + 3 )*A( J + 3, J + 3 )
  140          CONTINUE
            ELSE
C**** U n r o l l   t o   D e p t h   4 ********** L O O P  H **********
               JX = KX
               N4 = MOD( N, 4 ) + 1
C**** Take care of the beginning triangle ******************************
               IF( N4.EQ.4 )THEN
                  IND1 = JX + 2*INCX
                  IND2 = JX + INCX
                  IND3 = JX
                  JX = JX + 3*INCX
               ELSE IF( N4.EQ.3 )THEN
                  IND1 = JX + INCX
                  IND2 = JX
                  JX = JX + 2*INCX
               ELSE IF( N4.EQ.2 )THEN
                  IND1 = JX
                  JX = JX + INCX
               END IF
               IF( N4.GE.4 )THEN
                  IF( NOUNIT )
     $               X( IND3 ) = X( IND3 )*A( N4 - 3, N4 - 3 )
                  X( IND3 ) = X( IND3 ) + X( IND2 )*A( N4 - 2, N4 - 3 )
     $                         + X( IND1 )*A( N4 - 1, N4 - 3 )
               END IF
               IF( N4.GE.3 )THEN
                  IF( NOUNIT )
     $               X( IND2 ) = X( IND2 )*A( N4 - 2, N4 - 2 )
                  X( IND2 ) = X( IND2 ) + X( IND1 )*A( N4 - 1, N4 - 2 )
               END IF
               IF( N4.GE.2 )THEN
                  IF( NOUNIT )
     $               X( IND1 ) = X( IND1 )*A( N4 - 1, N4 - 1 )
               END IF
C**** Unrolled main loop ***********************************************
               DO 160, J = N4, N, 4
                  TEMP1 = X( JX )
                  TEMP2 = X( JX + INCX )
                  TEMP3 = X( JX + 2*INCX )
                  TEMP4 = X( JX + 3*INCX )
                  IF( TEMP1.NE.ZERO.OR.TEMP2.NE.ZERO.OR.TEMP3.NE.
     $                ZERO.OR.TEMP4.NE.ZERO )THEN
                     IX = KX
                     DO 150, I = 1, J - 1
                        X( IX ) = ( ( ( ( X( IX ) + TEMP1*A( J, I ) )
     $                            + TEMP2*A( J + 1, I ) )
     $                            + TEMP3*A( J + 2, I ) )
     $                            + TEMP4*A( J + 3, I ) )
                        IX = IX + INCX
  150                CONTINUE
                  END IF
C**** 4 by 4 triangle left over at end of block ************************
                  IF( NOUNIT )
     $               X( JX ) = X( JX )*A( J, J )
                  X( JX ) = X( JX ) + X( JX + INCX )*A( J + 1, J ) +
     $                      X( JX + 2*INCX )*A( J + 2, J ) +
     $                      X( JX + 3*INCX )*A( J + 3, J )
                  IF( NOUNIT )
     $               X( JX + INCX ) = X( JX + INCX )*A( J + 1, J + 1 )
                  X( JX + INCX ) = X( JX + INCX ) + X( JX + 2*INCX )
     $                             *A( J + 2, J + 1 ) + X( JX + 3*INCX )
     $                             *A( J + 3, J + 1 )
                  IF( NOUNIT )
     $               X( JX + 2*INCX ) = X( JX + 2*INCX )
     $                                  *A( J + 2, J + 2 )
                  X( JX + 2*INCX ) = X( JX + 2*INCX ) + X( JX + 3*INCX )
     $                               *A( J + 3, J + 2 )
                  IF( NOUNIT )
     $               X( JX + 3*INCX ) = X( JX + 3*INCX )
     $                                  *A( J + 3, J + 3 )
                  JX = JX + 4*INCX
  160          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06PFF (DTRMV ).
C
      END
      SUBROUTINE F06PJF( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     AXP4 VERSION FOR VECTOR MACHINES
C     .. Entry Points ..
      ENTRY      DTRSV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
C     .. Scalar Arguments ..
      INTEGER            INCX, LDA, N
      CHARACTER*1        DIAG, TRANS, UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * )
C     ..
C
C  Purpose
C  =======
C
C  DTRSV  solves one of the systems of equations
C
C     A*x = b,   or   A'*x = b,
C
C  where b and x are n element vectors and A is an n by n unit, or
C  non-unit, upper or lower triangular matrix.
C
C  No test for singularity or near-singularity is included in this
C  routine. Such tests must be performed before calling this routine.
C
C  Parameters
C  ==========
C
C  UPLO   - CHARACTER*1.
C           On entry, UPLO specifies whether the matrix is an upper or
C           lower triangular matrix as follows:
C
C              UPLO = 'U' or 'u'   A is an upper triangular matrix.
C
C              UPLO = 'L' or 'l'   A is a lower triangular matrix.
C
C           Unchanged on exit.
C
C  TRANS  - CHARACTER*1.
C           On entry, TRANS specifies the equations to be solved as
C           follows:
C
C              TRANS = 'N' or 'n'   A*x = b.
C
C              TRANS = 'T' or 't'   A'*x = b.
C
C              TRANS = 'C' or 'c'   A'*x = b.
C
C           Unchanged on exit.
C
C  DIAG   - CHARACTER*1.
C           On entry, DIAG specifies whether or not A is unit
C           triangular as follows:
C
C              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
C
C              DIAG = 'N' or 'n'   A is not assumed to be unit
C                                  triangular.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the order of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
C           Before entry with  UPLO = 'U' or 'u', the leading n by n
C           upper triangular part of the array A must contain the upper
C           triangular matrix and the strictly lower triangular part of
C           A is not referenced.
C           Before entry with UPLO = 'L' or 'l', the leading n by n
C           lower triangular part of the array A must contain the lower
C           triangular matrix and the strictly upper triangular part of
C           A is not referenced.
C           Note that when  DIAG = 'U' or 'u', the diagonal elements of
C           A are not referenced either, but are assumed to be unity.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, n ).
C           Unchanged on exit.
C
C  X      - DOUBLE PRECISION array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the n
C           element right-hand side vector b. On exit, X is overwritten
C           with the solution vector x.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C  -- DO-loops unrolled on 20-November-1986.
C     Peter Mayes, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   TEMP1, TEMP2, TEMP3, TEMP4
      INTEGER            I, INFO, IX, J, JX, KX, N4
      LOGICAL            NOUNIT
C     .. External Subroutines ..
      EXTERNAL           F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC          MAX, MOD
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( .NOT.(UPLO .EQ.'U' .OR. UPLO .EQ.'u').AND.
     $         .NOT.(UPLO .EQ.'L' .OR. UPLO .EQ.'l')      )THEN
         INFO = 1
      ELSE IF( .NOT.(TRANS.EQ.'N' .OR. TRANS.EQ.'n').AND.
     $         .NOT.(TRANS.EQ.'T' .OR. TRANS.EQ.'t').AND.
     $         .NOT.(TRANS.EQ.'C' .OR. TRANS.EQ.'c')      )THEN
         INFO = 2
      ELSE IF( .NOT.(DIAG .EQ.'U' .OR. DIAG .EQ.'u').AND.
     $         .NOT.(DIAG .EQ.'N' .OR. DIAG .EQ.'n')      )THEN
         INFO = 3
      ELSE IF( N.LT.0 )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      END IF
      IF( INFO.NE.0 )THEN
         CALL F06AAZ( 'F06PJF/DTRSV ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( N.EQ.0 )
     $   RETURN
C
      NOUNIT = (DIAG.EQ.'N' .OR. DIAG.EQ.'n')
C
C     Set up the start point in X if the increment is not unity. This
C     will be  ( N - 1 )*INCX  too small for descending loops.
C
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
C
C     Start the operations. In this version the inner loops are all
C     equivalent to AXPY operations.
C
      IF( (TRANS.EQ.'N' .OR. TRANS.EQ.'n') )THEN
C
C        Form  x := inv( A )*x.
C
         IF( (UPLO.EQ.'U' .OR. UPLO.EQ.'u') )THEN
            IF( INCX.EQ.1 )THEN
C**** U n r o l l   t o   D e p t h  4 ********** L O O P  A ***********
               N4 = MOD( N, 4 ) + 1
               DO 20, J = N, N4, -4
C**** 4 by 4 triangle at bottom of block *******************************
                  IF( NOUNIT )
     $               X( J ) = X( J )/A( J, J )
                  X( J - 1 ) = X( J - 1 ) - X( J )*A( J - 1, J )
                  IF( NOUNIT )
     $               X( J - 1 ) = X( J - 1 )/A( J - 1, J - 1 )
                  X( J - 2 ) = X( J - 2 ) - X( J )*A( J - 2, J ) -
     $                         X( J - 1 )*A( J - 2, J - 1 )
                  IF( NOUNIT )
     $               X( J - 2 ) = X( J - 2 )/A( J - 2, J - 2 )
                  X( J - 3 ) = X( J - 3 ) - X( J )*A( J - 3, J ) -
     $                         X( J - 1 )*A( J - 3, J - 1 ) - X( J - 2 )
     $                         *A( J - 3, J - 2 )
                  IF( NOUNIT )
     $               X( J - 3 ) = X( J - 3 )/A( J - 3, J - 3 )
C**** Unrolled main loop ***********************************************
                  TEMP1 = X( J )
                  TEMP2 = X( J - 1 )
                  TEMP3 = X( J - 2 )
                  TEMP4 = X( J - 3 )
                  IF( TEMP1.NE.ZERO.OR.TEMP2.NE.ZERO.OR.TEMP3.NE.
     $                ZERO.OR.TEMP4.NE.ZERO )THEN
                     DO 10, I = J - 4, 1, -1
                        X( I ) = ( ( ( ( X( I ) - TEMP1*A( I, J ) )
     $                           - TEMP2*A( I, J - 1 ) )
     $                           - TEMP3*A( I, J - 2 ) )
     $                           - TEMP4*A( I, J - 3 ) )
   10                CONTINUE
                  END IF
   20          CONTINUE
C**** Left-overs on top left corner ************************************
               IF( N4.GE.2 )THEN
                  IF( NOUNIT )
     $               X( N4 - 1 ) = X( N4 - 1 )/A( N4 - 1, N4 - 1 )
               END IF
               IF( N4.GE.3 )THEN
                  X( N4 - 2 ) = X( N4 - 2 ) - X( N4 - 1 )
     $                          *A( N4 - 2, N4 - 1 )
                  IF( NOUNIT )
     $               X( N4 - 2 ) = X( N4 - 2 )/A( N4 - 2, N4 - 2 )
               END IF
               IF( N4.GE.4 )THEN
                  X( N4 - 3 ) = X( N4 - 3 ) - X( N4 - 1 )
     $                          *A( N4 - 3, N4 - 1 ) - X( N4 - 2 )
     $                          *A( N4 - 3, N4 - 2 )
                  IF( NOUNIT )
     $               X( N4 - 3 ) = X( N4 - 3 )/A( N4 - 3, N4 - 3 )
               END IF
            ELSE
               JX = KX + ( N - 1 )*INCX
C**** U n r o l l   t o   D e p t h  4 ********** L O O P  B ***********
               N4 = MOD( N, 4 ) + 1
               DO 40, J = N, N4, -4
C**** 4 by 4 triangle at bottom of block *******************************
                  IF( NOUNIT )
     $               X( JX ) = X( JX )/A( J, J )
                  X( JX - INCX ) = X( JX - INCX ) - X( JX )
     $                             *A( J - 1, J )
                  IF( NOUNIT )
     $               X( JX - INCX ) = X( JX - INCX )/A( J - 1, J - 1 )
                  X( JX - 2*INCX ) = X( JX - 2*INCX ) - X( JX )
     $                               *A( J - 2, J ) - X( JX - INCX )
     $                               *A( J - 2, J - 1 )
                  IF( NOUNIT )
     $               X( JX - 2*INCX ) = X( JX - 2*INCX )
     $                                  /A( J - 2, J - 2 )
                  X( JX - 3*INCX ) = X( JX - 3*INCX ) - X( JX )
     $                               *A( J - 3, J ) - X( JX - INCX )
     $                               *A( J - 3, J - 1 ) -
     $                               X( JX - 2*INCX )*A( J - 3, J - 2 )
                  IF( NOUNIT )
     $               X( JX - 3*INCX ) = X( JX - 3*INCX )
     $                                  /A( J - 3, J - 3 )
C**** Unrolled main loop ***********************************************
                  TEMP1 = X( JX )
                  TEMP2 = X( JX - INCX )
                  TEMP3 = X( JX - 2*INCX )
                  TEMP4 = X( JX - 3*INCX )
                  IF( TEMP1.NE.ZERO.OR.TEMP2.NE.ZERO.OR.TEMP3.NE.
     $                ZERO.OR.TEMP4.NE.ZERO )THEN
                     IX = JX - 3*INCX
                     DO 30, I = J - 4, 1, -1
                        IX = IX - INCX
                        X( IX ) = ( ( ( ( X( IX ) - TEMP1*A( I, J ) )
     $                            - TEMP2*A( I, J - 1 ) )
     $                            - TEMP3*A( I, J - 2 ) )
     $                            - TEMP4*A( I, J - 3 ) )
   30                CONTINUE
                  END IF
                  JX = JX - 4*INCX
   40          CONTINUE
C**** Left-overs on top left corner ************************************
               IF( N4.GE.2 )THEN
                  IF( NOUNIT )
     $               X( JX ) = X( JX )/A( N4 - 1, N4 - 1 )
               END IF
               IF( N4.GE.3 )THEN
                  X( JX - INCX ) = X( JX - INCX ) - X( JX )
     $                             *A( N4 - 2, N4 - 1 )
                  IF( NOUNIT )
     $               X( JX - INCX ) = X( JX - INCX )/A( N4 - 2, N4 - 2 )
               END IF
               IF( N4.GE.4 )THEN
                  X( JX - 2*INCX ) = X( JX - 2*INCX ) - X( JX )
     $                               *A( N4 - 3, N4 - 1 ) -
     $                               X( JX - INCX )*A( N4 - 3, N4 - 2 )
                  IF( NOUNIT )
     $               X( JX - 2*INCX ) = X( JX - 2*INCX )
     $                                  /A( N4 - 3, N4 - 3 )
               END IF
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
C**** U n r o l l   t o   D e p t h  4 ********** L O O P  C ***********
               N4 = 4*( N/4 )
               DO 60, J = 1, N4, 4
C**** 4 by 4 triangle at bottom of block *******************************
                  IF( NOUNIT )
     $               X( J ) = X( J )/A( J, J )
                  X( J + 1 ) = X( J + 1 ) - X( J )*A( J + 1, J )
                  IF( NOUNIT )
     $               X( J + 1 ) = X( J + 1 )/A( J + 1, J + 1 )
                  X( J + 2 ) = X( J + 2 ) - X( J )*A( J + 2, J ) -
     $                         X( J + 1 )*A( J + 2, J + 1 )
                  IF( NOUNIT )
     $               X( J + 2 ) = X( J + 2 )/A( J + 2, J + 2 )
                  X( J + 3 ) = X( J + 3 ) - X( J )*A( J + 3, J ) -
     $                         X( J + 1 )*A( J + 3, J + 1 ) - X( J + 2 )
     $                         *A( J + 3, J + 2 )
                  IF( NOUNIT )
     $               X( J + 3 ) = X( J + 3 )/A( J + 3, J + 3 )
C**** Unrolled main loop ***********************************************
                  TEMP1 = X( J )
                  TEMP2 = X( J + 1 )
                  TEMP3 = X( J + 2 )
                  TEMP4 = X( J + 3 )
                  IF( TEMP1.NE.ZERO.OR.TEMP2.NE.ZERO.OR.TEMP3.NE.
     $                ZERO.OR.TEMP4.NE.ZERO )THEN
                     DO 50, I = J + 4, N
                        X( I ) = ( ( ( ( X( I ) - TEMP1*A( I, J ) )
     $                           - TEMP2*A( I, J + 1 ) )
     $                           - TEMP3*A( I, J + 2 ) )
     $                           - TEMP4*A( I, J + 3 ) )
   50                CONTINUE
                  END IF
   60          CONTINUE
C**** Left-overs on top left corner ************************************
               IF( N4 + 1.LE.N )THEN
                  IF( NOUNIT )
     $               X( N4 + 1 ) = X( N4 + 1 )/A( N4 + 1, N4 + 1 )
               END IF
               IF( N4 + 2.LE.N )THEN
                  X( N4 + 2 ) = X( N4 + 2 ) - X( N4 + 1 )
     $                          *A( N4 + 2, N4 + 1 )
                  IF( NOUNIT )
     $               X( N4 + 2 ) = X( N4 + 2 )/A( N4 + 2, N4 + 2 )
               END IF
               IF( N4 + 3.LE.N )THEN
                  X( N4 + 3 ) = X( N4 + 3 ) - X( N4 + 1 )
     $                          *A( N4 + 3, N4 + 1 ) - X( N4 + 2 )
     $                          *A( N4 + 3, N4 + 2 )
                  IF( NOUNIT )
     $               X( N4 + 3 ) = X( N4 + 3 )/A( N4 + 3, N4 + 3 )
               END IF
            ELSE
               JX = KX
C**** U n r o l l   t o   D e p t h  4 ********** L O O P  D ***********
               N4 = 4*( N/4 )
               DO 80, J = 1, N4, 4
C**** 4 by 4 triangle at bottom of block *******************************
                  IF( NOUNIT )
     $               X( JX ) = X( JX )/A( J, J )
                  X( JX + INCX ) = X( JX + INCX ) - X( JX )
     $                             *A( J + 1, J )
                  IF( NOUNIT )
     $               X( JX + INCX ) = X( JX + INCX )/A( J + 1, J + 1 )
                  X( JX + 2*INCX ) = X( JX + 2*INCX ) - X( JX )
     $                               *A( J + 2, J ) - X( JX + INCX )
     $                               *A( J + 2, J + 1 )
                  IF( NOUNIT )
     $               X( JX + 2*INCX ) = X( JX + 2*INCX )
     $                                  /A( J + 2, J + 2 )
                  X( JX + 3*INCX ) = X( JX + 3*INCX ) - X( JX )
     $                               *A( J + 3, J ) - X( JX + INCX )
     $                               *A( J + 3, J + 1 ) -
     $                               X( JX + 2*INCX )*A( J + 3, J + 2 )
                  IF( NOUNIT )
     $               X( JX + 3*INCX ) = X( JX + 3*INCX )
     $                                  /A( J + 3, J + 3 )
C**** Unrolled main loop ***********************************************
                  TEMP1 = X( JX )
                  TEMP2 = X( JX + INCX )
                  TEMP3 = X( JX + 2*INCX )
                  TEMP4 = X( JX + 3*INCX )
                  IF( TEMP1.NE.ZERO.OR.TEMP2.NE.ZERO.OR.TEMP3.NE.
     $                ZERO.OR.TEMP4.NE.ZERO )THEN
                     IX = JX + 3*INCX
                     DO 70, I = J + 4, N
                        IX = IX + INCX
                        X( IX ) = ( ( ( ( X( IX ) - TEMP1*A( I, J ) )
     $                            - TEMP2*A( I, J + 1 ) )
     $                            - TEMP3*A( I, J + 2 ) )
     $                            - TEMP4*A( I, J + 3 ) )
   70                CONTINUE
                  END IF
                  JX = JX + 4*INCX
   80          CONTINUE
C**** Left-overs on top left corner ************************************
               IF( N4 + 1.LE.N )THEN
                  IF( NOUNIT )
     $               X( JX ) = X( JX )/A( N4 + 1, N4 + 1 )
               END IF
               IF( N4 + 2.LE.N )THEN
                  X( JX + INCX ) = X( JX + INCX ) - X( JX )
     $                             *A( N4 + 2, N4 + 1 )
                  IF( NOUNIT )
     $               X( JX + INCX ) = X( JX + INCX )/A( N4 + 2, N4 + 2 )
               END IF
               IF( N4 + 3.LE.N )THEN
                  X( JX + 2*INCX ) = X( JX + 2*INCX ) - X( JX )
     $                               *A( N4 + 3, N4 + 1 ) -
     $                               X( JX + INCX )*A( N4 + 3, N4 + 2 )
                  IF( NOUNIT )
     $               X( JX + 2*INCX ) = X( JX + 2*INCX )
     $                                  /A( N4 + 3, N4 + 3 )
               END IF
            END IF
         END IF
      ELSE
C
C        Form  x := inv( A' )*x.
C
         IF( (UPLO.EQ.'U' .OR. UPLO.EQ.'u') )THEN
            IF( INCX.EQ.1 )THEN
C**** U n r o l l   t o   D e p t h  4 ********** L O O P  E ***********
               N4 = 4*( N/4 )
               DO 100, J = 1, N4, 4
C**** 4 by 4 triangle at bottom of block *******************************
                  IF( NOUNIT )
     $               X( J ) = X( J )/A( J, J )
                  X( J + 1 ) = X( J + 1 ) - X( J )*A( J, J + 1 )
                  IF( NOUNIT )
     $               X( J + 1 ) = X( J + 1 )/A( J + 1, J + 1 )
                  X( J + 2 ) = X( J + 2 ) - X( J )*A( J, J + 2 ) -
     $                         X( J + 1 )*A( J + 1, J + 2 )
                  IF( NOUNIT )
     $               X( J + 2 ) = X( J + 2 )/A( J + 2, J + 2 )
                  X( J + 3 ) = X( J + 3 ) - X( J )*A( J, J + 3 ) -
     $                         X( J + 1 )*A( J + 1, J + 3 ) - X( J + 2 )
     $                         *A( J + 2, J + 3 )
                  IF( NOUNIT )
     $               X( J + 3 ) = X( J + 3 )/A( J + 3, J + 3 )
C**** Unrolled main loop ***********************************************
                  TEMP1 = X( J )
                  TEMP2 = X( J + 1 )
                  TEMP3 = X( J + 2 )
                  TEMP4 = X( J + 3 )
                  IF( TEMP1.NE.ZERO.OR.TEMP2.NE.ZERO.OR.TEMP3.NE.
     $                ZERO.OR.TEMP4.NE.ZERO )THEN
                     DO 90, I = J + 4, N
                        X( I ) = ( ( ( ( X( I ) - TEMP1*A( J, I ) )
     $                           - TEMP2*A( J + 1, I ) )
     $                           - TEMP3*A( J + 2, I ) )
     $                           - TEMP4*A( J + 3, I ) )
   90                CONTINUE
                  END IF
  100          CONTINUE
C**** Left-overs on top left corner ************************************
               IF( N4 + 1.LE.N )THEN
                  IF( NOUNIT )
     $               X( N4 + 1 ) = X( N4 + 1 )/A( N4 + 1, N4 + 1 )
               END IF
               IF( N4 + 2.LE.N )THEN
                  X( N4 + 2 ) = X( N4 + 2 ) - X( N4 + 1 )
     $                          *A( N4 + 1, N4 + 2 )
                  IF( NOUNIT )
     $               X( N4 + 2 ) = X( N4 + 2 )/A( N4 + 2, N4 + 2 )
               END IF
               IF( N4 + 3.LE.N )THEN
                  X( N4 + 3 ) = X( N4 + 3 ) - X( N4 + 1 )
     $                          *A( N4 + 1, N4 + 3 ) - X( N4 + 2 )
     $                          *A( N4 + 2, N4 + 3 )
                  IF( NOUNIT )
     $               X( N4 + 3 ) = X( N4 + 3 )/A( N4 + 3, N4 + 3 )
               END IF
            ELSE
               JX = KX
C**** U n r o l l   t o   D e p t h  4 ********** L O O P  D ***********
               N4 = 4*( N/4 )
               DO 120, J = 1, N4, 4
C**** 4 by 4 triangle at bottom of block *******************************
                  IF( NOUNIT )
     $               X( JX ) = X( JX )/A( J, J )
                  X( JX + INCX ) = X( JX + INCX ) - X( JX )
     $                             *A( J, J + 1 )
                  IF( NOUNIT )
     $               X( JX + INCX ) = X( JX + INCX )/A( J + 1, J + 1 )
                  X( JX + 2*INCX ) = X( JX + 2*INCX ) - X( JX )
     $                               *A( J, J + 2 ) - X( JX + INCX )
     $                               *A( J + 1, J + 2 )
                  IF( NOUNIT )
     $               X( JX + 2*INCX ) = X( JX + 2*INCX )
     $                                  /A( J + 2, J + 2 )
                  X( JX + 3*INCX ) = X( JX + 3*INCX ) - X( JX )
     $                               *A( J, J + 3 ) - X( JX + INCX )
     $                               *A( J + 1, J + 3 ) -
     $                               X( JX + 2*INCX )*A( J + 2, J + 3 )
                  IF( NOUNIT )
     $               X( JX + 3*INCX ) = X( JX + 3*INCX )
     $                                  /A( J + 3, J + 3 )
C**** Unrolled main loop ***********************************************
                  TEMP1 = X( JX )
                  TEMP2 = X( JX + INCX )
                  TEMP3 = X( JX + 2*INCX )
                  TEMP4 = X( JX + 3*INCX )
                  IF( TEMP1.NE.ZERO.OR.TEMP2.NE.ZERO.OR.TEMP3.NE.
     $                ZERO.OR.TEMP4.NE.ZERO )THEN
                     IX = JX + 3*INCX
                     DO 110, I = J + 4, N
                        IX = IX + INCX
                        X( IX ) = ( ( ( ( X( IX ) - TEMP1*A( J, I ) )
     $                            - TEMP2*A( J + 1, I ) )
     $                            - TEMP3*A( J + 2, I ) )
     $                            - TEMP4*A( J + 3, I ) )
  110                CONTINUE
                  END IF
                  JX = JX + 4*INCX
  120          CONTINUE
C**** Left-overs on top left corner ************************************
               IF( N4 + 1.LE.N )THEN
                  IF( NOUNIT )
     $               X( JX ) = X( JX )/A( N4 + 1, N4 + 1 )
               END IF
               IF( N4 + 2.LE.N )THEN
                  X( JX + INCX ) = X( JX + INCX ) - X( JX )
     $                             *A( N4 + 1, N4 + 2 )
                  IF( NOUNIT )
     $               X( JX + INCX ) = X( JX + INCX )/A( N4 + 2, N4 + 2 )
               END IF
               IF( N4 + 3.LE.N )THEN
                  X( JX + 2*INCX ) = X( JX + 2*INCX ) - X( JX )
     $                               *A( N4 + 1, N4 + 3 ) -
     $                               X( JX + INCX )*A( N4 + 2, N4 + 3 )
                  IF( NOUNIT )
     $               X( JX + 2*INCX ) = X( JX + 2*INCX )
     $                                  /A( N4 + 3, N4 + 3 )
               END IF
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
C**** U n r o l l   t o   D e p t h  4 ********** L O O P  G ***********
               N4 = MOD( N, 4 ) + 1
               DO 140, J = N, N4, -4
C**** 4 by 4 triangle at bottom of block *******************************
                  IF( NOUNIT )
     $               X( J ) = X( J )/A( J, J )
                  X( J - 1 ) = X( J - 1 ) - X( J )*A( J, J - 1 )
                  IF( NOUNIT )
     $               X( J - 1 ) = X( J - 1 )/A( J - 1, J - 1 )
                  X( J - 2 ) = X( J - 2 ) - X( J )*A( J, J - 2 ) -
     $                         X( J - 1 )*A( J - 1, J - 2 )
                  IF( NOUNIT )
     $               X( J - 2 ) = X( J - 2 )/A( J - 2, J - 2 )
                  X( J - 3 ) = X( J - 3 ) - X( J )*A( J, J - 3 ) -
     $                         X( J - 1 )*A( J - 1, J - 3 ) - X( J - 2 )
     $                         *A( J - 2, J - 3 )
                  IF( NOUNIT )
     $               X( J - 3 ) = X( J - 3 )/A( J - 3, J - 3 )
C**** Unrolled main loop ***********************************************
                  TEMP1 = X( J )
                  TEMP2 = X( J - 1 )
                  TEMP3 = X( J - 2 )
                  TEMP4 = X( J - 3 )
                  IF( TEMP1.NE.ZERO.OR.TEMP2.NE.ZERO.OR.TEMP3.NE.
     $                ZERO.OR.TEMP4.NE.ZERO )THEN
                     DO 130, I = J - 4, 1, -1
                        X( I ) = ( ( ( ( X( I ) - TEMP1*A( J, I ) )
     $                           - TEMP2*A( J - 1, I ) )
     $                           - TEMP3*A( J - 2, I ) )
     $                           - TEMP4*A( J - 3, I ) )
  130                CONTINUE
                  END IF
  140          CONTINUE
C**** Left-overs on top left corner ************************************
               IF( N4.GE.2 )THEN
                  IF( NOUNIT )
     $               X( N4 - 1 ) = X( N4 - 1 )/A( N4 - 1, N4 - 1 )
               END IF
               IF( N4.GE.3 )THEN
                  X( N4 - 2 ) = X( N4 - 2 ) - X( N4 - 1 )
     $                          *A( N4 - 1, N4 - 2 )
                  IF( NOUNIT )
     $               X( N4 - 2 ) = X( N4 - 2 )/A( N4 - 2, N4 - 2 )
               END IF
               IF( N4.GE.4 )THEN
                  X( N4 - 3 ) = X( N4 - 3 ) - X( N4 - 1 )
     $                          *A( N4 - 1, N4 - 3 ) - X( N4 - 2 )
     $                          *A( N4 - 2, N4 - 3 )
                  IF( NOUNIT )
     $               X( N4 - 3 ) = X( N4 - 3 )/A( N4 - 3, N4 - 3 )
               END IF
            ELSE
               JX = KX + ( N - 1 )*INCX
C**** U n r o l l   t o   D e p t h  4 ********** L O O P  B ***********
               N4 = MOD( N, 4 ) + 1
               DO 160, J = N, N4, -4
C**** 4 by 4 triangle at bottom of block *******************************
                  IF( NOUNIT )
     $               X( JX ) = X( JX )/A( J, J )
                  X( JX - INCX ) = X( JX - INCX ) - X( JX )
     $                             *A( J, J - 1 )
                  IF( NOUNIT )
     $               X( JX - INCX ) = X( JX - INCX )/A( J - 1, J - 1 )
                  X( JX - 2*INCX ) = X( JX - 2*INCX ) - X( JX )
     $                               *A( J, J - 2 ) - X( JX - INCX )
     $                               *A( J - 1, J - 2 )
                  IF( NOUNIT )
     $               X( JX - 2*INCX ) = X( JX - 2*INCX )
     $                                  /A( J - 2, J - 2 )
                  X( JX - 3*INCX ) = X( JX - 3*INCX ) - X( JX )
     $                               *A( J, J - 3 ) - X( JX - INCX )
     $                               *A( J - 1, J - 3 ) -
     $                               X( JX - 2*INCX )*A( J - 2, J - 3 )
                  IF( NOUNIT )
     $               X( JX - 3*INCX ) = X( JX - 3*INCX )
     $                                  /A( J - 3, J - 3 )
C**** Unrolled main loop ***********************************************
                  TEMP1 = X( JX )
                  TEMP2 = X( JX - INCX )
                  TEMP3 = X( JX - 2*INCX )
                  TEMP4 = X( JX - 3*INCX )
                  IF( TEMP1.NE.ZERO.OR.TEMP2.NE.ZERO.OR.TEMP3.NE.
     $                ZERO.OR.TEMP4.NE.ZERO )THEN
                     IX = JX - 3*INCX
                     DO 150, I = J - 4, 1, -1
                        IX = IX - INCX
                        X( IX ) = ( ( ( ( X( IX ) - TEMP1*A( J, I ) )
     $                            - TEMP2*A( J - 1, I ) )
     $                            - TEMP3*A( J - 2, I ) )
     $                            - TEMP4*A( J - 3, I ) )
  150                CONTINUE
                  END IF
                  JX = JX - 4*INCX
  160          CONTINUE
C**** Left-overs on top left corner ************************************
               IF( N4.GE.2 )THEN
                  IF( NOUNIT )
     $               X( JX ) = X( JX )/A( N4 - 1, N4 - 1 )
               END IF
               IF( N4.GE.3 )THEN
                  X( JX - INCX ) = X( JX - INCX ) - X( JX )
     $                             *A( N4 - 1, N4 - 2 )
                  IF( NOUNIT )
     $               X( JX - INCX ) = X( JX - INCX )/A( N4 - 2, N4 - 2 )
               END IF
               IF( N4.GE.4 )THEN
                  X( JX - 2*INCX ) = X( JX - 2*INCX ) - X( JX )
     $                               *A( N4 - 1, N4 - 3 ) -
     $                               X( JX - INCX )*A( N4 - 2, N4 - 3 )
                  IF( NOUNIT )
     $               X( JX - 2*INCX ) = X( JX - 2*INCX )
     $                                  /A( N4 - 3, N4 - 3 )
               END IF
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06PJF (DTRSV ).
C
      END
      DOUBLE PRECISION FUNCTION X02AJF()
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     RETURNS  (1/2)*B**(1-P)  IF ROUNDS IS .TRUE.
C     RETURNS  B**(1-P)  OTHERWISE
C
      DOUBLE PRECISION X02CON
      DATA X02CON /4.6d-24/
C     .. Executable Statements ..
      X02AJF = X02CON
      RETURN
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
      X=0.0D0
      X01AAF = X+3.14159265358979323846264338328D0
      RETURN
      END

      SUBROUTINE CONVER(NR,FM,GRE,GIM)
      INTEGER J,I,NR,K
      DOUBLE PRECISION KOFF,GIM(36,NR),GRE(36,NR),FM,FRE(1670)
      DO J=1,36
         DO I=1,NR
            GIM(J,I)=(GIM(J,I)/1.0D26)
            GRE(J,I)=(GRE(J,I)/1.0D26)
         ENDDO
      ENDDO
      DO J=1,36
         DO I=1,NR/2
            GIM(J,I)=GIM(J,I*2)
         ENDDO
         DO I=NR/2+1,NR
            GIM(J,I)=0.0D0
         ENDDO
      ENDDO
      DO I=1,NR
         FRE(I)=DBLE(I)/DBLE(NR)*6.4D13
      ENDDO
      PRINT*,'OK'
      OPEN(2,FILE='DATA1.TXT')
      DO I=1,NR
         WRITE(2,2)FRE(I),GIM(1,I),GRE(1,I)
      ENDDO
      CLOSE(2)
      OPEN(3,FILE='DATA2.TXT')
      DO I=1,NR
         WRITE(3,2)FRE(I),GIM(2,I),GRE(2,I)
      ENDDO
      CLOSE(3)
      OPEN(4,FILE='DATA3.TXT')
      DO I=1,nr
         WRITE(4,2)FRE(I),GIM(3,I),GRE(3,I)
      ENDDO
      CLOSE(4)
      OPEN(5,FILE='DATA4.TXT')
      DO I=1,nr
         WRITE(5,2)FRE(I),GIM(4,I),GRE(4,I)
      ENDDO
      CLOSE(5)
      OPEN(6,FILE='DATA5.TXT')
      DO I=1,nr
         WRITE(6,2)FRE(I),GIM(5,I),GRE(5,I)
      ENDDO
      CLOSE(6)
      OPEN(7,FILE='DATA6.TXT')
      DO I=1,nr
         WRITE(7,2)FRE(I),GIM(6,I),GRE(6,I)
      ENDDO
      CLOSE(7)
      OPEN(8,FILE='DATA7.TXT')
      DO I=1,nr
         WRITE(8,2)FRE(I),GIM(7,I),GRE(7,I)
      ENDDO
      CLOSE(8)
      OPEN(9,FILE='DATA8.TXT')
      DO I=1,nr
         WRITE(9,2)FRE(I),GIM(8,I),GRE(8,I)
      ENDDO
      CLOSE(9)
      OPEN(10,FILE='DATA9.TXT')
      DO I=1,nr
         WRITE(10,2)FRE(I),GIM(9,I),GRE(9,I)
      ENDDO
      CLOSE(10)
      OPEN(11,FILE='DATA10.TXT')
      DO I=1,nr
         WRITE(11,2)FRE(I),GIM(10,I),GRE(10,I)
      ENDDO
      CLOSE(11)
      OPEN(12,FILE='DATA11.TXT')
      DO I=1,nr
         WRITE(12,2)FRE(I),GIM(11,I),GRE(11,I)
      ENDDO
      CLOSE(12)
      OPEN(13,FILE='DATA12.TXT')
      DO I=1,nr
         WRITE(13,2)FRE(I),GIM(12,I),GRE(12,I)
      ENDDO
      CLOSE(13)
      OPEN(14,FILE='DATA13.TXT')
      DO I=1,NR
         WRITE(14,2)FRE(I),GIM(13,I),GRE(13,I)
      ENDDO
      CLOSE(14)
      OPEN(15,FILE='DATA14.TXT')
      DO I=1,NR
         WRITE(15,2)FRE(I),GIM(14,I),GRE(14,I)
      ENDDO
      CLOSE(15)
      OPEN(16,FILE='DATA15.TXT')
      DO I=1,nr
         WRITE(16,2)FRE(I),GIM(15,I),GRE(15,I)
      ENDDO
      CLOSE(16)
      OPEN(17,FILE='DATA16.TXT')
      DO I=1,nr
         WRITE(17,2)FRE(I),GIM(16,I),GRE(16,I)
      ENDDO
      CLOSE(17)
      OPEN(18,FILE='DATA17.TXT')
      DO I=1,nr
         WRITE(18,2)FRE(I),GIM(17,I),GRE(17,I)
      ENDDO
      CLOSE(18)
      OPEN(19,FILE='DATA18.TXT')
      DO I=1,nr
         WRITE(19,2)FRE(I),GIM(18,I),GRE(18,I)
      ENDDO
      CLOSE(19)
      OPEN(20,FILE='DATA19.TXT')
      DO I=1,nr
         WRITE(20,2)FRE(I),GIM(19,I),GRE(19,I)
      ENDDO
      CLOSE(20)
      OPEN(21,FILE='DATA20.TXT')
      DO I=1,nr
         WRITE(21,2)FRE(I),GIM(20,I),GRE(20,I)
      ENDDO
      CLOSE(21)
      OPEN(22,FILE='DATA21.TXT')
      DO I=1,nr
         WRITE(22,2)FRE(I),GIM(21,I),GRE(21,I)
      ENDDO
      CLOSE(22)
      OPEN(23,FILE='DATA22.TXT')
      DO I=1,nr
         WRITE(23,2)FRE(I),GIM(22,I),GRE(22,I)
      ENDDO
      CLOSE(23)
      OPEN(24,FILE='DATA23.TXT')
      DO I=1,nr
         WRITE(24,2)FRE(I),GIM(23,I),GRE(23,I)
      ENDDO
      CLOSE(24)
      OPEN(25,FILE='DATA24.TXT')
      DO I=1,nr
         WRITE(25,2)FRE(I),GIM(24,I),GRE(24,I)
      ENDDO
      CLOSE(25)
      OPEN(26,FILE='DATA25.TXT')
      DO I=1,nr
         WRITE(26,2)FRE(I),GIM(25,I),GRE(25,I)
      ENDDO
      CLOSE(26)
      OPEN(27,FILE='DATA26.TXT')
      DO I=1,nr
         WRITE(27,2)FRE(I),GIM(26,I),GRE(26,I)
      ENDDO
      CLOSE(27)
      OPEN(28,FILE='DATA27.TXT')
      DO I=1,nr
         WRITE(28,2)FRE(I),GIM(27,I),GRE(27,I)
      ENDDO
      CLOSE(28)
      OPEN(29,FILE='DATA28.TXT')
      DO I=1,nr
         WRITE(29,2)FRE(I),GIM(28,I),GRE(28,I)
      ENDDO
      CLOSE(29)
      OPEN(30,FILE='DATA29.TXT')
      DO I=1,nr
         WRITE(30,2)FRE(I),GIM(29,I),GRE(29,I)
      ENDDO
      CLOSE(30)
      OPEN(31,FILE='DATA30.TXT')
      DO I=1,nr
         WRITE(31,2)FRE(I),GIM(30,I),GRE(30,I)
      ENDDO
      CLOSE(31)
      OPEN(32,FILE='DATA31.TXT')
      DO I=1,nr
         WRITE(32,2)FRE(I),GIM(31,I),GRE(31,I)
      ENDDO
      CLOSE(32)
      OPEN(33,FILE='DATA32.TXT')
      DO I=1,nr
         WRITE(33,2)FRE(I),GIM(32,I),GRE(32,I)
      ENDDO
      CLOSE(33)
      OPEN(34,FILE='DATA33.TXT')
      DO I=1,nr
         WRITE(34,2)FRE(I),GIM(33,I),GRE(33,I)
      ENDDO
      CLOSE(34)
      OPEN(35,FILE='DATA34.TXT')
      DO I=1,nr
         WRITE(35,2)FRE(I),GIM(34,I),GRE(34,I)
      ENDDO
      CLOSE(35)
      OPEN(36,FILE='DATA35.TXT')
      DO I=1,nr
         WRITE(36,2)FRE(I),GIM(35,I),GRE(35,I)
      ENDDO
      CLOSE(36)
      OPEN(37,FILE='DATA36.TXT')
      DO I=1,nr
         WRITE(37,2)FRE(I),GIM(36,I),GRE(36,I)
      ENDDO
      CLOSE(37)
      OPEN(38,FILE='DATA0.TXT')
      DO I=1,NR
         WRITE(38,2)FRE(I),GIM(37,I),GRE(37,I)
      ENDDO
      CLOSE(38)

 2    FORMAT(3E16.8)
      END


