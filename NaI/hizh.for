      PROGRAM HIZH
C     ===========================================================
      DOUBLE PRECISION ONE, TWO, THR
      PARAMETER (ONE=1.0D0, TWO=2.0D0, THR=ONE+TWO)
      DOUBLE PRECISION WM, V2, W, RHO, RHO11, RHO12, RHO22,
     *     BETA, A2, A0, PHI, ALPHA, G11, G22, G12, ETA1, ETA2,
     *     K4, B, OMEGAL1, OMEGAL2, OMEGAL, V4, SINPHI, COSPHI, 
     *     SIN2A, A, X
C     ===========================================================
      WM=TW0*DSQRT(V2)
      RHO=(W/WM-DSQRT(W**2/WM**2-ONE))**2
      RHO11=ONE-RHO
      RHO22=ONE-RHO**3
      RHO12=RHO*(ONE-RHO)
      BETA=A2/A0
      B=ONE/TWO**2*(ONE+BETA)**2
      SINPHI=(DSQRT(ONE+B**2)-ONE)/(B**2+(DSQRT(ONE+B**2)-ONE)**2)
      COSPHI=DSQRT(ONE-SINPHI**2)
      G11=RHO11*COSPHI**2+RHO22*SINPHI**2+RHO12*TWO*SINPHI*COSPHI
      G22=RHO11*SINPHI**2+RHO22*COSPHI**2+RHO12*TWO*SINPHI*COSPHI
      G12=RHO12*(COSPHI**2-SINPHI**2)+
     *     ONE/TWO*(RHO11-RGO22)**TWO*SINPHI*COSPHI
      ETA1=TWO*THR*V4*(ONE+B+DSQRT(ONE-B**2))
      ETA1=TWO*THR*V4*(ONE+B-DSQRT(ONE-B**2))
      SIN2A=TWO*DSQRT(ETA1*ETA2)*DABS(G12)/
     *     DSQRT((ETA1*G11-ETA2*G22)**2+TWO**2*ETA1*ETA2*G12**2)
      A=DASIN(SIN2A)/TWO
C     ===========================================================     
      X=K4*A0**2*(ONE+B+DSQRT(ONE-B**2))*
     *     (G11*DCOS(A)**2+G22*DSIN(A)**2*ETA2/ETA1+
     *     G12*DSQRT(ETA2/ETA1)*SIN2A)
C     ===========================================================
      OMEGAL1=(ONE/TWO)+DSQRT(ONE+TWO**2*X**2)/TWO
      OMEGAL1=(ONE/TWO)-DSQRT(ONE+TWO**2*X**2)/TWO
      
C     ===========================================================
      END
      
