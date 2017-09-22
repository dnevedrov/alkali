      PROGRAM SMOOTHING
      INTEGER N
      PARAMETER (N=3000)
      DOUBLE PRECISION X(N),Y1(N),Y2(N),Y3(N)
      INTEGER J
      OPEN(1,FILE='111')
      DO J=1,N
            READ(1,2)X(J),Y1(J),Y2(J),Y3(J)
      ENDDO
      CLOSE(1)

      OPEN(2,FILE='1kclna111.dat')
      DO J=1,N
         WRITE(2,1)X(J),Y1(J)/Y3(J)
      ENDDO
      CLOSE(2)

      OPEN(3,FILE='2kclna111.dat')
      DO J=1,N
         WRITE(3,1)X(J),Y2(J)*1.0D8
      ENDDO
      CLOSE(3)

      OPEN(4,FILE='3kclna111.dat')
      DO J=1,N
         WRITE(4,1)X(J),Y3(J)*1.0D12
      ENDDO
      CLOSE(4)

 1    FORMAT(2E16.8)
 2    FORMAT(4E16.8)
      END