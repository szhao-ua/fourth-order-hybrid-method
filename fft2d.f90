module mymodule
    implicit none
  
    contains

    
    SUBROUTINE FFT4DD(F,M,N,HX,HY)
      ! TO CALCULATE THE 4th ORDER PROBLEM
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER :: M,N
      REAL*8 :: F(M,N),F1(M,N),F2(M,N)
      REAL*8 :: X(M+1),Y(N+1)
      REAL*8 :: LAMBDAX(M+1),LAMBDAY(N+1),ERROR(M,N)
      REAL*8 :: HX,HY,LMAX,P
      INTEGER :: I,J
      
      F1(M,N)=0.0D0
      F2(M,N)=0.0D0
      
      PI=ACOS(-1.0D0)
      
      DO I=1,M
        Y=0.0D0
       DO J=1,N
          Y(J+1)=F(I,J)
       ENDDO
       CALL dsinft(y,n+1)
       DO J=1,N
       F1(I,J)=Y(J+1)
       END DO
      ENDDO
      
      DO J=1,N
        X=0.0D0
        DO I=1,M
          X(I+1)=F1(I,J)
        ENDDO
      
      CALL dsinft(X,M+1)
      DO I=1,M
      F2(I,J)=X(I+1)
      ENDDO
      ENDDO
      
      DO I=1,M
        LAMBDAX(I+1) =-1.0D0/3.0/HX/HX*(COS(I*PI/(M+1))-1.0D0)*(COS(I*PI/(M+1))-7.0D0)
      ENDDO
      
      DO I=1,N
      LAMBDAY(I+1)=-1.0D0/3.0/HY/HY*(COS(I*PI/(N+1))-1.0D0)*(COS(I*PI/(N+1))-7.0D0)
      ENDDO
      
      DO I=1,M
        DO J=1,N
         F2(I,J)=F2(I,J)*4.0D0/(M+1)/(N+1)*1.0D0/(LAMBDAX(I+1)+LAMBDAY(J+1)+P)
        ENDDO
      ENDDO
      
      DO J=1,N
        X=0.0D0
       DO I=1,M
          X(I+1)=F2(I,J)
       ENDDO
       CALL dsinft(X,M+1)
       DO I=1,M
       F1(I,J)=X(I+1)
       END DO
      ENDDO
      
      DO I=1,M
        Y=0.0D0
        DO J=1,N
          Y(J+1)=F1(I,J)
        ENDDO
      
      CALL dsinft(Y,N+1)
      DO J=1,N
      F(I,J)=Y(J+1)
      ENDDO
      ENDDO
      
      
      END SUBROUTINE

          
end module mymodule