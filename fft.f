      
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C   Calculates the Fourier transform of a set of n real-valued data points.
C   Replaces this data (which is stored in array data(1:n)) by the positive
C   frequency half of its complex Fourier transform. 
C   The real-valued first and last components of the complex transform are
C   returned as elements data(1) and data(2), respectively.
C   n must be a power of 2. This routine also calculates the inverse transform
C   of a complex data array if it is the transform of real data.
C   (Result in this case must be multiplied by 2/n.)
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE drealft(data,n,isign)
      IMPLICIT NONE
      INTEGER, INTENT(IN):: n,isign
      INTEGER i,i1,i2,i3,i4,n2p3
      REAL*8, DIMENSION(n), INTENT(INOUT) :: data
      REAL*8 c1,c2,h1r,h1i,h2r,h2i,wr,wi,wpr,wpi,wtemp,wrs,wis
      REAL*8 theta

C   USES dfour1
      
      theta=3.141592653589793d0/(n/2.0d0)   !Initialize the recurrence.
      c1=0.5d0
      if (isign.eq.1) then
        c2=-0.5d0
        call dfour1(data,n/2,+1)      !The forward transform is here.
      else
        c2=0.5d0         !Otherwise set up for an inverse transform.
        theta=-theta
      endif

      wpr=-2.0d0*sin(0.5d0*theta)**2
      wpi=sin(theta)
      wr=1.0d0+wpr
      wi=wpi

      n2p3=n+3
      do i=2,n/4      ! Case i=1 done separately below.
        i1=2*i-1
        i2=i1+1
        i3=n2p3-i2 
        i4=i3+1

        wrs=wr
        wis=wi
        h1r=c1*(data(i1)+data(i3))   ! The two separate transforms are separated out of
        h1i=c1*(data(i2)-data(i4))   ! data.
        h2r=-c2*(data(i2)+data(i4))
        h2i=c2*(data(i1)-data(i3))
        data(i1)=h1r+wrs*h2r-wis*h2i   !Here they are recombined to form the true 
        data(i2)=h1i+wrs*h2i+wis*h2r   !transform of the original real data.
        data(i3)=h1r-wrs*h2r+wis*h2i
        data(i4)=-h1i+wrs*h2i+wis*h2r
        wtemp=wr         !The recurrence.
        wr=wr*wpr-wi*wpi+wr
        wi=wi*wpr+wtemp*wpi+wi
      enddo

      if (isign.eq.1) then
        h1r=data(1)
        data(1)=h1r+data(2)
        data(2)=h1r-data(2)      
C   Squeeze the first and last data together to
C   get them all within the original array
      else
        h1r=data(1)
        data(1)=c1*(h1r+data(2))
        data(2)=c1*(h1r-data(2))
        call dfour1(data,n/2,-1)   ! This is the inverse transform for the case isign=-1.
      endif

      return
      END

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C   Calculates the sine transform of a set of n real-valued data points stored
C   in array y(1:n). The number n must be a power of 2. On exit y is replaced 
C   by its transform. This program, without changes, also calculates the inverse
C   sine transform, but in this case the output array should be multiplied by 2/n.
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE dsinft(y,n)
      IMPLICIT NONE
      REAL*8 theta,wr,wi,wpr,wpi,wtemp,y1,y2,sum
      REAL*8, DIMENSION(n), INTENT(INOUT) :: y
      INTEGER, INTENT(IN) :: n
      INTEGER j
C   USES drealft

      theta=3.141592653589793d0/dble(n)      !Initialize the recurrence.
      wr=1.0d0
      wi=0.0d0
      wpr=-2.0d0*sin(0.5d0*theta)**2
      wpi=sin(theta)
      y(1)=0.0d0

      do j=1,n/2
        wtemp=wr
        wr=wr*wpr-wi*wpi+wr      !Calculate the sine for the auxiliary array.
        wi=wi*wpr+wtemp*wpi+wi      !The cosine is needed to continue the recurrence.
        y1=wi*(y(j+1)+y(n-j+1))      !Construct the auxiliary array.
        y2=0.5d0*(y(j+1)-y(n-j+1))
        y(j+1)=y1+y2         !Terms j and N - j are related.
        y(n-j+1)=y1-y2
      enddo

      call drealft(y,n,+1)      !Transform the auxiliary array.
      sum=0.0
      y(1)=0.5d0*y(1)         !Initialize the sum used for odd terms below.
      y(2)=0.0d0

      do j=1,n-1,2
        sum=sum+y(j)
        y(j)=y(j+1)      !Even terms in the transform are determined directly.
        y(j+1)=sum      !Odd terms are determined by this running sum.
      enddo

      return
      END
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C   FAST FOURIER TRANSFORM
C   Replaces data(1:2*nn) by its discrete Fourier transform, if
C   isign is input as 1; or replaces data(1:2*nn) by nn times its
C   inverse discrete Fourier transform, if isign is input as -1.
C   data is a complex array of length nn or, equivalently, a real
C   array of length 2*nn. nn MUST be an integer power of 2 (this
C   is not checked for!).
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       SUBROUTINE DFOUR1(DATA,NN,ISIGN)
       IMPLICIT NONE
        INTEGER, INTENT(IN) :: NN,ISIGN
        REAL*8, DIMENSION(2*NN), INTENT(INOUT) :: DATA
        REAL*8 WR,WI,WPR,WPI,WTEMP,THETA
        REAL*8 TEMPR,TEMPI
        REAL*8 WRS, WIS
        REAL*8 PI
        INTEGER N,M,MMAX,ISTEP
        REAL*8, DIMENSION(NN) :: X,Y,LAMBDAX,LAMBDAY
        INTEGER I,J



       PI=ACOS(-1.D0)
       N=2*NN
       J=1
       
C   This is the bit-reversal section of the routine.
       DO I=1,N,2
          IF(J.GT.I) THEN
             TEMPR=DATA(J)
C   Exchange the two complex numbers.

             TEMPI=DATA(J+1)
             DATA(J)=DATA(I)
             DATA(J+1)=DATA(I+1)
             DATA(I)=TEMPR
             DATA(I+1)=TEMPI
          END IF
          M=N/2
1000      IF((M.GE.2).AND.(J.GT.M)) THEN
             J=J-M
             M=M/2
             GO TO  1000
          END IF  
          J=J+M
       END DO
       MMAX=2

C   Here begins the Danielson-Lanczos section of the routine.
1001   IF(N.GT.MMAX) THEN

C   Outer loop executed log_2(nn) times.
          ISTEP=2*MMAX
          THETA=2.D0*PI/(ISIGN*MMAX)
          WPR=-2.D0*SIN(0.5D0*THETA)**2
          WPI=SIN(THETA)
          WR=1.D0
          WI=0.D0

C   Here are the two nested inner loops.
          DO M=1,MMAX,2
             DO I=M,N,ISTEP
                J=I+MMAX

C   This is the Danielson-Lanczos formula:
                TEMPR=WR*DATA(J)-WI*DATA(J+1)
                TEMPI=WR*DATA(J+1)+WI*DATA(J)
                DATA(J)  =DATA(I)  -TEMPR
                DATA(J+1)=DATA(I+1)-TEMPI
                DATA(I)  =DATA(I)  +TEMPR
                DATA(I+1)=DATA(I+1)+TEMPI
             END DO 
             WTEMP=WR

C   Trigonometric recurrence.
             WR=WR*WPR-WI*WPI+WR
             WI=WI*WPR+WTEMP*WPI+WI
          END DO
          MMAX=ISTEP
       GOTO 1001
C   Not yet done.

       END IF

       return
       END          
