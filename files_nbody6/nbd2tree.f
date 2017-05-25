!     12/05/2002
      PARAMETER (nmax=32000)
      REAL AS(nmax),BODYS(nmax)
      REAL XS(3,nmax),VS(3,nmax)
      real rlagr(11),avmass(11)
      INTEGER NAME(nmax)
      character*4 str

      open (10,file='infile')
      read (10,*) i,t !,t,i,i
      read (10,*) NPAR,NFIX, i,i,i,i
      read (10,*) ai,ai,ai,ai, DELTAT, ai,ai,ai,ai
      close(10)
      OPEN (UNIT=3,STATUS='OLD',FORM='UNFORMATTED',FILE='OUT3')
    
      ilabel=0  
   1  READ (3,end=9)  NTOT, MODEL, NRUN, NK
        READ (3)  (AS(K),K=1,NK), (BODYS(J),J=1,NTOT),
     &           ((XS(K,J),K=1,3),J=1,NTOT), ((VS(K,J),K=1,3),J=1,NTOT),
     &           (NAME(J),J=1,NTOT)

        print *,as(1)
        Call i2str(ilabel,4,str)
        OPEN (4,FILE='c_'//str//'.dat')
        write (4,fmt='(A,1x,F7.1)')"# t(Nbody units)=",as(1)
c       write (4,fmt='(A)')"# x y z vx vy vz m label"
        write (4,fmt='(A)')"# x y z vx vy vz m"
        DO I=1,NTOT
c         WRITE (4,100)(XS(K,I),K=1,3),(VS(K,I),K=1,3),BODYS(I)
        if ((BODYS(I).GT.0.) .AND. (NAME(I).LT.(NPAR+1)) 
     &           .AND. (NAME(I).GT.0)) THEN

          WRITE (4,101)(XS(K,I),K=1,3),(VS(K,I),K=1,3),BODYS(I),NAME(I)
          endif
        END DO
        CLOSE(4)
        ilabel=ilabel+1
        go to 1
 9    CLOSE(3)

      open (7,file='fort.7')
!      open (7,file='fort.10')
      open (10,file='rlag6.dat')
 11         READ (7,110,end=99) TIME, (RLAGR(K),K=1,11),
     &           (AVMASS(K),K=1,11),ai,ai
            write (10,'(12(1x,g11.4))')time,(RLAGR(K),K=1,11)
!           READ (7,110,end=99) TIME, (RLAGR(K),K=1,11),
!     &           (AVMASS(K),K=1,11),ai,ai
!           READ (7,110,end=99) TIME, (RLAGR(K),K=1,11),
!     &           (AVMASS(K),K=1,11),ai,ai
!           READ (7,110,end=99) TIME, (RLAGR(K),K=1,11),
!     &           (AVMASS(K),K=1,11),ai,ai
            go to 11
 99        close(7)
           close(10)

 110         FORMAT(1X,1P,25(1X,D15.7))

 100  FORMAT (7(1X,G15.8))
 101  FORMAT (7(1X,G15.8),1X,I6)
      END

      Subroutine i2str(i,len,stringa)
        character*1 stringa(len)
        integer i,len
!       convert I into a string stringa with len leading 0s

        integer j
        do j=len,1,-1
c         i2str=trim(i2str)//char(48+I/10**(J-1)-I/10**J*10)
          stringa(len+1-j)=char(48+I/10**(J-1)-I/10**J*10)
        end do
      Return
      
      End


