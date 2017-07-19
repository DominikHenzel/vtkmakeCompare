c
c fft/p2s-fornberg.f
c (C) 2001 K. Fukagata 
c
      SUBROUTINE P2S(FP,FPT,NX,NY,NZ)
      implicit none
      COMMON /F2/ IFFLG,CTR
!$    include 'omp_lib.h'
      include 'fftw3.f'
      integer :: NX,NY,NZ,I,J,K,IFFLG,CTR
      double precision::FP(0:NX+1,0:NY+1,0:NZ+1)
      double complex::FPT(0:NX/2,0:NY+1,0:NZ+1)
      INTEGER*8::PLAN,IT_MAX,iret
      double precision,ALLOCATABLE :: FI(:,:)
      double complex, ALLOCATABLE :: FO(:,:)
      LOGICAL :: wisdom_exists
      character(len=80):: VarFileName,cntr
      character(len=8) :: fmt ! format descriptor

      ALLOCATE(FI(1:NX,1:NZ))
      ALLOCATE(FO(1:NX/2+1,1:NZ))


      IFFLG=1 ! Use wisdom? yes:1 no:0
      IT_MAX=1! Nomber of threads for OMP 
!      CTR=0   ! display ?
!$     IT_MAX=omp_get_max_threads()
!      write(6,*) IT_MAX
      do J=1,NY
        !$omp parallel private(I,K)
        !$omp do
        do K=0,NZ-1
        do I=0,NX-1
        FI(I+1,K+1)=FP(I,J,K)
        enddo
        enddo
        !$omp end do
        !$omp end parallel

c wisdom set up------------------------------------------------
!!$     tst = omp_get_wtime()
!$     call dfftw_init_threads(iret)
!$     if(iret.eq.0) then 
!$      write(6,*)'FFTW3 threads-initialization failed'
!$     endif
!$     call dfftw_plan_with_nthreads(IT_MAX)
!!$    ten = omp_get_wtime()
c      TIME_OMP = EN-ST
c      write(6,*)'INITIALIZATION :', TIME_OMP
        IF(IFFLG.EQ.0) THEN
        CALL dfftw_plan_dft_r2c_2d(PLAN,NX,NZ,FI,FO,FFTW_MEASURE)
        ELSEIF(IFFLG.EQ.1) THEN
  
         if(IT_MAX.ge.1)then
          fmt = '(I2.2)'
          write (cntr,fmt) IT_MAX
          VarFileName = 'wisdom_r2c_'//trim(cntr)//'.fftw'
         else
          VarFileName = 'wisdom_r2c.fftw'
         endif
         INQUIRE(FILE=VarFileName, EXIST=wisdom_exists)
         if(wisdom_exists)then
          IF(CTR.EQ.0) THEN
          write(6,*)'Loading widsom from '//VarFileName
          OPEN(66,FILE=VarFileName,STATUS='UNKNOWN')
          call import_wisdom_from_file(iret,66)
c          write(6,*)iret
c         stop
           if(iret.eq.0)then
c           write(6,*)'J = ',J
c           write(6,*)'iret = ',iret
           write(6,*)'FFTW3 wisdom-import failed!'
c           stop
           endif
           write(6,*)'Loading widsom finish'
           CTR=1
           ENDIF
C FOR 2D
          CALL dfftw_plan_dft_r2c_2d(PLAN,NX,NZ,FI,FO,FFTW_MEASURE)
         else
         write(6,*)'Creating plan'
C FOR 2D
         CALL dfftw_plan_dft_r2c_2d(PLAN,NX,NZ,FI,FO,FFTW_MEASURE)
         write(6,*)'Saving wisdom to '//VarFileName
         OPEN(66,FILE=VarFileName)
         call export_wisdom_to_file(66)
         CLOSE(66)
         endif
 
         ELSE
         WRITE(6,*)'ERROR - in PHY2SPC  IFLG_FFTW !',IFFLG
         ENDIF
c--------
         CALL dfftw_execute_dft_r2c(PLAN,FI,FO)
         CALL dfftw_destroy_plan(PLAN)
          
         !$omp parallel private(I,K)
         !$omp do
         do K=0,NZ-1
         do I=0,NX/2
         FPT(I,J,K)=FO(I+1,K+1)/DBLE(NX*NZ)
         enddo
         enddo
         !$omp end do 
         !$omp end parallel 
!$ call DFFTW_CLEANUP_THREADS
      enddo
cccccccccccccccccccccccc

c
      deallocate(FI)
      deallocate(FO)
      RETURN
      END

      SUBROUTINE S2P(FPT,FP,NX,NY,NZ)
      implicit none
      COMMON /F3/ IFFLGC,CTRC
      include 'fftw3.f'
!$    include 'omp_lib.h'
      integer :: NX,NY,NZ,I,J,K,IFFLGC,CTRC
      double precision::FP(0:NX+1,0:NY+1,0:NZ+1)
      double complex::FPT(0:NX/2,0:NY+1,0:NZ+1)
      INTEGER*8::PLAN,IT_MAX,iret
      double precision,ALLOCATABLE :: FI(:,:)
      double complex, ALLOCATABLE :: FO(:,:)
      LOGICAL :: wisdom_exists
      character(len=80):: VarFileName,cntr
      character(len=8) :: fmt ! format descriptor


      ALLOCATE(FI(1:NX,1:NZ))
      ALLOCATE(FO(1:NX/2+1,1:NZ))
      IFFLGC=1
      IT_MAX=1
!$     IT_MAX=omp_get_max_threads()
      do J=1,NY
        !$omp parallel private(I,K)
        !$omp do
        do K=1,NZ
        do I=1,NX/2+1
        FO(I,K)=FPT(I-1,J,K-1)
        enddo
        enddo
        !$omp end do
        !$omp end parallel

c wisdom set up------------------------------------------------
!!$     tst = omp_get_wtime()
!$     call dfftw_init_threads(iret)
!$     if(iret.eq.0) then 
!$      write(6,*)'FFTW3 threads-initialization failed'
!$     endif
!$     call dfftw_plan_with_nthreads(IT_MAX)
!!$    ten = omp_get_wtime()
c      TIME_OMP = EN-ST
c      write(6,*)'INITIALIZATION :', TIME_OMP
        IF(IFFLGC.EQ.0) THEN
        CALL dfftw_plan_dft_c2r_2d(PLAN,NX,NZ,FO,FI,FFTW_MEASURE)
        ELSEIF(IFFLGC.EQ.1) THEN
  
         if(IT_MAX.ge.1)then
          fmt = '(I2.2)'
          write (cntr,fmt) IT_MAX
          VarFileName = 'wisdom_c2r_'//trim(cntr)//'.fftw'
         else
          VarFileName = 'wisdom_c2r.fftw'
         endif
         INQUIRE(FILE=VarFileName, EXIST=wisdom_exists)
         if(wisdom_exists)then
          IF(CTRC.EQ.0) THEN
          write(6,*)'Loading widsom from '//VarFileName
          OPEN(66,FILE=VarFileName,STATUS='UNKNOWN')
          call import_wisdom_from_file(iret,66)
c          write(6,*)iret
c         stop
           if(iret.eq.0)then
c           write(6,*)'J = ',J
c           write(6,*)'iret = ',iret
           write(6,*)'FFTW3 wisdom-import failed!'
c           stop
           endif
           write(6,*)'Loading widsom finish'
           CTRC=1
           ENDIF
C FOR 2D
          CALL dfftw_plan_dft_c2r_2d(PLAN,NX,NZ,FO,FI,FFTW_MEASURE)
         else
         write(6,*)'Creating plan'
C FOR 2D
         CALL dfftw_plan_dft_c2r_2d(PLAN,NX,NZ,FO,FI,FFTW_MEASURE)
         write(6,*)'Saving wisdom to '//VarFileName
         OPEN(66,FILE=VarFileName)
         call export_wisdom_to_file(66)
         CLOSE(66)
         endif
 
         ELSE
         WRITE(6,*)'ERROR - in PHY2SPC  IFLG_FFTW !',IFFLGC
         ENDIF
c--------
        CALL dfftw_execute_dft_c2r(PLAN,FO,FI)
        CALL dfftw_destroy_plan(PLAN)

        do K=0,NZ-1
        do I=0,NX-1
        FP(I,J,K)=FI(I+1,K+1)
        enddo
        enddo 

        do K=0,NZ
        FP(NX,J,K)=FP(0,J,K)
        ENDDO  
        do I=0,NX
        FP(I,J,NZ)=FP(I,J,0)
        ENDDO  
 
!$ call DFFTW_CLEANUP_THREADS
      enddo
      deallocate(FI)
      deallocate(FO)
      RETURN
      END

cccccccccccccccccccccccccccccccccccccccccc
      subroutine write_char(c, iunit)
      character c
      integer iunit
      write(iunit,321) c
 321  format(a,$)
      end      

      subroutine export_wisdom_to_file(iunit)
      integer iunit
      external write_char
      call dfftw_export_wisdom(write_char, iunit)
      end

c     Fortran 77 does not have any portable way to read an arbitrary
c     file one character at a time.  The best alternative seems to be to
c     read a whole line into a buffer, since for fftw-exported wisdom we
c     can bound the line length.  (If the file contains longer lines,
c     then the lines will be truncated and the wisdom import should
c     simply fail.)  Ugh.
      subroutine read_char(ic, iunit)
      integer ic
      integer iunit
      character*256 buf
      save buf
      integer ibuf
      data ibuf/257/
      save ibuf
      if (ibuf .lt. 257) then
         ic = ichar(buf(ibuf:ibuf))
         ibuf = ibuf + 1
         return
      endif
      read(iunit,123,end=666) buf
      ic = ichar(buf(1:1))
      ibuf = 2
      return
 666  ic = -1
      ibuf = 257
 123  format(a256)
      end
      
      subroutine import_wisdom_from_file(isuccess, iunit)
      integer isuccess
      integer iunit
      external read_char
c      call dfftw_import_wisdom(isuccess, read_char, iunit)
      call dfftw_import_wisdom(isuccess, read_char, iunit)
      end
