c
c solver/mkp.f
c (C) 2004 K. Fukagata 
c
c      SUBROUTINE MKP(ALPHA)
!      SUBROUTINE MKP(p,u_temp,v_temp,w_temp,dt,dx,dy,dz,dx12,dy12,dz12,
!     &        Re,nx,ny,nz,res,res_max,res_int,res_var,ompnthreads,
!     &        u_ta,v_ta,w_ta,u_ep2,v_ep2,w_ep2,divu,u,v,w,xl,zl)!KAMETANI
      subroutine MKP(p,u_temp,v_temp,w_temp,dt,dx,dy,dz,dx12,dy12,dz12,
     &     Re,nx,ny,nz,ompnthreads,u_ta,v_ta,w_ta,u_ep2,v_ep2,w_ep2,
     &     divu,u,v,w,xl,zl)
      implicit none
!$    include 'omp_lib.h'
c      COMMON /PP/ PHI
c      COMMON /PT/ PHIT
c     COMMON /PG/ GXPHI,GYPHI,GZPHI
      double precision PI,Re
      double precision xl,zl,DT
      double precision ALPHA
      integer I,J,K,nx,ny,nz
      double precision PHI(0:nx+1,0:ny+1,0:nz+1)
      double precision GXPHI(0:nx+1,0:ny+1,0:nz+1)
      double precision GYPHI(0:nx+1,0:ny+1,0:nz+1)
      double precision GZPHI(0:nx+1,0:ny+1,0:nz+1)
      COMPLEX*16 PHIT(0:nx/2,0:ny+1,0:nz+1)
c MITSUHASHI
      double precision dx(0:nx+1),dy(0:ny+1),dz(0:nz+1)
      double precision x(0:nx+1),y(0:ny+1),z(0:nz+1)
      double precision x12(0:nx+1),y12(0:ny+1),z12(0:nz+1)
      double precision dx12(0:nx+1),dy12(0:ny+1),dz12(0:nz+1)
      double precision p(0:nx+1,0:ny+1,0:nz+1)
      double precision u_temp(0:nx+1,0:ny+1,0:nz+1)
      double precision v_temp(0:nx+1,0:ny+1,0:nz+1)
      double precision w_temp(0:nx+1,0:ny+1,0:nz+1)
      integer ompnthreads
      double precision u(0:nx+1,0:ny+1,0:nz+1)
      double precision v(0:nx+1,0:ny+1,0:nz+1)
      double precision w(0:nx+1,0:ny+1,0:nz+1)
      double precision u_ta(0:nx+1,0:ny+1,0:nz+1)
      double precision v_ta(0:nx+1,0:ny+1,0:nz+1)
      double precision w_ta(0:nx+1,0:ny+1,0:nz+1)
      double precision divu(0:nx+1,0:ny+1,0:nz+1)
      double precision u_ep2(0:nx+1,0:ny+1,0:nz+1)
      double precision v_ep2(0:nx+1,0:ny+1,0:nz+1)
      double precision w_ep2(0:nx+1,0:ny+1,0:nz+1)
c     FOR FFT&TDMA
      double precision AC(1:ny,3)
      double precision KXZ(0:nx/2,0:nz)
      double precision KX (0:nx/2)
      double precision KZ (0:nz)      
C------------------------------------------------------------------------------
C 2. Constants for FFT
C------------------------------------------------------------------------------
      PI=3.14159265358979312D0
c x direction
      DO 200 I=0,nx/2
         KX(I)=DBLE(I)/xl*PI*2.D0
 200  CONTINUE
c z direction
      DO 210 K=0,nz
         IF (K .LE. nz/2) THEN
            KZ(K)=DBLE(K)
         ELSE
            KZ(K)=DBLE(K-nz)
         ENDIF
 210  CONTINUE
      DO 220 K=0,nz-1
         KZ(K)=KZ(K)/zl*PI*2.0D0
 220  CONTINUE
c     Laplacian
      DO 230 K=0,nz
         DO 230 I=0,nx/2
         KXZ(I,K)=2.D0*(1.D0-dCOS(KX(I)*dx(I)))/dx(I)/dx(I)
     +          +2.D0*(1.D0-dCOS(KZ(K)*dz(K)))/dz(K)/dz(K)
 230  CONTINUE
C------------------------------------------------------------------------------
C 3. Matrix coefficients
C------------------------------------------------------------------------------
c For centered veriables (p, ux, uz)
      DO 300 J=1,ny
         AC(J,1)=1.D0/dy12(J-1)/dy(J)
         AC(J,3)=1.D0/dy12(J)  /dy(J)
         AC(J,2)=-AC(J,1)-AC(J,3)
 300  CONTINUE

c 1. RHS of Poisson eq.
c
!$omp parallel firstprivate(I,J,K,DX,DZ,DT,ALPHA)
!$omp do schedule(STATIC) collapse(2)
      DO 110 k=1,nz
      DO 110 j=1,ny         
      DO 110 i=1,nx
         PHI(i,j,k)=(
     +         (u_temp(i,j,k)-u_temp(i-1,j,k))/dx(i)
     +        +(v_temp(i,j,k)-v_temp(i,j-1,k))/dy(j)
     +        +(w_temp(i,j,k)-w_temp(i,j,k-1))/dz(k))/dt
 110  CONTINUE
!$omp end do
!$omp end parallel
 
!$omp parallel 
!$omp do schedule(STATIC) collapse(2)
      DO 120 j=1,ny
      DO 120 i=1,nx
         PHI(i,j,0)=PHI(i,j,nz)
 120  CONTINUE
!$omp end do
!$omp end parallel

!$omp parallel 
!$omp do schedule(STATIC) collapse(2)
      DO 130 k=1,nz
      DO 130 j=1,ny
         PHI(0,j,k)=PHI(nx,j,k)
 130   CONTINUE
!$omp end do
!$omp end parallel
 
      DO 140 j=1,ny
         PHI(0,j,0)=PHI(nx,j,nz)
 140   CONTINUE

 
c
c 2. Solving Poisson eq. using trigonometric expansion 
c
c P2S (external): 2D-FFT from physical to spectral spaces
c S2P (external): 2D-FFT from spectral to physical spaces
c SPN (below)   : Solver for non-zero mode
c SP0B (below)  : Solver for zero mode
c PHIBC (below) : Boundaary condition for PHI (Neumann)
c
      CALL P2S(PHI,PHIT,nx,ny,nz)
      CALL SPN(nx,ny,nz,AC,KXZ,PHIT)
      CALL SP0B(nx,ny,nz,AC,PHIT)
      CALL S2P(PHIT,PHI,nx,ny,nz)
      CALL PHIBC(nx,ny,nz,PHI) 
c
c 3. Compute Gradient Phi
c
!$omp parallel 
!$omp do schedule(STATIC) collapse(2)
      DO 310 k=0,nz
      DO 310 j=0,ny
      DO 310 i=0,nx
         GXPHI(i,j,k)=(PHI(i+1,j,k)-PHI(i,j,k))/dx12(i)
         GYPHI(i,j,k)=(PHI(i,j+1,k)-PHI(i,j,k))/dy12(j)
         GZPHI(i,j,k)=(PHI(i,j,k+1)-PHI(i,j,k))/dz12(k)
 310  CONTINUE
!$omp end do
!$omp end parallel 

c
c 4. Velocity correction
c
c!$omp parallel 
c!$omp do schedule(STATIC) collapse(2)
      DO 410 k=0,nz
      DO 410 j=0,ny
      DO 410 i=0,nx
         u(i,j,k)=u_temp(i,j,k)-dt*GXPHI(i,j,k) 
         v(i,j,k)=v_temp(i,j,k)-dt*GYPHI(i,j,k)
         w(i,j,k)=w_temp(i,j,k)-dt*GZPHI(i,j,k)
 410  CONTINUE
c!$omp end do
c!$omp end parallel 
c

c 5. Pressure 
c
!$omp parallel 
!$omp do schedule(STATIC) collapse(2)
      DO 510 K=1,nz
      DO 510 J=1,ny
      DO 510 I=1,nx
         p(I,J,K)=PHI(I,J,K)
 510  CONTINUE
!$omp end do
!$omp end parallel 

!$omp parallel 
!$omp do schedule(STATIC) collapse(2)
      do k = 1,nz
         do j = 1,ny
            do i = 1,nx
            u_ep2(i,j,k)= -dt*(p(i+1,j,k)-p(i,j,k))/dx12(i)+u_ta(i,j,k)
            v_ep2(i,j,k)= -dt*(p(i,j+1,k)-p(i,j,k))/dy12(j)+v_ta(i,j,k)
            w_ep2(i,j,k)= -dt*(p(i,j,k+1)-p(i,j,k))/dz12(k)+w_ta(i,j,k)
            end do
         end do
      end do
!$omp end do
!$omp end parallel 

!$omp parallel 
!$omp do schedule(STATIC) collapse(2)
      do k = 0,nz
         do j = 0,ny
            do i = 0,nx
               u_ta(i,j,k) = dt*GXPHI(i,j,k)
               v_ta(i,j,k) = dt*GYPHI(i,j,k)
               w_ta(i,j,k) = dt*GZPHI(i,j,k)
            end do
         end do
      end do
!$omp end do
!$omp end parallel 

      do k=1,nz
         do j=1,ny
            do i=1,nx
               divu(i,j,k)=(u(i,j,k)-u(i-1,j,k))/dx(i)
     &                 +(v(i,j,k)-v(i,j-1,k))/dy(j)
     &                 +(w(i,j,k)-w(i,j,k-1))/dz(k)               
            enddo
         enddo
      enddo

      RETURN
      END


      SUBROUTINE SPN(NX,NY,NZ,AC,KXZ,PHIT)
c
c CTRIDIAG (external): TDMA solver for complex matrix
      implicit double precision(a-h,o-z)
!$    include 'omp_lib.h' 
c      implicit double precision(a-h,o-z)
      REAL*8 AA(NY,3)
      COMPLEX*16 PHIT(0:NX/2,0:NY+1,0:NZ+1),B(NY),D(NY)
      double precision AC(1:NY,3)
      double precision KXZ(0:NX/2,0:NZ)

*poption parallel,tlocal(I,J,AA,B)
!$omp parallel private(I,J,AA,B)
!$omp do schedule(STATIC) collapse(2)
      DO 210 K=0,NZ-1
      DO 210 I=0,NX/2
         IF((I+1)*(K+1).GT.1)THEN
            DO 220 J=1,NY
               AA(J,1)=AC(J,1)
               AA(J,2)=AC(J,2)-KXZ(I,K)
               AA(J,3)=AC(J,3)
               B(J)=PHIT(I,J,K)
               D(J)=PHIT(I,J,K)
 220        CONTINUE
            AA(1,2)=AA(1,2)+AA(1,1)
            AA(1,1)=0.D0
            AA(NY,2)=AA(NY,2)+AA(NY,3)
            AA(NY,3)=0.D0
            CALL CTRIDIAG(AA,B,NY)
            DO 230 J=1,NY
               PHIT(I,J,K)=B(J)
 230        CONTINUE
         ENDIF
 210  CONTINUE
!$omp end do
!$omp end parallel
      RETURN
      END


      SUBROUTINE SP0B(NX,NY,NZ,AC,PHIT)
c      INCLUDE '../par.f'
c      INCLUDE '../common.f'
      implicit double precision(a-h,o-z)
c      COMMON /PT/ PHIT
      REAL*8 AA(NY,3)
      COMPLEX*16 PHIT(0:NX/2,0:NY+1,0:NZ+1),D(NY)
      double precision AC(1:NY,3)
      DO 200 J=1,NY
         AA(J,1)=AC(J,1)
         AA(J,2)=AC(J,2)
         AA(J,3)=AC(J,3)
         D(J)=PHIT(0,J,0)
 200  CONTINUE
c      AA(1,2)=AA(1,2)+AA(1,1)
c      AA(1,1)=0.D0
c      AA(NY,2)=AA(NY,2)+AA(NY,3)
c      AA(NY,3)=0.D0
c
      PHIT(0,0,0)=(0.D0,0.D0)
      PHIT(0,1,0)=(0.D0,0.D0)
      DO 230 J=1,NY-1
         PHIT(0,J+1,0)=(D(J)-AA(J,1)*PHIT(0,J-1,0)-AA(J,2)*PHIT(0,J,0))
     +        /AA(J,3)
 230  CONTINUE
      PHIT(0,NY+1,0)=PHIT(0,NY,0)
      RETURN
      END


      SUBROUTINE PHIBC(NX,NY,NZ,PHI)
c      INCLUDE '../par.f'
      implicit double precision(a-h,o-z)
c      COMMON /PP/ PHI
      REAL*8 PHI(0:NX+1,0:NY+1,0:NZ+1)
      DO 100 K=1,NZ
      DO 100 I=0,NX
         PHI(I,0,K)=PHI(I,1,K)
         PHI(I,NY+1,K)=PHI(I,NY,K)
 100  CONTINUE
      DO 110 J=1,NY
      DO 110 I=0,NX
         PHI(I,J,NZ+1)=PHI(I,J,1)
         PHI(I,J,0   )=PHI(I,J,NZ)
 110  CONTINUE
      DO 120 K=0,NZ
      DO 120 J=1,NY
         PHI(NX+1,J,K)=PHI(1, J,K)
         PHI(0,   J,K)=PHI(NX,J,K)
 120  CONTINUE
      DO 130 J=1,NY
         PHI(0,J,0)=PHI(NX,J,NZ) 
         PHI(0,J,NZ+1)=PHI(NX,J,1)
         PHI(NX+1,J,0)=PHI(1,J,NZ)
         PHI(NX+1,J,NZ+1)=PHI(1,J,1)
 130  CONTINUE
      RETURN
      END
