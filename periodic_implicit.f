      subroutine implicit_solver(SUX,nx,dx,dx12,Re,dt,ft,fb,eta,CUX)
      implicit double precision(a-h,o-z)
      double precision SUX(0:nx),u(0:nx+1)
      double precision dx(0:nx+1),dx12(0:nx+1)
      double precision CUX(1:nx)
      double precision alpha(1:nx)
      double precision aa(0:nx+1),bb(0:nx+1),cc(0:nx+1)

      beta = dt/Re/2.0d0
      
      do i=1,nx
         alpha(i) = dt*eta*CUX(i)+1.0d0
      enddo

      do i = 1,nx
         u(i) = SUX(i)
      end do
      
      u(0)=ft
      u(nx+1)=fb
      
      do k=0,nx+1
         aa(k)=0.0d0
         bb(k)=0.0d0
         cc(k)=0.0d0
      enddo

      do i = 1,nx             
         aa(i) = -beta/alpha(i)/dx(i)/dx12(i)
         bb(i) = 1.0d0
     &        +beta/alpha(i)*(1.0d0/dx(i+1)+1.0d0/dx(i))/dx12(i)
         cc(i) = -beta/alpha(i)/dx(i+1)/dx12(i)
      enddo
      aa(0)=0.0d0
      bb(0)=0.0d0
      cc(0)=1.0d0
      topright=-1.0d0
      bottomleft=-1.0d0
      aa(nx+1)=1.0d0
      bb(nx+1)=0.0d0
      cc(nx+1)=0.0d0
      call cyclic(aa,bb,cc,bottomleft,topright,u,SUX,nx)      
      return
      end

!---------------------------------Numerical recipes in Fortran pp67---------------      
      subroutine TDMA(a,b,c,r,u,n)
      INTEGER n,nmax
      INTEGER j
      PARAMETER (NMAX=500)
      DOUBLE PRECISION a(0:n+1),b(0:n+1),c(0:n+1),r(0:n+1),u(0:n+1)
      DOUBLE PRECISION bet,gam(NMAX)
      bet=b(0)
      u(0)=r(0)/bet
      do j=1,n+1
         gam(j)=c(j-1)/bet
         bet=b(j)-a(j)*gam(j)
         u(j)=(r(j)-a(j)*u(j-1))/bet
      enddo
      do j=n,0,-1
         u(j)=u(j)-gam(j+1)*u(j+1)
      enddo
      return
      END

!------------------------------------------------------
      subroutine cyclic(a,b,c,alpha,beta,r,x,n)
      INTEGER n,NMAX
      DOUBLE PRECISION alpha,beta,a(0:n+1),b(0:n+1),c(0:n+1),r(0:n+1),
     &     x(0:n+1)
      INTEGER i
      PARAMETER (NMAX=500)
      DOUBLE PRECISION fact,gamma,bb(0:NMAX),u(0:NMAX),z(0:NMAX)
!     gamma=-b(0)
      gamma=1.0d0
      bb(0)=b(0)-gamma
      bb(n+1)=b(n+1)-alpha*beta/gamma
      do i=1,n
         bb(i)=b(i)
      enddo
      call TDMA(a,bb,c,r,x,n)
      u(0)=gamma
      u(n+1)=alpha
      do i=1,n
         u(i)=0.0d0
      enddo
      call TDMA(a,bb,c,u,z,n)
      fact=(x(0)+beta*x(n+1)/gamma)/(1.0d0+z(0)+beta*z(n+1)/gamma)
      do i=0,n+1
         x(i)=x(i)-fact*z(i)
      enddo
      return
      END

