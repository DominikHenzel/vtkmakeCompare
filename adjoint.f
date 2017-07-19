      subroutine set_adj_convect
     &     (u,u_v,v,v_v,w,w_v,Re,dt,nx,ny,nz,dx,dy,dz,
     &     dx12,dy12,dz12,u_final,v_final,w_final)
      implicit double precision(a-h,o-z)
      double precision u(0:nx+1,0:ny+1,0:nz+1)
      double precision v(0:nx+1,0:ny+1,0:nz+1)
      double precision w(0:nx+1,0:ny+1,0:nz+1)
      double precision u_v(0:nx+1,0:ny+1,0:nz+1)
      double precision v_v(0:nx+1,0:ny+1,0:nz+1)
      double precision w_v(0:nx+1,0:ny+1,0:nz+1)
      double precision dx(0:nx+1),dy(0:ny+1),dz(0:nz+1)
      double precision dx12(0:nx+1),dy12(0:ny+1),dz12(0:nz+1)
      double precision u_final(0:nx+1,0:ny+1,0:nz+1)
      double precision v_final(0:nx+1,0:ny+1,0:nz+1)
      double precision w_final(0:nx+1,0:ny+1,0:nz+1)
!$omp parallel do schedule(STATIC) private(i,j,k)  
      do k = 1,nz
         do j = 1,ny
            do i = 1,nx
               u1=-((u_final(i,j,k)+u_final(i+1,j,k))*
     &              (u(i,j,k)+u(i+1,j,k))
     &              -(u_final(i,j,k)+u_final(i-1,j,k))*
     &              (u(i,j,k)+u(i-1,j,k))) /dx12(i)-
     &                          !additional term
     &              u_final(i,j,k)*(u(i-1,j,k)-u(i+1,j,k))/dx12(i)*2.0               

               u2=-((v_final(i+1,j,k)*dx(i+1)+v_final(i,j,k)*dx(i))
     &              *(u(i,j+1,k)+u(i,j,k))
     &              -(v_final(i+1,j-1,k)*dx(i+1)+v_final(i,j-1,k)*dx(i))
     &              *(u(i,j,k)+u(i,j-1,k)))/dx12(i)/dy(j)
     &                          !additional term
     &              -((v_final(i+1,j,k)+v_final(i+1,j-1,k))*dx(i+1)+
     &              (v_final(i,j,k)+v_final(i,j-1,k))*dx(i))/dx12(i)*
     &              (-v(i+1,j,k)-v(i+1,j-1,k)+v(i,j,k)+v(i,j-1,k))
     &              /dx12(i)/2.0d0 ! 7/3 should be diveided by dy??
     &              

               u3=-((w_final(i+1,j,k)*dx(i+1)+w_final(i,j,k)*dx(i))
     &              *(u(i,j,k+1)+u(i,j,k))
     &              -(w_final(i+1,j,k-1)*dx(i+1)+w_final(i,j,k-1)*dx(i))
     &              *(u(i,j,k)+u(i,j,k-1)))/dx12(i)/dz(k)
     &                          !additional term 
     &              -((w_final(i+1,j,k)+w_final(i+1,j,k-1))*dx(i+1)+
     &              (w_final(i,j,k)+w_final(i,j,k-1))*dx(i))/dx12(i)*
     &              (-w(i+1,j,k)-w(i+1,j,k-1)+w(i,j,k)+v(i,j,k-1))
     &              /dx12(i)/2.0d0

               u_v(i,j,k) = 0.25d0*(-u1-u2-u3)*dt+1.0d0*dt
            enddo
         enddo
      enddo
!$omp end parallel do
!$omp parallel do schedule(STATIC) private(i,j,k)  
      do k = 1,nz
         do j = 1,ny
            do i = 1,nx
               v1=((u_final(i,j+1,k)*dy(j+1)+u_final(i,j,k)*dy(j))
     &              *(v(i+1,j,k)+v(i,j,k))
     &              -(u_final(i-1,j+1,k)*dy(j+1)+u_final(i-1,j,k)*dy(j))
     &              *(v(i,j,k)+v(i-1,j,k)))/dy12(j)/dx(i)
     &                          !additional term
     &              -((u_final(i,j+1,k)+u_final(i-1,j+1,k))*dy(j+1)+
     &              (u_final(i,j,k)+u_final(i-1,j,k))*dy(j))/dy12(j)*
     &              (-u(i,j+1,k)-u(i-1,j+1,k)+u(i,j,k)+u(i-1,j,k))
     &              /dy12(j)/2.0d0

               v2=((v_final(i,j,k)+v_final(i,j+1,k))*
     &              (v(i,j,k)+v(i,j+1,k))
     &              -(v_final(i,j,k)+v_final(i,j-1,k))*
     &              (v(i,j,k)+v(i,j-1,k)))/dy12(j)-
     &                          !additional term
     &              v_final(i,j,k)*(v(i-1,j,k)-v(i+1,j,k))*2.0d0/dy12(j)

               v3=((w_final(i,j+1,k)*dy(j+1)+w_final(i,j,k)*dy(j))
     &              *(v(i,j,k+1)+v(i,j,k))
     &              -(w_final(i,j+1,k-1)*dy(j+1)+w_final(i,j,k-1)*dy(j))
     &              *(v(i,j,k)+v(i,j,k-1)))/dy12(j)/dz(k)
     &                          !additional term
     &              -((w_final(i,j,k+1)+w_final(i,j+1,k-1))*dy(j+1)+ 
     &              (w_final(i,j,k)+w_final(i,j,k-1))*dy(j))/dy12(j)*
     &              (-w(i,j+1,k)-w(i,j+1,k-1)+w(i,j,k)+w(i,j,k-1))
     &              /dy12(j)/2.0d0

               v_v(i,j,k)= 0.25d0*(-v1-v2-v3)*dt
            end do
         end do
      end do
!$omp end parallel do
!$omp parallel do schedule(STATIC) private(i,j,k)  
      do k = 1,nz
         do j = 1,ny
            do i = 1,nx
               w1=((u_final(i,j,k+1)*dz(k+1)+u_final(i,j,k)*dz(k))
     &              *(w(i+1,j,k)+w(i,j,k))
     &              -(u_final(i-1,j,k+1)*dz(k+1)+u_final(i-1,j,k)*dz(k))
     &              *(w(i,j,k)+w(i-1,j,k)))/dz12(k)/dx(i)
     &                          !additional term
     &              -((u_final(i,j,k+1)+u_final(i-1,j,k+1))*dz(k+1)+
     &              (u_final(i,j,k)+u_final(i-1,j,k))*dz(k))/dz12(k)*
     &              (-u(i,j,k+1)-u(i-1,j,k+1)+u(i,j,k)+u(i-1,j,k))
     &              /dz12(k)*2.0d0

               w2=((v_final(i,j,k+1)*dz(k+1)+v_final(i,j,k)*dz(k))
     &              *(w(i,j+1,k)+w(i,j,k))
     &              -(v_final(i,j-1,k+1)*dz(k+1)+v_final(i,j-1,k)*dz(k))
     &              *(w(i,j,k)+w(i,j-1,k)))/dz12(k)/dy(j)
     &                          !additional term
     &              -((v_final(i,j,k+1)+v_final(i,j-1,k+1))*dz(k+1)+
     &              (v_final(i,j,k)+v_final(i,j-1,k))*dz(k))/dz12(k)*
     &              (-v(i,j,k+1)-v(i,j-1,k+1)+v(i,j,k)+v(i,j-1,k))
     &              /dz12(k)*2.0d0

               w3=((w_final(i,j,k)+w_final(i,j,k+1))*
     &              (w(i,j,k)+w(i,j,k+1))
     &              -(w_final(i,j,k)+w_final(i,j,k-1))*
     &              (w(i,j,k)+w(i,j,k-1)))/dz12(k)-
     &                          !additional term
     &              w_final(i,j,k)*(w(i,j,k-1)-w(i,j,k+1))*2.0d0/dz12(k)

               w_v(i,j,k) = 0.25d0*(-w1-w2-w3)*dt
            end do
         end do
      end do
!$omp end parallel do
      return
      end

      subroutine set_adj_convection_for_scalar
     &     (u,v,w,Re,dt,nx,ny,nz,dx,dy,dz,
     &     dx12,dy12,dz12,T,u_vT,rhoCp,
     &     u_final,v_final,w_final)
      implicit double precision(a-h,o-z)
      double precision u(0:nx+1,0:ny+1,0:nz+1)
      double precision v(0:nx+1,0:ny+1,0:nz+1)
      double precision w(0:nx+1,0:ny+1,0:nz+1)
      double precision dx(0:nx+1),dy(0:ny+1),dz(0:nz+1)
      double precision dx12(0:nx+1),dy12(0:ny+1),dz12(0:nz+1)
      double precision T(0:nx+1,0:ny+1,0:nz+1)
      double precision u_vT(0:nx+1,0:ny+1,0:nz+1)
      double precision rhoCp(0:nx+1,0:ny+1,0:nz+1)
      double precision u_final(0:nx+1,0:ny+1,0:nz+1)      
      double precision v_final(0:nx+1,0:ny+1,0:nz+1)      
      double precision w_final(0:nx+1,0:ny+1,0:nz+1)      

!$omp parallel do schedule(STATIC) private(i,j,k)        
      do k=1,nz
         do j=1,ny
            do i=1,nx
               uT=(
     &              -((rhoCp(i-1,j,k)*T(i-1,j,k)*dx(i)
     &              +rhoCp(i,j,k)*T(i,j,k)*dx(i-1))*u_final(i-1,j,k))
     &              /2.0d0/dx12(i-1)
     &              +((rhoCp(i,j,k)*T(i,j,k)*dx(i+1)
     &              +rhoCp(i+1,j,k)*T(i+1,j,k)*dx(i))*u_final(i,j,k))
     &              /2.0d0/dx12(i)
     &              )/dx(i)

               vT=(
     &              -((rhoCp(i,j-1,k)*T(i,j-1,k)*dy(j)
     &              +rhoCp(i,j,k)*T(i,j,k)*dy(j-1))*v_final(i,j-1,k))
     &              /2.0d0/dy12(j-1)         
     &              +((rhoCp(i,j,k)*T(i,j,k)*dy(j+1)
     &              +rhoCp(i,j+1,k)*T(i,j+1,k)*dy(j))*v_final(i,j,k))
     &              /2.0d0/dy12(j)
     &              )/dy(j)

               wT=(
     &              -((rhoCp(i,j,k-1)*T(i,j,k-1)*dz(k)
     &              +rhoCp(i,j,k)*T(i,j,k)*dz(k-1))*w_final(i,j,k-1))
     &              /2.0d0/dz12(k-1)         
     &              +((rhoCp(i,j,k)*T(i,j,k)*dz(k+1)
     &              +rhoCp(i,j,k+1)*T(i,j,k+1)*dz(k))*w_final(i,j,k))
     &              /2.0d0/dz12(k)
     &              )/dz(k)

               u_vT(i,j,k)=(uT+vT+wT)*dt/rhoCp(i,j,k)+1.0d0*dt
            enddo
         enddo
      enddo
!$omp end parallel do
      end
      
!##########################################################################
c     set additive scalar for Navier-stokes equation
      subroutine set_adj_NS_scalar
     &     (dt,nx,ny,nz,dx12,dy12,dz12,T,u_s,v_s,w_s,
     &     u_final,v_final,w_final)
      implicit double precision(a-h,o-z)
      double precision dx12(0:nx+1),dy12(0:ny+1),dz12(0:nz+1)
      double precision T_final(0:nx+1,0:ny+1,0:nz+1)
      double precision T(0:nx+1,0:ny+1,0:nz+1)
      double precision u_s(0:nx+1,0:ny+1,0:nz+1)
      double precision v_s(0:nx+1,0:ny+1,0:nz+1)
      double precision w_s(0:nx+1,0:ny+1,0:nz+1)
      double precision u_final(0:nx+1,0:ny+1,0:nz+1)      
      double precision v_final(0:nx+1,0:ny+1,0:nz+1)      
      double precision w_final(0:nx+1,0:ny+1,0:nz+1)      
      do k=1,nz
         do j=1,ny
            do i=1,nx
               u_s(i,j,k)=dt*T_final(i,j,k)
     &              *(T(i+1,j,k)-T(i-1,j,k))/dx12(i)/2
               v_s(i,j,k)=dt*T_final(i,j,k)
     &              *(T(i,j+1,k)-T(i,j-1,k))/dy12(j)/2
               w_s(i,j,k)=dt*T_final(i,j,k)
     &              *(T(i,j,k+1)-T(i,j,k-1))/dz12(k)/2              
               
            enddo
         enddo
      enddo
      end
      
!#######################################################################
!      subroutine adjoint_output
