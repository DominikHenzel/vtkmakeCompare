      subroutine set_levelset(nx,ny,nz,dx,dx12,dy12,dz12,
     &     phi0,phi,x,y,z,cellvalue)
      implicit double precision(A-H,O-Z)
      PARAMETER(nnn=2000)
      PARAMETER(dt=1.0d-3)
      PARAMETER(epsilon=1.0d0-1)
      PARAMETER(pi=3.1415926535d0)
      DOUBLE PRECISION phi(0:nx+1,0:ny+1,0:nz+1)
      DOUBLE PRECISION phi_temp(0:nx+1,0:ny+1,0:nz+1)
      DOUBLE PRECISION phi0(0:nx+1,0:ny+1,0:nz+1)
      DOUBLE PRECISION x(0:nx+1),y(0:ny+1),z(0:nz+1)
      DOUBLE PRECISIOn dx(0:nx+1)
      double precision dx12(0:nx+1),dy12(0:ny+1),dz12(0:nz+1)
      double precision cellvalue(1:nx,1:ny,1:nz,1:4)
      double precision a,am,ap
      double precision b,bm,bp
      double precision c,cm,cp
      double precision d,dm,dp
      double precision e,em,ep
      double precision f,fm,fp
      do k=0,nz+1
         do j=0,ny+1
            do i=0,nx+1
               phi(i,j,k)=phi0(i,j,k)
               phi_temp(i,j,k)=0.0d0
            enddo
         enddo
      enddo
!###############################solve equation ###################
      do it=0,nnn
         updatesum=0.0d0
         do k=1,nz
            do j=1,ny
               do i=1,nx
                  a=(phi(i,j,k)-phi(i-1,j,k))/dx12(i-1)
                  b=(phi(i+1,j,k)-phi(i,j,k))/dx(i)
                  c=(phi(i,j,k)-phi(i,j-1,k))/dy12(j-1)
                  d=(phi(i,j+1,k)-phi(i,j,k))/dy12(j)
                  e=(phi(i,j,k)-phi(i,j,k-1))/dz12(k-1)
                  f=(phi(i,j,k+1)-phi(i,j,k))/dz12(k)
                  if(a .ge. 0.0d0)then
                     ap=a
                     am=0.0d0
                  else
                     ap=0.0d0
                     am=a
                  endif
                  if(b .ge. 0.0d0)then
                     bp=b
                     bm=0.0d0
                  else
                     bp=0.0d0
                     bm=b
                  endif
                  if(c .ge. 0.0d0)then
                     cp=c
                     cm=0.0d0
                  else
                     cp=0.0d0
                     cm=c
                  endif
                  if(d .ge. 0.0d0)then
                     dp=d
                     dm=0.0d0
                  else
                     dp=0.0d0
                     dm=d
                  endif
                  if(e .ge. 0.0d0)then
                     ep=e
                     em=0.0d0
                  else
                     ep=0.0d0
                     em=e
                  endif
                  if(f .ge. 0.0d0)then
                     fp=f
                     fm=0.0d0
                  else
                     fp=0.0d0
                     fm=f
                  endif
                  sss=phi0(i,j,k)/dsqrt(phi0(i,j,k)**2+epsilon**2)
                  if(phi0(i,j,k) .gt. 0.0d0)then
                     GGG=dsqrt(dmax1(ap*ap,bm*bm)+dmax1(cp*cp,dm*dm)+
     &                    dmax1(ep*ep,fm*fm))-1.0d0
                  elseif(phi0(i,j,k) .lt. 0.0d0)then
                     GGG=dsqrt(dmax1(am*am,bp*bp)+dmax1(cm*cm,dp*dp)+
     &                    dmax1(em*em,fp*fp))-1.0d0
                  else
                     GGG=0.0d0
                  endif
                  phi_temp(i,j,k)=phi(i,j,k)-dt*sss*GGG
               enddo
            enddo
         enddo

         do k=1,nz
            do j=1,ny
               do i=1,nx
                  updatesum=updatesum+
     &                 dabs((phi_temp(i,j,k)-phi(i,j,k))/phi(i,j,k))
                  phi(i,j,k)=phi_temp(i,j,k)
               enddo
            enddo
         enddo
         do j=0,ny+1
            do i=0,nx+1
               phi(i,j,0)=phi(i,j,1)
               phi(i,j,nz+1)=phi(i,j,nz)
            enddo
         enddo              
         do k=0,nz+1
            do i=0,nx+1
               phi(i,0,k)=phi(i,1,k)
               phi(i,ny+1,k)=phi(i,ny,k)
            enddo
         enddo              
         do k=0,nz+1
            do j=0,ny+1
               phi(0,j,k)=phi(1,j,k)
               phi(nx+1,j,k)=phi(nx,j,k)
            enddo
         enddo              
         write(6,*)'number of iteration:',it
         write(6,*)'difference=',updatesum/dble(nx*ny*nz)
      enddo


!##########################OUTPUT #############################
      open(300,file='levelset.vtk')
      write(300,"('# vtk DataFile Version 3.8.1')")
      write(300,"('3D structure of vessel')")
      write(300,"('ASCII ')")
      write(300,"('DATASET STRUCTURED_GRID')")
      write(300,"('DIMENSIONS ',3(1x,i3))") nx, ny, nz
      write(300,"('POINTS ',i9,' float')") nx*ny*nz
      do k=1,nz
         do j=1,ny
            do i=1,nx
               write(300,"(3e15.7)") x(i), y(j), z(k)
            enddo
         enddo
      enddo
      write(300,"('POINT_DATA ',i9)") nx*ny*nz      
      write(300,"('SCALARS phi float64')")
      write(300,"('LOOKUP_TABLE default')")
      do k=1,nz
         do j=1,ny
            do i=1,nx
               write(300,"(1e15.7)") phi(i,j,k) 
            enddo
         enddo
      enddo
      write(6,*)'data output'
      close(300)
      do k=1,nz
         do j=1,ny
            do i=1,nx
               if(phi(i,j,k).ge.dx(i))then
                  cellvalue(i,j,k,1)=1.0d0
               elseif(phi(i,j,k).le.-1.0d0*dx(i))then
                  cellvalue(i,j,k,1)=0.0d0
               else
                  cellvalue(i,j,k,1)=(dsin(pi*phi(i,j,k)/2.0d0/dx(i))
     &                 +1.0d0)/2.0d0
               endif
            enddo
         enddo
      enddo        
      open(1111,file="cellvalue1.vtk")
      write(1111,"('# vtk DataFile Version 3.8.1')")
      write(1111,"('3D structure of vessel')")
      write(1111,"('ASCII ')")         
      write(1111,"('DATASET STRUCTURED_GRID')")
      write(1111,"('DIMENSIONS ',3(1x,i3))") nx, ny, nz         
      write(1111,"('POINTS ',i9,' float')") nx*ny*nz
      do k=1,nz
         do j=1,ny
            do i=1,nx
               write(1111,"(3e15.7)") x(i), y(j), z(k)
            enddo
         enddo
      enddo
      
      write(1111,"('POINT_DATA ',i9)") nx*ny*nz
      
!     ! velocity vector (vector block)
      write(1111,"('SCALARS PHI float')")
      write(1111,"('LOOKUP_TABLE default')")
      do k=1,nz
         do j=1,ny
            do i=1,nx
               write(1111,"(e15.7)") cellvalue(i,j,k,1)
            enddo
         enddo
      enddo
      close(1111)           
      end
     
!###############################################read file###############
      subroutine readgrid3d(nx,ny,nz,phi0)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION phi0(0:nx+1,0:ny+1,0:nz+1)
      open(100,file='structure.dat',status='OLD')
      do k=1,nz
         do j=1,ny
            do i=1,nx
               read(100,*)a,b,c,phi0(i,j,k)
            enddo
         enddo
      enddo
      close(100)
      end
      
      
