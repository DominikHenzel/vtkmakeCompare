      program main
      integer nx,ny,nz
      double precision,allocatable,dimension(:)::x,y,z
      double precision,allocatable,dimension(:,:,:):: u,v,w,pres
     &     ,Temp,div,phi
      open(100,file='gridnumber.dat',form='unformatted')
      read(100)nx,ny,nz
      close(100)
      allocate(x(0:nx+1))
      allocate(y(0:nx+1))     
      allocate(z(0:nx+1))          
      allocate(u(0:nx+1,0:ny+1,0:nz+1))  
      allocate(v(0:nx+1,0:ny+1,0:nz+1))
      allocate(w(0:nx+1,0:ny+1,0:nz+1))
      allocate(pres(0:nx+1,0:ny+1,0:nz+1))
      allocate(Temp(0:nx+1,0:ny+1,0:nz+1))
      allocate(div(0:nx+1,0:ny+1,0:nz+1))
      allocate(phi(0:nx+1,0:ny+1,0:nz+1))
      open(300,file='velocity.dat',form='unformatted')
      open(400,file='pressure.dat',form='unformatted')
      open(500,file='temperature.dat',form='unformatted')
      open(600,file='divergence.dat',form='unformatted')
      open(700,file='phi.dat',form='unformatted')
      do k=1,nz
         do j=1,ny
            do i=1,nx
               read(300) u(i,j,k),v(i,j,k),w(i,j,k)
               read(400) pres(i,j,k)
               read(500) temp(i,j,k)
               read(600) div(i,j,k)
               read(700) phi(i,j,k)
            enddo
         enddo
      enddo
      close(300)
      close(400)
      close(500)
      close(600)
      close(700)
      open(200,file='gridx.dat',form='unformatted')   
      open(201,file='gridy.dat',form='unformatted')
      open(202,file='gridz.dat',form='unformatted')
      do i=1,nx
         read(200)x(i)
      enddo
      do j=1,ny
         read(201)y(j)
      enddo
      do k=1,nz
         read(202)z(k)
      enddo
      close(200)
      close(201)
      close(202)
      open(1000,file='maindata.vtk',form='formatted')
      write(1000,"('# vtk DataFile Version 3.8.1')")
      write(1000,"('3D structure of vessel')")
      write(1000,"('ASCII ')")
      write(1000,"('DATASET STRUCTURED_GRID')")
      write(1000,"('DIMENSIONS ',3(1x,i3))") nx, ny, nz
      write(1000,"('POINTS ',i9,' float')") nx*ny*nz
      do k=1,nz
         do j=1,ny
            do i=1,nx
               write(1000,"(3e15.7)") x(i), y(j), z(k)
            enddo
         enddo
      enddo
      write(1000,"('POINT_DATA ',i9)") nx*ny*nz      
      write(1000,"('VECTORS Velocity float64')")
      do k=1,nz
         do j=1,ny
            do i=1,nz
               write(1000,"(3e15.7)") u(i,j,k), v(i,j,k), w(i,j,k)
            enddo
         enddo
      enddo
      write(1000,"('SCALARS Pressure float64')")
      write(1000,"('LOOKUP_TABLE default')")
      do k=1,nz
         do j=1,ny
            do i=1,nx
               write(1000,"(e15.7)")pres(i,j,k)
            enddo
         enddo
      enddo
      write(1000,"('SCALARS Temperature float64')")
      write(1000,"('LOOKUP_TABLE default')")
      do k=1,nz
         do j=1,ny
            do i=1,nx
               write(1000,"(e15.7)")Temp(i,j,k)
            enddo
         enddo
      enddo
      write(1000,"('SCALARS Divergence float64')")
      write(1000,"('LOOKUP_TABLE default')")
      do k=1,nz
         do j=1,ny
            do i=1,nx
               write(1000,"(e15.7)")div(i,j,k)
            enddo
         enddo
      enddo
      write(1000,"('SCALARS phi float64')")
      write(1000,"('LOOKUP_TABLE default')")
      do k=1,nz
         do j=1,ny
            do i=1,nx
               write(1000,"(e15.7)")phi(i,j,k)
            enddo
         enddo
      enddo
      close(1000)
      open(700,file='exit_velocity.dat',form='formatted')
      i= nx
      k= nz/2
      errorsum=0
      do j=1,ny
         write(700,*)y(j),u(i,j,k)
         errorsum=errorsum + dabs(-1.0d1*y(j)*(y(j)-2)-u(i,j,k))
      enddo
      write(6,*)'error=',errorsum/dble(ny)
      close(700)
      open(800,file='exit_temperature.dat',form='formatted')
      i= nx
      k= nz/2
      do j=1,ny
         write(800,*)y(j),Temp(i,j,k)
      enddo
      close(800)
      open(900,file='exit_phi.dat',form='formatted')
      i=nx
      do k=1,nz
         do j=1,ny
            write(900,*)y(j),z(k),phi(i,j,k)
         enddo
      enddo
      close(900)
      deallocate(x)
      deallocate(y)
      deallocate(z)
      deallocate(u)
      deallocate(v)
      deallocate(w)
      deallocate(pres)
      deallocate(Temp)
      deallocate(div)
      deallocate(phi)
      end
