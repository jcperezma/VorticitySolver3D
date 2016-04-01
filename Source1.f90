    include 'mkl_dfti.f90'
    include 'mkl_poisson.f90'
    program Poisson_3D_double_precision
    ! Include modules defined by mkl_poisson.f90 and mkl_dfti.f90 header files
    use mkl_dfti
    use mkl_poisson
    use poissonSolveJacobi
    use tests
    implicit none
    integer nx,ny,nz, frame
    ! Note that the size of the transform nx must be even !!!
    parameter(nx=31, ny=31, nz=31)
    double precision pi
    parameter(pi=3.14159265358979324D0)
    integer ix, iy, iz, i, stat, j, maxIter, iter
    integer ipar(128)
    double precision ax, bx, ay, by, az, bz, lx, ly, lz, hx, hy, hz, xi, yi, zi, cx, cy, cz
    double precision dt, t , endTime, mu, Uwallbottom, Uwalltop, Uwallleft, Uwallright
    ! for 2D problems
    !double precision f(nx+1,ny+1), u(nx+1,ny+1), vt(nx+1,ny+1), vtnew(nx+1,ny+1) !for 2D
    !double precision bd_ax((ny+1)), bd_bx((ny+1)), bd_ay((nx+1)), bd_by((nx+1)) ! for 2D
    !character(4) BCtype
    !double precision dpar((5*nx/2+7))  
    double precision q, errorTol, error
    character(len=1024) :: filename
    logical     :: useMKL, printFrames
    type(DFTI_DESCRIPTOR), pointer :: xhandle, yhandle
    integer :: STATUS = 0
           
    ! For 3D Problems
     double precision, dimension(:,:,:,:),allocatable  :: f, u, vt, vtnew, err! for 3D
     !double precision bd_ax((ny+1),(nz+1)), bd_bx((ny+1),(nz+1)), bd_ay((nx+1),(nz+1)), bd_by((nx+1),(nz+1)), bd_az((nx+1),(ny+1)), bd_bz((nx+1),(ny+1)) for 3D   
     !character(6) BCtype 
     !double precision dpar(5*(nx+ny)/2+9)
     
     allocate(f(nx+1,ny+1,nz+1,3), u(nx+1,ny+1,nz+1,3), vt(nx+1,ny+1,nz+1,3), vtnew(nx+1,ny+1,nz+1,3), err(nx+1,ny+1,nz+1,3))
     
     
    dt = 0.0002
    t=0
    endTime = 2.4  !1.2

    useMKL = 0
    printFrames =0
    maxIter = 200
    errorTol = 0.001
    
    mu = 0.1
    Uwalltop = 1
    Uwallbottom = -1
    Uwallleft = 0
    Uwallright = 0
    ! Note that proper packing of data in right-hand side array f is
    ! automatically provided by the following declaration of the arrays

    


    ! Printing the header for the example
    print *, ''
    print *, ' Vorticity solver in 3D'
    print *, ' **********************************************'

    ax = 0
    bx = 1
    ay = 0
    by = 1
    az = 0
    bz = 1

    !*******************************************************************************
    q=0.0
    ! Computing the mesh size hx in x-direction
    lx=bx-ax
    hx=lx/nx
    ! Computing the mesh size hy in y-direction
    ly=by-ay
    hy=ly/ny
    ! Computing the mesh size hz in z-direction
    lz=bz-az
    hz=lz/nz
    
    
    !******************** Test Unit **********************************
    call testConvTermComputationSF(f, u, vt, nx, ny, nz, lx, ly, lz, hx, hy ,hz )
   
    print *, 'This program is going to exit.'
    call EXIT(STATUS)
      !******************** end of Test Unit **********************************  

    ! Initialize all variables
    do iy=1,ny+1
        do ix=1,nx+1
            do iz = 1,nz+1 !for 3D
                do i = 1,3
                xi=hx*(ix-1)/lx
                yi=hy*(iy-1)/ly
                zi=hz*(iz-1)/lz
                f(ix,iy,iz,i)=0 ! right hand side
                vt(ix,iy,iz,i)=0
                vtnew(ix,iy,iz,i)=0
                u(ix,iy,iz,i)=0
                enddo
            enddo
        enddo
    enddo

    ! Setting the type of the boundary conditions on each side of the rectangular
    ! domain:
    ! On the boundary laying on the line x=0(=ax) Newmann boundary condition
    ! will be used
    ! On the boundary laying on the line x=1(=bx) Newmann boundary condition
    ! will be used
    ! On the boundary laying on the line y=0(=ay) Dirichlet boundary condition
    ! will be used
    ! On the boundary laying on the line y=1(=by) Dirichlet boundary condition
    !BCtype = 'PPDD' !2D

    ! Setting the values at the boundary for the function f(x,y,z,i) 
    
    ! There are 6 faces, the stream function has 3 components. set dirichlet boundary conditions
    ! set boundary bd_bx at x=0
    do iy = 1,ny+1
        do iz =1 ,nz+1
            do i=1,3
            f(1,iy,iz,i)=0
            enddo
        enddo
    enddo

    ! set boundary bd_bx at x=1
    do iy = 1,ny+1
        do iz =1 ,nz+1
            do i=1,3
            f(nx+1,iy,iz,i)=0
            enddo
        enddo
    enddo

    ! set boundary bd_ay at y=0
    do ix = 1,nx+1
        do iz =1 ,nz+1
            do i=1,3
                f(ix,1,iz,i)=0 !-0.5 
            enddo
        enddo
    enddo
    
    ! set boundary bd_by at y=1
    do ix = 1,nx+1
        do iz =1 ,nz+1
            do i=1,3
                f(ix,ny+1,iz,i)=0 !-0.5 
            enddo
        enddo
    enddo
    
    ! set boundary bd_ay at z=0
    do ix = 1,nx+1
        do iy =1 ,ny+1
            do i=1,3
                f(ix,iy,1,i)=0 !-0.5 
            enddo
        enddo
    enddo
    
    ! set boundary bd_by at z=1
    do ix = 1,nx+1
        do iy =1 ,ny+1
            do i=1,3
                f(ix,iy,nz+1,i)=0 !-0.5 
            enddo
        enddo
    enddo
    
    frame =1


    ! Main loop
    DO WHILE ( t .LT. endTime ) 
      
        !1.  Solve the poisson equation for each component of the stream function    
            
        call solvePoisson3DVorticity(f, vt, vtnew, nx, ny, nz, hx, maxIter, 1, errorTol) ! first component of stream function
        call solvePoisson3DVorticity(f, vt, vtnew, nx, ny, nz, hx, maxIter, 2, errorTol) ! second component of stream function
        call solvePoisson3DVorticity(f, vt, vtnew, nx, ny, nz, hx, maxIter, 3, errorTol) ! third component of stream function

        vt = vtnew

        
        !2. Find the velocity u form the stram function to simplify computations afterwards
        
        call computeVelocities(f, u, nx, ny, nz, hx)
        
        !3. Compute the vorticity at the boundaries
        
        call updateBoundary(f, vt, vtnew, u, nx, ny, nz, hx)
        ! TODO update edges and corners 
        
        
        ! Wall x = 0    
        !vt(1:nx+1,1)= 2/hy**2*(f(1:nx+1,1)- f(1:nx+1,2)) + 2* Uwallbottom /hy 
        !print *, 'corner bottom ', 2/hy**2*(f(1,1)- f(1,2)) + 2* Uwallbottom /hy, '  bottom next node', 2/hy**2*(f(2,1)- f(2,2)) + 2* Uwallbottom /hy  

        ! Wall x = 1    
        !vt(1:nx+1,ny+1)= 2/hy**2*(f(1:nx+1,ny+1)- f(1:nx+1,ny))- 2* Uwalltop /hy 
        !print *, 'corner top ', 2/hy**2*(f(1,ny+1)- f(1,ny)) + 2* Uwallbottom /hy, '  bottom next node', 2/hy**2*(f(2,ny+1)- f(2,ny)) + 2* Uwallbottom /hy  

        vtnew = vt

        !4. advance the vorticity

        !! inner nodes
        !do ix = 2,nx
        !    do iy = 2,ny
        !        i = ix
        !        j = iy
        !        vtnew(ix,iy) = vt(ix,iy) + dt * ( -( f(i,j+1)-f(i,j-1) ) /(2*hy) * (vt(i+1,j)-vt(i-1,j))/(2*hx) &
        !        + ( f(i+1,j)-f(i-1,j) )/(2*hx) * (vt(i,j+1)-vt(i,j-1))/(2*hx)   &
        !        + mu * ( vt(i+1,j)+vt(i-1,j)+vt(i,j+1) + vt(i,j-1) - 4* vt(i,j))/(hx*hy) )
        !    enddo            
        !enddo
        !
        !! left Wall
        !do iy = 2,ny
        !    i = 1
        !    j = iy
        !    vtnew(1,iy) = vt(1,iy) + dt * ( -( f(i,j+1)-f(i,j-1) ) /(2*hy) * (vt(2,j)-vt(nx,j))/(2*hx) &
        !    + ( f(2,j)-f(nx,j) )/(2*hx) * (vt(i,j+1)-vt(i,j-1))/(2*hx)   &
        !    + mu * ( vt(2,j)+vt(nx,j)+vt(i,j+1) + vt(i,j-1) - 4* vt(i,j))/(hx*hy) )
        !enddo
        !
        !
        !! right wall
        !do ix = 2,nx
        !    do iy = 2,ny
        !        i = nx+1
        !        j = iy
        !        vtnew(nx+1,iy) = vt(ix,iy) + dt * ( -( f(i,j+1)-f(i,j-1) ) /(2*hy) * (vt(2,j)-vt(nx,j))/(2*hx) &
        !        + ( f(2,j)-f(nx,j) )/(2*hx) * (vt(i,j+1)-vt(i,j-1))/(2*hx)   &
        !        + mu * ( vt(2,j)+vt(nx,j)+vt(i,j+1) + vt(i,j-1) - 4* vt(i,j))/(hx*hy) )
        !    enddo            
        !enddo
        
        !print *, ' sf(2,2) ' , f(2,2)
    
    if (printFrames) then
       ! Print stream function for this time step
           write (filename, "(A7,I3,A4)") "SFframe", frame ,'.txt'
           print *, trim(filename)
           open(4,file=filename)
        write (4,*), nx+1,ny+1, nz+1
           write (4,*), hx,hy,hz
           do ix=1,nx+1
               do iy=1,ny+1
                   write(4,*) f(ix,iy,iz,i)
                   !   write(4,'(F6.2, A1)',advance="no") vtnew(ix,iy), ' '
               enddo
               !write(4,*) ' '
           enddo
          close(4)    
           write (filename, "(A6,I3,A4)") "Vframe", frame ,'.txt'
           print *, trim(filename)
           open(4,file=filename)
        write (4,*), nx+1,ny+1, nz+1
           write (4,*), hx,hy,hz
           do ix=1,nx+1
               do iy=1,ny+1
                   write(4,*) vtnew(ix,iy,iz,i)
                   !   write(4,'(F6.2, A1)',advance="no") vtnew(ix,iy), ' '
               enddo
               !write(4,*) ' '
           enddo
          close(4)      
    endif
    
        t = t+dt
        print *,' Time ' ,  t
        frame = frame+1
    END DO 


    print *, ipar
    ! Printing the final results
    write(*,10) nx
    write(*,11) ny
    print *, ''
    ! Watching the error along the line x=hx

    
    print *, ''
    ! Success message to print if everything is OK

    print *, nx, ny

    open(4,file='results.txt')
    !write (4,*), ubound(hinges,1)

    write (4,*), nx+1,ny+1, nz+1
    write (4,*), hx,hy,hz
    do ix=1,nx+1
        do iy=1,ny+1
            write(4,*) vtnew(ix,iy,iz,i)
        enddo
    enddo
    close(4)

    open(4,file='results2.txt')
    !write (4,*), ubound(hinges,1)

    write (4,*), nx+1,ny+1, nz+1
    write (4,*), hx,hy,hz
    do ix=1,nx+1
        do iy=1,ny+1
            write(4,*) f(ix,iy,iz,i)
        enddo
    enddo
    close(4)

    ! Jumping over failure message
    go to 1
    ! Failure message to print if something went wrong
999 print *, 'Double precision 2D Poisson example FAILED to compute the solution...'
1   continue
10  format(1x,'The number of mesh intervals in x-direction is nx=',I3)
11  format(1x,'The number of mesh intervals in y-direction is ny=',I3)
12  format(1x,'In the mesh point (',F5.3,',',F5.3,') the error between the computed and the true solution is equal to ', E10.3)



    end