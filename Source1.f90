  
    program Poisson_3D_double_precision
    ! Include modules defined by mkl_poisson.f90 and mkl_dfti.f90 header files
    
    use vorticitySolverUtilities
    use tests
    implicit none
    integer nx,ny,nz, frame
    ! Note that the size of the transform nx must be even !!!
    parameter(nx=63, ny=63, nz=63)
    double precision pi
    parameter(pi=3.14159265358979324D0)
    integer ix, iy, iz, i, stat, j, maxIter, iter
    integer ipar(128)
    double precision ax, bx, ay, by, az, bz, lx, ly, lz, hx, hy, hz, xi, yi, zi, cx, cy, cz
    double precision dt, t , endTime, mu, Uwallbottom, Uwalltop, Uwallleft, Uwallright, eta, rho
    
    ! for 2D problems
    !double precision f(nx+1,ny+1), u(nx+1,ny+1), vt(nx+1,ny+1), vtnew(nx+1,ny+1) !for 2D
    !double precision bd_ax((ny+1)), bd_bx((ny+1)), bd_ay((nx+1)), bd_by((nx+1)) ! for 2D
    !character(4) BCtype
    !double precision dpar((5*nx/2+7))  
    double precision q, errorTol, error, D
    character(len=1024) :: filename
    logical     :: useMKL, printFrames, isStokesFlow
    type(DFTI_DESCRIPTOR), pointer :: xhandle, yhandle
    integer :: STATUS = 0
           
    ! For 3D Problems
     double precision, dimension(:,:,:,:),allocatable  :: f, u, vt, vtnew, err, convTerm, viscousTerm, ub, fb, nablacrossfb ! for 3D
     double precision, dimension(3) :: position, force, normVector
     double precision bd_ax((ny+1),(nz+1)), bd_bx((ny+1),(nz+1)), bd_ay((nx+1),(nz+1)), bd_by((nx+1),(nz+1)), bd_az((nx+1),(ny+1)), bd_bz((nx+1),(ny+1)) !for 3D   
     !character(6) BCtype 
     !double precision dpar(5*(nx+ny)/2+9)
     
     allocate(f(nx+1,ny+1,nz+1,3), u(nx+1,ny+1,nz+1,3), vt(nx+1,ny+1,nz+1,3), vtnew(nx+1,ny+1,nz+1,3), err(nx+1,ny+1,nz+1,3), fb(nx+1,ny+1,nz+1,3), nablacrossfb(nx+1,ny+1,nz+1,3))
     allocate(convTerm(nx+1,ny+1,nz+1,3), viscousTerm(nx+1,ny+1,nz+1,3), ub(nx+1,ny+1,nz+1,3))
     
     
    dt = 0.0002
    t=0
    endTime = 0.02  !1.2
    
    isStokesFlow =1 !0 no, 1 yes 
    useMKL = 1 !0 no, 1 yes
    printFrames =0 ! 0 no, 1 yes
    maxIter = 200
    errorTol = 0.0001
    ub=0
    fb=0
    convTerm=0
    viscousTerm=0
    eta = 1.6
    rho = 1000
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
    !Fb =0
    !call testaddForce2Fb(Fb, hx, hy, hz )
    !
    !print *, 'This program is going to exit.'
    !call EXIT(STATUS)
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

    ! Set the boundary conditions for the stream function, and initialize the velocities ub at the boundaries
    ! There are 6 faces, the stream function has 3 components. set dirichlet boundary conditions
    ! set boundary bd_bx at x=0
    ix=1
    do iy = 1,ny+1
        do iz =1 ,nz+1
            xi=hx*(ix-1)/lx
            yi=hy*(iy-1)/ly
            zi=hz*(iz-1)/lz
            do i=1,3
            f(1,iy,iz,i)=0
            enddo
            ub(ix,iy,iz,1)=0
            ub(ix,iy,iz,2)=0
            ub(ix,iy,iz,3)=0
            !ub(ix,iy,iz,1)=pi*cos(pi*yi)*sin(pi*xi)-sin(pi*xi)*Pi*cos(pi*zi)
            !ub(ix,iy,iz,2)=sin(pi*yi)*pi*cos(pi*zi)-sin(pi*yi)*pi*cos(pi*xi)
            !ub(ix,iy,iz,3)=pi*cos(pi*xi)*sin(pi*zi)-Pi*cos(pi*yi)*sin(pi*zi)
        enddo
    enddo

    ! set boundary bd_bx at x=1
    ix=nx+1
    do iy = 1,ny+1
        do iz =1 ,nz+1
            xi=hx*(ix-1)/lx
            yi=hy*(iy-1)/ly
            zi=hz*(iz-1)/lz
            do i=1,3
            f(nx+1,iy,iz,i)=0
            enddo
            ub(ix,iy,iz,1)=0
            ub(ix,iy,iz,2)=0
            ub(ix,iy,iz,3)=0
            !ub(ix,iy,iz,1)=pi*cos(pi*yi)*sin(pi*xi)-sin(pi*xi)*Pi*cos(pi*zi)
            !ub(ix,iy,iz,2)=sin(pi*yi)*pi*cos(pi*zi)-sin(pi*yi)*pi*cos(pi*xi)
            !ub(ix,iy,iz,3)=pi*cos(pi*xi)*sin(pi*zi)-Pi*cos(pi*yi)*sin(pi*zi)
        enddo
    enddo

    
    
    ! set boundary bd_ay at z=0
    iz=1
    do ix = 1,nx+1
        do iy =1 ,ny+1
            do i=1,3
                f(ix,iy,1,i)=0 !-0.5 
            enddo
            xi=hx*(ix-1)/lx
            yi=hy*(iy-1)/ly
            zi=hz*(iz-1)/lz
            ub(ix,iy,iz,1)=0
            ub(ix,iy,iz,2)=0
            ub(ix,iy,iz,3)=0
            !ub(ix,iy,iz,1)=pi*cos(pi*yi)*sin(pi*xi)-sin(pi*xi)*Pi*cos(pi*zi)
            !ub(ix,iy,iz,2)=sin(pi*yi)*pi*cos(pi*zi)-sin(pi*yi)*pi*cos(pi*xi)
            !ub(ix,iy,iz,3)=pi*cos(pi*xi)*sin(pi*zi)-Pi*cos(pi*yi)*sin(pi*zi)
        enddo
    enddo
    
    ! set boundary bd_by at z=1
    iz = nz+1
    do ix = 1,nx+1
        do iy =1 ,ny+1
            do i=1,3
                f(ix,iy,nz+1,i)=0 !-0.5 
            enddo
            xi=hx*(ix-1)/lx
            yi=hy*(iy-1)/ly
            zi=hz*(iz-1)/lz
            ub(ix,iy,iz,1)=0
            ub(ix,iy,iz,2)=0
            ub(ix,iy,iz,3)=0
            !ub(ix,iy,iz,1)=pi*cos(pi*yi)*sin(pi*xi)-sin(pi*xi)*Pi*cos(pi*zi)
            !ub(ix,iy,iz,2)=sin(pi*yi)*pi*cos(pi*zi)-sin(pi*yi)*pi*cos(pi*xi)
            !ub(ix,iy,iz,3)=pi*cos(pi*xi)*sin(pi*zi)-Pi*cos(pi*yi)*sin(pi*zi)
        enddo
    enddo
    
    ! set boundary bd_ay at y=0
    iy=1
    do ix = 1,nx+1
        do iz =1 ,nz+1
            do i=1,3
                f(ix,1,iz,i)=0 !-0.5
            enddo
            f(ix,iy,iz,3)=-0.5
            xi=hx*(ix-1)/lx
            yi=hy*(iy-1)/ly
            zi=hz*(iz-1)/lz
            ub(ix,iy,iz,1)=1
            ub(ix,iy,iz,2)=0
            ub(ix,iy,iz,3)=0
            !ub(ix,iy,iz,1)=pi*cos(pi*yi)*sin(pi*xi)-sin(pi*xi)*Pi*cos(pi*zi)
            !ub(ix,iy,iz,2)=sin(pi*yi)*pi*cos(pi*zi)-sin(pi*yi)*pi*cos(pi*xi)
            !ub(ix,iy,iz,3)=pi*cos(pi*xi)*sin(pi*zi)-Pi*cos(pi*yi)*sin(pi*zi)
        enddo
    enddo
    
    ! set boundary bd_by at y=1
    iy =ny+1
    do ix = 1,nx+1
        do iz =1 ,nz+1
            do i=1,3
                f(ix,ny+1,iz,i)=0 !-0.5 
            enddo
            f(ix,ny+1,iz,3)=0.5
            xi=hx*(ix-1)/lx
            yi=hy*(iy-1)/ly
            zi=hz*(iz-1)/lz
            ub(ix,iy,iz,1)=1
            ub(ix,iy,iz,2)=0
            ub(ix,iy,iz,3)=0
            !ub(ix,iy,iz,1)=pi*cos(pi*yi)*sin(pi*xi)-sin(pi*xi)*Pi*cos(pi*zi)
            !ub(ix,iy,iz,2)=sin(pi*yi)*pi*cos(pi*zi)-sin(pi*yi)*pi*cos(pi*xi)
            !ub(ix,iy,iz,3)=pi*cos(pi*xi)*sin(pi*zi)-Pi*cos(pi*yi)*sin(pi*zi)
        enddo
    enddo
    
    
    ! Boundary conditions for MKL
    bd_ax = 0
    bd_bx = 0
    bd_ay = 0
    bd_by = 0
    bd_az = 0
    bd_bz = 0 
    
    frame =1
    vtnew =0
    vt = 0
    
    !Set point force in the center of the box in the  x direction
    !position(1) = 0.5
    !position(2) = 0.5
    !position(3) = 0.5
    !
    !force(1) = 0.5
    !force(2) = 0
    !force(3) = 0
    !
    !call addForce2Fb(position, Fb, force, hx, hy, hz )
    !add a torque as a couple
    !position(1) = 0.5
    !position(2) = 0.5
    !position(3) = 0.5
    !
    !force(1) = 0
    !force(2) = 0
    !force(3) = -1
    !
    !D= hx/10
    !normVector(1) =1
    !normVector(2) =1
    !normVector(3) = (-force(1) - force(2) )/force(3)
    !
    !normVector = normVector / sqrt(dot_product(normVector,normVector))
    !
    !print *, normVector
    !print *, ' '
    !!call addForce2Fb(position, Fb, force, hx, hy, hz )
    !
    !call addForce2Fb(position+D/2*normVector, Fb, -crossProd(force,normVector)/D, hx, hy, hz )
    !print *, ' '
    !print *, ' '
    !call addForce2Fb(position-D/2*normVector, Fb, crossProd(force,normVector)/D, hx, hy, hz )
 
    call nablaCrossPer(fb, nablacrossfb, nx, ny, nz, hx, 2)

    ! Main loop
    DO WHILE ( t .LT. endTime ) 
        !call updateVelBCs(ub, nx, ny, nz, hx, hy, hz, lx, ly, lz, t )
            
        !1.  Solve the poisson equation for each component of the stream function    
        if (useMKL ) then
            f= vtnew
            !print *, " component 1 "
            call solvePoisson3DVorticityPerMKL(f(:,:,:,1), bd_ax, bd_bx , bd_ay, bd_by, bd_az, bd_bz, nx, ny, nz, ax, bx, ay, by, az, bz, 1)
            !print *, " component 2 "
            call solvePoisson3DVorticityPerMKL(f(:,:,:,2), bd_ax, bd_bx , bd_ay, bd_by, bd_az, bd_bz, nx, ny, nz, ax, bx, ay, by, az, bz, 2)
            !print *, " component 3 "
            call solvePoisson3DVorticityPerMKL(f(:,:,:,3), bd_ax, bd_bx , bd_ay, bd_by, bd_az, bd_bz, nx, ny, nz, ax, bx, ay, by, az, bz, 3)
            
        else
        
            call solvePoisson3DVorticityPer(f, vt, vtnew, nx, ny, nz, hx, maxIter, 1, errorTol) ! first component of stream function
            call solvePoisson3DVorticityPer(f, vt, vtnew, nx, ny, nz, hx, maxIter, 2, errorTol) ! second component of stream function
            call solvePoisson3DVorticityPer(f, vt, vtnew, nx, ny, nz, hx, maxIter, 3, errorTol) ! third component of stream function
        
        end if
        vt = vtnew

        
        !2. Find the velocity u from the stream function
        
        call nablaCrossPer(f, u, nx, ny, nz, hx, 2)
        
        !3. Compute the vorticity at the boundaries
        
        call updateBoundaryPer(f, vtnew, ub, nx, ny, nz, hx)
       
        
        !vtnew = vt

        !4. advance the vorticity in inner nodes
        
        ! Computes convective term
        ! nabla X ( vt X u)
        if(isStokesFlow)then
            convTerm =0
        else
        call computeConvectiveTermSFPer(f, u, convTerm, nx, ny, nz, hx, hy, hz)
        end if
        ! Computes Viscous term 
        ! (nabla dot nabla) vt
        call computeLaplacianPer(vtnew, nx, ny, nz, hx, hy, hz, viscousTerm)
        
        ! explicit euler on inner nodes        
        vtnew(1:nx+1,2:ny,1:nz+1,:) = vtnew(1:nx+1,2:ny,1:nz+1,:) + dt * ( eta * viscousTerm(1:nx+1,2:ny,1:nz+1,:)/rho - convTerm(1:nx+1,2:ny,1:nz+1,:) + nablacrossfb(1:nx+1,2:ny,1:nz+1,:)/rho ) 
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

    !Compute Velocity
    call nablaCrossPer(f, u, nx, ny, nz, hx, 2)

    !print *, ipar
    ! Printing the final results
    write(*,10) nx
    write(*,11) ny
    print *, ''
    ! Watching the error along the line x=hx

    
    print *, ''
    ! Success message to print if everything is OK

    print *, nx, ny

    open(4,file='vt1.txt')
    !write (4,*), ubound(hinges,1)

    write (4,*), nx+1,ny+1, nz+1
    write (4,*), hx,hy,hz
    do ix=1,nx+1
        do iy=1,ny+1
            do iz=1,nz+1
            write(4,*) vtnew(ix,iy,iz,1)
            enddo
        enddo
    enddo
    close(4)
    open(4,file='vt2.txt')
    !write (4,*), ubound(hinges,1)

    write (4,*), nx+1,ny+1, nz+1
    write (4,*), hx,hy,hz
    do ix=1,nx+1
        do iy=1,ny+1
            do iz=1,nz+1
            write(4,*) vtnew(ix,iy,iz,2)
            enddo
        enddo
    enddo
    close(4)
    open(4,file='vt3.txt')
    !write (4,*), ubound(hinges,1)

    write (4,*), nx+1,ny+1, nz+1
    write (4,*), hx,hy,hz
    do ix=1,nx+1
        do iy=1,ny+1
            do iz=1,nz+1
            write(4,*) vtnew(ix,iy,iz,3)
            enddo
        enddo
    enddo
    close(4)

    open(4,file='sf1.txt')
    !write (4,*), ubound(hinges,1)

    write (4,*), nx+1,ny+1, nz+1
    write (4,*), hx,hy,hz
    do ix=1,nx+1
        do iy=1,ny+1
            do iz=1,nz+1
            write(4,*) f(ix,iy,iz,1)
            enddo
        enddo
    enddo
    close(4)
    
    open(4,file='sf2.txt')
    !write (4,*), ubound(hinges,1)

    write (4,*), nx+1,ny+1, nz+1
    write (4,*), hx,hy,hz
    do ix=1,nx+1
        do iy=1,ny+1
            do iz=1,nz+1
            write(4,*) f(ix,iy,iz,2)
            enddo
        enddo
    enddo
    close(4)
    
    open(4,file='sf3.txt')
    !write (4,*), ubound(hinges,1)

    write (4,*), nx+1,ny+1, nz+1
    write (4,*), hx,hy,hz
    do ix=1,nx+1
        do iy=1,ny+1
            do iz=1,nz+1
            write(4,*) f(ix,iy,iz,3)
            enddo
        enddo
    enddo
    close(4)
    
    open(4,file='results3.txt')
    !write (4,*), ubound(hinges,1)

    write (4,*), nx+1,ny+1, nz+1
    write (4,*), hx,hy,hz
    do ix=1,nx+1
        do iy=1,ny+1
            do iz=1,nz+1
            write(4,*) ub(ix,iy,iz,1)
            enddo
        enddo
    enddo
    close(4)
    
    open(4,file='u1.txt')
    !write (4,*), ubound(hinges,1)

    write (4,*), nx+1,ny+1, nz+1
    write (4,*), hx,hy,hz
    do ix=1,nx+1
        do iy=1,ny+1
            do iz=1,nz+1
            write(4,*) u(ix,iy,iz,1)
            enddo
        enddo
    enddo
    close(4)
    
    open(4,file='u2.txt')
    !write (4,*), ubound(hinges,1)

    write (4,*), nx+1,ny+1, nz+1
    write (4,*), hx,hy,hz
    do ix=1,nx+1
        do iy=1,ny+1
            do iz=1,nz+1
            write(4,*) u(ix,iy,iz,2)
            enddo
        enddo
    enddo
    close(4)
    
    open(4,file='u3.txt')
    !write (4,*), ubound(hinges,1)

    write (4,*), nx+1,ny+1, nz+1
    write (4,*), hx,hy,hz
    do ix=1,nx+1
        do iy=1,ny+1
            do iz=1,nz+1
            write(4,*) u(ix,iy,iz,3)
            enddo
        enddo
    enddo
    close(4)
    
    open(4,file='conv1.txt')
    !write (4,*), ubound(hinges,1)

    write (4,*), nx+1,ny+1, nz+1
    write (4,*), hx,hy,hz
    do ix=1,nx+1
        do iy=1,ny+1
            do iz=1,nz+1
            write(4,*) convTerm(ix,iy,iz,1)
            enddo
        enddo
    enddo
    close(4)
    
    open(4,file='conv2.txt')
    !write (4,*), ubound(hinges,1)

    write (4,*), nx+1,ny+1, nz+1
    write (4,*), hx,hy,hz
    do ix=1,nx+1
        do iy=1,ny+1
            do iz=1,nz+1
            write(4,*) convTerm(ix,iy,iz,2)
            enddo
        enddo
    enddo
    close(4)
    
    open(4,file='conv3.txt')
    !write (4,*), ubound(hinges,1)

    write (4,*), nx+1,ny+1, nz+1
    write (4,*), hx,hy,hz
    do ix=1,nx+1
        do iy=1,ny+1
            do iz=1,nz+1
            write(4,*) convTerm(ix,iy,iz,3)
            enddo
        enddo
    enddo
    close(4)
    
    open(4,file='visc1.txt')
    !write (4,*), ubound(hinges,1)

    write (4,*), nx+1,ny+1, nz+1
    write (4,*), hx,hy,hz
    do ix=1,nx+1
        do iy=1,ny+1
            do iz=1,nz+1
            write(4,*) viscousTerm(ix,iy,iz,1)
            enddo
        enddo
    enddo
    close(4)
    
    open(4,file='visc2.txt')
    !write (4,*), ubound(hinges,1)

    write (4,*), nx+1,ny+1, nz+1
    write (4,*), hx,hy,hz
    do ix=1,nx+1
        do iy=1,ny+1
            do iz=1,nz+1
            write(4,*) viscousTerm(ix,iy,iz,2)
            enddo
        enddo
    enddo
    close(4)
    
    open(4,file='visc3.txt')
    !write (4,*), ubound(hinges,1)

    write (4,*), nx+1,ny+1, nz+1
    write (4,*), hx,hy,hz
    do ix=1,nx+1
        do iy=1,ny+1
            do iz=1,nz+1
            write(4,*) viscousTerm(ix,iy,iz,3)
            enddo
        enddo
    enddo
    close(4)
    
    open(4,file='ncfb1.txt')
    !write (4,*), ubound(hinges,1)

    write (4,*), nx+1,ny+1, nz+1
    write (4,*), hx,hy,hz
    do ix=1,nx+1
        do iy=1,ny+1
            do iz=1,nz+1
            write(4,*) fb(ix,iy,iz,1)
            enddo
        enddo
    enddo
    close(4)
    
    open(4,file='ncfb2.txt')
    !write (4,*), ubound(hinges,1)

    write (4,*), nx+1,ny+1, nz+1
    write (4,*), hx,hy,hz
    do ix=1,nx+1
        do iy=1,ny+1
            do iz=1,nz+1
            write(4,*) fb(ix,iy,iz,2)
            enddo
        enddo
    enddo
    close(4)
    
    open(4,file='ncfb3.txt')
    !write (4,*), ubound(hinges,1)

    write (4,*), nx+1,ny+1, nz+1
    write (4,*), hx,hy,hz
    do ix=1,nx+1
        do iy=1,ny+1
            do iz=1,nz+1
            write(4,*) fb(ix,iy,iz,3)
            enddo
        enddo
    enddo
    close(4)
    
    
    open(4,file='assembledVelResults.txt')
    !write (4,*), ubound(hinges,1)
    write(4,"(A2, A2, A2, A3, A3, A2)") 'x ','y ', 'z ', 'u1 ', 'u2 ', 'u3 '  
    do iz=1,nz+1
        do iy=1,ny+1
    do ix=1,nx+1
        
            
            xi=hx*(ix-1)/lx
            yi=hy*(iy-1)/ly
            zi=hz*(iz-1)/lz
            write(4,"( f6.2, A1, f6.2, A1, f6.2, A1, f6.2, A1, f6.2, A1, f6.2 )") xi, ' ', yi, ' ', zi ,' ',u(ix,iy,iz,1),' ',u(ix,iy,iz,2),' ',u(ix,iy,iz,3)
            enddo
        enddo
    enddo
    close(4)
    
    do ix = 1,nx+1
    do iy = 1,ny+1
        do iz =1 ,nz+1
            xi=hx*(ix-1)/lx
            yi=hy*(iy-1)/ly
            zi=hz*(iz-1)/lz
            f(ix,iy,iz,1)=  sin(pi*yi)*sin(pi*zi) ! right hand side ! yi**2+zi**2  !
            f(ix,iy,iz,2)=  sin(pi*zi)*sin(pi*xi) ! right hand side ! xi**2+zi**2  !
            f(ix,iy,iz,3)=  sin(pi*xi)*sin(pi*yi) ! right hand side ! xi**2+yi**2  !
            u(ix,iy,iz,1)=exp(t)*pi*cos(pi*yi)*sin(pi*xi)-sin(pi*xi)*Pi*cos(pi*zi)
            u(ix,iy,iz,2)=exp(t)*sin(pi*yi)*pi*cos(pi*zi)-sin(pi*yi)*pi*cos(pi*xi)
            u(ix,iy,iz,3)=exp(t)*pi*cos(pi*xi)*sin(pi*zi)-Pi*cos(pi*yi)*sin(pi*zi)
            vt(ix,iy,iz,1) = 2*pi**2*f(ix,iy,iz,1)
            vt(ix,iy,iz,2) = 2*pi**2*f(ix,iy,iz,2)
            vt(ix,iy,iz,3) = 2*pi**2*f(ix,iy,iz,3)
        enddo
    enddo
    enddo
    
    open(4,file='results4.txt')
    !write (4,*), ubound(hinges,1)

    write (4,*), nx+1,ny+1, nz+1
    write (4,*), hx,hy,hz
    do ix=1,nx+1
        do iy=1,ny+1
            do iz=1,nz+1
            write(4,*) u(ix,iy,iz,2)
            enddo
        enddo
    enddo
    close(4)
    
    open(4,file='results5.txt')
    !write (4,*), ubound(hinges,1)

    write (4,*), nx+1,ny+1, nz+1
    write (4,*), hx,hy,hz
    do ix=1,nx+1
        do iy=1,ny+1
            do iz=1,nz+1
            write(4,*) vt(ix,iy,iz,1)
            enddo
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