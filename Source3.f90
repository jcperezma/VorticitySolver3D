    module tests
    use poissonSolveJacobi
    implicit none 
    contains

    subroutine testComputeVelocities(f, u, nx, ny, nz, lx, ly, lz, hx, hy ,hz )
    double precision pi
    parameter(pi=3.14159265358979324D0)
    double precision lx, ly, lz, hx, hy, hz
    integer ix, iy, iz, i
    integer nx, ny, nz
    double precision xi, yi, zi, error
    double precision, dimension(:,:,:,:),allocatable :: f, u, err
    !double precision f(nx+1,ny+1,nz+1,3), u(nx+1,ny+1,nz+1,3), err(nx+1,ny+1,nz+1,3)
    allocate(err(nx+1,ny+1,nz+1,3) )

    ! compute stream function components
    error = 0
    do iy=1,ny+1
        do ix=1,nx+1
            do iz = 1,nz+1 !for 3D
                !do i = 1,3
                xi=hx*(ix-1)/lx
                yi=hy*(iy-1)/ly
                zi=hz*(iz-1)/lz
                f(ix,iy,iz,1)=  sin(pi*yi)*sin(pi*zi) ! right hand side ! yi**2+zi**2  !
                f(ix,iy,iz,2)=  sin(pi*zi)*sin(pi*xi) ! right hand side ! xi**2+zi**2  !
                f(ix,iy,iz,3)=  sin(pi*xi)*sin(pi*yi) ! right hand side ! xi**2+yi**2  !

                !vt(ix,iy,iz,i)=0
                !vtnew(ix,iy,iz,i)=0
                !u(ix,iy,iz,i)=0
                !enddo
            enddo
        enddo
    enddo
    print*, error
    ! find the velocity from the stream function

    call computeVelocities(f, u, nx, ny, nz, hx)

    ! Compute error
    error =0
    do iy=1,ny+1
        do ix=1,nx+1
            do iz = 1,nz+1 !for 3D
                !do i = 1,3
                xi=hx*(ix-1)/lx
                yi=hy*(iy-1)/ly
                zi=hz*(iz-1)/lz
                err(ix,iy,iz,1) = abs( u(ix,iy,iz,1) -( pi*sin(pi*xi)*cos(pi*yi) - pi *sin(pi*xi) * cos(pi*zi) ) ) ! abs( u(ix,iy,iz,1)-(2*yi-2*zi) )
                err(ix,iy,iz,2) = abs( u(ix,iy,iz,2)-(2*zi-2*xi) ) ! abs( u(ix,iy,iz,2)-(2*zi-2*xi) )
                err(ix,iy,iz,3) = abs( u(ix,iy,iz,3)-(2*xi-2*yi) ) ! abs( u(ix,iy,iz,3)-(2*xi-2*yi) )

                !f(ix,iy,iz,1)=sin(pi*yi)*sin(pi*zi) ! right hand side
                !f(ix,iy,iz,2)=sin(pi*zi)*sin(pi*xi) ! right hand side
                !f(ix,iy,iz,3)=sin(pi*xi)*sin(pi*yi) ! right hand side

                !vt(ix,iy,iz,i)=0
                !vtnew(ix,iy,iz,i)=0
                !u(ix,iy,iz,i)=0
                !enddo
            enddo
        enddo
    enddo

    ix=15
    iy=15
    iz=30

    xi=hx*(ix-1)/lx
    yi=hy*(iy-1)/ly
    zi=hz*(iz-1)/lz
    print *, u(ix,iy,iz,1), ' ', ( pi*sin(pi*xi)*cos(pi*yi) - pi *sin(pi*xi) * cos(pi*zi) ), ' yi ', yi , ' zi ', zi 

    open(4,file='results.txt')
    !write (4,*), ubound(hinges,1)

    write (4,*), nx+1,ny+1, nz+1
    write (4,*), hx,hy,hz
    do ix=1,nx+1
        do iy=1,ny+1
            do iz=1,nz+1
                write(4,*) err(ix,iy,iz,1)
            enddo
        enddo
    enddo
    close(4)          

    open(4,file='results2.txt')
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

    print*, error
    deallocate(err)

    end subroutine testComputeVelocities


    subroutine testLaplacian(f, vt, nx, ny, nz, lx, ly, lz, hx, hy ,hz )
    double precision pi
    parameter(pi=3.14159265358979324D0)
    double precision lx, ly, lz, hx, hy, hz
    integer ix, iy, iz, i
    integer nx, ny, nz
    double precision xi, yi, zi, error
    double precision, dimension(:,:,:,:),allocatable :: f, vt, err, laplacian, laplacianComp
    !double precision f(nx+1,ny+1,nz+1,3), u(nx+1,ny+1,nz+1,3), err(nx+1,ny+1,nz+1,3)
    allocate(err(nx+1,ny+1,nz+1,3), laplacian(nx+1,ny+1,nz+1,3), laplacianComp(nx+1,ny+1,nz+1,3) )
    laplacian =0
    ! compute stream function components
    error = 0
    do iy=1,ny+1
        do ix=1,nx+1
            do iz = 1,nz+1 !for 3D
                !do i = 1,3 
                xi=hx*(ix-1)/lx
                yi=hy*(iy-1)/ly
                zi=hz*(iz-1)/lz
                f(ix,iy,iz,1)=  sin(pi*yi)*sin(pi*zi) ! right hand side ! yi**2+zi**2  !
                f(ix,iy,iz,2)=  sin(pi*zi)*sin(pi*xi) ! right hand side ! xi**2+zi**2  !
                f(ix,iy,iz,3)=  sin(pi*xi)*sin(pi*yi) ! right hand side ! xi**2+yi**2  !
                laplacianComp(ix,iy,iz,1) = -2*pi**2*f(ix,iy,iz,1)
                laplacianComp(ix,iy,iz,2) = -2*pi**2*f(ix,iy,iz,2)
                laplacianComp(ix,iy,iz,3) = -2*pi**2*f(ix,iy,iz,3)
                !vt(ix,iy,iz,i)=0
                !vtnew(ix,iy,iz,i)=0
                !u(ix,iy,iz,i)=0
                !enddo
            enddo
        enddo
    enddo


    ! Computes the laplacian of the 3D function

    call computeLaplacian(f, nx, ny, nz, hx, hy, hz, laplacian) !


    ! Compute error
    error =0
    do iy=1,ny+1
        do ix=1,nx+1
            do iz = 1,nz+1 !for 3D
                !do i = 1,3
                xi=hx*(ix-1)/lx
                yi=hy*(iy-1)/ly
                zi=hz*(iz-1)/lz
                err(ix,iy,iz,1) = abs( laplacian(ix,iy,iz,1)-( -2*pi**2*f(ix,iy,iz,1) ) ) ! abs( u(ix,iy,iz,1)-(2*yi-2*zi) )
                err(ix,iy,iz,2) = abs( laplacian(ix,iy,iz,2)-( -2*pi**2*f(ix,iy,iz,2) ) ) ! abs( u(ix,iy,iz,2)-(2*zi-2*xi) )
                err(ix,iy,iz,3) = abs( laplacian(ix,iy,iz,3)-( -2*pi**2*f(ix,iy,iz,3) ) ) ! abs( u(ix,iy,iz,3)-(2*xi-2*yi) )

                !f(ix,iy,iz,1)=sin(pi*yi)*sin(pi*zi) ! right hand side
                !f(ix,iy,iz,2)=sin(pi*zi)*sin(pi*xi) ! right hand side
                !f(ix,iy,iz,3)=sin(pi*xi)*sin(pi*yi) ! right hand side

                !vt(ix,iy,iz,i)=0
                !vtnew(ix,iy,iz,i)=0
                !u(ix,iy,iz,i)=0
                !enddo
            enddo
        enddo
    enddo

    ix=15
    iy=15
    iz=30

    xi=hx*(ix-1)/lx
    yi=hy*(iy-1)/ly
    zi=hz*(iz-1)/lz
    print *, err(ix,iy,iz,1), ' ', err(ix,iy,iz,2), ' ', err(ix,iy,iz,3)

    open(4,file='results.txt')
    !write (4,*), ubound(hinges,1)

    write (4,*), nx+1,ny+1, nz+1
    write (4,*), hx,hy,hz
    do ix=1,nx+1
        do iy=1,ny+1
            do iz=1,nz+1
                write(4,*) err(ix,iy,iz,1)
            enddo
        enddo
    enddo
    close(4)          

    open(4,file='results2.txt')
    !write (4,*), ubound(hinges,1)

    write (4,*), nx+1,ny+1, nz+1
    write (4,*), hx,hy,hz
    do ix=1,nx+1
        do iy=1,ny+1
            do iz=1,nz+1
                write(4,*) laplacian(ix,iy,iz,1)
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
                write(4,*) laplacianComp(ix,iy,iz,1)
            enddo
        enddo
    enddo
    close(4) 

    print*, error
    deallocate(err, laplacian, laplacianComp)

    end subroutine testLaplacian


    subroutine testConvectiveTermComputation(f, u, vt, nx, ny, nz, lx, ly, lz, hx, hy ,hz )
    double precision pi
    parameter(pi=3.14159265358979324D0)
    double precision lx, ly, lz, hx, hy, hz
    integer ix, iy, iz, i
    integer nx, ny, nz
    double precision xi, yi, zi, error
    double precision, dimension(:,:,:,:),allocatable :: f, vt, u, err, conv, convComp
    !double precision f(nx+1,ny+1,nz+1,3), u(nx+1,ny+1,nz+1,3), err(nx+1,ny+1,nz+1,3)
    allocate(err(nx+1,ny+1,nz+1,3), conv(nx+1,ny+1,nz+1,3), convComp(nx+1,ny+1,nz+1,3) )
    conv=0
    ! compute vorticity and velocity components
    error = 0
    do iy=1,ny+1
        do ix=1,nx+1
            do iz = 1,nz+1 !for 3D
                !do i = 1,3 
                xi=hx*(ix-1)/lx
                yi=hy*(iy-1)/ly
                zi=hz*(iz-1)/lz
                f(ix,iy,iz,1)=  sin(pi*yi)*sin(pi*zi) ! right hand side ! yi**2+zi**2  !
                f(ix,iy,iz,2)=  sin(pi*zi)*sin(pi*xi) ! right hand side ! xi**2+zi**2  !
                f(ix,iy,iz,3)=  sin(pi*xi)*sin(pi*yi) ! right hand side ! xi**2+yi**2  !
                u(ix,iy,iz,1) = pi*cos(pi*yi)*sin(pi*xi) - pi*cos(pi*zi)*sin(pi*xi)
                u(ix,iy,iz,2) = pi*sin(pi*yi)*cos(pi*zi) - pi*sin(pi*yi)*cos(pi*xi)
                u(ix,iy,iz,3) = pi*sin(pi*zi)*cos(pi*xi) - pi*cos(pi*yi)*sin(pi*zi)
                vt(ix,iy,iz,1) = 2*pi**2*f(ix,iy,iz,1)
                vt(ix,iy,iz,2) = 2*pi**2*f(ix,iy,iz,2)
                vt(ix,iy,iz,3) = 2*pi**2*f(ix,iy,iz,3)
                convComp(ix,iy,iz,1) = -4*pi**4*sin(pi*yi)*cos(pi*xi)*sin(pi*zi)*(cos(pi*yi)-cos(pi*zi))
                convComp(ix,iy,iz,2) =  4*pi**4*cos(pi*yi)*sin(pi*xi)*sin(pi*zi)*(cos(pi*xi)-cos(pi*zi))
                convComp(ix,iy,iz,3) = -4*pi**4*sin(pi*yi)*sin(pi*xi)*cos(pi*zi)*(cos(pi*xi)-cos(pi*yi))
                !vt(ix,iy,iz,i)=0
                !vtnew(ix,iy,iz,i)=0
                !u(ix,iy,iz,i)=0
                !enddo
            enddo
        enddo
    enddo


    ! Computes the laplacian of the 3D function

    call computeConvectiveTermDirect(vt, u, conv, nx, ny, nz, hx, hy, hz)



    ! Compute error
    error =0
    do iy=1,ny+1
        do ix=1,nx+1
            do iz = 1,nz+1 !for 3D
                !do i = 1,3
                xi=hx*(ix-1)/lx
                yi=hy*(iy-1)/ly
                zi=hz*(iz-1)/lz
                err(ix,iy,iz,1) = abs( conv(ix,iy,iz,1)- convComp(ix,iy,iz,1)) ! abs( u(ix,iy,iz,1)-(2*yi-2*zi) )
                err(ix,iy,iz,2) = abs( conv(ix,iy,iz,2)- convComp(ix,iy,iz,2) ) ! abs( u(ix,iy,iz,2)-(2*zi-2*xi) )
                err(ix,iy,iz,3) = abs( conv(ix,iy,iz,3)- convComp(ix,iy,iz,3) ) ! abs( u(ix,iy,iz,3)-(2*xi-2*yi) )

                !f(ix,iy,iz,1)=sin(pi*yi)*sin(pi*zi) ! right hand side
                !f(ix,iy,iz,2)=sin(pi*zi)*sin(pi*xi) ! right hand side
                !f(ix,iy,iz,3)=sin(pi*xi)*sin(pi*yi) ! right hand side

                !vt(ix,iy,iz,i)=0
                !vtnew(ix,iy,iz,i)=0
                !u(ix,iy,iz,i)=0
                !enddo
            enddo
        enddo
    enddo

    ix=15
    iy=15
    iz=30

    xi=hx*(ix-1)/lx
    yi=hy*(iy-1)/ly
    zi=hz*(iz-1)/lz
    print *, err(ix,iy,iz,1), ' ', err(ix,iy,iz,2), ' ', err(ix,iy,iz,3)

    open(4,file='results.txt')
    !write (4,*), ubound(hinges,1)

    write (4,*), nx+1,ny+1, nz+1
    write (4,*), hx,hy,hz
    do ix=1,nx+1
        do iy=1,ny+1
            do iz=1,nz+1
                write(4,*) err(ix,iy,iz,1)
            enddo
        enddo
    enddo
    close(4)          

    open(4,file='results2.txt')
    !write (4,*), ubound(hinges,1)

    write (4,*), nx+1,ny+1, nz+1
    write (4,*), hx,hy,hz
    do ix=1,nx+1
        do iy=1,ny+1
            do iz=1,nz+1
                write(4,*) conv(ix,iy,iz,1)
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
                write(4,*) convComp(ix,iy,iz,1)
            enddo
        enddo
    enddo
    close(4) 

    print*, error
    deallocate(err, conv, convComp)

    end subroutine testConvectiveTermComputation


    subroutine testConvTermComputationSF(f, u, vt, nx, ny, nz, lx, ly, lz, hx, hy ,hz )
    double precision pi
    parameter(pi=3.14159265358979324D0)
    double precision lx, ly, lz, hx, hy, hz
    integer ix, iy, iz, i
    integer nx, ny, nz
    double precision xi, yi, zi, error
    double precision, dimension(:,:,:,:),allocatable :: f, vt, u, err, conv, convComp,WxU
    !double precision f(nx+1,ny+1,nz+1,3), u(nx+1,ny+1,nz+1,3), err(nx+1,ny+1,nz+1,3)
    allocate(err(nx+1,ny+1,nz+1,3), conv(nx+1,ny+1,nz+1,3), convComp(nx+1,ny+1,nz+1,3) , WxU(nx+1,ny+1,nz+1,3))
    conv=0
    ! compute vorticity and velocity components
    error = 0

    do ix=1,nx+1
        do iy=1,ny+1
            do iz = 1,nz+1 !for 3D
                !do i = 1,3 
                xi=hx*(ix-1)/lx
                yi=hy*(iy-1)/ly
                zi=hz*(iz-1)/lz
                f(ix,iy,iz,1)=  sin(pi*yi)*sin(pi*zi) ! right hand side ! yi**2+zi**2  !
                f(ix,iy,iz,2)=  sin(pi*zi)*sin(pi*xi) ! right hand side ! xi**2+zi**2  !
                f(ix,iy,iz,3)=  sin(pi*xi)*sin(pi*yi) ! right hand side ! xi**2+yi**2  !
                u(ix,iy,iz,1) = pi*cos(pi*yi)*sin(pi*xi) - pi*cos(pi*zi)*sin(pi*xi)
                u(ix,iy,iz,2) = pi*sin(pi*yi)*cos(pi*zi) - pi*sin(pi*yi)*cos(pi*xi)
                u(ix,iy,iz,3) = pi*sin(pi*zi)*cos(pi*xi) - pi*cos(pi*yi)*sin(pi*zi)
                vt(ix,iy,iz,1) = 2*pi**2*f(ix,iy,iz,1)
                vt(ix,iy,iz,2) = 2*pi**2*f(ix,iy,iz,2)
                vt(ix,iy,iz,3) = 2*pi**2*f(ix,iy,iz,3)
                convComp(ix,iy,iz,1) = -4*pi**4*sin(pi*yi)*cos(pi*xi)*sin(pi*zi)*(cos(pi*yi)-cos(pi*zi))
                convComp(ix,iy,iz,2) =  4*pi**4*cos(pi*yi)*sin(pi*xi)*sin(pi*zi)*(cos(pi*xi)-cos(pi*zi))
                convComp(ix,iy,iz,3) = -4*pi**4*sin(pi*yi)*sin(pi*xi)*cos(pi*zi)*(cos(pi*xi)-cos(pi*yi))
                WxU(ix,iy,iz,1) = -2*pi**3*sin(pi*xi)*(cos(pi*zi)**2*cos(pi*xi)-cos(pi*zi)**2*cos(pi*yi)-cos(pi*zi)*cos(pi*yi)**2+cos(pi*xi)*cos(pi*yi)**2+cos(pi*zi)-2*cos(pi*xi)+cos(pi*yi))
                !vt(ix,iy,iz,i)=0
                !vtnew(ix,iy,iz,i)=0
                !u(ix,iy,iz,i)=0
                !enddo
            enddo
        enddo
    enddo
    ix=1
    iy=15
    iz=16

    xi=hx*(ix-1)/lx
    yi=hy*(iy-1)/ly
    zi=hz*(iz-1)/lz
   
    print *, 'u analytical ' , u(ix,iy,iz,1), ' ', u(ix,iy,iz,2), ' ', u(ix,iy,iz,3)
    print *, 'vt analytical ' , vt(ix,iy,iz,1), ' ', vt(ix,iy,iz,2), ' ', vt(ix,iy,iz,3)
    
    print *, 'WxU analytical ' , WxU(ix,iy,iz,1)
    print *, u(ix,iy,iz,3) , ' ', vt(ix,iy,iz,2),' ' , u(ix,iy,iz,2),' ',vt(ix,iy,iz,3)
    
    
    !convComp = WxU
    ! Compute velocity from the stream function 
    call nablaCross(f, u, nx, ny, nz, hx,1)
    !For some reason that I dont understand for the function Im testing the velocity has to be computed using first order FD at the borders, 
    !while the vorticity with a second order. This will lead to the smallest error. 

    ! Computes the laplacian of the 3D function

    call computeConvectiveTermSF(f, u, conv, nx, ny, nz, hx, hy, hz)

    

    ! Compute error
    error =0
    do iy=1,ny+1
        do ix=1,nx+1
            do iz = 1,nz+1 !for 3D
                !do i = 1,3
                xi=hx*(ix-1)/lx
                yi=hy*(iy-1)/ly
                zi=hz*(iz-1)/lz
                err(ix,iy,iz,1) = abs( conv(ix,iy,iz,1)- convComp(ix,iy,iz,1)) ! abs( u(ix,iy,iz,1)-(2*yi-2*zi) )
                err(ix,iy,iz,2) = abs( conv(ix,iy,iz,2)- convComp(ix,iy,iz,2) ) ! abs( u(ix,iy,iz,2)-(2*zi-2*xi) )
                err(ix,iy,iz,3) = abs( conv(ix,iy,iz,3)- convComp(ix,iy,iz,3) ) ! abs( u(ix,iy,iz,3)-(2*xi-2*yi) )

                !f(ix,iy,iz,1)=sin(pi*yi)*sin(pi*zi) ! right hand side
                !f(ix,iy,iz,2)=sin(pi*zi)*sin(pi*xi) ! right hand side
                !f(ix,iy,iz,3)=sin(pi*xi)*sin(pi*yi) ! right hand side

                !vt(ix,iy,iz,i)=0
                !vtnew(ix,iy,iz,i)=0
                !u(ix,iy,iz,i)=0
                !enddo
            enddo
        enddo
    enddo
     ix=1
    iy=15
    iz=16
 print *, 'error ' , err(ix,iy,iz,1), ' ', err(ix,iy,iz,2), ' ', err(ix,iy,iz,3)
    
    open(4,file='results.txt')
    !write (4,*), ubound(hinges,1)

    write (4,*), nx+1,ny+1, nz+1
    write (4,*), hx,hy,hz
    do ix=1,nx+1
        do iy=1,ny+1
            do iz=1,nz+1
                write(4,*) err(ix,iy,iz,1)
            enddo
        enddo
    enddo
    close(4)          

    open(4,file='results2.txt')
    !write (4,*), ubound(hinges,1)

    write (4,*), nx+1,ny+1, nz+1
    write (4,*), hx,hy,hz
    do ix=1,nx+1
        do iy=1,ny+1
            do iz=1,nz+1
                write(4,*) conv(ix,iy,iz,1)
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
                write(4,*) convComp(ix,iy,iz,1)
            enddo
        enddo
    enddo
    close(4) 

    print*, error
    deallocate(err, conv, convComp,WxU)

    end subroutine testConvTermComputationSF

    end module tests