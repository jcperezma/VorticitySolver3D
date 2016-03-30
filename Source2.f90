    module poissonSolveJacobi
    implicit none 
    contains
    subroutine solvePoisson3DVorticity(f, vt, vtnew, nx, ny, nz, hx, maxIter, i, errorTol)

    double precision f(nx+1,ny+1,nz+1,3), vt(nx+1,ny+1,nz+1,3), vtnew(nx+1,ny+1,nz+1,3) ! 
    double precision errorTol, error, hx
    integer nx,ny,nz, maxIter, i, iter, ix, iy, iz
    ! i is the component

    !Solves the poisson equation in 3d using SOR. the SOR factor is 1.5
    !assumes a uniform grid, that is hx=hy=hz


    do iter = 1 , maxIter
        vt=f;
        ! update inner nodes
        do ix=2,nx
            do iy=2,ny
                do iz=2,nz
                    f(ix,iy,iz,i)  = 1.5 * 0.16666 * ( hx*hx * vtnew(ix,iy,iz,i) + f(ix-1,iy,iz,i)+f(ix+1,iy,iz,i) + f(ix,iy-1,iz,i)+f(ix,iy+1,iz,i)  &
                    +  f(ix,iy,iz-1,i) +f(ix,iy,iz+1,i) )+ (1 - 1.5) * f(ix,iy,iz,i)
                enddo
            enddo
        enddo

        ! update boundary nodes for the normal component with a neumann boundary condition df/dn =0  f(1) - f(-1) = 0

        SELECT CASE (i) ! select component
        CASE (1)
            ! x component 
            ! x=0 
            ix =1
            do iy=2,ny
                do iz=2,nz
                    f(ix,iy,iz,i)  = 1.5 * 0.16666 * ( hx*hx * vtnew(ix,iy,iz,i) + 2*f(ix+1,iy,iz,i) + f(ix,iy-1,iz,i)+f(ix,iy+1,iz,i)  &
                    +  f(ix,iy,iz-1,i) +f(ix,iy,iz+1,i) )+ (1 - 1.5) * f(ix,iy,iz,i)
                enddo
            enddo

            ! x=1
            ix = nx+1
            do iy=2,ny
                do iz=2,nz
                    f(ix,iy,iz,i)  = 1.5 * 0.16666 * ( hx*hx * vtnew(ix,iy,iz,i) + 2*f(ix-1,iy,iz,i) + f(ix,iy-1,iz,i)+f(ix,iy+1,iz,i)  &
                    +  f(ix,iy,iz-1,i) +f(ix,iy,iz+1,i) )+ (1 - 1.5) * f(ix,iy,iz,i)
                enddo
            enddo


        CASE (2)
            ! y component 
            ! y = 0
            iy =1
            do ix=2,nx
                do iz=2,nz    
                    f(ix,iy,iz,i)  = 1.5 * 0.16666 * ( hx*hx * vtnew(ix,iy,iz,i) + f(ix-1,iy,iz,i)+f(ix+1,iy,iz,i) + 2 * f(ix,iy+1,iz,i)  &
                    + f(ix,iy,iz-1,i) +f(ix,iy,iz+1,i) )+ (1 - 1.5) * f(ix,iy,iz,i)
                enddo
            enddo 
            ! y = 1
            iy =ny+1
            do ix=2,nx
                do iz=2,nz                   
                    f(ix,iy,iz,i)  = 1.5 * 0.16666 * ( hx*hx * vtnew(ix,iy,iz,i) + f(ix-1,iy,iz,i)+f(ix+1,iy,iz,i) + 2 *f(ix,iy-1,iz,i)  &
                    + f(ix,iy,iz-1,i) +f(ix,iy,iz+1,i) )+ (1 - 1.5) * f(ix,iy,iz,i)
                enddo
            enddo  

        CASE (3)
            ! z component
            ! z = 0
            iz =1
            do ix=2,nx
                do iy=2,ny                        
                    f(ix,iy,iz,i)  = 1.5 * 0.16666 * ( hx*hx * vtnew(ix,iy,iz,i) + f(ix-1,iy,iz,i)+f(ix+1,iy,iz,i) + f(ix,iy-1,iz,i)+f(ix,iy+1,iz,i)  &
                    + 2 *f(ix,iy,iz+1,i) )+ (1 - 1.5) * f(ix,iy,iz,i)
                enddo
            enddo 
            ! z = 1
            iz =nz+1
            do ix=2,nx
                do iy=2,ny                        
                    f(ix,iy,iz,i)  = 1.5 * 0.16666 * ( hx*hx * vtnew(ix,iy,iz,i) + f(ix-1,iy,iz,i)+f(ix+1,iy,iz,i) + f(ix,iy-1,iz,i)+f(ix,iy+1,iz,i)  &
                    + 2 *f(ix,iy,iz-1,i) )+ (1 - 1.5) * f(ix,iy,iz,i)
                enddo
            enddo     

        END SELECT

        error = 0 

        do iy=1,ny+1
            do ix=1,nx+1
                do iz=1,nz+1 
                    error = error + abs(vt(ix,iy,iz,i) - f(ix,iy,iz,i))
                enddo
            enddo
        enddo

        !print *, ' error ' , error 
        if (error .lt. errorTol) EXIT 


    enddo



    end subroutine solvePoisson3DVorticity


    subroutine updateBoundary(f, vt, vtnew, u, nx, ny, nz, hx)
    double precision f(nx+1,ny+1,nz+1,3), vt(nx+1,ny+1,nz+1,3), vtnew(nx+1,ny+1,nz+1,3),u(nx+1,ny+1,nz+1,3) ! 
    double precision errorTol, error, hx
    integer nx,ny,nz, maxIter, i, iter, ix, iy, iz
    ! use U and f to compute the vorticity at the boundaries
    !!!!!!!!!!!!!!!!!!!!!!!!!!!! For faces Z =const !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! for face Z =0
    ! normal component for this face is 3
    ! for inner nodes   
    iz = 1
    do iy=2,ny
        do ix=2,nx
            vt(ix,iy,iz,3) = 0.5/hx* ( u(ix+1,iy,iz,2) - u(ix-1,iy,iz,2) - u(ix,iy+1,iz,1) + u(ix,iy-1,iz,1) )
        enddo
    enddo

    ! tangential components are 1 and 2
    ! Component 1
    ! inner nodes
    do iy=2,ny
        do ix=2,nx
            vt(ix,iy,iz,1) = -2* f(ix,iy,iz+1,1) / hx**2 + 2/hx* ( (f(ix+1,iy,iz,3) - f(ix-1,iy,iz,3))/(2*hx) + u(ix,iy,iz,2))
        enddo
    enddo

    ! Component 2
    ! inner nodes
    do iy=2,ny
        do ix=2,nx
            vt(ix,iy,iz,1) = -2* f(ix,iy,iz+1,2) / hx**2 + 2/hx* ( (f(ix,iy+1,iz,3) - f(ix,iy-1,iz,3))/(2*hx) - u(ix,iy,iz,1))
        enddo
    enddo


    ! for face Z =1
    ! normal component for this face is 3
    ! for inner nodes   
    iz = nz+1
    do iy=2,ny
        do ix=2,nx
            vt(ix,iy,iz,3) = 0.5/hx* ( u(ix+1,iy,iz,2) - u(ix-1,iy,iz,2) - u(ix,iy+1,iz,1) + u(ix,iy-1,iz,1) )
        enddo
    enddo

    ! tangential components are 1 and 2
    ! Component 1
    ! inner nodes
    do iy=2,ny
        do ix=2,nx
            vt(ix,iy,iz,1) = -2* f(ix,iy,iz-1,1) / hx**2 - 2/hx* ( (f(ix+1,iy,iz,3) - f(ix-1,iy,iz,3))/(2*hx) + u(ix,iy,iz,2))
        enddo
    enddo

    ! Component 2
    ! inner nodes
    do iy=2,ny
        do ix=2,nx
            vt(ix,iy,iz,1) = -2* f(ix,iy,iz-1,2) / hx**2 - 2/hx* ( (f(ix,iy+1,iz,3) - f(ix,iy-1,iz,3))/(2*hx) - u(ix,iy,iz,1))
        enddo
    enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!! For faces X =const !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! for face X =0
    ! normal component for this face is 1
    ! for inner nodes   
    ix = 1
    do iy=2,ny
        do iz=2,nz
            vt(ix,iy,iz,1) = 0.5/hx* ( u(ix,iy+1,iz,3) - u(ix,iy-1,iz,3) - u(ix,iy,iz+1,2) + u(ix,iy,iz-1,2) )
        enddo
    enddo

    ! tangential components are 2 and 3
    ! Component 2
    ! inner nodes
    do iy=2,ny
        do iz=2,nz
            vt(ix,iy,iz,2) = -2* f(ix+1,iy,iz,2) / hx**2 + 2/hx* ( (f(ix,iy+1,iz,1) - f(ix,iy-1,iz,1))/(2*hx) + u(ix,iy,iz,3))
        enddo
    enddo

    ! Component 3
    ! inner nodes
    do iy=2,ny
        do iz=2,nz
            vt(ix,iy,iz,3) = -2* f(ix+1,iy,iz,3) / hx**2 + 2/hx* ( (f(ix,iy,iz+1,1) - f(ix,iy,iz-1,1))/(2*hx) - u(ix,iy,iz,2))
        enddo
    enddo


    ! for face X =1
    ! normal component for this face is 3
    ! for inner nodes   
    ix = nx+1
    do iy=2,ny
        do iz=2,nz
            vt(ix,iy,iz,1) = 0.5/hx* ( u(ix,iy+1,iz,3) - u(ix,iy-1,iz,3) - u(ix,iy,iz+1,2) + u(ix,iy,iz-1,2) )
        enddo
    enddo

    ! tangential components are 2 and 3
    ! Component 2
    ! inner nodes
    do iy=2,ny
        do iz=2,nz
            vt(ix,iy,iz,2) = -2* f(ix-1,iy,iz,2) / hx**2 - 2/hx* ( (f(ix,iy+1,iz,1) - f(ix,iy-1,iz,1))/(2*hx) + u(ix,iy,iz,3))
        enddo
    enddo

    ! Component 3
    ! inner nodes
    do iy=2,ny
        do iz=2,nz
            vt(ix,iy,iz,3) = -2* f(ix-1,iy,iz,3) / hx**2 - 2/hx* ( (f(ix,iy,iz+1,1) - f(ix,iy,iz-1,1))/(2*hx) - u(ix,iy,iz,2))
        enddo
    enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!! For faces Y =const !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! for face Y =0
    ! normal component for this face is 2
    ! for inner nodes   
    iy = 1
    do ix=2,nx
        do iz=2,nz
            vt(ix,iy,iz,2) = 0.5/hx* ( u(ix,iy,iz+1,1) - u(ix,iy,iz-1,1) - u(ix+1,iy,iz,3) + u(ix-1,iy,iz,3) )
        enddo
    enddo

    ! tangential components are 1 and 3
    ! Component 1
    ! inner nodes
    do ix=2,nx
        do iz=2,nz
            vt(ix,iy,iz,1) = -2* f(ix,iy+1,iz,2) / hx**2 + 2/hx* ( (f(ix+1,iy,iz,2) - f(ix-1,iy,iz,2))/(2*hx) - u(ix,iy,iz,3))
        enddo
    enddo

    ! Component 3
    ! inner nodes
    do ix=2,nx
        do iz=2,nz
            vt(ix,iy,iz,3) = -2* f(ix,iy+1,iz,3) / hx**2 + 2/hx* ( (f(ix,iy,iz+1,1) - f(ix,iy,iz-1,1))/(2*hx) + u(ix,iy,iz,1))
        enddo
    enddo

    ! for face Y =1
    ! normal component for this face is 2
    ! for inner nodes   
    iy = ny+1
    do ix=2,nx
        do iz=2,nz
            vt(ix,iy,iz,2) = 0.5/hx* ( u(ix,iy,iz+1,1) - u(ix,iy,iz-1,1) - u(ix+1,iy,iz,3) + u(ix-1,iy,iz,3) )
        enddo
    enddo

    ! tangential components are 1 and 3
    ! Component 1
    ! inner nodes
    do iy=2,ny
        do iz=2,nz
            vt(ix,iy,iz,1) = -2* f(ix,iy-1,iz,2) / hx**2 - 2/hx* ( (f(ix+1,iy,iz,2) - f(ix-1,iy,iz,2))/(2*hx) - u(ix,iy,iz,3))
        enddo
    enddo

    ! Component 3
    ! inner nodes
    do iy=2,ny
        do iz=2,nz
            vt(ix,iy,iz,3) = -2* f(ix,iy-1,iz,3) / hx**2 - 2/hx* ( (f(ix,iy,iz+1,1) - f(ix,iy,iz-1,1))/(2*hx) + u(ix,iy,iz,1))
        enddo
    enddo

    end subroutine updateBoundary


    subroutine computeVelocities(f, u, nx, ny, nz, hx)
    implicit none
    double precision, dimension(:,:,:,:),allocatable :: f, u
    !double precision f(nx+1,ny+1,nz+1,3), u(nx+1,ny+1,nz+1,3)! 
    double precision hx, xCoeff(3), yCoeff(3), zCoeff(3), error
    integer nx,ny,nz, i, j, k, ix, ix2, iy, iy2, iz, iz2, indexX, indexY, indexZ
    ! Computes the components of the velocity from the stream function derivatives
    ! assumes a uniform grid; hx=hy=hz, though to change this is really easy.
    ! I tested it with simple functions and with the functions used by liu the max error is 0.25%

    do i=1,nx+1
        if(i .EQ. 1) then
            indexX =1   ! 0
            !Forward difference
            ix = 2      ! 1
            ix2 = 3     ! 2
            xCoeff(1) = -3.0/2.0
            xCoeff(2) = 2.0
            xCoeff(3) = -0.50
        elseif(i .EQ. nx+1) then
            !Backward difference
            indexX =nx+1 ! 0
            ix = nx      ! -1   
            ix2 = nx-1   ! -2
            xCoeff(1) = 3.0/2.0
            xCoeff(2) = -2.0
            xCoeff(3) = 0.50
        else
            !Central difference 
            indexX = i  !0   1
            ix =i-1     !-1  2
            ix2 = i+1   !1   3
            xCoeff(1) = 0
            xCoeff(2) = -0.5
            xCoeff(3) = 0.5
        endif
        do j=1,ny+1
            if(j .EQ. 1) then
                indexY=1   ! 0
                !Forward difference
                iy = 2      ! 1
                iy2 = 3     ! 2
                yCoeff(1) = -3.0/2.0
                yCoeff(2) = 2.0
                yCoeff(3) = -0.50
            elseif(j .EQ. ny+1) then
                !Backward difference
                indexY =ny+1 ! 0
                iy = ny      ! -1   
                iy2 = ny-1   ! -2
                yCoeff(1) = 3.0/2.0
                yCoeff(2) = -2.0
                yCoeff(3) = 0.50
            else
                !Central difference
                indexY = j  !0   1
                iy =j-1     !-1  2
                iy2 = j+1   !1   3
                yCoeff(1) = 0
                yCoeff(2) = -0.5
                yCoeff(3) = 0.5
            endif
            do k=1,nz+1
                if(k .EQ. 1) then
                    indexZ=1   ! 0
                    !Forward difference
                    iz = 2      ! 1
                    iz2 = 3     ! 2
                    zCoeff(1) = -3.0/2.0
                    zCoeff(2) = 2.0
                    zCoeff(3) = -0.50
                elseif (k .EQ. nz+1) then
                    !Backward difference
                    indexZ =nz+1 ! 0
                    iz = nz      ! -1   
                    iz2 = nz-1   ! -2
                    zCoeff(1) = 3.0/2.0
                    zCoeff(2) = -2.0
                    zCoeff(3) = 0.50
                else
                    !central difference
                    indexZ = k  !0  1
                    iz =k-1     !-1 2
                    iz2 = k+1   !1  3
                    zCoeff(1) = 0
                    zCoeff(2) = -0.5
                    zCoeff(3) = 0.5
                endif
                u(indexX, indexY, indexZ,1) = ( yCoeff(1)*f(indexX,indexY,indexZ,3)+ yCoeff(2)*f(indexX,iy,indexZ,3) + yCoeff(3)*f(indexX,iy2,indexZ,3) &
                - ( zCoeff(1)*f(indexX,indexY,indexZ,2)+ zCoeff(2)*f(indexX,indexY,iz,2) + zCoeff(3)*f(indexX,indexY,iz2,2) ) )/hx
                u(indexX, indexY, indexZ,2) = ( zCoeff(1)*f(indexX,indexY,indexZ,1)+ zCoeff(2)*f(indexX,indexY,iz,1) + zCoeff(3)*f(indexX,indexY,iz2,1)&
                - ( xCoeff(1)*f(indexX,indexY,indexZ,3)+ xCoeff(2)*f(ix,indexY,indexZ,3) + xCoeff(3)*f(ix2,indexY,indexZ,3) ) )/hx
                u(indexX, indexY, indexZ,3) = ( xCoeff(1)*f(indexX,indexY,indexZ,2)+ xCoeff(2)*f(ix,indexY,indexZ,2) + xCoeff(3)*f(ix2,indexY,indexZ,2) &
                - ( yCoeff(1)*f(indexX,indexY,indexZ,1)+ yCoeff(2)*f(indexX,iy,indexZ,1) + yCoeff(3)*f(indexX,iy2,indexZ,1) ) )/hx

            enddo
        enddo
    enddo


    end subroutine computeVelocities


    subroutine computeLaplacian(vt, nx, ny, nz, hx, hy, hz, laplacian)
    double precision, dimension(:,:,:,:),allocatable :: vt, laplacian
    double precision hx, hy, hz
    integer nx,ny,nz, i, j, k, ix, ix2, iy, iy2, iz, iz2, indexX, indexY, indexZ


    ! Need to compute the laplacian just for inner nodes, and eventually if we have periodic conditions
    ! need to evaluate it in the periodic boundaries. 
    do ix=2,nx
        do iy=2,ny
            do iz=2,nz
                do i=1,3
                    laplacian(ix,iy,iz,i)= ( -2*vt(ix,iy,iz,i)+vt(ix-1,iy,iz,i)+vt(ix+1,iy,iz,i)  ) / hx**2 &
                                          +( -2*vt(ix,iy,iz,i)+vt(ix,iy-1,iz,i)+vt(ix,iy+1,iz,i)  ) / hx**2 &
                                          +( -2*vt(ix,iy,iz,i)+vt(ix,iy,iz-1,i)+vt(ix,iy,iz+1,i)  ) / hx**2 
                enddo
            enddo
        enddo
    enddo
    
    ! the computation over the inner nodes seems ok. I tested it with the function from Liu.

    end subroutine computeLaplacian
    
    subroutine computeConvectiveTerm(vt, u, conv, nx, ny, nz, hx, hy, hz)
    double precision, dimension(:,:,:,:),allocatable :: vt, u, conv
    double precision hx, hy, hz
    integer nx,ny,nz, i, j, k, ix, ix2, iy, iy2, iz, iz2, indexX, indexY, indexZ


    ! Compute the convective term for the inner nodes, eventually I need to do this for the periodic 
    ! boundaries
    do ix=2,nx
        do iy=2,ny
            do iz=2,nz
                !X component
                conv(ix,iy,iz,1) =  u(ix,iy,iz,2)  * 0.5 *( vt(ix,iy+1,iz,1) - vt(ix,iy-1,iz,1)  )/hx  &   !     U_2 * dW_1/dy 
								 + vt(ix,iy,iz,1)  * 0.5 *(  u(ix,iy+1,iz,2) -  u(ix,iy-1,iz,2)  )/hx  &   !   + W_1 * dU_2/dy 
                                - ( u(ix,iy,iz,1)  * 0.5 *( vt(ix,iy+1,iz,2) - vt(ix,iy-1,iz,2)  )/hx  &   !   -(U_1 * dW_2/dy 
                                +  vt(ix,iy,iz,2)  * 0.5 *(  u(ix,iy+1,iz,1) -  u(ix,iy-1,iz,1)  )/hx) &   !   + W_2 * dU_1/dy)
                                -(  u(ix,iy,iz,1)  * 0.5 *( vt(ix,iy,iz+1,3) - vt(ix,iy,iz-1,3)  )/hx  &   !   -(U_1 * dW_3/dz
				                +  vt(ix,iy,iz,3)  * 0.5 *(  u(ix,iy,iz+1,1) -  u(ix,iy,iz-1,1)  )/hx  &   !   + W_3 * dU_1/dz
                                - ( u(ix,iy,iz,3)  * 0.5 *( vt(ix,iy,iz+1,1) - vt(ix,iy,iz-1,1)  )/hx  &   !   -(U_3 * dW_1/dz 
                                +  vt(ix,iy,iz,1)  * 0.5 *(  u(ix,iy,iz+1,3) -  u(ix,iy,iz-1,3)  )/hx) )   !   + W_1 * dU_3/dz))
                
                !Y component
                conv(ix,iy,iz,2) =  u(ix,iy,iz,3)  * 0.5 *( vt(ix,iy,iz+1,2) - vt(ix,iy,iz-1,2)  )/hx  &   !     U_3 * dW_2/dz 
								 + vt(ix,iy,iz,2)  * 0.5 *(  u(ix,iy,iz+1,3) -  u(ix,iy,iz-1,3)  )/hx  &   !   + W_2 * dU_3/dz 
                                - ( u(ix,iy,iz,2)  * 0.5 *( vt(ix,iy,iz+1,3) - vt(ix,iy,iz-1,3)  )/hx  &   !   -(U_2 * dW_3/dz 
                                +  vt(ix,iy,iz,3)  * 0.5 *(  u(ix,iy,iz+1,2) -  u(ix,iy,iz-1,2)  )/hx) &   !   + W_3 * dU_2/dz)
                                -(  u(ix,iy,iz,2)  * 0.5 *( vt(ix+1,iy,iz,1) - vt(ix-1,iy,iz,1)  )/hx  &   !   -(U_2 * dW_1/dx
				                +  vt(ix,iy,iz,1)  * 0.5 *(  u(ix+1,iy,iz,2) -  u(ix-1,iy,iz,2)  )/hx  &   !   + W_1 * dU_2/dx
                                - ( u(ix,iy,iz,1)  * 0.5 *( vt(ix+1,iy,iz,2) - vt(ix-1,iy,iz,2)  )/hx  &   !   -(U_1 * dW_2/dx 
                                +  vt(ix,iy,iz,2)  * 0.5 *(  u(ix+1,iy,iz,1) -  u(ix-1,iy,iz,1)  )/hx) )   !   + W_2 * dU_1/dx))
                
                !Z component
                conv(ix,iy,iz,3) =  u(ix,iy,iz,1)  * 0.5 *( vt(ix+1,iy,iz,3) - vt(ix-1,iy,iz,3)  )/hx  &   !     U_1 * dW_3/dx 
								 + vt(ix,iy,iz,3)  * 0.5 *(  u(ix+1,iy,iz,1) -  u(ix-1,iy,iz,1)  )/hx  &   !   + W_3 * dU_1/dx 
                                - ( u(ix,iy,iz,3)  * 0.5 *( vt(ix+1,iy,iz,1) - vt(ix-1,iy,iz,1)  )/hx  &   !   -(U_3 * dW_1/dx 
                                +  vt(ix,iy,iz,1)  * 0.5 *(  u(ix+1,iy,iz,3) -  u(ix-1,iy,iz,3)  )/hx) &   !   + W_1 * dU_3/dx)
                                -(  u(ix,iy,iz,3)  * 0.5 *( vt(ix,iy+1,iz,2) - vt(ix,iy-1,iz,2)  )/hx  &   !   -(U_3 * dW_2/dy
				                +  vt(ix,iy,iz,2)  * 0.5 *(  u(ix,iy+1,iz,3) -  u(ix,iy-1,iz,3)  )/hx  &   !   + W_2 * dU_3/dy
                                - ( u(ix,iy,iz,2)  * 0.5 *( vt(ix,iy+1,iz,3) - vt(ix,iy-1,iz,3)  )/hx  &   !   -(U_2 * dW_3/dy 
                                +  vt(ix,iy,iz,3)  * 0.5 *(  u(ix,iy+1,iz,2) -  u(ix,iy-1,iz,2)  )/hx) )   !   + W_3 * dU_2/dy))
            
            enddo
        enddo
    enddo
    
    ! the computation over the inner nodes seems ok. I tested it with the function from Liu.

    end subroutine computeConvectiveTerm
    
    end module poissonSolveJacobi