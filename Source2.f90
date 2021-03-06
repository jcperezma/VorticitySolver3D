    include 'mkl_dfti.f90'
    include 'mkl_poisson.f90'
    module vorticitySolverUtilities
    use mkl_dfti
    use mkl_poisson
    implicit none 
    contains
    
     function crossProd(vec1, vec2)
    real(8), dimension (3):: crossProd, vec1, vec2
    crossProd(1)=  vec1(2)*vec2(3)-vec1(3)*vec2(2)
    crossProd(2)=-(vec1(1)*vec2(3)-vec1(3)*vec2(1))
    crossProd(3)=  vec1(1)*vec2(2)-vec1(2)*vec2(1)
     end function crossProd
    
     
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

    subroutine solvePoisson3DVorticityPer(f, vt, vtnew, nx, ny, nz, hx, maxIter, i, errorTol)

    double precision f(nx+1,ny+1,nz+1,3), vt(nx+1,ny+1,nz+1,3), vtnew(nx+1,ny+1,nz+1,3) ! 
    double precision errorTol, error, hx, fxterm , fzterm 
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

        ! update nodes at the periodic walls

        ! x=0 
        ! forward difference with respect of X
        ix =1
        do iy=2,ny ! the stream function is defined at the Y faces
            do iz=1,nz+1

                if (iz .eq. 1) then
                    fzterm =  f(ix,iy,2,i) +f(ix,iy,nz,i) 

                elseif(iz .eq. nz+1) then
                    fzterm =  f(ix,iy,2,i) +f(ix,iy,nz,i)  
                else
                    fzterm =  f(ix,iy,iz-1,i) +f(ix,iy,iz+1,i)
                endif

                f(ix,iy,iz,i)  = 1.5 * 0.16666 * ( hx*hx * vtnew(ix,iy,iz,i) + f(nx,iy,iz,i)+f(2,iy,iz,i) + f(ix,iy-1,iz,i)+f(ix,iy+1,iz,i)  &
                +  fzterm )+ (1 - 1.5) * f(ix,iy,iz,i)
            enddo
        enddo

        ! x=1 
        ix =nx+1
        do iy=2,ny ! the stream function is defined at the Y faces
            do iz=1,nz+1
                if (iz .eq. 1) then
                    fzterm =  f(ix,iy,2,i) +f(ix,iy,nz,i) 

                elseif(iz .eq. nz+1) then
                    fzterm =  f(ix,iy,2,i) +f(ix,iy,nz,i)  
                else
                    fzterm =  f(ix,iy,iz-1,i) +f(ix,iy,iz+1,i)
                endif

                f(ix,iy,iz,i)  = 1.5 * 0.16666 * ( hx*hx * vtnew(ix,iy,iz,i) + f(nx,iy,iz,i)+f(2,iy,iz,i) + f(ix,iy-1,iz,i)+f(ix,iy+1,iz,i)  &
                +  fzterm )+ (1 - 1.5) * f(ix,iy,iz,i)
            enddo
        enddo

        ! z = 0
        iz =1
        do ix=1,nx+1
            do iy=2,ny    
                if (ix .eq. 1) then
                    fxterm =  f(2,iy,iz,i) +f(nx,iy,iz,i) 

                elseif(ix .eq. nz+1) then
                    fxterm =  f(2,iy,iz,i) +f(nx,iy,iz,i)  
                else
                    fxterm =  f(ix-1,iy,iz,i) +f(ix+1,iy,iz,i)
                endif
                f(ix,iy,iz,i)  = 1.5 * 0.16666 * ( hx*hx * vtnew(ix,iy,iz,i) + fxterm + f(ix,iy-1,iz,i)+f(ix,iy+1,iz,i)  &
                + f(ix,iy,nz,i) +f(ix,iy,2,i)  )+ (1 - 1.5) * f(ix,iy,iz,i)
            enddo
        enddo 

        ! z = 1
        iz =nz+1
        do ix=1,nx+1
            do iy=2,ny    
                if (ix .eq. 1) then
                    fxterm =  f(2,iy,iz,i) +f(nx,iy,iz,i) 

                elseif(ix .eq. nx+1) then
                    fxterm =  f(2,iy,iz,i) +f(nx,iy,iz,i)  
                else
                    fxterm =  f(ix-1,iy,iz,i) +f(ix+1,iy,iz,i)
                endif
                f(ix,iy,iz,i)  = 1.5 * 0.16666 * ( hx*hx * vtnew(ix,iy,iz,i) + fxterm + f(ix,iy-1,iz,i)+f(ix,iy+1,iz,i)  &
                + f(ix,iy,nz,i) +f(ix,iy,2,i)  )+ (1 - 1.5) * f(ix,iy,iz,i)
            enddo
        enddo 



        ! update boundary nodes for the normal component with a neumann boundary condition df/dn =0  f(1) - f(-1) = 0





        SELECT CASE (i) ! select component
        CASE (5)
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
            do ix=1,nx+1
                do iz=1,nz+1    
                    if (ix .eq. 1) then
                    fxterm =  f(2,iy,iz,i) +f(nx,iy,iz,i) 

                elseif(ix .eq. nx+1) then
                    fxterm =  f(2,iy,iz,i) +f(nx,iy,iz,i)  
                else
                    fxterm =  f(ix-1,iy,iz,i) +f(ix+1,iy,iz,i)
                endif
                
                if (iz .eq. 1) then
                    fzterm =  f(ix,iy,nz,i) +f(ix,iy,2,i) 

                elseif(iz .eq. nz+1) then
                    fzterm =  f(ix,iy,nz,i) +f(ix,iy,2,i)  
                else
                    fzterm =  f(ix,iy,iz+1,i) +f(ix,iy,iz-1,i)
                endif
                
                    f(ix,iy,iz,i)  = 1.5 * 0.16666 * ( hx*hx * vtnew(ix,iy,iz,i) + fxterm + 2 * f(ix,iy+1,iz,i)  &
                    + fzterm )+ (1 - 1.5) * f(ix,iy,iz,i)
                enddo
            enddo 
            ! y = 1
            iy =ny+1
            do ix=1,nx+1
                do iz=1,nz+1
                    if (ix .eq. 1) then
                    fxterm =  f(2,iy,iz,i) +f(nx,iy,iz,i) 

                elseif(ix .eq. nx+1) then
                    fxterm =  f(2,iy,iz,i) +f(nx,iy,iz,i)  
                else
                    fxterm =  f(ix-1,iy,iz,i) +f(ix+1,iy,iz,i)
                endif
                
                if (iz .eq. 1) then
                    fzterm =  f(ix,iy,nz,i) +f(ix,iy,2,i) 

                elseif(iz .eq. nz+1) then
                    fzterm =  f(ix,iy,nz,i) +f(ix,iy,2,i)  
                else
                    fzterm =  f(ix,iy,iz-1,i) +f(ix,iy,iz+1,i)
                endif
                
                    f(ix,iy,iz,i)  = 1.5 * 0.16666 * ( hx*hx * vtnew(ix,iy,iz,i) + fxterm + 2 *f(ix,iy-1,iz,i)  &
                    + fzterm )+ (1 - 1.5) * f(ix,iy,iz,i)
                enddo
            enddo  

        CASE (6)
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



    end subroutine solvePoisson3DVorticityPer

    
    subroutine solvePoisson3DVorticityPerMKL(f, bd_ax, bd_bx , bd_ay, bd_by, bd_az, bd_bz, nx, ny, nz, ax, bx, ay, by, az, bz, i)
    double precision dpar(5*(nx+ny)/2+9)
    integer ipar(128), i, stat, nx, ny, nz
    double precision :: bd_ax(:,:), bd_bx(:,:) , bd_ay(:,:), bd_by(:,:), bd_az(:,:), bd_bz(:,:)
    double precision :: f(:,:,:)
    character(6) BCtype
    double precision q, ax, bx, ay, by, az, bz
    type(DFTI_DESCRIPTOR), pointer :: xhandle, yhandle
    
    SELECT CASE (i) ! select component
    CASE (1)
        BCtype = 'PPDDPP'
        bd_ay  =0
        bd_by  =0
    case (2)
        BCtype = 'PPNNPP'
        
    case (3)
        BCtype = 'PPDDPP'
        bd_ay  =-0.5 ! for unidirectional flow
        bd_by  =0.5
    end select
    q=0.0
    ipar =0
    ! solve the poisson eq 
    
    call D_INIT_HELMHOLTZ_3D(ax, bx, ay, by, az, bz, nx, ny, nz, BCtype , q, ipar, dpar, stat)
    ipar(3)=0
    !print *, " warning 1 "
    !if (stat.ne.0) goto 999
    call D_COMMIT_HELMHOLTZ_3D (f, bd_ax, bd_bx, bd_ay, bd_by, bd_az , bd_bz, xhandle, yhandle, ipar, dpar, stat)
    !print *, " warning 2 "
    !if (stat.ne.0) goto 999
    call D_HELMHOLTZ_3D(f, bd_ax, bd_bx, bd_ay, bd_by, bd_az , bd_bz, xhandle, yhandle, ipar, dpar, stat)
    !print *, " warning 3 "
    !if (stat.ne.0) goto 999
    call free_Helmholtz_3D(xhandle, yhandle, ipar, stat)
    !print *, " warning 4 "
    !if (stat.ne.0) goto 999
    
    ! need boundary conditions and right hand side, those are inputs
    
  !999 print *, 'Double precision 2D Poisson example FAILED to compute the solution...'  
    
    end subroutine solvePoisson3DVorticityPerMKL
    
    subroutine updateBoundary(f, vt, u, nx, ny, nz, hx)
    double precision f(nx+1,ny+1,nz+1,3), vt(nx+1,ny+1,nz+1,3),u(nx+1,ny+1,nz+1,3) ! 
    double precision errorTol, error, hx, dfdx, dfdy, dfdz
    integer nx,ny,nz, maxIter, i, iter, ix, iy, iz
    ! use U and f to compute the vorticity at the boundaries, 
    ! U here is the velocity at the boundaries, I am wasting tons of memory but I shouldnt worry about memory for now.

    !!!!!!!!!!!!!!!!!!!!!!!!!!!! For faces Z =const !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! for face Z =0
    ! normal component for this face is 3
    ! for inner nodes   
    iz = 1
    do iy=1,ny+1
        do ix=1,nx+1
            if (iy .eq. 1)then
                !Forward difference
                dfdy = (u(ix,iy+1,iz,1) -  u(ix,iy,iz,1))/hx
            elseif (iy .eq. ny+1) then
                ! Backward Difference
                dfdy = (u(ix,iy,iz,1) -  u(ix,iy-1,iz,1))/hx
            else
                ! Central Difference  
                dfdy = 0.5*(u(ix,iy+1,iz,1) -  u(ix,iy-1,iz,1))/hx
            endif


            if (ix .eq. 1)then
                !Forward difference
                dfdx = (u(ix+1,iy,iz,2) -  u(ix,iy,iz,2))/hx
            elseif (ix .eq. nx+1) then
                ! Backward Difference
                dfdx = (u(ix,iy,iz,2) -  u(ix-1,iy,iz,2))/hx
            else
                ! Central Difference  
                dfdx = 0.5*(u(ix+1,iy,iz,2) -  u(ix-1,iy,iz,2))/hx
            endif

            vt(ix,iy,iz,3) = dfdx-dfdy ! du_2/dx - du_1/dy
        enddo
    enddo

    ! tangential components are 1 and 2
    ! Component 1
    ! inner nodes
    do iy=1,ny+1
        do ix=1,nx+1
            if (ix .eq. 1)then
                !Forward difference
                dfdx = (f(ix+1,iy,iz,3) -  f(ix,iy,iz,3))/hx
            elseif (ix .eq. nx+1) then
                ! Backward Difference
                dfdx = (f(ix,iy,iz,3) -  f(ix-1,iy,iz,3))/hx
            else
                ! Central Difference  
                dfdx = 0.5*(f(ix+1,iy,iz,3) -  f(ix-1,iy,iz,3))/hx
            endif

            !vt(ix,iy,iz,1) = -2* f(ix,iy,iz+1,1) / hx**2 + 2/hx* ( (f(ix+1,iy,iz,3) - f(ix-1,iy,iz,3))/(2*hx) + u(ix,iy,iz,2))
            vt(ix,iy,iz,1) = -2* f(ix,iy,iz+1,1) / hx**2 + 2/hx* ( dfdx + u(ix,iy,iz,2))
        enddo
    enddo

    ! Component 2
    ! all nodes
    do iy=1,ny+1
        do ix=1,nx+1
            if (iy .eq. 1)then
                !Forward difference
                dfdy = (f(ix,iy+1,iz,3) -  f(ix,iy,iz,3))/hx
            elseif (iy .eq. ny+1) then
                ! Backward Difference
                dfdy = (f(ix,iy,iz,3) -  f(ix,iy-1,iz,3))/hx
            else
                ! Central Difference  
                dfdy = 0.5*(f(ix,iy+1,iz,3) -  f(ix,iy-1,iz,3))/hx
            endif

            !vt(ix,iy,iz,1) = -2* f(ix,iy,iz+1,2) / hx**2 + 2/hx* ( (f(ix,iy+1,iz,3) - f(ix,iy-1,iz,3))/(2*hx) - u(ix,iy,iz,1))
            vt(ix,iy,iz,1) = -2* f(ix,iy,iz+1,2) / hx**2 + 2/hx* ( dfdy - u(ix,iy,iz,1))
        enddo
    enddo


    ! for face Z =1
    ! normal component for this face is 3
    ! for inner nodes   
    iz = nz+1
    do iy=1,ny+1
        do ix=1,nx+1
            if (iy .eq. 1)then
                !Forward difference
                dfdy = (u(ix,iy+1,iz,1) -  u(ix,iy,iz,1))/hx
            elseif (iy .eq. ny+1) then
                ! Backward Difference
                dfdy = (u(ix,iy,iz,1) -  u(ix,iy-1,iz,1))/hx
            else
                ! Central Difference  
                dfdy = 0.5*(u(ix,iy+1,iz,1) -  u(ix,iy-1,iz,1))/hx
            endif


            if (ix .eq. 1)then
                !Forward difference
                dfdx = (u(ix+1,iy,iz,2) -  u(ix,iy,iz,2))/hx
            elseif (ix .eq. nx+1) then
                ! Backward Difference
                dfdx = (u(ix,iy,iz,2) -  u(ix-1,iy,iz,2))/hx
            else
                ! Central Difference  
                dfdx = 0.5*(u(ix+1,iy,iz,2) -  u(ix-1,iy,iz,2))/hx
            endif

            vt(ix,iy,iz,3) = dfdx-dfdy ! du_2/dx - du_1/dy
        enddo
    enddo

    ! tangential components are 1 and 2
    ! Component 1
    ! inner nodes
    do iy=1,ny+1
        do ix=1,nx+1
            if (ix .eq. 1)then
                !Forward difference
                dfdx = (f(ix+1,iy,iz,3) -  f(ix,iy,iz,3))/hx
            elseif (ix .eq. nx+1) then
                ! Backward Difference
                dfdx = (f(ix,iy,iz,3) -  f(ix-1,iy,iz,3))/hx
            else
                ! Central Difference  
                dfdx = 0.5*(f(ix+1,iy,iz,3) -  f(ix-1,iy,iz,3))/hx
            endif

            !vt(ix,iy,iz,1) = -2* f(ix,iy,iz-1,1) / hx**2 - 2/hx* ( (f(ix+1,iy,iz,3) - f(ix-1,iy,iz,3))/(2*hx) + u(ix,iy,iz,2))
            vt(ix,iy,iz,1) = -2* f(ix,iy,iz-1,1) / hx**2 - 2/hx* ( dfdx + u(ix,iy,iz,2))
        enddo
    enddo

    ! Component 2
    ! inner nodes
    do iy=1,ny+1
        do ix=1,nx+1
            if (iy .eq. 1)then
                !Forward difference
                dfdy = (f(ix,iy+1,iz,3) -  f(ix,iy,iz,3))/hx
            elseif (iy .eq. ny+1) then
                ! Backward Difference
                dfdy = (f(ix,iy,iz,3) -  f(ix,iy-1,iz,3))/hx
            else
                ! Central Difference  
                dfdy = 0.5*(f(ix,iy+1,iz,3) -  f(ix,iy-1,iz,3))/hx
            endif
            !vt(ix,iy,iz,1) = -2* f(ix,iy,iz-1,2) / hx**2 - 2/hx* ( (f(ix,iy+1,iz,3) - f(ix,iy-1,iz,3))/(2*hx) - u(ix,iy,iz,1))
            vt(ix,iy,iz,1) = -2* f(ix,iy,iz-1,2) / hx**2 - 2/hx* ( dfdy - u(ix,iy,iz,1))       
        enddo
    enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!! For faces X =const !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! for face X =0
    ! normal component for this face is 1
    ! for all nodes   
    ix = 1
    do iy=1,ny+1
        do iz=1,nz+1
            if (iy .eq. 1)then
                !Forward difference
                dfdy = (u(ix,iy+1,iz,3) -  u(ix,iy,iz,3))/hx
            elseif (iy .eq. ny+1) then
                ! Backward Difference
                dfdy = (u(ix,iy,iz,3) -  u(ix,iy-1,iz,3))/hx
            else
                ! Central Difference  
                dfdy = 0.5*(u(ix,iy+1,iz,3) -  u(ix,iy-1,iz,3))/hx
            endif


            if (iz .eq. 1)then
                !Forward difference
                dfdz = (u(ix,iy,iz+1,2) -  u(ix,iy,iz,2))/hx
            elseif (iz .eq. nz+1) then
                ! Backward Difference
                dfdz = (u(ix,iy,iz,2) -  u(ix,iy,iz-1,2))/hx
            else
                ! Central Difference  
                dfdz = 0.5*(u(ix,iy,iz+1,2) -  u(ix,iy,iz-1,2))/hx
            endif

            vt(ix,iy,iz,1) = dfdy -dfdz ! du_3/dy - du_2/dz 
            !print '(A5, f6.2, A5, f6.2, A5, f6.2 )', 'dfdy', dfdy, 'dfdz', dfdz, ' vt ', vt(ix,iy,iz,1)

        enddo
    enddo

    ! tangential components are 2 and 3
    ! Component 2
    ! inner nodes
    do iy=1,ny+1
        do iz=1,nz+1
            if (iy .eq. 1)then
                !Forward difference
                dfdy = (f(ix,iy+1,iz,1) -  f(ix,iy,iz,1))/hx
            elseif (iy .eq. ny+1) then
                ! Backward Difference
                dfdy = (f(ix,iy,iz,1) -  f(ix,iy-1,iz,1))/hx
            else
                ! Central Difference  
                dfdy = 0.5*(f(ix,iy+1,iz,1) -  f(ix,iy-1,iz,1))/hx
            endif
            !vt(ix,iy,iz,2) = -2* f(ix+1,iy,iz,2) / hx**2 + 2/hx* ( (f(ix,iy+1,iz,1) - f(ix,iy-1,iz,1))/(2*hx) + u(ix,iy,iz,3))
            vt(ix,iy,iz,2) = -2* f(ix+1,iy,iz,2) / hx**2 + 2/hx* ( dfdy + u(ix,iy,iz,3))
        enddo
    enddo

    ! Component 3
    ! inner nodes
    do iy=1,ny+1
        do iz=1,nz+1
            if (iz .eq. 1)then
                !Forward difference
                dfdz = (f(ix,iy,iz+1,1) -  f(ix,iy,iz,1))/hx
            elseif (iz .eq. nz+1) then
                ! Backward Difference
                dfdz = (f(ix,iy,iz,1) -  f(ix,iy,iz-1,1))/hx
            else
                ! Central Difference  
                dfdz = 0.5*(f(ix,iy,iz+1,1) -  f(ix,iy,iz-1,1))/hx
            endif

            !vt(ix,iy,iz,3) = -2* f(ix+1,iy,iz,3) / hx**2 + 2/hx* ( (f(ix,iy,iz+1,1) - f(ix,iy,iz-1,1))/(2*hx) - u(ix,iy,iz,2))
            vt(ix,iy,iz,3) = -2* f(ix+1,iy,iz,3) / hx**2 + 2/hx* ( dfdz - u(ix,iy,iz,2))
        enddo
    enddo


    ! for face X =1
    ! normal component for this face is 3
    ! for inner nodes   
    ix = nx+1
    do iy=1,ny+1
        do iz=1,nz+1
            if (iy .eq. 1)then
                !Forward difference
                dfdy = (u(ix,iy+1,iz,3) -  u(ix,iy,iz,3))/hx
            elseif (iy .eq. ny+1) then
                ! Backward Difference
                dfdy = (u(ix,iy,iz,3) -  u(ix,iy-1,iz,3))/hx
            else
                ! Central Difference  
                dfdy = 0.5*(u(ix,iy+1,iz,3) -  u(ix,iy-1,iz,3))/hx
            endif


            if (iz .eq. 1)then
                !Forward difference
                dfdz = (u(ix,iy,iz+1,2) -  u(ix,iy,iz,2))/hx
            elseif (iz .eq. nz+1) then
                ! Backward Difference
                dfdz = (u(ix,iy,iz,2) -  u(ix,iy,iz-1,2))/hx
            else
                ! Central Difference  
                dfdz = 0.5*(u(ix,iy,iz+1,2) -  u(ix,iy,iz-1,2))/hx
            endif

            vt(ix,iy,iz,1) = dfdy -dfdz ! du_3/dy - du_2/dz 
        enddo
    enddo

    ! tangential components are 2 and 3
    ! Component 2
    ! inner nodes
    do iy=1,ny+1
        do iz=1,nz+1
            if (iy .eq. 1)then
                !Forward difference
                dfdy = (f(ix,iy+1,iz,1) -  f(ix,iy,iz,1))/hx
            elseif (iy .eq. ny+1) then
                ! Backward Difference
                dfdy = (f(ix,iy,iz,1) -  f(ix,iy-1,iz,1))/hx
            else
                ! Central Difference  
                dfdy = 0.5*(f(ix,iy+1,iz,1) -  f(ix,iy-1,iz,1))/hx
            endif

            !vt(ix,iy,iz,2) = -2* f(ix-1,iy,iz,2) / hx**2 - 2/hx* ( (f(ix,iy+1,iz,1) - f(ix,iy-1,iz,1))/(2*hx) + u(ix,iy,iz,3))
            vt(ix,iy,iz,2) = -2* f(ix-1,iy,iz,2) / hx**2 - 2/hx* ( dfdy + u(ix,iy,iz,3))
        enddo
    enddo

    ! Component 3
    ! inner nodes
    do iy=1,ny+1
        do iz=1,nz+1
            if (iz .eq. 1)then
                !Forward difference
                dfdz = (f(ix,iy,iz+1,1) -  f(ix,iy,iz,1))/hx
            elseif (iz .eq. nz+1) then
                ! Backward Difference
                dfdz = (f(ix,iy,iz,1) -  f(ix,iy,iz-1,1))/hx
            else
                ! Central Difference  
                dfdz = 0.5*(f(ix,iy,iz+1,1) -  f(ix,iy,iz-1,1))/hx
            endif

            !vt(ix,iy,iz,3) = -2* f(ix-1,iy,iz,3) / hx**2 - 2/hx* ( (f(ix,iy,iz+1,1) - f(ix,iy,iz-1,1))/(2*hx) - u(ix,iy,iz,2))
            vt(ix,iy,iz,3) = -2* f(ix-1,iy,iz,3) / hx**2 - 2/hx* ( dfdz - u(ix,iy,iz,2))
        enddo
    enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!! For faces Y =const !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! for face Y =0
    ! normal component for this face is 2
    ! for inner nodes   
    iy = 1
    do ix=1,nx+1
        do iz=1,nz+1
            if (ix .eq. 1)then
                !Forward difference
                dfdx = (u(ix+1,iy,iz,3) -  u(ix,iy,iz,3))/hx
            elseif (ix .eq. nx+1) then
                ! Backward Difference
                dfdx = (u(ix,iy,iz,3) -  u(ix-1,iy,iz,3))/hx
            else
                ! Central Difference  
                dfdx = 0.5*(u(ix+1,iy,iz,3) -  u(ix-1,iy,iz,3))/hx
            endif


            if (iz .eq. 1)then
                !Forward difference
                dfdz = (u(ix,iy,iz+1,1) -  u(ix,iy,iz,1))/hx
            elseif (iz .eq. nz+1) then
                ! Backward Difference
                dfdz = (u(ix,iy,iz,1) -  u(ix,iy,iz-1,1))/hx
            else
                ! Central Difference  
                dfdz = 0.5*(u(ix,iy,iz+1,1) -  u(ix,iy,iz-1,1))/hx
            endif

            vt(ix,iy,iz,2) = dfdz - dfdx !du_1/dz -du_3/dx 0.5/hx* ( u(ix,iy,iz+1,1) - u(ix,iy,iz-1,1) - u(ix+1,iy,iz,3) + u(ix-1,iy,iz,3) )

        enddo
    enddo

    ! tangential components are 1 and 3
    ! Component 1
    ! inner nodes
    do ix=1,nx+1
        do iz=1,nz+1
            if (ix .eq. 1)then
                !Forward difference
                dfdx = (f(ix+1,iy,iz,2) -  f(ix,iy,iz,2))/hx
            elseif (ix .eq. nx+1) then
                ! Backward Difference
                dfdx = (f(ix,iy,iz,2) -  f(ix-1,iy,iz,2))/hx
            else
                ! Central Difference  
                dfdx = 0.5*(f(ix+1,iy,iz,2) -  f(ix-1,iy,iz,2))/hx
            endif

            !vt(ix,iy,iz,1) = -2* f(ix,iy+1,iz,1) / hx**2 + 2/hx* ( (f(ix+1,iy,iz,2) - f(ix-1,iy,iz,2))/(2*hx) - u(ix,iy,iz,3))
            vt(ix,iy,iz,1) = -2* f(ix,iy+1,iz,1) / hx**2 + 2/hx* ( dfdx - u(ix,iy,iz,3))
            !print '(A5, f6.2, A5, f6.2, A5, f6.2 )', 'u ', u(ix,iy,iz,3), 'dfdx', hx, ' vt ', vt(ix,iy,iz,1)
        enddo
    enddo

    ! Component 3
    ! all nodes
    do ix=1,nx+1
        do iz=1,nz+1
            if (iz .eq. 1)then
                !Forward difference
                dfdz = (f(ix,iy,iz+1,1) -  f(ix,iy,iz,1))/hx
            elseif (iz .eq. nz+1) then
                ! Backward Difference
                dfdz = (f(ix,iy,iz,1) -  f(ix,iy,iz-1,1))/hx
            else
                ! Central Difference  
                dfdz = 0.5*(f(ix,iy,iz+1,1) -  f(ix,iy,iz-1,1))/hx
            endif

            !vt(ix,iy,iz,3) = -2* f(ix,iy+1,iz,3) / hx**2 + 2/hx* ( (f(ix,iy,iz+1,1) - f(ix,iy,iz-1,1))/(2*hx) + u(ix,iy,iz,1))
            vt(ix,iy,iz,3) = -2* f(ix,iy+1,iz,3) / hx**2 + 2/hx* ( dfdz + u(ix,iy,iz,1))
        enddo
    enddo

    ! for face Y =1
    ! normal component for this face is 2
    ! for inner nodes   
    iy = ny+1
    do ix=1,nx+1
        do iz=1,nz+1
            if (ix .eq. 1)then
                !Forward difference
                dfdx = (u(ix+1,iy,iz,3) -  u(ix,iy,iz,3))/hx
            elseif (ix .eq. nx+1) then
                ! Backward Difference
                dfdx = (u(ix,iy,iz,3) -  u(ix-1,iy,iz,3))/hx
            else
                ! Central Difference  
                dfdx = 0.5*(u(ix+1,iy,iz,3) -  u(ix-1,iy,iz,3))/hx
            endif
            if (iz .eq. 1)then
                !Forward difference
                dfdz = (u(ix,iy,iz+1,1) -  u(ix,iy,iz,1))/hx
            elseif (iz .eq. nz+1) then
                ! Backward Difference
                dfdz = (u(ix,iy,iz,1) -  u(ix,iy,iz-1,1))/hx
            else
                ! Central Difference  
                dfdz = 0.5*(u(ix,iy,iz+1,1) -  u(ix,iy,iz-1,1))/hx
            endif

            vt(ix,iy,iz,2) = dfdz - dfdx !du_1/dz -du_3/dx 0.5/hx* ( u(ix,iy,iz+1,1) - u(ix,iy,iz-1,1) - u(ix+1,iy,iz,3) + u(ix-1,iy,iz,3) )
        enddo
    enddo

    ! tangential components are 1 and 3
    ! Component 1
    ! inner nodes
    do ix=1,nx+1
        do iz=1,nz+1
            if (ix .eq. 1)then
                !Forward difference
                dfdx = (f(ix+1,iy,iz,2) -  f(ix,iy,iz,2))/hx
            elseif (ix .eq. nx+1) then
                ! Backward Difference
                dfdx = (f(ix,iy,iz,2) -  f(ix-1,iy,iz,2))/hx
            else
                ! Central Difference  
                dfdx = 0.5*(f(ix+1,iy,iz,2) -  f(ix-1,iy,iz,2))/hx
            endif

            !vt(ix,iy,iz,1) = -2* f(ix,iy-1,iz,1) / hx**2 - 2/hx* ( (f(ix+1,iy,iz,2) - f(ix-1,iy,iz,2))/(2*hx) - u(ix,iy,iz,3))
            vt(ix,iy,iz,1) = -2* f(ix,iy-1,iz,1) / hx**2 - 2/hx* ( dfdx - u(ix,iy,iz,3))
        enddo
    enddo

    ! Component 3
    ! inner nodes
    do ix=1,nx+1
        do iz=1,nz+1
            if (iz .eq. 1)then
                !Forward difference
                dfdz = (f(ix,iy,iz+1,2) -  f(ix,iy,iz,2))/hx
            elseif (iz .eq. nz+1) then
                ! Backward Difference
                dfdz = (f(ix,iy,iz,2) -  f(ix,iy,iz-1,2))/hx
            else
                ! Central Difference  
                dfdz = 0.5*(f(ix,iy,iz+1,2) -  f(ix,iy,iz-1,2))/hx
            endif

            !vt(ix,iy,iz,3) = -2* f(ix,iy-1,iz,3) / hx**2 - 2/hx* ( (f(ix,iy,iz+1,2) - f(ix,iy,iz-1,2))/(2*hx) + u(ix,iy,iz,1))
            vt(ix,iy,iz,3) = -2* f(ix,iy-1,iz,3) / hx**2 - 2/hx* ( dfdz + u(ix,iy,iz,1))
        enddo
    enddo

    end subroutine updateBoundary

    subroutine updateBoundaryPer(f, vt, u, nx, ny, nz, hx)
    double precision f(nx+1,ny+1,nz+1,3), vt(nx+1,ny+1,nz+1,3),u(nx+1,ny+1,nz+1,3) ! 
    double precision errorTol, error, hx, dfdx, dfdy, dfdz
    integer nx,ny,nz, maxIter, i, iter, ix, iy, iz
    ! use U and f to compute the vorticity at the boundaries, 
    ! U here is the velocity at the boundaries, I am wasting tons of memory but I shouldnt worry about memory for now.

    ! This is for the case where the only boundary is in the Y faces. There is periodicity in the X and Z direction. 




    !!!!!!!!!!!!!!!!!!!!!!!!!!!! For faces Y =const !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! for face Y =0
    ! normal component for this face is 2
    ! for inner nodes   
    iy = 1
    do ix=1,nx+1
        do iz=1,nz+1
            if (ix .eq. 1)then
                !Periodic
                dfdx = 0.5*(u(2,iy,iz,3) -  u(nx,iy,iz,3))/hx
            elseif (ix .eq. nx+1) then
                ! Periodic
                dfdx = 0.5*(u(2,iy,iz,3) -  u(nx,iy,iz,3))/hx
            else
                ! Central Difference  
                dfdx = 0.5*(u(ix+1,iy,iz,3) -  u(ix-1,iy,iz,3))/hx
            endif


            if (iz .eq. 1)then
                !Periodic
                dfdz = 0.5*(u(ix,iy,2,1) -  u(ix,iy,nz,1))/hx
            elseif (iz .eq. nz+1) then
                !Periodic
                dfdz = 0.5*(u(ix,iy,2,1) -  u(ix,iy,nz,1))/hx
            else
                ! Central Difference  
                dfdz = 0.5*(u(ix,iy,iz+1,1) -  u(ix,iy,iz-1,1))/hx
            endif

            vt(ix,iy,iz,2) = dfdz - dfdx !du_1/dz -du_3/dx 0.5/hx* ( u(ix,iy,iz+1,1) - u(ix,iy,iz-1,1) - u(ix+1,iy,iz,3) + u(ix-1,iy,iz,3) )

        enddo
    enddo

    ! tangential components are 1 and 3
    ! Component 1
    ! inner nodes
    do ix=1,nx+1
        do iz=1,nz+1
            if (ix .eq. 1)then
                !Periodic
                dfdx = 0.5*(f(2,iy,iz,2) -  f(nx,iy,iz,2))/hx
            elseif (ix .eq. nx+1) then
                !Periodic
                dfdx = 0.5*(f(2,iy,iz,2) -  f(nx,iy,iz,2))/hx
            else
                ! Central Difference  
                dfdx = 0.5*(f(ix+1,iy,iz,2) -  f(ix-1,iy,iz,2))/hx
            endif

            !vt(ix,iy,iz,1) = -2* f(ix,iy+1,iz,1) / hx**2 + 2/hx* ( (f(ix+1,iy,iz,2) - f(ix-1,iy,iz,2))/(2*hx) - u(ix,iy,iz,3))
            vt(ix,iy,iz,1) = -2* f(ix,iy+1,iz,1) / hx**2 + 2/hx* ( dfdx - u(ix,iy,iz,3))
            !print '(A5, f6.2, A5, f6.2, A5, f6.2 )', 'u ', u(ix,iy,iz,3), 'dfdx', hx, ' vt ', vt(ix,iy,iz,1)
        enddo
    enddo

    ! Component 3
    ! all nodes
    do ix=1,nx+1
        do iz=1,nz+1
            if (iz .eq. 1)then
                !Periodic
                dfdz =  0.5*(f(ix,iy,2,1) -  f(ix,iy,nz,1))/hx
            elseif (iz .eq. nz+1) then
                !Periodic
                dfdz =  0.5*(f(ix,iy,2,1) -  f(ix,iy,nz,1))/hx
            else
                ! Central Difference  
                dfdz = 0.5*(f(ix,iy,iz+1,1) -  f(ix,iy,iz-1,1))/hx
            endif

            !vt(ix,iy,iz,3) = -2* f(ix,iy+1,iz,3) / hx**2 + 2/hx* ( (f(ix,iy,iz+1,1) - f(ix,iy,iz-1,1))/(2*hx) + u(ix,iy,iz,1))
            vt(ix,iy,iz,3) = -2* (f(ix,iy+1,iz,3) -f(ix,iy,iz,3) ) / hx**2 + 2/hx* ( dfdz + u(ix,iy,iz,1))
        enddo
    enddo

    ! for face Y =1
    ! normal component for this face is 2
    ! for inner nodes   
    iy = ny+1
    do ix=1,nx+1
        do iz=1,nz+1
            if (ix .eq. 1)then
                !Periodic
                dfdx = 0.5*(u(2,iy,iz,3) -  u(nx,iy,iz,3))/hx
            elseif (ix .eq. nx+1) then
                ! Periodic
                dfdx = 0.5*(u(2,iy,iz,3) -  u(nx,iy,iz,3))/hx
            else
                ! Central Difference  
                dfdx = 0.5*(u(ix+1,iy,iz,3) -  u(ix-1,iy,iz,3))/hx
            endif


            if (iz .eq. 1)then
                !Periodic
                dfdz = 0.5*(u(ix,iy,2,1) -  u(ix,iy,nz,1))/hx
            elseif (iz .eq. nz+1) then
                !Periodic
                dfdz = 0.5*(u(ix,iy,2,1) -  u(ix,iy,nz,1))/hx
            else
                ! Central Difference  
                dfdz = 0.5*(u(ix,iy,iz+1,1) -  u(ix,iy,iz-1,1))/hx
            endif

            vt(ix,iy,iz,2) = dfdz - dfdx !du_1/dz -du_3/dx 0.5/hx* ( u(ix,iy,iz+1,1) - u(ix,iy,iz-1,1) - u(ix+1,iy,iz,3) + u(ix-1,iy,iz,3) )
        enddo
    enddo

    ! tangential components are 1 and 3
    ! Component 1
    ! inner nodes
    do ix=1,nx+1
        do iz=1,nz+1
            if (ix .eq. 1)then
                !Periodic
                dfdx = 0.5*(f(2,iy,iz,2) -  f(nx,iy,iz,2))/hx
            elseif (ix .eq. nx+1) then
                !Periodic
                dfdx = 0.5*(f(2,iy,iz,2) -  f(nx,iy,iz,2))/hx
            else
                ! Central Difference  
                dfdx = 0.5*(f(ix+1,iy,iz,2) -  f(ix-1,iy,iz,2))/hx
            endif

            !vt(ix,iy,iz,1) = -2* f(ix,iy-1,iz,1) / hx**2 - 2/hx* ( (f(ix+1,iy,iz,2) - f(ix-1,iy,iz,2))/(2*hx) - u(ix,iy,iz,3))
            vt(ix,iy,iz,1) = -2* (f(ix,iy-1,iz,1) ) / hx**2 - 2/hx* ( dfdx - u(ix,iy,iz,3))
        enddo
    enddo

    ! Component 3
    ! inner nodes
    do ix=1,nx+1
        do iz=1,nz+1
            if (iz .eq. 1)then
                !Periodic
                dfdz = 0.5*(f(ix,iy,2,2) -  f(ix,iy,nz,2))/hx
            elseif (iz .eq. nz+1) then
                !Periodic
                dfdz = 0.5*(f(ix,iy,2,2) -  f(ix,iy,nz,2))/hx
            else
                ! Central Difference  
                dfdz = 0.5*(f(ix,iy,iz+1,2) -  f(ix,iy,iz-1,2))/hx
            endif

            !vt(ix,iy,iz,3) = -2* f(ix,iy-1,iz,3) / hx**2 - 2/hx* ( (f(ix,iy,iz+1,2) - f(ix,iy,iz-1,2))/(2*hx) + u(ix,iy,iz,1))
            vt(ix,iy,iz,3) = -2* (f(ix,iy-1,iz,3) - f(ix,iy,iz,3)   ) / hx**2 - 2/hx* ( dfdz + u(ix,iy,iz,1))
        enddo
    enddo

    end subroutine updateBoundaryPer

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

    subroutine updateVelBCs(ub, nx, ny, nz, hx, hy, hz, lx, ly, lz, t )
    implicit none
    double precision, dimension(:,:,:,:),allocatable ::  ub
    !double precision f(nx+1,ny+1,nz+1,3), u(nx+1,ny+1,nz+1,3)! 
    double precision hx, hy, hz, lx, ly, lz, t, xi, yi, zi,pi
    integer nx,ny,nz, i, j, k, ix, ix2, iy, iy2, iz, iz2, indexX, indexY, indexZ
    ! Computes the components of the velocity from the stream function derivatives
    ! assumes a uniform grid; hx=hy=hz, though to change this is really easy.
    ! I tested it with simple functions and with the functions used by liu the max error is 0.25%
    pi=3.14159265358979324D0
    ix=1
    do iy = 1,ny+1
        do iz =1 ,nz+1
            xi=hx*(ix-1)/lx
            yi=hy*(iy-1)/ly
            zi=hz*(iz-1)/lz
            ub(ix,iy,iz,1)=exp(t)*pi*cos(pi*yi)*sin(pi*xi)-sin(pi*xi)*Pi*cos(pi*zi)
            ub(ix,iy,iz,2)=exp(t)*sin(pi*yi)*pi*cos(pi*zi)-sin(pi*yi)*pi*cos(pi*xi)
            ub(ix,iy,iz,3)=exp(t)*pi*cos(pi*xi)*sin(pi*zi)-Pi*cos(pi*yi)*sin(pi*zi)
        enddo
    enddo

    ! set boundary bd_bx at x=1
    ix=nx+1
    do iy = 1,ny+1
        do iz =1 ,nz+1
            xi=hx*(ix-1)/lx
            yi=hy*(iy-1)/ly
            zi=hz*(iz-1)/lz
            ub(ix,iy,iz,1)=exp(t)*pi*cos(pi*yi)*sin(pi*xi)-sin(pi*xi)*Pi*cos(pi*zi)
            ub(ix,iy,iz,2)=exp(t)*sin(pi*yi)*pi*cos(pi*zi)-sin(pi*yi)*pi*cos(pi*xi)
            ub(ix,iy,iz,3)=exp(t)*pi*cos(pi*xi)*sin(pi*zi)-Pi*cos(pi*yi)*sin(pi*zi)
        enddo
    enddo

    ! set boundary bd_ay at y=0
    iy=1
    do ix = 1,nx+1
        do iz =1 ,nz+1
            xi=hx*(ix-1)/lx
            yi=hy*(iy-1)/ly
            zi=hz*(iz-1)/lz
            ub(ix,iy,iz,1)=exp(t)*pi*cos(pi*yi)*sin(pi*xi)-sin(pi*xi)*Pi*cos(pi*zi)
            ub(ix,iy,iz,2)=exp(t)*sin(pi*yi)*pi*cos(pi*zi)-sin(pi*yi)*pi*cos(pi*xi)
            ub(ix,iy,iz,3)=exp(t)*pi*cos(pi*xi)*sin(pi*zi)-Pi*cos(pi*yi)*sin(pi*zi)
        enddo
    enddo

    ! set boundary bd_by at y=1
    iy =ny+1
    do ix = 1,nx+1
        do iz =1 ,nz+1

            xi=hx*(ix-1)/lx
            yi=hy*(iy-1)/ly
            zi=hz*(iz-1)/lz
            ub(ix,iy,iz,1)=exp(t)*pi*cos(pi*yi)*sin(pi*xi)-sin(pi*xi)*Pi*cos(pi*zi)
            ub(ix,iy,iz,2)=exp(t)*sin(pi*yi)*pi*cos(pi*zi)-sin(pi*yi)*pi*cos(pi*xi)
            ub(ix,iy,iz,3)=exp(t)*pi*cos(pi*xi)*sin(pi*zi)-Pi*cos(pi*yi)*sin(pi*zi)
        enddo
    enddo

    ! set boundary bd_ay at z=0
    iz=1
    do ix = 1,nx+1
        do iy =1 ,ny+1

            xi=hx*(ix-1)/lx
            yi=hy*(iy-1)/ly
            zi=hz*(iz-1)/lz
            ub(ix,iy,iz,1)=exp(t)*pi*cos(pi*yi)*sin(pi*xi)-sin(pi*xi)*Pi*cos(pi*zi)
            ub(ix,iy,iz,2)=exp(t)*sin(pi*yi)*pi*cos(pi*zi)-sin(pi*yi)*pi*cos(pi*xi)
            ub(ix,iy,iz,3)=exp(t)*pi*cos(pi*xi)*sin(pi*zi)-Pi*cos(pi*yi)*sin(pi*zi)
        enddo
    enddo

    ! set boundary bd_by at z=1
    iz = nz+1
    do ix = 1,nx+1
        do iy =1 ,ny+1

            xi=hx*(ix-1)/lx
            yi=hy*(iy-1)/ly
            zi=hz*(iz-1)/lz
            ub(ix,iy,iz,1)=exp(t)*pi*cos(pi*yi)*sin(pi*xi)-sin(pi*xi)*Pi*cos(pi*zi)
            ub(ix,iy,iz,2)=exp(t)*sin(pi*yi)*pi*cos(pi*zi)-sin(pi*yi)*pi*cos(pi*xi)
            ub(ix,iy,iz,3)=exp(t)*pi*cos(pi*xi)*sin(pi*zi)-Pi*cos(pi*yi)*sin(pi*zi)
        enddo
    enddo


    end subroutine updateVelBCs

    subroutine nablaCross(f, u, nx, ny, nz, hx, order)
    implicit none
    double precision, dimension(:,:,:,:),allocatable :: f, u
    !double precision f(nx+1,ny+1,nz+1,3), u(nx+1,ny+1,nz+1,3)! 
    double precision hx, xCoeff(3), yCoeff(3), zCoeff(3), error
    double precision df2_dx, df3_dx, df1_dy, df3_dy, df1_dz, df2_dz, ax,bx,ay,by,az,bz,lx,ly,lz,hy,hz,pi,xi,yi,zi
    integer nx,ny,nz, i, j, k, ix, ix2, iy, iy2, iz, iz2, indexX, indexY, indexZ, order
    ! Computes the cross product (nabla X f) = u
    ! assumes a uniform grid; hx=hy=hz, though to change this is really easy.
    ! I tested it with simple functions and with the functions used by liu the max error is 0.25%
    pi=3.14159265358979324D0
    ax = 0
    bx = 1
    ay = 0
    by = 1
    az = 0
    bz = 1

    !*******************************************************************************

    ! Computing the mesh size hx in x-direction
    lx=bx-ax
    hx=lx/nx
    ! Computing the mesh size hy in y-direction
    ly=by-ay
    hy=ly/ny
    ! Computing the mesh size hz in z-direction
    lz=bz-az
    hz=lz/nz

    do j=1,ny+1
        do k=1,nz+1
            do i=1,nx+1
                ! Derivatives with respect of X
                xi=hx*(i-1)/lx
                yi=hy*(j-1)/ly
                zi=hz*(k-1)/lz
                if(i .EQ. 1) then
                    !Forward difference
                    if (order .eq. 1) then
                        df2_dx = (f(2,j,k,2) - f(1,j,k,2)) /hx 
                        df3_dx = (f(2,j,k,3) - f(1,j,k,3)) /hx 
                    endif
                    !second order finite difference
                    if (order .eq. 2) then
                        df2_dx = (-3.0/2.0*f(1,j,k,2) + 2*f(2,j,k,2) - 0.5 * f(3,j,k,2)) /hx 
                        df3_dx = (-3.0/2.0*f(1,j,k,3) + 2*f(2,j,k,3) - 0.5 * f(3,j,k,3)) /hx
                    endif
                    ! third order finite difference
                    !df2_dx = (-11.0/6.0*f(1,j,k,2) + 3*f(2,j,k,2) - 1.5 * f(3,j,k,2) + 1.0/3.0 * f(4,j,k,2)) /hx 
                    !df3_dx = (-11.0/6.0*f(1,j,k,3) + 3*f(2,j,k,3) - 1.5 * f(3,j,k,3) + 1.0/3.0 * f(4,j,k,3)) /hx

                elseif(i .EQ. nx+1) then
                    !Backward difference
                    if (order .eq. 1) then
                        df2_dx = (f(i,j,k,2) - f(i-1,j,k,2)) /hx 
                        df3_dx = (f(i,j,k,3) - f(i-1,j,k,3)) /hx
                    endif
                    !second order finite difference
                    if (order .eq. 2) then
                        df2_dx = (3.0/2.0*f(i,j,k,2) - 2*f(i-1,j,k,2) + 0.5 * f(i-2,j,k,2)) /hx 
                        df3_dx = (3.0/2.0*f(i,j,k,3) - 2*f(i-1,j,k,3) + 0.5 * f(i-2,j,k,3)) /hx
                    endif
                    ! third order finite difference
                    !df2_dx = (11.0/6.0*f(i,j,k,2) - 3.0*f(i-1,j,k,2) + 1.5 * f(i-2,j,k,2) - 1.0/3.0 * f(i-3,j,k,2)) /hx 
                    !df3_dx = (11.0/6.0*f(i,j,k,3) - 3.0*f(i-1,j,k,3) + 1.5 * f(i-2,j,k,3) - 1.0/3.0 * f(i-3,j,k,3)) /hx

                else
                    !Central difference 
                    df2_dx = 0.5* ( f(i+1,j,k,2) - f(i-1,j,k,2))  /hx
                    df3_dx = 0.5* ( f(i+1,j,k,3) - f(i-1,j,k,3))  /hx
                endif

                ! Derivatives with respect of Y
                if(j .EQ. 1) then
                    ! Forward difference

                    if ( k .eq. 140) then

                        print *, i, ' ',j, ' ',k
                        print *, ' analytical ', ' df3_dy ', pi*cos(pi*yi)*sin(pi*xi) !, ' yi ', yi, ' zi ', zi
                        df1_dy = (f(i,2,k,1) - f(i,1,k,1)) /hx 
                        df3_dy = (f(i,2,k,3) - f(i,1,k,3)) /hx 
                        print *, ' 1st order ', ' df3_dy ' ,df3_dy


                        !df1_dy = (-3.0/2.0*f(i,1,k,1) + 2*f(i,2,k,1) - 0.5 * f(i,3,k,1)) /hx 
                        !df3_dy = (-3.0/2.0*f(i,1,k,3) + 2*f(i,2,k,3) - 0.5 * f(i,3,k,3)) /hx 

                        print *, ' 2nd order ', ' df3_dy ' ,df3_dy
                        ! fourth order finite difference

                        !df1_dy = (-25.0/12.0*f(i,1,k,1) + 4*f(i,2,k,1) - 3 * f(i,3,k,1) + 4.0/3.0 * f(i,4,k,1) - 0.25 * f(i,5,k,1)) /hx 
                        !df3_dy = (-25.0/12.0*f(i,1,k,3) + 4*f(i,2,k,3) - 3 * f(i,3,k,3) + 4.0/3.0 * f(i,4,k,3) - 0.25 * f(i,5,k,1)) /hx
                        print *, ' 4th order ', ' df3_dy ' ,df3_dy
                    endif
                    !first order finite diference
                    if (order .eq. 1) then
                        df1_dy = (f(i,2,k,1) - f(i,1,k,1)) /hx 
                        df3_dy = (f(i,2,k,3) - f(i,1,k,3)) /hx 
                    endif
                    !second order finite diference
                    if (order .eq. 2) then
                        df1_dy = (-3.0/2.0*f(i,1,k,1) + 2*f(i,2,k,1) - 0.5 * f(i,3,k,1)) /hx 
                        df3_dy = (-3.0/2.0*f(i,1,k,3) + 2*f(i,2,k,3) - 0.5 * f(i,3,k,3)) /hx 
                    endif
                    !third order finite difference
                    !df1_dy = (-11.0/6.0*f(i,1,k,1) + 3.0*f(i,2,k,1) - 1.5 * f(i,3,k,1) + 1.0/3.0 * f(i,4,k,1)) /hx 
                    !df3_dy = (-11.0/6.0*f(i,1,k,3) + 3.0*f(i,2,k,3) - 1.5 * f(i,3,k,3) + 1.0/3.0 * f(i,4,k,3)) /hx
                    ! fourth order finite difference

                    !df1_dy = (-25.0/12.0*f(i,1,k,1) + 4.0*f(i,2,k,1) - 3.0 * f(i,3,k,1) + 4.0/3.0 * f(i,4,k,1) - 0.25 * f(i,5,k,1)) /hx 
                    !df3_dy = (-25.0/12.0*f(i,1,k,3) + 4.0*f(i,2,k,3) - 3.0 * f(i,3,k,3) + 4.0/3.0 * f(i,4,k,3) - 0.25 * f(i,5,k,1)) /hx
                    if ( k .eq. 140) print *, ' 3rd order ', ' df3_dy ' ,df3_dy

                elseif(j .EQ. ny+1) then                
                    !Backward difference
                    if (order .eq. 1) then
                        df1_dy = (f(i,j,k,1) - f(i,j-1,k,1)) /hx 
                        df3_dy = (f(i,j,k,3) - f(i,j-1,k,3)) /hx 
                    endif
                    !second order finite diference
                    if (order .eq. 2) then
                        df1_dy = (3.0/2.0*f(i,j,k,1) - 2*f(i,j-1,k,1) + 0.5 * f(i,j-2,k,1)) /hx
                        df3_dy = (3.0/2.0*f(i,j,k,3) - 2*f(i,j-1,k,3) + 0.5 * f(i,j-2,k,3)) /hx 
                    endif

                    !third order finite difference
                    !df1_dy = (11.0/6.0*f(i,j,k,1) - 3*f(i,j-1,k,1) + 1.5 * f(i,j-2,k,1) - 1.0/3.0 * f(i,j-3,k,1)) /hx
                    !df3_dy = (11.0/6.0*f(i,j,k,3) - 3*f(i,j-1,k,3) + 1.5 * f(i,j-2,k,3) - 1.0/3.0 * f(i,j-3,k,3)) /hx 
                    ! fourth order finite difference

                    !df1_dy = (25.0/12.0*f(i,j,k,1) - 4.0*f(i,j-1,k,1) + 3.0 * f(i,j-2,k,1) - 4.0/3.0 * f(i,j-3,k,1) - 0.25 * f(i,j-4,k,1) ) /hx
                    !df3_dy = (25.0/12.0*f(i,j,k,3) - 4.0*f(i,j-1,k,3) + 3.0 * f(i,j-2,k,3) - 4.0/3.0 * f(i,j-3,k,3) - 0.25 * f(i,j-4,k,3)) /hx 

                else
                    !Central difference
                    df1_dy = 0.5* ( f(i,j+1,k,1) - f(i,j-1,k,1))  /hx
                    df3_dy = 0.5* ( f(i,j+1,k,3) - f(i,j-1,k,3))  /hx
                endif

                ! Derivatives with respect of Y
                if(k .EQ. 1) then
                    ! Forward difference
                    if (order .eq. 1) then
                        df1_dz = (f(i,j,2,1) - f(i,j,1,1) ) /hx 
                        df2_dz = (f(i,j,2,2) - f(i,j,1,2) ) /hx
                    endif

                    if (order .eq. 2) then
                        df1_dz = (-3.0/2.0*f(i,j,1,1) + 2*f(i,j,2,1) - 0.5 * f(i,j,3,1) ) /hx 
                        df2_dz = (-3.0/2.0*f(i,j,1,2) + 2*f(i,j,2,2) - 0.5 * f(i,j,3,2) ) /hx
                    endif

                    ! third order finite difference
                    !df1_dz = (-11.0/6.0*f(i,j,1,1) + 3*f(i,j,2,1) - 1.5 * f(i,j,3,1) + 1.0/3.0 * f(i,j,4,1)) /hx 
                    !df2_dz = (-11.0/6.0*f(i,j,1,2) + 3*f(i,j,2,2) - 1.5 * f(i,j,3,2) + 1.0/3.0 * f(i,j,4,2)) /hx

                elseif (k .EQ. nz+1) then
                    !Backward difference
                    if (order .eq. 1) then
                        df1_dz = (f(i,j,k,1) - f(i,j,k-1,1) ) /hx 
                        df2_dz = (f(i,j,k,2) - f(i,j,k-1,2) ) /hx
                    end if
                    if (order .eq. 2) then
                        df1_dz = (3.0/2.0*f(i,j,k,1) - 2.0*f(i,j,k-1,1) + 0.5 * f(i,j,k-2,1) ) /hx 
                        df2_dz = (3.0/2.0*f(i,j,k,2) - 2.0*f(i,j,k-1,2) + 0.5 * f(i,j,k-2,2) ) /hx 
                    endif
                    !df1_dz = (11.0/6.0*f(i,j,k,1) - 3*f(i,j,k-1,1) + 1.5 * f(i,j,k-2,1) - 1.0/3.0 * f(i,j,k-3,1)) /hx 
                    !df2_dz = (11.0/6.0*f(i,j,k,2) - 3*f(i,j,k-1,2) + 1.5 * f(i,j,k-2,2) - 1.0/3.0 * f(i,j,k-3,2)) /hx 


                else
                    !central difference
                    df1_dz = 0.5* ( f(i,j,k+1,1) - f(i,j,k-1,1))  /hx
                    df2_dz = 0.5* ( f(i,j,k+1,2) - f(i,j,k-1,2))  /hx
                    if((j .EQ. 1) .and. (k.eq. 140)) then

                        print *, ' analytical ', ' df2_dz ', pi*cos(pi*zi)*sin(pi*xi) !, ' yi ', yi, ' zi ', zi

                        print *, ' 2nd order ', ' df2_dz ' ,df2_dz

                    endif
                endif

                u(i, j, k,1) = df3_dy - df2_dz  ! df3/dy - df2/dz 
                u(i, j, k,2) = df1_dz - df3_dx  ! df1/dz - df3/dx
                u(i, j, k,3) = df2_dx - df1_dy  ! df2/dx - df1/dy

                if((j .EQ. 1) .and. (k.eq. 140) ) then
                    print *, ' result ', u(i, j, k,1)
                    print *, ' error ' , u(i, j, k,1) - (pi*cos(pi*yi)*sin(pi*xi) -pi*cos(pi*zi)*sin(pi*xi))
                endif
            enddo
        enddo
    enddo


    end subroutine nablaCross


    subroutine nablaCrossPer(f, u, nx, ny, nz, hx, order)
    implicit none
    double precision, dimension(:,:,:,:),allocatable :: f, u
    !double precision f(nx+1,ny+1,nz+1,3), u(nx+1,ny+1,nz+1,3)! 
    double precision hx, xCoeff(3), yCoeff(3), zCoeff(3), error
    double precision df2_dx, df3_dx, df1_dy, df3_dy, df1_dz, df2_dz, ax,bx,ay,by,az,bz,lx,ly,lz,hy,hz,pi,xi,yi,zi
    integer nx,ny,nz, i, j, k, ix, ix2, iy, iy2, iz, iz2, indexX, indexY, indexZ, order
    ! Computes the cross product (nabla X f) = u
    ! assumes a uniform grid; hx=hy=hz, though to change this is really easy.
    ! For a peridodic system in X and Z
    pi=3.14159265358979324D0
    ax = 0
    bx = 1
    ay = 0
    by = 1
    az = 0
    bz = 1

    !*******************************************************************************

    ! Computing the mesh size hx in x-direction
    lx=bx-ax
    hx=lx/nx
    ! Computing the mesh size hy in y-direction
    ly=by-ay
    hy=ly/ny
    ! Computing the mesh size hz in z-direction
    lz=bz-az
    hz=lz/nz

    do j=1,ny+1
        do k=1,nz+1
            do i=1,nx+1
                ! Derivatives with respect of X
                xi=hx*(i-1)/lx
                yi=hy*(j-1)/ly
                zi=hz*(k-1)/lz
                if(i .EQ. 1) then
                    !Peridodic
                    df2_dx = 0.5* ( f(2,j,k,2) - f(nx,j,k,2))  /hx
                    df3_dx = 0.5* ( f(2,j,k,3) - f(nx,j,k,3))  /hx             
                    
                elseif(i .EQ. nx+1) then
                    !Peridodic
                    df2_dx = 0.5* ( f(2,j,k,2) - f(nx,j,k,2))  /hx
                    df3_dx = 0.5* ( f(2,j,k,3) - f(nx,j,k,3))  /hx   
                else
                    !Central difference 
                    df2_dx = 0.5* ( f(i+1,j,k,2) - f(i-1,j,k,2))  /hx
                    df3_dx = 0.5* ( f(i+1,j,k,3) - f(i-1,j,k,3))  /hx
                endif

                ! Derivatives with respect of Y
                if(j .EQ. 1) then
                    ! Forward difference
                    !first order finite diference
                    if (order .eq. 1) then
                        df1_dy = (f(i,2,k,1) - f(i,1,k,1)) /hx 
                        df3_dy = (f(i,2,k,3) - f(i,1,k,3)) /hx 
                    endif
                    !second order finite diference
                    if (order .eq. 2) then
                        df1_dy = (-3.0/2.0*f(i,1,k,1) + 2*f(i,2,k,1) - 0.5 * f(i,3,k,1)) /hx 
                        df3_dy = (-3.0/2.0*f(i,1,k,3) + 2*f(i,2,k,3) - 0.5 * f(i,3,k,3)) /hx 
                    endif
                    !third order finite difference
                    !df1_dy = (-11.0/6.0*f(i,1,k,1) + 3.0*f(i,2,k,1) - 1.5 * f(i,3,k,1) + 1.0/3.0 * f(i,4,k,1)) /hx 
                    !df3_dy = (-11.0/6.0*f(i,1,k,3) + 3.0*f(i,2,k,3) - 1.5 * f(i,3,k,3) + 1.0/3.0 * f(i,4,k,3)) /hx
                    ! fourth order finite difference

                    !df1_dy = (-25.0/12.0*f(i,1,k,1) + 4.0*f(i,2,k,1) - 3.0 * f(i,3,k,1) + 4.0/3.0 * f(i,4,k,1) - 0.25 * f(i,5,k,1)) /hx 
                    !df3_dy = (-25.0/12.0*f(i,1,k,3) + 4.0*f(i,2,k,3) - 3.0 * f(i,3,k,3) + 4.0/3.0 * f(i,4,k,3) - 0.25 * f(i,5,k,1)) /hx
                    if ( k .eq. 140) print *, ' 3rd order ', ' df3_dy ' ,df3_dy

                elseif(j .EQ. ny+1) then                
                    !Backward difference
                    if (order .eq. 1) then
                        df1_dy = (f(i,j,k,1) - f(i,j-1,k,1)) /hx 
                        df3_dy = (f(i,j,k,3) - f(i,j-1,k,3)) /hx 
                    endif
                    !second order finite diference
                    if (order .eq. 2) then
                        df1_dy = (3.0/2.0*f(i,j,k,1) - 2*f(i,j-1,k,1) + 0.5 * f(i,j-2,k,1)) /hx
                        df3_dy = (3.0/2.0*f(i,j,k,3) - 2*f(i,j-1,k,3) + 0.5 * f(i,j-2,k,3)) /hx 
                    endif

                    !third order finite difference
                    !df1_dy = (11.0/6.0*f(i,j,k,1) - 3*f(i,j-1,k,1) + 1.5 * f(i,j-2,k,1) - 1.0/3.0 * f(i,j-3,k,1)) /hx
                    !df3_dy = (11.0/6.0*f(i,j,k,3) - 3*f(i,j-1,k,3) + 1.5 * f(i,j-2,k,3) - 1.0/3.0 * f(i,j-3,k,3)) /hx 
                    ! fourth order finite difference

                    !df1_dy = (25.0/12.0*f(i,j,k,1) - 4.0*f(i,j-1,k,1) + 3.0 * f(i,j-2,k,1) - 4.0/3.0 * f(i,j-3,k,1) - 0.25 * f(i,j-4,k,1) ) /hx
                    !df3_dy = (25.0/12.0*f(i,j,k,3) - 4.0*f(i,j-1,k,3) + 3.0 * f(i,j-2,k,3) - 4.0/3.0 * f(i,j-3,k,3) - 0.25 * f(i,j-4,k,3)) /hx 

                else
                    !Central difference
                    df1_dy = 0.5* ( f(i,j+1,k,1) - f(i,j-1,k,1))  /hx
                    df3_dy = 0.5* ( f(i,j+1,k,3) - f(i,j-1,k,3))  /hx
                endif

                ! Derivatives with respect of Z
                if(k .EQ. 1) then
                    !Periodic
                    df1_dz = 0.5* ( f(i,j,2,1) - f(i,j,nz,1))  /hx
                    df2_dz = 0.5* ( f(i,j,2,2) - f(i,j,nz,2))  /hx

                elseif (k .EQ. nz+1) then
                    !Periodic
                    df1_dz = 0.5* ( f(i,j,2,1) - f(i,j,nz,1))  /hx
                    df2_dz = 0.5* ( f(i,j,2,2) - f(i,j,nz,2))  /hx
                else
                    !central difference
                    df1_dz = 0.5* ( f(i,j,k+1,1) - f(i,j,k-1,1))  /hx
                    df2_dz = 0.5* ( f(i,j,k+1,2) - f(i,j,k-1,2))  /hx
                endif

                u(i, j, k,1) = df3_dy - df2_dz  ! df3/dy - df2/dz 
                u(i, j, k,2) = df1_dz - df3_dx  ! df1/dz - df3/dx
                u(i, j, k,3) = df2_dx - df1_dy  ! df2/dx - df1/dy

                if((j .EQ. 1) .and. (k.eq. 140) ) then
                    print *, ' result ', u(i, j, k,1)
                    print *, ' error ' , u(i, j, k,1) - (pi*cos(pi*yi)*sin(pi*xi) -pi*cos(pi*zi)*sin(pi*xi))
                endif
            enddo
        enddo
    enddo


    end subroutine nablaCrossPer
    
    
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

    subroutine computeLaplacianPer(vt, nx, ny, nz, hx, hy, hz, laplacian)
    double precision, dimension(:,:,:,:),allocatable :: vt, laplacian
    double precision hx, hy, hz, dvtdx, dvtdz
    integer nx,ny,nz, i, j, k, ix, ix2, iy, iy2, iz, iz2, indexX, indexY, indexZ
    !Computes the laplacian for a periodic function in the X and Z direction

    ! Need to compute the laplacian just for inner nodes, and eventually if we have periodic conditions
    ! need to evaluate it in the periodic boundaries. 
    do ix=1,nx+1
        do iy=2,ny
            do iz=1,nz+1
                do i=1,3
                    if( ix .eq. 1 )then
                        ! apply periodicity in X
                        dvtdx = -2*vt(ix,iy,iz,i)+vt(2,iy,iz,i)+vt(nx,iy,iz,i)

                    elseif(ix .eq. nx+1) then    
                        ! apply periodicity in X
                        dvtdx = -2*vt(ix,iy,iz,i)+vt(2,iy,iz,i)+vt(nx,iy,iz,i)
                    else
                        ! inner node
                        dvtdx = -2*vt(ix,iy,iz,i)+vt(ix-1,iy,iz,i)+vt(ix+1,iy,iz,i)                       
                    endif

                    if( iz .eq. 1 )then
                        ! apply periodicity in Z
                        dvtdz = -2*vt(ix,iy,iz,i)+vt(ix,iy,2,i)+vt(ix,iy,nz,i)

                    elseif(iz .eq. nz+1) then    
                        ! apply periodicity in Z
                        dvtdz = -2*vt(ix,iy,iz,i)+vt(ix,iy,2,i)+vt(nx,iy,nz,i)
                    else
                        ! inner node
                        dvtdz = -2*vt(ix,iy,iz,i)+vt(ix,iy,iz-1,i)+vt(ix,iy,iz+1,i)                       
                    endif


                    laplacian(ix,iy,iz,i)= ( dvtdx  ) / hx**2 &
                    +( -2*vt(ix,iy,iz,i)+vt(ix,iy-1,iz,i)+vt(ix,iy+1,iz,i)  ) / hx**2 &
                    +( dvtdz  ) / hx**2 

                enddo
            enddo
        enddo
    enddo
    end subroutine computeLaplacianPer


    subroutine computeConvectiveTermDirect(vt, u, conv, nx, ny, nz, hx, hy, hz)
    double precision, dimension(:,:,:,:),allocatable :: vt, u, conv
    double precision hx, hy, hz
    integer nx,ny,nz, i, j, k, ix, ix2, iy, iy2, iz, iz2, indexX, indexY, indexZ
    !I tested it with the function used by Liu the results had less than 0.1% of error

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

    end subroutine computeConvectiveTermDirect

    subroutine computeConvectiveTermSF(f, u, conv, nx, ny, nz, hx, hy, hz)
    double precision, dimension(:,:,:,:),allocatable :: f, u, conv, vt, WxU
    double precision hx, hy, hz
    integer nx,ny,nz, i, j, k, ix, ix2, iy, iy2, iz, iz2, indexX, indexY, indexZ
    allocate(WxU(nx+1,ny+1,nz+1,3), vt(nx+1,ny+1,nz+1,3) )
    WxU=0
    vt=0
    ! Computes the convective term just using the stream functions
    ! as in (nabla X ( nabla X sf) ) X nabla X sf 
    !is eq 3.4 in Lius paper
    ! by now we already have the velocity computed from the stream function, 

    ! we need to compute the vorticity from the velocity with w = nabla x u

    call nablaCross(u, vt, nx, ny, nz, hx,2) ! here we get the vorticity
    !call nablaCross(f, vt, nx, ny, nz, hx) ! here we get the velocity

    ! Now we need to compute the cross product WxU

    do ix=1,nx+1
        do iy=1,ny+1
            do iz=1,nz+1
                WxU(ix,iy,iz,1) = u(ix,iy,iz,3)*vt(ix,iy,iz,2) - u(ix,iy,iz,2)*vt(ix,iy,iz,3)
                WxU(ix,iy,iz,2) = u(ix,iy,iz,1)*vt(ix,iy,iz,3) - u(ix,iy,iz,3)*vt(ix,iy,iz,1)
                WxU(ix,iy,iz,3) = u(ix,iy,iz,2)*vt(ix,iy,iz,1) - u(ix,iy,iz,1)*vt(ix,iy,iz,2)
            enddo
        enddo
    enddo
    ix=1
    iy=15
    iz=16
    print *, 'u numerical ' , u(ix,iy,iz,1), ' ', u(ix,iy,iz,2), ' ', u(ix,iy,iz,3)
    print *, 'vt numerical ' , vt(ix,iy,iz,1), ' ', vt(ix,iy,iz,2), ' ', vt(ix,iy,iz,3)
    print *, 'WxU numerical ' , WxU(ix,iy,iz,1)
    print *, u(ix,iy,iz,3) , ' ', vt(ix,iy,iz,2),' ' , u(ix,iy,iz,2),' ',vt(ix,iy,iz,3)
    ! the convective term is nabla x W x U
    !conv = WxU 

    call nablaCross(WxU, conv , nx, ny, nz, hx,2)
    !conv(1,:,:,:) =0
    deallocate(WxU,vt)
    end subroutine computeConvectiveTermSF

    subroutine computeConvectiveTermSFPer(f, u, conv, nx, ny, nz, hx, hy, hz)
    double precision, dimension(:,:,:,:),allocatable :: f, u, conv, vt, WxU
    double precision hx, hy, hz
    integer nx,ny,nz, i, j, k, ix, ix2, iy, iy2, iz, iz2, indexX, indexY, indexZ
    allocate(WxU(nx+1,ny+1,nz+1,3), vt(nx+1,ny+1,nz+1,3) )
    WxU=0
    vt=0
    ! Computes the convective term just using the stream functions
    ! as in (nabla X ( nabla X sf) ) X nabla X sf 
    !is eq 3.4 in Lius paper
    ! by now we already have the velocity computed from the stream function, 

    ! we need to compute the vorticity from the velocity with w = nabla x u

    call nablaCrossPer(u, vt, nx, ny, nz, hx,2) ! here we get the vorticity
    !call nablaCross(f, vt, nx, ny, nz, hx) ! here we get the velocity

    ! Now we need to compute the cross product WxU

    do ix=1,nx+1
        do iy=1,ny+1
            do iz=1,nz+1
                WxU(ix,iy,iz,1) = u(ix,iy,iz,3)*vt(ix,iy,iz,2) - u(ix,iy,iz,2)*vt(ix,iy,iz,3)
                WxU(ix,iy,iz,2) = u(ix,iy,iz,1)*vt(ix,iy,iz,3) - u(ix,iy,iz,3)*vt(ix,iy,iz,1)
                WxU(ix,iy,iz,3) = u(ix,iy,iz,2)*vt(ix,iy,iz,1) - u(ix,iy,iz,1)*vt(ix,iy,iz,2)
            enddo
        enddo
    enddo
    ix=16
    iy=16
    iz=16
    print *, 'u numerical ' , u(ix,iy,iz,1), ' ', u(ix,iy,iz,2), ' ', u(ix,iy,iz,3)
    print *, 'vt numerical ' , vt(ix,iy,iz,1), ' ', vt(ix,iy,iz,2), ' ', vt(ix,iy,iz,3)
    print *, 'WxU numerical ' , WxU(ix,iy,iz,1)
    print *, u(ix,iy,iz,3) , ' ', vt(ix,iy,iz,2),' ' , u(ix,iy,iz,2),' ',vt(ix,iy,iz,3)
    ! the convective term is nabla x W x U
    !conv = WxU 

    call nablaCrossPer(WxU, conv , nx, ny, nz, hx,2)
    !conv(1,:,:,:) =0
    deallocate(WxU,vt)
    end subroutine computeConvectiveTermSFPer
    
    
    
    subroutine interpolateVelVort(position, velocity, vorticity, u, vt, hx, hy, hz, nx, ny, nz )
    ! takes a position and interpolates the velocity and the vorticity 
    ! assumes that the mesh origin is (0,0,0)
    ! uses a trilinear interpolation
    integer :: nx, ny, nz
    integer:: idx, idy, idz
    double precision :: hx, hy, hz
    double precision, dimension(3) :: position, velocity, vorticity
    double precision, dimension(3) :: xd, c00, c01, c10, c11, c0, c1
    double precision, dimension(:,:,:,:),allocatable:: u, vt
    
    !find index of lower most index of cell containing the point
    idx = FLOOR(position(1)/hx)
    idy = FLOOR(position(2)/hy)
    idz = FLOOR(position(3)/hz)
    
    if (idx<nx+1 .and. idx>-1 .and.  idy<ny+1 .and. idy>-1 .and.  idz<nz+1 .and. idz>-1) then   ! if the point is inside the domain.
        
    xd(1) = (position(1) - (idx-1)*hx )/hx
    xd(2) = (position(2) - (idy-1)*hy )/hy
    xd(3) = (position(3) - (idz-1)*hz )/hz
    
    !For Velocity
     c00 = u(idx,idy,idz,:)    *(1-xd(1))+ u(idx+1,idy,idz,:)    *xd(1)
     c01 = u(idx,idy,idz+1,:)  *(1-xd(1))+ u(idx+1,idy,idz+1,:)  *xd(1)
     c10 = u(idx,idy+1,idz,:)  *(1-xd(1))+ u(idx+1,idy+1,idz,:)  *xd(1)
     c01 = u(idx,idy+1,idz+1,:)*(1-xd(1))+ u(idx+1,idy+1,idz+1,:)*xd(1)
     
     c0 = c00*(1-xd(2)) + c10*xd(2)
     c1 = c01*(1-xd(2)) + c11*xd(2)
     
     velocity = c0*(1-xd(3))+ c1*xd(3)
     
     !For Vorticity
     c00 = vt(idx,idy,idz,:)    *(1-xd(1))+ vt(idx+1,idy,idz,:)    *xd(1)
     c01 = vt(idx,idy,idz+1,:)  *(1-xd(1))+ vt(idx+1,idy,idz+1,:)  *xd(1)
     c10 = vt(idx,idy+1,idz,:)  *(1-xd(1))+ vt(idx+1,idy+1,idz,:)  *xd(1)
     c01 = vt(idx,idy+1,idz+1,:)*(1-xd(1))+ vt(idx+1,idy+1,idz+1,:)*xd(1)
     
     c0 = c00*(1-xd(2)) + c10*xd(2)
     c1 = c01*(1-xd(2)) + c11*xd(2)
     
     vorticity = c0*(1-xd(3))+ c1*xd(3)
    else
        ! the point is outside of the domain
        velocity =0
        vorticity =0
    end if
    
    end subroutine interpolateVelVort
    
    
    subroutine addForce2Fb(position, Fb, force, hx, hy, hz )
    double precision, dimension(:,:,:,:),allocatable:: Fb
    integer:: idx, idy, idz
    double precision :: hx, hy, hz, x0, x1, y0, y1, z0, z1
    double precision, dimension(3) :: position, force, xd
    
    !find index of lower most index of cell containing the point
    idx = FLOOR(position(1)/hx) +1
    idy = FLOOR(position(2)/hy) +1
    idz = FLOOR(position(3)/hz) +1
    
    ! convert to force per unit volume
    force = force/ (hx*hy*hz)
    print '(A11, e10.2,A1,e10.2,A1,e10.2 )', 'force2 ', force(1), ' ', force(2), ' ',force(3)
    x0 = (position(1) - (idx-1)*hx )
    y0 =  (position(2) - (idy-1)*hy )
    z0 = (position(3) - (idz-1)*hz )
    x1 = (position(1) - (idx)*hx )
    y1 = (position(2) - (idy)*hy )
    z1 = (position(3) - (idz)*hz )
    
    print *, 'x0 ', x0/hx, ' y0 ', y0/hx, ' z0 ', z0/hx, 'x1 ', x1/hx, ' y1 ', y1/hx, ' z1 ', z1 /hx
    ! For each point
    ! point 1 0,0,0

    Fb(idx,idy,idz,:) = Fb(idx,idy,idz,:) + force * ( x0 )/hx * ( y0 )/hy * ( z0 )/hz
    print *, force * ( x0 )/hx * ( y0 )/hy * ( z0 )/hz
    
    ! point 2 1,0,0
    Fb(idx+1,idy,idz,:) = Fb(idx+1,idy,idz,:) + force * ( -x1 )/hx * ( y0 )/hy * ( z0 )/hz
    print *, force * ( -x1 )/hx * ( y0 )/hy * ( z0 )/hz
    ! point 3 1,1,0
    Fb(idx+1,idy+1,idz,:) = Fb(idx+1,idy+1,idz,:) + force * ( -x1 )/hx * ( -y1 )/hy * ( z0 )/hz
 print *, force * ( -x1 )/hx * ( -y1 )/hy * ( z0 )/hz
    ! point 4 1,1.1
    Fb(idx+1,idy+1,idz+1,:) = Fb(idx+1,idy+1,idz+1,:) + force * ( -x1 )/hx * ( -y1 )/hy * ( -z1 )/hz
     print *, force * ( -x1 )/hx * ( -y1 )/hy * ( -z1 )/hz
    ! point 5 1,0,1
    Fb(idx+1,idy,idz+1,:) = Fb(idx+1,idy,idz+1,:) + force * ( -x1 )/hx * ( y0 )/hy * ( -z1 )/hz
     print *, force * ( -x1 )/hx * ( y0 )/hy * ( -z1 )/hz
    ! point 6 0,1,1
    Fb(idx,idy+1,idz+1,:) = Fb(idx,idy+1,idz+1,:) + force * ( x0 )/hx * ( -y1 )/hy * (-z1 )/hz
     print *, force * ( x0 )/hx * ( -y1 )/hy * (-z1 )/hz
    ! point  7 0,1,0
    Fb(idx,idy+1,idz,:) = Fb(idx,idy+1,idz,:) + force * ( x0 )/hx * ( -y1 )/hy * ( z0 )/hz
     print *, force * ( x0 )/hx * ( -y1 )/hy * ( z0 )/hz
    ! point 8 0,0,1
    Fb(idx,idy,idz+1,:) = Fb(idx,idy,idz+1,:) + force * ( x0 )/hx * ( y0 )/hy * ( -z1 )/hz
     print *, force * ( x0 )/hx * ( y0 )/hy * ( -z1 )/hz
    
     
     print *, ' '
    print *, Fb(idx,idy,idz,:),Fb(idx+1,idy,idz,:), Fb(idx+1,idy+1,idz,:), Fb(idx+1,idy+1,idz+1,:), Fb(idx+1,idy,idz+1,:),Fb(idx,idy+1,idz+1,:), Fb(idx,idy+1,idz,:), Fb(idx,idy,idz+1,:)
    ! Tested it, the force has the same sign in all 8 pints and their sum is equal to the force
    print *, 'sum ' ,Fb(idx,idy,idz,:)+Fb(idx+1,idy,idz,:)+ Fb(idx+1,idy+1,idz,:)+ Fb(idx+1,idy+1,idz+1,:)+ Fb(idx+1,idy,idz+1,:)+Fb(idx,idy+1,idz+1,:)+ Fb(idx,idy+1,idz,:)+ Fb(idx,idy,idz+1,:)
    
    end subroutine addForce2Fb

    
    
    end module vorticitySolverUtilities