!!==============================================================================
!! This file is an implementation of eq(S53) in: 10.1002/wcms.1619
!! @author He X <xshinhe@pku.edu.cn>.
!!
!! to cite: He X, Wu B, Shang Y, Li B, Cheng X, Liu J. New phase space
!!          formulations and quantum dynamics approaches. WIREs Comput.
!!          Mol. Sci. 2022. e1619. https://doi.org/10.1002/wcms.1619
!!
!! We would like to share following facts:
!!
!! 1) Dynamics in Stratonovich-Weyl(SW) phase space inevitably face the
!!      singularity problems,especially when population is focused on a pure state.
!!      However, this is the very commonly used initial condition.
!!      In simulations, one has to choose initial condition little deviated
!!      from the singular point to make use of dynamics in SW phase space.
!!
!! 2) 2*F-2 dynamics variables are not enough to describe all properties
!!      whose dynamics firstly proposed in eq(S48) of 10.1002/wcms.1619.
!!      It means SW phase space is intrinsically not a proper representation for 
!!      properties with global phase. We will detail it in another EOMs_Liu_CI.f90
!!      file).
!!
!! 3) We firstly construct correct one-to-one mapping from SW to Meyer-Miller
!!      variables, by introducing an EVOLVING global phase factor, in eq(S53) of
!!      10.1002/wcms.1619. It utilizes 2*F-1 variables other than 2*F-2.
!!      We again emphasize that global phase should not be ignored.
!!
!! We provide basic interface for everyone who is interested in eq(S53)
!! in 10.1002/wcms.1619 (it contains global phase), though we recommend
!! people to use the linear Meyer-Miller EOMs other than 2*F-1 angular variables.
!! The former is easy to use and free of singularity problems.
!!
!!==============================================================================

!!==============================================================================
!! matrix inverse. [ only used for double check of analytical
!! symplectic structure. ]
!! :    Ainv = A ^ -1
subroutine calc_inv(Ainv, A, M)
implicit none
    integer, intent(in) :: M
    real(kind=8),intent(in) :: A(M,M)
    real(kind=8),intent(inout):: Ainv(M,M)
    real(kind=8) :: work(size(A,1))
    integer :: n,info,ipiv(size(A,1))
    Ainv = A
    n = size(A,1)
    call DGETRF(n,n,Ainv,n,ipiv,info) !< LU factorization
    if (info.ne.0) stop 'Matrix is numerically singular!'
    call DGETRI(n,Ainv,n,ipiv,work,n,info) !< inverse of a matrix using the LU
    if (info.ne.0) stop 'Matrix inversion failed!'
end subroutine

!!==============================================================================
!! Eigenproblem; only used for test (for exact solution of TDSE).
subroutine lin_ev(Es, Vs, A, M)
    real(8), intent(in) :: A(M,M)
    real(8), intent(inout) :: Vs(M,M)
    real(8), intent(inout) :: Es(M)
    integer :: LDA, n, lwork, info
    real(8) :: tmp_work(1)
    real(8), allocatable :: work(:)
    ! External procedures defined in LAPACK
    external DSYEV
    ! Store A
    Vs = A

    n = size(A,1)
    LDA = n

    ! Query the optimal workspace.
    lwork = -1
    call DSYEV( 'V', 'U', n, Vs, LDA, Es, tmp_work, lwork, info )
    lwork = max( int(tmp_work(1)), 3*n-1 )
    allocate( work(lwork ) )

    ! Solve eigenproblem.
    CALL DSYEV( 'V', 'U', n, Vs, LDA, Es, work, lwork, info )

    deallocate( work )
    ! Check for convergence.
    if( info /= 0 ) then
        stop 'The algorithm failed to compute eigenvalues.'
    endif
end subroutine

!!==============================================================================
!! Definition of coherent state.
!! Liu's work uses Tilma's definition, which gives a very sparse symplectic
!! matrix.
subroutine calc_cs(c, theta, varphi, F)
implicit none
    integer, intent(in) :: F
    real(kind=8), intent(in) :: theta(F-1), varphi(F-1)
    complex(kind=8), intent(inout) :: c(F)
    integer :: n,l
    
    ! c(1) special case
    c(1) = (1.0d0,0.0d0)
    do l=1,F-2
        c(1) = c(1) * cos(theta(l)) * exp((0.0d0,1.0d0)*varphi(l))
    enddo
    c(1) = c(1) * sin(theta(F-1)) * exp((0.0d0,1.0d0)*varphi(F-1))

    ! c(F) special case
    c(F) = cos(theta(F-1))

    if (F==2) return
    ! c(2) special case
    c(2) = - sin(theta(1)) * exp(-(0.0d0,1.0d0)*varphi(1))
    do l=2,F-2
        c(2) = c(2) * cos(theta(l)) * exp((0.0d0,1.0d0)*varphi(l))
    enddo
    c(2) = c(2) * sin(theta(F-1)) * exp((0.0d0,1.0d0)*varphi(F-1))

    ! other cases
    do n=3, F-1
        c(n) = - sin(theta(n-1))
        do l=n,F-2
            c(n) = c(n) * cos(theta(l)) * exp((0.0d0,1.0d0)*varphi(l))
        enddo
        c(n) = c(n) * sin(theta(F-1)) * exp((0.0d0,1.0d0)*varphi(F-1))
    enddo
end subroutine

!!==============================================================================
!! Derivatives of coherent state with parameteried angles (theta, varphi).
!! dcdtheta = dc / dtheta, dcdvarphi = dc / dvarphi
subroutine calc_dcsdangles(dcdtheta, dcdvarphi, c, theta, varphi, F)
implicit none 
    integer, intent(in) :: F
    real(kind=8), intent(in) :: theta(F-1), varphi(F-1)
    complex(kind=8), intent(in) :: c(F)
    complex(kind=8), intent(inout) :: dcdtheta(F, F-1), dcdvarphi(F, F-1)
    integer :: n,l
    
    dcdtheta = (0.0d0,0.0d0)
    dcdvarphi = (0.0d0,0.0d0)

    ! dc(1,:) special case
    dcdtheta(1,1) = c(1) * (-tan(theta(1)))
    dcdvarphi(1,1) = c(1) * (0.0d0, 1.0d0)
    do l=1,F-2
        dcdtheta(1,l) = c(1) * (-tan(theta(l)))
        dcdvarphi(1,l) = c(1) * (0.0d0, 1.0d0)
    enddo
    dcdtheta(1,F-1) = c(1) / (tan(theta(F-1)))
    dcdvarphi(1,F-1) = c(1) * (0.0d0, 1.0d0)

    ! dc(F,:) special case
    dcdtheta(F,F-1) = - sin(theta(F-1))

    if (F==2) return
    ! dc(2,:) special case
    dcdtheta(2,1) = c(2) / (tan(theta(1)))
    dcdvarphi(2,1) = - c(2) * (0.0d0, 1.0d0)
    do l=2,F-2
        dcdtheta(2,l) = c(2) * (-tan(theta(l)))
        dcdvarphi(2,l) =  c(2) * (0.0d0, 1.0d0)
    enddo
    dcdtheta(2,F-1) = c(2) / (tan(theta(F-1)))
    dcdvarphi(2,F-1) =  c(2) * (0.0d0, 1.0d0)

    ! other cases
    do n=3, F-1
        dcdtheta(n,n-1) = c(n) / (tan(theta(n-1)))
        do l=n,F-2
            dcdtheta(n,l) = c(n) * (-tan(theta(l)))
            dcdvarphi(n,l) =  c(n) * (0.0d0, 1.0d0)
        enddo
        dcdtheta(n,F-1) = c(n) / (tan(theta(F-1)))
        dcdvarphi(n,F-1) =  c(n) * (0.0d0, 1.0d0)
    enddo

    ! print *, dcdtheta(F,:)
    ! stop
end subroutine

!!==============================================================================
!! Symplectic structure for Tilma's coherent state.
!! It is an implementation of eq(S49) in 10.1002/wcms.1619.
subroutine calc_symplectic_structure_analytical(sigma, theta, varphi, F)
implicit none
    integer, intent(in) :: F
    real(kind=8), intent(in) :: theta(F-1), varphi(F-1)
    real(kind=8), intent(inout) :: sigma(2*F-2, 2*F-2)
    integer :: i,j,l

    sigma = 0.0d0
    ! j = 1
    sigma(F,1) = 0.5d0 / sin(2*theta(1)) / sin(theta(F-1))**2
    sigma(F+1,1) = -0.5d0 / tan(2*theta(1)) / sin(theta(F-1))**2
    do l=2,F-2
        sigma(F,1) = sigma(F,1) / cos(theta(l))**2
        sigma(F+1,1) = sigma(F+1,1) / cos(theta(l))**2
    enddo
    ! j = F-1
    sigma(2*F-2,F-1) = -1.0d0/sin(2*theta(F-1))
    ! others
    do j=2,F-2
        ! i=j
        i = j
        sigma(F-1+i,j) = 1.0d0 / sin(2*theta(i)) / sin(theta(F-1))**2
        do l=i+1,F-2
            sigma(F-1+i,j) = sigma(F-1+i,j) / cos(theta(l))**2
        enddo
        ! i=j+1
        i = j+1
        sigma(F-1+i,j) = -0.5d0 / tan(theta(i)) / sin(theta(F-1))**2
        do l=i+1,F-2
            sigma(F-1+i,j) = sigma(F-1+i,j) / cos(theta(l))**2
        enddo
    enddo
    sigma(1:F-1,F:2*F-1) = - transpose(sigma(F:2*F-1,1:F-1))
end subroutine

!!==============================================================================
!! Symplectic structure can be numerical calculated from definition
!! [see more on Provost' method in Commun. Math. Phys. 76, 289(1980)].
!! Do not use this subroutine. Only for test.
!! We have already given the analytical formula above.
subroutine calc_symplectic_structure_numerical(sigma, dcdtheta, dcdvarphi, F)
implicit none
    integer, intent(in) :: F 
    complex(kind=8), intent(in) :: dcdtheta(F,F-1), dcdvarphi(F,F-1)
    real(kind=8), intent(inout) :: sigma(2*F-2, 2*F-2)
    real(kind=8) :: sigmainv(2*F-2, 2*F-2)

    sigmainv = 0.0d0
    sigmainv(F:2*F-2,1:F-1) = -2*aimag( &
                                matmul(transpose(conjg(dcdvarphi)),dcdtheta))
    sigmainv(1:F-1,F:2*F-2) = -2*aimag( &
                                matmul(transpose(conjg(dcdtheta)), dcdvarphi))
    call calc_inv(sigma, sigmainv, 2*F-2)
end subroutine

!!==============================================================================
!! Derivatives of mapping Hamiltonian with angles. It is easily obtained with
!! derivatives of coherent states.
subroutine calc_dHdangles(dHdtheta, dHdvarphi, Hmat, c, dcdtheta, dcdvarphi, F)
implicit none
    integer, intent(in) :: F 
    real(kind=8), intent(in) :: Hmat(F,F)
    complex(kind=8), intent(in) :: c(F), dcdtheta(F,F-1), dcdvarphi(F,F-1)
    real(kind=8), intent(out) :: dHdtheta(F-1), dHdvarphi(F-1)

    ! let variables lambda=1, that it's a conserved variable
    dHdtheta = real( matmul(transpose(conjg(dcdtheta)), matmul(Hmat, c)) &
                + matmul(conjg(c), matmul(Hmat, dcdtheta)) )
    dHdvarphi = real( matmul(transpose(conjg(dcdvarphi)), matmul(Hmat, c)) &
                + matmul(conjg(c), matmul(Hmat, dcdvarphi)) )
end subroutine

program main
implicit none
    integer :: F

    real(kind=8), allocatable :: Hmat(:,:), theta(:), varphi(:), &
                                dHdtheta(:), dHdvarphi(:), sigma(:,:), &
                                Vs(:,:), Es(:)
    complex(kind=8), allocatable :: Uprop(:,:), c(:), c_corr(:), c_exact(:), &
                                dcdtheta(:,:), dcdvarphi(:,:), expmiHdt(:,:)
    real(kind=8) :: dt, global_phase = 0.0d0, gp0
    real(kind=8), parameter :: pi = 3.1415926535897932384626d0
    integer :: i,j,icount,Ntime,Ns
    character(len=32) :: inputf, outputf

    icount = command_argument_count()
    if (icount.ne.1) then
        stop "need one input file"
    endif
    call get_command_argument(1, inputf)
    outputf = trim(inputf) // "-L.out"

    open(20, file=inputf)
    read(20,*)
    read(20,*) F, Ntime, Ns, dt

    ! allocate array
    allocate(Hmat(F,F), theta(F-1), varphi(F-1), dHdtheta(F-1), dHdvarphi(F-1))
    allocate(sigma(2*F-2,2*F-2), Vs(F,F), Es(F))
    allocate(Uprop(F,F),c(F), c_corr(F), c_exact(F))
    allocate(dcdtheta(F,F-1),dcdvarphi(F,F-1),expmiHdt(F,F))

    read(20,*)
    ! read(20,*) theta(:) ! read initial angles
    ! read(20,*) varphi(:)
    read(20,*) c_exact(:) ! read c alternatively
    read(20,*)
    read(20,*) Hmat(:,:) ! read Hamiltonian
    close(20)

    c_exact = c_exact / sqrt(sum(abs(c_exact)**2))
    c(:) = c_exact(:)

    global_phase = atan2(aimag(c(3)), real(c(3)))
    gp0 = global_phase ! save "global phase" at zero time

    c(:) = c(:) / exp((0.0d0,1.0d0)*gp0) ! extract gp0 (and add in the end)

    !! conversion c(:) to theta & varphi
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    if(F.ne.3) stop "this conversion is specific in F=3 case"
    theta(2) = acos(real(c(3)))
    theta(1) = asin(abs(c(2)) / sin(theta(2)))
    if(theta(1) > pi) theta(1) = 2*pi - theta(1)

    c(1) = c(1) / cos(theta(1)) / sin(theta(2))
    c(2) = -c(2) / sin(theta(1)) / sin(theta(2))

    varphi(1) = 0.5d0 * (atan2(aimag(c(1)), real(c(1))) &
        - atan2(aimag(c(2)), real(c(2))) )
    varphi(2) = 0.5d0 * (atan2(aimag(c(1)), real(c(1))) &
        + atan2(aimag(c(2)), real(c(2))) )

    call calc_cs(c, theta, varphi, F)
    c_exact(:) = c(:) * exp((0.0d0,1.0d0)*gp0)

    !! conversion report
    print*, "initial theta:", theta
    print*, "initial varphi:", varphi
    print*, "initial c:", c_exact
    ! stop
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    call lin_ev(Es, Vs, Hmat, F)
    expmiHdt = 0
    do i=1,F
        expmiHdt(i,i) = exp(-(0.0d0,1.0d0)*Es(i)*dt)
    enddo
    Uprop = matmul(Vs, matmul(expmiHdt, transpose(Vs)))

    open(29, file=outputf)
    do i=0,Ntime
        if(mod(i,Ns) == 0) then
            c_corr = c
            c_corr = c_corr * exp((0.0d0,1.0d0)*global_phase)
            write(29,*) i*dt, theta, varphi, global_phase,   &
                real(c*exp((0.0d0,1.0d0)*gp0)),              &
                aimag(c*exp((0.0d0,1.0d0)*gp0)), abs(c)**2,  &
                real(c_corr), aimag(c_corr), abs(c_corr)**2, &
                real(c_exact), aimag(c_exact), abs(c_exact)**2
        endif

        !! this part evolved by Meyer-Miller EOMs in VV (real part x, imag part p)
        c_exact = c_exact - (0.0d0,1.0d0)*matmul(Hmat, real(c_exact)) *0.5d0*dt
        c_exact = c_exact + matmul(Hmat, aimag(c_exact)) * dt
        c_exact = c_exact - (0.0d0,1.0d0)*matmul(Hmat, real(c_exact)) *0.5d0*dt
        !! by exact solution
        ! c_exact = matmul(Uprop, c_exact)

        !!!<<<<<<<< examples using SW EOMs corrected with global phase
        !! first step in VV
        !call calc_cs(c, theta, varphi, F) !< update it if you have not
        call calc_dcsdangles(dcdtheta, dcdvarphi, c, theta, varphi, F)
        call calc_symplectic_structure_analytical(sigma, theta, varphi, F)
        call calc_dHdangles(dHdtheta,dHdvarphi,Hmat,c,dcdtheta,dcdvarphi,F)
        
        theta = theta + matmul(sigma(1:F-1,F:2*F-2), dHdvarphi) * 0.5d0*dt
        call calc_cs(c, theta, varphi, F) !< update current coherent state

        !! middle step in VV
        call calc_dcsdangles(dcdtheta, dcdvarphi, c, theta, varphi, F)
        call calc_symplectic_structure_analytical(sigma, theta, varphi, F)
        call calc_dHdangles(dHdtheta,dHdvarphi,Hmat,c,dcdtheta,dcdvarphi,F)

        varphi = varphi + matmul(sigma(F:2*F-2,1:F-1), dHdtheta) * dt
        global_phase = global_phase - dot_product(c, matmul(Hmat, c)) * dt &
                        + 0.5d0 * tan(theta(F-1)) * dHdtheta(F-1) * dt
        call calc_cs(c, theta, varphi, F) !< update current coherent state

        !! last step in VV
        call calc_dcsdangles(dcdtheta, dcdvarphi, c, theta, varphi, F)
        call calc_symplectic_structure_analytical(sigma, theta, varphi, F)
        call calc_dHdangles(dHdtheta,dHdvarphi,Hmat,c,dcdtheta,dcdvarphi,F)

        theta = theta + matmul(sigma(1:F-1,F:2*F-2), dHdvarphi) * 0.5d0*dt
        call calc_cs(c, theta, varphi, F) !< update current coherent state
        !!!>>>>>>>>
    enddo
    close(29)
    stop
    deallocate(Hmat, theta, varphi, dHdtheta, dHdvarphi, sigma, Vs, Es)
    deallocate(Uprop, c, c_corr, c_exact, dcdtheta, dcdvarphi, expmiHdt)
end program