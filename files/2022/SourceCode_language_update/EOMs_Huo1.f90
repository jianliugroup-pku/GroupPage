!!==============================================================================
!! Implement of the first version of equations of motion(EOMs) of paper we
!! comment on.
!!
!! Note: in this version, the specific definition of coherent state used by Huo 
!! et. al. totally mismatched with EOMs at all! They did not give the correct  
!! implementation of EOMs on Stratonovich-Weyl phase space at this time, as the 
!! population dynamics is far from exact results in following tests.
!!
!!==============================================================================

!!==============================================================================
!! eigen value problem. only used for test (for exact solution of TDSE).
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
!! symmetrization for Hamiltonian matrix
!! Asym-type (alpha labeled): saved in L-block
!! Bsym-type (beta labeled): saved in U-block
!! Hmat generally can be complex-valued. But in this file, real-value is enough
subroutine calc_Hsym(Hsym, Hmat, F) ! EQ (11)
implicit none
    integer, intent(in) :: F 
    real(kind=8), intent(in) :: Hmat(F,F)
    real(kind=8), intent(inout) :: Hsym(F,F)
    integer :: n,m,l 

    do m=1,F-1
    do n=m+1,F ! n>m, so (n,m) saves Asym-type, (m,n) saves Bsym-type
    Hsym(n,m) = real(Hmat(m,n) + Hmat(n,m))
    Hsym(m,n) = real((0.0d0, 1.0d0)*(Hmat(m,n) - Hmat(n,m)))
    enddo 
    enddo

    do n=2,F
    Hsym(n,n) = -sqrt(2.0d0*(n-1)/n) * real(Hmat(n,n))
    do l=1,n-1 
        Hsym(n,n) = Hsym(n,n) + sqrt(2.0d0/(n*(n-1))) * real(Hmat(l,l))
    enddo
    enddo
end subroutine

!!==============================================================================
!! derivatives of mapping Hamiltonian with varphi angles
subroutine calc_dHdvarphi(dHdvarphi, Hsym, Omega, F) ! EQ (101)
implicit none
    integer, intent(in) :: F 
    real(kind=8), intent(in) :: Omega(F,F), Hsym(F,F)
    real(kind=8), intent(out) :: dHdvarphi(F-1)
    integer :: n,j,k
    do n=1,F-1
    dHdvarphi(n) = 0
    do j=n+1,F
    do k=1,n ! k <= n < j, so (j,k) saves Asym-type, (k,j) saves B
        dHdvarphi(n) = dHdvarphi(n)-(Hsym(j,k)*Omega(k,j)-Hsym(k,j)*Omega(j,k))
    enddo
    enddo
    enddo
end subroutine

!!==============================================================================
!! definition of coherent state used by Huo et. al., which mismatches the EOMs.
subroutine calc_cs(c, theta, varphi, F) ! EQ (17)
implicit none
    integer, intent(in) :: F
    real(kind=8), intent(in) :: theta(F-1), varphi(F-1)
    complex(kind=8), intent(inout) :: c(F)
    integer :: n,l
    c(1) = cos(theta(1)/2) * exp(-(0.0d0,1.0d0)*varphi(1) /2)
    do n=2, F-1
        c(n) = cos(theta(n)/2) * exp(-(0.0d0,1.0d0)*varphi(n)/2)
        do l=1,n-1
            c(n) = c(n) * sin(theta(l)/2) * exp((0.0d0,1.0d0)*varphi(l)/2)
        enddo
    enddo
    c(F) = (1.0d0,0.0d0)
    do l=1,F-1
        c(F) = c(F) * sin(theta(l)/2) * exp((0.0d0,1.0d0)*varphi(l)/2)
    enddo
end subroutine

!!==============================================================================
!! conjugated variables with varphi
subroutine calc_Omega(Omega, theta, varphi, F) ! Eq (B2-B4)
implicit none
    integer, intent(in) :: F
    real(kind=8), intent(in) :: theta(F-1), varphi(F-1)
    real(kind=8), intent(out) :: Omega(F, F)

    real(kind=8) :: sumvarphis, u
    integer :: j,k,l,n,m
    
    !! OmegaAsym: saved left > right
    !! OmegaBsym: saved left < right
    !! OmegaCsym: saved left = right !! from 2 to count!
    do m=1,F-1
    do n=m+1,F !! always n > m, so (n,m) save Asym, (m,n) save Bsym

    Omega(n,m) = 1.0d0
    !if(m > 1) then
    do j=1, m-1
        Omega(n,m) = Omega(n,m) * sin(theta(j)/2) ** 2
    enddo
    !endif

    Omega(n,m) = Omega(n,m) * cos(theta(m)/2) 
    do k=m, n-1
        Omega(n,m) = Omega(n,m) * sin(theta(k)/2)
    enddo
    if(n .ne. F ) then
        Omega(n,m) = Omega(n,m) * cos(theta(n)/2)
    endif

    sumvarphis = 0.0d0
    do l=m, n-1
        sumvarphis = sumvarphis + varphi(l)
    enddo
    Omega(m,n) = Omega(n,m) * sin(sumvarphis) !! OmegaBsym
    Omega(n,m) = Omega(n,m) * cos(sumvarphis) !! OmegaAsym

    enddo
    enddo

    !! OmegaCsym
    do n=2, F
        Omega(n,n) = 0.0
        do j=1,n-1
            u = cos(theta(j)/2)**2
            do k=1,j-1
                u = u * sin(theta(k)/2)**2
            enddo
            Omega(n,n) = Omega(n,n) + u
        enddo

        u = (1-n)
        if(n .ne. F) then 
            u = u * cos(theta(n)/2)**2
        endif
        do j=1,n-1
            u = u * sin(theta(j)/2)**2
        enddo
        Omega(n,n) = Omega(n,n) + u
        Omega(n,n) = Omega(n,n) / sqrt(2.0d0 * n * (n-1))
    enddo
end subroutine

!!==============================================================================
!! derivatives of conjugated variables Omega with time
!! (only necessary values are calculated)
subroutine calc_dOmegadt(dOmegadt, Hsym, Omega, F) ! Eq (D5-D6)
implicit none
    integer, intent(in) :: F
    real(kind=8), intent(inout) :: Hsym(F,F), Omega(F,F), dOmegadt(F,F)
    integer :: n,k,j,l

    do n=1,F
        k=n+1  ! n<k, so (k,n) Asym, (n,k) Bsym. And k for Csym
        !! Asym
        dOmegadt(k,n) = sqrt((n-1.0d0)/(2*n)) * &
            (Hsym(n,n) * Omega(n,k) - Hsym(n,k)*Omega(n,n)) &
            - sqrt((n+1.0d0)/(2*n)) * &
            (Hsym(k,k) * Omega(n,k) - Hsym(n,k)*Omega(k,k))
        !! Bsym
        dOmegadt(n,k) = sqrt((n+1.0d0)/(2*n)) * &
            (Hsym(k,k) * Omega(k,n) - Hsym(k,n)*Omega(k,k)) &
            - sqrt((n-1.0d0)/(2*n)) * &
            (Hsym(n,n) * Omega(k,n) - Hsym(k,n)*Omega(n,n))

        do j=1,n-1 ! j < n < k, so (n,j) Asym, (j,n) Bsym
        dOmegadt(k,n) = dOmegadt(k,n) + 0.5d0 * (&
            Hsym(j,n)*Omega(k,j) - Hsym(k,j)*Omega(j,n) &
          - Hsym(n,j)*Omega(j,k) + Hsym(j,k)*Omega(n,j) &
        ) 
        dOmegadt(n,k) = dOmegadt(n,k) + 0.5d0 * (&
            Hsym(n,j)*Omega(k,j) - Hsym(k,j)*Omega(n,j) &
          + Hsym(j,n)*Omega(j,k) - Hsym(j,k)*Omega(j,n) &
        ) 
        enddo
        do l=n+2,F ! only case k=n+1, so n<k<l
        dOmegadt(k,n) = dOmegadt(k,n) + 0.5d0 * (&
            Hsym(l,n)*Omega(k,l) - Hsym(k,l)*Omega(l,n) &
          - Hsym(n,l)*Omega(l,k) + Hsym(l,k)*Omega(n,l) &
        ) 
        dOmegadt(n,k) = dOmegadt(n,k) + 0.5d0 * (&
            Hsym(l,n)*Omega(l,k) - Hsym(l,k)*Omega(l,n) &
          + Hsym(n,l)*Omega(k,l) - Hsym(k,l)*Omega(n,l) &
        ) 
        enddo
    enddo
end subroutine

!!==============================================================================
!! derivatives of angle variables with time ! Eq (105)
subroutine calc_dangledt(dthetadt, dvarphidt, theta, varphi, dHdvarphi,&
 Omega, dOmegadt, F)
implicit none 
    integer, intent(in) :: F 
    real(kind=8), intent(in) :: theta(F-1), varphi(F-1), dOmegadt(F, F), &
                        Omega(F, F), dHdvarphi(F-1)
    real(kind=8), intent(inout) :: dthetadt(F-1), dvarphidt(F-1)
    integer :: n,j
    real(kind=8) :: prods2

    do n=1,F-1
        !! dthetadt
        prods2 = 1.0d0 
        do j=1,n-1
            prods2 = prods2 * sin(theta(j)/2)**2
        enddo
        if (n .eq. 1) then
            dthetadt(n) = dHdvarphi(n) * 2.0d0 / sin(theta(n)) / prods2
        else 
            dthetadt(n) = ( dHdvarphi(n) * 2.0d0 / sin(theta(n)) &
                          - dHdvarphi(n-1) * tan(theta(n)/2) ) / prods2
        endif
        !! dvarphidt
        dvarphidt(n) = ( dOmegadt(n,n+1)*Omega(n+1,n) - &
            Omega(n,n+1)*dOmegadt(n+1,n) ) / (Omega(n+1,n)**2 + Omega(n,n+1)**2)
    enddo
end subroutine

program main
implicit none
    integer :: F
    real(kind=8), allocatable :: Hsym(:,:), Hmat(:,:), theta(:), varphi(:), &
                        dHdvarphi(:), Omega(:,:), dOmegadt(:,:), &
                        dthetadt(:), dvarphidt(:), &
                        Vs(:,:), Es(:)
    complex(kind=8), allocatable :: Uprop(:,:), expmiHdt(:,:),&
                        c(:), c_corr(:), c_exact(:)
    real(kind=8) :: dt, global_phase = 0.0d0
    real(kind=8), parameter :: pi = 3.1415926535897932384626d0

    integer :: i,j,icount,Ntime,Ns
    character(len=32) :: inputf, outputf

    icount = command_argument_count()
    if (icount.ne.1) then
        stop "need one input file"
    endif
    call get_command_argument(1, inputf)
    outputf = trim(inputf) // "-H1.out"

    open(20, file=inputf)
    read(20,*)
    read(20,*) F, Ntime, Ns, dt

    ! allocate array
    allocate(Hmat(F,F),Hsym(F,F),theta(F-1),varphi(F-1),&
        dthetadt(F-1),dvarphidt(F-1))
    allocate(dHdvarphi(F-1),Omega(F,F), dOmegadt(F,F), Vs(F,F), Es(F))
    allocate(Uprop(F,F),c(F), c_corr(F), c_exact(F), expmiHdt(F,F))

    read(20,*)
    read(20,*) theta(:)
    read(20,*) varphi(:)
    ! read(20,*) c_exact(:) ! read c
    read(20,*)
    read(20,*) Hmat(:,:) ! read Hamiltonian
    close(20)

    call calc_Hsym(Hsym, Hmat, F) ! EQ (11)
    call calc_cs(c, theta, varphi, F)
    c_exact(:) = c(:)

    print*, "initial theta:", theta
    print*, "initial varphi:", varphi
    print*, "initial c:", c_exact

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
            write(29,*) i*dt, theta, varphi, global_phase, &
                real(c), aimag(c), abs(c)**2, &
                real(c_corr), aimag(c_corr), abs(c_corr)**2, &
                real(c_exact), aimag(c_exact), abs(c_exact)**2
        endif

        !! this part evolved by Meyer-Miller EOMs in VV (real part x, imag part p)
        c_exact = c_exact - (0.0d0,1.0d0)*matmul(Hmat, real(c_exact)) *0.5d0*dt
        c_exact = c_exact + matmul(Hmat, aimag(c_exact)) * dt
        c_exact = c_exact - (0.0d0,1.0d0)*matmul(Hmat, real(c_exact)) *0.5d0*dt
        !! by exact solution (only for debug)
        ! c_exact = matmul(Uprop, c_exact)

        !!!<<<<<<<< examples using SW EOMs in this version
        !! first step in VV
        !call calc_cs(c, theta, varphi, F) !< update it if you have not
        call calc_Omega(Omega, theta, varphi, F)
        call calc_dOmegadt(dOmegadt, Hsym, Omega, F)
        call calc_dHdvarphi(dHdvarphi, Hsym, Omega, F)
        call calc_dangledt(dthetadt,dvarphidt, theta,varphi,dHdvarphi,&
            Omega,dOmegadt,F)
        theta = theta + dthetadt * 0.5d0*dt

        !! middle step in VV
        call calc_Omega(Omega, theta, varphi, F)
        call calc_dOmegadt(dOmegadt, Hsym, Omega, F)
        call calc_dHdvarphi(dHdvarphi, Hsym, Omega, F)
        call calc_dangledt(dthetadt,dvarphidt, theta,varphi,dHdvarphi,&
            Omega,dOmegadt,F)
        varphi = varphi + dvarphidt * dt

        !! last step in VV
        call calc_Omega(Omega, theta, varphi, F)
        call calc_dOmegadt(dOmegadt, Hsym, Omega, F)
        call calc_dHdvarphi(dHdvarphi, Hsym, Omega, F)
        call calc_dangledt(dthetadt,dvarphidt, theta,varphi,dHdvarphi,&
            Omega,dOmegadt,F)
        theta = theta + dthetadt * 0.5d0*dt
        call calc_cs(c, theta, varphi, F)  ! EQ (17) & (B2)
        !!!>>>>>>>>
    enddo
    close(29)
    stop
    deallocate(Hmat, Hsym, theta, varphi, dthetadt, dHdvarphi)
    deallocate(dHdvarphi, Omega, dOmegadt, Vs, Es)
    deallocate(Uprop, c, c_corr, c_exact, expmiHdt)
end program
