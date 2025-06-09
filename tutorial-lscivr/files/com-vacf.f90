program watervacf
    implicit none
    integer,parameter :: molecule = 216                 ! molecule number
    integer,parameter :: atom = 3                       ! atom number within single molecule
    integer,parameter :: atoms = molecule*atom          ! total atom number
    integer,parameter :: total_step = 40000             ! time steps for nve simu
    integer,parameter :: tau_step = 4000                ! correlation steps
    real(8)           :: velocity(total_step, 3*atoms)  ! velocity of each atom in each dimension
    real(8)           :: v_center(total_step, molecule, 3)
    real(8)           :: vv_acf(0:tau_step,molecule), nacf(0:tau_step,molecule), avg_vv_acf(0:tau_step)
    real(8)           :: mass(atom), sum_mass=0d0       ! atom mass; molecule mass
    integer           :: i, j, k, l                     ! loop variables
    integer           :: iT, iTau, tau, taumax          ! loop varibales used in acf calc
    character         :: headline                       ! first line of mdvel

    ! read velocity from mdvel file ans assign atom masses
    print *, 'start reading velocity'
    open(14,file='mdvel')
    read(14,*) headline
    do i=1,total_step
        read(14,*) velocity(i,:)
    enddo
    close(14)
    print *, 'finish writing velocity'
    mass(1)=16.01
    mass(2)=1.008
    mass(3)=1.008
    do l=1,atom
        sum_mass = sum_mass+mass(l)
    enddo

    ! calculate center-of-mass velocity
    v_center(:,:,:)=0d0
    do i=1,total_step       ! loop over total steps
        do j=1,molecule     ! loop over molecules
            do k=1,atom     ! loop over atoms within each single molecule
                do l=1,3    ! loop over three dimensions
                    v_center(i,j,l)=v_center(i,j,l)+velocity(i,(j-1)*atom*3+(k-1)*3+l)*mass(k)/sum_mass
                enddo
            enddo
        enddo
    enddo
    print *, 'finish COM velocity calculation'

    ! calculate vv acf
    vv_acf(:,:)=0d0
    nacf(:,:)=0d0
    do j=1,molecule
        do iT=1,total_step
            taumax=min(total_step, iT+tau_step)
            do iTau=iT,taumax
                tau=iTau-iT
                vv_acf(tau,j)=vv_acf(tau,j)+dot_product(v_center(iT,j,:),v_center(iTau,j,:))
                nacf(tau,j)=nacf(tau,j)+1
            enddo
        enddo
    enddo
    ! divide acf by the sum number
    do iTau=0,tau_step
        do i=1,molecule
            vv_acf(iTau,i)=vv_acf(iTau,i)/nacf(iTau,i)
        enddo
    enddo
    print *, 'finish vv_acf calculation'

    ! write output files
    do tau=0,tau_step
        avg_vv_acf=sum(vv_acf, dim=2)/molecule  ! average acf over different molecules
    enddo
    open(15,file='vv-correlation-mine.dat')
    do iTau=0,tau_step
        write(15,*) avg_vv_acf(iTau)
    enddo
    close(15)
    print *, 'finish writing output file'

end program watervacf