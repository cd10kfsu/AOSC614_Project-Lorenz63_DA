program test_inc4dvar
!$$$  program documentation block
!         .           .            .
!  program name: test_3dvar
!    programmer: da,cheng        org: umd      date: 2015-Jan-11
!
!  purpose:
!
!  revision history:
!    2015-Jan-11     da    - creator
!
!  file dependencies:
!
!  attributes: 
!    language: fortran 90
!    machine : 
!
!
!$$$ end documentation block
  use mod_type,           only : rdef
  use mod_lorenz63_fwd,   only : lorenz63_rk4, t_lorenz63
  use mod_rnorm,          only : rnorm, set_random_seed
  use mod_lorenz63_inc4dvar, only : lorenz63_inc4dvar
  implicit none

  integer,parameter :: nx = 3

  integer          :: ncycle, nspin
  integer          :: kitermax_in, kitermax_out
  type(t_lorenz63) :: xt, xb, xa
  real(rdef)       :: dx(nx)
  real(rdef)       :: errb(nx,nx)
  real(rdef)       :: infl
  real(rdef)       :: erro(nx)
  real(rdef),allocatable :: yo(:,:) ! (ncycle,nx)

  real(rdef),allocatable :: erro4d(:,:) ! (nx,nsteps_da_window)
  real(rdef),allocatable :: yo4d(:,:) ! (nx,nsteps_da_window)
  logical,   allocatable :: luse4d(:,:) ! (nx,nsteps_da_window)
  integer :: nsteps_per_da, nsteps_da_window

  type(t_lorenz63) :: tmpx
  integer :: i, j, k


!----------------------------------------------------------------------------------
! 0. configuration for lorenz63 model
!----------------------------------------------------------------------------------
  xt%dt   = 0.01d0   ! length for one time step
  xt%nt   = 8       ! generate nt-step forecast
  xt%para = (/ 10.0d0, 28.0d0, 8.0d0/3.0d0 /)  

  nspin   = 4000     ! spin up the model with nspin steps, then the final state 
                    ! is used as the initial condition
  ncycle  = 4000     ! total n-step forecast cycles
  erro    = (/ 1.414d0, 1.414d0, 1.414d0 /)
  dx      = (/ 0.1d0, 0.1d0, 0.1d0 /) ! percent of error for the intial condition
  infl    = 1.  ! inflation factor in front of the raw B
  kitermax_in  = 100 ! maximum number of inner loop
  kitermax_out = 3  ! maximum number of outer loop

  nsteps_per_da = 2    ! perform 4d-var every nsteps_per_da steps
  nsteps_da_window = 2 ! for each 4d-var cycle, use obs up to  nsteps_da_window steps


!----------------------------------------------------------------------------------
! 1. spin up the lorenz 63 model
!----------------------------------------------------------------------------------
  xt%x0   = (/  0.1d0,  0.1d0, 0.1d0 /)
  Write(6,*) "spin up lorenz63 model: dt =", xt%dt, ", nspin =", nspin
  Call lorenz63_rk4( xt%x0(:), xt%xn(:), xt%para, nspin, xt%dt )
  xt%x0(:) = xt%xn(:)
  Write(6,*) "initial condition: x=", xt%x0(:)


!----------------------------------------------------------------------------------
! 2. generate simulated observations
!----------------------------------------------------------------------------------
  Write(6,*) "generate simulated obs: err(yo) = ", erro(:)
  Allocate( yo(ncycle,nx) )
  tmpx = xt
  !Call set_random_seed( 0 )
  Write(6,*) "write simulated obs (yo) into fort.10010"
  Do i = 1, ncycle
     Call lorenz63_rk4( tmpx%x0(:), tmpx%xn(:), tmpx%para, tmpx%nt, tmpx%dt )
     Do j = 1, nx
        yo(i,j) = tmpx%xn(j) + erro(j)*rnorm()
     Enddo
     Write(10000,*) ( tmpx%xn(j), j = 1, nx )
     Write(10010,*) ( yo(i,j),    j = 1, nx )
     tmpx%x0(:) = tmpx%xn(:)
  Enddo


!----------------------------------------------------------------------------------
! 3. run 4d-var at every cycle
!----------------------------------------------------------------------------------
  allocate(erro4d(nx,nsteps_da_window))
  allocate(yo4d(nx,nsteps_da_window))
  allocate(luse4d(nx,nsteps_da_window))

  Write(6,*) "load background error covariance matrix B(3x3) from fort.1040"
  Do i = 1, nx
     Read(1040,*) ( errb(i,j), j = 1, nx )
     Write(6,*)   ( errb(i,j), j = 1, nx )
  Enddo
  Write(6,*) "rawB="
  Do i = 1, nx
     Write(6,*)   ( errb(i,j), j = 1, nx )
  Enddo

  Write(6,*) "infl=",infl
  errb = infl*errb
  Write(6,*) "finalB="
  Do i = 1, nx
     Write(6,*)   ( errb(i,j), j = 1, nx )
  Enddo

  Write(6,*) "run 3d-var at every cycle"
  xb = xt
  Do i = 1, nx
     xb%x0(i) = xt%x0(i) + dx(i)*xt%x0(i)*rnorm()  ! perturb the initial condition
  Enddo
  print*, xt%x0(:)
  print*, xb%x0(:)
  Write(6,*) "write truth (xt)     into fort.10020"
  Write(6,*) "write background(xb) into fort.10030"
  Write(6,*) "write analysis (xa)  into fort.10040"
  Do i = 1, ncycle

     ! true trajectory
     Call lorenz63_rk4( xt%x0(:), xt%xn(:), xt%para, xt%nt, xt%dt )
     ! background
     Call lorenz63_rk4( xb%x0(:), xb%xn(:), xt%para, xt%nt, xt%dt )

     if ( mod(i, nsteps_per_da) == 0 .and. i+nsteps_da_window-1 <= ncycle ) then
         If ( MOD(i,1) == 0 ) Write(6,*)  "run 4d-var for cycle", i
         ! run inc3d-var

         do k = 1, nsteps_da_window
            do j = 1, nx
               erro4d(j,k) = erro(j)
               yo4d(j,k)   = yo(i+k-1,j)
               luse4d(j,k) = .true.
            enddo
         enddo
         Call lorenz63_inc4dvar( nx            = nx, &
                                 nyo           = nx, &
                                 nda           = nsteps_da_window, &
                                 xb            = xb%xn(:), &
                                 B             = errb, &
                                 lyo           = luse4d, &
                                 yo            = yo4d, &
                                 erro          = erro4d, &
                                 xa            = xa%xn(:), & 
                                 kitermax_in   = kitermax_in, &
                                 kitermax_out  = kitermax_out, &
                                 dt            = xt%dt, &
                                 ndt           = xt%nt, &
                                 para          = xt%para )

     else

         xa%xn(:) = xb%xn(:) 

     endif

     Write(10020,*) ( xt%xn(j), j = 1, nx )
     Write(10030,*) ( xb%xn(j), j = 1, nx )
     Write(10040,*) ( xa%xn(j), j = 1, nx )

     xt%x0(:) = xt%xn(:)
     xb%x0(:) = xa%xn(:)

     !print*, "i=", i
     !print*, "bak  =", xb%xn(:)
     !print*, "truth=", xt%xn(:)
     !print*, "ana  =", xa%xn(:)
     !pause(1)


     
  Enddo

  Deallocate( yo ) 

endprogram

