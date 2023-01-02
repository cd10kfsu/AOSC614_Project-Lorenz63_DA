program test_lv
!$$$  program documentation block
!         .           .            .
!  program name: test_lv
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
  use mod_lorenz63_fwd,   only : lorenz63_rk4, t_lorenz63, step, get_tlm
  use mod_lorenz63_lv,    only : le_blv, le_flv 
  use mod_rnorm,          only : rnorm1d, rnorm
  implicit none

  integer,parameter :: nx = 3

  integer          :: ncycle, nspin
  type(t_lorenz63) :: xt
  real(rdef),allocatable :: M(:,:,:), Q(:,:)
  real(rdef),allocatable :: le(:)
  real(rdef) :: buft
  integer :: i, j


!----------------------------------------------------------------------------------
! 0. configuration for lorenz63 model
!----------------------------------------------------------------------------------
  xt%dt   = 0.008d0   ! length for one time step
  xt%nt   = 1       ! generate nt-step forecast
  xt%para = (/ 10.0d0, 28.0d0, 8.0d0/3.0d0 /)  

  nspin   = 10000     ! spin up the model with nspin steps, then the final state 
                    ! is used as the initial condition
  ncycle  = 400000     ! total n-step forecast cycles


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
  allocate(Q(nx,nx))
  do i = 1, nx
     do j = 1, nx
        Q(i,j) = rnorm()
     enddo
  enddo
  allocate(M(nx,nx,ncycle))
  allocate(le(nx))
  do i = 1, ncycle
     call get_tlm(xt%x0(:),xt%para,xt%dt,M(:,:,i))
     call step(xt%x0(:), xt%para, xt%dt, buft, xt%xn(:))
     xt%x0(:) = xt%xn(:)
  enddo


  !call le_blv(nx,ncycle,xt%dt,M,le, Q)
  call le_blv(nx,ncycle-1,xt%dt,M,le)
  print*, "bwd: le=", le, sum(le)

  do i = 1, ncycle
     M(:,:,i) = transpose(M(:,:,i))
  enddo
  !call le_flv(nx,ncycle,xt%dt,M,le, Q)
  call le_flv(nx,ncycle-1,xt%dt,M,le)
  print*, "fwd: le=", le, sum(le)




endprogram

