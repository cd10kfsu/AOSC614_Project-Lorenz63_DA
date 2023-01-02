program test_fwd
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
  implicit none

  integer,parameter :: nx = 3

  integer          :: ncycle, nspin
  type(t_lorenz63) :: xt
  integer :: i, j


!----------------------------------------------------------------------------------
! 0. configuration for lorenz63 model
!----------------------------------------------------------------------------------
  xt%dt   = 0.01d0   ! length for one time step
  xt%nt   = 1       ! generate nt-step forecast
  xt%para = (/ 10.0d0, 28.0d0, 8.0d0/3.0d0 /)  
  !xt%para = (/ 10.0d0, 400.d0, 8.0d0/3.0d0 /)  
  nspin   = 4000     ! spin up the model with nspin steps, then the final state 
                    ! is used as the initial condition
  ncycle  = 16000     ! total n-step forecast cycles

!----------------------------------------------------------------------------------
! 1. spin up the lorenz 63 model
!----------------------------------------------------------------------------------
  xt%x0   = (/  0.1d0,  0.1d0, 0.1d0 /)
  Write(6,*) "spin up lorenz63 model: dt =", xt%dt, ", nspin =", nspin
  Call lorenz63_rk4( xt%x0(:), xt%xn(:), xt%para, nspin, xt%dt )
  xt%x0(:) = xt%xn(:)
  Write(6,*) "initial condition: x=", xt%x0(:)

!----------------------------------------------------------------------------------
! 2. free run
!----------------------------------------------------------------------------------
  Write(6,*) "write truth (xt)     into fort.10020"
  Do i = 1, ncycle
     ! true trajectory
     Call lorenz63_rk4( xt%x0(:), xt%xn(:), xt%para, xt%nt, xt%dt )
     Write(10020,*) ( xt%xn(j), j = 1, nx )
     xt%x0(:) = xt%xn(:)
  Enddo

endprogram

