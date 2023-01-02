module mod_inflation
!$$$  program documentation block
!         .           .            .
!  program name: mod_lorenz63_ensrf
!    programmer: da,cheng        org: umd      date: 2016-Mar-23
!
!  purpose:
!
!  revision history:
!    2016-Mar-23     da    - creator
!
!  file dependencies:
!
!  attributes: 
!    language: fortran 90
!    machine : 
!
!
!$$$ end documentation block
  use mod_type, only : rdef, rdp
  implicit none

  private 
  public :: inflate_rtps, inflate_rtpp


  integer,parameter :: nx = 3


contains

subroutine inflate_rtpp( nx, nn, xb, xa, alpha )
  implicit none
! passed args
  integer,   intent(in   ) :: nx
  integer,   intent(in   ) :: nn
  real(rdef),intent(in   ) :: xb(nx,nn)
  real(rdef),intent(inout) :: xa(nx,nn)
  real(rdef),intent(in   ) :: alpha
! local vars
  real(rdef) :: xam
  real(rdef),allocatable :: dxb(:), dxa(:)
  integer :: i

  allocate( dxb(nn) )
  allocate( dxa(nn) )
  do i = 1, nx
     dxb(:) = xb(i,:) - SUM(xb(i,:))/nn
     xam = SUM(xa(i,:))/nn
     dxa(:) = xa(i,:) - xam
     ! inflate xa perturbation
     dxa(:) = (1-alpha)*dxa(:) + alpha*dxb(:)
     ! add to xa mean state
     xa(i,:) = xam + dxa(:)
  enddo
  deallocate( dxa, dxb )

endsubroutine


subroutine inflate_rtps( nx, nn, xb, xa, alpha )
  implicit none
! passed args
  integer,   intent(in   ) :: nx
  integer,   intent(in   ) :: nn
  real(rdef),intent(in   ) :: xb(nx,nn)
  real(rdef),intent(inout) :: xa(nx,nn)
  real(rdef),intent(in   ) :: alpha
! local vars
  real(rdef) :: xam
  real(rdef) :: bsprd
  real(rdef) :: asprd
  real(rdef) :: coef
  real(rdef),allocatable :: dx(:)
  integer :: i

  allocate( dx(nn) )
  do i = 1, nx
     ! calculate prior spread
     dx(:) = xb(i,:) - SUM(xb(i,:))/nn
     bsprd = SQRT( SUM(dx(:)**2)/(nn-1) )
     ! calculate posterior spread
     xam = SUM(xa(i,:))/nn
     dx(:) = xa(i,:) - xam
     asprd = SQRT( SUM(dx(:)**2)/(nn-1) )
     ! inflate xa perturbation
     coef = ( alpha*( bsprd - asprd )/asprd + 1.0_rdef )
     if ( coef < 1.0_rdef ) then
        print*, "[msg] inflate_rtps: coef < 1 for x(i), i, coef=", i, coef
        coef = 1.0_rdef
     endif
     dx(:) = dx(:) * coef
     ! add to xa mean state
     xa(i,:) = xam + dx(:)
  enddo
  deallocate( dx )

endsubroutine


endmodule
