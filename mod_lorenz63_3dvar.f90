module mod_lorenz63_3dvar
!$$$  program documentation block
!         .           .            .
!  program name: mod_lorenz63_3dvar
!    programmer: da,cheng        org: umd      date: 2015-Jan-11
!
!  purpose:
!    3D-Var system for the Lorenz63 system
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
  use mod_type, only : rdef, rdp
  use mod_math, only : invmtx
  Use mod_optimization, only : run_lbfgs
  implicit none

  private 
  public :: lorenz63_3dvar


  integer,parameter :: nx = 3
  integer,parameter :: nyo = 3
  integer,parameter :: kitermax = 100
  real(rdp),parameter :: epsl   = 1.d-6


contains

subroutine lorenz63_3dvar( xb, B, lyo, yo, erro, xa )
  implicit none
! passed args
  real(rdef),intent(in   ) :: xb(nx)
  real(rdef),intent(in   ) :: B(nx,nx)
  logical,   intent(in   ) :: lyo(nyo)
  real(rdef),intent(in   ) :: yo(nyo)
  real(rdef),intent(in   ) :: erro(nyo)
  real(rdef),intent(  out) :: xa(nx)
! local vars
  real(rdef) :: invB(nx,nx), invR(nyo)
  real(rdp) :: x(nx)
  real(rdp) :: Jc, Jo, Jb
  real(rdp) :: dJc(nx), dJo(nx), dJb(nx)

  integer :: kiter

  integer :: ierr
  integer :: i, j

! get B^-1, and R^-1
  Call invmtx( nx, B, invB, ierr )

  Write(6,*) "B ="
  Do i = 1, nx
     Write(6,*)   ( B(i,j), j = 1, nx )
  Enddo

  Write(6,*) "inv(B) ="
  Do i = 1, nx
     Write(6,*)   ( invB(i,j), j = 1, nx )
  Enddo

  If ( ierr /= 0 ) Stop "[lorenz63_3dvar] error: fail to get inverse of B"
  Do i = 1, nyo
     If ( lyo(i) ) Then
        invR(i) = 1.d0/( erro(i)*erro(i) )
     Endif
  Enddo

! Set xa=xb at the inital step
  x = xb

  kiter = 0
! Minimization loop
  lp_mini: Do 

    ! Calcualte Jb = 0.5 * ( x - xb) * B^-1 * ( x - xb )
    Jb = 0.0d0
    Do i = 1, nx
       Do j = 1, nx
          Jb = Jb + 0.5d0*( x(i)-xb(i) )*( x(j)-xb(j) )*invB(i,j) 
       Enddo
    Enddo
    ! Calculate Jo = 0.5 * ( h[x] - yo ) * R^-1 * ( h[x] - yo )
    Jo = 0.0d0
    Do i = 1, nyo
       If ( lyo(i) ) Then
          Jo = Jo + 0.5d0*( x(i)-yo(i) )*( x(i)-yo(i) )*invR(i)
       Endif
    Enddo
    ! Calculate total cost function Jc = Jb + Jo
    Jc = Jb + Jo

    ! Calculate Gradient dJb = B^-1 * ( x - xb )
    dJb = 0.0d0
    Do i = 1, nx
       Do j = 1, nx
          dJb(i) = dJb(i) + invB(i,j) * ( x(j)-xb(j) ) 
       Enddo
    Enddo
    ! Calculate Gradient dJo = H^T * R^-1 * ( h[x] - yo )
    dJo = 0.0d0
    Do i = 1, nx
       If ( lyo(i) ) Then
          dJo(i) = invR(i) * ( x(i) - yo(i) ) 
       Endif
    Enddo
    ! Calculate Total gradient dJc = dJo + dJb
    dJc = dJo + dJb

    ! 
    Call run_lbfgs( nx, x, Jc, dJc, epsl, ierr )
    If ( ierr < 0 ) Stop "[lorenz63_3dvar] error:fail to run lbfgs."
    kiter = kiter + 1
    If ( kiter > kitermax ) Then
       Write(6,*) "[lorenz63_3dvar] warning: maximum iter reached. "
       exit lp_mini
    Endif
    If ( ierr == 0 ) exit lp_mini

  Enddo lp_mini

  xa = x

endsubroutine


endmodule
