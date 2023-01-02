module mod_lorenz63_oi
!$$$  program documentation block
!         .           .            .
!  program name: mod_lorenz63_oi
!    programmer: da,cheng        org: umd      date: 2015-Jan-18
!
!  purpose:
!    Optimal interpolation system for the Lorenz63 system. Note this
!    system directly solve
!     Xa = Xb + B * H^T * (HBH^T+R)^-1 * ( Yobs - h(x) ) 
!
!    This system is intended to verify the lorenz63 3d-var system
!
!
!  revision history:
!    2015-Jan-18     da    - creator
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
  !use mod_math, only : invmtx, inv
  use mod_math, only : inv
  implicit none

  private 
  public :: nx, nyo
  public :: lorenz63_oi

  integer,parameter :: nx = 3
  integer,parameter :: nyo = 3


contains

subroutine lorenz63_oi( xb, B, yo, erro, xa )
  implicit none
!
! Assume obs operator is the unit array
!
! passed args
  real(rdef),intent(in   ) :: xb(nx)
  real(rdef),intent(in   ) :: B(nx,nx)
  real(rdef),intent(in   ) :: yo(nyo)
  real(rdef),intent(in   ) :: erro(nyo)
  real(rdef),intent(  out) :: xa(nx)
! local args
  real(rdef) :: R(nx,nx)
  real(rdef),allocatable :: r1d(:)
  real(rdef),allocatable :: r2d(:,:), r2d2(:,:)
  integer    :: ierr
  integer    :: i

! get Yobs - h(x)
  Allocate( r1d(nyo) )
  r1d(:) = yo(:) - xb(:) 
! get HBH^T + R where H is eye(3)
  R(:,:) = 0.0d0
  Do i = 1, nyo
     R(i,i) = erro(i)**2
  Enddo
  Allocate( r2d(nyo,nyo) )
  r2d(:,:) = B(:,:) + R(:,:)
! get (HBH^T + R)^-1
  Allocate( r2d2(nyo,nyo) )
  !Call invmtx( nyo, r2d, r2d2, ierr )
  Call inv( nyo, r2d, r2d2)
  ierr=0
  If ( ierr /= 0 ) Then 
     Write(6,*) "[warning] lorenz63_oi: fail to invert (HBH^T+R)"
     xa = xb
  Else 
!    get BH^T(HBH^T + R)^-1
     r2d = Matmul( B, r2d2 )
!    get analysis BH^T(HBH^T + R)^-1 * (Yobs-h(x))
     xa = xb + Matmul( r2d, r1d )
  Endif

  Deallocate( r1d, r2d, r2d2 )

Endsubroutine


Endmodule
