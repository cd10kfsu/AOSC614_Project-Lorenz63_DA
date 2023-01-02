module mod_lorenz63_letkf
!$$$  program documentation block
!         .           .            .
!  program name: mod_lorenz63_ensrf
!    programmer: da,cheng        org: umd      date: 2016-Mar-23
!
!  purpose:
!    Local ensemble transform kalman filter for lorenz63 system
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
  use mod_math, only : eigenmtx
  implicit none

  private 
  public :: lorenz63_letkf

  integer,parameter :: nx = 3
  integer,parameter :: nyo = 3


contains

subroutine lorenz63_letkf( nn, xb, lyo, yo, erro, xa, xam )
  implicit none
! passed args
  integer,    intent(in   ) :: nn
  real(rdef), intent(in   ) :: xb(nx,nn)
  logical,    intent(in   ) :: lyo(nyo)
  real(rdef), intent(in   ) :: yo(nyo)
  real(rdef), intent(in   ) :: erro(nyo)
  real(rdef), intent(  out) :: xa(nx,nn)
  real(rdef), intent(  out) :: xam(nx)
! local vars
  integer :: nyo_use
  real(rdef),allocatable :: dyb(:,:)
  real(rdef),allocatable :: ybm(:)
  real(rdef),allocatable :: erro_use(:)
  real(rdef),allocatable :: yo_ybm(:)
  real(rdef) :: xbm(nx)
  real(rdef) :: invSPa(nn,nn), SPa(nn,nn) ! Pa in ensemble space
  real(rdef),allocatable :: C(:,:)
  real(rdef) :: eigval(nn), eigvect(nn,nn)
  integer :: np
  real(rdef) :: wam(nn), Wa(nn,nn)
  real(rdef) :: dxb(nn)
  integer :: ierr
  integer :: i, j, k

  nyo_use = 0
! calculate mean yb and perturbation matrix Yb
  do i = 1, nyo
     if ( lyo(i) ) then
        nyo_use = nyo_use + 1
     else
        write(6,*) "skip obs: iob =", i, ", yo(iob)=", yo(i)
     endif
  enddo
  allocate( dyb(nyo_use,nn) )
  allocate( ybm(nyo_use), yo_ybm(nyo_use), erro_use(nyo_use) )

  k = 0
  do i = 1, nyo
     if ( lyo(i) ) then
        k           = k + 1
        dyb(k,:)    = xb(i,:)
        ybm(k)      = SUM(dyb(k,:))/nn
        dyb(k,:)    = dyb(k,:) - ybm(k)
        yo_ybm(k)   = yo(i) - ybm(k)
        erro_use(k) = erro(i)
        !print*, "k, erro, yo, ybm=", k, erro_use(k), yo(i), ybm(k)
        !pause
     endif
  enddo

! calculate mean xb
  do i = 1, nx
     xbm(i) = SUM(xb(i,:))/nn
  enddo
 
  allocate( C(nn,nyo_use) )

! C=(Yb)^T * R^-1
  do j = 1, nyo_use
     C(:,j) = dyb(j,:)/(erro_use(j)**2)
  enddo
! invSPa = [(k-1)I+C*Yb]
  invSPa = MATMUL( C, dyb )
  do k = 1, nn
     invSPa(k,k) = invSPa(k,k)+(nn-1)
  enddo
! let A=invSPa, where A=Q*V*Q^T
! so A^-1=Q*V^-1*Q^T
  call eigenmtx( nn, invSPa, eigvect, eigval, np, ierr )
  if ( ierr/=0 .or. np/=nn ) then
     write(6,*) "[warning] lorenz63_letkf: fail to find nn (+) eigenvetor"
     xam = xbm
     xa = xb
     return
  endif
  do k = 1, nn
     SPa(:,k) = eigvect(:,k)/eigval(k)
  enddo
  SPa = MATMUL( SPa, TRANSPOSE(eigvect) )
  ! Wa = [(k-1)SPa]^1/2
  !    = [(k-1)A^-1]^1/2
  !    = [ Q*(k-1)V^-1*Q^T ]^1/2
  !    = [ Q*sqrt((k-1)/V)*Q^T ]
  do k = 1, nn
     Wa(:,k) = eigvect(:,k)*SQRT( (nn-1)/eigval(k) )
  enddo
  Wa = MATMUL( Wa, TRANSPOSE(eigvect) )
  ! wam = SPa * C *( yo -ybm )
  wam = MATMUL( C, yo_ybm )
  wam = MATMUL( SPa, wam )  !!!
  ! Wa = Wa + wam
  do k = 1, nn
     Wa(:,k) = Wa(:,k) + wam(:)
  enddo

  do i = 1, nx
     ! xam = xbm + Xb*wam
     dxb = xb(i,:)-xbm(i)
     xam(i) = xbm(i) + SUM(dxb(:)*wam(:))

     do k = 1, nn
        xa(i,k) = xbm(i) + SUM(dxb(:)*Wa(:,k))
     enddo
  enddo

endsubroutine

endmodule
