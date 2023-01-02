module mod_lorenz63_enkf
!$$$  program documentation block
!         .           .            .
!  program name: mod_lorenz63_enkf
!    programmer: da,cheng        org: umd      date: 2016-Mar-23
!
!  purpose:
!    ENKF system for the Lorenz63 system. Here we use perturbed obs
!    In addition, we assimilate obs serially
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
  use mod_rnorm, only : rnorm, set_random_seed
  use mod_math, only : invmtx
  implicit none

  private 
  public :: lorenz63_enkf, lorenz63_enkf_m, lorenz63_enkf_mupdate


  integer,parameter :: nx = 3
  integer,parameter :: nyo = 3


contains
subroutine lorenz63_enkf_mupdate( xb, B, yo, erro, xa )
! update one ensemble member with matrix inversion
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
  Call invmtx( nyo, r2d, r2d2, ierr )
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



subroutine lorenz63_enkf_m( nn, xb, Pb, yo, erro, xa )
! this is a quick check if my lorenz63_enkf works correctly
! the matrix inversion is used here to get kalman gain
  implicit none

  integer,intent(in) :: nn
  real(rdef),intent(in) :: xb(nx,nn)
  real(rdef),intent(in) :: Pb(nx,nx)
  real(rdef),intent(in) :: yo(nyo)
  real(rdef),intent(in) :: erro(nyo)
  real(rdef),intent(out) :: xa(nx,nn)
  
  real(rdef),allocatable :: dyo(:)
  real(rdef),allocatable :: yoens(:,:)
  real(rdef) :: R(nyo,nyo),dhxb(nyo)
  real(rdef) :: r2d(nyo,nyo), r2d2(nyo,nyo)
  real(rdef) :: kfgain(nx,nyo)
  integer :: i, k
  integer :: ierr

  allocate( yoens(nyo,nn) )
  !yoens = yo
  allocate( dyo(nn) )
  yoens = 0.0d0
  do i = 1, nyo
     do k = 1, nn
        dyo(k) = erro(i)*rnorm()
     enddo
     dyo = dyo - SUM(dyo)/nn
     yoens(i,:) = yo(i) +dyo
     !print*, "[lorenz63_enkf_m] msg: E<yp(i,:)>, yo=", &
     !        SUM( yoens(i,:) )/nn, yo(i)
  enddo

  R(:,:) = 0.0_rdef
  do i = 1, nyo
     R(i,i) = erro(i)**2
  enddo
  r2d(:,:) = Pb(:,:) + R(:,:)
  call invmtx( nyo, r2d, r2d2, ierr )
  if ( ierr/=0 ) then
     write(6,*) "[lorenz63_enkf_m] warning: fail to invert matrix"
     xa(:,:) = xb(:,:)
     return
  else
     kfgain = MATMUL( Pb, r2d2 )
     do k = 1, nn
        dhxb(:) = yoens(:,k) - xb(:,k)
        xa(:,k) = xb(:,k) + MATMUL( kfgain, dhxb )
     enddo
  endif

endsubroutine
   
    
subroutine lorenz63_enkf( nn, xb, lyo, yo, erro, xa, xam )
! assimilate obs one at a time
  implicit none
! passed args
  integer,   intent(in   ) :: nn        ! number of ensembles
  real(rdef),intent(in   ) :: xb(nx,nn) 
  logical,   intent(in   ) :: lyo(nyo)  
  real(rdef),intent(in   ) :: yo(nyo)
  real(rdef),intent(in   ) :: erro(nyo) ! assume R is diagonal
  real(rdef),intent(  out) :: xa(nx,nn)
  real(rdef),intent(  out),optional:: xam(nx)
! local vars
  real(rdef),allocatable :: dyo(:)      ! obs perturbation
  real(rdef),allocatable :: dxb(:,:), hdxf(:)
  real(rdef),allocatable :: xb2(:,:)
  real(rdef) :: hpht
  real(rdef) :: pht(nx)
  real(rdef) :: kfgain(nx)
  integer :: iob, j, k

  allocate( xb2(nx,nn) )
  allocate( dyo(nn) )
  allocate( dxb(nx,nn) )
  allocate( hdxf(nn) )

  xb2(:,:) = xb(:,:)
  xa(:,:) = xb(:,:)

  ! in case all obs are skipped
  xa(:,:) = xb2(:,:)
  if (present(xam)) then
  do k =1 , nx
     xam(k) = SUM(xb2(k,:))/nn
  enddo
  endif

  do iob = 1, nyo
     if ( .not. lyo(iob) ) then
        write(6,*) "skip obs: iob =", iob, ", yo(iob)=", yo(iob)
        cycle  ! cycle if obs for xi is not available
     endif

     ! impose simple QC
     !if ( ABS( SUM(xb2(iob,:))/nn-yo(iob) ) > 5*erro(iob) ) then
     !   write(6,*) "reject obs: "
     !   write(6,*) "iob, yo, E<Hx>=", iob, yo(iob), SUM(xb2(iob,:))/nn
     !   write(6,*) "|yo-E<Hx>|, 5*std=", &
     !              ABS( SUM(xb2(iob,:))/nn-yo(iob) ), &
     !              5*erro(iob)
     !   cycle
     !endif

     ! perturb obs
     do j = 1, nn
        dyo(j) = erro(iob)*rnorm()
     enddo
     ! ensure the mean is zero
     dyo(:) = dyo(:) - SUM( dyo(:) )/nn

     ! HPH^T = E<(Hx-E<Hx>)(Hx-E<Hx>)>
     hdxf = xb2(iob,:) - SUM(xb2(iob,:))/nn
     hpht = SUM(hdxf*hdxf)/(nn-1)
     ! PH^T = E<(x-Ex)(Hx-E<hx>)>
     do k = 1, nx
        dxb(k,:) = xb2(k,:) - SUM(xb2(k,:))/nn
     enddo
     pht = MATMUL(dxb,hdxf)/(nn-1)
     
     ! K=PH^T (HPH^T+R)^-1
     kfgain(:) = pht/( hpht + erro(iob)**2 )
     
     ! get the analysis xa for each ensemble member
     ! Xa = Xb + K(Yo - HXb)
     do j = 1, nn
        xa(:,j) = xb2(:,j) + kfgain*( yo(iob)+dyo(j)-xb2(iob,j) )
     enddo

     ! save analysis mean
     if (present(xam)) then
     do k = 1, nx
        xam(k) = sum(xa(k,:))/nn
     enddo
     endif
     
     xb2(:,:) = xa(:,:)
  enddo
  deallocate( dyo, dxb, hdxf,xb2 )

endsubroutine


endmodule
