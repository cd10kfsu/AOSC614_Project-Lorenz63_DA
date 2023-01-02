module mod_math
!$$$  program documentation block
!         .           .            .
!  program name: mod_math
!    programmer: da,cheng        org: umd      date: 2015-Jan-08
!
!  purpose:
!    operations related to matrix manipulations. Based on Lapack
!
!
!  revision history:
!    2015-Jan-08     da    - creator
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
  public covar, invmtx, eigenmtx, sqrtmtx
  public eye, qr
  public :: inv

contains

!
! QR decomposition of a square matrix
!
subroutine QR(n,A,Q,R)
  implicit none

  integer,intent(in) :: n
  real(rdp),intent(in) :: A(n,n)
  real(rdp),intent(out) :: Q(n,n), R(n,n)

  integer :: ierr
  real(rdp) :: work(n), bufA(n,n), tau(n)
  integer :: lda
  integer :: i

  ierr=0
  lda = n
  bufA = A
  call dgeqrf(n,n,bufA,lda,tau,work,n,ierr)
  if (ierr/=0) stop "ERROR in QR: dgetqrf"
  R=0.d0
  do i = 1, n
     R(i,i:n) = bufA(i,i:n)
  enddo
  Q=bufA
  call dorgqr(n,n,n,Q,lda,tau,work,n,ierr)
  if (ierr/=0) stop "ERROR in QR: qorgqr"

endsubroutine

!
! identity matrix
!
subroutine EYE(n,A)
  implicit none

  integer,intent(in) :: n
  real(rdp),intent(out) :: A(n,n)

  integer :: i

  A=0.d0
  do i = 1, n
     A(i,i) = 1.d0
  enddo

endsubroutine


function covar( n, x, y ) result ( cov )
  implicit none
! passed args
  integer,  intent(in) :: n
  real(rdef),intent(in) :: x(n)
  real(rdef),intent(in) :: y(n)
  real(rdef)            :: cov
! local vars
  real(rdef) :: xm, ym

  xm = SUM(x(:))/n
  ym = SUM(y(:))/n
  cov = SUM( ( x(:)-xm )*( y(:) - ym ) )/(n-1)

endfunction

subroutine INV(n,A,invA)
  implicit none

  integer,intent(in) :: n
  real(rdef),intent(in) :: A(n,n)
  real(rdef),intent(out) :: invA(n,n)

  integer :: info, ipiv(n)
  integer :: lwork, lda 
  real(rdef) :: work(n)

  invA = A 
  lda  = n; lwork = n 
  call dgetrf(n,n,invA,lda,ipiv,info)
  if (info/=0) stop "ERROR in INV: dgetrf"
  call dgetri(n,invA,lda,ipiv,work,lwork,info)
  if (info/=0) stop "ERROR in INV: dgetri"
endsubroutine

subroutine invmtx( n, A, B, ierr )
  implicit none
! passed args
  integer,    intent(in   ) :: n
  real(rdef), intent(in   ) :: A(n,n)
  real(rdef), intent(  out) :: B(n,n)
  integer,    intent(  out) :: ierr
! local vars
  real(rdp),allocatable :: tmpA(:,:)
  integer :: i,j 

  ierr = 0
  Allocate( tmpA(n,n) )
  tmpA = REAL(A, rdp )

! get inverse with Cholesky factorization. subroutines from 
! lapack
  Call dpotrf( 'U', n, tmpA, n, ierr )
  if ( ierr /= 0 ) return
  Call dpotri( 'U', n, tmpA, n, ierr )
  if ( ierr /= 0 ) return

  B = tmpA
  Do i = 2, n
     Do j = 1, i-1
        B(i,j) = B(j,i) 
     Enddo
  Enddo

  Deallocate( tmpA )

endsubroutine


subroutine sqrtmtx( n, A, L )
  implicit none

  integer,intent(in) :: n
  real(8),intent(in) :: A(n,n)
  real(8),intent(out) :: L(n,n)

  real(8) :: S(n,n), D1d(n)
  integer :: i, np, ierr

  call eigenmtx(n, A, S, D1d, np, ierr)
  do i = 1, n
     L(:,i) = S(:,i)*sqrt(D1d(i))
  enddo

endsubroutine


subroutine eigenmtx( n, A, eigvect, eigval, np, ierr )
!
! check "dsyev" in the LAPACK
!
implicit none
! passed args
  integer,   intent(in   ) :: n
  real(rdef),intent(in   ) :: A(n,n)
  real(rdef),intent(  out) :: eigvect(n,n)
  real(rdef),intent(  out) :: eigval(n)
  integer,   intent(  out) :: np
  integer,   intent(  out) :: ierr
! local vars
  real(rdp)         :: r8A(n,n)
  real(rdp)         :: r8eigvect(n,n)
  real(rdp)         :: r8eigval(n)
  integer           :: lwk
  integer,parameter :: lwkmax = 3000 
  real(rdp)         :: wk(lwkmax)
  integer           :: i

  ierr = 0
  r8A      = REAL(A, rdp)
  r8eigvect = 0.0d0
  r8eigval = 0.0d0
  wk       = 0.0d0

  if ( 3*n-1 > lwkmax ) then
      write(6,*) "[error] eigenmtx: increase lwkmax to", 3*n-1 
      ierr = -1
      stop
  endif
! get the optimal work space
  lwk = -1
  call dsyev( 'V', 'U', n, r8A, n, r8eigval, wk, lwk, ierr )
  lwk = MIN( lwkmax, INT(wk(1)) )
! solve eigenproblem
  call dsyev( 'V', 'U', n, r8A, n, r8eigval, wk, lwk, ierr )
  if ( ierr /= 0 ) then
     write(6,*) "[warning] eigenmtx: fail to solve eigenvalue/vector"
     return
  endif
  np = 0
  do i = 1, n
     if ( r8eigval(i)>0 ) np = np + 1
  enddo
  !if ( np<n ) then
  !   write(6,*) "[warning] eigenmtx: have nonpositive eigenvalues"
  !endif

! put the largest eigenvalue as the 1st column, which brings the largest variance
  do i = 1, n
     eigval(i)    = r8eigval(n+1-i)
     eigvect(:,i) = r8A(:,n+1-i)
  enddo

endsubroutine


endmodule

