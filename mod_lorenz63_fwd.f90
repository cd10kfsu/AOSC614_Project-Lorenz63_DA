module mod_lorenz63_fwd
!$$$  program documentation block
!         .           .            .
!  program name: mod_lorenz63_fwd
!    programmer: da,cheng        org: umd      date: 2015-Jan-05
!
!  purpose:
!    forward model for 3-variable lorenz model.
!    The module structure follows Miyoshi's lorenz96 module
!
!  revision history:
!    2015-Jan-05     da    - creator
!    2018-oct-12     Da,Cheng     - merge all schemes together
!
!  file dependencies:
!
!  attributes: 
!    language: fortran 90
!    machine : 
!
!
!$$$ end documentation block

  use mod_type, only: r_kind => rdef
  implicit none

  private 
!------public vars & subs
  integer,     parameter,public :: nx              = 3
  real(r_kind),parameter,public :: para_default(3) = (/ 10.d0, 28.d0, 8.d0/3 /)

  public :: t_lorenz63
  public :: step, tl_step, adj_step                
  public :: lorenz63_rk4, tl_lorenz63_rk4, adj_lorenz63_rk4
  !public :: step_rk4, tl_step_rk4, adj_step_rk4   ! RK-4 scheme
  !public :: step_ft,  tl_step_ft,  adj_step_ft    ! Forward in time scheme
  public :: get_tlm
  !public :: get_tlm_1                             ! exact TLM
  !public :: get_tlm_2                             ! approximated TLM with perturbation method
  !public :: get_tlm_3                             ! appromxiated TLM directly from analytical jacobian
  !public :: jacobian
  !public :: tendency, tl_tendency, adj_tendency 

!------private vars & subs
  real(r_kind),private :: dx = 1.D-10

!------solver selection
  interface get_tlm 
    module procedure get_tlm_1 
    !module procedure get_tlm_2 
    !module procedure get_tlm_3
  endinterface

  interface step
    module procedure step_rk4
    !module procedure step_ft
  endinterface

  interface tl_step
    module procedure tl_step_rk4
    !module procedure tl_step_ft
  endinterface

  interface adj_step
    module procedure adj_step_rk4
    !module procedure adj_step_ft
  endinterface

  type t_lorenz63
    integer   :: nt
    real(r_kind) :: dt
    real(r_kind) :: para(nx)
    real(r_kind) :: x0(nx)
    real(r_kind) :: xn(nx)
  endtype


!subroutine lorenz63_rk4( xin, xout, para, nt, dt )

contains

!-----------------------------------------------------------------------------------------
!  output explicit TLM matrix
!-----------------------------------------------------------------------------------------
!
! exact TLM
!
subroutine get_tlm_1(xin,para,dt,M)
  implicit none

  real(r_kind),intent(in)  :: xin(nx)
  real(r_kind),intent(in)  :: para(nx)
  real(r_kind),intent(in)  :: dt
  real(r_kind),intent(out) :: M(nx,nx)

  real(r_kind) :: xout(nx), xinp(nx), xoutp(nx)
  real(r_kind) :: buft
  integer :: n

  call step(xin, para, dt, buft, xout)
  do n = 1, nx
     xinp=0.d0
     xinp(n) = 1.d0
     call tl_step( xin,  xinp, para, dt, buft, xout, xoutp )
     M(:,n) = xoutp
  enddo

endsubroutine

!
! use small perturbed run to approximate TLM
!
subroutine get_tlm_2(xin,para,dt,M)
  implicit none

  real(r_kind),intent(in)  :: xin(nx)
  real(r_kind),intent(in)  :: para(nx)
  real(r_kind),intent(in)  :: dt
  real(r_kind),intent(out) :: M(nx,nx)

  real(r_kind) :: xout(nx), xinp(nx), xoutp(nx)
  real(r_kind) :: buft
  integer :: n

  call step(xin, para, dt, buft, xout)
  do n = 1, nx
     xinp = xin
     xinp(n) = xinp(n) + dx
     call step(xinp, para, dt, buft, xoutp)
     M(:,n) = (xoutp-xout)/dx
  enddo

endsubroutine

!
! use jacobian & fwd scheme to approximate TLM
!
subroutine get_tlm_3(xin,para,dt,M)
  implicit none

  real(r_kind),intent(in)  :: xin(nx)
  real(r_kind),intent(in)  :: para(nx)
  real(r_kind),intent(in)  :: dt
  real(r_kind),intent(out) :: M(nx,nx)

  real(r_kind) :: DJ(nx,nx), I(nx,nx)
  integer :: n

  I=0.d0
  do n = 1, nx
     I(n,n) = 1.d0
  enddo
  call jacobian( xin, para, DJ)
  M(:,:) = I+ DJ*dt

endsubroutine


subroutine jacobian(xin, para, Jac)
  implicit none

  real(r_kind),intent(in) :: xin(nx)
  real(r_kind),intent(in) :: para(nx)
  real(r_kind),intent(out) :: Jac(nx,nx)

  real(r_kind) :: sigma, rho, beta

  sigma = para(1)
  rho   = para(2)
  beta  = para(3)

  Jac(1,:) = (/-sigma, sigma, 0.d0/)
  Jac(2,:) = (/rho-xin(3), -1.d0, -xin(1)/)
  Jac(3,:) = (/xin(2), xin(1), -beta/)

endsubroutine

!-----------------------------------------------------------------------------------------
!  FWD, TL, ADJ of lorenz63 with Runge-Kutta 4th order method
!-----------------------------------------------------------------------------------------
! multi-step RK-4 FWD
!subroutine lorenz63_rk4( xin,  para, nt, dt, t, xout )
subroutine lorenz63_rk4( xin, xout, para, nt, dt )
  implicit none
! passed args
  real(r_kind),intent(in   ) :: xin(nx)
  real(r_kind),intent(in   ) :: para(nx)
  integer,     intent(in   ) :: nt
  real(r_kind),intent(in   ) :: dt
  !real(r_kind),intent(inout) :: t
  real(r_kind),intent(  out) :: xout(nx)
! local vars
  real(r_kind) :: t
  real(r_kind) :: cp_xin(nx)
  real(r_kind) :: x2(nx), x3(nx), x4(nx)
  real(r_kind) :: tend1(nx), tend2(nx), tend3(nx), tend4(nx)
  integer :: k

  cp_xin(:) = xin(:)
  do k = 1, nt
    ! RK4 1st stage 
      Call tendency( cp_xin, tend1, para )
      x2(:) = cp_xin(:) + 0.5d0 * dt * tend1(:)
    ! RK4 2nd stage
      Call tendency( x2, tend2, para )
      x3(:) = cp_xin(:) + 0.5d0 * dt * tend2(:)
    ! RK4 3rd stage
      Call tendency( x3, tend3, para )
      x4(:) = cp_xin(:) +         dt * tend3(:)
    ! RK4 4th stage
      Call tendency( x4, tend4, para )
      xout(:) = cp_xin(:) + &
            dt * ( tend1(:) + 2.0d0*tend2(:) + 2.0d0*tend3(:) + tend4(:) )/6.0d0
      t = t+dt

      cp_xin(:) = xout(:)
  enddo

endsubroutine


! multi-step RK-4 TL
!subroutine tl_lorenz63_rk4( xin9,  xin, para, nt, dt, t, xout9, xout )
subroutine tl_lorenz63_rk4( xin9,  xin, xout9, xout, para, nt, dt)
  implicit none
! passed args
  real(r_kind),intent(in   ) :: xin9(nx)
  real(r_kind),intent(in   ) :: xin(nx)
  real(r_kind),intent(in   ) :: para(nx)
  integer,     intent(in   ) :: nt
  real(r_kind),intent(in   ) :: dt
  !real(r_kind),intent(inout) :: t
  real(r_kind),intent(  out) :: xout9(nx)
  real(r_kind),intent(  out) :: xout(nx)
! local vars
  real(r_kind) :: t
  real(r_kind) :: cp_xin(nx)
  real(r_kind) :: cp_xin9(nx)
  real(r_kind) :: x2(nx), x3(nx), x4(nx)
  real(r_kind) :: x29(nx), x39(nx), x49(nx)
  real(r_kind) :: tend1(nx), tend2(nx), tend3(nx), tend4(nx)
  real(r_kind) :: tend19(nx), tend29(nx), tend39(nx), tend49(nx)
  integer :: k

  cp_xin(:) = xin(:)
  cp_xin9(:) = xin9(:)
  do k = 1, nt
    ! RK4 1st stage 
      call tl_tendency( cp_xin9, cp_xin, tend19, tend1, para )
      x2(:) = cp_xin(:) + 0.5d0 * dt * tend1(:)
      x29(:) = cp_xin9(:) + 0.5d0 * dt * tend19(:)

    ! RK4 2nd stage
      call tl_tendency( x29, x2, tend29, tend2, para )
      x3(:) = cp_xin(:) + 0.5d0 * dt * tend2(:)
      x39(:) = cp_xin9(:) + 0.5d0 * dt * tend29(:)

    ! RK4 3rd stage
      call tl_tendency( x39, x3, tend39, tend3, para )
      x4(:) = cp_xin(:) +         dt * tend3(:)
      x49(:) = cp_xin9(:) +         dt * tend39(:)

    ! RK4 4th stage
      call tl_tendency( x49, x4, tend49, tend4, para )
      xout(:) = cp_xin(:) + &
            dt * ( tend1(:) + 2.0d0*tend2(:) + 2.0d0*tend3(:) + tend4(:) )/6.0d0
      xout9(:) = cp_xin9(:) + &
            dt * ( tend19(:) + 2.0d0*tend29(:) + 2.0d0*tend39(:) + tend49(:) )/6.0d0
      t = t+dt

      if (k<nt) then
          cp_xin(:) = xout(:)
          cp_xin9(:) = xout9(:)
      endif

  enddo

endsubroutine


! multi-step RK-4 ADJ
!subroutine adj_lorenz63_rk4( xin9,  xin, para, nt, dt, t, xout9, xout )
subroutine adj_lorenz63_rk4( xin9,  xin, xout9, xout, para, nt, dt )
  implicit none
! passed args
  real(r_kind),intent(in   ) :: xin9(nx)
  real(r_kind),intent(inout) :: xin(nx)
  real(r_kind),intent(in   ) :: para(nx)
  integer,     intent(in   ) :: nt
  real(r_kind),intent(in   ) :: dt
  !real(r_kind),intent(inout) :: t
  real(r_kind),intent(inout) :: xout9(nx)
  real(r_kind),intent(inout) :: xout(nx)
! local vars
  real(r_kind) :: t
  real(r_kind) :: cp_xin(nx)
  real(r_kind) :: cp_xin9(nx,nt)
  real(r_kind) :: x2(nx), x3(nx), x4(nx)
  real(r_kind) :: x29(nx,nt), x39(nx,nt), x49(nx,nt)
  real(r_kind) :: tend1(nx), tend2(nx), tend3(nx), tend4(nx)
  real(r_kind) :: tend19(nx,nt), tend29(nx,nt), tend39(nx,nt), tend49(nx,nt)
  integer :: k

! base trajectory calculation
  cp_xin9(:,1) = xin9(:)
  do k = 1, nt
      Call tendency( cp_xin9(:,k), tend19(:,k), para )
      x29(:,k) = cp_xin9(:,k) + 0.5d0 * dt * tend19(:,k)
      Call tendency( x29(:,k), tend29(:,k), para )
      x39(:,k) = cp_xin9(:,k) + 0.5d0 * dt * tend29(:,k)
      Call tendency( x39(:,k), tend39(:,k), para )
      x49(:,k) = cp_xin9(:,k) +         dt * tend39(:,k)
      Call tendency( x49(:,k), tend49(:,k), para )
      xout9(:) = cp_xin9(:,k) + &
            dt * ( tend19(:,k) + 2.0d0*tend29(:,k) + 2.0d0*tend39(:,k) + tend49(:,k) )/6.0d0
      if (k<nt) cp_xin9(:,k+1) = xout9(:)
  enddo

  do k = nt, 1, -1
      t = t - dt
      xin = 0.d0; x2 = 0.d0; x3 = 0.d0; x4 = 0.d0
      tend1 = 0.d0; tend2=0.d0; tend3=0.d0; tend4=0.d0
    ! ADJ: RK4 4th stage
      !xout(:) = xin(:) + &
      !      dt * ( tend1(:) + 2.0d0*tend2(:) + 2.0d0*tend3(:) + tend4(:) )/6.0d0
      xin(:) = xout(:) 
      tend1(:) = dt*xout(:)/6.d0
      tend2(:) = dt*2.0d0*xout(:)/6.d0
      tend3(:) = dt*2.0d0*xout(:)/6.d0
      tend4(:) = dt*xout(:)/6.d0
      !call tl_tendency( x49, x4, tend49, tend4, para )
      call adj_tendency( x49(:,k), x4, tend49(:,k), tend4, para )

    ! ADJ: RK4 3rd stage
      !x4(:) = xin(:) +         dt * tend3(:)
      xin(:) = x4(:) + xin(:)
      tend3(:) = dt * x4(:) + tend3(:)
      !call tl_tendency( x39, x3, tend39, tend3, para )
      call adj_tendency( x39(:,k), x3, tend39(:,k), tend3, para )

    ! ADJ: RK4 2nd stage
      !x3(:) = xin(:) + 0.5d0 * dt * tend2(:)
      xin(:) = x3(:) + xin(:)
      tend2(:) =  0.5d0 * dt * x3(:) + tend2(:)
      !call tl_tendency( x29, x2, tend29, tend2, para )
      call adj_tendency( x29(:,k), x2, tend29(:,k), tend2, para )

    ! ADJ: RK4 1st stage 
      !x2(:) = xin(:) + 0.5d0 * dt * tend1(:)
      xin(:) = x2(:) + xin(:)
      tend1(:) = 0.5d0 * dt * x2(:) + tend1(:)
      !call tl_tendency( xin9, xin, tend19, tend1, para )
      call adj_tendency( cp_xin9(:,k), xin, tend19(:,k), tend1, para )

    ! go back to the previous step
      if (k>1) xout(:) = xin(:)
  enddo

endsubroutine


! RK-4 FWD
subroutine step_rk4( xin,  para, dt, t, xout )
  implicit none
! passed args
  real(r_kind),intent(in   ) :: xin(nx)
  real(r_kind),intent(in   ) :: para(nx)
  real(r_kind),intent(in   ) :: dt
  real(r_kind),intent(inout) :: t
  real(r_kind),intent(  out) :: xout(nx)
! local vars
  real(r_kind) :: x2(nx), x3(nx), x4(nx)
  real(r_kind) :: tend1(nx), tend2(nx), tend3(nx), tend4(nx)

! RK4 1st stage 
  Call tendency( xin, tend1, para )
  x2(:) = xin(:) + 0.5d0 * dt * tend1(:)
! RK4 2nd stage
  Call tendency( x2, tend2, para )
  x3(:) = xin(:) + 0.5d0 * dt * tend2(:)
! RK4 3rd stage
  Call tendency( x3, tend3, para )
  x4(:) = xin(:) +         dt * tend3(:)
! RK4 4th stage
  Call tendency( x4, tend4, para )
  xout(:) = xin(:) + &
        dt * ( tend1(:) + 2.0d0*tend2(:) + 2.0d0*tend3(:) + tend4(:) )/6.0d0
  t = t+dt

endsubroutine


! RK-4 TL
subroutine tl_step_rk4( xin9,  xin, para, dt, t, xout9, xout )
  implicit none
! passed args
  real(r_kind),intent(in   ) :: xin9(nx)
  real(r_kind),intent(in   ) :: xin(nx)
  real(r_kind),intent(in   ) :: para(nx)
  real(r_kind),intent(in   ) :: dt
  real(r_kind),intent(inout) :: t
  real(r_kind),intent(  out) :: xout9(nx)
  real(r_kind),intent(  out) :: xout(nx)
! local vars
  real(r_kind) :: x2(nx), x3(nx), x4(nx)
  real(r_kind) :: x29(nx), x39(nx), x49(nx)
  real(r_kind) :: tend1(nx), tend2(nx), tend3(nx), tend4(nx)
  real(r_kind) :: tend19(nx), tend29(nx), tend39(nx), tend49(nx)

! RK4 1st stage 
  call tl_tendency( xin9, xin, tend19, tend1, para )
  x2(:) = xin(:) + 0.5d0 * dt * tend1(:)
  x29(:) = xin9(:) + 0.5d0 * dt * tend19(:)

! RK4 2nd stage
  call tl_tendency( x29, x2, tend29, tend2, para )
  x3(:) = xin(:) + 0.5d0 * dt * tend2(:)
  x39(:) = xin9(:) + 0.5d0 * dt * tend29(:)

! RK4 3rd stage
  call tl_tendency( x39, x3, tend39, tend3, para )
  x4(:) = xin(:) +         dt * tend3(:)
  x49(:) = xin9(:) +         dt * tend39(:)

! RK4 4th stage
  call tl_tendency( x49, x4, tend49, tend4, para )
  xout(:) = xin(:) + &
        dt * ( tend1(:) + 2.0d0*tend2(:) + 2.0d0*tend3(:) + tend4(:) )/6.0d0
  xout9(:) = xin9(:) + &
        dt * ( tend19(:) + 2.0d0*tend29(:) + 2.0d0*tend39(:) + tend49(:) )/6.0d0
  t = t+dt

endsubroutine


! RK-4 ADJ
subroutine adj_step_rk4( xin9,  xin, para, dt, t, xout9, xout )
  implicit none
! passed args
  real(r_kind),intent(in   ) :: xin9(nx)
  real(r_kind),intent(inout) :: xin(nx)
  real(r_kind),intent(in   ) :: para(nx)
  real(r_kind),intent(in   ) :: dt
  real(r_kind),intent(inout) :: t
  real(r_kind),intent(inout) :: xout9(nx)
  real(r_kind),intent(inout) :: xout(nx)
! local vars
  real(r_kind) :: x2(nx), x3(nx), x4(nx)
  real(r_kind) :: x29(nx), x39(nx), x49(nx)
  real(r_kind) :: tend1(nx), tend2(nx), tend3(nx), tend4(nx)
  real(r_kind) :: tend19(nx), tend29(nx), tend39(nx), tend49(nx)

! base trajectory calculation
  t = t - dt
  Call tendency( xin9, tend19, para )
  x29(:) = xin9(:) + 0.5d0 * dt * tend19(:)
  Call tendency( x29, tend29, para )
  x39(:) = xin9(:) + 0.5d0 * dt * tend29(:)
  Call tendency( x39, tend39, para )
  x49(:) = xin9(:) +         dt * tend39(:)
  Call tendency( x49, tend49, para )
  xout9(:) = xin9(:) + &
        dt * ( tend19(:) + 2.0d0*tend29(:) + 2.0d0*tend39(:) + tend49(:) )/6.0d0

  xin = 0.d0; x2 = 0.d0; x3 = 0.d0; x4 = 0.d0
  tend1 = 0.d0; tend2=0.d0; tend3=0.d0; tend4=0.d0
! ADJ: RK4 4th stage
  !xout(:) = xin(:) + &
  !      dt * ( tend1(:) + 2.0d0*tend2(:) + 2.0d0*tend3(:) + tend4(:) )/6.0d0
  xin(:) = xout(:) 
  tend1(:) = dt*xout(:)/6.d0
  tend2(:) = dt*2.0d0*xout(:)/6.d0
  tend3(:) = dt*2.0d0*xout(:)/6.d0
  tend4(:) = dt*xout(:)/6.d0
  !call tl_tendency( x49, x4, tend49, tend4, para )
  call adj_tendency( x49, x4, tend49, tend4, para )

! ADJ: RK4 3rd stage
  !x4(:) = xin(:) +         dt * tend3(:)
  xin(:) = x4(:) + xin(:)
  tend3(:) = dt * x4(:) + tend3(:)
  !call tl_tendency( x39, x3, tend39, tend3, para )
  call adj_tendency( x39, x3, tend39, tend3, para )

! ADJ: RK4 2nd stage
  !x3(:) = xin(:) + 0.5d0 * dt * tend2(:)
  xin(:) = x3(:) + xin(:)
  tend2(:) =  0.5d0 * dt * x3(:) + tend2(:)
  !call tl_tendency( x29, x2, tend29, tend2, para )
  call adj_tendency( x29, x2, tend29, tend2, para )

! ADJ: RK4 1st stage 
  !x2(:) = xin(:) + 0.5d0 * dt * tend1(:)
  xin(:) = x2(:) + xin(:)
  tend1(:) = 0.5d0 * dt * x2(:) + tend1(:)
  !call tl_tendency( xin9, xin, tend19, tend1, para )
  call adj_tendency( xin9, xin, tend19, tend1, para )

endsubroutine

!-----------------------------------------------------------------------------------------
!  FWD, TL, ADJ of lorenz63 with forward in time method
!-----------------------------------------------------------------------------------------
subroutine step_ft( xin,  para, dt, t, xout )
  implicit none
! passed args
  real(r_kind),intent(in   ) :: xin(nx)
  real(r_kind),intent(in   ) :: para(nx)
  real(r_kind),intent(in   ) :: dt
  real(r_kind),intent(inout) :: t
  real(r_kind),intent(  out) :: xout(nx)
! local vars
  real(r_kind) :: tend1(nx)

  Call tendency( xin, tend1, para )
  xout(:) = xin(:) + dt * tend1(:)
  t = t+dt

endsubroutine 


subroutine tl_step_ft( xin9, xin, para, dt, t, xout9, xout )
  implicit none
! passed args
  real(r_kind),intent(in   ) :: xin(nx)
  real(r_kind),intent(in   ) :: xin9(nx)
  real(r_kind),intent(in   ) :: para(nx)
  real(r_kind),intent(in   ) :: dt
  real(r_kind),intent(inout) :: t
  real(r_kind),intent(  out) :: xout(nx)
  real(r_kind),intent(  out) :: xout9(nx)
! local vars
  real(r_kind) :: tend1(nx)
  real(r_kind) :: tend19(nx)

  call tl_tendency( xin9, xin, tend19, tend1, para )
  xout(:)  = xin(:)  + dt * tend1(:)
  xout9(:) = xin9(:) + dt * tend19(:)
  t = t+dt

endsubroutine 


subroutine adj_step_ft( xin9, xin,  para, dt, t, xout9, xout )
  implicit none
! passed args
  real(r_kind),intent(  out) :: xin(nx)
  real(r_kind),intent(in   ) :: xin9(nx)
  real(r_kind),intent(in   ) :: para(nx)
  real(r_kind),intent(in   ) :: dt
  real(r_kind),intent(inout) :: t
  real(r_kind),intent(inout) :: xout(nx)
  real(r_kind),intent(inout) :: xout9(nx)
! local vars
  real(r_kind) :: tend1(nx)
  real(r_kind) :: tend19(nx)

  ! base trajectory
  call step_ft( xin9,  para, dt, t, xout9 )

  t = t-dt
  !TL: xout(:)  = xin(:)  + dt * tend1(:)
  xin(:) = xout(:)
  tend1(:) = dt*xout(:)
  !TL: call tl_tendency( xin9, xin, tend1, tend19, tend1, para )
  call adj_tendency( xin9, xin, tend19, tend1, para )
  
endsubroutine 

!-----------------------------------------------------------------------------------------
!  FWD, TL, ADJ of tendency term
!-----------------------------------------------------------------------------------------
subroutine tendency( xin, tend, para )
  implicit none

  real(r_kind),intent(in   ) :: xin(nx)
  real(r_kind),intent(in   ) :: para(nx)
  real(r_kind),intent(  out) :: tend(nx)

  tend(1) = -para(1)*( xin(1) - xin(2) )
  tend(2) = para(2)*xin(1) - xin(2) - xin(1)*xin(3)
  tend(3) = xin(1)*xin(2) - para(3)*xin(3)

endsubroutine


subroutine tl_tendency( xin9, xin, tend9, tend, para )
  implicit none

  real(r_kind),intent(in   ) :: xin(nx)
  real(r_kind),intent(in   ) :: xin9(nx)
  real(r_kind),intent(in   ) :: para(nx)
  real(r_kind),intent(  out) :: tend(nx)
  real(r_kind),intent(  out) :: tend9(nx)

  ! first copy the original code
  ! then based on the original code, insert the trajectory var directly
  ! modify fwd codes by adding 9

  tend(1) = -para(1)*( xin(1) - xin(2) )
  tend9(1) = -para(1)*( xin9(1) - xin9(2) )

  tend(2) = para(2)*xin(1) - xin(2) - xin(1)*xin9(3) - xin9(1)*xin(3)
  tend9(2) = para(2)*xin9(1) - xin9(2) - xin9(1)*xin9(3)

  tend(3) = xin9(1)*xin(2) + xin(1)*xin9(2) - para(3)*xin(3)
  tend9(3) = xin9(1)*xin9(2) - para(3)*xin9(3)

endsubroutine


subroutine adj_tendency( xin9, xin, tend9, tend, para )
  implicit none

  real(r_kind),intent(inout) :: xin(nx)
  real(r_kind),intent(in   ) :: xin9(nx)
  real(r_kind),intent(in   ) :: para(nx)
  real(r_kind),intent(in   ) :: tend(nx)
  real(r_kind),intent(inout) :: tend9(nx)

  !tend(1) = -para(1)*( xin(1) - xin(2) )
  !tend(2) = para(2)*xin(1) - xin(2) - xin(1)*xin9(3) - xin9(1)*xin(3)
  !tend(3) = xin9(1)*xin(2) + xin(1)*xin9(2) - para(3)*xin(3)

  !tend(3) = xin9(1)*xin(2) + xin(1)*xin9(2) - para(3)*xin(3)
  xin(1) = tend(3)*xin9(2)  + xin(1) 
  xin(2) = xin9(1)*tend(3)  + xin(2)
  xin(3) = -para(3)*tend(3) + xin(3)

  !tend(2) = para(2)*xin(1) - xin(2) - xin(1)*xin9(3) - xin9(1)*xin(3)
  xin(1) = para(2)*tend(2)  + xin(1)
  xin(2) = -tend(2)         + xin(2)
  xin(1) = -tend(2)*xin9(3) + xin(1)
  xin(3) = -xin9(1)*tend(2) + xin(3)

  !tend(1) = -para(1)*( xin(1) - xin(2) ) = -para(1)*xin(1)+para(1)*xin(2)
  xin(1) = -para(1)*tend(1) + xin(1)
  xin(2) =  para(1)*tend(1) + xin(2)

endsubroutine


endmodule 

