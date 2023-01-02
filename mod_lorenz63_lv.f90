!===============================================================================
! program name: m_lv        
!
! programmer: da,cheng        org: umd (cda@umd.edu)     date: 2018-03-23
!
! Description:
!   a module related to Lyapunov vector calculations
!
!-------------------------------------------------------------------------------
module mod_lorenz63_lv
  use mod_math, only : eye, qr
  implicit none

  private

  integer,parameter :: r_kind = kind(0.d0)


  public :: le_blv, le_flv

contains

!
! Lyapunov exponents calculated from backward Lyapunov vectors
!
subroutine le_blv(nx,nt,dt,M,LE,Q0)
  implicit none

  integer,intent(in) :: nx
  integer,intent(in) :: nt
  real(r_kind),intent(in)  :: dt
  real(r_kind),intent(in)  :: M(nx,nx,nt)   ! M(:,:,1)
  real(r_kind),intent(out) :: LE(nx)
  real(r_kind),intent(in),optional :: Q0(nx,nx)

  real(r_kind) :: Q(nx,nx,nt), R(nx,nx,nt)
  integer :: n

  call eye(nx,Q(:,:,1))
  if(present(Q0)) Q(:,:,1) = Q0
  call fwd_qr(nx,nt,M,Q,R,Q0=Q(:,:,1))
  do n = 1, nx
     LE(n) = sum(log(abs( R(n,n,1:nt) )))
  enddo
  LE(:)=LE(:)/(nt*dt)

endsubroutine
!
! Lyapunov exponents calculated from forward Lyapunov vectors
!
subroutine le_flv(nx,nt,dt,Mt,LE,Qn)
  implicit none

  integer,intent(in) :: nx
  integer,intent(in) :: nt
  real(r_kind),intent(in)  :: dt
  real(r_kind),intent(in)  :: Mt(nx,nx,nt)   ! M(:,:,i)   for TLM from step i to step i+1
  real(r_kind),intent(out) :: LE(nx)
  real(r_kind),intent(in),optional :: Qn(nx,nx)

  real(r_kind) :: Q(nx,nx,nt)  ! Q(:,:,i) for step i
  real(r_kind) :: R(nx,nx,nt)  ! R(:,:,i) from step i to step i+1
  real(r_kind) :: bufQ(nx,nx)
  integer :: n

  call eye(nx,bufQ)
  if(present(Qn)) bufQ = Qn
  call bwd_qr(nx,nt,Mt,Q,R,Qn=bufQ)
  do n = 1, nx
     LE(n) = sum(log(abs( R(n,n,1:nt) )))
  enddo
  LE(:)=LE(:)/(nt*dt)

endsubroutine
!
!
!
subroutine fwd_qr(nx,nt,M,Q,R,Q0)
  implicit none

  integer,intent(in) :: nx
  integer,intent(in) :: nt
  real(r_kind),intent(in) :: M(nx,nx,nt)
  real(r_kind),intent(out) :: Q(nx,nx,nt)
  real(r_kind),intent(out) :: R(nx,nx,nt)
  real(r_kind),intent(in),optional :: Q0(nx,nx)

  real(r_kind) :: bufQ(nx,nx), bufM(nx,nx)
  integer :: n

  call eye(nx,Q(:,:,1))  ! identiy matrix
  if (present(Q0)) Q(:,:,1)=Q0
  
  do n = 1, nt
     bufM = matmul(M(:,:,n),Q(:,:,n))  ! evolve from step i -> step i+1
     call qr(nx,bufM,bufQ,R(:,:,n))    ! get Q at step i+1, & R from step i->step i+1
     if (n<nt) Q(:,:,n+1) = bufQ
  enddo

endsubroutine
!
!
!
subroutine bwd_qr(nx,nt,Mt,Q,R,Qn)
  implicit none

  integer,intent(in) :: nx
  integer,intent(in) :: nt
  real(r_kind),intent(in) :: Mt(nx,nx,nt)
  real(r_kind),intent(out) :: Q(nx,nx,nt)
  real(r_kind),intent(out) :: R(nx,nx,nt)
  real(r_kind),intent(in),optional :: Qn(nx,nx)

  real(r_kind) :: bufQ(nx,nx), bufMt(nx,nx)
  integer :: n

  call eye(nx,bufQ)  ! identiy matrix
  if (present(Qn)) bufQ=Qn
  
  do n = nt, 1, -1
     bufMt = matmul(Mt(:,:,n),bufQ)  ! evolve from step i+1 -> step i
     call qr(nx,bufMt,bufQ,R(:,:,n))       ! get Q at step i, & R from step i+1 ->step i
     Q(:,:,n) = bufQ
  enddo


endsubroutine


endmodule
