module mod_lorenz63_inc4dvar
!$$$  program documentation block
!         .           .            .
!  program name: mod_lorenz63_3dvar
!    programmer: da,cheng        org: umd      date: 2015-Jan-11
!
!  purpose:
!    incremental 4D-Var system with multi-outerloop
!    CVT is used
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
  !use mod_math, only : inv=>invmtx, sqrtm=>sqrtmtx
  use mod_math, only : inv, sqrtm=>sqrtmtx
  use mod_optimization, only : run_lbfgs, destroy_lbfgs
  implicit none

  private 
  public :: lorenz63_inc4dvar

  !real(rdp),parameter :: epsl   = 1.d-7
  real(rdp),parameter :: epsl   = 1.d-5
contains

!
! Incremental 4D-Var with control variable transform
!
subroutine lorenz63_inc4dvar( nx, nyo, nda, xb, B, lyo, yo, erro, xa, kitermax_in, kitermax_out, dt, ndt, para)
  use mod_lorenz63_fwd, only : lorenz63_rk4, tl_lorenz63_rk4, adj_lorenz63_rk4 
  implicit none

! passed args
  integer,intent(in   ) :: nx
  integer,intent(in   ) :: nyo
  integer,intent(in   ) :: nda
  real(8),intent(in   ) :: xb(nx)
  real(8),intent(in   ) :: B(nx,nx)
  logical,intent(in   ) :: lyo(nyo,nda)
  real(8),intent(in   ) :: yo(nyo,nda)
  real(8),intent(in   ) :: erro(nyo,nda)
  real(8),intent(  out) :: xa(nx)
  integer,intent(in   ) :: kitermax_in
  integer,intent(in   ) :: kitermax_out
  real(8),intent(in   ) :: dt
  integer,intent(in   ) :: ndt
  real(8),intent(in   ) :: para(nx)
! local vars
  real(8) :: invB(nx,nx)
  real(8) :: invR(nyo,nda)
  real(8) :: Jc, Jo, Jb
  real(8) :: dJc(nx), dJo(nx), dJb(nx)

  real(8) :: ygt(nyo)       ! H(xg)
  real(8) :: dyg(nyo,nda)   ! yobs - h(xg)
  real(8) :: dy(nyo,nda)    ! yobs- h(xg)-Hdx
  real(8) :: xg(nx)         ! reference state at first step
  real(8) :: xgt0(nx,nda) ! reference state at all time slots
  real(8) :: dx(nx)         ! x-xg
  real(8) :: dxt0(nx,nda) 
  real(8) :: tmpx(nx)

  real(8) :: forcing0(nx)
  real(8) :: tmpforcing0(nx)

  real(8) :: dub(nx)
  real(8) :: du(nx), dv(nx)
  real(8) :: L(nx,nx), invL(nx,nx)

  real(8) :: gnorm, xnorm

  real(8) :: t
  integer :: ierr
  integer :: i, j, m, n
  logical :: lfound


  call sqrtm(nx,B,L)
  !call inv(nx,L,invL,ierr)
  call inv(nx,L,invL)

  do j = 1, nda
     do i = 1, nyo
        invR(i,j) = 1.d0/( erro(i,j)*erro(i,j) )
     enddo
  enddo

! we transform with Ldu = dx
  xg = xb
  outer_loop: do m = 1, kitermax_out

     print*, "outerloop=", m

     ! advance NWP & calculate nonlinear H(x)
     xgt0(1:nx,1) = xg
     do j = 1, nda

        if ( j>1 ) then
           !call step( xgt0(:,j-1), t, dt, xgt0(:,j) )
           call lorenz63_rk4( xgt0(:,j-1), xgt0(:,j), para, ndt, dt )
        endif
        
        do i = 1, nyo
           ! <FIXME>: y=h(xg)
           ygt(i) = xgt0(i,j)  
           dyg(i,j) = yo(i,j) - ygt(i)   ! (yo_i - h_i(M(xg)))
        enddo

     enddo

     ! calculate dub = L^-1 * dxb = L^-1 * (xb - xg )
     dub = matmul(invL,xb-xg)


     ! innter loop to get increment on the reference state
     lfound = .false.
     du = 0.d0
     inner_loop: do n = 1, kitermax_in
     
        ! calculate Jb = 0.5* (du-dub)^t * (du-dub)
        !           du = L^-1 * dx
        Jb = 0.d0
        do i = 1, nx
           Jb = Jb + 0.5d0*(du(i)-dub(i))**2
        enddo

        ! calculate Jo  =  sum { 0.5 * ( yo_i-h_i(x_i) )          * R^-1 * ( yo_i-h_i(x_i)          )  } 
        !              ~=  sum { 0.5 * ( yo_i-h_i(N(xg))- HMdx )  * R^-1 * ( yo_i-h_i(N(xg))- HMdx  )  }
        !              ~=  sum { 0.5 * ( yo_i-h_i(N(xg))- HMLdu ) * R^-1 * ( yo_i-h_i(N(xg))- HMLdu )  } used here
        dx  = matmul(L,du)
        Jo  = 0.d0
        dxt0(1:nx,1) = dx
        do j = 1, nda
           ! advance perturb with TLM
           if (j>1) then
              !call tl_step(dxt0(:,j-1),xgt0(:,j-1),t,dt,dxt0(:,j) )
              call tl_lorenz63_rk4( xin9  = xgt0(:,j-1),  & 
                                    xin   = dxt0(:,j-1),  &
                                    xout9 = tmpx(:),  &
                                    xout  = dxt0(:,j),  &
                                    para  = para,  &
                                    nt    = ndt,  &
                                    dt    = dt )
           endif
           ! HMdx 
           do i = 1, nyo
              if ( lyo(i,j) ) then
              ! ----------
              ! add Hd(x) here
              ! ----------
                 dy(i,j) = dyg(i,j) - dxt0(i,j)  ! modify here
                 Jo = Jo + 0.5d0*dy(i,j)*dy(i,j)*invR(i,j) 
              endif
           enddo
        enddo

        ! calculate total cost function Jc = Jb + Jo
        Jc = Jb + Jo

        ! calculate gradient dJb = I^-1 * (du-dub)
        dJb = 0.0d0
        do i = 1, nx
           dJb(i) = du(i) - dub(i)
        enddo

        ! calculate gradient dJo  = -L^t * sum { M^t * H^t * R^-1 * ( yo - h(xg+dx)     ) }
        !                        ~= -L^t * sum { M^t * H^t * R^-1 * ( yo - h(xg) - Hdx  ) }
        !                        ~= -L^t * sum { M^t * H^t * R^-1 * ( yo - h(xg) - HLdu ) } used here
        dJo = 0.d0
        dv  = 0.d0   ! dv = sum { M^t * H^t * R^-1 * ( yo - h(xg) - HLdu ) }
                     ! dJo = -L^t * dv
!        do j = nda, 1, -1
!
!           forcing0 = 0.d0
!           do i = 1, nyo
!              if ( lyo(i,j) ) then
!                 forcing0(i) = dy(i,j)*invR(i,j)
!              endif
!           enddo
!           if (j>1) then
!              do i = j, 2, -1
!                 call ad_step(forcing0, xgt0(:,i-1),t,dt,forcing0)
!              enddo
!           endif
!           dv = dv + forcing0(1:nx)
!
!        enddo
        
        ! add obs misfit advected from previous step with adjoint model.
        ! this approach saves computational cost
        forcing0 = 0.d0
        do j = nda, 1, -1

           do i = 1, nyo
              if ( lyo(i,j) ) then
                 forcing0(i) = forcing0(i)+dy(i,j)*invR(i,j)
              endif
           enddo
           if (j>1) then
              !call ad_step(forcing0, xgt0(:,j-1),t,dt,forcing0)
              call adj_lorenz63_rk4( xin9  = xgt0(:,j-1), &
                                     xin   = tmpforcing0, &
                                     xout9 = tmpx(:), &
                                     xout  = forcing0, &
                                     para  = para, &
                                     nt    = ndt, &
                                     dt    = dt)
              forcing0 = tmpforcing0
           endif

        enddo
        dv = dv + forcing0(1:nx)


        do i = 1, nx
           do j = 1, nx
              dJo(i) = dJo(i) + L(j,i)*dv(j)   ! L(j,i)*dv(j) = L^t(i,j)*dv(j)
           enddo
        enddo
        dJo(:) = -dJo(:)

        ! calculate total gradient
        dJc = dJo + dJb

        call run_lbfgs( nx, du, Jc, dJc, epsl, ierr )
        !gnorm = sqrt(dot_product(dJc,dJc))
        !xnorm = max(sqrt(dot_product(du,du)),1.d0)
        !if (gnorm/xnorm<epsl) then
        !   print*, "[msg] gnorm/xnorm<epsl, skip this outerloop: g/x,epsl=", gnorm, 1.d-10
        !   ierr = 0
        !endif
        if ( ierr < 0 ) Stop "[err] inc4dvar: fail to run lbfgs."
        if ( ierr == 0 ) then
           call destroy_lbfgs()
           lfound = .true.
           exit inner_loop
        endif

     enddo inner_loop

     if (.not.lfound) then
        print*, "[warning] inc4dvar: reach the maximum iteration=", kitermax_in
        call destroy_lbfgs()
     endif
     xg = xg + matmul(L,du)
 
  enddo outer_loop

  xa = xg 

  print*, "finished"

endsubroutine 

endmodule


