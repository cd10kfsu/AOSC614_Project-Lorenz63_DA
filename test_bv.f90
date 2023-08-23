program test_bv
!$$$  program documentation block
!         .           .            .
!  program name: test_bv
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
  use mod_lorenz63_fwd,   only : lorenz63_rk4, t_lorenz63, step, get_tlm
  use mod_rnorm,          only : rnorm, set_random_seed
  implicit none

  integer,parameter :: nx = 3
  integer          :: ncycle, nspin
  type(t_lorenz63) :: xb, x, x2, xInfy, xInfyRun2
  real(rdef) :: buft
  integer :: i, j, itry
  real(rdef) :: dx(nx), dx2(nx), dxInfy(nx), dxInfyRun2(nx)
  real(rdef) :: dx0(nx), dx0Run2(nx)
  real(rdef) :: le1, le2, leInfy
  real(rdef) :: sum1, sum2, sumInfy
  real(rdef) :: fracForLLECalc
  integer :: nsum
  integer :: ntry

  real(rdef),allocatable :: lle_table(:,:) ! ntry*3 |lle1|lle2|lleInfy|



!----------------------------------------------------------------------------------
! 0. configuration for lorenz63 model
!----------------------------------------------------------------------------------
  !xb%dt   = 0.008d0   ! length for one time step
  xb%dt   = 0.001d0   ! length for one time step
  xb%nt   = 1       ! generate nt-step forecast
  xb%para = (/ 10.0d0, 28.0d0, 8.0d0/3.0d0 /)  

  nspin   = 200000     ! spin up the model with nspin steps, then the final state 
                    ! is used as the initial condition
  ncycle  = 80000000     ! total n-step forecast cycles

  fracForLLECalc = 0.18 ! start to calculate LLE after cylce > ncyclc&fracForLLECalc

  ntry = 10
  allocate(lle_table(ntry,4))


  do itry = 1, ntry
     call set_random_seed( iseed=itry+62 )

      do i = 1, nx
         dx0(i) = rnorm()
         dx0Run2(i) = rnorm()
      enddo

      print*, "-----------------------------------------------"
      print*, "try #", itry, "random_seed=", itry+62
      print*, "START: random initial perturbations before rescaling, dx0, dx0Run2 = ", dx0, dx0Run2

     
    !----------------------------------------------------------------------------------
    ! 1. spin up the lorenz 63 model
    !----------------------------------------------------------------------------------
      !xb%x0   = [ 0.1d0,  0.1d0, 0.1d0 ]
      do i = 1, nx
         xb%x0(i) = rnorm()
      enddo
      Write(6,*) "spin up lorenz63 model: dt =", xb%dt, ", nspin =", nspin
      Call lorenz63_rk4( xb%x0(:), xb%xn(:), xb%para, nspin, xb%dt )
      xb%x0(:) = xb%xn(:)
      Write(6,*) ""
      Write(6,*) "INITIAL CONDITION: x=", xb%x0(:)


    !----------------------------------------------------------------------------------
    ! 2. rescale initial perturbations
    !----------------------------------------------------------------------------------
      print*, ""

      dx = dx0 ![ 0.6, 0.2, 0.5]
      !dx = dx/sqrt(sum(dx*dx))
      dx = dx/sum(abs(dx))
      print*, "l1-norm:",dx,sum(abs(dx))

      dx2 = dx0 ![ 0.05, 0.1, 0.8]
      dx2 = dx2/sqrt(sum(dx2*dx2))
      print*, "l2-norm:", dx2,sum(dx2*dx2)

      dxInfy = dx0 ![ 0.01, 2.0, 0.01]
      dxInfy = dxInfy/maxval(abs(dxInfy))
      print*, "lInfy-norm for member 1:", dxInfy,maxval(abs(dxInfy))

      dxInfyRun2 = dx0Run2 ![ 0.01, 2.0, 0.01]
      dxInfyRun2 = dxInfyRun2/maxval(abs(dxInfyRun2))
      print*, "lInfy-norm for member 2:", dxInfy,maxval(abs(dxInfyRun2))

      x%x0(:)     = xb%x0(:) + dx 
      x2%x0(:)    = xb%x0(:) + dx2 
      xInfy%x0(:) = xb%x0(:) + dxInfy 
      xInfyRun2%x0(:) = xb%x0(:) + dxInfyRun2

      sum1 = 0.d0; sum2 = 0.d0; sumInfy = 0.d0
      nsum = 0

    !----------------------------------------------------------------------------------
    ! 3. breeding
    !----------------------------------------------------------------------------------
      do i = 1, ncycle
         call step(xb%x0(:), xb%para, xb%dt, buft, xb%xn(:))
         call step(x%x0(:), xb%para, xb%dt, buft, x%xn(:))
         call step(x2%x0(:), xb%para, xb%dt, buft, x2%xn(:))
         call step(xInfy%x0(:), xb%para, xb%dt, buft, xInfy%xn(:))
         call step(xInfyRun2%x0(:), xb%para, xb%dt, buft, xInfyRun2%xn(:))

         ! BV now
         dx(:) = x%xn(:) - xb%xn(:)
         dx2(:) = x2%xn(:) - xb%xn(:)
         dxInfy(:) = xInfy%xn(:) - xb%xn(:)
         dxInfyRun2(:) = xInfyRun2%xn(:) - xb%xn(:)
         if (i > ncycle*fracForLLECalc) then
            le1    = log( sum(abs(dx)) / sum(abs(x%x0-xb%x0)) ) / xb%dt
            le2    = log( sqrt(sum(dx2*dx2)) / sqrt(sum((x2%x0-xb%x0)*(x2%x0-xb%x0))) )/xb%dt
            leInfy = log( maxval(abs(dxInfy)) / maxval(abs(xInfy%x0-xb%x0))    )/xb%dt

            sum1 = sum1 + le1; sum2 = sum2 + le2; sumInfy = sumInfy + leInfy
            nsum = nsum + 1
         end if

         if (i<=5 .or. i== int(ncycle*fracForLLECalc) .or. mod(i,int(ncycle/10.))==0) then  ! check local LLV
            print*, ""
            print*, "CHECK Local LLV at time step: ", i
            print*, "llv1=", dx
            print*, "llv2=", dx2
            print*, "llvInfy for member 1=", dxInfy
            print*, "llvInfy for member 2=", dxInfyRun2
         endif
  

         ! rescale BV 
         dx =  dx/sum(abs(dx))
         dx2 = dx2/sqrt(sum(dx2*dx2))
         dxInfy = dxInfy / maxval(abs(dxInfy))
         dxInfyRun2 = dxInfyRun2 / maxval(abs(dxInfyRun2))

         xb%x0(:) = xb%xn(:)
         x%x0(:)  = xb%x0(:) + dx(:)
         x2%x0(:)  = xb%x0(:) + dx2(:)
         xInfy%x0(:) = xb%x0(:) + dxInfy(:)
         xInfyRun2%x0(:) = xb%x0(:) + dxInfyRun2(:)
      enddo

      print*, ""
      print*, "SECTION_RESULTS: lle1, lle2, lleInfy, n=", sum1/nsum, sum2/nsum, sumInfy/nsum, nsum

      lle_table(itry,1) = sum1/nsum
      lle_table(itry,2) = sum2/nsum
      lle_table(itry,3) = sumInfy/nsum

  enddo

  print*, "RESULTS SUMMARY:"
  print*, "| try# | lle1 | lle2 | lleInfy|"
  print*, "-------------------------------"
  do itry = 1, ntry
     print*, itry, (lle_table(itry,j), j=1,3)
  enddo
  print*, "-------------------------------"
  print*, "avg", (sum(lle_table(:,j))/ntry, j=1,3)





endprogram

