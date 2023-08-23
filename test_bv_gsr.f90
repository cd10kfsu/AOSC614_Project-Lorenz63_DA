program test_bv_gsr

  use mod_type,           only : rdef
  use mod_lorenz63_fwd,   only : lorenz63_rk4, t_lorenz63, step, get_tlm
  use mod_rnorm,          only : rnorm, set_random_seed
  implicit none

  integer,parameter :: nx = 3
  integer :: nbv
  integer          :: ncycle, nspin
  type(t_lorenz63),allocatable :: x(:), x2(:), xInfy(:)  ! nbv
  type(t_lorenz63) :: xb
  real(rdef) :: buft
  integer :: i, j, itry, ibv, inorm
  !real(rdef) :: dx(nx), dx2(nx), dxInfy(nx)
  !real(rdef) :: dx0(nx)
  real(rdef),allocatable :: dx(:,:), dx2(:,:), dxInfy(:,:), dx0(:,:) ! nx*nbv
  real(rdef),allocatable :: bufdx(:,:), buf2dx(:,:) ! nx*nbv
  real(rdef) :: le1, le2, leInfy
  real(rdef),allocatable :: sum1(:), sum2(:), sumInfy(:) ! nbv
  real(rdef) :: fracForLLECalc
  integer :: nsum
  integer :: ntry

  real(rdef),allocatable :: lle_table(:,:,:) ! nbv*ntry*3 |lle1|lle2|lleInfy|



!----------------------------------------------------------------------------------
! 0. configuration for lorenz63 model
!----------------------------------------------------------------------------------
  !xb%dt   = 0.008d0   ! length for one time step
  xb%dt   = 0.001d0   ! length for one time step
  xb%nt   = 1       ! generate nt-step forecast
  xb%para = (/ 10.0d0, 28.0d0, 8.0d0/3.0d0 /)  

  nspin   = 200000     ! spin up the model with nspin steps, then the final state 
                    ! is used as the initial condition
  ncycle  = 20000000     ! total n-step forecast cycles

  fracForLLECalc = 0.18 ! start to calculate LLE after cylce > ncyclc&fracForLLECalc

  ntry = 10
  nbv = 3

  allocate(lle_table(nbv,ntry,3))

  allocate(x(nbv),x2(nbv),xInfy(nbv))
  allocate(dx(nx,nbv), dx2(nx,nbv), dxInfy(nx,nbv),dx0(nx,nbv))
  allocate(bufdx(nx,nbv),buf2dx(nx,nbv))
  allocate(sum1(nbv), sum2(nbv), sumInfy(nbv))


  do itry = 1, ntry
     call set_random_seed( iseed=itry+62 )

      do ibv = 1, nbv; do i = 1, nx
         dx0(i,ibv) = rnorm()
      enddo; enddo

      print*, "-----------------------------------------------"
      print*, "try #", itry, "random_seed=", itry+62
      print*, "START: random initial perturbations before rescaling, dx0 = ", dx0

     
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

      do ibv = 1, nbv
          dx(:,ibv) = dx0(:,ibv)
          !dx = dx/sqrt(sum(dx*dx))
          dx(:,ibv) = dx(:,ibv)/sum(abs(dx(:,ibv)))

          dx2(:,ibv) = dx0(:,ibv)
          dx2(:,ibv) = dx2(:,ibv)/sqrt(sum(dx2(:,ibv)*dx2(:,ibv)))

          dxInfy(:,ibv) = dx0(:,ibv) ![ 0.01, 2.0, 0.01]
          dxInfy(:,ibv) = dxInfy(:,ibv)/maxval(abs(dxInfy(:,ibv)))

          x(ibv)%x0(:)     = xb%x0(:) + dx(:,ibv) 
          x2(ibv)%x0(:)    = xb%x0(:) + dx2(:,ibv) 
          xInfy(ibv)%x0(:) = xb%x0(:) + dxInfy(:,ibv) 

      enddo

      sum1(:) = 0.d0; sum2(:) = 0.d0; sumInfy(:) = 0.d0
      nsum = 0

    !----------------------------------------------------------------------------------
    ! 3. breeding
    !----------------------------------------------------------------------------------
      do i = 1, ncycle

         do ibv = 1, nbv
             call step(xb%x0(:), xb%para, xb%dt, buft, xb%xn(:))
             call step(x(ibv)%x0(:), xb%para, xb%dt, buft, x(ibv)%xn(:))
             call step(x2(ibv)%x0(:), xb%para, xb%dt, buft, x2(ibv)%xn(:))
             call step(xInfy(ibv)%x0(:), xb%para, xb%dt, buft, xInfy(ibv)%xn(:))

             ! BV now
             dx(:,ibv) = x(ibv)%xn(:) - xb%xn(:)
             dx2(:,ibv) = x2(ibv)%xn(:) - xb%xn(:)
             dxInfy(:,ibv) = xInfy(ibv)%xn(:) - xb%xn(:)
         enddo

         ! hard-coded GS
         bufdx(:,:) = dx(:,:); buf2dx(:,:) = dx(:,:)
         buf2dx(:,1) = bufdx(:,1) ! 1st Local LV
         buf2dx(:,2) = bufdx(:,2) - buf2dx(:,1) * sum(buf2dx(:,1)*bufdx(:,2))/sum(buf2dx(:,1)*buf2dx(:,1))
         buf2dx(:,3) = bufdx(:,3) - buf2dx(:,1) * sum(buf2dx(:,1)*bufdx(:,3))/sum(buf2dx(:,1)*buf2dx(:,1)) &
                              - buf2dx(:,2) * sum(buf2dx(:,2)*bufdx(:,3))/sum(buf2dx(:,2)*buf2dx(:,2))
         dx(:,:) = buf2dx(:,:)
         
         bufdx(:,:) = dx2(:,:); buf2dx(:,:) = dx2(:,:)
         buf2dx(:,1) = bufdx(:,1) ! 1st Local LV
         buf2dx(:,2) = bufdx(:,2) - buf2dx(:,1) * sum(buf2dx(:,1)*bufdx(:,2))/sum(buf2dx(:,1)*buf2dx(:,1))
         buf2dx(:,3) = bufdx(:,3) - buf2dx(:,1) * sum(buf2dx(:,1)*bufdx(:,3))/sum(buf2dx(:,1)*buf2dx(:,1)) &
                               - buf2dx(:,2) * sum(buf2dx(:,2)*bufdx(:,3))/sum(buf2dx(:,2)*buf2dx(:,2))
         dx2(:,:) = buf2dx(:,:)

         bufdx(:,:) = dxInfy(:,:); buf2dx(:,:) = dxInfy(:,:)
         buf2dx(:,1) = bufdx(:,1) ! 1st Local LV
         buf2dx(:,2) = bufdx(:,2) - buf2dx(:,1) * sum(buf2dx(:,1)*bufdx(:,2))/sum(buf2dx(:,1)*buf2dx(:,1))
         buf2dx(:,3) = bufdx(:,3) - buf2dx(:,1) * sum(buf2dx(:,1)*bufdx(:,3))/sum(buf2dx(:,1)*buf2dx(:,1)) &
                               - buf2dx(:,2) * sum(buf2dx(:,2)*bufdx(:,3))/sum(buf2dx(:,2)*buf2dx(:,2))
         dxInfy(:,:) = buf2dx(:,:)

        do ibv = 1, nbv
             if (i > ncycle*fracForLLECalc) then
                le1    = log( sum(abs(dx(:,ibv))) / sum(abs(x(ibv)%x0-xb%x0)) ) / xb%dt
                le2    = log( sqrt(sum(dx2(:,ibv)*dx2(:,ibv))) / sqrt(sum((x2(ibv)%x0-xb%x0)*(x2(ibv)%x0-xb%x0))) )/xb%dt
                leInfy = log( maxval(abs(dxInfy(:,ibv))) / maxval(abs(xInfy(ibv)%x0-xb%x0))    )/xb%dt

                sum1(ibv) = sum1(ibv) + le1; sum2(ibv) = sum2(ibv) + le2; sumInfy(ibv) = sumInfy(ibv) + leInfy
                if (ibv==1) nsum = nsum + 1
             end if

             !if (i<=5 .or. i== int(ncycle*fracForLLECalc) .or. mod(i,int(ncycle/10.))==0) then  ! check local LLV
             !   print*, ""
             !   print*, "CHECK Local LLV at time step: ", i
             !   print*, "llv1=", dx
             !   print*, "llv2=", dx2
             !   print*, "llvInfy for member 1=", dxInfy
             !endif
        enddo
      
        do ibv = 1, nbv
             ! rescale BV 
             dxInfy(:,ibv) = dxInfy(:,ibv) / maxval(abs(dxInfy(:,ibv)))
             dx(:,ibv) =  dx(:,ibv)/sum(abs(dx(:,ibv)))
             dx2(:,ibv) = dx2(:,ibv)/sqrt(sum(dx2(:,ibv)*dx2(:,ibv)))

             xb%x0(:)         = xb%xn(:)
             x(ibv)%x0(:)     = xb%x0(:) + dx(:,ibv)
             x2(ibv)%x0(:)    = xb%x0(:) + dx2(:,ibv)
             xInfy(ibv)%x0(:) = xb%x0(:) + dxInfy(:,ibv)


         enddo


      enddo


      lle_table(:,itry,1) = sum1(:)/nsum ! nbv*ntry*3
      lle_table(:,itry,2) = sum2(:)/nsum
      lle_table(:,itry,3) = sumInfy(:)/nsum

      print*, "SECTION_RESULTS: lle1,    n=", sum1/nsum
      print*, "SECTION_RESULTS: lle2,    n=", sum2/nsum
      print*, "SECTION_RESULTS: lleInfy, n=", sumInfy/nsum

  enddo

  do ibv = 1, nbv
      print*, "RESULTS SUMMARY for LLE #:",ibv
      print*, "| try# | lle1 | lle2 | lleInfy|"
      print*, "-------------------------------"
      do itry = 1, ntry
         print*, itry, (lle_table(ibv,itry,j), j=1,3)
      enddo
      print*, "-------------------------------"
      print*, "avg", (sum(lle_table(ibv,:,j))/ntry, j=1,3)
  enddo





endprogram

