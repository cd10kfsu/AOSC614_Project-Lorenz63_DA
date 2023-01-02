Module mod_optimization
!$$$  program documentation block
!         .           .            .
!  program name: mod_optimization
!    programmer: da,cheng        org: umd      date: 2015-Dec-07
!
!  purpose:
!    minimization module 
!
!  revision history:
!    2015-Dec-07     da    - creator
!
!  file dependencies:
!
!  attributes: 
!    language: fortran 90
!    machine : 
!
!
!$$$ end documentation block
  Implicit none

  Private 
  Public :: run_lbfgs
  Public :: destroy_lbfgs


  integer,parameter :: r8 = 8

  Type :: tLBGFS
    logical :: linit   = .false.
    logical :: ldiagco = .false.
    integer :: ndim    = 0      ! dimensions of the control variables (n)
    integer :: ncorr   = 0      ! number of correction used in the bgfs update (m)
    integer :: iprint(2)    = 0  
    integer :: iflag        = 0
    real(r8) :: eps_machine = 0.d0 ! (xtol)
    real(r8),allocatable :: wk(:) ! workspace (w)
    real(r8),allocatable :: rdiag(:)  ! (diag)
  Endtype

  Type(tLBGFS), save :: config_lbfgs


Contains

Subroutine initialize_lbfgs( ndim )
!$$$ subprogram documentation block
!             .                   .                    . 
!  Program name:
!    programmer: Da,Cheng       org: umd            date: 2015-MON-DD
!
!  Purpose:
!    
!
!  Revision history:
!    2015-MON-DD     Da,Cheng     -Creator
!
!  Input args:
!
!  Output args:
!
!
!
!$$$
  Implicit None
! Passed args
  Integer,intent(in) :: ndim
! Local vars
  Integer :: nwk

! default settings
  config_lbfgs%ncorr     = 5  ! should be >0. 3~7 recommended
  config_lbfgs%ldiagco   = .FALSE.
  config_lbfgs%iprint(1) = 1
  config_lbfgs%iprint(2) = 0
  config_lbfgs%iflag     = 0

! set dims of control vars
  config_lbfgs%ndim      = ndim
! get machine precision 
  config_lbfgs%eps_machine = EPSILON(0.0d0)
! allocate arrays
  allocate( config_lbfgs%rdiag(config_lbfgs%ndim) )
  config_lbfgs%rdiag = 0.0d0
  nwk = config_lbfgs%ndim*(2*config_lbfgs%ncorr+1) + &
        2*config_lbfgs%ncorr
  allocate( config_lbfgs%wk(nwk) )

  config_lbfgs%linit = .TRUE.
  
Endsubroutine


Subroutine run_lbfgs( ndim, x, func, grad, eps_opt, istatus )
!$$$ subprogram documentation block
!             .                   .                    . 
!  Program name:
!    programmer: Da,Cheng       org: umd            date: 2015-MON-DD
!
!  Purpose:
!    
!
!  Revision history:
!    2015-MON-DD     Da,Cheng     -Creator
!
!  Input args:
!
!  Output args:
!
!
!
!$$$
  Implicit None
! Passed args
  Integer,intent(in) :: ndim
  Real(r8),intent(inout) :: x(ndim)
  Real(r8),intent(inout) :: func
  Real(r8),intent(inout) :: grad(ndim)
  Real(r8),intent(in) :: eps_opt
  Integer,intent(out) :: istatus

! vars required by subroutine LBFGS
  external LB2 
  integer :: MP, LP
  real(r8) :: GTOL, STPMIN, STPMAX
  COMMON /LB3/MP,LP,GTOL,STPMIN,STPMAX
  

! peform some check
  If ( .not.(config_lbfgs%linit) ) Then
     Write(6,*), "1st run of LBFGS algorithm: initialize vars"
     Call initialize_lbfgs( ndim )
  Endif

     
  Call LBFGS( config_lbfgs%ndim, config_lbfgs%ncorr, x, func, grad, &
              config_lbfgs%ldiagco, config_lbfgs%rdiag, &
              config_lbfgs%iprint, eps_opt, config_lbfgs%eps_machine, &
              config_lbfgs%wk, config_lbfgs%iflag )
  istatus = config_lbfgs%iflag
  If ( config_lbfgs%iflag == 0 ) Then
     Write(6,*) "finish minimization procedure"
     Call Destroy_lbfgs()
  Endif

Endsubroutine


Subroutine destroy_lbfgs()
!$$$ subprogram documentation block
!             .                   .                    . 
!  Program name:
!    programmer: Da,Cheng       org: umd            date: 2015-MON-DD
!
!  Purpose:
!    
!
!  Revision history:
!    2015-MON-DD     Da,Cheng     -Creator
!
!  Input args:
!
!  Output args:
!
!
!
!$$$
  Implicit None
  
  If ( .not.( config_lbfgs%linit ) ) Return
  config_lbfgs%ndim  = 0
  config_lbfgs%ncorr = 0
  config_lbfgs%eps_machine = 0.d0
  If ( Allocated( config_lbfgs%wk ) ) Deallocate( config_lbfgs%wk )
  If ( Allocated( config_lbfgs%rdiag ) ) Deallocate( config_lbfgs%rdiag )

  config_lbfgs%linit = .FALSE.

Endsubroutine

Endmodule
