module mod_type
!$$$  program documentation block
!         .           .            .
!  program name: mod_type
!    programmer: da,cheng        org: umd      date: 2015-Jan-05
!
!  purpose:
!
!  revision history:
!    2015-Jan-05     da    - creator
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

  public

  integer,parameter :: isp = 4
  integer,parameter :: idp = 8

  integer,parameter :: rsp = kind(0.0e0)
  integer,parameter :: rdp = kind(0.0d0)
  !integer,parameter :: rdef = rsp
  integer,parameter :: rdef = rdp

endmodule

