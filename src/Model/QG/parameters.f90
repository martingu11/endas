!! Copyright (C) 2008 Pavel Sakov
!! 
!! This file is part of EnKF-Matlab. EnKF-Matlab is a free software. See 
!! LICENSE for details.

! File:           parameters.f90
!
! Created:        31/08/2007
!
! Last modified:  08/02/2008
!
! Author:         Pavel Sakov
!                 CSIRO Marine and Atmospheric Research
!                 NERSC
!
! Purpose:        Fortran code for QG model. Reads and stores the model
!                 parameters.
!
! Description:    
!
! Revisions:

module parameters_mod

  !use utils_mod
  implicit none

  save

  integer, parameter, private :: FID_PRM = 111

  ! This is the place to change the resolution of the model, by
  ! changing the parameter MREFIN:
  ! To make the ocean a rectangular (at your own risk!) change nx1,
  ! ny1 to (2,4), and adjust domain size (Lx, Ly)...
  !
  integer, parameter, public :: MREFIN = 7
  integer, parameter, public :: NX1 = 2
  integer, parameter, public :: NY1 = 2
  integer, parameter, public :: M = NX1 * 2 ** (MREFIN - 1) + 1
  integer, parameter, public :: N = NY1 * 2 ** (MREFIN - 1) + 1

  real(8), parameter, public :: PI = 3.14159265358979323d0

  real(8), parameter, public :: rf_coeff = 0.1 ! Roberts filter coefficient
                                               ! for the leap-frog scheme

  !real(8), public :: dt            ! time step
  !real(8), public :: rkb           ! bottom friction
  !real(8), public :: rkh           ! horizontal friction
  !real(8), public :: rkh2          ! biharmonic friction
  !real(8), public :: F             ! Froud number, (L / Rd)^2
  !real(8), public :: r             ! factor in front of nonlinear
                                   ! advection term
  integer, public :: rstart = 0
  !real(8), public :: tend
  !real(8), public :: dtout
  !integer, public :: verbose = VERBOSE_DEF
  !character(STRLEN), public :: scheme
  !character(STRLEN), public :: restartfname
  !character(STRLEN), public :: outfname


  real(8), dimension(N, M), public :: CURLT
  real(8), public :: dx
  real(8), public :: dy

  contains

  subroutine params_init() bind (c, name="qg_params_init")
    integer :: j
    integer :: m1 = M - 1

    dx = 1.0d0 / dble(N - 1)
    dy = 1.0d0 / dble(M - 1)

    do j = 1, M
       CURLT(:, j) = - 2.0d0 * PI * sin(2.0d0 * PI * (j - 1) / m1)
    end do

  end subroutine params_init

       
end module parameters_mod
