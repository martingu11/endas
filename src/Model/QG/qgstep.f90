!! Copyright (C) 2008 Pavel Sakov
!! 
!! This file is part of EnKF-Matlab. EnKF-Matlab is a free software. See 
!! LICENSE for details.

! File:           qgstep.f90
!
! Created:        31/08/2007
!
! Last modified:  08/02/2008
!
! Author:         Pavel Sakov
!                 CSIRO Marine and Atmospheric Research
!                 NERSC
!
! Purpose:        Fortran code for QG model. Integrators.
!
! Description:    
!
! Revisions:

module qgstep_mod
  use qgflux_mod

contains

  subroutine qg_step_rk4(t, dt, rkb, rkh, rkh2, F, r, PSI, Q) bind (c, name="qg_step_rk4")
    use iso_c_binding, only: c_int, c_double
    real(c_double), intent(in), value :: t, dt, rkb, rkh, rkh2, F, r
    real(c_double), dimension(N, M), intent(inout) :: Q, PSI

    real(8), dimension(N, M) :: QFLUX1, QFLUX2, QFLUX3, QFLUX4
    real(8), dimension(N, M) :: PP
    real(8), dimension(N, M) :: Q2, Q3, Q4
    real(8) :: tt

    ! Given vorticity Q, this call calculates its flux QFLUX1. 
    ! Solves for PSI1 as a by-product, using PSI as the first guess
    !
    call qg_flux(t, rkb, rkh, rkh2, F, r, Q, PSI, PP, QFLUX1)
    tt = t + 0.5d0
    Q2 = Q + (0.5d0 * dt) * QFLUX1
    call qg_flux(tt, rkb, rkh, rkh2, F, r, Q2, PP, PSI, QFLUX2)
    Q3 = Q + (0.5d0 * dt) * QFLUX2
    call qg_flux(tt, rkb, rkh, rkh2, F, r, Q3, PSI, PP, QFLUX3)
    Q4 = Q + dt * QFLUX3
    tt = t + dt
    call qg_flux(tt, rkb, rkh, rkh2, F, r, Q4, PP, PSI, QFLUX4)

    !t = t + dt
    Q = Q + (QFLUX1 + 2.0d0 * (QFLUX2 + QFLUX3) + QFLUX4) * (dt / 6.0d0)
  end subroutine qg_step_rk4

 

end module qgstep_mod
