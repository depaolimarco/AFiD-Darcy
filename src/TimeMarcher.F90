!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: TimeMarcher.F90                                !
!    CONTAINS: subroutine TimeMarcher                     !
!                                                         !
!    PURPOSE: Main time integrating routine, which calls  !
!     other subroutines for calculating the Navier-Stokes !
!     equations and advancing velocity and temperature in !
!     time                                                !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine TimeMarcher

    use param
    use local_arrays
    use mpih
    use decomp_2d

    implicit none

    beta=dt/ren*0.5d0

    do ns=1,nsst

        ! Coefficients for time marching integration (alpha, gamma, rho)

        al=alm(ns)
        ga=gam(ns)
        ro=rom(ns)

        call ExplicitTermsTemp

        call ImplicitAndUpdateTemp
        call update_halo(temp,lvlhalo)

        call SetVelocityField
        call update_halo(vy,lvlhalo)
        call update_halo(vz,lvlhalo)

        call CalcLocalDivergence
        call SolvePressureCorrection
        call CorrectPressure
        call update_halo(dphhalo,lvlhalo)
        call update_halo(pr,lvlhalo)

        call CorrectVelocity
        call update_halo(vx,lvlhalo)
        call update_halo(vy,lvlhalo)
        call update_halo(vz,lvlhalo)

        time = time + (al*dt)

    enddo

    return
    
end subroutine TimeMarcher
