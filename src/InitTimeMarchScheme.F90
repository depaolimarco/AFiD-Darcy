!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: InitTimeMarchScheme.F90                        !
!    CONTAINS: subroutine InitTimeMarchScheme             !
!                                                         !
!    PURPOSE: Initialize the time-marching constants for  !
!     the integrator                                      !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine InitTimeMarchScheme

    use param

    implicit none

    ! integer ns

    if(nsst.gt.1) then
        gam(1)=8.d0/15.d0
        gam(2)=5.d0/12.d0
        gam(3)=3.d0/4.d0
        rom(1)=0.d0
        rom(2)=-17.d0/60.d0
        rom(3)=-5.d0/12.d0

        if(ismaster) then
            write(6,'(A61)') '     Time-stepping scheme             = III order Runge-Kutta'
            write(6,'(A40,3F8.3)') '     gam                              = ',(gam(ns),ns=1,nsst)
            write(6,'(A40,3F8.3)') '     ro                               = ',(rom(ns),ns=1,nsst)
        endif

    else
        gam(1)=1.5d0
        gam(2)=0.d0
        gam(3)=0.d0
        rom(1)=-0.5d0
        rom(2)=0.d0
        rom(3)=0.d0

        if(ismaster) then
            write(6,'(A54)') '     Time-stepping scheme             = Adams-Bashfort'
            write(6,'(A40,F8.3)') '     gam                              = ',gam(1)
            write(6,'(A40,F8.3)') '     ro                               = ',rom(1)
        endif
    endif

    do ns=1,nsst
        alm(ns)=(gam(ns)+rom(ns))
    end do

    return

end subroutine InitTimeMarchScheme