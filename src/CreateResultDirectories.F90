!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: CreateResultDirectories.F90                    !
!    CONTAINS: subroutine CreateResultDirectories         !
!                                                         !
!    PURPOSE: Create directories for saving results       !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CreateResultDirectories

    implicit none

    call system('mkdir -p Results')
    call system('mkdir -p Results/Temp')
    call system('mkdir -p Results/Vx')
    call system('mkdir -p Results/Vy')
    call system('mkdir -p Results/Vz')
    call system('mkdir -p Results/Stats')
    call system('mkdir -p Results/Grid')

end subroutine CreateResultDirectories