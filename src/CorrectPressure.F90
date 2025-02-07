!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: CorrectPressure.F90                            !
!    CONTAINS: subroutine CorrectPressure                 !
!                                                         !
!    PURPOSE: Copy pressure from dph without halo cells   !
!    into dphhalo with halo cells for computing gradients !
!    and update/correct the physical pressure pr with the !
!    correction dph.                                      !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CorrectPressure

    use param
    use local_arrays, only: pr,dph,dphhalo
    use decomp_2d, only: xstart,xend

    implicit none

    integer :: j,k,i

    !EP this copy can be avoided by changing transpose_x_to_y_real and
    !transpose_y_to_x_real so these routines can handle arrays with
    !halo. This copy is a defacto array temporary. Using inferred size
    !arrays in the transpose calls results in 5 more of these, and more
    !memory usage.  Time spent on this copy is 0.1% for 65^3 grid.

    do i=xstart(3),xend(3)
        do j=xstart(2),xend(2)
            do k=1,nxm
                dphhalo(k,j,i) = dph(k,j,i)
                pr(k,j,i) = pr(k,j,i) + dph(k,j,i)
            enddo
        enddo
    enddo

end subroutine CorrectPressure