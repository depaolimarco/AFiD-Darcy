!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: CorrectVelocity.F90                            !
!    CONTAINS: subroutine CorrectVelocity                 !
!                                                         !
!    PURPOSE: Update velocities with the pressure         !
!     correction to enforce incompresibility              !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CorrectVelocity
    use param
    use local_arrays, only: vy,vx,dphhalo,vz
    use decomp_2d, only: xstart,xend
    implicit none
    integer :: jc,jm,kc,km,ic,im
    real    :: locdph

    do ic=xstart(3),xend(3)
        im=ic-1
        do jc=xstart(2),xend(2)
            jm=jc-1
            do kc=1,nxm
                km=kmv(kc)
                locdph=dphhalo(kc,jc,ic)
                vx(kc,jc,ic) = vx(kc,jc,ic) - ((locdph-dphhalo(km,jc,ic))*udx3c(kc))
                vy(kc,jc,ic) = vy(kc,jc,ic) - ((locdph-dphhalo(kc,jm,ic))*dy)
                vz(kc,jc,ic) = vz(kc,jc,ic) - ((locdph-dphhalo(kc,jc,im))*dz)
            enddo
        enddo
    enddo

    return
end subroutine CorrectVelocity

