!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: CalcLocalDivergence.F90                        !
!    CONTAINS: subroutine CalcLocalDivergence             !
!                                                         !
!    PURPOSE: Compute the divergence of the intermediate  !
!     velocity at every point for the pressure            !
!     correction step                                     !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CalcLocalDivergence

    use param
    use local_arrays, only: vz,vy,vx,dph
    use decomp_2d, only: xstart,xend

    implicit none
    
    integer :: jc,jp,kc,kp,ic,ip
    real    :: dqcap

    do ic=xstart(3),xend(3)
        ip=ic+1
        do jc=xstart(2),xend(2)
            jp=jc+1
            do kc=1,nxm
                kp=kc+1
                dph(kc,jc,ic) = ((vz(kc,jc,ip)-vz(kc,jc,ic))*dz) + ((vy(kc,jp,ic)-vy(kc,jc,ic))*dy) + ((vx(kp,jc,ic)-vx(kc,jc,ic))*udx3m(kc))
            enddo
        enddo
    enddo

    return
    
end subroutine CalcLocalDivergence
