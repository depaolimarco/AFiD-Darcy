!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: ImplicitAndUpdateTemp.F90                      !
!    CONTAINS: subroutine ImplicitAndUpdateTemp           !
!                                                         !
!    PURPOSE: Compute the linear terms associated to      !
!     the temperature and call the implicit solver.       !
!     After this routine, the temperature has been        !
!     updated to the new timestep                         !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ImplicitAndUpdateTemp

    use param
    use local_arrays, only: temp,hro,rutemp,rhsx,dten
    use decomp_2d, only: xstart,xend

    implicit none

    integer :: im,jm,km
    integer :: ic,jc,kc
    integer :: ip,jp,kp
    real    :: alpec
    real    :: dxxt,dyyt,dzzt
    real    :: tplo,tpcc,tphi

    alpec=al/pec

    do ic=xstart(3),xend(3)
        im=ic-1
        ip=ic+1
        do jc=xstart(2),xend(2)
            jm=jc-1
            jp=jc+1
            do kc=1,nxm
                km=kc-1
                kp=kc+1

                tpcc = temp(kc,jc,ic)

                if (kc.eq.1) then
                    tplo = (tpcc*real(1-inslws)) + (((2.0d0*tempbp)-tpcc)*real(inslws))
                    tphi = temp(kp,jc,ic)
                else if (kc.eq.nxm) then
                    tplo = temp(km,jc,ic)
                    tphi = (tpcc*real(1-inslwn)) + (((2.0d0*temptp)-tpcc)*real(inslwn))
                else
                    tplo = temp(km,jc,ic)
                    tphi = temp(kp,jc,ic)
                end if

                ! Calculate second derivatives of temperature.
                ! These are computed for current timestep part of implicit Crank–Nicolson scheme.

                dzzt = (temp(kc,jc,ip) -2.0d0*temp(kc,jc,ic) + temp(kc,jc,im))*dzq
                dyyt = (temp(kc,jp,ic) -2.0d0*temp(kc,jc,ic) + temp(kc,jm,ic))*dyq
                dxxt = tphi*ap3sk(kc) + tpcc*ac3sk(kc) + tplo*am3sk(kc)

                ! Calculate right hand side of Eq. 5 (VO96)

                rhsx(kc,jc,ic) = (ga*hro(kc,jc,ic) + ro*rutemp(kc,jc,ic) + alpec*(dxxt + (real(flagimp)*(dyyt + dzzt))))*dt

                ! Store the non-linear terms for the calculation of the next timestep

                rutemp(kc,jc,ic) = hro(kc,jc,ic)

            enddo
        enddo
    enddo

    !  Solve equation and update temperature

    call SolveImpEqnUpdate_Temp
    
    return

end subroutine ImplicitAndUpdateTemp