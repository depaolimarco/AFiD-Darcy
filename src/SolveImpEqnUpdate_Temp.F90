!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: SolveImpEqnUpdate_Temp.F90                     !
!    CONTAINS: subroutine SolveImpEqnUpdate_Temp          !
!                                                         !
!    PURPOSE: Inverts the implicit equation for           !
!     temperature, and updates it to time t+dt            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SolveImpEqnUpdate_Temp

    use param
    use local_arrays, only : temp,rhsx,rhsy,rhsz
    use decomp_2d

    implicit none

    integer :: ic,jc,kc,info
    real    :: betadx,gvalue

    real    :: acil_b,acjl_b,ackl_b

    integer :: ipiv(nzm),ipjv(nym),ipkv(nxm)

    real    :: amil(nzm),acil(nzm),apil(nzm)
    real    :: amjl(nym),acjl(nym),apjl(nym)
    real    :: amkl(nxm),ackl(nxm),apkl(nxm)

    real    :: amiT(nzm-1),aciT(nzm),apiT(nzm-1),appi(nzm-2)
    real    :: amjT(nym-1),acjT(nym),apjT(nym-1),appj(nym-2)
    real    :: amkT(nxm-1),ackT(nxm),apkT(nxm-1),appk(nxm-2)

    real    :: rz1d(nzm),ry1d(nym),rx1d(nxm)

    real    :: fluz(nzm),flvz(nzm),flqz
    real    :: fluy(nym),flvy(nym),flqy

    betadx   =  0.5d0*al*dt/pec

    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)

            do kc=1,nxm
                ackl_b   =  1.0d0/(1.0d0 - ac3sk(kc)*betadx)
                amkl(kc) = -am3sk(kc)*betadx*ackl_b
                ackl(kc) =  1.0d0
                apkl(kc) = -ap3sk(kc)*betadx*ackl_b
                rx1d(kc) =  rhsx(kc,jc,ic)*ackl_b
            end do

            ackl(1)   = ackl(1)   + (real(1-inslws)*amkl(1))   - (real(inslws)*amkl(1))
            ackl(nxm) = ackl(nxm) + (real(1-inslwn)*apkl(nxm)) - (real(inslwn)*apkl(nxm))

            amkT=amkl(2:nxm)
            apkT=apkl(1:nxm-1)
            ackT=ackl(1:nxm)

            call dgttrf(nxm,amkT,ackT,apkT,appk,ipkv,info)
            call dgttrs('N',nxm,1,amkT,ackT,apkT,appk,ipkv,rx1d,nxm,info)

            do kc=1,nxm
                rhsx(kc,jc,ic) = rx1d(kc)
            end do

        end do
    end do

    if (flagimp.eq.1) then

        call transpose_x_to_y(rhsx,rhsy)

        do ic=ystart(3),yend(3)
            do kc=ystart(1),yend(1)

                do jc=1,nym
                    acjl_b   =  1.0d0/(1.0d0 - ac2sj*betadx)
                    amjl(jc) = -am2sj*betadx*acjl_b
                    acjl(jc) =  1.0d0
                    apjl(jc) = -ap2sj*betadx*acjl_b
                    ry1d(jc) =  rhsy(kc,jc,ic)*acjl_b
                    fluy(jc) =  0.0d0
                    flvy(jc) =  0.0d0
                end do

                gvalue    = (amjl(1)*apjl(nym))**0.5d0  ! This makes the matrix symmetric

                acjl(1)   = acjl(1)   - gvalue
                acjl(nym) = acjl(nym) - ((amjl(1)*apjl(nym))/gvalue)

                fluy(1)   = gvalue
                fluy(nym) = apjl(nym)

                flvy(1)   = 1.0d0
                flvy(nym) = amjl(1)/gvalue

                amjT = amjl(2:nym)
                acjT = acjl(1:nym)
                apjT = apjl(1:nym-1)

                call dgttrf(nym,amjT,acjT,apjT,appj,ipjv,info)

                call dgttrs('N',nym,1,amjT,acjT,apjT,appj,ipjv,ry1d,nym,info)
                call dgttrs('N',nym,1,amjT,acjT,apjT,appj,ipjv,fluy,nym,info)

                flqy = (flvy(1)*ry1d(1) + flvy(nym)*ry1d(nym))/(1.0d0 + (flvy(1)*fluy(1) + flvy(nym)*fluy(nym)))

                do jc=1,nym                
                    rhsy(kc,jc,ic) = ry1d(jc) - (flqy*fluy(jc))
                enddo
                
            end do
        end do

        if (nzm.gt.1) then
        
            call transpose_y_to_z(rhsy,rhsz)

            do jc=zstart(2),zend(2)
                do kc=zstart(1),zend(1)

                    do ic=1,nzm
                        acil_b   =  1.0d0/(1.0d0 - ac1si*betadx)
                        amil(ic) = -am1si*betadx*acil_b
                        acil(ic) =  1.0d0
                        apil(ic) = -ap1si*betadx*acil_b
                        rz1d(ic) =  rhsz(kc,jc,ic)*acil_b
                        fluz(ic) =  0.0d0
                        flvz(ic) =  0.0d0
                    end do

                    gvalue    = (amil(1)*apil(nzm))**0.5d0  ! This makes the matrix symmetric

                    acil(1)   =  acil(1)   - gvalue
                    acil(nzm) =  acil(nzm) - ((amil(1)*apil(nzm))/gvalue)

                    fluz(1)   = gvalue
                    fluz(nzm) = apil(nzm)

                    flvz(1)   = 1.0d0
                    flvz(nzm) = amil(1)/gvalue

                    amiT = amil(2:nzm)
                    aciT = acil(1:nzm)
                    apiT = apil(1:nzm-1)

                    call dgttrf(nzm,amiT,aciT,apiT,appi,ipiv,info)

                    call dgttrs('N',nzm,1,amiT,aciT,apiT,appi,ipiv,rz1d,nzm,info)
                    call dgttrs('N',nzm,1,amiT,aciT,apiT,appi,ipiv,fluz,nzm,info)

                    flqz = (flvz(1)*rz1d(1) + flvz(nzm)*rz1d(nzm))/(1.0d0 + (flvz(1)*fluz(1) + flvz(nzm)*fluz(nzm)))

                    do ic=1,nzm                
                        rhsz(kc,jc,ic) = rz1d(ic) - (flqz*fluz(ic))
                    enddo

                end do
            end do

            call transpose_z_to_x(rhsz,rhsx)

        else

            call transpose_y_to_x(rhsy,rhsx)

        end if

    end if

    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
            do kc=1,nxm
                temp(kc,jc,ic) = temp(kc,jc,ic) + rhsx(kc,jc,ic)
            end do
        end do
    end do

    return

end subroutine SolveImpEqnUpdate_Temp
