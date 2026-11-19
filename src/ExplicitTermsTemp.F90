!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: ExplicitTermsTemp.F90                          !
!    CONTAINS: subroutine ExplicitTermsTemp               !
!                                                         !
!    PURPOSE: Compute the non-linear terms associated to  !
!     the temperature.                                    !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ExplicitTermsTemp

    use param
    use local_arrays, only: vy,vx,temp,vz,hro
    use decomp_2d, only: xstart,xend

    implicit none

    integer :: im,ic,ip
    integer :: jm,jc,jp
    integer :: km,kc,kp
    real    :: htx,hty,htz
    real    :: dyyt,dzzt
    real    :: tplo,tphi

    do ic=xstart(3),xend(3)
        im=ic-1
        ip=ic+1
        do jc=xstart(2),xend(2)
            jm=jc-1
            jp=jc+1
            do kc=1,nxm
                km=kc-1
                kp=kc+1

                if (kc.eq.1) then
                    tplo = (temp(kc,jc,ic)*real(1-inslws)) + (tempbp*real(inslws))
                    tphi = (dx3c(kp)*temp(kc,jc,ic) + dx3c(kc)*temp(kp,jc,ic))/(2.0d0*g3rc(kp))
                else if (kc.eq.nxm) then
                    tplo = (dx3c(km)*temp(kc,jc,ic) + dx3c(kc)*temp(km,jc,ic))/(2.0d0*g3rc(kc))
                    tphi = (temp(kc,jc,ic)*real(1-inslwn)) + (temptp*real(inslwn))
                else
                    tplo = (dx3c(km)*temp(kc,jc,ic) + dx3c(kc)*temp(km,jc,ic))/(2.0d0*g3rc(kc))
                    tphi = (dx3c(kp)*temp(kc,jc,ic) + dx3c(kc)*temp(kp,jc,ic))/(2.0d0*g3rc(kp))
                end if

                ! d(vz.temp)/dz
                htz = ((vz(kc,jc,ip)*(temp(kc,jc,ip)+temp(kc,jc,ic))) - (vz(kc,jc,ic)*(temp(kc,jc,ic)+temp(kc,jc,im))))*dz*0.5d0

                ! d(vy.temp)/dy
                hty = ((vy(kc,jp,ic)*(temp(kc,jp,ic)+temp(kc,jc,ic))) - (vy(kc,jc,ic)*(temp(kc,jc,ic)+temp(kc,jm,ic))))*dy*0.5d0

                ! d(vx.temp)/dx
                htx = ((vx(kp,jc,ic)*tphi) - (vx(kc,jc,ic)*tplo))*udx3m(kc)

                dzzt = (temp(kc,jc,ip) -2.0d0*temp(kc,jc,ic) + temp(kc,jc,im))*dzq/pec
                dyyt = (temp(kc,jp,ic) -2.0d0*temp(kc,jc,ic) + temp(kc,jm,ic))*dyq/pec

                hro(kc,jc,ic) = -(htx+hty+htz) + (real(1-flagimp)*(dyyt+dzzt))

            enddo
        enddo
    enddo

    return

end subroutine ExplicitTermsTemp
