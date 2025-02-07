!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: CalcOutput.F90                                 !
!    CONTAINS: subroutine CalcOutput                      !
!                                                         !
!    PURPOSE: Calculate the Nusselt number and other      !
!    useful quantities at the top and bottom plates       !
!    and output to a file.                                !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CalcOutput

    use param
    use local_arrays,only: vz,vy,vx,temp
    use stat_arrays, only: dissth,disste
    use mpih
    use decomp_2d, only: xstart,xend

    implicit none

    integer         :: ic,jc,kc
    integer         :: im,jm,km
    integer         :: ip,jp,kp
    integer         :: kh

    integer         :: isamples

    real            :: nu_th,volt
    real            :: dissipte,dissipth
    real            :: idisste(nxm),idissth(nxm),idarray(nxm)

    real            :: tplo,tpcc,tphi
    real            :: vellow,velupp

    real            :: txm,txp,ty,tz
    real            :: vxx,vxy1,vxz1,vxy2,vxz2,vyx1,vyx2,vyy,vyz,vzx1,vzx2,vzy,vzz

    real            :: c1,dc1dt
    real            :: c2,dc2dt
    real            :: qlow,qupp
    real            :: res1,res2

    real            :: res_dummy

    character*50    :: filename,dsetname
    character*8     :: citime

    logical         :: fexist

    nu_th       = 0.0d0
    c1          = 0.0d0
    c2          = 0.0d0

    nuslow      = 0.0d0
    nusupp      = 0.0d0


    idisste(:)  = 0.0d0
    idissth(:)  = 0.0d0

    ! Computation of Nusselt Numbers
    do ic=xstart(3),xend(3)
        im = ic - 1
        ip = ic + 1
        do jc=xstart(2),xend(2)
            jm = jc - 1
            jp = jc + 1

            vellow      = dsqrt(0.5d0*((vy(1,  jc,ic)**2.0d0) + (vy(1,  jp,ic)**2.0d0) + (vz(1,  jc,ic)**2.0d0) + (vz(1,  jc,ip)**2.0d0)))
            velupp      = dsqrt(0.5d0*((vy(nxm,jc,ic)**2.0d0) + (vy(nxm,jp,ic)**2.0d0) + (vz(nxm,jc,ic)**2.0d0) + (vz(nxm,jc,ip)**2.0d0)))
            nuslow      = nuslow      - real(inslws)*(temp(1,  jc,ic) - tempbp)*(2.0d0*udx3m(1)  )
            nusupp      = nusupp      - real(inslwn)*(temptp - temp(nxm,jc,ic))*(2.0d0*udx3m(nxm))

        enddo
    end do

    nuslow      = nuslow        / real(nzm*nym)
    nusupp      = nusupp        / real(nzm*nym)

    call MpiSumRealScalar(nuslow,res_dummy)
    nuslow      = res_dummy
    call MpiSumRealScalar(nusupp,res_dummy)
    nusupp      = res_dummy


    ! Computation of kinetic and thermal dissipation rates
    do ic=xstart(3),xend(3)
        ip = ic+1
        im = ic-1
        do jc=xstart(2),xend(2)
            jp = jc+1
            jm = jc-1
            do kc=1,nxm
                kp = kpv(kc)
                kh = kc+1
                km = kmv(kc)

                if (kc.eq.1) then
                    tplo = (temp(kc,jc,ic)*real(1-inslws)) + (((2.0d0*tempbp)-temp(kc,jc,ic))*real(inslws))
                    tphi = temp(kp,jc,ic)
                else if (kc.eq.nxm) then
                    tplo = temp(km,jc,ic)
                    tphi = (temp(kc,jc,ic)*real(1-inslwn)) + (((2.0d0*temptp)-temp(kc,jc,ic))*real(inslwn))
                else
                    tplo = temp(km,jc,ic)
                    tphi = temp(kp,jc,ic)
                end if

                tpcc     = temp(kc,jc,ic)

                txm      = (tpcc - tplo)*udx3c(kc)
                txp      = (tphi - tpcc)*udx3c(kh)

                ty       = (temp(kc,jp,ic)-temp(kc,jc,ic))*dy
                tz       = (temp(kc,jc,ip)-temp(kc,jc,ic))*dz

                vxx      = (vx(kp,jc,ic)-vx(kc,jc,ic))*udx3m(kc)

                vxy1     = (vx(kc,jp,ic)-vx(kc,jc,ic))*dy
                vxz1     = (vx(kc,jc,ip)-vx(kc,jc,ic))*dz
                vxy2     = (vx(kp,jp,ic)-vx(kp,jc,ic))*dy
                vxz2     = (vx(kp,jc,ip)-vx(kp,jc,ic))*dz

                vyy      = (vy(kc,jp,ic)-vy(kc,jc,ic))*dy
                vyz      = (vy(kc,jc,ip)-vy(kc,jc,ic))*dz

                vzy      = (vz(kc,jp,ic)-vz(kc,jc,ic))*dy
                vzz      = (vz(kc,jc,ip)-vz(kc,jc,ic))*dz

                vyx1     = (vy(kc,jc,ic)-vy(km,jc,ic))*udx3c(kc)
                vyx2     = (vy(kp,jc,ic)-vy(kc,jc,ic))*udx3c(kh)
                vzx1     = (vz(kc,jc,ic)-vz(km,jc,ic))*udx3c(kc)
                vzx2     = (vz(kp,jc,ic)-vz(kc,jc,ic))*udx3c(kh)

                dissipte = (vxx*vxx)+(0.5*((vxy1*vxy1)+(vxz1*vxz1)+(vxy2*vxy2)+(vxz2*vxz2)+(vyx1*vyx1)+(vyx2*vyx2)+(vzx1*vzx1)+(vzx2*vzx2)))+(vyy*vyy)+(vyz*vyz)+(vzy*vzy)+(vzz*vzz)
                dissipth = (0.5d0*((txm**2.0d0) + (txp**2.0d0))) + (ty**2.0d0) + (tz**2.0d0)
                nu_th    = nu_th + dissipth*dx3c(kc)
                c1       = c1    + temp(kc,jc,ic)*dx3c(kc)
                c2       = c2    + (temp(kc,jc,ic)**2.0d0)*dx3c(kc)

                if (disscal) then
                    idisste(kc) = idisste(kc) + (dissipte*dt) / (real(nym)*real(nzm))
                    idissth(kc) = idissth(kc) + (dissipth*dt) / (real(nym)*real(nzm))
                    if (time.gt.tsta) then
                        disste(kc) =  disste(kc) + (dissipte*dt) / (real(nym)*real(nzm))
                        dissth(kc) =  dissth(kc) + (dissipth*dt) / (real(nym)*real(nzm))
                    end if
                end if

            end do
        end do
    end do

    call MpiSumRealScalar(nu_th,res_dummy)
    nu_th = res_dummy
    call MpiSumRealScalar(c1,res_dummy)
    c1    = res_dummy
    call MpiSumRealScalar(c2,res_dummy)
    c2    = res_dummy

    volt = 1.d0/(real(nxm)*real(nzm)*real(nym))

    if(ismaster) then

        nu_th = nu_th*volt
        c1    = c1   *volt
        dc1dt = (c1 - c1prev)/dt
        c2    = c2   *volt
        dc2dt = (c2 - c2prev)/dt

        qupp  = nusupp 
        qlow  = nuslow 
        res1  = dc1dt + ((qupp - qlow)/ray)

        qupp  = real(inslwn)*temptp*(nusupp)
        qlow  = real(inslws)*tempbp*(nuslow)
        res2  = (0.5d0*dc2dt) + ((qupp - qlow)/ray) + (nu_th/ray)

    endif

    ! Writing to file
    if(ismaster) then
        inquire(file="Results/nu_darcy.out",exist=fexist)
        open(97,file="Results/nu_darcy.out",status='unknown',access='sequential',position='append')
        if (.not.fexist) then
            write(97,'(13(A16,X))') 'Time', 'DeltaTime', 'NuB(dT/dy)', 'NuT(dT/dy)', 'chi', 'c1', 'dc1dt', 'res1', 'c2', 'dc2dt', 'res2'
        end if
        write(97,'(13(E16.8,X))') time, dt, nuslow, nusupp, nu_th, c1, dc1dt, res1, c2, dc2dt, res2
        close(97)
    endif

    if (disscal .and. instprof) then

        call MpiSumReal1D(idisste,idarray,nxm)
        idisste(:) = idarray(:)
        call MpiSumReal1D(idissth,idarray,nxm)
        idissth(:) = idarray(:)

        if (ismaster) then

            filename = trim('Results/instprof_master.h5')
            dsetname = trim('samples')
            call HdfSerialReadIntScalar(dsetname,filename,isamples)
            isamples = isamples + 1
            write (citime,"(i8.8)") isamples
            dsetname = trim("disste/"//citime)
            call HdfSerialWriteReal1D(dsetname,filename,idisste,nxm)
            dsetname = trim("dissth/"//citime)
            call HdfSerialWriteReal1D(dsetname,filename,idissth,nxm)

        end if
        
    end if

    return

end subroutine CalcOutput

subroutine CalcCNPrev

    use param
    use local_arrays,only: temp
    use mpih
    use decomp_2d, only: xstart,xend    

    implicit none

    integer         :: ic,jc,kc
    real            :: c1,c2,res_dummy

    c1 = 0.0d0
    c2 = 0.0d0
    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
            do kc=1,nxm
                c1 = c1 + temp(kc,jc,ic)*dx3c(kc)
                c2 = c2 + (temp(kc,jc,ic)**2.0d0)*dx3c(kc)
            end do
        end do
    end do

    call MpiSumRealScalar(c1,res_dummy)
    c1 = res_dummy

    call MpiSumRealScalar(c2,res_dummy)
    c2 = res_dummy

    if (ismaster) then
        c1prev = c1/(real(nxm)*real(nzm)*real(nym))
        c2prev = c2/(real(nxm)*real(nzm)*real(nym))
    end if

end subroutine CalcCNPrev