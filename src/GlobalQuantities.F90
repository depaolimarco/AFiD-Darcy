!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: GlobalQuantities.F90                           !
!    CONTAINS: subroutine GlobalQuantities                !
!                                                         !
!    PURPOSE: Calculate maximum velocity and temperature. !
!     volume averaged Nusselt number and Reynolds number  !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine GlobalQuantities

    use param
    use local_arrays, only: vy,vx,vz,temp
    use decomp_2d, only: xstart,xend
    use mpih

    implicit none

    integer :: jc,kc,kp,ic
    real    :: anusin,vol
    real    :: vy_rms_vol,vz_rms_vol
    real    :: vx_rms_vol,vzvyvx_rms_vol
    real    :: res_dummy
    logical :: fexist

    vmax(1)         = -huge(0.0d0)
    vmax(2)         = -huge(0.0d0)
    vmax(3)         = -huge(0.0d0)
    tempmax         = -huge(0.0d0)
    tempmin         =  huge(0.0d0)
    tempm           =  0.0d0
    anusin          =  0.d0
    vx_rms_vol      =  0.0d0
    vy_rms_vol      =  0.0d0
    vz_rms_vol      =  0.0d0
    vzvyvx_rms_vol  =  0.0d0
    vmax            =  0.0d0
    vol             =  1.d0/(alx3*dx*real(nzm)*real(nym))

    do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
            do kc=1,nxm
                kp             = kc + 1
                vmax(1)        = max(vmax(1),abs(vz(kc,jc,ic)))
                vmax(2)        = max(vmax(2),abs(vy(kc,jc,ic)))
                vmax(3)        = max(vmax(3),abs(vx(kc,jc,ic)))
                tempmax        = max(tempmax,temp(kc,jc,ic))
                tempmin        = min(tempmin,temp(kc,jc,ic))
                anusin         = anusin + 0.5d0*temp(kc,jc,ic)*(vx(kc,jc,ic)+vx(kp,jc,ic))*g3rc(kc)
                tempm          = tempm  + temp(kc,jc,ic)*g3rc(kc)
                vx_rms_vol     = vx_rms_vol + g3rc(kc)*vx(kc,jc,ic)**2
                vy_rms_vol     = vy_rms_vol + g3rc(kc)*vy(kc,jc,ic)**2
                vz_rms_vol     = vz_rms_vol + g3rc(kc)*vz(kc,jc,ic)**2
                vzvyvx_rms_vol = vzvyvx_rms_vol + g3rc(kc)*(vz(kc,jc,ic)**2 + vy(kc,jc,ic)**2 + 0.5d0*(vx(kc,jc,ic)**2 + vx(kp,jc,ic)**2))
            enddo
        enddo
    enddo

    call MpiSumRealScalar(tempm,res_dummy)
    tempm = res_dummy
    call MpiSumRealScalar(anusin,res_dummy)
    anusin = res_dummy
    call MpiSumRealScalar(vx_rms_vol,res_dummy)
    vx_rms_vol = res_dummy
    call MpiSumRealScalar(vy_rms_vol,res_dummy)
    vy_rms_vol = res_dummy
    call MpiSumRealScalar(vz_rms_vol,res_dummy)
    vz_rms_vol = res_dummy
    call MpiSumRealScalar(vzvyvx_rms_vol,res_dummy)
    vzvyvx_rms_vol = res_dummy
    call MpiMinRealScalar(tempmin,res_dummy)
    tempmin = res_dummy
    call MpiMaxRealScalar(tempmax,res_dummy)
    tempmax = res_dummy
    call MpiMaxRealScalar(vmax(1),res_dummy)
    vmax(1) = res_dummy
    call MpiMaxRealScalar(vmax(2),res_dummy)
    vmax(2) = res_dummy
    call MpiMaxRealScalar(vmax(3),res_dummy)
    vmax(3) = res_dummy

    if(ismaster) then

        anusin = 1.d0 + pec*anusin*vol

        inquire(file='Results/nu_vol.out',exist=fexist)
        open(95,file='Results/nu_vol.out',status='unknown',access='sequential',position='append')
        if (.not.fexist) then
            write(95,'(2(A16,X))') 'Time', 'Nusselt(VolAvg)'
        end if
        write(95,'(2(E16.8,X))') time, anusin
        close(95)

        tempm           = tempm*vol
        vx_rms_vol      = dsqrt(vx_rms_vol*vol)*ren
        vy_rms_vol      = dsqrt(vy_rms_vol*vol)*ren
        vz_rms_vol      = dsqrt(vz_rms_vol*vol)*ren
        vzvyvx_rms_vol  = dsqrt(vzvyvx_rms_vol*vol)*ren

        inquire(file='Results/rms_vel.out',exist=fexist)
        open(94,file='Results/rms_vel.out',status='unknown',position='append',access='sequential')
        if (.not.fexist) then
            write(94,'(5(A16,X))') 'Time', 'Re_RMS_Vz', 'Re_RMS_Vy', 'Re_RMS_Vx', 'Re_RMS_VxVyVz'
        end if
        write(94,'(5(E16.8,X))') time, vz_rms_vol, vy_rms_vol, vx_rms_vol, vzvyvx_rms_vol
        close(94)

    endif

    return

end subroutine GlobalQuantities
