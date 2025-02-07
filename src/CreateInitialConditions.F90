!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: CreateInitialConditions.F90                    !
!    CONTAINS: subroutine CreateInitialConditions         !
!                                                         !
!    PURPOSE: Initialization routine. Sets initial        !
!     conditions for velocity and temperature             !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CreateInitialConditionsOld

    use param
    use local_arrays, only: vy,vx,temp,vz
    use decomp_2d, only: xstart,xend
    use mpih

    implicit none

    integer :: i,j,k
    real :: xxx,yyy,zzz,val
    real :: lvly,lvlz
    character*50 :: filname

    ! x -> wall normal
    ! y -> spanwise
    ! z -> streamwise

    lvly = real(floor(ylen/alx3))
    lvlz = real(floor(zlen/alx3))

    vx(:,:,:)   = 0.0d0
    vy(:,:,:)   = 0.0d0
    vz(:,:,:)   = 0.0d0

    do i=xstart(3),xend(3)
        do j=xstart(2),xend(2)
            do k=1,nxm
                xxx = xm(k)/alx3 
                yyy = ym(j)/ylen
                zzz = zm(i)/zlen
                val = xxx + (eps*(dcos(2.0d0*pi*lvlz*zzz)*dcos(2.0d0*pi*lvly*yyy)*dsin(pi*xxx))/pi)
                temp(k,j,i) = tempbp + (temptp-tempbp)*val
            enddo
        enddo
    enddo

    filname = trim('Results/Temp/continua_temp00000000.h5')
    call HdfWriteRealHalo3D(filname,temp)
    filname = trim('Results/Vx/continua_vx00000000.h5')
    call HdfWriteRealHalo3D(filname,vx)
    filname = trim('Results/Vy/continua_vy00000000.h5')
    call HdfWriteRealHalo3D(filname,vy)
    filname = trim('Results/Vz/continua_vz00000000.h5')
    call HdfWriteRealHalo3D(filname,vz)

    return

end subroutine CreateInitialConditionsOld

subroutine CreateInitialConditions

    use param
    use local_arrays, only: vy,vx,temp,vz
    use decomp_2d, only: xstart,xend
    use mpih

    implicit none

    integer         :: i,j,k
    real            :: xxx,ttt
    real            :: rnum
    character*50    :: filname

    ! x -> wall normal
    ! y -> spanwise
    ! z -> streamwise

    vx(:,:,:)   = 0.0d0
    vy(:,:,:)   = 0.0d0
    vz(:,:,:)   = 0.0d0

    if (initime.le.0.0d0) then
        time = dt
    else
        time = initime/ray
    end if
    ttt  = time*ray

    do i=xstart(3),xend(3)
        do j=xstart(2),xend(2)
            do k=1,nxm
                call random_number(rnum)
                rnum = 2.0d0*eps*(rnum - 0.5d0)
                if (inslws.le.0) then
                    if (inslwn.le.0) then
                        ! This is the Rayleigh–Taylor case
                        xxx = (xm(k)/alx3 - 0.5d0)*ray
                        temp(k,j,i) = tempbp + (temptp-tempbp)*(0.5d0*(1+erf(xxx/(2.0d0*dsqrt(ttt)))) + eps*rnum)
                    else
                        ! This is the one-sided case
                        xxx = (xm(k)/alx3 - 1.0d0)*ray
                        temp(k,j,i) = tempbp + (temptp-tempbp)*(1.0d0+erf(xxx/(2.0d0*dsqrt(ttt))))*(1.0d0 + eps*rnum)
                    end if
                else
                    if (inslwn.le.0) then
                        ! This is the one-sided case
                        xxx = (xm(k)/alx3)*ray
                        temp(k,j,i) = temptp - (temptp-tempbp)*(1.0d0-erf(xxx/(2.0d0*dsqrt(ttt))))*(1.0d0 + eps*rnum)
                    else
                        ! This is the Rayleigh–Bénard case
                        xxx = (xm(k)/alx3)
                        temp(k,j,i) = tempbp + ((temptp-tempbp)*xxx) + (eps*rnum)
                    end if
                end if
            enddo
        enddo
    enddo

    filname = trim('Results/Temp/continua_temp00000000.h5')
    call HdfWriteRealHalo3D(filname,temp)
    filname = trim('Results/Vx/continua_vx00000000.h5')
    call HdfWriteRealHalo3D(filname,vx)
    filname = trim('Results/Vy/continua_vy00000000.h5')
    call HdfWriteRealHalo3D(filname,vy)
    filname = trim('Results/Vz/continua_vz00000000.h5')
    call HdfWriteRealHalo3D(filname,vz)

    return

end subroutine CreateInitialConditions