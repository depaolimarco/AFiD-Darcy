!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: StatRoutines.F90                               !
!    CONTAINS: subroutine CalcStats,WriteStats            !
!                                                         !
!    PURPOSE: Calculates and writes out statistics for    !
!     the flow field. All quantities are averaged in   	  !
!	  the two horizontal (homogeneous) directions.   	  !
!														  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CalcStats

    use param
    use local_arrays, only: vz,vy,vx,temp
    use decomp_2d, only: xstart,xend
    use stat_arrays
    use mpih

    implicit none

    integer :: ic,jc,kc
    integer :: ip,jp,kp
    integer :: im,jm,km

    real :: factor
    !! Cell center values calculated from linear interpolation of face center values
    real :: vxcc,vycc,vzcc,tempcc			
    real :: delta_y,delta_z

    nstatsamples = nstatsamples + 1
    tstat        = tstat        + dt

    factor = dt/(nym*nzm)

    delta_y = ylen/nym
    delta_z = zlen/nzm

    do ic=xstart(3),xend(3)

        ip = ic+1
        im = ic-1

        do jc=xstart(2),xend(2)

            jp = jc+1
            jm = jc-1

            do kc=1,nxm

                kp = kc+1
                km = kc-1

                !! Quantities based on linear interpolation at cell center

                vxcc 	= 0.5*(vx(kc,jc,ic) 	+ vx(kp,jc,ic))
                vycc 	= 0.5*(vy(kc,jc,ic) 	+ vy(kc,jp,ic))
                vzcc 	= 0.5*(vz(kc,jc,ic) 	+ vz(kc,jc,ip))
                tempcc 	= temp(kc,jc,ic)

                vxcc_m1(kc) 	= vxcc_m1(kc) 	+ vxcc*factor
                vycc_m1(kc) 	= vycc_m1(kc) 	+ vycc*factor
                vzcc_m1(kc) 	= vzcc_m1(kc) 	+ vzcc*factor
                tempcc_m1(kc)   = tempcc_m1(kc) + tempcc*factor

                vxcc_m2(kc) 	= vxcc_m2(kc) 	+ vxcc**2*factor
                vycc_m2(kc) 	= vycc_m2(kc) 	+ vycc**2*factor
                vzcc_m2(kc) 	= vzcc_m2(kc) 	+ vzcc**2*factor
                tempcc_m2(kc)   = tempcc_m2(kc) + tempcc**2*factor

                vxcc_vycc_m1(kc) 	= vxcc_vycc_m1(kc) 	 + (vxcc*vycc)*factor
                vxcc_vzcc_m1(kc)	= vxcc_vzcc_m1(kc) 	 + (vxcc*vzcc)*factor
                vxcc_tempcc_m1(kc) 	= vxcc_tempcc_m1(kc) + (vxcc*tempcc)*factor
                vycc_vzcc_m1(kc) 	= vycc_vzcc_m1(kc) 	 + (vycc*vzcc)*factor
                vycc_tempcc_m1(kc) 	= vycc_tempcc_m1(kc) + (vycc*tempcc)*factor
                vzcc_tempcc_m1(kc) 	= vzcc_tempcc_m1(kc) + (vzcc*tempcc)*factor

                vxcc_vycc_m2(kc) 	= vxcc_vycc_m1(kc) 	 + (vxcc*vycc)**2*factor
                vxcc_vzcc_m2(kc)	= vxcc_vzcc_m1(kc) 	 + (vxcc*vzcc)**2*factor
                vxcc_tempcc_m2(kc) 	= vxcc_tempcc_m1(kc) + (vxcc*tempcc)**2*factor
                vycc_vzcc_m2(kc) 	= vycc_vzcc_m1(kc) 	 + (vycc*vzcc)**2*factor
                vycc_tempcc_m2(kc) 	= vycc_tempcc_m1(kc) + (vycc*tempcc)**2*factor
                vzcc_tempcc_m2(kc) 	= vzcc_tempcc_m1(kc) + (vzcc*tempcc)**2*factor

                vz_vz_vx_m1(kc) = vz_vz_vx_m1(kc) + (vzcc*vzcc*vxcc)*factor
                vy_vy_vx_m1(kc) = vy_vy_vx_m1(kc) + (vycc*vycc*vxcc)*factor
                vx_vx_vx_m1(kc) = vx_vx_vx_m1(kc) + (vxcc*vxcc*vxcc)*factor

                vz_vz_vx_m2(kc) = vz_vz_vx_m2(kc) + (vzcc*vzcc*vxcc)**2*factor
                vy_vy_vx_m2(kc) = vy_vy_vx_m2(kc) + (vycc*vycc*vxcc)**2*factor
                vx_vx_vx_m2(kc) = vx_vx_vx_m2(kc) + (vxcc*vxcc*vxcc)**2*factor

                !! Quantities based on older statroutines file

                ! RS: Mean quantities
                vx_m1(kc)   = vx_m1(kc)   +   vx(kc,jc,ic)*factor
                vy_m1(kc)   = vy_m1(kc)   +   vy(kc,jc,ic)*factor
                vz_m1(kc)   = vz_m1(kc)   +   vz(kc,jc,ic)*factor
                temp_m1(kc) = temp_m1(kc) + temp(kc,jc,ic)*factor

                ! RS: Second moments
                vx_m2(kc)   = vx_m2(kc)   +   vx(kc,jc,ic)**2*factor
                vy_m2(kc)   = vy_m2(kc)   +   vy(kc,jc,ic)**2*factor
                vz_m2(kc)   = vz_m2(kc)   +   vz(kc,jc,ic)**2*factor
                temp_m2(kc) = temp_m2(kc) + temp(kc,jc,ic)**2*factor

                ! RS: Third moments
                vx_m3(kc)   = vx_m3(kc)   +   vx(kc,jc,ic)**3*factor
                vy_m3(kc)   = vy_m3(kc)   +   vy(kc,jc,ic)**3*factor
                vz_m3(kc)   = vz_m3(kc)   +   vz(kc,jc,ic)**3*factor
                temp_m3(kc) = temp_m3(kc) + temp(kc,jc,ic)**3*factor

                ! RS: Fourth moments
                vx_m4(kc)   = vx_m4(kc)   +   vx(kc,jc,ic)**4*factor
                vy_m4(kc)   = vy_m4(kc)   +   vy(kc,jc,ic)**4*factor
                vz_m4(kc)   = vz_m4(kc)   +   vz(kc,jc,ic)**4*factor
                temp_m4(kc) = temp_m4(kc) + temp(kc,jc,ic)**4*factor

                ! RS: heat flux
                tempvx_m1(kc) = tempvx_m1(kc) + (temp(kc,jc,ic)*0.5d0*(vx(kc,jc,ic) + vx(kp,jc,ic)))*factor

            end do
        end do
    end do

    return

end subroutine CalcStats

!******************************************************************************************************

subroutine WriteStatsEnd
    implicit none

    CHARACTER*50 :: fname

    fname = trim('Results/Stats/stafield_master.h5')
    call WriteStats(fname)

    return
end subroutine WriteStatsEnd

!******************************************************************************************************

subroutine WriteStatsSnap
    use param, only: time

    implicit none

    CHARACTER*50 :: fname
    CHARACTER*8 :: citime

    write(citime,"(I8.8)") nint(time)
    fname = trim('Results/Stats/stafield_master'//trim(citime)//'.h5')
    call WriteStats(fname)

    return
end subroutine WriteStatsSnap

!******************************************************************************************************

subroutine WriteStats(filnam)

    use mpih
    use param
    use stat_arrays
    use hdf5

    implicit none

    CHARACTER*50, INTENT(inout) :: filnam

    INTEGER                  :: nstatsamples_old, nstatsamples_new
    REAL                     :: tstat_old, tstat_new, tout_old, tout_new
    LOGICAL                  :: fexist
    CHARACTER*50             :: filnam_old,dsetname

    filnam     = trim(filnam)
    filnam_old = trim('Results/Stats/stafield_master.h5')

    inquire(file=filnam_old,exist=fexist)
    if (.not.fexist) then
        if(ismaster) write(6,*) 'Unable to read statistical files'
        if(ismaster) write(6,*) 'Restarting statistics from zero'
        readstats=.false.
    end if

    ! Adjustments writing the averaging time
    if (ismaster) then
        if (readstats) then
            dsetname = trim('averaging_samples')
            call HdfSerialReadIntScalar(dsetname,filnam_old,nstatsamples_old)
            dsetname = trim('averaging_time')
            call HdfSerialReadRealScalar(dsetname,filnam_old,tstat_old)
            dsetname = trim('averaging_interval')
            call HdfSerialReadRealScalar(dsetname,filnam_old,tout_old)
            nstatsamples_new = nstatsamples + nstatsamples_old
            tstat_new        = tstat        + tstat_old
            tout_new         = time-tsta    + tout_old
            if(filnam .ne. filnam_old) call HdfCreateBlankFile(filnam)
        else
            nstatsamples_new = nstatsamples
            tstat_new        = tstat
            tout_new         = time-tsta
            call HdfCreateBlankFile(filnam)
        end if

        dsetname = trim('averaging_samples')
        call HdfSerialWriteIntScalar(dsetname,filnam,nstatsamples_new)
        dsetname = trim('averaging_time')
        call HdfSerialWriteRealScalar(dsetname,filnam,tstat_new)
        dsetname = trim('averaging_interval')
        call HdfSerialWriteRealScalar(dsetname,filnam,tout_new)
        dsetname = trim('sample_interval')
        call HdfSerialWriteIntScalar(dsetname,filnam,nout)

        dsetname = trim('xc')
        call HdfSerialWriteReal1D(dsetname,filnam,xc,nx)
        dsetname = trim('yc')
        call HdfSerialWriteReal1D(dsetname,filnam,yc,ny)
        dsetname = trim('zc')
        call HdfSerialWriteReal1D(dsetname,filnam,zc,nz)

        dsetname = trim('xm')
        call HdfSerialWriteReal1D(dsetname,filnam,xm,nxm)
        dsetname = trim('ym')
        call HdfSerialWriteReal1D(dsetname,filnam,ym,nym)
        dsetname = trim('zm')
        call HdfSerialWriteReal1D(dsetname,filnam,zm,nzm)

        dsetname = trim('Rayleigh Number')
        call HdfSerialWriteRealScalar(dsetname,filnam,ray)

    end if

    !! Quantities based on linear interpolation at cell center

    dsetname = trim('vxcc_m1')
    call StatReadReduceWrite(vxcc_m1,filnam_old,filnam,dsetname)
    dsetname = trim('vycc_m1')
    call StatReadReduceWrite(vycc_m1,filnam_old,filnam,dsetname)
    dsetname = trim('vzcc_m1')
    call StatReadReduceWrite(vzcc_m1,filnam_old,filnam,dsetname)

    dsetname = trim('vxcc_m2')
    call StatReadReduceWrite(vxcc_m2,filnam_old,filnam,dsetname)
    dsetname = trim('vycc_m2')
    call StatReadReduceWrite(vycc_m2,filnam_old,filnam,dsetname)
    dsetname = trim('vzcc_m2')
    call StatReadReduceWrite(vzcc_m2,filnam_old,filnam,dsetname)

    dsetname = trim('vxcc_vycc_m1')
    call StatReadReduceWrite(vxcc_vycc_m1,filnam_old,filnam,dsetname)
    dsetname = trim('vxcc_vzcc_m1')
    call StatReadReduceWrite(vxcc_vzcc_m1,filnam_old,filnam,dsetname)
    dsetname = trim('vxcc_tempcc_m1')
    call StatReadReduceWrite(vxcc_tempcc_m1,filnam_old,filnam,dsetname)
    dsetname = trim('vycc_vzcc_m1')
    call StatReadReduceWrite(vycc_vzcc_m1,filnam_old,filnam,dsetname)
    dsetname = trim('vycc_tempcc_m1')
    call StatReadReduceWrite(vycc_tempcc_m1,filnam_old,filnam,dsetname)
    dsetname = trim('vzcc_tempcc_m1')
    call StatReadReduceWrite(vzcc_tempcc_m1,filnam_old,filnam,dsetname)

    dsetname = trim('vxcc_vycc_m2')
    call StatReadReduceWrite(vxcc_vycc_m2,filnam_old,filnam,dsetname)
    dsetname = trim('vxcc_vzcc_m2')
    call StatReadReduceWrite(vxcc_vzcc_m2,filnam_old,filnam,dsetname)
    dsetname = trim('vxcc_tempcc_m2')
    call StatReadReduceWrite(vxcc_tempcc_m2,filnam_old,filnam,dsetname)
    dsetname = trim('vycc_vzcc_m2')
    call StatReadReduceWrite(vycc_vzcc_m2,filnam_old,filnam,dsetname)
    dsetname = trim('vycc_tempcc_m2')
    call StatReadReduceWrite(vycc_tempcc_m2,filnam_old,filnam,dsetname)
    dsetname = trim('vzcc_tempcc_m2')
    call StatReadReduceWrite(vzcc_tempcc_m2,filnam_old,filnam,dsetname)

    dsetname = trim('vz_vz_vx_m1')
    call StatReadReduceWrite(vz_vz_vx_m1,filnam_old,filnam,dsetname)
    dsetname = trim('vy_vy_vx_m1')
    call StatReadReduceWrite(vy_vy_vx_m1,filnam_old,filnam,dsetname)
    dsetname = trim('vx_vx_vx_m1')
    call StatReadReduceWrite(vx_vx_vx_m1,filnam_old,filnam,dsetname)

    dsetname = trim('vz_vz_vx_m2')
    call StatReadReduceWrite(vz_vz_vx_m2,filnam_old,filnam,dsetname)
    dsetname = trim('vy_vy_vx_m2')
    call StatReadReduceWrite(vy_vy_vx_m2,filnam_old,filnam,dsetname)
    dsetname = trim('vx_vx_vx_m2')
    call StatReadReduceWrite(vx_vx_vx_m2,filnam_old,filnam,dsetname)

    !! Quantities based on older statroutines file

    dsetname = trim('vx_m1')
    call StatReadReduceWrite(vx_m1,filnam_old,filnam,dsetname)
    dsetname = trim('vy_m1')
    call StatReadReduceWrite(vy_m1,filnam_old,filnam,dsetname)
    dsetname = trim('vz_m1')
    call StatReadReduceWrite(vz_m1,filnam_old,filnam,dsetname)

    dsetname = trim('vx_m2')
    call StatReadReduceWrite(vx_m2,filnam_old,filnam,dsetname)
    dsetname = trim('vy_m2')
    call StatReadReduceWrite(vy_m2,filnam_old,filnam,dsetname)
    dsetname = trim('vz_m2')
    call StatReadReduceWrite(vz_m2,filnam_old,filnam,dsetname)

    dsetname = trim('vx_m3')
    call StatReadReduceWrite(vx_m3,filnam_old,filnam,dsetname)
    dsetname = trim('vy_m3')
    call StatReadReduceWrite(vy_m3,filnam_old,filnam,dsetname)
    dsetname = trim('vz_m3')
    call StatReadReduceWrite(vz_m3,filnam_old,filnam,dsetname)

    dsetname = trim('vx_m4')
    call StatReadReduceWrite(vx_m4,filnam_old,filnam,dsetname)
    dsetname = trim('vy_m4')
    call StatReadReduceWrite(vy_m4,filnam_old,filnam,dsetname)
    dsetname = trim('vz_m4')
    call StatReadReduceWrite(vz_m4,filnam_old,filnam,dsetname)

    dsetname = trim('temp_m1')
    call StatReadReduceWrite(temp_m1,filnam_old,filnam,dsetname)
    dsetname = trim('temp_m2')
    call StatReadReduceWrite(temp_m2,filnam_old,filnam,dsetname)
    dsetname = trim('temp_m3')
    call StatReadReduceWrite(temp_m3,filnam_old,filnam,dsetname)
    dsetname = trim('temp_m4')
    call StatReadReduceWrite(temp_m4,filnam_old,filnam,dsetname)

    dsetname = trim('tempvx_m1')
    call StatReadReduceWrite(tempvx_m1,filnam_old,filnam,dsetname)

    if(disscal) then
        dsetname = trim('dissth')
        call StatReadReduceWrite(dissth,filnam_old,filnam,dsetname)
        dsetname = trim('disste')
        call StatReadReduceWrite(disste,filnam_old,filnam,dsetname)
    end if

    return

end subroutine WriteStats

!******************************************************************************************************

subroutine InitInstantaneousProfile

    use param
    use hdf5

    implicit none

    character*50                :: filename,dsetname
    integer(HID_T)              :: file_id,dset_id
    integer(HID_T)              :: filespace
    integer(HSIZE_T)            :: dims(1)
    logical                     :: fexist,lexist
    character*10, dimension(21) :: statgroups
    integer                     :: hdf_error,i,ngroups

    statgroups = [character(len=10) :: &
    'time',&
    'deltatime',&
    'vx_m1', &
    'vx_m2', &
    'vx_m3', &
    'vx_m4', &
    'vy_m1', &
    'vy_m2', &
    'vy_m3', &
    'vy_m4', &
    'vz_m1', &
    'vz_m2', &
    'vz_m3', &
    'vz_m4', &
    'temp_m1', &
    'temp_m2', &
    'temp_m3', &
    'temp_m4', &
    'tempvx_m1', &
    'disste', &
    'dissth' &
    ]

    if (disscal) then 
        ngroups = 21
    else
        ngroups = 19
    end if

    filename = trim('Results/instprof_master.h5')
    
    if (ismaster) then
        inquire(file=filename,exist=fexist)
        if (.not.fexist) then
            call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,hdf_error)
            do i = 1,ngroups
                dsetname = trim(statgroups(i))
                call h5lexists_f(file_id,dsetname,lexist,hdf_error)
                if (.not.lexist) then
                    call h5gcreate_f(file_id,dsetname,dset_id,hdf_error)
                    call h5gclose_f(dset_id,hdf_error)
                end if
            end do
            
            dims(1) = 1
            dsetname = trim('samples')
            call h5screate_simple_f(1,dims,filespace,hdf_error)
            call h5dcreate_f(file_id,dsetname,H5T_NATIVE_INTEGER,filespace,dset_id,hdf_error)
            call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,0,dims,hdf_error)
            call h5dclose_f(dset_id,hdf_error)
            call h5sclose_f(filespace,hdf_error)

            call h5fclose_f(file_id,hdf_error)
        
        end if
    end if

end subroutine InitInstantaneousProfile

subroutine WriteInstantaneousProfile

    use param
    use local_arrays, only: vz,vy,vx,temp
    use decomp_2d, only: xstart,xend
    use stat_arrays
    use mpih
    
    implicit none

    integer         :: ic,jc,kc
    integer         :: ip,jp,kp
    integer         :: im,jm,km
    
    integer         :: isamples

    real            :: factor
    real            :: ivx_m1(nxm),ivx_m2(nxm),ivx_m3(nxm),ivx_m4(nxm)
    real            :: ivy_m1(nxm),ivy_m2(nxm),ivy_m3(nxm),ivy_m4(nxm)
    real            :: ivz_m1(nxm),ivz_m2(nxm),ivz_m3(nxm),ivz_m4(nxm)
    real            :: itemp_m1(nxm),itemp_m2(nxm),itemp_m3(nxm),itemp_m4(nxm)
    real            :: itempvx_m1(nxm)
    real            :: res_dummy(nxm)

    character*50    :: filename,dsetname
    character*8     :: citime

    factor          = 1.0d0/(nym*nzm)

    ivx_m1(:)       = 0.0d0
    ivx_m2(:)       = 0.0d0
    ivx_m3(:)       = 0.0d0
    ivx_m4(:)       = 0.0d0
    
    ivy_m1(:)       = 0.0d0
    ivy_m2(:)       = 0.0d0
    ivy_m3(:)       = 0.0d0
    ivy_m4(:)       = 0.0d0

    ivz_m1(:)       = 0.0d0
    ivz_m2(:)       = 0.0d0
    ivz_m3(:)       = 0.0d0
    ivz_m4(:)       = 0.0d0

    itemp_m1(:)     = 0.0d0
    itemp_m2(:)     = 0.0d0
    itemp_m3(:)     = 0.0d0
    itemp_m4(:)     = 0.0d0

    itempvx_m1(:)   = 0.0d0

    do ic=xstart(3),xend(3)

        ip = ic+1
        im = ic-1

        do jc=xstart(2),xend(2)

            jp = jc+1
            jm = jc-1

            do kc=1,nxm

                kp = kc+1
                km = kc-1

                ivx_m1(kc)     = ivx_m1(kc)     + (vx(kc,jc,ic)**1)*factor
                ivx_m2(kc)     = ivx_m2(kc)     + (vx(kc,jc,ic)**2)*factor
                ivx_m3(kc)     = ivx_m3(kc)     + (vx(kc,jc,ic)**3)*factor
                ivx_m4(kc)     = ivx_m4(kc)     + (vx(kc,jc,ic)**4)*factor
                
                ivy_m1(kc)     = ivy_m1(kc)     + (vy(kc,jc,ic)**1)*factor
                ivy_m2(kc)     = ivy_m2(kc)     + (vy(kc,jc,ic)**2)*factor
                ivy_m3(kc)     = ivy_m3(kc)     + (vy(kc,jc,ic)**3)*factor
                ivy_m4(kc)     = ivy_m4(kc)     + (vy(kc,jc,ic)**4)*factor
                
                ivz_m1(kc)     = ivz_m1(kc)     + (vz(kc,jc,ic)**1)*factor
                ivz_m2(kc)     = ivz_m2(kc)     + (vz(kc,jc,ic)**2)*factor
                ivz_m3(kc)     = ivz_m3(kc)     + (vz(kc,jc,ic)**3)*factor
                ivz_m4(kc)     = ivz_m4(kc)     + (vz(kc,jc,ic)**4)*factor
                
                itemp_m1(kc)   = itemp_m1(kc)   + (temp(kc,jc,ic)**1)*factor
                itemp_m2(kc)   = itemp_m2(kc)   + (temp(kc,jc,ic)**2)*factor
                itemp_m3(kc)   = itemp_m3(kc)   + (temp(kc,jc,ic)**3)*factor
                itemp_m4(kc)   = itemp_m4(kc)   + (temp(kc,jc,ic)**4)*factor
                
                itempvx_m1(kc) = itempvx_m1(kc) + (temp(kc,jc,ic)*0.5d0*(vx(kc,jc,ic) + vx(kp,jc,ic)))*factor

            end do
        end do
    end do

    call MpiSumReal1D(ivx_m1,res_dummy,nxm)
    ivx_m1(:) = res_dummy(:)
    call MpiSumReal1D(ivx_m2,res_dummy,nxm)
    ivx_m2(:) = res_dummy(:)
    call MpiSumReal1D(ivx_m3,res_dummy,nxm)
    ivx_m3(:) = res_dummy(:)
    call MpiSumReal1D(ivx_m4,res_dummy,nxm)
    ivx_m4(:) = res_dummy(:)

    call MpiSumReal1D(ivy_m1,res_dummy,nxm)
    ivy_m1(:) = res_dummy(:)
    call MpiSumReal1D(ivy_m2,res_dummy,nxm)
    ivy_m2(:) = res_dummy(:)
    call MpiSumReal1D(ivy_m3,res_dummy,nxm)
    ivy_m3(:) = res_dummy(:)
    call MpiSumReal1D(ivy_m4,res_dummy,nxm)
    ivy_m4(:) = res_dummy(:)

    call MpiSumReal1D(ivz_m1,res_dummy,nxm)
    ivz_m1(:) = res_dummy(:)
    call MpiSumReal1D(ivz_m2,res_dummy,nxm)
    ivz_m2(:) = res_dummy(:)
    call MpiSumReal1D(ivz_m3,res_dummy,nxm)
    ivz_m3(:) = res_dummy(:)
    call MpiSumReal1D(ivz_m4,res_dummy,nxm)
    ivz_m4(:) = res_dummy(:)

    call MpiSumReal1D(itemp_m1,res_dummy,nxm)
    itemp_m1(:) = res_dummy(:)
    call MpiSumReal1D(itemp_m2,res_dummy,nxm)
    itemp_m2(:) = res_dummy(:)
    call MpiSumReal1D(itemp_m3,res_dummy,nxm)
    itemp_m3(:) = res_dummy(:)
    call MpiSumReal1D(itemp_m4,res_dummy,nxm)
    itemp_m4(:) = res_dummy(:)

    call MpiSumReal1D(itempvx_m1,res_dummy,nxm)
    itempvx_m1(:) = res_dummy(:)

    if (ismaster) then

        filename = trim('Results/instprof_master.h5')
        dsetname = trim('samples')
        call HdfSerialReadIntScalar(dsetname,filename,isamples)
        isamples = isamples + 1
        write (citime,"(i8.8)") isamples
        call HdfSerialWriteIntScalar(dsetname,filename,isamples)

        dsetname = trim("time/"//citime)
        call HdfSerialWriteRealScalar(dsetname,filename,time)

        dsetname = trim("deltatime/"//citime)
        call HdfSerialWriteRealScalar(dsetname,filename,dt)

        dsetname = trim("vx_m1/"//citime)
        call HdfSerialWriteReal1D(dsetname,filename,ivx_m1,nxm)
        dsetname = trim("vx_m2/"//citime)
        call HdfSerialWriteReal1D(dsetname,filename,ivx_m2,nxm)
        dsetname = trim("vx_m3/"//citime)
        call HdfSerialWriteReal1D(dsetname,filename,ivx_m3,nxm)
        dsetname = trim("vx_m4/"//citime)
        call HdfSerialWriteReal1D(dsetname,filename,ivx_m4,nxm)

        dsetname = trim("vy_m1/"//citime)
        call HdfSerialWriteReal1D(dsetname,filename,ivy_m1,nxm)
        dsetname = trim("vy_m2/"//citime)
        call HdfSerialWriteReal1D(dsetname,filename,ivy_m2,nxm)
        dsetname = trim("vy_m3/"//citime)
        call HdfSerialWriteReal1D(dsetname,filename,ivy_m3,nxm)
        dsetname = trim("vy_m4/"//citime)
        call HdfSerialWriteReal1D(dsetname,filename,ivy_m4,nxm)

        dsetname = trim("vz_m1/"//citime)
        call HdfSerialWriteReal1D(dsetname,filename,ivz_m1,nxm)
        dsetname = trim("vz_m2/"//citime)
        call HdfSerialWriteReal1D(dsetname,filename,ivz_m2,nxm)
        dsetname = trim("vz_m3/"//citime)
        call HdfSerialWriteReal1D(dsetname,filename,ivz_m3,nxm)
        dsetname = trim("vz_m4/"//citime)
        call HdfSerialWriteReal1D(dsetname,filename,ivz_m4,nxm)

        dsetname = trim("temp_m1/"//citime)
        call HdfSerialWriteReal1D(dsetname,filename,itemp_m1,nxm)
        dsetname = trim("temp_m2/"//citime)
        call HdfSerialWriteReal1D(dsetname,filename,itemp_m2,nxm)
        dsetname = trim("temp_m3/"//citime)
        call HdfSerialWriteReal1D(dsetname,filename,itemp_m3,nxm)
        dsetname = trim("temp_m4/"//citime)
        call HdfSerialWriteReal1D(dsetname,filename,itemp_m4,nxm)

        dsetname = trim("tempvx_m1/"//citime)
        call HdfSerialWriteReal1D(dsetname,filename,itempvx_m1,nxm)

    end if

end subroutine WriteInstantaneousProfile