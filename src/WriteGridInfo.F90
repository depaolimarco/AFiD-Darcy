!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: WriteGridInfo.F90                              !
!    CONTAINS: subroutine WriteGridInfo                   !
!                                                         !
!    PURPOSE: Write the grid information in               !
!     cordin_info.h5                                      !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine WriteGridInfo

    use mpih
    use param
    use hdf5

    implicit none

    character*70 :: namfile
    character*50 :: dsetname

    if (ismaster) then

        namfile='Results/cordin_info.h5'
        call HdfCreateBlankFile(namfile)

        dsetname = trim('xc')
        call HdfSerialWriteReal1D(dsetname,namfile,xc,nx)
        dsetname = trim('yc')
        call HdfSerialWriteReal1D(dsetname,namfile,yc,ny)
        dsetname = trim('zc')
        call HdfSerialWriteReal1D(dsetname,namfile,zc,nz)

        dsetname = trim('xm')
        call HdfSerialWriteReal1D(dsetname,namfile,xm,nxm)
        dsetname = trim('ym')
        call HdfSerialWriteReal1D(dsetname,namfile,ym,nym)
        dsetname = trim('zm')
        call HdfSerialWriteReal1D(dsetname,namfile,zm,nzm)

    endif

    return
    
end subroutine WriteGridInfo