!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: StatReadReduceWrite.F90                        !
!    CONTAINS: subroutine StatReadReduceWrite             !
!                                                         !
!    PURPOSE: Reduce statistics across processors, read   !
!     old statistics if required and write them out       !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine StatReadReduceWrite(var,filename_old,filename,dsetname)

    use param
    use mpih

    implicit none

    character*50, intent(in) :: filename_old, filename, dsetname
    real, dimension(1:nxm)   :: var, var_sum, var_old

    call MpiSumReal1D(var,var_sum,nxm)

    if (ismaster) then
        if(readstats) then
            call HdfSerialReadReal1D(dsetname,filename_old,var_old,nxm)
            var_sum = var_sum + var_old
        endif
        call HdfSerialWriteReal1D(dsetname,filename,var_sum,nxm)
    end if

    return

end subroutine StatReadReduceWrite



