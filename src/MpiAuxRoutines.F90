!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: HdfRoutines.F90                                !
!    CONTAINS: subroutines MPI*                           !
!                                                         !
!    PURPOSE: Wrappers for MPI Routines                   !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine MpiBcastInt(n)
    use mpih
    implicit none
    integer, intent(in) :: n
    call MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    return
end subroutine MpiBcastInt
!==============================================================================

subroutine MpiBcastReal(n)
    use mpih
    implicit none
    real, intent(in) :: n
    call MPI_BCAST(n,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    return
end subroutine MpiBcastReal
!==============================================================================

subroutine MpiBarrier
    use mpih
    implicit none
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    return
end subroutine MpiBarrier
!==============================================================================

subroutine MpiSumRealScalar(var,res)
    use mpih
    implicit none
    real, intent(in)  :: var
    real, intent(out) :: res
    call MPI_REDUCE(var,res,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    return
end subroutine MpiSumRealScalar
!==============================================================================

subroutine MpiMaxRealScalar(var,res)
    use mpih
    implicit none
    real, intent(in)  :: var
    real, intent(out) :: res
    call MPI_REDUCE(var,res,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
    return
end subroutine MpiMaxRealScalar
!==============================================================================

subroutine MpiAllSumRealScalar(var,res)
    use mpih
    implicit none
    real, intent(in)  :: var
    real, intent(out) :: res
    call MPI_ALLREDUCE(var,res,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    return
end subroutine MpiAllSumRealScalar
!==============================================================================

subroutine MpiAllMaxRealScalar(var,res)
    use mpih
    implicit none
    real, intent(in)  :: var
    real, intent(out) :: res
    call MPI_ALLREDUCE(var,res,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
    return
end subroutine MpiAllMaxRealScalar
!==============================================================================

subroutine MpiMinRealScalar(var,res)
    use mpih
    implicit none
    real, intent(in)  :: var
    real, intent(out) :: res
    call MPI_REDUCE(var,res,1,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,ierr)
    return
end subroutine MpiMinRealScalar
!==============================================================================

subroutine MpiSumReal1D(var,res,sz)
    use mpih
    implicit none
    integer, intent(in) :: sz
    real, intent(in),  dimension(1:sz) :: var
    real, intent(out), dimension(1:sz) :: res
    call MPI_REDUCE(var,res,sz,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    return
end subroutine MpiSumReal1D
! ==============================================================================

subroutine MpiAbort
    use mpih
    implicit none
    call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
    return
end subroutine MpiAbort
