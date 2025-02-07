!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: CalcFluxes.F90                                 !
!    CONTAINS: subroutine CalcFluxes                      !
!                                                         !
!    PURPOSE: Calculate the instantaneous volume          !
!    averaged velocities                                  !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CalcFluxes

    use param
    use local_arrays,only: vx,vy,vz,temp
    use decomp_2d, only: xstart,xend,nrank
    use stat_arrays
    use mpih

    implicit none

    integer :: i,j,k
    real    :: vx_vol_avg,vy_vol_avg,vz_vol_avg,temp_vol_avg
    real    :: factor
    real    :: res_dummy
    logical :: fexist

    factor = 1.d0/(dx*real(nym*nzm))

    vx_vol_avg   = 0.d0
    vy_vol_avg   = 0.d0
    vz_vol_avg   = 0.d0
    temp_vol_avg = 0.d0

    do i = xstart(3),xend(3)
        do j = xstart(2),xend(2)
            do k = 1,nxm

                temp_vol_avg = temp_vol_avg + (factor*dx3c(k) * temp(k,j,i))
                vx_vol_avg   = vx_vol_avg   + (factor*dx3c(k) * 0.5d0*(vx(k,j,i)  +vx(k+1,j,i)  ))
                vy_vol_avg   = vy_vol_avg   + (factor*dx3c(k) * vy(k,j,i))
                vz_vol_avg   = vz_vol_avg   + (factor*dx3c(k) * vz(k,j,i))

            end do
        end do
    end do

    call MpiSumRealScalar(vx_vol_avg,res_dummy)
    vx_vol_avg = res_dummy
    call MpiSumRealScalar(vy_vol_avg,res_dummy)
    vy_vol_avg = res_dummy
    call MpiSumRealScalar(vz_vol_avg,res_dummy)
    vz_vol_avg = res_dummy
    call MpiSumRealScalar(temp_vol_avg,res_dummy)
    temp_vol_avg = res_dummy

    if (ismaster) then

        inquire(file='Results/vol_avg.out',exist=fexist)
        open(99,file='Results/vol_avg.out',status='unknown',access='sequential',position='append')
        if (.not.fexist) then
            write(99,'(5(A16,X))') 'Time', 'Vx_vol_avg', 'Vy_vol_avg', 'Vz_vol_avg', 'Temp_vol_avg'
        end if
        write(99,'(5(E16.8,X))') time,vx_vol_avg,vy_vol_avg,vz_vol_avg,temp_vol_avg
        close(99)

    end if

end subroutine CalcFluxes