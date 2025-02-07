!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: DeallocateVariables.F90                        !
!    CONTAINS: subroutine DeallocateVariables             !
!                                                         !
!    PURPOSE: Finalization routine. Deallocates all       !
!     variables used in the code                          !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine DeallocateVariables
    
    use param
    use local_arrays
    use stat_arrays
    use coefficient_arrays
    use AuxiliaryRoutines
    
    implicit none

    call DestroyReal1DArray(xc)
    call DestroyReal1DArray(xm)

    call DestroyReal1DArray(zc)
    call DestroyReal1DArray(zm)
    call DestroyReal1DArray(ak1)
    call DestroyReal1DArray(ao)

    call DestroyReal1DArray(yc)
    call DestroyReal1DArray(ym)
    call DestroyReal1DArray(ak2)
    call DestroyReal1DArray(ap)

    call DestroyReal1DArray(g3rc)
    call DestroyReal1DArray(g3rm)

    call DestroyReal1DArray(dx3c)
    call DestroyReal1DArray(dx3m)

    call DestroyReal1DArray(ixcm)
    call DestroyReal1DArray(ixcp)

    call DestroyReal1DArray(udx3c)
    call DestroyReal1DArray(udx3m)

    call DestroyReal1DArray(ap3sk)
    call DestroyReal1DArray(ac3sk)
    call DestroyReal1DArray(am3sk)

    call DestroyReal1DArray(amphk)
    call DestroyReal1DArray(acphk)
    call DestroyReal1DArray(apphk)

    call DestroyInt1dArray(kmc)
    call DestroyInt1dArray(kpc)
    call DestroyInt1dArray(kmv)
    call DestroyInt1dArray(kpv)

    if (statcal) then

        call DestroyReal1DArray(vxcc_m1)
        call DestroyReal1DArray(vycc_m1)
        call DestroyReal1DArray(vzcc_m1)
        call DestroyReal1DArray(tempcc_m1)

        call DestroyReal1DArray(vxcc_m2)
        call DestroyReal1DArray(vycc_m2)
        call DestroyReal1DArray(vzcc_m2)
        call DestroyReal1DArray(tempcc_m2)

        call DestroyReal1DArray(vxcc_vycc_m1)
        call DestroyReal1DArray(vxcc_vzcc_m1)
        call DestroyReal1DArray(vxcc_tempcc_m1)
        call DestroyReal1DArray(vycc_vzcc_m1)
        call DestroyReal1DArray(vycc_tempcc_m1)
        call DestroyReal1DArray(vzcc_tempcc_m1)

        call DestroyReal1DArray(vxcc_vycc_m2)
        call DestroyReal1DArray(vxcc_vzcc_m2)
        call DestroyReal1DArray(vxcc_tempcc_m2)
        call DestroyReal1DArray(vycc_vzcc_m2)
        call DestroyReal1DArray(vycc_tempcc_m2)
        call DestroyReal1DArray(vzcc_tempcc_m2)

        call DestroyReal1DArray(cvx_m1)
        call DestroyReal1DArray(cvy_m1)
        call DestroyReal1DArray(cvz_m1)
        call DestroyReal1DArray(ctemp_m1)

        call DestroyReal1DArray(cvx_m2)
        call DestroyReal1DArray(cvy_m2)
        call DestroyReal1DArray(cvz_m2)
        call DestroyReal1DArray(ctemp_m2)

        call DestroyReal1DArray(vz_vz_vx_m1)
        call DestroyReal1DArray(vy_vy_vx_m1)
        call DestroyReal1DArray(vx_vx_vx_m1)

        call DestroyReal1DArray(vz_vz_vx_m2)
        call DestroyReal1DArray(vy_vy_vx_m2)
        call DestroyReal1DArray(vx_vx_vx_m2)

        call DestroyReal1DArray(vx_m1)
        call DestroyReal1DArray(vy_m1)
        call DestroyReal1DArray(vz_m1)
        call DestroyReal1DArray(temp_m1)

        call DestroyReal1DArray(vx_m2)
        call DestroyReal1DArray(vy_m2)
        call DestroyReal1DArray(vz_m2)
        call DestroyReal1DArray(temp_m2)

        call DestroyReal1DArray(vx_m3)
        call DestroyReal1DArray(vy_m3)
        call DestroyReal1DArray(vz_m3)
        call DestroyReal1DArray(temp_m3)

        call DestroyReal1DArray(vx_m4)
        call DestroyReal1DArray(vy_m4)
        call DestroyReal1DArray(vz_m4)
        call DestroyReal1DArray(temp_m4)

        call DestroyReal1DArray(tempvx_m1)

        if (disscal) then
            call DestroyReal1DArray(disste)
            call DestroyReal1DArray(dissth)
        end if

    end if

    call DestroyReal3DArray(vx)
    call DestroyReal3DArray(vy)
    call DestroyReal3DArray(vz)
    call DestroyReal3DArray(temp)
    call DestroyReal3DArray(pr)
    call DestroyReal3DArray(dphhalo)
    call DestroyReal3DArray(dten)

    call DestroyReal3DArray(dph)
    call DestroyReal3DArray(rutemp)
    call DestroyReal3DArray(hro)

    call DestroyReal3DArray(rhsx)
    call DestroyReal3DArray(rhsy)
    call DestroyReal3DArray(rhsz)

    return
    
end subroutine DeallocateVariables

