!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: InitVariables.F90                              !
!    CONTAINS: subroutine InitVariables                   !
!                                                         !
!    PURPOSE: Initialization routine. Sets to zero all    !
!     variables used in the code                          !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine InitVariables

    use param
    use local_arrays
    use stat_arrays
    use coefficient_arrays
    use decomp_2d
    use AuxiliaryRoutines

    implicit none

    !-------------------------------------------------
    ! Arrays for grid making
    !-------------------------------------------------

    call AllocateReal1DArray(xc,1,nx)
    call AllocateReal1DArray(xm,1,nx)

    call AllocateReal1DArray(zc,1,nz)
    call AllocateReal1DArray(zm,1,nz)
    call AllocateReal1DArray(ak1,1,nz)
    call AllocateReal1DArray(ao,1,nz)

    call AllocateReal1DArray(yc,1,ny)
    call AllocateReal1DArray(ym,1,ny)
    call AllocateReal1DArray(ak2,1,ny)
    call AllocateReal1DArray(ap,1,ny)
    
    call AllocateReal1DArray(g3rc,1,nx)
    call AllocateReal1DArray(g3rm,1,nx)

    call AllocateReal1DArray(dx3c,1,nx)
    call AllocateReal1DArray(dx3m,1,nx)

    call AllocateReal1DArray(ixcm,1,nx)
    call AllocateReal1DArray(ixcp,1,nx)

    call AllocateReal1DArray(udx3c,1,nx)
    call AllocateReal1DArray(udx3m,1,nx)

    call AllocateReal1DArray(ap3sk,1,nx)
    call AllocateReal1DArray(ac3sk,1,nx)
    call AllocateReal1DArray(am3sk,1,nx)

    call AllocateReal1DArray(amphk,1,nx)
    call AllocateReal1DArray(acphk,1,nx)
    call AllocateReal1DArray(apphk,1,nx)

    call AllocateInt1dArray(kmc,1,nx)
    call AllocateInt1dArray(kpc,1,nx)
    call AllocateInt1dArray(kmv,1,nx)
    call AllocateInt1dArray(kpv,1,nx)

    !-------------------------------------------------
    ! Arrays for statistics
    !-------------------------------------------------

    if (statcal) then

        call AllocateReal1DArray(vxcc_m1,1,nxm)
        call AllocateReal1DArray(vycc_m1,1,nxm)
        call AllocateReal1DArray(vzcc_m1,1,nxm)
        call AllocateReal1DArray(tempcc_m1,1,nxm)

        call AllocateReal1DArray(vxcc_m2,1,nxm)
        call AllocateReal1DArray(vycc_m2,1,nxm)
        call AllocateReal1DArray(vzcc_m2,1,nxm)
        call AllocateReal1DArray(tempcc_m2,1,nxm)

        call AllocateReal1DArray(vxcc_vycc_m1,1,nxm)
        call AllocateReal1DArray(vxcc_vzcc_m1,1,nxm)
        call AllocateReal1DArray(vxcc_tempcc_m1,1,nxm)
        call AllocateReal1DArray(vycc_vzcc_m1,1,nxm)
        call AllocateReal1DArray(vycc_tempcc_m1,1,nxm)
        call AllocateReal1DArray(vzcc_tempcc_m1,1,nxm)

        call AllocateReal1DArray(vxcc_vycc_m2,1,nxm)
        call AllocateReal1DArray(vxcc_vzcc_m2,1,nxm)
        call AllocateReal1DArray(vxcc_tempcc_m2,1,nxm)
        call AllocateReal1DArray(vycc_vzcc_m2,1,nxm)
        call AllocateReal1DArray(vycc_tempcc_m2,1,nxm)
        call AllocateReal1DArray(vzcc_tempcc_m2,1,nxm)

        call AllocateReal1DArray(vz_vz_vx_m1,1,nxm)
        call AllocateReal1DArray(vy_vy_vx_m1,1,nxm)
        call AllocateReal1DArray(vx_vx_vx_m1,1,nxm)

        call AllocateReal1DArray(vz_vz_vx_m2,1,nxm)
        call AllocateReal1DArray(vy_vy_vx_m2,1,nxm)
        call AllocateReal1DArray(vx_vx_vx_m2,1,nxm)

        !! For sanity check with SpecRoutines.F90

        call AllocateReal1DArray(cvx_m1,1,nxm)
        call AllocateReal1DArray(cvy_m1,1,nxm)
        call AllocateReal1DArray(cvz_m1,1,nxm)
        call AllocateReal1DArray(ctemp_m1,1,nxm)

        call AllocateReal1DArray(cvx_m2,1,nxm)
        call AllocateReal1DArray(cvy_m2,1,nxm)
        call AllocateReal1DArray(cvz_m2,1,nxm)
        call AllocateReal1DArray(ctemp_m2,1,nxm)

        !! Quantities based on older statroutines file

        call AllocateReal1DArray(vx_m1,1,nxm)
        call AllocateReal1DArray(vy_m1,1,nxm)
        call AllocateReal1DArray(vz_m1,1,nxm)

        call AllocateReal1DArray(vx_m2,1,nxm)
        call AllocateReal1DArray(vy_m2,1,nxm)
        call AllocateReal1DArray(vz_m2,1,nxm)

        call AllocateReal1DArray(vx_m3,1,nxm)
        call AllocateReal1DArray(vy_m3,1,nxm)
        call AllocateReal1DArray(vz_m3,1,nxm)

        call AllocateReal1DArray(vx_m4,1,nxm)
        call AllocateReal1DArray(vy_m4,1,nxm)
        call AllocateReal1DArray(vz_m4,1,nxm)

        call AllocateReal1DArray(temp_m1,1,nxm)
        call AllocateReal1DArray(temp_m2,1,nxm)
        call AllocateReal1DArray(temp_m3,1,nxm)
        call AllocateReal1DArray(temp_m4,1,nxm)

        call AllocateReal1DArray(tempvx_m1,1,nxm)

        if (disscal) then
            call AllocateReal1DArray(disste,1,nxm)
            call AllocateReal1DArray(dissth,1,nxm)
        end if
    end if

    !-------------------------------------------------
    ! Arrays with ghost cells
    !-------------------------------------------------

    call AllocateReal3DArray(vy,1,nx,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)
    call AllocateReal3DArray(vz,1,nx,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)
    call AllocateReal3DArray(vx,1,nx,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)
    call AllocateReal3DArray(temp,1,nx,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)
    call AllocateReal3DArray(pr,1,nxm,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)
    call AllocateReal3DArray(dphhalo,1,nxm,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)
    call AllocateReal3DArray(dten,1,nx,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)

    !-----------------------------------------------
    ! Arrays without ghost cells
    !-----------------------------------------------
    
    call AllocateReal3DArray(dph,1,nxm,xstart(2),xend(2),xstart(3),xend(3))
    call AllocateReal3DArray(hro,1,nxm,xstart(2),xend(2),xstart(3),xend(3))
    call AllocateReal3DArray(rutemp,1,nxm,xstart(2),xend(2),xstart(3),xend(3))
    
    call AllocateReal3DArray(rhsx,1,nxm,xstart(2),xend(2),xstart(3),xend(3))
    call AllocateReal3DArray(rhsy,ystart(1),yend(1),1,nym,ystart(3),yend(3))
    call AllocateReal3DArray(rhsz,zstart(1),zend(1),zstart(2),zend(2),1,nzm)


    return

end subroutine InitVariables
