!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: SetVelocityField.F90                           !
!    CONTAINS: subroutine SetVelocityField                !
!                                                         !
!    PURPOSE: Set the velocity fields as given by         !
!    Rayleigh–Darcy equations                             !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SetVelocityField

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays, only: vx,vy,vz,temp,pr, permxc, perm, kfx, kfy, kfz

    implicit none

    integer :: k,j,i

    do i = xstart(3),xend(3)
        do j = xstart(2),xend(2)

            vx(1,j,i) = 0.0d0
            do k = 2,nxm
                vx(k,j,i) = kfx(k,j,i)*(-gravity_x*(dx3c(k-1)*temp(k,j,i) + dx3c(k)*temp(k-1,j,i))/(2.0d0*g3rc(k)) - ((pr(k,j,i)-pr(k-1,j,i))*udx3c(k)))
            end do
            vx(nx,j,i) = 0.0d0

            do k = 1,nxm
                vy(k,j,i) = kfy(k,j,i)*(-gravity_y*0.5d0*(temp(k,j,i) + temp(k,j-1,i)) - ((pr(k,j,i)-pr(k,j-1,i))*dy))


                vz(k,j,i) = kfz(k,j,i)*(-gravity_z*0.5d0*(temp(k,j,i) + temp(k,j,i-1)) - ((pr(k,j,i)-pr(k,j,i-1))*dz))
            end do

        end do
    end do

end subroutine SetVelocityField
