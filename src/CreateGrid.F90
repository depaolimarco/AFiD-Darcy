!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: CreateGrid.F90                                 !
!    CONTAINS: subroutine CreateGrid                      !
!                                                         !
!    PURPOSE: Compute the indices, grid, grid metrics     !
!     and coefficients for differentiation                !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! FOR UNIFORM CLUSTERING
subroutine GridType0

    use param
    implicit none
    integer :: kc

    do kc = 1,nx
        xc(kc) = real(kc - 1)*alx3/real(nxm)
    end do

end subroutine GridType0
    
! FOR HYPERBOLIC TANGENT CLUSTERING
subroutine GridType4

    use param
    implicit none
    integer :: kc

    do kc = 1,nx
        xc(kc) = (1 + tanh(str3*(real(2*kc-nx-1)/real(nxm)))/tanh(str3))*0.5*alx3
        if (xc(kc) .lt. 0 .or. xc(kc) .gt. alx3) then
            write (*, *) 'Grid is too streched: ', 'zc(', kc, ')=', xc(kc)
            stop
        end if
    end do

end subroutine GridType4

! FOR CLIPPED CHEBYSHEV CLUSTERING
subroutine GridType6

    use param
    use AuxiliaryRoutines

    implicit none

    integer :: kc
    integer :: nxmo, nclip
    real    :: delet,etain
    real, allocatable, dimension(:) :: etaz, etazm

    call AllocateReal1DArray(etaz,1,nx+500)
    call AllocateReal1DArray(etazm,1,nx+500)

    nclip = int(str3)
    nxmo = nx + nclip + nclip
    do kc = 1, nxmo
        etazm(kc) = +cos(pi*(float(kc) - 0.5)/float(nxmo))
    end do
    do kc = 1, nx
        etaz(kc) = etazm(kc + nclip)
    end do
    delet = etaz(1) - etaz(nx)
    etain = etaz(1)
    do kc = 1, nx
        etaz(kc) = etaz(kc)/(0.5*delet)
    end do
    xc(1) = 0.
    do kc = 2, nxm
        xc(kc) = alx3*(1.-etaz(kc))*0.5
    end do
    xc(nx) = alx3

    call DestroyReal1DArray(etaz)
    call DestroyReal1DArray(etazm)

end subroutine GridType6

! FOR ERROR FUNCTION TYPE CLUSTERING
subroutine GridType7

    use param
    implicit none
    integer :: kc
    real    :: x_temporary

    do kc = 1, nx
        x_temporary = (real(kc - 1)/real(nxm)) - 0.5
        x_temporary = erf(str3*x_temporary)/erf(str3/2)
        xc(kc) = 0.5*(1 + x_temporary)*alx3
    end do

end subroutine GridType7

! TO READ IN GRID FROM INPUT FILE
subroutine GridType8

    use param
    implicit none
    integer :: kc
    logical :: fexist

    inquire (file='./wall_normal_grid.in', exist=fexist)
    if (fexist) then
        if (ismaster) write (6, *) 'Reading custom grid from wall_normal_grid.in'
    else
        if (ismaster) write (6, *) 'Warning: wall_normal_grid.in not found!'
        call MpiAbort
    end if
    open (unit=78, file='wall_normal_grid.in', status='old')
    do kc = 1, nx
        read (78, *) xc(kc)
    end do
    close (78)

end subroutine GridType8

!=============================================================================!

subroutine CreateGrid

    use param
    use AuxiliaryRoutines

    implicit none

    real    :: x1, x2, x3
    real    :: a33, a33m, a33p
    real    :: delet, etain, tstr3
    real    :: z2dp
    integer :: ic, jc, kc, km, kp

    ! INDICES IN VERTICAL DIRECTION

    do kc = 1, nxm
        kmv(kc) = kc - 1
        kpv(kc) = kc + 1
        if (kc .eq. 1) kmv(kc) = kc
        if (kc .eq. nxm) kpv(kc) = kc
    end do

    do kc = 1, nxm
        kpc(kc) = kpv(kc) - kc
        kmc(kc) = kc - kmv(kc)
    end do

    ! UNIFORM (HORIZONTAL DIRECTIONS) GRID

    do ic = 1, nz
        zc(ic) = real(ic - 1)*zlen/real(nzm)
    end do

    do ic = 1, nzm
        zm(ic) = (zc(ic) + zc(ic + 1))*0.5d0
    end do

    do jc = 1, ny
        yc(jc) = real(jc - 1)*ylen/real(nym)
    end do

    do jc = 1, nym
        ym(jc) = (yc(jc) + yc(jc + 1))*0.5d0
    end do

    ! VERTICAL COORDINATE DEFINITION

    if (istr3 .eq. 0) call GridType0
    if (istr3 .eq. 4) call GridType4
    if (istr3 .eq. 6) call GridType6
    if (istr3 .eq. 7) call GridType7
    if (istr3 .eq. 8) call GridType8

    ! METRIC FOR UNIFORM DIRECTIONS

    dx = real(nxm)/alx3
    dy = real(nym)/ylen
    dz = real(nzm)/zlen

    dxq = dx*dx
    dyq = dy*dy
    dzq = dz*dz

    ! STAGGERED COORDINATES AND
    ! METRIC QUANTITIES FOR NON-UNIFORM
    ! DIRECTIONS

    do kc = 1, nxm
        xm(kc) = (xc(kc) + xc(kc + 1))*0.5d0
        dx3c(kc) = (xc(kc + 1) - xc(kc))*dx
    end do

    do kc = 1, nxm - 1
        dx3m(kc) = (xm(kc + 1) - xm(kc))*dx
    end do

    g3rc(1) = (xc(2) - xc(1))*dx
    do kc = 2, nxm
        g3rc(kc) = (xc(kc + 1) - xc(kc - 1))*dx*0.5d0
    end do
    g3rc(nx) = (xc(nx) - xc(nxm))*dx

    g3rm(1) = (xm(2) - xc(1))*dx*0.5d0
    do kc = 2, nxm - 1
        g3rm(kc) = (xm(kc + 1) - xm(kc - 1))*dx*0.5d0
    end do
    g3rm(nxm) = (xc(nx) - xm(nxm - 1))*dx*0.5d0

    ixcm(1) = 0.5d0
    ixcp(1) = 0.5d0
    do kc = 2, nxm
        ixcm(kc) = (xc(kc) - xc(kc - 1))/(xc(kc + 1) - xc(kc - 1))
        ixcp(kc) = (xc(kc + 1) - xc(kc))/(xc(kc + 1) - xc(kc - 1))
    end do 
    ixcm(nx) = 0.5d0
    ixcp(nx) = 0.5d0

    do kc = 1, nxm
        udx3m(kc) = dx/dx3c(kc)
    end do

    do kc = 1, nx
        udx3c(kc) = dx/g3rc(kc)
    end do

    ! WRITE GRID INFORMATION

    if (ismaster) then
        
        open (unit=78, file='Results/axicor.out', status='unknown')
        write (78, '(A8,2X,4(A24,2X))') 'kc', 'xc(kc)', 'xm(kc)', 'g3rc(kc)', 'dx3c(kc)'
        do kc = 1, nx
            write (78, '(I8,2X,4(E24.16,2X))') kc, xc(kc), xm(kc), g3rc(kc), dx3c(kc)
        end do
        close (78)

        ! QUANTITIES FOR DERIVATIVES
        
        open (unit=78, file='Results/fact3.out', status='unknown')
        write (78, '(A8,2X,2(A24,2X))') 'kc', 'udx3m(kc)', 'udx3c(kc)'
        do kc = 1, nxm
            write (78, '(I8,2X,2(E24.16,2X))') kc, udx3m(kc), udx3c(kc)
        end do
        write (78, '(I8,2X,2(E24.16,2X))') nx, udx3m(nxm), udx3c(nx)
        close (78)

    end if

    ! COEFFICIENTS FOR D2_/DZ2

    am1si = dzq
    ac1si = -2.0d0*dzq
    ap1si = dzq

    ! COEFFICIENTS FOR D2_/DY2

    am2sj = dyq
    ac2sj = -2.0d0*dyq
    ap2sj = dyq

    ! COEFFICIENTS FOR D2_/DX2 (XM VARIABLE)

    ap3sk(1) = dxq/g3rm(1)/dx3m(1)
    am3sk(1) = dxq/g3rm(1)/g3rc(1)
    ac3sk(1) = -(ap3sk(1) + am3sk(1))

    do kc = 2, nxm - 1
        kp = kc + 1
        km = kc - 1
        ap3sk(kc) = dxq/g3rm(kc)/dx3m(kc)
        am3sk(kc) = dxq/g3rm(kc)/dx3m(km)
        ac3sk(kc) = -(ap3sk(kc) + am3sk(kc))
    end do

    ap3sk(nxm) = dxq/g3rm(nxm)/g3rc(nx)
    am3sk(nxm) = dxq/g3rm(nxm)/dx3m(nxm - 1)
    ac3sk(nxm) = -(ap3sk(nxm) + am3sk(nxm))

    return

end subroutine CreateGrid

