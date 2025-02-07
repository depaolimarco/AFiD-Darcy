!***********************************************************
! Declaration of global variables
!***********************************************************

module param

    implicit none

    !==========================================================
    !       read from input file bou.in
    !==========================================================
    ! 1XX -----------------------
    integer     :: nread 
    integer     :: starea
    integer     :: ireset
    ! 2XX -----------------------
    integer     :: nxm, nym, nzm
    real        :: alx3,ylen,zlen
    integer     :: istr3
    real        :: str3
    ! 3XX -----------------------
    integer     :: ntst
    real        :: tmax
    real        :: walltimemax
    ! 4XX -----------------------
    integer     :: nsst 
    integer     :: idtv
    real        :: dt
    real        :: dtmin,dtmax
    real        :: limitCFL
    real        :: limitDSP
    real        :: resid
    real        :: limitVel
    real        :: eps
    real        :: initime
    integer     :: flagimp
    ! 5XX -----------------------
    real        :: ray
    integer     :: inslws
    integer     :: inslwn
    real        :: tempbp
    real        :: temptp
    real        :: gravity_x
    real        :: gravity_y
    real        :: gravity_z
    ! 6XX -----------------------
    integer     :: nout
    integer     :: flagstat
    integer     :: flaginst
    integer     :: flagbal
    real        :: tsta
    real        :: startsnap
    real        :: savesnap
    integer     :: stst3flag
    integer     :: movieflag
    real        :: tframe

    !=================================================
    !       end of input file
    !=================================================

    !******* Grid parameters **********************************
    integer                             :: nx, ny, nz
    real                                :: dx,dy,dz,dxq,dyq,dzq
    real, allocatable, dimension(:)     :: xc,xm
    real, allocatable, dimension(:)     :: yc,ym
    real, allocatable, dimension(:)     :: zc,zm
    real, allocatable, dimension(:)     :: g3rc,g3rm
    real, allocatable, dimension(:)     :: dx3c,dx3m
    real, allocatable, dimension(:)     :: ixcm,ixcp

    !******* Quantities for derivatives ***********************
    real, allocatable, dimension(:)     :: udx3c,udx3m

    !******* Grid indices *************************************
    integer, allocatable, dimension(:)  :: kmc,kpc,kmv,kpv

    !******* Metric coefficients ******************************
    real, allocatable, dimension(:)     :: ap3sk,ac3sk,am3sk
    real                                :: ap2sj,ac2sj,am2sj
    real                                :: ap1si,ac1si,am1si

    !******* Variables for FFTW and Poisson solver*************
    real, allocatable, dimension(:)     :: ak2,ap
    real, allocatable, dimension(:)     :: ak1,ao
    ! real, allocatable, dimension(:)     :: ak0,an
    real, allocatable, dimension(:)     :: amphk,acphk,apphk

    !******* Other variables ***********************************
    
    real                                :: ren, pec
    real                                :: pi
    real                                :: qqmax,qqtot
    real                                :: re
    real                                :: tempmax,tempmin,tempm
    real, dimension(1:3)                :: vmax

    real                                :: time
    integer                             :: ntime,ns
    real                                :: beta
    real                                :: al,ga,ro
    real, dimension(1:3)                :: gam,rom,alm

    logical                             :: dumpslabs        = .false.
    logical                             :: statcal          = .false.
    logical                             :: disscal          = .false.
    logical                             :: readflow         = .false.
    logical                             :: readstats        = .false.
    logical                             :: ismaster         = .false.
    logical                             :: resetlogstime    = .false.
    logical                             :: variabletstep    = .true.
    logical                             :: instprof         = .false.
    logical                             :: iterative        = .false.

    integer                             :: lvlhalo          = 2

    real                                :: nusupp
    real                                :: nuslow

    real                                :: c1prev,c2prev
    real                                :: currentTol

end module param

!************* End of param module******************************
!===============================================================

!******* 3D arrays, dynamically allocated by each process*******
module local_arrays

    use param

    implicit none

    real,allocatable,dimension(:,:,:)   :: vx,vy,vz,temp
    real,allocatable,dimension(:,:,:)   :: rutemp,hro,dten
    real,allocatable,dimension(:,:,:)   :: rhsx,rhsy,rhsz
    real,allocatable,dimension(:,:,:)   :: dph,dphhalo,pr
    real,allocatable,dimension(:,:,:)   :: advc,advm

end module local_arrays
!===============================================================

!******* 3D arrays, dynamically allocated by each process*******
module coefficient_arrays

    use param

    implicit none

    real,allocatable,dimension(:,:,:)   :: fmmc,fmcm,fmcc,fmcp,fmpc
    real,allocatable,dimension(:,:,:)   :: fcmm,fcmc,fcmp
    real,allocatable,dimension(:,:,:)   :: fccm,fccc,fccp
    real,allocatable,dimension(:,:,:)   :: fcpm,fcpc,fcpp
    real,allocatable,dimension(:,:,:)   :: fpmc,fpcm,fpcc,fpcp,fppc

end module coefficient_arrays
!===============================================================

module stat_arrays
    implicit none

    real,allocatable,dimension(:)       :: vx_m1,vy_m1,vz_m1,temp_m1
    real,allocatable,dimension(:)       :: vx_m2,vy_m2,vz_m2,temp_m2
    real,allocatable,dimension(:)       :: vx_m3,vy_m3,vz_m3,temp_m3
    real,allocatable,dimension(:)       :: vx_m4,vy_m4,vz_m4,temp_m4
    real,allocatable,dimension(:)       :: disste,dissth,tempvx_m1

    real,allocatable,dimension(:)       :: vxcc_m1,vycc_m1,vzcc_m1,tempcc_m1
    real,allocatable,dimension(:)       :: vxcc_m2,vycc_m2,vzcc_m2,tempcc_m2
    real,allocatable,dimension(:)       :: vxcc_vycc_m1,vxcc_vzcc_m1,vxcc_tempcc_m1,vycc_vzcc_m1,vycc_tempcc_m1,vzcc_tempcc_m1
    real,allocatable,dimension(:)       :: vxcc_vycc_m2,vxcc_vzcc_m2,vxcc_tempcc_m2,vycc_vzcc_m2,vycc_tempcc_m2,vzcc_tempcc_m2

    real,allocatable,dimension(:)       :: cvx_m1,cvy_m1,cvz_m1,ctemp_m1
    real,allocatable,dimension(:)       :: cvx_m2,cvy_m2,cvz_m2,ctemp_m2

    real,allocatable,dimension(:)       :: vz_vz_vx_m1,vy_vy_vx_m1,vx_vx_vx_m1
    real,allocatable,dimension(:)       :: vz_vz_vx_m2,vy_vy_vx_m2,vx_vx_vx_m2

    integer                             :: nstatsamples
    real                                :: tstat

end module stat_arrays
!=====================================================

module stat3_param
    implicit none
    integer :: kslab(1:9)
    real    :: xslab(1:9)
end module stat3_param
!=====================================================

module mpih
    implicit none
    include 'mpif.h'
    integer :: ierr
    integer, parameter :: master=0
end module mpih
!====================================================

module fftw_params
    ! use param, only: m2m,m2mh,m1m
    use iso_c_binding

    type, bind(C) :: fftw_iodim
        integer(C_INT) n, is, os
    end type fftw_iodim

    interface

        type(C_PTR) function fftw_plan_guru_dft(rank,dims, &
            howmany_rank,howmany_dims,in,out,sign,flags) &
            bind(C, name='fftw_plan_guru_dft')
            import
            integer(C_INT), value :: rank
            type(fftw_iodim), dimension(*), intent(in) :: dims
            integer(C_INT), value :: howmany_rank
            type(fftw_iodim), dimension(*), intent(in) :: howmany_dims
            complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
            complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
            integer(C_INT), value :: sign
            integer(C_INT), value :: flags
        end function fftw_plan_guru_dft

        type(C_PTR) function fftw_plan_guru_dft_r2c(rank,dims, &
            howmany_rank,howmany_dims,in,out,flags) &
            bind(C, name='fftw_plan_guru_dft_r2c')
            import
            integer(C_INT), value :: rank
            type(fftw_iodim), dimension(*), intent(in) :: dims
            integer(C_INT), value :: howmany_rank
            type(fftw_iodim), dimension(*), intent(in) :: howmany_dims
            real(C_DOUBLE), dimension(*), intent(out) :: in
            complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
            integer(C_INT), value :: flags
        end function fftw_plan_guru_dft_r2c

        type(C_PTR) function fftw_plan_guru_dft_c2r(rank,dims, &
            howmany_rank,howmany_dims,in,out,flags)  &
            bind(C, name='fftw_plan_guru_dft_c2r')
            import
            integer(C_INT), value :: rank
            type(fftw_iodim), dimension(*), intent(in) :: dims
            integer(C_INT), value :: howmany_rank
            type(fftw_iodim), dimension(*), intent(in) :: howmany_dims
            complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
            real(C_DOUBLE), dimension(*), intent(out) :: out
            integer(C_INT), value :: flags
        end function fftw_plan_guru_dft_c2r

    end interface

    integer FFTW_PATIENT, FFTW_FORWARD, FFTW_BACKWARD,FFTW_ESTIMATE
    parameter (FFTW_PATIENT=32)
    parameter (FFTW_MEASURE=0)
    parameter (FFTW_ESTIMATE=64)
    parameter (FFTW_FORWARD=-1)
    parameter (FFTW_BACKWARD=1)
    type(C_PTR) :: fwd_guruplan_y,bwd_guruplan_y
    type(C_PTR) :: fwd_guruplan_z,bwd_guruplan_z
    logical :: planned=.false.

    real,allocatable,dimension(:,:,:) :: ry1,rz1
    complex,allocatable,dimension(:,:,:) :: cy1,cz1,dphc

end module fftw_params
