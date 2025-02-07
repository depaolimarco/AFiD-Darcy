!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: ReadInputFile.F90                              !
!    CONTAINS: subroutine ReadInputFile                   !
!                                                         !
!    PURPOSE: Read parameters from bou.in file            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ReadInputFile
    
    use param
    
    implicit none
    
    real                :: gravity_m
    character(len=4)    :: dummy
    logical             :: fexist

    call read_from_bouin

    nx=nxm+1
    ny=nym+1
    nz=nzm+1

    ! DEFINITIONS FOR THE NATURAL CONVECTION

    ren = dsqrt(ray)
    pec = ray
    pi  = 2.d0*dasin(1.d0)

    gravity_m = dsqrt((gravity_x*gravity_x) + (gravity_y*gravity_y) + (gravity_z*gravity_z))
    gravity_x = gravity_x/gravity_m
    gravity_y = gravity_y/gravity_m
    gravity_z = gravity_z/gravity_m

    if(flagstat.ne.0) statcal = .true.
    if(idtv.eq.0) variabletstep = .false.
    if(flagbal.ne.0) disscal = .true.
    if(flaginst.ne.0) instprof = .true.
    if(nread.ne.0) readflow = .true.
    if(ireset.ne.0) resetlogstime = .true.


    if(stst3flag.ne.0) then
        inquire(file='./stst3.in', exist=fexist)
        if(fexist) then
            dumpslabs = .true.
        else
            write(6,*) "stst3.in not found, turning off slab dump"
        end if
    endif

    if(starea.ne.0) then
        readstats = .true.
        if (.not. readflow) write(6,*) 'Warning: Restarting flowfield with statistics read'
    endif

    return

end subroutine ReadInputFile


SUBROUTINE read_from_bouin
    
    USE param, ONLY: nread, starea, ireset, nxm, nym, nzm, alx3, ylen, zlen, &
        istr3, str3, ntst, tmax, walltimemax, nsst, idtv, dt, &
        dtmin, dtmax, limitCFL, limitDSP, resid, limitVel, eps, initime, flagimp, &
        ray, &
        inslws, inslwn, tempbp, temptp, gravity_x, gravity_y, gravity_z, &
        nout, tsta, startsnap, savesnap, &
        flagstat, flaginst, flagbal, stst3flag, movieflag, tframe
   
    IMPLICIT NONE

    INTEGER i, narg
    INTEGER, PARAMETER :: LARGE_INTEGER=200, nimp=100
    CHARACTER*200 line
    CHARACTER*100 ss(nimp)
    CHARACTER*1   stringdummy1

    open(unit=15,file='bou_newStyle.in',status='old')
    DO 222 i=1,LARGE_INTEGER
        read(15,'(a200)',end=444) line

        IF(line(1:3)=='101') THEN
            ! ###### READ DATA ######
            call scan_string (line, 1, ss, narg)
            stringdummy1=ss(1)
            if('y'==stringdummy1) then
                nread=1
            elseif('n'==stringdummy1) then
                nread=0
            else
                Write(*,*)"ERROR: Input value of parameter NREAD not valid"
                Write(*,*)"       ===> valid values 'n' or 'y'            "
                call stop_config
            endif

        ELSEIF(line(1:3)=='102') THEN
            ! ###### READ STATS ######
            call scan_string (line, 1, ss, narg)
            stringdummy1=ss(1)
            if('y'==stringdummy1) then
                starea=1
            elseif('n'==stringdummy1) then
                starea=0
            else
                Write(*,*)"ERROR: Input value of parameter STAREAD not valid"
                Write(*,*)"       ===> valid values 'n' or 'y'             "
                call stop_config
            endif

        ELSEIF(line(1:3)=='103') THEN
            ! ###### RESET COUNTERS ######
            call scan_string (line, 1, ss, narg)
            stringdummy1=ss(1)
            if('y'==stringdummy1) then
                ireset=1
            elseif('n'==stringdummy1) then
                ireset=0
            else
                Write(*,*)"ERROR: Input value of parameter IRESET not valid"
                Write(*,*)"       ===> valid values 'n' or 'y'             "
                call stop_config
            endif

        ELSEIF(line(1:3)=='201') THEN
            ! ###### GRID DIMENSION ######
            call scan_string (line, 3, ss, narg)
            read(ss(1),*) nxm
            read(ss(2),*) nym
            read(ss(3),*) nzm

        ELSEIF(line(1:3)=='202') THEN
            ! ###### AXIS EXTENDS ######
            call scan_string (line, 3, ss, narg)
            read(ss(1),*) alx3
            read(ss(2),*) ylen
            read(ss(3),*) zlen

        ELSEIF(line(1:3)=='203') THEN
            ! ###### AXIS PARAMETER 1 ######
            call scan_string (line, 1, ss, narg)
            read(ss(1),*) istr3

        ELSEIF(line(1:3)=='204') THEN
            ! ###### AXIS PARAMETER 2 ######
            call scan_string (line, 1, ss, narg)
            read(ss(1),*) str3

        ELSEIF(line(1:3)=='301') THEN
            ! ###### MAX STEPS ######
            call scan_string (line, 1, ss, narg)
            read(ss(1),*) ntst

        ELSEIF(line(1:3)=='302') THEN
            ! ###### MAX TIME ######
            call scan_string (line, 1, ss, narg)
            read(ss(1),*) tmax

        ELSEIF(line(1:3)=='303') THEN
            ! ###### MAX WALLCLOCK TIME ######
            call scan_string (line, 1, ss, narg)
            read(ss(1),*) walltimemax

        ELSEIF(line(1:3)=='401') THEN
            ! ###### NUMERICAL SCHEME ######
            call scan_string (line, 1, ss, narg)
            read(ss(1),*) nsst

        ELSEIF(line(1:3)=='402') THEN
            ! ###### VARY TIME STEP ######
            call scan_string (line, 1, ss, narg)
            stringdummy1=ss(1)
            if('y'==stringdummy1) then
                idtv=1
            elseif('n'==stringdummy1) then
                idtv=0
            else
                Write(*,*)"ERROR: Input value of parameter IDTV not valid"
                Write(*,*)"       ===> valid values 'n' or 'y'             "
                call stop_config
            endif

        ELSEIF(line(1:3)=='403') THEN
            ! ###### FIRST TIME STEP ######
            call scan_string (line, 1, ss, narg)
            read(ss(1),*) dt

        ELSEIF(line(1:3)=='404') THEN
            ! ###### MIN MAX STEP SIZE ######
            call scan_string (line, 2, ss, narg)
            read(ss(1),*) dtmin
            read(ss(2),*) dtmax

        ELSEIF(line(1:3)=='405') THEN
            ! ###### TIMESTEP STABILITY ######
            call scan_string (line, 1, ss, narg)
            read(ss(1),*) limitCFL

        ELSEIF(line(1:3)=='406') THEN
            ! ###### RESID ######
            call scan_string (line, 1, ss, narg)
            read(ss(1),*) resid

        ELSEIF(line(1:3)=='407') THEN
            ! ###### MAX VELOCITY ######
            call scan_string (line, 1, ss, narg)
            read(ss(1),*) limitVel

        ELSEIF(line(1:3)=='408') THEN
            ! ###### RESID ######
            call scan_string (line, 1, ss, narg)
            read(ss(1),*) eps

        ELSEIF(line(1:3)=='409') THEN
            ! ###### NUMERICAL SCHEME ######
            call scan_string (line, 1, ss, narg)
            read(ss(1),*) initime

        ELSEIF(line(1:3)=='410') THEN
            ! ###### IMPLICIT TREATMENT OF DIFFUSIVE TERMS IN PERIODIC DIRECTIONS ######
            call scan_string (line, 1, ss, narg)
            stringdummy1=ss(1)
            if('y'==stringdummy1) then
                flagimp=1
            elseif('n'==stringdummy1) then
                flagimp=0
            else
                Write(*,*)"ERROR: Input value of parameter IMPLICIT not valid"
                Write(*,*)"       ===> valid values 'n' or 'y'            "
                call stop_config
            endif

        ELSEIF(line(1:3)=='501') THEN
            ! ###### THERM OPTION & RAYLEIGH NUMBER ######
            call scan_string (line, 1, ss, narg)
            read(ss(1),*) ray

        ELSEIF(line(1:3)=='504') THEN
            ! ###### BOUNDARY CONDITION TYPE ######
            call scan_string (line, 2, ss, narg)
            read(ss(1),*) inslws
            read(ss(2),*) inslwn

        ELSEIF(line(1:3)=='505') THEN
            ! ###### BOUNDARY CONDITION VALUE ######
            call scan_string (line, 2, ss, narg)
            read(ss(1),*) tempbp
            read(ss(2),*) temptp

        ELSEIF(line(1:3)=='506') THEN
            ! ###### DIRECTION OF GRAVITY VECTOR ######
            call scan_string (line, 3, ss, narg)
            read(ss(1),*) gravity_x
            read(ss(2),*) gravity_y
            read(ss(3),*) gravity_z

        ELSEIF(line(1:3)=='601') THEN
            ! ###### OUTPUT INTERVAL ######
            call scan_string (line, 1, ss, narg)
            read(ss(1),*) nout

        ELSEIF(line(1:3)=='602') THEN
            ! ###### STAT CALC ######
            call scan_string (line, 2, ss, narg)
            stringdummy1=ss(1)
            if('y'==stringdummy1) then
                flagstat=1
            elseif('n'==stringdummy1) then
                flagstat=0
            else
                Write(*,*)"ERROR: Input value of parameter STATON not valid"
                Write(*,*)"       ===> valid values 'n' or 'y'             "
                call stop_config
            endif
            stringdummy1=ss(2)
            if('y'==stringdummy1) then
                flaginst=1
            elseif('n'==stringdummy1) then
                flaginst=0
            else
                Write(*,*)"ERROR: Input value of parameter STATON not valid"
                Write(*,*)"       ===> valid values 'n' or 'y'             "
                call stop_config
            endif

        ELSEIF(line(1:3)=='603') THEN
            ! ###### DISSIPATION CALC ######
            call scan_string (line, 1, ss, narg)
            stringdummy1=ss(1)
            if('y'==stringdummy1) then
                flagbal=1
            elseif('n'==stringdummy1) then
                flagbal=0
            else
                Write(*,*)"ERROR: Input value of parameter BLANCEON not valid"
                Write(*,*)"       ===> valid values 'n' or 'y'             "
                call stop_config
            endif

        ELSEIF(line(1:3)=='604') THEN
            ! ###### START STAT WRITING ######
            call scan_string (line, 1, ss, narg)
            read(ss(1),*) tsta

        ELSEIF(line(1:3)=='605') THEN
            ! ###### SNAPSHOTS ######
            call scan_string (line, 2, ss, narg)
            read(ss(1),*) startsnap
            read(ss(2),*) savesnap

        ELSEIF(line(1:3)=='606') THEN
            ! ###### SLAB OUTPUT ######
            call scan_string (line, 1, ss, narg)
            stringdummy1=ss(1)
            if('y'==stringdummy1) then
                stst3flag=1
            elseif('n'==stringdummy1) then
                stst3flag=0
            else
                Write(*,*)"ERROR: Input value of parameter SLABDUMP(STST3) not valid"
                Write(*,*)"       ===> valid values 'n' or 'y'             "
                call stop_config
            endif

        ELSEIF(line(1:3)=='607') THEN
            ! ###### MOVIE OUTPUT ######
            call scan_string (line, 2, ss, narg)
            stringdummy1=ss(1)
            if('y'==stringdummy1) then
                movieflag=1
                read(ss(2),*) tframe
            elseif('n'==stringdummy1) then
                movieflag=0
                read(ss(2),*) tframe
            else
                Write(*,*)"ERROR: Input value of parameter MOVIEON not valid"
                Write(*,*)"       ===> valid values 'n' or 'y'             "
                call stop_config
            endif

        ENDIF

    222 CONTINUE

    444 close(15)

end subroutine read_from_bouin

SUBROUTINE SCAN_STRING (S, N, SS, M)
    IMPLICIT NONE
    INTEGER I, J, N, M, APO(2*N)
    CHARACTER*200 S ; CHARACTER*100 SS(N)
    CHARACTER*1, PARAMETER :: PRIME="'"
    !
    ! ----- Exits for n=0
    if(n==0) then
        Write(*, *) "SCAN_STRING has nothing to do"
        call stop_config
    endif
    
    ! ----- Looking for primes (', not ") within the input string
    i=0
    do j=1, 200
        if(s(j:j)==prime) then
            i=i+1 ; apo(i)=j
        endif
    enddo
    m=i/2

    if(i==0) then
        ! ----- Exits if no primes are found
        Write(*, *) "SCAN_STRING has found NO primes"
        Write(*,'(a100)') s
        call stop_config

    elseif(mod(i,2)/=0) then
        ! ----- Exits if an odd number of primes is found
        Write(*, *) "SCAN_STRING has found an even number of primes:", i
        call stop_config

    elseif(m/=n) then
        ! ----- Exits if m/=n, otherwise determines the substrings
        Write(*, *) "SCAN_STRING has found ", m, " substrings"
        Write(*, *) "SCAN_STRING  expected ", n, " substrings"
        call stop_config
    else
        do i=1, n
            ss(i)=s(apo(2*i-1)+1:apo(2*i)-1)
        enddo
    endif

end subroutine SCAN_STRING

SUBROUTINE STOP_CONFIG
    use mpih
    IMPLICIT NONE
    write(* ,*) 'STOP_CONFIG: The program will STOP! -----------------------------'
    call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    STOP 1
end subroutine STOP_CONFIG
