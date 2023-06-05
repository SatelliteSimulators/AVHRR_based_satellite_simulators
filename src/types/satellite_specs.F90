MODULE satellite_specs

  ! A module that gathers all satellite related info into one structure
  !
  ! Salomon Eliasson@smhi.se

  USE namelist_input, ONLY: name_list
  USE netcdf,         ONLY: nf90_open,             &
       nf90_nowrite,          &
       nf90_inquire,          &
       nf90_inquire_dimension,&
       nf90_inq_varid,        &
       nf90_get_var,          &
       nf90_close

  USE my_netcdftools, ONLY: check
  USE cosp_kinds,     ONLY: wp

  IMPLICIT NONE

  PUBLIC  :: assign_satellite_specs
  PRIVATE :: is_ch3b,overpass_time

  TYPE satellite

     ! string name (e.g. noaa14)
     CHARACTER (len=20)   :: name

     ! ascending or descending ('asc','dec','all')
     CHARACTER (len=3)    :: node

     ! average satellite zenith angle as a function of latitude

     ! local equatorial overpass time for the middle of the year
     REAL(wp), ALLOCATABLE    :: equatorial_overpass_time(:,:)

     REAL(wp), ALLOCATABLE    :: years(:), latitudes(:), months(:)

     REAL(wp)                 :: overpass

     LOGICAL :: is_ch3b ! flag is channel 3b is used

     INTEGER :: nlat
  END TYPE satellite

CONTAINS

  SUBROUTINE assign_satellite_specs(sat,options)

    IMPLICIT NONE

    TYPE(name_list), INTENT(in)    :: options
    TYPE(satellite), INTENT(out)   :: sat

    CHARACTER (len = 1000) :: filename
    INTEGER :: ncid, varid, dimid, len, nDimensions
    INTEGER              :: nyrs, nmon
    CHARACTER(len= 30)   :: str, field
    LOGICAL              :: exists
    INTEGER              :: year, month
    sat%node = options%L2b%node
    year  = options%epoch%year
    month = options%epoch%month

    ! ----------------------------
    ! Get the overpass times and average satellite zenith angles
    ! from file
    ! Go to the auxiliary files to find the appropriate local
    ! overpass times and satellite zenith angles

    WRITE(filename, '("data/satellite_specs/"A,"_L2b_specs.nc")')&
         TRIM(options%L2b%satellite)
    INQUIRE(FILE=TRIM(filename), EXIST=exists)
    IF (.NOT. exists) THEN
       WRITE (*,*) "satellite specs file is missing. ", TRIM(filename)
       STOP "stop in satellite_specs.F90"
    END IF

    ncid=0
    dimid=0
    CALL check( nf90_open( filename, nf90_nowrite, ncid ) )
    CALL check( nf90_inquire( ncid, nDimensions ) )

    ! Get the dimensions that span the look up tables
    DO dimid = 1, nDimensions
       CALL check (nf90_inquire_dimension( ncid, dimid, str, len ) )
       CALL check( nf90_inq_varid(ncid, str, varid) )
       IF (str .EQ. 'year') THEN
          ALLOCATE( sat%years(len) )
          CALL check( nf90_get_var(ncid, varid, sat%years),options%dbg,'reading '//str);
          nyrs = len
       ELSE IF (str .EQ. 'month') THEN
          ALLOCATE( sat%months(len) )
          CALL check( nf90_get_var(ncid, varid, sat%months),options%dbg,'reading '//str);
          nmon = len
          END IF
    END DO

    IF ( (sat%node.EQ."asc") .OR. (sat%node.EQ."all")) THEN
       field="eqtr_crossing_asc"
    ELSEIF (sat%node.EQ."dec") THEN
       field="eqtr_crossing_dec"
    END IF
    ALLOCATE ( sat%equatorial_overpass_time(nmon,nyrs) )
    CALL check( nf90_inq_varid(ncid, field, varid),options%dbg,'reading '//field)
    CALL check( nf90_get_var(ncid, varid, sat%equatorial_overpass_time) )

    CALL check( nf90_close(ncid))

    IF (options%L2b%doL2bSampling) THEN
       sat%overpass = OVERPASS_TIME(sat,year,month)
       IF (ISNAN(sat%overpass)) THEN
          WRITE (*,*) "overpass time is NaN. ",&
               "Probably the desired satellite didn't exist at this time"
          STOP "stopped in satellite_specs.F90"
       ELSE
          PRINT '(A)', " --- Subsampling the data to level2b"
          PRINT '(A,1x,A,1x,A,1x,F5.2)', &
               " --- Equatorial overpass time for",&
               TRIM(options%L2b%satellite),"=",sat%overpass
       END IF
    ELSE
       PRINT '(A)', " ---- Processing the data in a legacy fashion (same utc)"
    END IF
    ! Whether or not AVHRR channel 3B is being used
    sat%is_ch3b = IS_CH3B(sat,year,month)

  END SUBROUTINE assign_satellite_specs

  SUBROUTINE deallocate_satellite_specs(sat)

    IMPLICIT NONE

    TYPE(satellite), INTENT(inout)   :: sat

    IF (ALLOCATED(sat%years)) THEN
      DEALLOCATE(sat%equatorial_overpass_time,&
                 sat%years                   ,&
                 sat%months                   )
   END IF

  END SUBROUTINE deallocate_satellite_specs

  ELEMENTAL FUNCTION overpass_time(sat,year,month) RESULT(ot)
    !  Get the equatorial overpass time of a given satellite
    IMPLICIT NONE

    TYPE(satellite), INTENT(in) :: sat
    INTEGER, INTENT(in)         :: year,month
    INTEGER :: y_ind, m_ind

    REAL(wp) :: ot

    y_ind = 1
    m_ind = 1
    DO WHILE (sat%years(y_ind) < year)
       y_ind = y_ind+1
    END DO
    DO WHILE (sat%months(m_ind) < month)
       m_ind = m_ind+1
    END DO

    ot = sat%equatorial_overpass_time(m_ind,y_ind)

  END FUNCTION overpass_time

  ELEMENTAL FUNCTION is_ch3B(sat,iyr,imon)

    ! Different satellites have used different channels during the day time

    IMPLICIT NONE

    TYPE(satellite), INTENT(in) :: sat
    INTEGER, INTENT(in)         :: imon, iyr

    ! OUT
    LOGICAL :: is_ch3B

    is_ch3B = .TRUE.

    SELECT CASE (TRIM(sat%name))
    CASE ('noaa16')
       IF (iyr .LE. 2002 .OR. (iyr .EQ. 2003 .AND. imon .LE. 4)) THEN
          is_ch3B = .FALSE.
       END IF
    CASE ('noaa17')
       is_ch3b = .FALSE.
    CASE ('metop-a')
       is_ch3B = .FALSE.
    END SELECT
  END FUNCTION is_ch3B
END MODULE satellite_specs
