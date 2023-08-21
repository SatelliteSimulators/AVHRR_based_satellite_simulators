MODULE simulator_variables
  ! Structure types for the simulated products
  !
  ! Salomon.Eliasson@smhi.se

  USE cosp_kinds,     ONLY: wp
  USE handy,          ONLY: &
       build_filename
  USE model_input,    ONLY: &
       get_model_aux,       &
       model_aux
  USE my_maths,       ONLY: &
       day_since_day,       &
       day_since_year
  USE namelist_input, ONLY: name_list


  IMPLICIT NONE

  PUBLIC :: allocate_simulator,  &
            deallocate_simulator,&
            initialise_simulator,&
            satellite_simulator, &
            set_timeStep

  !-------------------
  ! SIMULATOR OUTPUT
  ! 'netcdf_file' = name of file
  ! 'epoch'     = days since 1970/1/1

  TYPE satellite_simulator
     CHARACTER(len=200)    :: netcdf_file
     INTEGER               :: epoch
     INTEGER               :: date,dtg
     INTEGER               :: fvi = -9
     REAL(wp)              :: fvr = -9._wp ! 1.e+20
     INTEGER               :: ncols
     REAL(wp), ALLOCATABLE :: time_of_day  (:)
  END TYPE satellite_simulator

CONTAINS

 SUBROUTINE allocate_simulator(sim,ngrids)

    IMPLICIT NONE
    TYPE(satellite_simulator), INTENT(out) :: sim
    INTEGER, INTENT(in)           :: ngrids

    ALLOCATE (sim%time_of_day(ngrids) )

  END SUBROUTINE allocate_simulator
  SUBROUTINE initialise_simulator(sim,ngrids)

    IMPLICIT NONE
    TYPE(satellite_simulator), INTENT(inout) :: sim
    INTEGER, INTENT(in)           :: ngrids

    sim%time_of_day(1:ngrids) = -9._wp

  END SUBROUTINE initialise_simulator
  SUBROUTINE deallocate_simulator(sim)

    IMPLICIT NONE
    TYPE(satellite_simulator), INTENT(inout) :: sim

    DEALLOCATE (sim%time_of_day )

  END SUBROUTINE deallocate_simulator

  SUBROUTINE set_timeStep(sim,aux,options,iday,t1,t2)
    ! Set time step in loop. This snippet of code is shared by all
    ! simulators, so I moved it here


    IMPLICIT NONE

    TYPE(satellite_simulator), INTENT(inout) :: sim
    TYPE(model_aux),           INTENT(inout) :: aux
    TYPE(name_list), INTENT(inout)           :: options
    INTEGER, INTENT(in)  :: iday
    INTEGER, INTENT(OUT) :: t1,t2
    INTEGER  :: year,month
    REAL(wp) :: deltaTime

    year  = options%epoch%year
    month = options%epoch%month

    ! ---------
    ! TIME
    ! set date and dtg in daily_product; assumed UTC noon for dtg
    sim%date = year*10000 + month*100 + iday
    sim%dtg  = ( sim%date ) * 100 + 12

    ! get the start day since 1970 1 1
    sim%epoch = DAY_SINCE_YEAR(1970,year,month,iday)

    IF (options%paths%dailyFiles) THEN
       aux%netcdf_file = BUILD_FILENAME(&
            options%paths%model_input,y=year,m=month,d=iday)
       IF (ALLOCATED(aux%time)) DEALLOCATE(aux%time)
       CALL get_model_aux(aux, options, .TRUE.) ! .true. means to get 'time' only
       t1 = 1
       t2 = aux%ntlen
       !EvMb ... in case of dailyFiles=True and racmo operate internal
       !     time-loop t from t1 to t2-1 ... simply because the
       !     interval [00:01] or [00:03] UTC is not going to be
       !     samples
       if (options%model.eq.'racmo') t2=t2-1
       !EvMe

    ELSE ! if monthly files expected
       deltaTime=aux%time(2)-aux%time(1)
       t1 = INT((iday-1)*24._wp/deltaTime +1)
       t2 = INT(iday*24._wp/deltaTime)

    END IF
    ! this brings day in the right frame. Number of days since
    ! reference time

    ! --------- end time
  END SUBROUTINE set_timeStep

END MODULE simulator_variables
