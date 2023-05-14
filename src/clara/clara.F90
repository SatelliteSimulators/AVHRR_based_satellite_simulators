PROGRAM CLARA_SIMULATOR
  ! Wrapper for clara satellite simulator
  !
  ! Code authors: Salomon Eliasson, Jan Fokke Meirink
  !
  ! salomon.eliasson@smhi.se


  USE AUXILIARY_FUNCTIONS,       ONLY:&
       FIND_TEMPERATURE_INVERSIONS,   &
       SOLAR_ZENITH_ANGLE
  USE BRIGHTNESS_TEMPERATURES,   ONLY:&
       CORRECTED_TEMPERATURE_PROFILE, &
       TB_ICARUS
  USE CALC_FROM_MODEL,           ONLY:&
       CALC_MODEL_VERTICAL_PROPERTIES,&
       GET_MODEL_SSA_AND_G
  USE CLARA_FUNCTIONS,           ONLY:&
       CHECK_GRID_AVERAGES,           &
       CHECK_VARIABLES_SUBGRID,       &
       CTTH,                          &
       GRID_AVERAGE,                  &
       DAY_ADD,                       &
       DAY_AVERAGE
  USE CLARA_M,                   ONLY:&
       ALLOCATE_CLARA,                &
       CLARA_TYPE,                    &
       DEALLOCATE_CLARA,              &
       GET_NAMELIST_CLARA,            &
       INITIALISE_CLARA
  USE CLARA_NETCDF,              ONLY:&
       MAKE_NETCDF
  USE COSP_KINDS,                ONLY:&
       WP
  USE DATA_CHECK,                ONLY:&
       CHECK_VARIABLES,               &
       MISSING
  USE FROM_COSP2,                ONLY:&
       NUM_TRIAL_RES,                 &
       TRIAL_G_AND_W0
  USE HANDY,                     ONLY:&
       BUILD_FILENAME,                &
       CHECK_FILE,                    &
       TIME_KEEPER
  USE INTERNAL_SIMULATOR,        ONLY:&
       ALLOCATE_INTERNAL_SIMULATOR,   &
       DEALLOCATE_INTERNAL_SIMULATOR, &
       INITIALISE_INTERNAL_SIMULATOR, &
       INTERNAL
  USE MODEL_INPUT,               ONLY:&
       ALLOCATE_MODEL_MATRIX,         &
       DEALLOCATE_MODEL_AUX,          &
       DEALLOCATE_MODEL_MATRIX,       &
       GET_MODEL_AUX,                 &
       INITIALISE_MODEL_MATRIX,       &
       MODEL_TYPE,                    &
       READ_MODEL
  USE MY_MATHS,                  ONLY:&
       DAYOFYEAR
  USE NAMELIST_INPUT,            ONLY:&
       DEALLOCATE_NAMELIST,           &
       NAME_LIST
  USE OPTICS_M,                  ONLY:&
       DEALLOCATE_OPTICS,             &
       POPULATE_EFFECTIVE_RADIUS_LUT
  USE SATELLITE_SPECS,           ONLY:&
       ASSIGN_SATELLITE_SPECS,        &
       DEALLOCATE_SATELLITE_SPECS,    &
       SATELLITE
  USE SIMULATE_CLOUD_MICROPHYS,  ONLY:&
       CLOUD_EFFECTIVE_RADIUS,        &
       GET_CLOUD_MICROPHYSICS,        &
       GET_PTAU
  USE SIMULATOR_INPUT_VARIABLES, ONLY:&
       ALLOCATE_SIM_INPUT,            &
       DEALLOCATE_SIM_INPUT,          &
       INITIALISE_SIM_INPUT,          &
       SUBSET
  USE SIMULATOR_VARIABLES,       ONLY:&
       ALLOCATE_SIMULATOR,            &
       DEALLOCATE_SIMULATOR,          &
       INITIALISE_SIMULATOR,          &
       SATELLITE_SIMULATOR,           &
       SET_TIMESTEP
  USE SUBCOLUMNS,                ONLY:&
       GET_SUBCOLUMNS,                &
       CORRECT_MODEL_CLOUD_FRACTION

  IMPLICIT NONE

  TYPE(internal)                       :: inter
  TYPE(model_type)                     :: model,previous
  TYPE(name_list)                      :: options
  TYPE(satellite)                      :: sat
  TYPE(satellite_simulator)            :: sim
  TYPE(subset)                         :: sub
  CHARACTER(len=1000)                  :: namelist_file
  TYPE(clara_type)                     :: clara
  CHARACTER(3),PARAMETER               :: simVersionNumber='1.1'
  REAL(wp)                             :: utc,day_of_year
  REAL(wp)                             :: elapsed
  INTEGER                              :: startTime,endTime,clock_rate
  REAL(wp),PARAMETER                   :: rho_w = 1._wp!      [10^3 kg/m^3]
  REAL(wp),PARAMETER                   :: rho_i = 0.93_wp !   [10^3 kg/m^3]
  REAL(wp),PARAMETER, DIMENSION(2)     :: rho = [rho_w,rho_i]
  REAL(wp),ALLOCATABLE,DIMENSION(:)    :: LST,TOD
  REAL(wp),ALLOCATABLE,DIMENSION(:,:)  :: frac_out2
  REAL(wp),ALLOCATABLE,DIMENSION(:,:,:):: frac_out

  INTEGER :: year,month,day
  INTEGER :: ins,i,itime,d1,t1,t2
  INTEGER :: nc,nlon,nlat,nl,ng,n_tbins,n_pbins
  INTEGER, PARAMETER :: phaseIsLiquid = 1,phaseIsIce = 2
  LOGICAL :: atLeastOne,L2b,need2average

  CALL INITIALIZE_LOCAL_SCALARS()

  CALL GET_COMMAND_ARGUMENT(1,namelist_file)

  ! Get the options from the namelist
  CALL GET_NAMELIST_CLARA(options,namelist_file)

  year  = options%epoch%year
  month = options%epoch%month
  day   = options%epoch%day
  L2b   = options%L2b%doL2bSampling
  options%simVersionNumber=simVersionNumber

  need2Average = (.NOT.options%L2b%doL2bSampling)
  ! Check if I am going to do anything at all, and immediately leave if not

  sim%netcdf_file=trim(options%paths%sim_output)
  IF (.NOT. options%overwrite_existing) THEN
     atLeastOne=.FALSE.
     IF ( .NOT.CHECK_FILE(sim%netcdf_file) ) THEN
        atLeastOne=.TRUE.
     END IF
     IF (.NOT.atLeastOne) THEN
        STOP " --- all files exist. set OVERWITE=.TRUE. to overwrite"
     END IF
  END IF

  ! ------------------------
  ! GET MODEL DIMENSIONS etc
  ! ------------------------

  CALL ASSIGN_SATELLITE_SPECS(sat,options)
  ! ---------------
  ! --- Auxiliary
  model%aux%netcdf_file = BUILD_FILENAME(options%paths%model_input,&
       y=year,m=month,d=day,model=options%model)
  PRINT '(A,A)'," --- Reading from file:",trim(model%aux%netcdf_file)

  CALL GET_MODEL_AUX(model%aux,options)
  ! get time here if using monthly files
  IF (.NOT.options%paths%dailyFiles) THEN
     CALL GET_MODEL_AUX(model%aux,options,.TRUE.)
  END IF

  nlat = model%aux%nlat
  nlon = model%aux%nlon
  nl   = model%aux%nlev
  ng   = model%aux%ngrids
  nc   = options%ncols
  n_tbins = options%ctp_tau%n_tbins
  n_pbins = options%ctp_tau%n_pbins
  ! --------------- aux

  ALLOCATE(frac_out(ng,nc,nl),&
       frac_out2   (   nc,nl),&
       LST         (ng      ),&
       TOD         (ng      ))

  frac_out = missing
  frac_out2= missing
  LST      = missing
  TOD      = missing

  CALL ALLOCATE_MODEL_MATRIX(model,ng,nl)
  CALL INITIALISE_MODEL_MATRIX(model,ng,nl,mv=0._wp)
  CALL ALLOCATE_SIM_INPUT(sub,options,ng,nl)
  CALL ALLOCATE_SIMULATOR(sim,ng)
  CALL ALLOCATE_INTERNAL_SIMULATOR(inter,nc,nl,options)


  ! ----------------------------
  !          CLARA LUT's
  ! ------------------------------
  ! 1) Get the PODS per valid for the entire period
  ! 2) Get the right testing G AND W0 dependent on satellite and period.

  CALL ALLOCATE_CLARA(clara,options,model%aux)

  ! read g0 and w0 look up tables
  options%sim_aux%LUT%ice%optics%re  = &
       POPULATE_EFFECTIVE_RADIUS_LUT(options%CDR,phaseIsIce,sat%is_ch3b)
  options%sim_aux%LUT%water%optics%re= &
       POPULATE_EFFECTIVE_RADIUS_LUT(options%CDR,phaseIsLiquid,sat%is_ch3b)

  CALL TRIAL_G_AND_W0(options%sim_aux%LUT)
  ! --------- LUT


  ! ------------
  ! Check if daily netcdf file is already there
  sim%netcdf_file = BUILD_FILENAME(options%paths%sim_output,&
       CDR=options%CDR,model=options%model,y=year,m=month,d=day,&
       sat=options%L2b%satellite,node=options%L2b%node)
  !
  ! ------
  ! Restart these at the end of every full day
  CALL INITIALISE_CLARA(clara,options,ng)
  CALL INITIALISE_SIM_INPUT(sub,options,ng,nl,model%aux%lat,sat)
  CALL INITIALISE_SIMULATOR(sim,ng)
  CALL SET_TIMESTEP(sim,model%aux,options,day,t1,t2)

  ! start measuring the time it takes
  CALL SYSTEM_CLOCK (startTime,clock_rate)

  utc = MOD(model%aux%time(itime),24._wp) + model%aux%ref%hour
  day_of_year = DAYOFYEAR(year,month,day)
  TOD = sim%fvr

  ! ---------------------
  ! -- Sample Level 2B --
  IF (L2b) THEN
     !The following routine makes interpolated model fields, and
     ! a isL2b If the local time is 'x', what is that time in
     ! UTC at that longitude?

     CALL MAKE_LEVEL2B(year,month,day,t1,t2,&
          sat%overpass,options,previous,model,TOD)
     ! next time skip to the next day since I use t1 and t2 in the above routine

  ELSE
     write(*,'(2(a,1x),i4,2(a,i2),a,1x,f3.0)') &
          "reading model input for:","yr =",year,', mn =',month,&
          ', day =',day,&
          ', utc =',utc
     CALL READ_MODEL(model,itime,options)
  END IF
  ! --------- end sampling model

  ! --------------
  ! local solar time
  !
  sim%time_of_day=TOD

  LST(1:ng) = &
       & MERGE(sat%overpass,&
       MOD(utc+24._wp/360._wp*RESHAPE(model%aux%lon,(/ng/)),24._wp),L2b)

  sub%solzen(1:ng) = &
       SOLAR_ZENITH_ANGLE(model%aux%lat_v(1:ng),day_of_year,&
       LST(1:ng),ng)

  sub%sunlit(1:ng)=MERGE(1,0,sub%solzen .LE. options%daynight%daylim)

  ! Need to immediately clear away bad values
  CALL CORRECT_MODEL_CLOUD_FRACTION(model,ng,nl,nc)

  CALL CALC_MODEL_VERTICAL_PROPERTIES(ng,nl,model,sub,options)

  CALL GET_MODEL_SSA_AND_G(sub,options%sim_aux%LUT,ng,nl)

  sub%Tcorr(1:ng,1:nl+1) = CORRECTED_TEMPERATURE_PROFILE(model,ng,nl)

  sub%inv_layers(1:ng,1:nl) = FIND_TEMPERATURE_INVERSIONS(sub,options,ng,nl)

  IF (options%dbg>0) THEN
     CALL CHECK_VARIABLES(ng,nl,options,sub%data_mask,model,sub)
  END IF

  CALL GET_SUBCOLUMNS(ng,nc,nl,model%PSURF,&
       model%CC,model%CV,frac_out(1:ng,1:nc,1:nl),&
       options%subsampler)

  ! ---------------------------
  ! Loop over grid points
  DO d1 =1,ng

     ! I cannot move these 3 to outside the grid loop,
     ! because that uses too much memory

     ! empty internal
     CALL INITIALISE_INTERNAL_SIMULATOR(inter,nc,nl,options)

     frac_out2=frac_out(d1,:,:)
     CALL GET_CLOUD_MICROPHYSICS(d1,nc,nl,frac_out2,&
          inter,sub,options)

     inter%Tb(1:nc)=TB_ICARUS(d1,nc,nl,model,sub,inter,frac_out2)

     DO ins = 1,nc
        ! --- loop over sub-columns
        IF (inter%cflag(ins) .EQ. 0) CYCLE ! i.e., I should also
        ! simulate very thin clouds

        ! the CTTH routines need some rewriting if I want to
        ! avoid the loop
        CALL CTTH(d1,ins,options%CDR,model,sub,inter)
        IF (sub%sunlit(d1).EQ.1) THEN
           ! consider moving water and ice together to avoid if-statement."
           IF (inter%cph(ins) .EQ. 1) THEN
              inter%reff(ins) = CLOUD_EFFECTIVE_RADIUS(d1,ins,nl,&
                   sub,inter,options%sim_aux%LUT%water%optics)
           ELSEIF (inter%cph(ins) .EQ. 2) THEN
              inter%reff(ins) = CLOUD_EFFECTIVE_RADIUS(d1,ins,nl,&
                   sub,inter,options%sim_aux%LUT%ice%optics)
           END IF
        END IF
     END DO
     IF (sub%sunlit(d1).EQ.1) THEN
        WHERE(inter%cflag(1:nc).GT.1)
           inter%cwp      (1:nc) = 0.667_wp*1.e-3*&
                inter%reff(1:nc)*&
                inter%tau       (1:nc)*&
                rho(inter%cph   (1:nc))
        ELSEWHERE
           inter%cwp(1:nc) = 0._wp
        END WHERE
     END IF

     IF (options%dbg>0) THEN
        CALL CHECK_VARIABLES_SUBGRID(d1,nc,model,sub,inter,options,&
             frac_out2)
     END IF
     ! ------- post processing
     clara%av%hist2d_cot_ctp(d1,1:n_tbins,1:n_pbins,1:2) = &
          GET_PTAU(nc,n_tbins,n_pbins,options,&
          inter%tau,inter%ctp,inter%cph)

     CALL GRID_AVERAGE(d1,sub,inter,nc,clara%av)
  END DO ! end ng
  IF (options%dbg>0) CALL CHECK_GRID_AVERAGES(model,sub,options,clara%av)

  IF (need2Average) THEN
     ! I need to save the sum and number of elements for all
     ! variables and empty them after each time step
     ! make_grid average is really made for only one value per grid per day
     CALL DAY_ADD(clara,sub,ng,n_pbins,n_tbins)
  END IF

  CALL SYSTEM_CLOCK(endTime)
  elapsed = elapsed+REAL(endTime-startTime)/REAL(clock_rate)
  CALL TIME_KEEPER(elapsed)

  IF (need2Average) THEN
     ! I need to save the sum and number of elements for all
     ! variables and empty them after each time step
     ! make_grid average is really made for only one value per grid per day
     CALL DAY_AVERAGE(clara,ng)
     CALL CHECK_GRID_AVERAGES(model,sub,options,clara%av)
  END IF

  PRINT '(A,A)',"Writing simulator results to NetCDF file:&
       & ",TRIM(sim%netcdf_file)

  CALL MAKE_NETCDF(model,sim,options,day,S=sub,sat=sat,clara=clara)

  ! ---------------
  ! DEALLOCATE everything
  CALL DEALLOCATE_INTERNAL_SIMULATOR(inter,options)
  CALL DEALLOCATE_MODEL_AUX(model%aux)
  CALL DEALLOCATE_MODEL_MATRIX(model)
  IF (ALLOCATED(previous%T)) CALL DEALLOCATE_MODEL_MATRIX(previous)
  CALL DEALLOCATE_NAMELIST(options)
  CALL DEALLOCATE_SIMULATOR(sim)
  CALL DEALLOCATE_SIM_INPUT(sub,options)
  CALL DEALLOCATE_SATELLITE_SPECS(sat)
  CALL DEALLOCATE_CLARA(clara,options)

  DEALLOCATE (frac_out,frac_out2,LST,TOD)

  SELECT CASE (options%cloudMicrophys%cf_method)
     ! ONLY use POD method for cloud fraction
  CASE (1)
     DEALLOCATE(options%sim_aux%POD_layers   ,&
          options%sim_aux%random_numbers     ,&
          options%sim_aux%POD_tau_bin_centers,&
          options%sim_aux%POD_tau_bin_edges   )
  END SELECT

  PRINT *,"Finished !!!"
CONTAINS

  !
  !------
  !

  SUBROUTINE INITIALIZE_LOCAL_SCALARS()
    ! Immediately initialize everything
    utc                      = 0._wp
    elapsed                  = 0._wp
    startTime                = 0
    endTime                  = 0
    day_of_year              = 0._wp
    t1                       = 0
    t2                       = 0
    itime                    = 1
    i                        = 0

  END SUBROUTINE INITIALIZE_LOCAL_SCALARS

END PROGRAM CLARA_SIMULATOR
