PROGRAM CLAAS_SIMULATOR
   ! Wrapper for claas satellite simulator
   !
   ! Code authors: Salomon Eliasson
   !
   !  salomon.eliasson@smhi.se

   USE AUXILIARY_FUNCTIONS, ONLY: &
      SOLAR_ZENITH_ANGLE
   USE CALC_FROM_MODEL, ONLY: &
      CALC_MODEL_VERTICAL_PROPERTIES, &
      GET_MODEL_SSA_AND_G
   USE CLAAS_FUNCTIONS, ONLY: &
      CHECK_GRID_AVERAGES, &
      CHECK_VARIABLES_SUBGRID, &
      CTTH, &
      GRID_AVERAGE, &
      DAY_ADD, &
      DAY_AVERAGE
   USE CLAAS_M, ONLY: &
      ALLOCATE_CLAAS, &
      CLAAS_TYPE, &
      DEALLOCATE_CLAAS, &
      GET_NAMELIST_CLAAS, &
      INITIALISE_CLAAS
   USE CLAAS_NETCDF, ONLY: &
      MAKE_NETCDF
   USE COSP_KINDS, ONLY: &
      WP
   USE DATA_CHECK, ONLY: &
      CHECK_VARIABLES, &
      MISSING
   USE FROM_COSP2, ONLY: &
      NUM_TRIAL_RES, &
      TRIAL_G_AND_W0
   USE HANDY, ONLY: &
      BUILD_FILENAME, &
      CHECK_FILE, &
      TIME_KEEPER
   USE INTERNAL_SIMULATOR, ONLY: &
      ALLOCATE_INTERNAL_SIMULATOR, &
      DEALLOCATE_INTERNAL_SIMULATOR, &
      INITIALISE_INTERNAL_SIMULATOR, &
      INTERNAL
   USE MODEL_INPUT, ONLY: &
      ALLOCATE_MODEL_MATRIX, &
      DEALLOCATE_MODEL_AUX, &
      DEALLOCATE_MODEL_MATRIX, &
      GET_MODEL_AUX, &
      INITIALISE_MODEL_MATRIX, &
      MODEL_TYPE, &
      READ_MODEL
   USE MY_MATHS, ONLY: &
      DAYOFYEAR
   USE NAMELIST_INPUT, ONLY: &
      DEALLOCATE_NAMELIST, &
      NAME_LIST
   USE OPTICS_M, ONLY: &
      DEALLOCATE_OPTICS, &
      POPULATE_EFFECTIVE_RADIUS_LUT
   USE SIMULATE_CLOUD_MICROPHYS, ONLY: &
      CLOUD_EFFECTIVE_RADIUS, &
      GET_CLOUD_MICROPHYSICS, &
      GET_PTAU
   USE SIMULATOR_INPUT_VARIABLES, ONLY: &
      ALLOCATE_SIM_INPUT, &
      DEALLOCATE_SIM_INPUT, &
      INITIALISE_SIM_INPUT, &
      SUBSET
   USE SIMULATOR_VARIABLES, ONLY: &
      ALLOCATE_SIMULATOR, &
      DEALLOCATE_SIMULATOR, &
      INITIALISE_SIMULATOR, &
      SATELLITE_SIMULATOR, &
      SET_TIMESTEP
   USE SUBCOLUMNS, ONLY: &
      GET_SUBCOLUMNS, &
      CORRECT_MODEL_CLOUD_FRACTION

   IMPLICIT NONE

   TYPE(internal)                           :: inter
   TYPE(model_type)                         :: model, previous
   TYPE(name_list)                          :: options
   TYPE(satellite_simulator)                :: sim
   TYPE(subset)                             :: sub
   CHARACTER(len=1000)                      :: namelist_file
   TYPE(claas_type)                         :: claas
   CHARACTER(3), PARAMETER                  :: simVersionNumber = '1'
   REAL(wp)                                 :: utc, day_of_year
   REAL(wp)                                 :: elapsed
   INTEGER                                  :: startTime, endTime, clock_rate
   REAL(wp), PARAMETER                      :: rho_w = 1._wp!      [10^3 kg/m^3]
   REAL(wp), PARAMETER                      :: rho_i = 0.93_wp !   [10^3 kg/m^3]
   REAL(wp), PARAMETER, DIMENSION(2)        :: rho = [rho_w, rho_i]
   REAL(wp), ALLOCATABLE, DIMENSION(:)      :: LST
   LOGICAL, ALLOCATABLE, DIMENSION(:)       :: isWithinGeoStationaryFOV
   REAL(wp), ALLOCATABLE, DIMENSION(:, :)   :: frac_out
   REAL(wp), ALLOCATABLE, DIMENSION(:)       :: tmpCC, tmpCV, tmpLon
   REAL(wp) :: fvr = -9._wp ! 1.e+20
   REAL(wp) :: clearFactor
   INTEGER  :: year, month, day
   INTEGER  :: ins, i, itime, d1, t1, t2
   INTEGER  :: nc, nlon, nlat, nl, ng, n_tbins, n_pbins
   INTEGER, PARAMETER :: phaseIsLiquid = 1, phaseIsIce = 2
   LOGICAL  :: atLeastOne

   CALL INITIALIZE_LOCAL_SCALARS()

   CALL GET_COMMAND_ARGUMENT(1, namelist_file)

   ! Get the options from the namelist
   CALL GET_NAMELIST_CLAAS(options, namelist_file)

   year = options%epoch%year
   month = options%epoch%month
   day = options%epoch%day
   options%simVersionNumber = simVersionNumber

   sim%netcdf_file = trim(options%paths%sim_output)
   IF (.NOT. options%overwrite_existing) THEN
      atLeastOne = .FALSE.
      IF (.NOT. CHECK_FILE(sim%netcdf_file)) THEN
         atLeastOne = .TRUE.
      END IF
      IF (.NOT. atLeastOne) THEN
         STOP " --- all files exist. set OVERWRITE=.TRUE. to overwrite"
      END IF
   END IF

   ! ------------------------
   ! GET MODEL DIMENSIONS etc
   ! ------------------------

   ! ---------------
   ! --- Auxiliary
   model%aux%netcdf_file = BUILD_FILENAME(options%paths%model_input, &
                                          y=year, m=month, d=day, model=options%model)
   PRINT '(A,A)', " --- Reading from file:", trim(model%aux%netcdf_file)

   CALL GET_MODEL_AUX(model%aux, options)
   ! get time here if using monthly files
   IF (.NOT. options%paths%dailyFiles) THEN
      CALL GET_MODEL_AUX(model%aux, options, .TRUE.)
   END IF

   nlat = model%aux%nlat
   nlon = model%aux%nlon
   nl = model%aux%nlev
   ng = model%aux%ngrids
   nc = options%ncols
   n_tbins = options%ctp_tau%n_tbins
   n_pbins = options%ctp_tau%n_pbins
   ! --------------- aux

   ALLOCATE (frac_out(nc, nl), &
             LST(ng), &
             isWithinGeoStationaryFOV(ng))

   ALLOCATE (tmpCC(nl))
   ALLOCATE (tmpCV(nl))
   ALLOCATE (tmpLon(ng))
   tmpCC = fvr
   tmpCV = fvr

   frac_out = missing
   LST = missing
   isWithinGeoStationaryFOV = .FALSE.

   CALL ALLOCATE_MODEL_MATRIX(model, ng, nl)
   CALL INITIALISE_MODEL_MATRIX(model, ng, nl, mv=0._wp)
   CALL ALLOCATE_SIM_INPUT(sub, options, ng, nl)
   CALL ALLOCATE_SIMULATOR(sim, ng)
   CALL ALLOCATE_INTERNAL_SIMULATOR(inter, nc, nl, options)

   ! ---------------
   ! Get valid geostationary lat lons.
   ! This should mask everything away from more than options%daynight%daylim (84 deg) of origo (lat=0,lon=0)
   tmpLon=model%aux%lon_v
   IF (MAXVAL(tmpLon) .gt. 180) tmpLon=MERGE(tmpLon-360,tmpLon,tmpLon.gt.180)
   isWithinGeoStationaryFOV(1:ng) = MERGE(.TRUE., .FALSE., &
                                          SQRT(model%aux%lat_v**2 + tmpLon**2)/options%daynight%daylim .LE. 1)
   sub%data_mask(1:ng, 1:nl) = SPREAD(isWithinGeoStationaryFOV, 1, nl)

   !---------------------

   ! ----------------------------
   !          CLAAS LUT's
   ! ------------------------------
   ! 1) Get the PODS per valid for the entire period
   ! 2) Get the right testing G AND W0 dependent on satellite and period.

   CALL ALLOCATE_CLAAS(claas, options, model%aux)

   ! read g0 and w0 look up tables
   options%sim_aux%LUT%ice%optics%re = &
      POPULATE_EFFECTIVE_RADIUS_LUT(options%CDR, phaseIsIce, .TRUE.)
   options%sim_aux%LUT%water%optics%re = &
      POPULATE_EFFECTIVE_RADIUS_LUT(options%CDR, phaseIsLiquid, .TRUE.)

   CALL TRIAL_G_AND_W0(options%sim_aux%LUT)

   ! --------- LUT

   ! ------------
   ! Check if daily netcdf file is already there
   sim%netcdf_file = BUILD_FILENAME(options%paths%sim_output, &
                                    CDR=options%CDR, model=options%model, y=year, m=month, d=day)

   !
   ! ------
   ! Restart these at the end of every full day
   CALL INITIALISE_CLAAS(claas, options, ng)
   CALL INITIALISE_SIM_INPUT(sub, options, ng, nl)
   CALL INITIALISE_SIMULATOR(sim, ng)
   CALL SET_TIMESTEP(sim, model%aux, options, day, t1, t2)

   DO itime = t1,t2
      ! start measuring the time it takes
      CALL SYSTEM_CLOCK(startTime, clock_rate)

      utc = MOD(model%aux%time(itime), 24._wp) + model%aux%ref%hour

      ! ---------------------
      ! -- get model data --

      write (*, '(2(a,1x),i4,2(a,i2),a,1x,f3.0)') &
         "reading model input for:", "yr =", year, ', mn =', month, &
         ', day =', day, &
         ', utc =', utc
      CALL READ_MODEL(model, itime, options)
      ! --------- end sampling model

      ! --------------
      ! local solar time
      !
      day_of_year = DAYOFYEAR(year, month, day)

      WHERE (isWithinGeoStationaryFOV(1:ng))
         LST(1:ng) = MOD(utc + 24._wp/360._wp*RESHAPE(model%aux%lon, (/ng/)), 24._wp)

         sub%solzen(1:ng) = &
                           SOLAR_ZENITH_ANGLE(model%aux%lat_v(1:ng), day_of_year, &
                                                LST(1:ng), ng)

         sub%sunlit(1:ng) = MERGE(1, 0, sub%solzen .LE. options%daynight%daylim)
      END WHERE

      PRINT *, "--- Clearing away bad values"
      clearFactor = 1._wp/(2*nc)
      WHERE (.NOT. sub%data_mask(1:ng, 1:nl))
         model%CC = MERGE(model%CC, 0._wp, model%CC .GE. clearFactor)
         model%CIWC = MERGE(model%CIWC, 0._wp, model%CC .GE. clearFactor)
         model%CLWC = MERGE(model%CLWC, 0._wp, model%CC .GE. clearFactor)
      ELSEWHERE
         model%CC = fvr
         model%CIWC = fvr
         model%CLWC = fvr
      END WHERE

      CALL CALC_MODEL_VERTICAL_PROPERTIES(ng, nl, model, sub, options)

      CALL GET_MODEL_SSA_AND_G(sub, options%sim_aux%LUT, ng, nl)

      IF (options%dbg > 0) THEN
         CALL CHECK_VARIABLES(ng, nl, options, sub%data_mask, model, sub)
      END IF

      ! ---------------------------
      ! Loop over grid points

      DO d1 = 1, ng
         IF (.NOT. isWithinGeoStationaryFOV(d1)) CYCLE
         tmpCC = model%CC(d1, 1:nl)
         tmpCV = model%CV(d1, 1:nl)
         CALL GET_SUBCOLUMNS(1, nc, nl, model%PSURF(d1), &
                           tmpCC, tmpCC, &
                           frac_out, options%subsampler)

         ! empty internal
         CALL INITIALISE_INTERNAL_SIMULATOR(inter, nc, nl, options)

         CALL GET_CLOUD_MICROPHYSICS(d1, nc, nl, frac_out, &
                                     inter, sub, options)

         DO ins = 1, nc
            ! --- loop over sub-columns
            IF (inter%cflag(ins) == 0 .OR. & ! cloud free
                  inter%cflag(ins) == 4) CYCLE ! cloudy due to FAR

            ! i.e., Simulate very thin (1), semi-transparent (2), and opaque (3) clouds

            CALL CTTH(d1, ins, model, sub, inter)
            IF (sub%sunlit(d1) .EQ. 1) THEN
               ! consider moving water and ice together to avoid if-statement."
               IF (inter%cph(ins) .EQ. 1) THEN
                  inter%reff(ins) = CLOUD_EFFECTIVE_RADIUS(d1, ins, nl, &
                                                           sub, inter, options%sim_aux%LUT%water%optics)
               ELSE
                  inter%reff(ins) = CLOUD_EFFECTIVE_RADIUS(d1, ins, nl, &
                                                           sub, inter, options%sim_aux%LUT%ice%optics)
               END IF
               inter%cwp(ins) = 0.667_wp*1.e-3* &
                                 inter%reff(ins)* &
                                 inter%tau(ins)* &
                                 rho(inter%cph(ins))

            END IF
         END DO
         IF (options%dbg > 0) THEN
            CALL CHECK_VARIABLES_SUBGRID(d1, nc, model, sub, inter, options, frac_out)
         END IF
         ! ------- post processing
         claas%av%hist2d_cot_ctp(d1, 1:n_tbins, 1:n_pbins, 1:2) = GET_PTAU(nc, n_tbins, n_pbins, options, &
                                                                           inter%tau, inter%ctp, inter%cph)

         CALL GRID_AVERAGE(d1, sub, inter, nc, claas%av)

         ! Keep track of optical depth occurances
         claas%av%cflag_tot(d1, 1) = claas%av%cflag_tot(d1, 1) + COUNT(inter%cflag .EQ. 0)
         claas%av%cflag_tot(d1, 2) = claas%av%cflag_tot(d1, 2) + COUNT(inter%cflag .EQ. 1)
         claas%av%cflag_tot(d1, 3) = claas%av%cflag_tot(d1, 3) + COUNT(inter%cflag .EQ. 2)
         claas%av%cflag_tot(d1, 4) = claas%av%cflag_tot(d1, 4) + COUNT(inter%cflag .EQ. 3)
         claas%av%cflag_tot(d1, 5) = claas%av%cflag_tot(d1, 5) + COUNT(inter%cflag .EQ. 4)
      END DO ! end grid

      ! I need to save the sum and number of elements for all
      ! variables and empty them after each time step
      ! make_grid average is really made for only one value per grid per day
      CALL DAY_ADD(claas, sub, ng, n_pbins, n_tbins)

      CALL SYSTEM_CLOCK(endTime)
      elapsed = elapsed + REAL(endTime - startTime)/REAL(clock_rate)
      CALL TIME_KEEPER(elapsed)
   END DO
   ! I need to save the sum and number of elements for all
   ! variables and empty them after each time step
   ! make_grid average is really made for only one value per grid per day
   CALL DAY_AVERAGE(claas, ng, isWithinGeoStationaryFOV)

   PRINT '(A,A)', "Writing simulator results to NetCDF file:&
     & ", TRIM(sim%netcdf_file)

   CALL MAKE_NETCDF(model, sim, options, day, S=sub, claas=claas)

   ! ---------------
   ! DEALLOCATE everything
   CALL DEALLOCATE_INTERNAL_SIMULATOR(inter, options)
   CALL DEALLOCATE_MODEL_AUX(model%aux)
   CALL DEALLOCATE_MODEL_MATRIX(model)
   IF (ALLOCATED(previous%T)) CALL DEALLOCATE_MODEL_MATRIX(previous)
   CALL DEALLOCATE_NAMELIST(options)
   CALL DEALLOCATE_SIMULATOR(sim)
   CALL DEALLOCATE_SIM_INPUT(sub, options)
   CALL DEALLOCATE_claas(claas, options)

   DEALLOCATE (frac_out, LST)
   DEALLOCATE (isWithinGeoStationaryFOV)
   DEALLOCATE (tmpCC, tmpCV, tmpLon)

   SELECT CASE (options%cloudMicrophys%cf_method)
      ! ONLY use POD method for cloud fraction
   CASE (1)
      DEALLOCATE (options%sim_aux%POD_layers, &
                  options%sim_aux%random_numbers, &
                  options%sim_aux%POD_tau_bin_centers, &
                  options%sim_aux%POD_tau_bin_edges)
   END SELECT

   PRINT *, "Finished !!!"
CONTAINS

   !
   !------
   !

   SUBROUTINE INITIALIZE_LOCAL_SCALARS()
      ! Immediately initialize everything
      utc = 0._wp
      elapsed = 0._wp
      startTime = 0
      endTime = 0
      day_of_year = 0._wp
      t1 = 0
      t2 = 0
      itime = 1
      i = 0

   END SUBROUTINE INITIALIZE_LOCAL_SCALARS

END PROGRAM claas_SIMULATOR
