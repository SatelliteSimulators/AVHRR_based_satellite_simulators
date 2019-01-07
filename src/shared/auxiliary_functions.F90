MODULE AUXILIARY_FUNCTIONS
  !
  !
  ! SALOMON.ELIASSON@SMHI.SE

  USE DATA_CHECK,                 ONLY: MINIMUMTROPOPAUSEPRESSURE
  USE COSP_KINDS,                 ONLY: WP
  USE COSP_PHYS_CONSTANTS,        ONLY: GRAV ! [M/S^2]
  USE COSP_MATH_CONSTANTS,        ONLY: PI
  USE HANDY,                      ONLY: CHECK_FILE
  USE INTERNAL_SIMULATOR,         ONLY: INTERNAL
  USE MODEL_INPUT,                ONLY: MODEL_TYPE,MODEL_AUX
  USE MY_NETCDFTOOLS,             ONLY: CHECK
  USE NETCDF,                     ONLY:&
       NF90_CLOSE,                     &
       NF90_GET_VAR,                   &
       NF90_INQUIRE,                   &
       NF90_INQUIRE_DIMENSION,         &
       NF90_INQ_VARID,                 &
       NF90_OPEN,                      &
       NF90_NOWRITE              
  USE NAMELIST_INPUT,            ONLY: NAME_LIST
  USE SIMULATOR_INPUT_VARIABLES, ONLY: SUBSET

  IMPLICIT NONE

  PUBLIC :: AUX_DATA,              &
       CTP_TAU,                    &
       FIND_TEMPERATURE_INVERSIONS,&
       SOLAR_ZENITH_ANGLE


  TYPE FLAG
     LOGICAL :: model_lat_descend,data_lat_descend
     LOGICAL :: model_lon_180,data_lon_180
     LOGICAL :: modelIsFiner, modelIsCoarser
  END TYPE FLAG

CONTAINS

  SUBROUTINE AUX_DATA(file,variable,strDim1,&
       strDim2,data2D,conform2model,aux,&
       data1D,strDim3,data3D,nX,nY,nZ,dbg)

    IMPLICIT NONE
    

    CHARACTER(len=1000), INTENT(in)       :: file ! the full path
    CHARACTER(len=*), INTENT(in)          :: variable,strDim1
    CHARACTER(len=*), INTENT(in),OPTIONAL :: strDim2,strDim3
    LOGICAL, INTENT(in),OPTIONAL          :: conform2model
    TYPE(model_aux), INTENT(in),OPTIONAL  :: aux
    INTEGER, INTENT(in), OPTIONAL         :: dbg 

    
    ! output
    REAL(wp), ALLOCATABLE, INTENT(out), OPTIONAL :: data1D(:)
    REAL(wp), ALLOCATABLE, INTENT(out), OPTIONAL :: data2D(:,:)
    REAL(wp), ALLOCATABLE, INTENT(out), OPTIONAL :: data3D(:,:,:)
    INTEGER, INTENT(out), OPTIONAL ::nX,nY,nZ
    
    TYPE(flag)            :: flags
    CHARACTER(len = 100)  :: str
    LOGICAL               :: matchModel
    INTEGER               :: ndim1,ndim2, ndim3,ndims
    INTEGER:: ncid, varid, dimid, nDimensions, len
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: dim1, dim2, dim3

    ! putting the mask together
    INTEGER :: lt, ln, iln, ilt
    REAL(wp) :: dLon,dLat

    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: netcdf_data3D,tmp_data3D
    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: netcdf_data2D,tot,tmp_data2D

    ! =====
    ! READ 
    !
    ndims=1
    matchModel=.FALSE.
    IF (PRESENT(strDim2)) ndims=ndims+1
    IF (PRESENT(strDim3)) ndims=ndims+1
    IF (PRESENT(conform2model)) matchModel=.TRUE.

    ncid=0
    varid=0
    dimid=0
    ndim1=0
    ndim2=0
    ndim3=0
    CALL CHECK( nf90_open( file, nf90_nowrite, ncid ) )
    CALL CHECK( nf90_inquire( ncid, nDimensions ) ) 

    ! Get the dimensions that span the look up tables
    DO dimid = 1, nDimensions
       
       CALL CHECK( nf90_inquire_dimension( ncid, dimid, str, len ) )
       CALL CHECK( nf90_inq_varid(ncid, str, varid) )

       IF (str.EQ.strDim1) THEN
          ALLOCATE( dim1(len) )
          CALL CHECK( nf90_get_var(ncid, varid, dim1) ,dbg, 'reading '//trim(str))
          ndim1 = len
       ELSEIF (str.EQ.strDim2) THEN
          ALLOCATE( dim2(len) )
          CALL CHECK( nf90_get_var(ncid, varid, dim2) ,dbg, 'reading '//trim(str))
          ndim2 = len
       ELSEIF (str.EQ.strDim3) THEN
          ALLOCATE( dim3(len) )
          CALL CHECK( nf90_get_var(ncid, varid, dim3) ,dbg, 'reading '//trim(str))
          ndim3 = len
       END IF
    END DO

    CALL CHECK( nf90_inq_varid(ncid, TRIM(variable), varid))

    SELECT CASE(ndims)
    CASE(1)
       ALLOCATE ( data1D(ndim1) )       
       IF (ndims.GT.1) THEN
          CALL CHECK( nf90_get_var(ncid, varid, data1D),dbg,'reading '//TRIM(variable))
       ELSE
          data1D=dim1
       END IF
    CASE(2)
       ALLOCATE ( netcdf_data2D(ndim1, ndim2) )
       CALL CHECK( nf90_get_var(ncid, varid, netcdf_data2D),dbg,'reading '//TRIM(variable))
    CASE(3)
       ALLOCATE ( netcdf_data3D(ndim1, ndim2, ndim3) )
       CALL CHECK(nf90_get_var(ncid,varid,netcdf_data3D),dbg,'reading '//TRIM(variable))
    END SELECT

    CALL CHECK( nf90_close(ncid))

    ! allocate
    SELECT CASE(ndims)
    CASE(2)
       IF (matchModel) THEN
          ALLOCATE(data2D(aux%nlon,aux%nlat),&
               tmp_data2D(aux%nlon,aux%nlat))
       ELSE
          ALLOCATE(data2D(ndim1,ndim2))
       END IF
    CASE(3)
       IF (matchModel) THEN
          ALLOCATE(data3D(aux%nlon,aux%nlat,ndim3),&
               tmp_data3D(aux%nlon,aux%nlat,ndim3))
       ELSE
          ALLOCATE(data3D(ndim1,ndim2,ndim3))
       END IF
    END SELECT
    !
    ! ========================

    IF (matchModel) THEN
       SELECT CASE(ndims)
       CASE(2)

          CALL CONFORM_GRID_DATA_TO_MODEL(aux,flags,nlon=ndim1,nlat=ndim2,lon=dim1,lat=dim2,&
               data2=netcdf_data2D)

       CASE(3)
          CALL CONFORM_GRID_DATA_TO_MODEL(aux,flags,nlon=ndim1,nlat=ndim2,lon=dim1,lat=dim2,&
               ndim3=ndim3,data3=netcdf_data3D)
       END SELECT
       
    ! --------------------
    ! REGRID
    !
       ALLOCATE(tot(aux%nlon,aux%nlat))
       
       IF (flags%modelIsCoarser) THEN
       ! assuming the netcdf data has a FINER resolution than the model

          ilt = 1
          DO lt = 1, ndim2
             iln = 1
             DO ln = 1, ndim1
                IF (dim1(ln) > aux%lon(iln,1) .AND. iln<aux%nlon) iln=iln+1
                IF (flags%model_lat_descend) THEN
                   IF (dim2(lt) < aux%lat(1,ilt) .AND. ilt<aux%nlat) ilt=ilt+1
                ELSE
                   IF (dim2(lt) > aux%lat(1,ilt) .AND. ilt<aux%nlat) ilt=ilt+1
                END IF
                
                tot(iln,ilt) = tot(iln,ilt) + 1
                SELECT CASE(ndims)
                CASE(2)
                   tmp_data2D(iln,ilt) = tmp_data2D(iln,ilt) +&
                        & netcdf_data2D(ln,lt)
                CASE(3)
                   tmp_data3D(iln,ilt,1:ndim3) = tmp_data3D(iln,ilt,:) +&
                        & netcdf_data3D(ln,lt,1:ndim3)
                END SELECT
             END DO
          END DO
          SELECT CASE(ndims)
          CASE(2)
             data2D(1:aux%nlon,1:aux%nlat) = tmp_data2D/tot
          CASE(3)
             data3D(1:aux%nlon,1:aux%nlat,1:ndim3) = tmp_data3D/SPREAD(tot,3,ndim3)
          END SELECT
       ELSEIF (flags%modelIsFiner) THEN
       ! assuming the netcdf data has a COARSER resolution than the model

          dLon=ABS((dim1(1)-dim1(2))/2)
          dLat=ABS((dim2(1)-dim2(2))/2)
          ilt = 1
          DO lt = 1,aux%nlat
             iln = 1
             IF (flags%model_lat_descend) THEN
                IF (dim2(ilt)-dLat > aux%lat(1,lt) .AND. ilt<ndim2) ilt=ilt+1
             ELSE
                IF (dim2(ilt)+dLat < aux%lat(1,lt) .AND. ilt<ndim2) ilt=ilt+1
             END IF
!             PRINT '("lat(",i3,")=",f5.1,"   dim2(",i3,")=",f5.1)',lt,aux%lat(1,lt),ilt,dim2(ilt)

             DO ln=1,aux%nlon
                IF (dim1(iln)+dLon < aux%lon(ln,1) .AND. iln<ndim1) iln=iln+1
                SELECT CASE(ndims)
                CASE(2)
                   data2D(ln,lt) = netcdf_data2D(iln,ilt)
                CASE(3)
                   data3D(ln,lt,1:ndim3) = netcdf_data3D(iln,ilt,1:ndim3)
                END SELECT
             END DO
          END DO

       ELSEIF (ndim1.EQ.aux%nlon .AND. ndim2.EQ.aux%nlat) THEN
          SELECT CASE(ndims)
          CASE(2)
             data2D = netcdf_data2D
          CASE(3)
             data3D = netcdf_data3D
          END SELECT
       ELSE
!stop "I need to do something here for unfortunate grids"
       END IF

       DEALLOCATE(tot)
    ELSE
       SELECT CASE(ndims)
       CASE(2)
          data2D=netcdf_data2D
       CASE(3)
          data3D=netcdf_data3D
       END SELECT
    END IF
    
    DEALLOCATE (dim1)
    SELECT CASE(ndims)
    CASE(2)
       DEALLOCATE(netcdf_data2D,tmp_data2D,dim2)
    CASE(3)
       DEALLOCATE(netcdf_data3D,tmp_data3D,dim2,dim3)
    END SELECT

    IF (PRESENT(nX)) nX=ndim1
    IF (PRESENT(nY)) nX=ndim2
    IF (PRESENT(nZ)) nX=ndim3
    
  END SUBROUTINE AUX_DATA
  
  FUNCTION SOLAR_ZENITH_ANGLE(lat,day_of_year,local_solar_time,ngrids,onlyThese)
    !
    ! This subroutine calculated the solar zenith angle

    ! Source:
    ! The equations to calculate the solar zenith angle are copied from 
    ! 1) http://en.wikipedia.org/wiki/Equation_of_time (alternative calculation)
    ! and 
    ! 2) http://en.wikipedia.org/wiki/Solar_zenith_angle
    ! 
    ! Most comments at each calculation in the code are copy-pasted from
    ! the above wikipedia pages
    !
    ! Salomon.Eliasson@smhi.se

    IMPLICIT NONE

    INTEGER, INTENT(in)                   :: ngrids
    REAL(wp),INTENT(in)                   :: day_of_year
    REAL(wp),INTENT(in),DIMENSION(ngrids) :: lat
    REAL(wp),INTENT(in),DIMENSION(ngrids) :: local_solar_time
    REAL(wp),           DIMENSION(ngrids) :: solar_zenith_angle
    LOGICAL, INTENT(in),DIMENSION(ngrids), OPTIONAL :: onlyThese



    ! the Earth's orbirtal eccentricity
    REAL(wp), PARAMETER                       :: eoe  = 0.167_wp      
    ! The Earth's mean angular orbital velocity in degrees per day
    REAL(wp), PARAMETER                       :: W    = 360/365.24_wp 
    ! the Earth's obliquity angle in degrees
    REAL(wp), PARAMETER                       :: tilt = 23.44      
    REAL(wp), PARAMETER                       :: d180 = 180._wp
    ! convert radians to degrees
    REAL(wp), PARAMETER                       :: R_D  = d180/pi    
    ! convert degrees to radians
    REAL(wp), PARAMETER                       :: D_R  = pi/d180    

    REAL(wp) :: A, B, C, D
    REAL(wp) :: H(ngrids)
    REAL(wp) :: EoT                         ! equation of time
    REAL(wp) :: dcl                         ! declination
    REAL(wp) :: LST(ngrids)
    LOGICAL, DIMENSION(ngrids) :: data_mask

    IF (PRESENT(onlyThese)) THEN
       data_mask(1:ngrids) = onlyThese(1:ngrids)
    ELSE
       data_mask(1:ngrids) = .TRUE.
    END IF

    solar_zenith_angle(1:ngrids) = -999._wp

    ! ----------------
    ! EQUATION OF TIME
    ! http://en.wikipedia.org/wiki/Equation_of_time

    ! Angle the earth would move on its orbit at its average speed from
    ! the December solstice to date D.  D is the date, in days starting at
    ! zero on 1 January. 10 is the approximate number of days from the
    ! December solstice to 1 January.
    D = day_of_year - 1
    A = W * (D+10) !deg

    !Angle the Earth moves from the solstice to date D, including a first-order
    !correction for the Earth's orbital eccentricity, 0.167. The number 2 is
    !the number of days from January 1 to the date of the Earth's perihelion.

    B = A+ R_D* ( (360/pi)* (D_R*eoe)* SIN( D_R* (W*(D-2) ) ) )


    !C is the difference between the angles moved at mean speed, and at the
    !corrected speed projected onto the equatorial plane, and divided by 180 to
    !get the difference in "half turns". The subtraction gives the
    !conventional sign to the equation of time.

    C = (A - R_D* ATAN(TAN( D_R* B )/COS( D_R* tilt ) ) )/180

    !EoT is the equation of time in minutes. Subtracting nint(C) leaves a
    !small positive or negative fractional number of half turns, which is
    !multiplied by 720, the number of minutes (12 hours) that the Earth
    !takes to rotate one half turn relative to the Sun, to get the
    !equation of time.

    EoT = 750*(C-NINT(C))


    ! Declination
    !
    ! The value of B in the above calculation is an accurate value for the
    ! Sun's ecliptic longitude (shifted by 90 degrees), so the solar
    ! declination becomes readily available and is accurate to within a
    ! fraction of a degree.
    dcl = R_D* (-1._wp * ASIN(SIN(D_R* tilt)* COS(D_R* B) ))


    WHERE (data_mask)
       ! ------------
       ! Local solar time 
       ! add correction to local solar time
       LST(1:ngrids) = local_solar_time(1:ngrids) + EoT/60 

       !The Hour Angle converts the local solar time (LST) into the number of
       !degrees which the sun moves across the sky. By definition, the Hour Angle
       !is 0° at solar noon. Since the Earth rotates 15° per hour, each hour away
       !from solar noon corresponds to an angular motion of the sun in the sky of
       !15°. In the morning the hour angle is negative, in the afternoon the hour
       !angle is positive.
       H(1:ngrids) = 15*(LST(1:ngrids)-12);

       !% Solar Zenith Angle
       !
       ! The solar zenith angle, sza is estimated using results from spherical
       ! trigonometry. (from http://en.wikipedia.org/wiki/Solar_zenith_angle)
    END WHERE
    WHERE (data_mask)

       solar_zenith_angle(1:ngrids) = &
            R_D* ACOS(&
            SIN(D_R* lat(1:ngrids))*SIN(D_R* dcl)+&
            COS(D_R* lat(1:ngrids))*COS(D_R* dcl)*&
            COS(D_R*   H(1:ngrids))&
            )
    END WHERE

  END FUNCTION SOLAR_ZENITH_ANGLE

  FUNCTION CTP_TAU(ctp,tau,options,OD_lim,ncol,n_t,n_p)

    ! finds the number of hits for every CTP-tau bin.
    !
    ! salomon.eliasson@smhi.se

    IMPLICIT NONE

    ! IN
    !
    ! 'OD_lim' = Cloud optical depth threshold. This should be listed in
    !                  the namelist. All profiles with less than this
    !                  value are considered cloud free. The optical
    !                  depth is reduced by this amount too to take into
    !                  account the decreased sensitivity
    !
    ! 'tau'         = Cloud optical thickness for every subcolumn
    ! 'ctp'         = Cloud top pressure for every subcolumn [Pa]
    !
    ! CTP-TAU diagrams. 
    !
    ! 'n_t'         = number of optical thickness bins for P_tau histograms
    !                 (length(tbins)-1)
    ! 'n_p'         = number of pressure bins for P_tau histograms
    !                 (length(pbins)-1)

    TYPE(name_list), INTENT(in) :: options

    INTEGER,  INTENT(in)                  :: ncol, n_t, n_p
    REAL(wp), INTENT(in)                  :: OD_lim
    REAL(wp), INTENT(in), DIMENSION(ncol) :: ctp, tau

    INTEGER, DIMENSION(n_t+1,n_p) :: ctp_tau !number of boxes

    ! internal
    INTEGER :: ins, ixt, ixp, itau, ipres

    ctp_tau(1:n_t+1,1:n_p) = 0

    DO ins = 1,ncol
       ixt = 0
       ixp = 0

       IF (tau(ins) >= OD_lim ) THEN
          ! CHECK if there is enough optical depth to bin the data

          ! --- find correct cloud optical thickness bin
          IF (tau(ins) >= options%ctp_tau%tbin_edges(n_t+1) ) THEN
             ixt=n_t+1
          ELSE              
             DO itau = 1,n_t
                IF (tau(ins) >= options%ctp_tau%tbin_edges(itau) .AND. &
                     & tau(ins) < options%ctp_tau%tbin_edges(itau+1)) THEN

                   ixt = itau
                   EXIT
                END IF
             END DO
          END IF

          ! --- find correct pressure bin
          DO ipres = 1,n_p
             IF (ipres .EQ. n_p) THEN
                ixp = ipres
             ELSE IF (ctp(ins) > options%ctp_tau%pbin_edges(ipres) .AND. &
                  & ctp(ins) <= options%ctp_tau%pbin_edges(ipres+1)) THEN

                ixp = ipres
                EXIT
             END IF
          END DO
          IF (ixp .EQ. 0 .OR. ixt .EQ. 0) THEN
             PRINT *,"OD_lim=",OD_lim
             PRINT *,"COT = ",tau(ins)
             PRINT *,"CTP = ",ctp(ins)/100
             PRINT *, "ixp=",ixp,"p_bins=",options%ctp_tau%pbin_edges(1:n_p)
             PRINT *, "ixt=",ixp,"t_bins=",options%ctp_tau%tbin_edges(1:n_t+1)
             STOP "this shouldn't be"
          END IF

          ! fill ctp-tau by one
          ctp_tau(ixt,ixp) = ctp_tau(ixt,ixp) + 1
       END IF !--- there is a cloud that is thick enough
    END DO

    IF (options%dbg > 1 ) PRINT *, '--- Finished 2d hists'

  END FUNCTION CTP_TAU

  FUNCTION FIND_TEMPERATURE_INVERSIONS(sub,options,ngrids,nlev)&
       RESULT(inversions)
    !
    ! --- FIND THE TEMPERATURE INVERSIONS IN EACH PROFILE -----------------------
    ! 
    ! The purpose of this is to find all of the temperature inversions
    ! in each model grid. These are used in some of the cloud tests and
    ! the tropopause is also detected using this function. 
    !
    ! - The inversions are placed at the local minima in the profile. 
    ! - As a last step the highest inversion, assumed to be the tropopause is
    !   bumped up a level.
    ! - The tropopause must exist between 500-50hPa if none is found,
    !   the tropopause is placed at the level closest to the 50hPa level
    !
    ! Salomon.Eliasson@smhi.se

    IMPLICIT NONE

    TYPE(name_list), INTENT(in) :: options
    INTEGER, INTENT(in) :: ngrids,nlev


    ! "Inversions" = indices where there are temperature inversions
    !                (where, from the surface up, the temperature does
    !                not decrease compared to the temperature at i-1 level)

    TYPE(subset), INTENT(inout) :: sub

    ! internal variables
    INTEGER :: inl

    ! Relaxation for clouds in the tropopause region
    REAL(wp), PARAMETER :: MaximumTropopausePressure = 50000 ![Pa]
    REAL(wp), PARAMETER :: MinimumTropopausePressure = 5000  ![Pa]

    INTEGER, DIMENSION(ngrids) :: invInd, trop_lev,tropopause_count
    REAL(wp),DIMENSION(ngrids) :: inversion_pressure
    LOGICAL, DIMENSION(ngrids,nlev) :: previousNotInversion
    INTEGER :: clara_tropopauseThickness
    INTEGER :: d1
    ! out
    INTEGER, DIMENSION(ngrids,nlev) :: inversions

    trop_lev                    = 1 ! need to put it at the top for if no tropopause is found
    inversion_pressure          = -999._wp
    inversions(1:ngrids,1:nlev) = 0
    clara_tropopauseThickness   = MAX(2,nlev/30)
    invInd = 0

    previousNotInversion  = .TRUE.

    ! Loop from the surface up to the tropopause
    tropopause_count = 0
    
    DO d1 = 1,ngrids
       DO inl = nlev,2,-1
          
          !IF (sub%data_mask(d1,inl)) CYCLE
          IF ( sub%Tcorr(d1,inl) .GE. sub%Tcorr(d1,inl+1) ) THEN
             ! 
             ! You are in an inversion
             !
             IF (previousNotInversion(d1,inl)) THEN
                ! You have found the first or another
                ! inversion. Increase the index for the inversion
                ! array if the level below this one was not an
                ! inversion (as long as we are still in the
                ! tropopause), save the model level as an inversion level
                
                invInd(d1) = invInd(d1)+1
                inversions(d1,invInd(d1)) = inl
                inversion_pressure(d1) = sub%p_mid(d1,inl)
             END IF
             ! Since this layer is an inversion, the following statement should
             ! be false for the NEXT iteration (higher altitude)
             previousNotInversion(d1,inl-1) = .FALSE.
             ! start looking for the tropopause higher than 500hPa
             IF (sub%p_mid(d1,inl).LE.MaximumTropopausePressure) THEN !500 hPa
                tropopause_count(d1) = tropopause_count(d1)+1
             END IF
          ELSE
             tropopause_count(d1) = 0
             previousNotInversion(d1,inl-1) = .TRUE.
          END IF
          ! Find the tropopause: CLARA-A2's reasoning for the
          ! tropopause is: to find the first 2 or 3 consecutive
          ! inversion layers (depends on the vertical resolution of
          ! the model) with pressures less than 500hPa
          IF (sub%p_mid(d1,inl) .LT. MinimumTropopausePressure) THEN
             ! stop looking for inversions. Call the highest
             ! inversion the tropopause, or if none where found
             ! set the tropopause to this level
             trop_lev(d1) = MAXVAL([inl,MINVAL(inversions(d1,1:nlev))])
             
             ! exit this level loop
             EXIT
          ENDIF
          IF (tropopause_count(d1) .EQ. clara_tropopauseThickness) THEN
             ! go back to the first of the inversion layers that made up this tropopause
             trop_lev(d1) = inl+clara_tropopauseThickness-1
          END IF
       END DO ! nlevel
    END DO ! ngrids
    
!!$    WHERE( .NOT.sub%data_mask(1:ngrids,1) .AND.&
!!$         ( invInd(1:ngrids) .EQ. 0 .OR. & !no inversions
!!$         inversion_pressure(1:ngrids).GT.MaximumTropopausePressure))!too low
       
    WHERE( invInd(1:ngrids) .EQ. 0 .OR. & !no inversions
         inversion_pressure(1:ngrids).GT.MaximumTropopausePressure)!too low

       inversions(1:ngrids,1)=trop_lev(1:ngrids)

    END WHERE

    IF (options%dbg > 0) PRINT *, '--- Finished model temperature inversions'

  END FUNCTION FIND_TEMPERATURE_INVERSIONS

  SUBROUTINE CONFORM_GRID_DATA_TO_MODEL(aux,flags,nlon,nlat,lon,lat,data2,&
       ndim3,data3)

    ! reorder data to match the model data if need be

    IMPLICIT NONE

    TYPE(model_aux), INTENT(in)      :: aux ! model data
    TYPE(flag), INTENT(out)          :: flags

    ! aux data
    INTEGER,  INTENT(in)             :: nlon,nlat
    INTEGER,  INTENT(in),OPTIONAL    :: ndim3
    REAL(wp), INTENT(inout)          :: lon(nlon),lat(nlat)
    REAL(wp), INTENT(inout),OPTIONAL :: data2(nlon,nlat)
    REAL(wp), INTENT(inout),OPTIONAL :: data3(:,:,:)

    ! temporary arrays
    REAL(wp) :: out_lon(nlon),out_lat(nlat)

    REAL(wp), ALLOCATABLE ::out_data2(:,:),  tmp_out_data2(:,:)
    REAL(wp), ALLOCATABLE ::out_data3(:,:,:),tmp_out_data3(:,:,:)    
    REAL(wp) :: tmp_out_lon(nlon) 
    LOGICAL  :: model_lat_descend,data_lat_descend,model_lon_180,data_lon_180
    LOGICAL  :: loncond1,loncond2,latcond1,latcond2,workOnLat,workOnLon

    REAL(wp), PARAMETER :: BigNum= 999._wp
    REAL(wp), PARAMETER :: fill  =-999._wp
    INTEGER             :: i
    INTEGER             :: ind(1)
    LOGICAL             :: dim2

    
    ! initialise all
    out_lon      = fill
    out_lat      = fill

    IF (PRESENT(data2)) THEN
       dim2=.TRUE.
       ALLOCATE(out_data2(nlon,nlat),&
            tmp_out_data2(nlon,nlat))
       out_data2     = fill
       tmp_out_data2 = fill
    ELSE
       dim2=.FALSE.
       ALLOCATE(out_data3(nlon,nlat,ndim3),&
            tmp_out_data3(nlon,nlat,ndim3))
       out_data3     = fill
       tmp_out_data3 = fill
    END IF
    
    tmp_out_lon  = fill

    latcond1=.FALSE.
    latcond2=.FALSE.
    workOnLat=.FALSE.
    loncond1=.FALSE.
    loncond2=.FALSE.
    workOnLon=.FALSE.
    model_lat_descend=.FALSE.
    data_lat_descend=.FALSE.
    model_lon_180=.FALSE.
    data_lon_180=.FALSE.

    ! -------
    ! figure out if the model has ascending or decending lats/lons
    ! figure out if the detection limit -||-
    model_lat_descend = aux%lat(1,2) < aux%lat(1,1)
    data_lat_descend  = lat(2)       < lat(1)
    model_lon_180     = aux%lon(1,1) < 0
    data_lon_180      = lon(1)       < 0

    flags%model_lat_descend = model_lat_descend
    flags%data_lat_descend  = data_lat_descend 
    flags%model_lon_180     = model_lon_180    
    flags%data_lon_180      = data_lon_180     
    flags%modelIsCoarser    = (nlon.GT.aux%nlon .AND. nlat.GT.aux%nlat)
    flags%modelIsFiner      = (nlon.LT.aux%nlon .AND. nlat.LT.aux%nlat)

    ! LATITUDES
    latcond1  = (model_lat_descend .AND. .NOT. data_lat_descend)
    latcond2  = (.NOT.model_lat_descend .AND. data_lat_descend)
    workOnLat = latcond1 .OR. latcond2

    IF (workOnLat) THEN
       PRINT *, "--- reordering aux data for latitudes"
       out_lat  = lat(nlat:1:-1)
       IF (dim2) THEN
          out_data2 = data2(:,nlat:1:-1)
       ELSE
          out_data3 = data3(:,nlat:1:-1,1:ndim3)
       END IF
    ELSE
       out_lat=lat
       IF (dim2) THEN
          out_data2 = data2
       ELSE
          out_data3 = data3
       END IF
    END IF

    ! LONGITUDES
    loncond1  = model_lon_180 .AND. .NOT.data_lon_180
    loncond2  = .NOT.model_lon_180 .AND. data_lon_180
    workOnLon = loncond1 .OR. loncond2

    IF (workOnLon) THEN
       IF (loncond1) tmp_out_lon = MERGE(lon-360, lon, lon .GT. 180)
       IF (loncond2) tmp_out_lon = MERGE(lon+360, lon, lon .LT. 0  )

       ! my funky way of sorting
!
       IF (dim2) THEN
          tmp_out_data2=out_data2(1:nlon,1:nlat)
       ELSE
          tmp_out_data3=out_data3(1:nlon,1:nlat,1:ndim3)
       END IF
       
       DO i = 1,nlon
          ind = MINLOC(tmp_out_lon)
          out_lon(i)=tmp_out_lon(ind(1))
          tmp_out_lon(ind) = BigNum ! overwrite the smallest number

          IF (dim2) THEN
             out_data2(i,1:nlat)=tmp_out_data2(ind(1),1:nlat)
          ELSE
             out_data3(i,1:nlat,1:ndim3)=tmp_out_data3(ind(1),1:nlat,1:ndim3)
          END IF
       END DO
       !
       ! -----

       PRINT *, "--- changing longitude regime to match the model"

    ELSE
       out_lon=lon
    END IF

    IF (dim2) THEN
       data2 = out_data2
       DEALLOCATE(out_data2,tmp_out_data2)
    ELSE
       data3 = out_data3
       DEALLOCATE(out_data3,tmp_out_data3)
    END IF
    
    lon  = out_lon
    lat  = out_lat

  END SUBROUTINE CONFORM_GRID_DATA_TO_MODEL
END MODULE auxiliary_functions
