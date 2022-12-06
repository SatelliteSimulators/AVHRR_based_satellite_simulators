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
  USE MY_NETCDFTOOLS,             ONLY: &
       ADD_VARIABLE_ATTRIBUTES,         &
       CHECK
  USE NETCDF,                     ONLY:&
       NF90_CHAR,                      &
       NF90_CLOSE,                     &
       NF90_DEF_VAR,                   &
       NF90_ENDDEF,                    &
       NF90_GET_VAR,                   &
       NF90_INQUIRE,                   &
       NF90_INQUIRE_DIMENSION,         &
       NF90_INQ_VARID,                 &
       NF90_OPEN,                      &
       NF90_PUT_ATT,                   &
       NF90_PUT_VAR,                   &
       NF90_REDEF,                     &              
       NF90_NOWRITE              
  USE NAMELIST_INPUT,            ONLY: NAME_LIST
  USE SIMULATOR_INPUT_VARIABLES, ONLY: SUBSET

  IMPLICIT NONE

  PUBLIC :: AUX_DATA,              &
       CTP_TAU,                    &
       FIND_TEMPERATURE_INVERSIONS,&
       HANDLE_VARIABLE_ATTRIBUTES, &
       PUTVAR,                     &
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

  SUBROUTINE HANDLE_VARIABLE_ATTRIBUTES(ncid,variable,ndim,id,varid,O,M)

    ! Repository of variable names for use of saving to netcdf file
    !
    !
    !
    ! Salomon.Eliasson@smhi.se
    !

    USE DATA_CHECK,       ONLY: height_max,p_tau_max,P_max,tau_max,T_min,T_max
    USE MODEL_INPUT,      ONLY: MODEL_TYPE,MODEL_AUX
    USE NAMELIST_INPUT,   ONLY: NAME_LIST


    IMPLICIT NONE

    INTEGER, INTENT(in)                  :: ncid
    CHARACTER(len=*),INTENT(in)          :: variable
    INTEGER, INTENT(out)                 :: varid
    TYPE(name_list),INTENT(in)           :: O
    TYPE(model_type),OPTIONAL,INTENT(in) :: M
    INTEGER, INTENT(in)                  :: ndim
    INTEGER, INTENT(in)                  :: id(ndim)

    !local
    TYPE(model_aux)                      :: A
    REAL(4)                              :: fvr = -999.0 
    INTEGER, PARAMETER                   :: fvi = -999
    LOGICAL                              :: isinteger
    INTEGER                              :: varid_rotpole,i,nchan
    CHARACTER(len=1000)                  :: units,description,long_name,&
         string1,string2
    CHARACTER(len=20)                    :: coordinates
    CHARACTER(len=12)                    :: grid_mapping
    CHARACTER(len=20)                    :: calendar
    CHARACTER(len=1)                     :: axis
    CHARACTER(len=20)                    :: standard_name
    REAL(4) :: scale_factor,add_offset,valid_min,valid_max 

4   FORMAT(a,1x,i3,a,f5.2)
    IF (PRESENT(M)) A=M%aux
    varid       = 0
    isinteger   = .FALSE.
    grid_mapping= ''

    coordinates = ''

    ! so far everything is stored to scale and is not offseted
    scale_factor= 1._4
    add_offset  = 0._4
    nchan       = O%sim%nchannels
    valid_min   = fvr
    valid_max   = fvr 
    IF (O%model.EQ.'racmo') THEN
       IF (ndim.GT.1) THEN
          IF (variable.NE.'time_bnds'.AND.variable.NE.'lon'.AND.variable.NE.'lat') THEN
             grid_mapping= 'rotated_pole'
             coordinates = 'lon lat'   
          END IF
       END IF
    END IF

    IF (O%sim%doRTTOV) THEN
       WRITE(string1,'(I2)') O%sim%sensor
       DO i = 1,nchan
          IF (i==1) THEN
             WRITE(string2,'(I2)') O%sim%channels(i)
          ELSE
             WRITE(string2,'(A,",",I2)') TRIM(string2),O%sim%channels(i)
          END IF
       END DO
    END IF

    CALL check( nf90_redef(ncid) )

    SELECT CASE (variable)
    CASE ('areacella')
       units = "m2"
       WRITE(description,'(A)') "" 
       WRITE(long_name,'(A)') "Atmosphere Grid-Cell Area"
    CASE ('block1')
       units = ""
       isinteger = .TRUE.
       WRITE(description,'(A)') ""
       WRITE(long_name,'(A)') "GRIB Definition Block1"
    CASE ('block2')
       units = ""
       isinteger = .TRUE.
       WRITE(description,'(A)') ""
       WRITE(long_name,'(A)') "GRIB Definition Block2"
    CASE ('dtg')
       units = "yyyymmddhh"
       isinteger = .TRUE.
       WRITE(description,'(A)') "dtg"
       WRITE(long_name,'(A)') "Verifying Datum-Time Group"
    CASE ('lat')
       units = "degrees_north"
       WRITE(description,'(A)') "North of the equator" 
       WRITE(long_name,'(A)') "latitude"
    CASE ('lon')
       units = "degrees_east"
       WRITE(description,'(A)') "East from Greenwich"
       WRITE(long_name,'(A)') "longitude"
    CASE ('lsm')
       units =  "fraction"
       WRITE(description,'(A)') "Grid fractional cover of land"
       WRITE(long_name,'(A)') "Grid fractional cover of land"
    CASE ('rlon')
       units = "degrees"
       description = ''
       WRITE(long_name,'(A)') "longitude in rotated pole grid"
    CASE ('rlat')
       units = "degrees"
       description = ''
       WRITE(long_name,'(A)') "latitude in rotated pole grid"
    CASE ('hist2d_ctp_bin_border')
       units = "Pa"
       WRITE(description,'(A)') "Box boundaries of the vertical profile of atmospheric &
            &pressure, used to make CTP -tau histograms"
       WRITE(long_name,'(A)') "Box boundaries of the vertical profile of atmospheric pressure"
    CASE ('hist2d_cot_bin_border')
       units = "-"
       WRITE(description,'(A)') &
            "optical thickness box boundaries used to make CTP -tau histograms"
       WRITE(long_name,'(A)') "optical thickness box boundaries" 
    CASE ('hist2d_ctp_bin_centre')
       units = "Pa"
       WRITE(description,'(A)') "Box centres of the vertical profile of atmospheric &
            &pressure, used to make CTP -tau histograms"
       WRITE(long_name,'(A)') "Box centres of the vertical profile of atmospheric pressure"
    CASE ('hist2d_cot_bin_centre')
       units = "-"
       WRITE(description,'(A)') &
            "optical thickness box centres used to make CTP -tau histograms"
       WRITE(long_name,'(A)') "optical thickness box centres"
    CASE ('hist_phase')
       units = "-"
       WRITE(description,'(A)') &
            "cloud phase dimension for hist2d_cot_ctp. 0=ice,1=liq"
       WRITE(long_name,'(A)') "cloud phase dimension"
    CASE ('POD_layers')
       units = ""
       WRITE(description,'(A)') "Probability of detection in optical depth bins"
       WRITE(long_name,'(A)') "probability of detection"
       valid_min=0
       valid_max=1
    CASE ('POD_tau_bin_centers')
       units = ""
       WRITE(description,'(A)') "Optical depth bin centers of bins between which POD and FAR are calculated"
       WRITE(long_name,'(A)') "COT bin centers used for POD layers"
       valid_min=0
       valid_max=5
    CASE ('POD_tau_bin_edges')
       units = ""
       WRITE(description,'(A)') "Optical depth bin edges between which POD and FAR are calculated"
       WRITE(long_name,'(A)') "COT bin edges used for POD layers"
       valid_min=0
       valid_max=9999
    CASE ('solzen')
       units =  "deg"
       WRITE(description,'(A)') "Grid Solar Zenith Angle"
       WRITE(long_name,'(A)') "Grid Solar Zenith Angle"
    CASE ('time_of_day')
       units =  "hr"
       WRITE(description,'(A)') "Overpass Time of Day"
       WRITE(long_name,'(A)') "overpass time of day"
    CASE ('time')
       units = "days since 1970-01-01 00:00:00.0"
       WRITE(description,'(A)') "time"
       WRITE(long_name,'(A)') "time"
       WRITE(calendar,'(A)') "standard"
    CASE ('time_bnds')
       units = "days since 1970-01-01 00:00:00.0"
       WRITE(description,'(A)') "time_bnds"
       WRITE(long_name,'(A)') "time_bnds"

       !
       ! AUXILIARY
       ! ------

       ! -------------------------------
       ! SIMULATED VARIABLES
       ! 
    CASE ('albedo')
       units='fraction'
       WRITE(description,'(A)') "simulated cloudy albedo"
       WRITE(long_name,'(A)') "simulated cloudy albedo during daytime" 
       valid_min=0
       valid_max=1
    CASE ('cfc')
       units = "fraction"
       WRITE(description,'(A)') "simulated total cloud cover"
       WRITE(long_name,'(A)') "total cloud cover" 
       valid_min=0
       valid_max=1
    CASE ('cfc_day')
       units = "fraction"
       WRITE(description,'(A)') "simulated total cloud cover during daytime conditions only"
       WRITE(long_name,'(A)') "total cloud cover day clouds"
       valid_min=0
       valid_max=1
    CASE ('cfc_low')
       units = "fraction"
       WRITE(description,'(A)') "simulated total cloud cover for&
            &clouds with CTP greater than 680 hPa"
       WRITE(long_name,'(A)') "total cloud cover low clouds"
       valid_min=0
       valid_max=1
    CASE ('cfc_mid')
       units = "fraction"
       WRITE(description,'(A)') "simulated total cloud cover for&
            &clouds with CTP less than 680 hPa and greater than 440 hPa"
       WRITE(long_name,'(A)') "total cloud cover mid clouds"
       valid_min=0
       valid_max=1
    CASE ('cfc_high')
       units = "fraction"
       WRITE(description,'(A)') "simulated total cloud cover for&
            &clouds with CTP less than 440 hPa"
       WRITE(long_name,'(A)') "total cloud cover high clouds"
       valid_min=0
       valid_max=1
    CASE ('cla_vis006')
       units = "%"
       WRITE(description,'(A)')&
            "Simulated cloud albedo at 0.6 micron"
       WRITE(long_name,'(A)') "Simulated grid average cloud albedo" 
       valid_min=0
       valid_max=100
    CASE ('cth')
       units = 'm'
       WRITE(description,'(A)') "simulated cloud top height"
       WRITE(long_name,'(A)') "cloud top height" 
       valid_min=1
       valid_max=height_max
    CASE ('cth_corrected')
       units = 'm'
       WRITE(description,'(A,1x,F5.2,A,F5.2)')&
            "simulated CTH based on the CLOUD_CCI method (height corrected) of finding the cloud&
            & top (where tau=1), but also where tau >",O%cloudMicrophys%tau_min
       WRITE(long_name,'(A)') "corrected cloud top equivalent cloud top height" 
       valid_min=1
       valid_max=height_max
    CASE ('ctp')
       units = 'Pa'
       WRITE(description,'(A)') "simulated cloud top pressure. &
            &Derived from linear-averaging of the sub-grid CTP"
       WRITE(long_name,'(A)') "cloud top pressure"
       valid_min=0
       valid_max=P_max
    CASE ('ctp_log')
       units = 'Pa'
       WRITE(description,'(A)') "simulated cloud top pressure. &
            &Derived from log-averaging of the sub-grid CTP"
       WRITE(long_name,'(A)') "cloud top pressure" 
       valid_min=1
       valid_max=P_max
    CASE ('ctp_corrected')
       units = 'Pa'
       WRITE(description,'(A,1x,F5.2,A,F5.2)')&
            "simulated CTP based on the CLOUD_CCI method (height corrected) &
            &of finding the cloud top (where tau=0.3), but also where tau >",O%cloudMicrophys%tau_min
       WRITE(long_name,'(A)') "corrected cloud top equivalent cloud top pressure" 
       valid_min=1
       valid_max=P_max
    CASE ('ctt')
       units = 'K'
       WRITE(description,'(A)') "simulated cloud top temperature"
       WRITE(long_name,'(A)') "cloud top temperature"
       valid_min=T_min
       valid_max=T_max
    CASE ('ctt_corrected')
       units = 'K'
       WRITE(description,'(A,1x,F5.2,A,F5.2)')&
            "simulated CTT based on the CLOUD_CCI method (height corrected) of finding the cloud&
            & top (where tau=0.3), but also where tau >",O%cloudMicrophys%tau_min
       WRITE(long_name,'(A)') "corrected cloud top equivalent cloud top temperature"
       valid_min=T_min
       valid_max=T_max
    CASE ('icf')
       units = "fraction"
       WRITE(description,'(A)') "simulated fraction of clouds that contain more ice&
            &than liquid particles"
       WRITE(long_name,'(A)') "fraction of grid covered by ice clouds"
       valid_min=0
       valid_max=1
    CASE ('ireff','cer_ice','ref_ice')
       units = "micron"
       WRITE(description,'(A)') "simulated ice effective radius"
       WRITE(long_name,'(A)') "ice effective radius"
       valid_min=1
       valid_max=155
    CASE ('itau','cot_ice')
       units = "-"
       WRITE(description,'(A)') "simulated grid average cloud optical thickness for ice clouds"
       WRITE(long_name,'(A)') "Grid average optical thickness ice cloud" 
       valid_min=0
       valid_max=tau_max
    CASE ('iwp')
       units = "kg/m^2"
       WRITE(description,'(A)') "simulated grid average ice water path"
       WRITE(long_name,'(A)') "ice water path" 
       valid_min=0
       valid_max=1e6
    CASE ('lcf')
       units = "fraction"
       WRITE(description,'(A)') "simulated fraction of clouds that contain more liquid&
            &than ice particles"
       WRITE(long_name,'(A)') "fraction of grid covered by liquid clouds"
       valid_min=0
       valid_max=1
    CASE ('lreff','cer_liq','ref_liq')
       units = "micron"
       WRITE(description,'(A)') "simulated liquid effective radiusd"
       WRITE(long_name,'(A)') "liquid effective radius" 
       valid_min=1
       valid_max=155
    CASE ('ltau','cot_liq')
       units = "-"
       WRITE(description,'(A)') "simulated grid average cloud optical thickness for liquid clouds"
       WRITE(long_name,'(A)') "Grid average optical thickness liquid cloud" 
       valid_min=1
       valid_max=tau_max
    CASE ('lwp')
       units = "kg/m^2"
       WRITE(description,'(A)') "grid average liquid water path"
       WRITE(long_name,'(A)') "liquid water path" 
       valid_min=1
       valid_max=1e6
    CASE ('hist2d_cot_ctp')
       units = "unitless"
       isinteger = .TRUE.
       WRITE(description,4) "CTP--tau hits (lon,lat,n*tau bins,n*pressure bins) &
            & based on ",O%ncols,"sub-grids in each grid box and only includ&
            &ing clouds that have tau > ",O%cloudMicrophys%tau_min
       WRITE(long_name,'(A)') "Cloud top pressure- cloud optical thickness histograms"
       valid_min=0
       valid_max=1e6
    CASE ('tau','cot')
       units = "-"
       WRITE(description,'(A)') "grid average cloud optical thickness"
       WRITE(long_name,'(A)') "optical thickness" 
       valid_min=0
       valid_max=tau_max
    CASE ('tau_subcolumn')
       units = "-"
       WRITE(description,'(A)') "tau for individual subcolumns"
       WRITE(long_name,'(A)') "cloud optical thickness for individual subcolumns"
       valid_min=0
       valid_max=tau_max
       ! This has changed its meaning after the introduction of the ISCCP simulator
       !CASE ('Tb')
       !   units = "K"
       !   WRITE(description,'(5(A,1x))') "Tb from RTTOV, sensor=",TRIM(string1),&
       !        &"channels =",TRIM(string2),". See RTTOV user guide for details"
       !   WRITE(long_name,'(A)') "Brightness temperature"
       !   valid_min=T_min
       !   valid_max=T_max
       !CASE ('Tb_clr')
       !   units = "K"
       !   WRITE(description,'(5(A,1x))') "Box average clear sky Tb from RTTOV, &
       !        &sensor=",TRIM(string1),"channels =",TRIM(string2),&
       !        ". See RTTOV user guide for details"
       !   WRITE(long_name,'(A)') "Brightness temperature from clear sky"
       !   valid_min=T_min
       !   valid_max=T_max
       !CASE ('Tb_subcolumn')
       !   units = "K"
       !   WRITE(description,'(5(A,1x))') "Tb for individual subcolumns. &
       !        &Tb from RTTOV, sensor=",TRIM(string1),&
       !        &"channels =",TRIM(string2),". See RTTOV user guide for details"
       !   WRITE(long_name,'(A)') "Brightness temperature  for individual subcolumns"
       !   valid_min=T_min
       !   valid_max=T_max
       !
       ! END Simulated fields
       !------------------------

       !------------------------
       ! DIRECT MODEL FIELDS
       !
    CASE ('CI')
       units = "-"
       WRITE(description,'(A)') "Sea ice fraction [0-1]. code = 31"
       WRITE(long_name,'(A)') "Sea ice fraction"
       valid_min=0
       valid_max=1
    CASE ('ireff3D')
       units = "micron"
       WRITE(description,'(A)') "model layered ice effective radius"
       WRITE(long_name,'(A)') "model ice effective radius"
       valid_min=1
       valid_max=155
    CASE ('lreff3D')
       units = "micron"
       WRITE(description,'(A)') "model layerd liquid effective radius"
       WRITE(long_name,'(A)') "liquid effective radius"
       valid_min=1
       valid_max=155
    CASE ('SKT')
       units = "K"
       WRITE(description,'(A)') "Model Skin temperature. Temperature of the surface skin &
            & (radiative surface temperature). Before 01/10/2008, the &
            & skin temperature was equal to the bulk SST over the ocean. paramId=235"
       WRITE(long_name,'(A)') "model skin temperature" 
       valid_min=T_min
       valid_max=T_max
    CASE ('TCC')
       units = "fraction"
       WRITE(description,'(A)') "model total cloud cover"
       WRITE(long_name,'(A)') description
       valid_min=0
       valid_max=1
    CASE ('TCWV')
       units = "kg m**-2"
       WRITE(description,'(A)') "model total column water vapour. paramid = 137"
       WRITE(long_name,'(A)') "model total column water vapour"
       valid_min=0
       valid_max=100
       !
       ! END Direct model fields
       !------------------------

    CASE DEFAULT
       PRINT '(3(a,1x))',"variable:",variable,"is not setup"
       STOP "stopped in addVariable"
    END SELECT

    CALL add_variable_attributes(ncid,varid,variable,id(1:ndim), ndim,&
         add_offset=add_offset,coordinates=coordinates,description=description,&
         fillvalue_i=fvi,fillvalue_r=fvr,grid_mapping=grid_mapping,&
         isinteger=isinteger,long_name=long_name,scale_factor=scale_factor,&
         unit=units,valid_min=valid_min,valid_max=valid_max,dbg=O%dbg)

    ! And a little bit more....
    CALL check( nf90_redef(ncid) )

    SELECT CASE (variable)
    CASE ('time')
       WRITE(calendar,'(A)') "standard"
       CALL check( nf90_put_att(ncid,varid, "bounds", "time_bnds"))
       CALL check( nf90_put_att(ncid,varid, "calendar", calendar))
    CASE ('rlat')
       WRITE(axis,'(A)')"Y"
       WRITE(standard_name,'(A)')"grid_latitude"
       CALL check( nf90_put_att(ncid,varid, "axis", axis))
       CALL check( nf90_put_att(ncid,varid, "standard_name", standard_name))
       ! ... abuse invoking 'rlat' to specify 'rotated_pole' which contains only attributes
       CALL check( nf90_def_var(ncid=ncid,name='rotated_pole',xtype=NF90_CHAR,varID=varid_rotpole))
       CALL check( nf90_put_att(ncid,varid_rotpole,"grid_mapping_name",A%rotgrid))
       CALL check( nf90_put_att(ncid,varid_rotpole,"proj",A%rotproj))
       CALL check( nf90_put_att(ncid,varid_rotpole,"grid_north_pole_latitude",A%float_rotpolat))
       CALL check( nf90_put_att(ncid,varid_rotpole,"grid_north_pole_longitude",A%float_rotpolon))
    CASE ('rlon')
       WRITE(axis,'(A)')"X"
       WRITE(standard_name,'(A)')"grid_longitude"
       CALL check( nf90_put_att(ncid,varid, "axis", axis))
       CALL check( nf90_put_att(ncid,varid, "standard_name", standard_name))
    CASE ('lon')
       WRITE(standard_name,'(A)')"longitude"
       CALL check( nf90_put_att(ncid,varid, "standard_name", standard_name))
    CASE ('lat')
       WRITE(standard_name,'(A)')"latitude"
       CALL check( nf90_put_att(ncid,varid, "standard_name", standard_name))
    END SELECT

    CALL check( nf90_enddef(ncid) )

  END SUBROUTINE HANDLE_VARIABLE_ATTRIBUTES

  SUBROUTINE putVar(ncid,varid,A,data1,data2,data3,data4)

    IMPLICIT NONE

    INTEGER,         INTENT(in) :: ncid,varid
    TYPE(model_aux), INTENT(in) :: A
    REAL(wp),        INTENT(in),OPTIONAL :: data1(:),data2(:,:)
    REAL(wp),        INTENT(in),OPTIONAL :: data3(:,:,:),data4(:,:,:,:)
    INTEGER,ALLOCATABLE                  :: sz(:)
    INTEGER :: ndim

    IF (PRESENT(data1)) ndim=1
    IF (PRESENT(data2)) ndim=2
    IF (PRESENT(data3)) ndim=3
    IF (PRESENT(data4)) ndim=4

    ALLOCATE(sz(ndim))

    SELECT CASE (ndim)
    CASE(1)
       sz(1:ndim)=SHAPE(data1)
       IF (sz(1)==A%ngrids) THEN
          CALL check(nf90_put_var(ncid,varid,&
               RESHAPE(data1,(/A%nlon,A%nlat/))))
          ndim=ndim+1
       ELSE
          CALL check(nf90_put_var(ncid,varid,data1))
       END IF
    CASE(2)
       sz(1:ndim)=SHAPE(data2)
       IF (sz(1)==A%ngrids) THEN
          CALL check(nf90_put_var(ncid,varid,&
               RESHAPE(data2,(/A%nlon,A%nlat,sz(2)/))))
          ndim=ndim+1
       ELSE
          CALL check(nf90_put_var(ncid,varid,data2))
       END IF
    CASE(3)
       sz(1:ndim)=SHAPE(data3)
       IF (sz(1)==A%ngrids) THEN
          CALL check(nf90_put_var(ncid,varid,&
               RESHAPE(data3,(/A%nlon,A%nlat,sz(2),sz(3)/))))
          ndim=ndim+1
       ELSE
          CALL check(nf90_put_var(ncid,varid,data3))
       END IF
    CASE(4)
       sz(1:ndim)=SHAPE(data4)
       IF (sz(1)==A%ngrids) THEN
          CALL check(nf90_put_var(ncid,varid,&
               RESHAPE(data4,(/A%nlon,A%nlat,sz(2),sz(3),sz(4)/))))
          ndim=ndim+1
       ELSE
          CALL check(nf90_put_var(ncid,varid,data4))
       END IF
    END SELECT

    DEALLOCATE(sz)    

99  FORMAT('Added ',i1,' dimensional data')
    !PRINT 99,ndim

  END SUBROUTINE putVar

  FUNCTION SOLAR_ZENITH_ANGLE(lat,day_of_year,local_solar_time,ngrids)
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

    ! The Earth's mean angular orbital velocity in degrees per day
    REAL(wp) :: W
    REAL(wp) :: A, B, C, D
    REAL(wp) :: H(ngrids)
    REAL(wp) :: EoT                         ! equation of time
    REAL(wp) :: dcl                         ! declination
    REAL(wp) :: LST(ngrids)

    LST                          = -999._wp
    solar_zenith_angle(1:ngrids) = -999._wp

    ! ----------------
    ! EQUATION OF TIME
    ! http://en.wikipedia.org/wiki/Equation_of_time

    W = 360/365.24_wp 
    ! W is the Earth's mean angular orbital velocity in degrees per day.


    D = day_of_year - 1
    A = W * (D+10)

    ! D is the date, in days starting at zero on 1 January (i.e. the
    ! days part of the ordinal date minus 1). 10 is the approximate
    ! number of days from the December solstice to 1 January. A is the
    ! angle the earth would move on its orbit at its average speed
    ! from the December solstice to date D.

    B = A + 360/pi * 0.0167 * SIND(W*(D-2))

    !B is the angle the Earth moves from the solstice to date D, including a first-order
    !correction for the Earth's orbital eccentricity, 0.0167. The number 2 is
    !the number of days from January 1 to the date of the Earth's perihelion.

    C = (A - ATAND(TAND(B)/COSD(23.44)))/180
    ! C is the difference between the angles moved at mean speed, and at
    ! the corrected speed projected onto the equatorial plane, and divided
    ! by 180 to get the difference in "half turns". The value 23.44° is
    ! the obliquity (tilt) of the Earth's axis. The subtraction gives the
    ! conventional sign to the equation of time. For any given value of x,
    ! arctan x (sometimes written as tan−1 x) has multiple values,
    ! differing from each other by integer numbers of half turns. The
    ! value generated by a calculator or computer may not be the
    ! appropriate one for this calculation. This may cause C to be wrong
    ! by an integer number of half turns. The excess half turns are
    ! removed in the next step of the calculation to give the equation of
    ! time:

    EoT = 750*(C-NINT(C))

    !EoT is the equation of time in minutes. Subtracting nint(C) leaves a
    !small positive or negative fractional number of half turns, which is
    !multiplied by 720, the number of minutes (12 hours) that the Earth
    !takes to rotate one half turn relative to the Sun, to get the
    !equation of time.

    ! Declination
    !
    ! The value of B in the above calculation is an accurate value for the
    ! Sun's ecliptic longitude (shifted by 90 degrees), so the solar
    ! declination becomes readily available and is accurate to within a
    ! fraction of a degree.
    dcl =-1._wp * ASIND(SIND(23.44)*COSD(B))

    ! ------------
    ! Local solar time 
    ! add correction to local solar time
    
    LST(1:ngrids) = local_solar_time(1:ngrids) + EoT/60 

    ! The Hour Angle converts the local solar time (LST) into the
    ! number of degrees which the sun moves across the sky. By
    ! definition, the Hour Angle is 0° at solar noon (hence LST-12h).
    ! Since the Earth rotates 15° per hour, each hour away
    ! from solar noon corresponds to an angular motion of the sun in
    ! the sky of 15°. In the morning the hour angle is negative, in
    ! the afternoon the hour angle is positive.
    H(1:ngrids) = 15*(LST(1:ngrids)-12);

    !% Solar Zenith Angle
    !
    ! The solar zenith angle, sza is estimated using results from spherical
    ! trigonometry. (from http://en.wikipedia.org/wiki/Solar_zenith_angle)

    solar_zenith_angle(1:ngrids) = &
         ACOSD(&
         SIND(lat(1:ngrids))*SIND(dcl)+&
         COSD(lat(1:ngrids))*COSD(dcl)*COSD(H(1:ngrids))&
         )

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

    PRINT *, "--- Finding temperature inversions"
    
    trop_lev                    = 1 ! put at the top if no tropopause is found
    inversion_pressure          = -999._wp
    inversions(1:ngrids,1:nlev) = 0
    clara_tropopauseThickness   = MAX(2,nlev/30)
    invInd = 0

    previousNotInversion  = .TRUE.

    ! Loop from the surface up to the tropopause
    tropopause_count = 0

    DO d1 = 1,ngrids
       DO inl = nlev,2,-1
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
