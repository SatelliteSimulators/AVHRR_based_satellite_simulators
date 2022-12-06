MODULE MODEL_INPUT

  !  A tool to handle model data. So far this is set up for EC Earth and RACMO.
  !
  !  Salomon.Eliasson@smhi.se
  ! adapted for RACMO by Erik van Meijgaard (KNMI)
  !

  USE COSP_KINDS,          ONLY:&
       WP
  USE HANDY,               ONLY:&
       CHECK_FILE,              &
       REPLACE_TEXT,            &
       TALLY
  USE MY_MATHS,            ONLY:&
       DAYSINMONTH
  USE MY_NETCDFTOOLS,      ONLY:&
       CHECK
  USE NAMELIST_INPUT,      ONLY:&
       NAME_LIST                      
  USE NETCDF,              ONLY:&
       NF90_CLOSE,              &
       NF90_GET_VAR,            &
       NF90_GET_ATT,            &
       NF90_INQUIRE,            &
       NF90_INQUIRE_DIMENSION,  &
       NF90_INQ_VARID,          &
       NF90_OPEN,               &
       NF90_NOWRITE

  PUBLIC :: ALLOCATE_MODEL_MATRIX,&
       DEALLOCATE_MODEL_AUX,      &
       DEALLOCATE_MODEL_MATRIX,   &
       GET_MODEL_AUX,             &
       INITIALISE_MODEL_MATRIX,   &
       MODEL_TYPE,                &
       READ_MODEL

  PRIVATE :: GET_REF_TIME, GET_LAND_SEA_MASK

  ! These variables definitions are from
  ! http://old.ecmwf.int/publications/manuals/d/gribapi/param/filter=grib1/order=paramId/order_type=asc/p=1/table=128/
  !
  ! Auxilliary
  ! 
  ! 'areacella' = area of gridcell 
  ! 'hyai' = hybrid A coefficient at layer interfaces [Pa].
  ! 'hyam' = hybrid A coefficient at layer midpoints [Pa].
  ! 'hybi' = hybrid B coefficient at layer interfaces [1].
  ! 'hybm' = hybrid B coefficient at layer midpoints [1].
  ! 'rlat' = latitude [degrees_north] in Rotated Frame (or one dimensional lat)
  ! 'lat'  = latitude [degrees_north] 2-Dimensional Array
  ! 'rlon' = Longitude [Degrees East] in Rotated Frame (or one dimensional lon)
  ! 'lon'  = Longitude [Degrees East] 2-Dimensional Array
  ! 'lsm'  = land/sea mask (fraction land cover)
  ! 'time' = time step (e.g. hours since the start of the file)
  ! 'modelRes'= resolution of the model in degress
  ! 'fvr,fvi'= fill values
  ! 'nlon, nlat, nlev, ntlen' = dimensions of the data
  ! 'netcdf_file' = filename of model input
  ! 'time_units'= units for time
  ! 'calendar' = 'calendar used' 
  !
  ! MODEL FIELDS
  !
  ! 'CC'   = Horizontal fraction of the grid box covered by cloud [0-1]. code=248
  ! 'CI'   = Sea Ice Cover [0-1]. code=31
  ! 'CIWC' = Specific cloud ice water content. Grid-box
  !          mean specific cloud ice water content (mass of condensate
  !          / mass of moist air).  [kg kg**-1]. code=247
  ! 'CLWC' = Specific cloud liquid water content. Grid-box mean
  !          specific cloud liquid water content (mass of condensate /
  !          mass of moist air) [kg kg**-1]. code=246
  ! 'PSURF'= Surface pressure [(Pa)]. code=1
  ! 'Q'    = Specific humidity. Grid box mean (mass of water vapour /
  !          mass of moist air). [kg kg**-1]. code=133
  ! 'SKT'  = Skin temperature. Temperature of the surface skin
  !          (radiative surface temperature). Before 01/10/2008, the skin
  !          temperature was equal to the bulk SST over the
  !          ocean. [K]. code=235
  ! 'T'    = Temperature [K]. code=130
  ! 'T2M'  = 2 metre temperature [K]. code=167
  ! 'TCC'  = Total cloud cover. Total cloud cover derived from model
  !          levels using the model's overlap assumption [0-1]. code=164
  ! 'TCWV' = Total column water vapour (Vertically integrated water
  !          vapour) [kg m**-2]. code=137
  !       Missing but created for now
  ! 'CV'   = convective cloud cover [0-1]; doesn't exist in RACMO

  ! 'Q2M'  = specific humidity at 2m [kg kg**-1]
  ! 'surf_type' = surface type (0=land, 1=sea, 2=sea-ice)

  TYPE REF_TIME_TYPE
     ! time of reference
     CHARACTER(len =  100) :: time_units
     CHARACTER(len =  100) :: calendar
     INTEGER :: year,month,day,hour,idtg
     REAL(wp) :: deltaTime
     
  END TYPE REF_TIME_TYPE
  TYPE MODEL_AUX

     REAL(wp), ALLOCATABLE, DIMENSION(:)   :: areacella
     REAL(wp), ALLOCATABLE, DIMENSION(:)   :: hyai, hyam, hybi, hybm
     REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: lat, lon
     REAL(wp), ALLOCATABLE, DIMENSION(:)   :: lon_v, lat_v
     REAL(wp), ALLOCATABLE, DIMENSION(:)   :: lsm
     REAL(wp), ALLOCATABLE, DIMENSION(:)   :: time
     TYPE(ref_time_type) :: ref 
     REAL(wp) :: modelRes
     REAL(wp) :: fvr = -999._wp ! 1.e+20
     INTEGER  :: fvi = -999
     INTEGER  :: nlon, nlat, nlev, ntlen, ngrids
     CHARACTER(len = 1000) :: netcdf_file

     ! RACMO specific
     REAL(wp), ALLOCATABLE, DIMENSION(:)   :: rlat, rlon
     CHARACTER(len = 128) :: rotgrid,rotproj
     REAL                 :: float_rotpolat,float_rotpolon
  END TYPE MODEL_AUX

  TYPE MODEL_TYPE
     TYPE(model_aux) :: aux
     REAL(wp), ALLOCATABLE, DIMENSION(:,:)  :: CC   
     REAL(wp), ALLOCATABLE, DIMENSION(:)    :: CI   
     REAL(wp), ALLOCATABLE, DIMENSION(:,:)  :: CIWC 
     REAL(wp), ALLOCATABLE, DIMENSION(:,:)  :: CLWC
     REAL(wp), ALLOCATABLE, DIMENSION(:)    :: PSURF
     REAL(wp), ALLOCATABLE, DIMENSION(:,:)  :: Q    
     REAL(wp), ALLOCATABLE, DIMENSION(:)    :: SKT  
     REAL(wp), ALLOCATABLE, DIMENSION(:,:)  :: T    
     REAL(wp), ALLOCATABLE, DIMENSION(:)    :: T2M  
     REAL(wp), ALLOCATABLE, DIMENSION(:)    :: TCC  
     REAL(wp), ALLOCATABLE, DIMENSION(:)    :: TCWV 
     ! this is not actually in the model
     REAL(wp), ALLOCATABLE, DIMENSION(:,:)  :: CV
  END TYPE MODEL_TYPE

CONTAINS

  ! -----------------
  ! aux
  ! -----------------

  SUBROUTINE GET_MODEL_AUX(aux,options,tim)
    ! This function is only for variables that are constant for every
    ! model file
    !
    !
    IMPLICIT NONE 

    TYPE(model_aux), INTENT(inout) :: aux
    TYPE(name_list), INTENT(inout) :: options
    LOGICAL, OPTIONAL, INTENT(in)  :: tim ! optional time flag

    !local
    LOGICAL              :: time,lorlon,lorlat
    INTEGER              :: ncid,ncid2,varid,dimid,len,nDimensions
    INTEGER              :: nlon,nlat,ln,lt
    REAL(wp)             :: drlonrad,drlatrad
    REAL(wp),ALLOCATABLE ::lsm(:,:),tmp_lon(:),tmp_lat(:)
    CHARACTER (len=64)   :: str
    CHARACTER (len=64)   :: ctstr
    CHARACTER (len=64)   :: ccstr
    CHARACTER(len = 128) :: crotgrid,crotpolat,crotpolon,crotproj
    CHARACTER(len = 1000):: file
    REAL(wp), PARAMETER  :: pi = 3.1415926536
    REAL(wp), PARAMETER  :: d2r = pi/180.
    REAL(wp), PARAMETER  :: earth_radius = 6.3781e+6
    INTEGER, PARAMETER   :: step = 10
    INTEGER              :: dbg

    dbg = options%dbg

    time = .FALSE.
    IF (PRESENT(tim)) time = tim

    ncid = 0
    varid = 0
    dimid = 0
    CALL CHECK( nf90_open(aux%netcdf_file, nf90_nowrite, ncid) )
    CALL CHECK( nf90_inquire(ncid, nDimensions) ) 
    DO dimid = 1, nDimensions

       ! Get 'len' and 'str' for every variable
       CALL CHECK (nf90_inquire_dimension(ncid, dimid, str, len) )

       ! skip if str = time if you don't want it or vise versa
       IF (.NOT. time .AND. str .EQ. 'time') CYCLE
       IF (time .AND. str .NE. 'time') CYCLE
       IF (str .EQ. 'time') THEN 
          ALLOCATE( aux%time(len) )
          aux%time=-999._wp
          CALL CHECK( nf90_inq_varid(ncid, str, varid) )
          CALL CHECK( nf90_get_var(ncid, varid, aux%time) ,dbg,'reading '//str);
          aux%ntlen = len
          CALL CHECK( nf90_get_att(ncid,varid,"units",ctstr) )
          CALL CHECK( nf90_get_att(ncid,varid,"calendar",ccstr) )
        aux%ref%time_units = ctstr
          aux%ref%calendar   = ccstr
          aux%ref = GET_REF_TIME(aux,dbg)
       ELSE IF (str .EQ. 'nhym') THEN
          ALLOCATE(aux%hyam(len),aux%hybm(len) )
          CALL CHECK( nf90_inq_varid(ncid, 'hyam', varid) )
          CALL CHECK( nf90_get_var(ncid, varid, aux%hyam) ,dbg,'reading '//str);
          CALL CHECK( nf90_inq_varid(ncid, 'hybm', varid) )
          CALL CHECK( nf90_get_var(ncid, varid, aux%hybm) ,dbg,'reading '//str);
          aux%nlev = len
       ELSE IF (str .EQ. 'nhyi') THEN
          ALLOCATE(aux%hyai(len), aux%hybi(len) ) 
          CALL CHECK( nf90_inq_varid(ncid, 'hyai', varid) )
          CALL CHECK( nf90_get_var(ncid, varid, aux%hyai) ,dbg,'reading '//str);
          CALL CHECK( nf90_inq_varid(ncid, 'hybi', varid) )
          CALL CHECK( nf90_get_var(ncid, varid, aux%hybi) ,dbg,'reading '//str);
       END IF

       ! EC EARTH-specifc
       IF (options%model .EQ. 'ec_earth') THEN
          IF (str .EQ. 'lon') THEN
             ALLOCATE(tmp_lon(len))
             CALL CHECK( nf90_inq_varid(ncid, str, varid) )
             CALL CHECK( nf90_get_var(ncid, varid, tmp_lon) ,dbg, 'reading '//str);
             aux%nlon = len
             lorlon   =.TRUE.
          ELSE IF (str .EQ. 'lat') THEN
             ALLOCATE(tmp_lat(len))
             CALL CHECK( nf90_inq_varid(ncid, str, varid) )
             CALL CHECK( nf90_get_var(ncid, varid, tmp_lat) ,dbg,'reading '//str);
             aux%nlat = len
             lorlat   =.TRUE.
          END IF
       END IF

       ! RACMO-specific
       IF (options%model .EQ. 'racmo') THEN
          IF (str .EQ. 'rlon') THEN
             ALLOCATE( aux%rlon(len) ) 
             CALL CHECK( nf90_inq_varid(ncid, str, varid) )
             CALL CHECK( nf90_get_var(ncid, varid, aux%rlon) ,dbg, 'reading '//str);
             aux%nlon = len
             lorlon = .TRUE.
          ELSE IF (str .EQ. 'rlat') THEN
             ALLOCATE( aux%rlat(len) ) 
             CALL CHECK( nf90_inq_varid(ncid, str, varid) )
             CALL CHECK( nf90_get_var(ncid, varid, aux%rlat) ,dbg,'reading '//str);
             aux%nlat = len
             lorlat = .TRUE.
          END IF
       END IF
    END DO

    IF (.NOT.time) THEN
       IF (lorlon .AND. lorlat) THEN

          ! --------------------
          ! LON LAT
          !
          nlon = aux%nlon
          nlat = aux%nlat 
          IF (options%model .EQ. 'ec_earth') THEN
             ! ... lon lat
             ALLOCATE(aux%lon(nlon,nlat),&
                  aux%lat(nlon,nlat))
             aux%lon = SPREAD(tmp_lon,2,nlat)
             aux%lat = SPREAD(tmp_lat,1,nlon)
          END IF
          IF (options%model .EQ. 'racmo') THEN
             ! horizontal spatial dimensions and corresponding coordinate variables have been imported
             ! two-dimensional arrays lon and lat can be read here 
             ! ... 2-dimensional longitude
             str='lon'
             ALLOCATE( aux%lon(nlon,nlat) ) 
             CALL CHECK( nf90_inq_varid(ncid, str, varid) )
             CALL CHECK( nf90_get_var(ncid, varid, aux%lon) ,dbg, 'reading '//str);
             ! ... 2-dimensional latitude
             str='lat'
             ALLOCATE( aux%lat(nlon,nlat) ) 
             CALL CHECK( nf90_inq_varid(ncid, str, varid) )
             CALL CHECK( nf90_get_var(ncid, varid, aux%lat) ,dbg, 'reading '//str);
          ENDIF

          aux%ngrids=nlon*nlat

          ! make a reshaped vector of lon lat
          ALLOCATE(aux%lon_v(aux%ngrids),&
               aux%lat_v(aux%ngrids))
          aux%lon_v = RESHAPE(aux%lon,(/ aux%ngrids /))
          aux%lat_v = RESHAPE(aux%lat,(/ aux%ngrids /)) 
          ! -------------------- end LON LAT

          ! --------------------
          ! Land Sea mask
          !
          ALLOCATE(aux%lsm(nlon*nlat))
          IF (options%model .EQ. 'ec_earth') THEN
             ALLOCATE( lsm(nlon,nlat) ) 
             WRITE(file, '("data/land_sea_mask.nc")')
             ncid2=0
             dimid=1
             IF (CHECK_FILE(TRIM(file))) THEN
                PRINT '(2A)'," --- getting land mask from:",TRIM(file) 
                CALL CHECK( nf90_open     (TRIM(file),nf90_nowrite,ncid2 ) )
                ! check if this land sea mask fits 
                str='lon'
                CALL CHECK(nf90_inquire_dimension( ncid2,dimid,str,len ) )
                IF  (len.NE.nlon) THEN 
                   PRINT '(100A)', "land sea mask dimensions do not match model.&
                        Will make a land sea mask that matches the model"
                   aux%lsm = GET_LAND_SEA_MASK(aux,options)
                ELSE
                   str='LSM'
                   CALL CHECK( nf90_inq_varid(ncid2,str,varid) )
                   CALL CHECK( nf90_get_var  (ncid2,varid,lsm) ,dbg, 'reading '//str);
                   aux%lsm=RESHAPE(lsm, (/ nlon*nlat /))
                END IF
                CALL CHECK( nf90_close    (ncid2))             
             ELSE
                aux%lsm = GET_LAND_SEA_MASK(aux,options)
             END IF
          ELSEIF (options%model .EQ. 'racmo') THEN
             ! ... two-dimensional array lsm can be read here and stored in lsm2
             str='lsm'
             ALLOCATE( lsm(nlon,nlat) ) 
             CALL CHECK( nf90_inq_varid(ncid, str, varid) )
             CALL CHECK( nf90_get_var(ncid, varid, lsm) ,dbg, 'reading '//str);
             aux%lsm=RESHAPE(lsm, (/ nlon*nlat /))

             IF (dbg>1) THEN
                WRITE (0,'(a)') 'EvM: land sea mask'
                WRITE (0,'(3x,100I7)') (ln,ln=1,nlon,nlat/step)
                DO lt=nlat,1,-nlat/step
                   WRITE (0,'(I3,100F7.3)') lt,(NINT(lsm(ln,lt)*1000.)/1000.,ln=1,nlon,nlat/step)
                ENDDO
             ENDIF
          ELSE
             ! ... land sea mask
             aux%lsm = GET_LAND_SEA_MASK(aux, options)
          END IF
          ! -------------------- end LAND SEA MASK

          ! --------------------
          ! extras
          IF (options%model .EQ. 'racmo') THEN
             ! ... rotated pole attributes
             str='rotated_pole'
             CALL CHECK( nf90_inq_varid(ncid, str, varid) )
             CALL CHECK( nf90_get_att(ncid,varid,"grid_mapping_name",crotgrid) )
             CALL CHECK( nf90_get_att(ncid,varid,"grid_north_pole_latitude",crotpolat) )
             CALL CHECK( nf90_get_att(ncid,varid,"grid_north_pole_longitude",crotpolon) )
             CALL CHECK( nf90_get_att(ncid,varid,"proj",crotproj) )
             aux %  rotgrid  = crotgrid
             !        aux %  rotpolat = crotpolat
             !        aux %  rotpolon = crotpolon
             aux %  rotproj  = crotproj
             READ (crotpolat,'(f6.2)') aux %  float_rotpolat
             READ (crotpolon,'(f6.2)') aux %  float_rotpolon
             WRITE (0,'(a)') ' imported rotated pole attributes'
          ENDIF
          ! ------ end EXTRAS

          ! --------------------
          ! Grid area model resolution etc
          !
          IF (options%model .EQ. 'ec_earth') THEN
             ! ... grid area
             ALLOCATE( aux%areacella(nlat) ) 
             drlonrad=ABS(aux%lon(2,1)-aux%lon(1,1))*d2r
             drlatrad=ABS(aux%lat(1,2)-aux%lat(1,1))*d2r
             aux%areacella=2.*drlonrad*COS(aux%lat_v*d2r)*SIN(drlatrad/2.)*earth_radius**2

             ! ... model resolution
             aux%modelRes = (ABS(tmp_lon(2)-tmp_lon(1))+ABS(tmp_lat(2)-tmp_lat(1)))/2._wp
          END IF
          IF (options%model .EQ. 'racmo') THEN
             ! ... grid area
             ALLOCATE( aux%areacella(nlat) ) 
             drlonrad=ABS(aux%rlon(2)-aux%rlon(1))*d2r
             drlatrad=ABS(aux%rlat(2)-aux%rlat(1))*d2r
             aux%areacella=2.*drlonrad*COS(aux%rlat*d2r)*SIN(drlatrad/2.)*earth_radius**2

             ! ... model resolution
             aux%modelRes = (ABS(aux%rlon(2)-aux%rlon(1))+ABS(aux%rlat(2)-aux%rlat(1)))/2._wp
          END IF

          ! ---------
          ! Allocating ncols
          !
          ! Recommended 100xgrid degree. There is hardly any difference if
          ! you use 20 subgrids per model grid column or the number that
          ! would give a 5km footprint
          !
          
          options%ncols = CEILING(aux%modelRes*100)

          PRINT '(4(a,1x,i3))', " --- Using a static number of columns for a",&
               aux%nlon," x",aux%nlat," x",aux%nlev," grid: ncol =",options%ncols

          !-------------------- end AREA RESOLUTION ETC

       ELSE
          WRITE (0,'(a)') ' racmo_model_input: dimensions rlon and rlat not encountered '
          STOP
       ENDIF

    END IF

    IF (ALLOCATED(tmp_lon)) DEALLOCATE(tmp_lon,tmp_lat)

    CALL CHECK( nf90_close (ncid))

  END SUBROUTINE get_model_aux

  SUBROUTINE deallocate_model_aux(aux)

    IMPLICIT NONE
    TYPE(model_aux), INTENT(inout) :: aux

    DEALLOCATE (&
         aux%areacella,&
         aux%hyam     ,&   
         aux%hybm     ,&
         aux%hyai     ,&
         aux%hybi     ,&
         aux%lat      ,&
         aux%lon      ,&
         aux%lat_v    ,&
         aux%lon_v    ,&
         aux%lsm      ,&
         aux%time)

    IF (ALLOCATED(aux%rlon)) DEALLOCATE(aux%rlat,aux%rlon)

  END SUBROUTINE deallocate_model_aux

  ! -------
  ! model fields
  ! -------

  SUBROUTINE allocate_model_matrix(M,ngrids,nlev)
    ! include ngrids, nlev, etc in the input as M%aux is not always allocated

    IMPLICIT NONE

    TYPE(model_type), INTENT(inout) :: M
    INTEGER, INTENT(in)             :: ngrids,nlev

    ALLOCATE ( M%CC     (ngrids, nlev) ,&  
               M%CI     (ngrids      ) ,&              
               M%CIWC   (ngrids, nlev) ,&  
               M%CLWC   (ngrids, nlev) ,&  
               M%CV     (ngrids, nlev) ,&
               M%PSURF  (ngrids      ) ,&           
               M%Q      (ngrids, nlev) ,&  
               M%SKT    (ngrids      ) ,&              
               M%T      (ngrids, nlev) ,&      
               M%T2M    (ngrids      ) ,&              
               M%TCC    (ngrids      ) ,&              
               M%TCWV   (ngrids      ) )              

  END SUBROUTINE allocate_model_matrix

  SUBROUTINE initialise_model_matrix(M,ngrids,nlev,mv)
    ! include nlon, etc in the input as M%aux is not always allocated

    IMPLICIT NONE

    TYPE(model_type), INTENT(inout) :: M
    INTEGER, INTENT(in)             :: ngrids,nlev
    REAL(wp),INTENT(in),OPTIONAL    :: mv

    REAL(wp)                        :: init

    IF (.NOT.PRESENT(mv)) THEN
       init=-999._wp
    ELSE
       init=mv
    END IF
    M%CC   (1:ngrids, 1:nlev)  = init
    M%CI   (1:ngrids)          = init 
    M%CIWC (1:ngrids, 1:nlev)  = init
    M%CLWC (1:ngrids, 1:nlev)  = init
    M%CV   (1:ngrids, 1:nlev)  = init
    M%PSURF(1:ngrids)          = init
    M%Q    (1:ngrids, 1:nlev)  = init  
    M%SKT  (1:ngrids)          = init
    M%T    (1:ngrids, 1:nlev)  = init
    M%T2M  (1:ngrids)          = init
    M%TCC  (1:ngrids)          = init   
    M%TCWV (1:ngrids)          = init

  END SUBROUTINE initialise_model_matrix

  SUBROUTINE deallocate_model_matrix(M)

    IMPLICIT NONE
    TYPE(model_type), INTENT(inout) :: M

    DEALLOCATE ( M%CC,&
         M%CI   ,& 
         M%CIWC ,&
         M%CLWC ,&
         M%CV   ,&
         M%PSURF,&
         M%Q    ,&  
         M%SKT  ,&   
         M%T    ,&
         M%T2M  ,&
         M%TCC  ,&   
         M%TCWV)

  END SUBROUTINE deallocate_model_matrix

  SUBROUTINE READ_MODEL(M,itime,options,verbose)
    ! Function that reads in the needed model fields and returns them
    ! in a vectorised format which scops.f and the simulator wants it in

    IMPLICIT NONE

    TYPE(model_type), INTENT(inout) :: M
    INTEGER, INTENT(in)             :: itime
    INTEGER, OPTIONAL, INTENT(in)   :: verbose 
    TYPE(name_list), INTENT(in)     :: options

    INTEGER               :: ncid, varid
    INTEGER, DIMENSION(3) :: start3, count3
    INTEGER, DIMENSION(4) :: start4, count4
    INTEGER :: dbg
    INTEGER :: ngrids,nlon,nlat,nlev
    CHARACTER(64):: &
         strCC     ,&
         strCI     ,&
         strCIWC   ,&
         strCLWC   ,&
         strQ      ,&
         strSKT    ,&
         strT      ,&
         strT2M    ,&
         strTCC    ,&
         strTCWV   ,&
         var

    ! temporary variables before regridding
    REAL(wp), DIMENSION(M%aux%nlon,M%aux%nlat) :: CI,SKT,PSURF,T2M,TCC,TCWV 
    REAL(wp), DIMENSION(M%aux%nlon,M%aux%nlat,M%aux%nlev) :: CC,CIWC,CLWC,Q,T
    REAL(wp), DIMENSION(M%aux%nlon,M%aux%nlat,1) :: tmpPsurf

    nlon  =M%aux%nlon
    nlat  =M%aux%nlat
    nlev  =M%aux%nlev
    ngrids=M%aux%ngrids
    IF (.NOT.PRESENT(verbose)) THEN
       dbg = options%dbg
    ELSE
       dbg = verbose
    END IF
    ncid = 0
    varid = 0

    ! open file
    CALL CHECK( NF90_OPEN(M%aux%netcdf_file, NF90_NOWRITE, ncid))

    ! I know the names and dimensions of what I want, so I'll just read
    !them directly

    start3 = [1,1,itime];   count3 = [nlon,nlat,1]
    start4 = [1,1,1,itime]; count4 = [nlon,nlat,nlev,1]
    IF (dbg>1) THEN
       WRITE (0,*) 'EvM-read_model: start3 = ',start3
       WRITE (0,*) 'EvM-read_model: count3 = ',count3
    ENDIF

    ! Fields where only the names are different between models
    SELECT CASE(options%model)
    CASE('ec_earth')
       strCC='CC'
       strCI='CI'
       strCIWC='CIWC'
       strCLWC='CLWC'
       strQ='Q'
       strSKT='SKT'
       strT='T'
       strT2M='T2M'
       strTCC='TCC'
       strTCWV='TCWV'
    CASE('racmo')
       strCC='cfrac'
       strCI='ci'
       strCIWC='qi'
       strCLWC='ql'
       strQ='q'
       strSKT='tsurf'
       strT='t'
       strT2M='t2m'
       strTCC='tcc'
       strTCWV='qvi'
    CASE default
       WRITE(0,*) "variable string names not set up for ",options%model
    END SELECT

    ! cloud cover
!    PRINT '(A,A)',"Reading ",strCC 
    CALL CHECK( NF90_INQ_VARID(ncid, TRIM(strCC) , varid),dbg)
    CALL CHECK( NF90_GET_VAR(ncid, varid, CC, start4, count4) )

    ! sea ice
!    PRINT '(A,A)',"Reading ",strCI 
    CALL CHECK( NF90_INQ_VARID(ncid, trim(strCI) , varid),dbg)
    CALL CHECK( NF90_GET_VAR(ncid, varid, CI, start3, count3) )

    ! cloud ice water content
!    PRINT '(A,A)',"Reading ",strCIWC
    CALL CHECK( NF90_INQ_VARID(ncid, trim(strCIWC) , varid),dbg)
    CALL CHECK( NF90_GET_VAR(ncid, varid, CIWC, start4, count4) )

    ! cloud liquid water content
!    PRINT '(A,A)',"Reading ",strCLWC 
    CALL CHECK( NF90_INQ_VARID(ncid, trim(strCLWC) , varid),dbg)
    CALL CHECK( NF90_GET_VAR(ncid, varid, CLWC, start4, count4) )

    ! water vapour profile
!    PRINT '(A,A)',"Reading ",strQ 
    CALL CHECK( NF90_INQ_VARID(ncid, trim(strQ) , varid),dbg)
    CALL CHECK( NF90_GET_VAR(ncid, varid, Q, start4, count4) )

    ! surface skin temperature
!    PRINT '(A,A)',"Reading ",strSKT 
    CALL CHECK( NF90_INQ_VARID(ncid, trim(strSKT) , varid),dbg)
    CALL CHECK( NF90_GET_VAR(ncid, varid, SKT, start3, count3) )

    ! Temperature profile
!    PRINT '(A,A)',"Reading ",strT 
    CALL CHECK( NF90_INQ_VARID(ncid, trim(strT) , varid),dbg)
    CALL CHECK( NF90_GET_VAR(ncid, varid, T, start4, count4) ) 

    ! 2m temperature
!    PRINT '(A,A)',"Reading ",strT2M 
    CALL CHECK( NF90_INQ_VARID(ncid, trim(strT2M) , varid),dbg)
    CALL CHECK( NF90_GET_VAR(ncid, varid, T2M, start3, count3) )

    ! Total cloud cover
!    PRINT '(A,A)',"Reading ",strTCC 
    CALL CHECK( NF90_INQ_VARID(ncid, trim(strTCC) , varid),dbg)
    CALL CHECK( NF90_GET_VAR(ncid, varid, TCC, start3, count3) )

   ! Total column water vapour
!    PRINT '(A,A)',"Reading ",strTCWV 
    CALL CHECK( NF90_INQ_VARID(ncid, TRIM(strTCWV) , varid),dbg)
    CALL CHECK( NF90_GET_VAR(ncid, varid, TCWV, start3, count3) )

    ! Surface pressure
    IF (options%model.EQ.'ec_earth') THEN 
       var='LNSP'
!       PRINT '(A,A)',"Reading ",var 
       CALL CHECK( NF90_INQ_VARID(ncid, TRIM(var) , varid),dbg)
       CALL CHECK( NF90_GET_VAR(ncid, varid, tmpPsurf, start4, [nlon,nlat,1,1]) )
       PSURF(1:nlon,1:nlat)=EXP(tmpPsurf(1:nlon,1:nlat,1))
    ELSEIF (options%model.EQ.'racmo') THEN
!       PRINT '(A,A)',"Reading ",var 
       var='ps'
       CALL CHECK( NF90_INQ_VARID(ncid, TRIM(var) , varid),dbg)
       CALL CHECK( NF90_GET_VAR(ncid, varid, PSURF, start3, count3) ) 
    END IF

    CALL CHECK( NF90_CLOSE(ncid))

    ! Now do the regridding
    M%CC    = RESHAPE(CC   , (/ ngrids,nlev /))
    M%CI    = RESHAPE(CI   , (/ ngrids      /))
    M%CIWC  = RESHAPE(CIWC , (/ ngrids,nlev /))
    M%CLWC  = RESHAPE(CLWC , (/ ngrids,nlev /))
    M%PSURF = RESHAPE(PSURF, (/ ngrids      /))
    M%Q     = RESHAPE(Q    , (/ ngrids,nlev /))
    M%SKT   = RESHAPE(SKT  , (/ ngrids      /))
    M%T     = RESHAPE(T    , (/ ngrids,nlev /))
    M%T2M   = RESHAPE(T2M  , (/ ngrids      /))
    M%TCC   = RESHAPE(TCC  , (/ ngrids      /))
    M%TCWV  = RESHAPE(TCWV , (/ ngrids      /))

  END SUBROUTINE READ_MODEL

  FUNCTION get_ref_time(aux,dbg) RESULT(ref)

    IMPLICIT NONE

    TYPE(model_aux),INTENT(in)   :: aux
    INTEGER, INTENT(in)          :: dbg
    TYPE(ref_time_type)          :: ref
    CHARACTER(len=124)           :: time_units
    CHARACTER(len=124)           :: calendar
    CHARACTER(len=124)           :: time_ref_str
    INTEGER                      :: itime
    time_units= aux%ref%time_units
    calendar  = aux%ref%calendar
    ref=aux%ref

    IF (TRIM(calendar).NE.'standard'.AND.TRIM(calendar).NE.'gregorian'.AND.TRIM(calendar).NE.'proleptic_gregorian') THEN
       WRITE (0,'(3a)') ' calendar attribute ',TRIM(calendar),' not anticipated ... stop'
    ENDIF
    IF (time_units(1:11).EQ.'hours since') THEN
       time_ref_str=time_units(13:31)
    ELSE IF (time_units(1:11).EQ.'days since') THEN
       time_ref_str=time_units(12:30)
    ELSE
       WRITE (0,'(3a)') ' time units attribute ',TRIM(time_units),' not anticipated ... stop'
    ENDIF
    READ(time_ref_str( 1: 4),'(i4)') ref%year
    READ(time_ref_str( 6: 7),'(i2)') ref%month
    READ(time_ref_str( 9:10),'(i2)') ref%day
    READ(time_ref_str(12:13),'(i2)') ref%hour
    ref%idtg=ref%year*1000000+ref%month*10000+ref%day*100+ref%hour
    ref%deltaTime=aux%time(2)-aux%time(1)
    IF (dbg>1) THEN
       WRITE (0,'(a,i6)') 'EvM: ntlen = ',aux%ntlen
       WRITE (0,'(a,a)')  'EvM: time_units ',TRIM(aux%ref%time_units)
       WRITE (0,'(a,a)')  'EvM: calendar   ',TRIM(aux%ref%calendar)
       WRITE (0,'(a,i6)') 'EvM: times '
       WRITE (0,'(8f15.5)') (aux%time(itime),itime=1,aux%ntlen)
       WRITE (0,'(a,i12)')  'EvM: reference dtg ',ref%idtg
    END IF

  END FUNCTION get_ref_time

  FUNCTION get_land_sea_mask(aux, options) RESULT(lsm)
    !
    ! Internally the model knows if it is over land or not, but running
    ! the simulator offline means that I do not have this information. I
    ! import a high resolution land mask and calculate the land fraction
    ! from it
    !
    IMPLICIT NONE

    TYPE(model_aux), INTENT(in) :: aux
    TYPE(name_list), INTENT(in) :: options

    ! output
    REAL(wp), DIMENSION(aux%nlon*aux%nlat) :: lsm
    ! netcdf
    CHARACTER(len = 1000) :: file
    CHARACTER(len = 10) :: str
    INTEGER:: ncid, varid, dimid, nDimensions, len
    INTEGER ::nlon, nlat
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: lon, lat
    INTEGER(kind = 4), ALLOCATABLE, DIMENSION(:,:) :: HR_lsm
    ! putting the mask together
    INTEGER :: lt, ln, iln, ilt

    REAL(wp), DIMENSION(aux%nlon,aux%nlat) :: land, tot
    ! ====================================
    ! READ the land sea mask 
    !

    WRITE(file, '("data/land_sea_mask/land_sea_mask_1min.nc")')

    IF (.NOT. CHECK_FILE (file)) THEN
       STOP "where's the file?"
    END IF

    ncid=0
    dimid=0
    CALL CHECK( nf90_open( file, nf90_nowrite, ncid ) )
    CALL CHECK( nf90_inquire( ncid, nDimensions ) ) 

    ! Get the dimensions that span the look up tables
    DO dimid = 1, nDimensions
       CALL CHECK (nf90_inquire_dimension( ncid, dimid, str, len ) )
       CALL CHECK( nf90_inq_varid(ncid, str, varid) )

       IF (str .EQ. 'lon') THEN
          ALLOCATE( lon(len) )
          CALL CHECK( nf90_get_var(ncid, varid, lon) ,options%dbg, 'reading '//str);
          nlon = len
       ELSE IF (str .EQ. 'lat') THEN
          ALLOCATE( lat(len) )
          CALL CHECK( nf90_get_var(ncid, varid, lat) ,options%dbg, 'reading '//str);
          nlat = len
       END IF
    END DO

    ALLOCATE ( HR_lsm(nlon, nlat) )
    CALL CHECK( nf90_inq_varid(ncid, 'land_sea_mask', varid),&
         options%dbg, "reading land_sea_mask")
    CALL CHECK( nf90_get_var(ncid, varid, HR_lsm) ) 

    CALL CHECK( nf90_close(ncid))
    !
    ! ========================

    land(1:aux%nlon,1:aux%nlat) = 0._wp
    tot(1:aux%nlon,1:aux%nlat)  = 0._wp
    
    ilt = 1
    DO lt = 1, nlat
       iln = 1
       IF (lat(lt) < aux%lat(1,ilt) .AND. ilt<aux%nlat) ilt=ilt+1 ! descending latitudes
       DO ln = 1, nlon
          IF (lon(ln) > aux%lon(iln,1) .AND. iln<aux%nlon) iln=iln+1

          tot(iln,ilt) = tot(iln,ilt) + 1
          ! 1 = land, 0 = sea
          land(iln,ilt) = land(iln,ilt) + REAL(HR_lsm(ln,lt),wp)

       END DO
    END DO

    lsm = RESHAPE(land/tot, (/ aux%nlon*aux%nlat /))
    DEALLOCATE (lon, lat, HR_lsm)

  END FUNCTION get_land_sea_mask

END MODULE model_input

