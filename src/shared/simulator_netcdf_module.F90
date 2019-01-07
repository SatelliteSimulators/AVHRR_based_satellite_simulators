MODULE simulator_netcdf_module
  !
  ! Salomon.Eliasson@smhi.se
  !
  ! JFM, 15/07/2014: renamed + added some output fields
  ! Cloud_cci only

  USE CLOUD_CCI_M,          ONLY: CLOUD_CCI_TYPE
  USE COSP_KINDS,           ONLY: WP
  USE MODEL_INPUT,          ONLY: MODEL_TYPE,MODEL_AUX
  USE MY_NETCDFTOOLS,       ONLY: ADD_VARIABLE_ATTRIBUTES, CHECK
  USE NETCDF, ONLY :  &
       NF90_CLOBBER,  &
       NF90_CLOSE,    &
       NF90_CREATE,   &
       NF90_CHAR,     &
       NF90_UNLIMITED,&
       NF90_DEF_DIM,  &
       NF90_DEF_VAR,  &
       NF90_DEF_VAR_DEFLATE,&
       NF90_ENDDEF,   &
       NF90_GLOBAL,   &
       NF90_INT,      &
       NF90_NOERR,    &
       NF90_OPEN,     &
       NF90_PUT_ATT,  &
       NF90_PUT_VAR,  &
       NF90_REAL,     &
       NF90_REDEF,    &
       NF90_STRERROR, &
       NF90_WRITE
  USE NAMELIST_INPUT,            ONLY: NAME_LIST,VARIABLESCONTAINER
  USE SATELLITE_SPECS,           ONLY: SATELLITE
  USE SIMULATOR_VARIABLES,       ONLY: SATELLITE_SIMULATOR
  USE SIMULATOR_INPUT_VARIABLES, ONLY: SUBSET

  IMPLICIT NONE

  PUBLIC  :: MAKE_NETCDF
  PRIVATE :: ADDVARIABLE

CONTAINS

  SUBROUTINE MAKE_NETCDF(model,sim,O,iday,&
       S,sat,cloud_cci)

    !
    ! tailor-made module to save model and simulator output to a netcdffile
    !
    ! This function is constantly being updated depending on what I
    ! feel like saving (even though most inputs are optional).

    IMPLICIT NONE


    ! --- IN ---
    TYPE(name_list), INTENT(in)                     :: O
    INTEGER, INTENT(in)                             :: iday
    TYPE(model_type),                    INTENT(in) :: model
    TYPE(satellite_simulator),           INTENT(in) :: sim
    TYPE(satellite),           OPTIONAL, INTENT(in) :: sat
    TYPE(subset),              OPTIONAL, INTENT(in) :: S
    TYPE(cloud_cci_type),      OPTIONAL, INTENT(in) :: cloud_cci

    ! --- LOCAL ---
    TYPE(variablesContainer) :: V
    TYPE(model_aux)          :: A
    INTEGER                  :: ncid
    INTEGER                  :: lx,ly,lz,lt,n_tbins,n_pbins,nchan
    CHARACTER(len=1000)      :: string
    INTEGER, DIMENSION(8)    :: date_time
    INTEGER, DIMENSION(2)    :: utc
    CHARACTER(len=1000)      :: creation
    INTEGER :: year,month,date
    INTEGER :: block1id,block2id,bndsid
    INTEGER :: lnid,ltid,tauid,prsid,tauid2,prsid2,colid,tid,phaseid
    INTEGER, PARAMETER :: indef = 40,indes = 400
    INTEGER block1(indef),block2(indes)
    INTEGER, PARAMETER :: bnds = 2

    ! make netcdf file with dimensions and global attributes. I could
    ! make this more general.
    A  = model%aux
    lx = A%nlon
    ly = A%nlat
    lz = A%nlev
    lt = A%ntlen
    year    = O%epoch%year
    month   = O%epoch%month
    n_tbins = O%ctp_tau%n_tbins
    n_pbins = O%ctp_tau%n_pbins
    nchan   = O%sim%nchannels
    V       = O%vars 
    ! restart these or I'll run out
    ncid = 0

    ! Create the file. 
     CALL check( nf90_create(sim%netcdf_file,nf90_clobber,ncid) )
    ! ========================
    !  GLOBAL ATTRIBUTES

    ! institute etc
    CALL check( nf90_put_att(ncid,nf90_global,"Contact" ,"Salomon.Eliasson@smhi.se"))
    CALL check( nf90_put_att(ncid,nf90_global,"Institute" ,"SMHI"))

    ! creation
    CALL DATE_AND_TIME(values=date_time)
    utc(1) = FLOOR((REAL(60*date_time(5)+date_time(6)-date_time(4)))/60)
    utc(2) = 60*date_time(5)+date_time(6)-date_time(4)  - 60*utc(1)
    WRITE(creation,'(I4,A1,I2,A1,I2,3x,I2,A1,I2,2x,A3)')&
         date_time(1),'-',date_time(2),'-',date_time(3),&
         utc(1),':',utc(2),'UTC'
    CALL check( nf90_put_att(ncid,nf90_global,"Creation date",creation))

    ! technical specs
    string='Cloud_cci'
    CALL check( nf90_put_att(ncid,nf90_global,'simulated dataset',TRIM(string)))
    CALL check( nf90_put_att(ncid,nf90_global,'based_on_model',O%model))

    SELECT CASE (O%subsampler)
    CASE (0)
       string = "maximum random scops from COSP"
    CASE (1)
       string = "maximum random from SIMFERA"
    CASE (2)
       string = "Mc ICA"
    END SELECT
    CALL check( nf90_put_att(ncid,nf90_global,'subsampler_used',TRIM(string)))
    CALL check( nf90_put_att(ncid,nf90_global,'day_since_1_Jan_1970',sim%epoch))
    CALL check( nf90_put_att(ncid,nf90_global,'number of subcolumns',O%ncols))
    SELECT CASE (O%cloudMicrophys%cf_method)
    CASE (0)
       WRITE(string,'(A,1x,f3.1)') &
            "globally constant optical depth cloud treshold at tau=",O%cloudMicrophys%tau_min
    CASE (1)
       string = "Cloud mask created using global gridded detection limits"
    CASE (2)
       string = "probability of detection based on optical depth bins"
    END SELECT
    CALL check( nf90_put_att(ncid,nf90_global,'cloud_mask_method',TRIM(string)))
    IF (O%L2b%doL2bSampling) THEN
       WRITE(string,'(F4.1)') sat%overpass
       CALL check( nf90_put_att(ncid,nf90_global,'equatorial_overpass_time',TRIM(string)))
       CALL check( nf90_put_att(ncid,nf90_global,'satellite_node',O%L2b%node))
       CALL check( nf90_put_att(ncid,nf90_global,'satellite',O%L2b%satellite))
    END IF

    CALL check( nf90_put_att(ncid,nf90_global,'year',O%epoch%year))
    CALL check( nf90_put_att(ncid,nf90_global,'month',O%epoch%month))
    CALL check( nf90_put_att(ncid,nf90_global,'day',iday))

    ! ------------------------
    !  DIMENSIONS
    CALL check( nf90_def_dim(ncid,'time',nf90_unlimited,tid) )
    CALL check( nf90_def_dim(ncid,'bnds',bnds,bndsid) )
    IF (O%model.EQ.'racmo') THEN
       CALL check( nf90_def_dim(ncid,'rlon',lx,lnid) )
       CALL check( nf90_def_dim(ncid,'rlat',ly,ltid) )
    ELSE
       CALL check( nf90_def_dim(ncid,'lons',lx,lnid) )
       CALL check( nf90_def_dim(ncid,'lats',ly,ltid) )
    ENDIF
    CALL check( nf90_def_dim(ncid,'ncols',O%ncols,colid) )

    ! The number of interfaces
    CALL check( nf90_def_dim(ncid,'tau_bins',n_tbins,tauid) )
    CALL check( nf90_def_dim(ncid,'pres_bins',n_pbins,prsid) )
    CALL check( nf90_def_dim(ncid,'tau_bnds',n_tbins+1,tauid2) )
    CALL check( nf90_def_dim(ncid,'pres_bnds',n_pbins+1,prsid2) )
    CALL check( nf90_def_dim(ncid,'phase',2,phaseid) )

    IF (O%model.EQ.'racmo') THEN
       CALL check( nf90_def_dim(ncid,'nblock1',indef,block1id) )
       CALL check( nf90_def_dim(ncid,'nblock2',indes,block2id) )
    ENDIF

    CALL check( nf90_enddef(ncid) )

    ! Note: In addVariable,make sure to list the dimension id's in the
    !       right order!  

    ! ===============================
    ! Auxiliary data
    CALL addVariable(ncid,'time',     O,A,tid,       sim=sim)
    CALL addVariable(ncid,'time_bnds',O,A,bndsid,tid,sim=sim)
    IF (O%model.EQ.'racmo') THEN
       ! ... add fixed block1-info plotting purposes and gtx
       block1=0
       block1(1:24)=[24,0,0,1,96,1,255,128,254,105,0,0,2008,9,1,12,1,0,0,1,0,0,0,0]
       date=sim%date
       block1(13)=date/10000
       block1(14)=date/100-block1(13)*100
       block1(15)=date-block1(13)*10000-block1(14)*100
       block1(16)=12
       ! ... add fixed block2-info plotting purposes and gtx
       block2=0
       block2( 1)=50
       block2( 6)=10  !  rotated grid 
       block2( 7)=lx
       block2( 9)=ly
       block2(14)=NINT((A%rlon( 1)) * 1000.)
       block2(21)=NINT((A%rlon(lx)) * 1000.)
       block2(11)=NINT((A%rlat( 1)) * 1000.)
       block2(17)=136
       block2(18)=NINT((A%rlat(ly)) * 1000.)
       block2(24)=NINT((A%rlon(2)-A%rlon(1))*1000.)
       block2(26)=NINT((A%rlat(2)-A%rlat(1))*1000.)
       block2(28)=64  ! latitude from south to north
       block2(33)=NINT((A%float_rotpolat * (-1.)) * 1000.)  ! convert NP to SP latitude
       block2(36)=NINT((A%float_rotpolon -180.0 ) * 1000.)  ! convert NP to SP longitude
       block2(39)=1
       block2(40)=1

       CALL addVariable(ncid,'dtg',      O,A,tid,       sim=sim)
       CALL addVariable(ncid,'block1',   O,A,block1id,  block1=block1)
       CALL addVariable(ncid,'block2',   O,A,block2id,  block2=block2)
       CALL addVariable(ncid,'areacella',O,A,ltid,      M=model)
       CALL addVariable(ncid,'rlon',     O,A,lnid,      M=model)
       CALL addVariable(ncid,'rlat',     O,A,ltid,      M=model)
       CALL addVariable(ncid,'lon',      O,A,lnid,ltid, M=model)
       CALL addVariable(ncid,'lat',      O,A,lnid,ltid, M=model)
    ELSE
       CALL addVariable(ncid,'lon',      O,A,lnid,      M=model)
       CALL addVariable(ncid,'lat',      O,A,ltid,      M=model)
    END IF

    IF (V%land_sea)        CALL addVariable(ncid,'lsm',O,A,lnid,ltid,M=model)
    IF (V%solzen)          CALL addVariable(ncid,'solzen',O,A,lnid,ltid,S=S)
    IF (V%time_of_day)     CALL addVariable(ncid,'time_of_day',O,A,lnid,ltid,tid,sim=sim)
    IF (V%hist2d_cot_ctp) THEN 
       CALL addVariable(ncid,'hist_phase',           O,A,phaseid)
       CALL addVariable(ncid,'hist2d_ctp_bin_centre',O,A,prsid)
       CALL addVariable(ncid,'hist2d_cot_bin_centre',O,A,tauid)
       CALL addVariable(ncid,'hist2d_ctp_bin_border',O,A,prsid2)
       CALL addVariable(ncid,'hist2d_cot_bin_border',O,A,tauid2)
    END IF
    ! END auxiliary
    ! --------------------------------

    ! --------------------
    ! SIMULATORS
    ! --------------------
    
    
    IF (O%sim%doCloud_cci) THEN
       IF (V%cer_ice       ) CALL addVariable(ncid,'cer_ice'       ,O,A,lnid,ltid,tid,cloud_cci=cloud_cci)
       IF (V%cer_liq       ) CALL addVariable(ncid,'cer_liq'       ,O,A,lnid,ltid,tid,cloud_cci=cloud_cci)
       IF (V%cfc           ) CALL addVariable(ncid,'cfc'           ,O,A,lnid,ltid,tid,cloud_cci=cloud_cci)
       IF (V%cfc_low       ) CALL addVariable(ncid,'cfc_low'       ,O,A,lnid,ltid,tid,cloud_cci=cloud_cci)
       IF (V%cfc_mid       ) CALL addVariable(ncid,'cfc_mid'       ,O,A,lnid,ltid,tid,cloud_cci=cloud_cci)
       IF (V%cfc_high      ) CALL addVariable(ncid,'cfc_high'      ,O,A,lnid,ltid,tid,cloud_cci=cloud_cci)
       IF (V%cla_vis006    ) CALL addVariable(ncid,'cla_vis006'    ,O,A,lnid,ltid,tid,cloud_cci=cloud_cci)
       IF (V%cot           ) CALL addVariable(ncid,'cot'           ,O,A,lnid,ltid,tid,cloud_cci=cloud_cci)
       IF (V%cot_ice       ) CALL addVariable(ncid,'cot_ice'       ,O,A,lnid,ltid,tid,cloud_cci=cloud_cci)
       IF (V%cot_liq       ) CALL addVariable(ncid,'cot_liq'       ,O,A,lnid,ltid,tid,cloud_cci=cloud_cci)
       IF (V%cth           ) CALL addVariable(ncid,'cth'           ,O,A,lnid,ltid,tid,cloud_cci=cloud_cci)
       IF (V%ctp           ) CALL addVariable(ncid,'ctp'           ,O,A,lnid,ltid,tid,cloud_cci=cloud_cci)
       IF (V%ctp_log       ) CALL addVariable(ncid,'ctp_log'       ,O,A,lnid,ltid,tid,cloud_cci=cloud_cci)
       IF (V%ctt           ) CALL addVariable(ncid,'ctt'           ,O,A,lnid,ltid,tid,cloud_cci=cloud_cci)
       IF (V%cth_corrected ) CALL addVariable(ncid,'cth_corrected' ,O,A,lnid,ltid,tid,cloud_cci=cloud_cci)
       IF (V%ctp_corrected ) CALL addVariable(ncid,'ctp_corrected' ,O,A,lnid,ltid,tid,cloud_cci=cloud_cci)
       IF (V%ctt_corrected ) CALL addVariable(ncid,'ctt_corrected' ,O,A,lnid,ltid,tid,cloud_cci=cloud_cci)
       IF (V%hist2d_cot_ctp) &
            CALL addVariable(ncid,'hist2d_cot_ctp',O,A,lnid,ltid,tauid,prsid,phaseid,tid,cloud_cci=cloud_cci)
       IF (V%iwp           ) CALL addVariable(ncid,'iwp'           ,O,A,lnid,ltid,tid,cloud_cci=cloud_cci)
       IF (V%lwp           ) CALL addVariable(ncid,'lwp'           ,O,A,lnid,ltid,tid,cloud_cci=cloud_cci)
    END IF

    CALL check (nf90_close(ncid) )

  END SUBROUTINE make_netcdf

  SUBROUTINE addVariable(ncid,variable,O,A,&
       id1,id2,id3,id4,id5,id6,&
       M,S,sim,cloud_cci,&
       block1,block2)

    IMPLICIT NONE

    INTEGER, PARAMETER                              :: indef = 40,indes = 400

    INTEGER,                             INTENT(in) :: ncid
    CHARACTER(len=*),                    INTENT(in) :: variable
    TYPE(name_list),                     INTENT(in) :: O
    TYPE(model_aux),                     INTENT(in) :: A
    INTEGER,                             INTENT(in) :: id1
    TYPE(model_type),          OPTIONAL, INTENT(in) :: M
    TYPE(subset),              OPTIONAL, INTENT(in) :: S
    TYPE(satellite_simulator), OPTIONAL, INTENT(in) :: sim
    TYPE(cloud_cci_type),      OPTIONAL, INTENT(in) :: cloud_cci
    INTEGER,                   OPTIONAL, INTENT(in) :: block1(indef),block2(indes)
    INTEGER,                   OPTIONAL, INTENT(in) :: id2,id3,id4,id5,id6
    INTEGER                                         :: ndim,id(6)
    REAL(wp)                                        :: ztime,ztime_bnds(2)
    INTEGER                                         :: varid
    REAL(wp), ALLOCATABLE                           :: lon(:),lat(:)

    ! count the number of dimensions we are dealing with
    ndim=1
    id(1)=id1
    IF (PRESENT(id2)) THEN
       id(2)=id2
       ndim=ndim+1
    END IF
    IF (PRESENT(id3)) THEN
       id(3)=id3
       ndim=ndim+1
    END IF
    IF (PRESENT(id4)) THEN
       id(4)=id4
       ndim=ndim+1
    END IF
    IF (PRESENT(id5)) THEN
       id(5)=id5
       ndim=ndim+1
    END IF
    IF (PRESENT(id6)) THEN
       id(6)=id6
       ndim=ndim+1
    END IF

    ! get the varid
    IF (PRESENT(M)) THEN
       CALL handle_variable_attributes(ncid,variable,ndim,id,varid,O,M)
    ELSE
       CALL handle_variable_attributes(ncid,variable,ndim,id,varid,O)
    END IF

    !----------------------------
    ! NOW THE DATA
    ! 

    ! Reshape the data if I have to

    SELECT CASE (variable)
       ! ----------------
       ! auxiliary data
       !
    CASE ('areacella')
       CALL putVar(ncid,varid,A,A%areacella)
    CASE ('block1')
       CALL check(nf90_put_var(ncid,varid,block1))
    CASE ('block2')
       CALL check(nf90_put_var(ncid,varid,block2))
    CASE ('dtg')
       CALL check(nf90_put_var(ncid,varid,sim%dtg))   
    CASE ('lon')
       IF (O%model.EQ.'racmo') THEN
          CALL putVar(ncid,varid,A,data2=A%lon)
       ELSE
          ALLOCATE(lon(A%nlon))
          lon = A%lon(:,1)
          CALL putVar(ncid,varid,A,lon)
          DEALLOCATE(lon)
       ENDIF
    CASE ('lat')
       IF (O%model.EQ.'racmo') THEN
          CALL putVar(ncid,varid,A,data2=A%lat)
       ELSE
          ALLOCATE(lat(A%nlat))
          lat = A%lat(1,:)
          CALL putVar(ncid,varid,A,lat)
          DEALLOCATE(lat)
       ENDIF
    CASE ('lsm')
       CALL putVar(ncid,varid,A,A%lsm)
    CASE ('hist2d_ctp_bin_border')
       CALL putVar(ncid,varid,A,O%ctp_tau%pbin_edges)
    CASE ('hist2d_cot_bin_border')
       CALL putVar(ncid,varid,A,O%ctp_tau%tbin_edges)
    CASE ('hist2d_ctp_bin_centre')
       CALL putVar(ncid,varid,A,O%ctp_tau%pbin_centre)
    CASE ('hist2d_cot_bin_centre')
       CALL putVar(ncid,varid,A,O%ctp_tau%tbin_centre)
    CASE ('hist_phase')
       CALL putVar(ncid,varid,A,O%ctp_tau%hist_phase)
    CASE ('POD_tau_bin_centers')
       CALL putVar(ncid,varid,A,O%sim_aux%POD_tau_bin_centers)
    CASE ('POD_tau_bin_edges')
       CALL putVar(ncid,varid,A,O%sim_aux%POD_tau_bin_edges)
    CASE ('rlon')
       CALL putVar(ncid,varid,A,A%rlon)
    CASE ('rlat')
       CALL putVar(ncid,varid,A,A%rlat)
    CASE ('solzen')
       IF (PRESENT(S)) THEN
          CALL putVar(ncid,varid,A,S%solzen)
       END IF
    CASE ('time')
       ztime = float( sim%epoch ) + 0.5 ! center in the middle of the day to comply with dtg
       CALL check(nf90_put_var(ncid,varid,ztime))
    CASE ('time_bnds')
       ztime_bnds(1) = float( sim%epoch )       ! start of day at 00 UTC
       ztime_bnds(2) = float( sim%epoch ) + 1.0 ! end of day at 24 UTC
       CALL putVar(ncid,varid,A,ztime_bnds)
    CASE ('time_of_day')
       CALL putVar(ncid,varid,A,sim%time_of_day)

    END SELECT

    ! ---------------
    ! Cloud_cci
    IF (O%sim%doCloud_cci) THEN
       SELECT CASE (variable)
       CASE ('cer_ice')
          CALL putVar(ncid,varid,A,cloud_cci%av%cer_ice)
       CASE ('cer_liq')
          CALL putVar(ncid,varid,A,cloud_cci%av%cer_liq)
       CASE ('cfc')
          CALL putVar(ncid,varid,A,cloud_cci%av%cfc)
       CASE ('cfc_low')
          CALL putVar(ncid,varid,A,cloud_cci%av%cfc_low)
       CASE ('cfc_mid')
          CALL putVar(ncid,varid,A,cloud_cci%av%cfc_mid)
       CASE ('cfc_high')
          CALL putVar(ncid,varid,A,cloud_cci%av%cfc_high)
       CASE ('cla_vis006')
          CALL putVar(ncid,varid,A,cloud_cci%av%cla_vis006)
       CASE ('cot')
          CALL putVar(ncid,varid,A,cloud_cci%av%tau)
       CASE ('cot_ice')
          CALL putVar(ncid,varid,A,cloud_cci%av%cot_ice)
       CASE ('cot_liq')
          CALL putVar(ncid,varid,A,cloud_cci%av%cot_liq)
       CASE ('cth')
          CALL putVar(ncid,varid,A,cloud_cci%av%cth)
       CASE ('ctp')
          CALL putVar(ncid,varid,A,cloud_cci%av%ctp)
       CASE ('ctp_log')
          CALL putVar(ncid,varid,A,cloud_cci%av%ctp_log)
       CASE ('ctt')
          CALL putVar(ncid,varid,A,cloud_cci%av%ctt)
       CASE ('cth_corrected')
          CALL putVar(ncid,varid,A,cloud_cci%av%cth_c)
       CASE ('ctp_corrected')
          CALL putVar(ncid,varid,A,cloud_cci%av%ctp_c)
       CASE ('ctt_corrected')
          CALL putVar(ncid,varid,A,cloud_cci%av%ctt_c)
       CASE ('iwp')
          CALL putVar(ncid,varid,A,cloud_cci%av%iwp)
       CASE ('lwp')
          CALL putVar(ncid,varid,A,cloud_cci%av%lwp)
       CASE ('hist2d_cot_ctp')
          CALL putVar(ncid,varid,A,data4=cloud_cci%av%hist2d_cot_ctp)
       END SELECT
    END IF

  END SUBROUTINE addVariable

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
    PRINT 99,ndim

  END SUBROUTINE putVar


  SUBROUTINE handle_variable_attributes(ncid,variable,ndim,id,varid,O,M)

    USE data_check, ONLY: height_max,p_tau_max,P_max,tau_max,T_min,T_max
    
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
    CASE ('detection_limit')
       units = ""
       WRITE(description,'(A)') "optical depth values above which 50%&
            & or more of clouds are detected"
       WRITE(long_name,'(A)') "detection_limit"
       valid_min=0
       valid_max=100
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
            & top (where tau=",O%cloudMicrophys%tau_equivRadCldTop,&
            "), but also where tau >",O%cloudMicrophys%tau_min
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
            &of finding the cloud top (where tau=",&
            O%cloudMicrophys%tau_equivRadCldTop,&
            "), but also where tau >",O%cloudMicrophys%tau_min
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
            & top (where tau=",O%cloudMicrophys%tau_equivRadCldTop,&
            "), but also where tau >",O%cloudMicrophys%tau_min
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

  END SUBROUTINE handle_variable_attributes

END MODULE simulator_netcdf_module
