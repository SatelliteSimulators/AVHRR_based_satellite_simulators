MODULE clara_netcdf
  !
  ! Salomon.Eliasson@smhi.se
  !
  ! JFM, 15/07/2014: renamed + added some output fields
  !

  USE AUXILIARY_FUNCTIONS,  ONLY: &
       HANDLE_VARIABLE_ATTRIBUTES,&
       PUTVAR
  USE CLARA_M,              ONLY: CLARA_TYPE
  USE COSP_KINDS,           ONLY: WP
  USE MODEL_INPUT,          ONLY: MODEL_TYPE, MODEL_AUX
  USE MY_NETCDFTOOLS,       ONLY: CHECK
  USE NETCDF,               ONLY: &
       NF90_CLOBBER,              &
       NF90_CLOSE,                &
       NF90_CREATE,               &
       NF90_UNLIMITED,            &
       NF90_DEF_DIM,              &
       NF90_ENDDEF,               &
       NF90_GLOBAL,               &
       NF90_PUT_ATT,              &
       NF90_PUT_VAR
  USE NAMELIST_INPUT,            ONLY: NAME_LIST,VARIABLESCONTAINER
  USE SATELLITE_SPECS,           ONLY: SATELLITE
  USE SIMULATOR_VARIABLES,       ONLY: SATELLITE_SIMULATOR
  USE SIMULATOR_INPUT_VARIABLES, ONLY: SUBSET

  IMPLICIT NONE

  PUBLIC  :: MAKE_NETCDF
  PRIVATE :: ADDVARIABLE

CONTAINS

  SUBROUTINE MAKE_NETCDF(model,sim,O,iday,S,sat,clara)

    !
    ! tailor-made module to save model and simulator output to a netcdf file
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
    TYPE(clara_type),                    INTENT(in) :: clara

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
    INTEGER :: lnid,ltid,lvlid,tauid,prsid,tauid2,prsid2,&
         colid,tid,phaseid,cldthickid,sunid,podid,podid2
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
    CALL check( nf90_put_att(ncid,nf90_global,"Simulator version" ,O%simVersionNumber))
    CALL check( nf90_put_att(ncid,nf90_global,"Simulated Climate data record" ,O%CDR))
    CALL check( nf90_put_att(ncid,nf90_global,'model',O%model))

    ! creation
    CALL DATE_AND_TIME(values=date_time)
    utc(1) = FLOOR((REAL(60*date_time(5)+date_time(6)-date_time(4)))/60)
    utc(2) = 60*date_time(5)+date_time(6)-date_time(4)  - 60*utc(1)
    WRITE(creation,'(I4,A1,I2,A1,I2,3x,I2,A1,I2,2x,A3)')&
         date_time(1),'-',date_time(2),'-',date_time(3),&
         utc(1),':',utc(2),'UTC'
    CALL check( nf90_put_att(ncid,nf90_global,"Creation date",creation))


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
            "globally constant optical depth cloud threshold at tau=",O%cloudMicrophys%tau_min
    CASE (1)
       string = "probability of detection based on optical depth bins"
    END SELECT
    CALL check( nf90_put_att(ncid,nf90_global,'cloud_mask_method',TRIM(string)))

    CALL check( nf90_put_att(ncid,nf90_global,'satellite',O%L2b%satellite))
    IF (O%L2b%doL2bSampling) THEN
       WRITE(string,'(F4.1)') sat%overpass
       CALL check( nf90_put_att(ncid,nf90_global,'equatorial_overpass_time',TRIM(string)))
       CALL check( nf90_put_att(ncid,nf90_global,'satellite_node',O%L2b%node))
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
    CALL check( nf90_def_dim(ncid,'cloud_categories',4,cldthickid) )

    ! The number of interfaces
    CALL check( nf90_def_dim(ncid,'tau_bins',n_tbins,tauid) )
    CALL check( nf90_def_dim(ncid,'pres_bins',n_pbins,prsid) )
    CALL check( nf90_def_dim(ncid,'tau_bnds',n_tbins+1,tauid2) )
    CALL check( nf90_def_dim(ncid,'pres_bnds',n_pbins+1,prsid2) )
    CALL check( nf90_def_dim(ncid,'phase',2,phaseid) )
    IF (O%cloudMicrophys%cf_method.EQ.1) THEN
       CALL check( nf90_def_dim(ncid,'pod_bins',&
            SIZE(O%sim_aux%POD_tau_bin_centers),podid) )
       CALL check( nf90_def_dim(ncid,'pod_bnds',&
            SIZE(O%sim_aux%POD_tau_bin_edges),podid2) )
    END IF
    IF (O%cloudMicrophys%cf_method.GT.0) THEN
       CALL check( nf90_def_dim(ncid,'sunlit',2,sunid) )
    END IF
    IF  (O%sim%doModel) THEN
       CALL check( nf90_def_dim(ncid,'levels',lz,lvlid) )
    END IF
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
    IF (ALLOCATED(O%sim_aux%POD_layers)) THEN
       CALL addVariable(ncid,'POD_tau_bin_centers',O,A,podid ,clara=clara)
       CALL addVariable(ncid,'POD_tau_bin_edges',  O,A,podid2,clara=clara)
    END IF
    ! END auxiliary
    ! --------------------------------

    ! --------------------
    ! SIMULATORS
    ! --------------------

    IF (V%cfc      ) CALL addVariable(ncid,'cfc'      ,O,A,lnid,ltid,tid,clara=clara)
    IF (V%cfc_day  ) CALL addVariable(ncid,'cfc_day'  ,O,A,lnid,ltid,tid,clara=clara)
    IF (V%cfc_low  ) CALL addVariable(ncid,'cfc_low'  ,O,A,lnid,ltid,tid,clara=clara)
    IF (V%cfc_mid  ) CALL addVariable(ncid,'cfc_mid'  ,O,A,lnid,ltid,tid,clara=clara)
    IF (V%cfc_high ) CALL addVariable(ncid,'cfc_high' ,O,A,lnid,ltid,tid,clara=clara)
    IF (V%cflag_tot) CALL addVariable(ncid,'cflag_tot',O,A,lnid,ltid,cldthickid,tid,clara=clara)
    IF (V%cot_ice  ) CALL addVariable(ncid,'cot_ice'  ,O,A,lnid,ltid,tid,clara=clara)
    IF (V%cot_liq  ) CALL addVariable(ncid,'cot_liq'  ,O,A,lnid,ltid,tid,clara=clara)
    IF (V%cth      ) CALL addVariable(ncid,'cth'      ,O,A,lnid,ltid,tid,clara=clara)
    IF (V%ctp      ) CALL addVariable(ncid,'ctp'      ,O,A,lnid,ltid,tid,clara=clara)
    IF (V%ctp_log  ) CALL addVariable(ncid,'ctp_log'  ,O,A,lnid,ltid,tid,clara=clara)
    IF (V%ctt      ) CALL addVariable(ncid,'ctt'      ,O,A,lnid,ltid,tid,clara=clara)
    IF (V%iwp      ) CALL addVariable(ncid,'iwp'      ,O,A,lnid,ltid,tid,clara=clara)
    IF (V%lwp      ) CALL addVariable(ncid,'lwp'      ,O,A,lnid,ltid,tid,clara=clara)
    IF (V%ref_ice  ) CALL addVariable(ncid,'ref_ice'  ,O,A,lnid,ltid,tid,clara=clara)
    IF (V%ref_liq  ) CALL addVariable(ncid,'ref_liq'  ,O,A,lnid,ltid,tid,clara=clara)

    IF (V%hist2d_cot_ctp) &
         CALL addVariable(ncid,'hist2d_cot_ctp',O,A,lnid,ltid,tauid,prsid,phaseid,tid,clara=clara)
    IF (ALLOCATED(O%sim_aux%POD_layers)) &
         CALL addVariable(ncid,'POD_layers',O,A,lnid,ltid,podid,sunid,tid,clara=clara)

    CALL check (nf90_close(ncid) )

  END SUBROUTINE make_netcdf

  SUBROUTINE addVariable(ncid,variable,O,A,&
       id1,id2,id3,id4,id5,id6,&
       M,S,sim,clara,block1,block2)

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
    TYPE(clara_type),          OPTIONAL,  INTENT(in) :: clara
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
       CALL putVar(ncid,varid,A,S%solzen)
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
    ! CLARA
    !
    SELECT CASE (variable)
    CASE ('cfc')
       CALL putVar(ncid,varid,A,clara%av%cfc)
    CASE ('cfc_day')
       CALL putVar(ncid,varid,A,clara%av%cfc_day)
    CASE ('cfc_low')
       CALL putVar(ncid,varid,A,clara%av%cfc_low)
    CASE ('cfc_mid')
       CALL putVar(ncid,varid,A,clara%av%cfc_mid)
    CASE ('cfc_high')
       CALL putVar(ncid,varid,A,clara%av%cfc_high)
    CASE ('cflag_tot')
       CALL putVar(ncid,varid,A,data2i=clara%av%cflag_tot)
    CASE ('cot_ice')
       CALL putVar(ncid,varid,A,clara%av%cot_ice)
    CASE ('cot_liq')
       CALL putVar(ncid,varid,A,clara%av%cot_liq)
    CASE ('cth')
       CALL putVar(ncid,varid,A,clara%av%cth)
    CASE ('ctp')
       CALL putVar(ncid,varid,A,clara%av%ctp)
    CASE ('ctp_log')
       CALL putVar(ncid,varid,A,clara%av%ctp_log)
    CASE ('ctt')
       CALL putVar(ncid,varid,A,clara%av%ctt)
    CASE ('iwp')
       CALL putVar(ncid,varid,A,clara%av%iwp)
    CASE ('lwp')
       CALL putVar(ncid,varid,A,clara%av%lwp)
    CASE ('POD_layers')
       CALL putVar(ncid,varid,A,data3=O%sim_aux%POD_layers)
    CASE ('hist2d_cot_ctp')
       CALL putVar(ncid,varid,A,data4=clara%av%hist2d_cot_ctp)
    CASE ('ref_liq')
       CALL putVar(ncid,varid,A,clara%av%ref_liq)
    CASE ('ref_ice')
       CALL putVar(ncid,varid,A,clara%av%ref_ice)
    END SELECT

  END SUBROUTINE addVariable

END MODULE clara_netcdf
