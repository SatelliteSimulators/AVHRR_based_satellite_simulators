MODULE namelist_input

  ! putting the namelist into a structure instead
  !
  ! salomon.eliasson@smhi.se

  USE cosp_kinds, ONLY: wp
  USE handy,      ONLY: check_file
  USE my_maths,   ONLY: daysinmonth
  USE optics_m,   ONLY: simulator_aux

  IMPLICIT NONE

  PUBLIC :: &
       common_namelist,       &
       deallocate_namelist,   &
       namelist_ctp_tau,      &
       namelist_daynight,     &
       namelist_microphys,    &
       name_list

  TYPE paths
     ! ------- Paths -------
     !
     !  dailyFiles : Whether or not the input files are daily or monthly
     !  model_input:
     !             Is the full path to the model input files. Use
     !             the following keywords in the regular expression
     !             #SIM='sim',#MODEL=options%model(A),#Y4=year(i4),#M2=month(i2),
     !             #D2=day(i2),#Y2=year(i2),#M1=month(i1),#D1=day(i1),#UTC=utc(i2),
     !             #NODE=string,#VERSION=string,#STRING=string
     !             (e.g. "/some/where/#Y4/#M2/#D2/file_#Y4#M2#D2.nc", where
     !               Y4=year (i0.4),M2=month (i0.2), and D2=day (i0.2)
     !  sim_output :
     !             Same as story as model_input_regexp, but for the output NETCDF files

     CHARACTER(len=1000) :: sim_output, model_input
     LOGICAL             :: dailyFiles
  END TYPE paths
  TYPE epoch
     ! ------ EPOCH -------
     !  year	: year to run
     !  month	: month to run
     !  day	: day to loop over

     INTEGER :: year, month, day
  END TYPE epoch
  TYPE L2b

     ! ----- L2b ----------
     !

     ! This information is used to create a level2b-like output. The
     ! satellite, node, and
     ! year,month (provided to options%epoch) information is used to
     ! create the equator crossing times. The out put data will then
     ! have the same local time everywhere, the same as the satellites
     ! overpass times.
     !
     ! satellite     : name of satellite, e.g., 'noaa15'
     ! node          : 'asc','dec', 'all',''
     !
     ! doL2bsampling : .true. = do the above
     !
     ! interpolate   : .true. = don't just sample the equatorial overpass time,
     !                          linearly interpolate it so that the
     !                          local time matches this time everywhere
     !
     !
     LOGICAL           :: doL2bSampling
     LOGICAL           :: interpolate
     CHARACTER(len=20) :: satellite
     CHARACTER(len=3)  :: node
  END TYPE L2b
  TYPE cloudMicrophysics
     ! ------ cloudMicrophys ----------
     !
     ! 'cf_method'         technique to determine cloud fraction.
     !                     0 = global constant (tau_min)
     !                     1 = probability of detection based on optical depth bins
     !
     ! 'tau_min'            = Cloud optical depth threshold for CLARA-A1. This
     !                        should be listed in the namelist


     INTEGER  :: cf_method
     REAL(wp) :: tau_min

  END TYPE cloudMicrophysics
  TYPE illumination
     ! ------ daynight ----------
     !
     !   daylim   = maximum solar zenith angle for daylight
     !   nightlim = minimum solar zenith angle for nighttime
     real(wp):: daylim,nightlim
  END TYPE illumination
  TYPE ctp_tau
     ! -------- ctp_tau_hist_dim ----------
     !
     ! CTP-TAU diagrams. Dimensions and bin edges must be consistent
     ! (dim=length(edges)-1)
     !
     ! n_tbins	: number of optical thickness bins for P_tau
     ! 		 histograms
     ! n_pbins 	: number of pressure bins for P_tau histograms
     !
     ! --------- ctp_tau_hist_vec ---------
     !
     !  Bin edges values for cloud optical depth and cloud top pressure
     !  for CPT-tau histograms
     !
     ! tbin_edges      : vector for tau bin edges
     ! pbins_edges     : vector for CTP bin edges [Pa]
     ! hist_phase      : 0 (ice), 1 (liquid)
     !
     !
     REAL(wp), ALLOCATABLE :: tbin_edges(:),tbin_centre(:)
     REAL(wp), ALLOCATABLE :: pbin_edges(:),pbin_centre(:)
     REAL(WP) :: hist_phase(2)
     INTEGER  :: n_tbins, n_pbins
  END TYPE ctp_tau
  TYPE namelist_sim
     LOGICAL :: doClara
     LOGICAL :: doCloud_cci
     LOGICAL :: doISCCP
     LOGICAL :: doRTTOV
     LOGICAL :: doModel

     !    ! number of frequencies (only for RTTOV)
     INTEGER :: nchannels, sensor
     INTEGER, ALLOCATABLE :: channels(:)

     ! 0 = CDK model (ISCCP Tb), 1 = rttov Tb's
     INTEGER :: Tb
  END TYPE namelist_sim
  TYPE variablesContainer
     !
     ! These are all the available variables that can be output from
     ! the simulator
     !
     ! Cloud cover
     LOGICAL :: cfc,cfc_day,cfc_low,cfc_mid,cfc_high,icf,lcf
     ! Cloud top products
     LOGICAL :: cth,ctp,ctp_log,ctt,cth_corrected,ctp_corrected,ctt_corrected,cflag_tot
     ! Cloud effective radius
     LOGICAL :: ireff,lreff,ireff3D,lreff3D,ref_ice,ref_liq,cer_ice,cer_liq
     ! Cloud visible optical thickness
     LOGICAL :: albedo,itau,ltau,tau,cot_ice,cot_liq, cla_vis006, cot
     ! Cloud water path
     LOGICAL :: iwp,lwp
     ! histograms
     LOGICAL :: hist2d_cot_ctp
     ! ancillary data
     LOGICAL :: solzen,land_sea,time_of_day
     ! direct model variables
     LOGICAL :: CI,SKT,TCC,TCWV
     ! brightness temperature
     LOGICAL :: Tb,Tb_clr,tau_subcolumn,Tb_subcolumn

  END TYPE variablesContainer
  TYPE name_list
     TYPE(paths)             :: paths
     TYPE(epoch)             :: epoch
     TYPE(L2b)               :: L2b
     TYPE(cloudMicrophysics) :: cloudMicrophys
     TYPE(illumination)      :: daynight
     TYPE(ctp_tau)           :: ctp_tau
     TYPE(namelist_sim)      :: sim
     TYPE(variablesContainer):: vars
     TYPE(simulator_aux)     :: sim_aux

     ! ---- debug --------
     !
     ! dbg		: debug flag, 0 = no debug; -1 = use earlier
     ! 		: version of sim_clara; >0 = debug
     !		: (integer specifies npts to print)
     !
     ! 'ncols'              = number of subcolumns used in the model grid
     ! 'namelist_file'
     ! 'model'              = string name of the model
     ! 'overwrite_existing' = If .TRUE. it will overwite existing simulated files
     ! 'subsampler'         = 0 (scops (COSP default)), 1 (DWD SIMFERA), 2 McICA
     ! 'simVersionNumber'   = version of this software
     ! 'CDR'                = Name of the climate data record
     !                        (cloud_cci, clara_a2, clara_a3)

     INTEGER                 :: dbg
     CHARACTER(len=100)      :: model
     CHARACTER(len=1000)     :: namelist_file
     INTEGER                 :: ncols
     LOGICAL                 :: overwrite_existing
     INTEGER                 :: subsampler
     CHARACTER(len=3)        :: simVersionNumber
     CHARACTER(len=10)       :: CDR
  END TYPE name_list

CONTAINS

  ! --------------------
  ! SETUP SIMULATORS
  ! --------------------

  !
  ! COMMON namelist
  !
  SUBROUTINE common_namelist(x,file)

    CHARACTER(len=*),INTENT(in)    :: file
    TYPE(name_list), INTENT(inout) :: x
    INTEGER                        :: year, month, day
    CHARACTER(len=10)              :: model_name
    CHARACTER(len=1000)            :: model_input,sim_output
    CHARACTER(len=10)              :: CDR
    INTEGER                        :: dbg,subsampler
    LOGICAL                        :: exists,dailyFiles,overwrite_existing,use_satellite

    NAMELIST/epoch/year,month,day
    NAMELIST/paths/CDR,model_name,model_input,sim_output,dailyFiles
    NAMELIST/other/dbg,overwrite_existing,subsampler,use_satellite

    day=-999
    exists = .TRUE.
    INQUIRE(FILE=TRIM(file), EXIST=exists)
    IF (.NOT. exists) THEN
       WRITE (*,*) "namelist file is missing. It should be the first input &
            &argument to the executable (e.g. >> ./simulator.x namelist)"
       STOP "stopped in namelist_input.F90"
    END IF

    OPEN(10,file=file,status='old')

    READ(10,epoch)
    REWIND(10)
    READ(10,paths)
    REWIND(10)
    READ(10,other)
    REWIND(10)
    CLOSE(10)

    X%epoch%year          = year
    X%epoch%month         = month
    X%epoch%day           = day
    X%namelist_file       = TRIM(file)
    X%CDR                 = TRIM(CDR)
    X%model               = TRIM(model_name)
    X%paths%model_input   = TRIM(model_input)
    X%paths%sim_output    = TRIM(sim_output)
    X%paths%dailyFiles    = dailyFiles
    X%dbg                 = dbg
    X%subsampler          = subsampler
    X%overwrite_existing  = overwrite_existing

    IF (use_satellite) THEN
       CALL namelist_satellite(x,file)
    ELSE
        X%L2b%doL2bSampling    = .FALSE.
        X%L2b%interpolate      = .FALSE.
        X%L2b%node             = ''
        X%L2b%satellite        = 'none'
    END IF

  END SUBROUTINE common_namelist

  SUBROUTINE deallocate_namelist(options)

    IMPLICIT NONE

    TYPE(name_list),INTENT(inout):: options

    IF (ALLOCATED(options%ctp_tau%tbin_edges)) THEN
       DEALLOCATE(options%ctp_tau%tbin_edges,&
            options%ctp_tau%pbin_edges)
    END IF

  END SUBROUTINE deallocate_namelist

  SUBROUTINE initialise_variable_flag(var1,var2,var3,var4,var5,var6,&
       var7,var8,var9,var10,var11,var12,var13,var14,var15,var16,var17,&
       var18,var19,var20,var21,var22,var23,var24,var25,var26,var27,&
       var28,var29,var30,var31,var32)


    LOGICAL,INTENT(inout),OPTIONAL :: var1,var2,var3,var4,var5,var6,&
       var7,var8,var9,var10,var11,var12,var13,var14,var15,var16,var17,&
       var18,var19,var20,var21,var22,var23,var24,var25,var26,var27,&
       var28,var29,var30,var31,var32

    IF (PRESENT(var1 )) var1 =.FALSE.
    IF (PRESENT(var2 )) var2 =.FALSE.
    IF (PRESENT(var3 )) var3 =.FALSE.
    IF (PRESENT(var4 )) var4 =.FALSE.
    IF (PRESENT(var5 )) var5 =.FALSE.
    IF (PRESENT(var6 )) var6 =.FALSE.
    IF (PRESENT(var7 )) var7 =.FALSE.
    IF (PRESENT(var8 )) var8 =.FALSE.
    IF (PRESENT(var9 )) var9 =.FALSE.
    IF (PRESENT(var10)) var10=.FALSE.
    IF (PRESENT(var11)) var11=.FALSE.
    IF (PRESENT(var12)) var12=.FALSE.
    IF (PRESENT(var13)) var13=.FALSE.
    IF (PRESENT(var14)) var14=.FALSE.
    IF (PRESENT(var15)) var15=.FALSE.
    IF (PRESENT(var16)) var16=.FALSE.
    IF (PRESENT(var17)) var17=.FALSE.
    IF (PRESENT(var18)) var18=.FALSE.
    IF (PRESENT(var19)) var19=.FALSE.
    IF (PRESENT(var20)) var20=.FALSE.
    IF (PRESENT(var21)) var21=.FALSE.
    IF (PRESENT(var22)) var22=.FALSE.
    IF (PRESENT(var23)) var23=.FALSE.
    IF (PRESENT(var24)) var24=.FALSE.
    IF (PRESENT(var25)) var25=.FALSE.
    IF (PRESENT(var26)) var26=.FALSE.
    IF (PRESENT(var27)) var27=.FALSE.
    IF (PRESENT(var28)) var28=.FALSE.
    IF (PRESENT(var29)) var29=.FALSE.
    IF (PRESENT(var30)) var30=.FALSE.
    IF (PRESENT(var31)) var31=.FALSE.
    IF (PRESENT(var32)) var32=.FALSE.
  END SUBROUTINE initialise_variable_flag

  !
  ! Satellite
  !
  SUBROUTINE namelist_satellite(x,file)
    CHARACTER(len=*),INTENT(in)    :: file
    TYPE(name_list), INTENT(inout) :: x
    LOGICAL                        :: doL2bSampling, interpolate
    CHARACTER(len=20)              :: sat
    CHARACTER(len=3)               :: node
    NAMELIST/satellite/sat,node,doL2bSampling,interpolate

    OPEN(10,file=file,status='old')
    READ(10,satellite)
    CLOSE(10)

    X%L2b%node                 = node
    X%L2b%doL2bSampling        = doL2bSampling
    X%L2b%interpolate          = interpolate
    X%L2b%satellite            = sat

  END SUBROUTINE namelist_satellite

  !
  ! CTP_Tau
  !
  SUBROUTINE namelist_ctp_tau(x,file)

    CHARACTER(len=*),INTENT(in)    :: file
    TYPE(name_list), INTENT(inout) :: x
    REAL(wp), ALLOCATABLE          :: tbin_edges(:)
    REAL(wp), ALLOCATABLE          :: pbin_edges(:)
    INTEGER                        :: n_tbins, n_pbins
    INTEGER, PARAMETER             :: tmpbinsize = 100    ! for over-allocating ctp-tau bins
    REAL(wp), PARAMETER            :: tmpbinval = -3200.0 ! some unphysical value
    INTEGER :: i

    NAMELIST/ctp_tau_hist_vec/tbin_edges,pbin_edges

    OPEN(10,file=file,status='old')

    ! I don't know the size yet, so define a much larger array, I put
    ! initialise with unphysical values, and read the smaller array
    ! from the namelist into this array. I'll use a while loop to find
    ! the length of the bin vectors

    ALLOCATE ( tbin_edges(tmpbinsize) ,&
               pbin_edges(tmpbinsize) )
    tbin_edges(1:tmpbinsize) = tmpbinval
    pbin_edges(1:tmpbinsize) = tmpbinval
    READ(10,ctp_tau_hist_vec)
    CLOSE(10)

    ! ----------------
    ! CTP-TAU histograms
    ! Get the length of the arrays
    n_pbins = 1;
    DO WHILE (pbin_edges(n_pbins) .NE. tmpbinval)
       n_pbins = n_pbins+1
    END DO
    n_pbins  = n_pbins-2 ! -2 since length(pbin_edges)=n_pbins+1
    n_tbins = 1;
    DO WHILE (tbin_edges(n_tbins) .NE. tmpbinval)
       n_tbins = n_tbins+1
    END DO
    n_tbins  = n_tbins-2

    ALLOCATE (X%ctp_tau%tbin_edges (n_tbins+1) ,&
              X%ctp_tau%pbin_edges (n_pbins+1) )
    ALLOCATE (X%ctp_tau%tbin_centre(n_tbins  ) ,&
              X%ctp_tau%pbin_centre(n_pbins  ) )


    X%ctp_tau%tbin_edges(1:n_tbins+1) = tbin_edges(1:n_tbins+1)
    ! In the namelist the unit is [hPa], but we need it as [Pa] in the
    ! simulator
    X%ctp_tau%pbin_edges(1:n_pbins+1) = pbin_edges(1:n_pbins+1)*100

    DO i = 1, n_tbins
       X%ctp_tau%tbin_centre(i) = (tbin_edges(i)+tbin_edges(i+1))/2
    END DO
    DO i = 1, n_pbins
       X%ctp_tau%pbin_centre(i) = (pbin_edges(i)+pbin_edges(i+1))*50
    END DO

    X%ctp_tau%n_tbins      = n_tbins
    X%ctp_tau%n_pbins      = n_pbins
    X%ctp_tau%hist_phase   = (/0._wp,1._wp/)

    DEALLOCATE(tbin_edges,pbin_edges)

  END SUBROUTINE namelist_ctp_tau

  !
  ! Day night
  !
  SUBROUTINE namelist_daynight(x,file)

    CHARACTER(len=*),INTENT(in)    :: file
    TYPE(name_list), INTENT(inout) :: x
    REAL(wp)                       :: daylim,nightlim

    NAMELIST/daynight/daylim,nightlim

    OPEN(10,file=file,status='old')
    READ(10,daynight)
    CLOSE(10)

    X%daynight%daylim      = daylim
    X%daynight%nightlim    = nightlim
  END SUBROUTINE namelist_daynight

  !
  ! Microphysical
  !
  SUBROUTINE namelist_microphys(x,file)

    CHARACTER(len=*),INTENT(in)    :: file
    TYPE(name_list), INTENT(inout) :: x
    REAL(wp)                       :: tau_min
    INTEGER                        :: cf_method

    NAMELIST/cloudMicrophys/cf_method,tau_min

    OPEN(10,file=file,status='old')
    READ(10,cloudMicrophys)
    CLOSE(10)

    X%cloudMicrophys%cf_method          = cf_method
    X%cloudMicrophys%tau_min            = tau_min
  END SUBROUTINE namelist_microphys

END MODULE namelist_input
