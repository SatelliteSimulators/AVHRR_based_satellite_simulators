MODULE CLAAS_M
  !
  ! Derived type for claas
  !
  !
  ! Salomon.Eliasson@smhi.se

  USE AUXILIARY_FUNCTIONS, ONLY: AUX_DATA
  USE COSP_KINDS,          ONLY: WP
  USE MODEL_INPUT,         ONLY: MODEL_AUX
  USE NAMELIST_INPUT,      ONLY: &
       COMMON_NAMELIST,          &
       INITIALISE_VARIABLE_FLAG, &
       NAMELIST_CTP_TAU,         &
       NAMELIST_DAYNIGHT,        &
       NAMELIST_MICROPHYS,       &
       NAME_LIST,                &
       VARIABLESCONTAINER
  USE OPTICS_M,            ONLY: &
       DEALLOCATE_OPTICS,        &
       NUM_TRIAL_RES,            &
       N_POD_EDGES,              &
       SIMULATOR_AUX

  IMPLICIT NONE

  PUBLIC :: &
       GET_NAMELIST_CLAAS,  &
       ALLOCATE_CLAAS,      &
       CLAAS_TYPE,          &
       DEALLOCATE_CLAAS,    &
       INITIALISE_CLAAS,    &
       INITIALISE_CLAAS_SIM

  PRIVATE :: ALLOCATE_CLAAS_SIM, &
       DEALLOCATE_CLAAS_SIM,     &
       VARIABLES_CLAAS

  TYPE claas_fields

     REAL(wp), ALLOCATABLE :: cfc           (:)
     REAL(wp), ALLOCATABLE :: cfc_day       (:)
     REAL(wp), ALLOCATABLE :: cfc_low       (:)
     REAL(wp), ALLOCATABLE :: cfc_mid       (:)
     REAL(wp), ALLOCATABLE :: cfc_high      (:)
     INTEGER,  ALLOCATABLE :: cflag_tot     (:,:)
     REAL(wp), ALLOCATABLE :: cot_ice       (:)
     REAL(wp), ALLOCATABLE :: cot_liq       (:)
     REAL(wp), ALLOCATABLE :: cth           (:)
     REAL(wp), ALLOCATABLE :: ctp_log       (:)
     REAL(wp), ALLOCATABLE :: ctp           (:)
     REAL(wp), ALLOCATABLE :: ctt           (:)
     REAL(wp), ALLOCATABLE :: iwp           (:)
     REAL(wp), ALLOCATABLE :: hist2d_cot_ctp(:,:,:,:)
     REAL(wp), ALLOCATABLE :: lwp           (:)
     REAL(wp), ALLOCATABLE :: ref_ice       (:)
     REAL(wp), ALLOCATABLE :: ref_liq       (:)
     REAL(wp), ALLOCATABLE :: tau           (:)

  END TYPE claas_fields

  TYPE claas_numel
     REAL(wp), ALLOCATABLE :: grd  (:)
     REAL(wp), ALLOCATABLE :: cld  (:)
     REAL(wp), ALLOCATABLE :: icld (:)
     REAL(wp), ALLOCATABLE :: lcld (:)
     REAL(wp), ALLOCATABLE :: day  (:)
  END TYPE claas_numel

  TYPE claas_type
     TYPE(claas_fields) :: av
     TYPE(claas_fields) :: sum
     TYPE(claas_numel)  :: numel
  END TYPE claas_type

CONTAINS

    SUBROUTINE GET_NAMELIST_CLAAS(x,file)

    IMPLICIT NONE

    CHARACTER(len=*),INTENT(in)    :: file
    TYPE(name_list), INTENT(inout) :: x

    CALL common_namelist   (x,file)
    CALL namelist_ctp_tau  (x,file)
    CALL namelist_daynight (x,file)
    CALL namelist_microphys(x,file)
    CALL variables_claas   (x,file)

    X%sim%doClaas     = .TRUE.
    X%sim%doClara     = .FALSE.
    X%sim%doCloud_cci = .FALSE.
    X%sim%doISCCP     = .FALSE.
    X%sim%doModel     = .FALSE.
    X%sim%doRTTOV     = .FALSE.

  END SUBROUTINE GET_NAMELIST_CLAAS

  SUBROUTINE ALLOCATE_CLAAS(claas,options,aux)

    IMPLICIT NONE

    TYPE(claas_type), INTENT(out) :: claas
    TYPE(name_list), INTENT(inout):: options
    TYPE(model_aux),  INTENT(in)  :: aux
    INTEGER                       :: ngrids,n_tbins,n_pbins,ncol
    LOGICAL                       :: need2Average

    need2Average = .TRUE.

    ngrids  = aux%ngrids
    n_tbins = options%ctp_tau%n_tbins
    n_pbins = options%ctp_tau%n_pbins
    ncol    = options%ncols

    CALL ALLOCATE_CLAAS_SIM(claas%av,options,ngrids,n_pbins,n_tbins)
    CALL ALLOCATE_CLAAS_SIM(claas%sum,options,ngrids,n_pbins,n_tbins)

    ALLOCATE (claas%numel%grd  (ngrids),&
          claas%numel%cld  (ngrids),&
          claas%numel%icld (ngrids),&
          claas%numel%lcld (ngrids),&
          claas%numel%day  (ngrids))

    CALL get_POD_FAR(options, aux)

    ALLOCATE( &
         options%sim_aux%LUT%ice%optics%g0  (num_trial_res), &
         options%sim_aux%LUT%ice%optics%w0  (num_trial_res), &
         options%sim_aux%LUT%water%optics%g0(num_trial_res), &
         options%sim_aux%LUT%water%optics%w0(num_trial_res)  )

    ! This has to be done earlier (outside dayloop)
    options%sim_aux%LUT%ice%optics%g0  (1:num_trial_res) = 0._wp
    options%sim_aux%LUT%ice%optics%w0  (1:num_trial_res) = 0._wp
    options%sim_aux%LUT%water%optics%g0(1:num_trial_res) = 0._wp
    options%sim_aux%LUT%water%optics%w0(1:num_trial_res) = 0._wp

  END SUBROUTINE ALLOCATE_claas

  SUBROUTINE ALLOCATE_claas_SIM(IN,options,ngrids,n_pbins,n_tbins)

    IMPLICIT NONE
    TYPE(claas_fields), INTENT(inout) :: IN
    TYPE(name_list), INTENT(in)       :: options
    TYPE(variablesContainer)          :: V
    INTEGER, INTENT(in)               :: ngrids
    INTEGER, INTENT(in),OPTIONAL      :: n_tbins,n_pbins

    V = options%vars

    ALLOCATE(IN%cfc     (ngrids),&
             IN%cfc_day (ngrids),&
             IN%cfc_low (ngrids),&
             IN%cfc_mid (ngrids),&
             IN%cfc_high(ngrids),&
             IN%cot_ice (ngrids),&
             IN%cot_liq (ngrids),&
             IN%cth     (ngrids),&
             IN%ctp_log (ngrids),&
             IN%ctp     (ngrids),&
             IN%ctt     (ngrids),&
             IN%iwp     (ngrids),&
             IN%lwp     (ngrids),&
             IN%ref_ice (ngrids),&
             IN%ref_liq (ngrids),&
             IN%tau     (ngrids))
     ALLOCATE(IN%hist2d_cot_ctp(ngrids,n_tbins,n_pbins,2))
     ALLOCATE(IN%cflag_tot(ngrids,5))

  END SUBROUTINE ALLOCATE_claas_SIM

  SUBROUTINE INITIALISE_claas(claas,options,ngrids)

    IMPLICIT NONE
    TYPE(claas_type), INTENT(inout)     :: claas
    TYPE(name_list), INTENT(inout)      :: options
    INTEGER, INTENT(in)                 :: ngrids

    INTEGER :: n_tbins,n_pbins

    n_tbins = options%ctp_tau%n_tbins
    n_pbins = options%ctp_tau%n_pbins

    CALL INITIALISE_claas_SIM(claas%av,ngrids,-9._wp,n_pbins,n_tbins)

    IF (ALLOCATED(claas%sum%ctp)) THEN
       CALL INITIALISE_claas_SIM(claas%sum,ngrids,0._wp,n_pbins,n_tbins)
       claas%numel%grd  (1:ngrids) = 0._wp
       claas%numel%cld  (1:ngrids) = 0._wp
       claas%numel%icld (1:ngrids) = 0._wp
       claas%numel%lcld (1:ngrids) = 0._wp
       claas%numel%day  (1:ngrids) = 0._wp
    END IF

  END SUBROUTINE INITIALISE_claas

  SUBROUTINE INITIALISE_claas_SIM(IN,ngrids,fill,n_pbins,n_tbins)

    IMPLICIT NONE
    TYPE(claas_fields), INTENT(inout) :: IN
    INTEGER, INTENT(in)          :: ngrids
    REAL(wp), INTENT(in)         :: fill
    INTEGER, INTENT(in),OPTIONAL :: n_pbins,n_tbins

    IN%cfc      (1:ngrids) = fill
    IN%cfc_day  (1:ngrids) = fill
    IN%cfc_low  (1:ngrids) = fill
    IN%cfc_mid  (1:ngrids) = fill
    IN%cfc_high (1:ngrids) = fill
    IN%cot_ice  (1:ngrids) = fill
    IN%cot_liq  (1:ngrids) = fill
    IN%cth      (1:ngrids) = fill
    IN%ctp_log  (1:ngrids) = fill
    IN%ctp      (1:ngrids) = fill
    IN%ctt      (1:ngrids) = fill
    IN%iwp      (1:ngrids) = fill
    IN%lwp      (1:ngrids) = fill
    IN%ref_ice  (1:ngrids) = fill
    IN%ref_liq  (1:ngrids) = fill
    IN%tau      (1:ngrids) = fill
    IN%cflag_tot(1:ngrids,1:5) = 0
    IN%hist2d_cot_ctp(1:ngrids,1:n_tbins,1:n_pbins,1:2) = 0._wp

  END SUBROUTINE INITIALISE_claas_SIM

  SUBROUTINE DEALLOCATE_claas(claas,options)

    IMPLICIT NONE
    TYPE(claas_type), INTENT(inout)  :: claas
    TYPE(name_list),  INTENT(inout):: options

    CALL DEALLOCATE_claas_SIM(claas%av)

    IF (ALLOCATED(claas%sum%ctp)) THEN
       CALL DEALLOCATE_claas_SIM(claas%sum)
       DEALLOCATE (claas%numel%grd ,&
            claas%numel%cld ,&
            claas%numel%icld,&
            claas%numel%lcld,&
            claas%numel%day)
    END IF


    DEALLOCATE( &
         options%sim_aux%LUT%ice%optics%g0, &
         options%sim_aux%LUT%ice%optics%w0, &
         options%sim_aux%LUT%water%optics%g0, &
         options%sim_aux%LUT%water%optics%w0  )

  END SUBROUTINE DEALLOCATE_claas

  SUBROUTINE DEALLOCATE_claas_SIM(IN)

    IMPLICIT NONE
    TYPE(claas_fields), INTENT(inout) :: IN

    DEALLOCATE(IN%cfc           ,&
               IN%cfc_day       ,&
               IN%cfc_low       ,&
               IN%cfc_mid       ,&
               IN%cfc_high      ,&
               IN%cot_ice       ,&
               IN%cot_liq       ,&
               IN%cth           ,&
               IN%ctp_log       ,&
               IN%ctp           ,&
               IN%ctt           ,&
               IN%iwp           ,&
               IN%lwp           ,&
               IN%ref_ice       ,&
               IN%ref_liq       ,&
               IN%tau           )
    DEALLOCATE(IN%cflag_tot)
    DEALLOCATE(IN%hist2d_cot_ctp)
  END SUBROUTINE DEALLOCATE_claas_SIM

  SUBROUTINE VARIABLES_claas(x,file)

    IMPLICIT NONE
    CHARACTER(len=*),INTENT(in)    :: file
    TYPE(name_list), INTENT(inout) :: x
    LOGICAL :: land_sea,solzen,time_of_day,cfc,cfc_day,cfc_low,&
         cfc_mid,cfc_high,cot_ice,cflag_tot,cot_liq,cth,ctp,ctp_log,ctt,iwp,&
         lwp,hist2d_cot_ctp,ref_ice,ref_liq

    NAMELIST/variables2run/land_sea,solzen,time_of_day,cfc,cfc_day,cfc_low,&
         &cfc_mid,cfc_high,cot_ice,cflag_tot,cot_liq,cth,ctp,ctp_log,ctt,iwp,&
         &lwp,hist2d_cot_ctp,ref_ice,ref_liq

    CALL initialise_variable_flag(land_sea,solzen,time_of_day,&
         cfc,cfc_day,cfc_low,cfc_mid,cfc_high,cot_ice,cflag_tot,&
         cot_liq,cth,ctp,ctp_log,ctt,iwp,lwp,hist2d_cot_ctp,ref_ice,ref_liq)

    OPEN(10,file=file,status='old')
    READ(10,variables2run)
    CLOSE(10)

    x%vars%land_sea       = land_sea
    x%vars%solzen         = solzen
    x%vars%time_of_day    = time_of_day
    x%vars%cfc            = cfc
    x%vars%cfc_day        = cfc_day
    x%vars%cfc_low        = cfc_low
    x%vars%cfc_mid        = cfc_mid
    x%vars%cfc_high       = cfc_high
    x%vars%cot_ice        = cot_ice
    x%vars%cflag_tot      = cflag_tot
    x%vars%cot_liq        = cot_liq
    x%vars%cth            = cth
    x%vars%ctp            = ctp
    x%vars%ctp_log        = ctp_log
    x%vars%ctt            = ctt
    x%vars%iwp            = iwp
    x%vars%lwp            = lwp
    x%vars%hist2d_cot_ctp = hist2d_cot_ctp
    x%vars%ref_ice        = ref_ice
    x%vars%ref_liq        = ref_liq

  END SUBROUTINE VARIABLES_claas

  SUBROUTINE get_POD_FAR(options, aux)

     IMPLICIT NONE

     TYPE(name_list),  INTENT(inout):: options
     TYPE(model_aux),  INTENT(in)  :: aux
     REAL(wp), ALLOCATABLE         :: tmp3(:,:,:)
     REAL(wp), ALLOCATABLE         :: tmp2(:,:)
     INTEGER                       :: len
     CHARACTER(1000)               :: filename
     CHARACTER(100)                :: varstr
     CHARACTER(10)                 :: season
     INTEGER                       :: ngrids,ncol

     ngrids  = aux%ngrids
     ncol    = options%ncols
     WRITE(filename, '(A,A,A)') &
          'data/microphysics/CLAAS/cloud_properties/cloud_mask_limits_',&
          TRIM(options%CDR),'.nc'

     IF (options%epoch%month.GE.4 .AND. options%epoch%month.LE.9) THEN
          season='NHsummer'
     ELSE
          season='NHwinter'
     END IF

     !--------------------------------------
     ! POD
     !
     ! claas-3: Based on half-yearly statistics
     !           separated for day and night based on matchups
     !           between noaa18 and noaa19 and CALIPSO (5km collocated)
     CALL AUX_DATA(filename,'COT_edges',strDim1='COT_edges',&
          data1D=options%sim_aux%POD_tau_bin_edges,dbg=options%dbg)

     ! Making the last box also valid for all clouds greater than tau=5
     options%sim_aux%POD_tau_bin_edges(SIZE(options%sim_aux%POD_tau_bin_edges))=9999._wp

     CALL AUX_DATA(filename,'COT_centers',strDim1='COT_centers',&
            data1D=options%sim_aux%POD_tau_bin_centers,nX=len,dbg=options%dbg)

     ALLOCATE (options%sim_aux%POD_layers (ngrids,len,2))

     ! Årstiden summer är här definierat som april till september
     ! och vinter oktober till mars.

     WRITE (*,'(5(A,x))') &
            " --- Getting POD values based on",TRIM(season),"season for ",TRIM(options%CDR),"CDR"

     ! NIGHT values
     WRITE(varstr,'(A,A)') 'POD_NIGHT_',season
     CALL AUX_DATA(filename,TRIM(varstr),'lon','lat',strDim3='COT_centers',&
                    data3D=tmp3,conform2model=.TRUE.,aux=aux,dbg=options%dbg)

     options%sim_aux%POD_layers(:,:,1) = RESHAPE(tmp3,(/ngrids,len/))

     ! Day values
     WRITE(varstr,'(A,A)') 'POD_DAY_',season
     CALL AUX_DATA(filename,TRIM(varstr),'lon','lat',strDim3='COT_centers',&
            data3D=tmp3,conform2model=.TRUE.,aux=aux,dbg=options%dbg)

     options%sim_aux%POD_layers(:,:,2) = RESHAPE(tmp3,(/ngrids,len/))

     DEALLOCATE(tmp3)

     !---------------------------------------
     ! FAR
     !
     ALLOCATE (options%sim_aux%FAR (ngrids,2))

     ! Årstiden summer är här definierat som april till september
     ! och vinter oktober till mars.

     WRITE (*,'(5(A,x))') &
            " --- Getting FAR values based on",TRIM(season),"season for ",TRIM(options%CDR),"CDR"

     ! NIGHT values
     WRITE(varstr,'(A,A)') 'false_alarm_rate_NIGHT_',season
     CALL AUX_DATA(filename,TRIM(varstr),'lon','lat',&
                    data2D=tmp2,conform2model=.TRUE.,aux=aux,dbg=options%dbg)

     options%sim_aux%FAR(:,1) = RESHAPE(tmp2,(/ngrids/))

     ! Daytime values
     WRITE(varstr,'(A,A)') 'false_alarm_rate_DAY_',season
     CALL AUX_DATA(filename,TRIM(varstr),'lon','lat',&
            data2D=tmp2,conform2model=.TRUE.,aux=aux,dbg=options%dbg)

     options%sim_aux%FAR(:,2) = RESHAPE(tmp2,(/ngrids/))

     DEALLOCATE(tmp2)

     !---------------------------------------
     ! Random numbers
     !
     ALLOCATE (options%sim_aux%random_numbers(2,ngrids,ncol))
     CALL RANDOM_NUMBER(harvest = options%sim_aux%random_numbers)

END SUBROUTINE get_POD_FAR

END MODULE claas_M
