MODULE CLOUD_CCI_M
  !
  ! Derived type for Cloud_cci
  !
  !
  ! Salomon.Eliasson@smhi.se

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
  USE OPTICS_M,            ONLY: NUM_TRIAL_RES

  IMPLICIT NONE

  PUBLIC :: ALLOCATE_CLOUD_CCI, &
       CLOUD_CCI_TYPE,          &
       DEALLOCATE_CLOUD_CCI,    &
       GET_NAMELIST_CLOUD_CCI,  &
       INITIALISE_CLOUD_CCI,    &
       INITIALISE_CLOUD_CCI_SIM

  PRIVATE :: ALLOCATE_CLOUD_CCI_SIM, &
       DEALLOCATE_CLOUD_CCI_SIM,     &
       VARIABLES_CLOUD_CCI

  TYPE cloud_cci_fields

     REAL(wp), ALLOCATABLE :: cer_ice   (:)
     REAL(wp), ALLOCATABLE :: cer_liq   (:)
     REAL(wp), ALLOCATABLE :: cfc       (:)
     REAL(wp), ALLOCATABLE :: cfc_low   (:)
     REAL(wp), ALLOCATABLE :: cfc_mid   (:)
     REAL(wp), ALLOCATABLE :: cfc_high  (:)
     REAL(wp), ALLOCATABLE :: cla_vis006(:)
     REAL(wp), ALLOCATABLE :: cot_ice   (:)
     REAL(wp), ALLOCATABLE :: cot_liq   (:)
     REAL(wp), ALLOCATABLE :: cth       (:)
     REAL(wp), ALLOCATABLE :: ctp       (:)
     REAL(wp), ALLOCATABLE :: ctp_log   (:)
     REAL(wp), ALLOCATABLE :: ctt       (:)
     REAL(wp), ALLOCATABLE :: cth_c     (:)
     REAL(wp), ALLOCATABLE :: ctp_c     (:)
     REAL(wp), ALLOCATABLE :: ctt_c     (:)
     REAL(wp), ALLOCATABLE :: iwp       (:)
     REAL(wp), ALLOCATABLE :: lwp       (:)
     REAL(wp), ALLOCATABLE :: tau       (:)

     REAL(wp), ALLOCATABLE :: hist2d_cot_ctp(:,:,:,:)

  END TYPE cloud_cci_fields

  TYPE cloud_cci_numel
     REAL(wp), ALLOCATABLE :: grd  (:)
     REAL(wp), ALLOCATABLE :: alb  (:)
     REAL(wp), ALLOCATABLE :: cld  (:)
     REAL(wp), ALLOCATABLE :: icld (:)
     REAL(wp), ALLOCATABLE :: lcld (:)
     REAL(wp), ALLOCATABLE :: day  (:)
  END TYPE cloud_cci_numel

  TYPE cloud_cci_type
     TYPE(cloud_cci_fields) :: av
     TYPE(cloud_cci_fields) :: sum
     TYPE(cloud_cci_numel)  :: numel
  END TYPE cloud_cci_type

CONTAINS

  SUBROUTINE ALLOCATE_CLOUD_CCI(cloud_cci,options,aux)

    IMPLICIT NONE
    TYPE(cloud_cci_type),INTENT(out)  :: cloud_cci
    TYPE(name_list),     INTENT(inout):: options
    TYPE(model_aux),     INTENT(in)   :: aux
    INTEGER                           :: ngrids,n_tbins,n_pbins
    LOGICAL                           :: need2Average

    need2Average = (.NOT.options%L2b%doL2bSampling .OR. options%L2b%node.EQ."all")

    ngrids  = aux%ngrids
    n_tbins = options%ctp_tau%n_tbins
    n_pbins = options%ctp_tau%n_pbins

    CALL ALLOCATE_CLOUD_CCI_SIM(cloud_cci%av,options,ngrids,n_pbins,n_tbins)

    IF (need2Average) THEN
       CALL ALLOCATE_CLOUD_CCI_SIM(cloud_cci%sum,options,ngrids,n_pbins,n_tbins)

       ALLOCATE (cloud_cci%numel%grd  (ngrids),&
                 cloud_cci%numel%alb  (ngrids),&
                 cloud_cci%numel%cld  (ngrids),&
                 cloud_cci%numel%icld (ngrids),&
                 cloud_cci%numel%lcld (ngrids),&
                 cloud_cci%numel%day  (ngrids))

    END IF

    ALLOCATE( &
         options%sim_aux%LUT%ice%optics%g0  (num_trial_res), &
         options%sim_aux%LUT%ice%optics%w0  (num_trial_res), &
         options%sim_aux%LUT%water%optics%g0(num_trial_res), &
         options%sim_aux%LUT%water%optics%w0(num_trial_res)  )

  END SUBROUTINE ALLOCATE_CLOUD_CCI

  SUBROUTINE ALLOCATE_CLOUD_CCI_SIM(IN,options,ngrids,n_pbins,n_tbins)

    IMPLICIT NONE
    TYPE(cloud_cci_fields), INTENT(inout):: IN
    TYPE(name_list), INTENT(in)          :: options
    TYPE(variablesContainer)             :: V
    INTEGER, INTENT(in)                  :: ngrids
    INTEGER, INTENT(in),OPTIONAL         :: n_tbins,n_pbins

    V = options%vars

    ALLOCATE(IN%cer_ice   (ngrids),&
             IN%cer_liq   (ngrids),&
             IN%cfc       (ngrids),&
             IN%cfc_low   (ngrids),&
             IN%cfc_mid   (ngrids),&
             IN%cfc_high  (ngrids),&
             IN%cla_vis006(ngrids),&
             IN%cot_ice   (ngrids),&
             IN%cot_liq   (ngrids),&
             IN%cth       (ngrids),&
             IN%ctp       (ngrids),&
             IN%ctp_log   (ngrids),&
             IN%ctt       (ngrids),&
             IN%cth_c     (ngrids),&
             IN%ctp_c     (ngrids),&
             IN%ctt_c     (ngrids),&
             IN%iwp       (ngrids),&
             IN%lwp       (ngrids),&
             IN%tau       (ngrids))

!    IF (V%hist2d_cot_ctp) &
         ALLOCATE(IN%hist2d_cot_ctp(ngrids,n_tbins,n_pbins,2))

    ! The albedo lookup tables are allocated in read_cloud_cci_albedo_LUT()
  END SUBROUTINE ALLOCATE_CLOUD_CCI_SIM

  SUBROUTINE INITIALISE_CLOUD_CCI(cloud_cci,ngrids,n_pbins,n_tbins)

    IMPLICIT NONE
    TYPE(cloud_cci_type), INTENT(inout)  :: cloud_cci
    INTEGER, INTENT(in) :: ngrids,n_tbins,n_pbins

    CALL INITIALISE_CLOUD_CCI_SIM(cloud_cci%av,ngrids,-9._wp,n_pbins,n_tbins)

    IF (ALLOCATED(cloud_cci%sum%ctp)) THEN
       CALL INITIALISE_CLOUD_CCI_SIM(cloud_cci%sum,ngrids,0._wp,n_pbins,n_tbins)
       cloud_cci%numel%grd  (1:ngrids) = 0._wp
       cloud_cci%numel%alb  (1:ngrids) = 0._wp
       cloud_cci%numel%cld  (1:ngrids) = 0._wp
       cloud_cci%numel%icld (1:ngrids) = 0._wp
       cloud_cci%numel%lcld (1:ngrids) = 0._wp
       cloud_cci%numel%day  (1:ngrids) = 0._wp
    END IF

  END SUBROUTINE INITIALISE_CLOUD_CCI

  SUBROUTINE INITIALISE_CLOUD_CCI_SIM(IN,ngrids,fill,n_pbins,n_tbins)

    IMPLICIT NONE
    TYPE(cloud_cci_fields), INTENT(inout) :: IN
    INTEGER, INTENT(in)          :: ngrids
    REAL(wp), INTENT(in)         :: fill
    INTEGER, INTENT(in),OPTIONAL :: n_pbins,n_tbins

    IN%cer_ice   (1:ngrids) = fill
    IN%cer_liq   (1:ngrids) = fill
    IN%cfc       (1:ngrids) = fill
    IN%cfc_low   (1:ngrids) = fill
    IN%cfc_mid   (1:ngrids) = fill
    IN%cfc_high  (1:ngrids) = fill
    IN%cla_vis006(1:ngrids) = fill
    IN%cot_ice   (1:ngrids) = fill
    IN%cot_liq   (1:ngrids) = fill
    IN%cth       (1:ngrids) = fill
    IN%ctp       (1:ngrids) = fill
    IN%ctp_log   (1:ngrids) = fill
    IN%ctt       (1:ngrids) = fill
    IN%cth_c     (1:ngrids) = fill
    IN%ctp_c     (1:ngrids) = fill
    IN%ctt_c     (1:ngrids) = fill
    IN%iwp       (1:ngrids) = fill
    IN%lwp       (1:ngrids) = fill
    IN%tau       (1:ngrids) = fill

    IF (ALLOCATED(IN%hist2d_cot_ctp)) &
         IN%hist2d_cot_ctp(1:ngrids,1:n_tbins,1:n_pbins,1:2) = 0._wp

    ! The albedo lookup tables are initialized in read_cloud_cci_albedo_LUT()
  END SUBROUTINE INITIALISE_CLOUD_CCI_SIM

  SUBROUTINE DEALLOCATE_CLOUD_CCI(cloud_cci)

    IMPLICIT NONE
    TYPE(cloud_cci_type), INTENT(inout)  :: cloud_cci

    CALL DEALLOCATE_CLOUD_CCI_SIM(cloud_cci%av)

    IF (ALLOCATED(cloud_cci%sum%ctp)) THEN
       CALL DEALLOCATE_CLOUD_CCI_SIM(cloud_cci%sum)
       DEALLOCATE (cloud_cci%numel%grd  ,&
                   cloud_cci%numel%alb  ,&
                   cloud_cci%numel%cld  ,&
                   cloud_cci%numel%icld ,&
                   cloud_cci%numel%lcld ,&
                   cloud_cci%numel%day   )
    END IF

  END SUBROUTINE DEALLOCATE_CLOUD_CCI

  SUBROUTINE DEALLOCATE_CLOUD_CCI_SIM(IN)

    IMPLICIT NONE
    TYPE(cloud_cci_fields), INTENT(inout) :: IN

    DEALLOCATE ( IN%cer_ice       ,&
                 IN%cer_liq       ,&
                 IN%cfc           ,&
                 IN%cfc_low       ,&
                 IN%cfc_mid       ,&
                 IN%cfc_high      ,&
                 IN%cla_vis006    ,&
                 IN%cot_ice       ,&
                 IN%cot_liq       ,&
                 IN%cth           ,&
                 IN%ctp           ,&
                 IN%ctp_log       ,&
                 IN%ctt           ,&
                 IN%cth_c         ,&
                 IN%ctp_c         ,&
                 IN%ctt_c         ,&
                 IN%hist2d_cot_ctp,&
                 IN%iwp           ,&
                 IN%lwp           ,&
                 IN%tau           )

  END SUBROUTINE DEALLOCATE_CLOUD_CCI_SIM

  SUBROUTINE get_namelist_cloud_cci(x,file)

    CHARACTER(len=*),INTENT(in)    :: file
    TYPE(name_list), INTENT(inout) :: x

    CALL common_namelist    (x,file)
    CALL namelist_ctp_tau   (x,file)
    CALL namelist_daynight  (x,file)
    CALL namelist_microphys (x,file)
    CALL variables_cloud_cci(x,file)

    X%sim%doClara     = .FALSE.
    X%sim%doCloud_cci = .TRUE.
    X%sim%doISCCP     = .FALSE.
    X%sim%doModel     = .FALSE.
    X%sim%doRTTOV     = .FALSE.

  END SUBROUTINE get_namelist_cloud_cci

  SUBROUTINE variables_cloud_cci(x,file)

    IMPLICIT NONE
    CHARACTER(len=*),INTENT(in)    :: file
    TYPE(name_list), INTENT(inout) :: x
    LOGICAL :: land_sea,solzen,time_of_day,cer_ice,cer_liq,&
         cfc,cfc_low,cfc_mid,cfc_high,cla_vis006,cot,cot_ice,cot_liq,&
         cth,ctp,ctp_log,ctt,cth_corrected,ctp_corrected,ctt_corrected,&
         hist2d_cot_ctp,iwp,lwp

    NAMELIST/variables2run/land_sea,solzen,time_of_day,cer_ice,cer_liq,&
         cfc,cfc_low,cfc_mid,cfc_high,cla_vis006,cot,cot_ice,cot_liq,&
         cth,ctp,ctp_log,ctt,cth_corrected,ctp_corrected,ctt_corrected,&
         hist2d_cot_ctp,iwp,lwp

    CALL initialise_variable_flag(land_sea,solzen,time_of_day,cer_ice,cer_liq,&
         cfc,cfc_low,cfc_mid,cfc_high,cla_vis006,cot,cot_ice,cot_liq,&
         cth,ctp,ctp_log,ctt,cth_corrected,ctp_corrected,ctt_corrected,&
         hist2d_cot_ctp,iwp,lwp)


    OPEN(10,file=file,status='old')
    READ(10,variables2run)
    CLOSE(10)

    x%vars%land_sea       = land_sea
    x%vars%solzen         = solzen
    x%vars%time_of_day    = time_of_day
    x%vars%cer_ice        = cer_ice
    x%vars%cer_liq        = cer_liq
    x%vars%cfc            = cfc
    x%vars%cfc_low        = cfc_low
    x%vars%cfc_mid        = cfc_mid
    x%vars%cfc_high       = cfc_high
    x%vars%cla_vis006     = cla_vis006
    x%vars%cot            = cot
    x%vars%cot_ice        = cot_ice
    x%vars%cot_liq        = cot_liq
    x%vars%cth            = cth
    x%vars%ctp            = ctp
    x%vars%ctp_log        = ctp_log
    x%vars%ctt            = ctt
    x%vars%cth_corrected  = cth_corrected
    x%vars%ctp_corrected  = ctp_corrected
    x%vars%ctt_corrected  = ctt_corrected
    x%vars%hist2d_cot_ctp = hist2d_cot_ctp
    x%vars%iwp            = iwp
    x%vars%lwp            = lwp

  END SUBROUTINE variables_cloud_cci

END MODULE CLOUD_CCI_M
