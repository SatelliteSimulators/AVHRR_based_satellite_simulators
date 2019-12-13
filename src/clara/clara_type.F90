MODULE CLARA_M
  !
  ! Derived type for clara
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
       SIM_AUX

  IMPLICIT NONE

  PUBLIC :: ALLOCATE_CLARA, &
       CLARA_TYPE,          &
       DEALLOCATE_CLARA,    &
       GET_NAMELIST_CLARA,  &
       INITIALISE_CLARA,    &
       INITIALISE_CLARA_SIM

  PRIVATE :: ALLOCATE_CLARA_SIM, &
       DEALLOCATE_CLARA_SIM,     &
       VARIABLES_CLARA


  TYPE clara_fields

     REAL(wp), ALLOCATABLE :: cfc           (:)
     REAL(wp), ALLOCATABLE :: cfc_day       (:)
     REAL(wp), ALLOCATABLE :: cfc_low       (:)
     REAL(wp), ALLOCATABLE :: cfc_mid       (:)
     REAL(wp), ALLOCATABLE :: cfc_high      (:)
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

  END TYPE clara_fields

  TYPE clara_numel
     REAL(wp), ALLOCATABLE :: grd  (:)
     REAL(wp), ALLOCATABLE :: cld  (:)
     REAL(wp), ALLOCATABLE :: icld (:)
     REAL(wp), ALLOCATABLE :: lcld (:)
     REAL(wp), ALLOCATABLE :: day  (:)
  END TYPE clara_numel

  TYPE clara_type
     TYPE(clara_fields) :: av
     TYPE(clara_fields) :: sum
     TYPE(clara_numel)  :: numel
  END TYPE clara_type

CONTAINS

  SUBROUTINE ALLOCATE_CLARA(clara,options,aux)

    IMPLICIT NONE

    TYPE(clara_type), INTENT(out) :: clara
    TYPE(name_list), INTENT(inout):: options
    TYPE(model_aux),  INTENT(in)  :: aux
    INTEGER                       :: ngrids,n_tbins,n_pbins,ncol
    
    REAL(wp), ALLOCATABLE         :: tmp(:,:)
    REAL(wp), ALLOCATABLE         :: tmp3(:,:,:)
    INTEGER                       :: len
    CHARACTER(1000)               :: filename
    LOGICAL                       :: need2Average

    need2Average = (.NOT.options%L2b%doL2bSampling .OR. options%L2b%node.EQ."all")

    ngrids  = aux%ngrids 
    n_tbins = options%ctp_tau%n_tbins
    n_pbins = options%ctp_tau%n_pbins
    ncol    = options%ncols

    WRITE(filename, '(A,A)') TRIM(options%paths%data_dir),&
         'microphysics/CLARA/cloud_properties/cloud_mask_limits.nc'

    CALL ALLOCATE_CLARA_SIM(clara%av,options,ngrids,n_pbins,n_tbins)

    IF (need2Average) THEN
       CALL ALLOCATE_CLARA_SIM(clara%sum,options,ngrids,n_pbins,n_tbins)

       ALLOCATE (clara%numel%grd  (ngrids),&
            clara%numel%cld  (ngrids),&
            clara%numel%icld (ngrids),&
            clara%numel%lcld (ngrids),&
            clara%numel%day  (ngrids))

    END IF

    SELECT CASE (options%cloudMicrophys%cf_method)

    CASE(1)
       ALLOCATE (options%sim_aux%detection_limit(ngrids,2))
       ! use a map of optical depths

       ! Night values
       CALL AUX_DATA(filename,'detection_limit_night',&
            'lon','lat',tmp,.TRUE.,aux,dbg=options%dbg)
       
       options%sim_aux%detection_limit(:,1)=RESHAPE(tmp,(/ngrids/))

       ! Day values
       CALL AUX_DATA(filename,'detection_limit_day',&
            'lon','lat',tmp,.TRUE.,aux,dbg=options%dbg)
       
       options%sim_aux%detection_limit(:,2)=RESHAPE(tmp,(/ngrids/))

       DEALLOCATE(tmp)
    CASE(2)

       CALL AUX_DATA(filename,'COT_edges',strDim1='COT_edges',&
            data1D=options%sim_aux%POD_tau_bin_edges,dbg=options%dbg)

       ! Making the last box also valid for all clouds greater than tau=5
       options%sim_aux%POD_tau_bin_edges(SIZE(options%sim_aux%POD_tau_bin_edges))=9999._wp

       CALL AUX_DATA(filename,'COT_centers',strDim1='COT_centers',&
            data1D=options%sim_aux%POD_tau_bin_centers,nX=len,dbg=options%dbg)

       ALLOCATE (options%sim_aux%POD_layers (ngrids,len,2))

       ! NIGHT values
       CALL AUX_DATA(filename,'POD_night','lon','lat',strDim3='COT_centers',&
            data3D=tmp3,conform2model=.TRUE.,aux=aux,dbg=options%dbg)

       options%sim_aux%POD_layers(:,:,1) = RESHAPE(tmp3,(/ngrids,len/))

       ! Day values
       CALL AUX_DATA(filename,'POD_day','lon','lat',strDim3='COT_centers',&
            data3D=tmp3,conform2model=.TRUE.,aux=aux,dbg=options%dbg)

       options%sim_aux%POD_layers(:,:,2) = RESHAPE(tmp3,(/ngrids,len/))

       DEALLOCATE(tmp3)

       ! Get the random numbers I need (sunlit,:,:)
       ALLOCATE (options%sim_aux%random_numbers(2,ngrids,ncol))
       CALL RANDOM_NUMBER(harvest = options%sim_aux%random_numbers)

    CASE (3)
       ! same as case 2 yet for combined night and day

       CALL AUX_DATA(filename,'COT_edges',strDim1='COT_edges',&
            data1D=options%sim_aux%POD_tau_bin_edges,dbg=options%dbg)

       ! Making the last box also valid for all clouds greater than tau=5
       options%sim_aux%POD_tau_bin_edges(SIZE(options%sim_aux%POD_tau_bin_edges))=9999._wp

       CALL AUX_DATA(filename,'COT_centers',strDim1='COT_centers',&
            data1D=options%sim_aux%POD_tau_bin_centers,nX=len,dbg=options%dbg)

       ALLOCATE (options%sim_aux%POD_layers (ngrids,len,1))

       ! Combined values
       CALL AUX_DATA(filename,'POD_combined','lon','lat',strDim3='COT_centers',&
            data3D=tmp3,conform2model=.TRUE.,aux=aux,dbg=options%dbg)

       options%sim_aux%POD_layers(:,:,1) = RESHAPE(tmp3,(/ngrids,len/))

       DEALLOCATE(tmp3)

       ALLOCATE (options%sim_aux%random_numbers(1,ngrids,ncol))
       CALL RANDOM_NUMBER(harvest = options%sim_aux%random_numbers)

    END SELECT

    ALLOCATE( &
         options%sim_aux%LUT%ice%optics%g0  (num_trial_res), &
         options%sim_aux%LUT%ice%optics%w0  (num_trial_res), &
         options%sim_aux%LUT%water%optics%g0(num_trial_res), &
         options%sim_aux%LUT%water%optics%w0(num_trial_res)  )

  END SUBROUTINE ALLOCATE_CLARA

  SUBROUTINE ALLOCATE_CLARA_SIM(IN,options,ngrids,n_pbins,n_tbins)

    IMPLICIT NONE
    TYPE(clara_fields), INTENT(inout) :: IN
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

  END SUBROUTINE ALLOCATE_CLARA_SIM

  SUBROUTINE INITIALISE_CLARA(clara,ngrids,n_pbins,n_tbins)

    IMPLICIT NONE
    TYPE(clara_type), INTENT(inout)     :: clara
    INTEGER, INTENT(in) :: ngrids,n_tbins,n_pbins

    CALL INITIALISE_CLARA_SIM(clara%av,ngrids,-999._wp,n_pbins,n_tbins)

    IF (ALLOCATED(clara%sum%ctp)) THEN
       CALL INITIALISE_CLARA_SIM(clara%sum,ngrids,0._wp,n_pbins,n_tbins)
       clara%numel%grd  (1:ngrids) = 0._wp
       clara%numel%cld  (1:ngrids) = 0._wp
       clara%numel%icld (1:ngrids) = 0._wp
       clara%numel%lcld (1:ngrids) = 0._wp
       clara%numel%day  (1:ngrids) = 0._wp
    END IF

  END SUBROUTINE INITIALISE_CLARA

  SUBROUTINE INITIALISE_CLARA_SIM(IN,ngrids,fill,n_pbins,n_tbins)

    IMPLICIT NONE
    TYPE(clara_fields), INTENT(inout) :: IN
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

    IF (ALLOCATED(IN%hist2d_cot_ctp)) &
         IN%hist2d_cot_ctp(1:ngrids,1:n_tbins,1:n_pbins,1:2) = 0._wp

  END SUBROUTINE INITIALISE_CLARA_SIM

  SUBROUTINE DEALLOCATE_CLARA(clara)

    IMPLICIT NONE
    TYPE(clara_type), INTENT(inout)  :: clara

    CALL DEALLOCATE_CLARA_SIM(clara%av)

    IF (ALLOCATED(clara%sum%ctp)) THEN
       CALL DEALLOCATE_CLARA_SIM(clara%sum)
       DEALLOCATE (clara%numel%grd ,&
            clara%numel%cld ,&
            clara%numel%icld,&
            clara%numel%lcld,&
            clara%numel%day)
    END IF

  END SUBROUTINE DEALLOCATE_CLARA

  SUBROUTINE DEALLOCATE_CLARA_SIM(IN)

    IMPLICIT NONE
    TYPE(clara_fields), INTENT(inout) :: IN

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

    IF (ALLOCATED(IN%hist2d_cot_ctp)) DEALLOCATE(IN%hist2d_cot_ctp)
  END SUBROUTINE DEALLOCATE_CLARA_SIM

  SUBROUTINE GET_NAMELIST_CLARA(x,file)

    IMPLICIT NONE
    
    CHARACTER(len=*),INTENT(in)    :: file
    TYPE(name_list), INTENT(inout) :: x
    
    CALL common_namelist   (x,file)
    CALL namelist_ctp_tau  (x,file)
    CALL namelist_daynight (x,file)
    CALL namelist_microphys(x,file)
    CALL variables_clara   (x,file)

    X%sim%doClara     = .TRUE.
    X%sim%doCloud_cci = .FALSE.
    X%sim%doISCCP     = .FALSE.
    X%sim%doModel     = .FALSE.
    X%sim%doRTTOV     = .FALSE.

  END SUBROUTINE GET_NAMELIST_CLARA
  
  SUBROUTINE VARIABLES_CLARA(x,file)

    IMPLICIT NONE
    CHARACTER(len=*),INTENT(in)    :: file
    TYPE(name_list), INTENT(inout) :: x
    LOGICAL :: land_sea,solzen,time_of_day,cfc,cfc_day,cfc_low,&
         cfc_mid,cfc_high,cot_ice,cot_liq,cth,ctp,ctp_log,ctt,iwp,&
         lwp,hist2d_cot_ctp,ref_ice,ref_liq
    
    NAMELIST/variables2run/land_sea,solzen,time_of_day,cfc,cfc_day,cfc_low,&
         &cfc_mid,cfc_high,cot_ice,cot_liq,cth,ctp,ctp_log,ctt,iwp,&
         &lwp,hist2d_cot_ctp,ref_ice,ref_liq
    
    CALL initialise_variable_flag(land_sea,solzen,time_of_day,&
         cfc,cfc_day,cfc_low,cfc_mid,cfc_high,cot_ice,&
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
    
  END SUBROUTINE VARIABLES_CLARA
END MODULE CLARA_M
