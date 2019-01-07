MODULE INTERNAL_SIMULATOR
  ! Structure types for the simulated products
  !
  ! Salomon.Eliasson@smhi.se

  USE NAMELIST_INPUT, ONLY: NAME_LIST
  USE COSP_KINDS,     ONLY: WP

  IMPLICIT NONE

  PUBLIC :: ALLOCATE_INTERNAL_SIMULATOR,&
       INITIALISE_INTERNAL_SIMULATOR,&
       DEALLOCATE_INTERNAL_SIMULATOR


  ! ----------------
  ! INTERNAL
  ! ------------

  ! ----------------
  ! VARIABLES computed based on model output vars
  ! 
  ! 'cflag'       = cloud flag per column (cfree, subvisible, semi-transparent, or opaque)
  ! 'cfrac'       = cloud fraction at layer. Layer n is the layer between layers n and n+1
  ! 'cloud'       = RTTOV cloud type (6 types). contains LWC and IWC
  ! 'cph'         = cloud-top phase for each sub grid (1=liq,2=ice)
  ! 'cth'         = cloud top pressure per subgrid
  ! 'ctp'         = cloud top pressure per subgrid
  ! 'ctt'         = cloud top temperature per subgrid
  ! 'cth_c'       = corrected cloud top pressure per subgrid
  ! 'ctp_c'       = corrected cloud top pressure per subgrid
  ! 'ctt_c'       = corrected cloud top temperature per subgrid
  ! 'reff'        = effective radius (either ice or liquid)
  ! 'tau'         = simulated ice+liquid cloud optical depth for each sub gri
  ! 'tau_profile' = vertical profile of optical depth with cloud
  !                 fraction from scops.f applied
  ! 'Tb'          = Brightness temperature from RTTOV
  ! 'Tb_cld_CDK'  = Tb using the CDK approach

  TYPE internal_clara
     INTEGER, ALLOCATABLE, DIMENSION(:)     :: flagged
  END TYPE internal_clara

  TYPE internal_cloud_cci
     LOGICAL,  ALLOCATABLE, DIMENSION(:)    :: albedoIsDefined
     REAL(wp), ALLOCATABLE, DIMENSION(:)    :: cla_vis006
     REAL(wp), ALLOCATABLE, DIMENSION(:)    :: cth_c
     REAL(wp), ALLOCATABLE, DIMENSION(:)    :: ctp_c
     REAL(wp), ALLOCATABLE, DIMENSION(:)    :: ctt_c
  END TYPE internal_cloud_cci

  TYPE internal
     INTEGER,  ALLOCATABLE, DIMENSION(:)    :: cflag
     REAL(wp), ALLOCATABLE, DIMENSION(:,:)  :: cfrac
     REAL(wp), ALLOCATABLE, DIMENSION(:,:,:):: cloud
     INTEGER,  ALLOCATABLE, DIMENSION(:)    :: cph
     REAL(wp), ALLOCATABLE, DIMENSION(:)    :: cth
     REAL(wp), ALLOCATABLE, DIMENSION(:)    :: ctp
     REAL(wp), ALLOCATABLE, DIMENSION(:)    :: ctt
     REAL(wp), ALLOCATABLE, DIMENSION(:)    :: cwp
     REAL(wp), ALLOCATABLE, DIMENSION(:,:)  :: frac_in_up_cld
     REAL(wp), ALLOCATABLE, DIMENSION(:)    :: reff
     REAL(wp), ALLOCATABLE, DIMENSION(:)    :: tau
     REAL(wp), ALLOCATABLE, DIMENSION(:,:)  :: tau_profile
     REAL(wp), ALLOCATABLE, DIMENSION(:,:)  :: Tb
     REAL(wp), ALLOCATABLE, DIMENSION(:)    :: Tb_cld_CDK
     TYPE(internal_clara)                   :: clara
     TYPE(internal_cloud_cci)               :: cloud_cci
  END TYPE internal

CONTAINS

  SUBROUTINE ALLOCATE_INTERNAL_SIMULATOR(inter,ncol,nlev,options)

    IMPLICIT NONE

    ! Allocate local arrays

    TYPE(internal), INTENT(inout) :: inter
    TYPE(name_list), INTENT(in) :: options

    INTEGER, INTENT(in) :: ncol,nlev
    INTEGER :: nchan

    nchan = options%sim%nchannels

    ALLOCATE( inter%cflag          (ncol       ),&
              inter%cph            (ncol       ),&
              inter%cth            (ncol       ),&
              inter%ctp            (ncol       ),&
              inter%ctt            (ncol       ),&
              inter%cwp            (ncol       ),&
              inter%frac_in_up_cld (ncol, nlev ),&
              inter%reff           (ncol       ),&
              inter%tau            (ncol       ),&
              inter%tau_profile    (ncol, nlev ))

    IF (options%sim%doClara) THEN
       ALLOCATE(inter%clara%flagged(ncol),&
            inter%Tb_cld_CDK       (ncol),&
            inter%Tb               (ncol,  nchan  ))
       IF (options%sim%Tb.EQ.1) THEN
          ALLOCATE(inter%cfrac     (ncol,  nlev-1),&
               inter%cloud         (ncol,6,nlev-1))
       END IF
    ELSEIF (options%sim%doCloud_cci) THEN
       ALLOCATE( inter%cloud_cci%albedoIsDefined (ncol),&
                 inter%cloud_cci%cla_vis006      (ncol),&
                 inter%cloud_cci%cth_c           (ncol),&
                 inter%cloud_cci%ctp_c           (ncol),&
                 inter%cloud_cci%ctt_c           (ncol))
    ELSEIF (options%sim%doRTTOV) THEN
       ALLOCATE( inter%cfrac         (ncol,  nlev-1),&
                 inter%cloud         (ncol,6,nlev-1),&
                 inter%Tb            (ncol,  nchan  ))
    END IF
  END SUBROUTINE ALLOCATE_INTERNAL_SIMULATOR
  SUBROUTINE INITIALISE_INTERNAL_SIMULATOR(inter,ncol,nlev,options)

    IMPLICIT NONE

    INTEGER,INTENT(in) :: ncol,nlev
    TYPE(internal), INTENT(inout) :: inter
    TYPE(name_list), INTENT(in)   :: options
    INTEGER :: nchan

    nchan = options%sim%nchannels

    inter%cflag         (1:ncol       ) = -999
    inter%cph           (1:ncol       ) = -999
    inter%cth           (1:ncol       ) = -999._wp
    inter%ctp           (1:ncol       ) = -999._wp
    inter%ctt           (1:ncol       ) = -999._wp
    inter%cwp           (1:ncol       ) = -999._wp
    inter%frac_in_up_cld(1:ncol,1:nlev) = 0._wp
    inter%reff          (1:ncol       ) = -999._wp
    inter%tau           (1:ncol       ) = -999._wp
    inter%tau_profile   (1:ncol,1:nlev) = -999._wp

    IF (options%sim%doClara) THEN
       inter%clara%flagged (1:ncol) = 0
       inter%Tb_cld_CDK    (1:ncol) = -999._wp
       IF (options%sim%Tb.EQ.1) THEN
          inter%cfrac (1:ncol,    1:nlev-1) = -999._wp
          inter%cloud (1:ncol,1:6,1:nlev-1) = -999._wp
          inter%Tb    (1:ncol,    1:nchan ) = -999._wp
       END IF
    ELSEIF (options%sim%doCloud_cci) THEN
       inter%cloud_cci%albedoIsDefined(1:ncol) = .FALSE.
       inter%cloud_cci%cla_vis006     (1:ncol) = -999._wp
       inter%cloud_cci%cth_c          (1:ncol) = -999._wp
       inter%cloud_cci%ctp_c          (1:ncol) = -999._wp
       inter%cloud_cci%ctt_c          (1:ncol) = -999._wp
    ELSEIF (options%sim%doRTTOV) THEN
       inter%cfrac (1:ncol,    1:nlev-1) = -999._wp
       inter%cloud (1:ncol,1:6,1:nlev-1) = -999._wp
       inter%Tb    (1:ncol,    1:nchan ) = -999._wp
    END IF

  END SUBROUTINE INITIALISE_INTERNAL_SIMULATOR
  SUBROUTINE DEALLOCATE_INTERNAL_SIMULATOR(inter,options)

    IMPLICIT NONE
    TYPE(internal), INTENT(inout) :: inter
    TYPE(name_list), INTENT(in)   :: options

        DEALLOCATE(inter%cflag    ,&
              inter%cph           ,&
              inter%cth           ,&
              inter%ctp           ,&
              inter%ctt           ,&
              inter%cwp           ,&
              inter%frac_in_up_cld,&
              inter%reff          ,&
              inter%tau           ,&
              inter%tau_profile    )

    IF (options%sim%doClara) THEN
       DEALLOCATE(inter%clara%flagged,&
            inter%Tb_cld_CDK         ,&
            inter%Tb)
       IF (options%sim%Tb.EQ.1) THEN
          DEALLOCATE(inter%cfrac,&
               inter%cloud) 
       END IF
    ELSEIF (options%sim%doCloud_cci) THEN
       DEALLOCATE(inter%cloud_cci%albedoIsDefined,&
                 inter%cloud_cci%cla_vis006      ,&
                 inter%cloud_cci%cth_c           ,&
                 inter%cloud_cci%ctp_c           ,&
                 inter%cloud_cci%ctt_c           )
    ELSEIF (options%sim%doRTTOV) THEN
       DEALLOCATE(inter%cfrac,&
                 inter%cloud ,&
                 inter%Tb     )
    END IF

  END SUBROUTINE deallocate_internal_simulator
END MODULE internal_simulator
