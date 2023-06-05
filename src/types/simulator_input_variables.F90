MODULE simulator_input_variables
  ! All the input varibles for the satellite simulators are gathered in
  ! this function
  !
  ! Salomon.Eliasson@smhi.se

  USE namelist_input,  ONLY: name_list
  USE satellite_specs, ONLY: satellite
  USE cosp_kinds,      ONLY: wp

  IMPLICIT NONE

  PUBLIC :: allocate_sim_input, &
       initialise_sim_input, &
       deallocate_sim_input

  ! -----------------
  ! Derived
  !
  ! Variables derived from the model and on the model grid
  ! 'cloud_emis' = layer in-cloud broadband IR emissivity
  ! 'data_mask'  = combined mask (e.g. isL2b .and. sunlit)
  ! 'g'          = Combined liquid and ice asymmetry parameter
  ! 'height'     = height at model layer interfaces
  ! 'sunlit'     = 0 if nighttime, 1 if sunlit
  ! 'inv_layers' = levels that contain inversions (including
  !                    tropopause)
  ! 'irradiance' = Black body emission from model vertical layer
  !                temperature
  ! 'ireff'    = Layer cloud effective radius for ice [micron]
  ! 'itau'     = Layer cloud optical thickness for ice clouds []
  ! 'iwc'      = Average IWC between 2 layers (for RTTOV) [g/m^3]
  ! 'iwp'      = Layer ice water path [g/m^2]
  ! 'isL2b'    = sampling mask for level 2b data format
  ! 'lreff'    = Layer cloud effective radius for liquid [micron]
  ! 'ltau'     = Layer cloud optical thickness for liquid clouds []
  ! 'lwc'      = Average LWC between 2 layers (for RTTOV) [g/m^3]
  ! 'lwp'      = Layer liquid water path [g/m^2]
  ! 'rho'      = density of the dry atmosphere [kg/m^3]
  ! 'p_int'    = pressure at interfaces [Pa]
  ! 'p_mid'    = pressure at midpoints [Pa]
  ! 'Q_kgm2'   = Column water vapour (for corrected temperature profile)
  ! 'rho'      = density of dry air [kg/m^3]
  ! 'satzen'   = The latitudinal average satellite viewing angle (based
  !              on level 2b datasets)
  ! 'solzen'   = The solar zenith angle
  ! 'surfType' = RTTOV surface type (0=land,1=open ocean,2=sea ice)
  ! 'surf_irradiance' = Emission from the Earth's surface
  ! 'surf_IR_em'= broadband IR surface emission
  ! 'tau'      = Layer cloud optical thickness []
  ! 'Tb_clr_CDK' = Broadband IR Tb without clouds
  ! 'Tcorr'     = Simulated layer brightness temperature from model
  !               temperature and absorption by model layer water
  !               vapour
  ! 'wv_emis'   = Water vapour emissivity per vertical layer
  ! 'w0'        = Combined liquid and ice single-scattering albedo

  TYPE subset
     ! model
     REAL(wp),ALLOCATABLE, DIMENSION(:,:) :: cloud_emis
     LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: data_mask
     REAL(wp),ALLOCATABLE, DIMENSION(:,:) :: g0
     REAL(wp),ALLOCATABLE, DIMENSION(:,:) :: height
     INTEGER, ALLOCATABLE, DIMENSION(:)   :: sunlit
     INTEGER, ALLOCATABLE, DIMENSION(:,:) :: inv_layers
     REAL(wp),ALLOCATABLE, DIMENSION(:,:) :: irradiance
     REAL(wp),ALLOCATABLE, DIMENSION(:,:) :: ireff
     REAL(wp),ALLOCATABLE, DIMENSION(:,:) :: itau
     REAL(wp),ALLOCATABLE, DIMENSION(:,:) :: iwc
     REAL(wp),ALLOCATABLE, DIMENSION(:,:) :: iwp
     LOGICAL, ALLOCATABLE, DIMENSION(:)   :: isL2b
     REAL(wp),ALLOCATABLE, DIMENSION(:,:) :: lreff
     REAL(wp),ALLOCATABLE, DIMENSION(:,:) :: ltau
     REAL(wp),ALLOCATABLE, DIMENSION(:,:) :: lwc
     REAL(wp),ALLOCATABLE, DIMENSION(:,:) :: lwp
     REAL(wp),ALLOCATABLE, DIMENSION(:,:) :: p_int
     REAL(wp),ALLOCATABLE, DIMENSION(:,:) :: p_mid
     REAL(wp),ALLOCATABLE, DIMENSION(:,:) :: Q_kgm2
     REAL(wp),ALLOCATABLE, DIMENSION(:)   :: solzen
     INTEGER, ALLOCATABLE, DIMENSION(:)   :: surfType
     REAL(wp),ALLOCATABLE, DIMENSION(:)   :: surf_irradiance
     REAL(wp),ALLOCATABLE, DIMENSION(:)   :: surf_IR_em
     REAL(wp),ALLOCATABLE, DIMENSION(:,:) :: tau
     REAL(wp),ALLOCATABLE, DIMENSION(:)   :: Tb_clr_CDK
     REAL(wp),ALLOCATABLE, DIMENSION(:,:) :: Tcorr
     REAL(wp),ALLOCATABLE, DIMENSION(:,:) :: w0
     INTEGER, ALLOCATABLE, DIMENSION(:)   :: waterType
     REAL(wp),ALLOCATABLE, DIMENSION(:,:) :: wv_emis

  END TYPE subset

contains

  ! -------------
  ! vectorised model
  ! -------------

  SUBROUTINE allocate_sim_input(sub,options,ngrids,nlev)

    IMPLICIT NONE

    TYPE(subset), INTENT(out)   :: sub
    INTEGER, INTENT(in)         :: ngrids,nlev
    TYPE(name_list), INTENT(in) :: options

    ALLOCATE( sub%data_mask        (ngrids,nlev  ),&
         sub%height           (ngrids,nlev+1),&
         sub%sunlit           (ngrids       ),&
         sub%inv_layers (ngrids,nlev  ),&
         sub%ireff            (ngrids,nlev  ),&
         sub%itau             (ngrids,nlev  ),&
         sub%iwc              (ngrids,nlev-1),&
         sub%iwp              (ngrids,nlev  ),&
         sub%isL2b            (ngrids       ),&
         sub%lreff            (ngrids,nlev  ),&
         sub%ltau             (ngrids,nlev  ),&
         sub%lwc              (ngrids,nlev-1),&
         sub%lwp              (ngrids,nlev  ),&
         sub%p_int            (ngrids,nlev+1),&
         sub%p_mid            (ngrids,nlev  ),&
         sub%tau              (ngrids,nlev  ),&
         sub%Tcorr            (ngrids,nlev+1) ) ! include the surface
    IF (.NOT. options%sim%doRTTOV) THEN
       ALLOCATE(sub%solzen(ngrids))
    END IF
    IF (options%sim%doClara .OR. options%sim%doISCCP) THEN
       ALLOCATE( sub%cloud_emis     (ngrids,nlev ),&
            sub%irradiance     (ngrids,nlev ),&
            sub%Q_kgm2         (ngrids,nlev ),&
            sub%surf_irradiance(ngrids      ),&
            sub%surf_IR_em     (ngrids      ) )
    END IF
    IF (options%sim%doCLARA .OR. options%sim%doCloud_cci) THEN
       ALLOCATE( sub%g0          (ngrids,nlev    ),&
            sub%w0          (ngrids,nlev    ) )
    END IF
   IF (options%sim%doCLARA .OR. options%sim%doRTTOV) THEN
      ALLOCATE(sub%surfType    (ngrids         ),&
               sub%waterType   (ngrids         ) )
   END IF

  END SUBROUTINE allocate_sim_input
  SUBROUTINE initialise_sim_input(sub,options,ngrids,nlev)

    IMPLICIT NONE

    TYPE(subset), INTENT(inout)   :: sub
    TYPE(name_list),INTENT(in)    :: options
    INTEGER, INTENT(in)           :: ngrids,nlev

    sub%data_mask       (1:ngrids,1:nlev   )= .FALSE.
    sub%height          (1:ngrids,1:nlev+1 )= 0._wp
    sub%sunlit          (1:ngrids          )= 0
    sub%inv_layers(1:ngrids,1:nlev   )= 0
    sub%ireff           (1:ngrids,1:nlev   )= -999._wp
    sub%itau            (1:ngrids,1:nlev   )= 0._wp
    sub%iwc             (1:ngrids,1:nlev-1 )= 0._wp
    sub%iwp             (1:ngrids,1:nlev   )= 0._wp
    sub%isL2b           (1:ngrids          )= .FALSE.
    sub%lreff           (1:ngrids,1:nlev   )= -999._wp
    sub%ltau            (1:ngrids,1:nlev   )= 0._wp
    sub%lwc             (1:ngrids,1:nlev-1 )= 0._wp
    sub%lwp             (1:ngrids,1:nlev   )= 0._wp
    sub%p_int           (1:ngrids,1:nlev+1 )= -999._wp
    sub%p_mid           (1:ngrids,1:nlev   )= -999._wp
    sub%tau             (1:ngrids,1:nlev   )= 0._wp
    sub%Tcorr           (1:ngrids,1:nlev+1 )= -999._wp
    IF (.NOT. options%sim%doRTTOV) &
         sub%solzen     (1:ngrids          )= -999._wp

    IF (options%sim%doClara .OR. options%sim%doISCCP) THEN
       sub%cloud_emis     (1:ngrids,1:nlev   )= -999._wp
       sub%irradiance     (1:ngrids,1:nlev   )= -999._wp
       sub%Q_kgm2         (1:ngrids,1:nlev   )= -999._wp
       sub%surf_irradiance(1:ngrids          )= -999._wp
       sub%surf_IR_em     (1:ngrids          )= 0.98_wp
    END IF
    IF (options%sim%doCLARA .OR. options%sim%doCloud_cci) THEN
       sub%g0             (1:ngrids,1:nlev   )= -999._wp
       sub%w0             (1:ngrids,1:nlev   )= -999._wp
    END IF
    IF (options%sim%doCLARA .OR. options%sim%doRTTOV) THEN
       sub%surfType       (1:ngrids          )= -999
       sub%waterType      (1:ngrids          )= -999
    END IF

  END SUBROUTINE initialise_sim_input

  SUBROUTINE deallocate_sim_input(sub,options)

    IMPLICIT NONE
    TYPE(subset), INTENT(inout) :: sub
    TYPE(name_list), INTENT(in) :: options

    ! model
    DEALLOCATE(sub%data_mask  ,&
         sub%height           ,&
         sub%sunlit           ,&
         sub%inv_layers ,&
         sub%ireff            ,&
         sub%itau             ,&
         sub%iwc              ,&
         sub%iwp              ,&
         sub%isL2b            ,&
         sub%lreff            ,&
         sub%ltau             ,&
         sub%lwc              ,&
         sub%lwp              ,&
         sub%p_int            ,&
         sub%p_mid            ,&
         sub%tau              ,&
         sub%Tcorr            )
    IF (.NOT. options%sim%doRTTOV) DEALLOCATE(sub%solzen)
    IF (options%sim%doClara .or. options%sim%doISCCP) THEN
       DEALLOCATE(sub%cloud_emis     ,&
            sub%irradiance     ,&
            sub%Q_kgm2         ,&
            sub%surf_irradiance,&
            sub%surf_IR_em       )

    END IF
    IF (options%sim%doCLARA .OR. options%sim%doCloud_cci) THEN
       DEALLOCATE(sub%g0 ,&
            sub%w0 )
    END IF
    IF (options%sim%doCLARA .OR. options%sim%doRTTOV) THEN
       DEALLOCATE(Sub%surfType  ,&
            sub%waterType )
    END IF

  END SUBROUTINE deallocate_sim_input

END MODULE simulator_input_variables
