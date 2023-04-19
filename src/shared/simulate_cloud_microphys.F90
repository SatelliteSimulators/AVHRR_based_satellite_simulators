MODULE SIMULATE_CLOUD_MICROPHYS

  ! 
  ! Salomon.Eliasson@smhi.se

  USE COSP_KINDS,                ONLY: &
       WP
  USE DATA_CHECK,                ONLY: &
       PRINT_SUB_POINT
  USE FROM_COSP2,                ONLY: &
       COMPUTE_TOA_REFLECTANCE,        &
       INTERPOLATE_TO_MIN,             &
       RE_FILL,                        &
       TWO_STREAM_REFLECTANCE
  USE INTERNAL_SIMULATOR,        ONLY: &
       INTERNAL
  USE MODEL_INPUT,               ONLY: &
       MODEL_TYPE
  USE MOD_COSP_STATS,            ONLY:&
       HIST2D
  USE NAMELIST_INPUT,            ONLY: &
       NAME_LIST
  USE OPTICS_M,                  ONLY: &
       NUM_TRIAL_RES,                  &
       OPTICS_LUT,                     &
       SCATTERING_PROPERTIES,          &
       SIMULATOR_AUX
  USE SIMULATOR_INPUT_VARIABLES, ONLY: &
       SUBSET

  IMPLICIT NONE

  PUBLIC :: CLOUD_WATER_PHASE,        &
       UPPER_CLOUD,                   &
       CLOUD_EFFECTIVE_RADIUS,        &
       CLOUD_EFFECTIVE_RADIUS_NATIVE, &
       GET_CLOUDTYPE,                 &
       GET_CLOUD_MICROPHYSICS,        &
       GET_PTAU

CONTAINS

  SUBROUTINE GET_CLOUD_MICROPHYSICS(d1,nc,nlev,frac_out,inter,sub,options)

    ! A combined routine for cloud microphysics that all our simulators are using 
    ! Salomon.eliasson@smhi.se

    IMPLICIT NONE

    INTEGER, INTENT(in)           :: d1,nc,nlev
    REAL(wp), INTENT(in)          :: frac_out(nc,nlev)
    TYPE(internal), INTENT(inout) :: inter
    TYPE(subset), INTENT(in)      :: sub
    TYPE(name_list), INTENT(in)   :: options
    INTEGER                       :: sunlit

    inter%tau_profile (1:nc,1:nlev) = &
         MERGE(SPREAD((sub%itau(d1,1:nlev)+sub%ltau(d1,1:nlev)),1,nc),&
         0._wp,&
         frac_out(1:nc,1:nlev) .GT. 0)
    
    inter%tau(1:nc)   = SUM(inter%tau_profile(1:nc,1:nlev),DIM=2)

    sunlit=sub%sunlit(d1)

    inter%cflag(1:nc) = GET_CLOUDTYPE(d1,nc,inter,options,sunlit)
 
    inter%frac_in_up_cld(1:nc,1:nlev) = UPPER_CLOUD(nc,nlev,inter)
    
    inter%cph(1:nc)                   = CLOUD_WATER_PHASE(d1,nc,nlev,inter,sub)
    
  END SUBROUTINE GET_CLOUD_MICROPHYSICS

  FUNCTION GET_PTAU(nc,n_tbins,n_pbins,options,tau,ctp,cph) RESULT(ptau)
    
    IMPLICIT NONE

    INTEGER, INTENT(in)                  :: nc,n_tbins,n_pbins
    REAL(wp), INTENT(in)                 :: tau(nc),ctp(nc)
    INTEGER, INTENT(in)                  :: cph(nc)
    TYPE(name_list), INTENT(in)          :: options
    REAL(wp), DIMENSION(nc)              :: tmptau,tmpctp
    REAL(wp), DIMENSION(n_tbins,n_pbins) :: tmpptau

    REAL(wp)                             :: ptau    (n_tbins,n_pbins,2)

    ! make all tau less than 100 to match the clara, and
    ! cloud_cci observations

    tmptau = MERGE(tau,99.999_wp,tau.LT.99.999_wp)
    tmpctp = ctp

    ! ptau for ice clouds only
    tmptau=MERGE(tmptau,-999._wp,cph.EQ.2)
    tmpctp=MERGE(tmpctp,-999._wp,cph.EQ.2)

    tmpptau = 0
    CALL HIST2D(tmptau,tmpctp,nc,&
         options%ctp_tau%tbin_edges,n_tbins,&
         options%ctp_tau%pbin_edges,n_pbins,&
         tmpptau)

    ptau(1:n_tbins,1:n_pbins,1) = tmpptau

    tmptau = MERGE(tau,99.999_wp,tau.LT.99.999_wp)
    tmpctp = ctp
    tmpptau = 0

    ! ptau for liquid clouds only
    tmptau=MERGE(tmptau,-999._wp,cph.EQ.1)
    tmpctp=MERGE(tmpctp,-999._wp,cph.EQ.1)
    CALL HIST2D(tmptau,tmpctp,nc,&
            options%ctp_tau%tbin_edges,n_tbins,&
            options%ctp_tau%pbin_edges,n_pbins,&
            tmpptau)

    ptau(1:n_tbins,1:n_pbins,2) = tmpptau
    
  END FUNCTION GET_PTAU
  
  FUNCTION UPPER_CLOUD(nc,nlev,inter)

    !
    ! Find the upper part of the cloud. This is used for the cloud
    ! phase retrieval and the native model effective radius
    ! retrieval. upper_cloud(layer) = fraction of the layer that is in
    ! the upper part of the cloud


    IMPLICIT NONE

    INTEGER, INTENT(in)        :: nc, nlev
    TYPE(internal), INTENT(in) :: inter


    REAL(wp) :: upper_cloud(nc,nlev)

    INTEGER :: inl
    REAL(wp)          :: ucLim!  [optical depth]
    REAL(wp), ALLOCATABLE :: ac_tau(:,:)

    ! 'tau_equivRadCldTop' = Cloud optical depth where the radiantly equivalent cloud
    !                        top height is located. The MODIS satellite
    !                        simulator assumes tau_equivRadCldTop = 1(pincus et. al., 2012)
    ! Now hard coded
    ucLim = 1 !options%cloudMicrophys%tau_equivRadCldTop
    ALLOCATE( ac_tau (nc,nlev+1) )
    ac_tau (1:nc,1:nlev+1) = 0._wp ! padding the top row with 0s,
                                   ! since a cloud might be in the
                                   ! second highest layer
    upper_cloud (1:nc,1:nlev) = 0._wp

    DO inl = 1, nlev  ! loop over levels from TOA to surface
       ac_tau(1:nc,inl+1) = SUM(inter%tau_profile(1:nc,2:inl),DIM=2)

       WHERE (inter%tau_profile(1:nc,inl).GT.0) 
          WHERE (ac_tau(1:nc,inl+1).LE.ucLim) 
             ! The whole layer is in the upper part of the cloud
             upper_cloud(1:nc,inl) = 1._wp
          ELSEWHERE
             ! passed the threshold of upper cloud. find fraction of
             ! "by-how-much" 
             WHERE (ac_tau(1:nc,inl).LT.ucLim)
                ! part of the upper cloud is in this model layer
                upper_cloud(1:nc,inl) = &
                     (ucLim-ac_tau(1:nc,inl))/ac_tau(1:nc,inl+1)
             END WHERE
          END WHERE
       END WHERE
    END DO

    DEALLOCATE( ac_tau )

  END FUNCTION UPPER_CLOUD

  FUNCTION CLOUD_WATER_PHASE(d1,nc,nlev,inter,sub) RESULT(phase)
    !
    ! Calculate the cloud phase.
    ! Determines the cloud phase from the top optical depth of the
    ! cloud provided by UPPER_CLOUD()
    !
    IMPLICIT NONE

    INTEGER, INTENT(in)        :: d1,nc,nlev
    TYPE(internal), INTENT(in) :: inter
    TYPE(subset), INTENT(in)   :: sub

    ! OUT
    INTEGER, DIMENSION(nc)     :: phase

    ! internal
    REAL(wp),DIMENSION(nc)     :: num, den
    REAL(wp),DIMENSION(nc,nlev):: ltau,itau

    ltau = SPREAD(sub%ltau(d1,1:nlev),1,nc)
    itau = SPREAD(sub%itau(d1,1:nlev),1,nc)

    num   = 0._wp
    den   = 0._wp
    phase = 0

    WHERE(inter%cflag(1:nc).GT.0)
       num(1:nc) = SUM(inter%frac_in_up_cld(1:nc,1:nlev) * &
            ( 1._wp*ltau(1:nc,1:nlev)+2._wp*itau(1:nc,1:nlev) ),&
            DIM=2)

       den(1:nc) = SUM(inter%frac_in_up_cld(1:nc,1:nlev) * &
            ( ltau(1:nc,1:nlev)+itau(1:nc,1:nlev) ),&
            DIM=2)

       phase(1:nc) =NINT( num(1:nc)/den(1:nc) )
    END WHERE
    
  END FUNCTION CLOUD_WATER_PHASE

  FUNCTION CLOUD_EFFECTIVE_RADIUS(d1,ins,nlev,sub,inter,optics)
    ! 
    ! Find the effective radius using the MODIS approach:

    !1) Compute model TOA reflectance
    !    
    !   a) Using polynomial fits, compute asymmetry parameter (g) and
    !   single-scattering albedo (ω0) using the model radii for ice and
    !   liquid, we look up single values of g_i, g_w, ω0_i and ω0_w.
    !
    !   b) Combine with optical depth profiles of ice and liquid to get
    !   total g and ω0.
    !
    !     g = (g_i * τ_i + g_w * τ_w)/ (τ_i + τ_w)
    !     ω0 = (ω0_i * τ_i + ω0_w * τ_w)/ (τ_i + τ_w)
    !
    !   c) Use the two-stream method to obtain the model reflectance
    !
    !2) Compute the predicted TOA reflectance for a range of radii .
    !
    !   a) Same as 1a, but using a range of radii for both ice and
    !   liquid. We get now have vectors for g_i, g_w, ω0_i and ω0_w.
    !
    !   b) USE two stream method to obtain reflectance, but only for cloud
    !   phase, which was determined earlier.  3) Minimize #2 to #1. The
    !   simulated r_e is found from the g(r_e) ω0(r_e) used to calculate
    !   minimized reflectance.
    !
    !
    ! Note: This should only be called if there is a cloud
    !
    !
    ! Salomon.Eliasson@smhi.se

    IMPLICIT NONE

    TYPE(scattering_properties), INTENT(in) :: optics ! either water or ice (decided before call)
    INTEGER, INTENT(in)          :: d1,ins,nlev
    TYPE(subset), INTENT(in)     :: sub
    TYPE(internal), INTENT(in)   :: inter

    REAL(wp), DIMENSION(nlev)    :: g0_column,w0_column,tau_column
    REAL(wp)                     :: obs_Refl_nir
    REAL(wp), DIMENSION(num_trial_res) :: predicted_Refl_nir, trial_re
    REAL(wp) :: re_max,re_min
    INTEGER  :: diffloc(1)
    ! OUT
    REAL(wp) :: cloud_effective_radius

    predicted_Refl_nir(1:num_trial_res) = 0._wp
    g0_column         (1:nlev         ) = 0._wp
    w0_column         (1:nlev         ) = 0._wp
    tau_column        (1:nlev         ) = 0._wp
    cloud_effective_radius              = 0._wp
    obs_Refl_nir                        = 0._wp

    ! --------------------------------
    ! COMPUTE THE MODEL TOA REFLECTANCE
    !
    WHERE(inter%tau_profile(ins,1:nlev) .GT. 0)
       ! get the combined liquid and ice parameters from where
       ! scops says it's cloudy
       g0_column (1:nlev) = sub%g0(d1,1:nlev)
       w0_column (1:nlev) = sub%w0(d1,1:nlev)
       tau_column(1:nlev) = inter%tau_profile(ins,1:nlev) 
    END WHERE

    obs_Refl_nir = COMPUTE_TOA_REFLECTANCE(nlev, &
         tau_column(1:nlev), g0_column(1:nlev), w0_column(1:nlev))

    ! --------------------------------

    ! --------------------------------
    ! GET THE MINIMIZED SOLUTION FOR A GIVEN PHASE
    ! 1) Calculate some predicted reflectances based on the
    ! trial effective radius's 
    ! 2) Minimize the trial refl vs. calculated from the model
    ! to find the best fitting effective radius for a pre-decided
    ! cloud phase
    ! 3) I'm keeping it in microns

    predicted_Refl_nir(1:num_trial_res)=TWO_STREAM_REFLECTANCE(inter%tau(ins),&
         optics%g0(1:num_trial_res),optics%w0(1:num_trial_res))
    
    re_min   = optics%re%min
    re_max   = optics%re%max
    trial_re = optics%re%trial
    
    cloud_effective_radius = INTERPOLATE_TO_MIN(trial_re(1:num_trial_res), &
         predicted_Refl_nir(1:num_trial_res), obs_Refl_nir)

    IF (cloud_effective_radius .EQ. re_fill ) THEN

       ! sometimes the closest value does not 'bracket the root',
       ! which INTERPOLATE_TO_MIN requires (so basically there are two
       ! solutions), and in this case pick the closest solution. If
       ! obs_Refl_nir is outside the range, picking the closest
       ! solution also works. Cloud_effective_radius will be set
       ! to re_min if obs_refl.ge.prediceted_refl(1) since the
       ! reflectivity is largest for the smallest particles, and visa
       ! versa for large particles

       diffloc = MINLOC(ABS(predicted_Refl_nir-obs_Refl_nir))
       cloud_effective_radius = trial_re(diffloc(1))
    END IF
  END FUNCTION CLOUD_EFFECTIVE_RADIUS

  FUNCTION CLOUD_EFFECTIVE_RADIUS_NATIVE(d1,ins,nlev,sub,inter) RESULT(ref)
    ! 
    ! Calculate the model native bulk effective radius. This is done
    ! by finding the effective radius from the upper 1 optical depth of
    ! the cloud.
    !
    !
    !
    IMPLICIT NONE

    INTEGER, INTENT(in)          :: d1,ins
    TYPE(subset), INTENT(in)     :: sub
    TYPE(internal), INTENT(in)   :: inter
    INTEGER, INTENT(in)          :: nlev
    REAL(wp) :: num, den

    ! OUT
    REAL(wp) :: ref

    num = 0._wp
    den = 0._wp
    ref = 0._wp

    IF (inter%tau(ins) .GE. 0) THEN
       ! 0 = cloud top weighted-effective radius simulation
       IF ( inter%cph(ins) == 1 ) THEN ! 1 == liquid
          num = SUM(inter%frac_in_up_cld(ins,1:nlev)*&
               sub%lreff(d1,1:nlev)*                 &
               sub%ltau (d1,1:nlev))
          den = SUM(inter%frac_in_up_cld(ins,1:nlev)*&
               sub%ltau (d1,1:nlev))
       ELSEIF ( inter%cph(ins) == 2 ) THEN ! 2 == ice
          num = SUM(inter%frac_in_up_cld(ins,1:nlev)*&
               sub%ireff(d1,1:nlev)*                 &
               sub%itau (d1,1:nlev))
          den = SUM(inter%frac_in_up_cld(ins,1:nlev)*&
               sub%itau (d1,1:nlev))
       END IF
       ref = num/den

    ELSE
       ref = re_fill
    END IF

  END FUNCTION CLOUD_EFFECTIVE_RADIUS_NATIVE
  
  FUNCTION GET_CLOUDTYPE(d1,nc,inter,options,sunlit) RESULT(cflag)
    
    ! Make a variable that saves the cloud type based purely on optical depth
    !
    ! Salomon.Eliasson@smhi.se

    IMPLICIT NONE

    INTEGER, INTENT(in)          :: d1,nc
    TYPE(internal), INTENT(in)   :: inter
    TYPE(name_list), INTENT(in)  :: options
    INTEGER, INTENT(in)          :: sunlit 

    ! OUT
    integer :: cflag(nc)

    ! internal
    REAL(wp)                      :: tau_min
    REAL(wp)                      :: opaque,ClFree
    INTEGER                       :: lvl

    ClFree=0._wp
    opaque  = 4._wp ! approximate limit for accurate CALIPSO optical depth retrievals

    !---------
    ! CLOUD TYPE
    !
    ! clear=0, sub-visible=1, semi-transparent=2, opaque=3

    SELECT CASE(options%cloudMicrophys%cf_method)
    CASE(0)
       tau_min = options%cloudMicrophys%tau_min

       WHERE(inter%tau(1:nc).LE.ClFree)
          cflag(1:nc) = 0 ! cloud free
       ELSEWHERE((inter%tau(1:nc).GT.ClFree).AND.(inter%tau(1:nc).LT.tau_min))
          cflag(1:nc) = 1 ! sub-visible
       ELSEWHERE((inter%tau(1:nc).GE.tau_min).AND.(inter%tau(1:nc).LT.opaque))
          cflag(1:nc) = 2 ! semi-transparent
       ELSEWHERE(inter%tau(1:nc).GE.opaque)
          cflag(1:nc) = 3 ! opaque
       END WHERE

    CASE(1)

       ! Combine the information on the probability of detection for a
       ! geographical location and for a optical depth
       ! interval. Mostly, the POD is higher the higher the optical
       ! depth. The column is considered cloudy if a random number, x, 
       ! assigned to the (lon,lat,column) before the simulation is
       ! higher than x>(1-POD). I.e., the higher the POD the more
       ! likely this cloud is classified as cloudy. The POD is provided
       ! based on comprehensive comparisons between the satellite
       ! dataset and Calipso data which is tightly collocated with the
       ! data (Karlsson & Håkansson 2018)

       ! Cycle through the optical depth bins to find the POD to
       ! assume, singel out the columns that fall into that optical
       ! depth bin, and compare their predetermined random number to
       ! the reciprocal of the POD

       DO lvl=1,SIZE(options%sim_aux%POD_tau_bin_edges)-1
          WHERE(inter%tau(1:nc).GT.options%sim_aux%POD_tau_bin_edges(lvl)&
               & .AND.inter%tau.LE.options%sim_aux%POD_tau_bin_edges(lvl+1))
             
             cflag(1:nc)=MERGE(3,0,&
                  options%sim_aux%random_numbers(sunlit+1,d1,1:nc).GE.&
                  (1-SPREAD(options%sim_aux%POD_layers(d1,lvl,sunlit+1),1,nc)))
          END WHERE
       END DO
       ! reclassify the cloudy columns to semi-transparent if the optical
       ! depth is less than opaque (I should be using RTTOV channel
       ! differences for CLARA)

       WHERE(inter%tau(1:nc).EQ.0)
          cflag(1:nc)=0
       END WHERE
       WHERE((cflag(1:nc).EQ.3).AND.(inter%tau(1:nc).LT.opaque))
          cflag(1:nc)=2
       END WHERE
    END SELECT

  END FUNCTION GET_CLOUDTYPE

END MODULE SIMULATE_CLOUD_MICROPHYS
