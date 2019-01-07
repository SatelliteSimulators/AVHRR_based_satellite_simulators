MODULE optics_m
  ! Look up tables for the microphysics
  !
  ! Salomon.Eliasson@smhi.se
  !

  USE cosp_kinds, ONLY: wp


  INTEGER, PARAMETER :: num_trial_res = 20
  INTEGER, PARAMETER :: n_POD_edges = 20

  TYPE effective_radius


     ! single scattering albedo and assymetry parameter at 1.6 micron
     ! and 3.7 micron

     REAL(wp), ALLOCATABLE :: g0_16(:)
     REAL(wp), ALLOCATABLE :: w0_16(:)
     REAL(wp), ALLOCATABLE :: g0_37(:)
     REAL(wp), ALLOCATABLE :: w0_37(:)

     ! effective radius in LUT
     REAL(wp), ALLOCATABLE :: re(:)

     ! trial effective radius's
     INTEGER :: num_trial_res
     REAL(wp), ALLOCATABLE :: trial(:)

     INTEGER  :: nRe ! number of effective radius's
     REAL(wp) :: min,max
  END TYPE effective_radius

  TYPE scattering_properties
     REAL(wp), ALLOCATABLE  :: g0(:)
     REAL(wp), ALLOCATABLE  :: w0(:)
     TYPE(effective_radius) :: re
  END TYPE scattering_properties

  TYPE cloud_albedo
     REAL(wp), ALLOCATABLE :: Albedo(:,:,:)
     REAL(wp), ALLOCATABLE :: Re    (:)
     REAL(wp), ALLOCATABLE :: Solzen(:)
     REAL(wp), ALLOCATABLE :: Tau   (:)
     REAL(wp)              :: dRe
     REAL(wp)              :: dSolzen
     REAL(wp)              :: dTau
     INTEGER               :: nRe
     INTEGER               :: nSolzen
     INTEGER               :: nTau
  END TYPE cloud_albedo
  TYPE cloud_phase
     ! Put the different look up tables here
     !
     TYPE(scattering_properties) :: optics
     TYPE(cloud_albedo)          :: albedo
  END TYPE cloud_phase

  TYPE GZero_t
     ! Defines GZero_m module used to hold information about the nearest
     ! neighbour grid point array indices when interpolating Look-Up Tables
     ! in the ECP. (Used for cloud albedo)

     INTEGER :: iT0, iT1     ! Nearest neighbour indices for Tau value
     INTEGER :: iTm1,iTp1    ! Next-nearest neighbour indices for Tau value
     INTEGER :: iR0, iR1     ! Nearest neighbour indices for Re value
     INTEGER :: iRm1, iRp1   ! Next-nearest neighbour indices for Re value
     INTEGER :: iSoZ0, iSoZ1 ! Nearest neighbour indices for Sun zen value

     REAL(wp):: dT           ! Fraction of grid step in Tau from
     ! zero'th point to current Tau value
     REAL(wp):: dR           ! Fraction of grid step to current Re
     REAL(wp):: dSoZ         ! Fraction of grid step to current Sun zen

     REAL(wp):: T1           ! 1.0 - dT (stored for frequent use)
     REAL(wp):: R1           ! 1.0 - dR (stored for frequent use)
     REAL(wp):: So1          ! 1.0 - dSuZ (stored for frequent use)
  END TYPE GZero_t

  TYPE optics_LUT
     TYPE(cloud_phase) :: ice
     TYPE(cloud_phase) :: water
     TYPE(GZero_t)     :: GZero
  END TYPE optics_LUT

  ! for cloud fraction
  !
  TYPE sim_aux

     TYPE(optics_LUT)      :: LUT
     REAL(wp), ALLOCATABLE :: detection_limit    (:,:)
     REAL(wp), ALLOCATABLE :: POD_layers         (:,:,:)
     REAL(wp), ALLOCATABLE :: POD_tau_bin_edges  (:)
     REAL(wp), ALLOCATABLE :: POD_tau_bin_centers(:)
     REAL(wp), ALLOCATABLE :: random_numbers     (:,:,:)
  END TYPE sim_aux
  !
  ! ------

  PUBLIC :: populate_effective_radius_LUT, &
       deallocate_optics

CONTAINS

  ELEMENTAL FUNCTION populate_effective_radius_LUT(sim,phase) RESULT(Re)

    ! ----------
    ! CLARA LUT's
    ! based on:
    !
    ! Mie theory for two-PARAMETER gamma distributions (effective variance
    ! 0.15) of spherical droplets (W. De Rooij and C. Stap, van der,
    ! ‘Expansion of Mie scattering matrices in generalized spherical
    ! functions’. Astron. Astrophys. 131, 237–248 (1984)
    !
    ! .and.
    !
    ! Raytracing using geometric optics approximation for imperfect,
    ! randomly oriented, hexagonal, 'mono-size' ice crystals (Hess, H,
    ! R. B. A. Koelemeijer, and P. Stammes, 1998: Scattering matrices
    ! of imperfect hexagonal
    ! crystals. J. Quant. Spectrosc. Radiat. Transfer, 60, 301–308)
    ! ----------
    ! Cloud_cci LUT's
    ! based on:
    !
    ! fixme: more information...
    !
    ! Mie Theory for liquid
    !
    ! .and.  for ice
    !
    ! Baran, A. J., and S. Havemann. "The dependence of
    ! retrieved cirrus ice‐crystal effective dimension on assumed
    ! ice‐crystal geometry and size‐distribution function at solar
    ! wavelengths." Quarterly Journal of the Royal Meteorological
    ! Society 130.601 (2004): 2153-2167.
    !


    IMPLICIT NONE

    CHARACTER(*),INTENT(in) :: sim
    INTEGER, INTENT(in) :: phase ! liquid=1, ice=2
    TYPE(effective_radius):: re
    INTEGER :: i

    Re%num_trial_res=num_trial_res

    SELECT CASE(sim)
    CASE('clara')
       IF (phase.EQ.1) Re%nRe = 8
       IF (phase.EQ.2) Re%nRe = 9
    CASE('cloud_cci')
       IF (phase.EQ.1) Re%nRe = 12
       IF (phase.EQ.2) Re%nRe = 23
    END SELECT

    ALLOCATE(Re%re   (Re%nRe)      ,&
             Re%g0_16(Re%nRe)      ,&
             Re%w0_16(Re%nRe)      ,&
             Re%g0_37(Re%nRe)      ,&
             Re%w0_37(Re%nRe)      ,&
             Re%trial(num_trial_res))

    SELECT CASE(sim)

    CASE('clara')
       ! --------------
       !    CLARA
       ! -------------- 

       IF (phase.EQ.1) THEN! liquid

          Re%re    = (/3.00,4.25,6.00,8.50,12.00,17.00,24.00,34.00/)

          Re%g0_37 = (/0.806269,0.789406,0.767975,0.779810,0.813295,&
               0.843947,0.868019,0.889109/)

          Re%g0_16 = (/0.794294,0.794515,0.815903,0.836549,0.850984,&
               0.862017,0.870501,0.877463/)

          Re%w0_37 = (/0.976040,0.965643,0.945708,0.919468,0.890786,&
               0.857294,0.817537,0.770678/)

          Re%w0_16 = (/0.998404,0.997500,0.996338,0.994817,0.992806,&
               0.990160,0.986504,0.981392/)

       ELSEIF (phase.EQ.2) THEN !ice

          Re%re    = (/5.00,7.07,10.00,14.15,20.00,28.28,40.00,56.58,&
               80.00/)

          Re%g0_37 = (/0.756000,0.784000,0.814000,0.846000,0.878000,&
               0.907000,0.932000,0.948000,0.956000/)

          Re%g0_16 = (/0.779000,0.784000,0.789000,0.794000,0.801000,&
               0.809000,0.824000,0.848000,0.873000/)

          Re%w0_37 = (/0.882300,0.845000,0.800300,0.750300,0.698000,&
               0.647700,0.606200,0.577700,0.558300/)

          Re%w0_16 = (/0.987200,0.981900,0.974700,0.964800,0.951200,&
               0.932800,0.908100,0.875500,0.836500/)

       END IF

    CASE('cloud_cci')
       ! --------------
       !    Cloud_cci
       ! -------------- 

       IF (phase.EQ.1) THEN! liquid

          Re%re = (/1.00000,3.00000,5.00000,7.00000,9.00000,11.0000,&
               13.0000,15.0000,17.0000,19.0000,21.0000,23.0000/)

          Re%g0_37 = (/0.626871,0.819327,0.783379,0.766986,0.788780,&
               0.811922,0.829410,0.842919,0.852329,0.860446,0.866794&
               ,0.872355/)

          Re%g0_16 = (/0.827948,0.785390,0.802836,0.828255,0.842266,&
               0.849973,0.855060,0.859668,0.863104,0.866163,0.868552&
               ,0.870674/)

          Re%w0_37 = (/0.966989,0.976486,0.957191,0.932168,0.910892,&
               0.894659,0.880104,0.865600,0.852920,0.839949,0.828469&
               ,0.817487/)

          Re%w0_16 = (/0.991526,0.990282,0.989145,0.987935,0.986817,&
               0.985718,0.966989,0.976486,0.957191,0.932168,0.910892&
               ,0.894659/)


       ELSEIF (phase.EQ.2) THEN !ice

          Re%re = (/4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64,68,&
               72,76,80,84,88,92/)

          Re%g0_37 = (/ 0.581337,0.635218,0.656466,0.670513,0.682161,&
               0.696564,0.714014,0.731463,0.805473,0.832767,0.860062,&
               0.870130,0.861158,0.852185,0.843213,0.834240,0.825268,&
               0.816295,0.807323,0.810883,0.821191,0.831500,0.841809/)

          Re%g0_16 = (/ 0.726038,0.758328,0.769519,0.776531,0.782151,&
               0.788821,0.796652,0.804483,0.836676,0.847229,0.857781,&
               0.861654,0.858143,0.854632,0.851121,0.847610,0.844099,&
               0.840589,0.837078,0.838613,0.842866,0.847119,0.851372/)

          Re%w0_37 = (/ 0.921272,0.887802,0.862080,0.838227,0.814996,&
               0.790173,0.763588,0.737004,0.672871,0.653529,0.634187,&
               0.622179,0.618279,0.614378,0.610477,0.606577,0.602676,&
               0.598775,0.594875,0.589266,0.582739,0.576211,0.569683/)

          Re%w0_16 = (/ 0.987856,0.977456,0.967581,0.957769,0.947977,&
               0.938279,0.928684,0.919088,0.903135,0.893414,0.883692,&
               0.875283,0.868325,0.861367,0.854409,0.847451,0.840494,&
               0.833536,0.826578,0.819295,0.811838,0.804381,0.796923/)

       END IF
    END SELECT

    Re%min = MINVAL(Re%re)
    Re%max = MAXVAL(Re%re)

    Re%trial = Re%min + (Re%max - Re%min)/ &
         (num_trial_res-1) * (/ (i - 1, i = 1, num_trial_res) /)

  END FUNCTION populate_effective_radius_LUT

  SUBROUTINE deallocate_optics(LUT)

    IMPLICIT NONE 

    TYPE(optics_LUT), INTENT(inout) :: LUT

    DEALLOCATE(&
         LUT%water%optics%g0,&
         LUT%water%optics%w0,&
         LUT%ice%optics%g0,&
         LUT%ice%optics%w0)

    DEALLOCATE(&
         LUT%water%optics%re%re,   &
         LUT%water%optics%re%trial,&
         LUT%water%optics%re%g0_16,&
         LUT%water%optics%re%w0_16,&
         LUT%water%optics%re%g0_37,&
         LUT%water%optics%re%w0_37)

    DEALLOCATE(&
         LUT%ice%optics%re%re,   &
         LUT%ice%optics%re%trial,&
         LUT%ice%optics%re%g0_16,&
         LUT%ice%optics%re%w0_16,&
         LUT%ice%optics%re%g0_37,&
         LUT%ice%optics%re%w0_37)

    IF (ALLOCATED(LUT%ice%albedo%Albedo)) THEN
       DEALLOCATE ( &
            LUT%ice%albedo%Albedo  ,&
            LUT%ice%albedo%Re      ,&
            LUT%ice%albedo%Solzen  ,&
            LUT%ice%albedo%Tau     ,&
            LUT%water%albedo%Albedo,&
            LUT%water%albedo%Re    ,&
            LUT%water%albedo%Solzen,&
            LUT%water%albedo%Tau)
    END IF

  END SUBROUTINE deallocate_optics

END MODULE optics_m
