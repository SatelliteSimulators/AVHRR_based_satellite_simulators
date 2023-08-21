MODULE optics_m
  ! Look up tables for the microphysics
  !
  ! Salomon.Eliasson@smhi.se
  !

  USE cosp_kinds, ONLY: wp


  INTEGER, PARAMETER :: num_trial_res = 20
  INTEGER, PARAMETER :: n_POD_edges = 20
  INTEGER, PARAMETER :: small_liq=10
  INTEGER, PARAMETER :: small_ice=30

  TYPE polyfit
     REAL(wp), ALLOCATABLE :: g0_small(:)
     REAL(wp), ALLOCATABLE :: w0_small(:)
     REAL(wp), ALLOCATABLE :: g0_large(:)
     REAL(wp), ALLOCATABLE :: w0_large(:)
  END type polyfit

  TYPE effective_radius
     ! single scattering albedo and assymetry parameter at 1.6 micron
     ! and 3.7 micron

     ! fit functions
     TYPE(polyfit) :: fit

     ! effective radius in LUT
     REAL(wp), ALLOCATABLE :: Re(:)

     ! trial effective radius's
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
  TYPE simulator_aux
     TYPE(optics_LUT)      :: LUT
     REAL(wp), ALLOCATABLE :: POD_layers         (:,:,:)
     REAL(wp), ALLOCATABLE :: POD_tau_bin_edges  (:)
     REAL(wp), ALLOCATABLE :: POD_tau_bin_centers(:)
     REAL(wp), ALLOCATABLE :: random_numbers     (:,:,:)
     REAL(wp), ALLOCATABLE :: FAR                (:,:)
  END TYPE simulator_aux
  !
  ! ------

  PUBLIC :: populate_effective_radius_LUT, &
       deallocate_optics

CONTAINS

  ELEMENTAL FUNCTION populate_effective_radius_LUT(CDR,phase,is_ch3B)&
       & RESULT(re)

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

   CHARACTER(*),INTENT(in) :: CDR
   INTEGER, INTENT(in) :: phase ! liquid=1, ice=2
   LOGICAL, INTENT(in) :: is_ch3B
   TYPE(effective_radius):: re
   INTEGER :: i

   SELECT CASE(CDR)
   CASE('clara-a2')
      IF (phase.EQ.1) re%nRe = 8
      IF (phase.EQ.2) re%nRe = 9
   CASE('clara-a3','claas-3')
      IF (phase.EQ.1) re%nRe = 8
      IF (phase.EQ.2) re%nRe = 11
   CASE('cloud_cci')
      IF (phase.EQ.1) re%nRe = 12
      IF (phase.EQ.2) re%nRe = 23
   END SELECT

   ALLOCATE(re%Re(re%nRe))
   ALLOCATE(re%trial(num_trial_res))

   ! In COSP y=ax^2+bx+c is represented as y=c+x*(b+x*a) and
   SELECT CASE(CDR)

   CASE('clara-a2')

      ! For liquid, small means less than 10 microns, large is greater
      IF (phase.EQ.1) THEN! liquid
         re%Re    = (/3.00,4.25,6.00,8.50,12.00,17.00,24.00,34.00/)
         IF (is_ch3B) THEN
            ALLOCATE(re%fit%g0_small(4),re%fit%g0_large(4))
            re%fit%g0_small = (/0.917106,-0.051528, 0.005521,-0.000162/)
            re%fit%g0_large = (/0.649317, 0.020399,-0.000658, 0.000008/)
            ALLOCATE(re%fit%w0_small(4),re%fit%w0_large(3))
            re%fit%w0_small = (/0.992618,-0.001896,-0.001386, 0.000070/)
            re%fit%w0_large = (/0.990696,-0.009177, 0.000080/)
         ELSE
            ALLOCATE(re%fit%g0_small(4),re%fit%g0_large(4))
            re%fit%g0_small = (/0.754852, 0.013001,-0.000500,0.000007/)
            re%fit%g0_large = (/0.754852, 0.013001,-0.000500,0.000007/)
            ALLOCATE(re%fit%w0_small(3),re%fit%w0_large(3))
            re%fit%w0_small = (/1.000119,-0.000630, 0.000002/)
            re%fit%w0_large = (/1.000119,-0.000630, 0.000002/)
         END IF
      ELSEIF (phase.EQ.2) THEN !ice
         ! For ice, small means less than 30 microns, large is greater
         re%Re    = (/5.00,7.07,10.00,14.15,20.00,28.28,40.00,56.58,80.00/)
         IF (is_ch3B) THEN
            ALLOCATE(re%fit%g0_small(4),re%fit%g0_large(3))
            re%fit%g0_small = (/0.682793, 0.017103,-0.000456,0.000005/)
            re%fit%g0_large = (/0.829044, 0.003454,-0.000023/)
            ALLOCATE(re%fit%w0_small(4),re%fit%w0_large(3))
            re%fit%w0_small = (/0.984572,-0.023210, 0.000538,-0.000005/)
            re%fit%w0_large = (/0.773026,-0.005482, 0.000035/)
         ELSE
            ALLOCATE(re%fit%g0_small(3),re%fit%g0_large(3))
            re%fit%g0_small = (/0.774316, 0.001323,-0.000001/)
            re%fit%g0_large = (/0.774316, 0.001323,-0.000001/)
            ALLOCATE(re%fit%w0_small(3),re%fit%w0_large(3))
            re%fit%w0_small = (/0.999659,-0.002548, 0.000006/)
            re%fit%w0_large = (/0.999659,-0.002548, 0.000006/)
         END IF
      END IF

   CASE('clara-a3','claas-3')

      ! For liquid, small means less than 10 microns, large is greater
      IF (phase.EQ.1) THEN! liquid
         re%Re    = (/3.00,4.25,6.00,8.50,12.00,17.00,24.00,34.00/)
         IF (is_ch3B) THEN
            ALLOCATE(re%fit%g0_small(3),re%fit%g0_large(4))
            re%fit%g0_small = (/0.919734, -0.038016, 0.002498/)
            re%fit%g0_large = (/0.611819, 0.027453, -0.000983, 0.000012/)
            ALLOCATE(re%fit%w0_small(4),re%fit%w0_large(3))
            re%fit%w0_small = (/0.986414, 0.002607, -0.002204, 0.000108/)
            re%fit%w0_large = (/0.988073, -0.009212, 0.000078/)
         ELSE
            ALLOCATE(re%fit%g0_small(3),re%fit%g0_large(4))
            re%fit%g0_small = (/0.729135, 0.018741, -0.000696/)
            re%fit%g0_large = (/0.788729, 0.008174, -0.000279, 0.000003/)
            ALLOCATE(re%fit%w0_small(3),re%fit%w0_large(3))
            re%fit%w0_small = (/1.000046, -0.000625, 0.000002/)
            re%fit%w0_large = (/1.000046, -0.000625, 0.000002/)
         END IF
      ELSEIF (phase.EQ.2) THEN !ice
         ! For ice, small means less than 30 microns, large is greater
         re%Re   = (/5.00,7.50,10.00,12.50,15.00,20.00,25.00,30.00,40.00,50.00,60.00/)
         IF (is_ch3B) THEN
            ALLOCATE(re%fit%g0_small(4),re%fit%g0_large(3))
            re%fit%g0_small = (/0.784704, -0.000992, 0.000357, -0.000007/)
            re%fit%g0_large = (/0.774883, 0.004755, -0.000035/)
            ALLOCATE(re%fit%w0_small(4),re%fit%w0_large(3))
            re%fit%w0_small = (/0.970288, -0.017169, 0.000313, -0.000002/)
            re%fit%w0_large = (/0.853475, -0.007290, 0.000053/)
         ELSE
            ALLOCATE(re%fit%g0_small(4),re%fit%g0_large(3))
            re%fit%g0_small = (/0.782730, -0.002480, 0.000173, -0.000003/)
            re%fit%g0_large = (/0.761282, 0.001005, -0.000004/)
            re%fit%w0_small = (/1.000128, -0.001906, 0.000005/)
            re%fit%w0_large = (/1.000128, -0.001906, 0.000005/)
         END IF
      END IF

   CASE('cloud_cci')
      ! --------------
      !    Cloud_cci
      ! --------------

      IF (phase.EQ.1) THEN! liquid
         re%Re = (/1,3,5,7,9,11,13,15,17,19,21,23/)
         IF (is_ch3B) THEN
            ! For liquid, small means less than 10 microns, large is greater
            ALLOCATE(re%fit%g0_small(5),re%fit%g0_large(3))
            re%fit%g0_small = (/0.348833, 0.366684,-0.097879,0.010294,-0.000371/)
            re%fit%g0_large = (/0.665379, 0.017114,-0.000356/)
            ALLOCATE(re%fit%w0_small(4),re%fit%w0_large(3))
            re%fit%w0_small = (/0.953224, 0.019062,-0.004726,0.000229/)
            re%fit%w0_large = (/0.991867,-0.009941, 0.000103/)
         ELSE
            ALLOCATE(re%fit%g0_small(5),re%fit%g0_large(3))
            re%fit%g0_small = (/0.894848,-0.087999, 0.023229,-0.002236,0.000074/)
            re%fit%g0_large = (/0.804219, 0.005214,-0.000102/)
            ALLOCATE(re%fit%w0_small(3),re%fit%w0_large(5))
            re%fit%w0_small = (/0.992141,-0.000624, 0.000004/)
            re%fit%w0_large = (/1.608582,-0.177541, 0.018427,-0.000824,0.000013/)
         END IF
      ELSEIF (phase.EQ.2) THEN !ice

         re%Re = (/4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64,68,&
                  72,76,80,84,88,92/)

         ! For ice, small means less than 30 microns, large is greater
         IF (is_ch3B) THEN
             ALLOCATE(re%fit%g0_small(4),re%fit%g0_large(4))
            re%fit%g0_small = (/ 0.518512, 0.019587,-0.000833, 0.000013/)
            re%fit%g0_large = (/-0.272164, 0.056633,-0.000913, 0.000005/)
             ALLOCATE(re%fit%w0_small(4),re%fit%w0_large(4))
            re%fit%w0_small = (/ 0.957291,-0.010000, 0.000218,-0.000004/)
            re%fit%w0_large = (/ 1.374435,-0.033709, 0.000493,-0.000002/)
         ELSE
            ALLOCATE(re%fit%g0_small(4),re%fit%g0_large(4))
            re%fit%g0_small = (/0.688700, 0.011752,-0.000512,0.000008/)
            re%fit%g0_large = (/0.389050, 0.023442,-0.000378,0.000002/)
            ALLOCATE(re%fit%w0_small(3),re%fit%w0_large(3))
            re%fit%w0_small = (/1.001719,-0.002956, 0.000008/)
            re%fit%w0_large = (/1.001719,-0.002956, 0.000008/)
         END IF
      END IF
   END SELECT
   re%min = MINVAL(re%Re)
   re%max = MAXVAL(re%Re)
   re%trial = Re%min + (Re%max - Re%min)/ &
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
         LUT%water%optics%re%Re      ,&
         LUT%water%optics%re%trial   ,&
         LUT%water%optics%re%fit%g0_small,&
         LUT%water%optics%re%fit%w0_small,&
         LUT%water%optics%re%fit%g0_large,&
         LUT%water%optics%re%fit%w0_large)
    DEALLOCATE(&
         LUT%ice%optics%re%Re      ,&
         LUT%ice%optics%re%trial   ,&
         LUT%ice%optics%re%fit%g0_small,&
         LUT%ice%optics%re%fit%w0_small,&
         LUT%ice%optics%re%fit%g0_large,&
         LUT%ice%optics%re%fit%w0_large)

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
