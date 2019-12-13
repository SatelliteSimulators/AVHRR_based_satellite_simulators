MODULE brightness_temperatures
  ! Different types of brightness temperature calculations
  !
  ! Salomon.Eliasson@smhi.se

  USE COSP_KINDS,                ONLY: WP
  USE DATA_CHECK,                ONLY:&
       MISSING
  USE INTERNAL_SIMULATOR,        ONLY: INTERNAL
  USE MOD_COSP_ISCCP_INTERFACE,  ONLY:&
       COSP_ISCCP_INIT
  USE MOD_ICARUS,                ONLY:&
       ICARUS_SUBCOLUMN
  USE MODEL_INPUT,               ONLY: MODEL_TYPE
  USE NAMELIST_INPUT,            ONLY: NAME_LIST 
  USE SIMULATOR_INPUT_VARIABLES, ONLY: SUBSET

  IMPLICIT NONE

  PUBLIC :: CORRECTED_TEMPERATURE_PROFILE, &
       GET_TB ,       &
       TB_ICARUS

CONTAINS

  FUNCTION CORRECTED_TEMPERATURE_PROFILE(M,S,ngrids,nlev) RESULT(Tcorr)

    ! Get the corrected temperature profile using the 10.8 micron
    ! channel. Currently this is only for the CLARA simulator, so I
    ! can assume that I'm using the AVHRR.
    !
    ! The model temperature profile is corrected using the difference
    ! between the simulated Tb and the model skin temperature (SKT), ∆T
    ! tot = Tb − SKT, together with the vertical profile of water vapour
    ! amount (Q) [kg m−2].
    !
    ! The assumption is that it is sufficiently accurate to assume that
    ! the temperature correction at a certain layer is equal to ∆T times
    ! the fraction of the Total Column Water Vapour (TCWV) in and above
    ! model vertical grid number i. This approach implies that the
    ! difference in temperature is due to absorption of water vapour, and
    ! all other trace gases have a negligible effect on Tb.
    !
    ! It is also assumed that only a small error is introduced by
    ! assuming that the fraction of layer temperature correction is equal
    ! to the fraction of the layer particle water vapour path to the total
    ! column water vapour. By extension this assumes that a linear change
    ! in TCWV causes a linear change in delta T. This is technically not
    ! correct since the change is linear w.r.t. radiance but non-linear
    ! w.r.t. Tb. However, for such small temperatures (delta T) at 10.8
    ! micron, the Rayleigh-Jeans is close enough (<0.1 K).
    !
    ! One caveat to this approach is that since we are only comparing
    ! to ∆T tot the effect of temperature inversions is not taken into
    ! account at each level, but also this error per level is only a small
    ! fraction of a degree.
    !
    ! Salomon.Eliasson@smhi.se

    IMPLICIT NONE

    ! in
    TYPE(model_type), INTENT(in) :: M
    TYPE(subset),     INTENT(in) :: S
    INTEGER,          INTENT(in) :: ngrids,nlev

    ! internal
    REAL(wp), ALLOCATABLE, DIMENSION(:)   :: delT,TCWV
    REAL(wp)                              :: maxT(ngrids,2)
    ! out
    REAL(wp)                              :: Tcorr(ngrids,nlev+1)

    Tcorr(1:ngrids,1:nlev)   = M%T
    maxT (1:ngrids,1)        = M%T2M
    maxT (1:ngrids,2)        = M%SKT
    Tcorr(1:ngrids,nlev+1)   = MAXVAL(maxT,DIM=2)

  END FUNCTION CORRECTED_TEMPERATURE_PROFILE

  SUBROUTINE GET_TB(d1,ncol,model,sub,inter,options,frac_out)

    IMPLICIT NONE 

    ! in
    INTEGER,          INTENT(in)   :: ncol,d1
    TYPE(model_type), INTENT(in)   :: model
    TYPE(subset),     INTENT(in)   :: sub
    TYPE(internal),   INTENT(inout):: inter !cflag may change, and Tb
    TYPE(name_list),  INTENT(in)   :: options
    REAL(wp),         INTENT(in)   :: frac_out(ncol,model%aux%nlev)

    ! internal
    INTEGER                        :: nlev

    nlev = model%aux%nlev

    inter%Tb(1:ncol) = TB_ICARUS(d1,ncol,nlev,model,sub,inter,frac_out)
    
    !IF (options%dbg > 1) PRINT *, '--- Finished simulated Tb' 
 
  END SUBROUTINE GET_TB

  FUNCTION TB_ICARUS(d1,ncol,nlev,model,sub,inter,frac_out) RESULT(Tb)

    IMPLICIT NONE

    INTEGER,          INTENT(in) :: d1,ncol,nlev
    TYPE(model_type), INTENT(in) :: model
    TYPE(subset),     INTENT(in) :: sub
    TYPE(internal),   INTENT(in) :: inter
    REAL(wp),         INTENT(in) :: frac_out(ncol,nlev)

    ! out
    REAL(wp), DIMENSION(ncol)    :: Tb

    ! internal
    REAL(wp)                     :: frac_out2(1,ncol,nlev)
    REAL(WP), PARAMETER          :: emsfc_lw = 0.98       ! 10.5 micron emissivity of surface (fraction) 
    REAL(WP),DIMENSION(ncol,nlev):: demIN       ! Subcolumn emissivity
    REAL(WP)                     :: meantbclr(1)! Mean clear-sky 10.5 micron brightness temperature
    REAL(WP),DIMENSION(ncol)     :: &
         boxtau,                    & ! Optical thickness in each column
         boxptop,                   & ! Cloud top pressure (mb) in each column
         boxttop                    ! Cloud top temperature in each column
    INTEGER, DIMENSION(ncol)     :: levmatch ! ngrids,ncols
    REAL(wp), DIMENSION(nlev)    :: Q,T,p_mid
    REAL(wp), DIMENSION(nlev+1)  :: p_int

    boxtau   (1:ncol)        = missing
    boxptop  (1:ncol)        = missing
    boxttop  (1:ncol)        = missing
    meantbclr                = missing
    levmatch (1:ncol)        = -999
    demIN    (1:ncol,1:nlev) = missing

    ! arg1=isccp_top_height, arg2=height_direction
    CALL COSP_ISCCP_INIT(3,2)

    frac_out2(1,1:ncol,1:nlev)=frac_out

    demIN (1:ncol,1:nlev) = &
         MERGE(SPREAD(sub%cloud_emis(d1,1:nlev),1,ncol),&
         0._wp,&
         frac_out2(1,:,:) .GT. 0)

    Q=model%Q(d1,:)
    T=model%T(d1,:)
    p_mid=sub%p_mid(d1,:)
    p_int=sub%p_int(d1,:)

    CALL ICARUS_SUBCOLUMN(1,ncol,nlev,&
         sub%sunlit(d1)   ,&
         inter%tau_profile,&
         demiN            ,&
         model%skt(d1)    ,&
         emsfc_lw         ,&
         Q                ,&
         T                ,&
         p_mid            ,&
         p_int            ,&
         frac_out2        ,&
         levmatch,boxtau,boxptop,boxttop,meantbclr)

    Tb(1:ncol) = boxttop(1:ncol)

  END FUNCTION TB_ICARUS

END MODULE BRIGHTNESS_TEMPERATURES

