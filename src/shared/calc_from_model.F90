MODULE CALC_FROM_MODEL

  ! A suite of routines to derive some model quantities not directly provided
  !
  ! Salomon.Eliasson@smhi.se

  USE COSP_KINDS,                ONLY: WP
  USE COSP_MATH_CONSTANTS,       ONLY: PI
  USE FROM_COSP2,                ONLY: &
       GET_G_NIR,                      &
       GET_SSA_NIR
  USE MODEL_INPUT,               ONLY: MODEL_TYPE
  USE OPTICS_M,                  ONLY: OPTICS_LUT
  USE NAMELIST_INPUT,            ONLY: NAME_LIST
  USE SATELLITE_SPECS,           ONLY: SATELLITE
  USE SIMULATOR_INPUT_VARIABLES, ONLY: SUBSET


  IMPLICIT NONE

  PUBLIC :: calc_model_vertical_properties,&
       GET_MODEL_SSA_AND_G

  REAL(WP), PARAMETER :: g=9.80665, &
       Rd=287.04, &
       Rv=461
CONTAINS

  SUBROUTINE CALC_MODEL_VERTICAL_PROPERTIES(ngrids,nlev,M,S,options)

    IMPLICIT NONE

    INTEGER, INTENT(in)          :: ngrids,nlev
    TYPE(model_type), INTENT(in) :: M
    TYPE(subset), INTENT(inout)  :: S
    TYPE(name_list), INTENT(in)  :: options

    REAL(wp), DIMENSION(ngrids,nlev)  :: hyam_mat,hybm_mat
    REAL(wp), DIMENSION(ngrids,nlev+1):: hyai_mat,hybi_mat
    REAL(wp), DIMENSION(ngrids,nlev)  :: Tv                   ! virtual temperature
    REAL(wp), DIMENSION(ngrids,nlev)  :: CIWC_gm3, CLWC_gm3   ! kg/kg->g/m3
    REAL(wp), DIMENSION(ngrids,nlev)  :: DGEOZ                ! layer width (geopot) [m]
    REAL(wp), DIMENSION(ngrids,nlev)  :: rho                  ! density of dry air [kg/m^3]
    REAL(wp), DIMENSION(ngrids,nlev)  :: zminice,znum,zdesr,zd,zntot,zden,lsm_2D

    REAL(wp), PARAMETER :: rho_w = 1.0       ! density of water [10^3*kg/m^3]
    REAL(wp), PARAMETER :: rho_i = 0.93      ! density of ice   [10^3*kg/m^3]
    REAL(wp), PARAMETER :: replog = 1.e-12   ! security parameter
    REAL(wp), PARAMETER :: RTT = 273.15      ! Celsius to Kelvin
    REAL(wp), PARAMETER :: ZRefDe = 0.64952  ! Ratio reff to diameter = 3*sqrt(3)/8
    REAL(wp), PARAMETER :: a_conc_sea  = 50  ! aerosol concentration over sea [1/cm^3]
    REAL(wp), PARAMETER :: a_conc_land = 900 !              -||-         land [1/cm^3]

    ! fitting parmeters
    REAL(wp), PARAMETER :: a_sea   = -1.15e-3
    REAL(wp), PARAMETER :: b_sea   = 0.963
    REAL(wp), PARAMETER :: c_sea   = 5.3
    REAL(wp), PARAMETER :: zd_sea  = 0.33
    REAL(wp), PARAMETER :: a_land  = -2.1e-4
    REAL(wp), PARAMETER :: b_land  = 0.568
    REAL(wp), PARAMETER :: c_land  = -27.9
    REAL(wp), PARAMETER :: zd_land = 0.43
    REAL(wp), PARAMETER :: zminliq = 4, zmaxliq=16
    REAL(wp), PARAMETER :: zmaxice = 155

    INTEGER :: i

    !==================================================
    ! initialise some variables
    !==================================================
    hyam_mat = SPREAD(M%aux%hyam(1:nlev),  1,ngrids)
    hybm_mat = SPREAD(M%aux%hybm(1:nlev),  1,ngrids)
    hyai_mat = SPREAD(M%aux%hyai(1:nlev+1),1,ngrids)
    hybi_mat = SPREAD(M%aux%hybi(1:nlev+1),1,ngrids)
    lsm_2d   = SPREAD(M%aux%lsm,2,nlev) ! fraction land

    zdesr(1:ngrids,1:nlev) = 0.0
    DGEOZ=0

    ! LAYER GEOPOTENTIAL THICKNESS.
    ! From the hydrostatic equation.
    ! (e.g. ch1.6 Holton 2004; an introduction to dynamical meteorology)
    ! Z=R/g0 *<T>*ln(p1/p2). where <T> is the layer mean temperature

    !==================================================
    !    Virtual temperature
    !    Tv = T ( 1+ (Rv/Rd-1)*qv)
    !==================================================

    Tv(1:ngrids,1:nlev) = M%T(1:ngrids,1:nlev)*&
         (1+ (Rv/Rd-1)*M%Q(1:ngrids,1:nlev))

    !==================================================
    !    Pressure levels
    !    p_mid=hyam+hybm*Psurf
    !    p_int=hyai+hybi*Psurf
    !==================================================

    S%p_mid(1:ngrids,1:nlev) = &
         hyam_mat(1:ngrids,1:nlev) + &
         hybm_mat(1:ngrids,1:nlev)*SPREAD(M%PSURF(1:ngrids),2,nlev)

    S%p_int(1:ngrids,1:nlev+1) = &
         hyai_mat(1:ngrids,1:nlev+1) + &
         hybi_mat(1:ngrids,1:nlev+1)* SPREAD(M%PSURF(1:ngrids),2,nlev+1)

    !==================================================
    !    Density
    !    rho =  p /(Rd * Tv)
    !==================================================

    rho(1:ngrids,1:nlev) = S%p_mid(1:ngrids,1:nlev)/&
         (Rd*Tv(1:ngrids,1:nlev))

    !==================================================
    ! Height from surface
    ! z(i-1)=z(i) + (Rd*Tv)/g * ln(p_int(i)/p_int(i-1))
    !==================================================

    DO i=nlev,2,-1
          DGEOZ(:,i)=(Rd*Tv(:,i))/g*LOG(S%p_int(:,i+1)/S%p_int(:,i))
          S%height(:,i)=S%height(:,i+1) + DGEOZ(:,i)
    END DO
    ! to avoid infinity at the top of the atmosphere...
    DGEOZ(:,1)=(Rd*Tv(:,1))/g*LOG(S%p_int(:,2)/(S%p_int(:,2)/3))
    S%height(:,1)=S%height(:,2) + DGEOZ(:,1)
    !==================================================
    ! Cumulated water vapour
    !==================================================
    IF (ALLOCATED(S%Q_kgm2)) THEN
       S%Q_kgm2(1:ngrids,2:nlev) = &
            M%Q(1:ngrids,2:nlev)*&
            rho(1:ngrids,2:nlev)*&
            DGEOZ(1:ngrids,2:nlev)
    END IF

    !==================================================
    ! convert cwc to cwc_gm3
    ! CLWC_gm3=CLWC/CC*1000*rho
    !==================================================

    CLWC_gm3(1:ngrids,1:nlev) =&
         MERGE(M%CLWC(1:ngrids,1:nlev)/M%CC(1:ngrids,1:nlev)&
         *1e3*rho(1:ngrids,1:nlev),&
         0._wp,&
         M%CC(1:ngrids,1:nlev)>0)

    CIWC_gm3(1:ngrids,1:nlev) =&
         MERGE(M%CIWC(1:ngrids,1:nlev)/M%CC(1:ngrids,1:nlev)&
         *1e3*rho(1:ngrids,1:nlev),&
         0._wp,&
         M%CC(1:ngrids,1:nlev)>0)

    !==================================================
    ! Effective radius
    !==================================================

    ! Function that calculates effective radius according
    ! to EC-Earth parameterization
    !
    ! J.F. Meirink, 2014-07-09
    !
    ! Parameterizations for:
    ! - liquid clouds: Martin et al., JAS, 1994.
    ! - ice clouds:
    !  f(T,IWC): Sun and Rikus (1999) revised by Sun (2001)

    ! -----------------------
    !   LIQUID CLOUDS
    ! -----------------------

    ! Get a linear combination of the land and sea parametrisations
    ! Equations are from Martin et. al., 1995

    ! zntot: exactly equation 12 & 13 in paper
    ! lreff[micron] = (3*lwc_gm3 / (4*pi*rho_w*k*Ntot) )^ 1/3 (eq 11), but then revised
    ! where L=LWC
  
    WHERE( CLWC_gm3(1:ngrids,1:nlev) > 0) 
       zntot(1:ngrids,1:nlev) =                           &
            (1-lsm_2d(1:ngrids,1:nlev))*                  &
            (a_sea*a_conc_sea**2 + b_sea*a_conc_sea + c_sea) + &
            lsm_2d(1:ngrids,1:nlev)*                      &
            (a_land*a_conc_land**2 + b_land*a_conc_land + c_land)

       zd(1:ngrids,1:nlev) = zd_sea*(1-lsm_2d(1:ngrids,1:nlev))+&
            zd_land*lsm_2d(1:ngrids,1:nlev)

       znum(1:ngrids,1:nlev) = 3.0*CLWC_gm3(1:ngrids,1:nlev)*&
            (1.0_wp+3.0_wp*zd(1:ngrids,1:nlev)*zd(1:ngrids,1:nlev))**2

       zden(1:ngrids,1:nlev) = 4.0*pi*zntot(1:ngrids,1:nlev)*&
            (1.0 + zd(1:ngrids,1:nlev)*zd(1:ngrids,1:nlev))**3

       ! if znum/zden<replog then lreff=4, otherwise 4<=lreff<=16
       S%lreff(1:ngrids,1:nlev) = &
            MERGE( &
            MIN( &
            MAX(100._wp*EXP(0.333_wp*&
            LOG(znum(1:ngrids,1:nlev)/zden(1:ngrids,1:nlev))),&
            zminliq),& !max
            zmaxliq),& !min
            zminliq, & !merge=false
            znum(1:ngrids,1:nlev)/zden(1:ngrids,1:nlev) > replog) ! condition

    END WHERE

    ! -----------------------
    !   ICE CLOUDS
    ! -----------------------
    ! Is ultimately only a function of ciwc, temperature, and latitude

    zminice(1:ngrids,1:nlev)   = SPREAD(&
         20.0 + 40.0*COS(M%aux%lat_v(1:ngrids)*pi/180.0 ),2,nlev)

    WHERE ( CIWC_gm3(1:ngrids,1:nlev) > 0)

       zdesr(1:ngrids,1:nlev) = &
            (1.2351 + 0.0105*(M%T(1:ngrids,1:nlev) - RTT))*&
            ( (45.8966*CIWC_gm3(1:ngrids,1:nlev)**0.2214) +       &
            (0.7957*CIWC_gm3(1:ngrids,1:nlev)**0.2535) *          &
            (M%T(1:ngrids,1:nlev) - 83.15))

       ! ZRefDe* ( [20 60] <= ireff <= 155) , depending on latitude
       S%ireff(1:ngrids,1:nlev) = ZRefDe*&
            MIN( &
            MAX(zdesr(1:ngrids,1:nlev),&
            zminice(1:ngrids,1:nlev)),&
            zmaxice )

    END WHERE

    !==================================================
    ! Optical depth
    ! iwp = iwc*dz
    ! itau = 1e3*   3*iwp/(2*rho_i*ireff)
    ! emis=1-EXP(tau/x), (x_i=2.13 and x_l=2.56 (Rossow et al. 1996))
    !==================================================

    S%lwp(1:ngrids,1:nlev) = &
         CLWC_gm3(1:ngrids,1:nlev)*DGEOZ(1:ngrids,1:nlev)/1.e3 ![kg/m2] 
    S%ltau(1:ngrids,1:nlev) = &
         (1e3*S%lwp(1:ngrids,1:nlev)*3)/(2*rho_w*S%lreff(1:ngrids,1:nlev))
    S%iwp(1:ngrids,1:nlev) = &
         CIWC_gm3(1:ngrids,1:nlev)*DGEOZ(1:ngrids,1:nlev)/1.e3 ![kg/m2] 
    S%itau(1:ngrids,1:nlev) = &
         (1e3*S%iwp(1:ngrids,1:nlev)*3)/(2*rho_i*S%ireff(1:ngrids,1:nlev))

    S%tau=S%ltau+S%itau

    ! cloud_emis
    IF (ALLOCATED(S%cloud_emis)) THEN

       S%cloud_emis(1:ngrids,1:nlev) = 1. - &
            EXP( -(S%ltau(1:ngrids,1:nlev)/2.56+S%itau(1:ngrids,1:nlev)/2.13) )
    END IF

    IF ( (options%sim%doRTTOV).OR.&
         (options%sim%doClara.AND.options%sim%Tb .EQ. 1)) THEN
       ! ---------------------------------------
       ! Also need the layer averages for RTTOV
       ! ---------------------------------------
       S%LWC(1:ngrids,1:nlev-1) = &
            (CLWC_gm3(1:ngrids,1:nlev-1)+CLWC_gm3(1:ngrids,2:nlev))/2
       S%IWC(1:ngrids,1:nlev-1) = &
            (CIWC_gm3(1:ngrids,1:nlev-1)+CIWC_gm3(1:ngrids,2:nlev))/2
    END IF
  END SUBROUTINE CALC_MODEL_VERTICAL_PROPERTIES

  SUBROUTINE GET_MODEL_SSA_AND_G(sub,sat,ngrids,nlev,r_eff)

    !Calculate the combined liquid and ice single
    !scattering albedo, g, and the asymmetry parameter, w0 for the model. This is
    !done using look up tables of g and w0 as a function cloud
    !effective radius. At each model layer the g and w0 are derived
    !for ice and liquid separately and are then combined
    ! This is based on code from COSP (cosp_optics.F90)

    IMPLICIT NONE

    TYPE(satellite), INTENT(in)  :: sat
    TYPE(optics_LUT), INTENT(in) :: r_eff
    TYPE(subset), INTENT(inout)  :: sub
    INTEGER, INTENT(in)          :: ngrids,nlev

    REAL(wp), DIMENSION(ngrids,nlev) :: water_g,ice_g,water_w0,ice_w0

    water_g (1:ngrids,1:nlev) = 0.0   
    water_w0(1:ngrids,1:nlev) = 0.0
    ice_g   (1:ngrids,1:nlev) = 0.0
    ice_w0  (1:ngrids,1:nlev) = 0.0

    WHERE (.NOT.sub%data_mask .AND. SPREAD(sub%sunlit,2,nlev).EQ.1 )
       WHERE (sub%tau > 0) 
          ! there is a cloud
          WHERE (sub%ltau(1:ngrids,1:nlev) > 0)
             !there is a liquid cloud
             water_g (1:ngrids,1:nlev) = &
                  GET_G_NIR  (sub%lreff(1:ngrids,1:nlev),r_eff%water,sat%is_ch3b)
             water_w0(1:ngrids,1:nlev) = &
                  GET_SSA_NIR(sub%lreff(1:ngrids,1:nlev),r_eff%water,sat%is_ch3b)
          END WHERE
          WHERE (sub%itau(1:ngrids,1:nlev) > 0)
             ! there is an ice cloud
             ice_g   (1:ngrids,1:nlev) = &
                  GET_G_NIR  (sub%ireff(1:ngrids,1:nlev),r_eff%ice,sat%is_ch3b)
             ice_w0  (1:ngrids,1:nlev) = &
                  GET_SSA_NIR(sub%ireff(1:ngrids,1:nlev),r_eff%ice,sat%is_ch3b)
          END WHERE

          sub%g0 (1:ngrids,1:nlev) = &
               (sub%ltau(1:ngrids,1:nlev)*water_g(1:ngrids,1:nlev) + &
               sub%itau(1:ngrids,1:nlev)*ice_g(1:ngrids,1:nlev)) / & 
               sub%tau 
          sub%w0(1:ngrids,1:nlev) = ( &
               sub%ltau(1:ngrids,1:nlev)*&
               water_g(1:ngrids,1:nlev)*water_w0(1:ngrids,1:nlev) + &
               sub%itau(1:ngrids,1:nlev)*&
               ice_g(1:ngrids,1:nlev)*ice_w0(1:ngrids,1:nlev) ) / &
               (sub%g0(1:ngrids,1:nlev)*sub%tau)
       END WHERE
    END WHERE

  END SUBROUTINE GET_MODEL_SSA_AND_G

END MODULE CALC_FROM_MODEL
