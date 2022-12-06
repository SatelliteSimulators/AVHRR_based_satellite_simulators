MODULE CALC_FROM_MODEL

  ! A suite of routines to derive some model quantities not directly provided
  !
  ! Salomon.Eliasson@smhi.se

  USE COSP_KINDS,                ONLY: WP
  USE COSP_MATH_CONSTANTS,       ONLY: PI
  USE FROM_COSP2,                ONLY: &
       GET_G_NIR,                      &
       GET_SSA_NIR,                    &
       phaseIsIce, phaseIsLiquid
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

  SUBROUTINE CALC_MODEL_VERTICAL_PROPERTIES(ng,nl,M,S,options)

    IMPLICIT NONE

    INTEGER, INTENT(in)          :: ng,nl
    TYPE(model_type), INTENT(in) :: M
    TYPE(subset), INTENT(inout)  :: S
    TYPE(name_list), INTENT(in)  :: options

    REAL(wp), DIMENSION(ng,nl)  :: hyam_mat,hybm_mat
    REAL(wp), DIMENSION(ng,nl+1):: hyai_mat,hybi_mat
    REAL(wp), DIMENSION(ng,nl)  :: Tv                   ! virtual temperature
    REAL(wp), DIMENSION(ng,nl)  :: CIWC_gm3, CLWC_gm3   ! kg/kg->g/m3
    REAL(wp), DIMENSION(ng,nl)  :: DGEOZ                ! layer width (geopot) [m]
    REAL(wp), DIMENSION(ng,nl)  :: rho                  ! density of dry air [kg/m^3]
    REAL(wp), DIMENSION(ng,nl)  :: zminice,znum,zdesr,zd,zntot,zden,lsm_2D

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

    PRINT *, "--- Calculating model vertical properties"
    
    !==================================================
    ! initialise some variables
    !==================================================
    hyam_mat = SPREAD(M%aux%hyam(1:nl),  1,ng)
    hybm_mat = SPREAD(M%aux%hybm(1:nl),  1,ng)
    hyai_mat = SPREAD(M%aux%hyai(1:nl+1),1,ng)
    hybi_mat = SPREAD(M%aux%hybi(1:nl+1),1,ng)
    lsm_2d   = SPREAD(M%aux%lsm,2,nl) ! fraction land

    zdesr(1:ng,1:nl) = 0.0
    DGEOZ=0

    ! LAYER GEOPOTENTIAL THICKNESS.
    ! From the hydrostatic equation.
    ! (e.g. ch1.6 Holton 2004; an introduction to dynamical meteorology)
    ! Z=R/g0 *<T>*ln(p1/p2). where <T> is the layer mean temperature

    !==================================================
    !    Virtual temperature
    !    Tv = T ( 1+ (Rv/Rd-1)*qv)
    !==================================================

    Tv(1:ng,1:nl) = M%T(1:ng,1:nl)*&
         (1+ (Rv/Rd-1)*M%Q(1:ng,1:nl))

    !==================================================
    !    Pressure levels
    !    p_mid=hyam+hybm*Psurf
    !    p_int=hyai+hybi*Psurf
    !==================================================

    S%p_mid(1:ng,1:nl) = &
         hyam_mat(1:ng,1:nl) + &
         hybm_mat(1:ng,1:nl)*SPREAD(M%PSURF(1:ng),2,nl)

    S%p_int(1:ng,1:nl+1) = &
         hyai_mat(1:ng,1:nl+1) + &
         hybi_mat(1:ng,1:nl+1)* SPREAD(M%PSURF(1:ng),2,nl+1)

    !==================================================
    !    Density
    !    rho =  p /(Rd * Tv)
    !==================================================

    rho(1:ng,1:nl) = S%p_mid(1:ng,1:nl)/&
         (Rd*Tv(1:ng,1:nl))

    !==================================================
    ! Height from surface
    ! z(i-1)=z(i) + (Rd*Tv)/g * ln(p_int(i)/p_int(i-1))
    !==================================================

    DO i=nl,2,-1
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
       S%Q_kgm2(1:ng,2:nl) = &
            M%Q(1:ng,2:nl)*&
            rho(1:ng,2:nl)*&
            DGEOZ(1:ng,2:nl)
    END IF

    !==================================================
    ! convert cwc to cwc_gm3
    ! CLWC_gm3=CLWC/CC*1000*rho
    !==================================================

    CLWC_gm3(1:ng,1:nl) =&
         MERGE(M%CLWC(1:ng,1:nl)/M%CC(1:ng,1:nl)&
         *1e3*rho(1:ng,1:nl),&
         0._wp,&
         M%CC(1:ng,1:nl)>0)

    CIWC_gm3(1:ng,1:nl) =&
         MERGE(M%CIWC(1:ng,1:nl)/M%CC(1:ng,1:nl)&
         *1e3*rho(1:ng,1:nl),&
         0._wp,&
         M%CC(1:ng,1:nl)>0)

    !==================================================
    ! Effective radius
    !==================================================

    ! Function that calculates effective radius according
    ! to EC-Earth parameterization
    !
    ! J.F. Meirink, 2014-07-09
    !
    ! Parameterizations for: 
    ! - liquid clouds: Martin et al., JAS,
    !
    ! 1994.  - ice clouds: f(T,IWC): Sun and Rikus (1999) revised by
    ! Sun (2001). Note: The latitude dependence is based on the observations
    ! that tropical clouds typically have larger particles than extra
    ! tropical clouds, e.g., Field et. al., 2007, JAS

    ! -----------------------
    !   LIQUID CLOUDS
    ! -----------------------

    ! Get a linear combination of the land and sea parametrisations
    ! Equations are from Martin et. al., 1995

    ! zntot: exactly equation 12 & 13 in paper
    ! lreff[micron] = (3*lwc_gm3 / (4*pi*rho_w*k*Ntot) )^ 1/3 (eq 11), but then revised
    ! where L=LWC
  
    WHERE( CLWC_gm3(1:ng,1:nl) > 0) 
       zntot(1:ng,1:nl) =                           &
            (1-lsm_2d(1:ng,1:nl))*                  &
            (a_sea*a_conc_sea**2 + b_sea*a_conc_sea + c_sea) + &
            lsm_2d(1:ng,1:nl)*                      &
            (a_land*a_conc_land**2 + b_land*a_conc_land + c_land)

       zd(1:ng,1:nl) = zd_sea*(1-lsm_2d(1:ng,1:nl))+&
            zd_land*lsm_2d(1:ng,1:nl)

       znum(1:ng,1:nl) = 3.0*CLWC_gm3(1:ng,1:nl)*&
            (1.0_wp+3.0_wp*zd(1:ng,1:nl)*zd(1:ng,1:nl))**2

       zden(1:ng,1:nl) = 4.0*pi*zntot(1:ng,1:nl)*&
            (1.0 + zd(1:ng,1:nl)*zd(1:ng,1:nl))**3

       ! if znum/zden<replog then lreff=4, otherwise 4<=lreff<=16
       S%lreff(1:ng,1:nl) = &
            MERGE( &
            MIN( &
            MAX(100._wp*EXP(0.333_wp*&
            LOG(znum(1:ng,1:nl)/zden(1:ng,1:nl))),&
            zminliq),& !max
            zmaxliq),& !min
            zminliq, & !merge=false
            znum(1:ng,1:nl)/zden(1:ng,1:nl) > replog) ! condition

    END WHERE

    ! -----------------------
    !   ICE CLOUDS
    ! -----------------------
    ! Is ultimately only a function of ciwc, temperature, and latitude

    zminice(1:ng,1:nl)   = SPREAD(&
         20.0 + 40.0*COS(M%aux%lat_v(1:ng)*pi/180.0 ),2,nl)

    WHERE ( CIWC_gm3(1:ng,1:nl) > 0)

       zdesr(1:ng,1:nl) = &
            (1.2351 + 0.0105*(M%T(1:ng,1:nl) - RTT))*&
            ( (45.8966*CIWC_gm3(1:ng,1:nl)**0.2214) +       &
            (0.7957*CIWC_gm3(1:ng,1:nl)**0.2535) *          &
            (M%T(1:ng,1:nl) - 83.15))

       ! ZRefDe* ( [20 60] <= ireff <= 155) , depending on latitude
       S%ireff(1:ng,1:nl) = ZRefDe*&
            MIN( &
            MAX(zdesr(1:ng,1:nl),&
            zminice(1:ng,1:nl)),&
            zmaxice )

    END WHERE

    !==================================================
    ! Optical depth
    ! iwp = iwc*dz
    ! itau = 1e3*   3*iwp/(2*rho_i*ireff)
    ! emis=1-EXP(tau/x), (x_i=2.13 and x_l=2.56 (Rossow et al. 1996))
    !==================================================

    S%lwp(1:ng,1:nl) = &
         CLWC_gm3(1:ng,1:nl)*DGEOZ(1:ng,1:nl)/1.e3 ![kg/m2] 
    S%ltau(1:ng,1:nl) = &
         (1e3*S%lwp(1:ng,1:nl)*3)/(2*rho_w*S%lreff(1:ng,1:nl))
    S%iwp(1:ng,1:nl) = &
         CIWC_gm3(1:ng,1:nl)*DGEOZ(1:ng,1:nl)/1.e3 ![kg/m2] 
    S%itau(1:ng,1:nl) = &
         (1e3*S%iwp(1:ng,1:nl)*3)/(2*rho_i*S%ireff(1:ng,1:nl))

    S%tau=S%ltau+S%itau

    ! cloud_emis
    IF (ALLOCATED(S%cloud_emis)) THEN

       S%cloud_emis(1:ng,1:nl) = 1. - &
            EXP( -(S%ltau(1:ng,1:nl)/2.56+S%itau(1:ng,1:nl)/2.13) )
    END IF

    IF (options%sim%doRTTOV) THEN
       ! ---------------------------------------
       ! Also need the layer averages for RTTOV
       ! ---------------------------------------
       S%LWC(1:ng,1:nl-1) = &
            (CLWC_gm3(1:ng,1:nl-1)+CLWC_gm3(1:ng,2:nl))/2
       S%IWC(1:ng,1:nl-1) = &
            (CIWC_gm3(1:ng,1:nl-1)+CIWC_gm3(1:ng,2:nl))/2
    END IF
    
  END SUBROUTINE CALC_MODEL_VERTICAL_PROPERTIES

  SUBROUTINE GET_MODEL_SSA_AND_G(sub,LUT,ng,nl)

    !Calculate the combined liquid and ice single scattering albedo,
    !g, and the asymmetry parameter, w0 for the model. This is done
    !using look up tables of g and w0 as a function cloud effective
    !radius. At each model layer the g and w0 are derived for ice and
    !liquid separately and are then combined This is based on code
    !from COSP (cosp_optics.F90)

    IMPLICIT NONE

    TYPE(subset),     INTENT(inout) :: sub
    TYPE(optics_LUT), INTENT(in)    :: LUT
    INTEGER,          INTENT(in)    :: ng,nl

    ! internal
    REAL(wp), DIMENSION(ng,nl) :: water_g,ice_g,water_w0,ice_w0
    integer :: i,j !debugging
    
    PRINT *, "--- Calculating model cloud optical properties"

    water_g (1:ng,1:nl) = 0.0   
    water_w0(1:ng,1:nl) = 0.0
    ice_g   (1:ng,1:nl) = 0.0
    ice_w0  (1:ng,1:nl) = 0.0

!    do i=1,ng
!       do j =1,nl
!    do i=57303,57303
!       do j =41,nl
!    if (sub%tau(i,j) > 0) then
!             ! there is a cloud
!             if (sub%ltau(i,j) > 0) then
!                ! there is a liquid cloud
!
!                water_g (i,j) = &
!                     GET_G_NIR  (sub%lreff(i,j),LUT%water,phaseIsLiquid)
!                water_w0(1:ng,1:nl) = &
!                     GET_SSA_NIR(sub%lreff(i,j),LUT%water,phaseIsLiquid)
!
!
!             end if
!             if (sub%itau(i,j) > 0) then
!                ! there is an ice cloud
!
!                ice_g   (i,j) = &
!                     GET_G_NIR  (sub%ireff(i,j),LUT%ice,phaseIsIce)
!                ice_w0  (i,j) = &
!                     GET_SSA_NIR(sub%ireff(i,j),LUT%ice,phaseIsIce)
!                
!             end if
!
!             ! combined cloud scattering albedo
!             sub%g0 (i,j) = &
!                  (sub%ltau(i,j)*water_g(i,j) + sub%itau(i,j)*ice_g(i,j)) / & 
!                  sub%tau(i,j)
!
!          end if
!       end do
!    end do

    WHERE (.NOT.sub%data_mask .AND. SPREAD(sub%sunlit,2,nl).EQ.1 )
       WHERE (sub%tau(1:ng,1:nl) > 0) 
          ! there is a cloud

          water_g (1:ng,1:nl) = &
               GET_G_NIR  (sub%lreff(1:ng,1:nl),LUT%water,phaseIsLiquid)
          water_w0(1:ng,1:nl) = &
               GET_SSA_NIR(sub%lreff(1:ng,1:nl),LUT%water,phaseIsLiquid)
          ice_g   (1:ng,1:nl) = &
               GET_G_NIR  (sub%ireff(1:ng,1:nl),LUT%ice,phaseIsIce)
          ice_w0  (1:ng,1:nl) = &
                  GET_SSA_NIR(sub%ireff(1:ng,1:nl),LUT%ice,phaseIsIce)

          sub%w0 (1:ng,1:nl) = &
               (sub%ltau(1:ng,1:nl)*water_w0(1:ng,1:nl) + &
               sub%itau(1:ng,1:nl)*ice_w0(1:ng,1:nl)) / & 
               sub%tau(1:ng,1:nl)
          
          sub%g0(1:ng,1:nl) = ( &
               sub%ltau(1:ng,1:nl)*water_g(1:ng,1:nl)*water_w0(1:ng,1:nl) + &
               sub%itau(1:ng,1:nl)*ice_g(1:ng,1:nl)*ice_w0(1:ng,1:nl) ) / &
               (sub%w0(1:ng,1:nl)*sub%tau(1:ng,1:nl))
       END WHERE
    END WHERE

!    write(*,'(a5,x,i5,a,x,i2)') "grd =",57303,"top = 41 bottom =",nl
!    write(*,'(a6,x,51(f7.2,x))') "ltau =",sub%ltau(57303,41:nl)
!    write(*,'(a6,x,51(f7.2,x))') "lref =",sub%lreff(57303,41:nl)
!    write(*,'(a6,x,51(f7.2,x))') "g0_w =",water_g(57303,41:nl)
!    write(*,'(a6,x,51(f7.2,x))') "w0_w =",water_w0(57303,41:nl)
!    write(*,'(a6,x,51(f7.2,x))') "itau =",sub%itau(57303,41:nl)
!    write(*,'(a6,x,51(f7.2,x))') "iref =",sub%ireff(57303,41:nl)
!    write(*,'(a6,x,51(f7.2,x))') "g0_i =",ice_g(57303,41:nl)
!    write(*,'(a6,x,51(f7.2,x))') "w0_i =",ice_w0(57303,41:nl)
!    write(*,'(a6,x,51(f7.2,x))') "tau =",sub%tau(57303,41:nl)
!    write(*,'(a6,x,51(f7.2,x))') "g0 =",sub%g0(57303,41:nl)
!    write(*,'(a6,x,51(f7.2,x))') "w0 =",sub%w0(57303,41:nl)
  END SUBROUTINE GET_MODEL_SSA_AND_G

END MODULE CALC_FROM_MODEL
