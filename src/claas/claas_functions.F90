MODULE CLAAS_FUNCTIONS
  ! Simulate CLAAS properties
  !
  ! Salomon.Eliasson@smhi.se
   USE CLAAS_M,                   ONLY: &
      CLAAS_FIELDS,                   &
      CLAAS_TYPE,                     &
      INITIALISE_CLAAS_SIM
   USE COSP_KINDS,                ONLY: WP
   USE DATA_CHECK,                ONLY: &
      CLOUD_LOW, CLOUD_MID,           &
      CHECK_VARIABLES,                &
      CHECK_VAR,                      &
      COMMON_CHECK_SUBGRID,           &
      CTH_MAX,CTH_MIN,CTP_MAX,CTP_MIN,&
      CTT_MAX,CTT_MIN,EPSR,           &
      FRAC_MIN,FRAC_MAX,MISSING,      &
      P_TAU_MIN, P_TAU_MAX,           &
      TAU_MIN,TAU_MAX,T_MAX,T_MIN,    &
      PRINT_SUB_POINT ! FOR DEBUGGING
   USE INTERNAL_SIMULATOR,        ONLY: INTERNAL
   USE MODEL_INPUT,               ONLY: MODEL_TYPE
   USE NAMELIST_INPUT,            ONLY: NAME_LIST
   USE SIMULATOR_INPUT_VARIABLES, ONLY: SUBSET
   USE SIMPLE_CTTH_SIMULATOR,     ONLY: CTP_SIMPLE, CTH_SIMPLE

   IMPLICIT NONE

   PUBLIC :: CHECK_GRID_AVERAGES,&
      CHECK_PROFILES,          &
      CHECK_VARIABLES_SUBGRID, &
      CTTH,                    &
      GRID_AVERAGE,            &
      DAY_ADD,                 &
      DAY_AVERAGE
   PRIVATE :: CLAAS_CTTH_OPAQUE, &
      CLAAS_CTTH_ST_SIMPLE

   REAL(wp),PARAMETER :: ref_liq_min = 3,ref_liq_max = 34,&
      lwp_min = ref_liq_min*tau_min,&
      lwp_max = ref_liq_max*tau_max
   REAL(wp),PARAMETER :: ref_ice_min = 5,ref_ice_max = 80,&
      iwp_min = ref_ice_min*tau_min,&
      iwp_max = ref_ice_max*tau_max

   ! these are the ISSCP limits (http://isccp.giss.nasa.gov/cloudtypes.html)
   REAL(wp), PARAMETER :: lw_cl_lim  = 68000  ![Pa]
   REAL(wp), PARAMETER :: hi_cl_lim  = 44000  ![Pa]

   CONTAINS

SUBROUTINE CTTH(d1,ins,CDR,model,S,inter)
   ! The PPS CTTH algorithms
   !
   ! For thick clouds, use a simple approach
   ! For semi-transparent clouds use a special approach
   !
   IMPLICIT NONE

   TYPE(model_type), intent(in) :: model
   TYPE(subset), INTENT(in)     :: S
   TYPE(internal), INTENT(inout):: inter
   CHARACTER(*),INTENT(in)      :: CDR
   INTEGER, INTENT(in)          :: d1,ins
   INTEGER :: nlev

   ! debug
   REAL(wp) :: cth
   cth=-999._wp
   nlev = model%aux%nlev
   ! debug

   ! A very conservative limit of an optical thickness limit of 10
   ! for clouds at 3.7 micron  (Cattani et. al., 2007)
   IF (inter%cflag(ins) .EQ. 3) THEN
      ! OPAQUE CLOUDS
      CALL CLAAS_CTTH_OPAQUE(d1,ins,model,S,inter)
   ELSE
      CALL CLAAS_CTTH_ST(d1,ins,model,S,inter)
   END IF

END SUBROUTINE CTTH

SUBROUTINE CLAAS_CTTH_ST(d1,ins,model,S,inter)

   ! Simple approach using only the documented biases for
   ! semi-transparent clouds.
   !
   ! Salomon.Eliasson@smhi.se

   IMPLICIT NONE

   INTEGER, INTENT(in)          :: d1,ins
   TYPE(subset), INTENT(in)     :: S
   TYPE(model_type), INTENT(in) :: model
   TYPE(internal), intent(inout):: inter

   ! low,middle,high
   REAL(wp), ALLOCATABLE :: PPS_offset(:)
   REAL(wp) :: height,ctp,cth
   INTEGER :: inl,nlev

   height = 0._wp; ctp = missing; cth = missing
   nlev = model%aux%nlev

   ALLOCATE(PPS_offset(3))
   PPS_offset(1:3) = (/ 244, 266, -554/) ![m]

   ! Base it on the CTTH from the MODIS approach
   ctp = CTP_SIMPLE(d1,ins,nlev,S,inter,1._wp)
   cth = CTH_SIMPLE(d1,ins,nlev,S,inter,1._wp)
   IF (ctp .GT. lw_cl_lim) THEN
      ! low cloud
      height = cth + PPS_offset(1)
   ELSEIF ((ctp.LE.lw_cl_lim).AND.(ctp.GT.hi_cl_lim)) THEN
      ! middle cloud
      height = cth + PPS_offset(2)
   ELSE
      ! it is a high cloud
      height = cth + PPS_offset(3)
   END IF
   inl = nlev
   DO WHILE ((S%height(d1,inl)+S%height(d1,inl+1))/2.LT.height)
      inl = inl-1
   END DO
   inter%cth(ins) = MAXVAL([height,0._wp]) ! sometimes below the surface
   inter%ctp(ins) = S%p_mid(d1,inl)
   inter%ctt(ins) = S%Tcorr(d1,inl)

   DEALLOCATE(PPS_offset)
END SUBROUTINE CLAAS_CTTH_ST

SUBROUTINE CLAAS_CTTH_OPAQUE(d1,ins,model,S,inter)

   ! Simple approach using only the documented biases for
   ! opaque clouds for CLAAS-3

   IMPLICIT NONE

   INTEGER, INTENT(in)          :: d1,ins
   TYPE(subset), INTENT(in)     :: S
   TYPE(model_type), INTENT(in) :: model
   TYPE(internal), intent(inout):: inter

   REAL(wp), ALLOCATABLE :: PPS_offset(:)
   REAL(wp) :: height,ctp,cth
   INTEGER :: inl,nlev

   height = 0._wp; ctp = missing; cth = missing
   nlev = model%aux%nlev
   ALLOCATE(PPS_offset(3))
   ! --------------------
   ! Based on table 6, Karlsson et al., 2023
   ! CLARA-A3 compared with CALIPSO-CALIOP for exclusively
   ! single-layered clouds with cloud optical thicknesses above 0.4.
   ! For now, we use the same threshold offsets also for CLAAS-3

   PPS_offset(1:3) = (/ 177, 290, 212 /) ![m]

   ! Base it on the CTTH from the MODIS approach
   ctp = CTP_SIMPLE(d1,ins,nlev,S,inter,1._wp)
   cth = CTH_SIMPLE(d1,ins,nlev,S,inter,1._wp)
   IF (ctp .GT. lw_cl_lim) THEN
      ! low cloud
      height = cth + PPS_offset(1)
   ELSEIF ((ctp.LE.lw_cl_lim).AND.(ctp.GT.hi_cl_lim)) THEN
      ! middle cloud
      height = cth + PPS_offset(2)
   ELSE
      ! it is a high cloud
      height = cth + PPS_offset(3)
   END IF
   inl = nlev
   DO WHILE ((S%height(d1,inl)+S%height(d1,inl+1))/2.LT.height)
      inl = inl-1
   END DO
   inter%cth(ins) = MAXVAL([height,0._wp]) ! sometimes below the surface
   inter%ctp(ins) = S%p_mid(d1,inl)
   inter%ctt(ins) = S%Tcorr(d1,inl)

   DEALLOCATE(PPS_offset)
END SUBROUTINE CLAAS_CTTH_OPAQUE

  SUBROUTINE GRID_AVERAGE(d1,S,inter,ncol,IN)

    ! Calculate simulator averages on the model grid scale.

    IMPLICIT NONE

    INTEGER, INTENT(in)               :: d1,ncol
    TYPE(subset), INTENT(in)          :: S
    TYPE(internal), INTENT(in)        :: inter
    TYPE(CLAAS_fields), INTENT(inout) :: IN

    ! internal
    INTEGER :: ncld,nliq,nice
    LOGICAL,DIMENSION(ncol) :: iscloud,isliquid,isice,islow,ismid,ishigh
    REAL(wp):: tmpCTP_log(ncol)
    REAL(wp):: re_fill=-999

    iscloud  (1:ncol) = inter%cflag(1:ncol).GE.2 ! i.e., tau>tau_min
    isliquid (1:ncol) = (inter%cph (1:ncol).EQ.1).AND.iscloud(1:ncol)
    isice    (1:ncol) = (inter%cph (1:ncol).EQ.2).AND.iscloud(1:ncol)

    islow  = inter%ctp.GE.cloud_low
    ismid  = inter%ctp.LT.cloud_low.AND.inter%ctp.GE.cloud_mid
    ishigh = inter%ctp.LT.cloud_mid.AND.iscloud

    ncld = COUNT(iscloud)
    nliq = COUNT(isliquid)
    nice = COUNT(isice)

    ! initialise values that are defined even when no cloud is present
    IN%cot_ice(d1) = MERGE(0._wp,re_fill,S%sunlit(d1).EQ.1)
    IN%iwp    (d1) = MERGE(0._wp,re_fill,S%sunlit(d1).EQ.1)
    IN%cot_liq(d1) = MERGE(0._wp,re_fill,S%sunlit(d1).EQ.1)
    IN%lwp    (d1) = MERGE(0._wp,re_fill,S%sunlit(d1).EQ.1)
    IN%tau    (d1) = MERGE(0._wp,re_fill,S%sunlit(d1).EQ.1)

    tmpCTP_log(1:ncol) = 0._wp

    IF (ncld .GT. 0) THEN

       ! ----
       ! variables only defined if a cloud is present
       !
       WHERE (iscloud)
          tmpCTP_log  = LOG(inter%ctp)
       END WHERE
       IN%ctp_log(d1) = EXP(SUM(tmpCTP_log)  /ncld)
       IN%ctp    (d1) = SUM(inter%ctp,MASK=iscloud) /ncld
       IN%cth    (d1) = SUM(inter%cth,MASK=iscloud) /ncld
       IN%ctt    (d1) = SUM(inter%ctt,MASK=iscloud) /ncld

       IF (S%sunlit(d1).EQ.1) THEN
          ! Microphysical simulations
          ! -- require sunlight
          IF (nliq .GT. 0) &
               IN%ref_liq(d1)= SUM(inter%reff,MASK=isliquid)/nliq
          IF (nice .GT. 0) &
               IN%ref_ice(d1)= SUM(inter%reff,MASK=isice)   /nice

          !
          !--------

          ! ---------
          ! Variables that are always defined during sunlight
          !
          IN%cfc_day(d1) = (1._wp * ncld)/ncol
          IN%cot_ice(d1) = SUM(inter%tau,MASK=isice)   /ncol
          IN%iwp    (d1) = SUM(inter%cwp,MASK=isice)   /ncol
          IN%lwp    (d1) = SUM(inter%cwp,MASK=isliquid)/ncol
          IN%cot_liq(d1) = SUM(inter%tau,MASK=isliquid)/ncol
          IN%tau    (d1) = SUM(inter%tau,MASK=iscloud) /ncol
       END IF
    END IF

    ! Grid averaged quantities
    IN%cfc  (d1) = (1._wp * ncld)/ncol

    ! cloud height categories
    IN%cfc_low (d1) = COUNT(islow) /REAL(ncol,wp)
    IN%cfc_mid (d1) = COUNT(ismid) /REAL(ncol,wp)
    IN%cfc_high(d1) = COUNT(ishigh)/REAL(ncol,wp)

  END SUBROUTINE GRID_AVERAGE

  SUBROUTINE DAY_ADD(IN,S,ngrids,n_pbins,n_tbins)

    IMPLICIT NONE

    TYPE(CLAAS_type), INTENT(inout) :: IN
    TYPE(subset), INTENT(in)        :: S
    INTEGER, INTENT(in)             :: ngrids,n_tbins,n_pbins

    REAL(wp) :: tmpctp_log(ngrids)

    WHERE(.NOT.S%data_mask(:,1))

       ! Get the sum
       tmpctp_log     = MERGE(IN%av%ctp_log,0.1_wp,IN%av%ctp_log>0)
       IN%sum%cfc     = MERGE(IN%sum%cfc     +IN%av%cfc,      IN%sum%cfc,     IN%av%cfc   >0)
       IN%sum%cfc_day = MERGE(IN%sum%cfc_day +IN%av%cfc_day,  IN%sum%cfc_day, IN%av%cfc_day >0)
       IN%sum%cfc_low = MERGE(IN%sum%cfc_low +IN%av%cfc_low,  IN%sum%cfc_low, IN%av%cfc_low >0)
       IN%sum%cfc_mid = MERGE(IN%sum%cfc_mid +IN%av%cfc_mid,  IN%sum%cfc_mid, IN%av%cfc_mid >0)
       IN%sum%cfc_high= MERGE(IN%sum%cfc_high+IN%av%cfc_high, IN%sum%cfc_high,IN%av%cfc_high>0)
       IN%sum%cth     = MERGE(IN%sum%cth     +IN%av%cth,      IN%sum%cth,     IN%av%cth  >0)
       IN%sum%ctp_log = MERGE(IN%sum%ctp_log +LOG(tmpctp_log),IN%sum%ctp_log, IN%av%ctp_log>0)
       IN%sum%ctp     = MERGE(IN%sum%ctp     +IN%av%ctp,      IN%sum%ctp,     IN%av%ctp  >0)
       IN%sum%ctt     = MERGE(IN%sum%ctt     +IN%av%ctt,      IN%sum%ctt,     IN%av%ctt  >0)
       IN%sum%ref_ice = MERGE(IN%sum%ref_ice +IN%av%ref_ice,  IN%sum%ref_ice, IN%av%ref_ice>0)
       IN%sum%cot_ice = MERGE(IN%sum%cot_ice +IN%av%cot_ice,  IN%sum%cot_ice, IN%av%cot_ice >0)
       IN%sum%iwp     = MERGE(IN%sum%iwp     +IN%av%iwp,      IN%sum%iwp,     IN%av%iwp  >0)
       IN%sum%ref_liq = MERGE(IN%sum%ref_liq +IN%av%ref_liq,  IN%sum%ref_liq, IN%av%ref_liq>0)
       IN%sum%cot_liq = MERGE(IN%sum%cot_liq +IN%av%cot_liq,  IN%sum%cot_liq, IN%av%cot_liq >0)
       IN%sum%lwp     = MERGE(IN%sum%lwp     +IN%av%lwp,      IN%sum%lwp,     IN%av%lwp  >0)
       IN%sum%tau     = MERGE(IN%sum%tau     +IN%av%tau,      IN%sum%tau,     IN%av%tau  >0)

       ! Get the numel
       IN%numel%cld   = MERGE(IN%numel%cld+1, IN%numel%cld, IN%av%cfc>0 )
       IN%numel%day   = MERGE(IN%numel%day+1, IN%numel%day, S%sunlit.EQ.1)
       IN%numel%grd   = IN%numel%grd+1
    END WHERE
    IN%sum%hist2d_cot_ctp = IN%sum%hist2d_cot_ctp+IN%av%hist2d_cot_ctp

    ! reset these values
    CALL INITIALISE_CLAAS_SIM(IN%av,ngrids,-999._wp,n_pbins,n_tbins)

  END SUBROUTINE DAY_ADD

  SUBROUTINE DAY_AVERAGE(IN,ngrids)

    IMPLICIT NONE

    TYPE(CLAAS_type), INTENT(inout) :: IN
    INTEGER, INTENT(in) :: ngrids

    REAL(wp), DIMENSION(ngrids) :: &
         numel_grd, & ! number of grids
         numel_cld, & ! number of cloudy grids
         numel_icld,& ! number of ice cloudy grids
         numel_lcld,& ! number of liquid cloudy grids
         numel_day    ! number of daylit

    ! need temporary variables to avoid divide by 0
    numel_grd  = IN%numel%grd
    numel_cld  = MERGE(IN%numel%cld ,0.1_wp,IN%numel%cld>0)
    numel_icld = MERGE(IN%numel%icld,0.1_wp,IN%numel%icld>0)
    numel_lcld = MERGE(IN%numel%lcld,0.1_wp,IN%numel%lcld>0)
    numel_day  = MERGE(IN%numel%day, 0.1_wp,IN%numel%day>0)

    ! always defined
    IN%av%cfc     = MERGE(IN%sum%cfc     /numel_grd,-999._wp,numel_grd>0.1_wp)
    IN%av%cfc_low = MERGE(IN%sum%cfc_low /numel_grd,-999._wp,numel_grd>0.1_wp)
    IN%av%cfc_mid = MERGE(IN%sum%cfc_mid /numel_grd,-999._wp,numel_grd>0.1_wp)
    IN%av%cfc_high= MERGE(IN%sum%cfc_high/numel_grd,-999._wp,numel_grd>0.1_wp)

    ! always defined if cloud is present
    IN%av%cth    = MERGE(IN%sum%cth    /numel_cld,-999._wp,numel_cld>0.1_wp)
    IN%av%ctp    = MERGE(IN%sum%ctp    /numel_cld,-999._wp,numel_cld>0.1_wp)
    IN%av%ctp_log= MERGE(IN%sum%ctp_log/numel_cld,-999._wp,numel_cld>0.1_wp)
    IN%av%ctt    = MERGE(IN%sum%ctt    /numel_cld,-999._wp,numel_cld>0.1_wp)

    ! only defined during daytime
    IN%av%cfc_day= MERGE(IN%sum%cfc_day/numel_day,-999._wp,numel_day>0.1_wp)
    IN%av%cot_ice= MERGE(IN%sum%cot_ice/numel_day,-999._wp,numel_day>0.1_wp)
    IN%av%iwp    = MERGE(IN%sum%iwp    /numel_day,-999._wp,numel_day>0.1_wp)
    IN%av%cot_liq= MERGE(IN%sum%cot_liq/numel_day,-999._wp,numel_day>0.1_wp)
    IN%av%lwp    = MERGE(IN%sum%lwp    /numel_day,-999._wp,numel_day>0.1_wp)
    IN%av%tau    = MERGE(IN%sum%tau    /numel_day,-999._wp,numel_day>0.1_wp)

    ! only if correct phase/retrievable (also daytime)
    IN%av%ref_ice= MERGE(IN%sum%ref_ice/numel_icld,-999._wp,numel_icld>0.1_wp)
    IN%av%ref_liq= MERGE(IN%sum%ref_liq/numel_lcld,-999._wp,numel_lcld>0.1_wp)

    ! postprocessing
    IN%av%ctp_log  = MERGE(EXP(IN%av%ctp_log),-999._wp,IN%av%ctp_log>0  )
    IN%av%hist2d_cot_ctp = IN%sum%hist2d_cot_ctp

  END SUBROUTINE DAY_AVERAGE

  SUBROUTINE CHECK_VARIABLES_SUBGRID(d1,ncol,M,S,I,options,frac_out,print_profile)
    ! this is for valid values if there is a cloud

    IMPLICIT NONE

    INTEGER,          INTENT(in) :: d1,ncol
    TYPE(model_type), INTENT(in) :: M
    TYPE(subset),     INTENT(in) :: S
    TYPE(internal),   INTENT(in) :: I
    REAL(wp),         INTENT(in) :: frac_out(ncol,M%aux%nlev)
    TYPE(name_list),  INTENT(in) :: options
    LOGICAL, OPTIONAL,INTENT(in) :: print_profile
   INTEGER                      :: ins,ngrids,nlev
    LOGICAL                      :: whammy(ncol,M%aux%nlev)
    CHARACTER(len=1000)          :: cloudstr
    LOGICAL                      :: data_mask(M%aux%ngrids,M%aux%nlev)
    CHARACTER(5)                 :: form
    CHARACTER(30)                :: var
    LOGICAL                      :: mask(ncol)

    ngrids= M%aux%ngrids
    nlev  = M%aux%nlev

    ! keep track of which columns fail any variable check
    IF (PRESENT(print_profile)) THEN
       whammy = .TRUE.
    ELSE
       whammy = .FALSE.
    END IF

    ! start checking ...
    whammy = whammy .OR. COMMON_CHECK_SUBGRID(I,ncol,nlev)

    mask=I%cflag.LE.1

    IF (ANY(.NOT.mask)) THEN

       var='I%cth'
       form='f9.2'
       whammy = whammy.OR.CHECK_VAR(var,form,ncol,nlev,&
            cth_min,cth_max,I%cth,mask)

       var='I%ctp'
       form='f9.2'
       whammy = whammy.OR.CHECK_VAR(var,form,ncol,nlev,&
            ctp_min,ctp_max,I%ctp,mask,scale=0.01_wp)

       var='I%ctt'
       form='f9.2'
       whammy = whammy.OR.CHECK_VAR(var,form,ncol,nlev,&
            ctt_min,ctt_max,I%ctt,mask)

       IF (S%sunlit(d1).EQ.1) THEN
          var='I%reff'
          form='f7.2'
          whammy = whammy.OR.CHECK_VAR(var,form,ncol,nlev,&
               ref_liq_min,ref_ice_max,I%reff,mask)

          var='I%cwp'
          form='f9.2'
          whammy = whammy.OR.CHECK_VAR(var,form,ncol,nlev,&
               lwp_min,iwp_max,I%cwp,scale=1000._wp)

       END IF
    END IF

    IF (ANY(whammy).OR.PRESENT(print_profile)) THEN
       data_mask = .TRUE. !initialise to ignore every data point
       data_mask(d1,1:nlev) = .FALSE. ! except this point
       DO ins = 1,ncol

          IF (ANY(whammy(ins,:)).OR.PRESENT(print_profile)) THEN
             CALL CHECK_VARIABLES(ngrids,nlev,options,data_mask,M,S,&
                  scops_in=frac_out(ins,1:nlev))
             PRINT *," "
             PRINT '(2(a,1x,I4,1x))',&
                  "These are all the internal variables at d1=",d1,"ins=",ins
             SELECT CASE(I%cflag(ins))
             CASE (0)
                cloudstr = 'cloudfree'
             CASE (1)
                cloudstr = 'sub-visible'
             CASE (2)
                cloudstr = 'semi-transparent'
             CASE (3)
                cloudstr = 'opaque'
             END SELECT

             PRINT '(A19,I9,1x,"(",a,")")',"I%cflag = ",I%cflag (ins),TRIM(cloudstr)
             PRINT '(A19,I9)',"I%cph = ",I%cph(ins)
             IF (I%cflag(ins) .GT. 1) THEN
                IF (S%sunlit(d1).EQ.1) THEN
                   PRINT '(A19,F9.2," [micron]")',"I%reff = ",&
                        I%reff(ins)
                   PRINT '(A19,F9.2," [g/m2]")',"I%cwp = ",&
                        1000*I%cwp (ins)
                END IF
                PRINT '(A19,F9.2)',"I%cth     = ",I%cth    (ins)
                PRINT '(A19,F9.2)',"I%ctp     = ",I%ctp    (ins)/100
                PRINT '(A19,F9.2)',"I%ctt     = ",I%ctt    (ins)
             END IF
             IF (S%sunlit(d1).EQ.1) PRINT '(A19,F9.3)',"tau =",I%tau(ins)

             PRINT '(A19,1x,9x,F8.3,1x,9x," [K]")',"I%Tb=",I%Tb(ins)
             CALL CHECK_PROFILES(d1,ins,nlev,S,I,frac_out(ins,:))
             PRINT '(A13,I9)',"ins = ",ins
             IF (.NOT. PRESENT(print_profile)) &
                  STOP "The internal simulator values are outside the valid range"
          END IF
       END DO
    END IF
  END SUBROUTINE CHECK_VARIABLES_SUBGRID

  SUBROUTINE CHECK_GRID_AVERAGES(M,S,options,IN,print_profile)

    ! function for checking the consistency of simulated variables
    ! after grid averaging

    IMPLICIT NONE

    TYPE(name_list), INTENT(in)    :: options
    TYPE(model_type), INTENT(in)   :: M
    TYPE(subset), INTENT(in)       :: S
    TYPE(CLAAS_fields), INTENT(in) :: IN
    LOGICAL, OPTIONAL, INTENT(in)  :: print_profile

    INTEGER  :: n_tbins,n_pbins,ngrids,nlev,pb,tb,d1
    CHARACTER(len=100)   :: str,str2
    LOGICAL              :: whammy(M%aux%ngrids,M%aux%nlev)
    LOGICAL              :: dm(M%aux%ngrids,M%aux%nlev),dm1(M%aux%ngrids)
    CHARACTER(5)         :: form
    CHARACTER(30)        :: var

    ngrids  = M%aux%ngrids
    nlev    = M%aux%nlev
    n_tbins = options%ctp_tau%n_tbins
    n_pbins = options%ctp_tau%n_pbins

    dm  = S%data_mask
    dm1 = dm(1:ngrids,1)

    IF (PRESENT(print_profile)) THEN
       whammy = .TRUE.
    ELSE
       whammy = .FALSE.
    END IF

    var='AV%cfc'
    form='f8.3'
    whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,&
         frac_min,frac_max,IN%cfc,dm)

    var='AV%cfc_day'
    form='f8.3'
    whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,&
         frac_min,frac_max,IN%cfc_day,dm)

    var='AV%cfc_low'
    form='f8.3'
    whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,&
         frac_min,frac_max,IN%cfc_low,dm)

    var='AV%cfc_mid'
    form='f8.3'
    whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,&
         frac_min,frac_max,IN%cfc_mid,dm)

    var='AV%cfc_high'
    form='f8.3'
    whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,&
         frac_min,frac_max,IN%cfc_high,dm)

    var='AV%cth'
    form='f9.2'
    whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,&
         cth_min,cth_max,IN%cth,dm,except=missing)

    var='AV%ctp_log'
    form='f9.2'
    whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,&
         ctp_min,ctp_max,IN%ctp_log,dm,scale=0.01_wp,except=missing)

    var='AV%ctp'
    form='f9.2'
    whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,&
         ctp_min,ctp_max,IN%ctp,dm,scale=0.01_wp,except=missing)

    var='AV%ctt'
    form='f9.2'
    whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,&
         ctt_min,ctt_max,IN%ctt,dm,except=missing)

    var='AV%cot_ice'
    form='f9.2'
    whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,&
         tau_min,tau_max,IN%cot_ice,dm,except=missing)

    var='AV%cot_liq'
    form='f9.2'
    whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,&
         tau_min,tau_max,IN%cot_liq,dm,except=missing)

    var='AV%iwp'
    form='f9.2'
    whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,&
         iwp_min,iwp_max,IN%iwp,dm,scale=1000._wp,except=missing)

    var='AV%lwp'
    form='f9.2'
    whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,&
         lwp_min,lwp_max,IN%lwp,dm,scale=1000._wp,except=missing)

    var='AV%ref_ice'
    form='f9.2'
    whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,&
         ref_ice_min,ref_ice_max,IN%ref_ice,dm,except=missing)

    var='AV%ref_liq'
    form='f9.2'
    whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,&
         ref_liq_min,ref_liq_max,IN%ref_liq,dm,except=missing)

    var='AV%tau'
    form='f9.2'
    whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,&
         tau_min,tau_max,IN%tau,dm,except=missing)

    WRITE(str,'(a,I2,a)')        "(5x,",n_tbins,"(F7.2,1x))"
    WRITE(str2,'(a,":",1x,I2,a)') "(I4",n_tbins,"(I7,1x))"

    IF (ANY(whammy)) THEN
       dm = .TRUE. !initialise to ignore every data point
       FAIL:DO d1 = 1,ngrids
          IF (ANY(whammy(d1,:))) THEN
             dm(d1,1:nlev) = .FALSE. ! except this point
             PRINT *," "
             CALL CHECK_VARIABLES(ngrids,nlev,options,dm,M,S,.TRUE.)
             PRINT *," "
             PRINT '(a,1x,I4)',&
                  "These are all the gridded averages at d1=",d1

             PRINT '(a)',"CLAAS"
             PRINT '(A14,F9.2," [%]")',     "IN%cfc      = ",IN%cfc(d1)*100
             PRINT '(A14,F9.2," [m]")',     "IN%cth      = ",IN%cth(d1)
             PRINT '(A14,F9.2," [hPa]")',   "IN%ctp_log  = ",IN%ctp_log(d1)/100
             PRINT '(A14,F9.2," [hPa]")',   "IN%ctp      = ",IN%ctp(d1)/100
             PRINT '(A14,F9.2," [K]")',     "IN%ctt      = ",IN%ctt(d1)
             PRINT '(A14,F9.2," [micron]")',"IN%ref_ice  = ",IN%ref_ice(d1)
             PRINT '(A14,F9.2," [-]")',     "IN%cot_ice  = ",IN%cot_ice(d1)
             PRINT '(A14,F9.2," [g/m2]")',  "IN%iwp      = ",&
                  MERGE(IN%iwp(d1)*1000,IN%iwp(d1),IN%iwp(d1).GE.0)
             PRINT '(A14,F9.2," [g/m2]")',  "IN%iwp      = ",IN%iwp(d1)*1000
             PRINT '(A14,F9.2," [micron]")',"IN%ref_liq  = ",IN%ref_liq(d1)
             PRINT '(A14,F9.2," [-]")',     "IN%cot_liq  = ",IN%cot_liq(d1)
             PRINT '(A14,F9.2," [g/m2]")',  "IN%lwp      = ",&
                  MERGE(IN%lwp(d1)*1000,IN%lwp(d1),IN%lwp(d1).GE.0)
             PRINT '(A14,F9.2," [deg]")',   "S%solzen       = ",S%solzen(d1)
             PRINT '(A14,F9.2," [-]")',     "IN%tau      = ",IN%tau(d1)

             PRINT str,(options%ctp_tau%tbin_edges(tb),tb=1,n_tbins)
             DO pb = 1,n_pbins
                PRINT str2,&
                     INT(options%ctp_tau%pbin_edges(pb)/100),&
                     (INT(IN%hist2d_cot_ctp(d1,tb,pb,1)),tb=1,n_tbins)
             END DO
             PRINT str,(options%ctp_tau%tbin_edges(tb),tb=1,n_tbins)

             PRINT str,(options%ctp_tau%tbin_edges(tb),tb=1,n_tbins)
             DO pb = 1,n_pbins
                PRINT str2,&
                     INT(options%ctp_tau%pbin_edges(pb)/100),&
                     (INT(IN%hist2d_cot_ctp(d1,tb,pb,2)),tb=1,n_tbins)
             END DO
             PRINT str,(options%ctp_tau%tbin_edges(tb),tb=1,n_tbins)
             EXIT FAIL
          END IF
       END DO FAIL
    END IF
  END SUBROUTINE CHECK_GRID_AVERAGES

  SUBROUTINE CHECK_PROFILES(d1,ins,nlev,S,I,frac_out)

    ! function that produces useful output to the screen to visualise if
    ! the cloud tests are doing what we expect

    IMPLICIT NONE

    INTEGER,        INTENT(in) :: d1,ins,nlev
    TYPE(internal), INTENT(in) :: I
    TYPE(subset),   INTENT(in) :: S
    REAL(wp),       INTENT(in) :: frac_out(nlev)
    INTEGER                    :: inl,inv_ind
    INTEGER                    :: inversion(nlev)
    CHARACTER(len=1000)        :: arrow
    CHARACTER(len=1000)        :: description
    LOGICAL                    :: event
    CHARACTER(len=1000)        :: str
    INTEGER                    :: ii
    LOGICAL                    :: CLAAS_found

1   FORMAT ((A5,1x),2(A11,1x),(A9,1x),(A10,1x),(A8,1x),(A5,1x),(A4,1x),3(A7,1x),(A6,1x),A26)
    ! levels
2   FORMAT ((I5,1x),(F6.1,6x),(12x),(F6.1,4x),(F8.2,3x),(F8.2,1x),(2x,I1,3x),(5x),3(8x),(6x),(A6,1x),A26)
    ! layers
3   FORMAT ((F4.1,2x),(12x),(F6.1,6x),(10x),(11x),(9x),(6x),(F4.3,1x),3(F7.2,1x),(7x),26x)

    PRINT *," --------------- CTTH --------------"
    PRINT 1,"level","P_mid [hPa]","P_int [hPa]","Tcorr [K]","height [m]",&
         "tau","scops","cfra","St:land","St:sea","Cirrus","WHERE","Explanation"

    CLAAS_found    = .FALSE.
    inv_ind = 1
    inversion(1:nlev) = -999

    DO inl = nlev,1,-1
       IF (S%inv_layers(d1,inl) .GT. 0) THEN
          inversion(inv_ind) = S%inv_layers(d1,inl)
          inv_ind = inv_ind+1
       END IF
    END DO

    inv_ind = 1
    DO inl = 1,nlev

       ! from the top of the atmosphere down
       event       = .FALSE.
       arrow       = ""
       description = ""

       ! inversion (listed from the ground up)
       IF (inversion(inv_ind) .EQ. inl) THEN
          IF (inversion(inv_ind) .EQ. &
               MINVAL(inversion(1:nlev),MASK=inversion(1:nlev).GT. 0)) THEN
             WRITE(description,'(a,a)') TRIM(description),",Tropopause"
          ELSE
             WRITE(description,'(a,a)') TRIM(description),",Inversion"
          END IF
          event = .TRUE.
          inv_ind = inv_ind+1
       END IF
       IF ( (I%ctp(ins).LE.S%p_mid(d1,inl)+epsR).AND. &
            .NOT.CLAAS_found.AND. &
            .NOT.(I%ctp(ins).EQ.missing) ) THEN
          WRITE(description,'(a,a)') TRIM(description),",CLAAS"
          CLAAS_found = .TRUE.
          event = .TRUE.
       END IF

       IF (event) arrow = '<-----'

       PRINT 2,&
            inl,&
            S%p_mid       (d1 ,inl)/100,&
            S%Tcorr       (d1 ,inl),&
            (S%height(d1,inl)+S%height(d1,inl+1))/2,&
            I%tau_profile (ins,inl),&
            INT(frac_out  (    inl)),&
            arrow,&
            TRIM(description(2:200))

       IF (inl .LT. nlev) THEN
          IF (I%cfrac     (ins,inl).GT. 0) THEN
             PRINT 3,&
                  REAL(inl+.5),&
                  S%p_int(d1    ,inl+1)/100,&
                  I%cfrac(ins   ,inl)      ,&
                  I%cloud(ins,1,inl)       ,&
                  I%cloud(ins,2,inl)       ,&
                  I%cloud(ins,6,inl)
          END IF
       END IF
    END DO
    PRINT 1,"level","P_mid [hPa]","P_int [hPa]","Tcorr [K]","height [m]",&
         "tau","scops","cfrac","St:land","St:sea","Cirrus","WHERE","Explanation"
    PRINT *,""
    WRITE(str,'(a,I2,a)') "(a,",inv_ind-1,"(I2,1x))"
    PRINT TRIM(str),"inv_layers = ",(S%inv_layers(d1,ii),ii=1,inv_ind-1)
  END SUBROUTINE CHECK_PROFILES

END MODULE CLAAS_FUNCTIONS
