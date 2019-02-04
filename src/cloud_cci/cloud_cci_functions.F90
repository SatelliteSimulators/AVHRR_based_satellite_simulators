MODULE CLOUD_CCI_FUNCTIONS
  ! Cloud_cci
  !
  ! Salomon.Eliasson@smhi.se

  USE CLOUD_CCI_M,               ONLY: &
       CLOUD_CCI_TYPE,                 &
       CLOUD_CCI_FIELDS,               &
       INITIALISE_CLOUD_CCI_SIM
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
  USE FROM_CC4CL,                ONLY: &
       GRID_DIMENSION_READ,            &
       READ_VALUES_3D
  USE HANDY,                     ONLY: &
       BUILD_FILENAME,                 &
       GET_LUN,                        &
       TO_UPPER
  USE INTERNAL_SIMULATOR,        ONLY: INTERNAL
  USE MODEL_INPUT,               ONLY: MODEL_TYPE
  USE NAMELIST_INPUT,            ONLY: NAME_LIST
  USE OPTICS_M,                  ONLY: &
       CLOUD_ALBEDO,                   &
       GZERO_T
  USE SIMULATOR_INPUT_VARIABLES, ONLY: SUBSET
  USE SIMPLE_CTTH_SIMULATOR,     ONLY: &
       CTH_SIMPLE,                     &
       CTP_SIMPLE,                     &
       CTT_SIMPLE

  IMPLICIT NONE

  PUBLIC :: CHECK_GRID_AVERAGES,&
       CHECK_PROFILES,          &
       CHECK_VARIABLES_SUBGRID, &
       ALBEDO,                  &
       CTTH,                    &
       GRID_AVERAGE,            &
       DAY_ADD,                 &
       DAY_AVERAGE,             &
       READ_ALBEDO_LUT

  INTEGER, PARAMETER :: IXM1 = -1
  INTEGER, PARAMETER :: iX0  =  0
  INTEGER, PARAMETER :: iX1  =  1
  INTEGER, PARAMETER :: iXp1 =  2

  REAL(wp),PARAMETER :: cer_liq_min = 1,cer_liq_max = 23,&
       lwp_min = cer_liq_min*tau_min,&
       lwp_max = cer_liq_max*tau_max
  REAL(wp),PARAMETER :: cer_ice_min = 4,cer_ice_max = 92,&
       iwp_min = cer_ice_min*tau_min,&
       iwp_max = cer_ice_max*tau_max

CONTAINS

  SUBROUTINE CTTH(d1,ins,nlev,S,inter)
    !
    ! Purpose: Get the CTX products
    !
    ! NOTE: According to validation reports from Caroline Poulsen the
    ! CTX products are sensitive to the cloud height where tau = 1
    ! (uncorrected) from the cloud top (like the modis simulator). For
    ! the corrected algorithm the limit is set somewhat ambiguously to tau = 0.3
    !
    !
    IMPLICIT NONE

    TYPE(subset), INTENT(in)     :: S
    INTEGER, INTENT(in)          :: d1,ins,nlev
    TYPE(internal),INTENT(inout) :: inter
    REAL(wp), PARAMETER          :: uncorr_ctx_lim = 1
    REAL(wp), PARAMETER          :: corr_ctx_lim   = 0.3
    REAL(wp)                     :: T(nlev)

    T = S%Tcorr(d1,1:nlev) 
    IF (inter%cflag(ins) .GT. 1) THEN

       ! i.e. not subvisible (only semi-transparent or opaque)

       inter%cth(ins)             = CTH_SIMPLE(d1,ins,nlev  ,S,inter,uncorr_ctx_lim)
       inter%ctp(ins)             = CTP_SIMPLE(d1,ins,nlev  ,S,inter,uncorr_ctx_lim)
       inter%ctt(ins)             = CTT_SIMPLE(d1,ins,nlev,T,S,inter,uncorr_ctx_lim)
       inter%cloud_cci%cth_c(ins) = CTH_SIMPLE(d1,ins,nlev  ,S,inter,corr_ctx_lim)
       inter%cloud_cci%ctp_c(ins) = CTP_SIMPLE(d1,ins,nlev  ,S,inter,corr_ctx_lim)
       inter%cloud_cci%ctt_c(ins) = CTT_SIMPLE(d1,ins,nlev,T,S,inter,corr_ctx_lim)
    END IF

  END SUBROUTINE CTTH

  SUBROUTINE GRID_AVERAGE(d1,S,inter,ncol,IN)

    !
    ! Calculate simulator averages on the model grid scale.

    IMPLICIT NONE

    INTEGER, INTENT(in)                   :: d1,ncol
    TYPE(subset), INTENT(in)              :: S
    TYPE(internal), INTENT(in)            :: inter
    TYPE(cloud_cci_fields), INTENT(inout) :: IN

    ! internal
    INTEGER :: ncld,nliq,nice,nalb
    LOGICAL,DIMENSION(ncol) :: iscloud,isliquid,isice,isalb,islow,ismid,ishigh
    REAL(wp):: tmpCTP_log(ncol)
    REAL(wp):: re_fill=-999

    iscloud  (1:ncol) = inter%cflag(1:ncol).GE.2 ! i.e., tau>tau_min
    isliquid (1:ncol) = (inter%cph (1:ncol).EQ.1).AND.iscloud(1:ncol)
    isice    (1:ncol) = (inter%cph (1:ncol).EQ.2).AND.iscloud(1:ncol)
    isalb    (1:ncol) = inter%cloud_cci%albedoIsDefined(1:ncol)

    islow  = inter%ctp.GE.cloud_low
    ismid  = inter%ctp.LT.cloud_low.AND.inter%ctp.GE.cloud_mid
    ishigh = inter%ctp.LT.cloud_mid.AND.iscloud
    
    ncld = COUNT(iscloud)
    nliq = COUNT(isliquid)
    nice = COUNT(isice)
    nalb = COUNT(isalb)

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
       IN%ctp    (d1) = SUM(inter%ctp,            MASK=iscloud) /ncld
       IN%cth    (d1) = SUM(inter%cth,            MASK=iscloud) /ncld
       IN%ctt    (d1) = SUM(inter%ctt,            MASK=iscloud) /ncld
       IN%cth_c  (d1) = SUM(inter%cloud_cci%cth_c,MASK=iscloud) /ncld
       IN%ctp_c  (d1) = SUM(inter%cloud_cci%ctp_c,MASK=iscloud) /ncld
       IN%ctt_c  (d1) = SUM(inter%cloud_cci%ctt_c,MASK=iscloud) /ncld

       IF (S%sunlit(d1).EQ.1) THEN
          ! variables that depend on the sunlight

          IF (nliq .GT. 0) &
               IN%cer_liq(d1)   = SUM(inter%reff, MASK=isliquid)/nliq
          IF (nice .GT. 0) &
               IN%cer_ice(d1)   = SUM(inter%reff, MASK=isice)/nice
          IF (nalb .GT. 0) &
               IN%cla_vis006(d1)= SUM(inter%cloud_cci%cla_vis006,MASK=isalb)/nalb

          !
          !--------

          ! ---------
          ! Variables that are always defined during sunlight
          !
          IN%cot_ice(d1)  = SUM(inter%tau,MASK=isice)   /ncol
          IN%iwp    (d1)  = SUM(inter%cwp,MASK=isice)   /ncol
          IN%lwp    (d1)  = SUM(inter%cwp,MASK=isliquid)/ncol
          IN%cot_liq(d1)  = SUM(inter%tau,MASK=isliquid)/ncol
          IN%tau    (d1)  = SUM(inter%tau,MASK=iscloud) /ncol
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

    TYPE(cloud_cci_type), INTENT(inout) :: IN
    TYPE(subset), INTENT(in)            :: S
    INTEGER, INTENT(in)                 :: ngrids,n_tbins,n_pbins

    REAL(wp) :: tmpctp_log(ngrids)

    WHERE(.NOT.S%data_mask(:,1))

       ! Get the sum
       tmpctp_log       = MERGE(IN%av%ctp_log, 0.1_wp,IN%av%ctp_log >0)
       IN%sum%cfc       = MERGE(IN%sum%cfc     +IN%av%cfc,     IN%sum%cfc,     IN%av%cfc     >0) 
       IN%sum%cfc_low   = MERGE(IN%sum%cfc_low +IN%av%cfc_low, IN%sum%cfc_low, IN%av%cfc_low >0) 
       IN%sum%cfc_mid   = MERGE(IN%sum%cfc_mid +IN%av%cfc_mid, IN%sum%cfc_mid, IN%av%cfc_mid >0) 
       IN%sum%cfc_high  = MERGE(IN%sum%cfc_high+IN%av%cfc_high,IN%sum%cfc_high,IN%av%cfc_high>0) 
       IN%sum%cla_vis006= MERGE(IN%sum%cla_vis006+IN%av%cla_vis006,IN%sum%cla_vis006,IN%av%cla_vis006>0)
       IN%sum%cth       = MERGE(IN%sum%cth    +IN%av%cth,      IN%sum%cth,     IN%av%cth    >0)
       IN%sum%ctp_log   = MERGE(IN%sum%ctp_log+LOG(tmpctp_log),IN%sum%ctp_log, IN%av%ctp_log>0) 
       IN%sum%ctp       = MERGE(IN%sum%ctp    +IN%av%ctp,      IN%sum%ctp,     IN%av%ctp    >0) 
       IN%sum%ctt       = MERGE(IN%sum%ctt    +IN%av%ctt,      IN%sum%ctt,     IN%av%ctt    >0) 
       IN%sum%cth_c     = MERGE(IN%sum%cth_c  +IN%av%cth_c,    IN%sum%cth_c,   IN%av%cth_c  >0)
       IN%sum%ctp_c     = MERGE(IN%sum%ctp_c  +IN%av%ctp_c,    IN%sum%ctp_c,   IN%av%ctp_c  >0)
       IN%sum%ctt_c     = MERGE(IN%sum%ctt_c  +IN%av%ctt_c,    IN%sum%ctt_c,   IN%av%ctt_c  >0)
       IN%sum%cer_ice   = MERGE(IN%sum%cer_ice+IN%av%cer_ice,  IN%sum%cer_ice, IN%av%cer_ice>0) 
       IN%sum%cot_ice   = MERGE(IN%sum%cot_ice+IN%av%cot_ice,  IN%sum%cot_ice, IN%av%cot_ice>0) 
       IN%sum%iwp       = MERGE(IN%sum%iwp    +IN%av%iwp,      IN%sum%iwp,     IN%av%iwp    >0) 
       IN%sum%cer_liq   = MERGE(IN%sum%cer_liq+IN%av%cer_liq,  IN%sum%cer_liq, IN%av%cer_liq>0) 
       IN%sum%cot_liq   = MERGE(IN%sum%cot_liq+IN%av%cot_liq,  IN%sum%cot_liq, IN%av%cot_liq>0)
       IN%sum%lwp       = MERGE(IN%sum%lwp    +IN%av%lwp,      IN%sum%lwp,    IN%av%lwp     >0) 
       IN%sum%tau   = MERGE(IN%sum%tau  +IN%av%tau,    IN%sum%tau,  IN%av%tau   >0) 

       ! Get the numel
       IN%numel%alb   = MERGE(IN%numel%alb+1, IN%numel%alb, IN%av%cla_vis006>0.AND.S%sunlit.EQ.1)
       IN%numel%cld   = MERGE(IN%numel%cld+1, IN%numel%cld, IN%av%cfc>0 )
       IN%numel%day   = MERGE(IN%numel%day+1, IN%numel%day, S%sunlit.EQ.1)
       IN%numel%grd   = IN%numel%grd+1
    END WHERE
    IN%sum%hist2d_cot_ctp = IN%sum%hist2d_cot_ctp+IN%av%hist2d_cot_ctp

    ! reset these values
    CALL INITIALISE_CLOUD_CCI_SIM(IN%av,ngrids,-999._wp,n_pbins,n_tbins)

  END SUBROUTINE DAY_ADD

  SUBROUTINE DAY_AVERAGE(IN,ngrids)

    IMPLICIT NONE

    TYPE(cloud_cci_type), INTENT(inout) :: IN
    INTEGER, INTENT(in) :: ngrids

    REAL(wp), DIMENSION(ngrids) :: &
         numel_alb, & ! number of valid albedos
         numel_grd, & ! number of grids
         numel_cld, & ! number of cloudy grids
         numel_icld,& ! number of ice cloudy grids
         numel_lcld,& ! number of liquid cloudy grids
         numel_day    ! number of daylit

    ! need temporary variables to avoid divide by 0
    numel_alb  = MERGE(IN%numel%alb ,0.1_wp,IN%numel%alb>0)
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
    IN%av%cth_c  = MERGE(IN%sum%cth_c  /numel_cld,-999._wp,numel_cld>0.1_wp)
    IN%av%ctp_c  = MERGE(IN%sum%ctp_c  /numel_cld,-999._wp,numel_cld>0.1_wp)
    IN%av%ctt_c  = MERGE(IN%sum%ctt_c  /numel_cld,-999._wp,numel_cld>0.1_wp) 

    ! only defined during daytime  
    IN%av%cot_ice= MERGE(IN%sum%cot_ice/numel_day,-999._wp,numel_day>0.1_wp) 
    IN%av%iwp    = MERGE(IN%sum%iwp    /numel_day,-999._wp,numel_day>0.1_wp)
    IN%av%cot_liq= MERGE(IN%sum%cot_liq/numel_day,-999._wp,numel_day>0.1_wp) 
    IN%av%lwp    = MERGE(IN%sum%lwp    /numel_day,-999._wp,numel_day>0.1_wp) 
    IN%av%tau    = MERGE(IN%sum%tau    /numel_day,-999._wp,numel_day>0.1_wp)

    ! only if correct phase/retrievable (also daytime)
    IN%av%cla_vis006= MERGE(IN%sum%cla_vis006/numel_alb ,-999._wp,numel_alb >0.1_wp)
    IN%av%cer_ice   = MERGE(IN%sum%cer_ice   /numel_icld,-999._wp,numel_icld>0.1_wp) 
    IN%av%cer_liq   = MERGE(IN%sum%cer_liq   /numel_lcld,-999._wp,numel_lcld>0.1_wp) 

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

    INTEGER :: ins,ngrids,nlev
    LOGICAL               :: whammy(ncol,M%aux%nlev)
    CHARACTER(len=1000)   :: cloudstr
    LOGICAL               :: data_mask(M%aux%ngrids,M%aux%nlev)
    CHARACTER(5)          :: form
    CHARACTER(30)         :: var
    LOGICAL               :: mask(ncol)

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

       var='I%cth_c'
       form='f9.2'
       whammy = whammy.OR.CHECK_VAR(var,form,ncol,nlev,&
            cth_min,cth_max,I%cloud_cci%cth_c,mask)

       var='I%ctp_c'
       form='f9.2'
       whammy = whammy.OR.CHECK_VAR(var,form,ncol,nlev,&
            ctp_min,ctp_max,I%cloud_cci%ctp_c,mask,scale=0.01_wp)

       var='I%ctt_c'
       form='f9.2'
       whammy = whammy.OR.CHECK_VAR(var,form,ncol,nlev,&
            ctt_min,ctt_max,I%cloud_cci%ctt_c,mask)

       IF (S%sunlit(d1).EQ.1) THEN
          var='I%reff'
          form='f7.2'
          whammy = whammy.OR.CHECK_VAR(var,form,ncol,nlev,&
               cer_liq_min,cer_ice_max,I%reff,mask)

          var='I%cwp'
          form='f9.2'
          whammy = whammy.OR.CHECK_VAR(var,form,ncol,nlev,&
               lwp_min,iwp_max,I%cwp,scale=1000._wp)

       END IF
    END IF

    IF (S%sunlit(d1).EQ.1) THEN
       
       
       mask=.NOT.I%cloud_cci%albedoIsDefined
       
       IF (ANY(.NOT.mask)) THEN
          
          ! cloud albedo is undefined if outside the range of the
          ! lookup tables (no extrapolation is done)
          var='I%cloud_cci%cla_vis006'
          form='f7.2'
          whammy = whammy.OR.CHECK_VAR(var,form,ncol,nlev,&
               0._wp,100._wp,I%cloud_cci%cla_vis006,mask)
          
       END IF
    END IF
    
    IF (ANY(whammy).OR.PRESENT(print_profile)) THEN
       data_mask = .TRUE. !initialise to ignore every data point
       data_mask(d1,1:nlev) = .FALSE. ! except this point
       DO ins = 1, ncol

          IF (ANY(whammy(ins,:)).OR.PRESENT(print_profile)) THEN
             CALL CHECK_VARIABLES(ngrids,nlev,options,data_mask,M,S,&
                  scops_in=frac_out(ins,1:nlev))
             PRINT *, " "
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
             PRINT '(A19,I9)',"I%cph = ", I%cph(ins)
             IF (I%cflag(ins) .GT. 1) THEN
                IF (S%sunlit(d1).EQ.1) THEN
                   PRINT '(A19,F9.2," [micron]")',"I%reff = ",&
                        I%reff(ins)
                   PRINT '(A19,F9.2," [g/m2]")',"I%cwp = ", &
                        1000*I%cwp (ins)
                END IF
                PRINT '(A19,F9.2)',"I%cth   = ",I%cth(ins)
                PRINT '(A19,F9.2)',"I%ctp   = ",I%ctp(ins)/100
                PRINT '(A19,F9.2)',"I%ctt   = ",I%ctt(ins)
                PRINT '(A19,F9.2)',"I%cth_c = ",I%cloud_cci%cth_c(ins)
                PRINT '(A19,F9.2)',"I%ctp_c = ",I%cloud_cci%ctp_c(ins)/100
                PRINT '(A19,F9.2)',"I%ctt_c = ",I%cloud_cci%ctt_c(ins)
                PRINT '(A19,F9.2," [%]")',"I%albedo = ",I%cloud_cci%cla_vis006(ins)
             END IF
             IF (S%sunlit(d1).EQ.1) PRINT '(A19,F9.3)',"tau =", I%tau(ins)
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

    TYPE(name_list), INTENT(in)           :: options
    TYPE(model_type), INTENT(in)          :: M
    TYPE(subset), INTENT(in)              :: S
    TYPE(cloud_cci_fields), INTENT(in)    :: IN
    LOGICAL, OPTIONAL, INTENT(in)         :: print_profile

    INTEGER  :: n_tbins,n_pbins,ngrids,nlev,pb,tb,d1
    CHARACTER(len=100)   :: str,str2
    LOGICAL              :: whammy(M%aux%ngrids,M%aux%nlev)
    LOGICAL              :: dm    (M%aux%ngrids,M%aux%nlev),dm1(M%aux%ngrids)
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
    
    var='AV%cla_vis006'
    form='f9.2' 
    whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,&
         0._wp,100._wp,IN%cla_vis006,dm,except=missing)

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

    var='AV%cer_ice'
    form='f9.2' 
    whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,&
         cer_ice_min,cer_ice_max,IN%cer_ice,dm,except=missing)

    var='AV%cer_liq'
    form='f9.2' 
    whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,&
         cer_liq_min,cer_liq_max,IN%cer_liq,dm,except=missing)

    var='AV%iwp'
    form='f9.2' 
    whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,&
         iwp_min,iwp_max,IN%iwp,dm,scale=1000._wp,except=missing)

    var='AV%lwp'
    form='f9.2' 
    whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,&
         lwp_min,lwp_max,IN%lwp,dm,scale=1000._wp,except=missing)

    var='AV%tau'
    form='f9.2' 
    whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,&
         tau_min,tau_max,IN%tau,dm,except=missing)

    WRITE(str, '(a,I2,a)')        "(5x,",n_tbins,"(F7.2,1x))"
    WRITE(str2,'(a,":",1x,I2,a)') "(I4",n_tbins,"(I7,1x))"

    IF (ANY(whammy)) THEN
       dm = .TRUE. !initialise to ignore every data point
       FAIL:DO d1 = 1,ngrids
          IF (ANY(whammy(d1,:))) THEN
             dm(d1,1:nlev) = .FALSE. ! except this point
             PRINT *, " "
             CALL CHECK_VARIABLES(ngrids,nlev,options,dm,M,S,.TRUE.)
             PRINT *, " "
             PRINT '(a,1x,I4)',&
                  "These are all the gridded averages at d1=",d1
             
             PRINT '(a)', "Cloud_cci"
             PRINT '(A14,F9.2," [%]")',     "IN%cfc      = ",IN%cfc(d1)*100
             PRINT '(A19,F9.2," [%]")',     "IN%cla_vis006 = ",IN%cla_vis006(d1)
             PRINT '(A14,F9.2," [m]")',     "IN%cth      = ",IN%cth(d1)
             PRINT '(A14,F9.2," [hPa]")',   "IN%ctp_log  = ",IN%ctp_log(d1)/100
             PRINT '(A14,F9.2," [hPa]")',   "IN%ctp      = ",IN%ctp(d1)/100
             PRINT '(A14,F9.2," [K]")',     "IN%ctt      = ",IN%ctt(d1)
             PRINT '(A14,F9.2," [micron]")',"IN%cer_ice  = ",IN%cer_ice(d1)
             PRINT '(A14,F9.2," [-]")',     "IN%cot_ice  = ",IN%cot_ice(d1)
             PRINT '(A14,F9.2," [g/m2]")',  "IN%iwp      = ",&
                  MERGE(IN%iwp(d1)*1000,IN%iwp(d1),IN%iwp(d1).GE.0)
             PRINT '(A14,F9.2," [g/m2]")',  "IN%iwp      = ",IN%iwp(d1)*1000
             PRINT '(A14,F9.2," [micron]")',"IN%cer_liq  = ",IN%cer_liq(d1)
             PRINT '(A14,F9.2," [-]")',     "IN%cot_liq  = ",IN%cot_liq(d1)
             PRINT '(A14,F9.2," [g/m2]")',  "IN%lwp      = ",&
                  MERGE(IN%lwp(d1)*1000,IN%lwp(d1),IN%lwp(d1).GE.0)
             PRINT '(A14,F9.2," [deg]")',   "S%solzen    = ",S%solzen(d1)
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
    LOGICAL                    :: cloud_cci_c_found,cloud_cci_found

1   FORMAT ((A5,1x),2(A11,1x),(A9,1x),(A10,1x),(A8,1x),(A5,1x),(A6,1x),A26)
    ! levels
2   FORMAT ((I5,1x),(F6.1,6x),(12x),(F6.1,4x),(F8.2,3x),(F8.2,1x),(2x,I1,3x),(A6,1x),A26)
    ! layers

    PRINT *," --------------- CTTH --------------"
    PRINT 1,"level","P_mid [hPa]","P_int [hPa]","Tcorr [K]","height [m]",&
         "tau","scops","WHERE","Explanation" 

    cloud_cci_found = .FALSE.
    cloud_cci_c_found = .FALSE.
    inv_ind = 1
    inversion(1:nlev) = -999

    DO inl = nlev,1,-1
       IF (S%inversion_layers(d1,inl) .GT. 0) THEN
          inversion(inv_ind) = S%inversion_layers(d1,inl)
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
       IF (I%ctp(ins) .LE. S%p_mid(d1,inl)+epsR .AND. &
            .NOT. cloud_cci_found .AND. &
            .NOT. I%ctp(ins) .EQ. missing) THEN
          WRITE(description,'(a,a)') TRIM(description),",cloud_cci"
          cloud_cci_found = .TRUE.
          event = .TRUE.
       END IF
       IF (I%cloud_cci%ctp_c(ins) .LE. S%p_mid(d1,inl)+epsR .AND. &
            .NOT. cloud_cci_c_found .AND. &
            .NOT. I%cloud_cci%ctp_c(ins) .EQ. missing) THEN
          WRITE(description,'(a,a)') TRIM(description),",cloud_cci-corrected"
          cloud_cci_c_found = .TRUE.
          event = .TRUE.
       END IF

       IF (event) arrow = '<-----'

       PRINT 2,&
            inl,&
            S%p_mid       (d1 ,inl)/100,&
            S%Tcorr       (d1 ,inl),&
            (S%height(d1,inl)+S%height(d1,inl))/2,&
            I%tau_profile (ins,inl),&
            frac_out      (    inl),&
            arrow,&
            TRIM(description(2:200)) 

    END DO
    PRINT 1,"level","P_mid [hPa]","P_int [hPa]","Tcorr [K]","height [m]",&
         "tau","scops","WHERE","Explanation" 
    PRINT *, ""
    WRITE(str,'(a,I2,a)') "(a,",inv_ind-1,"(I2,1x))"
    PRINT TRIM(str),"Inversion_layers = ",(S%inversion_layers(d1,ii),ii=1,inv_ind-1)
  END SUBROUTINE CHECK_PROFILES

  FUNCTION ALBEDO(Grid,GZero) RESULT(cloud_cci_albedo)

    ! Purpose: Retrieve the spectral cloud albedo using cloud_cci
    ! LUTs. modified from  Int_LUT_Common() in the CC4CL codes

    IMPLICIT NONE

    TYPE(cloud_albedo), INTENT(in) :: Grid
    TYPE(GZero_t), INTENT(in)      :: GZero

    ! OUT
    REAL(wp) :: cloud_cci_albedo

    !internal
    REAL(wp), DIMENSION(Grid%nTau,Grid%nSolzen,Grid%nRe)  :: F ! LUT,
    ! i.e. the array to be
    ! interpolated 

    INTEGER                        :: j, jj, k, kk
    INTEGER, DIMENSION(-1:2)       :: T_index
    INTEGER, DIMENSION(-1:2)       :: R_index
    REAL(wp),DIMENSION(-1:2,-1:2)  :: G  ! A Matrix of dimension
    ! NTau,Nre,used to store array only
    ! interpolated to current viewing
    ! geometry

    REAL(wp), DIMENSION(4) :: Y          ! A vector to contain the values of F at
    ! (iT0,iR0), (iT0,iR1), (iT1,iR1) and
    ! (iT1,iR0) respectively (i.e. anticlockwise
    ! from the bottom left)
    REAL(wp), DIMENSION(4) :: dYdTau     ! Gradients of F wrt Tau at the same points
    ! as Y
    REAL(wp), DIMENSION(4) :: dYdRe      ! Gradients of F wrt Re  at the same points
    ! as Y
    REAL(wp), DIMENSION(4) :: ddY        ! 2nd order cross derivatives of ! F wrt Tau
    ! and Re at the same points
    REAL(wp)               :: a1, a2, a3 ! Temporary store for output of BiCubic


    F = Grid%Albedo

    T_index(-1) = GZero%iTm1
    T_index( 0) = GZero%iT0 
    T_index( 1) = GZero%iT1 
    T_index( 2) = GZero%iTp1
    R_index(-1) = GZero%iRm1
    R_index( 0) = GZero%iR0 
    R_index( 1) = GZero%iR1 
    R_index( 2) = GZero%iRp1

    DO j = iXm1, iXp1
       jj = T_index(j)
       DO k = iXm1, iXp1
          kk = R_index(k)
          G(j,k) = (GZero%So1*F(jj,GZero%iSoZ0,kk))+&
               (GZero%dSoZ*F(jj,GZero%iSoZ1,kk))
       END DO
    END DO

    ! Calculte the function derivatives at four LUT points around our X
    Y(1) = G(iX0,iX0)
    Y(4) = G(iX0,iX1)
    Y(3) = G(iX1,iX1)
    Y(2) = G(iX1,iX0)

    ! Now call the adapted Numerical Recipes BCuInt subroutine to
    ! perform the interpolation to our desired state vector
    ! WRT to Tau

    dYdTau(1) = (G(iX1,iX0)-G(iXm1,iX0))/(Grid%Tau(GZero%iT1)-Grid%Tau(GZero%iTm1))
    dYdTau(2) = (G(iXp1,iX0)-G(iX0,iX0))/(Grid%Tau(GZero%iTp1)-Grid%Tau(GZero%iT0))
    dYdTau(3) = (G(iXp1,iX1)-G(iX0,iX1))/(Grid%Tau(GZero%iTp1)-Grid%Tau(GZero%iT0))
    dYdTau(4) = (G(iX1,iX1)-G(iXm1,iX1))/(Grid%Tau(GZero%iT1)-Grid%Tau(GZero%iTm1))

    ! WRT to Re
    dYDRe(1)  = (G(iX0,iX1)-G(iX0,iXm1))/(Grid%Re (GZero%iR1 )-Grid%Re (GZero%iRm1))
    dYDRe(2)  = (G(iX1,iX1)-G(iX1,iXm1))/(Grid%Re (GZero%iR1 )-Grid%Re (GZero%iRm1))
    dYDRe(3)  = (G(iX1,iXp1)-G(iX1,iX0))/(Grid%Re (GZero%iRp1)-Grid%Re (GZero%iR0 ))
    dYDRe(4)  = (G(iX0,iXp1)-G(iX0,iX0))/(Grid%Re (GZero%iRp1)-Grid%Re (GZero%iR0 ))

    ! Cross derivatives (dY^2/dTaudRe)
    ddY(1) = (G(iX1,iX1)-G(iX1,iXm1)-G(iXm1,iX1)+G(iXm1,iXm1))/&
         ((Grid%Tau(GZero%iT1 )-Grid%Tau(GZero%iTm1))*&
         (Grid%Re (GZero%iR1 )-Grid%Re (GZero%iRm1)))
    ddY(2) = (G(iXp1,iX1)-G(iXp1,iXm1)-G(iX0,iX1)+G(iX0,iXm1))/&
         ((Grid%Tau(GZero%iTp1)-Grid%Tau(GZero%iT0 ))*&
         (Grid%Re (GZero%iR1 )-Grid%Re (GZero%iRm1)))
    ddY(3) = (G(iXp1,iXp1)-G(iXp1,iX0)-G(iX0,iXp1)+G(iX0,iX0))/&
         ((Grid%Tau(GZero%iTp1)-Grid%Tau(GZero%iT0 ))*&
         (Grid%Re (GZero%iRp1)-Grid%Re (GZero%iR0 )))
    ddY(4) = (G(iX1,iXp1)-G(iX1,iX0)-G(iXm1,iXp1)+G(iXm1,iX0))/&
         ((Grid%Tau(GZero%iT1 )-Grid%Tau(GZero%iTm1))*&
         (Grid%Re (GZero%iRp1)-Grid%Re (GZero%iR0 )))

    CALL BCUINT(Y,dYdTau,dYdRe,ddY, &
         Grid%Tau(GZero%iT0), &
         Grid%Tau(GZero%iT1), &
         Grid%Re(GZero%iR0), &
         Grid%Re(GZero%iR1), &
         GZero%dT,GZero%dR,a1,a2,a3)

    cloud_cci_albedo = a1

  END FUNCTION ALBEDO

  SUBROUTINE READ_ALBEDO_LUT(cloud_cci_albedo,phase,options)

    ! Get the cloud albedo
    !
    ! The spectral Cloud albedo look up tables are from:
    !http://proj.badc.rl.ac.uk/orac/browser#sad_dir
    !
    !RD Luts have the following format
    !Channel
    !18 (optical depth dimension)
    !10 solar zenith dimension
    !23 radius dimension
    !Rd(18,10,23)
    !
    !

    IMPLICIT NONE

    TYPE(name_list), INTENT(in) :: options
    CHARACTER(len=3), INTENT(in):: phase

    ! out
    TYPE(cloud_albedo), INTENT(out) :: cloud_cci_albedo

    ! internal
    CHARACTER(len=1000):: pathLUT
    CHARACTER(len=1000):: sat
    INTEGER            :: lun, iostat
    REAL(wp)           :: dTau,dSolzen,dRe
    INTEGER            :: nTau,nSolzen,nRe 
    INTEGER, PARAMETER :: lrgDim = 40 ! some large dimension length for allocating
    REAL(wp)           :: tmp_val(lrgDim)
    REAL(wp)           :: wavelength

    CHARACTER(1000)    :: dir
    ! -------------------------------
    !        FILE NAME
    !
    ! fixme: We need to get enable this for other sensors too.

    SELECT CASE (TRIM(options%L2b%satellite))
    CASE ('noaa6','noaa7','noaa8','noaa9','noaa10','noaa11','noaa12',&
         'noaa13','noaa14','noaa15','noaa16','noaa17','noaa18','noaa19')
       sat = TO_UPPER(TRIM(options%L2b%satellite))
       IF (sat.EQ.'NOAA6') sat="NOAA6 doesn't have any albedo files"
       IF (sat.EQ.'NOAA7') sat='NOAA07'
       IF (sat.EQ.'NOAA8') sat='NOAA08'
       IF (sat.EQ.'NOAA9') sat='NOAA09'
    CASE ('metopa')
       sat='Metop1'
    CASE ('metopb')
       sat = 'Metop2'
    CASE DEFAULT
       PRINT '(a,a,a)', "satellite: '",sat,&
            "' is not in the list"
    END SELECT

    WRITE(pathLUT,'(a,"microphysics/Cloud_cci/cloud_albedo/AVHRR-",a,"_",a3,"_RD_Ch1.sad.txt")') &
         TRIM(options%paths%data_dir),&
         TRIM(sat), &
         phase

    WRITE(dir,'(a,a)') TRIM(options%paths%data_dir),"microphysics/Cloud_cci/cloud_albedo/"
    pathLUT = BUILD_FILENAME(dir=dir,&
         formstr='#DS-#SAT_#STRING_RD_Ch1.sad.txt',&
         dataset='AVHRR',sat=sat,string=phase)

    ! -------------------------------
    !       READ THE DATA
    !

    lun = GET_LUN(100)

    OPEN(lun, file = pathLUT, status='old', iostat = iostat)

    IF (iostat .NE. 0) &
       WRITE(*,*) 'ERROR: cloud_cci_functions(): Error opening file: ',TRIM(pathLUT)

    ! I'm not using this variable
    READ(lun, *, iostat=iostat) wavelength

    ! ---------------
    ! Read the grid dimensions
    !

    tmp_val(1:lrgDim) = missing
    CALL GRID_DIMENSION_READ(pathLUT, lun, 'nTau', 'dTau', 'Tau', &
         nTau, dTau, tmp_val)
    ALLOCATE( cloud_cci_albedo%Tau(nTau) )
    cloud_cci_albedo%Tau(1:nTau) = tmp_val(1:nTau)

    ! Read the Solzen dimension
    tmp_val(1:lrgDim) = missing
    CALL GRID_DIMENSION_READ(pathLUT, lun, 'nSolzen', 'dSolzen', 'Solzen', &
         nSolzen, dSolzen, tmp_val)
    ALLOCATE( cloud_cci_albedo%Solzen(nSolzen) )
    cloud_cci_albedo%Solzen(1:nSolzen) = tmp_val(1:nSolzen)

    ! Read the Re dimension
    tmp_val(1:lrgDim) = missing
    CALL GRID_DIMENSION_READ(pathLUT, lun, 'nRe', 'dRe', 'Re', &
         nRe, dRe, tmp_val)
    ALLOCATE( cloud_cci_albedo%Re(nRe) )
    cloud_cci_albedo%Re(1:nRe) = tmp_val(1:nRe)

    ! Read in the LUT array
    ALLOCATE(cloud_cci_albedo%Albedo (nTau,nSolzen,nRe) )

    cloud_cci_albedo%Albedo = missing
    CALL READ_VALUES_3D(pathLUT, lun, &
         nTau, nSolzen, nRe,    &
         cloud_cci_albedo%Albedo)

    CLOSE(lun)

    cloud_cci_albedo%nTau   = nTau
    cloud_cci_albedo%nSolzen= nSolzen
    cloud_cci_albedo%nRe    = nRe
    cloud_cci_albedo%dTau   = dTau
    cloud_cci_albedo%dSolzen= dSolzen
    cloud_cci_albedo%dRe    = dRe

  END SUBROUTINE READ_ALBEDO_LUT

END MODULE CLOUD_CCI_FUNCTIONS
