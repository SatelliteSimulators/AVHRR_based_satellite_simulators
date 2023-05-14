MODULE data_check
  ! various routines dedicated to checking the validity of data. Good
  ! for debugging and runtime error checking
  !
  ! Salomon.Eliasson@smhi.se

  USE cosp_kinds,                ONLY: wp
  USE internal_simulator,        ONLY: internal
  USE model_input,               ONLY: model_type
  USE namelist_input,            ONLY: name_list
  USE optics_m,                  ONLY: optics_LUT
  USE simulator_variables,       ONLY: satellite_simulator
  USE simulator_input_variables, ONLY: subset


  IMPLICIT NONE

  PUBLIC  :: check_variables,  &
       common_check_subgrid,   &
       check_var,              &
       print_sub_point

  PRIVATE :: check_model,&
       check_derived

  REAL(wp), PARAMETER :: MinimumTropopausePressure = 50 ![Pa] CLARA-A2 has no limit on the tropopause
  REAL(wp), PARAMETER :: height_min = 0,        height_max  = 40e4
  REAL(wp), PARAMETER :: P_min = 0,             P_max    = 1100e2
  REAL(wp), PARAMETER :: T_min = 123.15,        T_max    = 353.15
  REAL(wp), PARAMETER :: cloud_low = 680e2,     cloud_mid = 440e2 ![hPa]
  REAL(wp), PARAMETER :: cth_min  = height_min, cth_max  = height_max
  REAL(wp), PARAMETER :: ctp_min  = MinimumTropopausePressure, ctp_max  = P_max
  REAL(wp), PARAMETER :: ctt_min  = T_min,      ctt_max  = T_max
  REAL(wp), PARAMETER :: epsR = 1000*EPSILON(1._wp)
  REAL(wp), PARAMETER :: frac_min = 0,          frac_max = 1
  REAL(wp), PARAMETER :: missing = -999._wp
  REAL(wp), PARAMETER :: p_tau_min=0,           p_tau_max=4.e7
  REAL(wp), PARAMETER :: tau_min = 0,           tau_max  = 4000 ! 200 would be reasonable...

CONTAINS

  SUBROUTINE CHECK_VARIABLES(ngrids,nlev,options,data_mask,&
       M,S,print_profile,scops_in)
    !
    ! Checks that the model and other input data to the simulator are
    ! within limits
    !
    IMPLICIT NONE

    TYPE(name_list), INTENT(in) :: options
    INTEGER, INTENT(in)         :: ngrids,nlev

    ! only take data where the data is not masked, i.e.,
    ! data_mask=.false., into consideration
    LOGICAL,INTENT(in)          :: data_mask(ngrids,nlev)

    TYPE(model_type), INTENT(in),OPTIONAL :: M
    TYPE(subset), INTENT(in),OPTIONAL     :: S
    LOGICAL, OPTIONAL, INTENT(in)         :: print_profile
    REAL(wp), OPTIONAL, INTENT(in)        :: scops_in(nlev)

    INTEGER :: d1,inl
    REAL(wp), ALLOCATABLE :: tmp_long(:)
    REAL(wp) :: scops(nlev)
    LOGICAL :: whammy(ngrids,nlev)

    IF (PRESENT(print_profile)) THEN
       whammy = .NOT.data_mask
    ELSE
       whammy = .FALSE.
    END IF
    IF (PRESENT(scops_in)) THEN
       scops = scops_in
    ELSE
       scops(1:nlev) = 1._wp
    END IF

    IF (PRESENT(M)) THEN
       whammy = whammy .OR. CHECK_MODEL(M,data_mask,ngrids,nlev)
    END IF

    IF (PRESENT(S)) THEN
       whammy = whammy .OR. CHECK_DERIVED(S,data_mask,options,ngrids,nlev)
    END IF
    IF (ANY(whammy)) THEN
       DO d1=1,ngrids
          IF (ANY(whammy(d1,:))) THEN
             PRINT *, " "
             PRINT '(a,1x,I4)',&
                  "These are all the variables at ngrid=",d1
             IF (PRESENT(scops_in)) PRINT *, "Note: SUBSETTED USING SCOPS!!"
             IF (PRESENT(M)) THEN
                PRINT *, "Model variables:"
14              FORMAT (A3,1x,(a2,1x),3(A11,1x),A13,1x,A8,1x,A3)
                WRITE(*,14) "lev","CV","CC","CIWC [g/kg]","CLWC [g/kg]",&
                     "Q [g/kg]","T","lev"
                DO inl=1,nlev
                   PRINT '(I3,1x,(i2,1x),3(F11.3,1x),F13.5,1x,F8.3,1x,I3)',&
                        inl,                           &
                        0,                             & ! CV is set to zero
                        scops(inl)*M%CC  (d1,inl),     &
                        scops(inl)*M%CIWC(d1,inl)*1000,&
                        scops(inl)*M%CLWC(d1,inl)*1000,&
                        M%Q(d1,inl)*1000,              &
                        M%T(d1,inl),                   &
                        inl
                END DO
                WRITE(*,14) "lev","CV","CC","CIWC [g/kg]","CLWC [g/kg]",&
                     "Q [g/kg]","T","lev"
                PRINT *,""

                PRINT '(A13,F8.3)',"M%CI      = ",M%CI   (d1)
                PRINT '(A13,F7.2)',"M%PSURF   = ",M%PSURF(d1)/100
                PRINT '(A13,F7.2)',"M%SKT     = ",M%SKT  (d1)
                PRINT '(A13,F7.2)',"M%T2M     = ",M%T2M  (d1)
                PRINT '(A13,F7.2)',"M%TCC     = ",M%TCC  (d1)
                PRINT '(A13,F8.3)',"M%TCWV    = ",M%TCWV (d1)
             END IF

             IF (PRESENT(S)) THEN
                CALL PRINT_SUB_POINT(S,d1,nlev,scops(1:nlev))

             END IF
             IF (.NOT. PRESENT(print_profile)) &
                  STOP "Model or simulator input fields are outside the valid range"
          END IF
       END DO
    END IF
    IF (ALLOCATED(tmp_long)) DEALLOCATE(tmp_long)
  END SUBROUTINE CHECK_VARIABLES
  !
  ! -----
  !
  FUNCTION CHECK_MODEL(M,data_mask,ngrids,nlev) RESULT(whammy)

    IMPLICIT NONE

    TYPE(model_type), INTENT(in) :: M
    ! only take data where the data is not masked, i.e.,
    ! data_mask=.false., into consideration
    INTEGER, INTENT(in) :: ngrids,nlev
    LOGICAL,INTENT(in)  :: data_mask(ngrids,nlev)
    REAL(wp), PARAMETER :: ciwc_min= 0,ciwc_max =1
    REAL(wp), PARAMETER :: clwc_min= 0,clwc_max =1
    REAL(wp), PARAMETER :: psurf_min=30000, psurf_max=110000
    REAL(wp), PARAMETER :: q_min   = 0,q_max   = 0.4129_wp
    REAL(wp), PARAMETER :: tcwv_min= 0,tcwv_max=100
    REAL(wp), PARAMETER :: lsm_min = 0,lsm_max = 1

    LOGICAL :: whammy(ngrids,nlev),dm(ngrids,nlev),dm1(ngrids)
    CHARACTER(30) :: var
    CHARACTER(5) :: form

    dm = data_mask
    dm1=dm(:,1)

    ! local check of one variable
    whammy   = .FALSE.

    var='M%CC'
    form='f5.3'
    whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,&
         frac_min,frac_max,mask2=dm,data2=M%cc)

    var='M%CI'
    form='f5.3'
    whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,&
         frac_min,frac_max,mask1=dm1,data1=M%ci)

    var='M%CIWC'
    form='f5.3'
    whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,&
         ciwc_min,ciwc_max,mask2=dm,data2=M%ciwc)

    var='M%CLWC'
    form='f5.3'
    whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,&
         clwc_min,clwc_max,mask2=dm,data2=M%clwc)

    var='M%PSURF'
    form='f7.2'
    whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,&
         psurf_min,psurf_max,mask1=dm1,data1=M%PSURF,scale=0.01_wp)

    var='M%Q'
    form='f6.4'
    whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,&
         q_min,q_max,mask2=dm,data2=M%Q)

    var='M%SKT'
    form='f7.2'
    whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,&
         T_min,T_max,mask1=dm1,data1=M%SKT)

    var='M%T'
    form='f7.2'
    whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,&
         T_min,T_max,mask2=dm,data2=M%T)

    var='M%T2M'
    form='f7.2'
    whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,&
         T_min,T_max,mask1=dm1,data1=M%T2M)

    var='M%TCC'
    form='f5.3'
    whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,&
         frac_min,frac_max,mask1=dm1,data1=M%TCC)

    var='M%TCWV'
    form='f7.3'
    whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,&
         tcwv_min,tcwv_max,mask1=dm1,data1=M%TCWV)

    var='M%CV'
    form='f5.3'
    whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,&
         frac_min,frac_max,mask2=dm,data2=M%CV)


  END FUNCTION CHECK_MODEL
  !
  ! -----
  !
  FUNCTION CHECK_DERIVED(S,data_mask,options,ngrids,nlev) RESULT(whammy)

    IMPLICIT NONE

    TYPE(subset), INTENT(in)   :: S
    TYPE(name_list),INTENT(in) :: options
    ! only take data where the data is not masked, i.e.,
    ! data_mask=.false., into consideration
    LOGICAL,INTENT(in)       :: data_mask(ngrids,nlev)
    INTEGER, INTENT(in)      :: ngrids,nlev

    REAL(wp), PARAMETER :: solzen_min  = 0, solzen_max  = 180
    REAL(wp), PARAMETER :: ciwc_gm3_min= 0, ciwc_gm3_max= 3e3
    REAL(wp), PARAMETER :: clwc_gm3_min= 0, clwc_gm3_max= 3e3
    REAL(wp)            :: lreff_min,lreff_max,lwp_min,lwp_max
    REAL(wp)            :: ireff_min,ireff_max,iwp_min,iwp_max
    REAL(wp),PARAMETER  :: ZRefDe = 0.64952!Ratio of effective radius to diameter
    LOGICAL :: whammy(ngrids,nlev)
    LOGICAL :: dm(ngrids,nlev),dm1(ngrids)
    CHARACTER(30) :: var
    CHARACTER(5)  :: form

    ! local check of one variable
    whammy   = .FALSE.

    ! -------
    ! This is the models' own effective radius
    lreff_min = 4
    lreff_max = 16
    lwp_min   = tau_min*lreff_min
    lwp_max   = tau_max*lreff_max
    ireff_min = 20*ZrefDe
    ireff_max = 155*ZrefDe
    iwp_min   = tau_min*ireff_min
    iwp_max   = tau_max*ireff_max

    dm=data_mask
    dm1=dm(:,1)

    IF (ALLOCATED(S%cloud_emis)) THEN
       var='S%cloud_emis'
       form='f5.3'
       whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,&
            0._wp,1._wp,mask2=dm,data2=S%cloud_emis)
    END IF

    IF (ALLOCATED(S%g0)) THEN

       var='S%g0'
       form='f6.4'
       whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,&
            0._wp,1._wp,mask2=dm,data2=S%g0,except=missing)
    END IF
    IF (ALLOCATED(S%height)) THEN
       var='S%height'
       form='f6.0'
       ! valid for model interfaces
       whammy=whammy.OR.SPREAD(ANY(CHECK_VAR(var,form,ngrids,nlev+1,&
            height_min,height_max,mask2=dm,data2=S%height),dim=2),2,nlev)
    END IF
    var='S%ireff'
    form='f7.2'
    whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,ireff_min,ireff_max,&
         mask2=dm,data2=S%ireff,except=missing)

    var='S%itau'
    form='f8.2'
    whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,tau_min,tau_max,&
         mask2=dm,data2=S%itau)

    var='S%iwp'
    form='f8.2'
    whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,iwp_min,iwp_max,&
         mask2=dm,data2=S%iwp,scale=1000._wp)

    var='S%lreff'
    form='f7.2'
    whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,lreff_min,lreff_max,&
         mask2=dm,data2=S%lreff,except=missing)

    var='S%ltau'
    form='f9.2'
    whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,tau_min,tau_max,&
         mask2=dm,data2=S%ltau)

    var='S%lwp'
    form='f9.2'
    whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,lwp_min,lwp_max,&
         mask2=dm,data2=S%lwp,scale=1000._wp)

    var='S%p_mid'
    form='f7.1'
    whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,p_min,p_max,&
         mask2=dm,data2=S%p_mid,scale=0.01_wp)

    var='S%solzen'
    form='f6.1'
    whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,solzen_min,solzen_max,&
         mask1=dm1,data1=S%solzen)

    IF (options%sim%doClara .OR. options%sim%doCloud_cci) THEN
       var='S%Tcorr'
       form='f7.2'
       whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,T_min,T_max,&
            mask2=dm,data2=S%Tcorr)
    END IF

    IF (ALLOCATED(S%w0)) THEN

       var='S%w0'
       form='f6.4'
       whammy=whammy.OR.CHECK_VAR(var,form,ngrids,nlev,0._wp,1._wp,&
            mask2=dm,data2=S%w0,except=missing)
    END IF


  END FUNCTION CHECK_DERIVED

  SUBROUTINE PRINT_SUB_POINT(S,d1,nlev,sc)

    ! because sometimes I want to call this function by itself when
    ! debugging somewhere else

    IMPLICIT NONE
    TYPE(subset), INTENT(in) :: S
    INTEGER,      INTENT(in) :: d1,nlev
    REAL(wp),OPTIONAL,INTENT(in) :: sc(nlev)
    REAL(wp), DIMENSION(nlev):: inversion, scops
    CHARACTER(9)             :: arrow
    INTEGER :: inl,inv_ind

    IF (PRESENT(sc)) THEN
       scops=sc
    ELSE
       scops = 1._wp
    END IF

15  FORMAT (A4,1x, (A6),   9(A8,1x),  a9)
16  FORMAT (I4,1x, (I5,1x),9(F8.3,1x),a9)
    PRINT *, "Subset:"
    IF (.NOT.ALL(scops.GT.0))&
         PRINT *,"extra simulator input variables (weighted by var*scops): "

    ! deal with inversions
    inv_ind = 1
    inversion(1:nlev) = -999
    DO inl = nlev,1,-1
       IF (S%inv_layers(d1,inl) .GT. 0) THEN
          inversion(inv_ind) = S%inv_layers(d1,inl)
          inv_ind = inv_ind+1
       END IF
    END DO
    inv_ind = 1

    WRITE(*,15) "lev:","ciwc","clwc",&
         "height","ireff","itau",&
         "iwp","lreff","ltau",&
         "lwp","p_int","p_mid","Tcorr","Inversion"
    WRITE(*,15) "[-]","[g/m3]","[g/m3]","[m]","[micron]","[-]",&
         "[g/m2]","[micron]","[-]","[g/m2]","[hPa]","[hPa]","[K]","[-]"
    DO inl = 1, nlev
       arrow       = ""
       IF (inversion(inv_ind) .EQ. inl) THEN
          arrow = '<-----'
          inv_ind = inv_ind+1
       END IF


       WRITE(*,16) inl,&
            INT((S%height(d1,inl)+S%height(d1,inl))/2),& !m
            S%ireff     (d1,inl)*scops(inl)    ,& !micron
            S%itau      (d1,inl)*scops(inl)    ,&
            S%iwp       (d1,inl)*1e3*scops(inl),& !g/m2
            S%lreff     (d1,inl)*scops(inl)    ,&
            S%ltau      (d1,inl)*scops(inl)    ,&
            S%lwp       (d1,inl)*1e3*scops(inl),& !g/m2
            S%p_int     (d1,inl)/100           ,&
            S%p_mid     (d1,inl)/100           ,&
            S%Tcorr     (d1,inl)               ,&
            arrow
    END DO
    WRITE(*,15) "lev:",&
         "height","ireff","itau",&
         "iwp","lreff","ltau",&
         "lwp","p_int","p_mid",&
         "Tcorr","inversion"
    WRITE(*,15) "[-]","[m]","[micron]","[-]",&
         "[g/m2]","[micron]","[-]","[g/m2]","[hPa]","[hPa]","[K]","[-]"

    PRINT '(A13,1x,F7.1," [deg]")',"S%solzen =",S%solzen(d1)
  END SUBROUTINE PRINT_SUB_POINT

  FUNCTION COMMON_CHECK_SUBGRID(I,ncol,nlev) RESULT(whammy)

    TYPE(internal), INTENT(in)  :: I
    INTEGER, INTENT(in)         :: ncol,nlev
    LOGICAL,DIMENSION(ncol,nlev):: whammy
    CHARACTER(5)                :: form
    CHARACTER(20)               :: var

    whammy = .FALSE.

    var='I%cflag'
    form='i1'
    whammy = whammy.OR.CHECK_VAR(var,form,ncol,nlev,&
         mini=0,maxi=3,data1I=I%cflag)

    var='I%tau_profile'
    form='f9.2'
    whammy = whammy.OR.CHECK_VAR(var,form,ncol,nlev,&
         min=tau_min,max=tau_max,data2=I%tau_profile)

    var='I%cph'
    form='i1'
    whammy = whammy.OR.CHECK_VAR(var,form,ncol,nlev,&
         mini=0,maxi=2,data1I=I%cph)

    var='I%tau'
    form='f9.2'
    whammy = whammy.OR.CHECK_VAR(var,form,ncol,nlev,&
         min=tau_min,max=tau_max,data1=I%tau)

  END FUNCTION COMMON_CHECK_SUBGRID
  !
  ! -------
  !
  FUNCTION CHECK_VAR(var,form,ndim1,ndim2,&
       min,max,data1,mask1,data2,mask2,&
       mini,maxi,data1I,data2I,&
       scale,except) RESULT(whammy)

    IMPLICIT NONE

    CHARACTER(*),INTENT(in)       :: var
    CHARACTER(*),INTENT(in)       :: form ! string format
    INTEGER,  INTENT(in)          :: ndim1,ndim2 ! I can't avoid the second dimension
    REAL(wp), INTENT(in), OPTIONAL:: min,max
    INTEGER,  INTENT(in), OPTIONAL:: mini,maxi
    REAL(wp), INTENT(in), OPTIONAL:: data1(ndim1)
    INTEGER,  INTENT(in), OPTIONAL:: data1I(ndim1)
    REAL(wp), INTENT(in), OPTIONAL:: data2(ndim1,ndim2)
    INTEGER,  INTENT(in), OPTIONAL:: data2I(ndim1,ndim2)
    LOGICAL,  INTENT(in), OPTIONAL:: mask1(ndim1)
    LOGICAL,  INTENT(in), OPTIONAL:: mask2(ndim1,ndim2)
    REAL(wp), INTENT(in), OPTIONAL:: scale
    REAL(wp), INTENT(in), OPTIONAL:: except

    ! internal
    LOGICAL                       :: m1(ndim1)
    LOGICAL                       :: m2(ndim1,ndim2)
    INTEGER                       :: d1,d
    REAL(wp)                      :: tmp(1) ,tmp2 (ndim2)
    INTEGER                       :: tmpi(1),tmp2i(ndim2)
    LOGICAL                       :: whammy(ndim1,ndim2),isint
    CHARACTER(1000) :: strF,loopF
    ! Quiet NAN, double precision.
    REAL(wp), PARAMETER :: D_QNAN = &
         TRANSFER((/ Z'00000000', Z'7FF80000' /),1.0_wp)
    REAL(wp):: ex,sc

    isint=.FALSE.
    IF (PRESENT(maxi )) THEN
       isint=.TRUE.
    END IF
    IF (PRESENT(data1).OR.PRESENT(data1I))d=1
    IF (PRESENT(data2).OR.PRESENT(data2I))d=2

    ! mask
    m1=.FALSE.
    IF (PRESENT(mask1)) m1=mask1
    m2=.FALSE.
    IF (PRESENT(mask2)) m2=mask2

    ! scale
    sc = 1._wp
    IF (PRESENT(scale)) sc=scale

    ! except this value
    ex = D_QNAN
    IF (PRESENT(except)) ex = except

    whammy = .FALSE.

    WRITE(strF,*) '("variable: ''',TRIM(var),&
         &''' should be between min = ",',TRIM(form),&
         &'," and max = ",',TRIM(form),&
         &',". ''',TRIM(var),''' =",)'

    IF (ANY(whammy)) THEN

       SELECT CASE(d)
       CASE(1)
          IF (isint) THEN
             whammy(:,1) = .NOT.m1 .AND.&
                  (data1i.LT.mini-epsR.OR.data1i.GT.maxi+epsR).AND.&
                  data1i.NE.ex
          ELSE
             whammy(:,1) = .NOT.m1 .AND.&
                  (data1.LT.min-epsR.OR.data1.GT.max+epsR).AND.&
                  data1.NE.ex
          END IF
          WRITE(loopF,*)'("',TRIM(var),'(",I4,",",I4,")=",',TRIM(form),',1x)'
       CASE(2)
          IF (isint) THEN
             whammy = .NOT.m2 .AND.&
                  (data2i.LT.mini-epsR.OR.data2i.GT.maxi+epsR).AND.&
                  data2i.NE.ex
          ELSE
             whammy = .NOT.m2 .AND.&
                  (data2.LT.min-epsR.OR.data2.GT.max+epsR).AND.&
                  data2.NE.ex
          END IF
          WRITE(loopF,*)'("',TRIM(var),'(",I4,")=",',TRIM(form),',1x)'
       END SELECT

       IF (isint) THEN
          PRINT strF,mini*sc,maxi*sc
       ELSE
          PRINT strF,min*sc,max*sc
       END IF

       IF (isint) THEN
          FAIL:DO d1 = 1,ndim1
             IF (ANY(whammy(d1,:))) THEN
                SELECT CASE(d)
                CASE (1)
                   tmpi = data1i(d1)
                CASE (2)
                   tmp2i = data2i(d1,1:ndim2)
                   PRINT loopF,d1,tmp2i*sc
                END SELECT
                EXIT FAIL
             END IF
          END DO FAIL
       ELSE
          FAIL2:DO d1 = 1,ndim1
             IF (ANY(whammy(d1,:))) THEN
                SELECT CASE(d)
                CASE (1)
                   tmp = data1(d1)
                   PRINT loopF,d1,tmp*sc
                CASE (2)
                   tmp2  = data2(d1,1:ndim2)
                   PRINT loopF,d1,tmp2*sc
                END SELECT
                EXIT FAIL2
             END IF
          END DO FAIL2
       END IF

    END IF ! whammy

  END FUNCTION CHECK_VAR
END MODULE DATA_CHECK
