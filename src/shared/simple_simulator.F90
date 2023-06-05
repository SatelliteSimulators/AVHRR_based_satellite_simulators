MODULE simple_ctth_simulator
  !
  ! A simple simulator the finds the top of the cloud by going some
  ! optical depth down into the cloud. This does not check the cloud
  ! mask (e.g. based on tau_min). It is assumed that this function is
  ! only called if there is a cloud
  !
  ! Salomon.Eliasson@smhi.se

  USE simulator_input_variables, ONLY: subset
  USE internal_simulator,        ONLY: internal
  USE cosp_kinds,                ONLY: wp

  IMPLICIT NONE

  PUBLIC :: cth_simple,&
       ctp_simple,&
       ctt_simple

  PRIVATE :: simple_ctth

CONTAINS

  FUNCTION cth_simple(d1,ins,nlev,sub,inter,tau_eqCT)

    IMPLICIT NONE
    INTEGER, INTENT(in)       :: d1,ins,nlev
    TYPE(subset), INTENT(in)  :: sub
    TYPE(internal),INTENT(in) :: inter
    REAL(wp), INTENT(in)      :: tau_eqCT ! limiters

    ! out
    REAL(wp) :: cth_simple

    CALL simple_ctth(d1,ins,nlev,sub,inter,tau_eqCT, &
         cth = cth_simple)

  END FUNCTION cth_simple

  FUNCTION ctp_simple(d1,ins,nlev,sub,inter,tau_eqCT)

    IMPLICIT NONE
    INTEGER, INTENT(in)       :: d1,ins,nlev
    TYPE(subset), INTENT(in)  :: sub
    TYPE(internal),INTENT(in) :: inter
    REAL(wp), INTENT(in)      :: tau_eqCT ! limiters

    ! out
    REAL(wp) :: ctp_simple

    CALL simple_ctth(d1,ins,nlev,sub,inter,tau_eqCT, &
         ctp = ctp_simple)

  END FUNCTION ctp_simple

  FUNCTION ctt_simple(d1,ins,nlev,T,sub,inter,tau_eqCT)

    IMPLICIT NONE
    INTEGER, INTENT(in)       :: d1,ins,nlev
    TYPE(subset), INTENT(in)  :: sub
    TYPE(internal),INTENT(in) :: inter
    REAL(wp), INTENT(in)      :: tau_eqCT ! limiters
    ! If you are looking for the cloud top temperature you need to
    ! supply a temperature profile (since it may need to be a corrected
    ! profile, or not)
    REAL(wp), INTENT(in)      :: T(nlev)
    ! out
    REAL(wp) :: ctt_simple

    CALL simple_ctth(d1,ins,nlev,sub,inter,tau_eqCT, &
         T = T,ctt = ctt_simple)

  END FUNCTION ctt_simple

  SUBROUTINE simple_ctth(d1,ins,nlev,sub,inter,tau_eqCT, &
       T, ctp, ctt, cth)

    !
    !
    ! Find the cloud top heights per grid box by finding the level
    ! where e.g., tau=1 (Holz et. al., 2006) integrated from the TOP
    ! DOWN (troposphere -> surface).
    !
    ! Salomon.Eliasson@smhi.se

    IMPLICIT NONE

    ! IN
    ! 'nlev'     = number of model layers
    ! 'tau_eqCT' = Cloud optical depth where the radiantly equivalent cloud
    !              top height is located. The MODIS satellite
    ! 'tau_prof'      = Calculated optical depth at 0.67 micron (in-cloud)
    ! 'P'        = pressure at level centers
    ! 'T'        = Simulated layer brightness temperature from model
    !              temperature and absorption by model layer water
    !              vapour at the model layer centers

    INTEGER, INTENT(in)       :: d1,ins,nlev
    TYPE(subset), INTENT(in)  :: sub
    TYPE(internal),INTENT(in) :: inter
    REAL(wp), INTENT(in)      :: tau_eqCT

    ! If you are looking for the cloud top temperature you need to
    ! supply a temperature profile (since it may need to be a corrected
    ! profile, or not)
    REAL(wp), INTENT(in), OPTIONAL :: T(nlev)

    ! OUT
    !
    ! 'ctp'     = Cloud top pressure
    ! 'ctt'     = Cloud top temperature

    REAL(wp),INTENT(out),OPTIONAL :: cth
    REAL(wp),INTENT(out),OPTIONAL :: ctp
    REAL(wp),INTENT(out),OPTIONAL :: ctt

    ! internal
    ! L2 = model level interface above the model layer where the
    !      total integrated optical depth is greater than tau_eqCT
    ! L1 = -||- below
    ! tau_CT_ratio = where in the layer the integrated optical depth =
    !              tau_eqCT (e.g. = (tau_eqCT-tau(L2)) / (tau(L1)-tau2(L2))
    !
    !

    INTEGER :: L2, L1, inl, trop_lev
    REAL(wp) :: tau2, tau1, tau_ct_ratio, tau_tot
    REAL(wp) :: epsR, h1,h2

    IF (.NOT. PRESENT(T) .AND. PRESENT(ctt)) &
         STOP "The temperature profile is needed to retrieve the CTT"
    epsR = EPSILON(1._wp)
    ! initialise variables
    tau2         = 0._wp
    tau_CT_ratio = 0._wp
    tau_tot      = SUM(inter%tau_profile(ins,1:nlev))
    L2 = 0
    L1 = 0
    trop_lev = -999

    ! Use the tropopause level to limit clouds to the troposphere.
    ! FIND_TEMPERATUREINVERSIONS() ensures that the tropopause is
    ! between 500-50hPa. But if the function was never called, use nlev

    IF (ANY(sub%inv_layers(d1,1:nlev) > 0)) THEN
       trop_lev = MINVAL(sub%inv_layers(d1,1:nlev),&
            mask=sub%inv_layers(d1,1:nlev) > 0)
    ELSE
       trop_lev = 1 ! top of the atmosphere
    END IF

    inl  = trop_lev
    tau1 = SUM(inter%tau_profile(ins,1:inl-1)) ! i.e., above the tropopause

    ! +epsR incase tau_eqCT is set to 0, i.e. catch the slightest of
    ! cloud tops (native model)

    IF (SUM(inter%tau_profile(ins,trop_lev:nlev)) .GE. tau_eqCT+epsR) THEN
       ! if the cloud is THICKER than the radiative cloud top height

       ! -----------
       ! Find the layer where the cloud is thicker than tau_eq_CT
       !
       DO WHILE (tau1 .LT. tau_eqCT+epsR)
          ! this level is the upper bounds of tau
          inl  = inl+1
          tau1 = SUM(inter%tau_profile(ins,trop_lev:inl) )
       END DO
       L1   = inl
       ! --------------
       ! Go back up the model layers to find the top tau_eqCT layer
       !
       tau2 = tau1
       DO WHILE (tau2 .GE. tau_eqCT+epsR)
          inl  = inl-1
          tau2 = SUM(inter%tau_profile(ins,trop_lev:inl))
       END DO
       L2 = MAXVAL([inl,trop_lev])

       tau_CT_ratio = (tau_eqCT-tau2)/(tau1-tau2)

       ! The pressures and Temperatures are at the layer centers, but
       ! the height is on interfaces
       IF (PRESENT(cth)) THEN
          h1=(sub%height(d1,L1)+sub%height(d1,L1+1))/2
          h2=(sub%height(d1,L2)+sub%height(d1,L2+1))/2
          cth = h2+tau_CT_ratio*(h1-h2)
       END IF
       IF (PRESENT(ctp)) THEN
          ctp = &
               EXP( LOG(sub%p_mid(d1,L2))+&
               tau_CT_ratio*(LOG(sub%p_mid(d1,L1)/sub%p_mid(d1,L2))))
       END IF
       IF (PRESENT(ctt)) THEN
          ctt = T(L2)+tau_CT_ratio*(T(L1)-T(L2))
       END IF

    ELSEIF (SUM(inter%tau_profile(ins,trop_lev:nlev)) .LT. tau_eqCT+epsR) THEN

       ! Then there is a cloud that is GREATER than the detection
       ! limit, yet THINNER than the level at what the radiative cloud top
       ! height is found (for MODIS it is tau=1)

       ! find the bottom of the cloud by looping up from
       !  the surface
       inl = nlev
       DO WHILE(inter%tau_profile(ins,inl) .EQ. 0)
          inl=inl-1
       END DO

       h1=(sub%height(d1,inl)+sub%height(d1,inl+1))/2
       IF (PRESENT(cth)) cth = h1
       IF (PRESENT(ctp)) ctp = sub%p_mid(d1,inl)
       IF (PRESENT(ctt)) ctt = T(inl)

    END IF ! if cloud is transparent

  END SUBROUTINE simple_ctth

END MODULE simple_ctth_simulator
