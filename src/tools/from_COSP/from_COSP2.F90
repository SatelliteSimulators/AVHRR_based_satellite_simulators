MODULE from_COSP2
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Copyright (c) 2015, Regents of the University of Colorado
  ! All rights reserved.
  !
  ! Redistribution and use in source and binary forms, with or without modification, are 
  ! permitted provided that the following conditions are met:
  !
  ! 1. Redistributions of source code must retain the above copyright notice, this list of 
  !    conditions and the following disclaimer.
  !
  ! 2. Redistributions in binary form must reproduce the above copyright notice, this list
  !    of conditions and the following disclaimer in the documentation and/or other 
  !    materials provided with the distribution.
  !
  ! 3. Neither the name of the copyright holder nor the names of its contributors may be 
  !    used to endorse or promote products derived from this software without specific prior
  !    written permission.
  !
  ! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
  ! EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
  ! MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL 
  ! THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
  ! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT 
  ! OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
  ! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
  ! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  ! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  !
  ! History
  ! May 2009:      Robert Pincus - Initial version
  ! June 2009:     Steve Platnick and Robert Pincus - Simple radiative transfer for size 
  !                retrievals
  ! August 2009:   Robert Pincus - Consistency and bug fixes suggested by Rick Hemler (GFDL) 
  ! November 2009: Robert Pincus - Bux fixes and speed-ups after experience with Rick Hemler 
  !                using AM2 (GFDL) 
  ! January 2010:  Robert Pincus - Added high, middle, low cloud fractions
  ! May 2015:      Dustin Swales - Modified for COSPv2.0
  ! Feb 2016:      Salomon.Eliasson@smhi.se - borrowed from https://github.com/CFMIP/COSPv2.0
  !
  ! Salomon.Eliasson@smhi.se

  USE COSP_KINDS,      ONLY: wp
  USE optics_m,        ONLY: &
       num_trial_res,        &
       optics_LUT,           &
       cloud_phase

  IMPLICIT NONE

  PUBLIC :: adding_doubling,    &
       compute_toa_reflectance, &
       get_g_nir,               &
       get_ssa_nir,             &
       interpolate_to_min,      &
       trial_g_and_w0,          & 
       two_stream,              &
       two_stream_reflectance

  !  INTEGER, PARAMETER ::               num_trial_res = 20,     & ! Increase to make the linear pseudo-retrieval of size 
  INTEGER, PARAMETER ::    phaseIsLiquid = 1,      & !
       phaseIsIce    = 2         !

  REAL(wp), PARAMETER :: re_fill = -999.0

CONTAINS

  ! #######################################################
  ! Make trial datasets of g and ssa
  ! #######################################################
  SUBROUTINE trial_g_and_w0(r_eff,is_ch3B)

    IMPLICIT NONE

    !EvMb
    INTEGER                         :: i
    !EvMe
    LOGICAL, INTENT(in)             :: is_ch3B ! if ch3B was on (else ch3a was on)
    TYPE(optics_LUT), INTENT(inout) :: r_eff

    REAL(wp), DIMENSION(num_trial_res) :: trial_re_w, trial_re_i

    trial_re_w = r_eff%water%optics%re%trial
    trial_re_i = r_eff%ice%optics%re%trial

    ! Water
    !EvMb ... next four modifications needed to overcome a bug in the gfortran compiler installed at KNMI 
    !   r_eff%water%optics%g0(1:num_trial_res) = get_g_nir(trial_re_w(1:num_trial_res),&
    !        r_eff%water,is_ch3B)
    do i=1,num_trial_res
       r_eff%water%optics%g0(i) = get_g_nir(trial_re_w(i),r_eff%water,is_ch3B)
    enddo

    !   r_eff%water%optics%w0(1:num_trial_res) = get_ssa_nir(trial_re_w(1:num_trial_res),&
    !        r_eff%water,is_ch3B)
    do i=1,num_trial_res
       r_eff%water%optics%w0(i) = get_ssa_nir(trial_re_w(i),r_eff%water,is_ch3B)
    enddo

    ! Ice
    !   r_eff%ice%optics%g0(1:num_trial_res)   = get_g_nir(trial_re_i(1:num_trial_res),&
    !         r_eff%ice,is_ch3B)
    do i=1,num_trial_res
       r_eff%ice%optics%g0(i)   = get_g_nir(trial_re_i(i),r_eff%ice,is_ch3B)
    enddo

    !   r_eff%ice%optics%w0(1:num_trial_res)   = get_ssa_nir(trial_re_i(1:num_trial_res),&
    !        r_eff%ice,is_ch3B)
    do i=1,num_trial_res
       r_eff%ice%optics%w0(i)   = get_ssa_nir(trial_re_i(i),r_eff%ice,is_ch3B)
    enddo
    !EvMe

  END SUBROUTINE trial_g_and_w0

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! FUNCTION get_g_nir
  ! Compute asymmetry parameter using provided radius and phase.
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ELEMENTAL FUNCTION get_g_nir(re_trial,phase,is_ch3B)
    ! INPUTS
    TYPE(cloud_phase), INTENT(in) :: phase
    LOGICAL, INTENT(in) :: is_ch3b

    ! Outputs
    REAL(wp) :: get_g_nir

    ! Local variables
    REAL(wp),INTENT(in) :: re_trial
    INTEGER  :: xi(2)
    INTEGER  :: nRe
    REAL(wp), ALLOCATABLE :: g0(:), re(:)
    REAL(wp) :: re_max,re_min

    ! list of trial effective radius's
    re_max    = phase%optics%re%max
    re_min    = phase%optics%re%min
    nRe       = phase%optics%re%nRe

    ALLOCATE(re(nRe))
    ALLOCATE(g0(nRe))
    re        = phase%optics%re%re
    IF (is_ch3B) THEN
       g0(1:nRe) = phase%optics%re%g0_37(1:nRe)
    ELSE
       g0(1:nRe) = phase%optics%re%g0_16(1:nRe)
    END IF

    ! Find interpolation bounds works for either cloud phase or wavelength
    xi = [MAXLOC(re-re_trial,re-re_trial .LE. 0),MAXLOC(re-re_trial,re-re_trial .LE. 0)+1]
    IF (MINVAL(ABS(re-re_trial)) .EQ. 0) xi=[xi(1),xi(1)]
    ! Interpolate
    IF (re_trial .GT. re_min .AND. re_trial .LT. re_max) THEN
       get_g_nir = g0(xi(1))+(g0(xi(2))-g0(xi(1)))*(re_trial-re(xi(1)))/    &
            (re(xi(2))-re(xi(1))) 
    ENDIF
    IF (xi(1) .EQ. xi(2)) get_g_nir = g0(xi(1))       ! On table node
    IF (re_trial .LT. re_min)  get_g_nir = g0(1)      ! Re too small
    IF (re_trial .GT. re_max)  get_g_nir = g0(nRe) ! Re too big

    DEALLOCATE(re)
    DEALLOCATE(g0)

  END FUNCTION get_g_nir
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! END FUNCTION get_g_nir
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! FUNCTION get_ssa_nirnir
  ! Compute single-scattering albedo for provided radius and phase.
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ELEMENTAL FUNCTION get_ssa_nir(re_trial,phase,is_ch3B)

    ! INPUTS
    TYPE(cloud_phase), INTENT(in) :: phase
    LOGICAL, INTENT(in) :: is_ch3b

    ! Outputs
    real(wp) :: get_ssa_nir

    ! Local variables
    REAL(wp),INTENT(in) :: re_trial
    INTEGER  :: xi(2)
    INTEGER  :: nRe
    REAL(wp), ALLOCATABLE :: w0(:), re(:)
    REAL(wp) :: re_max,re_min

    ! list of trial effective radius's
    re_max    = phase%optics%re%max
    re_min    = phase%optics%re%min
    nRe       = phase%optics%re%nRe
    ALLOCATE(re(nRe))
    ALLOCATE(w0(nRe))
    re        = phase%optics%re%re
    IF (is_ch3B) THEN
       w0(1:nRe) = phase%optics%re%w0_37(1:nRe)
    ELSE
       w0(1:nRe) = phase%optics%re%w0_16(1:nRe)
    END IF

    ! Find interpolation bounds works for either cloud phase or wavelength
    xi = [MAXLOC(re-re_trial,re-re_trial .LE. 0),MAXLOC(re-re_trial,re-re_trial .LE. 0)+1]
    IF (MINVAL(ABS(re-re_trial)) .EQ. 0) xi=[xi(1),xi(1)]
    ! Interpolate
    IF (re_trial .GT. re_min .AND. re_trial .LT. re_max) THEN
       get_ssa_nir = w0(xi(1))+(w0(xi(2))-w0(xi(1)))*(re_trial-re(xi(1)))/    &
            (re(xi(2))-re(xi(1))) 
    ENDIF
    IF (xi(1) .EQ. xi(2)) get_ssa_nir = w0(xi(1))      ! On table node
    IF (re_trial .LT. re_min)  get_ssa_nir = w0(1)   ! Re too small
    IF (re_trial .GT. re_max)  get_ssa_nir = w0(nRe) ! Re too big

    DEALLOCATE(re)
    DEALLOCATE(w0)

  END FUNCTION get_ssa_nir

  ! ########################################################################################
  ! Radiative transfer
  ! ########################################################################################
  PURE FUNCTION compute_toa_reflectance(nLevels,tau, g, w0)
    ! This wrapper reports reflectance only and strips out non-cloudy elements from the 
    ! calculation


    ! INPUTS
    INTEGER,INTENT(in)                     :: nLevels
    REAL(wp),INTENT(in),DIMENSION(nLevels) :: tau, g, w0
    ! OUTPUTS
    REAL(wp)                               :: compute_toa_reflectance
    ! LOCAL VARIABLES
    LOGICAL, DIMENSION(nLevels)                   :: cloudMask
    INTEGER, DIMENSION(COUNT(tau(1:nLevels) > 0)) :: cloudIndicies
    REAL(wp),DIMENSION(COUNT(tau(1:nLevels) > 0)) :: Refl,Trans
    REAL(wp)                                      :: Refl_tot, Trans_tot
    INTEGER                                       :: i

    cloudMask(1:nLevels) = tau(1:nLevels) > 0. 
    cloudIndicies = PACK((/ (i, i = 1, nLevels) /), mask = cloudMask) 
    DO i = 1, SIZE(cloudIndicies)
       CALL two_stream(tau(cloudIndicies(i)), g(cloudIndicies(i)), w0(cloudIndicies(i)), Refl(i), Trans(i))
    END DO

    CALL adding_doubling(COUNT(tau(1:nLevels) > 0),Refl(:), Trans(:), Refl_tot, Trans_tot)  

    compute_toa_reflectance = Refl_tot

  END FUNCTION compute_toa_reflectance

  ! ########################################################################################
  PURE SUBROUTINE two_stream(tauint, gint, w0int, ref, tra) 
    ! Compute reflectance in a single layer using the two stream approximation 
    !   The code itself is from Lazaros Oreopoulos via Steve Platnick 
    ! INPUTS
    REAL(wp), INTENT(in)  :: tauint, gint, w0int
    ! OUTPUTS
    REAL(wp), INTENT(out) :: ref, tra
    ! LOCAL VARIABLES
    !   for delta Eddington code
    !   xmu, gamma3, and gamma4 only used for collimated beam approximation (i.e., beam=1)
    INTEGER, PARAMETER :: beam = 2
    REAL(wp),PARAMETER :: xmu = 0.866, minConservativeW0 = 0.9999999
    REAL(wp) :: tau, w0, g, f, gamma1, gamma2, gamma3, gamma4, &
         rh, a1, a2, rk, r1, r2, r3, r4, r5, t1, t2, t3, t4, t5, beta, e1, e2, ef1, ef2, den, th

    ! Compute reflectance and transmittance in a single layer using the two stream approximation 
    !   The code itself is from Lazaros Oreopoulos via Steve Platnick 
    f   = gint**2
    tau = (1._wp - w0int * f) * tauint
    w0  = (1._wp - f) * w0int / (1._wp - w0int * f)
    g   = (gint - f) / (1._wp - f)

    ! delta-Eddington (Joseph et al. 1976)
    gamma1 =  (7._wp - w0* (4._wp + 3._wp * g)) / 4._wp
    gamma2 = -(1._wp - w0* (4._wp - 3._wp * g)) / 4._wp
    gamma3 =  (2._wp - 3._wp*g*xmu) / 4._wp
    gamma4 =   1._wp - gamma3

    IF (w0int > minConservativeW0) THEN
       ! Conservative scattering
       IF (beam == 1) THEN
          rh = (gamma1*tau+(gamma3-gamma1*xmu)*(1-EXP(-tau/xmu)))
          ref = rh / (1._wp + gamma1 * tau)
          tra = 1._wp - ref       
       ELSE IF(beam == 2) THEN
          ref = gamma1*tau/(1._wp + gamma1*tau)
          tra = 1._wp - ref
       ENDIF
    ELSE
       ! Non-conservative scattering
       a1 = gamma1 * gamma4 + gamma2 * gamma3
       a2 = gamma1 * gamma3 + gamma2 * gamma4

       rk = SQRT(gamma1**2 - gamma2**2)

       r1 = (1._wp - rk * xmu) * (a2 + rk * gamma3)
       r2 = (1._wp + rk * xmu) * (a2 - rk * gamma3)
       r3 = 2._wp * rk *(gamma3 - a2 * xmu)
       r4 = (1._wp - (rk * xmu)**2) * (rk + gamma1)
       r5 = (1._wp - (rk * xmu)**2) * (rk - gamma1)

       t1 = (1._wp + rk * xmu) * (a1 + rk * gamma4)
       t2 = (1._wp - rk * xmu) * (a1 - rk * gamma4)
       t3 = 2._wp * rk * (gamma4 + a1 * xmu)
       t4 = r4
       t5 = r5

       beta = -r5 / r4         

       e1 = MIN(rk * tau, 500._wp) 
       e2 = MIN(tau / xmu, 500._wp) 

       IF (beam == 1) THEN
          den = r4 * EXP(e1) + r5 * EXP(-e1)
          ref  = w0*(r1*EXP(e1)-r2*EXP(-e1)-r3*EXP(-e2))/den
          den = t4 * EXP(e1) + t5 * EXP(-e1)
          th  = EXP(-e2)
          tra = th-th*w0*(t1*EXP(e1)-t2*EXP(-e1)-t3*EXP(e2))/den
       ELSEIF (beam == 2) THEN
          ef1 = EXP(-e1)
          ef2 = EXP(-2*e1)
          ref = (gamma2*(1._wp-ef2))/((rk+gamma1)*(1._wp-beta*ef2))
          tra = (2._wp*rk*ef1)/((rk+gamma1)*(1._wp-beta*ef2))
       ENDIF
    END IF
  END SUBROUTINE two_stream

  ! ########################################################################################
  PURE SUBROUTINE adding_doubling (npts,Refl, Tran, Refl_tot, Tran_tot)      
    ! Use adding/doubling formulas to compute total reflectance and transmittance from 
    ! layer values

    ! INPUTS
    INTEGER,INTENT(in)                  :: npts
    REAL(wp),INTENT(in),DIMENSION(npts) :: Refl,Tran
    ! OUTPUTS
    REAL(wp),INTENT(out)                :: Refl_tot, Tran_tot
    ! LOCAL VARIABLES
    INTEGER :: i
    REAL(wp), DIMENSION(npts) :: Refl_cumulative, Tran_cumulative

    Refl_cumulative(1) = Refl(1)
    Tran_cumulative(1) = Tran(1)    

    DO i=2, npts
       ! place (add) previous combined layer(s) reflectance on top of layer i, w/black surface (or ignoring surface):
       Refl_cumulative(i) = Refl_cumulative(i-1) + Refl(i)*(Tran_cumulative(i-1)**2)/(1._wp - Refl_cumulative(i-1) * Refl(i))
       Tran_cumulative(i) = (Tran_cumulative(i-1)*Tran(i)) / (1._wp - Refl_cumulative(i-1) * Refl(i))
    END DO

    Refl_tot = Refl_cumulative(SIZE(Refl))
    Tran_tot = Tran_cumulative(SIZE(Refl))

  END SUBROUTINE adding_doubling

  ! ########################################################################################
  ELEMENTAL FUNCTION two_stream_reflectance(tauint, gint, w0int)
    ! Compute reflectance in a single layer using the two stream approximation 
    !   The code itself is from Lazaros Oreopoulos via Steve Platnick 

    ! INPUTS
    REAL(wp), INTENT(in) :: tauint, gint, w0int
    ! OUTPUTS
    REAL(wp)             :: two_stream_reflectance
    ! LOCAL VARIABLES
    !   for delta Eddington code
    !   xmu, gamma3, and gamma4 only used for collimated beam approximation (i.e., beam=1)
    INTEGER, PARAMETER :: beam = 2
    REAL(wp),PARAMETER :: xmu = 0.866, minConservativeW0 = 0.9999999
    REAL(wp) :: tau, w0, g, f, gamma1, gamma2, gamma3, gamma4, &
         rh, a1, a2, rk, r1, r2, r3, r4, r5, t1, t2, t3, t4, t5, beta, e1, e2, ef1, ef2, den

    f   = gint**2
    tau = (1._wp - w0int * f) * tauint
    w0  = (1._wp - f) * w0int / (1._wp - w0int * f)
    g   = (gint - f) / (1._wp - f)

    ! delta-Eddington (Joseph et al. 1976)
    gamma1 =  (7._wp - w0* (4._wp + 3._wp * g)) / 4._wp
    gamma2 = -(1._wp - w0* (4._wp - 3._wp * g)) / 4._wp
    gamma3 =  (2._wp - 3._wp*g*xmu) / 4._wp
    gamma4 =   1._wp - gamma3

    IF (w0int > minConservativeW0) THEN
       ! Conservative scattering
       IF (beam == 1) THEN
          rh = (gamma1*tau+(gamma3-gamma1*xmu)*(1-EXP(-tau/xmu)))
          two_stream_reflectance = rh / (1._wp + gamma1 * tau)
       ELSEIF (beam == 2) THEN
          two_stream_reflectance = gamma1*tau/(1._wp + gamma1*tau)
       ENDIF

    ELSE    !

       ! Non-conservative scattering
       a1 = gamma1 * gamma4 + gamma2 * gamma3
       a2 = gamma1 * gamma3 + gamma2 * gamma4

       rk = SQRT(gamma1**2 - gamma2**2)

       r1 = (1._wp - rk * xmu) * (a2 + rk * gamma3)
       r2 = (1._wp + rk * xmu) * (a2 - rk * gamma3)
       r3 = 2._wp * rk *(gamma3 - a2 * xmu)
       r4 = (1._wp - (rk * xmu)**2) * (rk + gamma1)
       r5 = (1._wp - (rk * xmu)**2) * (rk - gamma1)

       t1 = (1._wp + rk * xmu) * (a1 + rk * gamma4)
       t2 = (1._wp - rk * xmu) * (a1 - rk * gamma4)
       t3 = 2._wp * rk * (gamma4 + a1 * xmu)
       t4 = r4
       t5 = r5

       beta = -r5 / r4         

       e1 = MIN(rk * tau, 500._wp) 
       e2 = MIN(tau / xmu, 500._wp) 

       IF (beam == 1) THEN
          den = r4 * EXP(e1) + r5 * EXP(-e1)
          two_stream_reflectance  = w0*(r1*EXP(e1)-r2*EXP(-e1)-r3*EXP(-e2))/den
       ELSEIF (beam == 2) THEN
          ef1 = EXP(-e1)
          ef2 = EXP(-2*e1)
          two_stream_reflectance = (gamma2*(1._wp-ef2))/((rk+gamma1)*(1._wp-beta*ef2))
       ENDIF

    END IF
  END FUNCTION two_stream_reflectance

  ! ########################################################################################
  PURE FUNCTION interpolate_to_min(x, y, yobs)
    ! INPUTS
    REAL(wp),INTENT(in),DIMENSION(num_trial_res) :: x, y 
    REAL(wp),INTENT(in)                          :: yobs
    ! OUTPUTS
    REAL(wp)                                     :: interpolate_to_min
    ! LOCAL VARIABLES
    REAL(wp), DIMENSION(num_trial_res)           :: diff
    INTEGER                                      :: nPoints, minDiffLoc, lowerBound, upperBound

    ! Given a set of values of y as y(x), find the value of x that minimizes abs(y - yobs)
    !   y must be monotonic in x

    nPoints = SIZE(y)
    diff(1:num_trial_res) = y(1:num_trial_res) - yobs
    minDiffLoc = MINLOC(ABS(diff), dim = 1) 

    IF(minDiffLoc == 1) THEN 
       lowerBound = minDiffLoc
       upperBound = minDiffLoc + 1
    ELSE IF(minDiffLoc == nPoints) THEN
       lowerBound = minDiffLoc - 1
       upperBound = minDiffLoc
    ELSE
       IF(diff(minDiffLoc-1) * diff(minDiffLoc) < 0) THEN
          lowerBound = minDiffLoc-1
          upperBound = minDiffLoc
       ELSE 
          lowerBound = minDiffLoc
          upperBound = minDiffLoc + 1
       END IF
    END IF

    IF(diff(lowerBound) * diff(upperBound) < 0) THEN     
       !
       ! Interpolate the root position linearly if we bracket the root
       !
       interpolate_to_min = x(upperBound) - & 
            diff(upperBound) * (x(upperBound) - x(lowerBound)) / (diff(upperBound) - diff(lowerBound))
    ELSE
       interpolate_to_min = re_fill
    END IF

  END FUNCTION interpolate_to_min

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! FUNCTION weight_by_extinction
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  FUNCTION weight_by_extinction(nLevels,tauIncrement, f, tauLimit) 
    ! INPUTS
    integer, intent(in)                    :: nLevels
    real(wp),intent(in),dimension(nLevels) :: tauIncrement, f
    real(wp),intent(in)                    :: tauLimit
    ! OUTPUTS
    real(wp)                               :: weight_by_extinction
    ! LOCAL VARIABLES
    real(wp)                               :: deltaX, totalTau, totalProduct
    integer                                :: i 

    ! Find the extinction-weighted value of f(tau), assuming constant f within each layer
    totalTau = 0._wp; totalProduct = 0._wp
    do i = 1, size(tauIncrement)
       if(totalTau + tauIncrement(i) > tauLimit) then 
          deltaX       = tauLimit - totalTau
          totalTau     = totalTau     + deltaX
          totalProduct = totalProduct + deltaX * f(i) 
       else
          totalTau     = totalTau     + tauIncrement(i) 
          totalProduct = totalProduct + tauIncrement(i) * f(i) 
       end if
       if(totalTau >= tauLimit) exit
    end do
    weight_by_extinction = totalProduct/totalTau
  END FUNCTION weight_by_extinction

END MODULE from_COSP2
