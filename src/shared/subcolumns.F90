MODULE SUBCOLUMNS

  ! Routines related to dealing with making sub-columns in the model grid
  !
  ! Salomon.Eliasson@smhi.se

  USE COSP_KINDS,           ONLY:&
       wp
  USE INTERNAL_SIMULATOR,   ONLY:&
       internal
  USE MODEL_INPUT,          ONLY:&
       model_type
  USE MOD_RNG,              ONLY:&
       init_rng,                 &
       rng_state
  USE MOD_SCOPS,            ONLY:&
       scops

  PUBLIC :: get_subcolumns,&
       correct_model_cloud_fraction

CONTAINS

  SUBROUTINE GET_SUBCOLUMNS(npoints,ncol,nlev,number,CC,CV,frac_out,subsampler)
    !
    ! Wrapper around scops
    !
    ! Salomon.Eliasson@smhi.se
    !

    IMPLICIT NONE

    ! INPUT
    INTEGER, INTENT(in) :: npoints
    REAL(wp),INTENT(in) :: number(npoints)
    INTEGER, INTENT(in) :: ncol
    INTEGER, INTENT(in) :: nlev
    REAL(wp),INTENT(in) :: CC(nlev) ! statiform cloud fraction
    REAL(wp),INTENT(in) :: CV(nlev) ! convective cloud fraction
    INTEGER, INTENT(in) :: subsampler 
                       ! scops (COSP default) = 0
                       ! DWD's subsampler     = 1
                       ! Mc ICA               = 2

    ! OUT
    REAL(wp),INTENT(out):: frac_out(npoints,ncol,nlev)

    ! SCOPS
    ! internal
    TYPE(rng_state)     :: rngs(npoints)  ! Seeds for random number generator
    INTEGER             :: seed(npoints)

    PRINT *, "--- Running scops"
    SELECT CASE(subsampler)
    CASE (0)
       ! SCOPS from COSP
       seed(:)=0
       frac_out=0._wp
       seed = INT(number)
       CALL INIT_RNG(rngs, seed)  
       
       CALL SCOPS(npoints,nlev,ncol,rngs, &
            cc,cv,3,  & ! 3 = maximum/random overlap
            frac_out,0)
       
       ! Sometimes scops produces '2's (convection) even though no
       ! convection-fraction was supplied. If you want to mask these with
       ! '1' (stratiform cloudy) supply noConvection =true. As long as count(CV>0)==0
       
       ! map 1's to 2's if don't expect convection
       IF (COUNT(CV>0)==0) THEN
          frac_out=MERGE(frac_out,1._wp,frac_out.LT.2)
       END IF
    CASE (1)
         ! Subsampler from DWD
       CALL DWD_SUBCOLS_SLIM(ncol,nlev,3,cc,frac_out)

    CASE (2)
         ! Subsampler for McICA
       PRINT*, "Mc ICA sampler not set up yet. Subsampler failed..."
    END SELECT
    
  END SUBROUTINE GET_SUBCOLUMNS

!=========================================================================

  SUBROUTINE DWD_SUBCOLS_SLIM(ncol,nlev,overlap,icc,frac_out)

    ! Note: This is a slimmed down version of similarly named subroutine in
    ! the SIMFERA satellite simulator from the DWD. This routine
    ! replicates the cloud overlap scheme in the IFS models much
    ! closer than if using the SCOPS.f90 routine from COSP

    IMPLICIT NONE

    INTEGER, INTENT(in)  :: ncol, nlev, overlap
    REAL(wp), INTENT(in) :: icc(nlev)
    REAL(wp), INTENT(out):: frac_out(ncol,nlev)

    ! local
    INTEGER :: lastcloud,lastnfb
    INTEGER :: zi,nfb, bis,von
    INTEGER, DIMENSION(:), ALLOCATABLE :: ci, ci2
    INTEGER, PARAMETER :: over_rand=2 ! random overlap
    INTEGER, PARAMETER :: over_max_rand = 3 ! maximum-random overlap
    
    lastcloud = -2
    lastnfb = -1
    frac_out = 0._wp
    
    DO zi=1, nlev !loop over model levels 

! number of filled columns
       nfb = NINT( ncol * icc(zi) )

! this is in the DWD code, but this will reduce the cloud cover by ncol/100
!            IF ( nfb > ncol ) nfb = ncol

       IF ( nfb > 0 ) THEN
          
          IF ( lastcloud .NE. (zi-1) .OR. overlap == over_rand ) THEN
             ! completely random if "air" in between cloud layers
             IF ( ALLOCATED(ci) ) DEALLOCATE( ci )
             ALLOCATE( ci(ncol) )
             CALL GET_RANDOMU( ncol, ci )
          ENDIF

          IF ( nfb .GT. lastnfb .AND. lastcloud .EQ. (zi-1) ) THEN
             ! if number of filled boxes is greater then in the 
             ! previous layer .AND. the pervious layer contained some
             ! cloud ..
             IF ( ALLOCATED(ci2) ) DEALLOCATE( ci2 )
             ALLOCATE( ci2(ncol-lastnfb) )
             CALL GET_RANDOMU( (ncol-lastnfb), ci2 )
             CALL SWAP_INDICES(ncol, lastnfb, ci2, ci)
          ENDIF

          von = 1
          bis = nfb

          frac_out(ci(von:bis),zi) = 1 ! cloudy
          lastcloud = zi
          lastnfb = nfb

       END IF !end of nfb > 0 if-loop
    END DO !end of zi-loop
         
  END SUBROUTINE DWD_SUBCOLS_SLIM

!=========================================================================

  SUBROUTINE GET_RANDOMU( ncol, indices )

! Note: This subroutine is borrowed from the SIMFERA satellite
! simulator from the DWD. This routine replicates the cloud overlap
! scheme in the IFS models much closer than if using the SCOPS.f90
! routine from COSP

    IMPLICIT NONE
    
    INTEGER,                  INTENT(IN)    :: ncol
    INTEGER, DIMENSION(ncol), INTENT(INOUT) :: indices
    
    ! local variables
    INTEGER                  :: i
    INTEGER, DIMENSION(ncol) :: xlindgen
    REAL(wp),DIMENSION(ncol) :: rndRealArr
    LOGICAL                  :: check

    indices = 0
    check = .TRUE.

    ! Create integer vector with ncol elements
    DO i = 1, ncol
       xlindgen(i) = i
    END DO
        
    DO WHILE ( check )
       ! Get an array of random numbers using RANDOM_NUMBER, which
       ! Returns a single pseudorandom number or an array of pseudorandom
       ! numbers from the uniform distribution over the range 0 < x < 1.
       rndRealArr = 0.0
       CALL RANDOM_NUMBER( rndRealArr )
       check = RANDOM_NUMBER_DUPLICATES( ncol, rndRealArr )
       !IF ( check ) PRINT*, "-- WARNING: dupliate random number -> try again!"
    END DO

    ! Sort rndRealArr in ascending order and apply it to indices
    indices = xlindgen
    CALL GET_INDICES( ncol, rndRealArr, indices )
  END SUBROUTINE GET_RANDOMU

!=========================================================================

  SUBROUTINE GET_INDICES ( n, arr, brr )

    ! Note: This subroutine is borrowed from the SIMFERA satellite
    ! simulator from the DWD. This routine replicates the cloud overlap
    ! scheme in the IFS models much closer than if using the SCOPS.f90
    ! routine from COSP
    
    IMPLICIT NONE
    
    INTEGER                              :: n   !ncols
    REAL(wp),DIMENSION(n), INTENT(INOUT) :: arr !real values
    INTEGER, DIMENSION(n), INTENT(INOUT) :: brr !integer values
    
    ! local variables
    INTEGER  :: i, j, b
    REAL(wp) :: a
    
    outer: DO j=2, n
       
       a = arr(j)
       b = brr(j)
       
       inner: DO i = j-1, 1, -1
          IF ( arr(i) .LE. a ) EXIT inner
          arr(i+1) = arr(i)
          brr(i+1) = brr(i)
       END DO inner
       
       arr(i+1) = a
       brr(i+1) = b
       
       i = 0
       
    END DO outer
    
  END SUBROUTINE GET_INDICES
  
!=========================================================================

    SUBROUTINE SWAP_INDICES(ncol, lastnfb, ci2, ci)

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: ncol
        INTEGER, INTENT(IN) :: lastnfb
        INTEGER, DIMENSION(ncol-lastnfb), INTENT(IN)    :: ci2
        INTEGER, DIMENSION(ncol),         INTENT(INOUT) :: ci

        !local 
        INTEGER                  :: i
        INTEGER, DIMENSION(ncol) :: tmpci

        tmpci = ci
        ci = 0
        ci(1:lastnfb) = tmpci(1:lastnfb)

        DO i=1, ncol-lastnfb
            ci(i+lastnfb) = tmpci(ci2(i)+lastnfb)
        END DO

    END SUBROUTINE SWAP_INDICES

!==========================================================================
    
  SUBROUTINE CORRECT_MODEL_CLOUD_FRACTION(M,ngrids,nlev,ncol)!,mask)

    !
    ! These adjustments are made:
    !
    ! 1. If neither CLWC or CIWC are non-zero, make sure CC == 0
    ! 2. Make cloud variables zero if the cloud fraction is zero
    ! !3. Make the TCC zero if none of the levels are non-zero. This can't be done here
    !

    IMPLICIT NONE

    TYPE(model_type), INTENT(inout) :: M
    INTEGER, INTENT(in)             :: ngrids,nlev,ncol

    ! debugging
    REAL(wp),dimension(ngrids,nlev)::CC_new
    INTEGER :: g,lv
    REAL(wp) :: c,f

    PRINT *, "--- Correcting cloud fraction"
    CC_new=-999._wp

    f=1._wp/ncol
    DO lv=1,nlev
       DO g=1,ngrids
          c=0
          ! loop over intervals
          DO WHILE (M%CC(g,lv) .GE. c)
             c = c+f
          END DO
          c=c-f
          ! round to nearest discrete scops value
          CC_new(g,lv)=MERGE(c+f,c,M%CC(g,lv).GT.c+f/2) 
       END DO
    END DO
        
    WHERE(M%CC.GT.0)
       M%CIWC=M%CIWC*(CC_new/M%CC)
       M%CLWC=M%CLWC*(CC_new/M%CC)
    END WHERE
    M%CC=CC_new

  END SUBROUTINE CORRECT_MODEL_CLOUD_FRACTION

!==========================================================================

  FUNCTION RANDOM_NUMBER_DUPLICATES( n, vector ) RESULT( ret )
    
    IMPLICIT NONE
    
    INTEGER  :: n, i, j, cnt
    REAL(wp) :: vector(n)
    LOGICAL  :: ret
    
    ret = .FALSE.
    
    DO i=1, n
       cnt = 0
       DO j=1, n
          IF ( vector(i) == vector(j) ) cnt = cnt + 1
          IF ( cnt > 1 ) ret = .TRUE.
       END DO
    END DO
    
  END FUNCTION RANDOM_NUMBER_DUPLICATES

!==========================================================================

END MODULE SUBCOLUMNS
