MODULE from_CC4CL
  !
  ! Borrowed from CC4CL codes and adapted to satellite simulator
  !
  ! Salomon.Eliasson@smhi.se

  USE COSP_KINDS,      ONLY: wp
  USE OPTICS_M,        ONLY: GZero_t,cloud_albedo
  IMPLICIT NONE

  PUBLIC :: Set_GZero,grid_dimension_read,read_values_3d
  PRIVATE :: locate

CONTAINS

  !-------------------------------------------------------------------------------
  ! Name: SetGZero.F90
  !
  ! Purpose:
  ! Sets up the "zero'th point" grid information used in interpolating
  ! the SAD_LUT CRP arrays.
  !
  ! Description and Algorithm details:
  ! Use the grid info in SAD_LUT and the current Tau, Re and Geom values to
  ! find the "zero'th" point for interpolation (see IntLUTTauRe comments).
  ! Populate GZero struct with zero'th point data, i.e.
  !  - index of closest (lower) point in LUT data grid, in Tau, Re, SatZen etc
  !    (i.e. the zero'th point) - ensure that 0'th point index is between 1 and
  !    npoints-1
  !  - index of next-nearest points (set to nearest point if we're at the edge
  !    of the LUT).
  !  - Denominator for calculating finite difference gradients. Will either be
  !    equal to LUT grid spacing (at edges if LUT), or twice this.
  !  - fractional grid step from 0'th point to supplied Tau, Re, Sat zen etc
  !  - 1 minus fractional step from above
  !
  ! Arguments:
  ! Name    Type   In/Out/Both Description
  ! ------------------------------------------------------------------------------
  ! Tau     real   In          Current value of optical depth
  ! Re      real   In          Current value of effective radius
  ! Ctrl    struct In          Control structure
  ! SAD_LUT struct In          Static Application Data structure  containing
  !                            the arrays of Look-Up Tables to be interpolated
  !                            and the grids on which they are stored.
  ! GZero   Struct Out         Holds "0'th point" information relating to the
  !                            grid on which the SAD_LUT CRP arrays are based.
  ! status  int    Out         Standard status code set by ECP routines
  !
  ! History:
  ! 2001/01/23, AS: original version
  ! 2007/05/31, AY: Consolidation of code for dual-view BRDF retrieval. This
  !    means added calculation of new geometry variables for second view (copied
  !    existing code, adding "_2" suffix).
  ! 2008/03/21, GT: Changes for spline interpolation: Added mT0, mT1, mR0, mR1
  !    calculations - gradients at iT0, iT1, iR0 and iR1 LUT points.
  ! 2001/09/05, CA: 'locate' routine is now used for domain searching
  ! 2011/09/05, CA: next nearest neighbour indices iTm1,iTp1,iRm1,iRp1 now
  !    evaluated
  ! 2012/12/13, CP: added in allocation of isaz0 parameter
  ! 2012/09/15, CP: remove deallocation test from this routine and instead
  !    deallocate at the end of FM.F90 this solved g95
  ! 2013/05/08, MJ: changes loop boundaries and parameter list contents memory
  !    problem
  ! 2013/12/03, MJ: makes LUTs more flexible wrt channel and properties
  ! 2014/01/16, GM: Fixed indexing operations on uninitialized (garbage) values
  !    introducing use of SPixel%spixel_y_to_ctrl_y_index and the
  !    SAD_LUT%table_use* arrays.
  ! 2014/01/23, GM: Cleaned up code.
  ! 2014/09/09, GM: Changes related to new BRDF support.
  !
  !
  ! Bugs:
  ! None known.
  !-------------------------------------------------------------------------------

  SUBROUTINE Set_GZero(Tau, Re, Solzen,LUT, GZero, isDefined)

    IMPLICIT NONE

    REAL(wp),           INTENT(in)  :: Tau
    REAL(wp),           INTENT(in)  :: Re
    REAL(wp),           INTENT(in)  :: Solzen
    TYPE(cloud_albedo), INTENT(in)  :: LUT 
    TYPE(GZero_t),      INTENT(out) :: GZero
    LOGICAL, INTENT(out)            :: isDefined
    
    ! Set the "zero'th" array indices for the interpolations, i.e. find the array
    ! indices of the nearest neighbour grid points in each dimension. Use the
    ! 'locate' function

! check to see that we will be looking within the bounds of the look
!up table
!

    isDefined = (tau .GE. 10**LUT%Tau(1) .AND. &
         tau .LE. 10**LUT%Tau(LUT%nTau) ) .AND. &
         (Re .GE. LUT%Re(1) .AND. &
         Re .LE. LUT%Re(LUT%nRe) ) .AND. &
         (Solzen .GE. LUT%Solzen(1) .AND. &
         Solzen .LE. LUT%Solzen(LUT%nSolzen) )

    IF (.NOT. isDefined) RETURN

! Salomon: This assumes that LUT%Tau is log10(Tau)
    GZero%iT0 =   MAX(MIN(locate(10**LUT%Tau(1:LUT%nTau)   ,Tau)   ,LUT%nTau   -1),1)
    GZero%iR0 =   MAX(MIN(locate(LUT%Re     (1:LUT%nRe)    ,Re )   ,LUT%nRe    -1),1)
    GZero%iSoZ0 = MAX(MIN(locate(lUT%Solzen (1:LUT%nSolzen),Solzen),LUT%nSolzen-1),1)

    ! This sets the upper bracketing index, the locate above set the lower
    ! index
    GZero%iT1      = GZero%iT0      + 1
    GZero%iR1      = GZero%iR0      + 1
    GZero%iSoZ1    = GZero%iSoZ0    + 1

! This sets the next pair of bracketing indices around the primary one
 IF (GZero%iT0 == 1) THEN
    GZero%iTm1 = GZero%iT0
    GZero%iTp1 = GZero%iT1+1
 ELSE IF (GZero%iT1 == LUT%nTau) THEN
    GZero%iTm1 = GZero%iT0-1
    GZero%iTp1 = GZero%iT1
 ELSE
    GZero%iTm1 = GZero%iT0-1
    GZero%iTp1 = GZero%iT1+1
 END IF

 IF (GZero%iR0 == 1) THEN
    GZero%iRm1 = GZero%iR0
    GZero%iRp1 = GZero%iR1+1
 ELSE IF (GZero%iR1 == LUT%nRe) THEN
    GZero%iRm1 = GZero%iR0-1
    GZero%iRp1 = GZero%iR1
 ELSE
    GZero%iRm1 = GZero%iR0-1
    GZero%iRp1 = GZero%iR1+1
 END IF

! Calculate dT, dR: these are the distances in T, R, etc from the grid
! point with indices (0, 0, ...) to (Tau, Re, ...) expressed as a
! fraction of the LUT grid steps.

 GZero%dT = (Tau - 10**LUT%Tau(GZero%iT0)) / &
      (10**LUT%Tau(GZero%iT1) - 10**LUT%Tau(GZero%iT0))

 GZero%dR = (Re - LUT%Re(GZero%iR0)) / &
      (LUT%Re(GZero%iR1) - LUT%Re(GZero%iR0))


 ! The other two angles only exist when we have solar channels
 ! (including mixed)
 GZero%dSoZ = (SolZen - LUT%Solzen(GZero%iSoZ0)) / &
      (LUT%Solzen(GZero%iSoZ1) - LUT%Solzen(GZero%iSoZ0))
 
! Calculate 1.0 minus each of the d values above - used several times by
! the interpolation routines.

 GZero%T1    = 1.0 - GZero%dT
 GZero%R1    = 1.0 - GZero%dR
 GZero%So1   = 1.0 - GZero%dSoZ
 
end subroutine Set_GZero

!-------------------------------------------------------------------------------
! Name: Locate.F90
!
! Purpose:
! Finds the location of the pair of values in the set xx that bound x
!
! Description and Algorithm details:
! From Numerical Recipes in Fortran 90 [Press, Flannery].
!
! Arguments:
! Name   Type       In/Out/Both Description
! ------------------------------------------------------------------------------
! xx     real array In          Sorted array to search
! x      real       In          Value to search for
! locate int        Out         Index of array element that bounds x
!
! History:
! 2009/04/22, CA: Original version
!
!
! Bugs:
! None known.
!---------------------------------------------------------------------

FUNCTION locate(xx,x)
  IMPLICIT NONE
  REAL(wp), DIMENSION(:), INTENT(in) :: xx
  REAL(wp),               INTENT(in) :: x
  INTEGER :: locate
  INTEGER :: n,jl,jm,ju
  LOGICAL :: ascnd
  n=SIZE(xx)
  ascnd = (xx(n) >= xx(1))
  jl=0
  ju=n+1
  DO
     IF (ju-jl <= 1) EXIT
     jm=(ju+jl)/2
     IF (ascnd .EQV. (x >= xx(jm))) THEN
        jl=jm
     ELSE
        ju=jm
     END IF
  END DO
  IF (x == xx(1)) THEN
     locate=1
  ELSE IF (x == xx(n)) THEN
     locate=n-1
  ELSE
     locate=jl
  END IF
END FUNCTION Locate

SUBROUTINE grid_dimension_read(filename, lun, n_name, d_name, v_name, &
nValues, dValues, Values)

  ! Borrowed from ReadSADLUT.F90 in the CC4CL code

IMPLICIT NONE

! Argument declarations
CHARACTER(*), INTENT(in)    :: filename
CHARACTER(*), INTENT(in)    :: n_name
CHARACTER(*), INTENT(in)    :: d_name
CHARACTER(*), INTENT(in)    :: v_name
INTEGER,      INTENT(in)    :: lun
INTEGER,      INTENT(out)   :: nValues
REAL(wp),     INTENT(out)   :: dValues
REAL(wp),     INTENT(inout) :: Values(:)

! Local variables
INTEGER :: i
INTEGER :: iostat

READ(lun, *, iostat=iostat) nValues, dValues
IF (iostat .NE. 0) THEN
WRITE(*,*) 'ERROR: grid_dimension_read(): Error reading ', TRIM(n_name), &
   ' and ', TRIM(d_name), ' from SAD LUT file: ', TRIM(filename)
STOP
END IF

READ(lun, *, iostat=iostat) (Values(i), i=1, nValues)
IF (iostat .NE. 0) THEN
WRITE(*,*) 'ERROR: grid_dimension_read(): Error reading ', TRIM(v_name), &
   ' from SAD LUT file: ', TRIM(filename)
STOP
END IF

END SUBROUTINE grid_dimension_read

!------------------------------
! Read 3d LUT values
!------------------------------
SUBROUTINE read_values_3d(filename,lun, &
n_i, n_j, n_k, values)

  ! Also borrowed from ReadSADLUT.F90 in the CC4CL code

IMPLICIT NONE

CHARACTER(*), INTENT(in) :: filename
INTEGER,      INTENT(in) :: lun
INTEGER,      INTENT(in) :: n_i
INTEGER,      INTENT(in) :: n_j
INTEGER,      INTENT(in) :: n_k
REAL(wp), INTENT(inout)  :: values(n_i,n_j,n_k)

! Local variables
INTEGER :: i, j, k
INTEGER :: iostat

READ(lun, *, iostat=iostat) (((values(i, j, k), &
i = 1, n_i), j = 1, n_j), k = 1, n_k)
IF (iostat .NE. 0) THEN
WRITE(*,*) 'ERROR: read_values_3d(): Error reading albedo &
   & from file: ', TRIM(filename)
STOP
END IF

END SUBROUTINE read_values_3d

END MODULE from_CC4CL
