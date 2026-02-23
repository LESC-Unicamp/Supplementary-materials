MODULE OverlapAlgorithms

! Uses two modules: global variables and quaternion operations
USE GlobalVariables
USE QuaternionOperations

CONTAINS

! *********************************************************************************************** !
!    This subroutine takes the relative orientations of two molecular spherocylinders i and j     !
!    and the unit vector joining their centers of mass and calculates their contact distance.     !
!           See Vega and Lago, Computers Chem. 18, 55-59 (1993), for more information.            !
! *********************************************************************************************** !
SUBROUTINE OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, ContactDistance, &
&                           ParallelSPC, OverlapSPC, FullParticle )

IMPLICIT NONE

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 ) :: SquaredDistance            ! Squared vector distance between the centers of mass of particles i and j
REAL( Kind= Real64 ) :: iVectorDistanceOrientation ! Dot product of vector distance and orientation of particle i
REAL( Kind= Real64 ) :: jVectorDistanceOrientation ! Dot product of vector distance and orientation of particle j
REAL( Kind= Real64 ) :: AngleCosine                ! Dot product of both orientations (particles i and j)
REAL( Kind= Real64 ) :: DDistanceDLambda           ! Value that minimize r² (∂r²/∂λ = 0)
REAL( Kind= Real64 ) :: DDistanceDmu               ! Value that minimize r² (∂r²/∂μ = 0)
REAL( Kind= Real64 ) :: OrthogonalityCheck         ! Orthogonality check between orientations of particles i and j
REAL( Kind= Real64 ) :: iAuxiliarVar, jAuxiliarVar ! Auxiliary variables
REAL( Kind= Real64 ) :: ContactDistance            ! Vega-Lago contact distance (variable)
REAL( Kind= Real64 ) :: ShortestDistance           ! Shortest distance
REAL( Kind= Real64 ) :: SquaredShortestDistance    ! Shortest distance (squared)

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 2 ) :: HalfLength     ! Half-length of spherocylinders i and j
REAL( Kind= Real64 ), DIMENSION( 3 ) :: iOrientation   ! Orientation of particle i
REAL( Kind= Real64 ), DIMENSION( 3 ) :: jOrientation   ! Orientation of particle j
REAL( Kind= Real64 ), DIMENSION( 3 ) :: VectorDistance ! Vector distance between the centers of mass of particles i and j

! LOGICAL VARIABLES
LOGICAL :: OverlapSPC   ! Detects overlap between two spherocylinders: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: ParallelSPC  ! Checks whether spherocylinders are parallel or not
LOGICAL :: FullParticle ! Checks whether the overlap is being computed for the whole particle or just the cylinders

! Initialization
OverlapSPC  = .FALSE.
ParallelSPC = .FALSE.

! Shortest distance between two spherocylinders
IF( FullParticle ) ShortestDistance = 0.5D0 * ( CutoffSpherocylinder + CutoffSpherocylinder )
IF( .NOT. FullParticle ) ShortestDistance = 0.5D0 * ( cDiameter + cDiameter )

SquaredShortestDistance = ShortestDistance * ShortestDistance

! Half length of spherocylinder i
HalfLength(1) = ( 0.5D0 * cLength )
! Half length of spherocylinder j
HalfLength(2) = ( 0.5D0 * cLength )

! Initial calculation
iVectorDistanceOrientation = DOT_PRODUCT( VectorDistance, iOrientation )
jVectorDistanceOrientation = DOT_PRODUCT( VectorDistance, jOrientation )
AngleCosine                = DOT_PRODUCT( iOrientation, jOrientation )
OrthogonalityCheck         = 1.D0 - ( AngleCosine * AngleCosine )

! Checking whether the spherocylinders are parallel or not
IF( DABS( OrthogonalityCheck ) < 1.D-10 ) THEN
  ParallelSPC = .TRUE.
  ! Checking whether the parallel spherocylinders are not perpendicular to the intermolecular axis (avoid the indeterminate form 0/0)
  IF( DABS( iVectorDistanceOrientation ) >= 1.D-10 .AND. DABS( jVectorDistanceOrientation ) >= 1.D-10 ) THEN
    ! Take the extreme side of particle i
    DDistanceDLambda = DSIGN( HalfLength(1), iVectorDistanceOrientation )
    ! Closest point between particle i and particle j
    DDistanceDmu = ( DDistanceDLambda * AngleCosine ) - jVectorDistanceOrientation
    ! Take the extreme side of particle j if λ' > L/2
    IF( DABS( DDistanceDmu ) > HalfLength(2) ) THEN
      DDistanceDmu = DSIGN( HalfLength(2), DDistanceDmu )
    END IF
    ! Shortest distance (squared)
    ContactDistance = SquaredDistance + ( DDistanceDLambda * DDistanceDLambda ) + ( DDistanceDmu * DDistanceDmu ) - ( 2.D0 * &
    &                 DDistanceDLambda * DDistanceDmu * AngleCosine ) + ( 2.D0 * DDistanceDmu * jVectorDistanceOrientation ) - &
    &                 ( 2.D0 * DDistanceDLambda * iVectorDistanceOrientation )
    ! Overlap criterion
    IF( ContactDistance <= SquaredShortestDistance ) THEN
      OverlapSPC = .TRUE.
      RETURN
    END IF
    ! Take the extreme side of particle j
    DDistanceDLambda = DSIGN( HalfLength(2), jVectorDistanceOrientation )
    ! Closest point between particle j and particle i
    DDistanceDmu = ( DDistanceDLambda * AngleCosine ) - iVectorDistanceOrientation
    ! Take the extreme side of particle i if λ' > L/2
    IF( DABS( DDistanceDmu ) > HalfLength(1) ) THEN
      DDistanceDmu = DSIGN( HalfLength(1), DDistanceDmu )
    END IF
    ! Shortest distance (squared)
    ContactDistance = SquaredDistance + ( DDistanceDLambda * DDistanceDLambda ) + ( DDistanceDmu * DDistanceDmu ) - ( 2.D0 * &
    &                 DDistanceDLambda * DDistanceDmu * AngleCosine ) + ( 2.D0 * DDistanceDmu * jVectorDistanceOrientation ) - &
    &                 ( 2.D0 * DDistanceDLambda * iVectorDistanceOrientation )
    IF( ContactDistance <= SquaredShortestDistance ) THEN
      OverlapSPC = .TRUE.
      RETURN
    END IF
    ! No overlaps
    RETURN
  ! Parallel spherocylinders almost orthogonal to the intermolecular axis (avoid the indeterminate form 0/0)
  ELSE
    DDistanceDLambda = 0.D0
    DDistanceDmu     = 0.D0
    ! Shortest distance (squared)
    ContactDistance = SquaredDistance
    ! Overlap criterion
    IF( ContactDistance <= SquaredShortestDistance ) THEN
      OverlapSPC = .TRUE.
    END IF
    ! Return immediately
    RETURN
  END IF
END IF

! *********************************************************************************************** !
! STEP 1: Evaluation of (λ’, μ’) according to Equations (3) and (4)                               !
! *********************************************************************************************** !
!  See Vega and Lago, Computers Chem. (1994) for more information                                 !
! *********************************************************************************************** !
DDistanceDLambda = ( iVectorDistanceOrientation - ( AngleCosine * jVectorDistanceOrientation ) ) / OrthogonalityCheck
DDistanceDmu     = ( -jVectorDistanceOrientation + ( AngleCosine * iVectorDistanceOrientation ) ) / OrthogonalityCheck

! *********************************************************************************************** !
! STEP 2: Check whether the point (λ’, μ’) is in the rectangle (λ, μ), with λ = μ = [-L/2, L/2]   !
! *********************************************************************************************** !
!  See Vega and Lago, Computers Chem. (1994) for more information                                 !
! *********************************************************************************************** !
! Point (λ’, μ’) in the rectangle (λ, μ)
IF( ( DABS( DDistanceDLambda ) <= HalfLength(1) ) .AND. ( DABS( DDistanceDmu ) <= HalfLength(2) ) ) THEN
  ContactDistance = SquaredDistance + ( DDistanceDLambda * DDistanceDLambda ) + ( DDistanceDmu * DDistanceDmu ) - ( 2.D0 * &
  &                 DDistanceDLambda * DDistanceDmu * AngleCosine ) + ( 2.D0 * DDistanceDmu * jVectorDistanceOrientation ) - &
  &                 ( 2.D0 * DDistanceDLambda * iVectorDistanceOrientation )
  ! Overlap criterion
  IF( ContactDistance <= SquaredShortestDistance ) THEN
    OverlapSPC = .TRUE.
  END IF
  ! Return immediately
  RETURN
! Point (λ’, μ’) not in the rectangle (λ, μ)
ELSE
  iAuxiliarVar = DABS( DDistanceDLambda ) - HalfLength(1) ! Used to determine the region (extreme)
  jAuxiliarVar = DABS( DDistanceDmu ) - HalfLength(2)     ! Used to determine the region (extreme)
END IF

! *********************************************************************************************** !
! STEP 3-7: Shortest distance of the considered extreme (regions 1, 2, 3, or 4) to the line where !
! the other spherocylinder is contained                                                           !
! *********************************************************************************************** !
!  See Vega and Lago, Computers Chem. (1994) for more information                                 !
! *********************************************************************************************** !
! Region 1 or 3
IF( iAuxiliarVar > jAuxiliarVar ) THEN
  DDistanceDLambda = DSIGN( HalfLength(1), DDistanceDLambda )
  DDistanceDmu = ( DDistanceDLambda * AngleCosine ) - jVectorDistanceOrientation
  IF( DABS( DDistanceDmu ) > HalfLength(2) ) THEN
    DDistanceDmu = DSIGN( HalfLength(2), DDistanceDmu )
  END IF
! Region 2 or 4
ELSE
  DDistanceDmu = DSIGN( HalfLength(2), DDistanceDmu )
  DDistanceDLambda = ( DDistanceDmu * AngleCosine ) + iVectorDistanceOrientation
  IF( DABS( DDistanceDLambda ) > HalfLength(1) ) THEN
    DDistanceDLambda = DSIGN( HalfLength(1), DDistanceDLambda )
  END IF
END IF

! *********************************************************************************************** !
! STEP 8: Evaluate the shortest distance (squared)                                                !
! *********************************************************************************************** !
!  See Vega and Lago, Computers Chem. (1994) for more information                                 !
! *********************************************************************************************** !
ContactDistance = SquaredDistance + ( DDistanceDLambda * DDistanceDLambda ) + ( DDistanceDmu * DDistanceDmu ) - ( 2.D0 * &
&                 DDistanceDLambda * DDistanceDmu * AngleCosine ) + ( 2.D0 * DDistanceDmu * jVectorDistanceOrientation ) - &
&                 ( 2.D0 * DDistanceDLambda * iVectorDistanceOrientation )
! Overlap criterion
IF( ContactDistance <= SquaredShortestDistance ) THEN
  OverlapSPC = .TRUE.
END IF

! Return immediately
RETURN

END SUBROUTINE OverlapCheckSPC

! *********************************************************************************************** !
! This subroutine takes the relative quaternions/orientations of two molecular cylinders i and j  !
!      as well as the position of their centers of mass and the unit vector joining them and      !
!                        calculates whether the cylinders overlap or not.                         !
!             See Lopes et al., Chem. Phys. 154, 104902 (2021), for more information.             !
!            See Orellana et al., Eur. Phys. J. E 41, 51 (2018), for more information.            !
!       See Ibarra-Avalos et al., Mol. Simul. 33, 6, 505–515 (2007), for more information.        !
! *********************************************************************************************** !
SUBROUTINE OverlapCheckCYL( iQuaternion, jQuaternion, iOrientation, jOrientation, VectorDistance, iPosition, jPosition, &
&                           ParallelSPC, OverlapCYL )

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 ) :: iDisk, jDisk, dDisk ! Counters (disk)

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 ) :: cSquaredLength                  ! Cylindrical length (squared)
REAL( Kind= Real64 ) :: cSquaredDiameter                ! Cylindrical diameter (squared)
REAL( Kind= Real64 ) :: iVectorDistanceOrientation      ! Dot product of vector distance and orientation of particle i
REAL( Kind= Real64 ) :: jVectorDistanceOrientation      ! Dot product of vector distance and orientation of particle j
REAL( Kind= Real64 ) :: SquaredVectorDistanceParallel   ! Squared vector distance between particles i and j (parallel)
REAL( Kind= Real64 ) :: SquaredVectorDistanceOrthogonal ! Squared vector distance between particles i and j (orthogonal)

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 2 )    :: dSquaredDiameter             ! Diameters of the cylindrical disks for the disk-rim configuration
REAL( Kind= Real64 ), DIMENSION( 2 )    :: HalfLength                   ! Half-lengths of cylinders i and j
REAL( Kind= Real64 ), DIMENSION( 2 )    :: HalfDiameter                 ! Half-diameters of cylinders i and j
REAL( Kind= Real64 ), DIMENSION( 3 )    :: iOrientation, jOrientation   ! Orientation of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )    :: iPosition, jPosition         ! Position of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )    :: VectorDistance               ! Vector distance between particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )    :: VectorDistanceParallel       ! Vector distance between particles i and j (parallel)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: VectorDistanceOrthogonal     ! Vector distance between particles i and j (orthogonal)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: pRimPosition                 ! Position of the cylindrical rim
REAL( Kind= Real64 ), DIMENSION( 3 )    :: pRimOrientation              ! Orientation of the cylindrical rim
REAL( Kind= Real64 ), DIMENSION( 0: 3 ) :: iQuaternion, jQuaternion     ! Quaternion of particles i and j
REAL( Kind= Real64 ), DIMENSION( 0: 3 ) :: pDiskQuaternion              ! Quaternion of the cylindrical disks
REAL( Kind= Real64 ), DIMENSION( 3, 2 ) :: iDiskPosition, jDiskPosition ! Position of cylindrical disks of particles i and j

! LOGICAL VARIABLES
LOGICAL :: OverlapCYL      ! Detects overlap between two cylindrical particles : TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: OverlapRimRim   ! Detects overlap between two particles (rim-rim configuration) : TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: OverlapDiskDisk ! Detects overlap between two particles (disk-disk configuration) : TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: OverlapDiskRim  ! Detects overlap between two particles (disk-rim configuration) : TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: ParallelSPC     ! Checks the relative orientation of two spherocylinders : TRUE = parallel orientation; FALSE = non-parallel orientation

! Half length of cylinder i
HalfLength(1) = 0.5D0 * cLength
! Half length of cylinder j
HalfLength(2) = 0.5D0 * cLength
! Half diameter of cylinder i
HalfDiameter(1) = 0.5D0 * cDiameter
! Half diameter of cylinder j
HalfDiameter(2) = 0.5D0 * cDiameter

! Cylindrical length (squared)
cSquaredLength = 0.5D0 * ( cLength + cLength )
cSquaredLength = cSquaredLength * cSquaredLength
! Cylindrical diameter (squared)
cSquaredDiameter = 0.5D0 * ( cDiameter + cDiameter )
cSquaredDiameter = cSquaredDiameter * cSquaredDiameter

! Diameter of the disks (squared)
dSquaredDiameter(1) = cDiameter * cDiameter
dSquaredDiameter(2) = cDiameter * cDiameter

! Initialization
OverlapCYL = .FALSE.

! Initial calculation
iVectorDistanceOrientation = DOT_PRODUCT( VectorDistance, iOrientation )
jVectorDistanceOrientation = DOT_PRODUCT( VectorDistance, jOrientation )

! *********************************************************************************************** !
! CASE 1: PARALLEL CYLINDERS                                                                      !
! *********************************************************************************************** !
!  See Lopes et al., Chem. Phys. 154, 104902 (2021), for more information                         !
! *********************************************************************************************** !
IF( ParallelSPC ) THEN
  ! Vector distance between cylinders i and j (parallel)
  VectorDistanceParallel(1) = iOrientation(1) * iVectorDistanceOrientation
  VectorDistanceParallel(2) = iOrientation(2) * iVectorDistanceOrientation
  VectorDistanceParallel(3) = iOrientation(3) * iVectorDistanceOrientation
  ! Vector distance between cylinders i and j (orthogonal)
  VectorDistanceOrthogonal(1) = VectorDistance(1) - VectorDistanceParallel(1)
  VectorDistanceOrthogonal(2) = VectorDistance(2) - VectorDistanceParallel(2)
  VectorDistanceOrthogonal(3) = VectorDistance(3) - VectorDistanceParallel(3)
  ! Magnitude of vector distance between cylinders i and j (parallel)
  SquaredVectorDistanceParallel = DOT_PRODUCT( VectorDistanceParallel, VectorDistanceParallel )
  ! Magnitude of vector distance between cylinders i and j (orthogonal)
  SquaredVectorDistanceOrthogonal = DOT_PRODUCT( VectorDistanceOrthogonal, VectorDistanceOrthogonal )
  ! Overlap criterion
  IF( ( SquaredVectorDistanceParallel <= cSquaredLength ) .AND. ( SquaredVectorDistanceOrthogonal <= cSquaredDiameter ) ) THEN
    OverlapCYL = .TRUE.
    RETURN ! Return immediately
  ELSE
    OverlapCYL = .FALSE.
    RETURN ! Return immediately
  END IF
END IF

! *********************************************************************************************** !
! CASE 2: RIM-RIM CONFIGURATION                                                                   !
! *********************************************************************************************** !
!  See Vega and Lago, Computers Chem. 18, 55-59 (1993), for more information                      !
! *********************************************************************************************** !
CALL RimRimConfiguration( iOrientation, jOrientation, iPosition, jPosition, VectorDistance, HalfLength, cSquaredDiameter, &
&                         OverlapRimRim )
IF( OverlapRimRim ) THEN
  OverlapCYL = .TRUE.
  RETURN ! Return immediately
END IF

! *********************************************************************************************** !
! CASE 3: DISK-DISK CONFIGURATION                                                                 !
! *********************************************************************************************** !
!  See Ibarra-Avalos et al., Mol. Simul. 33, 6, 505–515 (2007), for more information              !
! *********************************************************************************************** !
! Disk 1 of cylinder i
iDiskPosition(1,1) = iPosition(1) + ( iOrientation(1) * HalfLength(1) )
iDiskPosition(2,1) = iPosition(2) + ( iOrientation(2) * HalfLength(1) )
iDiskPosition(3,1) = iPosition(3) + ( iOrientation(3) * HalfLength(1) )
! Disk 2 of cylinder i
iDiskPosition(1,2) = iPosition(1) - ( iOrientation(1) * HalfLength(1) )
iDiskPosition(2,2) = iPosition(2) - ( iOrientation(2) * HalfLength(1) )
iDiskPosition(3,2) = iPosition(3) - ( iOrientation(3) * HalfLength(1) )
! Disk 1 of cylinder j
jDiskPosition(1,1) = jPosition(1) + ( jOrientation(1) * HalfLength(2) )
jDiskPosition(2,1) = jPosition(2) + ( jOrientation(2) * HalfLength(2) )
jDiskPosition(3,1) = jPosition(3) + ( jOrientation(3) * HalfLength(2) )
! Disk 2 of cylinder j
jDiskPosition(1,2) = jPosition(1) - ( jOrientation(1) * HalfLength(2) )
jDiskPosition(2,2) = jPosition(2) - ( jOrientation(2) * HalfLength(2) )
jDiskPosition(3,2) = jPosition(3) - ( jOrientation(3) * HalfLength(2) )
! Search for overlaps between all pairs of disks from both cylinders
DO iDisk = 1, 2
  DO jDisk = 1, 2
    CALL DiskDiskConfiguration( iOrientation, jOrientation, iDiskPosition(:,iDisk), jDiskPosition(:,jDisk), HalfDiameter, &
    &                           OverlapDiskDisk )
    IF( OverlapDiskDisk ) THEN
      OverlapCYL = .TRUE.
      RETURN ! Return immediately
    END IF
  END DO
END DO

! *********************************************************************************************** !
! CASE 4: DISK-RIM CONFIGURATION                                                                  !
! *********************************************************************************************** !
!  *MODIFIED VERSION* of the algorithm developed by Lopes et al., Chem. Phys. 154, 104902 (2021)  !
! *********************************************************************************************** !
! Rim of cylinder i
pRimPosition(1) = iPosition(1)
pRimPosition(2) = iPosition(2)
pRimPosition(3) = iPosition(3)
! Orientation of cylinder i
pRimOrientation(1) = iOrientation(1)
pRimOrientation(2) = iOrientation(2)
pRimOrientation(3) = iOrientation(3)
! Quaternion of the disks of cylinder j
pDiskQuaternion(0) = jQuaternion(0)
pDiskQuaternion(1) = jQuaternion(1)
pDiskQuaternion(2) = jQuaternion(2)
pDiskQuaternion(3) = jQuaternion(3)
! Search for overlaps between both disks of cylinder j and the rim of cylinder i
DO dDisk = 1, 2
  CALL DiskRimConfiguration( jDiskPosition(:,dDisk), pDiskQuaternion, HalfDiameter(2), pRimPosition, pRimOrientation, &
  &                          HalfDiameter(1), HalfLength(1), OverlapDiskRim )
  IF( OverlapDiskRim ) THEN
    OverlapCYL = .TRUE.
    RETURN ! Return immediately
  END IF
END DO
! Rim of cylinder j
pRimPosition(1) = jPosition(1)
pRimPosition(2) = jPosition(2)
pRimPosition(3) = jPosition(3)
! Orientation of cylinder j
pRimOrientation(1) = jOrientation(1)
pRimOrientation(2) = jOrientation(2)
pRimOrientation(3) = jOrientation(3)
! Quaternion of the disks of cylinder i
pDiskQuaternion(0) = iQuaternion(0)
pDiskQuaternion(1) = iQuaternion(1)
pDiskQuaternion(2) = iQuaternion(2)
pDiskQuaternion(3) = iQuaternion(3)
! Search for overlaps between both disks of cylinder i and the rim of cylinder j
DO dDisk = 1, 2
  CALL DiskRimConfiguration( iDiskPosition(:,dDisk), pDiskQuaternion, HalfDiameter(1), pRimPosition, pRimOrientation, &
  &                          HalfDiameter(2), HalfLength(2), OverlapDiskRim )
  IF( OverlapDiskRim ) THEN
    OverlapCYL = .TRUE.
    RETURN ! Return immediately
  END IF
END DO

! No overlaps
RETURN

END SUBROUTINE OverlapCheckCYL

! *********************************************************************************************** !
!                                   Rim-Rim Overlap Algorithm                                     !
! *********************************************************************************************** !
!            See Vega and Lago, Computers Chem. 18, 55-59 (1993), for more information            !
! *********************************************************************************************** !
SUBROUTINE RimRimConfiguration( iOrientation, jOrientation, iPosition, jPosition, VectorDistance, HalfLength, cSquaredDiameter, &
&                               OverlapRimRim )

IMPLICIT NONE

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 ) :: cSquaredDiameter             ! Cylinder geometrical diameter (squared)
REAL( Kind= Real64 ) :: SquaredClosestVectorDistance ! Cylinder geometrical diameter (squared)
REAL( Kind= Real64 ) :: DDistanceDLambda             ! Value that minimize r² (∂r²/∂λ = 0)
REAL( Kind= Real64 ) :: DDistanceDmu                 ! Value that minimize r² (∂r²/∂μ = 0)

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 2 ) :: HalfLength                 ! Half length of cylinders i and j
REAL( Kind= Real64 ), DIMENSION( 3 ) :: iClosestReferencePosition  ! Auxiliary vector (from spherocylinder overlap algorithm)
REAL( Kind= Real64 ), DIMENSION( 3 ) :: jClosestReferencePosition  ! Auxiliary vector (from spherocylinder overlap algorithm)
REAL( Kind= Real64 ), DIMENSION( 3 ) :: iClosestPosition           ! Position of points of closest approach on the axes of the particle i 
REAL( Kind= Real64 ), DIMENSION( 3 ) :: jClosestPosition           ! Position of points of closest approach on the axes of the particle j
REAL( Kind= Real64 ), DIMENSION( 3 ) :: iOrientation, jOrientation ! Orientation of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 ) :: iPosition, jPosition       ! Position of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 ) :: VectorDistance             ! Vector distance between particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 ) :: ClosestVectorDistance      ! Closest vector distance between particles i and j

! LOGICAL VARIABLES
LOGICAL :: OverlapRimRim ! Detects overlap between two particles (rim-rim configuration) : TRUE = overlap detected; FALSE = overlap not detected

! Initialize logical variable
OverlapRimRim = .FALSE.

! Values that minimize r² (from Vega-Lago algorithm)
DDistanceDLambda = 1.D0 / (1.D0 - DOT_PRODUCT( iOrientation, jOrientation ) * DOT_PRODUCT( iOrientation, jOrientation )  ) * &
&                  ( DOT_PRODUCT( VectorDistance, iOrientation ) - DOT_PRODUCT( iOrientation, jOrientation ) * &
&                  DOT_PRODUCT( VectorDistance, jOrientation ) )
DDistanceDmu     = 1.D0 / (1.D0 - DOT_PRODUCT( iOrientation, jOrientation ) * DOT_PRODUCT( iOrientation, jOrientation )  ) * &
&                  ( - DOT_PRODUCT( VectorDistance, jOrientation ) + DOT_PRODUCT( iOrientation, jOrientation ) * &
&                  DOT_PRODUCT( VectorDistance, iOrientation ) )

! Position of λ and μ on the orientational axes of cylinders i and j
iClosestReferencePosition = DDistanceDLambda * iOrientation
jClosestReferencePosition = DDistanceDmu * jOrientation

! Position of the point of closest approach on the axis of cylinder i
iClosestPosition = iPosition + iClosestReferencePosition

! Position of the point of closest approach on the axis of cylinder j
jClosestPosition = jPosition + jClosestReferencePosition

! Vector distance between the points of closest approach on the axes of cylinders i and j
ClosestVectorDistance = jClosestPosition - iClosestPosition

! Vector distance between the points of closest approach on the axes of cylinders i and j (squared)
SquaredClosestVectorDistance = DOT_PRODUCT( ClosestVectorDistance, ClosestVectorDistance )

! Overlap criterion
IF( ( SquaredClosestVectorDistance <= cSquaredDiameter ) .AND. ( DABS( DDistanceDLambda ) <= HalfLength(1) ) .AND. &
&   ( DABS( DDistanceDmu ) <= HalfLength(2) ) ) THEN
  OverlapRimRim = .TRUE.
  RETURN ! Return immediately
END IF

! *********************************************************************************************** !
! OBSERVATION                                                                                     !
! *********************************************************************************************** !
! The first condition is precisely the Vega-Lago criterion. Since the overlap between cylinders   !
! is only taken into account when the circumscribing spherocylinders overlap each other, this     !
! condition must always hold true. The second and third conditions are additional and necessary   !
! conditions for cylinders after removing the hemispherical caps.                                 !
! *********************************************************************************************** !

! No overlaps
RETURN 

END SUBROUTINE RimRimConfiguration

! *********************************************************************************************** !
!                                  Disk-Disk Overlap Algorithm                                    !
! *********************************************************************************************** !
!        See Ibarra-Avalos et al., Mol. Simul. 33, 6, 505–515 (2007), for more information        !
! *********************************************************************************************** !
SUBROUTINE DiskDiskConfiguration( iOrientation, jOrientation, iDiskPosition, jDiskPosition, HalfDiameter, OverlapDiskDisk )

IMPLICIT NONE

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 ) :: dSquaredDistance            ! Vector distance between disks i and j (squared)
REAL( Kind= Real64 ) :: SquaredCutoffDistance       ! Magnitude of the vector distance between the circumscribing spheres (squared)
REAL( Kind= Real64 ) :: SegmentSum                  ! Sum of the segments of disks of cylinders i and j
REAL( Kind= Real64 ) :: iSegment                    ! Magnitude of the segment between a point on the circumference of disk of cylinder i and the point over the intersection line
REAL( Kind= Real64 ) :: jSegment                    ! Magnitude of the segment between a point on the circumference of disk of cylinder j and the point over the intersection line
REAL( Kind= Real64 ) :: NearestDistance             ! Magnitude of the vector distance between the points over the intersection line defined by the plane of the disks
REAL( Kind= Real64 ) :: iSquaredNearestPositionDisk ! Magnitude of the vector distance between a point over the intersection line and the center of mass of disk of cylinder i (squared)
REAL( Kind= Real64 ) :: jSquaredNearestPositionDisk ! Magnitude of the vector distance between a point over the intersection line and the center of mass of disk of cylinder j (squared)

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 2 ) :: rSquaredRadius               ! Radii of cylinders i and j (squared)
REAL( Kind= Real64 ), DIMENSION( 2 ) :: HalfDiameter                 ! Half diameter of cylinders i and j
REAL( Kind= Real64 ), DIMENSION( 3 ) :: dVectorDistance              ! Vector distance between the centers of mass of disks of cylinders i and j
REAL( Kind= Real64 ), DIMENSION( 3 ) :: iOrientation, jOrientation   ! Orientation of cylinders i and j
REAL( Kind= Real64 ), DIMENSION( 3 ) :: IntersectionLineOrientation  ! Orientation of the intersection line
REAL( Kind= Real64 ), DIMENSION( 3 ) :: iCrossProductOrientation     ! Cross product between the orientation of a disk of cylinder j and orientation of the intersection line
REAL( Kind= Real64 ), DIMENSION( 3 ) :: jCrossProductOrientation     ! Cross product between the orientation of a disk of cylinder i and orientation of the intersection line
REAL( Kind= Real64 ), DIMENSION( 3 ) :: iNearestPosition             ! Position of the points over the intersection line that are nearest of the centers of mass of disks of cylinder i
REAL( Kind= Real64 ), DIMENSION( 3 ) :: jNearestPosition             ! Position of the points over the intersection line that are nearest of the centers of mass of disks of cylinder j
REAL( Kind= Real64 ), DIMENSION( 3 ) :: NearestVectorDistance        ! Vector distance between the points over the intersection line defined by the plane of the disks
REAL( Kind= Real64 ), DIMENSION( 3 ) :: iNearestPositionDisk         ! Vector distance between the points over the intersection line and the centers of mass of disks of cylinder i
REAL( Kind= Real64 ), DIMENSION( 3 ) :: jNearestPositionDisk         ! Vector distance between the points over the intersection line and the centers of mass of disks of cylinder j
REAL( Kind= Real64 ), DIMENSION( 3 ) :: iDiskPosition, jDiskPosition ! Position of both disks of cylinders i and j

! LOGICAL VARIABLES
LOGICAL :: OverlapDiskDisk ! Detects overlap between two particles (disk-disk configuration) : TRUE = overlap detected; FALSE = overlap not detected

! Initialize logical variable
OverlapDiskDisk = .FALSE.

! Orientation of the intersection line
IntersectionLineOrientation(1) = iOrientation(2) * jOrientation(3) - iOrientation(3) * jOrientation(2)
IntersectionLineOrientation(2) = iOrientation(3) * jOrientation(1) - iOrientation(1) * jOrientation(3)
IntersectionLineOrientation(3) = iOrientation(1) * jOrientation(2) - iOrientation(2) * jOrientation(1)

! Cross product between the orientation of a disk of cylinder j and orientation of the intersection line
iCrossProductOrientation(1) = iOrientation(2) * IntersectionLineOrientation(3) - iOrientation(3) * IntersectionLineOrientation(2)
iCrossProductOrientation(2) = iOrientation(3) * IntersectionLineOrientation(1) - iOrientation(1) * IntersectionLineOrientation(3)
iCrossProductOrientation(3) = iOrientation(1) * IntersectionLineOrientation(2) - iOrientation(2) * IntersectionLineOrientation(1)

! Cross product between the orientation of a disk of cylinder i and orientation of the intersection line
jCrossProductOrientation(1) = jOrientation(2) * IntersectionLineOrientation(3) - jOrientation(3) * IntersectionLineOrientation(2)
jCrossProductOrientation(2) = jOrientation(3) * IntersectionLineOrientation(1) - jOrientation(1) * IntersectionLineOrientation(3)
jCrossProductOrientation(3) = jOrientation(1) * IntersectionLineOrientation(2) - jOrientation(2) * IntersectionLineOrientation(1)

! Magnitude of the vector distance between the circumscribing spheres (squared)
SquaredCutoffDistance = ( HalfDiameter(1) + HalfDiameter(2) ) * ( HalfDiameter(1) + HalfDiameter(2) )

! Vector distance between both disks of cylinders i and j
dVectorDistance(1) = jDiskPosition(1) - iDiskPosition(1)
dVectorDistance(2) = jDiskPosition(2) - iDiskPosition(2)
dVectorDistance(3) = jDiskPosition(3) - iDiskPosition(3)

! Magnitude of the vector distance between both disks
dSquaredDistance = DOT_PRODUCT( dVectorDistance, dVectorDistance )

! Condition of parallel disks
IF( DABS( DOT_PRODUCT( dVectorDistance, iOrientation ) ) <= 1.D-10 .AND. dSquaredDistance <= SquaredCutoffDistance ) THEN
  OverlapDiskDisk = .TRUE.
  RETURN ! Return immediately
! Non-overlapping condition (spheres circumscribing the disks)
ELSE IF( dSquaredDistance > SquaredCutoffDistance ) THEN
  RETURN ! Return immediately
END IF

! Position of a point over the intersection line that is nearest of the center of mass of disk of cylinder i
iNearestPosition = ( DOT_PRODUCT( iDiskPosition, IntersectionLineOrientation ) * IntersectionLineOrientation + &
&                  DOT_PRODUCT( iDiskPosition, iOrientation ) * jCrossProductOrientation - &
&                  DOT_PRODUCT( jDiskPosition, jOrientation ) * iCrossProductOrientation ) / &
&                  DOT_PRODUCT( IntersectionLineOrientation, IntersectionLineOrientation )

! Position of a point over the intersection line that is nearest of the center of mass of disk of cylinder j
jNearestPosition = ( DOT_PRODUCT( jDiskPosition, IntersectionLineOrientation ) * IntersectionLineOrientation + &
&                  DOT_PRODUCT( iDiskPosition, iOrientation ) * jCrossProductOrientation - &
&                  DOT_PRODUCT( jDiskPosition, jOrientation ) * iCrossProductOrientation ) / &
&                  DOT_PRODUCT( IntersectionLineOrientation, IntersectionLineOrientation )

! Vector distance between a point over the intersection line and the center of mass of disk of cylinder i
iNearestPositionDisk = iNearestPosition - iDiskPosition
! Magnitude of the vector distance between a point over the intersection line and the center of mass of disk of cylinder i (squared)
iSquaredNearestPositionDisk = DOT_PRODUCT( iNearestPositionDisk, iNearestPositionDisk )

! Vector distance between a point over the intersection line and the center of mass of disk of cylinder i
jNearestPositionDisk = jNearestPosition - jDiskPosition
! Magnitude of the vector distance between a point over the intersection line and the center of mass of disk of cylinder j (squared)
jSquaredNearestPositionDisk = DOT_PRODUCT( jNearestPositionDisk, jNearestPositionDisk )

! Radii (squared)
rSquaredRadius(1) = ( HalfDiameter(1) * HalfDiameter(1) )
rSquaredRadius(2) = ( HalfDiameter(2) * HalfDiameter(2) )

! Preliminary overlap condition
IF( ( iSquaredNearestPositionDisk < rSquaredRadius(1) ) .AND. ( jSquaredNearestPositionDisk < rSquaredRadius(2) ) ) THEN
  ! Length of the segments between the points on the circumference of disks of cylinders i and j and the points over the intersection line
  iSegment = DSQRT( rSquaredRadius(1) - iSquaredNearestPositionDisk )
  jSegment = DSQRT( rSquaredRadius(2) - jSquaredNearestPositionDisk )
  ! Sum of the segments
  SegmentSum = iSegment + jSegment
  ! Vector distance between the points over the intersection line
  NearestVectorDistance = iNearestPosition - jNearestPosition
  ! Magnitude of the vector distance between the points over the intersection line
  NearestDistance = DSQRT( DOT_PRODUCT( NearestVectorDistance, NearestVectorDistance ) )
  ! Overlap criterion
  IF( NearestDistance <= SegmentSum ) THEN
    OverlapDiskDisk = .TRUE.
    RETURN ! Return immediately
  END IF
END IF

! No overlaps
RETURN

END SUBROUTINE DiskDiskConfiguration

! *********************************************************************************************** !
!                                   Disk-Rim Overlap Algorithm                                    !
! *********************************************************************************************** !
!  *MODIFIED VERSION* of the algorithm developed by Lopes et al., Chem. Phys. 154, 104902 (2021)  !
! *********************************************************************************************** !
SUBROUTINE DiskRimConfiguration( pDiskPosition, pDiskQuaternion, dHalfDiameter, pRimPosition, pRimOrientation, rHalfDiameter, &
&                                rHalfLength, OverlapDiskRim )

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 ) :: cPolynomial ! Counter (polynomial coefficient)
INTEGER( Kind= Int64 ) :: cBrent      ! Counter (numerical method)
INTEGER( Kind= Int64 ) :: rInterval   ! Loop over λ intervals
INTEGER( Kind= Int64 ) :: nIntervals  ! Maximum number of intervals

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 ) :: drVectorDistanceOrientationX       ! Dot product of the vector distance between a disk (j or i) and a rim (i or j) and the orientation of the cylindrical disk along x-direction (j or i)
REAL( Kind= Real64 ) :: drVectorDistanceOrientationY       ! Dot product of the vector distance between a disk (j or i) and a rim (i or j) and the orientation of the cylindrical disk along y-direction (j or i)
REAL( Kind= Real64 ) :: drVectorDistanceOrientationZ       ! Dot product of the vector distance between a disk (j or i) and a rim (i or j) and the orientation of the cylindrical rim (i or j)
REAL( Kind= Real64 ) :: drOrientationXZ                    ! Dot product of the orientation of the cylindrical disk (j or i) along x-direction and the orientation of the cylindrical rim (i or j)
REAL( Kind= Real64 ) :: drOrientationYZ                    ! Dot product of the orientation of the cylindrical disk (j or i) along y-direction and the orientation of the cylindrical rim (i or j)
REAL( Kind= Real64 ) :: rSquaredClosestVectorDistanceDisk  ! Magnitude of the vector distance between the closest point on the cylindrical rim (i or j) and the center of the cylindrical disk (j or i)
REAL( Kind= Real64 ) :: Tolerance                          ! Numerical tolerance of the Brent's method
REAL( Kind= Real64 ) :: CosinePhi, SinePhi                 ! Cosine and sine of the angle φ that defines a point in the circumference of the cylindrical disk (j or i)
REAL( Kind= Real64 ) :: OppositeCathetus                   ! Opposite cathetus
REAL( Kind= Real64 ) :: AdjacentCathetus                   ! Adjacent cathetus
REAL( Kind= Real64 ) :: Hypotenuse                         ! Hypothenuse
REAL( Kind= Real64 ) :: drSquaredClosestDistanceParallel   ! Magnitude of vector T parallel to the cylindrical rim (i or j) (squared)
REAL( Kind= Real64 ) :: drSquaredClosestDistanceOrthogonal ! Magnitude of vector T orthogonal to the cylindrical rim (i or j) (squared)
REAL( Kind= Real64 ) :: drSquaredClosestDistance           ! Magnitude of the vector T (squared)
REAL( Kind= Real64 ) :: FuncR                              ! Objective function (root)
REAL( Kind= Real64 ) :: LambdaR                            ! Minimization variable related to a point on the cylindrical rim (i or j) (root)
REAL( Kind= Real64 ) :: LambdaC                            ! Minimization variable related to a point on the cylindrical rim (i or j) (root guess) [c]
REAL( Kind= Real64 ) :: LambdaD                            ! Minimization variable related to a point on the cylindrical rim (i or j) (root guess) [d]
REAL( Kind= Real64 ) :: FuncC                              ! Objective function (midpoint)
REAL( Kind= Real64 ) :: rSquaredRadius                     ! Radius of the cylindrical rim (squared)
REAL( Kind= Real64 ) :: dSquaredRadius                     ! Radius of the cylindrical disk (squared)
REAL( Kind= Real64 ) :: rSquaredLength                     ! Half-length of the cylindrical rim (squared)
REAL( Kind= Real64 ) :: drSquaredDiameter                  ! Combined diameter of both cylinders (squared)
REAL( Kind= Real64 ) :: dHalfDiameter                      ! Half diameter of the cylindrical disk
REAL( Kind= Real64 ) :: rHalfDiameter                      ! Half diameter of the cylindrical rim
REAL( Kind= Real64 ) :: rHalfLength                        ! Half length of the cylindrical rim
REAL( Kind= Real64 ) :: dSquaredVectorDistanceIntersection ! Distance between the center of the cylindrical disk and the intersection point (squared)
REAL( Kind= Real64 ) :: rSquaredVectorDistanceIntersection ! Distance between the center of the cylindrical rim and the intersection point (squared)
REAL( Kind= Real64 ) :: Alpha, Beta, Gamma                 ! Constants
REAL( Kind= Real64 ) :: CardanoQ, CardanoR                 ! Cardano's parameters
REAL( Kind= Real64 ) :: CardanoS, CardanoT                 ! Cardano's parameters
REAL( Kind= Real64 ) :: Theta, CosineArgument              ! Cardano's parameters
REAL( Kind= Real64 ) :: Discriminant                       ! Cardano's discriminant
REAL( Kind= Real64 ) :: MaxCoefficient                     ! Largest coefficient (absolute value) of the quartic function
REAL( Kind= Real64 ) :: ScaleParameter                     ! Scale parameter of the line equation that demarks the intersection point between the line and the disk
REAL( Kind= Real64 ) :: TempVar                            ! Temporary variable

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 3 )   :: CubicRoots                 ! Roots of the cubic equation
REAL( Kind= Real64 ), DIMENSION( 3 )   :: Func                       ! Objective function
REAL( Kind= Real64 ), DIMENSION( 3 )   :: Lambda                     ! Minimization variable related to a point on the cylindrical rim (i or j)
REAL( Kind= Real64 ), DIMENSION( 3 )   :: pRimPosition               ! Position of the cylindrical rim (i or j)
REAL( Kind= Real64 ), DIMENSION( 3 )   :: pRimOrientation            ! Orientation of the cylindrical rim (i or j)
REAL( Kind= Real64 ), DIMENSION( 3 )   :: pDiskOrientationX          ! Orientation of the cylindrical disks along x-direction (j or i)
REAL( Kind= Real64 ), DIMENSION( 3 )   :: pDiskOrientationY          ! Orientation of the cylindrical disks along y-direction (j or i)
REAL( Kind= Real64 ), DIMENSION( 3 )   :: pDiskOrientationZ          ! Orientation of the cylindrical disks along z-direction (j or i)
REAL( Kind= Real64 ), DIMENSION( 3 )   :: pDiskPosition              ! Position of cylindrical disks (j or i)
REAL( Kind= Real64 ), DIMENSION( 3 )   :: drVectorDistance           ! Vector distance between a cylindrical disk (j or i) and a cylindrical rim (i or j)
REAL( Kind= Real64 ), DIMENSION( 3 )   :: rClosestVectorDistanceDisk ! Vector distance between a cylindrical disk (j or i) and the closest point on the cylindrical rim (i or j)
REAL( Kind= Real64 ), DIMENSION( 3 )   :: rClosestPointDisk          ! Closest point on the cylindrical rim (i or j) with respect to the center of the cylindrical disk (j or i)
REAL( Kind= Real64 ), DIMENSION( 3 )   :: rClosestPoint              ! Closest point on the axis of the cylindrical rim (i or j) with respect to the circumference of the cylindrical disk (j or i)
REAL( Kind= Real64 ), DIMENSION( 3 )   :: dClosestPoint              ! Closest point on the circumference of the cylindrical disk (j or i) with respect to the axis of the cylindrical rim (i or j)
REAL( Kind= Real64 ), DIMENSION( 3 )   :: drClosestDistance          ! Vector distance between the closest point on the circumference of the cylindrical disk (j or i) and the center of the cylindrical rim (i or j)
REAL( Kind= Real64 ), DIMENSION( 3 )   :: dPoint                     ! Point in a plane defined by the cylindrical disk
REAL( Kind= Real64 ), DIMENSION( 3 )   :: rPointA, rPointB           ! Points in a line defined by the cylindrical rim
REAL( Kind= Real64 ), DIMENSION( 3 )   :: Intersection               ! Intersection point between the cylindrical axis and cylindrical disk
REAL( Kind= Real64 ), DIMENSION( 5 )   :: cQuartic                   ! Coefficients of the quartic function
REAL( Kind= Real64 ), DIMENSION( 4 )   :: cCubic                     ! Coefficients of the cubic function
REAL( Kind= Real64 ), DIMENSION( 3 )   :: cRoots                     ! Points where the derivative of the objective function with respect to λ are 0
REAL( Kind= Real64 ), DIMENSION( 0:3 ) :: pDiskQuaternion            ! Quaternion of the cylindrical disks (j or i)

! LOGICAL VARIABLES
LOGICAL :: OverlapDiskRim   ! Detects overlap between two particles (disk-rim configuration) : TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: BisectionLogical ! Checks whether the root from the bisection method will be used or not (Brent's method)

! Initialization
OverlapDiskRim   = .FALSE.
BisectionLogical = .FALSE.

! Initialize arrays and variables
nIntervals = 0
cRoots     = 0.D0
Func       = 0.D0
Lambda     = 0.D0
LambdaR    = 0.D0

! Radius of the cylindrical rim (squared)
rSquaredRadius = ( rHalfDiameter * rHalfDiameter )

! Radius of the cylindrical disk (squared)
dSquaredRadius = ( dHalfDiameter * dHalfDiameter )

! Half length of the cylindrical rim (squared)
rSquaredLength = ( rHalfLength * rHalfLength )

! Combined diameter of both cylinders (squared)
drSquaredDiameter = ( dHalfDiameter + rHalfDiameter )
drSquaredDiameter = ( drSquaredDiameter * drSquaredDiameter )

! Vector distance between a cylindrical disk and a cylindrical rim
drVectorDistance(1) = pDiskPosition(1) - pRimPosition(1)
drVectorDistance(2) = pDiskPosition(2) - pRimPosition(2)
drVectorDistance(3) = pDiskPosition(3) - pRimPosition(3)

! Projection of the vector distance between a cylindrical disk and a cylindrical rim along the orientation of the cylindrical rim
drVectorDistanceOrientationZ = DOT_PRODUCT( drVectorDistance, pRimOrientation )

! Closest point on a cylindrical rim with respect to the center of a cylindrical disk
rClosestPointDisk(1) = pRimPosition(1) + ( drVectorDistanceOrientationZ * pRimOrientation(1) )
rClosestPointDisk(2) = pRimPosition(2) + ( drVectorDistanceOrientationZ * pRimOrientation(2) )
rClosestPointDisk(3) = pRimPosition(3) + ( drVectorDistanceOrientationZ * pRimOrientation(3) )

! Vector distance between a cylindrical disk and the closest point on a cylindrical rim
rClosestVectorDistanceDisk(1) = pDiskPosition(1) - rClosestPointDisk(1)
rClosestVectorDistanceDisk(2) = pDiskPosition(2) - rClosestPointDisk(2)
rClosestVectorDistanceDisk(3) = pDiskPosition(3) - rClosestPointDisk(3)

! Magnitude of the vector distance between the closest point on the cylindrical rim and the center of the cylindrical disk
rSquaredClosestVectorDistanceDisk = DOT_PRODUCT( rClosestVectorDistanceDisk, rClosestVectorDistanceDisk )

! *********************************************************************************************** !
! PRELIMINARY CONDITIONS                                                                          !
! *********************************************************************************************** !
! Spheres circumscribing the cylindrical disk and the cross-sectional area of the cylindrical rim
IF( rSquaredClosestVectorDistanceDisk > drSquaredDiameter ) THEN
  ! No overlap
  RETURN ! Return immediately
! Center of the cylindrical disk in the space zone that corresponds to the extension of the cylindrical rim
ELSE IF( ( rSquaredClosestVectorDistanceDisk <= rSquaredRadius ) .AND. ( DABS( drVectorDistanceOrientationZ ) > &
  &      rHalfLength ) ) THEN
  ! Disk-disk configuration (already tested)
  RETURN ! Return immediately
! Center of the cylindrical disk in the space zone that corresponds to the cylindrical rim
ELSE IF( ( rSquaredClosestVectorDistanceDisk <= rSquaredRadius ) .AND. ( DABS( drVectorDistanceOrientationZ ) <= &
  &      rHalfLength ) ) THEN
  OverlapDiskRim = .TRUE.
  RETURN ! Return immediately
END IF

! *********************************************************************************************** !
! OTHER CONFIGURATIONS (INITIALIZATION)                                                           !
! *********************************************************************************************** !
! Orientation of the cylindrical disk along the x-direction
CALL VectorRotation( xAxis, pDiskQuaternion, pDiskOrientationX )
! Orientation of the cylindrical disk along the y-direction
CALL VectorRotation( yAxis, pDiskQuaternion, pDiskOrientationY )
! Orientation of the cylindrical disk along the z-direction
CALL VectorRotation( zAxis, pDiskQuaternion, pDiskOrientationZ )
! Dot product of the orientation of the cylindrical disk along x-direction and the orientation of the cylindrical rim
drOrientationXZ = DOT_PRODUCT( pDiskOrientationX, pRimOrientation )
! Dot product of the orientation of the cylindrical disk along y-direction and the orientation of the cylindrical rim
drOrientationYZ = DOT_PRODUCT( pDiskOrientationY, pRimOrientation )
! Projection of the vector distance between a disk and a rim on the orientation of the cylindrical disk along x-direction
drVectorDistanceOrientationX = DOT_PRODUCT( drVectorDistance, pDiskOrientationX )
! Projection of the vector distance between a disk and a rim on the orientation of the cylindrical disk along y-direction
drVectorDistanceOrientationY = DOT_PRODUCT( drVectorDistance, pDiskOrientationY )

! *********************************************************************************************** !
! NEW SPECIAL CASE (MAINLY FOR ANISOMORPHIC CYLINDERS)                                            !
! *********************************************************************************************** !
! Any single point in the plane defined by the cylindrical disk (φ = 0)
dPoint = pDiskPosition + dHalfDiameter * pDiskOrientationX
! Any two points in the line defined by the cylindrical axis (λ = L/2 and λ = -L/2)
rPointA = pRimPosition + rHalfLength * pRimOrientation
rPointB = pRimPosition - rHalfLength * pRimOrientation
! Checking whether the cylindrical axis is orthogonal to the normal vector of the plane (avoid the indeterminate form 0/0)
IF( DABS( DOT_PRODUCT( pDiskOrientationZ, pRimOrientation ) * DOT_PRODUCT( pDiskOrientationZ, pRimOrientation ) - 0.D0 ) >= &
&   EPSILON( 1.D0 ) ) THEN
  ! Scale parameter of the line equation that demarks the intersection point between the line and the disk
  ScaleParameter = ( DOT_PRODUCT( pDiskOrientationZ, dPoint ) - DOT_PRODUCT( pDiskOrientationZ, rPointA ) ) / &
  &                DOT_PRODUCT( pDiskOrientationZ, (rPointB - rPointA) )
  ! Coordinates of the intersection point
  Intersection = rPointA + ScaleParameter * (rPointB - rPointA)
  ! Distance between the intersection point and the center of the cylindrical disk (squared)
  dSquaredVectorDistanceIntersection = DOT_PRODUCT( (pDiskPosition - Intersection), (pDiskPosition - Intersection) )
  ! Distance between the intersection point and the center of the cylindrical axis (squared)
  rSquaredVectorDistanceIntersection = DOT_PRODUCT( (pRimPosition - Intersection), (pRimPosition - Intersection) )
  ! Checks whether the intersection point is within the limits of both the cylindrical disk and cylindrical axis
  IF( dSquaredVectorDistanceIntersection <= dSquaredRadius .AND. rSquaredVectorDistanceIntersection <= rSquaredLength ) THEN
    OverlapDiskRim = .TRUE.
    RETURN ! Return immediately
  END IF
END IF

! *********************************************************************************************** !
! MORE COMPLEX CONFIGURATIONS                                                                     !
! *********************************************************************************************** !
!  Calculates a point on the circumference of the cylindrical disk and a point on the axis of the !
!  cylindrical rim such that the distance between them is a minimum.                              !
! *********************************************************************************************** !

! Constants
Alpha = (drOrientationXZ * drOrientationXZ) + (drOrientationYZ * drOrientationYZ)
Beta  = (drOrientationXZ * drVectorDistanceOrientationX) + (drOrientationYZ * drVectorDistanceOrientationY)
Gamma = (drVectorDistanceOrientationX * drVectorDistanceOrientationX) + &
&       (drVectorDistanceOrientationY * drVectorDistanceOrientationY)

! Particles not aligned through their molecular axis (end-to-end configuration => α = 0)
IF( Alpha > 0.D0 ) THEN

  ! Coefficients of the quartic function (objective function)
  cQuartic(1) = Alpha
  cQuartic(2) = - 2.D0 * (Beta + Alpha * drVectorDistanceOrientationZ)
  cQuartic(3) = (4.D0 * Beta * drVectorDistanceOrientationZ) + (Alpha * drVectorDistanceOrientationZ * &
  &             drVectorDistanceOrientationZ) - (Alpha * Alpha * dHalfDiameter * dHalfDiameter) + Gamma 
  cQuartic(4) = 2.D0 * ( (Alpha * Beta * dHalfDiameter * dHalfDiameter) - (Gamma * drVectorDistanceOrientationZ) - (Beta * &
  &             drVectorDistanceOrientationZ * drVectorDistanceOrientationZ) )
  cQuartic(5) = (drVectorDistanceOrientationZ * drVectorDistanceOrientationZ * Gamma) - (Beta * Beta * &
  &             dHalfDiameter * dHalfDiameter)

  ! Reduction of the quartic function (ensures a faster convergence of the numerical methods)
  MaxCoefficient = 0.D0
  DO cPolynomial = 1, 5
    IF( DABS( cQuartic(cPolynomial) ) >= MaxCoefficient ) THEN
      MaxCoefficient = cQuartic(cPolynomial)
    END IF
  END DO
  cQuartic = cQuartic / MaxCoefficient

  ! Coefficients of the cubic function (derivative of the objective function with respect to λ)
  cCubic(1) = 4.D0 * cQuartic(1)
  cCubic(2) = 3.D0 * cQuartic(2)
  cCubic(3) = 2.D0 * cQuartic(3)
  cCubic(4) = cQuartic(4)

  ! ********************************************************************************************* !
  ! Cardano's solution for the cubic function (finding the points where the derivative is zero)   !
  ! ********************************************************************************************* !
  
  ! Cardano's coefficients
  CardanoQ = ( ( 3.D0 * cCubic(1) * cCubic(3) ) - ( cCubic(2) * cCubic(2) ) ) / ( 9.D0 * cCubic(1) * cCubic(1) )
  CardanoR = ( ( 9.D0 * cCubic(1) * cCubic(2) * cCubic(3) ) - ( 27.D0 *cCubic(1) * cCubic(1) * cCubic(4) ) - &
  &          ( 2.D0 * cCubic(2) * cCubic(2) * cCubic(2) ) ) / ( 54.D0 * cCubic(1) * cCubic(1) * cCubic(1) )
  
  ! Cardano's discriminant
  Discriminant = CardanoQ * CardanoQ * CardanoQ + CardanoR * CardanoR
  
  ! If Δ < 0, all roots are real and unequal (Trigonometric Solution)
  IF( Discriminant < 0.D0 ) THEN
    ! Trigonometric function
    CosineArgument = CardanoR / DSQRT( - CardanoQ * CardanoQ * CardanoQ )
    IF( CosineArgument <= - 1.D0 ) THEN
      CosineArgument = - 1.D0
    ELSE IF( CosineArgument >= 1.D0 ) THEN
      CosineArgument = 1.D0
    END IF
    Theta = DACOS( CosineArgument )
    ! Roots of the cubic equation
    CubicRoots(1) = 2.D0 * DSQRT( - CardanoQ ) * DCOS( (Theta / 3.D0) ) - cCubic(2) / ( 3.D0 * cCubic(1) )
    CubicRoots(2) = 2.D0 * DSQRT( - CardanoQ ) * DCOS( (Theta / 3.D0) + (2.D0 * cPi / 3.D0) ) - cCubic(2) / ( 3.D0 * cCubic(1) )
    CubicRoots(3) = 2.D0 * DSQRT( - CardanoQ ) * DCOS( (Theta / 3.D0) + (4.D0 * cPi / 3.D0) ) - cCubic(2) / ( 3.D0 * cCubic(1) )
    ! Sort roots in ascending order
    cRoots(1) = MINVAL( CubicRoots )
    cRoots(3) = MAXVAL( CubicRoots )
    DO cPolynomial = 1, 3
      IF( CubicRoots(cPolynomial) > cRoots(1) .AND. CubicRoots(cPolynomial) < cRoots(3) ) THEN
        cRoots(2) = CubicRoots(cPolynomial)
      END IF
    END DO
  ! If Δ = 0, all roots are real and at least two are equal
  ELSE IF( DABS( Discriminant - 0.D0 ) < EPSILON( 1.D0 ) ) THEN
    ! Cardano's coefficients
    CardanoS = CardanoR ** (1.D0 / 3.D0)
    CardanoT = CardanoS
    ! Roots of the cubic equation
    CubicRoots(1) = CardanoS + CardanoT - cCubic(2) / ( 3.D0 * cCubic(1) )
    CubicRoots(2) = - 0.5D0 * (CardanoS + CardanoT) - cCubic(2) / ( 3.D0 * cCubic(1) )
    CubicRoots(3) = CubicRoots(2)
    ! Sort roots in ascending order
    cRoots(1) = MINVAL( CubicRoots )
    cRoots(2) = MAXVAL( CubicRoots )
    cRoots(3) = cRoots(2)
  ! If Δ > 0, one roots is real and two are complex conjugates
  ELSE IF( Discriminant > 0.D0 ) THEN
    ! Cardano's coefficients
    CardanoS = ( CardanoR + DSQRT( (CardanoQ * CardanoQ * CardanoQ) + (CardanoR * CardanoR) ) ) ** (1.D0 / 3.D0)
    CardanoT = ( CardanoR - DSQRT( (CardanoQ * CardanoQ * CardanoQ) + (CardanoR * CardanoR) ) ) ** (1.D0 / 3.D0)
    ! Roots of the cubic equation
    CubicRoots(1) = CardanoS + CardanoT - cCubic(2) / ( 3.D0 * cCubic(1) )
    CubicRoots(2) = CubicRoots(1) ! We are only interested in the real solution
    CubicRoots(3) = CubicRoots(1) ! We are only interested in the real solution
    ! Sort roots in ascending order
    cRoots(1) = CubicRoots(1)
    cRoots(2) = cRoots(1)
    cRoots(3) = cRoots(1)
  END IF

  ! Number of root-finding intervals (based on the value of Cardano's discriminant)
  IF( Discriminant < 0.D0 ) THEN
    nIntervals = 4
  ELSE IF( DABS( Discriminant - 0.D0 ) < EPSILON( 1.D0 ) ) THEN
    nIntervals = 3
  ELSE IF( Discriminant > 0.D0 ) THEN
    nIntervals = 2
  END IF
  ! Single interval
  IF( MINVAL( cRoots ) >= rHalfLength .OR. MAXVAL( cRoots ) <= - rHalfLength ) THEN
    nIntervals = 1
  END IF

  ! ********************************************************************************************* !
  ! Loops (Initial Guesses/Intervals)                                                             !
  ! ********************************************************************************************* !
  ! Number of intervals = 4                                                                       !
  !   First loop  :     -L/2      <= λ0 <= CUBIC_ROOT(1)                                          !
  !   Second loop : CUBIC_ROOT(1) <= λ0 <= CUBIC_ROOT(2)                                          !
  !   Third loop  : CUBIC_ROOT(2) <= λ0 <= CUBIC_ROOT(3)                                          !
  !   Fourth loop : CUBIC_ROOT(3) <= λ0 <=     L/2                                                !
  !                                                                                               !
  ! Number of intervals = 3                                                                       !
  !   First loop  :     -L/2      <= λ0 <= CUBIC_ROOT(1)                                          !
  !   Second loop : CUBIC_ROOT(1) <= λ0 <= CUBIC_ROOT(2)                                          !
  !   Third loop  : CUBIC_ROOT(2) <= λ0 <= L/2                                                    !
  !                                                                                               !
  ! Number of intervals = 2                                                                       !
  !   First loop  :     -L/2      <= λ0 <= CUBIC_ROOT(1)                                          !
  !   Second loop : CUBIC_ROOT(1) <= λ0 <= L/2                                                    !
  !                                                                                               !
  ! Number of intervals = 1, if MIN( ROOT ) >= L/2                                                !
  !   Single loop :     -L/2      <= λ0 <= MIN( ROOT )                                            !
  !                                                                                               !
  ! Number of intervals = 1, if MAX( ROOT ) <= -L/2                                               !
  !   Single loop :  MAX( ROOT )  <= λ0 <= L/2                                                    !
  ! ********************************************************************************************* !
  RootFindingIntervals: DO rInterval = 1, nIntervals

    ! Ignore intervals where the function has the same sign at the extremum points (no real roots)
    IF( nIntervals == 4 ) THEN
      IF( rInterval == 1 ) THEN
        Lambda(2) = - rHalfLength
        Func(2)   = ( cQuartic(1) * Lambda(2) * Lambda(2) * Lambda(2) * Lambda(2) ) + ( cQuartic(2) * Lambda(2) * Lambda(2) * &
        &           Lambda(2) ) + ( cQuartic(3) * Lambda(2) * Lambda(2) ) + ( cQuartic(4) * Lambda(2) ) + cQuartic(5)
        Lambda(3) = cRoots(1)
        Func(3)   = ( cQuartic(1) * Lambda(3) * Lambda(3) * Lambda(3) * Lambda(3) ) + ( cQuartic(2) * Lambda(3) * Lambda(3) * &
        &           Lambda(3) ) + ( cQuartic(3) * Lambda(3) * Lambda(3) ) + ( cQuartic(4) * Lambda(3) ) + cQuartic(5)
        IF( ( Func(2) * Func(3) ) > 0.D0 ) THEN
          CYCLE RootFindingIntervals
        END IF
      ELSE IF( rInterval == 2 ) THEN
        Lambda(2) = cRoots(1)
        Func(2)   = ( cQuartic(1) * Lambda(2) * Lambda(2) * Lambda(2) * Lambda(2) ) + ( cQuartic(2) * Lambda(2) * Lambda(2) * &
        &           Lambda(2) ) + ( cQuartic(3) * Lambda(2) * Lambda(2) ) + ( cQuartic(4) * Lambda(2) ) + cQuartic(5)
        Lambda(3) = cRoots(2)
        Func(3)   = ( cQuartic(1) * Lambda(3) * Lambda(3) * Lambda(3) * Lambda(3) ) + ( cQuartic(2) * Lambda(3) * Lambda(3) * &
        &           Lambda(3) ) + ( cQuartic(3) * Lambda(3) * Lambda(3) ) + ( cQuartic(4) * Lambda(3) ) + cQuartic(5)
        IF( ( Func(2) * Func(3) ) > 0.D0 ) THEN
          CYCLE RootFindingIntervals
        END IF
      ELSE IF( rInterval == 3 ) THEN
        Lambda(2) = cRoots(2)
        Func(2)   = ( cQuartic(1) * Lambda(2) * Lambda(2) * Lambda(2) * Lambda(2) ) + ( cQuartic(2) * Lambda(2) * Lambda(2) * &
        &           Lambda(2) ) + ( cQuartic(3) * Lambda(2) * Lambda(2) ) + ( cQuartic(4) * Lambda(2) ) + cQuartic(5)
        Lambda(3) = cRoots(3)
        Func(3)   = ( cQuartic(1) * Lambda(3) * Lambda(3) * Lambda(3) * Lambda(3) ) + ( cQuartic(2) * Lambda(3) * Lambda(3) * &
        &           Lambda(3) ) + ( cQuartic(3) * Lambda(3) * Lambda(3) ) + ( cQuartic(4) * Lambda(3) ) + cQuartic(5)
        IF( ( Func(2) * Func(3) ) > 0.D0 ) THEN
          CYCLE RootFindingIntervals
        END IF
      ELSE IF( rInterval == 4 ) THEN
        Lambda(2) = cRoots(3)
        Func(2)   = ( cQuartic(1) * Lambda(2) * Lambda(2) * Lambda(2) * Lambda(2) ) + ( cQuartic(2) * Lambda(2) * Lambda(2) * &
        &           Lambda(2) ) + ( cQuartic(3) * Lambda(2) * Lambda(2) ) + ( cQuartic(4) * Lambda(2) ) + cQuartic(5)
        Lambda(3) = rHalfLength
        Func(3)   = ( cQuartic(1) * Lambda(3) * Lambda(3) * Lambda(3) * Lambda(3) ) + ( cQuartic(2) * Lambda(3) * Lambda(3) * &
        &           Lambda(3) ) + ( cQuartic(3) * Lambda(3) * Lambda(3) ) + ( cQuartic(4) * Lambda(3) ) + cQuartic(5)
        IF( ( Func(2) * Func(3) ) > 0.D0 ) THEN
          CYCLE RootFindingIntervals
        END IF
      END IF
    ELSE IF( nIntervals == 3 ) THEN
      IF( rInterval == 1 ) THEN
        Lambda(2) = - rHalfLength
        Func(2)   = ( cQuartic(1) * Lambda(2) * Lambda(2) * Lambda(2) * Lambda(2) ) + ( cQuartic(2) * Lambda(2) * Lambda(2) * &
        &           Lambda(2) ) + ( cQuartic(3) * Lambda(2) * Lambda(2) ) + ( cQuartic(4) * Lambda(2) ) + cQuartic(5)
        Lambda(3) = cRoots(1)
        Func(3)   = ( cQuartic(1) * Lambda(3) * Lambda(3) * Lambda(3) * Lambda(3) ) + ( cQuartic(2) * Lambda(3) * Lambda(3) * &
        &           Lambda(3) ) + ( cQuartic(3) * Lambda(3) * Lambda(3) ) + ( cQuartic(4) * Lambda(3) ) + cQuartic(5)
        IF( ( Func(2) * Func(3) ) > 0.D0 ) THEN
          CYCLE RootFindingIntervals
        END IF
      ELSE IF( rInterval == 2 ) THEN
        Lambda(2) = cRoots(1)
        Func(2)   = ( cQuartic(1) * Lambda(2) * Lambda(2) * Lambda(2) * Lambda(2) ) + ( cQuartic(2) * Lambda(2) * Lambda(2) * &
        &           Lambda(2) ) + ( cQuartic(3) * Lambda(2) * Lambda(2) ) + ( cQuartic(4) * Lambda(2) ) + cQuartic(5)
        Lambda(3) = cRoots(2)
        Func(3)   = ( cQuartic(1) * Lambda(3) * Lambda(3) * Lambda(3) * Lambda(3) ) + ( cQuartic(2) * Lambda(3) * Lambda(3) * &
        &           Lambda(3) ) + ( cQuartic(3) * Lambda(3) * Lambda(3) ) + ( cQuartic(4) * Lambda(3) ) + cQuartic(5)
        IF( ( Func(2) * Func(3) ) > 0.D0 ) THEN
          CYCLE RootFindingIntervals
        END IF
      ELSE IF( rInterval == 3 ) THEN
        Lambda(2) = cRoots(2)
        Func(2)   = ( cQuartic(1) * Lambda(2) * Lambda(2) * Lambda(2) * Lambda(2) ) + ( cQuartic(2) * Lambda(2) * Lambda(2) * &
        &           Lambda(2) ) + ( cQuartic(3) * Lambda(2) * Lambda(2) ) + ( cQuartic(4) * Lambda(2) ) + cQuartic(5)
        Lambda(3) = rHalfLength
        Func(3)   = ( cQuartic(1) * Lambda(3) * Lambda(3) * Lambda(3) * Lambda(3) ) + ( cQuartic(2) * Lambda(3) * Lambda(3) * &
        &           Lambda(3) ) + ( cQuartic(3) * Lambda(3) * Lambda(3) ) + ( cQuartic(4) * Lambda(3) ) + cQuartic(5)
        IF( ( Func(2) * Func(3) ) > 0.D0 ) THEN
          CYCLE RootFindingIntervals
        END IF
      END IF
    ELSE IF( nIntervals == 2 ) THEN
      IF( rInterval == 1 ) THEN
        Lambda(2) = - rHalfLength
        Func(2)   = ( cQuartic(1) * Lambda(2) * Lambda(2) * Lambda(2) * Lambda(2) ) + ( cQuartic(2) * Lambda(2) * Lambda(2) * &
        &           Lambda(2) ) + ( cQuartic(3) * Lambda(2) * Lambda(2) ) + ( cQuartic(4) * Lambda(2) ) + cQuartic(5)
        Lambda(3) = cRoots(1)
        Func(3)   = ( cQuartic(1) * Lambda(3) * Lambda(3) * Lambda(3) * Lambda(3) ) + ( cQuartic(2) * Lambda(3) * Lambda(3) * &
        &           Lambda(3) ) + ( cQuartic(3) * Lambda(3) * Lambda(3) ) + ( cQuartic(4) * Lambda(3) ) + cQuartic(5)
        IF( ( Func(2) * Func(3) ) > 0.D0 ) THEN
          CYCLE RootFindingIntervals
        END IF
      ELSE IF( rInterval == 2 ) THEN
        Lambda(2) = cRoots(1)
        Func(2)   = ( cQuartic(1) * Lambda(2) * Lambda(2) * Lambda(2) * Lambda(2) ) + ( cQuartic(2) * Lambda(2) * Lambda(2) * &
        &           Lambda(2) ) + ( cQuartic(3) * Lambda(2) * Lambda(2) ) + ( cQuartic(4) * Lambda(2) ) + cQuartic(5)
        Lambda(3) = rHalfLength
        Func(3)   = ( cQuartic(1) * Lambda(3) * Lambda(3) * Lambda(3) * Lambda(3) ) + ( cQuartic(2) * Lambda(3) * Lambda(3) * &
        &           Lambda(3) ) + ( cQuartic(3) * Lambda(3) * Lambda(3) ) + ( cQuartic(4) * Lambda(3) ) + cQuartic(5)
        IF( ( Func(2) * Func(3) ) > 0.D0 ) THEN
          CYCLE RootFindingIntervals
        END IF
      END IF
    ELSE IF( nIntervals == 1 .AND. MINVAL( cRoots ) >= rHalfLength ) THEN
      Lambda(2) = - rHalfLength
      Func(2)   = ( cQuartic(1) * Lambda(2) * Lambda(2) * Lambda(2) * Lambda(2) ) + ( cQuartic(2) * Lambda(2) * Lambda(2) * &
      &           Lambda(2) ) + ( cQuartic(3) * Lambda(2) * Lambda(2) ) + ( cQuartic(4) * Lambda(2) ) + cQuartic(5)
      Lambda(3) = MINVAL( cRoots )
      Func(3)   = ( cQuartic(1) * Lambda(3) * Lambda(3) * Lambda(3) * Lambda(3) ) + ( cQuartic(2) * Lambda(3) * Lambda(3) * &
      &           Lambda(3) ) + ( cQuartic(3) * Lambda(3) * Lambda(3) ) + ( cQuartic(4) * Lambda(3) ) + cQuartic(5)
      IF( ( Func(2) * Func(3) ) > 0.D0 ) THEN
        CYCLE RootFindingIntervals
      END IF
    ELSE IF( nIntervals == 1 .AND. MAXVAL( cRoots ) <= - rHalfLength ) THEN
      Lambda(2) = MAXVAL( cRoots )
      Func(2)   = ( cQuartic(1) * Lambda(2) * Lambda(2) * Lambda(2) * Lambda(2) ) + ( cQuartic(2) * Lambda(2) * Lambda(2) * &
      &           Lambda(2) ) + ( cQuartic(3) * Lambda(2) * Lambda(2) ) + ( cQuartic(4) * Lambda(2) ) + cQuartic(5)
      Lambda(3) = rHalfLength
      Func(3)   = ( cQuartic(1) * Lambda(3) * Lambda(3) * Lambda(3) * Lambda(3) ) + ( cQuartic(2) * Lambda(3) * Lambda(3) * &
      &           Lambda(3) ) + ( cQuartic(3) * Lambda(3) * Lambda(3) ) + ( cQuartic(4) * Lambda(3) ) + cQuartic(5)
      IF( ( Func(2) * Func(3) ) > 0.D0 ) THEN
        CYCLE RootFindingIntervals
      END IF
    END IF

    ! ******************************************************************************************* !
    ! Brent's method                                                                              !
    ! ******************************************************************************************* !
    
    ! Initialization
    BisectionLogical = .FALSE. ! Bisection root
    cBrent           = 0       ! Iteration counter
    Tolerance        = 1.D-10  ! Numerical tolerance
    LambdaD          = 0.D0    ! Third-to-last root guess

    ! Brent's condition
    IF( ( Func(2) * Func(3) ) <= 0.D0 ) THEN
      ! Swap condition
      IF( DABS( Func(2) ) < DABS( Func(3) ) ) THEN
        ! Swap bounds
        TempVar   = Lambda(2)
        Lambda(2) = Lambda(3)
        Lambda(3) = TempVar
        ! Swap function values
        TempVar = Func(2)
        Func(2) = Func(3)
        Func(3) = TempVar
      END IF
      ! Initialization of previous bound
      LambdaC = Lambda(2)
      FuncC   = Func(2)
      ! Initialize value of objective function
      FuncR = - 1.D0
      ! Initialize 'BISECTION' flag as TRUE since the last iteration used was the bisection method
      BisectionLogical = .TRUE.
      ! Stop criterion
      DO WHILE( DABS( Func(3) ) >= Tolerance .AND. DABS( FuncR ) >= Tolerance .AND. &
        &       DABS( ( Lambda(2) - Lambda(3) ) / Lambda(2) ) >= Tolerance )
        ! Initialize root of the function
        IF( DABS( Func(2) - FuncC ) >= EPSILON( 1.D0 ) .AND. DABS( Func(3) - FuncC ) >= EPSILON( 1.D0 ) ) THEN
          ! Inverse quadratic interpolation root-finding procedure
          LambdaR = ( Lambda(2) * Func(3) * FuncC ) / ( ( Func(2) - Func(3) ) * (Func(2) - FuncC) ) + &
          &         ( Lambda(3) * Func(2) * FuncC ) / ( ( Func(3) - Func(2) ) * ( Func(3) - FuncC ) ) + &
          &         ( LambdaC * Func(2) * Func(3) ) / ( ( FuncC - Func(2) ) * ( FuncC - Func(3) ) )
        ELSE
          ! False position formula
          LambdaR = Lambda(3) - Func(3) * ( Lambda(3) - Lambda(2) ) / ( Func(3) - Func(2) )
        END IF
        ! Check whether the root obtained from the interpolation method or false position formula will be used; otherwise, use the midpoint from bisection method
        IF( ( (LambdaR - ( 3.D0 * Lambda(2) + Lambda(3) ) / 4.D0) * ( LambdaR - Lambda(3) ) >= 0.D0 ) .OR. &
        &   ( BisectionLogical .AND. ( DABS( LambdaR - Lambda(3) ) >= DABS( (Lambda(3) - LambdaC) / 2.D0 ) ) ) .OR. &
        &   ( .NOT. BisectionLogical .AND. ( DABS( LambdaR - Lambda(3) ) >= DABS( (LambdaC - LambdaD) / 2.D0 ) ) ) .OR. &
        &   ( BisectionLogical .AND. ( DABS( LambdaR - LambdaC ) < Tolerance ) ) .OR. ( .NOT. BisectionLogical .AND. &
        &   ( DABS( LambdaC - LambdaD ) < Tolerance ) ) ) THEN
          ! Root from bisection method
          BisectionLogical = .TRUE.
          LambdaR          = 0.5D0 * ( Lambda(2) + Lambda(3) )
        ELSE
          ! Root from inverse quadratic interpolation or false position formula
          BisectionLogical = .FALSE.
        END IF
        ! Calculate function at the midpoint
        FuncR  = ( cQuartic(1) * LambdaR * LambdaR * LambdaR * LambdaR ) + ( cQuartic(2) * LambdaR * LambdaR * LambdaR ) + &
        &        ( cQuartic(3) * LambdaR * LambdaR ) + ( cQuartic(4) * LambdaR ) + cQuartic(5)
        ! Set third-to-last root guess
        LambdaD = LambdaC
        ! Set second-to-last root guess
        LambdaC = Lambda(3)
        FuncC   = Func(3)
        ! Check signs of functions
        IF( Func(2) * FuncR < 0.D0 ) THEN
          Lambda(3) = LambdaR
          Func(3)   = FuncR
        ELSE
          Lambda(2) = LambdaR
          Func(2)   = FuncR
        END IF
        ! Swap condition
        IF( DABS( Func(2) ) < DABS( Func(3) ) ) THEN
          ! Swap bounds
          TempVar   = Lambda(2)
          Lambda(2) = Lambda(3)
          Lambda(3) = TempVar
          ! Swap function values
          TempVar = Func(2)
          Func(2) = Func(3)
          Func(3) = TempVar
        END IF
        ! Iteration
        cBrent = cBrent + 1
        IF( cBrent > 100 ) THEN
          WRITE( *, "(3G0)" ) "Brent's method failed to converge after ", cBrent, " iterations. Exiting..."
          STOP
        END IF
      END DO
      ! Root (iterative process or initial guess)
      IF( DABS( FuncR ) < Tolerance ) THEN
        Lambda(1) = LambdaR
      ELSE IF( DABS( Func(3) ) < Tolerance .OR. DABS( ( Lambda(2) - Lambda(3) ) / Lambda(2) ) < Tolerance ) THEN
        Lambda(1) = Lambda(3)
      END IF
    ! Algorithm should not reach this point
    ELSE
      ! No real roots within this interval
      CYCLE RootFindingIntervals
    END IF

    ! Point that minimizes the objective function is outside the cylindrical rim
    IF( Lambda(1) < -rHalfLength .OR. Lambda(1) > rHalfLength) THEN
      CYCLE RootFindingIntervals
    END IF

    ! Find closest points on the circumference of the disk and the axis of the cylinder
    OppositeCathetus = Lambda(1) * drOrientationYZ - drVectorDistanceOrientationY 
    AdjacentCathetus = Lambda(1) * drOrientationXZ - drVectorDistanceOrientationX 
    Hypotenuse       = DSQRT( (OppositeCathetus * OppositeCathetus) + (AdjacentCathetus * AdjacentCathetus) )
    ! Trigonometric relations in the rectangle triangle
    CosinePhi = AdjacentCathetus / Hypotenuse
    SinePhi   = OppositeCathetus / Hypotenuse
    ! Closest point on the circumference of the cylindrical disk
    dClosestPoint(1) = pDiskPosition(1) + ( dHalfDiameter * CosinePhi * pDiskOrientationX(1) ) + &
    &                  ( dHalfDiameter * SinePhi * pDiskOrientationY(1) )
    dClosestPoint(2) = pDiskPosition(2) + ( dHalfDiameter * CosinePhi * pDiskOrientationX(2) ) + &
    &                  ( dHalfDiameter * SinePhi * pDiskOrientationY(2) )
    dClosestPoint(3) = pDiskPosition(3) + ( dHalfDiameter * CosinePhi * pDiskOrientationX(3) ) + &
    &                  ( dHalfDiameter * SinePhi * pDiskOrientationY(3) )
    ! Closest point on the axis of the cylindrical rim
    rClosestPoint(1) = pRimPosition(1) + ( Lambda(1) * pRimOrientation(1) )
    rClosestPoint(2) = pRimPosition(2) + ( Lambda(1) * pRimOrientation(2) )
    rClosestPoint(3) = pRimPosition(3) + ( Lambda(1) * pRimOrientation(3) )
    ! Vector distance between the closest point on the circumference of a cylindrical disk and the center of a cylindrical rim
    drClosestDistance(1) = dClosestPoint(1) - pRimPosition(1)
    drClosestDistance(2) = dClosestPoint(2) - pRimPosition(2)
    drClosestDistance(3) = dClosestPoint(3) - pRimPosition(3)
    ! Magnitude of the vector T (squared)
    drSquaredClosestDistance = DOT_PRODUCT( drClosestDistance, drClosestDistance )
    ! Magnitude of vector T parallel to the cylindrical rim (squared)
    drSquaredClosestDistanceParallel = DOT_PRODUCT( pRimOrientation, drClosestDistance )
    drSquaredClosestDistanceParallel = drSquaredClosestDistanceParallel * drSquaredClosestDistanceParallel
    ! Magnitude of vector T orthogonal to the cylindrical rim (squared)
    drSquaredClosestDistanceOrthogonal = drSquaredClosestDistance - drSquaredClosestDistanceParallel
    ! Overlap condition
    IF( drSquaredClosestDistanceParallel <= rSquaredLength .AND. drSquaredClosestDistanceOrthogonal <= rSquaredRadius ) THEN
      OverlapDiskRim = .TRUE.
      RETURN ! Return immediately
    ELSE
      OverlapDiskRim = .FALSE.
    END IF

  END DO RootFindingIntervals

END IF ! The end-to-end configuration have already been tested by other methods (α = 0)

! *********************************************************************************************** !
! -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*  IMPORTANT NOTES  -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* !
! *********************************************************************************************** !
! IF THE ALGORITHM REACHES THIS POINT:                                                            !
!                                                                                                 !
! (1) it is because the root obtained from the numerical method returns a point in the            !
! circumference of the cylindrical disk (j or i) that is outside the cylindrical rim (i or j) but !
! within the radial extension of its length. RESULT = NO OVERLAP                                  !
!                                                                                                 !
! (2) it is because the root obtained from the numerical method returns a point in the            !
! circumference of the cylindrical disk (j or i) that is outside the cylindrical rim (i or j) but !
! within the axial extension of its length. RESULT = NO OVERLAP                                   !
!                                                                                                 !
! (3) it is because there are no real roots within the analyzed intervals. RESULT = NO OVERLAP    !
!                                                                                                 !
! (4) it is because no root was found in the considered intervals (-L/2 <= λ0 <= L/2), which      !
! means the root that minimizes the objective function is beyond the limits of the length of      !
! the cylindrical rim (i or j). This situation will never lead to molecular overlaps except when  !
! the disks of both cylinders intersect each other. However, this situation has already been      !
! tested in the disk-disk configuration. RESULT = NO OVERLAP                                      !
!                                                                                                 !
! (5) it is because α = 0, but this situation has already been tested by other methods.           !
! RESULT = NO OVERLAP                                                                             !
! *********************************************************************************************** !

! No overlaps
RETURN

END SUBROUTINE DiskRimConfiguration

END MODULE OverlapAlgorithms