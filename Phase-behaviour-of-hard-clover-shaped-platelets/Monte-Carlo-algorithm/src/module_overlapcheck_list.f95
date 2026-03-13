MODULE OverlapCheckLists

! Uses three modules: global variables, overlap algorithms, and vector operations
USE GlobalVariables
USE OverlapAlgorithms, ONLY: OverlapCheckSPC, OverlapCheckCYL
USE VectorOperations

IMPLICIT NONE

! Make subroutines public
PUBLIC :: ListOverlapCheck, FullListOverlapCheck, ListOverlapCheckInitialConfiguration, FullListOverlapCheckInitialConfiguration

! *********************************************************************************************** !
!                                      INITIAL CONFIGURATION                                      !
! *********************************************************************************************** !
CONTAINS

! *********************************************************************************************** !
!         This subroutine uses lists to check whether there are any overlapping particles         !
!                                  in the initial configuration                                   !
! *********************************************************************************************** !
SUBROUTINE FullListOverlapCheckInitialConfiguration( ContactDistance, CurrentBoxLength, CurrentBoxLengthInverse, HalfNeighbours )

! Uses two modules: linked lists and overlap check algorithms
USE LinkedLists, ONLY: MakeList

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 )           :: iParticle, jParticle ! Counter (particle)
INTEGER( Kind= Int64 ), OPTIONAL :: HalfNeighbours               ! Checks whether a cell and its 26 surrounding cells are searched (PRESENT) or just its 13 neighbour cells (NOT PRESENT)

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 ) :: ContactDistance ! Contact distance (Perram-Wertheim or Vega-Lago methods)

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: BoxCutoff               ! Box cutoff (x-, y-, and z-directions)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: iPosition               ! Position of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )    :: iOrientation            ! Orientation of particles i and j
REAL( Kind= Real64 ), DIMENSION( 9 )    :: CurrentBoxLength        ! Box length
REAL( Kind= Real64 ), DIMENSION( 9 )    :: CurrentBoxLengthInverse ! Box length (inverse)
REAL( Kind= Real64 ), DIMENSION( 0:3 )  :: iQuaternion             ! Quaternions of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3, 4 ) :: ciPosition              ! Position of particles i and j

! LOGICAL VARIABLES
LOGICAL :: Overlap         ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: CellHalfLogical ! Checks whether a cell and its 26 surrounding cells are searched (FALSE) or just its 13 neighbour cells (TRUE)

! Initialization
Overlap         = .FALSE.
CellHalfLogical = .FALSE.

! Re-initialize cell list (if necessary)
BoxCutoff(1) = cLargestSphereDiameter / CurrentBoxLength(1)
BoxCutoff(2) = cLargestSphereDiameter / CurrentBoxLength(5)
BoxCutoff(3) = cLargestSphereDiameter / CurrentBoxLength(9)
CALL MakeList( BoxCutoff, pPosition, CurrentBoxLengthInverse )

! Check whether the number of cells along any direction (x, y, or z) is less than 3
IF( .NOT. CellListControl ) RETURN

! Consider all 26 neighbours or only 13
CellHalfLogical = PRESENT( HalfNeighbours )

! Loop over all particles (consider half of neighbours)
DO iParticle = 1, nParticles
  ! Quaternion of particle i
  iQuaternion = pQuaternion(:,iParticle)
  ! Orientation of particle i
  iOrientation = pOrientation(:,iParticle)
  ! Position of particle i
  iPosition = pPosition(:,iParticle)
  ! Position of particle i
  ciPosition = cPosition(:,:,iParticle)
  ! Overlap check between particle i and its neighbours
  CALL ListOverlapCheckInitialConfiguration( iParticle, iQuaternion, iOrientation, iPosition, ciPosition, ContactDistance, &
  &                                          CurrentBoxLength, CurrentBoxLengthInverse, Overlap, jParticle, CellHalfLogical )
  ! Overlap detected
  IF( Overlap ) THEN
    ! Overlap detected
    WRITE( *, "(5G0)" ) "Overlap found in the initial configuration between particles ", iParticle, " and ", &
    &                   jParticle, ". Exiting..."
    STOP
  END IF
END DO

! No overlaps
RETURN

END SUBROUTINE FullListOverlapCheckInitialConfiguration

! *********************************************************************************************** !
!         This subroutine uses lists to check whether there are any overlapping particles         !
!                                  in the initial configuration                                   !
! *********************************************************************************************** !
SUBROUTINE ListOverlapCheckInitialConfiguration( iParticle, iQuaternion, iOrientation, iPosition, ciPosition, ContactDistance, &
&                                                CurrentBoxLength, CurrentBoxLengthInverse, Overlap, jParticle, CellHalfLogical )

! Uses two modules: linked lists and overlap check algorithms
USE LinkedLists, ONLY: Neighbours, CellIndex

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 ) :: iCylinder, jCylinder ! Counters (particle)
INTEGER( Kind= Int64 ) :: iParticle, jParticle ! Counters (particle)
INTEGER( Kind= Int64 ) :: jList                ! Counters (neighbour list)

! INTEGER VARIABLES (ARRAY)
INTEGER( Kind= Int64 ), DIMENSION( 3 )          :: iCellIndex     ! 3D cell index
INTEGER( Kind= Int64 ), DIMENSION( nParticles ) :: jNeighbourList ! List of neighbour j particles

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 ) :: SquaredDistance ! Magnitude of the vector distance between particles i and j (squared)
REAL( Kind= Real64 ) :: ContactDistance ! Contact distance (Perram-Wertheim or Vega-Lago methods)

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: VectorDistance             ! Vector distance between particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )    :: iPosition, jPosition       ! Position of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )    :: iOrientation, jOrientation ! Orientation of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )    :: ScalingDistanceUnitBox     ! Scaled position (unit box)
REAL( Kind= Real64 ), DIMENSION( 9 )    :: CurrentBoxLength           ! Box length
REAL( Kind= Real64 ), DIMENSION( 9 )    :: CurrentBoxLengthInverse    ! Box length (inverse)
REAL( Kind= Real64 ), DIMENSION( 0:3 )  :: iQuaternion, jQuaternion   ! Quaternions of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3, 4 ) :: ciPosition, cjPosition     ! Position of particles i and j

! LOGICAL VARIABLES
LOGICAL :: OverlapSPC      ! Detects overlap between two spherocylinders: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: OverlapCYL      ! Detects overlap between two cylinders: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: ParallelSPC     ! Checks the relative orientation of two spherocylinders : TRUE = parallel orientation; FALSE = non-parallel orientation
LOGICAL :: Overlap         ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: CellHalfLogical ! Checks whether a cell and its 26 surrounding cells are searched (FALSE) or just its 13 neighbour cells (TRUE)

! Initialization
jList       = 0
Overlap     = .FALSE.
OverlapCYL  = .FALSE.
ParallelSPC = .FALSE.
OverlapSPC  = .FALSE.

! Spatial transformation
CALL MatrixVectorMultiplication( CurrentBoxLengthInverse, iPosition, ScalingDistanceUnitBox )

! Minimum image convention
ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )

! Index function allocating particle i to cell
iCellIndex = CellIndex( ScalingDistanceUnitBox )

! Debug check
IF( CellHalfLogical .AND. ANY( iCellIndex(:) /= pCellIndex(:,iParticle) ) ) THEN ! Should never happen
  WRITE( *, "(G0,2(G0,', '),2G0,2(G0,', '),2G0)" ) "Cell mismatch: [", iCellIndex(:), "] and [", pCellIndex(:,iParticle), "]"
  STOP
END IF

! List of neighbours
jNeighbourList = Neighbours( iParticle, iCellIndex, CellHalfLogical )

! Loop until no more entries in the neighbour list
DO
  ! Next entry
  jList = jList + 1
  ! Neighbour particle index
  jParticle = jNeighbourList(jList)
  ! List exhausted
  IF( jParticle == 0 ) EXIT
  ! Skip self (should never happen)
  IF( iParticle == jParticle ) CYCLE
  ! Position of particle j
  jPosition(:) = pPosition(:,jParticle)
  ! Position of particle j
  cjPosition(:,:) = cPosition(:,:,jParticle)
  ! Orientation of particle j
  jOrientation(:) = pOrientation(:,jParticle)
  ! Quaternion of particle j
  jQuaternion(:) = pQuaternion(:,jParticle)
  ! Vector distance between particles i and j
  VectorDistance(:) = jPosition(:) - iPosition(:)
  ! Minimum image convention
  CALL MatrixVectorMultiplication( CurrentBoxLengthInverse, VectorDistance, ScalingDistanceUnitBox ) ! Spatial transformation
  ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
  CALL MatrixVectorMultiplication( CurrentBoxLength, ScalingDistanceUnitBox, VectorDistance ) ! Spatial transformation
  ! Magnitude of the vector distance (squared)
  SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
  ! Initialization
  OverlapSPC = .FALSE.
  OverlapCYL = .FALSE.
  ! Cutoff distance (sphere circumscribing the non-convex body)
  IF( SquaredDistance <= SquaredCutoffSphere ) THEN
    ! Spherocylinder circumscribing the non-convex body
    CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, ContactDistance, ParallelSPC, OverlapSPC, &
    &                     .TRUE. )
    ! Vega-Lago criterion
    IF( OverlapSPC ) THEN
      OverlapCYL = .TRUE. ! Check overlap between cylinders
    ELSE
      CYCLE
    END IF
  END IF
  ! Test overlap between cylinders
  IF( OverlapCYL ) THEN
    ! First loop takes one of the four cylinders from particle i
    DO iCylinder = 1, 4
      ! Second loop takes one of the four cylinders from particle j
      DO jCylinder = 1, 4
        ! Vector distance between particles i and j
        VectorDistance(:) = cjPosition(:,jCylinder) - ciPosition(:,iCylinder)
        ! Minimum Image Convention (see Allen and Tildesley, 2nd Edition (2017), pages 35-45)
        CALL MatrixVectorMultiplication( CurrentBoxLengthInverse, VectorDistance, ScalingDistanceUnitBox )
        ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
        CALL MatrixVectorMultiplication( CurrentBoxLength, ScalingDistanceUnitBox, VectorDistance )
        ! Magnitude of the vector distance (squared)
        SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
        ! Cutoff distance (sphere circumscribing the cylinder)
        IF( SquaredDistance <= cSquaredCutoffSphere ) THEN
          ! Spherocylinder circumscribing the cylinder
          OverlapSPC = .FALSE.
          CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, ContactDistance, ParallelSPC, &
          &                     OverlapSPC, .FALSE. )
          ! Vega-Lago Criterion
          IF( OverlapSPC ) THEN
            ! Remove periodic boundary conditions
            cjPosition(:,jCylinder) = ciPosition(:,iCylinder) + VectorDistance(:)
            ! Check overlap between cylinders
            CALL OverlapCheckCYL( iQuaternion, jQuaternion, iOrientation, jOrientation, VectorDistance, ciPosition(:,iCylinder), &
            &                     cjPosition(:,jCylinder), ParallelSPC, OverlapCYL )
            ! Overlap detected
            IF( OverlapCYL ) THEN
              Overlap = .TRUE.
              RETURN ! Return immediately if an overlap is detected
            END IF
          END IF
        END IF
      END DO
    END DO
  END IF
END DO

RETURN ! No overlaps detected

END SUBROUTINE ListOverlapCheckInitialConfiguration

! *********************************************************************************************** !
!   This subroutine uses lists to check if a random displacement (translation or rotation) of a   !
!                   fixed particle i causes any overlaps with other particles j                   !
! *********************************************************************************************** !
SUBROUTINE ListOverlapCheck( iParticle, iQuaternion, iOrientation, iPosition, ciPosition, ContactDistance, CurrentBoxLength, &
&                            CurrentBoxLengthInverse, Overlap, CellHalfLogical )

! Uses two modules: linked lists and overlap check algorithms
USE LinkedLists, ONLY: Neighbours, CellIndex

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 ) :: iParticle, jParticle ! Counters (particle)
INTEGER( Kind= Int64 ) :: iCylinder, jCylinder ! Counters (cylinder)
INTEGER( Kind= Int64 ) :: jList                ! Counters (neighbour list)

! INTEGER VARIABLES (ARRAY)
INTEGER( Kind= Int64 ), DIMENSION( 3 )          :: iCellIndex     ! 3D cell index
INTEGER( Kind= Int64 ), DIMENSION( nParticles ) :: jNeighbourList ! List of neighbour j particles

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 ) :: SquaredDistance ! Magnitude of the vector distance between particles i and j (squared)
REAL( Kind= Real64 ) :: ContactDistance ! Contact distance (Perram-Wertheim or Vega-Lago methods)

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: VectorDistance             ! Vector distance between particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )    :: iPosition, jPosition       ! Position of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )    :: iOrientation, jOrientation ! Orientation of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )    :: ScalingDistanceUnitBox     ! Scaled position (unit box)
REAL( Kind= Real64 ), DIMENSION( 9 )    :: CurrentBoxLength           ! Box length
REAL( Kind= Real64 ), DIMENSION( 9 )    :: CurrentBoxLengthInverse    ! Box length (inverse)
REAL( Kind= Real64 ), DIMENSION( 0:3 )  :: iQuaternion, jQuaternion   ! Quaternions of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3, 4 ) :: ciPosition, cjPosition     ! Position of particles i and j

! LOGICAL VARIABLES 
LOGICAL :: Overlap         ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: OverlapHER      ! Detects overlap between two ellipsoids of revolution: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: OverlapSPC      ! Detects overlap between two spherocylinders: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: OverlapCYL      ! Detects overlap between two cylinders: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: ParallelSPC     ! Checks the relative orientation of two spherocylinders : TRUE = parallel orientation; FALSE = non-parallel orientation
LOGICAL :: CellHalfLogical ! Checks whether a cell and its 26 surrounding cells are searched (FALSE) or just its 13 neighbour cells (TRUE)

! Initialization
jList       = 0
Overlap     = .FALSE.
OverlapHER  = .FALSE.
OverlapCYL  = .FALSE.
ParallelSPC = .FALSE.

! Spatial transformation
CALL MatrixVectorMultiplication( CurrentBoxLengthInverse, iPosition, ScalingDistanceUnitBox )

! Minimum image convention
ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )

! Index function allocating particle i to cell
iCellIndex = CellIndex( ScalingDistanceUnitBox )

! Debug check
IF( CellHalfLogical .AND. ANY( iCellIndex(:) /= pCellIndex(:,iParticle) ) ) THEN ! Should never happen
  WRITE( *, "(G0,2(G0,', '),2G0,2(G0,', '),2G0)" ) "Cell mismatch: [", iCellIndex(:), "] and [", pCellIndex(:,iParticle), "]"
  STOP
END IF

! List of neighbours
jNeighbourList = Neighbours( iParticle, iCellIndex, CellHalfLogical )

! Loop until no more entries in the neighbour list
DO
  ! Next entry
  jList = jList + 1
  ! Neighbour particle index
  jParticle = jNeighbourList(jList)
  ! List exhausted
  IF( jParticle == 0 ) EXIT
  ! Skip self (should never happen)
  IF( iParticle == jParticle ) CYCLE
  ! Position of particle j
  jPosition(:) = pPositionMC(:,jParticle)
  ! Position of particle j
  cjPosition(:,:) = cPositionMC(:,:,jParticle)
  ! Orientation of particle j
  jOrientation(:) = pOrientationMC(:,jParticle)
  ! Quaternion of particle j
  jQuaternion(:) = pQuaternionMC(:,jParticle)
  ! Vector distance between particles i and j
  VectorDistance(:) = jPosition(:) - iPosition(:)
  ! Minimum image convention
  CALL MatrixVectorMultiplication( CurrentBoxLengthInverse, VectorDistance, ScalingDistanceUnitBox ) ! Spatial transformation
  ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
  CALL MatrixVectorMultiplication( CurrentBoxLength, ScalingDistanceUnitBox, VectorDistance ) ! Spatial transformation
  ! Magnitude of the vector distance (squared)
  SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
  OverlapSPC = .FALSE.
  OverlapCYL = .FALSE.
  IF( SquaredDistance <= SquaredCutoffSphere ) THEN
    ! Spherocylinder circumscribing the non-convex body
    CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, ContactDistance, ParallelSPC, OverlapSPC, &
    &                     .TRUE. )
    ! Vega-Lago criterion
    IF( OverlapSPC ) THEN
      OverlapCYL = .TRUE. ! Check overlap between cylinders
    ELSE
      OverlapCYL = .FALSE.
    END IF
  END IF
  ! Test overlap between cylinders
  IF( OverlapCYL ) THEN
    ! First loop takes one of the four cylinders from particle i
    DO iCylinder = 1, 4
      ! Second loop takes one of the four cylinders from particle j
      DO jCylinder = 1, 4
        ! Vector distance between particles i and j
        VectorDistance(:) = cjPosition(:,jCylinder) - ciPosition(:,iCylinder)
        ! Minimum Image Convention (see Allen and Tildesley, 2nd Edition (2017), pages 35-45)
        CALL MatrixVectorMultiplication( CurrentBoxLengthInverse, VectorDistance, ScalingDistanceUnitBox )
        ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
        CALL MatrixVectorMultiplication( CurrentBoxLength, ScalingDistanceUnitBox, VectorDistance )
        ! Magnitude of the vector distance (squared)
        SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
        ! Cutoff distance (sphere circumscribing the cylinder)
        IF( SquaredDistance <= cSquaredCutoffSphere ) THEN
          ! Spherocylinder circumscribing the cylinder
          OverlapSPC = .FALSE.
          CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, ContactDistance, ParallelSPC, &
          &                     OverlapSPC, .FALSE. )
          ! Vega-Lago Criterion
          IF( OverlapSPC ) THEN
            cjPosition(:,jCylinder) = ciPosition(:,iCylinder) + VectorDistance(:)
            CALL OverlapCheckCYL( iQuaternion, jQuaternion, iOrientation, jOrientation, VectorDistance, ciPosition(:,iCylinder), &
            &                     cjPosition(:,jCylinder), ParallelSPC, OverlapCYL )
            ! Overlap detected
            IF( OverlapCYL ) THEN
              Overlap = .TRUE.
              RETURN ! Return immediately if an overlap is detected
            END IF
          END IF
        END IF
      END DO
    END DO
  END IF
END DO

RETURN ! No overlaps detected

END SUBROUTINE ListOverlapCheck

! *********************************************************************************************** !
!    This subroutine uses lists to check if a random volume scaling (isotropic or anisotropic)    !
!                               causes any overlaps among particles                               !
! *********************************************************************************************** !
SUBROUTINE FullListOverlapCheck( ContactDistance, CurrentBoxLength, CurrentBoxLengthInverse, Overlap, HalfNeighbours )

! Uses two modules: linked lists and overlap check algorithms
USE LinkedLists, ONLY: MakeList

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 )           :: iParticle      ! Counter (particle)
INTEGER( Kind= Int64 ), OPTIONAL :: HalfNeighbours ! Checks whether a cell and its 26 surrounding cells are searched (PRESENT) or just its 13 neighbour cells (NOT PRESENT)

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 )                   :: ContactDistance         ! Contact distance (Perram-Wertheim or Vega-Lago methods)

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 3 )   :: BoxCutoff               ! Box cutoff (x-, y-, and z-directions)
REAL( Kind= Real64 ), DIMENSION( 3 )   :: iPosition               ! Position of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )   :: iOrientation            ! Orientation of particles i and j
REAL( Kind= Real64 ), DIMENSION( 9 )   :: CurrentBoxLength        ! Box length
REAL( Kind= Real64 ), DIMENSION( 9 )   :: CurrentBoxLengthInverse ! Box length (inverse)
REAL( Kind= Real64 ), DIMENSION( 0:3 ) :: iQuaternion             ! Quaternions of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3,4 ) :: ciPosition              ! Position of particles i and j

! LOGICAL VARIABLES
LOGICAL :: Overlap         ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: CellHalfLogical ! Checks whether a cell and its 26 surrounding cells are searched (FALSE) or just its 13 neighbour cells (TRUE)
LOGICAL :: SharedOverlap   ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: PrivateOverlap  ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected

! Initialization
Overlap         = .FALSE.
CellHalfLogical = .FALSE.
SharedOverlap   = .FALSE.
PrivateOverlap  = .FALSE.

! Re-initialize cell list (if necessary)
BoxCutoff(1) = cLargestSphereDiameter / CurrentBoxLength(1)
BoxCutoff(2) = cLargestSphereDiameter / CurrentBoxLength(5)
BoxCutoff(3) = cLargestSphereDiameter / CurrentBoxLength(9)
CALL MakeList( BoxCutoff, pPositionMC, CurrentBoxLengthInverse )

! Check whether the number of cells along any direction (x, y, or z) is less than 3
IF( .NOT. CellListControl ) RETURN

! Consider all 26 neighbours or only 13
CellHalfLogical = PRESENT( HalfNeighbours )

! Loop over all particles (consider half of neighbours)
!##################################################################################################!
!$OMP PARALLEL DO DEFAULT( Shared ) &                                                              !
!$OMP PRIVATE( iParticle, iQuaternion, iOrientation, iPosition, ciPosition, ContactDistance ) &    !
!$OMP PRIVATE( CellHalfLogical ) &                                                                 !
!$OMP REDUCTION( .OR. : PrivateOverlap )                                                           !
!##################################################################################################!
DO iParticle = 1, nParticles
  IF( SharedOverlap ) CYCLE
  ! Quaternion of particle i
  iQuaternion = pQuaternionMC(:,iParticle)
  ! Orientation of particle i
  iOrientation = pOrientationMC(:,iParticle)
  ! Position of particle i
  iPosition = pPositionMC(:,iParticle)
  ! Position of cylinders i
  ciPosition = cPositionMC(:,:,iParticle)
  ! Overlap check between particle i and its neighbours
  CALL ListOverlapCheck( iParticle, iQuaternion, iOrientation, iPosition, ciPosition, ContactDistance, CurrentBoxLength, &
  &                      CurrentBoxLengthInverse, PrivateOverlap, CellHalfLogical )
  ! Overlap found
  IF( PrivateOverlap ) SharedOverlap = .TRUE.
END DO
!#################################################################################################!
!$OMP END PARALLEL DO                                                                             !
!#################################################################################################!

! Returns control to the calling program unit if an overlap is detected
IF( SharedOverlap ) THEN
  Overlap = .TRUE.
  RETURN
END IF

RETURN ! No overlaps detected

END SUBROUTINE FullListOverlapCheck

END MODULE OverlapCheckLists