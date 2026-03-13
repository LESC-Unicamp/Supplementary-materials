MODULE OverlapCheckSystem

! Uses three modules: global variables, overlap algorithms, and vector operations
USE GlobalVariables
USE OverlapAlgorithms, ONLY: OverlapCheckSPC, OverlapCheckCYL
USE VectorOperations

CONTAINS

! *********************************************************************************************** !
! This subroutine checks if a random displacement (translation or rotation) of a fixed particle i !
!                           causes any overlaps with other particles j                            !
! *********************************************************************************************** !
SUBROUTINE ParticleOverlapCheck( iParticle, iQuaternion, iOrientation, iPosition, ciPosition, CurrentBoxLength, &
&                                CurrentBoxLengthInverse, Overlap )

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 ) :: iParticle, jParticle ! Counters (particles)
INTEGER( Kind= Int64 ) :: iCylinder, jCylinder ! Counters (cylinders)

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 ) :: ContactDistance ! Contact distance (Vega-Lago method)
REAL( Kind= Real64 ) :: SquaredDistance ! Magnitude of the vector distance between particles i and j (squared)

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: ScalingDistanceUnitBox     ! Scaled position (unit box)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: iPosition, jPosition       ! Position of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )    :: iOrientation, jOrientation ! Orientation of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )    :: VectorDistance             ! Vector distance between particles i and j
REAL( Kind= Real64 ), DIMENSION( 9 )    :: CurrentBoxLength           ! Box length
REAL( Kind= Real64 ), DIMENSION( 9 )    :: CurrentBoxLengthInverse    ! Box length (inverse)
REAL( Kind= Real64 ), DIMENSION( 0:3 )  :: iQuaternion, jQuaternion   ! Quaternion of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3, 4 ) :: ciPosition, cjPosition     ! Position of cylinders of particles i and j

! LOGICAL VARIABLES
LOGICAL :: Overlap        ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: SharedOverlap  ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: PrivateOverlap ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: OverlapSPC     ! Detects overlap between two spherocylinders: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: OverlapCYL     ! Detects overlap between two cylinders: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: ParallelSPC    ! Checks the relative orientation of two spherocylinders : TRUE = parallel orientation; FALSE = non-parallel orientation

! Initialization of overlap detector
PrivateOverlap = .FALSE.
SharedOverlap  = .FALSE.
OverlapCYL     = .FALSE.
Overlap        = .FALSE.

!###############################################################################################!
!$OMP PARALLEL DO DEFAULT( Shared ) &                                                           !
!$OMP PRIVATE( jParticle, jPosition, jOrientation, jQuaternion, cjPosition, VectorDistance ) &  !
!$OMP PRIVATE( ScalingDistanceUnitBox, SquaredDistance, OverlapSPC, OverlapCYL ) &              !
!$OMP PRIVATE( ContactDistance, ParallelSPC, iCylinder, jCylinder ) &                           !
!$OMP REDUCTION( .OR. : PrivateOverlap )                                                        !
!###############################################################################################!
DO jParticle = 1, nParticles
  ! Cycle if a thread founds an overlap
  IF( SharedOverlap ) CYCLE
  ! Cycle if a thread founds an overlap
  IF( jParticle == iParticle ) CYCLE
  ! Configuration of particle j
  jPosition(:)    = pPositionMC(:,jParticle)
  jOrientation(:) = pOrientationMC(:,jParticle)
  jQuaternion(:)  = pQuaternionMC(:,jParticle)
  cjPosition(:,:) = cPositionMC(:,:,jParticle)
  ! Vector distance between particles i and j
  VectorDistance(:) = jPosition(:) - iPosition(:)
  ! Minimum image convention
  CALL MatrixVectorMultiplication( CurrentBoxLengthInverse, VectorDistance, ScalingDistanceUnitBox )
  ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
  CALL MatrixVectorMultiplication( CurrentBoxLength, ScalingDistanceUnitBox, VectorDistance )
  ! Magnitude of the vector distance (squared)
  SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
  ! Initialization
  OverlapSPC = .FALSE.
  OverlapCYL = .FALSE.
  ! Cutoff distance (sphere circumscribing the spherocylinder enclosing the nonconvex body)
  IF( SquaredDistance <= SquaredCutoffSphere ) THEN
    ! Spherocylinder circumscribing the nonconvex body
    CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, ContactDistance, ParallelSPC, OverlapSPC, &
    &                     .TRUE. )
    ! Vega-Lago criterion
    IF( OverlapSPC ) THEN
      OverlapCYL = .TRUE.
    ELSE
      OverlapCYL = .FALSE.
    END IF
  END IF
  ! Considering cylinders (only if preliminary tests fail)
  IF( OverlapCYL ) THEN
    ! First loop takes one of the four cylinders from particle i
    DO iCylinder = 1, 4
      ! Second loop takes one of the four cylinders from particle j
      DO jCylinder = 1, 4
        ! Cycle if a thread founds an overlap
        IF( SharedOverlap ) CYCLE
        ! Vector distance between particles i and j
        VectorDistance(:) = cjPosition(:,jCylinder) - ciPosition(:,iCylinder)
        ! Minimum image convention
        CALL MatrixVectorMultiplication( CurrentBoxLengthInverse, VectorDistance, ScalingDistanceUnitBox )
        ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
        CALL MatrixVectorMultiplication( CurrentBoxLength, ScalingDistanceUnitBox, VectorDistance )
        ! Magnitude of the vector distance (squared)
        SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
        ! Cutoff distance (sphere circumscribing the spherocylinder enclosing the cylinder)
        IF( SquaredDistance <= cSquaredCutoffSphere ) THEN
          ! Spherocylinder circumscribing the cylinder
          CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, ContactDistance, ParallelSPC, &
          &                     OverlapSPC, .FALSE. )
          ! Vega-Lago Criterion
          IF( OverlapSPC ) THEN
            ! Apply periodic boundary conditions on the position of particle j
            cjPosition(:,jCylinder) = ciPosition(:,iCylinder) + VectorDistance(:)
            ! Overlap test for cylinders
            CALL OverlapCheckCYL( iQuaternion, jQuaternion, iOrientation, jOrientation, VectorDistance, ciPosition(:,iCylinder), &
            &                     cjPosition(:,jCylinder), ParallelSPC, PrivateOverlap )
            ! Overlap found
            IF( PrivateOverlap ) SharedOverlap = .TRUE.
          END IF
        END IF
      END DO
    END DO
  END IF
END DO
!###############################################################################################!
!$OMP END PARALLEL DO                                                                           !
!###############################################################################################!

! Returns control to the calling program unit if an overlap is detected
IF( SharedOverlap ) THEN
  Overlap = .TRUE.
  RETURN
END IF

RETURN ! No overlaps detected

END SUBROUTINE ParticleOverlapCheck

! *********************************************************************************************** !
!          This subroutine checks if a random volume scaling (isotropic or anisotropic)           !
!                               causes any overlaps among particles                               !
! *********************************************************************************************** !
SUBROUTINE FullOverlapCheck( ContactDistance, CurrentBoxLength, CurrentBoxLengthInverse, Overlap )

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 ) :: iParticle, jParticle ! Counters (particles)
INTEGER( Kind= Int64 ) :: iCylinder, jCylinder ! Counters (cylinders)

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 ) :: ContactDistance ! Contact distance (Perram-Wertheim or Vega-Lago methods)
REAL( Kind= Real64 ) :: SquaredDistance ! Magnitude of the vector distance between particles i and j (squared)

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 9 )   :: CurrentBoxLength           ! Box length
REAL( Kind= Real64 ), DIMENSION( 9 )   :: CurrentBoxLengthInverse    ! Box length (inverse)
REAL( Kind= Real64 ), DIMENSION( 3 )   :: ScalingDistanceUnitBox     ! Scaled position (unit box)
REAL( Kind= Real64 ), DIMENSION( 3 )   :: iPosition, jPosition       ! Position of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )   :: iOrientation, jOrientation ! Orientation of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )   :: VectorDistance             ! Vector distance between particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )   :: ciPosition, cjPosition     ! Position of cylinders of particles i and j
REAL( Kind= Real64 ), DIMENSION( 0:3 ) :: iQuaternion, jQuaternion   ! Quaternion of particles i and j

! LOGICAL VARIABLES
LOGICAL :: Overlap ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: SharedOverlap  ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: PrivateOverlap ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: OverlapSPC     ! Detects overlap between two spherocylinders: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: OverlapCYL     ! Detects overlap between two cylinders: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: ParallelSPC    ! Checks the relative orientation of two spherocylinders : TRUE = parallel orientation; FALSE = non-parallel orientation

! Initialization of overlap detector
Overlap        = .FALSE.
OverlapSPC     = .FALSE.
ParallelSPC    = .FALSE.
SharedOverlap  = .FALSE.
PrivateOverlap = .FALSE.

!###############################################################################################!
!$OMP PARALLEL DO COLLAPSE( 2 ) DEFAULT( Shared ) &                                             !
!$OMP PRIVATE( iParticle, jParticle, iPosition, jPosition, iOrientation, jOrientation ) &       !
!$OMP PRIVATE( iQuaternion, jQuaternion, ciPosition, cjPosition, ScalingDistanceUnitBox ) &     !
!$OMP PRIVATE( VectorDistance, SquaredDistance, OverlapSPC, OverlapCYL, ContactDistance ) &     !
!$OMP PRIVATE( ParallelSPC, iCylinder, jCylinder ) &                                            !
!$OMP REDUCTION( .OR. : PrivateOverlap )                                                        !
!###############################################################################################!
DO iParticle = 1, nParticles
  DO jParticle = 1, nParticles
    ! Cycle if a thread founds an overlap
    IF( SharedOverlap ) CYCLE
    ! Cycle if index j <= i (prevents double counting)
    IF( jParticle <= iParticle ) CYCLE
    ! Initialization
    OverlapCYL = .FALSE.
    ! Position of particles i and j
    iPosition(:) = pPositionMC(:,iParticle)
    jPosition(:) = pPositionMC(:,jParticle)
    ! Orientation of particles i and j
    iOrientation(:)  = pOrientationMC(:,iParticle)
    jOrientation(:)  = pOrientationMC(:,jParticle)
    ! Quaternion of particles i and j
    iQuaternion(:)  = pQuaternionMC(:,iParticle)
    jQuaternion(:)  = pQuaternionMC(:,jParticle)
    ! Vector distance between particles i and j
    VectorDistance(:) = jPosition(:) - iPosition(:)
    ! Minimum image convention
    CALL MatrixVectorMultiplication( CurrentBoxLengthInverse, VectorDistance, ScalingDistanceUnitBox )
    ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
    CALL MatrixVectorMultiplication( CurrentBoxLength, ScalingDistanceUnitBox, VectorDistance )
    ! Magnitude of the vector distance (squared)
    SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
    ! Cutoff distance (sphere circumscribing the spherocylinder enclosing the nonconvex body)
    IF( SquaredDistance <= SquaredCutoffSphere ) THEN
      ! Spherocylinder circumscribing the nonconvex body
      CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, ContactDistance, ParallelSPC, OverlapSPC, &
      &                     .TRUE. )
      ! Vega-Lago Criterion
      IF( OverlapSPC ) THEN
        OverlapCYL = .TRUE. ! Check overlap between cylinders
      ELSE
        OverlapCYL = .FALSE.
      END IF
    END IF
    ! Considering cylinders (only if preliminary tests fail)
    IF( OverlapCYL ) THEN
      ! First loop takes one of the four cylinders from particle i
      DO iCylinder = 1, 4
        ! Position of cylinder of particle i
        ciPosition(:) = cPositionMC(:,iCylinder,iParticle)
        ! Second loop takes one of the four cylinders from particle j
        DO jCylinder = 1, 4
          ! Cycle if a thread founds an overlap
          IF( SharedOverlap ) CYCLE
          ! Position of cylinder of particle j
          cjPosition(:) = cPositionMC(:,jCylinder,jParticle)
          ! Vector distance between cylinders of particles i and j
          VectorDistance(:) = cjPosition(:) - ciPosition(:)
          ! Minimum image convention
          CALL MatrixVectorMultiplication( CurrentBoxLengthInverse, VectorDistance, ScalingDistanceUnitBox )
          ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
          CALL MatrixVectorMultiplication( CurrentBoxLength, ScalingDistanceUnitBox, VectorDistance )
          ! Magnitude of the vector distance (squared)
          SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
          ! Cutoff distance (sphere circumscribing the spherocylinder enclosing the cylinder)
          IF( SquaredDistance <= cSquaredCutoffSphere ) THEN
            ! Cutoff distance (spherocylinder circumscribing the cylinder)
            CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, ContactDistance, ParallelSPC, &
            &                     OverlapSPC, .FALSE. )
            ! Vega-Lago Criterion
            IF( OverlapSPC ) THEN
              ! Apply periodic boundary conditions on the position of particle j
              cjPosition(1) = ciPosition(1) + VectorDistance(1)
              cjPosition(2) = ciPosition(2) + VectorDistance(2)
              cjPosition(3) = ciPosition(3) + VectorDistance(3)
              ! Overlap test for cylinders
              CALL OverlapCheckCYL( iQuaternion, jQuaternion, iOrientation, jOrientation, VectorDistance, ciPosition, cjPosition, &
              &                     ParallelSPC, PrivateOverlap )
              ! Overlap found
              IF( PrivateOverlap ) SharedOverlap = .TRUE.
            END IF
          END IF
        END DO
      END DO
    END IF
  END DO ! Loop of particle j
END DO ! Loop of particle i
!#############################################################################################!
!$OMP END PARALLEL DO                                                                         !
!#############################################################################################!

! Returns control to the calling program unit if an overlap is detected
IF( SharedOverlap ) THEN
  Overlap = .TRUE.
  RETURN
END IF

RETURN

END SUBROUTINE FullOverlapCheck

END MODULE OverlapCheckSystem