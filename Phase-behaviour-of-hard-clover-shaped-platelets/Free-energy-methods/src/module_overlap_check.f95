! ############################################################################################### !
!                          FREE-ENERGY METHODS FOR NON-CONVEX MOLECULES                           !
!               This code contains overlap check routines used in the main program.               !
!                                                                                                 !
! Version number: 1.0.0                                                                           !
! ############################################################################################### !
!                                University of Campinas (Unicamp)                                 !
!                                 School of Chemical Engineering                                  !
!                                     Nathan Barros de Souza                                      !
!                             --------------------------------------                              !
!                             Supervisor: Luís Fernando Mercier Franco                            !
!                             --------------------------------------                              !
!                                        August 13th, 2026                                        !
! ############################################################################################### !
! Main References:                  M. P. Allen, D. J. Tildesley                                  !
!                           Oxford University Press, 2nd Edition (2017)                           !
!                             DOI: 10.1093/oso/9780198803195.001.0001                             !
!                             --------------------------------------                              !
!                                          G. Marsaglia                                           !
!                            Ann. Math. Statist. 43(2), 645-646 (1972)                            !
!                                  DOI: 10.1214/aoms/1177692644                                   !
!                             ---------------------------------------                             !
!                     J. Graaf, L. Filion, M. Marechal, R. Roij, M. Dijkstra                      !
!                                J. Chem. Phys. 137, 214101 (2012)                                !
!                                     DOI: 10.1063/1.4767529                                      !
! ############################################################################################### !
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !
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
SUBROUTINE ParticleCheckOverlap( iParticle, iQuaternion, iOrientation, iPosition, ciPosition, Overlap )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: iParticle, jParticle ! Counter (particles)
INTEGER( Kind= Int64 ) :: iCylinder, jCylinder ! Counter (cylinders)

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                    :: ContactDistance            ! Vega-Lago Contact Distance (variable)
REAL( Kind= Real64 )                    :: SquaredDistance            ! Vector distance between particles i and j (squared)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: ScalingDistanceUnitBox     ! Scaled coordinates (unit box)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: iPosition, jPosition       ! Position of particle i and j
REAL( Kind= Real64 ), DIMENSION( 3 )    :: iOrientation, jOrientation ! Orientation of particle i
REAL( Kind= Real64 ), DIMENSION( 3 )    :: VectorDistance             ! Vector distance between particles i and j
REAL( Kind= Real64 ), DIMENSION( 0: 3 ) :: iQuaternion, jQuaternion   ! Quaternion of particle i
REAL( Kind= Real64 ), DIMENSION( 3, 4 ) :: ciPosition, cjPosition     ! Position of cylinders of particles i and j

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: Overlap        ! Detects overlap between two particles (Vega-Lago) : TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: ParallelSPC    ! Checks the relative orientation of two spherocylinders : TRUE = parallel orientation; FALSE = non-parallel orientation
LOGICAL :: OverlapCYL     ! Detects overlap between two cylinders : TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: SharedOverlap  ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: PrivateOverlap ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected

! Initialization of overlap detector
PrivateOverlap = .FALSE.
SharedOverlap = .FALSE.
OverlapCYL = .FALSE.
Overlap = .FALSE.

! *********************************************************************************************** !
! Particle loops                                                                                  !
! *********************************************************************************************** !
!  The loops below analyze the neighboring particles j around a fixed particle i                  !
!  and detect possible overlaps between them.                                                     !
! *********************************************************************************************** !

!#################################################################################################!
!$OMP PARALLEL DO DEFAULT( Shared ) &                                                             !
!$OMP PRIVATE( jParticle, jPosition, jOrientation, jQuaternion, cjPosition, VectorDistance ) &    !
!$OMP PRIVATE( ScalingDistanceUnitBox, SquaredDistance, ContactDistance, OverlapCYL ) &           !
!$OMP PRIVATE( ParallelSPC, iCylinder, jCylinder ) &                                              !
!$OMP REDUCTION( .OR. : PrivateOverlap )                                                          !
!#################################################################################################!
! First loop takes only particles whose j-indexes are below the i-index of the fixed particle
DO jParticle = 1, nParticles
  ! Cycle if a thread founds an overlap
  IF( SharedOverlap ) CYCLE
  ! Cycle if a thread founds an overlap
  IF( jParticle == iParticle ) CYCLE
  ! Considering initially the center of mass (preliminary tests)
  jPosition(:)    = pPositionMC(:,jParticle)
  jQuaternion(:)  = pQuaternionMC(:,jParticle)
  jOrientation(:) = pOrientationMC(:,jParticle)
  cjPosition(:,:) = cPositionMC(:,:,jParticle)
  ! Vector distance between particles i and j
  VectorDistance(:) = jPosition(:) - iPosition(:)
  ! Minimum Image Convention (see Allen and Tildesley, 2nd Edition (2017), pages 35-45)
  CALL MatrixVectorMultiplication( BoxLengthInverse, VectorDistance, ScalingDistanceUnitBox )
  ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
  CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, VectorDistance )
  ! Magnitude of the vector distance (squared)
  SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
  ! Initialization
  OverlapCYL = .FALSE.
  ! Cutoff distance (sphere circumscribing the non-convex body)
  IF( SquaredDistance <= SquaredCutoffSphere ) THEN
    ! Spherocylinder circumscribing the non-convex body
    CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, ContactDistance, ParallelSPC )
    ! Vega-Lago criterion
    IF( ContactDistance <= SquaredCutoffSpherocylinder ) THEN
      OverlapCYL = .TRUE. ! Check overlap between cylinders
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
        ! Minimum Image Convention (see Allen and Tildesley, 2nd Edition (2017), pages 35-45)
        CALL MatrixVectorMultiplication( BoxLengthInverse, VectorDistance, ScalingDistanceUnitBox )
        ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
        CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, VectorDistance )
        ! Magnitude of the vector distance (squared)
        SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
        ! Cutoff distance (sphere circumscribing the cylinder)
        IF( SquaredDistance <= cSquaredCutoffSphere ) THEN
          ! Spherocylinder circumscribing the cylinder
          CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, ContactDistance, ParallelSPC )
          ! Vega-Lago Criterion
          IF( ContactDistance <= ( cDiameter * cDiameter ) ) THEN
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
!#################################################################################################!
!$OMP END PARALLEL DO                                                                             !
!#################################################################################################!

! Returns control to the calling program unit if an overlap is detected
IF( SharedOverlap ) THEN
  Overlap = .TRUE.
  RETURN
END IF

RETURN ! No overlaps detected

END SUBROUTINE ParticleCheckOverlap

! *********************************************************************************************** !
!              This subroutine checks for overlaps among the particles in the system              !
! *********************************************************************************************** !
SUBROUTINE FullOverlapCheck( Overlap )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: iParticle, jParticle ! Counter (particles)
INTEGER( Kind= Int64 ) :: iCylinder, jCylinder ! Counter (cylinders)

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                    :: ContactDistance            ! Vega-Lago Contact Distance (variable)
REAL( Kind= Real64 )                    :: SquaredDistance            ! Vector distance between particles i and j (squared)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: ScalingDistanceUnitBox     ! Periodic boundary factor
REAL( Kind= Real64 ), DIMENSION( 3 )    :: iPosition, jPosition       ! Position of particle i and j
REAL( Kind= Real64 ), DIMENSION( 3 )    :: iOrientation, jOrientation ! Orientation of particle i
REAL( Kind= Real64 ), DIMENSION( 3 )    :: ciPosition, cjPosition     ! Position of cylinders of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )    :: VectorDistance             ! Vector distance between particles i and j
REAL( Kind= Real64 ), DIMENSION( 0: 3 ) :: iQuaternion, jQuaternion   ! Quaternion of particle i

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: Overlap        ! Detects overlap between two particles (Vega-Lago) : TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: ParallelSPC    ! Checks the relative orientation of two spherocylinders : TRUE = parallel orientation; FALSE = non-parallel orientation
LOGICAL :: OverlapCYL     ! Detects overlap between two cylinders : TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: SharedOverlap  ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: PrivateOverlap ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected

! Initialization of overlap detector
Overlap        = .FALSE.
OverlapCYL     = .FALSE.
SharedOverlap  = .FALSE.
PrivateOverlap = .FALSE.

! *********************************************************************************************** !
! Particle loops                                                                                  !
! *********************************************************************************************** !
!  The loops below analyze the neighboring particles j around a fixed particle i                  !
!  and detect possible overlaps between them.                                                     !
! *********************************************************************************************** !

!#################################################################################################!
!$OMP PARALLEL DO COLLAPSE( 2 ) DEFAULT( Shared ) &                                               !
!$OMP PRIVATE( iParticle, jParticle, iPosition, jPosition, iOrientation, jOrientation ) &         !
!$OMP PRIVATE( iQuaternion, jQuaternion, ciPosition, cjPosition, ScalingDistanceUnitBox ) &       !
!$OMP PRIVATE( VectorDistance, SquaredDistance, OverlapCYL, ContactDistance ) &                   !
!$OMP PRIVATE( ParallelSPC, iCylinder, jCylinder ) &                                              !
!$OMP REDUCTION( .OR. : PrivateOverlap )                                                          !
!#################################################################################################!
! First loop represents a particle with a fixed index i
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
    ! Minimum Image Convention (see Allen and Tildesley, 2nd Edition (2017), pages 35-45)
    CALL MatrixVectorMultiplication( BoxLengthInverse, VectorDistance, ScalingDistanceUnitBox )
    ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
    CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, VectorDistance )
    ! Magnitude of the vector distance (squared)
    SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
    ! Cutoff distance (sphere circumscribing the non-convex body)
    IF( SquaredDistance <= SquaredCutoffSphere ) THEN
      ! Spherocylinder circumscribing the non-convex body
      CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, ContactDistance, ParallelSPC )
      ! Vega-Lago Criterion
      IF( ContactDistance <= SquaredCutoffSpherocylinder ) THEN
        OverlapCYL = .TRUE. ! Check overlap between cylinders
      ELSE
        OverlapCYL = .FALSE.
      END IF
    END IF
    ! Considering the cylinders (only if preliminary tests fail)
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
          ! Minimum Image Convention (see Allen and Tildesley, 2nd Edition (2017), pages 35-45)
          CALL MatrixVectorMultiplication( BoxLengthInverse, VectorDistance, ScalingDistanceUnitBox )
          ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
          CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, VectorDistance )
          ! Magnitude of the vector distance (squared)
          SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
          ! Cutoff distance (sphere circumscribing a spherocylinder circumscribing the cylinder)
          IF( SquaredDistance <= cSquaredCutoffSphere ) THEN
            ! Cutoff distance (spherocylinder circumscribing the cylinder)
            CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, ContactDistance, ParallelSPC )
            ! Vega-Lago Criterion
            IF( ContactDistance <= ( cDiameter * cDiameter ) ) THEN
              ! Apply periodic boundary conditions on the position of particle j
              cjPosition(1) = ciPosition(1) + VectorDistance(1)
              cjPosition(2) = ciPosition(2) + VectorDistance(2)
              cjPosition(3) = ciPosition(3) + VectorDistance(3)
              ! Time-consuming overlap check
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
!#################################################################################################!
!$OMP END PARALLEL DO                                                                             !
!#################################################################################################!

! Returns control to the calling program unit if an overlap is detected
IF( SharedOverlap ) THEN
  Overlap = .TRUE.
  RETURN
END IF

RETURN ! No overlaps

END SUBROUTINE FullOverlapCheck

END MODULE OverlapCheckSystem