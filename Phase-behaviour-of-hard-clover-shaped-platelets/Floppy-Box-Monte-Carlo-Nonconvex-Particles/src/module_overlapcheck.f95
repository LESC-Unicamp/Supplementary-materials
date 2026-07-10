! ############################################################################################### !
!                             FLOPPY-BOX MONTE CARLO ALGORITHM (FBMC)                             !
!                  This module contains all subroutines used in the main program                  !
!                            to search for overlapping configurations.                            !
!                                                                                                 !
! Version number: 1.0.0                                                                           !
! ############################################################################################### !
!                                University of Campinas (Unicamp)                                 !
!                                 School of Chemical Engineering                                  !
!                                     Nathan Barros de Souza                                      !
!                             --------------------------------------                              !
!                             Supervisor: Luís Fernando Mercier Franco                            !
!                             --------------------------------------                              !
!                                       February 15th, 2024                                       !
! ############################################################################################### !
! Main References:                  M. P. Allen, D. J. Tildesley                                  !
!                           Oxford University Press, 2nd Edition (2017)                           !
!                             DOI: 10.1093/oso/9780198803195.001.0001                             !
! ############################################################################################### !
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !
MODULE OverlapCheck

! Uses one module: global variables
USE GlobalVar

IMPLICIT NONE

CONTAINS

! *********************************************************************************************** !
!                    This subroutine checks there is any overlapping particles                    !
!                                  in the initial configuration                                   !
! *********************************************************************************************** !
SUBROUTINE OverlapCheckInitialConfiguration(  )

! Uses one module: overlap check algorithms
USE OverlapCheckAlgorithms, ONLY: OverlapCheckSPC, OverlapCheckCYL

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: iParticle, jParticle   ! Counters (particle)
INTEGER( Kind= Int64 ) :: iCylinder, jCylinder   ! Counters (cylinders)
INTEGER( Kind= Int64 ) :: pImage                 ! Counters (image)

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                    :: SquaredDistance            ! Magnitude of the vector distance between particles i and j (squared)
REAL( Kind= Real64 )                    :: ContactDistance            ! Contact distance (Perram-Wertheim or Vega-Lago methods)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: VectorDistance             ! Vector distance between particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )    :: iPosition, jPosition       ! Position of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )    :: iOrientation, jOrientation ! Orientation of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )    :: ciPosition, cjPosition     ! Position of cylinders of particles i and j
REAL( Kind= Real64 ), DIMENSION( 0:3 )  :: iQuaternion, jQuaternion   ! Quaternions of particles i and j

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: Overlap     ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: OverlapSPC  ! Detects overlap between two spherocylinders: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: ParallelSPC ! Checks the relative orientation of two spherocylinders : TRUE = parallel orientation; FALSE = non-parallel orientation

! Initialization
Overlap     = .FALSE.
OverlapSPC  = .FALSE.
ParallelSPC = .FALSE.

! First loop represents all particles with indexes i of component Ci
DO iParticle = 1, nParticles - 1
  ! Second loop represents all particles with indexes j of component Cj
  DO jParticle = iParticle + 1, nParticles
    ! Position of particle i
    iPosition(1) = pPosition(1,iParticle)
    iPosition(2) = pPosition(2,iParticle)
    iPosition(3) = pPosition(3,iParticle)
    ! Position of particle j
    jPosition(1) = pPosition(1,jParticle)
    jPosition(2) = pPosition(2,jParticle)
    jPosition(3) = pPosition(3,jParticle)
    ! Orientation of particle i
    iOrientation(1) = pOrientation(1,iParticle)
    iOrientation(2) = pOrientation(2,iParticle)
    iOrientation(3) = pOrientation(3,iParticle)
    ! Orientation of particle j
    jOrientation(1) = pOrientation(1,jParticle)
    jOrientation(2) = pOrientation(2,jParticle)
    jOrientation(3) = pOrientation(3,jParticle)
    ! Quaternion of particle i
    iQuaternion(0) = pQuaternion(0,iParticle)
    iQuaternion(1) = pQuaternion(1,iParticle)
    iQuaternion(2) = pQuaternion(2,iParticle)
    iQuaternion(3) = pQuaternion(3,iParticle)
    ! Quaternion of particle j
    jQuaternion(0) = pQuaternion(0,jParticle)
    jQuaternion(1) = pQuaternion(1,jParticle)
    jQuaternion(2) = pQuaternion(2,jParticle)
    jQuaternion(3) = pQuaternion(3,jParticle)
    ! Vector distance between particles i and j
    VectorDistance(1) = jPosition(1) - iPosition(1)
    VectorDistance(2) = jPosition(2) - iPosition(2)
    VectorDistance(3) = jPosition(3) - iPosition(3)
    ! Magnitude of the vector distance (squared)
    SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
    ! Preliminary test (circumscribing spheres)
    IF( SquaredDistance <= pSquaredSphereDistance ) THEN
      ! Initialization
      OverlapSPC = .FALSE.
      ! Preliminary test (circumscribing spherocylinders)
      CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, ContactDistance, ParallelSPC )
      ! Overlap criterion
      IF( ContactDistance <= pSquaredSpheroCylinderDistance ) THEN
        ! Initialization
        OverlapSPC = .TRUE.
        ! First loop takes one of the four cylinders from particle i
        DO iCylinder = 1, 4
          ! Position of cylinder of particle i
          ciPosition(1) = cPosition(1,iCylinder,iParticle)
          ciPosition(2) = cPosition(2,iCylinder,iParticle)
          ciPosition(3) = cPosition(3,iCylinder,iParticle)
          ! Second loop takes one of the four cylinders from particle j
          DO jCylinder = 1, 4
            ! Position of cylinder of particle j
            cjPosition(1) = cPosition(1,jCylinder,jParticle)
            cjPosition(2) = cPosition(2,jCylinder,jParticle)
            cjPosition(3) = cPosition(3,jCylinder,jParticle)
            ! Vector distance between cylinders of particles i and j
            VectorDistance(1) = cjPosition(1) - ciPosition(1)
            VectorDistance(2) = cjPosition(2) - ciPosition(2)
            VectorDistance(3) = cjPosition(3) - ciPosition(3)
            ! Magnitude of the vector distance (squared)
            SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
            ! Cutoff distance (sphere circumscribing a spherocylinder circumscribing the cylinder)
            IF( SquaredDistance <= cSquaredSphereDistance ) THEN
              ! Cutoff distance (spherocylinder circumscribing the cylinder)
              CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, ContactDistance, ParallelSPC )
              ! Vega-Lago Criterion
              IF( ContactDistance <= cDiameter * cDiameter ) THEN
                ! Apply periodic boundary conditions on the position of particle j
                cjPosition(1) = ciPosition(1) + VectorDistance(1)
                cjPosition(2) = ciPosition(2) + VectorDistance(2)
                cjPosition(3) = ciPosition(3) + VectorDistance(3)
                ! Time-consuming overlap check
                CALL OverlapCheckCYL( iQuaternion, jQuaternion, iOrientation, jOrientation, VectorDistance, ciPosition, &
                &                     cjPosition, ParallelSPC, Overlap )
                ! Overlap criterion
                IF( Overlap ) THEN
                  ! Overlap detected
                  WRITE( *, "(5G0)" ) "Overlap found in the initial configuration between particles ", iParticle, " and ", &
                  &                   jParticle, ". Exiting..."
                  ! Stop program if overlap is detected
                  CALL Sleep( 1 )
                  CALL Exit(  )
                END IF
              END IF
            END IF
          END DO
        END DO
      END IF
    END IF
  END DO
END DO

! First loop represents all particles with indexes i of component Ci
DO iParticle = 1, nParticles
  ! Second loop represents all particles with indexes j of component Cj
  DO jParticle = 1, nParticles
    ! Loop over images
    DO pImage = 1, nImages
      ! Position of particle i
      iPosition(1) = pPosition(1,iParticle)
      iPosition(2) = pPosition(2,iParticle)
      iPosition(3) = pPosition(3,iParticle)
      ! Position of images of particle j
      jPosition(1) = imPosition(1,pImage,jParticle)
      jPosition(2) = imPosition(2,pImage,jParticle)
      jPosition(3) = imPosition(3,pImage,jParticle)
      ! Orientation of particle i
      iOrientation(1) = pOrientation(1,iParticle)
      iOrientation(2) = pOrientation(2,iParticle)
      iOrientation(3) = pOrientation(3,iParticle)
      ! Orientation of images of particle j
      jOrientation(1) = imOrientation(1,pImage,jParticle)
      jOrientation(2) = imOrientation(2,pImage,jParticle)
      jOrientation(3) = imOrientation(3,pImage,jParticle)
      ! Quaternion of particle i
      iQuaternion(0) = pQuaternion(0,iParticle)
      iQuaternion(1) = pQuaternion(1,iParticle)
      iQuaternion(2) = pQuaternion(2,iParticle)
      iQuaternion(3) = pQuaternion(3,iParticle)
      ! Quaternion of particle j
      jQuaternion(0) = imQuaternion(0,pImage,jParticle)
      jQuaternion(1) = imQuaternion(1,pImage,jParticle)
      jQuaternion(2) = imQuaternion(2,pImage,jParticle)
      jQuaternion(3) = imQuaternion(3,pImage,jParticle)
      ! Vector distance between particles i and j
      VectorDistance(1) = jPosition(1) - iPosition(1)
      VectorDistance(2) = jPosition(2) - iPosition(2)
      VectorDistance(3) = jPosition(3) - iPosition(3)
      ! Magnitude of the vector distance (squared)
      SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
      ! Preliminary test (circumscribing spheres)
      IF( SquaredDistance <= pSquaredSphereDistance ) THEN
        ! Initialization
        OverlapSPC = .FALSE.
        ! Preliminary test (circumscribing spherocylinders)
        CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, ContactDistance, ParallelSPC )
        ! Overlap criterion
        IF( ContactDistance <= pSquaredSpheroCylinderDistance ) THEN
          ! Initialization
          OverlapSPC = .TRUE.
          ! First loop takes one of the four cylinders from particle i
          DO iCylinder = 1, 4
            ! Position of cylinder of particle i
            ciPosition(1) = cPosition(1,iCylinder,iParticle)
            ciPosition(2) = cPosition(2,iCylinder,iParticle)
            ciPosition(3) = cPosition(3,iCylinder,iParticle)
            ! Second loop takes one of the four cylinders from particle j
            DO jCylinder = 1, 4
              ! Position of cylinder of particle j
              cjPosition(1) = imcPosition(1,jCylinder,pImage,jParticle)
              cjPosition(2) = imcPosition(2,jCylinder,pImage,jParticle)
              cjPosition(3) = imcPosition(3,jCylinder,pImage,jParticle)
              ! Vector distance between cylinders of particles i and j
              VectorDistance(1) = cjPosition(1) - ciPosition(1)
              VectorDistance(2) = cjPosition(2) - ciPosition(2)
              VectorDistance(3) = cjPosition(3) - ciPosition(3)
              ! Magnitude of the vector distance (squared)
              SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
              ! Cutoff distance (sphere circumscribing a spherocylinder circumscribing the cylinder)
              IF( SquaredDistance <= cSquaredSphereDistance ) THEN
                ! Cutoff distance (spherocylinder circumscribing the cylinder)
                CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, ContactDistance, ParallelSPC )
                ! Vega-Lago Criterion
                IF( ContactDistance <= cDiameter * cDiameter ) THEN
                  ! Apply periodic boundary conditions on the position of particle j
                  cjPosition(1) = ciPosition(1) + VectorDistance(1)
                  cjPosition(2) = ciPosition(2) + VectorDistance(2)
                  cjPosition(3) = ciPosition(3) + VectorDistance(3)
                  ! Time-consuming overlap check
                  CALL OverlapCheckCYL( iQuaternion, jQuaternion, iOrientation, jOrientation, VectorDistance, ciPosition, &
                  &                     cjPosition, ParallelSPC, Overlap )
                  ! Overlap criterion
                  IF( Overlap ) THEN
                    ! Overlap detected
                    WRITE( *, "(5G0)" ) "Overlap found in the initial configuration between particles ", iParticle, " and ", &
                    &                   jParticle, ". Exiting..."
                    ! Stop program if overlap is detected
                    CALL Sleep( 1 )
                    CALL Exit(  )
                  END IF
                END IF
              END IF
            END DO
          END DO
        END IF
      END IF
    END DO
  END DO
END DO

! No overlaps
RETURN

END SUBROUTINE OverlapCheckInitialConfiguration

! *********************************************************************************************** !
!          This subroutine checks if a random volume scaling (isotropic or anisotropic)           !
!                               causes any overlaps among particles                               !
! *********************************************************************************************** !
SUBROUTINE FullOverlapCheck( ContactDistance, Overlap )

! Uses one module: overlap check algorithms
USE OverlapCheckAlgorithms, ONLY: OverlapCheckSPC, OverlapCheckCYL

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: iParticle, jParticle   ! Counters (particle)
INTEGER( Kind= Int64 ) :: iCylinder, jCylinder   ! Counters (cylinders)
INTEGER( Kind= Int64 ) :: pImage                 ! Counters (image)

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                    :: SquaredDistance            ! Magnitude of the vector distance between particles i and j (squared)
REAL( Kind= Real64 )                    :: ContactDistance            ! Contact distance (Perram-Wertheim or Vega-Lago methods)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: VectorDistance             ! Vector distance between particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )    :: iPosition, jPosition       ! Position of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )    :: iOrientation, jOrientation ! Orientation of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )    :: ciPosition, cjPosition     ! Position of cylinders of particles i and j
REAL( Kind= Real64 ), DIMENSION( 0:3 )  :: iQuaternion, jQuaternion   ! Quaternions of particles i and j

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: Overlap     ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: OverlapSPC  ! Detects overlap between two spherocylinders: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: ParallelSPC ! Checks the relative orientation of two spherocylinders : TRUE = parallel orientation; FALSE = non-parallel orientation

! Initialization
Overlap     = .FALSE.
OverlapSPC  = .FALSE.
ParallelSPC = .FALSE.

! First loop represents all particles with indexes i of component Ci
DO iParticle = 1, nParticles - 1
  ! Second loop represents all particles with indexes j of component Cj
  DO jParticle = iParticle + 1, nParticles
    ! Position of particle i
    iPosition(1) = pPositionMC(1,iParticle)
    iPosition(2) = pPositionMC(2,iParticle)
    iPosition(3) = pPositionMC(3,iParticle)
    ! Position of particle j
    jPosition(1) = pPositionMC(1,jParticle)
    jPosition(2) = pPositionMC(2,jParticle)
    jPosition(3) = pPositionMC(3,jParticle)
    ! Orientation of particle i
    iOrientation(1) = pOrientationMC(1,iParticle)
    iOrientation(2) = pOrientationMC(2,iParticle)
    iOrientation(3) = pOrientationMC(3,iParticle)
    ! Orientation of particle j
    jOrientation(1) = pOrientationMC(1,jParticle)
    jOrientation(2) = pOrientationMC(2,jParticle)
    jOrientation(3) = pOrientationMC(3,jParticle)
    ! Quaternion of particle i
    iQuaternion(0) = pQuaternionMC(0,iParticle)
    iQuaternion(1) = pQuaternionMC(1,iParticle)
    iQuaternion(2) = pQuaternionMC(2,iParticle)
    iQuaternion(3) = pQuaternionMC(3,iParticle)
    ! Quaternion of particle j
    jQuaternion(0) = pQuaternionMC(0,jParticle)
    jQuaternion(1) = pQuaternionMC(1,jParticle)
    jQuaternion(2) = pQuaternionMC(2,jParticle)
    jQuaternion(3) = pQuaternionMC(3,jParticle)
    ! Vector distance between particles i and j
    VectorDistance(1) = jPosition(1) - iPosition(1)
    VectorDistance(2) = jPosition(2) - iPosition(2)
    VectorDistance(3) = jPosition(3) - iPosition(3)
    ! Magnitude of the vector distance (squared)
    SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
    ! Preliminary test (circumscribing spheres)
    IF( SquaredDistance <= pSquaredSphereDistance ) THEN
      ! Initialization
      OverlapSPC = .FALSE.
      ! Preliminary test (circumscribing spherocylinders)
      CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, ContactDistance, ParallelSPC )
      ! Overlap criterion
      IF( ContactDistance <= pSquaredSpheroCylinderDistance ) THEN
        ! Initialization
        OverlapSPC = .TRUE.
        ! First loop takes one of the four cylinders from particle i
        DO iCylinder = 1, 4
          ! Position of cylinder of particle i
          ciPosition(1) = cPositionMC(1,iCylinder,iParticle)
          ciPosition(2) = cPositionMC(2,iCylinder,iParticle)
          ciPosition(3) = cPositionMC(3,iCylinder,iParticle)
          ! Second loop takes one of the four cylinders from particle j
          DO jCylinder = 1, 4
            ! Position of cylinder of particle j
            cjPosition(1) = cPositionMC(1,jCylinder,jParticle)
            cjPosition(2) = cPositionMC(2,jCylinder,jParticle)
            cjPosition(3) = cPositionMC(3,jCylinder,jParticle)
            ! Vector distance between cylinders of particles i and j
            VectorDistance(1) = cjPosition(1) - ciPosition(1)
            VectorDistance(2) = cjPosition(2) - ciPosition(2)
            VectorDistance(3) = cjPosition(3) - ciPosition(3)
            ! Magnitude of the vector distance (squared)
            SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
            ! Cutoff distance (sphere circumscribing a spherocylinder circumscribing the cylinder)
            IF( SquaredDistance <= cSquaredSphereDistance ) THEN
              ! Cutoff distance (spherocylinder circumscribing the cylinder)
              CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, ContactDistance, ParallelSPC )
              ! Vega-Lago Criterion
              IF( ContactDistance <= cDiameter * cDiameter ) THEN
                ! Apply periodic boundary conditions on the position of particle j
                cjPosition(1) = ciPosition(1) + VectorDistance(1)
                cjPosition(2) = ciPosition(2) + VectorDistance(2)
                cjPosition(3) = ciPosition(3) + VectorDistance(3)
                ! Time-consuming overlap check
                CALL OverlapCheckCYL( iQuaternion, jQuaternion, iOrientation, jOrientation, VectorDistance, ciPosition, &
                &                     cjPosition, ParallelSPC, Overlap )
                ! Overlap criterion
                IF( Overlap ) THEN
                  ! Overlap detected
                  RETURN
                END IF
              END IF
            END IF
          END DO
        END DO
      END IF
    END IF
  END DO
END DO

! First loop represents all particles with indexes i of component Ci
DO iParticle = 1, nParticles
  ! Second loop represents all particles with indexes j of component Cj
  DO jParticle = 1, nParticles
    ! Loop over images
    DO pImage = 1, nImages
      ! Position of particle i
      iPosition(1) = pPositionMC(1,iParticle)
      iPosition(2) = pPositionMC(2,iParticle)
      iPosition(3) = pPositionMC(3,iParticle)
      ! Position of images of particle j
      jPosition(1) = imPositionMC(1,pImage,jParticle)
      jPosition(2) = imPositionMC(2,pImage,jParticle)
      jPosition(3) = imPositionMC(3,pImage,jParticle)
      ! Orientation of particle i
      iOrientation(1) = pOrientationMC(1,iParticle)
      iOrientation(2) = pOrientationMC(2,iParticle)
      iOrientation(3) = pOrientationMC(3,iParticle)
      ! Orientation of images of particle j
      jOrientation(1) = imOrientationMC(1,pImage,jParticle)
      jOrientation(2) = imOrientationMC(2,pImage,jParticle)
      jOrientation(3) = imOrientationMC(3,pImage,jParticle)
      ! Quaternion of particle i
      iQuaternion(0) = pQuaternionMC(0,iParticle)
      iQuaternion(1) = pQuaternionMC(1,iParticle)
      iQuaternion(2) = pQuaternionMC(2,iParticle)
      iQuaternion(3) = pQuaternionMC(3,iParticle)
      ! Quaternion of particle j
      jQuaternion(0) = imQuaternionMC(0,pImage,jParticle)
      jQuaternion(1) = imQuaternionMC(1,pImage,jParticle)
      jQuaternion(2) = imQuaternionMC(2,pImage,jParticle)
      jQuaternion(3) = imQuaternionMC(3,pImage,jParticle)
      ! Vector distance between particles i and j
      VectorDistance(1) = jPosition(1) - iPosition(1)
      VectorDistance(2) = jPosition(2) - iPosition(2)
      VectorDistance(3) = jPosition(3) - iPosition(3)
      ! Magnitude of the vector distance (squared)
      SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
      ! Preliminary test (circumscribing spheres)
      IF( SquaredDistance <= pSquaredSphereDistance ) THEN
        ! Initialization
        OverlapSPC = .FALSE.
        ! Preliminary test (circumscribing spherocylinders)
        CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, ContactDistance, ParallelSPC )
        ! Overlap criterion
        IF( ContactDistance <= pSquaredSpheroCylinderDistance ) THEN
          ! Initialization
          OverlapSPC = .TRUE.
          ! First loop takes one of the four cylinders from particle i
          DO iCylinder = 1, 4
            ! Position of cylinder of particle i
            ciPosition(1) = cPositionMC(1,iCylinder,iParticle)
            ciPosition(2) = cPositionMC(2,iCylinder,iParticle)
            ciPosition(3) = cPositionMC(3,iCylinder,iParticle)
            ! Second loop takes one of the four cylinders from particle j
            DO jCylinder = 1, 4
              ! Position of cylinder of particle j
              cjPosition(1) = imcPositionMC(1,jCylinder,pImage,jParticle)
              cjPosition(2) = imcPositionMC(2,jCylinder,pImage,jParticle)
              cjPosition(3) = imcPositionMC(3,jCylinder,pImage,jParticle)
              ! Vector distance between cylinders of particles i and j
              VectorDistance(1) = cjPosition(1) - ciPosition(1)
              VectorDistance(2) = cjPosition(2) - ciPosition(2)
              VectorDistance(3) = cjPosition(3) - ciPosition(3)
              ! Magnitude of the vector distance (squared)
              SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
              ! Cutoff distance (sphere circumscribing a spherocylinder circumscribing the cylinder)
              IF( SquaredDistance <= cSquaredSphereDistance ) THEN
                ! Cutoff distance (spherocylinder circumscribing the cylinder)
                CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, ContactDistance, ParallelSPC )
                ! Vega-Lago Criterion
                IF( ContactDistance <= cDiameter * cDiameter ) THEN
                  ! Apply periodic boundary conditions on the position of particle j
                  cjPosition(1) = ciPosition(1) + VectorDistance(1)
                  cjPosition(2) = ciPosition(2) + VectorDistance(2)
                  cjPosition(3) = ciPosition(3) + VectorDistance(3)
                  ! Time-consuming overlap check
                  CALL OverlapCheckCYL( iQuaternion, jQuaternion, iOrientation, jOrientation, VectorDistance, ciPosition, &
                  &                     cjPosition, ParallelSPC, Overlap )
                  ! Overlap criterion
                  IF( Overlap ) THEN
                    ! Overlap detected
                    RETURN
                  END IF
                END IF
              END IF
            END DO
          END DO
        END IF
      END IF
    END DO
  END DO
END DO

! No overlaps
RETURN

END SUBROUTINE FullOverlapCheck

! *********************************************************************************************** !
! This subroutine checks if a random displacement (translation or rotation) of a fixed particle i !
!                           causes any overlaps with other particles j                            !
! *********************************************************************************************** !
SUBROUTINE ParticleOverlapCheck( iParticle, iQuaternion, iQuaternionIMG, iOrientation, iOrientationIMG, iPosition, iPositionIMG, &
&                                ciPosition, ciPositionIMG, ContactDistance, Overlap )

! Uses one module: overlap check algorithms
USE OverlapCheckAlgorithms, ONLY: OverlapCheckSPC, OverlapCheckCYL

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: iParticle, jParticle ! Counters (particles)
INTEGER( Kind= Int64 ) :: iCylinder, jCylinder ! Counters (cylinders)
INTEGER( Kind= Int64 ) :: pImage               ! Counters (images)

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                             :: ContactDistance            ! Vega-Lago Contact Distance (variable)
REAL( Kind= Real64 )                             :: SquaredDistance            ! Vector distance between particles i and j (squared)
REAL( Kind= Real64 ), DIMENSION( 3 )             :: iPosition, jPosition       ! Position of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )             :: iOrientation, jOrientation ! Orientation of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )             :: VectorDistance             ! Vector distance between particles i and j
REAL( Kind= Real64 ), DIMENSION( 0:3 )           :: iQuaternion, jQuaternion   ! Quaternion of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3, 4 )          :: ciPosition, cjPosition     ! Position of cylinders of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3, nImages )    :: iPositionIMG               ! Position of images of particle i
REAL( Kind= Real64 ), DIMENSION( 3, nImages )    :: iOrientationIMG            ! Orientation of images of particle i
REAL( Kind= Real64 ), DIMENSION( 0:3, nImages )  :: iQuaternionIMG             ! Quaternion of images of particle i
REAL( Kind= Real64 ), DIMENSION( 3, 4, nImages ) :: ciPositionIMG              ! Position of cylinder images of particle i

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: Overlap     ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: OverlapSPC  ! Detects overlap between two spherocylinders: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: ParallelSPC ! Checks the relative orientation of two spherocylinders : TRUE = parallel orientation; FALSE = non-parallel orientation

! Initialization
Overlap     = .FALSE.
OverlapSPC  = .FALSE.
ParallelSPC = .FALSE.

! *********************************************************************************************** !
! Overlap Check - Central Particles                                                               !
! *********************************************************************************************** !
! First loop takes only particles whose j-indexes are below the i-index of the fixed particle
DO jParticle = 1, iParticle - 1
  ! Position of particle j
  jPosition(1) = pPositionMC(1,jParticle)
  jPosition(2) = pPositionMC(2,jParticle)
  jPosition(3) = pPositionMC(3,jParticle)
  ! Orientation of particle j
  jOrientation(1) = pOrientationMC(1,jParticle)
  jOrientation(2) = pOrientationMC(2,jParticle)
  jOrientation(3) = pOrientationMC(3,jParticle)
  ! Quaternion of particle j
  jQuaternion(0) = pQuaternionMC(0,jParticle)
  jQuaternion(1) = pQuaternionMC(1,jParticle)
  jQuaternion(2) = pQuaternionMC(2,jParticle)
  jQuaternion(3) = pQuaternionMC(3,jParticle)
  ! Vector distance between particles i and j
  VectorDistance(1) = jPosition(1) - iPosition(1)
  VectorDistance(2) = jPosition(2) - iPosition(2)
  VectorDistance(3) = jPosition(3) - iPosition(3)
  ! Magnitude of the vector distance (squared)
  SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
  ! Cutoff distance (sphere circumscribing the non-convex body)
  IF( SquaredDistance <= pSquaredSphereDistance ) THEN
    ! Initialization
    OverlapSPC = .FALSE.
    ! Preliminary test (circumscribing spherocylinders)
    CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, ContactDistance, ParallelSPC )
    ! Overlap criterion
    IF( ContactDistance <= pSquaredSpheroCylinderDistance ) THEN
      ! Initialization
      OverlapSPC = .TRUE.
      ! First loop takes one of the four cylinders from particle i
      DO iCylinder = 1, 4
        ! Second loop takes one of the four cylinders from particle j
        DO jCylinder = 1, 4
          ! Position of cylinder of particle j
          cjPosition(1,jCylinder) = cPositionMC(1,jCylinder,jParticle)
          cjPosition(2,jCylinder) = cPositionMC(2,jCylinder,jParticle)
          cjPosition(3,jCylinder) = cPositionMC(3,jCylinder,jParticle)
          ! Vector distance between cylinders of particles i and j
          VectorDistance(1) = cjPosition(1,jCylinder) - ciPosition(1,iCylinder)
          VectorDistance(2) = cjPosition(2,jCylinder) - ciPosition(2,iCylinder)
          VectorDistance(3) = cjPosition(3,jCylinder) - ciPosition(3,iCylinder)
          ! Magnitude of the vector distance (squared)
          SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
          ! Cutoff distance (sphere circumscribing a spherocylinder circumscribing the cylinder)
          IF( SquaredDistance <= cSquaredSphereDistance ) THEN
            ! Cutoff distance (spherocylinder circumscribing the cylinder)
            CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, ContactDistance, ParallelSPC )
            ! Vega-Lago Criterion
            IF( ContactDistance <= cDiameter * cDiameter ) THEN
              ! Apply periodic boundary conditions on the position of particle j
              cjPosition(1,jCylinder) = ciPosition(1,iCylinder) + VectorDistance(1)
              cjPosition(2,jCylinder) = ciPosition(2,iCylinder) + VectorDistance(2)
              cjPosition(3,jCylinder) = ciPosition(3,iCylinder) + VectorDistance(3)
              ! Check overlap between cylinders
              CALL OverlapCheckCYL( iQuaternion, jQuaternion, iOrientation, jOrientation, VectorDistance, &
              &                     ciPosition(:,iCylinder), cjPosition(:,jCylinder), ParallelSPC, Overlap )
              ! Overlap criterion
              IF( Overlap ) THEN
                ! Return immediately
                RETURN
              END IF
            END IF
          END IF
        END DO
      END DO
    END IF
  END IF
END DO
! Second loop takes only particles whose j-indexes are above the i-index of the fixed particle
DO jParticle = iParticle + 1, nParticles
! Position of particle j
  jPosition(1) = pPositionMC(1,jParticle)
  jPosition(2) = pPositionMC(2,jParticle)
  jPosition(3) = pPositionMC(3,jParticle)
  ! Orientation of particle j
  jOrientation(1) = pOrientationMC(1,jParticle)
  jOrientation(2) = pOrientationMC(2,jParticle)
  jOrientation(3) = pOrientationMC(3,jParticle)
  ! Quaternion of particle j
  jQuaternion(0) = pQuaternionMC(0,jParticle)
  jQuaternion(1) = pQuaternionMC(1,jParticle)
  jQuaternion(2) = pQuaternionMC(2,jParticle)
  jQuaternion(3) = pQuaternionMC(3,jParticle)
  ! Vector distance between particles i and j
  VectorDistance(1) = jPosition(1) - iPosition(1)
  VectorDistance(2) = jPosition(2) - iPosition(2)
  VectorDistance(3) = jPosition(3) - iPosition(3)
  ! Magnitude of the vector distance (squared)
  SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
  ! Cutoff distance (sphere circumscribing the non-convex body)
  IF( SquaredDistance <= pSquaredSphereDistance ) THEN
    ! Initialization
    OverlapSPC = .FALSE.
    ! Preliminary test (circumscribing spherocylinders)
    CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, ContactDistance, ParallelSPC )
    ! Overlap criterion
    IF( ContactDistance <= pSquaredSpheroCylinderDistance ) THEN
      ! Initialization
      OverlapSPC = .TRUE.
      ! First loop takes one of the four cylinders from particle i
      DO iCylinder = 1, 4
        ! Second loop takes one of the four cylinders from particle j
        DO jCylinder = 1, 4
          ! Position of cylinder of particle j
          cjPosition(1,jCylinder) = cPositionMC(1,jCylinder,jParticle)
          cjPosition(2,jCylinder) = cPositionMC(2,jCylinder,jParticle)
          cjPosition(3,jCylinder) = cPositionMC(3,jCylinder,jParticle)
          ! Vector distance between cylinders of particles i and j
          VectorDistance(1) = cjPosition(1,jCylinder) - ciPosition(1,iCylinder)
          VectorDistance(2) = cjPosition(2,jCylinder) - ciPosition(2,iCylinder)
          VectorDistance(3) = cjPosition(3,jCylinder) - ciPosition(3,iCylinder)
          ! Magnitude of the vector distance (squared)
          SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
          ! Cutoff distance (sphere circumscribing a spherocylinder circumscribing the cylinder)
          IF( SquaredDistance <= cSquaredSphereDistance ) THEN
            ! Cutoff distance (spherocylinder circumscribing the cylinder)
            CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, ContactDistance, ParallelSPC )
            ! Vega-Lago Criterion
            IF( ContactDistance <= cDiameter * cDiameter ) THEN
              ! Apply periodic boundary conditions on the position of particle j
              cjPosition(1,jCylinder) = ciPosition(1,iCylinder) + VectorDistance(1)
              cjPosition(2,jCylinder) = ciPosition(2,iCylinder) + VectorDistance(2)
              cjPosition(3,jCylinder) = ciPosition(3,iCylinder) + VectorDistance(3)
              ! Check overlap between cylinders
              CALL OverlapCheckCYL( iQuaternion, jQuaternion, iOrientation, jOrientation, VectorDistance, &
              &                     ciPosition(:,iCylinder), cjPosition(:,jCylinder), ParallelSPC, Overlap )
              ! Overlap criterion
              IF( Overlap ) THEN
                ! Return immediately
                RETURN
              END IF
            END IF
          END IF
        END DO
      END DO
    END IF
  END IF
END DO

! *********************************************************************************************** !
! Overlap Check - Central Particles and Other Particles' Images                                   !
! *********************************************************************************************** !
! Third loop takes only images of particles whose j-indexes are below the i-index of the fixed particle
DO jParticle = 1, iParticle - 1
  ! Images
  DO pImage = 1, nImages
    ! Position of images of particle j
    jPosition(1) = imPositionMC(1,pImage,jParticle)
    jPosition(2) = imPositionMC(2,pImage,jParticle)
    jPosition(3) = imPositionMC(3,pImage,jParticle)
    ! Orientation of images of particle j
    jOrientation(1) = imOrientationMC(1,pImage,jParticle)
    jOrientation(2) = imOrientationMC(2,pImage,jParticle)
    jOrientation(3) = imOrientationMC(3,pImage,jParticle)
    ! Quaternion of particle j
    jQuaternion(0) = imQuaternionMC(0,pImage,jParticle)
    jQuaternion(1) = imQuaternionMC(1,pImage,jParticle)
    jQuaternion(2) = imQuaternionMC(2,pImage,jParticle)
    jQuaternion(3) = imQuaternionMC(3,pImage,jParticle)
    ! Vector distance between particles i and j
    VectorDistance(1) = jPosition(1) - iPosition(1)
    VectorDistance(2) = jPosition(2) - iPosition(2)
    VectorDistance(3) = jPosition(3) - iPosition(3)
    ! Magnitude of the vector distance (squared)
    SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
    ! Preliminary test (circumscribing spheres)
    IF( SquaredDistance <= pSquaredSphereDistance ) THEN
      ! Initialization
      OverlapSPC = .FALSE.
      ! Preliminary test (circumscribing spherocylinders)
      CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, ContactDistance, ParallelSPC )
      ! Overlap criterion
      IF( ContactDistance <= pSquaredSpheroCylinderDistance ) THEN
        ! Initialization
        OverlapSPC = .TRUE.
        ! First loop takes one of the four cylinders from particle i
        DO iCylinder = 1, 4
          ! Second loop takes one of the four cylinders from particle j
          DO jCylinder = 1, 4
            ! Position of cylinder of particle j
            cjPosition(1,jCylinder) = imcPositionMC(1,jCylinder,pImage,jParticle)
            cjPosition(2,jCylinder) = imcPositionMC(2,jCylinder,pImage,jParticle)
            cjPosition(3,jCylinder) = imcPositionMC(3,jCylinder,pImage,jParticle)
            ! Vector distance between cylinders of particles i and j
            VectorDistance(1) = cjPosition(1,jCylinder) - ciPosition(1,iCylinder)
            VectorDistance(2) = cjPosition(2,jCylinder) - ciPosition(2,iCylinder)
            VectorDistance(3) = cjPosition(3,jCylinder) - ciPosition(3,iCylinder)
            ! Magnitude of the vector distance (squared)
            SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
            ! Cutoff distance (sphere circumscribing a spherocylinder circumscribing the cylinder)
            IF( SquaredDistance <= cSquaredSphereDistance ) THEN
              ! Cutoff distance (spherocylinder circumscribing the cylinder)
              CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, ContactDistance, ParallelSPC )
              ! Vega-Lago Criterion
              IF( ContactDistance <= cDiameter * cDiameter ) THEN
                ! Apply periodic boundary conditions on the position of particle j
                cjPosition(1,jCylinder) = ciPosition(1,iCylinder) + VectorDistance(1)
                cjPosition(2,jCylinder) = ciPosition(2,iCylinder) + VectorDistance(2)
                cjPosition(3,jCylinder) = ciPosition(3,iCylinder) + VectorDistance(3)
                ! Check overlap between cylinders
                CALL OverlapCheckCYL( iQuaternion, jQuaternion, iOrientation, jOrientation, VectorDistance, &
                &                     ciPosition(:,iCylinder), cjPosition(:,jCylinder), ParallelSPC, Overlap )
                ! Overlap criterion
                IF( Overlap ) THEN
                  ! Return immediately
                  RETURN
                END IF
              END IF
            END IF
          END DO
        END DO
      END IF
    END IF
  END DO
END DO
! Fourth loop takes only images of particles whose j-indexes are above the i-index of the fixed particle
DO jParticle = iParticle + 1, nParticles
  ! Images
  DO pImage = 1, nImages
    ! Position of images of particle j
    jPosition(1) = imPositionMC(1,pImage,jParticle)
    jPosition(2) = imPositionMC(2,pImage,jParticle)
    jPosition(3) = imPositionMC(3,pImage,jParticle)
    ! Orientation of images of particle j
    jOrientation(1) = imOrientationMC(1,pImage,jParticle)
    jOrientation(2) = imOrientationMC(2,pImage,jParticle)
    jOrientation(3) = imOrientationMC(3,pImage,jParticle)
    ! Quaternion of particle j
    jQuaternion(0) = imQuaternionMC(0,pImage,jParticle)
    jQuaternion(1) = imQuaternionMC(1,pImage,jParticle)
    jQuaternion(2) = imQuaternionMC(2,pImage,jParticle)
    jQuaternion(3) = imQuaternionMC(3,pImage,jParticle)
    ! Vector distance between particles i and j
    VectorDistance(1) = jPosition(1) - iPosition(1)
    VectorDistance(2) = jPosition(2) - iPosition(2)
    VectorDistance(3) = jPosition(3) - iPosition(3)
    ! Magnitude of the vector distance (squared)
    SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
    ! Preliminary test (circumscribing spheres)
    IF( SquaredDistance <= pSquaredSphereDistance ) THEN
      ! Initialization
      OverlapSPC = .FALSE.
      ! Preliminary test (circumscribing spherocylinders)
      CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, ContactDistance, ParallelSPC )
      ! Overlap criterion
      IF( ContactDistance <= pSquaredSpheroCylinderDistance ) THEN
        ! Initialization
        OverlapSPC = .TRUE.
        ! First loop takes one of the four cylinders from particle i
        DO iCylinder = 1, 4
          ! Second loop takes one of the four cylinders from particle j
          DO jCylinder = 1, 4
            ! Position of cylinder of particle j
            cjPosition(1,jCylinder) = imcPositionMC(1,jCylinder,pImage,jParticle)
            cjPosition(2,jCylinder) = imcPositionMC(2,jCylinder,pImage,jParticle)
            cjPosition(3,jCylinder) = imcPositionMC(3,jCylinder,pImage,jParticle)
            ! Vector distance between cylinders of particles i and j
            VectorDistance(1) = cjPosition(1,jCylinder) - ciPosition(1,iCylinder)
            VectorDistance(2) = cjPosition(2,jCylinder) - ciPosition(2,iCylinder)
            VectorDistance(3) = cjPosition(3,jCylinder) - ciPosition(3,iCylinder)
            ! Magnitude of the vector distance (squared)
            SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
            ! Cutoff distance (sphere circumscribing a spherocylinder circumscribing the cylinder)
            IF( SquaredDistance <= cSquaredSphereDistance ) THEN
              ! Cutoff distance (spherocylinder circumscribing the cylinder)
              CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, ContactDistance, ParallelSPC )
              ! Vega-Lago Criterion
              IF( ContactDistance <= cDiameter * cDiameter ) THEN
                ! Apply periodic boundary conditions on the position of particle j
                cjPosition(1,jCylinder) = ciPosition(1,iCylinder) + VectorDistance(1)
                cjPosition(2,jCylinder) = ciPosition(2,iCylinder) + VectorDistance(2)
                cjPosition(3,jCylinder) = ciPosition(3,iCylinder) + VectorDistance(3)
                ! Check overlap between cylinders
                CALL OverlapCheckCYL( iQuaternion, jQuaternion, iOrientation, jOrientation, VectorDistance, &
                &                     ciPosition(:,iCylinder), cjPosition(:,jCylinder), ParallelSPC, Overlap )
                ! Overlap criterion
                IF( Overlap ) THEN
                  ! Return immediately
                  RETURN
                END IF
              END IF
            END IF
          END DO
        END DO
      END IF
    END IF
  END DO
END DO

! *********************************************************************************************** !
! Overlap Check - Central Particles and Self-Images                                               !
! *********************************************************************************************** !
! Fifth loop takes only the image of particle whose j-index is the same as the i-index of the fixed particle
DO pImage = 1, nImages
  ! Position of images of particle i
  jPosition(1) = iPositionIMG(1,pImage)
  jPosition(2) = iPositionIMG(2,pImage)
  jPosition(3) = iPositionIMG(3,pImage)
  ! Orientation of images of particle i
  jOrientation(1) = iOrientationIMG(1,pImage)
  jOrientation(2) = iOrientationIMG(2,pImage)
  jOrientation(3) = iOrientationIMG(3,pImage)
  ! Quaternion of particle i
  jQuaternion(0) = iQuaternionIMG(0,pImage)
  jQuaternion(1) = iQuaternionIMG(1,pImage)
  jQuaternion(2) = iQuaternionIMG(2,pImage)
  jQuaternion(3) = iQuaternionIMG(3,pImage)
  ! Vector distance between particles i and j
  VectorDistance(1) = iPosition(1) - jPosition(1)
  VectorDistance(2) = iPosition(2) - jPosition(2)
  VectorDistance(3) = iPosition(3) - jPosition(3)
  ! Magnitude of the vector distance (squared)
  SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
  ! Preliminary test (circumscribing spheres)
  IF( SquaredDistance <= pSquaredSphereDistance ) THEN
    ! Initialization
    OverlapSPC = .FALSE.
    ! Preliminary test (circumscribing spherocylinders)
    CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, ContactDistance, ParallelSPC )
    ! Overlap criterion
    IF( ContactDistance <= pSquaredSpheroCylinderDistance ) THEN
      ! Initialization
      OverlapSPC = .TRUE.
      ! First loop takes one of the four cylinders from particle i
      DO iCylinder = 1, 4
        ! Second loop takes one of the four cylinders from particle j
        DO jCylinder = 1, 4
          ! Position of cylinder of particle i
          cjPosition(1,jCylinder) = ciPositionIMG(1,jCylinder,pImage)
          cjPosition(2,jCylinder) = ciPositionIMG(2,jCylinder,pImage)
          cjPosition(3,jCylinder) = ciPositionIMG(3,jCylinder,pImage)
          ! Vector distance between cylinders of particles i and j
          VectorDistance(1) = cjPosition(1,jCylinder) - ciPosition(1,iCylinder)
          VectorDistance(2) = cjPosition(2,jCylinder) - ciPosition(2,iCylinder)
          VectorDistance(3) = cjPosition(3,jCylinder) - ciPosition(3,iCylinder)
          ! Magnitude of the vector distance (squared)
          SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
          ! Cutoff distance (sphere circumscribing a spherocylinder circumscribing the cylinder)
          IF( SquaredDistance <= cSquaredSphereDistance ) THEN
            ! Cutoff distance (spherocylinder circumscribing the cylinder)
            CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, ContactDistance, ParallelSPC )
            ! Vega-Lago Criterion
            IF( ContactDistance <= cDiameter * cDiameter ) THEN
              ! Apply periodic boundary conditions on the position of particle j
              cjPosition(1,jCylinder) = ciPosition(1,iCylinder) + VectorDistance(1)
              cjPosition(2,jCylinder) = ciPosition(2,iCylinder) + VectorDistance(2)
              cjPosition(3,jCylinder) = ciPosition(3,iCylinder) + VectorDistance(3)
              ! Check overlap between cylinders
              CALL OverlapCheckCYL( iQuaternion, jQuaternion, iOrientation, jOrientation, VectorDistance, &
              &                     ciPosition(:,iCylinder), cjPosition(:,jCylinder), ParallelSPC, Overlap )
              ! Overlap criterion
              IF( Overlap ) THEN
                ! Return immediately
                RETURN
              END IF
            END IF
          END IF
        END DO
      END DO
    END IF
  END IF
END DO

RETURN ! No overlaps detected

END SUBROUTINE ParticleOverlapCheck

END MODULE