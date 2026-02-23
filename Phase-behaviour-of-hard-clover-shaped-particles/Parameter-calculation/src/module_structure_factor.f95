MODULE Structure

! Uses 4 modules: global variables, vector operations, quaternions, and image lists
USE GlobalVariables
USE VectorOperations
USE Quaternions
USE ImageLists

CONTAINS

! ********************************************************************************************************************************* !
! This subroutine calculates the structure factor of a molecular configuration (2D).                                                !
!                                                                                                                                   !
! ● The first step is to rotate the box in order to align the phase director with the Z-axis. This will rotate the box vectors, as  !
!   well as the quaternions and positions of the particles, making the analysis easier since the structure factor is calculated in  !
!   2D.                                                                                                                             !
!                                                                                                                                   !
! ● The structure factor is calculated for each plane (layer) of the configuration. For that reason, it is recommended to enable    !
!   <ComputeImages> to include the images of the particles in the plane.                                                            !
!                                                                                                                                   !
! ● Before defining the planes, the closest neighbors of each particle are found. The closest neighbors are those particles (and    !
!   images) that lie within the same layer/plane. This is done by projecting the vector distance between two particles onto the     !
!   normal of the plane that now coincides with the Z-axis. The projected vector (<ProjectedVector>) must have a magnitude less     !
!   than half the length of the particle to be considered within the same plane/layer.                                              !
!                                                                                                                                   !
! ● To prevent double counting the same plane, one particle is selected as the head of the plane.                                   !
!                                                                                                                                   !
! ● The structure factor is calculated using the formula:                                                                           !
!                                                                                                                                   !
!     S(k) = Σ Sₚ(k) / Nₚ,                                                                                                             !
!                                                                                                                                   !
!   where Nₚ is the number of planes and Sₚ(k) is the structure factor of a plane. The structure factor of a plane is calculated     !
!   using the formula:                                                                                                              !
!                                                                                                                                   !
!     Sₚ(k) = | Σⱼ exp( - i × k · rⱼ ) |² / ( N + Nᵢ )ₚ,                                                                             !
!                                                                                                                                   !
!   where k is the wave vector, rⱼ is the position of the particle (or image) j, and ( N + Nᵢ )ₚ is the number of particles         !
!   (and their images) in the plane.                                                                                                !
!                                                                                                                                   !
! ● The wave vector is calculated for each point in the 2D grid and is scaled by the particle diameter and an arbitrary factor.     !
!                                                                                                                                   !
! ● The term | … |² is the square of the magnitude of the complex number, expressed as:                                             !
!                                                                                                                                   !
!     | a + i × b |² = ( a + i × b ) × ( a - i × b ) = a² + b².                                                                     !
!                                                                                                                                   !
! For a better performance, this subroutine is parallelized using OpenMP.                                                           !
! ********************************************************************************************************************************* !
! For more information, please read: Avendaño, C. and Escobedo, F. Soft Matter, 2012, 8, 4675–4681.                                 !
!                                    Pandey, R. B. and Farmer B. L. Struct. Chem., 2017, 28, 625–633.                               !
! ********************************************************************************************************************************* !
SUBROUTINE Structure_Factor( CurrentFrame, TotalFrames, nStruc, PhaseDirector, StructureFactor, gOrientationOrder )

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 )               :: cCylinder, iParticle, jParticle, iNeighbor ! Counters
INTEGER( Kind= Int64 )               :: iStruc, jStruc, iImage, iPlane, kLayer     ! Counters
INTEGER( Kind= Int64 )               :: iNeighborImage                             ! Counters
INTEGER( Kind= Int64 )               :: nImages                                    ! Number of images
INTEGER( Kind= Int64 )               :: OldLenMessage, MaxLenMessage               ! Length of the print message
INTEGER( Kind= Int64 ), INTENT( IN ) :: CurrentFrame, TotalFrames                  ! Frames of the configuration file
INTEGER( Kind= Int64 ), INTENT( IN ) :: nStruc                                     ! Size of the structure (grid)

! INTEGER VARIABLES (ARRAY,ALLOCATABLE)
INTEGER( Kind= Int64 ), DIMENSION( : ), ALLOCATABLE    :: PlaneHeading           ! Heading of the plane (index of the particle)
INTEGER( Kind= Int64 ), DIMENSION( : ), ALLOCATABLE    :: pPlane                 ! Particles inside a plane
INTEGER( Kind= Int64 ), DIMENSION( :, : ), ALLOCATABLE :: ClosestNeighbors       ! Closest neighbors
INTEGER( Kind= Int64 ), DIMENSION( :, : ), ALLOCATABLE :: ClosestNeighborsImages ! Closest neighbors (images)

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 ) :: RotationAxisMagnitude ! Magnitude of the rotation axis
REAL( Kind= Real64 ) :: RotationAngle         ! Rotation angle
REAL( Kind= Real64 ) :: BoxVolumeBoxRotation  ! Box volume (after rotating the box)
REAL( Kind= Real64 ) :: CossineTerm           ! Sum of cossine terms
REAL( Kind= Real64 ) :: SineTerm              ! Sum of sine terms
REAL( Kind= Real64 ) :: Angle                 ! Angle
REAL( Kind= Real64 ) :: gOrientationOrder     ! Angle

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 3 )                                               :: RotatedOrientationX         ! Wave vector
REAL( Kind= Real64 ), DIMENSION( 3 )                                               :: wVector                     ! Wave vector
REAL( Kind= Real64 ), DIMENSION( 3 )                                               :: VectorDistance              ! Vector distance
REAL( Kind= Real64 ), DIMENSION( 3 )                                               :: ProjectedVector             ! Vector distance (projected)
REAL( Kind= Real64 ), DIMENSION( 3 )                                               :: ScalingDistanceUnitBox      ! Scaled position (unit box)
REAL( Kind= Real64 ), DIMENSION( 3 )                                               :: RotatedVector               ! Rotated vector
REAL( Kind= Real64 ), DIMENSION( 9 )                                               :: BoxLengthBoxRotation        ! Box length (after rotating the box)
REAL( Kind= Real64 ), DIMENSION( 9 )                                               :: BoxLengthInverseBoxRotation ! Inverse of box length (after rotating the box)
REAL( Kind= Real64 ), DIMENSION( 3 )                                               :: BoxVectorX                  ! Box vector, X
REAL( Kind= Real64 ), DIMENSION( 3 )                                               :: BoxVectorY                  ! Box vector, Y
REAL( Kind= Real64 ), DIMENSION( 3 )                                               :: BoxVectorZ                  ! Box vector, Z
REAL( Kind= Real64 ), DIMENSION( 3 )                                               :: RotationVector              ! Rotation vector
REAL( Kind= Real64 ), DIMENSION( 0:3 )                                             :: RotationQuaternion          ! Rotation quaternion
REAL( Kind= Real64 ), DIMENSION( 4, 3 )                                            :: cRotatedPosition            ! Rotated reference position (cylinder or petals of the geometry)
REAL( Kind= Real64 ), DIMENSION( : ), INTENT( IN )                                 :: PhaseDirector               ! Phase director (nematic director)
REAL( Kind= Real64 ), DIMENSION( nParticles, 0:3 )                                 :: pQuaternionBoxRotation      ! Quaternion of the particles (after rotating the box)
REAL( Kind= Real64 ), DIMENSION( nParticles, 3 )                                   :: pPositionBoxRotation        ! Position of the particles (after rotating the box)
REAL( Kind= Real64 ), DIMENSION( nParticles, 4, 3 )                                :: cPosition                   ! Position of the cylinders (after rotating the box)
REAL( Kind= Real64 ), DIMENSION( -nStruc:nStruc, -nStruc:nStruc ), INTENT( INOUT ) :: StructureFactor             ! Structure factor (Sk)

! REAL VARIABLES (ARRAY,ALLOCATABLE)
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE          :: gOrientationOrderPlane      ! Angle
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE       :: ScalingDistanceImageUnitBox ! Scaled position of images (unit box)
REAL( Kind= Real64 ), DIMENSION( :, :, : ), ALLOCATABLE    :: StructureFactorPlane        ! Structure factor of a plane
REAL( Kind= Real64 ), DIMENSION( :, :, : ), ALLOCATABLE    :: imPosition                  ! Position of the images
REAL( Kind= Real64 ), DIMENSION( :, :, :, : ), ALLOCATABLE :: imcPosition                 ! Position of the cylinder images

! LOGICAL VARIABLES (ARRAY)
LOGICAL, DIMENSION( nParticles ) :: ParticlesAlreadyCounted ! Check if the particle has already been counted

! CHARACTER STRINGS
CHARACTER( Len= 1000 ) :: PrintStatus ! Print status

! Initialization
ParticlesAlreadyCounted = .FALSE.

! Number of images per particle in all layers
nImages = 0
DO kLayer = 1, nLayers
  nImages = nImages + ( 24 * kLayer * kLayer ) + 2
END DO

! Allocation
IF( ALLOCATED( ClosestNeighbors ) ) DEALLOCATE( ClosestNeighbors )
ALLOCATE( ClosestNeighbors( nParticles, nParticles ) )
ClosestNeighbors = 0
IF( ALLOCATED( ClosestNeighborsImages ) ) DEALLOCATE( ClosestNeighborsImages )
ALLOCATE( ClosestNeighborsImages( nParticles, nParticles * nImages ) )
ClosestNeighborsImages = 0
IF( ALLOCATED( ScalingDistanceImageUnitBox ) ) DEALLOCATE( ScalingDistanceImageUnitBox )
ALLOCATE( ScalingDistanceImageUnitBox( nImages, 3 ) )
ScalingDistanceImageUnitBox = 0.D0
IF( ALLOCATED( imPosition ) ) DEALLOCATE( imPosition )
ALLOCATE( imPosition( nParticles, nImages, 3 ) )
imPosition = 0.D0
IF( ALLOCATED( imcPosition ) ) DEALLOCATE( imcPosition )
ALLOCATE( imcPosition( nParticles, nImages, 4, 3 ) )
imcPosition = 0.D0
IF( ALLOCATED( PlaneHeading ) ) DEALLOCATE( PlaneHeading )
ALLOCATE( PlaneHeading( nParticles ) )
PlaneHeading = 0

! *********************************************************************************************** !
! Rotate the box in order to align the phase director with the Z-axis                             !
! *********************************************************************************************** !

! Rotation angle
RotationAngle = DACOS( DOT_PRODUCT( BodyFixedAxis, PhaseDirector ) / DSQRT( DOT_PRODUCT( BodyFixedAxis, BodyFixedAxis ) * &
&               DOT_PRODUCT( PhaseDirector, PhaseDirector ) ) )

! Check if the rotation angle is not zero or 180 degrees (already aligned along the z-axis)
IF( DABS( DABS( DCOS( RotationAngle ) ) - 1.D0 ) > 1.D-10 ) THEN
  CALL Cross_Product( PhaseDirector, BodyFixedAxis, RotationVector ) ! Rotation vector
  RotationVector = RotationVector / DSQRT( DOT_PRODUCT( RotationVector, RotationVector ) ) ! Normalization
  RotationAxisMagnitude = DSQRT( DOT_PRODUCT( RotationVector, RotationVector ) ) ! Magnitude of the rotation vector
ELSE
  RotationVector = [ 0.D0, 0.D0, 0.D0 ] ! No rotation
END IF

! Calculate rotation quaternion
RotationQuaternion(0)   = DCOS( RotationAngle / 2.D0 ) ! Real part
RotationQuaternion(1:3) = DSIN( RotationAngle / 2.D0 ) * RotationVector(1:3) ! Imaginary part (vector)

! Box vectors
BoxVectorX = BoxLength(1:3) ! Box vector, X-axis
BoxVectorY = BoxLength(4:6) ! Box vector, Y-axis
BoxVectorZ = BoxLength(7:9) ! Box vector, Z-axis

! Check magnitude of the rotation vector
IF( DABS( RotationAxisMagnitude ) >= 1.D-10 ) THEN
  ! Rotated x-vector
  CALL VectorRotation( BoxVectorX / ( DSQRT( DOT_PRODUCT( BoxVectorX, BoxVectorX ) ) ), RotationQuaternion, RotatedVector )
  ! New x-vector of the simulation box
  BoxLengthBoxRotation(1:3) = DSQRT( DOT_PRODUCT( BoxLength(1:3), BoxLength(1:3) ) ) * RotatedVector
  ! Rotated y-vector
  CALL VectorRotation( BoxVectorY / ( DSQRT( DOT_PRODUCT( BoxVectorY, BoxVectorY ) ) ), RotationQuaternion, RotatedVector )
  ! New y-vector of the simulation box
  BoxLengthBoxRotation(4:6) = DSQRT( DOT_PRODUCT( BoxLength(4:6), BoxLength(4:6) ) ) * RotatedVector
  ! Rotated z-vector
  CALL VectorRotation( BoxVectorZ / ( DSQRT( DOT_PRODUCT( BoxVectorZ, BoxVectorZ ) ) ), RotationQuaternion, RotatedVector )
  ! New z-vector of the simulation box
  BoxLengthBoxRotation(7:9) = DSQRT( DOT_PRODUCT( BoxLength(7:9), BoxLength(7:9) ) ) * RotatedVector
  ! Calculate the new reciprocal box basis vectors
  CALL InverseMatrixCofactorVec( BoxLengthBoxRotation, BoxLengthInverseBoxRotation, BoxVolumeBoxRotation )
  ! Rescale positions and orientations of the particles accordingly (the spatial distribution of particles remains unchanged)
  DO iParticle = 1, nParticles
    ! Transform spatial coordinates using old box dimensions
    CALL MatrixVectorMultiplication( BoxLengthInverse, pPosition(iParticle,:), ScalingDistanceUnitBox )
    ! New spatial coordinates using new box dimensions
    CALL MatrixVectorMultiplication( BoxLengthBoxRotation, ScalingDistanceUnitBox, pPositionBoxRotation(iParticle,:) )
    ! Reorient particles
    CALL QuaternionMultiplication( RotationQuaternion, pQuaternion(iParticle,:), pQuaternionBoxRotation(iParticle,:) )
    ! New position of the cylinders (petals)
    DO cCylinder = 1, 4
      CALL VectorRotation( cReferencePosition(cCylinder,:), pQuaternionBoxRotation(iParticle,0:3), cRotatedPosition(cCylinder,:) ) ! Rotate body-fixed reference position
      cPosition(iParticle,cCylinder,:) = pPositionBoxRotation(iParticle,1:3) + cRotatedPosition(cCylinder,:) ! Position of the cylinders
    END DO
  END DO
END IF

! Construct list of images
IF( ComputeImages ) THEN
  DO iParticle = 1, nParticles
    ! Transform spatial coordinates into coordinates of a unit box
    CALL MatrixVectorMultiplication( BoxLengthInverseBoxRotation, pPositionBoxRotation(iParticle,:), ScalingDistanceUnitBox )
    ! Build image list
    CALL ImageConstruction( ScalingDistanceUnitBox, ScalingDistanceImageUnitBox )
    ! Particle images
    DO iImage = 1, nImages
      ! Transform spatial coordinates back to the original box
      CALL MatrixVectorMultiplication( BoxLengthBoxRotation, ScalingDistanceImageUnitBox(iImage,:), imPosition(iParticle,iImage,:) )
    END DO
    ! Cylinder images
    DO cCylinder = 1, 4
      ! Transform spatial coordinates into coordinates of a unit box
      CALL MatrixVectorMultiplication( BoxLengthInverseBoxRotation, cPosition(iParticle,cCylinder,:), ScalingDistanceUnitBox )
      ! Build image list
      CALL ImageConstruction( ScalingDistanceUnitBox, ScalingDistanceImageUnitBox )
      ! Cylinder images
      DO iImage = 1, nImages
        ! Transform spatial coordinates back to the original box
        CALL MatrixVectorMultiplication( BoxLengthBoxRotation, ScalingDistanceImageUnitBox(iImage,:), &
        &                                imcPosition(iParticle,iImage,cCylinder,:) )
      END DO
    END DO
  END DO
END IF

! Find the neighbors of each particle that lie within the same layer/plane
DO iParticle = 1, nParticles
  ! Initialization
  iNeighbor = 1
  iNeighborImage = 1
  DO jParticle = 1, nParticles
    VectorDistance = pPositionBoxRotation(jParticle,1:3) - pPositionBoxRotation(iParticle,1:3) ! Vector distance between particle i and j
    ProjectedVector = DOT_PRODUCT( VectorDistance, BodyFixedAxis ) * BodyFixedAxis ! Projected vector onto the normal of the plane that now coincides with the Z-axis
    ! Add the neighbors within the same layer/plane
    IF( DOT_PRODUCT( ProjectedVector, ProjectedVector ) <= 0.25D0 * pLengthClover * pLengthClover ) THEN ! The projection onto the normal of the plane must be less than half the length of the particle to be considered within the same plane/layer
      ClosestNeighbors(iParticle,iNeighbor) = jParticle ! Closest neighbors within the same plane/layer
      iNeighbor = iNeighbor + 1 ! Increment the counter
    END IF
    ! Check if the images are included as potential neighbors (this option is particularly useful when a plane has only one particle)
    IF( ComputeImages ) THEN
      DO iImage = 1, nImages
        VectorDistance = imPosition(jParticle,iImage,1:3) - pPositionBoxRotation(iParticle,1:3) ! Vector distance between particle i and the images of particle j
        ProjectedVector = DOT_PRODUCT( VectorDistance, BodyFixedAxis ) * BodyFixedAxis ! Projected vector onto the normal of the plane that now coincides with the Z-axis
        ! Add the closest neighbors within the same layer/plane
        IF( DOT_PRODUCT( ProjectedVector, ProjectedVector ) <= 0.25D0 * pLengthClover * pLengthClover ) THEN ! The projection onto the normal of the plane must be less than half the length of the particle to be considered within the same plane/layer
          ClosestNeighborsImages(iParticle,iNeighborImage) = ( jParticle - 1 ) * nImages + iImage ! Closest neighbors within the same plane/layer
          iNeighborImage = iNeighborImage + 1 ! Increment the counter
        END IF
      END DO
    END IF
  END DO
END DO

! Find unique planes (one particle is selected as the heading of the plane to prevent double counting the same plane)
iPlane = 1
DO iParticle = 1, nParticles
  IF( ParticlesAlreadyCounted(iParticle) ) CYCLE ! Loop if the particle has already been counted
  ParticlesAlreadyCounted(iParticle) = .TRUE. ! Mark the particle as counted
  DO jParticle = 1, COUNT( ClosestNeighbors(iParticle,:) /= 0 ) ! Loop over the neighbors of the particle i
    iNeighbor = ClosestNeighbors(iParticle,jParticle) ! Index of neighbors of the particle i
    ParticlesAlreadyCounted(iNeighbor) = .TRUE. ! Mark the neighbor as counted
  END DO
  PlaneHeading(iPlane) = iParticle ! Heading of the plane
  iPlane = iPlane + 1 ! Increment the counter
END DO

! Initialization
MaxLenMessage = 0

! Allocation
IF( ALLOCATED( StructureFactorPlane ) ) DEALLOCATE( StructureFactorPlane )
ALLOCATE( StructureFactorPlane( COUNT( PlaneHeading /= 0 ), -nStruc:nStruc, -nStruc:nStruc ) )
StructureFactorPlane = 0.D0
ALLOCATE( pPlane( COUNT( PlaneHeading /= 0 ) ) )
pPlane = 0
IF( ALLOCATED( gOrientationOrderPlane ) ) DEALLOCATE( gOrientationOrderPlane )
ALLOCATE( gOrientationOrderPlane( COUNT( PlaneHeading /= 0 ) ) )
gOrientationOrderPlane = 0.D0

! Calculate the structure factor
IF( ComputeStructureFactor ) THEN
  StructureFactor = 0.D0
  ! Loop over all planes
  DO iPlane = 1, COUNT( PlaneHeading /= 0 )
    iParticle = PlaneHeading(iPlane) ! Heading of the plane
    ! Status
    WRITE( PrintStatus, "(9G0)" ) "Reading frame: ", CurrentFrame, " of ", TotalFrames, " | Evaluating plane ", iPlane, " of ", &
    &                             COUNT( PlaneHeading /= 0 ), "..."
    OldLenMessage = LEN( TRIM( PrintStatus ) )
    IF( OldLenMessage > MaxLenMessage ) MaxLenMessage = OldLenMessage
    WRITE( *, "(2G0)", Advance= "No" ) CHAR(13), TRIM( PrintStatus )
    !$OMP PARALLEL DO COLLAPSE( 2 ) DEFAULT( Shared ) &
    !$OMP PRIVATE( iStruc, jStruc, jParticle, iNeighbor, iNeighborImage, CossineTerm, SineTerm, wVector )
    DO iStruc = -nStruc, nStruc
      DO jStruc = 0, nStruc ! Symmetric about the horizontal axis
        ! Wave vector
        wVector(1) = DBLE( iStruc ) * pDiameterClover * 0.01D0 ! The wave vector is scaled by the particle diameter (0.004 is an arbitrary scaling factor)
        wVector(2) = DBLE( jStruc ) * pDiameterClover * 0.01D0 ! The wave vector is scaled by the particle diameter (0.004 is an arbitrary scaling factor)
        ! Initialization
        CossineTerm = 0.D0
        SineTerm = 0.D0
        ! Loop over neighbors within the same layer/plane
        DO jParticle = 1, COUNT( ClosestNeighbors(iParticle,:) /= 0 )
          iNeighbor = ClosestNeighbors(iParticle,jParticle) ! Plane neighbours
          CossineTerm = CossineTerm + DCOS( - DOT_PRODUCT( wVector(1:2), pPositionBoxRotation(iNeighbor,1:2) ) ) ! Real part (the analysis is performed in 2D)
          SineTerm = SineTerm + DSIN( - DOT_PRODUCT( wVector(1:2), pPositionBoxRotation(iNeighbor,1:2) ) ) ! Imaginary part (the analysis is performed in 2D)
        END DO
        ! Compute only if images are included as potential neighbors
        IF( ComputeImages ) THEN
          ! Loop over neighbors within the same layer/plane
          DO jParticle = 1, COUNT( ClosestNeighborsImages(iParticle,:) /= 0 )
            iNeighbor = INT( ClosestNeighborsImages(iParticle,jParticle) / nImages ) ! Index of the neighbor (particle)
            iNeighbor = MERGE( iNeighbor, iNeighbor + 1, MOD( ClosestNeighborsImages(iParticle, jParticle), nImages ) == 0 ) ! Correct the index
            iNeighborImage = ClosestNeighborsImages(iParticle,jParticle) - ( INT( ClosestNeighborsImages(iParticle,jParticle) / &
            &                nImages ) ) * nImages ! Index of the neighbor (image)
            iNeighborImage = MERGE( nImages, iNeighborImage, iNeighborImage == 0 ) ! Correct the index
            CossineTerm = CossineTerm + DCOS( - DOT_PRODUCT( wVector(1:2), imPosition(iNeighbor,iNeighborImage,1:2) ) ) ! Real part (the analysis is performed in 2D)
            SineTerm = SineTerm + DSIN( - DOT_PRODUCT( wVector(1:2), imPosition(iNeighbor,iNeighborImage,1:2) ) ) ! Imaginary part (the analysis is performed in 2D)
          END DO
        END IF
        CossineTerm = CossineTerm * CossineTerm ! Square of the real part
        SineTerm = SineTerm * SineTerm ! Square of the imaginary part
        StructureFactor(iStruc,jStruc) = CossineTerm + SineTerm ! Structure factor
      END DO
    END DO
    !$OMP END PARALLEL DO
    StructureFactorPlane(iPlane,:,:) = StructureFactor(:,:) / DBLE( COUNT( ClosestNeighbors(iParticle,:) /= 0 ) + &
    &                                  COUNT( ClosestNeighborsImages(iParticle,:) /= 0 ) ) ! Average the structure factor over the number of particles (and their images) in the plane
    pPlane(iPlane) = COUNT( ClosestNeighbors(iParticle,:) /= 0 ) + COUNT( ClosestNeighborsImages(iParticle,:) /= 0 ) ! Number of particles (and their images) in the plane
  END DO
  ! Reflect other structure factors horizontally
  DO iPlane = 1, COUNT( PlaneHeading /= 0 )
    DO iStruc = -nStruc, nStruc
      DO jStruc = -nStruc, - 1
        StructureFactorPlane(iPlane,iStruc,jStruc) = StructureFactorPlane(iPlane,-iStruc,-jStruc)
      END DO
    END DO
  END DO
END IF

! Calculate the orientational order of molecular axis
IF( ComputeGlobalOrientOrder ) THEN
  gOrientationOrder = 0.D0
  ! Initialization
  CossineTerm = 0.D0
  SineTerm = 0.D0
  ! Loop over neighbors within the same layer/plane
  DO iParticle = 1, nParticles
    !iNeighbor = ClosestNeighbors(iParticle,jParticle) ! Plane neighbours
    CALL VectorRotation( [1.D0, 0.D0, 0.D0], pQuaternionBoxRotation(iParticle,0:3), RotatedOrientationX ) ! Rotate the x-axis
    RotatedOrientationX(3) = 0.D0 ! Set the z-component to zero
    RotatedOrientationX = RotatedOrientationX / DSQRT( DOT_PRODUCT( RotatedOrientationX, RotatedOrientationX ) ) ! Normalize the vector
    Angle = DOT_PRODUCT( RotatedOrientationX, [1.D0, 0.D0, 0.D0] ) / DSQRT( DOT_PRODUCT( RotatedOrientationX, &
    &                    RotatedOrientationX ) ) ! Angle between the rotated x-axis and the x-axis
    IF( Angle <= -1.D0 ) Angle = -1.D0 ! Check if the angle is less than -1
    IF( Angle >= 1.D0 ) Angle = 1.D0 ! Check if the angle is greater than 1
    Angle = DACOS( Angle ) ! Calculate the angle
    CossineTerm = CossineTerm + DCOS( 4.D0 * MOD( Angle, 2.D0 * cPi / 4.D0 ) ) ! Real part (the analysis is performed in 2D)
    SineTerm = SineTerm + DSIN( 4.D0 * MOD( Angle, 2.D0 * cPi / 4.D0 ) ) ! Imaginary part (the analysis is performed in 2D)
  END DO
  CossineTerm = CossineTerm * CossineTerm ! Square of the real part
  SineTerm = SineTerm * SineTerm ! Square of the imaginary part
  gOrientationOrder = DSQRT( CossineTerm + SineTerm ) / nParticles ! Average the structure factor over the number of particles (and their images) in the plane
END IF

! Structure factor of all planes
StructureFactor = 0.D0
IF( ComputeStructureFactor ) THEN
  DO iPlane = 1, COUNT( PlaneHeading /= 0 ) ! Loop over all planes
    StructureFactor(:,:) = StructureFactor(:,:) + StructureFactorPlane(iPlane,:,:) * DBLE( pPlane(iPlane) ) ! Structure factor of all planes
  END DO
END IF
IF( ComputeStructureFactor ) StructureFactor = StructureFactor / DBLE( SUM( pPlane ) ) ! Average the structure factor over the number of particles (and their images) in the plane

! Status
IF( ComputeStructureFactor ) THEN
  WRITE( PrintStatus, "(4G0)" ) "Reading frame: ", CurrentFrame, " of ", TotalFrames
  OldLenMessage = LEN( TRIM( PrintStatus ) )
  WRITE( *, "(6G0)", Advance= "No" ) CHAR(13), "Reading frame: ", CurrentFrame, " of ", TotalFrames, &
  &                                  REPEAT( " ", MaxLenMessage - OldLenMessage )
END IF

RETURN

END SUBROUTINE Structure_Factor

END MODULE Structure