MODULE BondParameters

! Uses three modules: global variables, vector operations, quaternions, and spherical harmonics
USE GlobalVariables
USE VectorOperations
USE Quaternions
USE Spherical_Harmonics_LM
USE Neighbors

CONTAINS

! ********************************************************************************************************************************* !
! This subroutine calculates the bond parameters (<Psi>) of a system in two-dimensions.                                             !
!                                                                                                                                   !
! ● The bond parameters are calculated for a 2D plane (e.g., a monolayer of particles). The neighbours of each particle are         !
!   identified as the 'n' nearest neighbours. The nth-fold bond order parameter is defined as:                                      !
!                                                                                                                                   !
!     ψₙ = (1 / N) × | Σᵢ (1 / nᵢ) × Σⱼ exp( 𝑖 × n × θᵢⱼ ) |,                                                                        !
!                                                                                                                                   !
!   where N is the total number of particles, nⱼ is the number of neighbours of particle j, and θᵢⱼ is the angle made by the vector !
!   distance between particle i and its nearest neighbours j with the arbitrary axis perpendicular to the phase director.           !
!                                                                                                                                   !
! ● Here, we apply the Euler's formula to transform "exp(𝑖nθ)" into "cos(nθ) + 𝑖 × sin(nθ)". Within the absolute modulus,            !
!   we'll have a sum of real terms and a sum of imaginary terms. To solve it, we can apply the absolute squared magnitude and       !
!   take its square root.                                                                                                           !
!                                                                                                                                   !
! ● The absolute squared magnitude is calculated by multiplying the complex number by its complex conjugate:                        !
!                                                                                                                                   !
!     | Σ (a + 𝑖 × b) |² = [ ( Σ ( a ) + 𝑖 × Σₖ ( b ) ) × ( Σ ( a ) - 𝑖 × Σₖ ( b ) ) ] = { [ Σ ( a ) ]² + [ Σ ( b ) ]² }              !
! ********************************************************************************************************************************* !
SUBROUTINE LocalBondParameters( PhaseDirector, Psi, ChiralOrder, ChiralityParticle )

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 ) :: iParticle, jParticle, kParticle  ! Counters
INTEGER( Kind= Int64 ) :: iNeighbor, iOrder, iAxis         ! Counters
INTEGER( Kind= Int64 ) :: cCylinder                        ! Counters
INTEGER( Kind= Int64 ) :: nVertices, nEdges                ! Number of vertices, faces, and edges
INTEGER( Kind= Int64 ) :: nCandidatesVoronoi               ! Number of candidates (Voronoi tessellation)
INTEGER( Kind= Int64 ) :: iCandidate                       ! Counters (candidates)
INTEGER( Kind= Int64 ) :: iVertex                          ! Counters (vertices)
INTEGER( Kind= Int64 ) :: fVertex                          ! Counters (vertices of a face)
INTEGER( Kind= Int64 ) :: xChiralPositive, xChiralNegative ! Counters (positive and negative chiralities)

! INTEGER VARIABLES (ARRAY)
INTEGER( Kind= Int64 ), DIMENSION( nParticles )     :: ParticleIndices   ! Particle indices
INTEGER( Kind= Int64 ), DIMENSION( nParticles )     :: nFaceEdges        ! Number of edges per face
INTEGER( Kind= Int64 ), DIMENSION( MaxVertices )    :: fVertexIndices    ! Index of the vertices of a face (neighbor)
INTEGER( Kind= Int64 ), DIMENSION( MaxVertices, 3 ) :: vIndex            ! Index of the vertices

! INTEGER VARIABLES (ARRAY,ALLOCATABLE)
INTEGER( Kind= Int64 ), DIMENSION( :, : ), ALLOCATABLE :: ClosestNeighbors ! Closest neighbors of each particle

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 )                :: LargestProjection     ! Largest projection of the vector distance onto the reference vector
REAL( Kind= Real64 )                :: CossineArgument       ! Cossine argument
REAL( Kind= Real64 )                :: RotationAxisMagnitude ! Magnitude of the rotation axis
REAL( Kind= Real64 )                :: RotationAngle         ! Rotation angle
REAL( Kind= Real64 )                :: BoxVolumeBoxRotation  ! Box volume (after rotating the box)
REAL( Kind= Real64 )                :: Angle                 ! Angle between the vector distance and the parallel vector
REAL( Kind= Real64 )                :: SubtendedSolidAngle   ! Solid angle of the spherical triangle
REAL( Kind= Real64 )                :: ChiralityAccumulator  ! Sum of the chirality
REAL( Kind= Real64 ), INTENT( OUT ) :: ChiralOrder           ! Chiral order parameter

! REAL VARIABLES (ARRAY,CONSTANT)
REAL( Kind= Real64 ), DIMENSION( 4, 3 ), PARAMETER :: PlaneAxes = RESHAPE( [ 1.D0, -1.D0, 0.D0, 0.D0, &
&                                                                            0.D0, 0.D0, 1.D0, -1.D0, &
&                                                                            0.D0, 0.D0, 0.D0, 0.D0   &
&                                                                          ], [4, 3] ) ! Body-fixed axes in the plane

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 3 )                      :: NormalVector                ! Normal vector of the plane
REAL( Kind= Real64 ), DIMENSION( 3 )                      :: OrthogonalVector            ! Vector orthogonal to the the reference vector and the vector distance
REAL( Kind= Real64 ), DIMENSION( 3 )                      :: ReferenceVector             ! Reference vector in the plane
REAL( Kind= Real64 ), DIMENSION( 4 )                      :: ReferenceProjection         ! Projection of the vector distance onto the reference vector
REAL( Kind= Real64 ), DIMENSION( 3 )                      :: RotatedVector               ! Rotated vector
REAL( Kind= Real64 ), DIMENSION( 9 )                      :: BoxLengthBoxRotation        ! Box length (after rotating the box)
REAL( Kind= Real64 ), DIMENSION( 9 )                      :: BoxLengthInverseBoxRotation ! Inverse of box length (after rotating the box)
REAL( Kind= Real64 ), DIMENSION( 3 )                      :: BoxVectorX                  ! Box vector, X
REAL( Kind= Real64 ), DIMENSION( 3 )                      :: BoxVectorY                  ! Box vector, Y
REAL( Kind= Real64 ), DIMENSION( 3 )                      :: BoxVectorZ                  ! Box vector, Z
REAL( Kind= Real64 ), DIMENSION( 3 )                      :: RotationVector              ! Rotation vector
REAL( Kind= Real64 ), DIMENSION( 0:3 )                    :: RotationQuaternion          ! Rotation quaternion
REAL( Kind= Real64 ), DIMENSION( 3 )                      :: VectorDistance              ! Vector distance
REAL( Kind= Real64 ), DIMENSION( 3 )                      :: ScalingDistanceUnitBox      ! Scaled distance in the unit box
REAL( Kind= Real64 ), DIMENSION( 4, 3 )                   :: BasisVector                 ! Body-fixed axis written in the system of coordinates
REAL( Kind= Real64 ), DIMENSION( 2, 2 )                   :: EdgeDistances               ! Vector distance between the vertices and centroid of the Voronoi cell
REAL( Kind= Real64 ), DIMENSION( 4, 3 )                   :: cRotatedPosition            ! Rotated reference position (cylinder or petals of the geometry)
REAL( Kind= Real64 ), DIMENSION( nParticles )             :: ChiralityParticle           ! Chirality per particle
REAL( Kind= Real64 ), DIMENSION( nParticles )             :: SquaredDistance             ! Squared distance
REAL( Kind= Real64 ), DIMENSION( nParticles, 0:3 )        :: pQuaternionBoxRotation      ! Quaternion of the particles (after rotating the box)
REAL( Kind= Real64 ), DIMENSION( nParticles, 3 )          :: pPositionBoxRotation        ! Position of the particles (after rotating the box)
REAL( Kind= Real64 ), DIMENSION( nParticles, 4, 3 )       :: cPosition                   ! Position of the cylinders (after rotating the box)
REAL( Kind= Real64 ), DIMENSION( MaxVertices, 3 )         :: vRelativePosition           ! Relative position of the vertices
REAL( Kind= Real64 ), DIMENSION( nParticles, nParticles ) :: SolidAngleNeighbor          ! Subtended solid angle of the neighbors (faces of the Voronoi cell)
REAL( Kind= Real64 ), DIMENSION( : ), INTENT( IN )        :: PhaseDirector               ! Phase director (nematic director)
REAL( Kind= Real64 ), DIMENSION( : ), INTENT( OUT )       :: Psi                         ! Bond parameters (ψ)

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE    :: Cossine                ! Cossine term
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE    :: Sine                   ! Sine term
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE    :: CossineTerm            ! Sum of cossine terms
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE    :: SineTerm               ! Sum of sine terms
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE    :: SquaredDistanceVoronoi ! Squared distance (Voronoi tessellation)
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE :: cRelativeDistance      ! Relative distance of the candidates

! Allocation
IF( ALLOCATED( ClosestNeighbors ) ) DEALLOCATE( ClosestNeighbors )
ALLOCATE( ClosestNeighbors( nParticles, nParticles ) )
ClosestNeighbors = 0
IF( ALLOCATED( Cossine ) ) DEALLOCATE( Cossine )
ALLOCATE( Cossine( MaxNeighbors - 3 ) )
Cossine = 0.D0
IF( ALLOCATED( Sine ) ) DEALLOCATE( Sine )
ALLOCATE( Sine( MaxNeighbors - 3 ) )
Sine = 0.D0
IF( ALLOCATED( CossineTerm ) ) DEALLOCATE( CossineTerm )
ALLOCATE( CossineTerm( MaxNeighbors - 3 ) )
CossineTerm = 0.D0
IF( ALLOCATED( SineTerm ) ) DEALLOCATE( SineTerm )
ALLOCATE( SineTerm( MaxNeighbors - 3 ) )
SineTerm = 0.D0
IF( ALLOCATED( SquaredDistanceVoronoi ) ) DEALLOCATE( SquaredDistanceVoronoi )
ALLOCATE( SquaredDistanceVoronoi( nParticles ) )
SquaredDistanceVoronoi = 0.D0
IF( ALLOCATED( cRelativeDistance ) ) DEALLOCATE( cRelativeDistance )
ALLOCATE( cRelativeDistance( nParticles, 3 ) )
cRelativeDistance = 0.D0

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

! Initialization
SolidAngleNeighbor = 0.D0
ChiralityAccumulator = 0.D0
ChiralityParticle = 0.D0

! Find the neighbors of each particle
!$OMP PARALLEL DO DEFAULT( Shared ) &
!$OMP PRIVATE( iParticle, jParticle, kParticle, NormalVector, VectorDistance, ScalingDistanceUnitBox, SquaredDistance ) &
!$OMP PRIVATE( ParticleIndices, iNeighbor, Cossine, Sine, Angle, CossineArgument, cRelativeDistance, ReferenceVector ) &
!$OMP PRIVATE( SquaredDistanceVoronoi, nCandidatesVoronoi, nVertices, nEdges, nFaceEdges, vRelativePosition, vIndex, iCandidate ) &
!$OMP PRIVATE( iVertex, fVertex, fVertexIndices, EdgeDistances, SubtendedSolidAngle, xChiralPositive, xChiralNegative, iAxis ) &
!$OMP PRIVATE( ReferenceProjection, BasisVector, OrthogonalVector, LargestProjection ) &
!$OMP REDUCTION( +: CossineTerm, SineTerm, ChiralityAccumulator )
DO iParticle = 1, nParticles ! Loop over all particles
  NormalVector = BodyFixedAxis ! Assume the normal vector as the body-fixed axis (z-axis)
  NormalVector = NormalVector / DSQRT( DOT_PRODUCT( NormalVector, NormalVector ) ) ! Normalization
  ! Initialization of the closest neighbors
  IF( NeighborOption2D == 1 ) THEN ! Nearest neighbors [1]
    ! Initialization
    SquaredDistance = HUGE( 0.D0 ) ! Large value
    DO jParticle = 1, nParticles ! Loop over all other particles
      IF( iParticle == jParticle ) CYCLE ! Skip the same particle
      VectorDistance = pPositionBoxRotation(jParticle,1:3) - pPositionBoxRotation(iParticle,1:3) ! Vector distance between particle i and j
      CALL MatrixVectorMultiplication( BoxLengthInverseBoxRotation, VectorDistance, ScalingDistanceUnitBox ) ! Spatial transformation
      ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox ) ! Periodic boundary conditions
      CALL MatrixVectorMultiplication( BoxLengthBoxRotation, ScalingDistanceUnitBox, VectorDistance ) ! Spatial transformation
      SquaredDistance(jParticle) = DOT_PRODUCT( VectorDistance, VectorDistance ) ! Squared distance
    END DO
    ParticleIndices = (/ ( kParticle, kParticle = 1, nParticles ) /) ! Create an array of particle indices from 1 to N
    CALL Sort_Array( ParticleIndices, SquaredDistance ) ! Sort the indices of the particles according to the squared distance in ascending order
    ! Initialization
    iNeighbor = 1
    DO jParticle = 1, nParticles ! Loop over all other particles
      VectorDistance = pPositionBoxRotation(ParticleIndices(jParticle),1:3) - pPositionBoxRotation(iParticle,1:3) ! Vector distance between particle i and its closest neighbors j
      CALL MatrixVectorMultiplication( BoxLengthInverseBoxRotation, VectorDistance, ScalingDistanceUnitBox ) ! Spatial transformation
      ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox ) ! Periodic boundary conditions
      CALL MatrixVectorMultiplication( BoxLengthBoxRotation, ScalingDistanceUnitBox, VectorDistance ) ! Spatial transformation
      VectorDistance = DOT_PRODUCT( VectorDistance, NormalVector ) * NormalVector ! Project the vector distance onto the normal vector (phase director)
      ! Add the closest neighbors within the same layer/plane
      IF( DOT_PRODUCT( VectorDistance, VectorDistance ) <= 0.25D0 * pLengthClover * pLengthClover ) THEN
        ClosestNeighbors(iParticle,iNeighbor) = ParticleIndices(jParticle) ! Closest neighbors within the same layer/plane
        iNeighbor = iNeighbor + 1 ! Increment the counter
      END IF
      IF( iNeighbor > MaxNeighbors ) EXIT ! Maximum number of neighbors reached
    END DO
  ELSE IF( NeighborOption2D == 2 ) THEN ! Voronoi tessellation
    ! Initialization
    nCandidatesVoronoi = 0
    SquaredDistanceVoronoi = HUGE( 0.D0 ) ! Large value
    ! Number of candidates
    DO jParticle = 1, nParticles ! Loop over all other particles
      IF( iParticle == jParticle ) CYCLE ! Skip the same particle
      VectorDistance = pPositionBoxRotation(jParticle,1:3) - pPositionBoxRotation(iParticle,1:3) ! Vector distance between particle i and its closest neighbors j
      CALL MatrixVectorMultiplication( BoxLengthInverseBoxRotation, VectorDistance, ScalingDistanceUnitBox ) ! Spatial transformation
      ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox ) ! Periodic boundary conditions
      CALL MatrixVectorMultiplication( BoxLengthBoxRotation, ScalingDistanceUnitBox, VectorDistance ) ! Spatial transformation
      cRelativeDistance(jParticle,:) = VectorDistance ! Relative distance of the candidates
      SquaredDistanceVoronoi(jParticle) = DOT_PRODUCT( VectorDistance, VectorDistance ) ! Squared distance
      VectorDistance = DOT_PRODUCT( VectorDistance, NormalVector ) * NormalVector ! Project the vector distance onto the normal vector (phase director)
      ! Consider only the candidates within the same layer/plane
      IF( DOT_PRODUCT( VectorDistance, VectorDistance ) <= 0.25D0 * pLengthClover * pLengthClover ) THEN
        nCandidatesVoronoi = nCandidatesVoronoi + 1
      ELSE
        SquaredDistanceVoronoi(jParticle) = HUGE( 0.D0 ) ! Large value
      END IF
    END DO
    ParticleIndices = (/ ( kParticle, kParticle = 1, nParticles ) /) ! Create an array of particle indices from 1 to N
    ! Sort the indices of the particles according to the squared distance in ascending order (Voronoi tessellation)
    CALL Sort_Array_Voronoi( ParticleIndices, cRelativeDistance, SquaredDistanceVoronoi )
    ! Voronoi tessellation
    CALL VoronoiConstruction2D( nCandidatesVoronoi, nVertices, nEdges, cRelativeDistance, SquaredDistanceVoronoi, nFaceEdges, &
    &                           vRelativePosition, vIndex )
    ! Find the closest neighbors within the Voronoi tessellation
    iNeighbor = 1
    DO iCandidate = 1, nCandidatesVoronoi
      IF( nFaceEdges(iCandidate) /= 0 ) THEN
        ClosestNeighbors(iParticle,iNeighbor) = ParticleIndices(iCandidate)
        iNeighbor = iNeighbor + 1
      END IF
    END DO
    ! Find vertices of the Voronoi tessellation
    DO iNeighbor = 1, COUNT( ClosestNeighbors(iParticle,:) /= 0, Dim= 1 ) ! Loop over the neighbors (faces)
      ! Initialization
      fVertexIndices = 0
      fVertex = 1
      DO iVertex = 1, nVertices ! Loop over all vertices of the Voronoi cell
        ! Check if the index of the closest neighbor is associated with the vertex
        IF( ANY( ParticleIndices(vIndex(iVertex,1:2)) == ClosestNeighbors(iParticle,iNeighbor) ) ) THEN
          ! Store the index of the vertex associated with the common face (vertices of a particular face)
          fVertexIndices(fVertex) = iVertex
          fVertex = fVertex + 1
        END IF
      END DO
      ! Relative distances of the vertices to the center of the Voronoi cell
      EdgeDistances(1,1:2) = vRelativePosition(fVertexIndices(1),1:2)
      EdgeDistances(2,1:2) = vRelativePosition(fVertexIndices(2),1:2)
      ! Angle argument
      CossineArgument = DOT_PRODUCT( EdgeDistances(1,1:2), EdgeDistances(2,1:2) ) / ( DSQRT( DOT_PRODUCT( EdgeDistances(1,1:2), &
      &     EdgeDistances(1,1:2) ) ) * DSQRT( DOT_PRODUCT( EdgeDistances(2,1:2), EdgeDistances(2,1:2) ) ) )
      IF( CossineArgument < -1.D0 ) CossineArgument = -1.D0 ! Avoid numerical errors
      IF( CossineArgument > 1.D0 ) CossineArgument = 1.D0 ! Avoid numerical errors
      ! Subtended solid angle of the Voronoi cell
      SubtendedSolidAngle = DACOS( CossineArgument ) ! Angle between the two edges of the Voronoi cell
      SolidAngleNeighbor(iParticle,iNeighbor) = SubtendedSolidAngle
    END DO
  END IF
  ! Initialization
  Cossine = 0.D0
  Sine = 0.D0
  xChiralPositive = 0
  xChiralNegative = 0
  IF( NeighborOption2D == 1 ) THEN ! Nearest neighbors [1]
    DO iNeighbor = 1, COUNT( ClosestNeighbors(iParticle,:) /= 0, Dim= 1 ) ! Loop over closest neighbors within the same layer/plane
      VectorDistance = pPositionBoxRotation(ClosestNeighbors(iParticle,iNeighbor),1:3) - pPositionBoxRotation(iParticle,1:3) ! Vector distance between particle i and its closest neighbors j within the same layer/plane
      CALL MatrixVectorMultiplication( BoxLengthInverseBoxRotation, VectorDistance, ScalingDistanceUnitBox ) ! Spatial transformation
      ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox ) ! Periodic boundary conditions
      CALL MatrixVectorMultiplication( BoxLengthBoxRotation, ScalingDistanceUnitBox, VectorDistance ) ! Spatial transformation
      ! Angle between the vector distance and a reference direction
      Angle = DATAN2( VectorDistance(2) , VectorDistance(1) )
      IF( Angle < 0.0D0 ) Angle = 2.0D0 * cPi + Angle
      ! Cossine and sine functions (Euler's formula)
      DO iOrder = 4, MaxNeighbors
        IF( iNeighbor <= iOrder ) Cossine(iOrder-3) = Cossine(iOrder-3) + DCOS( DBLE( iOrder ) * Angle ) ! Real part
        IF( iNeighbor <= iOrder ) Sine(iOrder-3) = Sine(iOrder-3) + DSIN( DBLE( iOrder ) * Angle ) ! Imaginary part
      END DO
      ! Project the vector distance onto the xy-plane (normal vector)
      VectorDistance(3) = 0.D0
      ! Chiral order
      ReferenceProjection = - HUGE( 0.D0 ) ! Tiny value
      ! Reference vector in the plane associated with the largest projection of the vector distance
      CALL VectorRotation( PlaneAxes(1,:), pQuaternionBoxRotation(iParticle,:), BasisVector(1,:) ) ! Basis vector in the plane [x-axis]
      BasisVector(1,3) = 0.D0 ! Z-component is zero
      ReferenceProjection(1) = DOT_PRODUCT( VectorDistance, BasisVector(1,:) ) ! Projection of the vector distance onto the reference vector
      ReferenceVector = BasisVector(1,:) ! Reference vector in the plane
      LargestProjection = ReferenceProjection(1) ! Largest projection of the vector distance onto the reference vector
      DO iAxis = 2, 4
        CALL VectorRotation( PlaneAxes(iAxis,:), pQuaternionBoxRotation(iParticle,:), BasisVector(iAxis,:) ) ! Basis vector in the plane [x-axis]
        BasisVector(iAxis,3) = 0.D0 ! Z-component is zero
        ReferenceProjection(iAxis) = DOT_PRODUCT( VectorDistance, BasisVector(iAxis,:) ) ! Projection of the vector distance onto the reference vector
        IF( ReferenceProjection(iAxis) > LargestProjection ) THEN
          LargestProjection = ReferenceProjection(iAxis) ! Largest projection of the vector distance onto the reference vector
          ReferenceVector = BasisVector(iAxis,:) ! Reference vector in the plane
        END IF
      END DO
      VectorDistance = VectorDistance / DSQRT( DOT_PRODUCT( VectorDistance, VectorDistance ) ) ! Normalization
      ReferenceVector = ReferenceVector / DSQRT( DOT_PRODUCT( ReferenceVector, ReferenceVector ) ) ! Normalization
      CALL Cross_Product( ReferenceVector, VectorDistance, OrthogonalVector ) ! Orthogonal vector to the reference vector and the vector distance
      OrthogonalVector = OrthogonalVector / DSQRT( DOT_PRODUCT( OrthogonalVector, OrthogonalVector ) ) ! Normalization
      ! Signed angle between the vector distance and the reference vector
      Angle = DOT_PRODUCT( OrthogonalVector, NormalVector) / DABS( DOT_PRODUCT( OrthogonalVector, NormalVector) ) * &
      &       DACOS( DOT_PRODUCT( VectorDistance, ReferenceVector ) / ( &
      &       DSQRT( DOT_PRODUCT( VectorDistance, VectorDistance ) ) * DSQRT( DOT_PRODUCT( ReferenceVector, ReferenceVector ) ) ) ) ! This operation represents the ratio of the signed sine and cosine
      IF( Angle > 0.0D0 ) xChiralPositive = xChiralPositive + 1 ! Positive chirality (clockwise)
      IF( Angle < 0.0D0 ) xChiralNegative = xChiralNegative + 1 ! Negative chirality (counterclockwise)
    END DO
    ! Average the trigonometric functions over the number of neighbors
    DO iOrder = 4, MaxNeighbors
      Cossine(iOrder-3) = Cossine(iOrder-3) / DBLE( iOrder ) ! Real part
      Sine(iOrder-3) = Sine(iOrder-3) / DBLE( iOrder ) ! Imaginary part
    END DO
  ELSE IF( NeighborOption2D == 2 ) THEN ! ! Voronoi tessellation
    DO iNeighbor = 1, COUNT( ClosestNeighbors(iParticle,:) /= 0, Dim= 1 ) ! Loop over closest neighbors within the same layer/plane
      VectorDistance = pPositionBoxRotation(ClosestNeighbors(iParticle,iNeighbor),1:3) - pPositionBoxRotation(iParticle,1:3) ! Vector distance between particle i and its closest neighbors j within the same layer/plane
      CALL MatrixVectorMultiplication( BoxLengthInverseBoxRotation, VectorDistance, ScalingDistanceUnitBox ) ! Spatial transformation
      ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox ) ! Periodic boundary conditions
      CALL MatrixVectorMultiplication( BoxLengthBoxRotation, ScalingDistanceUnitBox, VectorDistance ) ! Spatial transformation
      ! Angle between the vector distance and a reference direction
      Angle = DATAN2( VectorDistance(2) , VectorDistance(1) )
      IF( Angle < 0.0D0 ) Angle = 2.0D0 * cPi + Angle
      ! Cossine and sine functions (Euler's formula) weighted by the solid angle of the Voronoi cell
      DO iOrder = 4, MaxNeighbors
        Cossine(iOrder-3) = Cossine(iOrder-3) + DCOS( DBLE( iOrder ) * Angle ) * SolidAngleNeighbor(iParticle,iNeighbor) ! Real part
        Sine(iOrder-3) = Sine(iOrder-3) + DSIN( DBLE( iOrder ) * Angle ) * SolidAngleNeighbor(iParticle,iNeighbor) ! Imaginary part
      END DO
      ! Project the vector distance onto the xy-plane (normal vector)
      VectorDistance(3) = 0.D0
      ! Chiral order
      ReferenceProjection = - HUGE( 0.D0 ) ! Tiny value
      ! Reference vector in the plane associated with the largest projection of the vector distance
      CALL VectorRotation( PlaneAxes(1,:), pQuaternionBoxRotation(iParticle,:), BasisVector(1,:) ) ! Basis vector in the plane [x-axis]
      BasisVector(1,3) = 0.D0 ! Z-component is zero
      ReferenceProjection(1) = DOT_PRODUCT( VectorDistance, BasisVector(1,:) ) ! Projection of the vector distance onto the reference vector
      ReferenceVector = BasisVector(1,:) ! Reference vector in the plane
      LargestProjection = ReferenceProjection(1) ! Largest projection of the vector distance onto the reference vector
      DO iAxis = 2, 4
        CALL VectorRotation( PlaneAxes(iAxis,:), pQuaternionBoxRotation(iParticle,:), BasisVector(iAxis,:) ) ! Basis vector in the plane [x-axis]
        BasisVector(iAxis,3) = 0.D0 ! Z-component is zero
        ReferenceProjection(iAxis) = DOT_PRODUCT( VectorDistance, BasisVector(iAxis,:) ) ! Projection of the vector distance onto the reference vector
        IF( ReferenceProjection(iAxis) > LargestProjection ) THEN
          LargestProjection = ReferenceProjection(iAxis) ! Largest projection of the vector distance onto the reference vector
          ReferenceVector = BasisVector(iAxis,:) ! Reference vector in the plane
        END IF
      END DO
      VectorDistance = VectorDistance / DSQRT( DOT_PRODUCT( VectorDistance, VectorDistance ) ) ! Normalization
      ReferenceVector = ReferenceVector / DSQRT( DOT_PRODUCT( ReferenceVector, ReferenceVector ) ) ! Normalization
      CALL Cross_Product( ReferenceVector, VectorDistance, OrthogonalVector ) ! Orthogonal vector to the reference vector and the vector distance
      OrthogonalVector = OrthogonalVector / DSQRT( DOT_PRODUCT( OrthogonalVector, OrthogonalVector ) ) ! Normalization
      ! Signed angle between the vector distance and the reference vector
      Angle = DOT_PRODUCT( OrthogonalVector, NormalVector) / DABS( DOT_PRODUCT( OrthogonalVector, NormalVector) ) * &
      &       DACOS( DOT_PRODUCT( VectorDistance, ReferenceVector ) / ( &
      &       DSQRT( DOT_PRODUCT( VectorDistance, VectorDistance ) ) * DSQRT( DOT_PRODUCT( ReferenceVector, ReferenceVector ) ) ) )
      IF( Angle > 0.0D0 ) xChiralPositive = xChiralPositive + 1 ! Positive chirality (clockwise)
      IF( Angle < 0.0D0 ) xChiralNegative = xChiralNegative + 1 ! Negative chirality (counterclockwise)
    END DO
  END IF
  ! Accumulate the chirality
  ChiralityAccumulator = ChiralityAccumulator + DBLE( xChiralPositive - xChiralNegative ) / &
  &                      DBLE( COUNT( ClosestNeighbors(iParticle,:) /= 0, Dim= 1 ) ) ! Chirality
  ! Chirality of the particle
  ChiralityParticle(iParticle) = DBLE( xChiralPositive - xChiralNegative ) / DBLE( COUNT( &
  &                              ClosestNeighbors(iParticle,:) /= 0, Dim= 1 ) ) ! Chirality of the particle
  ! Accumulate the trigonometric functions
  IF( NeighborOption2D == 1 ) THEN ! ! Voronoi tessellation
    CossineTerm = CossineTerm + Cossine ! Sum of cossine terms
    SineTerm = SineTerm + Sine ! Sum of sine terms
  ELSE IF( NeighborOption2D == 2 ) THEN ! ! Voronoi tessellation
    CossineTerm = CossineTerm + Cossine / 2.D0 / cPi ! Sum of cossine terms
    SineTerm = SineTerm + Sine / 2.D0 / cPi ! Sum of sine terms
  END IF
END DO
!$OMP END PARALLEL DO

! Local bond parameters
Psi = DSQRT( CossineTerm * CossineTerm + SineTerm * SineTerm ) / DBLE( nParticles )

! Chiral order parameter
ChiralOrder = ChiralityAccumulator / DBLE( nParticles ) ! Chiral order parameter

RETURN

END SUBROUTINE LocalBondParameters

! ********************************************************************************************************************************* !
! This subroutine sorts an array in ascending order.                                                                                !
! ********************************************************************************************************************************* !
SUBROUTINE Sort_Array( Indices, Values )

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 ) :: i, j, TempIndex

! INTEGER VARIABLES (ARRAY)
INTEGER( Kind= Int64 ), DIMENSION( : ), INTENT(INOUT) :: Indices

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( : ), INTENT(IN) :: Values

! Sort the array in ascending order
DO i = 1, SIZE( Indices ) - 1
  DO j = i + 1, SIZE( Indices )
    IF( Values(Indices(i)) > Values(Indices(j)) ) THEN
      TempIndex = Indices(i)
      Indices(i) = Indices(j)
      Indices(j) = TempIndex
    END IF
  END DO
END DO

RETURN

END SUBROUTINE Sort_Array

! ********************************************************************************************************************************* !
! This subroutine calculates the Steinhardt order parameters (<SteinhardtParameter>) of a configuration. The parameter is           !
! calculated for each particle. If the option <ComputeWignerParameter> is enabled, the Wigner order parameters (<WignerParameter>)  !
! are also calculated.                                                                                                              !
!                                                                                                                                   !
! ● The subroutine starts by finding the neighbours of each particle in three-dimensions. The neighbours can be identified using:   !
!   [1] the nearest neighbours, [2] a cut-off distance, and [3] the Voronoi tessellation. If option [3] is chosen, the subtended    !
!   solid angles associated with the faces (neighbours) of the Voronoi cells (<SolidAngleNeighbor>) are also computed using         !
!   L'Huilier's rule. This is done to weight the qℓm(i) terms in the calculation of the Steinhardt order parameters.                !
!                                                                                                                                   !
! ● The Steinhardt order parameters are defined as:                                                                                 !
!                                                                                                                                   !
!     Qℓ(i) = (4π / (2ℓ + 1)) × Σₘ | qℓm(i) |²,                                                                                     !
!                                                                                                                                   !
!   where m ranges from -ℓ to ℓ, and qℓm(i) is defined as:                                                                          !
!                                                                                                                                   !
!     qℓm(i) = 1 / N(i) × Σⱼ Yℓm(θᵢⱼ, φᵢⱼ),                                                                                         !
!                                                                                                                                   !
!   where N(i) is the number of neighbours of particle i and Yℓm(θᵢⱼ, φᵢⱼ) is the spherical harmonics associated with the polar     !
!   and azimuthal angles of the vector distance between particles i and j (spherical coordinates).                                  !
!                                                                                                                                   !
! ● If enabled, the Wigner order parameters are computed as:                                                                        !
!                                                                                                                                   !
!     wℓ(i) = Σm₁+m₂+m₃=0 [ W3j(ℓ, ℓ, ℓ, m₁, m₂, m₃) × qℓm₁(i) × qℓm₂(i) × qℓm₃(i) ] / √[ Σₘ | qℓm(i) |² ]³                         !
!                                                                                                                                   !
!   where "W3j(ℓ, ℓ, ℓ, m₁, m₂, m₃)" is the Wigner 3-j symbol, and Σm₁+m₂+m₃=0 is a summation over the magnetic quantum numbers m₁, !
!   m₂, and m₃ in such a way that m₁ + m₂ + m₃ = 0 (angular momentum conservation). The product qℓm₁(i) × qℓm₂(i) × qℓm₃(i) is a    !
!   multiplication of complex numbers. To satisfy the angular momentum conservation, the imaginary part of the product needs to be  !
!   equal to zero (or very close to it).                                                                                            !
!                                                                                                                                   !
! ● For a better behavior of the order parameters (less scattering), consider enabling the option <AveragedSteinhardt>. This option !
!   includes the neighbours of the neighbours of the central particle in the calculation.                                           !
!                                                                                                                                   !
! For a better performance, this subroutine is parallelized using OpenMP.                                                           !
! ********************************************************************************************************************************* !
SUBROUTINE Steinhardt_Parameters( CurrentFrame, TotalFrames, SteinhardtParameter, WignerParameter, gSteinhardtParameter, &
&                                 gWignerParameter )

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 )               :: m1Counter, m2Counter, m3Counter          ! Counter (magnetic quantum numbers)
INTEGER( Kind= Int64 )               :: iParticle, jParticle, kParticle          ! Counters (particles)
INTEGER( Kind= Int64 )               :: iNeighbor, jNeighbor, nNeighborNeighbors ! Counters (neighbors)
INTEGER( Kind= Int64 )               :: iTerm                                    ! Counters (Steinhardt order parameters)
INTEGER( Kind= Int64 )               :: nCandidatesVoronoi                       ! Number of candidates (Voronoi tessellation)
INTEGER( Kind= Int64 )               :: iTriangulation                           ! Number of candidates (Voronoi tessellation)
INTEGER( Kind= Int64 )               :: iVertex                                  ! Counters (vertices)
INTEGER( Kind= Int64 )               :: fVertex                                  ! Counters (vertices of a face)
INTEGER( Kind= Int64 )               :: vNeighbor                                ! Counters (vertex shared by a neighbor)
INTEGER( Kind= Int64 )               :: iCandidate                               ! Counters (candidates)
INTEGER( Kind= Int64 )               :: mQuantumNumber                           ! Magnetic quantum number (m)
INTEGER( Kind= Int64 )               :: lQuantumNumber                           ! Azimuthal quantum number (l)
INTEGER( Kind= Int64 )               :: minInteger, maxInteger                   ! Minimum and maximum integer values
INTEGER( Kind= Int64 )               :: zInteger                                 ! Integer numbers
INTEGER( Kind= Int64 )               :: SignInteger, SignExponent                ! Sign of the integer summation
INTEGER( Kind= Int64 )               :: TotalNumberLoops                         ! Total number of loops in the parallel region
INTEGER( Kind= Int64 )               :: Progress                                 ! Progress in the parallel region
INTEGER( Kind= Int64 )               :: OldLenMessage, MaxLenMessage             ! Length of the print message
INTEGER( Kind= Int64 )               :: nVertices, nFaces, nEdges                ! Number of vertices, faces, and edges
INTEGER( Kind= Int64 )               :: VoronoiProgress                          ! Progress in the Voronoi tessellation
INTEGER( Kind= Int64 )               :: NextVertex                               ! Next vertex in the face
INTEGER( Kind= Int64 )               :: NextVertexNeighbor                       ! Index of the indices of the vertices of the face (neighbor)
INTEGER( Kind= Int64 )               :: LastVertex                               ! Last vertex of the triangle (fan triangulation)
INTEGER( Kind= Int64 ), INTENT( IN ) :: CurrentFrame, TotalFrames                ! Frames of the configuration file

! INTEGER VARIABLES (ARRAY)
INTEGER( Kind= Int64 ), DIMENSION( nParticles )     :: nFaceEdges         ! Number of edges per face
INTEGER( Kind= Int64 ), DIMENSION( nParticles )     :: ParticleIndices    ! Particle indices
INTEGER( Kind= Int64 ), DIMENSION( MaxVertices, 3 ) :: vIndex             ! Index of the vertices
INTEGER( Kind= Int64 ), DIMENSION( MaxVertices )    :: fVertexIndices     ! Index of the vertices of a face (neighbor)
INTEGER( Kind= Int64 ), DIMENSION( MaxVertices )    :: VertexOrder        ! Order of the vertices of a face (neighbor)

! INTEGER VARIABLES (ARRAY,ALLOCATABLE)
INTEGER( Kind= Int64 ), DIMENSION( nParticles, nParticles ) :: ClosestNeighbors ! Closest neighbors of each particle

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 ) :: RealTerm              ! Real part of the spherical harmonics
REAL( Kind= Real64 ) :: ImaginaryTerm         ! Imaginary part of the spherical harmonics
REAL( Kind= Real64 ) :: RadialDistance        ! Radial distance
REAL( Kind= Real64 ) :: PolarAngle            ! Polar angle
REAL( Kind= Real64 ) :: AzimuthalAngle        ! Azimuthal angle
REAL( Kind= Real64 ) :: Prefactor             ! Prefactor of the rotationally invariant order parameters
REAL( Kind= Real64 ) :: WignerNormalization   ! Wigner normalization factor (denominator)
REAL( Kind= Real64 ) :: Wigner3j              ! Wigner 3-j symbol
REAL( Kind= Real64 ) :: ProductLocalReal      ! Product of the local bond parameters (real part)
REAL( Kind= Real64 ) :: ProductLocalImaginary ! Product of the local bond parameters (imaginary part)
REAL( Kind= Real64 ) :: Wigner3jInteger       ! Wigner 3-j symbol per integer
REAL( Kind= Real64 ) :: Semiperimeter         ! Semiperimeter of the spherical triangle
REAL( Kind= Real64 ) :: SubtendedSolidAngle   ! Solid angle of the spherical triangle

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 3 )                                            :: FactorialTerm           ! Factorial term
REAL( Kind= Real64 ), DIMENSION( 3 )                                            :: VectorDistance          ! Vector distance
REAL( Kind= Real64 ), DIMENSION( 3 )                                            :: ScalingDistanceUnitBox  ! Scaled distance in the unit box
REAL( Kind= Real64 ), DIMENSION( 3 )                                            :: ArcLength               ! Arc length of the side of the subtended triangle
REAL( Kind= Real64 ), DIMENSION( 3, 3 )                                         :: TriangleDistances       ! Relative distances of the triangle vertices to the center of the Voronoi cell
REAL( Kind= Real64 ), DIMENSION( nParticles )                                   :: SquaredDistance         ! Squared distance
REAL( Kind= Real64 ), DIMENSION( MaxVertices, 3 )                               :: vRelativePosition       ! Relative position of the vertices
REAL( Kind= Real64 ), DIMENSION( nParticles, nParticles )                       :: SolidAngleNeighbor      ! Subtended solid angle of the neighbors (faces of the Voronoi cell)
REAL( Kind= Real64 ), DIMENSION( : ), INTENT( OUT )                             :: gSteinhardtParameter    ! Global Steinhardt parameters (Qℓ)
REAL( Kind= Real64 ), DIMENSION( :, : ), INTENT( OUT )                          :: gWignerParameter        ! Global Wigner parameters (wℓ)
REAL( Kind= Real64 ), DIMENSION( :, : ), INTENT( OUT )                          :: SteinhardtParameter     ! Steinhardt parameters (Qℓ)
REAL( Kind= Real64 ), DIMENSION( :, :, : ), INTENT( OUT )                       :: WignerParameter         ! Wigner parameters (wℓ)
REAL( Kind= Real64 ), DIMENSION( MaximumSteinhardt - MinimumSteinhardt + 1 )    :: Cossine                 ! Cossine term
REAL( Kind= Real64 ), DIMENSION( MaximumSteinhardt - MinimumSteinhardt + 1 )    :: CossineNeighbor         ! Cossine term of the neighbors
REAL( Kind= Real64 ), DIMENSION( MaximumSteinhardt - MinimumSteinhardt + 1 )    :: Sine                    ! Sine term
REAL( Kind= Real64 ), DIMENSION( MaximumSteinhardt - MinimumSteinhardt + 1 )    :: SineNeighbor            ! Sine term of the neighbors
REAL( Kind= Real64 ), DIMENSION( MaximumSteinhardt - MinimumSteinhardt + 1 )    :: CossineTerm             ! Sum of cossine terms (absolute squared magnitude of the spherical harmonics)
REAL( Kind= Real64 ), DIMENSION( MaximumSteinhardt - MinimumSteinhardt + 1 )    :: SineTerm                ! Sum of sine terms (absolute squared magnitude of the spherical harmonics)
REAL( Kind= Real64 ), DIMENSION( MaximumSteinhardt - MinimumSteinhardt + 1, 2 ) :: WignerParameterParticle ! Wigner parameters (wℓ) per particle
REAL( Kind= Real64 ), DIMENSION( nParticles, MaximumSteinhardt - MinimumSteinhardt + 1, -MaximumSteinhardt:MaximumSteinhardt, 2 ) &
&                                                                               :: GlobalSteinhardt        ! Global Steinhardt order parameter (Qℓ)

! REAL VARIABLES (ARRAY,ALLOCATABLE)
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE       :: SquaredDistanceVoronoi ! Squared distance (Voronoi tessellation)
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE    :: cRelativeDistance      ! Relative distance of the candidates
REAL( Kind= Real64 ), DIMENSION( :, :, : ), ALLOCATABLE :: LocalBondWigner        ! Local bond parameters (qℓm[i]) for the calculation of the Wigner parameters

! LOGICAL VARIABLES
LOGICAL, DIMENSION( MaxVertices ) :: VertexChosen ! Flags the vertices that have already been chosen

! CHARACTER STRINGS
CHARACTER( Len= 1000 ) :: PrintStatus ! Print status

! Initialization
ClosestNeighbors = 0
SteinhardtParameter = 0.D0
WignerParameterParticle = 0.D0
SolidAngleNeighbor = 0.D0
VoronoiProgress = 0
GlobalSteinhardt = 0.D0

! Allocation
IF( ALLOCATED( SquaredDistanceVoronoi ) ) DEALLOCATE( SquaredDistanceVoronoi )
ALLOCATE( SquaredDistanceVoronoi( nParticles ) )
SquaredDistanceVoronoi = 0.D0
IF( ALLOCATED( cRelativeDistance ) ) DEALLOCATE( cRelativeDistance )
ALLOCATE( cRelativeDistance( nParticles, 3 ) )
cRelativeDistance = 0.D0

! Initialization
MaxLenMessage = 0

! Find the neighbors of each particle
!$OMP PARALLEL DO DEFAULT( Shared ) &
!$OMP PRIVATE( iParticle, jParticle, kParticle, VectorDistance, ScalingDistanceUnitBox, SquaredDistance, cRelativeDistance ) &
!$OMP PRIVATE( ParticleIndices, iNeighbor, nVertices, nEdges, nFaces, nFaceEdges, vRelativePosition, vIndex, NextVertex ) &
!$OMP PRIVATE( SquaredDistanceVoronoi, nCandidatesVoronoi, iVertex, iCandidate, fVertexIndices, fVertex, vNeighbor, VertexChosen ) &
!$OMP PRIVATE( VertexOrder, NextVertexNeighbor, iTriangulation, LastVertex, ArcLength, TriangleDistances ) &
!$OMP PRIVATE( Semiperimeter, SubtendedSolidAngle )
DO iParticle = 1, nParticles ! Loop over all particles
  ! Initialization
  SquaredDistance = HUGE( 0.D0 ) ! Large value
  SquaredDistanceVoronoi = HUGE( 0.D0 ) ! Large value
  DO jParticle = 1, nParticles ! Loop over all other particles
    IF( iParticle == jParticle ) CYCLE ! Skip the same particle
    VectorDistance = pPosition(jParticle,1:3) - pPosition(iParticle,1:3) ! Vector distance between particle i and j
    CALL MatrixVectorMultiplication( BoxLengthInverse, VectorDistance, ScalingDistanceUnitBox ) ! Spatial transformation
    ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox ) ! Periodic boundary conditions
    CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, VectorDistance ) ! Spatial transformation
    SquaredDistance(jParticle) = DOT_PRODUCT( VectorDistance, VectorDistance ) ! Squared distance
    IF( NeighborOption == 3 ) cRelativeDistance(jParticle,:) = VectorDistance
    IF( NeighborOption == 3 ) SquaredDistanceVoronoi(jParticle) = SquaredDistance(jParticle)
  END DO
  ParticleIndices = (/ ( kParticle, kParticle = 1, nParticles ) /) ! Create an array of particle indices from 1 to N
  IF( NeighborOption == 1 .OR. NeighborOption == 2 ) CALL Sort_Array( ParticleIndices, SquaredDistance ) ! Sort the indices of the particles according to the squared distance in ascending order
  IF( NeighborOption == 3 ) CALL Sort_Array_Voronoi( ParticleIndices, cRelativeDistance, SquaredDistanceVoronoi ) ! Sort the indices of the particles according to the squared distance in ascending order (Voronoi tessellation)
  ! Initialization of the closest neighbors
  IF( NeighborOption == 1 .OR. NeighborOption == 2 ) THEN ! Nearest neighbors [1] or Cut-off distance [2]
    ! Initialization
    iNeighbor = 1
    DO jParticle = 1, nParticles ! Loop over all other particles
      ! Add the closest neighbors to the list
      IF( NeighborOption == 1 ) THEN ! Nearest neighbors
        ClosestNeighbors(iParticle,iNeighbor) = ParticleIndices(jParticle) ! Closest neighbors
        iNeighbor = iNeighbor + 1 ! Increment the counter
        IF( iNeighbor > nNeighbors ) EXIT ! Maximum number of neighbors reached [OPTION 1]
      ELSE IF( NeighborOption == 2 ) THEN ! Cut-off distance
        IF( ParticleIndices(jParticle) == iParticle ) CYCLE ! Skip the same particle
        ! Vector distance between particle i and its nearest neighbors j
        VectorDistance = pPosition(ParticleIndices(jParticle),1:3) - pPosition(iParticle,1:3) ! Vector distance between particle i and its nearest neighbors j
        CALL MatrixVectorMultiplication( BoxLengthInverse, VectorDistance, ScalingDistanceUnitBox ) ! Spatial transformation
        ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox ) ! Periodic boundary conditions
        CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, VectorDistance ) ! Spatial transformation
        ! Check if the distance is within the cut-off distance
        IF( DOT_PRODUCT( VectorDistance, VectorDistance ) <= CutoffDistance * CutoffDistance ) THEN
          ClosestNeighbors(iParticle,iNeighbor) = ParticleIndices(jParticle) ! Closest neighbors within the cut-off distance [OPTION 2]
          iNeighbor = iNeighbor + 1 ! Increment the counter
        END IF
      END IF
    END DO
  ELSE IF( NeighborOption == 3 ) THEN ! Voronoi tessellation
    ! Number of candidates
    nCandidatesVoronoi = 0
    DO kParticle = 1, nParticles
      IF( SquaredDistanceVoronoi(kParticle) <= CutoffDistanceVoronoi * CutoffDistanceVoronoi ) THEN
        nCandidatesVoronoi = nCandidatesVoronoi + 1
      END IF
    END DO
    ! Voronoi tessellation
    CALL VoronoiConstruction3D( nCandidatesVoronoi, nVertices, nEdges, nFaces, cRelativeDistance, SquaredDistanceVoronoi, &
    &                           nFaceEdges, vRelativePosition, vIndex )
    ! Find the closest neighbors within the Voronoi tessellation
    iNeighbor = 1
    DO iCandidate = 1, nCandidatesVoronoi
      IF( nFaceEdges(iCandidate) /= 0 ) THEN
        ClosestNeighbors(iParticle,iNeighbor) = ParticleIndices(iCandidate)
        iNeighbor = iNeighbor + 1
      END IF
    END DO
    ! Locate the vertices of the Voronoi tessellation of every face (neighbor)
    DO iNeighbor = 1, COUNT( ClosestNeighbors(iParticle,:) /= 0, Dim= 1 ) ! Loop over the neighbors (faces)
      ! Initialization
      fVertexIndices = 0
      fVertex = 1
      DO iVertex = 1, nVertices ! Loop over all vertices of the Voronoi cell
        ! Check if the index of the closest neighbor is associated with the vertex
        IF( ANY( ParticleIndices(vIndex(iVertex,:)) == ClosestNeighbors(iParticle,iNeighbor) ) ) THEN
          ! Store the index of the vertex associated with the common face (vertices of a particular face)
          fVertexIndices(fVertex) = iVertex
          fVertex = fVertex + 1
        END IF
      END DO
      ! Find the order of the vertices of the face (edge connections)
      VertexChosen = .FALSE.
      iVertex = 1
      VertexOrder = 0
      NextVertexNeighbor = 1
      ! Pick first vertex to start the face construction
      DO vNeighbor = 1, 3
        IF( ParticleIndices(vIndex(fVertexIndices(iVertex),vNeighbor)) /= ClosestNeighbors(iParticle,iNeighbor) ) THEN ! Do not pick the index of the vertex related to the neighbor to be the next vertex
          VertexChosen(iVertex) = .TRUE. ! Mark this vertex of the face as chosen so it will not be picked again
          VertexOrder(NextVertexNeighbor) = fVertexIndices(iVertex) ! Start the face construction
          NextVertex = ParticleIndices(vIndex(fVertexIndices(iVertex),vNeighbor)) ! Next vertex to be picked (it must be the index of one of the other two neighbors that share this vertex)
          NextVertexNeighbor = NextVertexNeighbor + 1
          EXIT
        END IF
      END DO
      ! Loop through the vertex list of the Voronoi cell to pick the remaining vertices (face construction)
      DO WHILE( iVertex < COUNT( fVertexIndices(:) /= 0, Dim= 1 ) )
        iVertex = iVertex + 1
        IF( .NOT. VertexChosen(iVertex) ) THEN ! Check if the vertex has already been chosen
          DO vNeighbor = 1, 3 ! Loop through the indices of the vertices
            IF( ParticleIndices(vIndex(fVertexIndices(iVertex),vNeighbor)) /= ClosestNeighbors(iParticle,iNeighbor) ) THEN ! Do not pick the index of the vertex related to the neighbor to be the next vertex
              IF( ParticleIndices(vIndex(fVertexIndices(iVertex),vNeighbor)) == NextVertex ) THEN ! Check if the vertex is the next vertex to be picked
                VertexChosen(iVertex) = .TRUE. ! Mark this vertex of the face as chosen so it will not be picked again
                VertexOrder(NextVertexNeighbor) = fVertexIndices(iVertex)
                NextVertexNeighbor = NextVertexNeighbor + 1
                ! Select the next vertex to be picked
                IF( ParticleIndices(vIndex(fVertexIndices(iVertex),1)) == ClosestNeighbors(iParticle,iNeighbor) ) THEN
                  IF( ParticleIndices(vIndex(fVertexIndices(iVertex),2)) == NextVertex ) THEN
                    NextVertex = ParticleIndices(vIndex(fVertexIndices(iVertex),3))
                  ELSE IF( ParticleIndices(vIndex(fVertexIndices(iVertex),3)) == NextVertex ) THEN
                    NextVertex = ParticleIndices(vIndex(fVertexIndices(iVertex),2))
                  END IF
                ELSE IF( ParticleIndices(vIndex(fVertexIndices(iVertex),2)) == ClosestNeighbors(iParticle,iNeighbor) ) THEN
                  IF( ParticleIndices(vIndex(fVertexIndices(iVertex),1)) == NextVertex ) THEN
                    NextVertex = ParticleIndices(vIndex(fVertexIndices(iVertex),3))
                  ELSE IF( ParticleIndices(vIndex(fVertexIndices(iVertex),3)) == NextVertex ) THEN
                    NextVertex = ParticleIndices(vIndex(fVertexIndices(iVertex),1))
                  END IF
                ELSE IF( ParticleIndices(vIndex(fVertexIndices(iVertex),3)) == ClosestNeighbors(iParticle,iNeighbor) ) THEN
                  IF( ParticleIndices(vIndex(fVertexIndices(iVertex),1)) == NextVertex ) THEN
                    NextVertex = ParticleIndices(vIndex(fVertexIndices(iVertex),2))
                  ELSE IF( ParticleIndices(vIndex(fVertexIndices(iVertex),2)) == NextVertex ) THEN
                    NextVertex = ParticleIndices(vIndex(fVertexIndices(iVertex),1))
                  END IF
                END IF
                ! Restart the loop
                iVertex = 1
                ! Exit the loop if the next connector has been found
                EXIT
              END IF
            END IF
          END DO
        END IF
      END DO
      ! Triangulate the faces and calculate the subtended solid angles
      LastVertex = 2
      DO iTriangulation = 1, COUNT( VertexOrder(:) /= 0, Dim= 1 ) - 2 ! Number of times to triangulate the face
        ! Relative distances of the triangle vertices to the center of the Voronoi cell
        TriangleDistances(1,1:3) = vRelativePosition(VertexOrder(1),1:3)
        TriangleDistances(2,1:3) = vRelativePosition(VertexOrder(LastVertex),1:3)
        TriangleDistances(3,1:3) = vRelativePosition(VertexOrder(LastVertex+1),1:3)
        ! Calculate the arc length of the side of the subtended triangle
        ArcLength(1) = DOT_PRODUCT( TriangleDistances(2,1:3), TriangleDistances(3,1:3) ) / ( DSQRT( DOT_PRODUCT( &
        &     TriangleDistances(2,1:3), TriangleDistances(2,1:3) ) ) * DSQRT( DOT_PRODUCT( TriangleDistances(3,1:3), &
        &     TriangleDistances(3,1:3) ) ) )
        IF( ArcLength(1) > 1.D0 ) ArcLength(1) = 1.D0   ! Avoid numerical errors
        IF( ArcLength(1) < -1.D0 ) ArcLength(1) = -1.D0 ! Avoid numerical errors
        ArcLength(1) = DACOS( ArcLength(1) )
        ArcLength(2) = DOT_PRODUCT( TriangleDistances(1,1:3), TriangleDistances(3,1:3) ) / ( DSQRT( DOT_PRODUCT( &
        &     TriangleDistances(1,1:3), TriangleDistances(1,1:3) ) ) * DSQRT( DOT_PRODUCT( TriangleDistances(3,1:3), &
        &     TriangleDistances(3,1:3) ) ) )
        IF( ArcLength(2) > 1.D0 ) ArcLength(2) = 1.D0   ! Avoid numerical errors
        IF( ArcLength(2) < -1.D0 ) ArcLength(2) = -1.D0 ! Avoid numerical errors
        ArcLength(2) = DACOS( ArcLength(2) )
        ArcLength(3) = DOT_PRODUCT( TriangleDistances(1,1:3), TriangleDistances(2,1:3) ) / ( DSQRT( DOT_PRODUCT( &
        &     TriangleDistances(1,1:3), TriangleDistances(1,1:3) ) ) * DSQRT( DOT_PRODUCT( TriangleDistances(2,1:3), &
        &     TriangleDistances(2,1:3) ) ) )
        IF( ArcLength(3) > 1.D0 ) ArcLength(3) = 1.D0   ! Avoid numerical errors
        IF( ArcLength(3) < -1.D0 ) ArcLength(3) = -1.D0 ! Avoid numerical errors
        ArcLength(3) = DACOS( ArcLength(3) )
        ! Semi-perimeter of the spherical triangle
        Semiperimeter = 0.5D0 * SUM( ArcLength )
        ! L'Huilier's rule (also see Green, Geoexploration, 24, 61-69, 1986, https://doi.org/10.1016/0016-7142(86)90019-0 )
        SubtendedSolidAngle = DTAN( 0.5D0 * Semiperimeter ) * DTAN( 0.5D0 * ( Semiperimeter - ArcLength(1) ) ) * &
        &                     DTAN( 0.5D0 * ( Semiperimeter - ArcLength(2) ) ) * DTAN( 0.5D0 * ( Semiperimeter - ArcLength(3) ) )
        IF( SubtendedSolidAngle < 0.D0 ) SubtendedSolidAngle = 0.D0 ! Avoid numerical errors (sometimes it is a very small, negative number)
        SubtendedSolidAngle = 4.D0 * DATAN( DSQRT( SubtendedSolidAngle ) )
        ! For every triangle of the face, accumulate the solid angle subtended by the it
        SolidAngleNeighbor(iParticle,iNeighbor) = SolidAngleNeighbor(iParticle,iNeighbor) + SubtendedSolidAngle
        ! Loop over the triangles
        LastVertex = LastVertex + 1
      END DO
    END DO
    !$OMP ATOMIC
    VoronoiProgress = VoronoiProgress + 1
    !$OMP END ATOMIC
    ! Print the progress
    !$OMP CRITICAL
    ! Status
    WRITE( PrintStatus, "(5G0,F6.2,G0)" ) "Reading frame: ", CurrentFrame, " of ", TotalFrames, " | Computing Voronoi "// &
    &     "tesselation: ", 100.D0 * DBLE( MIN( VoronoiProgress, nParticles ) ) / DBLE( nParticles ), "%"
    OldLenMessage = LEN( TRIM( PrintStatus ) )
    IF( OldLenMessage > MaxLenMessage ) MaxLenMessage = OldLenMessage
    WRITE( *, "(2G0)", Advance= "No" ) CHAR(13), TRIM( PrintStatus )
    !$OMP END CRITICAL
  END IF
END DO
!$OMP END PARALLEL DO

! Initialization
TotalNumberLoops = nParticles * ( MaximumSteinhardt - MinimumSteinhardt + 1 )
Progress = 0

! Calculate the Steinhardt order parameters
!$OMP PARALLEL DO COLLAPSE( 2 ) DEFAULT( Shared ) &
!$OMP PRIVATE( iParticle, jParticle, iNeighbor, jNeighbor, nNeighborNeighbors, VectorDistance, ScalingDistanceUnitBox ) &
!$OMP PRIVATE( lQuantumNumber, mQuantumNumber, iTerm, Cossine, Sine, CossineTerm, SineTerm, CossineNeighbor, SineNeighbor ) &
!$OMP PRIVATE( RadialDistance, PolarAngle, AzimuthalAngle, RealTerm, ImaginaryTerm, Prefactor ) &
!$OMP PRIVATE( LocalBondWigner, m1Counter, m2Counter, m3Counter, Wigner3j, ProductLocalReal, ProductLocalImaginary ) &
!$OMP PRIVATE( WignerNormalization, WignerParameterParticle, FactorialTerm, minInteger, maxInteger, zInteger, SignInteger ) &
!$OMP PRIVATE( SignExponent, Wigner3jInteger )
DO iParticle = 1, nParticles ! Loop over all particles
  DO iTerm = 1, MaximumSteinhardt - MinimumSteinhardt + 1 ! Loop over the Steinhardt order parameters (l)
    ! Initialization
    lQuantumNumber = MinimumSteinhardt + iTerm - 1 ! Azimuthal quantum number (l)
    CossineTerm = 0.D0
    SineTerm = 0.D0
    ! Allocation
    IF( ALLOCATED( LocalBondWigner ) ) DEALLOCATE( LocalBondWigner )
    ALLOCATE( LocalBondWigner( iTerm, -lQuantumNumber:lQuantumNumber, 2 ) )
    LocalBondWigner = 0.D0
    ! Loop over the magnetic quantum numbers (m)
    DO mQuantumNumber = - lQuantumNumber, lQuantumNumber
      ! Initialization
      Cossine = 0.D0
      Sine = 0.D0
      DO iNeighbor = 1, COUNT( ClosestNeighbors(iParticle,:) /= 0, Dim= 1 ) ! Loop over closest neighbors of particle i
        VectorDistance = pPosition(ClosestNeighbors(iParticle,iNeighbor),1:3) - pPosition(iParticle,1:3) ! Vector distance between particle i and its closest neighbors j
        CALL MatrixVectorMultiplication( BoxLengthInverse, VectorDistance, ScalingDistanceUnitBox ) ! Spatial transformation
        ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox ) ! Periodic boundary conditions
        CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, VectorDistance ) ! Spatial transformation
        CALL ConvertCartesianToSpherical( VectorDistance, RadialDistance, PolarAngle, AzimuthalAngle ) ! Convert the vector distance to spherical coordinates
        CALL Spherical_Harmonics( lQuantumNumber, mQuantumNumber, PolarAngle, AzimuthalAngle, RealTerm, ImaginaryTerm ) ! Calculate the spherical harmonics
        IF( NeighborOption == 3 ) THEN ! Voronoi tessellation [3]
          ! Calculate the average weighted by the solid angle subtended by the faces of the Voronoi cell
          RealTerm = RealTerm * SolidAngleNeighbor(iParticle,iNeighbor) ! Solid angle subtended by the face of the Voronoi cell
          ImaginaryTerm = ImaginaryTerm * SolidAngleNeighbor(iParticle,iNeighbor) ! Solid angle subtended by the face of the Voronoi cell
        END IF
        Cossine(iTerm) = Cossine(iTerm) + RealTerm ! Real part
        Sine(iTerm) = Sine(iTerm) + ImaginaryTerm ! Imaginary part
      END DO
      ! Average the trigonometric functions over the number of neighbors of particle i
      IF( NeighborOption == 1 .OR. NeighborOption == 2 ) THEN ! Nearest neighbors [1] or Cut-off distance [2]
        Cossine(iTerm) = Cossine(iTerm) / DBLE( COUNT( ClosestNeighbors(iParticle,:) /= 0, Dim= 1 ) ) ! Real part
        Sine(iTerm) = Sine(iTerm) / DBLE( COUNT( ClosestNeighbors(iParticle,:) /= 0, Dim= 1 ) ) ! Imaginary part
      ELSE IF( NeighborOption == 3 ) THEN ! Voronoi tessellation [3]
        Cossine(iTerm) = Cossine(iTerm) / 4.D0 / cPi ! Real part
        Sine(iTerm) = Sine(iTerm) / 4.D0 / cPi ! Imaginary part
      END IF
      ! Store the local bond parameters (qℓm[i]) for the calculation of the global parameters
      GlobalSteinhardt(iParticle,iTerm,mQuantumNumber,1) = Cossine(iTerm) ! Real part
      GlobalSteinhardt(iParticle,iTerm,mQuantumNumber,2) = Sine(iTerm) ! Imaginary part
      ! Averaged Steinhardt order parameters
      IF( AveragedSteinhardt ) THEN
        DO iNeighbor = 1, COUNT( ClosestNeighbors(iParticle,:) /= 0, Dim= 1 ) ! Loop over closest neighbors of particle i
          ! Initialization
          nNeighborNeighbors = COUNT( ClosestNeighbors(ClosestNeighbors(iParticle,iNeighbor),:) /= 0, Dim= 1 ) ! Number of neighbors of the neighbors
          CossineNeighbor(iTerm) = 0.D0
          SineNeighbor(iTerm) = 0.D0
          DO jNeighbor = 1, nNeighborNeighbors ! Loop over the neighbors of the neighbors of particle i
            ! Vector distance between the neighbors of the neighbors of particle i
            VectorDistance = pPosition(ClosestNeighbors(ClosestNeighbors(iParticle,iNeighbor),jNeighbor),1:3) - &
            &                pPosition(ClosestNeighbors(iParticle,iNeighbor),1:3)
            CALL MatrixVectorMultiplication( BoxLengthInverse, VectorDistance, ScalingDistanceUnitBox ) ! Spatial transformation
            ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox ) ! Periodic boundary conditions
            CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, VectorDistance ) ! Spatial transformation
            CALL ConvertCartesianToSpherical( VectorDistance, RadialDistance, PolarAngle, AzimuthalAngle ) ! Convert the vector distance to spherical coordinates
            CALL Spherical_Harmonics( lQuantumNumber, mQuantumNumber, PolarAngle, AzimuthalAngle, RealTerm, ImaginaryTerm ) ! Calculate the spherical harmonics
            IF( NeighborOption == 3 ) THEN ! Voronoi tessellation [3]
              ! Calculate the average weighted by the solid angle subtended by the faces of the Voronoi cell
              RealTerm = RealTerm * SolidAngleNeighbor(ClosestNeighbors(iParticle,iNeighbor),jNeighbor) ! Solid angle subtended by the face of the Voronoi cell
              ImaginaryTerm = ImaginaryTerm * SolidAngleNeighbor(ClosestNeighbors(iParticle,iNeighbor),jNeighbor) ! Solid angle subtended by the face of the Voronoi cell
            END IF
            CossineNeighbor(iTerm) = CossineNeighbor(iTerm) + RealTerm ! Real part
            SineNeighbor(iTerm) = SineNeighbor(iTerm) + ImaginaryTerm ! Imaginary part
          END DO
          ! Average the trigonometric functions over the number of neighbors of the neighbors of particle i
          IF( NeighborOption == 1 .OR. NeighborOption == 2 ) THEN ! Nearest neighbors [1] or Cut-off distance [2]
            CossineNeighbor(iTerm) = CossineNeighbor(iTerm) / DBLE( nNeighborNeighbors ) ! Real part
            SineNeighbor(iTerm) = SineNeighbor(iTerm) / DBLE( nNeighborNeighbors ) ! Imaginary part
          ELSE IF( NeighborOption == 3 ) THEN ! Voronoi tessellation [3]
            CossineNeighbor(iTerm) = CossineNeighbor(iTerm) / 4.D0 / cPi ! Real part
            SineNeighbor(iTerm) = SineNeighbor(iTerm) / 4.D0 / cPi ! Imaginary part
          END IF
          ! Accumulate the trigonometric functions
          Cossine(iTerm) = Cossine(iTerm) + CossineNeighbor(iTerm) ! Real part
          Sine(iTerm) = Sine(iTerm) + SineNeighbor(iTerm) ! Imaginary part
        END DO
        ! Average the trigonometric functions over the neighbors of particle i
        Cossine(iTerm) = Cossine(iTerm) / DBLE( 1 + COUNT( ClosestNeighbors(iParticle,:) /= 0, Dim= 1 ) ) ! Real part
        Sine(iTerm) = Sine(iTerm) / DBLE( 1 + COUNT( ClosestNeighbors(iParticle,:) /= 0, Dim= 1 ) ) ! Imaginary part
      END IF
      ! Save the trigonometric functions if the Wigner parameters are computed
      IF( ComputeWignerParameter ) THEN
        LocalBondWigner(iTerm,mQuantumNumber,1) = Cossine(iTerm) ! Real part
        LocalBondWigner(iTerm,mQuantumNumber,2) = Sine(iTerm) ! Imaginary part
      END IF
      ! Absolute squared magnitude of the spherical harmonics (product between the complex number, a + i×b, and its complex conjugate, a - i×b, yielding: a² + b²)
      CossineTerm(iTerm) = CossineTerm(iTerm) + Cossine(iTerm) * Cossine(iTerm)
      SineTerm(iTerm) = SineTerm(iTerm) + Sine(iTerm) * Sine(iTerm)
    END DO
    ! Compute Wigner's parameters
    IF( ComputeWignerParameter ) THEN
      ! Initialization
      WignerParameterParticle(iTerm,:) = 0.D0
      ! Sum over the magnetic quantum numbers
      DO m1Counter = - lQuantumNumber, lQuantumNumber
        DO m2Counter = - lQuantumNumber, lQuantumNumber
          DO m3Counter = - lQuantumNumber, lQuantumNumber
            IF( m1Counter + m2Counter + m3Counter /= 0 ) CYCLE ! Condition to satisfy the angular momentum conservation
            ! Initialization
            FactorialTerm = 0.D0
            ! First factorial term
            FactorialTerm(1) = Factorial( lQuantumNumber ) * Factorial( lQuantumNumber ) * Factorial( lQuantumNumber )
            FactorialTerm(1) = FactorialTerm(1) / Factorial( 3 * lQuantumNumber + 1 )
            FactorialTerm(1) = DSQRT( FactorialTerm(1) )
            ! Second factorial term
            FactorialTerm(2) = Factorial( lQuantumNumber + m1Counter ) * Factorial( lQuantumNumber - m1Counter ) * &
            &                  Factorial( lQuantumNumber + m2Counter ) * Factorial( lQuantumNumber - m2Counter ) * &
            &                  Factorial( lQuantumNumber + m3Counter ) * Factorial( lQuantumNumber - m3Counter )
            FactorialTerm(2) = DSQRT( FactorialTerm(2) )
            ! Condition to calculate the sum over integer values
            minInteger = MAX( 0, m2Counter, - m1Counter )
            maxInteger = MIN( lQuantumNumber - m1Counter, lQuantumNumber + m2Counter, lQuantumNumber )
            ! Condition of only positive factorial terms (z ≥ 0)
            IF( maxInteger >= minInteger ) THEN
              DO zInteger = minInteger, maxInteger
                SignExponent = zInteger - m3Counter
                SignInteger = MERGE( 1, - 1, MOD( SignExponent, 2 ) == 0 .OR. SignExponent == 0 )
                Wigner3jInteger = Factorial( zInteger ) * Factorial( lQuantumNumber - zInteger ) * Factorial( lQuantumNumber - &
                &                 m1Counter - zInteger ) * Factorial( lQuantumNumber + m2Counter - zInteger ) * &
                &                 Factorial( m1Counter + zInteger ) * Factorial( - m2Counter + zInteger )
                Wigner3jInteger = 1.D0 / Wigner3jInteger
                ! Third factorial term
                FactorialTerm(3) = FactorialTerm(3) + DSIGN( Wigner3jInteger, DBLE( SignInteger ) )
              END DO
            END IF
            ! Wigner 3-j symbol (See Eq. 106.14 (pp. 437-438) of Quantum Mechanics by Landau and Lifshitz, 3rd Edition, 1977)
            Wigner3j = FactorialTerm(1) * FactorialTerm(2) * FactorialTerm(3)
            ! Product of the local bond parameters
            ProductLocalReal = LocalBondWigner(iTerm,m1Counter,1) * LocalBondWigner(iTerm,m2Counter,1) * &
            &     LocalBondWigner(iTerm,m3Counter,1) - LocalBondWigner(iTerm,m1Counter,2) * LocalBondWigner(iTerm,m2Counter,2) * &
            &     LocalBondWigner(iTerm,m3Counter,1) - LocalBondWigner(iTerm,m1Counter,1) * LocalBondWigner(iTerm,m2Counter,2) * &
            &     LocalBondWigner(iTerm,m3Counter,2) - LocalBondWigner(iTerm,m2Counter,1) * LocalBondWigner(iTerm,m1Counter,2) * &
            &     LocalBondWigner(iTerm,m3Counter,2)
            ProductLocalImaginary = LocalBondWigner(iTerm,m1Counter,1) * LocalBondWigner(iTerm,m3Counter,1) * &
            &     LocalBondWigner(iTerm,m2Counter,2) + LocalBondWigner(iTerm,m2Counter,1) * LocalBondWigner(iTerm,m3Counter,1) * &
            &     LocalBondWigner(iTerm,m1Counter,2) + LocalBondWigner(iTerm,m1Counter,1) * LocalBondWigner(iTerm,m2Counter,1) * &
            &     LocalBondWigner(iTerm,m3Counter,2) - LocalBondWigner(iTerm,m1Counter,2) * LocalBondWigner(iTerm,m2Counter,2) * &
            &     LocalBondWigner(iTerm,m3Counter,2) ! This should be 0 (or very close to 0)
            ! Accumulate the Wigner parameters
            WignerParameterParticle(iTerm,1) = WignerParameterParticle(iTerm,1) + Wigner3j * ProductLocalReal
            WignerParameterParticle(iTerm,2) = WignerParameterParticle(iTerm,2) + Wigner3j * ProductLocalImaginary
          END DO
        END DO
      END DO
      WignerNormalization = DSQRT( CossineTerm(iTerm) + SineTerm(iTerm) ) * DSQRT( CossineTerm(iTerm) + SineTerm(iTerm) ) * &
      &                     DSQRT( CossineTerm(iTerm) + SineTerm(iTerm) )
      ! Normalize the Wigner parameters
      WignerParameter(iParticle,iTerm,1) = WignerParameterParticle(iTerm,1) / WignerNormalization
      WignerParameter(iParticle,iTerm,2) = WignerParameterParticle(iTerm,2) / WignerNormalization
    END IF
    ! Steinhardt parameters (rotationally invariant order parameters)
    Prefactor = 4.D0 * cPi / DBLE( 2 * lQuantumNumber + 1 ) ! Prefactor
    SteinhardtParameter(iParticle,iTerm) = DSQRT( Prefactor * ( CossineTerm(iTerm) + SineTerm(iTerm) ) )
    ! Update counter
    !$OMP ATOMIC
    Progress = Progress + 1
    !$OMP END ATOMIC
    ! Print the progress
    !$OMP CRITICAL
    ! Status
    WRITE( PrintStatus, "(5G0,F6.2,G0)" ) "Reading frame: ", CurrentFrame, " of ", TotalFrames, " | Calculating Steinhardt "// &
    &     "parameters: ", 100.D0 * DBLE( MIN( Progress * nParticles / TotalNumberLoops, nParticles ) ) / DBLE( nParticles ), "%"
    OldLenMessage = LEN( TRIM( PrintStatus ) )
    IF( OldLenMessage > MaxLenMessage ) MaxLenMessage = OldLenMessage
    WRITE( *, "(2G0)", Advance= "No" ) CHAR(13), TRIM( PrintStatus )
    !$OMP END CRITICAL
  END DO
END DO
!$OMP END PARALLEL DO

! Global parameters
DO iTerm = 1, MaximumSteinhardt - MinimumSteinhardt + 1
  ! Initialization
  lQuantumNumber = MinimumSteinhardt + iTerm - 1 ! Azimuthal quantum number (l)
  CossineTerm = 0.D0
  SineTerm = 0.D0
  ! Allocation
  IF( ALLOCATED( LocalBondWigner ) ) DEALLOCATE( LocalBondWigner )
  ALLOCATE( LocalBondWigner( iTerm, -lQuantumNumber:lQuantumNumber, 2 ) )
  LocalBondWigner = 0.D0
  ! Loop over the magnetic quantum numbers (m)
  DO mQuantumNumber = - lQuantumNumber, lQuantumNumber
    ! Initialization
    Cossine = 0.D0
    Sine = 0.D0
    ! Loop over the particles
    DO iParticle = 1, nParticles
      IF( NeighborOption == 1 .OR. NeighborOption == 2 ) THEN ! Nearest neighbors [1] or Cut-off distance [2]
        Cossine(iTerm) = Cossine(iTerm) + DBLE( COUNT( ClosestNeighbors(iParticle,:) /= 0, Dim= 1 ) ) * &
        &     GlobalSteinhardt(iParticle,iTerm,mQuantumNumber,1)
        Sine(iTerm) = Sine(iTerm) + DBLE( COUNT( ClosestNeighbors(iParticle,:) /= 0, Dim= 1 ) ) * &
        &     GlobalSteinhardt(iParticle,iTerm,mQuantumNumber,2)
      ELSE IF( NeighborOption == 3 ) THEN ! Voronoi tessellation [3]
        Cossine(iTerm) = Cossine(iTerm) + GlobalSteinhardt(iParticle,iTerm,mQuantumNumber,1)
        Sine(iTerm) = Sine(iTerm) + GlobalSteinhardt(iParticle,iTerm,mQuantumNumber,2)
      END IF      
    END DO
    ! Average the trigonometric functions over the number of particles
    IF( NeighborOption == 1 .OR. NeighborOption == 2 ) THEN ! Nearest neighbors [1] or Cut-off distance [2]
      Cossine(iTerm) = Cossine(iTerm) / DBLE( COUNT( ClosestNeighbors /= 0 ) ) ! Real part
      Sine(iTerm) = Sine(iTerm) / DBLE( COUNT( ClosestNeighbors /= 0 ) ) ! Imaginary part
    ELSE IF( NeighborOption == 3 ) THEN ! Nearest neighbors [1] or Cut-off distance [2]
      Cossine(iTerm) = Cossine(iTerm) / nParticles ! Real part
      Sine(iTerm) = Sine(iTerm) / nParticles ! Imaginary part
    END IF
    ! Save the trigonometric functions if the Wigner parameters are computed
    IF( ComputeWignerParameter ) THEN
      LocalBondWigner(iTerm,mQuantumNumber,1) = Cossine(iTerm) ! Real part
      LocalBondWigner(iTerm,mQuantumNumber,2) = Sine(iTerm) ! Imaginary part
    END IF
    ! Absolute squared magnitude of the spherical harmonics (product between the complex number, a + i×b, and its complex conjugate, a - i×b, yielding: a² + b²)
    CossineTerm(iTerm) = CossineTerm(iTerm) + Cossine(iTerm) * Cossine(iTerm)
    SineTerm(iTerm) = SineTerm(iTerm) + Sine(iTerm) * Sine(iTerm)
  END DO
  ! Compute Wigner's parameters
  IF( ComputeWignerParameter ) THEN
    ! Initialization
    WignerParameterParticle(iTerm,:) = 0.D0
    ! Sum over the magnetic quantum numbers
    DO m1Counter = - lQuantumNumber, lQuantumNumber
      DO m2Counter = - lQuantumNumber, lQuantumNumber
        DO m3Counter = - lQuantumNumber, lQuantumNumber
          IF( m1Counter + m2Counter + m3Counter /= 0 ) CYCLE ! Condition to satisfy the angular momentum conservation
          ! Initialization
          FactorialTerm = 0.D0
          ! First factorial term
          FactorialTerm(1) = Factorial( lQuantumNumber ) * Factorial( lQuantumNumber ) * Factorial( lQuantumNumber )
          FactorialTerm(1) = FactorialTerm(1) / Factorial( 3 * lQuantumNumber + 1 )
          FactorialTerm(1) = DSQRT( FactorialTerm(1) )
          ! Second factorial term
          FactorialTerm(2) = Factorial( lQuantumNumber + m1Counter ) * Factorial( lQuantumNumber - m1Counter ) * &
          &                  Factorial( lQuantumNumber + m2Counter ) * Factorial( lQuantumNumber - m2Counter ) * &
          &                  Factorial( lQuantumNumber + m3Counter ) * Factorial( lQuantumNumber - m3Counter )
          FactorialTerm(2) = DSQRT( FactorialTerm(2) )
          ! Condition to calculate the sum over integer values
          minInteger = MAX( 0, m2Counter, - m1Counter )
          maxInteger = MIN( lQuantumNumber - m1Counter, lQuantumNumber + m2Counter, lQuantumNumber )
          ! Condition of only positive factorial terms (z ≥ 0)
          IF( maxInteger >= minInteger ) THEN
            DO zInteger = minInteger, maxInteger
              SignExponent = zInteger - m3Counter
              SignInteger = MERGE( 1, - 1, MOD( SignExponent, 2 ) == 0 .OR. SignExponent == 0 )
              Wigner3jInteger = Factorial( zInteger ) * Factorial( lQuantumNumber - zInteger ) * Factorial( lQuantumNumber - &
              &                 m1Counter - zInteger ) * Factorial( lQuantumNumber + m2Counter - zInteger ) * &
              &                 Factorial( m1Counter + zInteger ) * Factorial( - m2Counter + zInteger )
              Wigner3jInteger = 1.D0 / Wigner3jInteger
              ! Third factorial term
              FactorialTerm(3) = FactorialTerm(3) + DSIGN( Wigner3jInteger, DBLE( SignInteger ) )
            END DO
          END IF
          ! Wigner 3-j symbol (See Eq. 106.14 (pp. 437-438) of Quantum Mechanics by Landau and Lifshitz, 3rd Edition, 1977)
          Wigner3j = FactorialTerm(1) * FactorialTerm(2) * FactorialTerm(3)
          ! Product of the local bond parameters
          ProductLocalReal = LocalBondWigner(iTerm,m1Counter,1) * LocalBondWigner(iTerm,m2Counter,1) * &
          &     LocalBondWigner(iTerm,m3Counter,1) - LocalBondWigner(iTerm,m1Counter,2) * LocalBondWigner(iTerm,m2Counter,2) * &
          &     LocalBondWigner(iTerm,m3Counter,1) - LocalBondWigner(iTerm,m1Counter,1) * LocalBondWigner(iTerm,m2Counter,2) * &
          &     LocalBondWigner(iTerm,m3Counter,2) - LocalBondWigner(iTerm,m2Counter,1) * LocalBondWigner(iTerm,m1Counter,2) * &
          &     LocalBondWigner(iTerm,m3Counter,2)
          ProductLocalImaginary = LocalBondWigner(iTerm,m1Counter,1) * LocalBondWigner(iTerm,m3Counter,1) * &
          &     LocalBondWigner(iTerm,m2Counter,2) + LocalBondWigner(iTerm,m2Counter,1) * LocalBondWigner(iTerm,m3Counter,1) * &
          &     LocalBondWigner(iTerm,m1Counter,2) + LocalBondWigner(iTerm,m1Counter,1) * LocalBondWigner(iTerm,m2Counter,1) * &
          &     LocalBondWigner(iTerm,m3Counter,2) - LocalBondWigner(iTerm,m1Counter,2) * LocalBondWigner(iTerm,m2Counter,2) * &
          &     LocalBondWigner(iTerm,m3Counter,2) ! This should be 0 (or very close to 0)
          ! Accumulate the Wigner parameters
          WignerParameterParticle(iTerm,1) = WignerParameterParticle(iTerm,1) + Wigner3j * ProductLocalReal
          WignerParameterParticle(iTerm,2) = WignerParameterParticle(iTerm,2) + Wigner3j * ProductLocalImaginary
        END DO
      END DO
    END DO
    WignerNormalization = DSQRT( CossineTerm(iTerm) + SineTerm(iTerm) ) * DSQRT( CossineTerm(iTerm) + SineTerm(iTerm) ) * &
    &                     DSQRT( CossineTerm(iTerm) + SineTerm(iTerm) )
    ! Normalize the Wigner parameters
    gWignerParameter(iTerm,1) = WignerParameterParticle(iTerm,1) / WignerNormalization
    gWignerParameter(iTerm,2) = WignerParameterParticle(iTerm,2) / WignerNormalization
  END IF
  ! Steinhardt parameters (rotationally invariant order parameters)
  Prefactor = 4.D0 * cPi / DBLE( 2 * lQuantumNumber + 1 ) ! Prefactor
  gSteinhardtParameter(iTerm) = DSQRT( Prefactor * ( CossineTerm(iTerm) + SineTerm(iTerm) ) )
END DO

! Status
WRITE( PrintStatus, "(4G0)" ) "Reading frame: ", CurrentFrame, " of ", TotalFrames
OldLenMessage = LEN( TRIM( PrintStatus ) )
WRITE( *, "(6G0)", Advance= "No" ) CHAR(13), "Reading frame: ", CurrentFrame, " of ", TotalFrames, &
&                                  REPEAT( " ", MaxLenMessage - OldLenMessage )

RETURN

END SUBROUTINE Steinhardt_Parameters

! ********************************************************************************************************************************* !
! This subroutine sorts an array in ascending order.                                                                                !
! ********************************************************************************************************************************* !
SUBROUTINE Sort_Array_Voronoi( Indices, VectorDistance, Values )

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 ) :: i, j, TempIndex

! INTEGER VARIABLES (ARRAY)
INTEGER( Kind= Int64 ), DIMENSION( : ), INTENT(INOUT) :: Indices

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 ) :: TempSquaredDistance

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 3 )                   :: TempVector
REAL( Kind= Real64 ), DIMENSION( nParticles )          :: SquaredDistance
REAL( Kind= Real64 ), DIMENSION( : ), INTENT(INOUT)    :: Values
REAL( Kind= Real64 ), DIMENSION( :, : ), INTENT(INOUT) :: VectorDistance

! Initialization
SquaredDistance = Values

! Sort the array in ascending order
DO i = 1, SIZE( Indices ) - 1
  DO j = i + 1, SIZE( Indices )
    IF( Values(Indices(i)) > Values(Indices(j)) ) THEN
      TempIndex = Indices(i)
      Indices(i) = Indices(j)
      Indices(j) = TempIndex
      TempVector = VectorDistance(i,:)
      VectorDistance(i,:) = VectorDistance(j,:)
      VectorDistance(j,:) = TempVector
      TempSquaredDistance = SquaredDistance(i)
      SquaredDistance(i) = SquaredDistance(j)
      SquaredDistance(j) = TempSquaredDistance      
    END IF
  END DO
END DO

! Update the values
Values = SquaredDistance

RETURN

END SUBROUTINE Sort_Array_Voronoi

END MODULE BondParameters