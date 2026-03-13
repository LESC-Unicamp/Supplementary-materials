MODULE OrderParameters

! Uses two modules: global variables and quaternions
USE GlobalVariables
USE Quaternions

CONTAINS

! ********************************************************************************************************************************* !
! This subroutine calculates the uniaxial nematic order parameter of a liquid crystal.                                              !
!                                                                                                                                   !
! ● The calculation involves the construction of the Q-tensor and the determination of the eigenvector (<PhaseDirector>) associated !
!   with the largest eigenvalue (<OrderParameterS2>).                                                                               !
!                                                                                                                                   !
! ● The Q-tensor is constructed using the orientation of the particles (<pCurrentOrientation>) and the second-order Legendre        !
!   polynomial.                                                                                                                     !
! ********************************************************************************************************************************* !
! For more information, please read: M. P. Allen, D. J. Tildesley, Oxford University Press, 2nd Edition (2017).                     !
!                                    O. K. Smith, Communications of the ACM, 4(4), 168 (1961), DOI: 10.1145/355578.366316.          !
! ********************************************************************************************************************************* !
SUBROUTINE NematicOrderParameter( pCurrentOrientation, OrderParameterS2, PhaseDirector )

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 ) :: Alpha, Beta ! Unit vector specifiers (î, ĵ, k̂)
INTEGER( Kind= Int64 ) :: Particle    ! Counter (particle)

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 )                :: KroneckerDelta   ! Kronecker delta
REAL( Kind= Real64 )                :: M                ! One-third of the trace (matrix)
REAL( Kind= Real64 )                :: HQ               ! One-half of the determinant (matrix)
REAL( Kind= Real64 )                :: P                ! Sum of squares of elements of a matrix
REAL( Kind= Real64 )                :: Phi              ! Angle
REAL( Kind= Real64 )                :: TestCondition    ! Condition for real eigenvalues
REAL( Kind= Real64 ), INTENT( OUT ) :: OrderParameterS2 ! Nematic order parameter

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 3 )                           :: Eigenvalues          ! Eigenvalues of the order tensor Q
REAL( Kind= Real64 ), DIMENSION( 3, 3 )                        :: pTensorQ             ! Order tensor Q of particle i (3 x 3 Matrix)
REAL( Kind= Real64 ), DIMENSION( 3, 3 )                        :: TensorQ              ! Order tensor Q of particle i (Average)
REAL( Kind= Real64 ), DIMENSION( 3 ), INTENT( OUT )            :: PhaseDirector        ! Phase director (eigenvector of the largest eigenvalue)
REAL( Kind= Real64 ), DIMENSION( nParticles, 3 ), INTENT( IN ) :: pCurrentOrientation  ! Orientation of particles

! Initialization
pTensorQ = 0.D0

! Q-tensor construction
DO Particle = 1, nParticles
  DO Alpha = 1, 3
    DO Beta = 1, 3
      ! Dyadic product (Kronecker delta)
      IF ( Alpha == Beta ) THEN
          KroneckerDelta = 1.D0
      ELSE
          KroneckerDelta = 0.D0
      END IF
      ! Second-order Legendre polynomial (P₂)
      pTensorQ(Alpha,Beta) = pTensorQ(Alpha,Beta) + ( 1.5D0 * pCurrentOrientation(Particle,Alpha) * &
      &                      pCurrentOrientation(Particle,Beta) ) - (0.5D0 * KroneckerDelta)
    END DO
  END DO
END DO

! Averaged Q-tensor
TensorQ = pTensorQ / DBLE( nParticles )

! *********************************************************************************************** !
!                             Eigenvalues of a Symmetric 3 x 3 Matrix                             !
! *********************************************************************************************** !
! See O. K. Smith, Communications of the ACM, 4(4), 168 (1961), DOI: 10.1145/355578.366316        !
! for more information.                                                                           !
! *********************************************************************************************** !

! One-third of the trace of symmetric Q matrix
M = TensorQ(1,1) + TensorQ(2,2) + TensorQ(3,3)
M = M / 3.D0

! Half of the determinant of (Q - mI), where I is the identity matrix
HQ = ( TensorQ(1,1) - M ) * ( TensorQ(2,2) - M ) * ( TensorQ(3,3) - M ) + ( TensorQ(1,2) * TensorQ(2,3) * TensorQ(3,1) ) + &
&    ( TensorQ(1,3) * TensorQ(3,2) * TensorQ(2,1) ) - ( TensorQ(1,3) * ( TensorQ(2,2) - M ) * TensorQ(3,1) ) - &
&    ( TensorQ(2,3) * TensorQ(3,2) * ( TensorQ(1,1) - M ) ) - ( ( TensorQ(3,3) - M ) * TensorQ(2,1) * TensorQ(1,2) )
HQ = 0.5D0 * HQ

! One-sixth of the sum of squares of elements of (Q - mI)
P = ( TensorQ(1,1) - M ) * ( TensorQ(1,1) - M ) + ( TensorQ(1,2) * TensorQ(1,2) ) + ( TensorQ(1,3) * TensorQ(1,3) ) + &
&   ( TensorQ(2,2) - M ) * ( TensorQ(2,2) - M ) + ( TensorQ(2,1) * TensorQ(2,1) ) + ( TensorQ(2,3) * TensorQ(2,3) ) + &
&   ( TensorQ(3,3) - M ) * ( TensorQ(3,3) - M ) + ( TensorQ(3,1) * TensorQ(3,1) ) + ( TensorQ(3,2) * TensorQ(3,2) )
P = P / 6.D0

! Test condition
TestCondition = ( P * P * P ) - ( HQ * HQ )

! Real eigenvalues condition (P³ ≥ HQ²)
IF ( TestCondition >= 0.D0 ) THEN
  Phi = DATAN( DSQRT( TestCondition ) / HQ ) / 3.D0 ! 0 ≤ ϕ ≤ π
ELSE
  Phi = 0.D0
END IF

! Eigenvalues of Q | Cardano's solution
Eigenvalues(1) = M + ( 2.D0 * DSQRT( P ) * DCOS( Phi ) )
Eigenvalues(2) = M - DSQRT( P ) * ( DCOS( Phi ) + ( DSQRT( 3.D0 ) * DSIN( Phi ) ) )
Eigenvalues(3) = M - DSQRT( P ) * ( DCOS( Phi ) - ( DSQRT( 3.D0 ) * DSIN( Phi ) ) )

! Nematic order parameter
OrderParameterS2 = MAXVAL( Eigenvalues ) ! Largest eigenvalue

! Phase director
CALL NematicPhaseDirector( TensorQ, OrderParameterS2, PhaseDirector )

RETURN

END SUBROUTINE NematicOrderParameter

! ********************************************************************************************************************************* !
! This subroutine solves a system of linear equations using Cramer's rule method.                                                   !
!                                                                                                                                   !
! ● The solution is the <Eigenvector> associated with one of the eigenvalues (<OrderParameter>)                                     !
!   of the Q-tensor matrix (<Matrix>).                                                                                              !
! ********************************************************************************************************************************* !
! OBS.: This calculation does not differ between antiparallel nematic phase directors, n and -n., which are usually distinguished   !
!       for ferroelectric liquid crystals and related to the dipole moment.                                                         !
! ********************************************************************************************************************************* !
SUBROUTINE NematicPhaseDirector( Matrix, OrderParameter, Eigenvector )

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 ) :: i, j, k   ! Counters
INTEGER( Kind= Int64 ) :: ArraySize ! Size of the matrix
INTEGER( Kind= Int64 ) :: row, col  ! Rows and columns

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 )               :: Determinant, DeterminantCramer ! Determinants
REAL( Kind= Real64 ), INTENT( IN ) :: OrderParameter                 ! Nematic order parameter

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( : ), INTENT( OUT )   :: Eigenvector ! Eigenvector associated with the largest eigenvalue
REAL( Kind= Real64 ), DIMENSION( :, : ), INTENT( IN ) :: Matrix      ! Input tensor Q matrix

! REAL VARIABLES (ARRAY,ALLOCATABLE)
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE :: AuxMatrix     ! Auxiliary matrix
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE :: ReducedMatrix ! Reduced matrix
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE :: CramerMatrix  ! Cramer's matrix
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE    :: ResultVector  ! Result vector

! Check if the matrix is square
IF( SIZE( Matrix, Dim= 1 ) /= SIZE( Matrix, Dim= 2 ) ) THEN
  WRITE( *, '(A)' ) 'The matrix is not square!'
  STOP
END IF

! Check if the matrix is 3x3
IF( SIZE( Matrix, Dim= 1 ) /= 3 ) THEN
  WRITE( *, '(A)' ) 'The matrix is not 3x3!'
  STOP
END IF

! Get the size of the matrix
ArraySize = SIZE( Matrix, Dim= 1 )

! Allocation
IF( ALLOCATED( ReducedMatrix ) ) DEALLOCATE( ReducedMatrix )
IF( ALLOCATED( AuxMatrix ) ) DEALLOCATE( AuxMatrix )
IF( ALLOCATED( ResultVector ) ) DEALLOCATE( ResultVector )
IF( ALLOCATED( CramerMatrix ) ) DEALLOCATE( CramerMatrix )
ALLOCATE( ReducedMatrix( ArraySize - 1, ArraySize - 1 ) )
ReducedMatrix = 0.D0
ALLOCATE( CramerMatrix( ArraySize - 1, ArraySize - 1 ) )
CramerMatrix = 0.D0
ALLOCATE( AuxMatrix( ArraySize, ArraySize ) )
AuxMatrix = 0.D0
ALLOCATE( ResultVector( ArraySize - 1 ) )
ResultVector = 0.D0

! Initialize the eigenvector
Eigenvector = 1.D0

! Rewrite input matrix by subtracting the diagonal matrix of the eigenvalues (A - λI)
DO i = 1, ArraySize
  DO j = 1, ArraySize
    AuxMatrix(i,j) = Matrix(i,j) - OrderParameter * IdentityMatrix(i,j)
  END DO
END DO

! Reduced matrix
k = 1
DO
  ! Stop if all reduced matrices are singular
  IF( k > 3 ) THEN
    WRITE( *, '(A)' ) 'Cannot find a the eigenvector of the nematic phase!'
    STOP
  END IF
  ! Trim the k-th row and k-th column of the auxiliary matrix
  row = 1
  DO i = 1, ArraySize
    IF( i == k ) CYCLE
    col = 1
    DO j = 1, ArraySize
      IF( j == k ) CYCLE
      ReducedMatrix(row,col) = AuxMatrix(i,j)
      col = col + 1
    END DO
    row = row + 1
  END DO
  ! Determinant of the reduced 3x3 matrix
  Determinant = ReducedMatrix(1,1) * ReducedMatrix(2,2) - ReducedMatrix(1,2) * ReducedMatrix(2,1)
  IF( DABS( Determinant ) <= 1.D-12 ) THEN ! Singular matrix
    k = k + 1 ! Trim another row and column
    CYCLE
  END IF
  EXIT
END DO

! Result vector (last column of the matrix)
j = 1
DO i = 1, ArraySize
  IF( i == k ) CYCLE
  ResultVector(j) = - Matrix(i,k) * Eigenvector(k) ! Xk is assumed to be 1 (could be any value)
  j = j + 1
END DO

! Cramer's matrix (substitution of the result vector in ith column)
j = 1
DO i = 1, SIZE( Eigenvector )
  IF( i == k ) CYCLE
  CramerMatrix = ReducedMatrix
  CramerMatrix(:,j) = ResultVector(:)
  ! Determinant of the Cramer's matrix
  DeterminantCramer = CramerMatrix(1,1) * CramerMatrix(2,2) - CramerMatrix(1,2) * CramerMatrix(2,1)
  ! Eigenvector elements (solution of the system of linear equations)
  Eigenvector(i) = DeterminantCramer / Determinant
  ! Next column of the Cramer's matrix
  j = j + 1
END DO

! Normalize the eigenvector
Eigenvector = Eigenvector / DSQRT( DOT_PRODUCT( Eigenvector, Eigenvector ) )

RETURN

END SUBROUTINE NematicPhaseDirector

! ********************************************************************************************************************************* !
! This subroutine calculates the cubatic order parameter (<cCubaticOrder>) of a liquid crystal.                                     !
!                                                                                                                                   !
! ● It uses simulated annealing to find the maximum value of the cubatic order parameter.                                           !
! ********************************************************************************************************************************* !
! For more information, please read: Haji-Akbari and Glotzer. J. Phys. A: Math. Theor. 48 485201 (2015).                            !
!                                    Duncan et al. Phys. Rev. E 79, 031702 (2009).                                                  !
! Original code in C++ at https://github.com/glotzerlab/freud/blob/main/freud/order/Cubatic.cc.                                     !
! ********************************************************************************************************************************* !
SUBROUTINE CubaticOrderParameter( CurrentFrame, TotalFrames, pCurrentQuaternion, cCubaticOrder )

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 )               :: iCounter, jCounter, kCounter, lCounter     ! Counters
INTEGER( Kind= Int64 )               :: aComponent                                 ! Counter (axis component)
INTEGER( Kind= Int64 )               :: sAnnealingCounter, MaxCounter, iRepetition ! Counter (simulated annealing)
INTEGER( Kind= Int64 )               :: pParticle                                  ! Counter (particle)
INTEGER( Kind= Int64 )               :: OldLenMessage, MaxLenMessage               ! Length of the print message
INTEGER( Kind= Int64 ), INTENT( IN ) :: CurrentFrame, TotalFrames                  ! Frames of the configuration file

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 )                :: KroneckerDeltaTerm ! Kronecker delta term
REAL( Kind= Real64 )                :: NewCubaticOrder    ! Cubatic order parameter (random)
REAL( Kind= Real64 )                :: BoltzmannFactor    ! Boltzmann factor
REAL( Kind= Real64 )                :: CurrentTemp        ! Current temperature
REAL( Kind= Real64 ), INTENT( OUT ) :: cCubaticOrder      ! Cubatic order parameter

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 6 )                             :: KroneckerDelta       ! Kronecker delta
REAL( Kind= Real64 ), DIMENSION( 0:3 )                           :: NewQuaternion        ! Random quaternion
REAL( Kind= Real64 ), DIMENSION( 0:3 )                           :: cQuaternion          ! Quaternion of the cubatic phase
REAL( Kind= Real64 ), DIMENSION( 0:3 )                           :: RandomQuaternion     ! Random quaternion
REAL( Kind= Real64 ), DIMENSION( 3, 3, 3, 3 )                    :: gCubaticTensor       ! Global cubatic tensor
REAL( Kind= Real64 ), DIMENSION( 3, 3, 3, 3 )                    :: rCubaticTensor       ! Random cubatic tensor
REAL( Kind= Real64 ), DIMENSION( 3, 3, 3, 3 )                    :: rNewCubaticTensor    ! New random cubatic tensor
REAL( Kind= Real64 ), DIMENSION( 3, 3, 3, 3 )                    :: TensorDifference     ! Difference between the global and random cubatic tensors
REAL( Kind= Real64 ), DIMENSION( 3, 3 )                          :: xyzOrientationSingle ! Random orientation (along x, y, and z axes)
REAL( Kind= Real64 ), DIMENSION( nParticles, 3, 3 )              :: xyzOrientation       ! Orientation of particles (along x, y, and z axes)
REAL( Kind= Real64 ), DIMENSION( nParticles, 0:3 ), INTENT( IN ) :: pCurrentQuaternion   ! Quaternion of particles

! REAL VARIABLES (ARRAY,ALLOCATABLE)
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE :: CubaticOrderAnnealing ! Cubatic order parameter (simulated annealing)

! CHARACTER STRINGS
CHARACTER( Len= 200 ) :: PrintStatus ! Print status

! Allocation
IF( ALLOCATED( CubaticOrderAnnealing ) ) DEALLOCATE( CubaticOrderAnnealing )
ALLOCATE( CubaticOrderAnnealing( MaxRepetitions ) )
CubaticOrderAnnealing = 0.D0

! Convert quaternion to orientation
DO pParticle = 1, nParticles
  CALL VectorRotation( [ 1.D0, 0.D0, 0.D0 ], pCurrentQuaternion(pParticle,0:3), xyzOrientation(pParticle,1,1:3) ) ! Hamilton's quaternion product (x-axis)
  CALL VectorRotation( [ 0.D0, 1.D0, 0.D0 ], pCurrentQuaternion(pParticle,0:3), xyzOrientation(pParticle,2,1:3) ) ! Hamilton's quaternion product (y-axis)
  CALL VectorRotation( [ 0.D0, 0.D0, 1.D0 ], pCurrentQuaternion(pParticle,0:3), xyzOrientation(pParticle,3,1:3) ) ! Hamilton's quaternion product (z-axis)
END DO

! Construct the global cubatic tensor - Eq. 27
gCubaticTensor = 0.D0
DO iCounter = 1, 3
  DO jCounter = 1, 3
    DO kCounter = 1, 3
      DO lCounter = 1, 3
        ! Kronecker delta
        KroneckerDelta = 0.D0
        IF( iCounter == jCounter ) KroneckerDelta(1) = 1.D0 ! Kronecker delta
        IF( kCounter == lCounter ) KroneckerDelta(2) = 1.D0 ! Kronecker delta
        IF( iCounter == kCounter ) KroneckerDelta(3) = 1.D0 ! Kronecker delta
        IF( jCounter == lCounter ) KroneckerDelta(4) = 1.D0 ! Kronecker delta
        IF( iCounter == lCounter ) KroneckerDelta(5) = 1.D0 ! Kronecker delta
        IF( jCounter == kCounter ) KroneckerDelta(6) = 1.D0 ! Kronecker delta
        ! Kronecker delta term
        KroneckerDeltaTerm = KroneckerDelta(1) * KroneckerDelta(2) + KroneckerDelta(3) * KroneckerDelta(4) + &
        &                    KroneckerDelta(5) * KroneckerDelta(6)
        KroneckerDeltaTerm = 0.4D0 * KroneckerDeltaTerm
        DO aComponent = 1, 3 ! Loop over x, y, and z axes
          DO pParticle = 1, nParticles
            ! Global cubatic tensor
            gCubaticTensor(iCounter,jCounter,kCounter,lCounter) = gCubaticTensor(iCounter,jCounter,kCounter,lCounter) + 2.D0 * &
            &     xyzOrientation(pParticle,aComponent,iCounter) * xyzOrientation(pParticle,aComponent,jCounter) * &
            &     xyzOrientation(pParticle,aComponent,kCounter) * xyzOrientation(pParticle,aComponent,lCounter)
          END DO
        END DO
        ! Global cubatic tensor
        gCubaticTensor(iCounter,jCounter,kCounter,lCounter) = gCubaticTensor(iCounter,jCounter,kCounter,lCounter) / &
        &     DBLE( nParticles ) - KroneckerDeltaTerm
      END DO
    END DO
  END DO
END DO

! Initialization
MaxLenMessage = 0

! Loop over replicates and then choose the one that maximizes the cubatic order parameter
DO iRepetition = 1, MaxRepetitions
  ! Status
  WRITE( PrintStatus, "(4(A,I0),A)" ) "Reading frame: ", CurrentFrame, " of ", TotalFrames, " | Simulated annealing repetition ", &
  &                             iRepetition, " of ", MaxRepetitions, "..."
  OldLenMessage = LEN( TRIM( PrintStatus ) )
  IF( OldLenMessage > MaxLenMessage ) MaxLenMessage = OldLenMessage
  WRITE( *, "(A,A)", Advance= "No" ) CHAR(13), TRIM( PrintStatus )
  ! Generate a random quaternion
  CALL RandomQuaternionGenerator( cQuaternion, 0.D0 ) ! Start with a scalar quaternion
  NewQuaternion = cQuaternion
  ! Random orientation (along x, y, and z axes)
  CALL VectorRotation( [ 1.D0, 0.D0, 0.D0 ], cQuaternion, xyzOrientationSingle(1,1:3) ) ! Hamilton's quaternion product (x-axis)
  CALL VectorRotation( [ 0.D0, 1.D0, 0.D0 ], cQuaternion, xyzOrientationSingle(2,1:3) ) ! Hamilton's quaternion product (y-axis)
  CALL VectorRotation( [ 0.D0, 0.D0, 1.D0 ], cQuaternion, xyzOrientationSingle(3,1:3) ) ! Hamilton's quaternion product (z-axis)
  ! Construct a random cubatic tensor from the random orientation
  rCubaticTensor = 0.D0
  DO iCounter = 1, 3
    DO jCounter = 1, 3
      DO kCounter = 1, 3
        DO lCounter = 1, 3
          ! Kronecker delta
          KroneckerDelta = 0.D0
          IF( iCounter == jCounter ) KroneckerDelta(1) = 1.D0 ! Kronecker delta
          IF( kCounter == lCounter ) KroneckerDelta(2) = 1.D0 ! Kronecker delta
          IF( iCounter == kCounter ) KroneckerDelta(3) = 1.D0 ! Kronecker delta
          IF( jCounter == lCounter ) KroneckerDelta(4) = 1.D0 ! Kronecker delta
          IF( iCounter == lCounter ) KroneckerDelta(5) = 1.D0 ! Kronecker delta
          IF( jCounter == kCounter ) KroneckerDelta(6) = 1.D0 ! Kronecker delta
          ! Kronecker delta term
          KroneckerDeltaTerm = KroneckerDelta(1) * KroneckerDelta(2) + KroneckerDelta(3) * KroneckerDelta(4) + &
          &                    KroneckerDelta(5) * KroneckerDelta(6)
          KroneckerDeltaTerm = 0.4D0 * KroneckerDeltaTerm
          DO aComponent = 1, 3 ! Loop over x, y, and z axes
            ! Random cubatic tensor
            rCubaticTensor(iCounter,jCounter,kCounter,lCounter) = rCubaticTensor(iCounter,jCounter,kCounter,lCounter) + 2.D0 * &
            &     xyzOrientationSingle(aComponent,iCounter) * xyzOrientationSingle(aComponent,jCounter) * &
            &     xyzOrientationSingle(aComponent,kCounter) * xyzOrientationSingle(aComponent,lCounter)
          END DO
          ! Random cubatic tensor
          rCubaticTensor(iCounter,jCounter,kCounter,lCounter) = rCubaticTensor(iCounter,jCounter,kCounter,lCounter) - &
          &     KroneckerDeltaTerm
        END DO
      END DO
    END DO
  END DO
  ! Calculate cubatic order parameter
  TensorDifference = gCubaticTensor - rCubaticTensor ! Difference between the global and random cubatic tensors
  cCubaticOrder = 1.D0 - DotProductCubaticTensor( TensorDifference, TensorDifference ) / DotProductCubaticTensor( &
  &              rCubaticTensor, rCubaticTensor ) ! Objective function
  NewCubaticOrder = cCubaticOrder
  ! Initialization
  CurrentTemp = InitialTemp
  sAnnealingCounter = 0
  MaxCounter = MaxSimulatedAnnealingCounter
  ! Simulated annealing (optimization | maximization)
  DO WHILE( sAnnealingCounter < MaxCounter .AND. CurrentTemp > FinalTemp ) ! Stop criterion
    ! Update counter
    sAnnealingCounter = sAnnealingCounter + 1
    ! Find new random quaternion
    CALL RandomQuaternionGenerator( RandomQuaternion, MaxAngularDisplc )
    CALL QuaternionMultiplication( cQuaternion, RandomQuaternion, NewQuaternion )
    ! Random orientation (along x, y, and z axes)
    CALL VectorRotation( [ 1.D0, 0.D0, 0.D0 ], NewQuaternion, xyzOrientationSingle(1,1:3) ) ! Hamilton's quaternion product (x-axis)
    CALL VectorRotation( [ 0.D0, 1.D0, 0.D0 ], NewQuaternion, xyzOrientationSingle(2,1:3) ) ! Hamilton's quaternion product (y-axis)
    CALL VectorRotation( [ 0.D0, 0.D0, 1.D0 ], NewQuaternion, xyzOrientationSingle(3,1:3) ) ! Hamilton's quaternion product (z-axis)
    ! Construct a random cubatic tensor
    rNewCubaticTensor = 0.D0
    DO iCounter = 1, 3
      DO jCounter = 1, 3
        DO kCounter = 1, 3
          DO lCounter = 1, 3
            ! Kronecker delta
            KroneckerDelta = 0.D0
            IF( iCounter == jCounter ) KroneckerDelta(1) = 1.D0 ! Kronecker delta
            IF( kCounter == lCounter ) KroneckerDelta(2) = 1.D0 ! Kronecker delta
            IF( iCounter == kCounter ) KroneckerDelta(3) = 1.D0 ! Kronecker delta
            IF( jCounter == lCounter ) KroneckerDelta(4) = 1.D0 ! Kronecker delta
            IF( iCounter == lCounter ) KroneckerDelta(5) = 1.D0 ! Kronecker delta
            IF( jCounter == kCounter ) KroneckerDelta(6) = 1.D0 ! Kronecker delta
            ! Kronecker delta term
            KroneckerDeltaTerm = KroneckerDelta(1) * KroneckerDelta(2) + KroneckerDelta(3) * KroneckerDelta(4) + &
            &                    KroneckerDelta(5) * KroneckerDelta(6)
            KroneckerDeltaTerm = 0.4D0 * KroneckerDeltaTerm
            ! Random cubatic tensor
            DO aComponent = 1, 3 ! Loop over x, y, and z axes
              ! Random cubatic tensor
              rNewCubaticTensor(iCounter,jCounter,kCounter,lCounter) = rNewCubaticTensor(iCounter,jCounter,kCounter,lCounter) + &
              &     2.D0 * xyzOrientationSingle(aComponent,iCounter) * xyzOrientationSingle(aComponent,jCounter) * &
              &     xyzOrientationSingle(aComponent,kCounter) * xyzOrientationSingle(aComponent,lCounter)
            END DO
            ! Random cubatic tensor
            rNewCubaticTensor(iCounter,jCounter,kCounter,lCounter) = rNewCubaticTensor(iCounter,jCounter,kCounter,lCounter) - &
            &     KroneckerDeltaTerm
          END DO
        END DO
      END DO
    END DO
    ! Calculate cubatic order parameter
    TensorDifference = gCubaticTensor - rNewCubaticTensor
    NewCubaticOrder = 1.D0 - DotProductCubaticTensor( TensorDifference, TensorDifference ) / DotProductCubaticTensor( &
    &                 rNewCubaticTensor, rNewCubaticTensor ) ! Objective function
    ! Acceptance criterion
    IF( NewCubaticOrder > cCubaticOrder ) THEN
      ! Update the cubatic order parameter
      cCubaticOrder = NewCubaticOrder
      ! Update the cubatic tensor
      rCubaticTensor = rNewCubaticTensor
      ! Update the quaternion
      cQuaternion = NewQuaternion
    ELSE
      ! Boltzmann factor
      BoltzmannFactor = DEXP( - ( cCubaticOrder - NewCubaticOrder ) / CurrentTemp )
      ! Random number
      CALL Random_Number( RandomNumber )
      ! Accept worse solutions with a certain probability (allow for uphill moves which prevent the algorithm from getting stuck in local minima)
      IF( BoltzmannFactor >= RandomNumber ) THEN
        ! Update the cubatic order parameter
        cCubaticOrder = NewCubaticOrder
        ! Update the cubatic tensor
        rCubaticTensor = rNewCubaticTensor
        ! Update the quaternion
        cQuaternion = NewQuaternion
      END IF
    END IF
    ! Update the temperature ("cooling")
    CurrentTemp = CurrentTemp * TempScale
  END DO
  CubaticOrderAnnealing(iRepetition) = cCubaticOrder
END DO

! Status
WRITE( PrintStatus, "(A,I0,A,I0)" ) "Reading frame: ", CurrentFrame, " of ", TotalFrames
OldLenMessage = LEN( TRIM( PrintStatus ) )
WRITE( *, "(6G0)", Advance= "No" ) CHAR(13), "Reading frame: ", CurrentFrame, " of ", TotalFrames, &
&                                  REPEAT( " ", MaxLenMessage - OldLenMessage )

! Cubatic order parameter
cCubaticOrder = MAXVAL( CubaticOrderAnnealing )

RETURN

END SUBROUTINE CubaticOrderParameter

! ********************************************************************************************************************************* !
! This function calculates the objective function of the simulated annealing process as the dot product of two cubatic tensors.     !
! ********************************************************************************************************************************* !
DOUBLE PRECISION FUNCTION DotProductCubaticTensor( Tensor_A, Tensor_B )

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 ) :: iCounter, jCounter, kCounter, lCounter ! Counters

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 3, 3, 3, 3 ), INTENT( IN ) :: Tensor_A ! Cubatic tensor A
REAL( Kind= Real64 ), DIMENSION( 3, 3, 3, 3 ), INTENT( IN ) :: Tensor_B ! Cubatic tensor B

! Scalar product of two cubatic tensors
DotProductCubaticTensor = 0.D0
DO iCounter = 1, 3
  DO jCounter = 1, 3
    DO kCounter = 1, 3
      DO lCounter = 1, 3
        DotProductCubaticTensor = DotProductCubaticTensor + Tensor_A(iCounter,jCounter,kCounter,lCounter) * &
        &                         Tensor_B(iCounter,jCounter,kCounter,lCounter)
      END DO
    END DO
  END DO
END DO

RETURN

END FUNCTION DotProductCubaticTensor

! ********************************************************************************************************************************* !
! This subroutine calculates the smectic order parameter of a liquid crystal.                                                       !
!                                                                                                                                   !
! ● Because the smectic order parameter depends on the layer spacing, the optimal layer spacing (<OptimalLayerSpacing>) that        !
!   maximizes the smectic order parameter (<SmecticOrder>) is optimized using the Nelder-Mead simplex algorithm.                    !
! ********************************************************************************************************************************* !
SUBROUTINE SmecticOrderParameter( CurrentFrame, TotalFrames, pCurrentPosition, PhaseDirector, SmecticOrder, OptimalLayerSpacing )

! Uses one module: simplex algorithm
USE SimplexAlgorithm, ONLY: Nelder_Mead

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 ), INTENT( IN ) :: CurrentFrame, TotalFrames ! Frames of the configuration file

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 ), INTENT( OUT ) :: SmecticOrder        ! Smectic order parameter
REAL( Kind= Real64 ), INTENT( OUT ) :: OptimalLayerSpacing ! Optimized layer spacing

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( : ), INTENT( IN )    :: PhaseDirector    ! Phase director of the nematic phase
REAL( Kind= Real64 ), DIMENSION( :, : ), INTENT( IN ) :: pCurrentPosition ! Position of particles

! Optimize layer spacing and calculate the smectic order parameter
CALL Nelder_Mead( 1_INT64, CurrentFrame, TotalFrames, PhaseDirector, pCurrentPosition, SmecticOrder, OptimalLayerSpacing )

RETURN

END SUBROUTINE SmecticOrderParameter

END MODULE OrderParameters