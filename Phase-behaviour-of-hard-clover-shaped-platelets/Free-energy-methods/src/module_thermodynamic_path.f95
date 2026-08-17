! ############################################################################################### !
!                          FREE-ENERGY METHODS FOR NON-CONVEX MOLECULES                           !
!   This module contains the subroutines used to calculate free-energy changes by thermodynamic   !
! integration using the Einstein crystal (Frenkel-Ladd) and externally oriented fluid             !
!                              (Frenkel-Mulder) as reference states.                              !
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
! Main References:        C. Vega, E. Sanz, J. L. F. Abascal, E. G. Noya                          !
!                           J. Phys.: Condens. Matter 20, 153101 (2008)                           !
!                               DOI: 10.1088/0953-8984/20/15/153101                               !
!                             --------------------------------------                              !
!                                    D. Frenkel, B. M. Mulder                                     !
!                               J. Mol. Phys. 55, 1171-1192 (1985)                                !
!                                 DOI: 10.1080/00268978500101971                                  !
! ############################################################################################### !
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !
MODULE ThermodynamicStep

! Uses four modules: global variables, overlap algorithms, overlap checks, quaternion operations, vector operations, and symmetry-dependent orientational fields
USE GlobalVariables
USE OverlapAlgorithms
USE OverlapCheckSystem
USE QuaternionOperations
USE VectorOperations
USE SymmetryOrientationalFields, ONLY: PointGroup_Dnh, SymmetryRotationFunction

IMPLICIT NONE

CONTAINS

! *********************************************************************************************** !
!         This subroutine evaluates the free energy change between an unconstrained solid         !
!    and the solid with fixed CM and also includes the calculation of the ideal contributions     !
! *********************************************************************************************** !
SUBROUTINE Compute_A0(  )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES -*- THIS IS SINGLE PRECISION -*-                                              !
! *********************************************************************************************** !
INTEGER :: SeedSize ! Seed array size

! *********************************************************************************************** !
! INTEGER VARIABLES (ALLOCATABLE) -*- THIS IS SINGLE PRECISION -*-                                !
! *********************************************************************************************** !
INTEGER, DIMENSION( : ), ALLOCATABLE :: Seed ! Random seed

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: Iteration         ! Counter
INTEGER( Kind= Int64 ) :: IntegralDimension ! Integral dimension
INTEGER( Kind= Int64 ) :: MCIterations      ! Maximum simulation cycles
INTEGER( Kind= Int64 ) :: FrameLEFT         ! Box frame dimension
INTEGER( Kind= Int64 ) :: FrameRIGHT        ! Box frame dimension

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                   :: PhiAngle, ThetaAngle, GammaAngle         ! Random φ, θ, and γ angles within domain
REAL( Kind= Real64 )                   :: iFunction                                ! Function value
REAL( Kind= Real64 )                   :: FunctionMean                             ! Function mean
REAL( Kind= Real64 )                   :: IntegralResult                           ! Integral result
REAL( Kind= Real64 )                   :: PsiAngleA, PsiAngleB                     ! Angle between lattice vectors and rotated vectors
REAL( Kind= Real64 )                   :: qAngle                                   ! Angle of the quaternion
REAL( Kind= Real64 )                   :: FreeEnergyTranslationalIDL               ! Ideal translational free energy
REAL( Kind= Real64 )                   :: FreeEnergyRotationalIDL                  ! Ideal rotational free energy
REAL( Kind= Real64 ), DIMENSION( 3 )   :: xOrientation, yOrientation, zOrientation ! Molecular orientation
REAL( Kind= Real64 ), DIMENSION( 3 )   :: xOrientationReference                    ! Molecular orientation along x-direction (reference)
REAL( Kind= Real64 ), DIMENSION( 3 )   :: yOrientationReference                    ! Molecular orientation along y-direction (reference)
REAL( Kind= Real64 ), DIMENSION( 3 )   :: zOrientationReference                    ! Molecular orientation along z-direction (reference)
REAL( Kind= Real64 ), DIMENSION( 3 )   :: RotationAxis                             ! Axis of rotation
REAL( Kind= Real64 ), DIMENSION( 0:3 ) :: rQuaternionInitial                       ! Rotation quaternion (initial)
REAL( Kind= Real64 ), DIMENSION( 0:3 ) :: rQuaternion                              ! Rotation quaternion

! *********************************************************************************************** !
! REAL VARIABLES (ALLOCATABLE)                                                                    !
! *********************************************************************************************** !
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE :: LowerBound       ! Lower bound values
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE :: UpperBound       ! Upper bound values
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE :: FunctionVariable ! Function independent variables

! *********************************************************************************************** !
! CHARACTER STRINGS                                                                               !
! *********************************************************************************************** !
CHARACTER( LEN= 03 ) :: AngleUnit       ! Angle unit (radian or degree)
CHARACTER( LEN= 03 ) :: CalculationType ! Calculation type

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: FixedSeed ! Seed type

! Calculation type
CalculationType = "IDL"

! Initial parameters
OPEN( Unit= 10, File= "ini_parameters.ini", Action= "READ" )
READ( 10, * ) DummyText, IntegralDimension
ALLOCATE( LowerBound(IntegralDimension), UpperBound(IntegralDimension), FunctionVariable(IntegralDimension) )
READ( 10, * ) DummyText, AngleUnit
CALL ToUpper( AngleUnit, LEN_TRIM( AngleUnit ), AngleUnit )
READ( 10, * ) DummyText, LowerBound
READ( 10, * ) DummyText, UpperBound
IF( AngleUnit == "DEG" ) THEN
  LowerBound = LowerBound * cPi / 180.D0
  UpperBound = UpperBound * cPi / 180.D0
END IF
READ( 10, * ) DummyText, MCIterations
READ( 10, * ) DummyText, FixedSeed
CLOSE( 10 )

! Summary
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 70 )//CH_UR
FrameLEFT  = (35 - FLOOR( REAL( 32 ) / 2.D0 ) )
FrameRIGHT = (35 - CEILING( REAL( 32 ) / 2.D0 ) )
WRITE( *, "(2G0)" ) CH_VS//REPEAT( " ", FrameLEFT )//"EINSTEIN CRYSTAL PROPERTIES (A0)"//REPEAT( " ", FrameRIGHT )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 70 )//CH_BR
WRITE( *, "(2G0)" ) "Integral dimension: ", IntegralDimension
WRITE( *, "(2G0)" ) "Angle unit: ", AngleUnit
IF( IntegralDimension == 1 ) THEN
  WRITE( *, "(1(G0,G0.5),G0)" ) "Lower bounds: φ = ", LowerBound(1), "rad"
  WRITE( *, "(1(G0,G0.5),G0)" ) "Upper bounds: φ = ", UpperBound(1), "rad"
ELSE IF( IntegralDimension == 2 ) THEN
  WRITE( *, "(2(G0,G0.5),G0)" ) "Lower bounds: φ = ", LowerBound(1), "rad, θ = ", LowerBound(2), "rad"
  WRITE( *, "(2(G0,G0.5),G0)" ) "Upper bounds: φ = ", UpperBound(1), "rad, θ = ", UpperBound(2), "rad"
ELSE IF( IntegralDimension == 3 ) THEN
  WRITE( *, "(3(G0,G0.5),G0)" ) "Lower bounds: φ = ", LowerBound(1), "rad, θ = ", LowerBound(2), "rad, γ = ", LowerBound(3), &
  &                             "rad"
  WRITE( *, "(3(G0,G0.5),G0)" ) "Upper bounds: φ = ", UpperBound(1), "rad, θ = ", UpperBound(2), "rad, γ = ", UpperBound(3), &
  &                             "rad"
END IF

WRITE( *, "(G0,G0.5)" ) "Maximum number of MC iterations: ", MCIterations
WRITE( *, "(G0,G0.5,G0)" ) "Reduced spring constant: ", rSpringConstant
WRITE( *, "(2G0)" ) "Fixed seed? ", FixedSeed
WRITE( *, "(G0)" ) " "

! Translational free energy (assuming that de Broglie's wavelength is 1Å)
FreeEnergyTranslationalIDL = 1.5D0 * ( DBLE( nParticles ) - 1.D0 ) * DLOG( cPi / rSpringConstant ) + 1.5D0 * &
&                            DLOG( DBLE( nParticles ) ) - DLOG( nDensity )
FreeEnergyTranslationalIDL = - FreeEnergyTranslationalIDL / DBLE( nParticles )

! Quaternion angle
qAngle = 0.D0

! Initial quaternion
rQuaternionInitial(0) = DCOS( qAngle * 0.5D0 )            ! Real part
rQuaternionInitial(1) = DSIN( qAngle * 0.5D0 ) * zAxis(1) ! Imaginary part (Vector)
rQuaternionInitial(2) = DSIN( qAngle * 0.5D0 ) * zAxis(2) ! Imaginary part (Vector)
rQuaternionInitial(3) = DSIN( qAngle * 0.5D0 ) * zAxis(3) ! Imaginary part (Vector)

! Molecular orientation
CALL VectorRotation( xAxis, rQuaternionInitial, xOrientation )
xOrientation = xOrientation / DSQRT( DOT_PRODUCT( xOrientation, xOrientation ) )
xOrientationReference = xOrientation
CALL VectorRotation( yAxis, rQuaternionInitial, yOrientation )
yOrientation = yOrientation / DSQRT( DOT_PRODUCT( yOrientation, yOrientation ) )
yOrientationReference = yOrientation
CALL VectorRotation( zAxis, rQuaternionInitial, zOrientation )
zOrientation = zOrientation / DSQRT( DOT_PRODUCT( zOrientation, zOrientation ) )
zOrientationReference = zOrientation

! Function mean
FunctionMean = 0.D0

! Seed type
IF( FixedSeed ) THEN
  ! Random number generator seed
  CALL Random_Seed( Size= SeedSize )
  IF( ALLOCATED( Seed ) ) DEALLOCATE( Seed )
  ALLOCATE( Seed(SeedSize) )
  Seed = 123456789
  CALL Random_Seed( Put= Seed )
ELSE IF( .NOT. FixedSeed ) THEN
  ! Random pseudorandom number generator seed
  CALL Random_Seed(  )
END IF

! Monte Carlo integration method
DO Iteration = 1, MCIterations
  ! Random X value
  CALL Random_Number( RandomNumber )
  ! φ angle
  PhiAngle = RandomNumber * ( UpperBound(1) - LowerBound(1) ) + LowerBound(1)
  ! Axis of rotation (z-axis)
  RotationAxis = zOrientation
  ! Rotation quaternion
  rQuaternion(0) = DCOS( PhiAngle * 0.5D0 )
  rQuaternion(1) = DSIN( PhiAngle * 0.5D0 ) * RotationAxis(1)
  rQuaternion(2) = DSIN( PhiAngle * 0.5D0 ) * RotationAxis(2)
  rQuaternion(3) = DSIN( PhiAngle * 0.5D0 ) * RotationAxis(3)
  ! Reorient x-axis
  CALL VectorRotation( xOrientation, rQuaternion, xOrientation )
  xOrientation = xOrientation / DSQRT( DOT_PRODUCT( xOrientation, xOrientation ) )
  ! Reorient y-axis
  CALL VectorRotation( yOrientation, rQuaternion, yOrientation )
  yOrientation = yOrientation / DSQRT( DOT_PRODUCT( yOrientation, yOrientation ) )
  IF( IntegralDimension > 1 ) THEN
    ! Random Y value
    CALL Random_Number( RandomNumber )
    ! θ angle
    ThetaAngle = RandomNumber * ( UpperBound(2) - LowerBound(2) ) + LowerBound(2)
    ! Axis of rotation (y'-axis)
    RotationAxis = yOrientation
    ! Rotation quaternion
    rQuaternion(0) = DCOS( ThetaAngle * 0.5D0 )
    rQuaternion(1) = DSIN( ThetaAngle * 0.5D0 ) * RotationAxis(1)
    rQuaternion(2) = DSIN( ThetaAngle * 0.5D0 ) * RotationAxis(2)
    rQuaternion(3) = DSIN( ThetaAngle * 0.5D0 ) * RotationAxis(3)
    ! Reorient x-axis
    CALL VectorRotation( xOrientation, rQuaternion, xOrientation )
    xOrientation = xOrientation / DSQRT( DOT_PRODUCT( xOrientation, xOrientation ) )
    ! Reorient z-axis
    CALL VectorRotation( zOrientation, rQuaternion, zOrientation )
    zOrientation = zOrientation / DSQRT( DOT_PRODUCT( zOrientation, zOrientation ) )
    IF( IntegralDimension > 2 ) THEN
      ! Random Z value
      CALL Random_Number( RandomNumber )
      ! χ/γ angle
      GammaAngle = RandomNumber * ( UpperBound(3) - LowerBound(3) ) + LowerBound(3)
      ! Axis of rotation (z''-axis)
      RotationAxis = zOrientation
      ! Rotation quaternion
      rQuaternion(0) = DCOS( GammaAngle * 0.5D0 )
      rQuaternion(1) = DSIN( GammaAngle * 0.5D0 ) * RotationAxis(1)
      rQuaternion(2) = DSIN( GammaAngle * 0.5D0 ) * RotationAxis(2)
      rQuaternion(3) = DSIN( GammaAngle * 0.5D0 ) * RotationAxis(3)
      ! Reorient x-axis
      CALL VectorRotation( xOrientation, rQuaternion, xOrientation )
      xOrientation = xOrientation / DSQRT( DOT_PRODUCT( xOrientation, xOrientation ) )
      ! Reorient y-axis
      CALL VectorRotation( yOrientation, rQuaternion, yOrientation )
      yOrientation = yOrientation / DSQRT( DOT_PRODUCT( yOrientation, yOrientation ) )
    END IF
  END IF
  ! Angle between lattice vector (y-direction) and rotated vector (y-direction) in an instantaneous configuration
  PsiAngleA = DOT_PRODUCT( yOrientationReference, yOrientation ) / DSQRT( DOT_PRODUCT( yOrientationReference, &
  &           yOrientationReference ) ) / DSQRT( DOT_PRODUCT( yOrientation, yOrientation ) )
  ! Can somehow be greater than |1|
  IF( PsiAngleA > 1.D0 ) THEN
    PsiAngleA = 1.D0
  ELSE IF( PsiAngleA < - 1.D0 ) THEN
    PsiAngleA = - 1.D0
  END IF
  PsiAngleA = DACOS( PsiAngleA )
  ! Angle between lattice vector (z-direction) and rotated vector (z-direction) in an instantaneous configuration
  PsiAngleB = DOT_PRODUCT( zOrientationReference, zOrientation ) / DSQRT( DOT_PRODUCT( zOrientationReference, &
  &           zOrientationReference ) ) / DSQRT( DOT_PRODUCT( zOrientation, zOrientation ) )
  ! Can somehow be greater than |1|
  IF( PsiAngleB > 1.D0 ) THEN
    PsiAngleB = 1.D0
  ELSE IF( PsiAngleB < - 1.D0 ) THEN
    PsiAngleB = - 1.D0
  END IF
  PsiAngleB = DACOS( PsiAngleB )
  ! Function value
  IF( IntegralDimension == 1 ) THEN
    FunctionVariable(1) = ThetaAngle
    iFunction = SymmetryRotationFunction( IntegralDimension, FunctionVariable(1) )
  ELSE IF( IntegralDimension == 2 ) THEN
    FunctionVariable(1) = ThetaAngle
    FunctionVariable(2) = PsiAngleA
    iFunction = SymmetryRotationFunction( IntegralDimension, FunctionVariable(1:2) )
  ELSE IF( IntegralDimension == 3 ) THEN
    FunctionVariable(1) = ThetaAngle
    FunctionVariable(2) = PsiAngleA
    FunctionVariable(3) = PsiAngleB
    iFunction = SymmetryRotationFunction( IntegralDimension, FunctionVariable(1:3) )
  END IF
  ! Accumulate function
  FunctionMean = FunctionMean + iFunction
  ! Progress
  IF( MOD( Iteration, 10000 ) == 0 ) CALL ProgressBar( Iteration, MCIterations, CalculationType )
  ! Reset orientation of molecules
  xOrientation = xOrientationReference
  yOrientation = yOrientationReference
  zOrientation = zOrientationReference
END DO

! Status
WRITE( *, "(G0)" ) " "
WRITE( *, "(G0)" ) " "

! Mean function
FunctionMean = FunctionMean / MCIterations

! Integral result
IntegralResult = FunctionMean * ( ( UpperBound(1) - LowerBound(1) ) * ( UpperBound(2) - LowerBound(2) ) * ( UpperBound(3) - &
&                LowerBound(3) ) )

! Rotational free energy
FreeEnergyRotationalIDL = - DLOG( ( 0.125D0 / (cPi * cPi) ) * IntegralResult )

! Summary
WRITE( *, "(G0)" ) "Writing log..."
WRITE( *, "(G0)" ) " "
WRITE( *, "(G0,G0.5)" ) "Translational ideal free energy is: ", FreeEnergyTranslationalIDL
WRITE( *, "(G0,G0.5)" ) "Orientational ideal free energy is: ", FreeEnergyRotationalIDL
WRITE( *, "(G0,G0.5)" ) "Integral results is: ", IntegralResult
WRITE( *, "(G0,G0.5)" ) "Ideal free energy (A0): ", FreeEnergyTranslationalIDL + FreeEnergyRotationalIDL

! Results file
OPEN( Unit= 10, File= "Results/"//TRIM( DateDescriptor )//"/"//TRIM( HourDescriptor )//"_A0_D"//TRIM( FileDescriptor(1) )// &
&                     "_L"//TRIM( FileDescriptor(2) )//"_LD"//TRIM( FileDescriptor(3) )//".dat" )
WRITE( 10, "(G0,G0.15)" ) "Translational ideal free energy is: ", FreeEnergyTranslationalIDL
WRITE( 10, "(G0,G0.15)" ) "Orientational ideal free energy is: ", FreeEnergyRotationalIDL
WRITE( 10, "(G0,G0.15)" ) "Integral results is: ", IntegralResult
WRITE( 10, "(G0,G0.15)" ) "Total ideal free energy (A0): ", FreeEnergyTranslationalIDL + FreeEnergyRotationalIDL
FLUSH( 10 )
CLOSE( 10 )

RETURN

END SUBROUTINE Compute_A0

! *********************************************************************************************** !
!                        This subroutine evaluates the free energy change                         !
!                       between an interacting EC and a non-interacting EC                        !
! *********************************************************************************************** !
SUBROUTINE Compute_A1(  )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES -*- THIS IS SINGLE PRECISION -*-                                              !
! *********************************************************************************************** !
INTEGER :: SeedSize ! Seed array size

! *********************************************************************************************** !
! INTEGER VARIABLES (ALLOCATABLE) -*- THIS IS SINGLE PRECISION -*-                                !
! *********************************************************************************************** !
INTEGER, DIMENSION( : ), ALLOCATABLE :: Seed ! Random seed

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: pParticle, iParticle, jParticle ! Loop counter (particle)
INTEGER( Kind= Int64 ) :: cCylinder, iCylinder, jCylinder ! Loop counter (cylinder)
INTEGER( Kind= Int64 ) :: Cycles, nValidCycles            ! Loop counter (cycle)
INTEGER( Kind= Int64 ) :: Particle                        ! Random particle
INTEGER( Kind= Int64 ) :: FrameLEFT                       ! Box frame dimension
INTEGER( Kind= Int64 ) :: FrameRIGHT                      ! Box frame dimension
INTEGER( Kind= Int64 ) :: nAcceptanceTranslation          ! Move acceptance counter: Translation
INTEGER( Kind= Int64 ) :: nAcceptanceRotation             ! Move acceptance counter: Rotation
INTEGER( Kind= Int64 ) :: nMovementTranslationCounter     ! Move counter (translation)
INTEGER( Kind= Int64 ) :: nMovementRotationCounter        ! Move counter (rotation)
INTEGER( Kind= Int64 ) :: nNonOverlappingConfigurations   ! Number of non-overlapping configurations

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                    :: ContactDistance                          ! Vega-Lago shortest distance
REAL( Kind= Real64 )                    :: SquaredDistance                          ! Vector distance between particles i and j (squared)
REAL( Kind= Real64 )                    :: FreeEnergyChangeA1                       ! Free energy change between an interacting Einstein crystal and a non-interacting Einstein crystal
REAL( Kind= Real64 )                    :: cPotentialEnergyDifference               ! Energy difference between old and current trial moves (total)
REAL( Kind= Real64 )                    :: OldPotentialEnergy, NewPotentialEnergy   ! Einstein crystal potential energy (before/after a trial move)
REAL( Kind= Real64 )                    :: iOldPotentialEnergy, iNewPotentialEnergy ! Einstein crystal potential energy of particle i (before/after a trial move)
REAL( Kind= Real64 )                    :: jOldPotentialEnergy, jNewPotentialEnergy ! Einstein crystal potential energy of particles j (before/after a trial move)
REAL( Kind= Real64 )                    :: PsiAngleA                                ! Angle between lattice vectors and rotated vectors
REAL( Kind= Real64 )                    :: PsiAngleB                                ! Angle between lattice vectors and rotated vectors
REAL( Kind= Real64 )                    :: MaxTranslationalDisplacement             ! Maximum displacement [+/-] (Translation)
REAL( Kind= Real64 )                    :: MaxAngularDisplacement                   ! Maximum displacement [+/-] (Rotation)
REAL( Kind= Real64 )                    :: PotentialEnergy                          ! Potential energy
REAL( Kind= Real64 )                    :: Ratio                                    ! Acceptance ratio (simulation)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: WritePosition                            ! Position of particles (center of mass)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: bPosition                                ! Position of the center of mass of the system (simulation box)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: ScalingDistanceUnitBox                   ! Scaling factor (unit box)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: iOrientation, jOrientation               ! Orientation of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )    :: yOrientation, zOrientation               ! Orientation of particles along y- and z-directions
REAL( Kind= Real64 ), DIMENSION( 3 )    :: iPosition, jPosition                     ! Position of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )    :: ciPosition, cjPosition                   ! Position of cylinders of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )    :: VectorDistance                           ! Vector distance
REAL( Kind= Real64 ), DIMENSION( 3 )    :: iOldOrientation, iNewOrientation         ! Orientation (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: yOrientationReference                    ! Lattice orientation along y-direction
REAL( Kind= Real64 ), DIMENSION( 3 )    :: zOrientationReference                    ! Lattice orientation along z-direction
REAL( Kind= Real64 ), DIMENSION( 3 )    :: iOldPosition, iNewPosition               ! Position of the center of mass (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: bOldPosition, bNewPosition               ! Position of the center of mass of the system (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( 3, 4 ) :: cOldPosition, cNewPosition               ! Position of cylinders (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( 3, 4 ) :: cRotatedPositionBasis                    ! Relative distance (after a trial move)
REAL( Kind= Real64 ), DIMENSION( 0: 3 ) :: iOldQuaternion, iNewQuaternion           ! Quaternion (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( 0: 3 ) :: iQuaternion, jQuaternion                 ! Quaternion of particles i and j

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: MovementRotationLogical    ! Rotation move selection : TRUE = movement selected; FALSE = movement not selected
LOGICAL :: MovementTranslationLogical ! Translation movement selection : TRUE = movement selected; FALSE = movement not selected
LOGICAL :: Overlap                    ! Detects overlap between two particles (initial) : TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: ParallelSPC                ! Checks the relative orientation of two spherocylinders : TRUE = parallel orientation; FALSE = non-parallel orientation
LOGICAL :: OverlapCYL                 ! Detects overlap between two cylinders : TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: FixedSeed                  ! Seed type

! *********************************************************************************************** !
! CHARACTER STRINGS                                                                               !
! *********************************************************************************************** !
CHARACTER( LEN= 03 )  :: CalculationType  ! Calculation type
CHARACTER( LEN= 140 ) :: DescriptorString ! Descriptor for strings

! Summary
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 70 )//CH_UR
FrameLEFT  = (35 - FLOOR( REAL( 33 ) / 2.D0 ) )
FrameRIGHT = (35 - CEILING( REAL( 33 ) / 2.D0 ) )
WRITE( *, "(2G0)" ) CH_VS//REPEAT( " ", FrameLEFT )//"EINSTEIN CRYSTAL PROPERTIES (ΔA1)"//REPEAT( " ", FrameRIGHT )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 70 )//CH_BR

! Read seed type
OPEN( Unit= 10, File= "ini_montecarlo.ini", Action= "READ" )
READ( 10, * )
READ( 10, * )
READ( 10, * )
READ( 10, * )
READ( 10, * )
READ( 10, * )
READ( 10, * )
READ( 10, * )
READ( 10, * ) DummyText, FixedSeed
CLOSE( 10 )

! Seed type
IF( FixedSeed ) THEN
  ! Random number generator seed
  CALL Random_Seed( Size= SeedSize )
  IF( ALLOCATED( Seed ) ) DEALLOCATE( Seed )
  ALLOCATE( Seed(SeedSize) )
  Seed = 123456789
  CALL Random_Seed( Put= Seed )
ELSE IF( .NOT. FixedSeed ) THEN
  ! Random pseudorandom number generator seed
  CALL Random_Seed(  )
END IF

! Center of mass of system (simulation box)
bPosition(1) = SUM( pPosition(1,:) ) / nParticles
bPosition(2) = SUM( pPosition(2,:) ) / nParticles
bPosition(3) = SUM( pPosition(3,:) ) / nParticles

! Center of mass of system (lattice reference position)
bPositionEq = bPosition

! Active transformation (rotation)
DO pParticle = 1, nParticles
  CALL VectorRotation( zAxis, pQuaternion(:,pParticle), pOrientation(:,pParticle) )
END DO

! Overlap check (initial configuration)
DO iParticle = 1, nParticles - 1
  DO jParticle = iParticle + 1, nParticles
    ! Initialization
    OverlapCYL = .FALSE.
    ! Position of particles i and j
    iPosition(:) = pPosition(:,iParticle)
    jPosition(:) = pPosition(:,jParticle)
    ! Orientation of particles i and j
    iOrientation(:) = pOrientation(:,iParticle)
    jOrientation(:) = pOrientation(:,jParticle)
    ! Quaternion of particles i and j
    iQuaternion(:) = pQuaternion(:,iParticle)
    jQuaternion(:) = pQuaternion(:,jParticle)
    ! Vector distance between particles i and j
    VectorDistance(:) = jPosition(:) - iPosition(:)
    ! Minimum Image Convention (see Allen and Tildesley, 2nd Edition (2017), pages 35-45)
    CALL MatrixVectorMultiplication( BoxLengthInverse, VectorDistance, ScalingDistanceUnitBox )
    ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
    CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, VectorDistance )
    ! Magnitude of the vector distance (squared)
    SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
    ! Cutoff distance (sphere circumscribing a spherocylinder circumscribing the non-convex body)
    IF( SquaredDistance <= SquaredCutoffSphere ) THEN
      ! Cutoff distance (spherocylinder circumscribing the non-convex body)
      CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, ContactDistance, ParallelSPC )
      ! Vega-Lago Criterion
      IF( ContactDistance <= SquaredCutoffSpherocylinder ) THEN
        OverlapCYL = .TRUE.  ! Check overlap between cylinders
      ELSE
        OverlapCYL = .FALSE. ! Do not check overlap between cylinders
      END IF
    END IF
    ! Considering the cylinders (only if preliminary tests fail)
    IF( OverlapCYL ) THEN
      ! First loop takes one of the four cylinders from particle i
      DO iCylinder = 1, 4
        ! Position of cylinder of particle i
        ciPosition(:) = cPosition(:,iCylinder,iParticle)
        ! Second loop takes one of the four cylinders from particle j
        DO jCylinder = 1, 4
          ! Position of cylinder of particle j
          cjPosition(:) = cPosition(:,jCylinder,jParticle)
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
            IF( ContactDistance <= cDiameter * cDiameter ) THEN
              ! Apply periodic boundary conditions on the position of particle j
              cjPosition(1) = ciPosition(1) + VectorDistance(1)
              cjPosition(2) = ciPosition(2) + VectorDistance(2)
              cjPosition(3) = ciPosition(3) + VectorDistance(3)
              ! Time-consuming overlap check
              CALL OverlapCheckCYL( iQuaternion, jQuaternion, iOrientation, jOrientation, VectorDistance, ciPosition, cjPosition, &
              &                     ParallelSPC, Overlap )
              ! Overlap detected
              IF( Overlap ) THEN
                WRITE( *, "(5G0)" ) "Overlap detected in initial configuration between particles ", iParticle, " and ", &
                &                   jParticle, "! Exiting..."
                STOP
              END IF
            END IF
          END IF
        END DO
      END DO
    END IF
  END DO
END DO

! *********************************************************************************************** !
! Output file units                                                                               !
! *********************************************************************************************** !

! Trajectory file (depends on user's choice)
IF( TrajectoryLogical ) THEN
  OPEN( Unit= 20, File= "Trajectories/"//TRIM( DateDescriptor )//"/"//TRIM( HourDescriptor )//"_traj_D" &
  &                     //TRIM( FileDescriptor(1) )//"_L"//TRIM( FileDescriptor(2) )//"_LD"//TRIM( FileDescriptor(3) )//".xyz" )
  WRITE( 20, "(I5)" ) nParticles * 4
  DescriptorString = "(G0,8(G0.6,1X),G0.6,G0,2(G0.6,1X),G0.6,2G0)"
  WRITE( 20, DescriptorString ) 'Lattice="', BoxLength(1:9), '" Origin="', -0.5D0 * ( BoxLength(1) + BoxLength(4) + &
  &                             BoxLength(7) ), -0.5D0 * ( BoxLength(2) + BoxLength(5) + BoxLength(8) ), -0.5D0 * &
  &                             ( BoxLength(3) + BoxLength(6) + BoxLength(9) ), '" ', &
  &                             "Properties=species:S:1:pos:R:3:orientation:R:4:aspherical_shape:R:3"
  DO pParticle = 1, nParticles
    ! Initial configuration for OVITO (reduced units)
    DO cCylinder = 1, 4
      ! Position of cylinders
      WritePosition(1) = cPosition(1,cCylinder,pParticle)
      WritePosition(2) = cPosition(2,cCylinder,pParticle)
      WritePosition(3) = cPosition(3,cCylinder,pParticle)
      WRITE( 20, "(G0,10(' ',G0.9))") pParticle, WritePosition(1), WritePosition(2), WritePosition(3), pQuaternion(1,pParticle), &
      &                               pQuaternion(2,pParticle), pQuaternion(3,pParticle), pQuaternion(0,pParticle), &
      &                               0.5D0 * cDiameter, 0.5D0 * cDiameter, cLength
    END DO
  END DO
  FLUSH( 20 )
END IF

! Ratio file (translation)
OPEN( Unit= 30, File= "Ratio/Translation/"//TRIM( DateDescriptor )//"/"//TRIM( HourDescriptor )//"_ratio_D" &
&                     //TRIM( FileDescriptor(1) )//"_L"//TRIM( FileDescriptor(2) )//"_LD"//TRIM( FileDescriptor(3) )//".dat" )
WRITE( 30, "(G0)" ) "Cycles,Ratio,MaxTranslationalDisplacement,AcceptanceRatioTranslation"
FLUSH( 30 )

! Ratio file (rotation)
OPEN( Unit= 40, File= "Ratio/Rotation/"//TRIM( DateDescriptor )//"/"//TRIM( HourDescriptor )//"_ratio_D" &
&                     //TRIM( FileDescriptor(1) )//"_L"//TRIM( FileDescriptor(2) )//"_LD"//TRIM( FileDescriptor(3) )//".dat" )
WRITE( 40, "(G0)" ) "Cycles,Ratio,MaxAngularDisplacement,AcceptanceRatioRotation"
FLUSH( 40 )

! Results file
OPEN( Unit= 70, File= "Results/"//TRIM( DateDescriptor )//"/"//TRIM( HourDescriptor )//"_results_D" &
&                     //TRIM( FileDescriptor(1) )//"_L"//TRIM( FileDescriptor(2) )//"_LD"//TRIM( FileDescriptor(3) )//".dat" )
WRITE( 70, "(G0)" ) "Cycles,PotentialEnergy,NonOverlappingConfigurations,NonOverlappingFraction"
FLUSH( 70 )

! *********************************************************************************************** !
! Monte Carlo parameters                                                                          !
! *********************************************************************************************** !
MovementTranslationLogical    = .FALSE.                ! Translational move selector              (initial value)
MovementRotationLogical       = .FALSE.                ! Rotational move selector                 (initial value)
MaxTranslationalDisplacement  = maxTranslationalDisplc ! Maximum translational displacement       (initial value)
MaxAngularDisplacement        = maxRotationalDisplc    ! Maximum rotational displacement          (initial value)
nAcceptanceTranslation        = 0                      ! Translational move acceptance counter    (initial value)
nAcceptanceRotation           = 0                      ! Rotational move acceptance counter       (initial value)
nMovementTranslationCounter   = 0                      ! Translational move counter               (initial value)
nMovementRotationCounter      = 0                      ! Rotational move counter                  (initial value)
cPositionMC(:,:,:)            = cPosition(:,:,:)       ! Position (cylinders)                     (initial value)
pQuaternionMC(:,:)            = pQuaternion(:,:)       ! Quaternion algebra                       (initial value)
pPositionMC(:,:)              = pPosition(:,:)         ! Position of particles                    (initial value)
pOrientationMC(:,:)           = pOrientation(:,:)      ! Orientation of particles                 (initial value)
PotentialEnergy               = 0.D0                   ! Intermolecular potential                 (initial value)
nValidCycles                  = 0                      ! Valid Monte Carlo cycles counter         (initial value)
nNonOverlappingConfigurations = 0                      ! Number of overlapping configurations     (initial value)
CalculationType               = "NVT"                  ! Calculation type                         (initial value)

! Simulation cycles
DO Cycles = 1, nSimulationCycles

  ! Progress bar
  IF( MOD( Cycles, nSavingFrequency ) == 0 ) THEN
    CALL ProgressBar( Cycles, nSimulationCycles, CalculationType )
  END IF

  ! Random particle
  CALL Random_Number( RandomNumber )
  Particle = INT( RandomNumber * DBLE( nParticles ) ) + 1

  ! Assignment of previous configuration (microstate m)
  iOldPosition(:)    = pPositionMC(:,Particle)    ! Position
  cOldPosition(:,:)  = cPositionMC(:,:,Particle)  ! Position (cylinders)
  iOldQuaternion(:)  = pQuaternionMC(:,Particle)  ! Quaternion
  iOldOrientation(:) = pOrientationMC(:,Particle) ! Orientation
  bOldPosition(:)    = bPosition(:)               ! Position of the center of mass of the system (simulation box)

  ! Random move selection
  CALL Random_Number( RandomNumber )
  ! Translation criterion
  IF( RandomNumber <= TranslationalProbability ) THEN
    MovementTranslationLogical  = .TRUE.  ! Enable translation
    MovementRotationLogical     = .FALSE. ! Disable rotation
    nMovementTranslationCounter = nMovementTranslationCounter + 1
  ! Rotation criterion
  ELSE IF( RandomNumber > TranslationalProbability ) THEN
    MovementRotationLogical    = .TRUE.  ! Enable rotation
    MovementTranslationLogical = .FALSE. ! Disable translation
    nMovementRotationCounter   = nMovementRotationCounter + 1
  END IF

  ! Translational movement
  IF( MovementTranslationLogical ) THEN
    ! Distance between particle old position and lattice position
    VectorDistance(:) = iOldPosition(:) - pPositionEq(:,Particle)
    ! Minimum image convention
    CALL MatrixVectorMultiplication( BoxLengthInverse, VectorDistance, ScalingDistanceUnitBox )
    ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
    CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, VectorDistance )
    ! Old potential energy of particle i (EC)
    iOldPotentialEnergy = DOT_PRODUCT( ( VectorDistance - ( bOldPosition - bPositionEq ) ), &
    &                                  ( VectorDistance - ( bOldPosition - bPositionEq ) ) )
    ! Old potential energy of all the other particles (EC)
    jOldPotentialEnergy = 0.D0
    DO jParticle = 1, nParticles
      ! Skip the particle that is being moved
      IF( jParticle == Particle ) CYCLE
      ! Distance between particle old position and lattice position
      VectorDistance(:) = pPositionMC(:,jParticle) - pPositionEq(:,jParticle)
      ! Minimum image convention
      CALL MatrixVectorMultiplication( BoxLengthInverse, VectorDistance, ScalingDistanceUnitBox )
      ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
      CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, VectorDistance )
      ! Old potential energy of particle j (EC)
      jOldPotentialEnergy = jOldPotentialEnergy + DOT_PRODUCT( ( VectorDistance - ( bOldPosition - bPositionEq ) ), &
      &                                                        ( VectorDistance - ( bOldPosition - bPositionEq ) ) )
    END DO
    ! Old potential energy of the system
    OldPotentialEnergy = ( iOldPotentialEnergy + jOldPotentialEnergy ) * rSpringConstant
    ! Random translation along x-axis
    CALL Random_Number( RandomNumber )
    iNewPosition(1) = iOldPosition(1) + ( ( 2.D0 * RandomNumber ) - 1.D0 ) * MaxTranslationalDisplacement                  ! Range [-MaxTranslationalDisplacement,MaxTranslationalDisplacement]
    bNewPosition(1) = bOldPosition(1) + ( ( ( 2.D0 * RandomNumber ) - 1.D0 ) * MaxTranslationalDisplacement ) / nParticles ! Range [-MaxTranslationalDisplacement,MaxTranslationalDisplacement]
    ! Random translation along y-axis
    CALL Random_Number( RandomNumber )
    iNewPosition(2) = iOldPosition(2) + ( ( 2.D0 * RandomNumber ) - 1.D0 ) * MaxTranslationalDisplacement                  ! Range [-MaxTranslationalDisplacement,MaxTranslationalDisplacement]
    bNewPosition(2) = bOldPosition(2) + ( ( ( 2.D0 * RandomNumber ) - 1.D0 ) * MaxTranslationalDisplacement ) / nParticles ! Range [-MaxTranslationalDisplacement,MaxTranslationalDisplacement]
    ! Random translation along z-axis
    CALL Random_Number( RandomNumber )
    iNewPosition(3) = iOldPosition(3) + ( ( 2.D0 * RandomNumber ) - 1.D0 ) * MaxTranslationalDisplacement                  ! Range [-MaxTranslationalDisplacement,MaxTranslationalDisplacement]
    bNewPosition(3) = bOldPosition(3) + ( ( ( 2.D0 * RandomNumber ) - 1.D0 ) * MaxTranslationalDisplacement ) / nParticles ! Range [-MaxTranslationalDisplacement,MaxTranslationalDisplacement]
    ! Minimum Image Convention
    CALL MatrixVectorMultiplication( BoxLengthInverse, iNewPosition, ScalingDistanceUnitBox )
    ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
    CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, iNewPosition )
    ! Distance between particle new position and lattice position
    VectorDistance(:) = iNewPosition(:) - pPositionEq(:,Particle)
    ! Minimum image convention
    CALL MatrixVectorMultiplication( BoxLengthInverse, VectorDistance, ScalingDistanceUnitBox )
    ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
    CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, VectorDistance )
    ! New potential energy of particle i (EC)
    iNewPotentialEnergy = DOT_PRODUCT( ( VectorDistance - ( bNewPosition - bPositionEq ) ), &
    &                                  ( VectorDistance - ( bNewPosition - bPositionEq ) ) )
    ! New potential energy of all the other particles (EC)
    jNewPotentialEnergy = 0.D0
    DO jParticle = 1, nParticles
      ! Skip the particle that is being moved
      IF( jParticle == Particle ) CYCLE
      ! Distance between particle old position and lattice position
      VectorDistance(:) = pPositionMC(:,jParticle) - pPositionEq(:,jParticle)
      ! Minimum image convention
      CALL MatrixVectorMultiplication( BoxLengthInverse, VectorDistance, ScalingDistanceUnitBox )
      ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
      CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, VectorDistance )
      ! Old potential energy of particle j (EC)
      jNewPotentialEnergy = jNewPotentialEnergy + DOT_PRODUCT( ( VectorDistance - ( bNewPosition - bPositionEq ) ), &
      &                                                        ( VectorDistance - ( bNewPosition - bPositionEq ) ) )
    END DO
    ! New potential energy of the system
    NewPotentialEnergy = ( iNewPotentialEnergy + jNewPotentialEnergy ) * rSpringConstant
    ! Energy difference
    cPotentialEnergyDifference = NewPotentialEnergy - OldPotentialEnergy
  ! No translational movement
  ELSE IF( .NOT. MovementTranslationLogical ) THEN
    iNewPosition(:) = iOldPosition(:)
    bNewPosition(:) = bOldPosition(:)
  END IF

  ! Rotational movement
  IF( MovementRotationLogical ) THEN
    ! Active transformation (rotation) along y-direction using old quaternion
    CALL VectorRotation( yAxis, iOldQuaternion, yOrientation )
    ! Active transformation (rotation) along z-direction using old quaternion
    CALL VectorRotation( zAxis, iOldQuaternion, zOrientation )
    ! Reference orientations of particle i
    yOrientationReference = pOrientationEq(2,:,Particle)
    zOrientationReference = pOrientationEq(3,:,Particle)
    ! Angle between lattice vector (y-direction) and rotated vector (y-direction) in an instantaneous configuration
    PsiAngleA = DOT_PRODUCT( yOrientationReference, yOrientation ) / DSQRT( DOT_PRODUCT( yOrientationReference, &
    &           yOrientationReference ) ) / DSQRT( DOT_PRODUCT( yOrientation, yOrientation ) )
    IF( PsiAngleA > 1.D0 ) THEN
      PsiAngleA = 1.D0
    ELSE IF( PsiAngleA < - 1.D0 ) THEN
      PsiAngleA = - 1.D0
    END IF
    PsiAngleA = DACOS( PsiAngleA )
    ! Angle between lattice vector (z-direction) and rotated vector (z-direction) in an instantaneous configuration
    PsiAngleB = DOT_PRODUCT( zOrientationReference, zOrientation ) / DSQRT( DOT_PRODUCT( zOrientationReference, &
    &           zOrientationReference ) ) / DSQRT( DOT_PRODUCT( zOrientation, zOrientation ) )
    IF( PsiAngleB > 1.D0 ) THEN
      PsiAngleB = 1.D0
    ELSE IF( PsiAngleB < - 1.D0 ) THEN
      PsiAngleB = - 1.D0
    END IF
    PsiAngleB = DACOS( PsiAngleB )
    ! Old potential energy of particle i
    iOldPotentialEnergy = PointGroup_Dnh( PsiAngleA, PsiAngleB )
    ! Old potential energy of the system
    OldPotentialEnergy = iOldPotentialEnergy * rSpringConstant
    ! New quaternion
    CALL QuaternionCombination( iOldQuaternion, iNewQuaternion, MaxAngularDisplacement )
    ! Active transformation (rotation)
    CALL VectorRotation( zAxis, iNewQuaternion, iNewOrientation )
    ! Active transformation (rotation) along y-direction using new quaternion
    CALL VectorRotation( yAxis, iNewQuaternion, yOrientation )
    ! Active transformation (rotation) along z-direction using new quaternion
    CALL VectorRotation( zAxis, iNewQuaternion, zOrientation )
    ! Reference orientations of particle i
    yOrientationReference = pOrientationEq(2,:,Particle)
    zOrientationReference = pOrientationEq(3,:,Particle)
    ! Angle between lattice vector (y-direction) and rotated vector (y-direction) in an instantaneous configuration
    PsiAngleA = DOT_PRODUCT( yOrientationReference, yOrientation ) / DSQRT( DOT_PRODUCT( yOrientationReference, &
    &           yOrientationReference ) ) / DSQRT( DOT_PRODUCT( yOrientation, yOrientation ) )
    IF( PsiAngleA > 1.D0 ) THEN
      PsiAngleA = 1.D0
    ELSE IF( PsiAngleA < - 1.D0 ) THEN
      PsiAngleA = - 1.D0
    END IF
    PsiAngleA = DACOS( PsiAngleA )
    ! Angle between lattice vector (z-direction) and rotated vector (z-direction) in an instantaneous configuration
    PsiAngleB = DOT_PRODUCT( zOrientationReference, zOrientation ) / DSQRT( DOT_PRODUCT( zOrientationReference, &
    &           zOrientationReference ) ) / DSQRT( DOT_PRODUCT( zOrientation, zOrientation ) )
    IF( PsiAngleB > 1.D0 ) THEN
      PsiAngleB = 1.D0
    ELSE IF( PsiAngleB < - 1.D0 ) THEN
      PsiAngleB = - 1.D0
    END IF
    PsiAngleB = DACOS( PsiAngleB )
    ! New potential energy of particle i
    iNewPotentialEnergy = PointGroup_Dnh( PsiAngleA, PsiAngleB )
    ! New potential energy of the system
    NewPotentialEnergy = iNewPotentialEnergy * rSpringConstant
    ! Energy difference
    cPotentialEnergyDifference = NewPotentialEnergy - OldPotentialEnergy
  ! No rotational movement
  ELSE IF( .NOT. MovementRotationLogical ) THEN
    iNewQuaternion(:) = iOldQuaternion(:)
    iNewOrientation(:) = iOldOrientation(:)
  END IF

  ! Random position of cylinders (after translation or rotation)
  DO cCylinder = 1, 4
    ! Active transformation (translation)
    CALL VectorRotation( cPositionBasis(:,cCylinder), iNewQuaternion, cRotatedPositionBasis(:,cCylinder) )
    cNewPosition(:,cCylinder) = iNewPosition(:) + cRotatedPositionBasis(:,cCylinder)
  END DO

  ! Metropolis criterion
  CALL Random_Number( RandomNumber )
  IF( cPotentialEnergyDifference <= 0.D0 .OR. DEXP( - cPotentialEnergyDifference ) >= RandomNumber ) THEN
    ! Update system configuration
    pPositionMC(:,Particle)    = iNewPosition(:)    ! Update position
    cPositionMC(:,:,Particle)  = cNewPosition(:,:)  ! Update position (cylinders)
    pQuaternionMC(:,Particle)  = iNewQuaternion(:)  ! Update quaternion
    pOrientationMC(:,Particle) = iNewOrientation(:) ! Update orientation
    bPosition(:)               = bNewPosition(:)    ! Update position of system's center of mass (simulation box)
    ! Update energy of the system
    PotentialEnergy = PotentialEnergy + cPotentialEnergyDifference
    ! Update displacement counter
    IF( MovementTranslationLogical ) THEN
      nAcceptanceTranslation = nAcceptanceTranslation + 1 ! Translational move counter
    ELSE IF( MovementRotationLogical ) THEN
      nAcceptanceRotation = nAcceptanceRotation + 1 ! Rotational move counter
    END IF
  ! Move rejected
  ELSE
    ! Retrieve old configuration
    pPositionMC(:,Particle)    = iOldPosition(:)    ! Position
    cPositionMC(:,:,Particle)  = cOldPosition(:,:)  ! Position (cylinders)
    pQuaternionMC(:,Particle)  = iOldQuaternion(:)  ! Quaternion
    pOrientationMC(:,Particle) = iOldOrientation(:) ! Orientation
    bPosition(:)               = bOldPosition(:)    ! Position of system's center of mass (simulation box)
  END IF

  ! Update total energy of the solid
  IF( Cycles > nEquilCycles ) THEN
    CALL FullOverlapCheck( Overlap )
    IF( .NOT. Overlap ) THEN
      nNonOverlappingConfigurations = nNonOverlappingConfigurations + 1
    END IF
    nValidCycles = nValidCycles + 1
  END IF

  ! Adjustment of maximum displacement
  IF( Cycles <= nEquilCycles ) THEN ! During equilibration only
    ! Adjustment frequency (translation and rotation)
    IF( MOD( Cycles, nAdjustmentFrequency ) == 0 ) THEN
      ! Acceptance ratio (non-overlapping microstates over sampled microstates)
      IF( nMovementTranslationCounter > 0 ) THEN
        Ratio = DBLE( nAcceptanceTranslation ) / DBLE( nMovementTranslationCounter )
        ! Translational adjustment
        IF( Ratio <= AcceptanceRatioTranslation ) THEN
          MaxTranslationalDisplacement  = 0.95D0 * MaxTranslationalDisplacement
        ELSE
          MaxTranslationalDisplacement  = 1.05D0 * MaxTranslationalDisplacement
        END IF
        ! Ratio data
        WRITE( 30, "(G0,3(',',G0.9))" ) Cycles, Ratio, MaxTranslationalDisplacement, AcceptanceRatioTranslation
        FLUSH( 30 )
        ! Reset counter
        nAcceptanceTranslation = 0
        nMovementTranslationCounter = 0
        ! Avoid multiple turns
        IF( MaxTranslationalDisplacement >= 4.D0 * MAXVAL( BoxLength ) ) THEN
          MaxTranslationalDisplacement = MaxTranslationalDisplacement - 2.D0 * MAXVAL( BoxLength )
        END IF
      END IF
      ! Acceptance ratio (non-overlapping microstates over sampled microstates)
      IF( nMovementRotationCounter > 0 ) THEN
        Ratio = DBLE( nAcceptanceRotation ) / DBLE( nMovementRotationCounter )
        ! Rotation adjustment
        IF( Ratio <= AcceptanceRatioRotation ) THEN
          MaxAngularDisplacement = 0.95D0 * MaxAngularDisplacement
        ELSE
          MaxAngularDisplacement = 1.05D0 * MaxAngularDisplacement
        END IF
        ! Ratio data
        WRITE( 40, "(G0,3(',',G0.9))" ) Cycles, Ratio, MaxAngularDisplacement, AcceptanceRatioRotation
        FLUSH( 40 )
        ! Reset counter
        nAcceptanceRotation = 0
        nMovementRotationCounter = 0
        ! Avoid multiple turns
        IF( MaxAngularDisplacement >= 4.D0 * cPi ) THEN
          MaxAngularDisplacement = MaxAngularDisplacement - 2.D0 * cPi
        END IF
      END IF
    END IF
  END IF

  ! Trajectory data
  IF( TrajectoryLogical ) THEN
    IF( MOD( Cycles, nSavingFrequency ) == 0 ) THEN
      WRITE( 20, "(I5)" ) nParticles * 4
      DescriptorString = "(G0,8(G0.6,1X),G0.6,G0,2(G0.6,1X),G0.6,2G0)"
      WRITE( 20, DescriptorString ) 'Lattice="', BoxLength(1:9), '" Origin="', -0.5D0 * ( BoxLength(1) + BoxLength(4) + &
      &                             BoxLength(7) ), -0.5D0 * ( BoxLength(2) + BoxLength(5) + BoxLength(8) ), -0.5D0 * &
      &                             ( BoxLength(3) + BoxLength(6) + BoxLength(9) ), '" ', &
      &                             "Properties=species:S:1:pos:R:3:orientation:R:4:aspherical_shape:R:3"
      DO iParticle = 1, nParticles
        ! Position of cylinders
        DO cCylinder = 1, 4
          WritePosition(1) = cPositionMC(1,cCylinder,iParticle)
          WritePosition(2) = cPositionMC(2,cCylinder,iParticle)
          WritePosition(3) = cPositionMC(3,cCylinder,iParticle)
          WRITE( 20, "(G0,10(' ',G0.9))" ) iParticle, WritePosition(1), WritePosition(2), WritePosition(3), &
          &                                pQuaternionMC(1,iParticle), pQuaternionMC(2,iParticle), pQuaternionMC(3,iParticle), &
          &                                pQuaternionMC(0,iParticle), 0.5D0 * cDiameter, 0.5D0 * cDiameter, cLength
        END DO
      END DO
      FLUSH( 20 )
    END IF
  END IF

  ! Results data
  IF( MOD( Cycles, nSavingFrequency ) == 0 ) THEN
    ! Potential energy and number of non-overlapping configurations over sampled cycles
    IF( Cycles <= nEquilCycles ) THEN
      WRITE( 70, "(G0,',',G0.9,',',G0,',',G0.9)" ) Cycles, PotentialEnergy, 0_Int64, 0.D0
    ELSE
      WRITE( 70, "(G0,',',G0.9,',',G0,',',G0.9)" ) Cycles, PotentialEnergy, nNonOverlappingConfigurations, &
      &                                            DBLE( nNonOverlappingConfigurations ) / DBLE( nValidCycles )
    END IF
    FLUSH( 70 )
  END IF

  ! Reset non-overlapping configurations counter
  IF( Cycles == nEquilCycles ) THEN
    nNonOverlappingConfigurations = 0
    nValidCycles = 0
  END IF

END DO

WRITE( *, * ) " "
WRITE( *, * ) " "

! Free energy change between an interacting EC and a non-interacting EC
FreeEnergyChangeA1 = - DLOG( DBLE( nNonOverlappingConfigurations ) / DBLE( nValidCycles ) ) / DBLE( nParticles )

! Output units                                         
IF( TrajectoryLogical ) THEN
  CLOSE( 20 )
END IF
CLOSE( 30 )
CLOSE( 40 )
CLOSE( 70 )

! Summary
WRITE( *, "(G0)" ) "Writing log..."
WRITE( *, "(G0)" ) " "
WRITE( *, "(G0,G0.9)" ) "Free energy change between an interacting EC and a non-interacting EC is: ", FreeEnergyChangeA1

! Results file
OPEN( Unit= 70, File= "Results/"//TRIM( DateDescriptor )//"/"//TRIM( HourDescriptor )//"_A1_D"//TRIM( FileDescriptor(1) )// &
&                     "_L"//TRIM( FileDescriptor(2) )//"_LD"//TRIM( FileDescriptor(3) )//".dat" )
WRITE( 70, "(G0,G0.15)" ) "Free energy change between an interacting EC and a non-interacting EC is: ", FreeEnergyChangeA1
FLUSH( 70 )
CLOSE( 70 )

RETURN

END SUBROUTINE Compute_A1

! *********************************************************************************************** !
!                        This subroutine evaluates the free energy change                         !
!                            between the solid and the interacting EC                             !
! *********************************************************************************************** !
SUBROUTINE Compute_A2(  )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES -*- THIS IS SINGLE PRECISION -*-                                              !
! *********************************************************************************************** !
INTEGER :: SeedSize ! Seed array size

! *********************************************************************************************** !
! INTEGER VARIABLES (ALLOCATABLE) -*- THIS IS SINGLE PRECISION -*-                                !
! *********************************************************************************************** !
INTEGER, DIMENSION( : ), ALLOCATABLE :: Seed ! Random seed

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: pParticle, iParticle, jParticle ! Loop counter (particles)
INTEGER( Kind= Int64 ) :: cCylinder, iCylinder, jCylinder ! Loop counter (cylinders)
INTEGER( Kind= Int64 ) :: Cycles, nEnergySamples          ! Loop counter
INTEGER( Kind= Int64 ) :: qPoint                          ! Quadrature point
INTEGER( Kind= Int64 ) :: qPointInitial, qPointFinal      ! Quadrature point (initial and final)
INTEGER( Kind= Int64 ) :: Particle                        ! Random particle
INTEGER( Kind= Int64 ) :: FrameLEFT                       ! Box frame dimension
INTEGER( Kind= Int64 ) :: FrameRIGHT                      ! Box frame dimension
INTEGER( Kind= Int64 ) :: nAcceptanceTranslation          ! Move acceptance counter: Translation
INTEGER( Kind= Int64 ) :: nAcceptanceRotation             ! Move acceptance counter: Rotation
INTEGER( Kind= Int64 ) :: nMovementTranslationCounter     ! Move counter (Translation)
INTEGER( Kind= Int64 ) :: nMovementRotationCounter        ! Move counter (Rotation)
INTEGER( Kind= Int64 ) :: nNonOverlappingConfigurations   ! Number of non-overlapping configurations

! *********************************************************************************************** !
! REAL PARAMETERS (GAUSS-LEGENDRE QUADRATURE)                                                     !
! *********************************************************************************************** !
REAL( Kind= Real64 ), PARAMETER :: maxCouplingParameter = 1.D0 ! Hamiltonian thermodynamic integration parameter (maximum)
REAL( Kind= Real64 ), PARAMETER :: minCouplingParameter = 0.D0 ! Hamiltonian thermodynamic integration parameter (minimum)
REAL( Kind= Real64 ), PARAMETER :: CouplingParameterMidpoint = 0.5D0 * (maxCouplingParameter + minCouplingParameter) ! Interval
REAL( Kind= Real64 ), PARAMETER :: CouplingParameterHalfWidth = 0.5D0 * (maxCouplingParameter - minCouplingParameter) ! Interval
REAL( Kind= Real64 ), DIMENSION( 10 ), PARAMETER :: GaussLegendreWeights = [ 0.152753387131D0, 0.149172986473D0, 0.142096109318D0, &
&                                                                            0.131688638449D0, 0.118194531962D0, 0.101930119817D0, &
&                                                                            0.083276741577D0, 0.062672048334D0, 0.040601429800D0, &
&                                                                            0.017614007139D0 ] ! Gauss-Legendre Quadrature Weights
REAL( Kind= Real64 ), DIMENSION( 10 ), PARAMETER :: GaussLegendreNodes = [ 0.076526521133D0, 0.227785851142D0, 0.373706088715D0, &
&                                                                          0.510867001951D0, 0.636053680727D0, 0.746331906460D0, &
&                                                                          0.839116971822D0, 0.912234428251D0, 0.963971927278D0, &
&                                                                          0.993128599185D0 ] ! Gauss-Legendre Quadrature Polynomial Roots
REAL( Kind= Real64 ), DIMENSION( 20 ), PARAMETER :: GaussLegendrePoints = [ &
&                                             ( CouplingParameterMidpoint - CouplingParameterHalfWidth * GaussLegendreNodes(10) ), &
&                                             ( CouplingParameterMidpoint - CouplingParameterHalfWidth * GaussLegendreNodes(9) ), &
&                                             ( CouplingParameterMidpoint - CouplingParameterHalfWidth * GaussLegendreNodes(8) ), &
&                                             ( CouplingParameterMidpoint - CouplingParameterHalfWidth * GaussLegendreNodes(7) ), &
&                                             ( CouplingParameterMidpoint - CouplingParameterHalfWidth * GaussLegendreNodes(6) ), &
&                                             ( CouplingParameterMidpoint - CouplingParameterHalfWidth * GaussLegendreNodes(5) ), &
&                                             ( CouplingParameterMidpoint - CouplingParameterHalfWidth * GaussLegendreNodes(4) ), &
&                                             ( CouplingParameterMidpoint - CouplingParameterHalfWidth * GaussLegendreNodes(3) ), &
&                                             ( CouplingParameterMidpoint - CouplingParameterHalfWidth * GaussLegendreNodes(2) ), &
&                                             ( CouplingParameterMidpoint - CouplingParameterHalfWidth * GaussLegendreNodes(1) ), &
&                                             ( CouplingParameterMidpoint + CouplingParameterHalfWidth * GaussLegendreNodes(1) ), &
&                                             ( CouplingParameterMidpoint + CouplingParameterHalfWidth * GaussLegendreNodes(2) ), &
&                                             ( CouplingParameterMidpoint + CouplingParameterHalfWidth * GaussLegendreNodes(3) ), &
&                                             ( CouplingParameterMidpoint + CouplingParameterHalfWidth * GaussLegendreNodes(4) ), &
&                                             ( CouplingParameterMidpoint + CouplingParameterHalfWidth * GaussLegendreNodes(5) ), &
&                                             ( CouplingParameterMidpoint + CouplingParameterHalfWidth * GaussLegendreNodes(6) ), &
&                                             ( CouplingParameterMidpoint + CouplingParameterHalfWidth * GaussLegendreNodes(7) ), &
&                                             ( CouplingParameterMidpoint + CouplingParameterHalfWidth * GaussLegendreNodes(8) ), &
&                                             ( CouplingParameterMidpoint + CouplingParameterHalfWidth * GaussLegendreNodes(9) ), &
&                                             ( CouplingParameterMidpoint + CouplingParameterHalfWidth * GaussLegendreNodes(10) ) ] ! Lambda points

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                    :: ContactDistance                            ! Vega-Lago shortest distance
REAL( Kind= Real64 )                    :: SquaredDistance                            ! Vector distance between particles i and j (squared)
REAL( Kind= Real64 )                    :: FreeEnergyChangeA2                         ! Free energy change between the solid and the interacting Einstein crystal
REAL( Kind= Real64 )                    :: cPotentialEnergyDifference                 ! Coupled energy difference between old and current trial moves
REAL( Kind= Real64 )                    :: SquaredSineAngleDifference                ! Einstein crystal energy difference between old and current trial moves
REAL( Kind= Real64 )                    :: NewPotentialEnergyEC, OldPotentialEnergyEC ! Einstein crystal potential energy (before/after a trial move)
REAL( Kind= Real64 )                    :: OldPotentialEnergy, NewPotentialEnergy     ! Potential energy (before/after a trial move)
REAL( Kind= Real64 )                    :: iOldPotentialEnergy, iNewPotentialEnergy   ! Potential energy of particle i (before/after a trial move)
REAL( Kind= Real64 )                    :: jOldPotentialEnergy, jNewPotentialEnergy   ! Potential energy of particles j (before/after a trial move)
REAL( Kind= Real64 )                    :: PsiAngleA                                  ! Angle between lattice vectors and rotated vectors
REAL( Kind= Real64 )                    :: PsiAngleB                                  ! Angle between lattice vectors and rotated vectors
REAL( Kind= Real64 )                    :: AverageECPotentialEnergy                   ! Average Einstein crystal potential energy
REAL( Kind= Real64 )                    :: ECPotentialEnergy                          ! Einstein crystal potential energy
REAL( Kind= Real64 )                    :: MaxTranslationalDisplacement               ! Maximum displacement [+/-] (Translation)
REAL( Kind= Real64 )                    :: MaxAngularDisplacement                     ! Maximum displacement [+/-] (Rotation)
REAL( Kind= Real64 )                    :: CoupledPotentialEnergy                     ! Coupled potential energy
REAL( Kind= Real64 )                    :: Ratio                                      ! Acceptance ratio (simulation)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: bPosition                                  ! Position of the center of mass of the system (simulation box)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: WritePosition                              ! Position of particles (center of mass)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: ScalingDistanceUnitBox                     ! Scaled coordinates (unit box)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: iOrientation, jOrientation                 ! Orientation of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )    :: yOrientation, zOrientation                 ! Orientation of particles along y- and z-directions
REAL( Kind= Real64 ), DIMENSION( 3 )    :: iPosition, jPosition                       ! Position of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )    :: ciPosition, cjPosition                     ! Position of cylinders of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )    :: VectorDistance                             ! Vector distance between particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )    :: iOldOrientation, iNewOrientation           ! Orientation (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: yOrientationReference                      ! Lattice orientation along y-direction
REAL( Kind= Real64 ), DIMENSION( 3 )    :: zOrientationReference                      ! Lattice orientation along z-direction
REAL( Kind= Real64 ), DIMENSION( 3 )    :: iOldPosition, iNewPosition                 ! Position of the center of mass (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: bOldPosition                               ! Position of the center of mass of the system (before a trial move)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: bNewPosition                               ! Position of the center of mass of the system (after a trial move)
REAL( Kind= Real64 ), DIMENSION( 20 )   :: ECEnsembleAverage                          ! Potential energy at quadrature node
REAL( Kind= Real64 ), DIMENSION( 10 )   :: ECQuadratureContributions                  ! ΔA2 quadrature contributions 
REAL( Kind= Real64 ), DIMENSION( 3, 4 ) :: cOldPosition, cNewPosition                 ! Position of cylinders (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( 3, 4 ) :: cRotatedPositionBasis                      ! Relative distance (after a trial move)
REAL( Kind= Real64 ), DIMENSION( 0: 3 ) :: iOldQuaternion, iNewQuaternion             ! Quaternion (before/after a trial move) and auxiliar quaternion
REAL( Kind= Real64 ), DIMENSION( 0: 3 ) :: iQuaternion, jQuaternion                   ! Quaternion of particles i and j

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: MovementRotationLogical    ! Rotation move selection : TRUE = movement selected; FALSE = movement not selected
LOGICAL :: MovementTranslationLogical ! Translation movement selection : TRUE = movement selected; FALSE = movement not selected
LOGICAL :: Overlap                    ! Detects overlap between two particles (initial) : TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: ParallelSPC                ! Checks the relative orientation of two spherocylinders : TRUE = parallel orientation; FALSE = non-parallel orientation
LOGICAL :: OverlapCYL                 ! Detects overlap between two cylinders : TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: FixedSeed                  ! Seed type

! *********************************************************************************************** !
! CHARACTER STRINGS                                                                               !
! *********************************************************************************************** !
CHARACTER( LEN= 03 ) :: CalculationType   ! Calculation type
CHARACTER( LEN= 140 ) :: DescriptorString ! Descriptor for strings

! Summary
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 70 )//CH_UR
FrameLEFT  = (35 - FLOOR( REAL( 33 ) / 2.D0 ) )
FrameRIGHT = (35 - CEILING( REAL( 33 ) / 2.D0 ) )
WRITE( *, "(2G0)" ) CH_VS//REPEAT( " ", FrameLEFT )//"EINSTEIN CRYSTAL PROPERTIES (ΔA2)"//REPEAT( " ", FrameRIGHT )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 70 )//CH_BR

! Read seed type
OPEN( Unit= 10, File= "ini_montecarlo.ini", Action= "READ" )
READ( 10, * )
READ( 10, * )
READ( 10, * )
READ( 10, * )
READ( 10, * )
READ( 10, * )
READ( 10, * )
READ( 10, * )
READ( 10, * ) DummyText, FixedSeed
CLOSE( 10 )

! Seed type
IF( FixedSeed ) THEN
  ! Random number generator seed
  CALL Random_Seed( Size= SeedSize )
  IF( ALLOCATED( Seed ) ) DEALLOCATE( Seed )
  ALLOCATE( Seed(SeedSize) )
  Seed = 123456789
  CALL Random_Seed( Put= Seed )
ELSE IF( .NOT. FixedSeed ) THEN
  ! Random pseudorandom number generator seed
  CALL Random_Seed(  )
END IF

! Results file
OPEN( Unit= 80, File= "Results/"//TRIM( DateDescriptor )//"/"//TRIM( HourDescriptor )//"_results_D" &
&                     //TRIM( FileDescriptor(1) )//"_L"//TRIM( FileDescriptor(2) )//"_LD" &
&                     //TRIM( FileDescriptor(3) )//"_λpoints.dat" )
WRITE( 80, "(G0)" ) "QuadraturePoint,CouplingParameter,AverageEinsteinCrystalPotentialEnergy"
FLUSH( 80 )

! Select initial and final quadrature points
WRITE( *, "(G0)" ) " "
WRITE( *, "(G0)" ) "Initial and final λ values: "
READ( *, * ) qPointInitial, qPointFinal
WRITE( *, "(G0)" ) " "

! *********************************************************************************************** !
! Gauss-Legendre Quadrature Method                                                                !
! *********************************************************************************************** !
DO qPoint = qPointInitial, qPointFinal

  ! Status
  WRITE( *, "(G0,G0.5,3G0)" ) "Current λ value: ", GaussLegendrePoints(qPoint), " (", qPoint, ")"
  WRITE( *, "(G0)" ) " "

  ! Descriptor (lambda)
  WRITE( FileDescriptor(4), "(G0)" ) qPoint

  ! Center of mass of the system
  bPosition(1) = SUM( pPosition(1,:) ) / nParticles
  bPosition(2) = SUM( pPosition(2,:) ) / nParticles
  bPosition(3) = SUM( pPosition(3,:) ) / nParticles

  ! Center of mass of system (lattice reference position)
  bPositionEq = bPosition

  ! Active transformation (rotation)
  DO pParticle = 1, nParticles
    CALL VectorRotation( zAxis, pQuaternion(:,pParticle), pOrientation(:,pParticle) )
  END DO

  ! Overlap check (initial configuration)
  DO iParticle = 1, nParticles - 1 ! First loop represents a particle with a fixed index i
    DO jParticle = iParticle + 1, nParticles ! Second loop represents all other particles with indexes j
      ! Initialization
      OverlapCYL = .FALSE.
      ! Position of particles i and j
      iPosition(:) = pPosition(:,iParticle)
      jPosition(:) = pPosition(:,jParticle)
      ! Orientation of particles i and j
      iOrientation(:) = pOrientation(:,iParticle)
      jOrientation(:) = pOrientation(:,jParticle)
      ! Quaternion of particles i and j
      iQuaternion(:) = pQuaternion(:,iParticle)
      jQuaternion(:) = pQuaternion(:,jParticle)
      ! Vector distance between particles i and j
      VectorDistance(:) = jPosition(:) - iPosition(:)
      ! Minimum Image Convention (see Allen and Tildesley, 2nd Edition (2017), pages 35-45)
      CALL MatrixVectorMultiplication( BoxLengthInverse, VectorDistance, ScalingDistanceUnitBox )
      ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
      CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, VectorDistance )
      ! Magnitude of the vector distance (squared)
      SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
      ! Cutoff distance (sphere circumscribing a spherocylinder circumscribing the non-convex body)
      IF( SquaredDistance <= SquaredCutoffSphere ) THEN
        ! Cutoff distance (spherocylinder circumscribing the non-convex body)
        CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, ContactDistance, ParallelSPC )
        ! Vega-Lago Criterion
        IF( ContactDistance <= SquaredCutoffSpherocylinder ) THEN
          OverlapCYL = .TRUE.  ! Check overlap between cylinders
        ELSE
          OverlapCYL = .FALSE. ! Do not check overlap between cylinders
        END IF
      END IF
      ! Considering the cylinders (only if preliminary tests fail)
      IF( OverlapCYL ) THEN
        ! First loop takes one of the four cylinders from particle i
        DO iCylinder = 1, 4
          ! Position of cylinder of particle i
          ciPosition(:) = cPosition(:,iCylinder,iParticle)
          ! Second loop takes one of the four cylinders from particle j
          DO jCylinder = 1, 4
            ! Position of cylinder of particle j
            cjPosition(:) = cPosition(:,jCylinder,jParticle)
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
              IF( ContactDistance <= cDiameter * cDiameter ) THEN
                cjPosition(1) = ciPosition(1) + VectorDistance(1)
                cjPosition(2) = ciPosition(2) + VectorDistance(2)
                cjPosition(3) = ciPosition(3) + VectorDistance(3)
                ! Time-consuming overlap check
                CALL OverlapCheckCYL( iQuaternion, jQuaternion, iOrientation, jOrientation, VectorDistance, ciPosition, &
                &                     cjPosition, ParallelSPC, Overlap )
                ! Overlap detected
                IF( Overlap ) THEN
                  WRITE( *, "(5G0)" ) "Overlap detected in initial configuration between particles ", iParticle, " and ", &
                  &                   jParticle, "! Exiting..."
                  STOP
                END IF
              END IF
            END IF
          END DO
        END DO
      END IF
    END DO
  END DO

  ! ********************************************************************************************* !
  ! Output file units                                                                             !
  ! ********************************************************************************************* !

  ! Trajectory file
  IF( TrajectoryLogical ) THEN
    OPEN( Unit= 20, File= "Trajectories/"//TRIM( DateDescriptor )//"/"//TRIM( HourDescriptor )//"_traj_D" &
    &                     //TRIM( FileDescriptor(1) )//"_L"//TRIM( FileDescriptor(2) )//"_LD" &
    &                     //TRIM( FileDescriptor(3) )//"_λ"//TRIM( FileDescriptor(4) )//".xyz" )
    WRITE( 20, "(I5)" ) nParticles * 4
    DescriptorString = "(G0,8(G0.6,1X),G0.6,G0,2(G0.6,1X),G0.6,2G0)"
    WRITE( 20, DescriptorString ) 'Lattice="', BoxLength(1:9), '" Origin="', -0.5D0 * ( BoxLength(1) + BoxLength(4) + &
    &                             BoxLength(7) ), -0.5D0 * ( BoxLength(2) + BoxLength(5) + BoxLength(8) ), -0.5D0 * &
    &                             ( BoxLength(3) + BoxLength(6) + BoxLength(9) ), '" ', &
    &                             "Properties=species:S:1:pos:R:3:orientation:R:4:aspherical_shape:R:3"
    DO pParticle = 1, nParticles
      ! Initial configuration for OVITO
      DO cCylinder = 1, 4
        ! Position of cylinders
        WritePosition(1) = cPosition(1,cCylinder,pParticle)
        WritePosition(2) = cPosition(2,cCylinder,pParticle)
        WritePosition(3) = cPosition(3,cCylinder,pParticle)
        WRITE( 20, "(G0,10(' ',G0.9))") pParticle, WritePosition(1), WritePosition(2), WritePosition(3), pQuaternion(1,pParticle), &
        &                               pQuaternion(2,pParticle), pQuaternion(3,pParticle), pQuaternion(0,pParticle), &
        &                               0.5D0 * cDiameter, 0.5D0 * cDiameter, cLength
      END DO
    END DO
    FLUSH( 20 )
  END IF

  ! Ratio file (translation)
  OPEN( Unit= 30, File= "Ratio/Translation/"//TRIM( DateDescriptor )//"/"//TRIM( HourDescriptor )//"_ratio_D" &
  &                     //TRIM( FileDescriptor(1) )//"_L"//TRIM( FileDescriptor(2) )//"_LD" &
  &                     //TRIM( FileDescriptor(3) )//"_λ"//TRIM( FileDescriptor(4) )//".dat" )
  WRITE( 30, "(G0)" ) "Cycles,Ratio,MaxTranslationalDisplacement,AcceptanceRatioTranslation"
  FLUSH( 30 )

  ! Ratio file (rotation)
  OPEN( Unit= 40, File= "Ratio/Rotation/"//TRIM( DateDescriptor )//"/"//TRIM( HourDescriptor )//"_ratio_D" &
  &                     //TRIM( FileDescriptor(1) )//"_L"//TRIM( FileDescriptor(2) )//"_LD" &
  &                     //TRIM( FileDescriptor(3) )//"_λ"//TRIM( FileDescriptor(4) )//".dat" )
  WRITE( 40, "(G0)" ) "Cycles,Ratio,MaxAngularDisplacement,AcceptanceRatioRotation"
  FLUSH( 40 )

  ! Results file
  OPEN( Unit= 70, File= "Results/"//TRIM( DateDescriptor )//"/"//TRIM( HourDescriptor )//"_results_D" &
  &                     //TRIM( FileDescriptor(1) )//"_L"//TRIM( FileDescriptor(2) )//"_LD" &
  &                     //TRIM( FileDescriptor(3) )//"_λ"//TRIM( FileDescriptor(4) )//".dat" )
  WRITE( 70, "(G0)" ) "Cycles,CoupledPotentialEnergy,EinsteinCrystalPotentialEnergy,CumulativeEinsteinCrystalPotentialEnergy"
  FLUSH( 70 )

  ! ********************************************************************************************* !
  ! Monte Carlo parameters                                                                        !
  ! ********************************************************************************************* !
  MovementTranslationLogical    = .FALSE.                ! Translational move selector              (initial value)
  MovementRotationLogical       = .FALSE.                ! Rotational move selector                 (initial value)
  MaxTranslationalDisplacement  = maxTranslationalDisplc ! Maximum translational displacement       (initial value)
  MaxAngularDisplacement        = maxRotationalDisplc    ! Maximum rotational displacement          (initial value)
  nAcceptanceTranslation        = 0                      ! Translational move acceptance counter    (initial value)
  nAcceptanceRotation           = 0                      ! Rotational move acceptance counter       (initial value)
  nMovementTranslationCounter   = 0                      ! Translational move counter               (initial value)
  nMovementRotationCounter      = 0                      ! Rotational move counter                  (initial value)
  cPositionMC(:,:,:)            = cPosition(:,:,:)       ! Position (cylinders)                     (initial value)
  pQuaternionMC(:,:)            = pQuaternion(:,:)       ! Quaternion algebra                       (initial value)
  pPositionMC(:,:)              = pPosition(:,:)         ! Position of particles                    (initial value)
  pOrientationMC(:,:)           = pOrientation(:,:)      ! Orientation of particles                 (initial value)
  CoupledPotentialEnergy        = 0.D0                   ! Coupled potential                        (initial value)
  nNonOverlappingConfigurations = 0                      ! Number of overlapping configurations     (initial value)
  CalculationType               = "NVT"                  ! Calculation type                         (initial value)
  AverageECPotentialEnergy      = 0.D0                   ! Average energy                           (initial value)
  nEnergySamples                = 0                      ! Number of energy samples                 (initial value)

  ! Potential energy of initial configuration
  ECPotentialEnergy = 0.D0
  DO pParticle = 1, nParticles
    ! Distance between particle old position and lattice position
    VectorDistance(:) = pPositionMC(:,pParticle) - pPositionEq(:,pParticle)
    ! Minimum Image Convention
    CALL MatrixVectorMultiplication( BoxLengthInverse, VectorDistance, ScalingDistanceUnitBox )
    ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
    CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, VectorDistance )
    ! Potential energy of the system (translational contribution)
    ECPotentialEnergy = ECPotentialEnergy + DOT_PRODUCT( ( VectorDistance - ( bPosition - bPositionEq ) ), &
    &                                                    ( VectorDistance - ( bPosition - bPositionEq ) ) )
    ! Active transformation (rotation) along y-direction
    CALL VectorRotation( yAxis, pQuaternionMC(:,pParticle), yOrientation )
    ! Active transformation (rotation) along z-direction
    CALL VectorRotation( zAxis, pQuaternionMC(:,pParticle), zOrientation )
    ! Reference orientations of particle i
    yOrientationReference = pOrientationEq(2,:,pParticle)
    zOrientationReference = pOrientationEq(3,:,pParticle)
    ! Angle between lattice vector (y-direction) and rotated vector (y-direction) in an instantaneous configuration
    PsiAngleA = DOT_PRODUCT( yOrientationReference, yOrientation ) / DSQRT( DOT_PRODUCT( yOrientationReference, &
    &           yOrientationReference ) ) / DSQRT( DOT_PRODUCT( yOrientation, yOrientation ) )
    IF( PsiAngleA > 1.D0 ) THEN
      PsiAngleA = 1.D0
    ELSE IF( PsiAngleA < - 1.D0 ) THEN
      PsiAngleA = - 1.D0
    END IF
    PsiAngleA = DACOS( PsiAngleA )
    ! Angle between lattice vector (z-direction) and rotated vector (z-direction) in an instantaneous configuration
    PsiAngleB = DOT_PRODUCT( zOrientationReference, zOrientation ) / DSQRT( DOT_PRODUCT( zOrientationReference, &
    &           zOrientationReference ) ) / DSQRT( DOT_PRODUCT( zOrientation, zOrientation ) )
    IF( PsiAngleB > 1.D0 ) THEN
      PsiAngleB = 1.D0
    ELSE IF( PsiAngleB < - 1.D0 ) THEN
      PsiAngleB = - 1.D0
    END IF
    PsiAngleB = DACOS( PsiAngleB )
    ! Potential energy of the system (rotational contribution)
    ECPotentialEnergy = ECPotentialEnergy + PointGroup_Dnh( PsiAngleA, PsiAngleB )
  END DO
  ! Total potential energy of the system
  ECPotentialEnergy = ECPotentialEnergy * rSpringConstant
  CoupledPotentialEnergy = ECPotentialEnergy * ( 1.D0 - GaussLegendrePoints(qPoint) )

  ! Simulation cycles
  DO Cycles = 1, nSimulationCycles

    ! Progress bar
    IF( MOD( Cycles, nSavingFrequency ) == 0 ) THEN
      CALL ProgressBar( Cycles, nSimulationCycles, CalculationType )
    END IF

    ! Random particle
    CALL Random_Number( RandomNumber )
    Particle = INT( RandomNumber * DBLE( nParticles ) ) + 1

    ! Assignment of previous configuration (microstate m)
    iOldPosition(:)    = pPositionMC(:,Particle)    ! Position
    cOldPosition(:,:)  = cPositionMC(:,:,Particle)  ! Position (cylinders)
    iOldQuaternion(:)  = pQuaternionMC(:,Particle)  ! Quaternion
    iOldOrientation(:) = pOrientationMC(:,Particle) ! Orientation
    bOldPosition(:)    = bPosition(:)               ! Position of the center of mass of the system (simulation box)

    ! Random move selection
    CALL Random_Number( RandomNumber )
    ! Translation criterion
    IF( RandomNumber <= TranslationalProbability ) THEN
      MovementTranslationLogical  = .TRUE.  ! Enable translation
      MovementRotationLogical     = .FALSE. ! Disable rotation
      nMovementTranslationCounter = nMovementTranslationCounter + 1
    ! Rotation criterion
    ELSE IF( RandomNumber > TranslationalProbability ) THEN
      MovementRotationLogical    = .TRUE.  ! Enable rotation
      MovementTranslationLogical = .FALSE. ! Disable translation
      nMovementRotationCounter   = nMovementRotationCounter + 1
    END IF

    ! Translational movement
    IF( MovementTranslationLogical ) THEN
      ! Distance between particle old position and lattice position
      VectorDistance(:) = iOldPosition(:) - pPositionEq(:,Particle)
      ! Minimum Image Convention
      CALL MatrixVectorMultiplication( BoxLengthInverse, VectorDistance, ScalingDistanceUnitBox )
      ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
      CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, VectorDistance )
      ! Old potential energy of particle i (EC)
      iOldPotentialEnergy = DOT_PRODUCT( ( VectorDistance - ( bOldPosition - bPositionEq ) ), &
      &                                  ( VectorDistance - ( bOldPosition - bPositionEq ) ) )
      ! Old potential energy of all the other particles (EC)
      jOldPotentialEnergy = 0.D0
      DO jParticle = 1, nParticles
        ! Skip the particle that is being moved
        IF( jParticle == Particle ) CYCLE
        ! Distance between particle old position and lattice position
        VectorDistance(:) = pPositionMC(:,jParticle) - pPositionEq(:,jParticle)
        ! Minimum image convention
        CALL MatrixVectorMultiplication( BoxLengthInverse, VectorDistance, ScalingDistanceUnitBox )
        ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
        CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, VectorDistance )
        ! Old potential energy of particle j (EC)
        jOldPotentialEnergy = jOldPotentialEnergy + DOT_PRODUCT( ( VectorDistance - ( bOldPosition - bPositionEq ) ), &
        &                                                        ( VectorDistance - ( bOldPosition - bPositionEq ) ) )
      END DO
      ! Old potential energy of the system
      OldPotentialEnergy = ( iOldPotentialEnergy + jOldPotentialEnergy ) * rSpringConstant
      OldPotentialEnergyEC = OldPotentialEnergy
      OldPotentialEnergy = OldPotentialEnergy * ( 1.D0 - GaussLegendrePoints(qPoint) )
      ! Random translation along x-axis
      CALL Random_Number( RandomNumber )
      iNewPosition(1) = iOldPosition(1) + ( ( 2.D0 * RandomNumber ) - 1.D0 ) * MaxTranslationalDisplacement                  ! Range [-MaxTranslationalDisplacement,MaxTranslationalDisplacement]
      bNewPosition(1) = bOldPosition(1) + ( ( ( 2.D0 * RandomNumber ) - 1.D0 ) * MaxTranslationalDisplacement ) / nParticles ! Range [-MaxTranslationalDisplacement,MaxTranslationalDisplacement]
      ! Random translation along y-axis
      CALL Random_Number( RandomNumber )
      iNewPosition(2) = iOldPosition(2) + ( ( 2.D0 * RandomNumber ) - 1.D0 ) * MaxTranslationalDisplacement                  ! Range [-MaxTranslationalDisplacement,MaxTranslationalDisplacement]
      bNewPosition(2) = bOldPosition(2) + ( ( ( 2.D0 * RandomNumber ) - 1.D0 ) * MaxTranslationalDisplacement ) / nParticles ! Range [-MaxTranslationalDisplacement,MaxTranslationalDisplacement]
      ! Random translation along z-axis
      CALL Random_Number( RandomNumber )
      iNewPosition(3) = iOldPosition(3) + ( ( 2.D0 * RandomNumber ) - 1.D0 ) * MaxTranslationalDisplacement                  ! Range [-MaxTranslationalDisplacement,MaxTranslationalDisplacement]
      bNewPosition(3) = bOldPosition(3) + ( ( ( 2.D0 * RandomNumber ) - 1.D0 ) * MaxTranslationalDisplacement ) / nParticles ! Range [-MaxTranslationalDisplacement,MaxTranslationalDisplacement]
      ! Minimum Image Convention
      CALL MatrixVectorMultiplication( BoxLengthInverse, iNewPosition, ScalingDistanceUnitBox )
      ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
      CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, iNewPosition )
      ! Distance between particle new position and lattice position
      VectorDistance(:) = iNewPosition(:) - pPositionEq(:,Particle)
      ! Minimum Image Convention
      CALL MatrixVectorMultiplication( BoxLengthInverse, VectorDistance, ScalingDistanceUnitBox )
      ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
      CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, VectorDistance )
      ! New potential energy of particle i (EC)
      iNewPotentialEnergy = DOT_PRODUCT( ( VectorDistance - ( bNewPosition - bPositionEq ) ), &
      &                                  ( VectorDistance - ( bNewPosition - bPositionEq ) ) )
      ! New potential energy of all the other particles (EC)
      jNewPotentialEnergy = 0.D0
      DO jParticle = 1, nParticles
        ! Skip the particle that is being moved
        IF( jParticle == Particle ) CYCLE
        ! Distance between particle old position and lattice position
        VectorDistance(:) = pPositionMC(:,jParticle) - pPositionEq(:,jParticle)
        ! Minimum image convention
        CALL MatrixVectorMultiplication( BoxLengthInverse, VectorDistance, ScalingDistanceUnitBox )
        ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
        CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, VectorDistance )
        ! Old potential energy of particle j (EC)
        jNewPotentialEnergy = jNewPotentialEnergy + DOT_PRODUCT( ( VectorDistance - ( bNewPosition - bPositionEq ) ), &
        &                                                        ( VectorDistance - ( bNewPosition - bPositionEq ) ) )
      END DO
      ! New potential energy of the system
      NewPotentialEnergy = ( iNewPotentialEnergy + jNewPotentialEnergy ) * rSpringConstant
      NewPotentialEnergyEC = NewPotentialEnergy
      NewPotentialEnergy = NewPotentialEnergy * ( 1.D0 - GaussLegendrePoints(qPoint) )
      ! Energy difference
      cPotentialEnergyDifference = NewPotentialEnergy - OldPotentialEnergy
      SquaredSineAngleDifference = NewPotentialEnergyEC - OldPotentialEnergyEC
    ! No translational movement
    ELSE IF( .NOT. MovementTranslationLogical ) THEN
      iNewPosition(:) = iOldPosition(:)
      bNewPosition(:) = bOldPosition(:)
    END IF

    ! Rotational movement
    IF( MovementRotationLogical ) THEN
      ! Active transformation (rotation) along y-direction using old quaternion
      CALL VectorRotation( yAxis, iOldQuaternion, yOrientation )
      ! Active transformation (rotation) along z-direction using old quaternion
      CALL VectorRotation( zAxis, iOldQuaternion, zOrientation )
      ! Reference orientations of particle i
      yOrientationReference = pOrientationEq(2,:,Particle)
      zOrientationReference = pOrientationEq(3,:,Particle)
      ! Angle between lattice vector (y-direction) and rotated vector (y-direction) in an instantaneous configuration
      PsiAngleA = DOT_PRODUCT( yOrientationReference, yOrientation ) / DSQRT( DOT_PRODUCT( yOrientationReference, &
      &           yOrientationReference ) ) / DSQRT( DOT_PRODUCT( yOrientation, yOrientation ) )
      IF( PsiAngleA > 1.D0 ) THEN
        PsiAngleA = 1.D0
      ELSE IF( PsiAngleA < - 1.D0 ) THEN
        PsiAngleA = - 1.D0
      END IF
      PsiAngleA = DACOS( PsiAngleA )
      ! Angle between lattice vector (z-direction) and rotated vector (z-direction) in an instantaneous configuration
      PsiAngleB = DOT_PRODUCT( zOrientationReference, zOrientation ) / DSQRT( DOT_PRODUCT( zOrientationReference, &
      &           zOrientationReference ) ) / DSQRT( DOT_PRODUCT( zOrientation, zOrientation ) )
      IF( PsiAngleB > 1.D0 ) THEN
        PsiAngleB = 1.D0
      ELSE IF( PsiAngleB < - 1.D0 ) THEN
        PsiAngleB = - 1.D0
      END IF
      PsiAngleB = DACOS( PsiAngleB )
      ! Old potential energy of particle i
      iOldPotentialEnergy = PointGroup_Dnh( PsiAngleA, PsiAngleB )      
      ! Old potential energy of the system
      OldPotentialEnergy = iOldPotentialEnergy * rSpringConstant
      OldPotentialEnergyEC = OldPotentialEnergy
      OldPotentialEnergy = OldPotentialEnergy * ( 1.D0 - GaussLegendrePoints(qPoint) )
      ! New quaternion
      CALL QuaternionCombination( iOldQuaternion, iNewQuaternion, MaxAngularDisplacement )
      ! Active transformation (rotation)
      CALL VectorRotation( zAxis, iNewQuaternion, iNewOrientation )
      ! Active transformation (rotation) along y-direction using new quaternion
      CALL VectorRotation( yAxis, iNewQuaternion, yOrientation )
      ! Active transformation (rotation) along z-direction using new quaternion
      CALL VectorRotation( zAxis, iNewQuaternion, zOrientation )
      ! Reference orientations of particle i
      yOrientationReference = pOrientationEq(2,:,Particle)
      zOrientationReference = pOrientationEq(3,:,Particle)
      ! Angle between lattice vector (y-direction) and rotated vector (y-direction) in an instantaneous configuration
      PsiAngleA = DOT_PRODUCT( yOrientationReference, yOrientation ) / DSQRT( DOT_PRODUCT( yOrientationReference, &
      &           yOrientationReference ) ) / DSQRT( DOT_PRODUCT( yOrientation, yOrientation ) )
      IF( PsiAngleA > 1.D0 ) THEN
        PsiAngleA = 1.D0
      ELSE IF( PsiAngleA < - 1.D0 ) THEN
        PsiAngleA = - 1.D0
      END IF
      PsiAngleA = DACOS( PsiAngleA )
      ! Angle between lattice vector (z-direction) and rotated vector (z-direction) in an instantaneous configuration
      PsiAngleB = DOT_PRODUCT( zOrientationReference, zOrientation ) / DSQRT( DOT_PRODUCT( zOrientationReference, &
      &           zOrientationReference ) ) / DSQRT( DOT_PRODUCT( zOrientation, zOrientation ) )
      IF( PsiAngleB > 1.D0 ) THEN
        PsiAngleB = 1.D0
      ELSE IF( PsiAngleB < - 1.D0 ) THEN
        PsiAngleB = - 1.D0
      END IF
      PsiAngleB = DACOS( PsiAngleB )
      ! New potential energy of particle i
      iNewPotentialEnergy = PointGroup_Dnh( PsiAngleA, PsiAngleB )
      ! New potential energy of the system
      NewPotentialEnergy = iNewPotentialEnergy * rSpringConstant
      NewPotentialEnergyEC = NewPotentialEnergy
      NewPotentialEnergy = NewPotentialEnergy * ( 1.D0 - GaussLegendrePoints(qPoint) )
      ! Energy difference
      cPotentialEnergyDifference = NewPotentialEnergy - OldPotentialEnergy
      SquaredSineAngleDifference = NewPotentialEnergyEC - OldPotentialEnergyEC
    ! No Rotation
    ELSE IF( .NOT. MovementRotationLogical ) THEN
      iNewQuaternion(:) = iOldQuaternion(:)
      iNewOrientation(:) = iOldOrientation(:)
    END IF

    ! Random position of cylinders (after translation or rotation)
    DO cCylinder = 1, 4
      ! Active transformation (translation)
      CALL VectorRotation( cPositionBasis(:,cCylinder), iNewQuaternion, cRotatedPositionBasis(:,cCylinder) )
      cNewPosition(:,cCylinder) = iNewPosition(:) + cRotatedPositionBasis(:,cCylinder)
    END DO

    ! Purely repulsive hard-core potential
    CALL ParticleCheckOverlap( Particle, iNewQuaternion, iNewOrientation, iNewPosition, cNewPosition, Overlap )

    ! Metropolis criterion
    CALL Random_Number( RandomNumber )
    IF( ( cPotentialEnergyDifference <= 0.D0 .OR. DEXP( - cPotentialEnergyDifference ) >= RandomNumber ) .AND. &
    &   ( .NOT. Overlap ) ) THEN
      ! Update system configuration
      pPositionMC(:,Particle)    = iNewPosition(:)    ! Update position
      cPositionMC(:,:,Particle)  = cNewPosition(:,:)  ! Update position (cylinders)
      pQuaternionMC(:,Particle)  = iNewQuaternion(:)  ! Update quaternion
      pOrientationMC(:,Particle) = iNewOrientation(:) ! Update orientation
      bPosition(:)               = bNewPosition(:)    ! Update position of system's center of mass (simulation box)
      ! Update energy of the system
      CoupledPotentialEnergy = CoupledPotentialEnergy + cPotentialEnergyDifference
      ! Update energy of the Einstein crystal
      ECPotentialEnergy = ECPotentialEnergy + SquaredSineAngleDifference
      ! Update displacement counter
      IF( MovementTranslationLogical ) THEN
        nAcceptanceTranslation = nAcceptanceTranslation + 1 ! Translational move counter
      ELSE IF( MovementRotationLogical ) THEN
        nAcceptanceRotation = nAcceptanceRotation + 1 ! Rotational move counter
      END IF
    ! Move rejected
    ELSE
      ! Retrieve old configuration
      pPositionMC(:,Particle)    = iOldPosition(:)    ! Position
      cPositionMC(:,:,Particle)  = cOldPosition(:,:)  ! Position (cylinders)
      pQuaternionMC(:,Particle)  = iOldQuaternion(:)  ! Quaternion
      pOrientationMC(:,Particle) = iOldOrientation(:) ! Orientation
      bPosition(:)               = bOldPosition(:)    ! Position of system's center of mass (simulation box)
    END IF

    ! Adjustment of maximum displacement
    IF( Cycles <= nEquilCycles ) THEN ! During equilibration only
      ! Adjustment frequency (translation and rotation)
      IF( MOD( Cycles, nAdjustmentFrequency ) == 0 ) THEN
        ! Acceptance ratio (non-overlapping microstates over sampled microstates)
        IF( nMovementTranslationCounter > 0 ) THEN
          Ratio = DBLE( nAcceptanceTranslation ) / DBLE( nMovementTranslationCounter )
          ! Translational adjustment
          IF( Ratio <= AcceptanceRatioTranslation ) THEN
            MaxTranslationalDisplacement  = 0.95D0 * MaxTranslationalDisplacement
          ELSE
            MaxTranslationalDisplacement  = 1.05D0 * MaxTranslationalDisplacement
          END IF
          ! Ratio data
          WRITE( 30, "(G0,3(',',G0.9))" ) Cycles, Ratio, MaxTranslationalDisplacement, AcceptanceRatioTranslation
          FLUSH( 30 )
          ! Reset counter
          nAcceptanceTranslation = 0
          nMovementTranslationCounter = 0
          ! Avoid multiple turns
          IF( MaxTranslationalDisplacement >= 4.D0 * MAXVAL( BoxLength ) ) THEN
            MaxTranslationalDisplacement = MaxTranslationalDisplacement - 2.D0 * MAXVAL( BoxLength )
          END IF
        END IF
        ! Acceptance ratio (non-overlapping microstates over sampled microstates)
        IF( nMovementRotationCounter > 0 ) THEN
          Ratio = DBLE( nAcceptanceRotation ) / DBLE( nMovementRotationCounter )
          ! Rotation adjustment
          IF( Ratio <= AcceptanceRatioRotation ) THEN
            MaxAngularDisplacement = 0.95D0 * MaxAngularDisplacement
          ELSE
            MaxAngularDisplacement = 1.05D0 * MaxAngularDisplacement
          END IF
          ! Ratio data
          WRITE( 40, "(G0,3(',',G0.9))" ) Cycles, Ratio, MaxAngularDisplacement, AcceptanceRatioRotation
          FLUSH( 40 )
          ! Reset counter
          nAcceptanceRotation = 0
          nMovementRotationCounter = 0
          ! Avoid multiple turns
          IF( MaxAngularDisplacement >= 4.D0 * cPi ) THEN
            MaxAngularDisplacement = MaxAngularDisplacement - 2.D0 * cPi
          END IF
        END IF
      END IF
    END IF

    ! Trajectory data
    IF( TrajectoryLogical ) THEN
      IF( MOD( Cycles, nSavingFrequency ) == 0 ) THEN
        WRITE( 20, "(I5)" ) nParticles * 4
        DescriptorString = "(G0,8(G0.6,1X),G0.6,G0,2(G0.6,1X),G0.6,2G0)"
        WRITE( 20, DescriptorString ) 'Lattice="', BoxLength(1:9), '" Origin="', -0.5D0 * ( BoxLength(1) + BoxLength(4) + &
        &                             BoxLength(7) ), -0.5D0 * ( BoxLength(2) + BoxLength(5) + BoxLength(8) ), -0.5D0 * &
        &                             ( BoxLength(3) + BoxLength(6) + BoxLength(9) ), '" ', &
        &                             "Properties=species:S:1:pos:R:3:orientation:R:4:aspherical_shape:R:3"
        DO pParticle = 1, nParticles
          ! Position of cylinders
          DO cCylinder = 1, 4
            WritePosition(1) = cPositionMC(1,cCylinder,pParticle)
            WritePosition(2) = cPositionMC(2,cCylinder,pParticle)
            WritePosition(3) = cPositionMC(3,cCylinder,pParticle)
            WRITE( 20, "(G0,10(' ',G0.9))" ) pParticle, WritePosition(1), WritePosition(2), WritePosition(3), &
            &                                pQuaternionMC(1,pParticle), pQuaternionMC(2,pParticle), pQuaternionMC(3,pParticle), &
            &                                pQuaternionMC(0,pParticle), 0.5D0 * cDiameter, 0.5D0 * cDiameter, cLength
          END DO
        END DO
        FLUSH( 20 )
      END IF
    END IF

    ! Results data
    IF( MOD( Cycles, nSavingFrequency ) == 0 ) THEN
      AverageECPotentialEnergy = AverageECPotentialEnergy + ECPotentialEnergy
      nEnergySamples = nEnergySamples + 1
      ! Potential energy
      WRITE( 70, "(G0,3(',',G0.15))" ) Cycles, CoupledPotentialEnergy, ECPotentialEnergy, AverageECPotentialEnergy
      FLUSH( 70 )
    END IF

    ! Reset accumulators
    IF( Cycles == nEquilCycles ) THEN
      AverageECPotentialEnergy = 0.D0
      nEnergySamples = 0
    END IF

  END DO

  WRITE( *, * ) " "
  WRITE( *, * ) " "

  ! Average energy
  AverageECPotentialEnergy = AverageECPotentialEnergy / ( DBLE( nEnergySamples ) )

  ! Potential energy (ensemble average) of quadrature node
  ECEnsembleAverage(qPoint) = AverageECPotentialEnergy

  ! Potential energy of quadrature node
  WRITE( 80, "(G0,2(',',G0.15))" ) qPoint, GaussLegendrePoints(qPoint), ECEnsembleAverage(qPoint)
  FLUSH( 80 )

  ! Output units                                         
  IF( TrajectoryLogical ) THEN
    CLOSE( 20 )
  END IF
  CLOSE( 30 )
  CLOSE( 40 )
  CLOSE( 70 )

END DO

! Output units 
CLOSE( 80 )

! 20-point Gauss-Legendre quadrature
ECQuadratureContributions(1)  = GaussLegendreWeights(1)  * ( ECEnsembleAverage(10) + ECEnsembleAverage(11) )
ECQuadratureContributions(2)  = GaussLegendreWeights(2)  * ( ECEnsembleAverage(9)  + ECEnsembleAverage(12) )
ECQuadratureContributions(3)  = GaussLegendreWeights(3)  * ( ECEnsembleAverage(8)  + ECEnsembleAverage(13) )
ECQuadratureContributions(4)  = GaussLegendreWeights(4)  * ( ECEnsembleAverage(7)  + ECEnsembleAverage(14) )
ECQuadratureContributions(5)  = GaussLegendreWeights(5)  * ( ECEnsembleAverage(6)  + ECEnsembleAverage(15) )
ECQuadratureContributions(6)  = GaussLegendreWeights(6)  * ( ECEnsembleAverage(5)  + ECEnsembleAverage(16) )
ECQuadratureContributions(7)  = GaussLegendreWeights(7)  * ( ECEnsembleAverage(4)  + ECEnsembleAverage(17) )
ECQuadratureContributions(8)  = GaussLegendreWeights(8)  * ( ECEnsembleAverage(3)  + ECEnsembleAverage(18) )
ECQuadratureContributions(9)  = GaussLegendreWeights(9)  * ( ECEnsembleAverage(2)  + ECEnsembleAverage(19) )
ECQuadratureContributions(10) = GaussLegendreWeights(10) * ( ECEnsembleAverage(1)  + ECEnsembleAverage(20) )

! Free energy change between the solid and the interacting Einstein crystal
FreeEnergyChangeA2 = - CouplingParameterHalfWidth * SUM( ECQuadratureContributions ) / DBLE( nParticles )

! Summary
WRITE( *, "(G0)" ) "Writing log..."
WRITE( *, "(G0)" ) " "
WRITE( *, "(G0,G0.5)" ) "Free energy change between the solid and the interacting EC is: ", FreeEnergyChangeA2

! Results file
OPEN( Unit= 95, File= "Results/"//TRIM( DateDescriptor )//"/"//TRIM( HourDescriptor )//"_A2_D"//TRIM( FileDescriptor(1) )// &
&                     "_L"//TRIM( FileDescriptor(2) )//"_LD"//TRIM( FileDescriptor(3) )//"_λENERGY.dat" )
WRITE( 95, "(G0,G0.15)" ) "Free energy change between the solid and the interacting EC is: ", FreeEnergyChangeA2
FLUSH( 95 )
CLOSE( 95 )

RETURN

END SUBROUTINE Compute_A2

! *********************************************************************************************** !
!                        This subroutine evaluates the free energy change                         !
!                            between the solid and the interacting EC                             !
! *********************************************************************************************** !
SUBROUTINE Compute_AField(  )

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES -*- THIS IS SINGLE PRECISION -*-                                              !
! *********************************************************************************************** !
INTEGER :: SeedSize ! Seed array size

! *********************************************************************************************** !
! INTEGER VARIABLES (ALLOCATABLE) -*- THIS IS SINGLE PRECISION -*-                                !
! *********************************************************************************************** !
INTEGER, DIMENSION( : ), ALLOCATABLE :: Seed ! Random seed

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: pParticle, iParticle, jParticle ! Loop counter (particles)
INTEGER( Kind= Int64 ) :: cCylinder, iCylinder, jCylinder ! Loop counter (cylinders)
INTEGER( Kind= Int64 ) :: Cycles, nEnergySamples          ! Loop counter
INTEGER( Kind= Int64 ) :: qPoint                          ! Quadrature point
INTEGER( Kind= Int64 ) :: qPointInitial, qPointFinal      ! Quadrature point (initial and final)
INTEGER( Kind= Int64 ) :: Particle                        ! Random particle
INTEGER( Kind= Int64 ) :: FrameLEFT                       ! Box frame dimension
INTEGER( Kind= Int64 ) :: FrameRIGHT                      ! Box frame dimension
INTEGER( Kind= Int64 ) :: nAcceptanceTranslation          ! Move acceptance counter: Translation
INTEGER( Kind= Int64 ) :: nAcceptanceRotation             ! Move acceptance counter: Rotation
INTEGER( Kind= Int64 ) :: nMovementTranslationCounter     ! Move counter (Translation)
INTEGER( Kind= Int64 ) :: nMovementRotationCounter        ! Move counter (Rotation)

! *********************************************************************************************** !
! REAL PARAMETERS (GAUSS-LEGENDRE QUADRATURE CONSTANTS)                                           !
! *********************************************************************************************** !
REAL( Kind= Real64 ), DIMENSION( 10 ), PARAMETER :: GaussLegendreWeights = [ 0.152753387131D0, 0.149172986473D0, 0.142096109318D0, &
&                                                                            0.131688638449D0, 0.118194531962D0, 0.101930119817D0, &
&                                                                            0.083276741577D0, 0.062672048334D0, 0.040601429800D0, &
&                                                                            0.017614007139D0 ] ! Gauss-Legendre Quadrature Weights
REAL( Kind= Real64 ), DIMENSION( 10 ), PARAMETER :: GaussLegendreNodes = [ 0.076526521133D0, 0.227785851142D0, 0.373706088715D0, &
&                                                                          0.510867001951D0, 0.636053680727D0, 0.746331906460D0, &
&                                                                          0.839116971822D0, 0.912234428251D0, 0.963971927278D0, &
&                                                                          0.993128599185D0 ] ! Gauss-Legendre Quadrature Polynomial Roots

! *********************************************************************************************** !
! REAL VARIABLES (GAUSS-LEGENDRE QUADRATURE VARIABLES)                                            !
! *********************************************************************************************** !
REAL( Kind= Real64 )                  :: maxCouplingParameter       ! Hamiltonian thermodynamic integration parameter (maximum)
REAL( Kind= Real64 )                  :: minCouplingParameter       ! Hamiltonian thermodynamic integration parameter (minimum)
REAL( Kind= Real64 )                  :: CouplingParameterMidpoint  ! Interval
REAL( Kind= Real64 )                  :: CouplingParameterHalfWidth ! Interval
REAL( Kind= Real64 ), DIMENSION( 20 ) :: GaussLegendrePoints        ! Lambda points

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                    :: ContactDistance                          ! Vega-Lago shortest distance
REAL( Kind= Real64 )                    :: CossineOrientFieldAngle                  ! Angle cossine
REAL( Kind= Real64 )                    :: NewPotentialEnergy, OldPotentialEnergy   ! Potential energy (before/after a trial move)
REAL( Kind= Real64 )                    :: PotentialEnergyDifference                ! Potential energy difference between old and current trial moves
REAL( Kind= Real64 )                    :: SquaredDistance                          ! Vector distance between particles i and j (squared)
REAL( Kind= Real64 )                    :: FreeEnergyChangeAField                   ! Free energy change between the externally oriented fluid and the liquid-crystalline phase of interest 
REAL( Kind= Real64 )                    :: SquaredSineAngleDifference               ! Squared sine angle difference between old and current trial moves
REAL( Kind= Real64 )                    :: OldSquaredSineAngle, NewSquaredSineAngle ! Squared sine angle (before/after a trial move)
REAL( Kind= Real64 )                    :: AverageSquaredSineAngle                  ! Average squared sine angle
REAL( Kind= Real64 )                    :: SquaredSineAngle                         ! Squared sine angle
REAL( Kind= Real64 )                    :: MaxTranslationalDisplacement             ! Maximum displacement [+/-] (Translation)
REAL( Kind= Real64 )                    :: MaxAngularDisplacement                   ! Maximum displacement [+/-] (Rotation)
REAL( Kind= Real64 )                    :: OrientationalFieldEnergy                 ! Potential energy of the orientational field of the system
REAL( Kind= Real64 )                    :: Ratio                                    ! Acceptance ratio (simulation)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: ScalingDistanceUnitBox                   ! Scaled coordinates (unit box)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: WritePosition                            ! Position of particles (center of mass)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: iOrientation, jOrientation               ! Orientation of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )    :: iPosition, jPosition                     ! Position of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )    :: ciPosition, cjPosition                   ! Position of cylinders of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )    :: VectorDistance                           ! Vector distance between particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )    :: iOldOrientation, iNewOrientation         ! Orientation (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: iOldPosition, iNewPosition               ! Position of the center of mass (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( 20 )   :: sqSineAngleEnsembleAverage               ! Squared sine angle at quadrature node
REAL( Kind= Real64 ), DIMENSION( 10 )   :: QuadratureContributions                  ! Quadrature contributions 
REAL( Kind= Real64 ), DIMENSION( 3, 4 ) :: cOldPosition, cNewPosition               ! Position of cylinders (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( 3, 4 ) :: cRotatedPositionBasis                    ! Relative distance (after a trial move)
REAL( Kind= Real64 ), DIMENSION( 0: 3 ) :: iOldQuaternion, iNewQuaternion           ! Quaternion (before/after a trial move) and auxiliar quaternion
REAL( Kind= Real64 ), DIMENSION( 0: 3 ) :: iQuaternion, jQuaternion                 ! Quaternion of particles i and j

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: MovementRotationLogical    ! Rotation move selection : TRUE = movement selected; FALSE = movement not selected
LOGICAL :: MovementTranslationLogical ! Translation movement selection : TRUE = movement selected; FALSE = movement not selected
LOGICAL :: Overlap                    ! Detects overlap between two particles (initial) : TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: ParallelSPC                ! Checks the relative orientation of two spherocylinders : TRUE = parallel orientation; FALSE = non-parallel orientation
LOGICAL :: OverlapCYL                 ! Detects overlap between two cylinders : TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: FixedSeed                  ! Seed type

! *********************************************************************************************** !
! CHARACTER STRINGS                                                                               !
! *********************************************************************************************** !
CHARACTER( LEN= 03 ) :: CalculationType   ! Calculation type
CHARACTER( LEN= 140 ) :: DescriptorString ! Descriptor for strings

! GL parameters (initialization)
maxCouplingParameter = OrientFieldStrength ! Hamiltonian thermodynamic integration parameter (maximum)
minCouplingParameter = 0.D0                ! Hamiltonian thermodynamic integration parameter (minimum)
CouplingParameterMidpoint = 0.5D0 * (maxCouplingParameter + minCouplingParameter)  ! Interval
CouplingParameterHalfWidth = 0.5D0 * (maxCouplingParameter - minCouplingParameter) ! Interval
GaussLegendrePoints = [ ( CouplingParameterMidpoint - CouplingParameterHalfWidth * GaussLegendreNodes(10) ), &
&                       ( CouplingParameterMidpoint - CouplingParameterHalfWidth * GaussLegendreNodes(9) ), &
&                       ( CouplingParameterMidpoint - CouplingParameterHalfWidth * GaussLegendreNodes(8) ), &
&                       ( CouplingParameterMidpoint - CouplingParameterHalfWidth * GaussLegendreNodes(7) ), &
&                       ( CouplingParameterMidpoint - CouplingParameterHalfWidth * GaussLegendreNodes(6) ), &
&                       ( CouplingParameterMidpoint - CouplingParameterHalfWidth * GaussLegendreNodes(5) ), &
&                       ( CouplingParameterMidpoint - CouplingParameterHalfWidth * GaussLegendreNodes(4) ), &
&                       ( CouplingParameterMidpoint - CouplingParameterHalfWidth * GaussLegendreNodes(3) ), &
&                       ( CouplingParameterMidpoint - CouplingParameterHalfWidth * GaussLegendreNodes(2) ), &
&                       ( CouplingParameterMidpoint - CouplingParameterHalfWidth * GaussLegendreNodes(1) ), &
&                       ( CouplingParameterMidpoint + CouplingParameterHalfWidth * GaussLegendreNodes(1) ), &
&                       ( CouplingParameterMidpoint + CouplingParameterHalfWidth * GaussLegendreNodes(2) ), &
&                       ( CouplingParameterMidpoint + CouplingParameterHalfWidth * GaussLegendreNodes(3) ), &
&                       ( CouplingParameterMidpoint + CouplingParameterHalfWidth * GaussLegendreNodes(4) ), &
&                       ( CouplingParameterMidpoint + CouplingParameterHalfWidth * GaussLegendreNodes(5) ), &
&                       ( CouplingParameterMidpoint + CouplingParameterHalfWidth * GaussLegendreNodes(6) ), &
&                       ( CouplingParameterMidpoint + CouplingParameterHalfWidth * GaussLegendreNodes(7) ), &
&                       ( CouplingParameterMidpoint + CouplingParameterHalfWidth * GaussLegendreNodes(8) ), &
&                       ( CouplingParameterMidpoint + CouplingParameterHalfWidth * GaussLegendreNodes(9) ), &
&                       ( CouplingParameterMidpoint + CouplingParameterHalfWidth * GaussLegendreNodes(10) ) ] ! Quadrature points (field strength)

! Summary
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 70 )//CH_UR
FrameLEFT  = (35 - FLOOR( REAL( 35 ) / 2.D0 ) )
FrameRIGHT = (35 - CEILING( REAL( 35 ) / 2.D0 ) )
WRITE( *, "(2G0)" ) CH_VS//REPEAT( " ", FrameLEFT )//"FRENKEL-MULDER PROPERTIES (ΔAfield)"//REPEAT( " ", FrameRIGHT )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 70 )//CH_BR

! Read seed type
OPEN( Unit= 10, File= "ini_montecarlo.ini", Action= "READ" )
READ( 10, * )
READ( 10, * )
READ( 10, * )
READ( 10, * )
READ( 10, * )
READ( 10, * )
READ( 10, * )
READ( 10, * )
READ( 10, * ) DummyText, FixedSeed
CLOSE( 10 )

! Seed type
IF( FixedSeed ) THEN
  ! Random number generator seed
  CALL Random_Seed( Size= SeedSize )
  IF( ALLOCATED( Seed ) ) DEALLOCATE( Seed )
  ALLOCATE( Seed(SeedSize) )
  Seed = 123456789
  CALL Random_Seed( Put= Seed )
ELSE IF( .NOT. FixedSeed ) THEN
  ! Random pseudorandom number generator seed
  CALL Random_Seed(  )
END IF

! Results file
OPEN( Unit= 80, File= "Results/"//TRIM( DateDescriptor )//"/"//TRIM( HourDescriptor )//"_results_D" &
&                     //TRIM( FileDescriptor(1) )//"_L"//TRIM( FileDescriptor(2) )//"_LD" &
&                     //TRIM( FileDescriptor(3) )//"_λpoints.dat" )
WRITE( 80, "(G0)" ) "QuadraturePoint,OrientationalFieldStrength,AverageOrientationalPotentialEnergy"
FLUSH( 80 )

! Select initial and final quadrature points
WRITE( *, "(G0)" ) " "
WRITE( *, "(G0)" ) "Initial and final ξ values: "
READ( *, * ) qPointInitial, qPointFinal
WRITE( *, "(G0)" ) " "

! *********************************************************************************************** !
! Gauss-Legendre Quadrature Method                                                                !
! *********************************************************************************************** !
DO qPoint = qPointInitial, qPointFinal

  ! Status
  WRITE( *, "(G0,G0.5,3G0)" ) "Current λ value: ", GaussLegendrePoints(qPoint), " (", qPoint, ")"
  WRITE( *, "(G0)" ) " "

  ! Descriptor (lambda)
  WRITE( FileDescriptor(4), "(G0)" ) qPoint

  ! Active transformation (rotation)
  DO pParticle = 1, nParticles
    CALL VectorRotation( zAxis, pQuaternion(:,pParticle), pOrientation(:,pParticle) )
  END DO

  ! Overlap check (initial configuration)
  DO iParticle = 1, nParticles - 1 ! First loop represents a particle with a fixed index i
    DO jParticle = iParticle + 1, nParticles ! Second loop represents all other particles with indexes j
      ! Initialization
      OverlapCYL = .FALSE.
      ! Position of particles i and j
      iPosition(:) = pPosition(:,iParticle)
      jPosition(:) = pPosition(:,jParticle)
      ! Orientation of particles i and j
      iOrientation(:) = pOrientation(:,iParticle)
      jOrientation(:) = pOrientation(:,jParticle)
      ! Quaternion of particles i and j
      iQuaternion(:) = pQuaternion(:,iParticle)
      jQuaternion(:) = pQuaternion(:,jParticle)
      ! Vector distance between particles i and j
      VectorDistance(:) = jPosition(:) - iPosition(:)
      ! Minimum Image Convention (see Allen and Tildesley, 2nd Edition (2017), pages 35-45)
      CALL MatrixVectorMultiplication( BoxLengthInverse, VectorDistance, ScalingDistanceUnitBox )
      ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
      CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, VectorDistance )
      ! Magnitude of the vector distance (squared)
      SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
      ! Cutoff distance (sphere circumscribing a spherocylinder circumscribing the non-convex body)
      IF( SquaredDistance <= SquaredCutoffSphere ) THEN
        ! Cutoff distance (spherocylinder circumscribing the non-convex body)
        CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, ContactDistance, ParallelSPC )
        ! Vega-Lago Criterion
        IF( ContactDistance <= SquaredCutoffSpherocylinder ) THEN
          OverlapCYL = .TRUE.  ! Check overlap between cylinders
        ELSE
          OverlapCYL = .FALSE. ! Do not check overlap between cylinders
        END IF
      END IF
      ! Considering the cylinders (only if preliminary tests fail)
      IF( OverlapCYL ) THEN
        ! First loop takes one of the four cylinders from particle i
        DO iCylinder = 1, 4
          ! Position of cylinder of particle i
          ciPosition(:) = cPosition(:,iCylinder,iParticle)
          ! Second loop takes one of the four cylinders from particle j
          DO jCylinder = 1, 4
            ! Position of cylinder of particle j
            cjPosition(:) = cPosition(:,jCylinder,jParticle)
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
              IF( ContactDistance <= cDiameter * cDiameter ) THEN
                cjPosition(1) = ciPosition(1) - VectorDistance(1)
                cjPosition(2) = ciPosition(2) - VectorDistance(2)
                cjPosition(3) = ciPosition(3) - VectorDistance(3)
                ! Time-consuming overlap check
                CALL OverlapCheckCYL( iQuaternion, jQuaternion, iOrientation, jOrientation, VectorDistance, ciPosition, &
                &                     cjPosition, ParallelSPC, Overlap )
                ! Overlap detected
                IF( Overlap ) THEN
                  WRITE( *, "(5G0)" ) "Overlap detected in initial configuration between particles ", iParticle, " and ", &
                  &                   jParticle, "! Exiting..."
                  STOP
                END IF
              END IF
            END IF
          END DO
        END DO
      END IF
    END DO
  END DO

  ! Trajectory file
  IF( TrajectoryLogical ) THEN
    OPEN( Unit= 20, File= "Trajectories/"//TRIM( DateDescriptor )//"/"//TRIM( HourDescriptor )//"_traj_D" &
    &                     //TRIM( FileDescriptor(1) )//"_L"//TRIM( FileDescriptor(2) )//"_LD" &
    &                     //TRIM( FileDescriptor(3) )//"_λ"//TRIM( FileDescriptor(4) )//".xyz" )
    WRITE( 20, "(I5)" ) nParticles * 4
    DescriptorString = "(G0,8(G0.6,1X),G0.6,G0,2(G0.6,1X),G0.6,2G0)"
    WRITE( 20, DescriptorString ) 'Lattice="', BoxLength(1:9), '" Origin="', -0.5D0 * ( BoxLength(1) + BoxLength(4) + &
    &                             BoxLength(7) ), -0.5D0 * ( BoxLength(2) + BoxLength(5) + BoxLength(8) ), -0.5D0 * &
    &                             ( BoxLength(3) + BoxLength(6) + BoxLength(9) ), '" ', &
    &                             "Properties=species:S:1:pos:R:3:orientation:R:4:aspherical_shape:R:3"
    DO pParticle = 1, nParticles
      ! Initial configuration for OVITO
      DO cCylinder = 1, 4
        ! Position of cylinders
        WritePosition(1) = cPosition(1,cCylinder,pParticle)
        WritePosition(2) = cPosition(2,cCylinder,pParticle)
        WritePosition(3) = cPosition(3,cCylinder,pParticle)
        WRITE( 20, "(G0,10(' ',G0.9))") pParticle, WritePosition(1), WritePosition(2), WritePosition(3), pQuaternion(1,pParticle), &
        &                               pQuaternion(2,pParticle), pQuaternion(3,pParticle), pQuaternion(0,pParticle), &
        &                               0.5D0 * cDiameter, 0.5D0 * cDiameter, cLength
      END DO
    END DO
    FLUSH( 20 )
  END IF

  ! Ratio file (translation)
  OPEN( Unit= 30, File= "Ratio/Translation/"//TRIM( DateDescriptor )//"/"//TRIM( HourDescriptor )//"_ratio_D" &
  &                     //TRIM( FileDescriptor(1) )//"_L"//TRIM( FileDescriptor(2) )//"_LD" &
  &                     //TRIM( FileDescriptor(3) )//"_λ"//TRIM( FileDescriptor(4) )//".dat" )
  WRITE( 30, "(G0)" ) "Cycles,Ratio,MaxTranslationalDisplacement,AcceptanceRatioTranslation"
  FLUSH( 30 )

  ! Ratio file (rotation)
  OPEN( Unit= 40, File= "Ratio/Rotation/"//TRIM( DateDescriptor )//"/"//TRIM( HourDescriptor )//"_ratio_D" &
  &                     //TRIM( FileDescriptor(1) )//"_L"//TRIM( FileDescriptor(2) )//"_LD" &
  &                     //TRIM( FileDescriptor(3) )//"_λ"//TRIM( FileDescriptor(4) )//".dat" )
  WRITE( 40, "(G0)" ) "Cycles,Ratio,MaxAngularDisplacement,AcceptanceRatioRotation"
  FLUSH( 40 )

  ! Results file
  OPEN( Unit= 70, File= "Results/"//TRIM( DateDescriptor )//"/"//TRIM( HourDescriptor )//"_results_D" &
  &                     //TRIM( FileDescriptor(1) )//"_L"//TRIM( FileDescriptor(2) )//"_LD" &
  &                     //TRIM( FileDescriptor(3) )//"_λ"//TRIM( FileDescriptor(4) )//".dat" )
  WRITE( 70, "(G0)" ) "Cycles,OrientationalFieldEnergy,SquaredSineAngle,CumulativeSquaredSineAngle"
  FLUSH( 70 )

  ! ********************************************************************************************* !
  ! Monte Carlo parameters                                                                        !
  ! ********************************************************************************************* !
  MovementTranslationLogical   = .FALSE.                ! Translational move selector              (initial value)
  MovementRotationLogical      = .FALSE.                ! Rotational move selector                 (initial value)
  MaxTranslationalDisplacement = maxTranslationalDisplc ! Maximum translational displacement       (initial value)
  MaxAngularDisplacement       = maxRotationalDisplc    ! Maximum rotational displacement          (initial value)
  nAcceptanceTranslation       = 0                      ! Translational move acceptance counter    (initial value)
  nAcceptanceRotation          = 0                      ! Rotational move acceptance counter       (initial value)
  nMovementTranslationCounter  = 0                      ! Translational move counter               (initial value)
  nMovementRotationCounter     = 0                      ! Rotational move counter                  (initial value)
  cPositionMC(:,:,:)           = cPosition(:,:,:)       ! Position (cylinders)                     (initial value)
  pQuaternionMC(:,:)           = pQuaternion(:,:)       ! Quaternion algebra                       (initial value)
  pPositionMC(:,:)             = pPosition(:,:)         ! Position of particles                    (initial value)
  pOrientationMC(:,:)          = pOrientation(:,:)      ! Orientation of particles                 (initial value)
  OrientationalFieldEnergy     = 0.D0                   ! Orientational field potential            (initial value)
  CalculationType              = "NVT"                  ! Calculation type                         (initial value)
  AverageSquaredSineAngle      = 0.D0                   ! Average energy                           (initial value)
  nEnergySamples               = 0                      ! Number of energy samples                 (initial value)

  ! Potential energy of initial configuration
  SquaredSineAngle = 0.D0
  DO pParticle = 1, nParticles
    CossineOrientFieldAngle = DOT_PRODUCT( pOrientationMC(:,pParticle), OrientFieldVector(:) ) / ( DSQRT( DOT_PRODUCT( &
    &                         pOrientationMC(:,pParticle), pOrientationMC(:,pParticle) ) ) * DSQRT( DOT_PRODUCT( &
    &                         OrientFieldVector, OrientFieldVector ) ) )
    IF( CossineOrientFieldAngle >= 1.D0 ) CossineOrientFieldAngle = 1.D0
    IF( CossineOrientFieldAngle <= - 1.D0 ) CossineOrientFieldAngle = - 1.D0
    ! Increment potential energy
    SquaredSineAngle = SquaredSineAngle + DSIN( DACOS( CossineOrientFieldAngle ) ) * DSIN( DACOS( CossineOrientFieldAngle ) )
  END DO
  ! Total potential energy of the system
  OrientationalFieldEnergy = SquaredSineAngle * GaussLegendrePoints(qPoint)

  ! Simulation cycles
  DO Cycles = 1, nSimulationCycles

    ! Progress bar
    IF( MOD( Cycles, nSavingFrequency ) == 0 ) THEN
      CALL ProgressBar( Cycles, nSimulationCycles, CalculationType )
    END IF

    ! Random particle
    CALL Random_Number( RandomNumber )
    Particle = INT( RandomNumber * DBLE( nParticles ) ) + 1

    ! Assignment of previous configuration (microstate m)
    iOldPosition(:)    = pPositionMC(:,Particle)    ! Position
    cOldPosition(:,:)  = cPositionMC(:,:,Particle)  ! Position (cylinders)
    iOldQuaternion(:)  = pQuaternionMC(:,Particle)  ! Quaternion
    iOldOrientation(:) = pOrientationMC(:,Particle) ! Orientation

    ! Random move selection
    CALL Random_Number( RandomNumber )
    ! Translation criterion
    IF( RandomNumber <= TranslationalProbability ) THEN
      MovementTranslationLogical  = .TRUE.  ! Enable translation
      MovementRotationLogical     = .FALSE. ! Disable rotation
      nMovementTranslationCounter = nMovementTranslationCounter + 1
    ! Rotation criterion
    ELSE IF( RandomNumber > TranslationalProbability ) THEN
      MovementRotationLogical    = .TRUE.  ! Enable rotation
      MovementTranslationLogical = .FALSE. ! Disable translation
      nMovementRotationCounter   = nMovementRotationCounter + 1
    END IF

    ! Translational movement
    IF( MovementTranslationLogical ) THEN
      ! Random translation along x-axis
      CALL Random_Number( RandomNumber )
      iNewPosition(1) = iOldPosition(1) + ( ( 2.D0 * RandomNumber ) - 1.D0 ) * MaxTranslationalDisplacement ! Range [-MaxTranslationalDisplacement,MaxTranslationalDisplacement]
      ! Random translation along y-axis
      CALL Random_Number( RandomNumber )
      iNewPosition(2) = iOldPosition(2) + ( ( 2.D0 * RandomNumber ) - 1.D0 ) * MaxTranslationalDisplacement ! Range [-MaxTranslationalDisplacement,MaxTranslationalDisplacement]
      ! Random translation along z-axis
      CALL Random_Number( RandomNumber )
      iNewPosition(3) = iOldPosition(3) + ( ( 2.D0 * RandomNumber ) - 1.D0 ) * MaxTranslationalDisplacement ! Range [-MaxTranslationalDisplacement,MaxTranslationalDisplacement]
      ! Minimum Image Convention
      CALL MatrixVectorMultiplication( BoxLengthInverse, iNewPosition, ScalingDistanceUnitBox )
      ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
      CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, iNewPosition )
    ! No translational movement
    ELSE IF( .NOT. MovementTranslationLogical ) THEN
      iNewPosition(:) = iOldPosition(:)
    END IF

    ! Rotational movement
    IF( MovementRotationLogical ) THEN
      ! Random Composed Unit Quaternion
      CALL QuaternionCombination( iOldQuaternion, iNewQuaternion, MaxAngularDisplacement )
      ! Active transformation (rotation)
      CALL VectorRotation( zAxis, iNewQuaternion, iNewOrientation )
    ! No Rotation
    ELSE IF( .NOT. MovementRotationLogical ) THEN
      iNewQuaternion(:) = iOldQuaternion(:)
      iNewOrientation(:) = iOldOrientation(:)
    END IF

    ! Random position of cylinders (after translation or rotation)
    DO cCylinder = 1, 4
      ! Active transformation (translation)
      CALL VectorRotation( cPositionBasis(:,cCylinder), iNewQuaternion, cRotatedPositionBasis(:,cCylinder) )
      cNewPosition(:,cCylinder) = iNewPosition(:) + cRotatedPositionBasis(:,cCylinder)
    END DO

    ! Purely repulsive hard-core potential
    CALL ParticleCheckOverlap( Particle, iNewQuaternion, iNewOrientation, iNewPosition, cNewPosition, Overlap )

    ! Potential energy (before trial move)
    CossineOrientFieldAngle = DOT_PRODUCT( iOldOrientation(:), OrientFieldVector(:) ) / ( DSQRT( DOT_PRODUCT( iOldOrientation, &
    &                         iOldOrientation ) ) * DSQRT( DOT_PRODUCT( OrientFieldVector, OrientFieldVector ) ) )
    IF( CossineOrientFieldAngle >= 1.D0 ) CossineOrientFieldAngle = 1.D0
    IF( CossineOrientFieldAngle <= - 1.D0 ) CossineOrientFieldAngle = - 1.D0
    OldPotentialEnergy  = DSIN( DACOS( CossineOrientFieldAngle ) ) * DSIN( DACOS( CossineOrientFieldAngle ) )
    OldSquaredSineAngle = OldPotentialEnergy
    OldPotentialEnergy  = GaussLegendrePoints(qPoint) * OldPotentialEnergy
    ! Potential energy (after trial move)
    CossineOrientFieldAngle = DOT_PRODUCT( iNewOrientation(:), OrientFieldVector(:) ) / ( DSQRT( DOT_PRODUCT( iNewOrientation, &
    &                         iNewOrientation ) ) * DSQRT( DOT_PRODUCT( OrientFieldVector, OrientFieldVector ) ) )
    IF( CossineOrientFieldAngle >= 1.D0 ) CossineOrientFieldAngle = 1.D0
    IF( CossineOrientFieldAngle <= - 1.D0 ) CossineOrientFieldAngle = - 1.D0
    CossineOrientFieldAngle = DACOS( CossineOrientFieldAngle )
    NewPotentialEnergy  = DSIN( CossineOrientFieldAngle ) * DSIN( CossineOrientFieldAngle )
    NewSquaredSineAngle = NewPotentialEnergy
    NewPotentialEnergy  = GaussLegendrePoints(qPoint) * NewPotentialEnergy

    ! Difference in potential energy and squared sine angle
    PotentialEnergyDifference = NewPotentialEnergy - OldPotentialEnergy
    SquaredSineAngleDifference = NewSquaredSineAngle - OldSquaredSineAngle

    ! Metropolis criterion
    CALL Random_Number( RandomNumber )
    IF( ( PotentialEnergyDifference <= 0.D0 .OR. DEXP( - PotentialEnergyDifference ) >= RandomNumber ) .AND. &
    &   ( .NOT. Overlap ) ) THEN
      ! Update system configuration
      pPositionMC(:,Particle)    = iNewPosition(:)    ! Update position
      cPositionMC(:,:,Particle)  = cNewPosition(:,:)  ! Update position (cylinders)
      pQuaternionMC(:,Particle)  = iNewQuaternion(:)  ! Update quaternion
      pOrientationMC(:,Particle) = iNewOrientation(:) ! Update orientation
      ! Update energy of the system
      OrientationalFieldEnergy = OrientationalFieldEnergy + PotentialEnergyDifference
      ! Update squared sine angle (ensemble average)
      SquaredSineAngle = SquaredSineAngle + SquaredSineAngleDifference
      ! Update displacement counter
      IF( MovementTranslationLogical ) THEN
        nAcceptanceTranslation = nAcceptanceTranslation + 1 ! Translational move counter
      ELSE IF( MovementRotationLogical ) THEN
        nAcceptanceRotation = nAcceptanceRotation + 1 ! Rotational move counter
      END IF
    ! Move rejected
    ELSE
      ! Retrieve old configuration
      pPositionMC(:,Particle)    = iOldPosition(:)    ! Position
      cPositionMC(:,:,Particle)  = cOldPosition(:,:)  ! Position (cylinders)
      pQuaternionMC(:,Particle)  = iOldQuaternion(:)  ! Quaternion
      pOrientationMC(:,Particle) = iOldOrientation(:) ! Orientation
    END IF

    ! Adjustment of maximum displacement
    IF( Cycles <= nEquilCycles ) THEN ! During equilibration only
      ! Adjustment frequency (translation and rotation)
      IF( MOD( Cycles, nAdjustmentFrequency ) == 0 ) THEN
        ! Acceptance ratio (non-overlapping microstates over sampled microstates)
        IF( nMovementTranslationCounter > 0 ) THEN
          Ratio = DBLE( nAcceptanceTranslation ) / DBLE( nMovementTranslationCounter )
          ! Translational adjustment
          IF( Ratio <= AcceptanceRatioTranslation ) THEN
            MaxTranslationalDisplacement  = 0.95D0 * MaxTranslationalDisplacement
          ELSE
            MaxTranslationalDisplacement  = 1.05D0 * MaxTranslationalDisplacement
          END IF
          ! Ratio data
          WRITE( 30, "(G0,3(',',G0.9))" ) Cycles, Ratio, MaxTranslationalDisplacement, AcceptanceRatioTranslation
          FLUSH( 30 )
          ! Reset counter
          nAcceptanceTranslation = 0
          nMovementTranslationCounter = 0
          ! Avoid multiple turns
          IF( MaxTranslationalDisplacement >= 4.D0 * MAXVAL( BoxLength ) ) THEN
            MaxTranslationalDisplacement = MaxTranslationalDisplacement - 2.D0 * MAXVAL( BoxLength )
          END IF
        END IF
        ! Acceptance ratio (non-overlapping microstates over sampled microstates)
        IF( nMovementRotationCounter > 0 ) THEN
          Ratio = DBLE( nAcceptanceRotation ) / DBLE( nMovementRotationCounter )
          ! Rotation adjustment
          IF( Ratio <= AcceptanceRatioRotation ) THEN
            MaxAngularDisplacement = 0.95D0 * MaxAngularDisplacement
          ELSE
            MaxAngularDisplacement = 1.05D0 * MaxAngularDisplacement
          END IF
          ! Ratio data
          WRITE( 40, "(G0,3(',',G0.9))" ) Cycles, Ratio, MaxAngularDisplacement, AcceptanceRatioRotation
          FLUSH( 40 )
          ! Reset counter
          nAcceptanceRotation = 0
          nMovementRotationCounter  = 0
          ! Avoid multiple turns
          IF( MaxAngularDisplacement >= 4.D0 * cPi ) THEN
            MaxAngularDisplacement = MaxAngularDisplacement - 2.D0 * cPi
          END IF
        END IF
      END IF
    END IF

    ! Trajectory data
    IF( TrajectoryLogical ) THEN
      IF( MOD( Cycles, nSavingFrequency ) == 0 ) THEN
        WRITE( 20, "(I5)" ) nParticles * 4
        DescriptorString = "(G0,8(G0.6,1X),G0.6,G0,2(G0.6,1X),G0.6,2G0)"
        WRITE( 20, DescriptorString ) 'Lattice="', BoxLength(1:9), '" Origin="', -0.5D0 * ( BoxLength(1) + BoxLength(4) + &
        &                             BoxLength(7) ), -0.5D0 * ( BoxLength(2) + BoxLength(5) + BoxLength(8) ), -0.5D0 * &
        &                             ( BoxLength(3) + BoxLength(6) + BoxLength(9) ), '" ', &
        &                             "Properties=species:S:1:pos:R:3:orientation:R:4:aspherical_shape:R:3"
        DO pParticle = 1, nParticles
          ! Position of cylinders
          DO cCylinder = 1, 4
            WritePosition(1) = cPositionMC(1,cCylinder,pParticle)
            WritePosition(2) = cPositionMC(2,cCylinder,pParticle)
            WritePosition(3) = cPositionMC(3,cCylinder,pParticle)
            WRITE( 20, "(G0,10(' ',G0.9))" ) pParticle, WritePosition(1), WritePosition(2), WritePosition(3), &
            &                                pQuaternionMC(1,pParticle), pQuaternionMC(2,pParticle), pQuaternionMC(3,pParticle), &
            &                                pQuaternionMC(0,pParticle), 0.5D0 * cDiameter, 0.5D0 * cDiameter, cLength
          END DO
        END DO
        FLUSH( 20 )
      END IF
    END IF

    ! Results data
    IF( MOD( Cycles, nSavingFrequency ) == 0 ) THEN
      AverageSquaredSineAngle = AverageSquaredSineAngle + SquaredSineAngle
      nEnergySamples = nEnergySamples + 1
      ! Potential energy
      WRITE( 70, "(G0,3(',',G0.15))" ) Cycles, OrientationalFieldEnergy, SquaredSineAngle, AverageSquaredSineAngle
      FLUSH( 70 )
    END IF

    ! Reset accumulators
    IF( Cycles == nEquilCycles ) THEN
      AverageSquaredSineAngle = 0.D0
      nEnergySamples = 0
    END IF

  END DO

  WRITE( *, * ) " "
  WRITE( *, * ) " "

  ! Average squared sine angle
  AverageSquaredSineAngle = AverageSquaredSineAngle / ( DBLE( nEnergySamples ) )

  ! Total squared sine angle (ensemble average) of quadrature node
  sqSineAngleEnsembleAverage(qPoint) = AverageSquaredSineAngle

  ! Potential energy and number of non-overlapping configurations over sampled cycles
  WRITE( 80, "(G0,2(',',G0.15))" ) qPoint, GaussLegendrePoints(qPoint), sqSineAngleEnsembleAverage(qPoint)
  FLUSH( 80 )

  ! Output units                                         
  IF( TrajectoryLogical ) THEN
    CLOSE( 20 )
  END IF
  CLOSE( 30 )
  CLOSE( 40 )
  CLOSE( 70 )

END DO

! Output units 
CLOSE( 80 )

! Gauss-Legendre quadrature
QuadratureContributions(1)  = GaussLegendreWeights(1)  * (sqSineAngleEnsembleAverage(10) + sqSineAngleEnsembleAverage(11) )
QuadratureContributions(2)  = GaussLegendreWeights(2)  * (sqSineAngleEnsembleAverage(9)  + sqSineAngleEnsembleAverage(12) )
QuadratureContributions(3)  = GaussLegendreWeights(3)  * (sqSineAngleEnsembleAverage(8)  + sqSineAngleEnsembleAverage(13) )
QuadratureContributions(4)  = GaussLegendreWeights(4)  * (sqSineAngleEnsembleAverage(7)  + sqSineAngleEnsembleAverage(14) )
QuadratureContributions(5)  = GaussLegendreWeights(5)  * (sqSineAngleEnsembleAverage(6)  + sqSineAngleEnsembleAverage(15) )
QuadratureContributions(6)  = GaussLegendreWeights(6)  * (sqSineAngleEnsembleAverage(5)  + sqSineAngleEnsembleAverage(16) )
QuadratureContributions(7)  = GaussLegendreWeights(7)  * (sqSineAngleEnsembleAverage(4)  + sqSineAngleEnsembleAverage(17) )
QuadratureContributions(8)  = GaussLegendreWeights(8)  * (sqSineAngleEnsembleAverage(3)  + sqSineAngleEnsembleAverage(18) )
QuadratureContributions(9)  = GaussLegendreWeights(9)  * (sqSineAngleEnsembleAverage(2)  + sqSineAngleEnsembleAverage(19) )
QuadratureContributions(10) = GaussLegendreWeights(10) * (sqSineAngleEnsembleAverage(1)  + sqSineAngleEnsembleAverage(20) )

! Free energy change between the externally oriented fluid and the liquid-crystalline phase of interest 
FreeEnergyChangeAField = - CouplingParameterHalfWidth * SUM( QuadratureContributions ) / DBLE( nParticles )

! Summary
WRITE( *, "(G0)" ) "Writing log..."
WRITE( *, "(G0)" ) " "
WRITE( *, "(G0,G0.5)" ) "Free energy change between the externally oriented fluid and the LC phase of interest  is: ", &
&                       FreeEnergyChangeAField

! Results file
OPEN( Unit= 95, File= "Results/"//TRIM( DateDescriptor )//"/"//TRIM( HourDescriptor )//"_A2_D"//TRIM( FileDescriptor(1) )// &
&                     "_L"//TRIM( FileDescriptor(2) )//"_LD"//TRIM( FileDescriptor(3) )//"_λENERGY.dat" )
WRITE( 95, "(G0,G0.15)" ) "Free energy change between the solid and the interacting EC is: ", FreeEnergyChangeAField
FLUSH( 95 )
CLOSE( 95 )

RETURN

END SUBROUTINE Compute_AField

END MODULE ThermodynamicStep
