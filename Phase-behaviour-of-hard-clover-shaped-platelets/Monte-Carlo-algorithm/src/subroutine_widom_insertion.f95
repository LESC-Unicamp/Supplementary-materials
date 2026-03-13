! *********************************************************************************************** !
!   This function applies the Widom insertion method to calculate the chemical potential of the   !
!                                           fluid phase                                           !
! *********************************************************************************************** !
SUBROUTINE WidomInsertionMethod( CurrentBoxLength, CurrentBoxLengthInverse, BoltzmannEnergyFactor )

! Uses four modules: global variables, vector operations, quaternion operations, and overlap check
USE GlobalVariables
USE VectorOperations
USE QuaternionOperations
USE OverlapCheckSystem, ONLY: ParticleOverlapCheck

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 ) :: tParticle             ! Test particle
INTEGER( Kind= Int64 ) :: cCylinder             ! Counter (cylinder)
INTEGER( Kind= Int64 ) :: NonOverlappingCounter ! Counter (non-overlapping configurations)
INTEGER( Kind= Int64 ) :: wCycle                ! Counter (Widom cycles)

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 ) :: BoltzmannEnergyFactor ! Boltzmann factor
REAL( Kind= Real64 ) :: RandomAngle           ! Random angle

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: RandomVector            ! Random unit vector
REAL( Kind= Real64 ), DIMENSION( 3 )    :: ScalingDistanceUnitBox  ! Scaled position (unit box)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: wPosition               ! Random position (Widom insertion)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: wOrientation            ! Random orientation (Widom insertion)
REAL( Kind= Real64 ), DIMENSION( 9 )    :: CurrentBoxLength        ! Length (x,y,z) of simulation box
REAL( Kind= Real64 ), DIMENSION( 9 )    :: CurrentBoxLengthInverse ! Length (x,y,z) of simulation box (inverse)
REAL( Kind= Real64 ), DIMENSION( 0:3 )  :: wQuaternion             ! Random quaternion (Widom insertion)
REAL( Kind= Real64 ), DIMENSION( 3, 4 ) :: wcPosition              ! Random position of cylinders (Widom insertion)
REAL( Kind= Real64 ), DIMENSION( 3, 4 ) :: cRotatedPositionBasis   ! Relative distance

! LOGICAL VARIABLES
LOGICAL :: Overlap ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected

! Initialization
tParticle = 0
NonOverlappingCounter = 0

! Widom cycles
DO wCycle = 1, nWidomCycles
  ! Random x-position (unit box)
  CALL Random_Number( RandomNumber )
  ScalingDistanceUnitBox(1) = RandomNumber - 0.5D0
  ! Random y-position (unit box)
  CALL Random_Number( RandomNumber )
  ScalingDistanceUnitBox(2) = RandomNumber - 0.5D0
  ! Random z-position (unit box)
  CALL Random_Number( RandomNumber )
  ScalingDistanceUnitBox(3) = RandomNumber - 0.5D0
  ! Random position in the simulation box
  CALL MatrixVectorMultiplication( CurrentBoxLength, ScalingDistanceUnitBox, wPosition )
  ! Random rotation angle
  CALL Random_Number( RandomNumber )
  RandomAngle = ( 2.D0 * RandomNumber ) * cPi ! Range [0, 2π]
  ! Random vector generator
  CALL RandomVectorGenerator( RandomVector )
  ! Random quaternion
  wQuaternion(0) = DCOS( RandomAngle * 0.5D0 )                   ! Real part
  wQuaternion(1) = DSIN( RandomAngle * 0.5D0 ) * RandomVector(1) ! Imaginary part (Vector)
  wQuaternion(2) = DSIN( RandomAngle * 0.5D0 ) * RandomVector(2) ! Imaginary part (Vector)
  wQuaternion(3) = DSIN( RandomAngle * 0.5D0 ) * RandomVector(3) ! Imaginary part (Vector)
  ! Random orientation
  CALL VectorRotation( zAxis, wQuaternion, wOrientation )
  ! Random position of cylinders
  DO cCylinder = 1, 4
    ! Active transformation (translation)
    CALL VectorRotation( cPositionBasis(:,cCylinder), wQuaternion, cRotatedPositionBasis(:,cCylinder) )
    wcPosition(:,cCylinder) = wPosition(:) + cRotatedPositionBasis(:,cCylinder)
  END DO
  ! Check overlap
  CALL ParticleOverlapCheck( tParticle, wQuaternion, wOrientation, wPosition, wcPosition, CurrentBoxLength, &
  &                          CurrentBoxLengthInverse, Overlap )
  ! Counter of non-overlapping configurations
  IF( .NOT. Overlap ) THEN
    NonOverlappingCounter = NonOverlappingCounter + 1
  END IF
END DO

! Average Boltzmann factor
BoltzmannEnergyFactor = DBLE( NonOverlappingCounter ) / nWidomCycles

RETURN

END SUBROUTINE WidomInsertionMethod