MODULE InitialConfiguration

! Uses four modules: global variables, quaternion operations, vector operations, and overlap check
USE GlobalVariables
USE QuaternionOperations
USE VectorOperations
USE OverlapCheckSystem, ONLY: ParticleOverlapCheck

IMPLICIT NONE

CONTAINS

! *********************************************************************************************** !
!          This subroutine allows the user to choose the initial molecular configuration          !
! *********************************************************************************************** !
SUBROUTINE Configuration_Selection(  )

IMPLICIT NONE

! Initialization
ConfigurationLogical(:) = .FALSE.

! Molecular properties
OPEN( Unit= 100, File= "ini_initial_configuration.ini" )

! Molecular structure
READ( 100, * ) Dummy, ConfigurationInquiry
CALL ToUpper( ConfigurationInquiry, LEN_TRIM( ConfigurationInquiry ), ConfigurationInquiry )

! Initial configuration
IF( ConfigurationInquiry == "PB" ) THEN
  ! Packed box
  ConfigurationLogical(1) = .TRUE.
ELSE IF( ConfigurationInquiry == "FB" ) THEN
  ! Floppy box
  ConfigurationLogical(2) = .TRUE.
END IF

CLOSE( 100 )

RETURN

END SUBROUTINE Configuration_Selection

! *********************************************************************************************** !
!         This subroutine allows the user to choose the surface geometry of the molecules         !
! *********************************************************************************************** !
SUBROUTINE Geometry_Selection(  )

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 ) :: FrameLeft = 0  ! Box frame dimension
INTEGER( Kind= Int64 ) :: FrameRight = 0 ! Box frame dimension

! *********************************************************************************************** !
!  FOUR-LEAF CLOVER GEOMETRIES                                                                    !
!                                                                                                 !
!  The geometry of the cross sections consists of four circles and are constructed as follows:    !
!                                                                                                 !
!  (1) = RULES: A circle intersects the centers of the two nearest neighboring circles            !
!               A circle is equidistant from the two nearest neighboring circles                  !
!                                                                                                 !
!  (2) = RULES: All four circles intersect the center of mass of the composed geometry            !
!               A circle is equidistant from the two nearest neighboring circles                  !
! *********************************************************************************************** !
GeometrySelection(:) = .FALSE.

! Packed box configuration
IF( ConfigurationLogical(1) .OR. Path == 1 ) THEN
  OPEN( Unit= 100, File= "ini_initial_configuration.ini", Action= "READ" )
  READ( 100, * ) Dummy, Dummy
  READ( 100, * ) Dummy, pGeometry
  IF( pGeometry == 1 ) THEN
    GeometrySelection(1) = .TRUE.
  ELSE IF( pGeometry == 2 ) THEN
    GeometrySelection(2) = .TRUE.
  END IF
  CLOSE( 100 )
! Floppy box configuration
ELSE IF( ConfigurationLogical(2) .AND. Path /= 1 ) THEN
  INQUIRE( File= "Floppy-box/", Exist= FileExist )
  IF( .NOT. FileExist ) CALL SYSTEM( "mkdir Floppy-box" )
  INQUIRE( File= "Floppy-box/fbox.dat", Exist= FileExist )
  WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 70 )//CH_UR
  FrameLeft  = (35 - FLOOR( REAL( 15 ) / 2.D0 ) )
  FrameRight = (35 - CEILING( REAL( 15 ) / 2.D0 ) )
  WRITE( *, "(2G0)" ) CH_VS//REPEAT( " ", FrameLeft )//"FLOPPY-BOX DATA"//REPEAT( " ", FrameRight )//CH_VS
  WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 70 )//CH_BR
  IF( .NOT. FileExist ) THEN
    WRITE( *, "(G0)" ) "Floppy-box data file not found! Exiting..."
    WRITE( *, "(G0)") " "
    STOP
  ELSE IF( FileExist ) THEN
    WRITE( *, "(G0)" ) "Floppy-box data file found! Resuming..."
    WRITE( *, "(G0)" ) " "
  END IF
  ! Molecular properties
  OPEN( Unit= 10, File= "Floppy-box/fbox.dat", Action= "READ" )
  ! Cross-section geometry
  READ( 10, * ) Dummy, pGeometry
  GeometrySelection(:) = .FALSE.
  ! Rule (1)
  IF( pGeometry == 1 ) THEN
    GeometrySelection(1) = .TRUE.
  ! Rule (2)
  ELSE IF( pGeometry == 2 ) THEN
    GeometrySelection(2) = .TRUE.
  END IF
  ! Cylindrical diameter
  READ( 10, * ) Dummy, cDiameter ! Å
  ! Cylindrical length
  READ( 10, * ) Dummy, cLength ! Å
  ! Aspect ratio of cylinders
  cAspectRatio = cLength / cDiameter
  CLOSE( 10 )
END IF

RETURN

END SUBROUTINE Geometry_Selection

! *********************************************************************************************** !
!               Subroutine to calculate the surface area of the non-convex geometry               !
!       The non-convex geometry is composed of four identical cylinders, such that [1] the        !
! circumference of a cylinder passes through the centers of the two nearest cylinders or [2] the  !
!           circumference of all cylinders passes through the molecular center of mass            !
!                                  (common intersection point).                                   !
!   In the case [1], the distance between the centers of the circumferences of the two nearest    !
!     cylinders is the radius of one cylinder, and the width of the non-convex body is 1.5D,      !
!                            where D is the diameter of the cylinder.                             !
!   In the case [2], the distance between the centers of the circumferences of the two nearest    !
!      cylinders is the hypotenuse (hyp) of an isosceles triangle whose side is equal to the      !
!            radius of one cylinder, and the width of the non-convex body is D + hyp.             !
! *********************************************************************************************** !
SUBROUTINE CrossSectionArea( pSurfaceArea )

! Uses one module: global variables
USE GlobalVariables

IMPLICIT NONE

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 ) :: aCircle, aIntersection2Circles, aIntersection4Circles ! Area of the circles and their intersections
REAL( Kind= Real64 ) :: pSurfaceArea                                          ! Surface area

! Geometry [1]
IF( GeometrySelection(1) ) THEN
  ! Area of circles
  aCircle = ( ( cDiameter * cDiameter ) / 4.D0 ) * cPi
  ! Area of the intersection between two circles
  aIntersection2Circles = ( ( cDiameter * cDiameter ) / 4.D0 ) * ( ( 2.D0 * cPi / 3.D0 ) - ( DSQRT( 3.D0 ) / 2.D0 ) )
  ! Area of the intersection between four circles
  aIntersection4Circles = ( ( cDiameter * cDiameter ) / 4.D0 ) * ( ( cPi / 3.D0 ) - ( DSQRT( 3.D0 ) ) + 1.D0 )
  ! Surface area
  pSurfaceArea = (4.D0 * aCircle ) - ( 4.D0 * aIntersection2Circles ) + ( aIntersection4Circles )
! Geometry [2]
ELSE IF( GeometrySelection(2) ) THEN
  ! Area of circles
  aCircle = ( ( cDiameter * cDiameter ) / 4.D0 ) * cPi
  ! Area of the intersection between two circles
  aIntersection2Circles = ( ( cDiameter * cDiameter ) / 4.D0 ) * ( cPi / 2.D0 - 1.D0 )
  ! Surface area
  pSurfaceArea = 4.D0 * aCircle - 4.D0 * aIntersection2Circles
END IF

RETURN

END SUBROUTINE CrossSectionArea

! *********************************************************************************************** !
!            This subroutine allocates particles according to the PB configuration                !
! *********************************************************************************************** !
SUBROUTINE PackedBoxConfiguration(  )

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 ) :: iParticle = 0, jParticle = 0    ! Counters (particles)
INTEGER( Kind= Int64 ) :: kParticle = 0, pCounter = 0     ! Counters (particles)
INTEGER( Kind= Int64 ) :: iRow = 0                        ! Counter (row)
INTEGER( Kind= Int64 ) :: cCylinder = 0                   ! Counter (cylinders)
INTEGER( Kind= Int64 ) :: nAcceptanceTranslation = 0      ! Move acceptance counter: Translation
INTEGER( Kind= Int64 ) :: nAcceptanceRotation = 0         ! Move acceptance counter: Rotation
INTEGER( Kind= Int64 ) :: nMovementTranslationCounter = 0 ! Move counter (Translation)
INTEGER( Kind= Int64 ) :: nMovementRotationCounter = 0    ! Move counter (Rotation)
INTEGER( Kind= Int64 ) :: Cycles = 0                      ! Cycles
INTEGER( Kind= Int64 ) :: FrameLeft = 0                   ! Box frame dimension
INTEGER( Kind= Int64 ) :: FrameRight = 0                  ! Box frame dimension

! INTEGER VARIABLES (ARRAY)
INTEGER( Kind= Int64 ), DIMENSION( 3 ) :: nParticlesRow = 0 ! Number of particles along x, y, and z directions

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 ) :: pWidth = 0.D0                       ! Width of the molecular geometry
REAL( Kind= Real64 ) :: Ratio = 0.D0                        ! Acceptance ratio (simulation)
REAL( Kind= Real64 ) :: MaxTranslationalDisplacement = 0.D0 ! Maximum displacement [+/-] (Translation)
REAL( Kind= Real64 ) :: MaxAngularDisplacement = 0.D0       ! Maximum displacement [+/-] (Rotation)

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: ScalingDistanceUnitBox = 0.D0                  ! Scaling factor
REAL( Kind= Real64 ), DIMENSION( 3 )    :: LengthIncrement = 0.D0                         ! Length
REAL( Kind= Real64 ), DIMENSION( 3 )    :: iOldOrientation = 0.D0, iNewOrientation = 0.D0 ! Orientation (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: iOldPosition = 0.D0, iNewPosition = 0.D0       ! Position of the center of mass (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( 3, 4 ) :: cOldPosition = 0.D0, cNewPosition = 0.D0       ! Position of cylinders (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( 3, 4 ) :: cRotatedPositionBasis = 0.D0                   ! Relative distance (after a trial move)
REAL( Kind= Real64 ), DIMENSION( 0: 3 ) :: iOldQuaternion = 0.D0, iNewQuaternion = 0.D0   ! Quaternion (before/after a trial move)

! LOGICAL VARIABLES
LOGICAL :: Overlap = .FALSE.                    ! Detects overlap between two particles (Vega-Lago) : TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: MovementRotationLogical = .FALSE.    ! Rotation move selection                           : TRUE = movement selected; FALSE = movement not selected
LOGICAL :: MovementTranslationLogical = .FALSE. ! Translation movement selection                    : TRUE = movement selected; FALSE = movement not selected

! CHARACTER STRINGS
CHARACTER( LEN= 03 ) :: EnsembleType = ' ' ! Ensemble type

! Summary
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 70 )//CH_UR
FrameLeft  = (35 - FLOOR( REAL( 24 ) / 2.D0 ) )
FrameRight = (35 - CEILING( REAL( 24 ) / 2.D0 ) )
WRITE( *, "(2G0)" ) CH_VS//REPEAT( " ", FrameLeft )//"PACKED BOX CONFIGURATION"//REPEAT( " ", FrameRight )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 70 )//CH_BR

! Convert degrees to radians
QuaternionAngle = QuaternionAngle * cPi / 180.D0

! Quaternion algebra
pQuaternion(0,:) = DCOS( QuaternionAngle * 0.5D0 )            ! Real part
pQuaternion(1,:) = DSIN( QuaternionAngle * 0.5D0 ) * zAxis(1) ! Imaginary part (Vector)
pQuaternion(2,:) = DSIN( QuaternionAngle * 0.5D0 ) * zAxis(2) ! Imaginary part (Vector)
pQuaternion(3,:) = DSIN( QuaternionAngle * 0.5D0 ) * zAxis(3) ! Imaginary part (Vector)

! Volume of the simulation box (initial estimative)
BoxVolume = 1.D0

! Status
WRITE( *, "(G0)" ) " "
WRITE( *, "(G0,G0.5,G0)" ) "Initial volume estimative: ", BoxVolume, "Å³"
WRITE( *, "(G0)" ) " "

! Initialization
nParticlesRow(:) = 1

! Molecular width
IF( GeometrySelection(1) ) THEN ! Geometry [1]
  pWidth = 1.5D0 * cDiameter ! Å
ELSE IF( GeometrySelection(2) ) THEN ! Geometry [2]
  pWidth = ( 1.D0 + 0.5D0 * DSQRT( 2.D0 ) ) * cDiameter ! Å
END IF

! Optimization of the volume of the simulation box
OptimizationLoop: DO
  ! Number density
  Density  = ( DBLE( nParticles ) / BoxVolume ) ! Å⁻³
  ! Packing fraction
  PackingFraction  = Density * pVolume
  ! Cubic box (only diagonal elements)
  BoxLength(:) = 0.D0
  BoxLength(1) = BoxVolume ** ( 1.D0 / 3.D0 ) ! Å
  BoxLength(5) = BoxVolume ** ( 1.D0 / 3.D0 ) ! Å
  BoxLength(9) = BoxVolume ** ( 1.D0 / 3.D0 ) ! Å
  ! Structure Validation
  IF( ( nParticlesRow(1) * nParticlesRow(2) * nParticlesRow(3) ) < nParticles ) THEN
    ! Update maximum number of particles along x, y and z directions
    nParticlesRow(1) = FLOOR( BoxLength(1) / pWidth )
    nParticlesRow(2) = FLOOR( BoxLength(5) / pWidth )
    nParticlesRow(3) = FLOOR( BoxLength(9) / cLength )
    ! Attempt to increase box volume (new estimative)
    BoxVolume = BoxVolume * 1.01D0 ! Å³
    CYCLE OptimizationLoop
  END IF
  ! Length increment
  LengthIncrement(1) = ( BoxLength(1) - ( pWidth * nParticlesRow(1) ) ) / nParticlesRow(1)  ! Å
  LengthIncrement(2) = ( BoxLength(5) - ( pWidth * nParticlesRow(2) ) ) / nParticlesRow(2)  ! Å
  LengthIncrement(3) = ( BoxLength(9) - ( cLength * nParticlesRow(3) ) ) / nParticlesRow(3) ! Å
  ! Avoid super compact configurations
  DO iRow = 1, 3
    IF( cAspectRatio <= 1.D0 ) THEN
      IF( LengthIncrement(iRow) < 1.D-1 * cLength ) THEN
        BoxVolume = BoxVolume * 1.01D0
        CYCLE OptimizationLoop
      END IF
    ELSE IF( cAspectRatio > 1.D0 ) THEN
      IF( LengthIncrement(iRow) < 1.D-1 / cLength ) THEN
        BoxVolume = BoxVolume * 1.01D0
        CYCLE OptimizationLoop
      END IF
    END IF
  END DO
  EXIT OptimizationLoop
END DO OptimizationLoop

! Make it N times bigger
BoxVolume = BoxFactor * BoxVolume ! Å³
! New number density
Density  = ( DBLE( nParticles ) / BoxVolume ) ! Å⁻³
! New packing fraction
PackingFraction  = Density * pVolume
! New cubic box (only diagonal elements)
BoxLength(:) = 0.D0
BoxLength(1) = BoxVolume ** ( 1.D0 / 3.D0 ) ! Å
BoxLength(5) = BoxVolume ** ( 1.D0 / 3.D0 ) ! Å
BoxLength(9) = BoxVolume ** ( 1.D0 / 3.D0 ) ! Å
! New length increment
LengthIncrement(1) = ( BoxLength(1) - ( ( pWidth ) * nParticlesRow(1) ) ) / nParticlesRow(1) ! Å
LengthIncrement(2) = ( BoxLength(5) - ( ( pWidth ) * nParticlesRow(2) ) ) / nParticlesRow(2) ! Å
LengthIncrement(3) = ( BoxLength(9) - ( ( cLength ) * nParticlesRow(3) ) ) / nParticlesRow(3)     ! Å

! Status
WRITE( *, "(G0,G0.5,G0)" ) "Final volume estimative: ", BoxVolume, "Å³" ! Å³
WRITE( *, "(G0)" ) " "
WRITE( *, "(G0)", Advance= "NO" ) "Positioning particles..."

! Position of nonconvex particles (centers of mass)
pCounter = 1
ParticlePlacement: DO iParticle = 1, nParticlesRow(1)
  DO jParticle = 1, nParticlesRow(2)
    DO kParticle = 1, nParticlesRow(3)
      ! Positioning (centers of mass)
      pPosition(1,pCounter) = DBLE( iParticle - 1 ) * ( pWidth + LengthIncrement(1) ) + 1.D-1  ! Å
      pPosition(2,pCounter) = DBLE( jParticle - 1 ) * ( pWidth + LengthIncrement(2) ) + 1.D-1  ! Å
      pPosition(3,pCounter) = DBLE( kParticle - 1 ) * ( cLength + LengthIncrement(3) ) + 1.D-1 ! Å
      ! Cylinders
      CALL VectorRotation( cPositionBasis(:,1), pQuaternion(:,pCounter), cRotatedPositionBasis(:,1) )
      cPosition(:,1,pCounter) = pPosition(:,pCounter) + cRotatedPositionBasis(:,1) ! First quarter
      CALL VectorRotation( cPositionBasis(:,2), pQuaternion(:,pCounter), cRotatedPositionBasis(:,2) )
      cPosition(:,2,pCounter) = pPosition(:,pCounter) + cRotatedPositionBasis(:,2) ! Second quarter
      CALL VectorRotation( cPositionBasis(:,3), pQuaternion(:,pCounter), cRotatedPositionBasis(:,3) )
      cPosition(:,3,pCounter) = pPosition(:,pCounter) + cRotatedPositionBasis(:,3) ! Third quarter
      CALL VectorRotation( cPositionBasis(:,4), pQuaternion(:,pCounter), cRotatedPositionBasis(:,4) )
      cPosition(:,4,pCounter) = pPosition(:,pCounter) + cRotatedPositionBasis(:,4) ! Fourth quarter
      ! Increment
      pCounter = pCounter + 1
      ! Condition
      IF( pCounter > nParticles ) THEN
        EXIT ParticlePlacement
      END IF
    END DO
  END DO
END DO ParticlePlacement

! Inverse matrix of the box
CALL InverseMatrixCofactorVec( BoxLength, BoxLengthInverse, BoxVolume )

! Centralizing the simulation box at origin of the coordinate system (0, 0, 0)
DO iParticle = 1, nParticles
  CALL MatrixVectorMultiplication( BoxLengthInverse, pPosition(:,iParticle), ScalingDistanceUnitBox )
  ScalingDistanceUnitBox = ScalingDistanceUnitBox - 0.5D0
  CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, pPosition(:,iParticle) )
  DO cCylinder = 1, 4
    CALL MatrixVectorMultiplication( BoxLengthInverse, cPosition(:,cCylinder,iParticle), ScalingDistanceUnitBox )
    ScalingDistanceUnitBox = ScalingDistanceUnitBox - 0.5D0
    CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, cPosition(:,cCylinder,iParticle) )
  END DO
END DO

! Active transformation
DO iParticle = 1, nParticles
  CALL VectorRotation( zAxis, pQuaternion(:,iParticle), pOrientation(:,iParticle) )
END DO

! Status
WRITE( *, "(G0)", Advance= "YES" ) " Done!"
WRITE( *, "(G0)", Advance= "YES" ) " "

! ********************************************************************************************* !
! Initialization (Randomize System Configuration)                                               !
! ********************************************************************************************* !
MovementTranslationLogical   = .FALSE.                          ! Translational move selector           (initial value)
MovementRotationLogical      = .FALSE.                          ! Rotational move selector              (initial value)
MaxTranslationalDisplacement = UserMaxTranslationalDisplacement ! Maximum translational displacement    (initial value)
MaxAngularDisplacement       = UserMaxRotationalDisplacement    ! Maximum rotational displacement       (initial value)
nAcceptanceTranslation       = 0                                ! Translational move acceptance counter (initial value)
nAcceptanceRotation          = 0                                ! Rotational move acceptance counter    (initial value)
nMovementTranslationCounter  = 0                                ! Translational move counter            (initial value)
nMovementRotationCounter     = 0                                ! Rotational move counter               (initial value)
cPositionMC(:,:,:)           = cPosition(:,:,:)                 ! Position (cylinders)                  (initial value)
pQuaternionMC(:,:)           = pQuaternion(:,:)                 ! Quaternion algebra                    (initial value)
pPositionMC(:,:)             = pPosition(:,:)                   ! Position of particles                 (initial value)
pOrientationMC(:,:)          = pOrientation(:,:)                ! Orientation of particles              (initial value)
EnsembleType                 = "NVT"                            ! Ensemble type                         (initial value)
Cycles                       = 0                                ! Number of cycles                      (initial value)

! Summary
IF( MaxSimulationCyclesInit > 0 ) THEN
  WRITE( *, "(3G0)" ) "Starting an NVT-Monte Carlo simulation with ", MaxSimulationCyclesInit, &
  &                   " cycle(s) to randomly distribute particles..."
  WRITE( *, "(G0)" ) " "
END IF

! Displace particles as long as the number of cycles is lower than the specified value
DO WHILE( Cycles < MaxSimulationCyclesInit )

  ! Progress bar
  CALL Progress_Bar( Cycles, MaxSimulationCyclesInit, EnsembleType )

  ! Loop over all particles
  DO iParticle = 1, nParticles

    ! Assignment of previous configuration (Microstate m)
    iOldPosition(:)    = pPositionMC(:,iParticle)    ! Position
    iOldQuaternion(:)  = pQuaternionMC(:,iParticle)  ! Quaternion
    cOldPosition(:,:)  = cPositionMC(:,:,iParticle)  ! Position (cylinders)
    iOldOrientation(:) = pOrientationMC(:,iParticle) ! Orientation

    ! Pseudorandom number generator (uniform distribution) (see 'subroutines' code for more details)
    CALL Random_Number( RandomNumber )

    ! Translation criterion
    IF( RandomNumber < TranslationalProbability ) THEN
      MovementTranslationLogical  = .TRUE.  ! Enable translation
      MovementRotationLogical     = .FALSE. ! Disable rotation
      nMovementTranslationCounter = nMovementTranslationCounter + 1 ! Increment move counter
    ! Rotation criterion
    ELSE IF( RandomNumber >= TranslationalProbability ) THEN
      MovementRotationLogical    = .TRUE.  ! Enable rotation
      MovementTranslationLogical = .FALSE. ! Disable translation
      nMovementRotationCounter   = nMovementRotationCounter + 1 ! Increment move counter
    END IF

    ! Translation
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
      ! Minimum Image Convention (see Allen and Tildesley, 2nd Edition (2017), pages 35-45)
      CALL MatrixVectorMultiplication( BoxLengthInverse, iNewPosition, ScalingDistanceUnitBox )
      ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
      CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, iNewPosition )
    ! No Translation
    ELSE IF( .NOT. MovementTranslationLogical ) THEN
      iNewPosition(:) = iOldPosition(:)
    END IF

    ! Rotation
    IF( MovementRotationLogical ) THEN
      ! Random Composed Unit Quaternion (see 'Subroutines' code for more details)
      CALL QuaternionCombination( iOldQuaternion, iNewQuaternion, MaxAngularDisplacement )
      ! Active transformation (rotation)
      CALL VectorRotation( zAxis, iNewQuaternion, iNewOrientation )
    ! No Rotation
    ELSE IF( .NOT. MovementRotationLogical ) THEN
      iNewQuaternion(:)  = iOldQuaternion(:)
      iNewOrientation(:) = iOldOrientation(:)
    END IF

    ! Random position of cylinders (after translation or rotation)
    DO cCylinder = 1, 4
      ! Active transformation (translation)
      CALL VectorRotation( cPositionBasis(:,cCylinder), iNewQuaternion, cRotatedPositionBasis(:,cCylinder) )
      cNewPosition(:,cCylinder) = iNewPosition(:) + cRotatedPositionBasis(:,cCylinder)
    END DO

    ! Overlap Check
    CALL ParticleOverlapCheck( iParticle, iNewQuaternion, iNewOrientation, iNewPosition, cNewPosition, BoxLength, &
    &                          BoxLengthInverse, Overlap )

    ! Acceptance Criterion
    IF( .NOT. Overlap ) THEN
      ! Update system configuration
      pPositionMC(:,iParticle)    = iNewPosition(:)    ! Update position
      cPositionMC(:,:,iParticle)  = cNewPosition(:,:)  ! Update position (cylinders)
      pQuaternionMC(:,iParticle)  = iNewQuaternion(:)  ! Update quaternion
      pOrientationMC(:,iParticle) = iNewOrientation(:) ! Update orientation
      ! Update displacement counter
      IF( MovementTranslationLogical ) THEN
        nAcceptanceTranslation = nAcceptanceTranslation + 1 ! Translational move counter
      ELSE IF( MovementRotationLogical ) THEN
        nAcceptanceRotation = nAcceptanceRotation + 1 ! Rotational move counter
      END IF
    ELSE
      ! Retrieve old configuration
      pPositionMC(:,iParticle)    = iOldPosition(:)    ! Position
      cPositionMC(:,:,iParticle)  = cOldPosition(:,:)  ! Position (cylinders)
      pQuaternionMC(:,iParticle)  = iOldQuaternion(:)  ! Quaternion
      pOrientationMC(:,iParticle) = iOldOrientation(:) ! Orientation
    END IF

  END DO

  ! Adjustment of maximum translation
  IF( MOD( Cycles, nAdjustmentMovementFrequency ) == 0 .AND. nMovementTranslationCounter > 0 ) THEN
    ! Acceptance ratio (non-overlapping microstates over sampled microstates)
    Ratio = DBLE( nAcceptanceTranslation ) / DBLE( nMovementTranslationCounter )
    ! Translational adjustment
    IF( Ratio <= 0.5D0 ) THEN
      MaxTranslationalDisplacement = 0.95D0 * MaxTranslationalDisplacement
    ELSE
      MaxTranslationalDisplacement = 1.05D0 * MaxTranslationalDisplacement
    END IF
    ! Avoid large translations
    IF( ( MaxTranslationalDisplacement >= 2.D0 * MAXVAL( BoxLength ) ) ) THEN
      MaxTranslationalDisplacement = MaxTranslationalDisplacement - MAXVAL( BoxLength )
    END IF
    ! Reset counter
    nMovementTranslationCounter = 0
    nAcceptanceTranslation      = 0
  END IF

  ! Adjustment of maximum rotation
  IF( MOD( Cycles, nAdjustmentMovementFrequency ) == 0 .AND. nMovementRotationCounter > 0 ) THEN
    ! Acceptance ratio (non-overlapping microstates over sampled microstates)
    Ratio = DBLE( nAcceptanceRotation ) / DBLE( nMovementRotationCounter )
    ! Rotational adjustment
    IF( Ratio <= 0.5D0 ) THEN
      MaxAngularDisplacement = 0.95D0 * MaxAngularDisplacement
    ELSE
      MaxAngularDisplacement = 1.05D0 * MaxAngularDisplacement
    END IF
    ! Avoid 4π-rotations (multiple turns)
    IF( ( MaxAngularDisplacement >= 4.D0 * cPi ) ) THEN
      MaxAngularDisplacement = MaxAngularDisplacement - 2.0D0 * cPi
    END IF
    ! Reset counter
    nAcceptanceRotation = 0
    nMovementRotationCounter  = 0
  END IF

  ! Update number of cycles
  Cycles = Cycles + 1

  ! Progress bar
  CALL Progress_Bar( Cycles, MaxSimulationCyclesInit, EnsembleType )

END DO

IF( MaxSimulationCyclesInit > 0 ) THEN
  WRITE( *, * ) " "
  WRITE( *, * ) " "
END IF

! Configuration update
pPosition(:,:)    = pPositionMC(:,:)   ! Update position
cPosition(:,:,:)  = cPositionMC(:,:,:) ! Update position (cylinders)
pQuaternion(:,:)  = pQuaternionMC(:,:)   ! Update quaternion
pOrientation(:,:) = pOrientationMC(:,:)   ! Update orientation

RETURN

END SUBROUTINE PackedBoxConfiguration

! *********************************************************************************************** !
!            This subroutine allocates particles according to the FB configuration                !
! *********************************************************************************************** !
SUBROUTINE FloppyBoxConfiguration(  )

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 ) :: iCell = 0, jCell = 0, kCell = 0              ! Counters (unit cells)
INTEGER( Kind= Int64 ) :: ParticleCell = 0, Particle = 0, pCounter = 0 ! Counter (particles)
INTEGER( Kind= Int64 ) :: iParticle = 0                                ! Counters (particles)
INTEGER( Kind= Int64 ) :: pCell = 0                                    ! Number of particles in unit cell
INTEGER( Kind= Int64 ) :: nCells = 0                                   ! Number of unit cells
INTEGER( Kind= Int64 ) :: FrameLeft = 0                                ! Box frame dimension
INTEGER( Kind= Int64 ) :: FrameRight = 0                               ! Box frame dimension

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 )                    :: CubicRootCheck         ! Cube root verification
REAL( Kind= Real64 )                    :: CellLengthCorrection   ! Scalar density-dependent factor to correct length of unit cell
REAL( Kind= Real64 )                    :: CellVolume             ! Cell volume

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 9 )    :: CellLength             ! Length of unit cell (triclinic)
REAL( Kind= Real64 ), DIMENSION( 9 )    :: CellLengthInverse      ! Inverse of length of unit cell (triclinic)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: ScalingDistanceUnitBox ! Scaling factor (unit cubic box)
REAL( Kind= Real64 ), DIMENSION( 3, 4 ) :: cRotatedPositionBasis  ! Relative distance (cylinders)

! REAL VARIABLES (ARRAY,ALLOCATABLE)
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE    :: pPositionCell   ! Position of particles in the unit cell
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE    :: pQuaternionCell ! Orientation (quaternion) of particles in the unit cell

! Floppy-box data
OPEN( Unit= 10, File= "Floppy-box/fbox.dat", Action= "READ" )
! Skip
READ( 10, * )
READ( 10, * )
READ( 10, * )
READ( 10, * )
! Cell length
READ( 10, * ) Dummy, CellLength(1:9)
! Cell volume
READ( 10, * ) Dummy, CellVolume
! Skip
READ( 10, * )
READ( 10, * )
! Number of particles in unit cell
READ( 10, * ) Dummy, pCell
! Allocation
ALLOCATE( pPositionCell(3,pCell) )
ALLOCATE( pQuaternionCell(0:3,pCell) )
! Skip
READ( 10, * )
READ( 10, * )
READ( 10, * )
DO Particle = 1, pCell
  ! Equilibrium positions and orientations (quaternions) in the unit cell
  READ( 10, * ) Dummy, pPositionCell(1:3,Particle), pQuaternionCell(0:3,Particle)
END DO
CLOSE( 10 )

! Unrotated reference position (cylinders)
IF( GeometrySelection(1) ) THEN ! Geometry [1]
  ! First quarter
  cPositionBasis(1,1) = 0.25D0 * cDiameter
  cPositionBasis(2,1) = 0.25D0 * cDiameter
  cPositionBasis(3,1) = 0.D0
  ! Second quarter
  cPositionBasis(1,2) = -0.25D0 * cDiameter
  cPositionBasis(2,2) = 0.25D0 * cDiameter
  cPositionBasis(3,2) = 0.D0
  ! Third quarter
  cPositionBasis(1,3) = -0.25D0 * cDiameter
  cPositionBasis(2,3) = -0.25D0 * cDiameter
  cPositionBasis(3,3) = 0.D0
  ! Fourth quarter
  cPositionBasis(1,4) = 0.25D0 * cDiameter
  cPositionBasis(2,4) = -0.25D0 * cDiameter
  cPositionBasis(3,4) = 0.D0
ELSE IF( GeometrySelection(2) ) THEN ! Geometry [2]
  ! First quarter
  cPositionBasis(1,1) = 0.25D0 * DSQRT( 2.D0 ) * cDiameter
  cPositionBasis(2,1) = 0.25D0 * DSQRT( 2.D0 ) * cDiameter
  cPositionBasis(3,1) = 0.D0
  ! Second quarter
  cPositionBasis(1,2) = -0.25D0 * DSQRT( 2.D0 ) * cDiameter
  cPositionBasis(2,2) = 0.25D0 * DSQRT( 2.D0 ) * cDiameter
  cPositionBasis(3,2) = 0.D0
  ! Third quarter
  cPositionBasis(1,3) = -0.25D0 * DSQRT( 2.D0 ) * cDiameter
  cPositionBasis(2,3) = -0.25D0 * DSQRT( 2.D0 ) * cDiameter
  cPositionBasis(3,3) = 0.D0
  ! Fourth quarter
  cPositionBasis(1,4) = 0.25D0 * DSQRT( 2.D0 ) * cDiameter
  cPositionBasis(2,4) = -0.25D0 * DSQRT( 2.D0 ) * cDiameter
  cPositionBasis(3,4) = 0.D0
END IF

! Validation of the number of particles
ParticleValidation: DO
  CubicRootCheck = ( DBLE( nParticles ) / DBLE( pCell ) ) ** ( 1.D0 / 3.D0 )
  IF( DABS( CubicRootCheck - DNINT( CubicRootCheck ) ) <= 1.D-10 ) THEN
    EXIT ParticleValidation
  ELSE
    WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 70 )//CH_UR
    FrameLeft  = (35 - FLOOR( REAL( 23 ) / 2.D0 ) )
    FrameRight = (35 - CEILING( REAL( 23 ) / 2.D0 ) )
    WRITE( *, "(2G0)" ) CH_VS//REPEAT( " ", FrameLeft )//"VALIDATION OF STRUCTURE"//REPEAT( " ", FrameRight )//CH_VS
    WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 70 )//CH_BR
    WRITE( *, "(G0)") "Invalid number of particles for the selected floppy-box configuration!"
    WRITE( *, "(5G0)") "The total number of particles [", nParticles,"] divided by ", pCell, &
    &                  " (number of particles inside a unit cell) must be a perfect cube root."
    WRITE( *, "(G0)") " "
    WRITE( *, "(G0)") "Exiting..."
    WRITE( *, "(G0)") " "
    STOP
  END IF
END DO ParticleValidation

! Inverse of cell length
CALL InverseMatrixCofactorVec( CellLength, CellLengthInverse, CellVolume )

! Number of unit cells per axis
nCells = NINT( ( DBLE( nParticles ) / DBLE( pCell ) ) ** ( 1.D0 / 3.D0 ) )

! Cell length correction based on the number density of the unit cell
CellLengthCorrection = ( DBLE( pCell ) / ( Density * CellVolume ) ) ** ( 1.D0 / 3.D0 )

! New cell length based on the number density of the unit cell
CellLength(:) = CellLengthCorrection * CellLength(:)

! New positions based on the number density of the unit cell
DO Particle = 1, pCell
  ! Scaled coordinates using old cell length
  CALL MatrixVectorMultiplication( CellLengthInverse, pPositionCell(:,Particle), ScalingDistanceUnitBox )
  ! New real coordinates using new cell length
  CALL MatrixVectorMultiplication( CellLength, ScalingDistanceUnitBox, pPositionCell(:,Particle) )
END DO

! New inverse of cell length and volume
CALL InverseMatrixCofactorVec( CellLength, CellLengthInverse, CellVolume )

! Simulation box length
BoxLength(:) = CellLength(:) * DBLE( nCells )

! Inverse of simulation box length
CALL InverseMatrixCofactorVec( BoxLength, BoxLengthInverse, BoxVolume )

! Position of nonconvex particles (centers of mass)
pCounter = 1
ParticlePlacement: DO
  DO ParticleCell = 1, pCell
    ! Scaled equilibrium position of particle in unit cell
    CALL MatrixVectorMultiplication( CellLengthInverse, pPositionCell(:,ParticleCell), ScalingDistanceUnitBox )
    DO iCell = 1, nCells
      DO jCell = 1, nCells
        DO kCell = 1, nCells
          ! Position of particles in a unit cubic box
          pPosition(1,pCounter) = ScalingDistanceUnitBox(1) + DBLE( iCell - 1 )
          pPosition(2,pCounter) = ScalingDistanceUnitBox(2) + DBLE( jCell - 1 )
          pPosition(3,pCounter) = ScalingDistanceUnitBox(3) + DBLE( kCell - 1 )
          ! Quaternion of particles
          pQuaternion(0,pCounter) = pQuaternionCell(0,ParticleCell)
          pQuaternion(1,pCounter) = pQuaternionCell(1,ParticleCell)
          pQuaternion(2,pCounter) = pQuaternionCell(2,ParticleCell)
          pQuaternion(3,pCounter) = pQuaternionCell(3,ParticleCell)
          ! Increment
          pCounter = pCounter + 1
          ! Condition
          IF( pCounter > nParticles ) THEN
            EXIT ParticlePlacement
          END IF
        END DO
      END DO
    END DO
  END DO
END DO ParticlePlacement

! Correct coordinates inside the triclinic simulation box
DO iParticle = 1, nParticles
  CALL MatrixVectorMultiplication( CellLength, pPosition(:,iParticle), pPosition(:,iParticle) )
END DO

! Centralize unit box at origin of coordination system (0, 0, 0)
pPosition(1,:) = pPosition(1,:) - 0.5D0 * ( BoxLength(1) + BoxLength(4) + BoxLength(7) )
pPosition(2,:) = pPosition(2,:) - 0.5D0 * ( BoxLength(2) + BoxLength(5) + BoxLength(8) )
pPosition(3,:) = pPosition(3,:) - 0.5D0 * ( BoxLength(3) + BoxLength(6) + BoxLength(9) )

! Position of the cylindrical petals (centers of mass)
DO iParticle = 1, nParticles
  ! Apply periodic boundary conditions
  CALL MatrixVectorMultiplication( BoxLengthInverse, pPosition(:,iParticle), ScalingDistanceUnitBox )
  ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
  CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, pPosition(:,iParticle) )
  ! Reposition of cylinders
  CALL VectorRotation( cPositionBasis(:,1), pQuaternion(:,iParticle), cRotatedPositionBasis(:,1) )
  cPosition(:,1,iParticle) = pPosition(:,iParticle) + cRotatedPositionBasis(:,1) ! First quarter
  CALL VectorRotation( cPositionBasis(:,2), pQuaternion(:,iParticle), cRotatedPositionBasis(:,2) )
  cPosition(:,2,iParticle) = pPosition(:,iParticle) + cRotatedPositionBasis(:,2) ! Second quarter
  CALL VectorRotation( cPositionBasis(:,3), pQuaternion(:,iParticle), cRotatedPositionBasis(:,3) )
  cPosition(:,3,iParticle) = pPosition(:,iParticle) + cRotatedPositionBasis(:,3) ! Third quarter
  CALL VectorRotation( cPositionBasis(:,4), pQuaternion(:,iParticle), cRotatedPositionBasis(:,4) )
  cPosition(:,4,iParticle) = pPosition(:,iParticle) + cRotatedPositionBasis(:,4) ! Fourth quarter
END DO

RETURN

END SUBROUTINE FloppyBoxConfiguration

! *********************************************************************************************** !
!          This subroutine creates a file with all particles' positions and orientations          !
! *********************************************************************************************** !
SUBROUTINE ConfigurationOutput(  )

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 ) :: iParticle = 0, cCylinder = 0 ! Counter

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 3 ) :: WritePosition ! Position of particles (center of mass)

! CHARACTER STRINGS
CHARACTER( LEN= 140 ) :: DescriptorString = ' ' ! Descriptor for strings

! Packed-box structure
IF( ConfigurationLogical(1) ) THEN
  OPEN( Unit= 10, File= "Initial_Configuration/OVITO/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )//"_initconf_LD" &
  &                     //TRIM( DescriptorAR )//"_P"//TRIM( DescriptorPressure )//"_pb.xyz" )
  WRITE( 10, "(I5)" ) nParticles * 4
  DescriptorString = "(G0,8(G0.6,1X),G0.6,G0,2(G0.6,1X),G0.6,2G0)"
  WRITE( 10, DescriptorString ) 'Lattice="', BoxLength(1:9), '" Origin="', -0.5D0 * ( BoxLength(1) + BoxLength(4) + &
  &                             BoxLength(7) ), -0.5D0 * ( BoxLength(2) + BoxLength(5) + BoxLength(8) ), -0.5D0 * &
  &                             ( BoxLength(3) + BoxLength(6) + BoxLength(9) ), '" ', &
  &                             "Properties=species:S:1:pos:R:3:orientation:R:4:aspherical_shape:R:3"
  DO iParticle = 1, nParticles
    ! Initial configuration for OVITO (reduced units)
    DO cCylinder = 1, 4
      ! Position of cylinders
      WritePosition(1) = cPosition(1,cCylinder,iParticle)
      WritePosition(2) = cPosition(2,cCylinder,iParticle)
      WritePosition(3) = cPosition(3,cCylinder,iParticle)
      WRITE( 10, "(G0,10(' ',G0.6))" ) iParticle, WritePosition(1), WritePosition(2), WritePosition(3), pQuaternion(1,iParticle), &
      &                                pQuaternion(2,iParticle), pQuaternion(3,iParticle), pQuaternion(0,iParticle), &
      &                                0.5D0 * cDiameter, 0.5D0 * cDiameter, cLength
    END DO
  END DO
  CLOSE( 10 )
! Floppy-box structure
ELSE IF( ConfigurationLogical(2) ) THEN
  OPEN( Unit= 10, File= "Initial_Configuration/OVITO/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )//"_initconf_LD" &
  &                     //TRIM( DescriptorAR )//"_P"//TRIM( DescriptorPressure )//"_fb.xyz" )
  WRITE( 10, "(I5)" ) nParticles * 4
  DescriptorString = "(G0,8(G0.6,1X),G0.6,G0,2(G0.6,1X),G0.6,2G0)"
  WRITE( 10, DescriptorString ) 'Lattice="', BoxLength(1:9), '" Origin="', -0.5D0 * ( BoxLength(1) + BoxLength(4) + &
  &                             BoxLength(7) ), -0.5D0 * ( BoxLength(2) + BoxLength(5) + BoxLength(8) ), -0.5D0 * &
  &                             ( BoxLength(3) + BoxLength(6) + BoxLength(9) ), '" ', &
  &                             "Properties=species:S:1:pos:R:3:orientation:R:4:aspherical_shape:R:3"
  DO iParticle = 1, nParticles
    ! Initial configuration for OVITO (reduced units)
    DO cCylinder = 1, 4
      ! Position of cylinders
      WritePosition(1) = cPosition(1,cCylinder,iParticle)
      WritePosition(2) = cPosition(2,cCylinder,iParticle)
      WritePosition(3) = cPosition(3,cCylinder,iParticle)
      WRITE( 10, "(G0,10(' ',G0.6))" ) iParticle, WritePosition(1), WritePosition(2), WritePosition(3), pQuaternion(1,iParticle), &
      &                                pQuaternion(2,iParticle), pQuaternion(3,iParticle), pQuaternion(0,iParticle), &
      &                                0.5D0 * cDiameter, 0.5D0 * cDiameter, cLength
    END DO
  END DO
  CLOSE( 10 )
END IF

RETURN

END SUBROUTINE ConfigurationOutput

END MODULE InitialConfiguration