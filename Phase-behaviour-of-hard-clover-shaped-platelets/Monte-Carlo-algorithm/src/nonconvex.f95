PROGRAM NonconvexMonteCarlo

! Uses eleven modules: global variables, variable initialization, initial configuration, directory creator, quaternion operations, 
!                      linked lists, overlap check, overlap algorithms, vector operations, and lattice reduction methods
USE GlobalVariables
USE InitializeVariables
USE InitialConfiguration
USE Folders
USE QuaternionOperations
USE LinkedLists, ONLY: FinalizeList, MakeList, BoxCheckNPT, ParticleTranslationNVT
USE OverlapCheckLists
USE OverlapCheckSystem
USE OverlapAlgorithms, ONLY: OverlapCheckCYL, OverlapCheckSPC
USE VectorOperations
USE LatticeReductionMethods, ONLY: LatticeReduction

IMPLICIT NONE

! INTEGER VARIABLES -*- THIS IS SINGLE PRECISION -*-
INTEGER :: SeedSize ! Seed array size

! INTEGER VARIABLES (ALLOCATABLE) -*- THIS IS SINGLE PRECISION -*-
INTEGER, DIMENSION( : ), ALLOCATABLE :: Seed ! Random seed

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 ) :: pParticle, iParticle, jParticle    ! Counters (particles)
INTEGER( Kind= Int64 ) :: cCylinder, iCylinder, jCylinder    ! Counters (cylinders)
INTEGER( Kind= Int64 ) :: bEdge                              ! Counter (box edge)
INTEGER( Kind= Int64 ) :: DensityCounter                     ! Counter (density)
INTEGER( Kind= Int64 ) :: pStack                             ! Counter (particles)
INTEGER( Kind= Int64 ) :: InitialCycle                       ! Initial cycle of the Monte Carlo simulation
INTEGER( Kind= Int64 ) :: wCounter                           ! Counter (Widom insertion)
INTEGER( Kind= Int64 ) :: Particle                           ! Random particle
INTEGER( Kind= Int64 ) :: Cycles                             ! Counter (cycles)
INTEGER( Kind= Int64 ) :: BoxMatrixComponent                 ! Random component of the box matrix
INTEGER( Kind= Int64 ) :: nAcceptanceTranslation             ! Move acceptance counter: Translation
INTEGER( Kind= Int64 ) :: nAcceptanceRotation                ! Move acceptance counter: Rotation
INTEGER( Kind= Int64 ) :: nAcceptanceStack                   ! Move acceptance counter: Stack Rotation
INTEGER( Kind= Int64 ) :: nAcceptanceIsotropicVolumeChange   ! Move acceptance counter: Isotropic volume change
INTEGER( Kind= Int64 ) :: nAcceptanceAnisotropicVolumeChange ! Move acceptance counter: Anistropic volume change
INTEGER( Kind= Int64 ) :: nMovementTranslationCounter        ! Move counter (Translation)
INTEGER( Kind= Int64 ) :: nMovementRotationCounter           ! Move counter (Rotation)
INTEGER( Kind= Int64 ) :: nMovementStackCounter              ! Move counter (Stack rotation)
INTEGER( Kind= Int64 ) :: nMovementIsoVolumeChangeCounter    ! Move counter (Isotropic volume change)
INTEGER( Kind= Int64 ) :: nMovementAnisoVolumeChangeCounter  ! Move counter (Anisotropic volume change)
INTEGER( Kind= Int64 ) :: FrameLeft                          ! Box frame dimension
INTEGER( Kind= Int64 ) :: FrameRight                         ! Box frame dimension

! INTEGER VARIABLES (ARRAY,ALLOCATABLE)
INTEGER( Kind= Int64 ), DIMENSION( : ), ALLOCATABLE :: pStackedParticles, pStackedParticlesTemp ! Particles in a stack (column)

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 ) :: ContactDistance                ! Vega-Lago Contact Distance (variable)
REAL( Kind= Real64 ) :: SquaredDistance                ! Vector distance between particles i and j (squared)
REAL( Kind= Real64 ) :: OldBoxVolume, NewBoxVolume     ! Volume of simulation box (before/after a trial move)
REAL( Kind= Real64 ) :: EnthalpyChange                 ! Enthalpy criterion (reduced)
REAL( Kind= Real64 ) :: VolumeScalingFactor            ! Scaling factor of the volume (isotropic)
REAL( Kind= Real64 ) :: Ratio                          ! Acceptance ratio (simulation)
REAL( Kind= Real64 ) :: wBoltzmannFactor               ! Energy after adding one particle to the system
REAL( Kind= Real64 ) :: wAccumulatedBoltzmannFactor    ! Energy sum after adding one particle to the system
REAL( Kind= Real64 ) :: wAverageBoltzmannFactor        ! Average energy after adding one particle to the system
REAL( Kind= Real64 ) :: BoxDistortionMC                ! Box distortion
REAL( Kind= Real64 ) :: ResidualChemicalPotential      ! Residual chemical potential from the Widom insertion method
REAL( Kind= Real64 ) :: TotalPotential                 ! Total potential energy (orientaitonal field)
REAL( Kind= Real64 ) :: NewPotential, OldPotential     ! Potential energy (before/after a trial move)
REAL( Kind= Real64 ) :: PotentialDifference            ! Potential energy difference (before/after a trial move)
REAL( Kind= Real64 ) :: OrderParameterS2               ! Nematic order parameter
REAL( Kind= Real64 ) :: DensityStandardDeviation       ! Standard deviation of the density
REAL( Kind= Real64 ) :: MaxTranslationalDisplacement   ! Maximum displacement [+/-] (Translation)
REAL( Kind= Real64 ) :: MaxAngularDisplacement         ! Maximum displacement [+/-] (Rotation)
REAL( Kind= Real64 ) :: MaxIsoVolumetricDisplacement   ! Maximum displacement [+/-] (Isotropic volume change)
REAL( Kind= Real64 ) :: MaxAnisoVolumetricDisplacement ! Maximum displacement [+/-] (Anisotropic volume change)
REAL( Kind= Real64 ) :: AlignmentAngle                 ! Alignment angle between the orientation of a particle and the phase drector
REAL( Kind= Real64 ) :: BoxVolumeMC                    ! Volume of simulation box
REAL( Kind= Real64 ) :: StartTimer                     ! Start timer of Monte Carlo simulation
REAL( Kind= Real64 ) :: StopTimer                      ! Stop timer of Monte Carlo simulation
REAL( Kind= Real64 ) :: sRotationAngle                 ! Rotation angle (stack rotation)

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: BoxVectorAngle                           ! Cossine of angle between box vectors
REAL( Kind= Real64 ), DIMENSION( 3 )    :: BoxCutoff                                ! Box cutoff (x-, y-, and z-directions)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: BoxEdgeLength                            ! Length of box edges
REAL( Kind= Real64 ), DIMENSION( 3 )    :: BoxEdgeRatio                             ! Length ratio of box edges
REAL( Kind= Real64 ), DIMENSION( 3 )    :: cRepositionFactor                        ! Rescaling factor for the position of the cylinders
REAL( Kind= Real64 ), DIMENSION( 3 )    :: ScalingDistanceUnitBox                   ! Scaling factor (unit box)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: iOrientation, jOrientation               ! Orientation of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )    :: iPosition, jPosition                     ! Position of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )    :: ciPosition, cjPosition                   ! Position of cylinders of particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )    :: VectorDistance                           ! Vector distance between particles i and j
REAL( Kind= Real64 ), DIMENSION( 3 )    :: WritePosition                            ! Position of the center of mass
REAL( Kind= Real64 ), DIMENSION( 3 )    :: iOldOrientation, iNewOrientation         ! Orientation of a particle (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: iOldPosition, iNewPosition               ! Position of the center of mass of a particle (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( 9 )    :: OldBoxLength, NewBoxLength               ! Length of simulation box (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( 9 )    :: OldBoxLengthInverse, NewBoxLengthInverse ! Inverse of the length of simulation box (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( 9 )    :: BoxLengthMC                              ! Length (x,y,z) of simulation box
REAL( Kind= Real64 ), DIMENSION( 9 )    :: BoxLengthInverseMC                       ! Length (x,y,z) of simulation box (inverse)
REAL( KIND= Real64 ), DIMENSION( 3 )    :: sRandomRotationAxis                      ! Random rotation axis (stack rotation)
REAL( KIND= Real64 ), DIMENSION( 3 )    :: OrthogonalVector1, OrthogonalVector2     ! Unit vectors orthogonal to the phase director of the stack
REAL( KIND= Real64 ), DIMENSION( 3 )    :: PhaseDirector                            ! Phase director (nematic)
REAL( Kind= Real64 ), DIMENSION( 0:3 )  :: iQuaternion, jQuaternion                 ! Quaternion of particles i and j
REAL( Kind= Real64 ), DIMENSION( 0:3 )  :: iOldQuaternion, iNewQuaternion           ! Quaternion of a particle (before/after a trial move)
REAL( KIND= Real64 ), DIMENSION( 0:3 )  :: sRotationQuaternion                      ! Rotation quaternion (stack rotation)
REAL( Kind= Real64 ), DIMENSION( 3, 4 ) :: cOldPosition, cNewPosition               ! Position of cylinders (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( 3, 4 ) :: cRotatedPositionBasis                    ! Rotated body-fixed position (after a trial move)

! REAL VARIABLES (ARRAY,ALLOCATABLE)
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE       :: AverageDensity                           ! Average density
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE    :: pPositionSave                            ! Position array
REAL( KIND= Real64 ), DIMENSION( :, : ), ALLOCATABLE    :: pPositionStack, pPositionStackTemp       ! Position of the center of mass (stack rotation)
REAL( KIND= Real64 ), DIMENSION( :, : ), ALLOCATABLE    :: pQuaternionStack, pQuaternionStackTemp   ! Quaternion of the center of mass (stack rotation)
REAL( KIND= Real64 ), DIMENSION( :, : ), ALLOCATABLE    :: pOrientationStack, pOrientationStackTemp ! Orientation of the center of mass (stack rotation)
REAL( KIND= Real64 ), DIMENSION( :, :, : ), ALLOCATABLE :: cPositionStack, cPositionStackTemp       ! Position of cylinders (stack rotation)
REAL( Kind= Real64 ), DIMENSION( :, :, : ), ALLOCATABLE :: cPositionSave                            ! Position array (cylinders)

! LOGICAL VARIABLES
LOGICAL :: Overlap                          ! Detects overlap between two particles (Vega-Lago) : TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: OverlapSPC                       ! Detects overlap between two particles (preliminary) : TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: ParallelSPC                      ! Checks the relative orientation of two spherocylinders : TRUE = parallel orientation; FALSE = non-parallel orientation
LOGICAL :: OverlapCYL                       ! Detects overlap between two cylinders : TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: LatticeReductionLogical          ! Detects if a lattice reduction is necessary : TRUE = lattice reduction; FALSE = box shape preserved
LOGICAL :: CheckBoxDistortion               ! Detects if a box deformation is valid or not : TRUE = ignore box deformation; FALSE = consider box deformation
LOGICAL :: MovementRotationLogical          ! Rotation move selection : TRUE = movement selected; FALSE = movement not selected
LOGICAL :: MovementTranslationLogical       ! Translation movement selection : TRUE = movement selected; FALSE = movement not selected
LOGICAL :: MovementIsoVolumeChangeLogical   ! Isotropic volume change selection : TRUE = movement selected; FALSE = movement not selected
LOGICAL :: MovementAnisoVolumeChangeLogical ! Anisotropic volume change selection : TRUE = movement selected; FALSE = movement not selected

! CHARACTER STRINGS
CHARACTER( LEN= 03 )  :: EnsembleType         ! Ensemble type
CHARACTER( LEN= 14 )  :: DescriptorBackupFile ! Ensemble type
CHARACTER( LEN= 140 ) :: DescriptorString     ! Descriptor for strings

! Get the maximum number of threads
#ifdef _OPENMP
  nThreads = OMP_GET_MAX_THREADS(  )
#else
  nThreads = 1
#endif

! System properties
OPEN( Unit= 10, File= "ini_system.ini" )
! Skip
READ( 10, * ) Dummy, Dummy
READ( 10, * ) Dummy, Dummy
READ( 10, * ) Dummy, Dummy
READ( 10, * ) Dummy, Dummy
READ( 10, * ) Dummy, Dummy
READ( 10, * ) Dummy, Dummy
! Simulation path ([1]: NPT/NVT simulation ramp, [2]: Normal NPT/NVT simulation)
READ( 10, * ) Dummy, Path
CLOSE( 10 )

! Molecular configuration selection (see 'InitialConfiguration' module)
CALL Configuration_Selection(  )

! Surface geometry selection (see 'InitialConfiguration' module)
CALL Geometry_Selection(  )

! CPU Clock
CALL Date_And_Time( Values= DateTime )

! Initialization of common variables (see 'InitializeVariables' module)
CALL Common_Variables(  )

! Allocation
ALLOCATE( pQuaternion(0:3,nParticles) )
ALLOCATE( pPosition(3,nParticles) )
ALLOCATE( pOrientation(3,nParticles) )
ALLOCATE( cPosition(3,4,nParticles) )
ALLOCATE( pQuaternionMC(0:3,nParticles) )
ALLOCATE( pPositionMC(3,nParticles) )
ALLOCATE( pPositionSave(3,nParticles) )
ALLOCATE( pOrientationMC(3,nParticles) )
ALLOCATE( cPositionMC(3,4,nParticles) )
ALLOCATE( cPositionSave(3,4,nParticles) )

! Initialization of Monte Carlo parameters (see 'InitializeVariables' module)
CALL MonteCarlo_Variables(  )

! Initialization of Inquiry/Control variables (see 'InitializeVariables' module)
CALL Control_Variables(  )

! Allocation
ALLOCATE( AverageDensity( INT( DBLE( nSimulationCycles - nEquilibration ) / DBLE( nSavingFrequency ) ) ) )
AverageDensity = 0.D0

! Status
WRITE( *, "(G0)" ) "Do you wish to continue? [Y/N]"
READ( *, * ) Dummy
CALL ToUpper( Dummy, LEN_TRIM( Dummy ), Dummy )
IF( Dummy /= "Y" ) STOP
WRITE( *, "(G0)" ) " "

! Initialize directories
CALL Initialize_Folders(  )

! Unrotated reference position
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

! Fixed seed
IF( .NOT. RandomSeed ) THEN
  CALL RANDOM_SEED( Size= SeedSize )
  ALLOCATE( Seed(SeedSize) )
  Seed = 123456789
  CALL RANDOM_SEED( Put= Seed )
END IF

! NVT-NPT Monte Carlo simulation (ramp)
IF( Path == 1 .AND. .NOT. RestoreBackup ) THEN
  OPEN( 10, File= "config_old.xyz", Action= "READ" )
  READ( 10, * ) PotentialStrength
  READ( 10, * ) OrientationalFieldVector
  READ( 10, * ) MaxTranslationalDisplacement, MaxAngularDisplacement, MaxIsoVolumetricDisplacement, MaxAnisoVolumetricDisplacement
  READ( 10, * ) BoxLengthMC
  READ( 10, * ) BoxLengthInverseMC
  READ( 10, * ) BoxVolumeMC
  READ( 10, * ) pPositionMC
  READ( 10, * ) pQuaternionMC
  READ( 10, * ) pOrientationMC
  READ( 10, * ) cPositionMC
  CLOSE( 10 )
  ! Initial configuration
  pOrientation = pOrientationMC
  pQuaternion = pQuaternionMC
  pPosition = pPositionMC
  cPosition = cPositionMC
  BoxLength = BoxLengthMC
  BoxLengthInverse = BoxLengthInverseMC
  BoxVolume = BoxVolumeMC
END IF

! Initial configuration (see 'InitialConfiguration' module)
IF( .NOT. RestoreBackup ) THEN
  IF( Path /= 1 ) THEN
    IF( ConfigurationLogical(1) ) THEN
      ! Calls 'PackedBoxConfiguration' subroutine if the user chooses a packed box structure
      CALL PackedBoxConfiguration(  )
    ELSE IF( ConfigurationLogical(2) ) THEN
      ! Calls 'FloppyBoxConfiguration' subroutine if the user chooses a floppy-box structure
      CALL FloppyBoxConfiguration(  )
    END IF
  END IF
  ! Initial cycle
  InitialCycle = 0
  ! Initialization
  DensityCounter = 1
ELSE
  WRITE( *, "(G0)" ) "User has chosen to restore a backup file from a previous simulation. Type the 14-digit backup descriptor: "
  READ( *, * ) DescriptorBackupFile
  INQUIRE( File= "Backup/"//TRIM( DescriptorBackupFile )//"_backup.dat", Exist= FileExist )
  IF( .NOT. FileExist ) THEN
    WRITE( *, "(G0)" ) "File not found."
    STOP
  ELSE
    OPEN( Unit= 220, File= "Backup/"//TRIM( DescriptorBackupFile )//"_backup.dat", Action= "READ" )
    DO pParticle = 1, nParticles
      READ( 220, * ) pPosition(:,pParticle), pQuaternion(:,pParticle)
      DO cCylinder = 1, 4
        READ( 220, * ) cPosition(:,cCylinder,pParticle)
      END DO
    END DO
    READ( 220, * ) BoxLength
    READ( 220, * ) UserMaxTranslationalDisplacement, UserMaxRotationalDisplacement, UserMaxIsoVolumetricDisplacement, &
    &               UserMaxAnisoVolumetricDisplacement
    READ( 220, * ) Density
    READ( 220, * ) PackingFraction
    READ( 220, * ) AverageDensity
    READ( 220, * ) DensityCounter
    IF( WidomInsertion ) THEN
      READ( 220, * ) wAccumulatedBoltzmannFactor
      READ( 220, * ) wCounter
    END IF
    READ( 220, * ) InitialCycle
    READ( 220, * ) Seed
    CALL RANDOM_SEED( Put= Seed )
    CLOSE( 220 )
  END IF
  ! Inverse matrix of the box
  CALL InverseMatrixCofactorVec( BoxLength, BoxLengthInverse, BoxVolume )
END IF

! Initial configuration file (see 'InitialConfiguration' module)
CALL ConfigurationOutput(  )

! Active transformation (rotation)
DO pParticle = 1, nParticles
  CALL VectorRotation( zAxis, pQuaternion(:,pParticle), pOrientation(:,pParticle) )
END DO

! Summary
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 70 )//CH_UR
FrameLeft  = (35 - FLOOR( REAL( 22 ) / 2.D0 ) )
FrameRight = (35 - CEILING( REAL( 22 ) / 2.D0 ) )
WRITE( *, "(2G0)" ) CH_VS//REPEAT( " ", FrameLeft )//"MONTE CARLO SIMULATION"//REPEAT( " ", FrameRight )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 70 )//CH_BR

! Initialize the diameter of the largest sphere
cLargestSphereDiameter = 0.D0
IF( CellListLogical ) cLargestSphereDiameter = CutoffSphere

! Start simulation timer
CALL Cpu_Time( StartTimer )

! Overlap check (initial configuration)
IF( .NOT. RestoreBackup ) THEN
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
      ! Minimum image convention
      CALL MatrixVectorMultiplication( BoxLengthInverse, VectorDistance, ScalingDistanceUnitBox )
      ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
      CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, VectorDistance )
      ! Magnitude of the vector distance (squared)
      SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
      ! Cutoff distance (sphere circumscribing a spherocylinder circumscribing the nonconvex body)
      IF( SquaredDistance <= SquaredCutoffSphere ) THEN
        ! Cutoff distance (spherocylinder circumscribing the non-convex body)
        CALL OverlapCheckSPC( iOrientation, jOrientation, VectorDistance, SquaredDistance, ContactDistance, ParallelSPC, &
        &                     OverlapSPC, .TRUE. )
        ! Vega-Lago criterion (non-convex body)
        IF( OverlapSPC ) THEN
          OverlapCYL = .TRUE.  ! Check overlap between cylinders
        ELSE
          OverlapCYL = .FALSE. ! Do not check overlap between cylinders
        END IF
      END IF
      ! Considering cylinders (only if preliminary tests fail)
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
            ! Minimum image convention
            CALL MatrixVectorMultiplication( BoxLengthInverse, VectorDistance, ScalingDistanceUnitBox )
            ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
            CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, VectorDistance )
            ! Magnitude of the vector distance (squared)
            SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
            ! Cutoff distance (sphere circumscribing a spherocylinder circumscribing the cylinder)
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
END IF

! *********************************************************************************************** !
! Output file units                                                                               !
! *********************************************************************************************** !
IF( .NOT. RestoreBackup ) THEN
  ! Ratio file (translation)
  OPEN( Unit= 30, File= "Ratio/Translation/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )//"_ratio_LD" &
  &                     //TRIM( DescriptorAR )//"_P"//TRIM( DescriptorPressure )//".dat" )
  ! Ratio file (rotation)
  OPEN( Unit= 40, File= "Ratio/Rotation/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )//"_ratio_LD" &
  &                     //TRIM( DescriptorAR )//"_P"//TRIM( DescriptorPressure )//".dat" )
  ! Ratio file (volume)
  OPEN( Unit= 50, File= "Ratio/Volume/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )//"_ratio_LD" &
  &                     //TRIM( DescriptorAR )//"_P"//TRIM( DescriptorPressure )//".dat" )
  ! Ratio file (box length and volume)
  OPEN( Unit= 55, File= "Box/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )//"_box_LD"//TRIM( DescriptorAR )//"_P" &
  &                     //TRIM( DescriptorPressure )//".dat" )
  ! Results file
  OPEN( Unit= 70, File= "Results/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )//"_results_LD" &
  &                     //TRIM( DescriptorAR )//"_P"//TRIM( DescriptorPressure )//".dat" )
  ! Stack rotation file
  IF( StackRotationLogical ) THEN
    OPEN( Unit= 80, File= "Stack_Rotation/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )//"_stack_LD" &
    &                     //TRIM( DescriptorAR )//"_P"//TRIM( DescriptorPressure )//".dat" )
  END IF
  ! Potential file
  IF( PotentialLogical ) THEN
    OPEN( Unit= 90, File= "Potential/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )//"_potential_LD" &
    &                     //TRIM( DescriptorAR )//"_P"//TRIM( DescriptorPressure )//".dat" )
  END IF
ELSE
  ! Ratio file (translation)
  OPEN( Unit= 30, File= "Ratio/Translation/"//TRIM( DescriptorBackupFile(1:8) )//"/"//TRIM( DescriptorBackupFile(9:14) )// &
  &                     "_ratio_LD" //TRIM( DescriptorAR )//"_P"//TRIM( DescriptorPressure )//".dat", Position= "APPEND" )
  ! Ratio file (rotation)
  OPEN( Unit= 40, File= "Ratio/Rotation/"//TRIM( DescriptorBackupFile(1:8) )//"/"//TRIM( DescriptorBackupFile(9:14) )// &
  &                     "_ratio_LD" //TRIM( DescriptorAR )//"_P"//TRIM( DescriptorPressure )//".dat", Position= "APPEND" )
  ! Ratio file (volume)
  OPEN( Unit= 50, File= "Ratio/Volume/"//TRIM( DescriptorBackupFile(1:8) )//"/"//TRIM( DescriptorBackupFile(9:14) )// &
  &                     "_ratio_LD" //TRIM( DescriptorAR )//"_P"//TRIM( DescriptorPressure )//".dat", Position= "APPEND" )
  ! Ratio file (box length and volume)
  OPEN( Unit= 55, File= "Box/"//TRIM( DescriptorBackupFile(1:8) )//"/"//TRIM( DescriptorBackupFile(9:14) )//"_box_LD" &
  &                     //TRIM( DescriptorAR )//"_P"//TRIM( DescriptorPressure )//".dat", Position= "APPEND" )
  ! Results file
  OPEN( Unit= 70, File= "Results/"//TRIM( DescriptorBackupFile(1:8) )//"/"//TRIM( DescriptorBackupFile(9:14) )// &
  &                     "_results_LD" //TRIM( DescriptorAR )//"_P"//TRIM( DescriptorPressure )//".dat", Position= "APPEND" )
  ! Stack rotation file
  IF( StackRotationLogical ) THEN
    OPEN( Unit= 80, File= "Stack_Rotation/"//TRIM( DescriptorBackupFile(1:8) )//"/"//TRIM( DescriptorBackupFile(9:14) )// &
    &                     "_stack_LD" //TRIM( DescriptorAR )//"_P"//TRIM( DescriptorPressure )//".dat", Position= "APPEND" )
  END IF
  ! Potential file
  IF( PotentialLogical ) THEN
    OPEN( Unit= 90, File= "Potential/"//TRIM( DescriptorBackupFile(1:8) )//"/"//TRIM( DescriptorBackupFile(9:14) )// &
    &                     "_potential_LD" //TRIM( DescriptorAR )//"_P"//TRIM( DescriptorPressure )//".dat", Position= "APPEND" )
  END IF
END IF

! *********************************************************************************************** !
! Monte Carlo parameters                                                                          !
! *********************************************************************************************** !
MovementTranslationLogical         = .FALSE.                            ! Translational move selector                      (initial value)
MovementRotationLogical            = .FALSE.                            ! Rotational move selector                         (initial value)
MovementIsoVolumeChangeLogical     = .FALSE.                            ! Isotropic volume change selector                 (initial value)
MovementAnisoVolumeChangeLogical   = .FALSE.                            ! Anisotropic volume change selector               (initial value)
MaxTranslationalDisplacement       = UserMaxTranslationalDisplacement   ! Maximum translational displacement               (initial value)
MaxAngularDisplacement             = UserMaxRotationalDisplacement      ! Maximum rotational displacement                  (initial value)
MaxIsoVolumetricDisplacement       = UserMaxIsoVolumetricDisplacement   ! Maximum volumetric displacement (isotropic)      (initial value)
MaxAnisoVolumetricDisplacement     = UserMaxAnisoVolumetricDisplacement ! Maximum volumetric displacement (anisotropic)    (initial value)
nAcceptanceTranslation             = 0                                  ! Translational move acceptance counter            (initial value)
nAcceptanceRotation                = 0                                  ! Rotational move acceptance counter               (initial value)
nAcceptanceIsotropicVolumeChange   = 0                                  ! Volumetric move acceptance counter (isotropic)   (initial value)
nAcceptanceAnisotropicVolumeChange = 0                                  ! Volumetric move acceptance counter (anisotropic) (initial value)
nMovementTranslationCounter        = 0                                  ! Translational move counter                       (initial value)
nMovementRotationCounter           = 0                                  ! Rotational move counter                          (initial value)
nMovementIsoVolumeChangeCounter    = 0                                  ! Volume change counter (isotropic)                (initial value)
nMovementAnisoVolumeChangeCounter  = 0                                  ! Volume change counter (anisotropic)              (initial value)
IF( Path /= 1 ) THEN
  cPositionMC(:,:,:)                 = cPosition(:,:,:)                   ! Position (cylinders)                             (initial value)
  pQuaternionMC(:,:)                 = pQuaternion(:,:)                   ! Quaternion algebra                               (initial value)
  pPositionMC(:,:)                   = pPosition(:,:)                     ! Position of particles                            (initial value)
  pOrientationMC(:,:)                = pOrientation(:,:)                  ! Orientation of particles                         (initial value)
  BoxLengthMC(:)                     = BoxLength(:)                       ! Box length                                       (initial value)
  BoxLengthInverseMC(:)              = BoxLengthInverse(:)                ! Box length (inverse)                             (initial value)
  BoxVolumeMC                        = BoxVolume                          ! Box volume                                       (initial value)
END IF
IF( StackRotationLogical ) THEN
  nAcceptanceStack                 = 0                                  ! Stack rotation acceptance counter                (initial value)
  nMovementStackCounter            = 0                                  ! Stack rotation counter                           (initial value)
END IF

! Widom insertion parameters
IF( WidomInsertion .AND. .NOT. RestoreBackup ) THEN
  wAccumulatedBoltzmannFactor = 0.D0 ! Energy sum after adding one particle (Boltzmann factor)
  wCounter = 0 ! Counter (Widom insertion)
END IF

! Ensemble type
IF( DABS( MovementProbability - 1.D0 ) < EPSILON( 1.D0 ) ) THEN
  EnsembleType = "NVT"
ELSE
  EnsembleType = "NPT"
END IF

! Trajectory file (depends on user's choice)
IF( PrintTrajectory ) THEN
  IF( .NOT. RestoreBackup ) THEN
    OPEN( Unit= 20, File= "Trajectories/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )//"_traj_LD" &
    &                     //TRIM( DescriptorAR )//"_P"//TRIM( DescriptorPressure )//".xyz" )
    WRITE( 20, "(I5)" ) nParticles * 4
    DescriptorString = "(G0,8(G0.9,1X),G0.9,G0,2(G0.9,1X),G0.9,2G0)"
    WRITE( 20, DescriptorString ) 'Lattice="', BoxLengthMC(1:9), '" Origin="', -0.5D0 * ( BoxLengthMC(1) + BoxLengthMC(4) + &
    &                             BoxLengthMC(7) ), -0.5D0 * ( BoxLengthMC(2) + BoxLengthMC(5) + BoxLengthMC(8) ), -0.5D0 * &
    &                             ( BoxLengthMC(3) + BoxLengthMC(6) + BoxLengthMC(9) ), '" ', &
    &                             "Properties=species:S:1:pos:R:3:orientation:R:4:aspherical_shape:R:3"
    DO pParticle = 1, nParticles
      ! Initial configuration for OVITO (reduced units)
      DO cCylinder = 1, 4
        ! Position of cylinders
        WritePosition(1) = cPositionMC(1,cCylinder,pParticle)
        WritePosition(2) = cPositionMC(2,cCylinder,pParticle)
        WritePosition(3) = cPositionMC(3,cCylinder,pParticle)
        WRITE( 20, "(G0,10(' ',G0.9))") pParticle, WritePosition(1), WritePosition(2), WritePosition(3), &
        &                               pQuaternionMC(1,pParticle), pQuaternionMC(2,pParticle), pQuaternionMC(3,pParticle), &
        &                               pQuaternionMC(0,pParticle), 0.5D0 * cDiameter, 0.5D0 * cDiameter, cLength
      END DO
    END DO
    FLUSH( 20 )
    OPEN( Unit= 25, File= "Trajectories/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )//"_evaluate_bond_parameter_LD" &
    &                     //TRIM( DescriptorAR )//"_P"//TRIM( DescriptorPressure )//".xyz" )
    WRITE( 25, "(2G0)" ) "FRAME #01"
    WRITE( 25, * ) " "
    WRITE( 25, "(G0)", Advance="NO" ) "Box Length: "
    WRITE( 25, "(9(G0.9,1X))" ) BoxLengthMC
    WRITE( 25, * ) " "
    DO pParticle = 1, nParticles
      ! Position of cylinders
      WritePosition(1) = pPositionMC(1,pParticle)
      WritePosition(2) = pPositionMC(2,pParticle)
      WritePosition(3) = pPositionMC(3,pParticle)
      WRITE( 25, "(G0,7(' ',G0.9))" ) pParticle, WritePosition(1), WritePosition(2), WritePosition(3), pQuaternionMC(0,pParticle), &
      &                               pQuaternionMC(1,pParticle), pQuaternionMC(2,pParticle), pQuaternionMC(3,pParticle)
    END DO
    WRITE( 25, * ) " "
    FLUSH( 25 )
  ELSE
    OPEN( Unit= 20, File= "Trajectories/"//TRIM( DescriptorBackupFile(1:8) )//"/"//TRIM( DescriptorBackupFile(9:14) )//"_traj_LD" &
    &                     //TRIM( DescriptorAR )//"_P"//TRIM( DescriptorPressure )//".xyz", Position= "APPEND" )
    OPEN( Unit= 25, File= "Trajectories/"//TRIM( DescriptorBackupFile(1:8) )//"/"//TRIM( DescriptorBackupFile(9:14) )// &
    &                     "_evaluate_bond_parameter_LD" //TRIM( DescriptorAR )//"_P"//TRIM( DescriptorPressure )//".xyz", &
    &                     Position= "APPEND" )
  END IF
END IF

! Initialize cell list
IF( CellListLogical ) THEN
  BoxCutoff(1) = cLargestSphereDiameter / BoxLengthMC(1)
  BoxCutoff(2) = cLargestSphereDiameter / BoxLengthMC(5)
  BoxCutoff(3) = cLargestSphereDiameter / BoxLengthMC(9)
  CALL MakeList( BoxCutoff, pPositionMC, BoxLengthInverseMC )
  OldBoxLength = BoxLengthMC
  OldBoxLengthInverse = BoxLengthInverseMC
END IF

! Total potential calculation
IF( PotentialLogical ) THEN
  TotalPotential = 0.D0
  DO pParticle = 1, nParticles
    AlignmentAngle = DOT_PRODUCT( pOrientationMC(:,pParticle), OrientationalFieldVector(:) ) / ( DSQRT( &
    &                DOT_PRODUCT( pOrientationMC(:,pParticle), pOrientationMC(:,pParticle) ) ) * DSQRT( &
    &                DOT_PRODUCT( OrientationalFieldVector, OrientationalFieldVector ) ) )
    IF( AlignmentAngle >= 1.D0 ) AlignmentAngle = 1.D0
    IF( AlignmentAngle <= - 1.D0 ) AlignmentAngle = - 1.D0
    AlignmentAngle = DACOS( AlignmentAngle )
    TotalPotential = TotalPotential + DSIN( AlignmentAngle ) * DSIN( AlignmentAngle )
  END DO
  TotalPotential = TotalPotential * PotentialStrength
END IF

! Progress bar
CALL Progress_Bar( 0_INT64, nSimulationCycles, EnsembleType )

! Simulation cycles
DO Cycles = InitialCycle + 1, nSimulationCycles

  ! Progress bar
  IF( MOD( Cycles, nSavingFrequency ) == 0 ) THEN
    CALL Progress_Bar( Cycles, nSimulationCycles, EnsembleType )
  END IF

  ! Prepare simulation box for single-particle moves
  IF( CellListLogical .AND. ( MovementIsoVolumeChangeLogical .OR. MovementAnisoVolumeChangeLogical ) ) THEN ! Check cells only after a volume change
    CALL BoxCheckNPT( pPositionMC, OldBoxLength, BoxLengthMC, BoxLengthInverseMC )
  END IF

  ! Random sampling (molecular movement or volume change)
  CALL Random_Number( RandomNumber )

  ! Molecular movement
  IF( RandomNumber < MovementProbability ) THEN

    ! Disable volume movement
    MovementIsoVolumeChangeLogical = .FALSE.
    MovementAnisoVolumeChangeLogical = .FALSE.

    ! Random particle
    CALL Random_Number( RandomNumber )
    Particle = INT( RandomNumber * DBLE( nParticles ) ) + 1

    ! Assignment of previous configuration
    iOldPosition(:)    = pPositionMC(:,Particle)    ! Position
    cOldPosition(:,:)  = cPositionMC(:,:,Particle)  ! Position (cylinders)
    iOldQuaternion(:)  = pQuaternionMC(:,Particle)  ! Quaternion
    iOldOrientation(:) = pOrientationMC(:,Particle) ! Orientation

    ! Generates a random number
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
      CALL MatrixVectorMultiplication( BoxLengthInverseMC, iNewPosition, ScalingDistanceUnitBox )
      ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
      CALL MatrixVectorMultiplication( BoxLengthMC, ScalingDistanceUnitBox, iNewPosition )
    ! Disabled translation
    ELSE IF( .NOT. MovementTranslationLogical ) THEN
      iNewPosition(:) = iOldPosition(:)
    END IF

    ! Rotational movement
    IF( MovementRotationLogical ) THEN
      ! Random Composed Unit Quaternion
      CALL QuaternionCombination( iOldQuaternion, iNewQuaternion, MaxAngularDisplacement )
      ! Active transformation (rotation)
      CALL VectorRotation( zAxis, iNewQuaternion, iNewOrientation )
    ! Disabled rotation
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

    ! Overlap check after displacement of a particle
    IF( .NOT. CellListControl ) THEN
      ! Whole system
      CALL ParticleOverlapCheck( Particle, iNewQuaternion, iNewOrientation, iNewPosition, cNewPosition, BoxLengthMC, &
      &                          BoxLengthInverseMC, Overlap )
    ELSE
      ! Linked lists
      CALL ListOverlapCheck( Particle, iNewQuaternion, iNewOrientation, iNewPosition, cNewPosition, ContactDistance, BoxLengthMC, &
      &                      BoxLengthInverseMC, Overlap, .FALSE. )
    END IF

    ! Calculate potential energy difference
    IF( PotentialLogical ) THEN
      ! Potential of the old configuration
      AlignmentAngle = DOT_PRODUCT( iOldOrientation(:), OrientationalFieldVector(:) ) / ( DSQRT( DOT_PRODUCT( iOldOrientation, &
      &                iOldOrientation ) ) * DSQRT( DOT_PRODUCT( OrientationalFieldVector, OrientationalFieldVector ) ) )
      IF( AlignmentAngle >= 1.D0 ) AlignmentAngle = 1.D0
      IF( AlignmentAngle <= - 1.D0 ) AlignmentAngle = - 1.D0
      OldPotential = PotentialStrength * DSIN( DACOS( AlignmentAngle ) ) * DSIN( DACOS( AlignmentAngle ) )
      ! Potential of the new configuration
      AlignmentAngle = DOT_PRODUCT( iNewOrientation(:), OrientationalFieldVector(:) ) / ( DSQRT( DOT_PRODUCT( iNewOrientation, &
      &                iNewOrientation ) ) * DSQRT( DOT_PRODUCT( OrientationalFieldVector, OrientationalFieldVector ) ) )
      IF( AlignmentAngle >= 1.D0 ) AlignmentAngle = 1.D0
      IF( AlignmentAngle <= - 1.D0 ) AlignmentAngle = - 1.D0
      AlignmentAngle = DACOS( AlignmentAngle )
      NewPotential = PotentialStrength * DSIN( AlignmentAngle ) * DSIN( AlignmentAngle )
      ! Potential difference between the new and old configurations
      PotentialDifference = NewPotential - OldPotential
    END IF

    ! Random number
    IF( PotentialLogical ) CALL Random_Number( RandomNumber )

    ! Metropolis criterion
    IF( .NOT. PotentialLogical .OR. DEXP( - PotentialDifference ) >= RandomNumber ) THEN
      ! Acceptance condition
      IF( .NOT. Overlap ) THEN
        ! Update system configuration
        pPositionMC(:,Particle)    = iNewPosition(:)    ! Update position
        cPositionMC(:,:,Particle)  = cNewPosition(:,:)  ! Update position (cylinders)
        pQuaternionMC(:,Particle)  = iNewQuaternion(:)  ! Update quaternion
        pOrientationMC(:,Particle) = iNewOrientation(:) ! Update orientation
        ! Update displacement counter
        IF( MovementTranslationLogical ) THEN
          IF( CellListControl ) CALL ParticleTranslationNVT( Particle, ScalingDistanceUnitBox ) ! Update cell
          nAcceptanceTranslation = nAcceptanceTranslation + 1 ! Translational move counter
        ELSE IF( MovementRotationLogical ) THEN
          nAcceptanceRotation = nAcceptanceRotation + 1 ! Rotational move counter
        END IF
        ! Update potential energy
        TotalPotential = TotalPotential + PotentialDifference
      ELSE
        ! Retrieve old configuration
        pPositionMC(:,Particle)    = iOldPosition(:)    ! Position
        cPositionMC(:,:,Particle)  = cOldPosition(:,:)  ! Position (cylinders)
        pQuaternionMC(:,Particle)  = iOldQuaternion(:)  ! Quaternion
        pOrientationMC(:,Particle) = iOldOrientation(:) ! Orientation
      END IF
    ELSE
      ! Retrieve old configuration
      pPositionMC(:,Particle)    = iOldPosition(:)    ! Position
      cPositionMC(:,:,Particle)  = cOldPosition(:,:)  ! Position (cylinders)
      pQuaternionMC(:,Particle)  = iOldQuaternion(:)  ! Quaternion
      pOrientationMC(:,Particle) = iOldOrientation(:) ! Orientation
    END IF

  ! Volume change
  ELSE IF( RandomNumber >= MovementProbability ) THEN

    ! Disable translation and rotation
    MovementTranslationLogical = .FALSE.
    MovementRotationLogical   = .FALSE.

    ! Assignment of previous configuration (box)
    OldBoxLength        = BoxLengthMC        ! Box length
    OldBoxLengthInverse = BoxLengthInverseMC ! Box length (inverse)
    OldBoxVolume        = BoxVolumeMC        ! Box volume

    ! Expansion/compression type
    CALL Random_Number( RandomNumber )

    ! Isotropic volume change
    IF( RandomNumber < IsoVolumetricProbability ) THEN
      ! Random scaling factor
      CALL Random_Number( RandomNumber )
      VolumeScalingFactor = DLOG( OldBoxVolume ) + ( 2.D0 * RandomNumber - 1.0D0 ) * MaxIsoVolumetricDisplacement
      VolumeScalingFactor = DEXP( VolumeScalingFactor )
      VolumeScalingFactor = ( VolumeScalingFactor / OldBoxVolume ) ** ( 1.D0 / 3.D0 )
      ! Proportional box length
      NewBoxLength = OldBoxLength * VolumeScalingFactor
      CALL InverseMatrixCofactorVec( NewBoxLength, NewBoxLengthInverse, NewBoxVolume )
      ! Movement counter
      nMovementIsoVolumeChangeCounter = nMovementIsoVolumeChangeCounter + 1
      ! Movement type
      MovementIsoVolumeChangeLogical = .TRUE. ! Enable isotropic volume change
      MovementAnisoVolumeChangeLogical = .FALSE. ! Disable anisotropic volume change
    ! Anisotropic volume change
    ELSE IF( RandomNumber >= IsoVolumetricProbability ) THEN
      ! Random box component
      CALL Random_Number( RandomNumber )
      BoxMatrixComponent = INT( RandomNumber * 6.D0 ) + 1
      IF( BoxMatrixComponent == 1 ) THEN
        BoxMatrixComponent = 1 ! xx component
      ELSE IF( BoxMatrixComponent == 2 ) THEN
        BoxMatrixComponent = 4 ! xy component
      ELSE IF( BoxMatrixComponent == 3 ) THEN
        BoxMatrixComponent = 5 ! yy component
      ELSE IF( BoxMatrixComponent == 4 ) THEN
        BoxMatrixComponent = 7 ! xz component
      ELSE IF( BoxMatrixComponent == 5 ) THEN
        BoxMatrixComponent = 8 ! yz component
      ELSE IF( BoxMatrixComponent == 6 ) THEN
        BoxMatrixComponent = 9 ! zz component
      END IF
      NewBoxLength = OldBoxLength
      ! Random deformation of the box
      CALL Random_Number( RandomNumber )
      NewBoxLength(BoxMatrixComponent) = OldBoxLength(BoxMatrixComponent) + MaxAnisoVolumetricDisplacement * &
      &                                  ( 2.D0 * RandomNumber - 1.0D0 )
      ! Calculate the new reciprocal box basis vectors and the volume of the system
      CALL InverseMatrixCofactorVec( NewBoxLength, NewBoxLengthInverse, NewBoxVolume )
      ! Movement counter
      nMovementAnisoVolumeChangeCounter = nMovementAnisoVolumeChangeCounter + 1
      ! Movement type
      MovementIsoVolumeChangeLogical = .FALSE. ! Disable isotropic volume change
      MovementAnisoVolumeChangeLogical = .TRUE. ! Enable anisotropic volume change
    END IF

    ! Reset condition of anisotropic volume change
    CheckBoxDistortion = .FALSE.

    ! Condition of anisotropic volume change (box distortion)
    IF( MovementAnisoVolumeChangeLogical ) THEN
      ! Box length
      BoxEdgeLength(1) = DSQRT( DOT_PRODUCT( NewBoxLength(1:3), NewBoxLength(1:3) ) )
      BoxEdgeLength(2) = DSQRT( DOT_PRODUCT( NewBoxLength(4:6), NewBoxLength(4:6) ) )
      BoxEdgeLength(3) = DSQRT( DOT_PRODUCT( NewBoxLength(7:9), NewBoxLength(7:9) ) )
      ! Length ratio
      BoxEdgeRatio(1) = BoxEdgeLength(1) / BoxEdgeLength(2)
      BoxEdgeRatio(2) = BoxEdgeLength(1) / BoxEdgeLength(3)
      BoxEdgeRatio(3) = BoxEdgeLength(2) / BoxEdgeLength(3)
      ! Angle between box vectors
      BoxVectorAngle(1) = DOT_PRODUCT( NewBoxLength(1:3), NewBoxLength(4:6) ) / ( BoxEdgeLength(1) * BoxEdgeLength(2) )
      BoxVectorAngle(2) = DOT_PRODUCT( NewBoxLength(1:3), NewBoxLength(7:9) ) / ( BoxEdgeLength(1) * BoxEdgeLength(3) )
      BoxVectorAngle(3) = DOT_PRODUCT( NewBoxLength(4:6), NewBoxLength(7:9) ) / ( BoxEdgeLength(2) * BoxEdgeLength(3) )
      ! Avoid big distortions of the simulation box
      DO bEdge = 1, 3
        ! Angle distortion
        IF( BoxVectorAngle(bEdge) < DCOS( (cPi / 2.D0) + MaxBoxAngle ) .OR. &
        &   BoxVectorAngle(bEdge) > DCOS( (cPi / 2.D0) - MaxBoxAngle ) ) THEN
          BoxVolumeMC        = OldBoxVolume
          BoxLengthMC        = OldBoxLength
          BoxLengthInverseMC = OldBoxLengthInverse
          CheckBoxDistortion = .TRUE.
          EXIT
        END IF
        ! Length distortion
        IF( BoxEdgeRatio(bEdge) > MaxLengthRatio .OR. BoxEdgeRatio(bEdge) < 1.D0 / MaxLengthRatio ) THEN
          BoxVolumeMC        = OldBoxVolume
          BoxLengthMC        = OldBoxLength
          BoxLengthInverseMC = OldBoxLengthInverse
          CheckBoxDistortion = .TRUE.
          EXIT
        END IF
      END DO
    END IF

    ! Box not too distorted
    IF( .NOT. CheckBoxDistortion ) THEN

      ! Enthalpy criterion
      IF( MovementIsoVolumeChangeLogical ) EnthalpyChange = ( Pressure * ( NewBoxVolume - OldBoxVolume ) ) - &
      &                                                     ( DBLE( nParticles + 1 ) * DLOG( NewBoxVolume / OldBoxVolume ) )
      IF( MovementAnisoVolumeChangeLogical ) EnthalpyChange = ( Pressure * ( NewBoxVolume - OldBoxVolume ) ) - &
      &                                                       ( DBLE( nParticles ) * DLOG( NewBoxVolume / OldBoxVolume ) )

      ! Random number
      CALL Random_Number( RandomNumber )

      ! Enthalpy Criterion
      IF( DEXP( - EnthalpyChange ) >= RandomNumber ) THEN

        ! System configuration update
        cPositionSave(:,:,:) = cPositionMC(:,:,:) ! Old configuration
        pPositionSave(:,:) = pPositionMC(:,:) ! Old configuration

        ! Isotropic volume change
        IF( MovementIsoVolumeChangeLogical ) THEN
          ! Rescale positions of particles accordingly
          DO pParticle = 1, nParticles
            ! Rescale lattice position
            pPositionMC(:,pParticle) = pPositionMC(:,pParticle) * VolumeScalingFactor
            ! Rescale factor
            cRepositionFactor(:) = pPositionMC(:,pParticle) - pPositionSave(:,pParticle) ! For cylinders
            ! Cylinders
            DO cCylinder = 1, 4
              cPositionMC(:,cCylinder,pParticle) = cPositionMC(:,cCylinder,pParticle) + cRepositionFactor(:)
            END DO
          END DO
        ! Anisotropic volume change
        ELSE IF( MovementAnisoVolumeChangeLogical ) THEN
          ! Rescale positions of particles accordingly
          DO pParticle = 1, nParticles
            ! Scaling lattice coordinates using the old box length
            CALL MatrixVectorMultiplication( OldBoxLengthInverse, pPositionMC(:,pParticle), ScalingDistanceUnitBox )
            ! New real lattice coordinates using the new box length
            CALL MatrixVectorMultiplication( NewBoxLength, ScalingDistanceUnitBox, pPositionMC(:,pParticle) )
            ! Rescale factor
            cRepositionFactor(:) = pPositionMC(:,pParticle) - pPositionSave(:,pParticle) ! For cylinders
            ! Cylinders
            DO cCylinder = 1, 4
              cPositionMC(:,cCylinder,pParticle) = cPositionMC(:,cCylinder,pParticle) + cRepositionFactor(:)
            END DO
          END DO
        END IF

        ! Overlap check after expansion/compression of the simulation box
        IF( .NOT. CellListControl ) THEN
          ! Whole system
          CALL FullOverlapCheck( ContactDistance, NewBoxLength, NewBoxLengthInverse, Overlap )
        ELSE
          ! Linked lists
          CALL FullListOverlapCheck( ContactDistance, NewBoxLength, NewBoxLengthInverse, Overlap, HalfNeighboursControl )
          ! In case the number of cells in one direction (x, y, or z) becomes less than 3
          IF( .NOT. CellListControl ) CALL FullOverlapCheck( ContactDistance, NewBoxLength, NewBoxLengthInverse, Overlap )
        END IF

        ! Acceptance condition
        IF( .NOT. Overlap ) THEN
          ! Assigns the simulation box properties of a trial volume change to the system configuration.
          BoxVolumeMC           = NewBoxVolume           ! Update volume
          BoxLengthMC(:)        = NewBoxLength(:)        ! Update length
          BoxLengthInverseMC(:) = NewBoxLengthInverse(:) ! Update length (inverse)
          ! Displacement counter update
          IF( MovementIsoVolumeChangeLogical ) THEN
            nAcceptanceIsotropicVolumeChange = nAcceptanceIsotropicVolumeChange + 1 ! Isotropic move counter
          ELSE IF( MovementAnisoVolumeChangeLogical ) THEN
            nAcceptanceAnisotropicVolumeChange = nAcceptanceAnisotropicVolumeChange + 1 ! Anisotropic move counter
          END IF
          ! Update packing fraction and reduced number density
          PackingFraction = TotalParticleVolume / NewBoxVolume
          Density = DBLE( nParticles ) / NewBoxVolume
          ! Re-initialization
          CheckBoxDistortion = .FALSE.
          ! Lattice reduction
          LatticeReductionLogical = .FALSE.
          CALL LatticeReduction( BoxLengthMC, BoxDistortionMC, LatticeReductionLogical )
          IF( LatticeReductionLogical ) THEN
            ! Calculate the new reciprocal box basis vectors
            CALL InverseMatrixCofactorVec( BoxLengthMC, BoxLengthInverseMC, BoxVolumeMC )
            DO pParticle = 1, nParticles
              ! Minimum image convention
              CALL MatrixVectorMultiplication( BoxLengthInverseMC, pPositionMC(:,pParticle), ScalingDistanceUnitBox )
              ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
              CALL MatrixVectorMultiplication( BoxLengthMC, ScalingDistanceUnitBox, pPositionMC(:,pParticle) )
              ! Cylinders
              DO cCylinder = 1, 4
                ! Active transformation (translation)
                CALL VectorRotation( cPositionBasis(:,cCylinder), pQuaternionMC(:,pParticle), cRotatedPositionBasis(:,cCylinder) )
                cPositionMC(:,cCylinder,pParticle) = pPositionMC(:,pParticle) + cRotatedPositionBasis(:,cCylinder)
              END DO
            END DO
            ! Check orientation of the box (eliminate box rotations)
            IF( DABS( BoxLengthMC(2) - 0.D0 ) >= EPSILON( 1.D0 ) .OR. DABS( BoxLengthMC(3) - 0.D0 ) >= EPSILON( 1.D0 ) .OR. &
            &   DABS( BoxLengthMC(6) - 0.D0 ) >= EPSILON( 1.D0 ) ) THEN
              ! Undo box rotation
              CALL UndoBoxRotation( BoxLengthMC, BoxLengthInverseMC, BoxVolumeMC )
            END IF
          END IF
        ! Retrieve old properties of the system configuration
        ELSE
          BoxVolumeMC           = OldBoxVolume           ! Box volume
          BoxLengthMC(:)        = OldBoxLength(:)        ! Box length
          BoxLengthInverseMC(:) = OldBoxLengthInverse(:) ! Box length (inverse)
          pPositionMC(:,:)      = pPositionSave(:,:)     ! Position
          cPositionMC(:,:,:)    = cPositionSave(:,:,:)   ! Position (cylinders)
        END IF

      END IF ! Enthalpy criterion

    END IF ! Box distortion criterion

  END IF

  ! *********************************************************************************************** !
  ! Stack rotation moves                                                                            !
  ! *********************************************************************************************** !
  ! See Marechal M., Patti A., Dennison M., and Dijkstra M. Phys. Rev. Lett. 108, 206101 (2012)     !
  ! for more information.                                                                           !
  ! *********************************************************************************************** !

  ! Condition for stack rotation moves
  IF( Density >= StackRotationDensityThreshold .AND. StackRotationLogical ) THEN
    ! Random number
    CALL Random_Number( RandomNumber )
    IF( RandomNumber < StackRotationProbability ) THEN
      ! Update counter
      nMovementStackCounter = nMovementStackCounter + 1
      ! Initialize number of particles in stack
      pStack = 1
      ! Initial allocation
      IF( ALLOCATED( pPositionStack ) ) DEALLOCATE( pPositionStack )
      IF( ALLOCATED( pQuaternionStack ) ) DEALLOCATE( pQuaternionStack )
      IF( ALLOCATED( pOrientationStack ) ) DEALLOCATE( pOrientationStack )
      IF( ALLOCATED( cPositionStack ) ) DEALLOCATE( cPositionStack )
      IF( ALLOCATED( pStackedParticles ) ) DEALLOCATE( pStackedParticles )
      ALLOCATE( pPositionStack(3,pStack), pQuaternionStack(0:3,pStack), pOrientationStack(3,pStack), cPositionStack(3,4,pStack) )
      ALLOCATE( pStackedParticles(pStack) )
      ! Random particle
      CALL Random_Number( RandomNumber )
      Particle = INT( RandomNumber * DBLE( nParticles ) ) + 1
      ! Add particle to the stack list
      pStackedParticles(pStack) = Particle
      ! Position of the particle in the stack
      pPositionStack(:,pStack) = pPositionMC(:,Particle)
      ! Quaternion of the particle in the stack
      pQuaternionStack(:,pStack) = pQuaternionMC(:,Particle)
      ! Orientation of the particle in the stack
      pOrientationStack(:,pStack) = pOrientationMC(:,Particle)
      ! Position of the cylinders in the stack
      cPositionStack(:,:,pStack) = cPositionMC(:,:,Particle)
      ! Find and place other particles in the stack
      DO pParticle = 1, nParticles
        ! Skip the particle already in the stack
        IF( pParticle == Particle ) CYCLE
        ! Vector distance between particles i and j
        VectorDistance(:) = pPositionMC(:,pParticle) - pPositionMC(:,Particle)
        ! Minimum image convention
        CALL MatrixVectorMultiplication( BoxLengthInverseMC, VectorDistance, ScalingDistanceUnitBox )
        ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
        CALL MatrixVectorMultiplication( BoxLengthMC, ScalingDistanceUnitBox, VectorDistance )
        ! Magnitude of the vector distance (squared)
        SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
        ! Condition for stacking
        IF( SquaredDistance < 0.25D0 * CutoffSpherocylinder * CutoffSpherocylinder ) THEN
          ! Increment number of particles in the stack
          pStack = pStack + 1
          ! Initial allocation of temporary arrays
          IF( ALLOCATED( pPositionStackTemp ) ) DEALLOCATE( pPositionStackTemp )
          IF( ALLOCATED( pQuaternionStackTemp ) ) DEALLOCATE( pQuaternionStackTemp )
          IF( ALLOCATED( pOrientationStackTemp ) ) DEALLOCATE( pOrientationStackTemp )
          IF( ALLOCATED( cPositionStackTemp ) ) DEALLOCATE( cPositionStackTemp )
          IF( ALLOCATED( pStackedParticlesTemp ) ) DEALLOCATE( pStackedParticlesTemp )
          ALLOCATE( pPositionStackTemp(3,SIZE(pStackedParticles)), pQuaternionStackTemp(0:3,SIZE(pStackedParticles)) )
          ALLOCATE( pOrientationStackTemp(3,SIZE(pStackedParticles)), cPositionStackTemp(3,4,SIZE(pStackedParticles)) )
          ALLOCATE( pStackedParticlesTemp(SIZE(pStackedParticles)) )
          ! Assign temporary arrays
          pPositionStackTemp(:,:) = pPositionStack(:,:)
          pQuaternionStackTemp(:,:) = pQuaternionStack(:,:)
          pOrientationStackTemp(:,:) = pOrientationStack(:,:)
          cPositionStackTemp(:,:,:) = cPositionStack(:,:,:)
          pStackedParticlesTemp(:) = pStackedParticles(:)
          ! Re-allocation of the stack list
          IF( ALLOCATED( pPositionStack ) ) DEALLOCATE( pPositionStack )
          IF( ALLOCATED( pQuaternionStack ) ) DEALLOCATE( pQuaternionStack )
          IF( ALLOCATED( pOrientationStack ) ) DEALLOCATE( pOrientationStack )
          IF( ALLOCATED( cPositionStack ) ) DEALLOCATE( cPositionStack )
          IF( ALLOCATED( pStackedParticles ) ) DEALLOCATE( pStackedParticles )
          ALLOCATE( pPositionStack(3,SIZE(pStackedParticlesTemp)+1), pQuaternionStack(0:3,SIZE(pStackedParticlesTemp)+1) )
          ALLOCATE( pOrientationStack(3,SIZE(pStackedParticlesTemp)+1), cPositionStack(3,4,SIZE(pStackedParticlesTemp)+1) )
          ALLOCATE( pStackedParticles(SIZE(pStackedParticlesTemp)+1) )
          ! Reassign stack list
          pPositionStack(:,1:SIZE(pStackedParticlesTemp)) = pPositionStackTemp(:,1:SIZE(pStackedParticlesTemp))
          pQuaternionStack(:,1:SIZE(pStackedParticlesTemp)) = pQuaternionStackTemp(:,1:SIZE(pStackedParticlesTemp))
          pOrientationStack(:,1:SIZE(pStackedParticlesTemp)) = pOrientationStackTemp(:,1:SIZE(pStackedParticlesTemp))
          cPositionStack(:,:,1:SIZE(pStackedParticlesTemp)) = cPositionStackTemp(:,:,1:SIZE(pStackedParticlesTemp))
          pStackedParticles(1:SIZE(pStackedParticlesTemp)) = pStackedParticlesTemp(1:SIZE(pStackedParticlesTemp))
          ! Add particle to the stack list
          pStackedParticles(pStack) = pParticle
          ! Position of the particle in the stack
          pPositionStack(:,pStack) = pPositionMC(:,pParticle)
          ! Quaternion of the particle in the stack
          pQuaternionStack(:,pStack) = pQuaternionMC(:,pParticle)
          ! Orientation of the particle in the stack
          pOrientationStack(:,pStack) = pOrientationMC(:,pParticle)
          ! Position of the cylinders in the stack
          cPositionStack(:,:,pStack) = cPositionMC(:,:,pParticle)
        END IF
      END DO
      ! Find vector orthogonal to the director of the stack
      IF( SIZE( pStackedParticles ) > 1 ) THEN
        ! Phase director of the stack
        CALL UniaxialNematicOrderParameter( INT( SIZE( pStackedParticles ), Kind= Int64 ), OrderParameterS2, pOrientationStack, &
        &                                   PhaseDirector, .TRUE. )
        ! Normalize phase director
        PhaseDirector = PhaseDirector / DSQRT( DOT_PRODUCT( PhaseDirector, PhaseDirector ) )
        ! First trial vector (orthogonal to the phase director)
        OrthogonalVector1 = [ -PhaseDirector(2), PhaseDirector(1), 0.D0 ]
        ! Check if the phase director is aligned with the z-axis
        IF( DABS( PhaseDirector(3) ) >= (1.D0 - 1.D-8) ) THEN
          OrthogonalVector1 = [ 1.D0, 0.D0, 0.D0 ]
        END IF
        ! Second trial vector (orthogonal to the phase director and the first trial vector)
        CALL Cross_Product( PhaseDirector, OrthogonalVector1, OrthogonalVector2 )
        ! Random rotation angle
        CALL Random_Number( RandomNumber ) ! Range: [0, 1)
        sRotationAngle = 2.D0 * cPi * RandomNumber ! Range: [0, 2π)
        ! Random point in a unit circle
        sRandomRotationAxis = DCOS( sRotationAngle ) * OrthogonalVector1 + DSIN( sRotationAngle ) * OrthogonalVector2
      ELSE
        ! Single particle
        PhaseDirector = pOrientationMC(:,Particle)
        ! Normalize phase director
        PhaseDirector = PhaseDirector / DSQRT( DOT_PRODUCT( PhaseDirector, PhaseDirector ) )
        ! First trial vector (orthogonal to the phase director)
        OrthogonalVector1 = [ -PhaseDirector(2), PhaseDirector(1), 0.D0 ]
        ! Check if the phase director is aligned with the z-axis
        IF( DABS( PhaseDirector(3) ) >= (1.D0 - 1.D-8) ) THEN
          OrthogonalVector1 = [ 1.D0, 0.D0, 0.D0 ]
        END IF
        ! Second trial vector (orthogonal to the phase director and the first trial vector)
        CALL Cross_Product( PhaseDirector, OrthogonalVector1, OrthogonalVector2 )
        ! Random rotation angle
        CALL Random_Number( RandomNumber ) ! Range: [0, 1)
        sRotationAngle = 2.D0 * cPi * RandomNumber ! Range: [0, 2π)
        ! Random point in a unit circle
        sRandomRotationAxis = DCOS( sRotationAngle ) * OrthogonalVector1 + DSIN( sRotationAngle ) * OrthogonalVector2
      END IF
      ! Normalize vector
      sRandomRotationAxis = sRandomRotationAxis / DSQRT( DOT_PRODUCT( sRandomRotationAxis, sRandomRotationAxis ) )
      ! Rotation angle
      sRotationAngle = 0.5D0 * cPi ! 90 degrees
      ! Rotation quaternion
      sRotationQuaternion(0) = DCOS( sRotationAngle * 0.5D0 ) ! Real part
      sRotationQuaternion(1) = DSIN( sRotationAngle * 0.5D0 ) * sRandomRotationAxis(1) ! Imaginary part (Vector)
      sRotationQuaternion(2) = DSIN( sRotationAngle * 0.5D0 ) * sRandomRotationAxis(2) ! Imaginary part (Vector)
      sRotationQuaternion(3) = DSIN( sRotationAngle * 0.5D0 ) * sRandomRotationAxis(3) ! Imaginary part (Vector)
      ! Update positions and orientations of the particles after stack rotation
      DO pParticle = 1, SIZE( pStackedParticles )
        ! Head of the stack
        IF( pStackedParticles(pParticle) == Particle ) THEN
          ! Update quaternion
          CALL QuaternionMultiplication( sRotationQuaternion, pQuaternionMC(:,pStackedParticles(pParticle)), &
          &                              pQuaternionMC(:,pStackedParticles(pParticle)) )
          ! Update orientation
          CALL VectorRotation( zAxis, pQuaternionMC(:,pStackedParticles(pParticle)), &
          &                    pOrientationMC(:,pStackedParticles(pParticle)) )
          ! Update position of the cylinders
          DO cCylinder = 1, 4
            ! Active transformation (translation)
            CALL VectorRotation( cPositionBasis(:,cCylinder), pQuaternionMC(:,pStackedParticles(pParticle)), &
            &                    cRotatedPositionBasis(:,cCylinder) )
            cPositionMC(:,cCylinder,pStackedParticles(pParticle)) = pPositionMC(:,pStackedParticles(pParticle)) + &
            &                                                       cRotatedPositionBasis(:,cCylinder)
          END DO
        ! Apply periodic conditions to a particle in the stack because the rotation of the whole stack is performed with respect to the center of mass of the first particle in the stack
        ELSE
          ! Vector distance between the particle in the stack and the particle around which the stack is rotated
          VectorDistance(:) = pPositionMC(:,pStackedParticles(pParticle)) - pPositionMC(:,Particle)
          ! Minimum image convention
          CALL MatrixVectorMultiplication( BoxLengthInverseMC, VectorDistance, ScalingDistanceUnitBox )
          ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
          CALL MatrixVectorMultiplication( BoxLengthMC, ScalingDistanceUnitBox, VectorDistance )
          ! Rotate the vector distance
          CALL VectorRotation( VectorDistance, sRotationQuaternion, VectorDistance )
          ! Update position by adding the rotated vector distance to the position of the particle around which the stack is rotated
          pPositionMC(:,pStackedParticles(pParticle)) = pPositionMC(:,Particle) + VectorDistance
          ! Apply periodic conditions to the particle in the stack
          CALL MatrixVectorMultiplication( BoxLengthInverseMC, pPositionMC(:,pStackedParticles(pParticle)), ScalingDistanceUnitBox )
          ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
          CALL MatrixVectorMultiplication( BoxLengthMC, ScalingDistanceUnitBox, pPositionMC(:,pStackedParticles(pParticle)) )
          ! Update quaternion
          CALL QuaternionMultiplication( sRotationQuaternion, pQuaternionMC(:,pStackedParticles(pParticle)), &
          &                              pQuaternionMC(:,pStackedParticles(pParticle)) )
          ! Update orientation
          CALL VectorRotation( zAxis, pQuaternionMC(:,pStackedParticles(pParticle)), &
          &                    pOrientationMC(:,pStackedParticles(pParticle)) )
          ! Update position of the cylinders
          DO cCylinder = 1, 4
            ! Active transformation (translation)
            CALL VectorRotation( cPositionBasis(:,cCylinder), pQuaternionMC(:,pStackedParticles(pParticle)), &
            &                    cRotatedPositionBasis(:,cCylinder) )
            cPositionMC(:,cCylinder,pStackedParticles(pParticle)) = pPositionMC(:,pStackedParticles(pParticle)) + &
            &                                                       cRotatedPositionBasis(:,cCylinder)
          END DO
        END IF
      END DO
      ! Check if the stack is overlapping other particles in the system after rotation
      DO pParticle = 1, SIZE( pStackedParticles )
        CALL ParticleOverlapCheck( pStackedParticles(pParticle), pQuaternionMC(:,pStackedParticles(pParticle)), &
        &                          pOrientationMC(:,pStackedParticles(pParticle)), pPositionMC(:,pStackedParticles(pParticle)), &
        &                          cPositionMC(:,:,pStackedParticles(pParticle)), BoxLengthMC, BoxLengthInverseMC, Overlap )
        ! Exit loop if overlap is detected
        IF( Overlap ) EXIT
      END DO
      IF( .NOT. Overlap ) THEN
        ! Update acceptance counter for stack rotation moves
        nAcceptanceStack = nAcceptanceStack + 1
        ! Update cell list
        IF( CellListControl ) THEN
          DO pParticle = 1, SIZE( pStackedParticles )
            CALL MatrixVectorMultiplication( BoxLengthInverseMC, pPositionMC(:,pStackedParticles(pParticle)), &
            &                                ScalingDistanceUnitBox )
            ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
            CALL ParticleTranslationNVT( pStackedParticles(pParticle), ScalingDistanceUnitBox ) ! Update cell
          END DO
        END IF
      ELSE
        ! Retrieve old configuration
        DO pParticle = 1, SIZE( pStackedParticles )
          pPositionMC(:,pStackedParticles(pParticle)) = pPositionStack(:,pParticle)       ! Position
          pQuaternionMC(:,pStackedParticles(pParticle)) = pQuaternionStack(:,pParticle)   ! Quaternion
          pOrientationMC(:,pStackedParticles(pParticle)) = pOrientationStack(:,pParticle) ! Orientation
          cPositionMC(:,:,pStackedParticles(pParticle)) = cPositionStack(:,:,pParticle)   ! Position (cylinders)
        END DO
      END IF
    END IF
  END IF

  ! Adjustment of maximum displacements
  IF( Cycles <= nEquilibration ) THEN ! During equilibration only

    ! Translational adjustment
    IF( MOD( Cycles, nAdjustmentMovementFrequency ) == 0 ) THEN
      IF( nMovementTranslationCounter > 0 ) THEN
        ! Acceptance ratio (non-overlapping microstates over sampled microstates)
        Ratio = DBLE( nAcceptanceTranslation ) / DBLE( nMovementTranslationCounter )
        ! Translational adjustment
        IF( Ratio <= AcceptanceRatioTranslation ) THEN
          MaxTranslationalDisplacement  = 0.95D0 * MaxTranslationalDisplacement
        ELSE
          MaxTranslationalDisplacement  = 1.05D0 * MaxTranslationalDisplacement
        END IF
        ! Ratio data
        WRITE( 30, "(G0,3(' ',G0.6))" ) Cycles, Ratio, MaxTranslationalDisplacement, AcceptanceRatioTranslation
        FLUSH( 30 )
        ! Reset counter
        nAcceptanceTranslation = 0
        nMovementTranslationCounter  = 0
        ! Avoid multiple turns
        IF( MaxTranslationalDisplacement >= 4.D0 * MAXVAL( BoxLengthMC ) ) THEN
          MaxTranslationalDisplacement = MaxTranslationalDisplacement - 2.D0 * MAXVAL( BoxLengthMC )
        END IF
      END IF
    END IF

    ! Rotational adjustment
    IF( MOD( Cycles, nAdjustmentMovementFrequency ) == 0 ) THEN
      IF( nMovementRotationCounter > 0 ) THEN
        ! Acceptance ratio (non-overlapping microstates over sampled microstates)
        Ratio = DBLE( nAcceptanceRotation ) / DBLE( nMovementRotationCounter )
        ! Rotation adjustment
        IF( Ratio <= AcceptanceRatioRotation ) THEN
          MaxAngularDisplacement = 0.95D0 * MaxAngularDisplacement
        ELSE
          MaxAngularDisplacement = 1.05D0 * MaxAngularDisplacement
        END IF
        ! Ratio data
        WRITE( 40, "(G0,3(' ',G0.6))" ) Cycles, Ratio, MaxAngularDisplacement, AcceptanceRatioRotation
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

    ! Volumetric adjustment (isotropic)
    IF( MOD( Cycles, nAdjustmentVolumeFrequency ) == 0 ) THEN
      IF( nMovementIsoVolumeChangeCounter > 0 ) THEN
        ! Acceptance ratio (non-overlapping microstates over sampled microstates)
        Ratio = DBLE( nAcceptanceIsotropicVolumeChange ) / DBLE( nMovementIsoVolumeChangeCounter )
        ! Volumetric adjustment
        IF( Ratio <= AcceptanceRatioIsoVolumeChange ) THEN
          MaxIsoVolumetricDisplacement = 0.95D0 * MaxIsoVolumetricDisplacement
        ELSE
          MaxIsoVolumetricDisplacement = 1.05D0 * MaxIsoVolumetricDisplacement
        END IF
        ! Ratio data
        WRITE( 50, "(G0,3(' ',G0.6),G0)" ) Cycles, Ratio, MaxIsoVolumetricDisplacement, AcceptanceRatioIsoVolumeChange, "I"
        FLUSH( 50 )
        ! Reset counter
        nAcceptanceIsotropicVolumeChange = 0
        nMovementIsoVolumeChangeCounter  = 0
      END IF
    END IF

    ! Volumetric adjustment (anisotropic)
    IF( MOD( Cycles, nAdjustmentVolumeFrequency ) == 0 ) THEN
      IF( nMovementAnisoVolumeChangeCounter > 0 ) THEN
        ! Acceptance ratio (non-overlapping microstates over sampled microstates)
        Ratio = DBLE( nAcceptanceAnisotropicVolumeChange ) / DBLE( nMovementAnisoVolumeChangeCounter )
        ! Volumetric adjustment
        IF( Ratio <= AcceptanceRatioAnisoVolumeChange ) THEN
          MaxAnisoVolumetricDisplacement = 0.95D0 * MaxAnisoVolumetricDisplacement
        ELSE
          MaxAnisoVolumetricDisplacement = 1.05D0 * MaxAnisoVolumetricDisplacement
        END IF
        ! Ratio data
        WRITE( 50, "(G0,3(' ',G0.6),G0)" ) Cycles, Ratio, MaxAnisoVolumetricDisplacement, AcceptanceRatioAnisoVolumeChange, "A"
        FLUSH( 50 )
        ! Reset counter
        nAcceptanceAnisotropicVolumeChange = 0
        nMovementAnisoVolumeChangeCounter  = 0
      END IF
    END IF

  END IF

  ! Ratio data (box properties)
  IF( MOD( Cycles, nSavingFrequency ) == 0 .AND. EnsembleType == "NPT" ) THEN
    WRITE( 55, "(2G0,G0.3)" ) Cycles, ",", BoxDistortionMC
    FLUSH( 55 )
  END IF

  ! Trajectory data
  IF( PrintTrajectory ) THEN
    IF( MOD( Cycles, nSavingFrequency ) == 0 ) THEN
      WRITE( 20, "(I5)" ) nParticles * 4
      DescriptorString = "(G0,8(G0.9,1X),G0.9,G0,2(G0.9,1X),G0.9,2G0)"
      WRITE( 20, DescriptorString ) 'Lattice="', BoxLengthMC(1:9), '" Origin="', -0.5D0 * ( BoxLengthMC(1) + BoxLengthMC(4) + &
      &                             BoxLengthMC(7) ), -0.5D0 * ( BoxLengthMC(2) + BoxLengthMC(5) + BoxLengthMC(8) ), -0.5D0 * &
      &                             ( BoxLengthMC(3) + BoxLengthMC(6) + BoxLengthMC(9) ), '" ', &
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
      WRITE( 25, "(2G0)" ) "Frame #", Cycles
      WRITE( 25, * ) " "
      WRITE( 25, "(G0)", Advance="NO" ) "Box Length: "
      WRITE( 25, "(9(G0.9,1X))" ) BoxLengthMC
      WRITE( 25, * ) " "
      DO pParticle = 1, nParticles
        WritePosition(1) = pPositionMC(1,pParticle)
        WritePosition(2) = pPositionMC(2,pParticle)
        WritePosition(3) = pPositionMC(3,pParticle)
        WRITE( 25, "(G0,7(' ',G0.9))" ) pParticle, WritePosition(1), WritePosition(2), WritePosition(3), &
        &                               pQuaternionMC(0,pParticle), pQuaternionMC(1,pParticle), pQuaternionMC(2,pParticle), &
        &                               pQuaternionMC(3,pParticle)
      END DO
      WRITE( 25, * ) " "
      FLUSH( 25 )
    END IF
  END IF

  ! Results data
  IF( MOD( Cycles, nSavingFrequency ) == 0 .AND. EnsembleType == "NPT" ) THEN
    ! Packing fraction, density and volume
    CALL UniaxialNematicOrderParameter( nParticles, OrderParameterS2, pOrientationMC, PhaseDirector, .FALSE. )
    WRITE( 70, "(G0,4(' ',G0.6))" ) Cycles, PackingFraction, Density, Pressure, OrderParameterS2
    FLUSH( 70 )
  END IF

  ! Stack rotation moves
  IF( Density >= StackRotationDensityThreshold .AND. StackRotationLogical ) THEN
    ! Update stack rotation counter
    IF( MOD( Cycles, nSavingFrequency ) == 0 ) THEN
      WRITE( 80, "(G0,' ',G0.5)" ) Cycles, DBLE( nAcceptanceStack ) / DBLE( nMovementStackCounter )
      FLUSH( 80 )
    END IF
  END IF

  ! Potential energy
  IF( PotentialLogical ) THEN
    IF( MOD( Cycles, nSavingFrequency ) == 0 ) THEN
      WRITE( 90, "(G0,' ',G0.5)" ) Cycles, TotalPotential
      FLUSH( 90 )
    END IF
  END IF

  ! Widom insertion method to calculate the residual chemical potential
  IF( ( Cycles > nEquilibration ) .AND. WidomInsertion .AND. MOD( Cycles, nSavingFrequency ) == 0 ) THEN
    CALL WidomInsertionMethod( BoxLengthMC, BoxLengthInverseMC, wBoltzmannFactor )
    wAccumulatedBoltzmannFactor = wAccumulatedBoltzmannFactor + wBoltzmannFactor
    wCounter = wCounter + 1
  END IF

  ! Average density to calculate the DensityStandardDeviation and standard deviation
  IF( ( Cycles > nEquilibration ) .AND. MOD( Cycles, nSavingFrequency ) == 0 ) THEN
    AverageDensity(DensityCounter) = Density
    DensityCounter = DensityCounter + 1
  END IF

  ! Backup file
  IF( ( MOD( Cycles, nSavingFrequency ) == 0 ) ) THEN
    IF( .NOT. RestoreBackup ) THEN
      OPEN( Unit= 220, File= "Backup/"//TRIM( DescriptorDate )//""//TRIM( DescriptorHour )//"_backup.dat" )
    ELSE
      OPEN( Unit= 220, File= "Backup/"//TRIM( DescriptorBackupFile(1:8) )//""//TRIM( DescriptorBackupFile(9:14) )//"_backup.dat" )
    END IF
    DO pParticle = 1, nParticles
      WRITE( 220, * ) pPositionMC(:,pParticle), pQuaternionMC(:,pParticle)
      DO cCylinder = 1, 4
        WRITE( 220, * ) cPositionMC(:,cCylinder,pParticle)
      END DO
    END DO
    WRITE( 220, * ) BoxLengthMC
    WRITE( 220, * ) MaxTranslationalDisplacement, MaxAngularDisplacement, MaxIsoVolumetricDisplacement, &
    &               MaxAnisoVolumetricDisplacement
    WRITE( 220, * ) Density
    WRITE( 220, * ) PackingFraction
    WRITE( 220, * ) AverageDensity
    WRITE( 220, * ) DensityCounter
    IF( WidomInsertion ) THEN
      WRITE( 220, * ) wAccumulatedBoltzmannFactor
      WRITE( 220, * ) wCounter
    END IF
    WRITE( 220, * ) Cycles
    CALL RANDOM_SEED( Get= Seed )
    WRITE( 220, * ) Seed
    CLOSE( 220 )
  END IF

END DO

! Final results
CALL UniaxialNematicOrderParameter( nParticles, OrderParameterS2, pOrientationMC, PhaseDirector, .FALSE. )
WRITE( *, "(G0)" ) " "
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 70 )//CH_UR
FrameLeft  = (35 - FLOOR( REAL( 13 ) / 2.D0 ) )
FrameRight = (35 - CEILING( REAL( 13 ) / 2.D0 ) )
WRITE( *, "(2G0)" ) CH_VS//REPEAT( " ", FrameLeft )//"FINAL RESULTS"//REPEAT( " ", FrameRight )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 70 )//CH_BR
WRITE( *, "(G0,G0.5)" ) "Potential strength: ", PotentialStrength
WRITE( *, "(G0,G0.5)" ) "Potential energy: ", TotalPotential
WRITE( *, "(G0,G0.5)" ) "Order parameter: ", OrderParameterS2
WRITE( *, "(G0,G0.5)" ) "Packing fraction: ", PackingFraction
WRITE( *, "(G0,G0.5)" ) "Number density: ", Density
WRITE( *, "(G0)" ) " "

! Save configuration for pressure ramp
IF( Path == 1 ) THEN
  OPEN( Unit= 10, File= "config_new.xyz")
  WRITE( 10, * ) PotentialStrength
  WRITE( 10, * ) OrientationalFieldVector
  WRITE( 10, * ) MaxTranslationalDisplacement, MaxAngularDisplacement, MaxIsoVolumetricDisplacement, MaxAnisoVolumetricDisplacement
  WRITE( 10, * ) BoxLengthMC
  WRITE( 10, * ) BoxLengthInverseMC
  WRITE( 10, * ) BoxVolumeMC
  WRITE( 10, * ) pPositionMC
  WRITE( 10, * ) pQuaternionMC
  WRITE( 10, * ) pOrientationMC
  WRITE( 10, * ) cPositionMC
  CLOSE( 10 )
END IF

! Save last configuration
OPEN( Unit= 10, File= "Configurations/config"//TRIM( DescriptorDate )//TRIM( DescriptorHour )//"_P"//TRIM( DescriptorPressure )// &
&                     ".xyz")
WRITE( 10, * ) PotentialStrength
WRITE( 10, * ) OrientationalFieldVector
WRITE( 10, * ) MaxTranslationalDisplacement, MaxAngularDisplacement, MaxIsoVolumetricDisplacement, MaxAnisoVolumetricDisplacement
WRITE( 10, * ) BoxLengthMC
WRITE( 10, * ) BoxLengthInverseMC
WRITE( 10, * ) BoxVolumeMC
WRITE( 10, * ) pPositionMC
WRITE( 10, * ) pQuaternionMC
WRITE( 10, * ) pOrientationMC
WRITE( 10, * ) cPositionMC
CLOSE( 10 )

! Sample standard deviation
AverageDensity = AverageDensity - SUM( AverageDensity ) / DBLE( SIZE( AverageDensity ) )
DensityStandardDeviation = SUM( AverageDensity ** 2 )
DensityStandardDeviation = DSQRT( DensityStandardDeviation / DBLE( SIZE( AverageDensity ) - 1 ) )

! Pressure-density curve
IF( Path == 1 ) THEN
  OPEN( Unit= 10, File= "Curve/curve.dat", Position= "APPEND" )
  WRITE( 10, "(G0.9,3(' ',G0.9))" ) Pressure, Density, DensityStandardDeviation, OrderParameterS2
  CLOSE( 10 )
END IF

! Widom insertion
IF( WidomInsertion ) THEN
  WRITE( *, "(G0)" ) "Calculating chemical potential using Widom insertion method..."
  ! Average energy after adding one particle to the system
  wAverageBoltzmannFactor = wAccumulatedBoltzmannFactor / wCounter
  ! Chemical potential (in kT units)
  ResidualChemicalPotential = - DLOG( wAverageBoltzmannFactor )
  WRITE( *, "(G0)" ) " "
  WRITE( *, "(G0,G0.5)" ) "Calculated Residual Chemical Potential (kT units): ", ResidualChemicalPotential
  WRITE( *, "(G0)" ) " "
  WRITE( 70, "(G0)" ) " "
  WRITE( 70, "(G0,G0.5)" ) "Residual Chemical Potential (kT units) = ", ResidualChemicalPotential
  WRITE( 70, "(G0)" ) " "
  FLUSH( 70 )
END IF

! End of Metropolis algorithm
WRITE( *, "(G0)" ) "Monte Carlo simulation finished successfully! See directories for results."
WRITE( *, "(G0)" ) " "

! Output units                                         
IF( PrintTrajectory ) THEN
  CLOSE( 20 )
END IF
CLOSE( 30 )
CLOSE( 40 )
CLOSE( 50 )
CLOSE( 55 )
CLOSE( 70 )
IF( StackRotationLogical ) CLOSE( 80 )
IF( PotentialLogical ) CLOSE( 90 )
CLOSE( 110 )
IF( CellListLogical ) CALL FinalizeList(  )

! Status
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 55 )//CH_UR
WRITE( *, "(G0)" ) CH_VS//REPEAT( " ", 19 )//"SIMULATION LENGTH"//REPEAT( " ", 19 )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 55 )//CH_BR

! End simulation timer
CALL Cpu_Time( StopTimer )
WRITE( *, "(G0,G0.5,G0)" ) "Elapsed Time: ", (StopTimer - StartTimer), "s."
WRITE( *, "(G0)" ) " "

! Deallocation
DEALLOCATE( pQuaternion, pPosition, pOrientation, cPosition, pQuaternionMC, pPositionMC, pPositionSave, pOrientationMC, &
&           cPositionMC, cPositionSave )

END PROGRAM NonconvexMonteCarlo