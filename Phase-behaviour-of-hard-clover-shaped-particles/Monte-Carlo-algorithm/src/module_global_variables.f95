MODULE GlobalVariables

! Use kind Real64 and Int64
USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: Real64, Int64, Output_Unit

! OpenMP API
#ifdef _OPENMP
USE OMP_LIB
#endif

IMPLICIT NONE

! INTEGER VARIABLES -*- THIS IS SINGLE PRECISION -*-
INTEGER, DIMENSION( 8 ) :: DateTime ! Computer clock (date and time)

! INTEGER VARIABLES (SCALAR,CONSTANT)
INTEGER( Kind= Int64 ), PARAMETER :: HalfNeighboursControl = 1 ! Checks whether a cell and its 26 surrounding cells are searched or just its 13 neighbour cells

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 ) :: nParticles                   ! Number of particles
INTEGER( Kind= Int64 ) :: pGeometry                    ! Surface geometry
INTEGER( Kind= Int64 ) :: MaxSimulationCyclesInit      ! Maximum number of cycles (Optimization)
INTEGER( Kind= Int64 ) :: nThreads                     ! Number of threads (OpenMP)
INTEGER( Kind= Int64 ) :: Path                         ! Path ([1]: NVT/NPT Ramp Simulation, [2]: Conventional NVT/NPT Simulation)
INTEGER( Kind= Int64 ) :: nSimulationCycles            ! Total number of cycles
INTEGER( Kind= Int64 ) :: nEquilibration               ! Number of equilibration cycles
INTEGER( Kind= Int64 ) :: nSavingFrequency             ! Results saving frequency
INTEGER( Kind= Int64 ) :: nAdjustmentMovementFrequency ! Maximum molecular displacement adjustment frequency
INTEGER( Kind= Int64 ) :: nAdjustmentVolumeFrequency   ! Maximum volumetric change adjustment frequency
INTEGER( Kind= Int64 ) :: nWidomCycles                 ! Total number of cycles regarding the Widom insertion method

! INTEGER VARIABLES (ARRAY)
INTEGER( Kind= Int64 ), DIMENSION( 3 ) :: pCells ! Number of cells

! INTEGER VARIABLES (ARRAY,ALLOCATABLE)
INTEGER( Kind= Int64 ), DIMENSION( : ), ALLOCATABLE       :: pCellList  ! Cell list
INTEGER( Kind= Int64 ), DIMENSION( :, : ), ALLOCATABLE    :: pCellIndex ! 3D-cell index of each particle
INTEGER( Kind= Int64 ), DIMENSION( :, :, : ), ALLOCATABLE :: pCellHead  ! Cell head

! REAL VARIABLES (SCALAR,CONSTANT)
REAL( Kind= Real64 ), PARAMETER :: cPi = 4.D0 * DATAN( 1.D0 ) ! π
REAL( Kind= Real64 ), PARAMETER :: cBoltzmann = 1.380649D-23  ! m².kg/(s².K)

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 ) :: QuaternionAngle                    ! Quaternion angle [real part, W] (for initial configuration only)
REAL( Kind= Real64 ) :: BoxVolume                          ! Volume of simulation box
REAL( Kind= Real64 ) :: Density                            ! Reduced number density
REAL( Kind= Real64 ) :: Pressure                           ! Reduced pressure
REAL( Kind= Real64 ) :: PackingFraction                    ! Packing fraction
REAL( Kind= Real64 ) :: PotentialStrength                  ! Strength of the external potential (orientational field)
REAL( Kind= Real64 ) :: StackRotationDensityThreshold      ! Threshold density for stack rotation
REAL( Kind= Real64 ) :: Temperature                        ! Temperature
REAL( Kind= Real64 ) :: cDiameter                          ! Diameter (cylinder)
REAL( Kind= Real64 ) :: cLength                            ! Length (cylinder)
REAL( Kind= Real64 ) :: pVolume                            ! Particle volume
REAL( Kind= Real64 ) :: cAspectRatio                       ! Aspect ratio (cylinder)
REAL( Kind= Real64 ) :: TotalParticleVolume                ! Total particle volume
REAL( Kind= Real64 ) :: RandomNumber                       ! Random number from a pseudorandom number generator subroutine
REAL( Kind= Real64 ) :: MovementProbability                ! Probability of movement
REAL( Kind= Real64 ) :: TranslationalProbability           ! Probability of movement (translation)
REAL( Kind= Real64 ) :: RotationalProbability              ! Probability of movement (rotation)
REAL( Kind= Real64 ) :: VolumeChangeProbability            ! Probability of volume change
REAL( Kind= Real64 ) :: IsoVolumetricProbability           ! Probability of volume change (isotropic)
REAL( Kind= Real64 ) :: AnisoVolumetricProbability         ! Probability of volume change (anisotropic)
REAL( Kind= Real64 ) :: StackRotationProbability           ! Probability of stack rotation
REAL( Kind= Real64 ) :: CutoffSphere                       ! Cutoff diameter of the circumscribing sphere (particle)
REAL( Kind= Real64 ) :: cCutoffSphere                      ! Cutoff diameter of the circumscribing sphere (cylinder)
REAL( Kind= Real64 ) :: SquaredCutoffSphere                ! Squared cutoff diameter of the circumscribing sphere (particle)
REAL( Kind= Real64 ) :: cSquaredCutoffSphere               ! Squared cutoff diameter of the circumscribing sphere (cylinder)
REAL( Kind= Real64 ) :: CutoffSpherocylinder               ! Cutoff diameter of the circumscribing spherocylinder (particle)
REAL( Kind= Real64 ) :: SquaredCutoffSpherocylinder        ! Squared cutoff diameter of the circumscribing spherocylinder (particle)
REAL( Kind= Real64 ) :: MaxLengthRatio                     ! Maximum length ratio of simulation box during anisotropic volume changes
REAL( Kind= Real64 ) :: MaxBoxAngle                        ! Maximum angle between box vectors during anisotropic volume changes
REAL( Kind= Real64 ) :: BoxFactor                          ! Volume factor of the simulation box for the packed-box configuration
REAL( Kind= Real64 ) :: cLargestSphereDiameter             ! Diameter of the largest circumscribing sphere
REAL( Kind= Real64 ) :: UserMaxTranslationalDisplacement   ! User maximum displacement [+/-] (Translation)
REAL( Kind= Real64 ) :: UserMaxRotationalDisplacement      ! User maximum displacement [+/-] (Rotation)
REAL( Kind= Real64 ) :: UserMaxIsoVolumetricDisplacement   ! User maximum displacement [+/-] (Isotropic volume change)
REAL( Kind= Real64 ) :: UserMaxAnisoVolumetricDisplacement ! User maximum displacement [+/-] (Anisotropic volume change)
REAL( Kind= Real64 ) :: AcceptanceRatioTranslation         ! Acceptance ratio threshold (Translation)
REAL( Kind= Real64 ) :: AcceptanceRatioRotation            ! Acceptance ratio threshold (Rotation)
REAL( Kind= Real64 ) :: AcceptanceRatioIsoVolumeChange     ! Acceptance ratio threshold (Isotropic volume change)
REAL( Kind= Real64 ) :: AcceptanceRatioAnisoVolumeChange   ! Acceptance ratio threshold (Anisotropic volume change)
REAL( Kind= Real64 ) :: MaxBoxDistortion                   ! Maximum box distortion before lattice reduction

! REAL VARIABLES (ARRAY,CONSTANT)
REAL( Kind= Real64 ), DIMENSION( 3 ), PARAMETER    :: xAxis = [1.D0, 0.D0, 0.D0] ! Body-fixed axis of rotation (x)
REAL( Kind= Real64 ), DIMENSION( 3 ), PARAMETER    :: yAxis = [0.D0, 1.D0, 0.D0] ! Body-fixed axis of rotation (y)
REAL( Kind= Real64 ), DIMENSION( 3 ), PARAMETER    :: zAxis = [0.D0, 0.D0, 1.D0] ! Body-fixed axis of rotation (z)
REAL( Kind= Real64 ), DIMENSION( 3, 3 ), PARAMETER :: IdentityMatrix = RESHAPE( [ 1.D0, 0.D0, 0.D0, &
&                                                                                 0.D0, 1.D0, 0.D0, &
&                                                                                 0.D0, 0.D0, 1.D0 ], [ 3, 3 ] ) ! Identity matrix

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: OrientationalFieldVector ! Orientation of the external potential field
REAL( Kind= Real64 ), DIMENSION( 9 )    :: BoxLength                ! Length (x,y,z) of simulation box
REAL( Kind= Real64 ), DIMENSION( 9 )    :: BoxLengthInverse         ! Length (x,y,z) of simulation box (inverse)
REAL( Kind= Real64 ), DIMENSION( 3, 4 ) :: cPositionBasis           ! Body-fixed reference position

! REAL VARIABLES (ARRAY,ALLOCATABLE)
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE    :: pQuaternion, pQuaternionMC   ! Quaternion array
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE    :: pPosition, pPositionMC       ! Position array
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE    :: pOrientation, pOrientationMC ! Orientation array
REAL( Kind= Real64 ), DIMENSION( :, :, : ), ALLOCATABLE :: cPosition, cPositionMC       ! Position array (cylinders)

! CHARACTER VARIABLES
CHARACTER( Len= 10  ) :: DescriptorAR           ! Descriptor for output file
CHARACTER( Len= 10  ) :: DescriptorPressure     ! Descriptor for output file
CHARACTER( Len= 10  ) :: DescriptorDate         ! Descriptor for output folder
CHARACTER( Len= 10  ) :: DescriptorHour         ! Descriptor for output file
CHARACTER( Len= 32  ) :: FormatAR               ! String format for output file
CHARACTER( Len= 32  ) :: FormatPressure         ! String format for output file
CHARACTER( Len= 32  ) :: FormatDate             ! String format for output folder
CHARACTER( Len= 32  ) :: FormatHour             ! String format for output file
CHARACTER( Len= 03  ) :: ConfigurationInquiry   ! Initial configuration inquiry
CHARACTER( Len= 03  ) :: LatticeReductionMethod ! Lattice reduction type
CHARACTER( Len= 01  ) :: Dummy                  ! Dummy

! CHARACTER STRINGS (BOX FRAMES)
CHARACTER( Len= 3 ), PARAMETER :: CH_HS = "═" ! Box drawing symbol
CHARACTER( Len= 3 ), PARAMETER :: CH_VS = "║" ! Box drawing symbol
CHARACTER( Len= 3 ), PARAMETER :: CH_UL = "╔" ! Box drawing symbol
CHARACTER( Len= 3 ), PARAMETER :: CH_BL = "╚" ! Box drawing symbol
CHARACTER( Len= 3 ), PARAMETER :: CH_UR = "╗" ! Box drawing symbol
CHARACTER( Len= 3 ), PARAMETER :: CH_BR = "╝" ! Box drawing symbol
CHARACTER( Len= 3 ), PARAMETER :: CH_VL = "╠" ! Box drawing symbol
CHARACTER( Len= 3 ), PARAMETER :: CH_VR = "╣" ! Box drawing symbol
CHARACTER( Len= 3 ), PARAMETER :: SH_VL = "╟" ! Box drawing symbol
CHARACTER( Len= 3 ), PARAMETER :: SH_VR = "╢" ! Box drawing symbol
CHARACTER( Len= 3 ), PARAMETER :: SS_HS = "─" ! Box drawing symbol
CHARACTER( Len= 3 ), PARAMETER :: SS_VS = "│" ! Box drawing symbol
CHARACTER( Len= 3 ), PARAMETER :: SS_VR = "├" ! Box drawing symbol
CHARACTER( Len= 3 ), PARAMETER :: SS_VL = "┤" ! Box drawing symbol
CHARACTER( Len= 3 ), PARAMETER :: SS_UL = "┌" ! Box drawing symbol
CHARACTER( Len= 3 ), PARAMETER :: SS_BL = "└" ! Box drawing symbol
CHARACTER( Len= 3 ), PARAMETER :: SS_UR = "┐" ! Box drawing symbol
CHARACTER( Len= 3 ), PARAMETER :: SS_BR = "┘" ! Box drawing symbol
CHARACTER( Len= 3 ), PARAMETER :: C_FUL = "█" ! Box drawing symbol

! LOGICAL VARIABLES (SCALAR)
LOGICAL :: FileExist            ! Checks whether a file exists or not
LOGICAL :: FolderExist          ! Checks whether folder exists or not
LOGICAL :: PrintTrajectory      ! Checks whether the trajectory will be printed
LOGICAL :: RestoreBackup        ! Checks whether a previous simulation will be restored from a backup file
LOGICAL :: RandomSeed           ! Checks whether a random seed will be used
LOGICAL :: WidomInsertion       ! Checks whether the Widom insertion method will be used
LOGICAL :: StackRotationLogical ! Checks whether the stack rotation moves will be applied
LOGICAL :: PotentialLogical     ! Checks whether the external potential field will be applied
LOGICAL :: CellListLogical      ! Checks whether a cell list will be used to compute the potential (NVT only) or search for overlaps
LOGICAL :: CellListControl      ! Toggles control of cell list

! LOGICAL VARIABLES (ARRAY)
LOGICAL, DIMENSION( 2 ) :: ConfigurationLogical        ! Selected molecular configuration
LOGICAL, DIMENSION( 2 ) :: LatticeReductionTypeLogical ! Lattice reduction selection
LOGICAL, DIMENSION( 2 ) :: GeometrySelection           ! Selected surface geometry

END MODULE GlobalVariables