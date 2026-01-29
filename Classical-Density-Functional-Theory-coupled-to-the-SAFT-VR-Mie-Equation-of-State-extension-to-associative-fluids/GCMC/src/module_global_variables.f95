! ----------------------------------------------------------------------- !
! Global Variables Module                                                 !
! ----------------------------------------------------------------------- !
MODULE GlobalVariables

! Use kind Real64 and Int64
USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: Real64, Int64, Output_Unit

! OpenMP API
#ifdef _OPENMP
USE OMP_LIB
#endif

IMPLICIT NONE

! INTEGER VARIABLES
INTEGER( Kind= Int64 ) :: nParticles      ! Number of particles
INTEGER( Kind= Int64 ) :: nAdjustment     ! Adjustment frequency
INTEGER( Kind= Int64 ) :: nSaveConfig     ! Saving frequency (configuration)
INTEGER( Kind= Int64 ) :: nSave           ! Saving frequency
INTEGER( Kind= Int64 ) :: nCycles         ! Number of MC cycles
INTEGER( Kind= Int64 ) :: nEquilibration  ! Number of equilibration cycles
INTEGER( Kind= Int64 ) :: nEquilRandom    ! Number of equilibration cycles (random configuration)
INTEGER( Kind= Int64 ) :: nThreads        ! Number of threads
INTEGER( Kind= Int64 ) :: nSlabs          ! Number of slabs for density profile
INTEGER( Kind= Int64 ) :: maxAttemptsNPT  ! Maximum number of attempts in the NPT ensemble (random configuration)
INTEGER( Kind= Int64 ) :: nAssocSites     ! Number of associated sites
INTEGER( Kind= Int64 ) :: tAssocSites     ! Total number of possible associated sites
INTEGER( Kind= Int64 ) :: nCells(3)       ! Number of cells along X, Y, and Z
INTEGER( Kind= Int64 ) :: nWidomFrequency ! Frequency of Widom test particle insertions

! INTEGER VARIABLES (ALLOCATABLE)
INTEGER( Kind= Int64 ), ALLOCATABLE :: AssociatedSites( :, :, : ) ! Indexes of the associated site

! REAL CONSTANTS
REAL( Kind= Real64 ), PARAMETER :: cAvogadro = 6.02214076D23                     ! Avogadro's number
REAL( Kind= Real64 ), PARAMETER :: cBoltzmann = 1.38064852D-23                   ! Boltzmann constant in J/K
REAL( Kind= Real64 ), PARAMETER :: cPlanck = 6.62607004D-34                      ! Planck's constant in Joule . second
REAL( Kind= Real64 ), PARAMETER :: cPi = 4.D0 * ATAN( 1.D0 )                     ! Pi number
REAL( Kind= Real64 ), PARAMETER :: UniversalGasConstant = cAvogadro * cBoltzmann ! Ideal gas constant in J/(mol K)
REAL( Kind= Real64 ), PARAMETER :: xAxis(3) = [ 1.D0, 0.D0, 0.D0 ]               ! X axis unit vector
REAL( Kind= Real64 ), PARAMETER :: yAxis(3) = [ 0.D0, 1.D0, 0.D0 ]               ! Y axis unit vector
REAL( Kind= Real64 ), PARAMETER :: zAxis(3) = [ 0.D0, 0.D0, 1.D0 ]               ! Z axis unit vector
REAL( Kind= Real64 ), PARAMETER :: waterMolarMass = 18.01535D-3                  ! Molar mass of water in kg/mol
REAL( Kind= Real64 ), PARAMETER :: waterMass = waterMolarMass / cAvogadro        ! Mass of one water molecule in kg

! REAL VARIABLES
REAL( Kind= Real64 ) :: cellLength(3)                                                         ! Cell lengths
REAL( Kind= Real64 ) :: deBroglieThermalWavelength, sqBroglieWavelength, cbBroglieWavelength  ! De Broglie thermal wavelength
REAL( Kind= Real64 ) :: BoxLength(9), BoxLengthInverse(9), BoxVolume                          ! Box parameters
REAL( Kind= Real64 ) :: poreHeight, poreWidth, poreLength                                     ! Pore dimensions
REAL( Kind= Real64 ) :: SigmaSolid, SigmaFluid, squaredSigmaFluid, invSigmaFluid, &           ! Diameter and related parameters
&                       cubicSigmaFluid                                                       ! Diameter and related parameters
REAL( Kind= Real64 ) :: EpsilonSolid, EpsilonFluid, EpsilonWater, &                           ! Well depth and related parameters
&                       reducedEpsilonFluid, reducedEpsilonWater                              ! Well depth and related parameters
REAL( Kind= Real64 ) :: attractiveRange, repulsiveRange, SWRange, sqSWRange                   ! Potential range
REAL( Kind= Real64 ) :: EpsilonSolid_Fluid, SigmaSolid_Fluid                                  ! Combined parameters
REAL( Kind= Real64 ) :: MolarMass, Density, Temperature, InterlayerSpacing                    ! System properties
REAL( Kind= Real64 ) :: SolidDensity                                                          ! Solid density
REAL( Kind= Real64 ) :: longrangeCorrectionPotential, longrangeCorrectionPressure             ! Long-range corrections
REAL( Kind= Real64 ) :: cutoffRadiusWater, cutoffSquaredRadiusWater                           ! Cutoff parameters (water)
REAL( Kind= Real64 ) :: cutoffRadiusAssociation, cutoffSquaredRadiusAssociation               ! Cutoff parameters (association)
REAL( Kind= Real64 ) :: cutoffRadiusMie, cutoffSquaredRadiusMie, inverseCutoffRadiusMie, &    ! Cutoff parameters (Mie)
&                       inverseCubicCutoffRadiusMie, InverseNinthCutoffRadiusMie              ! Cutoff parameters (Mie)
REAL( Kind= Real64 ) :: NormalisationCoefficient, ratioRange, ratioRepulsiveRange, &          ! Potential parameters
&                       ratioAttractiveRange                                                  ! Potential parameters
REAL( Kind= Real64 ) :: maxTranslationalDisplc, maxRotationalDisplc                           ! Maximum displacements
REAL( Kind= Real64 ) :: maxTranslationalDisplcRandom                                          ! Maximum displacements (random configuration)
REAL( Kind= Real64 ) :: maxRotationalDisplcRandom                                             ! Maximum displacements (random configuration)
REAL( Kind= Real64 ) :: BoxEdgeMaxRatio                                                       ! Maximum box edge ratio (random configuration)
REAL( Kind= Real64 ) :: MaxIsoVolumetricDisplacementRandom                                    ! Maximum displacements (random configuration)
REAL( Kind= Real64 ) :: MaxAnisoVolumetricDisplacementRandom                                  ! Maximum displacements (random configuration)
REAL( Kind= Real64 ) :: MovementProbability, TranslationalProbability, RotationalProbability  ! Probabilities
REAL( Kind= Real64 ) :: VolumeChangeProbability, IsotropicProbability, AnisotropicProbability ! Probabilities
REAL( Kind= Real64 ) :: CreationProbability, DestructionProbability, PCDProbability           ! Probabilities
REAL( Kind= Real64 ) :: TranslationalProbabilityXY                                            ! Probabilities
REAL( Kind= Real64 ) :: AcceptanceRatioTranslation, AcceptanceRatioRotation, &                ! Acceptance ratio thresholds
&                       AcceptanceRatioIsotropic, AcceptanceRatioAnisotropic, &               ! Acceptance ratio thresholds
&                       AcceptanceRatioCreation, AcceptanceRatioDestruction                   ! Acceptance ratio thresholds
REAL( Kind= Real64 ) :: SiteRadius, AngleWater, ReferencePositions( 3, 4 ), &                 ! Four-site model parameters (OOHH)
&                       ReferenceQuaternions( 0:3, 4 )                                        ! Four-site model parameters (OOHH)
REAL( Kind= Real64 ) :: RandomNumber                                                          ! Random number
REAL( Kind= Real64 ) :: TargetPressureRandom, ReducedTargetPressure                           ! Target pressure (random configuration)
REAL( Kind= Real64 ) :: ReducedChemicalPotential                                              ! Chemical potential
REAL( Kind= Real64 ) :: calculatedReducedChemPot, excessReducedChemPot                        ! Chemical potential
REAL( Kind= Real64 ) :: GrandCanonicalConstant                                                ! Grand canonical constant
REAL( Kind= Real64 ) :: fractionAssociatedSites, fractionFreeSites                            ! Fraction of associated and free sites

! REAL VARIABLES (ALLOCATABLE)
REAL( Kind= Real64 ), ALLOCATABLE :: pPosition( :, : ), pQuaternion( :, : ), pOrientation( :, : ), sitePosition( :, :, : ), &
&                                    siteOrientation( :, :, : ), siteQuaternion( :, :, : ), pPositionRandom( :, : ), &
&                                    pQuaternionRandom( :, : ), pOrientationRandom( :, : ), sitePositionRandom( :, :, : ), &
&                                    siteOrientationRandom( :, :, : ), siteQuaternionRandom( :, :, : ) ! Particle and sites positions and orientations

! CHARACTER STRINGS
CHARACTER( Len= 1 ) :: Dummy        ! Dummy variable for reading
CHARACTER( Len= 1 ) :: SiteNames(4) ! Site names
CHARACTER( LEN= 1 ) :: LineIndex    ! First character (index) of a line

! LOGICAL VARIABLES
LOGICAL :: RandomSeed                    ! Check for random seed for the random number generator
LOGICAL :: RandomConf                    ! Check for random configuration generation
LOGICAL :: Confinement                   ! Check for confinement
LOGICAL :: ApplyNPTRandom                ! Apply NPT moves in random configuration
LOGICAL :: FileExists                    ! Check if file exists
LOGICAL :: isSolidPotentialEnabled       ! Check if solid potential is enabled
LOGICAL :: isMiePotentialEnabled         ! Check if attractive potential is enabled
LOGICAL :: isSWPotentialEnabled          ! Check if attractive potential is enabled
LOGICAL :: isAssociativePotentialEnabled ! Check if associative potential is enabled
LOGICAL :: isHardSpherePotentialEnabled  ! Check if hard sphere potential is enabled
LOGICAL :: isWidomEnabled                ! Check if Widom's insertion is enabled
LOGICAL :: OverwriteCombiningRule        ! Check if the combining rule for diameter and potential depth will be overwritten by the user
LOGICAL :: xyTranslationIndependency     ! Check if there is translation dependency in the XY plane

LOGICAL, ALLOCATABLE :: isSiteAssociated( :, : ) ! Checks if a site is associated

END MODULE GlobalVariables