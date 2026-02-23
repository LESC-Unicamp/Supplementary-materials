MODULE GlobalVariables

! Use intrinsic module: ISO_FORTRAN_ENV (for real and integer kinds)
USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : Real64, Int64

! OpenMP API
#ifdef _OPENMP
USE OMP_LIB
#endif

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 ) :: nParticles                   ! Number of particles
INTEGER( Kind= Int64 ) :: nSave                        ! Saving frequency (frames)
INTEGER( Kind= Int64 ) :: MaxCycles                    ! Maximum number of cycles
INTEGER( Kind= Int64 ) :: nEquilibration               ! Number of equilibration cycles
INTEGER( Kind= Int64 ) :: nThreads                     ! Number of threads (parallelization)
INTEGER( Kind= Int64 ) :: MaxNeighbors                 ! Maximum number of neighbours (for calculating the bond parameters, Ψ)
INTEGER( Kind= Int64 ) :: nLayers                      ! Number of layers of minimum images
INTEGER( Kind= Int64 ) :: nSkip                        ! Skipping frequency (frames)
INTEGER( Kind= Int64 ) :: nBins                        ! Number of bins (for the radial distribution functions)
INTEGER( Kind= Int64 ) :: MaxSimulatedAnnealingCounter ! Maximum number of iterations (simulated annealing)
INTEGER( Kind= Int64 ) :: MaxRepetitions               ! Maximum number of repetitions (simulated annealing)
INTEGER( Kind= Int64 ) :: nNeighbors                   ! Number of nearest neighbours (Steinhardt order parameters)
INTEGER( Kind= Int64 ) :: NeighborOption               ! Option to compute the neighbors (Steinhardt order parameters)
INTEGER( Kind= Int64 ) :: NeighborOption2D             ! Option to compute the neighbors (local bond parameters)
INTEGER( Kind= Int64 ) :: MaximumSteinhardt            ! Maximum value of the Steinhardt order parameter (ℓ)
INTEGER( Kind= Int64 ) :: MinimumSteinhardt            ! Minimum value of the Steinhardt order parameter (ℓ)
INTEGER( Kind= Int64 ) :: MaxVertices                  ! Maximum number of vertices (for the Voronoi construction)
INTEGER( Kind= Int64 ) :: nStruc                       ! Size of the structure factor grid
INTEGER( Kind= Int64 ) :: cGlobalSimplex               ! Number of cycles for the simplex optimization (Nelder-Mead)

! REAL VARIABLES (CONSTANTS)
REAL( Kind= Real64 ), PARAMETER :: cPi = 4.D0 * DATAN( 1.D0 ) ! π

! REAL VARIABLES (CONSTANTS,ARRAY)
REAL( Kind= Real64 ), DIMENSION( 3 ), PARAMETER    :: BodyFixedAxis = [ 0.D0, 0.D0, 1.D0 ] ! Body-fixed axis (Z-axis)
REAL( Kind= Real64 ), DIMENSION( 3, 3 ), PARAMETER :: IdentityMatrix = RESHAPE( [ 1.D0, 0.D0, 0.D0, &
&                                                                                 0.D0, 1.D0, 0.D0, &
&                                                                                 0.D0, 0.D0, 1.D0 ], [ 3, 3 ] ) ! Identity matrix

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 ) :: InitialLayerSpacing   ! Initial layer spacing for the smectic phase (estimation)
REAL( Kind= Real64 ) :: MinLayerSpacing       ! Minimum layer spacing for the smectic phase (optimization)
REAL( Kind= Real64 ) :: MaxLayerSpacing       ! Maximum layer spacing for the smectic phase (optimization)
REAL( Kind= Real64 ) :: BoxVolume             ! Volume of the simulation box
REAL( Kind= Real64 ) :: pDiameter             ! Diameter of the cylinders
REAL( Kind= Real64 ) :: pDiameterClover       ! Diameter of the clover (particle)
REAL( Kind= Real64 ) :: pLength               ! Length of the cylinders
REAL( Kind= Real64 ) :: pLengthClover         ! Length of the clover (particle)
REAL( Kind= Real64 ) :: CappingFactor         ! Capping factor of the structure factor
REAL( Kind= Real64 ) :: bWidth                ! Width of the bins
REAL( Kind= Real64 ) :: RandomNumber          ! Random number
REAL( Kind= Real64 ) :: InitialTemp           ! Initial temperature (simulated annealing)
REAL( Kind= Real64 ) :: FinalTemp             ! Final temperature (simulated annealing)
REAL( Kind= Real64 ) :: TempScale             ! Temperature scaling (simulated annealing)
REAL( Kind= Real64 ) :: MaxAngularDisplc      ! Maximum angular displacement (simulated annealing)
REAL( Kind= Real64 ) :: CutoffDistance        ! Cutoff distance for the neighbors (Steinhardt order parameters)
REAL( Kind= Real64 ) :: CutoffDistanceVoronoi ! Cutoff distance for the neighbors (Voronoi construction)
REAL( Kind= Real64 ) :: ParallelRange         ! Extension of the cylindrical region (radial) for the calculation of the parallel radial distribution function
REAL( Kind= Real64 ) :: OrthogonalRange       ! Extension of the cylindrical region (axial) for the calculation of the orthogonal radial distribution function

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 9 )    :: BoxLength          ! Box length
REAL( Kind= Real64 ), DIMENSION( 9 )    :: BoxLengthInverse   ! Inverse of the box length
REAL( Kind= Real64 ), DIMENSION( 4, 3 ) :: cReferencePosition ! Body-fixed position

! REAL VARIABLES (ARRAY,ALLOCATABLE)
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE :: pPosition    ! Position of the particles
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE :: pQuaternion  ! Quaternion of the particles
REAL( Kind= Real64 ), DIMENSION( :, : ), ALLOCATABLE :: pOrientation ! Orientation of the particles

! LOGICAL VARIABLES
LOGICAL :: ComputeLocalBondParameters ! Compute local bond parameters
LOGICAL :: ComputeParameterEquil      ! Include equilibration data in the calculation of the order parameters and other properties of the system
LOGICAL :: ComputeRDF                 ! Compute radial distribution function
LOGICAL :: ComputeStructureFactor     ! Compute structure factor
LOGICAL :: ComputeGlobalOrientOrder   ! Compute global orientational order parameter
LOGICAL :: ComputeImages              ! Include images of particles in the plane when calculating the structure factor
LOGICAL :: ComputeCubaticOrder        ! Compute cubatic order parameter
LOGICAL :: ComputeNematicOrder        ! Compute nematic order parameter
LOGICAL :: ComputeSmecticOrder        ! Compute smectic order parameter
LOGICAL :: ComputeSteinhardtOrder     ! Compute Steinhardt order parameters (Qℓ)
LOGICAL :: ApplyCapping               ! Apply capping factor to the structure factor
LOGICAL :: AveragedSteinhardt         ! Compute averaged Steinhardt order parameters (includes the neighbours of the neighbours of the central particle)
LOGICAL :: ComputeWignerParameter     ! Compute Wigner's parameters

! CHARACTER STRINGS
CHARACTER( Len= 100 ) :: FileName ! File name
CHARACTER( Len= 12 )  :: Dummy    ! Dummy

END MODULE GlobalVariables 