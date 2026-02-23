MODULE Initialization

! Uses one module: global variables
USE GlobalVariables

CONTAINS

! ********************************************************************************************************************************* !
! This subroutine reads the file "<initialization.ini>" and initializes the simulation parameters.                                  !
!                                                                                                                                   !
! These parameters include:                                                                                                         !
!                                                                                                                                   !
! ● Parameters necessary for reading the trajectory file: the file name <FileName>, the number of simulated particles <nParticles>, !
!        the number of simulation cycles <MaxCycles>, the number of equilibration cycles <nEquilibration>, the saving               !
!        frequency <nSave>, and the skipping frequency <nSkip>. The skipping frquency is used to skip lines in the trajectory file. !
!                                                                                                                                   !
! ● Parameters of the particles: the length <pLength> and diameter <pDiameter> of the cylinders (petals).                           !
!                                                                                                                                   !
! ● Logical parameters for the analysis: <ComputeParameterEquil> checks whether equilibration data is included in the analysis,     !
!        <ComputeNematicOrder> checks whether the nematic order parameter is computed, <ComputeCubaticOrder> checks whether the     !
!        cubatic order parameter is computed, <ComputeSmecticOrder> checks whether the smectic order parameter is computed,         !
!        <ComputeStructureFactor> checks whether the structure factor is computed, <ComputeLocalBondParameters> checks whether the  !
!        local bond parameters (Ψ) are computed, <ComputeSteinhardtOrder> checks whether the Steinhardt order parameters (Qℓ) are   !
!        computed, and <ComputeRDF> checks whether the radial distribution functions (RDF) are computed, including the parallel and !
!        perpendicular components and the distribution of the Legendre polynomials (P₂ and P₄).                                     !
!                                                                                                                                   !
! ● Parameter for the calculation of the local bond parameters (Ψ): the maximum number of neighbors <MaxNeighbors> in a plane (2D). !
!        The Voronoi diagram can also be used to calculate the neighbors by setting <NeighborOption2D> = 2.                         !
!                                                                                                                                   !
! ● Parameters for the calculation of the structure factor: <ComputeImages> is a logical variable that determines whether the images!
!        of the central particles are computed to generate a plane and <nLayers> defines how many periodic layers are considered.   !
!        <ApplyCapping> is another logical variable that determines whether a capping is applied to the structure factor, and       !
!        <CappingFactor> is the capping threshold. The structure factor is calculated on a two-dimensional grid with <nStruc>       !
!        points along ±x- and ±y-directions.                                                                                        !
!                                                                                                                                   !
! ● Parameters for the calculation of the radial distribution functions: <nBins> is the number of bins and <bWidth> is the width    !
!        of the bins. The parallel RDF is by default limited to a cylindrical region whose radial extension is equal to the         !
!        diameter of the particle. The orthogonal RDF is by default limited to a cylindrical region whose axial extension is        !
!        equal to the length of the particle. <ParallelRange> and <OrthogonalRange> are parameters that extend the volume of this   !
!        cylindrical region along the radial and axial directions, respectively. Use <ParallelRange> = 1 and <OrthogonalRange> = 1  !
!        to use the default volumes of the cylindrical region.                                                                      !
!                                                                                                                                   !
! ● Parameters for the calculation of the cubatic order parameter using simulated annealing: <InitialTemp> and <FinalTemp> are      !
!        the initial and final temperatures, <TempScale> is the temperature scaling factor, <MaxSimulatedAnnealingCounter> is the   !
!        number of steps, <MaxAngularDisplc> is the maximum angular displacement used to randomly generate a test unit              !
!        quaternion, and <MaxRepetitions> is the maximum number of repetitions.                                                     !
!                                                                                                                                   !
! ● Parameters for the calculation of the smectic order parameter: <InitialLayerSpacing> is the initial estimate of the layer       !
!        spacing, while <MaxLayerSpacing> and <MinLayerSpacing> are the maximum and minimum layer spacings, respectively, used      !
!        in the Nelder-Mead simplex algorithm to optimize the layer spacing. The optimization is performed <cGlobalSimplex> times.  !
!                                                                                                                                   !
! ● Parameters for the calculation of the Steinhardt order parameters: <NeighborOption> is the option for selecting the closest     !
!        neighbors in a 3D space [1: for the first <nNeighbors> nearest neighbors, 2: for nearest neighbors within a cutoff         !
!        distance <CutoffDistance>, and 3: for the nearest neighbors in a Voronoi diagram considering a cutoff distance             !
!        <CutoffDistanceVoronoi>]. For the option 3, <MaxVertices> defines the maximum number of vertices in a Voronoi cell. For a  !
!        more accurate calculation of the Steinhardt order parameters, consider enabling <AveragedSteinhardt> (q̄ℓ) which will also  !
!        include the neighbors of the neighbors of the central particle, and <ComputeWignerParameter> (wℓ or w̄ℓ) which will include !
!        the Wigner's parameters. The rotational invariants qℓ are computed for ℓ ranging from <MinimumSteinhardt> to               !
!        <MaximumSteinhardt>.                                                                                                       !
!                                                                                                                                   !
! ********************************************************************************************************************************* !
SUBROUTINE Initialization_Parameters(  )

IMPLICIT NONE

! CHARACTER STRINGS
CHARACTER( LEN= 1 ) :: LineIndex ! First character (index) of a line

! File name
OPEN( Unit= 10, File= "initialization.ini", Action= "READ" )

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! File name
READ( 10, * ) Dummy, FileName

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize the number of particles
READ( 10, * ) Dummy, nParticles
IF( nParticles <= 0 ) THEN
  WRITE( *, "(3G0)" ) "The number of particles [", nParticles, "] must be greater than 0!"
  STOP
END IF

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize the number of simulation cycles
READ( 10, * ) Dummy, MaxCycles
IF( MaxCycles <= 0 ) THEN
  WRITE( *, "(3G0)" ) "The number of cycles [", MaxCycles, "] must be greater than 0!"
  STOP
END IF

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize the number of equilibration cycles
READ( 10, * ) Dummy, nEquilibration
IF( nEquilibration < 0 ) THEN
  WRITE( *, "(3G0)" ) "The number of equilibration cycles [", nEquilibration, "] must be greater than or equal to 0!"
  STOP
END IF
IF( nEquilibration >= MaxCycles ) THEN
  WRITE( *, "(5G0)" ) "The number of equilibration cycles [", nEquilibration, "] must be less than the number of cycles [", &
  &                   MaxCycles, "]!"
  STOP
END IF

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize the saving frequency
READ( 10, * ) Dummy, nSave
IF( nSave < 0 ) THEN
  WRITE( *, "(3G0)" ) "The saving frequency [", nSave, "] must be greater than 0!"
  STOP
END IF

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize the length of the cylinders (petals)
READ( 10, * ) Dummy, pLength
IF( pLength < 0.D0 ) THEN
  WRITE( *, "(3G0)" ) "The length of the cylinders [", pLength, "] must be greater than 0!"
  STOP
END IF

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize the diameter of the cylinders (petals)
READ( 10, * ) Dummy, pDiameter
IF( pDiameter < 0.D0 ) THEN
  WRITE( *, "(3G0)" ) "The diameter of the cylinders [", pDiameter, "] must be greater than 0!"
  STOP
END IF

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize the logical variable that determines whether the equilibration data is included in the analysis
READ( 10, * ) Dummy, ComputeParameterEquil

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize the skipping frequency
READ( 10, * ) Dummy, nSkip
IF( nSkip < 0 ) THEN
  WRITE( *, "(3G0)" ) "The skipping frequency [", nSkip, "] must be greater than or equal to 0!"
  STOP
END IF

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize the logical variable that determines whether the nematic order parameter is computed
READ( 10, * ) Dummy, ComputeNematicOrder

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize the logical variable that determines whether the structure factor is computed
READ( 10, * ) Dummy, ComputeStructureFactor

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize the logical variable that determines whether the global orientational order parameter is computed
READ( 10, * ) Dummy, ComputeGlobalOrientOrder

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize the logical variable that determines whether the local bond parameters (Ψ) are computed
READ( 10, * ) Dummy, ComputeLocalBondParameters

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize the logical variable that determines whether the radial distribution functions (RDF) are computed
READ( 10, * ) Dummy, ComputeRDF

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize the logical variable that determines whether the cubatic order parameter is computed
READ( 10, * ) Dummy, ComputeCubaticOrder

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize the logical variable that determines whether the smectic order parameter is computed
READ( 10, * ) Dummy, ComputeSmecticOrder

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize the logical variable that determines whether the Steinhardt order parameters (Qℓ) are computed
READ( 10, * ) Dummy, ComputeSteinhardtOrder

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize the option for the calculation of neighbors in 2D (local bond parameters)
READ( 10, * ) Dummy, NeighborOption2D
IF( NeighborOption2D < 1 .OR. NeighborOption2D > 2 ) THEN
  WRITE( *, "(3G0)" ) "The option for the calculation of neighbors in 2D [", NeighborOption2D, "] must be 1 or 2!"
  STOP
END IF

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize the maximum number of neighbours for calculating the bond parameters (Ψ)
READ( 10, * ) Dummy, MaxNeighbors
IF( MaxNeighbors < 4 ) THEN
  WRITE( *, "(3G0)" ) "The number of neighbors [", MaxNeighbors, "] must be greater than 3!"
  STOP
END IF

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize the logical variable that determines whether the images of the central particles are computed
READ( 10, * ) Dummy, ComputeImages

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize number of layers of periodic images
READ( 10, * ) Dummy, nLayers
IF( nLayers < 1 ) THEN
  WRITE( *, "(3G0)" ) "The number of layers [", nLayers, "] must be greater than 0!"
  STOP
END IF

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize the logical variable that determines whether the capping is applied to the structure factor
READ( 10, * ) Dummy, ApplyCapping

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize the capping of the structure factor
READ( 10, * ) Dummy, CappingFactor
IF( CappingFactor < 0.D0 ) THEN
  WRITE( *, "(3G0)" ) "The capping factor [", CappingFactor, "] must be greater than 0!"
  STOP
END IF

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize the size of the two-dimensional grid for the structure factor
READ( 10, * ) Dummy, nStruc
IF( nStruc <= 0 ) THEN
  WRITE( *, "(3G0)" ) "The size of the structure factor grid [", nStruc, "] must be greater than 0!"
  STOP
END IF

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize number of bins (for the radial distribution functions)
READ( 10, * ) Dummy, nBins
IF( nBins <= 1 ) THEN
  WRITE( *, "(3G0)" ) "The number of bins [", nBins, "] must be greater than 1!"
  STOP
END IF

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize width of the bins
READ( 10, * ) Dummy, bWidth
IF( bWidth <= 0.D0 ) THEN
  WRITE( *, "(3G0)" ) "The bin width [", bWidth, "] must be greater than 0!"
  STOP
END IF

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize the parallel range for the radial distribution functions
READ( 10, * ) Dummy, ParallelRange
IF( ParallelRange <= 0.D0 ) THEN
  WRITE( *, "(3G0)" ) "The extension factor of the cylindrical region along the radial direction [", ParallelRange, &
  &                   "] must be greater than 0!"
  STOP
END IF

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize the orthogonal range for the radial distribution functions
READ( 10, * ) Dummy, OrthogonalRange
IF( OrthogonalRange <= 0.D0 ) THEN
  WRITE( *, "(3G0)" ) "The extension factor of the cylindrical region along the axial direction [", OrthogonalRange, &
  &                   "] must be greater than 0!"
  STOP
END IF

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize the initial temperature for the simulated annealing
READ( 10, * ) Dummy, InitialTemp
IF( InitialTemp < 0.D0 ) THEN
  WRITE( *, "(3G0)" ) "The initial temperature [", InitialTemp, "] must be greater than or equal to 0!"
  STOP
END IF

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize the final temperature for the simulated annealing
READ( 10, * ) Dummy, FinalTemp
IF( FinalTemp < 0.D0 ) THEN
  WRITE( *, "(3G0)" ) "The final temperature [", FinalTemp, "] must be greater than or equal to 0!"
  STOP
END IF
IF( FinalTemp >= InitialTemp ) THEN
  WRITE( *, "(5G0)" ) "The final temperature [", FinalTemp, "] must be less than the initial temperature [", InitialTemp, "]!"
  STOP
END IF

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize the temperature scaling factor for the simulated annealing
READ( 10, * ) Dummy, TempScale
IF( TempScale <= 0.D0 .OR. TempScale >= 1.D0 ) THEN
  WRITE( *, "(3G0)" ) "The temperature scaling factor [", TempScale, "] must be greater than 0 and less than 1!"
  STOP
END IF

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize the maximum number of simulated annealing steps (iterations)
READ( 10, * ) Dummy, MaxSimulatedAnnealingCounter
IF( MaxSimulatedAnnealingCounter < 0 ) THEN
  WRITE( *, "(3G0)" ) "The maximum number of simulated annealing steps [", MaxSimulatedAnnealingCounter, &
  &                   "] must be greater than or equal to 0!"
  STOP
END IF

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize the maximum angular displacement for the simulated annealing
READ( 10, * ) Dummy, MaxAngularDisplc
IF( MaxAngularDisplc <= 0.D0 ) THEN
  WRITE( *, "(3G0)" ) "The maximum angular displacement [", MaxAngularDisplc, "] must be greater than 0!"
  STOP
END IF

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize the maximum number of repetitions for the simulated annealing
READ( 10, * ) Dummy, MaxRepetitions
IF( MaxRepetitions <= 0 ) THEN
  WRITE( *, "(3G0)" ) "The maximum number of repetitions [", MaxRepetitions, "] must be greater than 0!"
  STOP
END IF

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize the (optimum) layer spacing for the smectic phase 
READ( 10, * ) Dummy, InitialLayerSpacing
IF( InitialLayerSpacing <= 0.D0 ) THEN
  WRITE( *, "(3G0)" ) "The initial layer spacing [", InitialLayerSpacing, "] must be greater than 0!"
  STOP
END IF

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize the maximum layer spacing for the smectic phase
READ( 10, * ) Dummy, MaxLayerSpacing
IF( MaxLayerSpacing <= 0.D0 ) THEN
  WRITE( *, "(3G0)" ) "The maximum layer spacing [", MaxLayerSpacing, "] must be greater than 0!"
  STOP
END IF

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize the minimum layer spacing for the smectic phase
READ( 10, * ) Dummy, MinLayerSpacing
IF( MinLayerSpacing <= 0.D0 ) THEN
  WRITE( *, "(3G0)" ) "The minimum layer spacing [", MinLayerSpacing, "] must be greater than 0!"
  STOP
END IF

! Check if the minimum layer spacing is less than the maximum layer spacing
IF( MinLayerSpacing >= MaxLayerSpacing ) THEN
  WRITE( *, "(5G0)" ) "The minimum layer spacing [", MinLayerSpacing, "] must be less than the maximum layer spacing [", &
  &                   MaxLayerSpacing, "]!"
  STOP
END IF

! Check if the initial layer spacing is within the limits
IF( InitialLayerSpacing < MinLayerSpacing .OR. InitialLayerSpacing > MaxLayerSpacing ) THEN
  WRITE( *, "(7G0)" ) "The initial layer spacing [", InitialLayerSpacing, "] must be between the minimum layer spacing [", &
  &                   MinLayerSpacing, "] and the maximum layer spacing [", MaxLayerSpacing, "]!"
  STOP
END IF

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize the number of cycles for the global simplex optimization
READ( 10, * ) Dummy, cGlobalSimplex
IF( cGlobalSimplex <= 0 ) THEN
  WRITE( *, "(3G0)" ) "The number of cycles for the simplex optimization algorithm [", cGlobalSimplex, "] must be greater than 0!"
  STOP
END IF

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize the option for the calculation of neighbors in 3D (only for Steinhardt order parameters)
READ( 10, * ) Dummy, NeighborOption
IF( NeighborOption < 1 .OR. NeighborOption > 3 ) THEN
  WRITE( *, "(3G0)" ) "The neighbor option [", NeighborOption, "] must be between 1 and 3!"
  STOP
END IF

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize the number of neighbors for the Steinhardt order parameters [option 1]
READ( 10, * ) Dummy, nNeighbors
IF( nNeighbors <= 0 ) THEN
  WRITE( *, "(3G0)" ) "The number of neighbors [", nNeighbors, "] must be greater than 0!"
  STOP
END IF
IF( nNeighbors > nParticles ) THEN
  WRITE( *, "(5G0)" ) "The number of neighbors [", nNeighbors, "] must be less than the number of particles [", nParticles, "]!"
  STOP
END IF

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize the cutoff distance for the Steinhardt order parameters [option 2]
READ( 10, * ) Dummy, CutoffDistance
IF( CutoffDistance <= 0.D0 ) THEN
  WRITE( *, "(3G0)" ) "The cutoff distance [", CutoffDistance, "] must be a positive value!"
  STOP
END IF

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize the maximum number of vertices for the Voronoi construction [option 3]
READ( 10, * ) Dummy, MaxVertices
IF( MaxVertices < 4 ) THEN
  WRITE( *, "(3G0)" ) "The maximum number of vertices [", MaxVertices, "] must be greater than 3!"
  STOP
END IF

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize the cutoff distance for the Voronoi construction [option 3]
READ( 10, * ) Dummy, CutoffDistanceVoronoi
IF( CutoffDistanceVoronoi <= 0.D0 ) THEN
  WRITE( *, "(3G0)" ) "The cutoff distance for the Voronoi construction [", CutoffDistanceVoronoi, "] must be a positive value!"
  STOP
END IF

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize the logical variable that computes the averaged Steinhardt order parameters (q̄ℓ)
READ( 10, * ) Dummy, AveragedSteinhardt

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize the minimum value of the Steinhardt order parameter (azimuthal quantum number, ℓ)
READ( 10, * ) Dummy, MinimumSteinhardt
IF( MinimumSteinhardt < 1 .OR. MinimumSteinhardt > 12 ) THEN
  WRITE( *, "(3G0)" ) "The minimum value of the Steinhardt order parameter [", MinimumSteinhardt, "] must be between 1 and 12!"
  STOP
END IF

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize the maximum value of the Steinhardt order parameter (azimuthal quantum number, ℓ)
READ( 10, * ) Dummy, MaximumSteinhardt

! Check if the maximum value of the Steinhardt order parameter is greater than the minimum value
IF( MaximumSteinhardt < MinimumSteinhardt ) THEN
  WRITE( *, "(5G0)" ) "The maximum value of the Steinhardt order parameter [", MaximumSteinhardt, &
  &                   "] must be greater than the minimum value [", MinimumSteinhardt, "]!"
  STOP
END IF

! Check if the maximum value of the Steinhardt order parameter is within the limits
IF( MinimumSteinhardt < 1 .OR. MaximumSteinhardt > 12 ) THEN
  WRITE( *, "(5G0)" ) "The values of the Steinhardt order parameter [", MinimumSteinhardt, " - ", MaximumSteinhardt, &
  &                   "] must be between 1 and 12!"
  STOP
END IF

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Initialize the logical variable that determines whether the Wigner's parameters are computed (wℓ or w̄ℓ depending whether the averaged Steinhardt order parameters are computed)
READ( 10, * ) Dummy, ComputeWignerParameter

! Close unit
CLOSE( 10 )

RETURN

END SUBROUTINE Initialization_Parameters

END MODULE Initialization