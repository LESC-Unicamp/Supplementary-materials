MODULE InitializeVariables

! Uses two modules: global variables and initial configuration
USE GlobalVariables
USE InitialConfiguration, ONLY: CrossSectionArea

IMPLICIT NONE

CONTAINS

! *********************************************************************************************** !
!                                Initialization of common variables                               !
! *********************************************************************************************** !
SUBROUTINE Common_Variables(  )

IMPLICIT NONE

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 ) :: pSurfaceArea ! Surface area of the considered molecular geometry

! System properties
OPEN( Unit= 10, File= "ini_system.ini" )

! Number of particles
READ( 10, * ) Dummy, nParticles
IF( nParticles <= 1 ) THEN
  WRITE( *, "(3G0)" ) "Error: The number of particles [", nParticles, "] is invalid."
  STOP
END IF
! Reduced number density
READ( 10, * ) Dummy, Density
IF( Density <= 0.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "Error: The reduced number density [", Density, "] is invalid."
  STOP
END IF
IF( .NOT. ConfigurationLogical(2) .OR. Path == 1 ) THEN
  ! Diameter
  READ( 10, * ) Dummy, cDiameter ! Å
  IF( cDiameter <= 0.D0 ) THEN
    WRITE( *, "(G0,G0.5,G0)" ) "Error: The diameter [", cDiameter, "] is invalid."
    STOP
  END IF
  ! Length
  READ( 10, * ) Dummy, cLength   ! Å
  IF( cLength <= 0.D0 ) THEN
    WRITE( *, "(G0,G0.5,G0)" ) "Error: The length [", cLength, "] is invalid."
    STOP
  END IF
ELSE
  READ( 10, * ) Dummy, Dummy
  READ( 10, * ) Dummy, Dummy
END IF
! Aspect ratio (cylinder)
cAspectRatio = cLength / cDiameter
! Temperature
READ( 10, * ) Dummy, Temperature ! K
! Reduced pressure
READ( 10, * ) Dummy, Pressure ! P* = Pσ₀³/(kT), σ₀ = 1Å
IF( Pressure <= 0.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "Error: The reduced pressure [", Pressure, "] is invalid."
  STOP
END IF
! Cross-sectional area
CALL CrossSectionArea( pSurfaceArea )
! Molecular volume
pVolume = pSurfaceArea * cLength ! Å³
! Total molecular volume
TotalParticleVolume = nParticles * pVolume ! Å³

! Cutoff diameter (non-convex body)
IF( GeometrySelection(1) ) THEN
  CutoffSphere = ( 0.5D0 * DSQRT( 2.D0 ) * cDiameter + cDiameter + cLength )
  SquaredCutoffSphere = CutoffSphere * CutoffSphere
  CutoffSpherocylinder = ( 0.5D0 * DSQRT( 2.D0 ) * cDiameter + cDiameter )
  SquaredCutoffSpherocylinder = CutoffSpherocylinder * CutoffSpherocylinder
ELSE IF( GeometrySelection(2) ) THEN
  CutoffSphere = ( 2.D0 * cDiameter + cLength )
  SquaredCutoffSphere = CutoffSphere * CutoffSphere
  CutoffSpherocylinder = 2.D0 * cDiameter
  SquaredCutoffSpherocylinder = CutoffSpherocylinder * CutoffSpherocylinder
END IF

! Cutoff diameter (cylinders)
cCutoffSphere = ( cDiameter + cLength )
cSquaredCutoffSphere = cCutoffSphere * cCutoffSphere

CLOSE( 10 )

! Initial configuration properties
OPEN( Unit= 10, File= "ini_initial_configuration.ini" )

! Skip
READ( 10 , * ) Dummy, Dummy
READ( 10 , * ) Dummy, Dummy
! Quaternion angle
READ( 10 , * ) Dummy, QuaternionAngle
! Maximum number of cycles (random configuration)
READ( 10 , * ) Dummy, MaxSimulationCyclesInit
IF( MaxSimulationCyclesInit < 0 ) THEN
  WRITE( *, "(4G0)" ) "Error: The maximum number of cycles [", MaxSimulationCyclesInit, "] for the initial configuration ", &
  &                   "is invalid."
  STOP
END IF
! Restore simulation from backup
READ( 10 , * ) Dummy, RestoreBackup
! Box volume factor (packed box)
READ( 10 , * ) Dummy, BoxFactor
IF( BoxFactor <= 0.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "Error: The box factor [", BoxFactor, "] is invalid."
  STOP
END IF

CLOSE( 10 )

! Output file descriptors (Aspect Ratio [1])
FormatAR = "(F0.5)"
WRITE( DescriptorAR, FormatAR ) cAspectRatio
! Output file descriptors (Pressure [2])
FormatPressure = "(F0.5)"
WRITE( DescriptorPressure, FormatPressure ) Pressure

RETURN

END SUBROUTINE Common_Variables

! *********************************************************************************************** !
!                            Initialization of Monte Carlo parameters                             !
! *********************************************************************************************** !
SUBROUTINE MonteCarlo_Variables(  )

IMPLICIT NONE

! Simulation parameters
OPEN( Unit= 10, File= "ini_montecarlo.ini" )

! Total number of cycles
READ( 10 , * ) Dummy, nSimulationCycles
IF( nSimulationCycles <= 0 ) THEN
  WRITE( *, "(3G0)" ) "Error: The total number of cycles [", nSimulationCycles, "] is invalid."
  STOP
END IF
! Number of equilibration cycles
READ( 10 , * ) Dummy, nEquilibration
IF( nEquilibration < 0 ) THEN
  WRITE( *, "(3G0)" ) "Error: The number of equilibration cycles [", nEquilibration, "] is invalid."
  STOP
END IF
IF(  nEquilibration >= nSimulationCycles ) THEN
  WRITE( *, "(6G0)" ) "Error: The number of equilibration cycles [", nEquilibration, "] is greater than the total number ", &
  &                   "of cycles [", nSimulationCycles, "]."
  STOP
END IF
! Saving frequency
READ( 10 , * ) Dummy, nSavingFrequency
IF( nSavingFrequency <= 0 ) THEN
  WRITE( *, "(3G0)" ) "Error: The saving frequency [", nSavingFrequency, "] is invalid."
  STOP
END IF
! Adjustment frequency (molecular movement)
READ( 10 , * ) Dummy, nAdjustmentMovementFrequency
IF( nAdjustmentMovementFrequency <= 0 ) THEN
  WRITE( *, "(3G0)" ) "Error: The adjustment frequency [", nAdjustmentMovementFrequency, "] is invalid."
  STOP
END IF
! Adjustment frequency (volume change)
READ( 10 , * ) Dummy, nAdjustmentVolumeFrequency
IF( nAdjustmentVolumeFrequency <= 0 ) THEN
  WRITE( *, "(3G0)" ) "Error: The adjustment frequency [", nAdjustmentVolumeFrequency, "] is invalid."
  STOP
END IF
! Maximum translational displacement
READ( 10 , * ) Dummy, UserMaxTranslationalDisplacement
! Maximum rotational displacement
READ( 10 , * ) Dummy, UserMaxRotationalDisplacement
! Maximum volume change (isotropic)
READ( 10 , * ) Dummy, UserMaxIsoVolumetricDisplacement
! Maximum volume change (anisotropic)
READ( 10 , * ) Dummy, UserMaxAnisoVolumetricDisplacement
! Maximum box distortion
READ( 10 , * ) Dummy, MaxBoxDistortion
IF( MaxBoxDistortion <= 1.D0 ) THEN
  WRITE( *, "(3G0)" ) "Error: The maximum box distortion [", MaxBoxDistortion, "] is invalid."
  STOP
END IF
! Maximum length ratio (box distortion)
READ( 10, * ) Dummy, MaxLengthRatio
IF( MaxLengthRatio <= 0.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "Error: The maximum length ratio [", MaxLengthRatio, "] is invalid."
  STOP
END IF
! Maximum angle (box distortion)
READ( 10, * ) Dummy, MaxBoxAngle
IF( MaxBoxAngle <= 0.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "Error: The maximum box angle [", MaxBoxAngle, "] is invalid."
  STOP
END IF
MaxBoxAngle = MaxBoxAngle * cPi / 180.D0
! Lattice reduction method
READ( 10 , * ) Dummy, LatticeReductionMethod
CALL ToUpper( LatticeReductionMethod, LEN_TRIM( LatticeReductionMethod ), LatticeReductionMethod )
LatticeReductionTypeLogical(:) = .FALSE.
! Lattice reduction: Gottwald method
IF( LatticeReductionMethod == "FBM" ) THEN
  LatticeReductionTypeLogical(1) = .TRUE.
! lattice reduction: Lenstra-Lenstra-Lovász method
ELSE IF( LatticeReductionMethod == "LLL" ) THEN
  LatticeReductionTypeLogical(2) = .TRUE.
END IF

CLOSE( 10 )

! Potential parameters
OPEN( Unit= 10, File= "ini_orientational_field.ini" )

! Orientation field
READ( 10 , * ) Dummy, PotentialLogical
! Strength of the orientational field
READ( 10 , * ) Dummy, PotentialStrength
! Orientational field vector
READ( 10 , * ) Dummy, OrientationalFieldVector
! Normalization of the orientational field
OrientationalFieldVector = OrientationalFieldVector / DSQRT( DOT_PRODUCT( OrientationalFieldVector, OrientationalFieldVector ) )

CLOSE( 10 )

! Simulation parameters
OPEN( Unit= 10, File= "ini_probabilities.ini" )

! Displacement probability
READ( 10, * ) Dummy, MovementProbability
IF( MovementProbability < 0.D0 ) MovementProbability = 0.D0
IF( MovementProbability > 1.D0 ) MovementProbability = 1.D0
! Volume change probability
VolumeChangeProbability = 1.D0 - MovementProbability
! Translational probability
READ( 10, * ) Dummy, TranslationalProbability
IF( TranslationalProbability < 0.D0 ) TranslationalProbability = 0.D0
IF( TranslationalProbability > 1.D0 ) TranslationalProbability = 1.D0
! Rotational probability
RotationalProbability = 1.D0 - TranslationalProbability
! Isotropic volume change probability
READ( 10, * ) Dummy, IsoVolumetricProbability
IF( IsoVolumetricProbability < 0.D0 ) IsoVolumetricProbability = 0.D0
IF( IsoVolumetricProbability > 1.D0 ) IsoVolumetricProbability = 1.D0
! Anisotropic volume change probability
AnisoVolumetricProbability = 1.D0 - IsoVolumetricProbability
! Stack rotation probability
READ( 10, * ) Dummy, StackRotationProbability
IF( StackRotationProbability < 0.D0 ) StackRotationProbability = 0.D0
IF( StackRotationProbability > 1.D0 ) StackRotationProbability = 1.D0

CLOSE( 10 )

! Simulation parameters
OPEN( Unit= 10, File= "ini_ratios.ini" )

! Translational acceptance ratio threshold
READ( 10, * ) Dummy, AcceptanceRatioTranslation
IF( AcceptanceRatioTranslation < 0.D0 ) AcceptanceRatioTranslation = 0.D0
IF( AcceptanceRatioTranslation > 1.D0 ) AcceptanceRatioTranslation = 1.D0
! Rotational acceptance ratio threshold
READ( 10, * ) Dummy, AcceptanceRatioRotation
IF( AcceptanceRatioRotation < 0.D0 ) AcceptanceRatioRotation = 0.D0
IF( AcceptanceRatioRotation > 1.D0 ) AcceptanceRatioRotation = 1.D0
! Volumetric acceptance ratio threshold (isotropic)
READ( 10, * ) Dummy, AcceptanceRatioIsoVolumeChange
IF( AcceptanceRatioIsoVolumeChange < 0.D0 ) AcceptanceRatioIsoVolumeChange = 0.D0
IF( AcceptanceRatioIsoVolumeChange > 1.D0 ) AcceptanceRatioIsoVolumeChange = 1.D0
! Volumetric acceptance ratio threshold (anisotropic)
READ( 10, * ) Dummy, AcceptanceRatioAnisoVolumeChange
IF( AcceptanceRatioAnisoVolumeChange < 0.D0 ) AcceptanceRatioAnisoVolumeChange = 0.D0
IF( AcceptanceRatioAnisoVolumeChange > 1.D0 ) AcceptanceRatioAnisoVolumeChange = 1.D0

CLOSE( 10 )

RETURN

END SUBROUTINE MonteCarlo_Variables

! *********************************************************************************************** !
!                           Initialization of Inquiry/Control variables                           !
! *********************************************************************************************** !
SUBROUTINE Control_Variables(  )

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 ) :: FrameLeft  ! Box frame dimension
INTEGER( Kind= Int64 ) :: FrameRight ! Box frame dimension

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 ) :: pSurfaceArea ! Surface area of the considered molecular geometry

! Cross-sectional area
CALL CrossSectionArea( pSurfaceArea )

! Control parameters
OPEN( Unit= 10, File= "ini_control.ini" )

! Trajectory inquiry
READ( 10 , * ) Dummy, PrintTrajectory
! Seed type
READ( 10, * ) Dummy, RandomSeed
! Widom insertion method
READ( 10, * ) Dummy, WidomInsertion
! Number of cycles per configuration (Widom insertion)
READ( 10 , * ) Dummy, nWidomCycles
IF( WidomInsertion .AND. nWidomCycles <= 0 ) THEN
  WRITE( *, "(3G0)" ) "Error: The number of Widom cycles [", nWidomCycles, "] is invalid."
  STOP
END IF
! Stack rotation moves
READ( 10 , * ) Dummy, StackRotationLogical
! Density after which stack rotations are performed
READ( 10 , * ) Dummy, StackRotationDensityThreshold
IF( StackRotationLogical .AND. StackRotationDensityThreshold < 0.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "Error: The stack rotation density threshold [", StackRotationDensityThreshold, "] is invalid."
  STOP
END IF
! Cell list inquiry
READ( 10, * ) Dummy, CellListLogical
IF( .NOT. CellListLogical ) CellListControl = .FALSE.

CLOSE( 10 )

! Summary
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 70 )//CH_UR
FrameLeft  = (35 - FLOOR( REAL( 17 ) / 2.D0 ) )
FrameRight = (35 - CEILING( REAL( 17 ) / 2.D0 ) )
WRITE( *, "(2G0)" ) CH_VS//REPEAT( " ", FrameLeft )//"SYSTEM PROPERTIES"//REPEAT( " ", FrameRight )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 70 )//CH_BR
WRITE( *, "(2G0)" ) "Number of particles: ", nParticles
WRITE( *, "(2G0)" ) "Number of cylinders: ", 4 * nParticles
WRITE( *, "(G0,G0.5)" ) "Reduced number density: ", Density
WRITE( *, "(G0,G0.5,G0)" ) "Diameter (cylinder): ", cDiameter, "Å"
WRITE( *, "(G0,G0.5,G0)" ) "Length (cylinder): ", cLength, "Å"
WRITE( *, "(G0,G0.5)" ) "Aspect ratio (cylinder): ", cAspectRatio
WRITE( *, "(G0,G0.5,G0)" ) "Temperature: ", Temperature, "K"
WRITE( *, "(G0,G0.5)" ) "Reduced pressure: ", Pressure
WRITE( *, "(G0,G0.5,G0)" ) "Real pressure: ", Pressure * cBoltzmann * Temperature / 1.D-30 / 1.D6, "MPa"
WRITE( *, "(G0,G0.5,G0)" ) "Surface area: ", pSurfaceArea, "Å²"
WRITE( *, "(G0,G0.5,G0)" ) "Particle volume: ", pVolume, "Å³"
WRITE( *, "(G0)" ) " "

! Summary
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 70 )//CH_UR
FrameLeft  = (35 - FLOOR( REAL( 32 ) / 2.D0 ) )
FrameRight = (35 - CEILING( REAL( 32 ) / 2.D0 ) )
WRITE( *, "(2G0)" ) CH_VS//REPEAT( " ", FrameLeft )//"INITIAL CONFIGURATION PROPERTIES"//REPEAT( " ", FrameRight )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 70 )//CH_BR
IF( ConfigurationLogical(1) ) WRITE( *, "(G0)" ) "Initial Configuration: Packed Box"
IF( ConfigurationLogical(2) ) WRITE( *, "(G0)" ) "Initial Configuration: Floppy Box"
WRITE( *, "(2G0)" ) "Arrangement Type: ", pGeometry
IF( .NOT. ConfigurationLogical(2) ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "Quaternion angle: ", QuaternionAngle, "°"
  WRITE( *, "(G0)" ) "Unrotated reference axis: [0, 0, 1]"
  WRITE( *, "(2G0)" ) "Maximum number of cycles (random configuration): ", MaxSimulationCyclesInit
  WRITE( *, "(G0,G0.5)" ) "Box stretching/contraction factor: ", BoxFactor
END IF
WRITE( *, "(G0)" ) " "

! Summary
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 70 )//CH_UR
FrameLeft  = (35 - FLOOR( REAL( 21 ) / 2.D0 ) )
FrameRight = (35 - CEILING( REAL( 21 ) / 2.D0 ) )
WRITE( *, "(2G0)" ) CH_VS//REPEAT( " ", FrameLeft )//"SIMULATION PROPERTIES"//REPEAT( " ", FrameRight )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 70 )//CH_BR
WRITE( *, "(2G0)" ) "Number of Cycles: ", nSimulationCycles
WRITE( *, "(2G0)" ) "Number of Equilibration Cycles: ", nEquilibration
WRITE( *, "(2G0)" ) "Number of Production Cycles: ", nSimulationCycles - nEquilibration
WRITE( *, "(3G0)" ) "Saving Frequency: Every ", nSavingFrequency, " Cycle(s)"
WRITE( *, "(3G0)" ) "Adjustment Frequency (Movement): Every ", nAdjustmentMovementFrequency, " Cycle(s)"
WRITE( *, "(3G0)" ) "Adjustment Frequency (Volume): Every ", nAdjustmentVolumeFrequency, " Cycle(s)"
WRITE( *, "(G0,G0.5,G0)" ) "Maximum Translational Displacement: ", UserMaxTranslationalDisplacement, "Å"
WRITE( *, "(G0,G0.5,G0)" ) "Maximum Rotational Displacement: ", UserMaxRotationalDisplacement, "rad"
WRITE( *, "(G0,G0.5,G0)" ) "Maximum Volumetric Displacement (Isotropic): ", UserMaxIsoVolumetricDisplacement, "Å"
WRITE( *, "(G0,G0.5,G0)" ) "Maximum Volumetric Displacement (Anisotropic): ", UserMaxAnisoVolumetricDisplacement, "Å"
WRITE( *, "(G0,G0.5)" ) "Translational Acceptance Ratio Threshold: ", AcceptanceRatioTranslation
WRITE( *, "(G0,G0.5)" ) "Rotational Acceptance Ratio Threshold: ", AcceptanceRatioRotation
WRITE( *, "(G0,G0.5)" ) "Volumetric Acceptance Ratio Threshold (Isotropic): ", AcceptanceRatioIsoVolumeChange
WRITE( *, "(G0,G0.5)" ) "Volumetric Acceptance Ratio Threshold (Anisotropic): ", AcceptanceRatioAnisoVolumeChange
WRITE( *, "(G0,G0.5,G0)" ) "Probability of Molecular Movement: ", MovementProbability * 100.D0, "%"
WRITE( *, "(G0,G0.5,G0)" ) "Probability of Translational Displacement: ", TranslationalProbability * 100.D0, "%"
WRITE( *, "(G0,G0.5,G0)" ) "Probability of Rotational Displacement: ", RotationalProbability * 100.D0, "%"
WRITE( *, "(G0,G0.5,G0)" ) "Probability of Volume Change: ", VolumeChangeProbability * 100.D0, "%"
WRITE( *, "(G0,G0.5,G0)" ) "Probability of Isotropic Volume Change: ", IsoVolumetricProbability * 100.D0, "%"
WRITE( *, "(G0,G0.5,G0)" ) "Probability of Anisotropic Volume Change: ", AnisoVolumetricProbability * 100.D0, "%"
WRITE( *, "(G0,G0.5)" ) "Maximum Box Distortion Factor: ", MaxBoxDistortion
IF( LatticeReductionTypeLogical(1) ) THEN
  WRITE( *, "(G0)" ) "Lattice Reduction Method: Gottwald"
ELSE IF( LatticeReductionTypeLogical(2) ) THEN
  WRITE( *, "(G0)" ) "Lattice Reduction Method: Lenstra-Lenstra-Lovász"
END IF
IF( WidomInsertion ) THEN
  WRITE( *, "(G0,G0)" ) "Maximum Number of Widom Cycles per Molecular Configuration: ", nWidomCycles
END IF
WRITE( *, "(G0)" ) " "

! Summary
WRITE( *, "(G0)" ) CH_UL//REPEAT( CH_HS, 70 )//CH_UR
FrameLeft  = (35 - FLOOR( REAL( 18 ) / 2.D0 ) )
FrameRight = (35 - CEILING( REAL( 18 ) / 2.D0 ) )
WRITE( *, "(2G0)" ) CH_VS//REPEAT( " ", FrameLeft )//"CONTROL PROPERTIES"//REPEAT( " ", FrameRight )//CH_VS
WRITE( *, "(G0)" ) CH_BL//REPEAT( CH_HS, 70 )//CH_BR
IF( PrintTrajectory ) WRITE( *, "(G0)" ) "Print Trajectory? [YES]"
IF( .NOT. PrintTrajectory ) WRITE( *, "(G0)" ) "Print Trajectory? [NO]"
IF( CellListLogical ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "Cell List: Enabled"
ELSE
  WRITE( *, "(G0,G0.5,G0)" ) "Cell List: Disabled"
END IF
IF( RandomSeed ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "Seed type: Random"
ELSE
  WRITE( *, "(G0,G0.5,G0)" ) "Seed type: Fixed"  
END IF

WRITE( *, "(G0)" ) " "
WRITE( *, "(2G0)" ) "Number of Threads: ", nThreads
WRITE( *, "(G0)" ) " "

RETURN

END SUBROUTINE Control_Variables

END MODULE InitializeVariables