! ----------------------------------------------------------------------- !
! Variable Initialization Subroutine                                      !
! ----------------------------------------------------------------------- !
SUBROUTINE Variable_Initialization(  )

USE GlobalVariables
USE VectorOperations, ONLY: InverseMatrixCofactorVec

IMPLICIT NONE

! Check for 'parameters.ini' file existence
INQUIRE( File= "parameters.ini", EXIST= FileExists )
IF( .NOT. FileExists ) THEN
  WRITE( *, "(A)" ) "Error: 'parameters.ini' file not found."
  STOP
END IF

! Read parameters from 'parameters.ini'
OPEN( Unit= 10, File= "parameters.ini", Action= 'READ' )

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Number of particles
READ( 10, * ) Dummy, nParticles
IF( nParticles <= 0 ) THEN
  WRITE( *, "(A)" ) "Error: Number of particles must be greater than zero."
  STOP
END IF

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Density
READ( 10, * ) Dummy, Density ! 1/Å³
IF( Density <= 0.D0 ) THEN
  WRITE( *, "(A)" ) "Error: Density must be greater than zero."
  STOP
END IF

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Temperature
READ( 10, * ) Dummy, Temperature ! K
IF( Temperature <= 0.D0 ) THEN
  WRITE( *, "(A)" ) "Error: Temperature must be greater than zero."
  STOP
END IF

! De Broglie thermal wavelength
deBroglieThermalWavelength = SQRT( ( cPlanck * cPlanck ) / ( 2.D0 * cPi * waterMass * cBoltzmann * Temperature ) ) * 1.D10 ! Å
sqBroglieWavelength = deBroglieThermalWavelength * deBroglieThermalWavelength ! Å²
cbBroglieWavelength = deBroglieThermalWavelength * deBroglieThermalWavelength * deBroglieThermalWavelength ! Å³

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Reduced chemical potential
READ( 10, * ) Dummy, ReducedChemicalPotential
READ( 10, * ) Dummy, excessReducedChemPot

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Pore dimensions
READ( 10, * ) Dummy, poreWidth ! Å
IF( poreWidth <= 0.D0 ) THEN
  WRITE( *, "(A)" ) "Error: Pore width must be greater than zero."
  STOP
END IF
READ( 10, * ) Dummy, poreHeight ! Å
IF( poreHeight <= 0.D0 ) THEN
  WRITE( *, "(A)" ) "Error: Pore height must be greater than zero."
  STOP
END IF
! Pore length [Å]
poreLength = nParticles / Density / poreHeight / poreWidth ! Å

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Solid properties
READ( 10, * ) Dummy, SolidDensity ! 1/Å³
IF( SolidDensity <= 0.D0 ) THEN
  WRITE( *, "(A)" ) "Error: Solid density must be greater than zero."
  STOP
END IF
READ( 10, * ) Dummy, InterlayerSpacing ! Å
IF( InterlayerSpacing <= 0.D0 ) THEN
  WRITE( *, "(A)" ) "Error: Interlayer spacing must be greater than zero."
  STOP
END IF
READ( 10, * ) Dummy, SigmaSolid ! Å
IF( SigmaSolid <= 0.D0 ) THEN
  WRITE( *, "(A)" ) "Error: Solid diameter must be greater than zero."
  STOP
END IF
READ( 10, * ) Dummy, EpsilonSolid ! K
IF( EpsilonSolid < 0.D0 ) THEN
  WRITE( *, "(A)" ) "Error: Solid well depth must be non-negative."
  STOP
END IF
READ( 10, * ) Dummy, isSolidPotentialEnabled

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Fluid properties (Mie/SW)
READ( 10, * ) Dummy, EpsilonFluid ! K
IF( EpsilonFluid < 0.D0 ) THEN
  WRITE( *, "(A)" ) "Error: Fluid well depth must be non-negative."
  STOP
END IF
READ( 10, * ) Dummy, SigmaFluid ! Å
IF( SigmaFluid <= 0.D0 ) THEN
  WRITE( *, "(A)" ) "Error: Fluid diameter must be greater than zero."
  STOP
END IF
READ( 10, * ) Dummy, cutoffRadiusMie ! Å
IF( cutoffRadiusMie <= 0.D0 ) THEN
  WRITE( *, "(A)" ) "Error: Cutoff radius (Mie) must be greater than zero."
  STOP
END IF
READ( 10, * ) Dummy, attractiveRange
IF( attractiveRange <= 0.D0 ) THEN
  WRITE( *, "(A)" ) "Error: Attractive exponent (Mie) must be greater than zero."
  STOP
END IF
READ( 10, * ) Dummy, repulsiveRange
IF( repulsiveRange <= attractiveRange ) THEN
  WRITE( *, "(A)" ) "Error: Repulsive exponent (Mie) must be greater than attractive exponent."
  STOP
END IF
READ( 10, * ) Dummy, SWRange
IF( SWRange <= 1.D0 ) THEN
  WRITE( *, "(A)" ) "Error: Square-well range must be greater than one."
  STOP
END IF
READ( 10, * ) Dummy, isMiePotentialEnabled
READ( 10, * ) Dummy, isHardSpherePotentialEnabled
READ( 10, * ) Dummy, isSWPotentialEnabled

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Fluid properties (Water)
READ( 10, * ) Dummy, EpsilonWater ! K
IF( EpsilonWater < 0.D0 ) THEN
  WRITE( *, "(A)" ) "Error: Association well depth must be non-negative."
  STOP
END IF
READ( 10, * ) Dummy, cutoffRadiusWater
IF( cutoffRadiusWater <= 0.D0 ) THEN
  WRITE( *, "(A)" ) "Error: Cutoff radius (association) must be greater than zero."
  STOP
END IF
cutoffRadiusWater = cutoffRadiusWater * SigmaFluid ! Å
READ( 10, * ) Dummy, SiteRadius
IF( SiteRadius <= 0.D0 ) THEN
  WRITE( *, "(A)" ) "Error: Site radius must be greater than zero."
  STOP
END IF
SiteRadius = SiteRadius * SigmaFluid ! Å
READ( 10, * ) Dummy, AngleWater ! °
IF( AngleWater <= 0.D0 .OR. AngleWater >= 180.D0 ) THEN
  WRITE( *, "(A)" ) "Error: Bond angle must be in the range (0, 180) degrees."
  STOP
END IF
AngleWater = AngleWater * ( cPi / 180.D0 ) ! Radians
READ( 10, * ) Dummy, isAssociativePotentialEnabled

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Combining rules
READ( 10, * ) Dummy, OverwriteCombiningRule
IF( OverwriteCombiningRule ) THEN
  READ( 10, * ) Dummy, SigmaSolid_Fluid   ! Å
  READ( 10, * ) Dummy, EpsilonSolid_Fluid ! K
ELSE
  READ( 10, * ) Dummy, Dummy
  READ( 10, * ) Dummy, Dummy
END IF

! Check if translation along x and y are intertwined to translations in z
READ( 10, * ) Dummy, xyTranslationIndependency
READ( 10, * ) Dummy, TranslationalProbabilityXY
IF( TranslationalProbabilityXY < 0.D0 .OR. TranslationalProbabilityXY > 1.D0 ) THEN
  WRITE( *, "(A)" ) "Error: Isotropic volume change probability must be in the range [0, 1]."
  STOP
END IF

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Monte Carlo simulation parameters
READ( 10, * ) Dummy, nCycles
IF( nCycles <= 0 ) THEN
  WRITE( *, "(A)" ) "Error: Number of MC cycles must be greater than zero."
  STOP
END IF
READ( 10, * ) Dummy, nEquilibration
IF( nEquilibration < 0 .OR. nEquilibration >= nCycles ) THEN
  WRITE( *, "(A)" ) "Error: Number of equilibration cycles must be non-negative and less than the number of MC cycles."
  STOP
END IF
READ( 10, * ) Dummy, nEquilRandom   ! Random configuration only
IF( nEquilRandom < 0 ) THEN
  WRITE( *, "(A)" ) "Error: Number of equilibration cycles (random configuration) must be non-negative."
  STOP
END IF
READ( 10, * ) Dummy, maxAttemptsNPT ! Random configuration only
IF( maxAttemptsNPT <= 0 ) THEN
  WRITE( *, "(A)" ) "Error: Maximum number of attempts in the NPT ensemble (random configuration) must be greater than zero."
  STOP
END IF
READ( 10, * ) Dummy, nAdjustment
IF( nAdjustment <= 0 ) THEN
  WRITE( *, "(A)" ) "Error: Adjustment frequency must be greater than zero."
  STOP
END IF
READ( 10, * ) Dummy, nSave
IF( nSave <= 0 ) THEN
  WRITE( *, "(A)" ) "Error: Saving frequency must be greater than zero."
  STOP
END IF
READ( 10, * ) Dummy, nSaveConfig
IF( nSaveConfig <= 0 ) THEN
  WRITE( *, "(A)" ) "Error: Configuration saving frequency must be greater than zero."
  STOP
END IF
READ( 10, * ) Dummy, maxTranslationalDisplc
IF( maxTranslationalDisplc <= 0.D0 ) THEN
  WRITE( *, "(A)" ) "Error: Maximum translational displacement must be greater than zero."
  STOP
END IF
maxTranslationalDisplc = maxTranslationalDisplc * SigmaFluid ! Å
READ( 10, * ) Dummy, maxTranslationalDisplcRandom
maxTranslationalDisplcRandom = maxTranslationalDisplcRandom * SigmaFluid ! Å
READ( 10, * ) Dummy, maxRotationalDisplc
READ( 10, * ) Dummy, maxRotationalDisplcRandom
READ( 10, * ) Dummy, MaxIsoVolumetricDisplacementRandom
READ( 10, * ) Dummy, MaxAnisoVolumetricDisplacementRandom
READ( 10, * ) Dummy, TargetPressureRandom
IF( TargetPressureRandom <= 0.D0 ) THEN
  WRITE( *, "(A)" ) "Error: Target pressure (random configuration) must be greater than zero."
  STOP
END IF
! Convert from MPa to Pa
TargetPressureRandom = TargetPressureRandom * 1.D6 ! Pa
! Reduced pressure
ReducedTargetPressure = TargetPressureRandom * 1.D-30 / cBoltzmann / Temperature
READ( 10, * ) Dummy, MovementProbability
IF( MovementProbability < 0.D0 .OR. MovementProbability > 1.D0 ) THEN
  WRITE( *, "(A)" ) "Error: Movement probability must be in the range [0, 1]."
  STOP
END IF
VolumeChangeProbability = 1.D0 - MovementProbability
PCDProbability = 1.D0 - MovementProbability
READ( 10, * ) Dummy, TranslationalProbability
IF( TranslationalProbability < 0.D0 .OR. TranslationalProbability > 1.D0 ) THEN
  WRITE( *, "(A)" ) "Error: Translational probability must be in the range [0, 1]."
  STOP
END IF
RotationalProbability = 1.D0 - TranslationalProbability
READ( 10, * ) Dummy, IsotropicProbability
IF( IsotropicProbability < 0.D0 .OR. IsotropicProbability > 1.D0 ) THEN
  WRITE( *, "(A)" ) "Error: Isotropic volume change probability must be in the range [0, 1]."
  STOP
END IF
AnisotropicProbability = 1.D0 - IsotropicProbability
READ( 10, * ) Dummy, CreationProbability
IF( CreationProbability < 0.D0 .OR. CreationProbability > 1.D0 ) THEN
  WRITE( *, "(A)" ) "Error: Creation probability must be in the range [0, 1]."
  STOP
END IF
DestructionProbability = 1.D0 - CreationProbability
READ( 10, * ) Dummy, AcceptanceRatioTranslation
IF( AcceptanceRatioTranslation < 0.D0 .OR. AcceptanceRatioTranslation > 1.D0 ) THEN
  WRITE( *, "(A)" ) "Error: Acceptance ratio for translation must be in the range [0, 1]."
  STOP
END IF
READ( 10, * ) Dummy, AcceptanceRatioRotation
IF( AcceptanceRatioRotation < 0.D0 .OR. AcceptanceRatioRotation > 1.D0 ) THEN
  WRITE( *, "(A)" ) "Error: Acceptance ratio for rotation must be in the range [0, 1]."
  STOP
END IF
READ( 10, * ) Dummy, AcceptanceRatioIsotropic
IF( AcceptanceRatioIsotropic < 0.D0 .OR. AcceptanceRatioIsotropic > 1.D0 ) THEN
  WRITE( *, "(A)" ) "Error: Acceptance ratio for isotropic volume change must be in the range [0, 1]."
  STOP
END IF
READ( 10, * ) Dummy, AcceptanceRatioAnisotropic
IF( AcceptanceRatioAnisotropic < 0.D0 .OR. AcceptanceRatioAnisotropic > 1.D0 ) THEN
  WRITE( *, "(A)" ) "Error: Acceptance ratio for anisotropic volume change must be in the range [0, 1]."
  STOP
END IF
READ( 10, * ) Dummy, AcceptanceRatioCreation
IF( AcceptanceRatioCreation < 0.D0 .OR. AcceptanceRatioCreation > 1.D0 ) THEN
  WRITE( *, "(A)" ) "Error: Acceptance ratio for creation must be in the range [0, 1]."
  STOP
END IF
READ( 10, * ) Dummy, AcceptanceRatioDestruction
IF( AcceptanceRatioDestruction < 0.D0 .OR. AcceptanceRatioDestruction > 1.D0 ) THEN
  WRITE( *, "(A)" ) "Error: Acceptance ratio for destruction must be in the range [0, 1]."
  STOP
END IF
READ( 10, * ) Dummy, BoxEdgeMaxRatio
IF( BoxEdgeMaxRatio <= 1.D0 ) THEN
  WRITE( *, "(A)" ) "Error: Maximum box edge ratio (random configuration) must be greater than one."
  STOP
END IF
READ( 10, * ) Dummy, nWidomFrequency
IF( nWidomFrequency < 0 ) THEN
  WRITE( *, "(A)" ) "Error: Widom insertion frequency must be non-negative."
  STOP
END IF

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! Logical variables
READ( 10, * ) Dummy, RandomConf
READ( 10, * ) Dummy, RandomSeed
READ( 10, * ) Dummy, ApplyNPTRandom
READ( 10, * ) Dummy, Confinement
READ( 10, * ) Dummy, isWidomEnabled

! Skip comments and blank lines
DO
  READ( 10, "(A)" ) LineIndex
  IF( INDEX( LineIndex, "#" ) >= 1 .OR.  TRIM( ADJUSTL( LineIndex ) ) == "" ) CYCLE
  IF( INDEX( LineIndex, "#" ) == 0 .AND. TRIM( ADJUSTL( LineIndex ) ) /= "" ) EXIT
END DO
BACKSPACE( 10 )

! RDF calculation parameters
READ( 10, * ) Dummy, nSlabs
IF( nSlabs <= 0 ) THEN
  WRITE( *, "(A)" ) "Error: Number of slabs must be greater than zero."
  STOP
END IF

CLOSE( Unit= 10 )

! Allocation
ALLOCATE( pPosition(3,nParticles), pQuaternion(0:3,nParticles), pOrientation(3,nParticles), sitePosition(3,4,nParticles) )
ALLOCATE( siteOrientation(3,4,nParticles), siteQuaternion(0:3,4,nParticles) )
ALLOCATE( pPositionRandom(3,nParticles), pQuaternionRandom(0:3,nParticles), pOrientationRandom(3,nParticles) )
ALLOCATE( sitePositionRandom(3,4,nParticles), siteOrientationRandom(3,4,nParticles), siteQuaternionRandom(0:3,4,nParticles) )
ALLOCATE( AssociatedSites( 2, 4, nParticles ), isSiteAssociated( 4, nParticles ) )

! Combined rules
IF( .NOT. OverwriteCombiningRule ) THEN
  SigmaSolid_Fluid   = 0.5D0 * ( SigmaSolid + SigmaFluid )
  EpsilonSolid_Fluid = SQRT( EpsilonSolid * EpsilonFluid )
END IF

! Associate box lengths [Å]
BoxLength = 0.D0
BoxLength(1) = poreHeight
BoxLength(5) = poreLength
BoxLength(9) = poreWidth

! Inverse diameter [1/Å]
invSigmaFluid = 1.D0 / SigmaFluid

! Number of cells
nCells(1) = FLOOR( poreHeight * invSigmaFluid)
nCells(2) = FLOOR( poreLength * invSigmaFluid) 
nCells(3) = FLOOR( poreWidth * invSigmaFluid)
IF( ANY( nCells == 0 ) ) THEN
  WRITE( *, "(A)" ) "Error: Number of cells in at least one dimension is zero. Increase pore size or decrease particle size."
  STOP
END IF

! Cell lengths
cellLength(1) = poreHeight / nCells(1)
cellLength(2) = poreLength / nCells(2)
cellLength(3) = poreWidth / nCells(3)

! Volume calculation
CALL InverseMatrixCofactorVec( BoxLength, BoxLengthInverse, BoxVolume )

! Grand-canonical constant
GrandCanonicalConstant = LOG( BoxVolume / cbBroglieWavelength )

! Well depth calculation
reducedEpsilonFluid = EpsilonFluid / Temperature
reducedEpsilonWater = EpsilonWater / Temperature

! Diameter calculations
squaredSigmaFluid = SigmaFluid * SigmaFluid ! Å²
cubicSigmaFluid = SigmaFluid * SigmaFluid * SigmaFluid ! Å³

! Cutoff radius calculations
cutoffRadiusAssociation = 2.D0 * SiteRadius + cutoffRadiusWater ! Å
cutoffSquaredRadiusAssociation = cutoffRadiusAssociation * cutoffRadiusAssociation ! Å²
cutoffSquaredRadiusWater = cutoffRadiusWater * cutoffRadiusWater ! Å²
cutoffSquaredRadiusMie = cutoffRadiusMie * cutoffRadiusMie ! Å²
inverseCutoffRadiusMie = SigmaFluid / cutoffRadiusMie
inverseCubicCutoffRadiusMie = inverseCutoffRadiusMie * inverseCutoffRadiusMie * inverseCutoffRadiusMie
InverseNinthCutoffRadiusMie = inverseCubicCutoffRadiusMie * inverseCubicCutoffRadiusMie * inverseCubicCutoffRadiusMie

! Potential range calculations
ratioRange = repulsiveRange / attractiveRange
ratioRepulsiveRange = repulsiveRange / ( repulsiveRange - attractiveRange )
ratioAttractiveRange = attractiveRange / ( repulsiveRange - attractiveRange )
NormalisationCoefficient = ratioRepulsiveRange * ratioRange ** ( ratioAttractiveRange )
sqSWRange = SWRange * SWRange * squaredSigmaFluid

! Long-range corrections
longrangeCorrectionPotential = 8.D-30 * cPi * Density * ( SigmaFluid ** 3.D0 ) * reducedEpsilonFluid * &
&                              ( InverseNinthCutoffRadiusMie / 3.D0 - inverseCubicCutoffRadiusMie ) / 3.D0
longrangeCorrectionPressure  = 16.D-30 * cPi * Density * Density * SigmaFluid ** 3.D0 * reducedEpsilonFluid * ( 2.D0 * &
&                              InverseNinthCutoffRadiusMie / 3.D0 - inverseCubicCutoffRadiusMie ) / cAvogadro / 3.D0

! Total number of association sites
tAssocSites = 4 * nParticles / 2 ! Removes redundant counting

! Write simulation properties
WRITE( *, "(A)" ) " "
WRITE( *, "(A)" ) "-------------------------------"
WRITE( *, "(A)" ) "------ System Parameters ------"
WRITE( *, "(A)" ) "-------------------------------"
WRITE( *, "(2G0)" ) "Number of particles: ", nParticles
WRITE( *, "(2G0.5)" ) "Fluid density (1/Å³): ", Density
WRITE( *, "(2G0.5)" ) "Solid density (1/Å³): ", SolidDensity
WRITE( *, "(2G0.5)" ) "Temperature (K): ", Temperature
WRITE( *, "(2G0.5)" ) "Pore height (Å): ", poreHeight
WRITE( *, "(2G0.5)" ) "Pore width (Å): ", poreWidth
WRITE( *, "(2G0.5)" ) "Pore length (Å): ", poreLength
WRITE( *, "(2G0.5)" ) "Pore volume (Å³): ", BoxVolume
WRITE( *, "(2G0.5)" ) "Interlayer spacing (Å): ", InterlayerSpacing
WRITE( *, "(2G0.5)" ) "Number of slabs: ", nSlabs
IF( isSolidPotentialEnabled ) THEN
  WRITE( *, "(A)" ) "Solid-fluid potential: [ENABLED]"
ELSE
  WRITE( *, "(A)" ) "Solid-fluid potential: [DISABLED]"
END IF
WRITE( *, "(A)" ) " "
WRITE( *, "(A)" ) "-------------------------------"
WRITE( *, "(A)" ) "------- Model Parameters ------"
WRITE( *, "(A)" ) "-------------------------------"
WRITE( *, "(2G0.5)" ) "Fluid diameter (Å): ", SigmaFluid
WRITE( *, "(2G0.5)" ) "Solid diameter (Å): ", SigmaSolid
WRITE( *, "(2G0.5)" ) "Fluid well depth (K): ", EpsilonFluid
WRITE( *, "(2G0.5)" ) "Solid well depth (K): ", EpsilonSolid
WRITE( *, "(2G0.5)" ) "Association well depth (K): ", EpsilonWater
WRITE( *, "(2G0.5)" ) "Site radius (Å): ", SiteRadius
WRITE( *, "(2G0.5)" ) "Bond angle (°): ", AngleWater * ( 180.D0 / cPi )
WRITE( *, "(2G0.5)" ) "Cutoff radius (association) (Å): ", cutoffRadiusWater
WRITE( *, "(2G0.5)" ) "Cutoff radius (Mie) (Å): ", cutoffRadiusMie
WRITE( *, "(2G0.5)" ) "Attractive exponent (Mie): ", attractiveRange
WRITE( *, "(2G0.5)" ) "Repulsive exponent (Mie): ", repulsiveRange
WRITE( *, "(2G0.5)" ) "Square-well range: ", SWRange
IF( isMiePotentialEnabled ) THEN
   WRITE( *, "(A)" ) "Attractive Mie potential: [ENABLED]"
ELSE
  WRITE( *, "(A)" ) "Attractive potential: [DISABLED]"
END IF
IF( isSWPotentialEnabled ) THEN
  WRITE( *, "(A)" ) "Square-well potential: [ENABLED]"
ELSE
  WRITE( *, "(A)" ) "Square-well potential: [DISABLED]"
END IF
IF( isHardSpherePotentialEnabled ) THEN
  WRITE( *, "(A)" ) "Hard-sphere potential: [ENABLED]"
ELSE
  WRITE( *, "(A)" ) "Hard-sphere potential: [DISABLED]"
END IF
IF( isAssociativePotentialEnabled ) THEN
  WRITE( *, "(A)" ) "Association potential: [ENABLED]"
ELSE
  WRITE( *, "(A)" ) "Association potential: [DISABLED]"
END IF
WRITE( *, "(A)" ) " "
WRITE( *, "(A)" ) "-------------------------------"
WRITE( *, "(A)" ) "---- Simulation Parameters ----"
WRITE( *, "(A)" ) "-------------------------------"
WRITE( *, "(2G0)" ) "Number of MC cycles: ", nCycles
WRITE( *, "(2G0)" ) "Number of equilibration cycles: ", nEquilibration
WRITE( *, "(2G0)" ) "Number of equilibration cycles (random configuration): ", nEquilRandom
WRITE( *, "(2G0)" ) "Adjustment frequency: ", nAdjustment
WRITE( *, "(2G0)" ) "Saving frequency: ", nSave
WRITE( *, "(2G0)" ) "Configuration saving frequency: ", nSaveConfig
WRITE( *, "(2G0.5)" ) "Movement probability: ", MovementProbability
WRITE( *, "(2G0.5)" ) "Volume change probability: ", VolumeChangeProbability
WRITE( *, "(2G0.5)" ) "Translational probability: ", TranslationalProbability
WRITE( *, "(2G0.5)" ) "Rotational probability: ", RotationalProbability
WRITE( *, "(2G0.5)" ) "Isotropic volume change probability: ", IsotropicProbability
WRITE( *, "(2G0.5)" ) "Anisotropic volume change probability: ", AnisotropicProbability
WRITE( *, "(2G0.5)" ) "Particle creation/destruction probability: ", PCDProbability
WRITE( *, "(2G0.5)" ) "Creation probability: ", CreationProbability
WRITE( *, "(2G0.5)" ) "Destruction probability: ", DestructionProbability
WRITE( *, "(2G0.5)" ) "Acceptance ratio threshold (translation): ", AcceptanceRatioTranslation
WRITE( *, "(2G0.5)" ) "Acceptance ratio threshold (rotation): ", AcceptanceRatioRotation
WRITE( *, "(2G0.5)" ) "Acceptance ratio threshold (isotropic volume change): ", AcceptanceRatioIsotropic
WRITE( *, "(2G0.5)" ) "Acceptance ratio threshold (anisotropic volume change): ", AcceptanceRatioAnisotropic
WRITE( *, "(2G0.5)" ) "Acceptance ratio threshold (creation): ", AcceptanceRatioCreation
WRITE( *, "(2G0.5)" ) "Acceptance ratio threshold (destruction): ", AcceptanceRatioDestruction
WRITE( *, "(2G0.5)" ) "Maximum box edge ratio (random configuration): ", BoxEdgeMaxRatio
IF( RandomConf ) THEN
  WRITE( *, "(A)" ) "Generate random configuration: [YES]"
ELSE
  WRITE( *, "(A)" ) "Generate random configuration: [NO]"
END IF
IF( RandomSeed ) THEN
  WRITE( *, "(A)" ) "Random seed: [YES]"
ELSE
  WRITE( *, "(A)" ) "Random seed: [NO]"
END IF
IF( ApplyNPTRandom ) THEN
  WRITE( *, "(A)" ) "Apply NPT moves (random configuration): [YES]"
ELSE
  WRITE( *, "(A)" ) "Apply NPT moves (random configuration): [NO]"
END IF
IF( Confinement ) THEN
  WRITE( *, "(A)" ) "Confinement enabled: [YES]"
ELSE
  WRITE( *, "(A)" ) "Confinement enabled: [NO]"
END IF
IF( isWidomEnabled ) THEN
  WRITE( *, "(A)" ) "Widom test particle insertion: [ENABLED]"
ELSE
  WRITE( *, "(A)" ) "Widom test particle insertion: [DISABLED]"
END IF
WRITE( *, "(A)" ) " "

RETURN

END SUBROUTINE Variable_Initialization