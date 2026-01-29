! ------------------------------------------------------------------------- !
! Random Configuration Subroutine                                           !
! ------------------------------------------------------------------------- !
SUBROUTINE RandomConfiguration(  )

USE GlobalVariables
USE VectorOperations
USE QuaternionOperations

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES (GENERAL)                                                                     !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: iParticle, pParticle ! Counters (particle)
INTEGER( Kind= Int64 ) :: OverlappingParticles ! Counter of overlapping particles
INTEGER( Kind= Int64 ) :: CounterEquilibration ! Counter of overlapping particles
INTEGER( Kind= Int64 ) :: iSite, bEdge

! *********************************************************************************************** !
! INTEGER VARIABLES (PARAMETER)                                                                   !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ), PARAMETER :: FixPFraction = 1 ! Control variable for the algorithm that fixes the packing fraction of the system

! *********************************************************************************************** !
! INTEGER VARIABLES (MONTE CARLO PARAMETERS)                                                      !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: pCycle, BoxMatrixComponent                      ! Counter of cycles
INTEGER( Kind= Int64 ) :: nAttempts                   ! Counter
INTEGER( Kind= Int64 ) :: nAcceptanceTranslation      ! Move acceptance counter: Translation
INTEGER( Kind= Int64 ) :: nAcceptanceRotation
INTEGER( Kind= Int64 ) :: nAcceptanceIsotropicVolumeChange
INTEGER( Kind= Int64 ) :: nAcceptanceAnisotropicVolumeChange
INTEGER( Kind= Int64 ) :: nMovementTranslationCounter
INTEGER( Kind= Int64 ) :: nMovementRotationCounter
INTEGER( Kind= Int64 ) :: nMovementIsoVolumeChangeCounter
INTEGER( Kind= Int64 ) :: nMovementAnisoVolumeChangeCounter


! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                     :: iOldPotential, iOldVirial, iOldAssociationTerm, iOldMieSWTerm, iOldSteeleTerm
REAL( Kind= Real64 )                     :: iNewPotential, iNewVirial, iNewAssociationTerm, iNewMieSWTerm, iNewSteeleTerm
REAL( Kind= Real64 )                     :: PotentialEnergyDifference, AssociationEnergyDifference, MieEnergyDifference
REAL( Kind= Real64 )                     :: tPotentialEnergy, VirialContribution, tAssociationTerm, tMieSWTerm, tSteeleTerm
REAL( Kind= Real64 )                     :: OldPotentialEnergy, NewPotentialEnergy, OldVirialContribution, NewVirialContribution
REAL( Kind= Real64 )                     :: OldAssociationTerm, NewAssociationTerm, OldMieTerm, NewMieSWTerm, VolumeScalingFactor
REAL( Kind= Real64 )                     :: OldSteeleTerm, NewSteeleTerm, DensityRandom
REAL( Kind= Real64 )                     :: SteeleEnergyDifference, VirialDifference
REAL( Kind= Real64 )                     :: OldBoxLength(9), OldBoxLengthInverse(9), OldBoxVolume
REAL( Kind= Real64 )                     :: NewBoxLength(9), NewBoxLengthInverse(9), NewBoxVolume
REAL( Kind= Real64 )                     :: BoxEdgeLength(9), BoxEdgeRatio(9)
REAL( Kind= Real64 )                     :: Ratio, EnthalpyChange, EnergyChange                        ! Acceptance ratio (simulation)
REAL( Kind= Real64 )                     :: MaxTranslationalDisplacement ! Maximum displacement [+/-] (Translation)
REAL( Kind= Real64 )                     :: MaxAngularDisplacement       ! Maximum displacement [+/-] (Rotation)
REAL( Kind= Real64 )                     :: MaxIsoVolumetricDisplacement    ! Maximum volume change
REAL( Kind= Real64 )                     :: MaxAnisoVolumetricDisplacement  ! Maximum volume change
REAL( Kind= Real64 )                     :: SquaredDistance              ! Squared distance
REAL( Kind= Real64 )                     :: BoxVolumeRandom              ! Box volume (random configuration)
REAL( Kind= Real64 ), DIMENSION( 3 )     :: VectorDistance               ! Vector distance
REAL( Kind= Real64 ), DIMENSION( 3 )     :: ScalingDistanceUnitBox       ! Scaled position (unit box)
REAL( Kind= Real64 ), DIMENSION( 3 )     :: iOldPosition, iNewPosition   ! Position of particle (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( 0:3 )   :: iOldQuaternion, iNewQuaternion   ! Position of particle (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( 3 )     :: iOldOrientation, iNewOrientation   ! Position of particle (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( 3 )     :: RotatedPosition
REAL( Kind= Real64 ), DIMENSION( 9 )     :: BoxLengthRandom              ! Box length (random configuration)
REAL( Kind= Real64 ), DIMENSION( 9 )     :: BoxLengthInverseRandom       ! Inverse of box length (random configuration)
REAL( Kind= Real64 ), DIMENSION( 3,4 )   :: iOldSitePosition, iNewSitePosition   ! Position of particle (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( 3,4 )   :: iOldSiteOrientation, iNewSiteOrientation   ! Position of particle (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( 0:3,4 ) :: iOldSiteQuaternion, iNewSiteQuaternion   ! Position of particle (before/after a trial move)
REAL( Kind= Real64 ), DIMENSION( 3, nParticles ) :: pPositionRandomSave
REAL( Kind= Real64 ), DIMENSION( 3, 4, nParticles ) :: sitePositionRandomSave


! *********************************************************************************************** !
! CHARACTER STRINGS                                                                               !
! *********************************************************************************************** !
CHARACTER( LEN= 140 ) :: DescriptorString ! Descriptor for strings

! *********************************************************************************************** !
! LOGICAL VARIABLES                                                                               !
! *********************************************************************************************** !
LOGICAL :: Overlap, PrivateOverlap, SharedOverlap ! Detects overlap between two particles: TRUE = overlap detected; FALSE = overlap not detected
LOGICAL :: MovementTranslationLogical             ! Translational movement selection : TRUE = movement selected; FALSE = movement not selected
LOGICAL :: MovementRotationLogical                ! Rotational movement selection    : TRUE = movement selected; FALSE = movement not selected
LOGICAL :: MovementIsoVolumeChangeLogical         ! Isotropic volume change selection: TRUE = movement selected; FALSE = movement not selected
LOGICAL :: MovementAnisoVolumeChangeLogical       ! Anisotropic volume change selection: TRUE = movement selected; FALSE = movement not selected
LOGICAL :: CheckBoxDistortion       ! Anisotropic volume change selection: TRUE = movement selected; FALSE = movement not selected

! *********************************************************************************************** !
! LOGICAL VARIABLES (ALLOCATABLE)                                                                 !
! *********************************************************************************************** !
LOGICAL, DIMENSION( : ), ALLOCATABLE :: OverlapCounterLogical ! Checks how many particles are overlapping each other during the hit-and-miss algorithm

! Allocation
ALLOCATE( OverlapCounterLogical(nParticles) )

! Initialization
OverlapCounterLogical = .TRUE.
OverlappingParticles  = COUNT( OverlapCounterLogical, Dim= 1 ) - 1
CounterEquilibration = 0

! Positioning of particles (centers of mass)
pPosition = 0.D0

! *********************************************************************************************** !
! Monte Carlo parameters (NVT simulation)                                                         !
! *********************************************************************************************** !
MovementTranslationLogical  = .FALSE.      ! Translational move selector           (initial value)
nAcceptanceTranslation      = 0            ! Translational move acceptance counter (initial value)
nMovementTranslationCounter = 0            ! Translational move counter            (initial value)
nAttempts                   = 0            ! Number of attempts                    (initial value)
pPositionRandom             = pPosition    ! Position of particles                 (initial value)

! Maximum translational displacement (initial value)
MaxTranslationalDisplacement = maxTranslationalDisplcRandom

! Recreate box
IF( .NOT. Confinement ) THEN
  BoxVolume = DBLE( nParticles ) / Density
  BoxLength = 0.D0
  BoxLength(1) = BoxVolume ** ( 1.D0 / 3.D0 )
  BoxLength(5) = BoxLength(1)
  BoxLength(9) = BoxLength(1)
  CALL InverseMatrixCofactorVec( BoxLength, BoxLengthInverse, BoxVolume )
END IF

! Summary
WRITE( *, "(A)" ) "-------------------------------"
WRITE( *, "(A)" ) "---- Initial Configuration ----"
WRITE( *, "(A)" ) "-------------------------------"
WRITE( *, "(G0)" ) "Attempting to randomly distribute particles inside a cubic box. It may take a while..."
WRITE( *, "(G0)" ) " "

! Open initial configuration file
OPEN( Unit= 55, File= "initial_configuration.xyz", Status= "REPLACE" )

! *********************************************************************************************** !
! Hit-and-miss (NVT simulation)                                                                   !
! *********************************************************************************************** !
HitAndMissNVT: DO

  ! Loop over particles
  DO pCycle = 1, nParticles

    ! Rewind output unit
    REWIND( 55 )

    ! Movement selection
    MovementTranslationLogical  = .TRUE. ! Enable translation
    nMovementTranslationCounter = nMovementTranslationCounter + 1 ! Increment move counter

    ! Pseudorandom number generator (uniform distribution)
    CALL Random_Number( RandomNumber )
    ! Random selection of particles of component i
    iParticle = INT( RandomNumber * DBLE( nParticles ) ) + 1

    ! Assignment of previous configuration (microstate m)
    iOldPosition(:) = pPositionRandom(:,iParticle) ! Old position

    ! Translational movement
    IF( MovementTranslationLogical ) THEN
      ! Random translation along x-axis
      IF( .NOT. xyTranslationIndependency ) THEN
        CALL Random_Number( RandomNumber )
        iNewPosition(1) = iOldPosition(1) + ( ( 2.D0 * RandomNumber ) - 1.D0 ) * MaxTranslationalDisplacement ! Range [-drmax,drmax]
        ! Random translation along y-axis
        CALL Random_Number( RandomNumber )
        iNewPosition(2) = iOldPosition(2) + ( ( 2.D0 * RandomNumber ) - 1.D0 ) * MaxTranslationalDisplacement ! Range [-drmax,drmax]
        ! Random translation along z-axis
        CALL Random_Number( RandomNumber )
        iNewPosition(3) = iOldPosition(3) + ( ( 2.D0 * RandomNumber ) - 1.D0 ) * MaxTranslationalDisplacement ! Range [-drmax,drmax]
      ELSE
        CALL Random_Number( RandomNumber )
        ! Random translation along the XY-plane
        IF( RandomNumber < TranslationalProbabilityXY ) THEN
          CALL Random_Number( RandomNumber )
          iNewPosition(1) = iOldPosition(1) + ( ( 2.D0 * RandomNumber ) - 1.D0 ) * MaxTranslationalDisplacement ! Range [-drmax,drmax]
          ! Random translation along y-axis
          CALL Random_Number( RandomNumber )
          iNewPosition(2) = iOldPosition(2) + ( ( 2.D0 * RandomNumber ) - 1.D0 ) * MaxTranslationalDisplacement ! Range [-drmax,drmax]
          ! Z is the same
          iNewPosition(3) = iOldPosition(3)
        ! Random translation along Z
        ELSE
          ! X and Y are the same
          iNewPosition(1) = iOldPosition(1)
          iNewPosition(2) = iOldPosition(2)
          ! Random translation along z-axis
          CALL Random_Number( RandomNumber )
          iNewPosition(3) = iOldPosition(3) + ( ( 2.D0 * RandomNumber ) - 1.D0 ) * MaxTranslationalDisplacement ! Range [-drmax,drmax]
        END IF
      END IF
      ! Prevent molecules from approximating too close to the pore walls
      IF( Confinement ) THEN
        IF( (iNewPosition(3) + 0.5D0 * poreWidth) < (0.5D0 * SigmaFluid) ) THEN
          iNewPosition(3) = 0.D0
        ELSE IF( (iNewPosition(3) + 0.5D0 * poreWidth) > (poreWidth - 0.5D0 * SigmaFluid) ) THEN
          iNewPosition(3) = 0.D0
        END IF
      END IF
      ! Minimum image convention
      CALL MatrixVectorMultiplication( BoxLengthInverse, iNewPosition, ScalingDistanceUnitBox ) ! Spatial transformation
      IF( Confinement ) THEN
        ScalingDistanceUnitBox(1:2) = ScalingDistanceUnitBox(1:2) - ANINT( ScalingDistanceUnitBox(1:2) ) ! Periodic boundary conditions (system is confined along z-axis)
      ELSE
        ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
      END IF
      CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, iNewPosition ) ! Spatial transformation
    END IF

    ! Overlap check
    Overlap = .FALSE.
    SharedOverlap = .FALSE.
    DO pParticle = 1, nParticles
      ! Cycle if a thread founds an overlap
      IF( SharedOverlap ) CYCLE
      ! Skip itself
      IF( iParticle == pParticle ) CYCLE
      ! Vector distance
      VectorDistance(:) = pPositionRandom(:,pParticle) - iNewPosition(:)
      ! Apply minimum image convention
      CALL MatrixVectorMultiplication( BoxLengthInverse, VectorDistance, ScalingDistanceUnitBox ) ! Spatial transformation
      IF( Confinement ) THEN
        ScalingDistanceUnitBox(1:2) = ScalingDistanceUnitBox(1:2) - ANINT( ScalingDistanceUnitBox(1:2) ) ! Periodic boundary conditions (system is confined along z-axis)
      ELSE
        ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
      END IF
      CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, VectorDistance ) ! Spatial transformation
      ! Squared distance
      SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
      ! Check for overlap
      IF( SquaredDistance <= squaredSigmaFluid ) THEN
        PrivateOverlap = .TRUE.
        ! Overlap found
        IF( PrivateOverlap ) SharedOverlap = .TRUE.
      END IF
    END DO

    ! Set overlap state
    IF( SharedOverlap ) THEN
      Overlap = .TRUE.
    END IF

    ! Acceptance criterion
    IF( .NOT. Overlap ) THEN
      ! System configuration update
      pPositionRandom(:,iParticle) = iNewPosition(:) ! Update position
      ! Update counter of overlapping configurations
      IF( OverlapCounterLogical(iParticle) .AND. CounterEquilibration < 1 ) THEN
        OverlapCounterLogical(iParticle) = .FALSE.
        OverlappingParticles = OverlappingParticles - 1
        IF( OverlappingParticles < 0 ) OverlappingParticles = 0
      END IF
      ! Displacement counter update
      IF( MovementTranslationLogical ) THEN
        nAcceptanceTranslation = nAcceptanceTranslation + 1 ! Translational move counter
      END IF
    ELSE
      ! Retrieve old configuration
      pPositionRandom(:,iParticle) = iOldPosition(:) ! Retrieve old position
    END IF

  END DO

  ! Adjustment of maximum displacement (translation and rotation)
  IF( MOD( nMovementTranslationCounter, nAdjustment ) == 0 .AND. CounterEquilibration < 1 ) THEN

    ! Acceptance ratio (translation)
    IF( nMovementTranslationCounter > 0 ) THEN
      Ratio = DBLE( nAcceptanceTranslation ) / DBLE( nMovementTranslationCounter )
      ! Translational adjustment
      IF( Ratio <= AcceptanceRatioTranslation ) THEN
        MaxTranslationalDisplacement = 0.95D0 * MaxTranslationalDisplacement
      ELSE
        MaxTranslationalDisplacement = 1.05D0 * MaxTranslationalDisplacement
      END IF
      ! Set minimum translational displacement (arbitrary)
      IF( MaxTranslationalDisplacement <= SigmaFluid ) MaxTranslationalDisplacement = SigmaFluid
    END IF

    ! Reset counters
    nAcceptanceTranslation      = 0
    nMovementTranslationCounter = 0

  END IF

  ! Iteration
  IF( CounterEquilibration < 1 ) THEN
    nAttempts = nAttempts + 1
  END IF

  ! Initial configuration (partial)
  WRITE( 55, "(G0)" ) nParticles
  DescriptorString = "(G0,8(G0,1X),G0,G0,2(G0,1X),G0,2G0)"
  WRITE( 55, DescriptorString ) 'Lattice="', BoxLength(1:9), '" Origin="', -0.5D0 * ( BoxLength(1) + BoxLength(4) + &
  &                             BoxLength(7) ), -0.5D0 * ( BoxLength(2) + BoxLength(5) + BoxLength(8) ), -0.5D0 * &
  &                             ( BoxLength(3) + BoxLength(6) + BoxLength(9) ), '" ', &
  &                             "Properties=id:I:1:species:S:1:pos:R:3:orientation:R:4:aspherical_shape:R:3"
  DO iParticle = 1, nParticles
    WRITE( 55, "(12(G0,1X))" ) iParticle, "C", pPositionRandom(1,iParticle), pPositionRandom(2,iParticle), &
    &                          pPositionRandom(3,iParticle), 0.D0, 0.D0, 0.D0, 1.D0, 0.5D0 * SigmaFluid, 0.5D0 * SigmaFluid, &
    &                          0.5D0 * SigmaFluid
  END DO
  FLUSH( 55 )

  ! Progress bar
  IF( CounterEquilibration < 1 ) THEN
    CALL ProgressBarHitAndMiss( nAttempts, OverlappingParticles )
  ELSE IF( CounterEquilibration < 2 ) THEN
    ! Summary
    WRITE( *, "(G0)" ) " "
    WRITE( *, "(G0)" ) " "
    IF( nAttempts == 1 ) THEN
      WRITE( *, "(3G0)" ) "Possible random initial configuration found by the hit-and-miss algorithm with ", nAttempts, " attempt."
    ELSE IF( nAttempts > 1 ) THEN
      WRITE( *, "(3G0)" ) "Possible random initial configuration found by the hit-and-miss algorithm with ", nAttempts, " attempts."
    END IF
    WRITE( *, "(G0)" ) " "
    WRITE( *, "(G0)" ) "Now running equilibration phase of the hit-and-miss algorithm..."
    WRITE( *, "(G0)" ) " "
    CALL ProgressBarEquilibration( CounterEquilibration, "NVT" )
  ELSE
    CALL ProgressBarEquilibration( CounterEquilibration, "NVT" )
  END IF

  ! Check number of overlapping particles
  IF( OverlappingParticles == 0 ) THEN
    CounterEquilibration = CounterEquilibration + 1
  END IF

  ! Exit condition
  IF( CounterEquilibration > nEquilRandom ) THEN
    EXIT HitAndMissNVT
  END IF

END DO HitAndMissNVT

! Reassign positions of all particles
pPosition = pPositionRandom ! Position of particles

! Generate quaternions and sites
DO iParticle = 1, nParticles
  ! Random quaternions
  CALL RandomQuaternionGenerator( pQuaternion(:,iParticle), 2.D0 * cPi )
  ! Orientations
  CALL VectorRotation( zAxis, pQuaternion(:,iParticle), pOrientation(:,iParticle) )
  DO iSite = 1, 4
    CALL VectorRotation( ReferencePositions(:,iSite), pQuaternion(:,iParticle), RotatedPosition(:) )
    ! Site positions
    sitePosition(:,iSite,iParticle) = pPosition(:,iParticle) + RotatedPosition(:)
    ! Site orientations
    siteOrientation(:,iSite,iParticle) = RotatedPosition(:) / SQRT( DOT_PRODUCT( RotatedPosition(:), RotatedPosition(:) ) )
    ! Site quaternions
    CALL QuaternionMultiplication( pQuaternion(:,iParticle), ReferenceQuaternions(:,iSite), siteQuaternion(:,iSite,iParticle) )
  END DO
END DO

! Initial configuration (partial)
REWIND( 55 )
WRITE( 55, "(G0)" ) nParticles * 5
DescriptorString = "(G0,8(G0.9,1X),G0.9,G0,2(G0.9,1X),G0.9,2G0)"
WRITE( 55, DescriptorString ) 'Lattice="', BoxLength(1:9), '" Origin="', -0.5D0 * ( BoxLength(1) + BoxLength(4) + BoxLength(7) ), &
&                             -0.5D0 * ( BoxLength(2) + BoxLength(5) + BoxLength(8) ), -0.5D0 * ( BoxLength(3) + BoxLength(6) + &
&                             BoxLength(9) ), '" ', "Properties=id:I:1:species:S:1:pos:R:3:orientation:R:4:aspherical_shape:R:3"
DO iParticle = 1, nParticles
  WRITE( 55, "(12(G0,1X))" ) iParticle, "C", pPosition(1:3,iParticle), pQuaternion(1:3,iParticle), pQuaternion(0,iParticle), &
  &                          0.5D0 * SigmaFluid, 0.5D0 * SigmaFluid, 0.5D0 * SigmaFluid
  DO iSite = 1, 4
    WRITE( 55, "(12(G0,1X))" ) iParticle, SiteNames(iSite), sitePosition(1:3,iSite,iParticle), &
    &                          siteQuaternion(1:3,iSite,iParticle), siteQuaternion(0,iSite,iParticle), &
    &                          0.1D0 * SigmaFluid, 0.1D0 * SigmaFluid, 0.4D0 * SigmaFluid
  END DO
END DO

! Close output unit
CLOSE( 55 )

! Deallocation
DEALLOCATE( OverlapCounterLogical )

! Summary
WRITE( *, "(G0)" ) " "
WRITE( *, "(G0)" ) " "

! Apply NPT Monte Carlo if required
IF( ApplyNPTRandom ) THEN

  ! *********************************************************************************************** !
  ! Monte Carlo parameters (NPT simulation)                                                         !
  ! *********************************************************************************************** !
  MovementTranslationLogical         = .FALSE.                              ! Translational move selector           (initial value)
  MovementRotationLogical            = .FALSE.                              ! Rotational move selector              (initial value)
  MovementIsoVolumeChangeLogical     = .FALSE.                              ! Volume move selector                  (initial value)
  MovementAnisoVolumeChangeLogical   = .FALSE.                              ! Volume move selector                  (initial value)
  nAcceptanceTranslation             = 0                                    ! Translational move acceptance counter (initial value)
  nAcceptanceRotation                = 0                                    ! Rotational move acceptance counter    (initial value)
  nAcceptanceIsotropicVolumeChange   = 0                                    ! Volumetric move acceptance counter    (initial value)
  nAcceptanceAnisotropicVolumeChange = 0                                    ! Volumetric move acceptance counter    (initial value)
  nMovementTranslationCounter        = 0                                    ! Translational move counter            (initial value)
  nMovementRotationCounter           = 0                                    ! Rotational move counter               (initial value)
  nMovementIsoVolumeChangeCounter    = 0                                    ! Volume scaling counter                (initial value)
  nMovementAnisoVolumeChangeCounter  = 0                                    ! Volume scaling counter                (initial value)
  nAttempts                          = 0                                    ! Number of attempts                    (initial value)
  pQuaternionRandom                  = pQuaternion                          ! Quaternion algebra                    (initial value)
  pPositionRandom                    = pPosition                            ! Position of particles                 (initial value)
  pOrientationRandom                 = pOrientation                         ! Orientation of particles              (initial value)
  sitePositionRandom                 = sitePosition                         ! Site positions                        (initial value)
  siteOrientationRandom              = siteOrientation                      ! Site orientations                     (initial value)
  siteQuaternionRandom               = siteQuaternion                       ! Site quaternions                      (initial value)
  BoxLengthRandom                    = BoxLength                            ! Box length                            (initial value)
  BoxLengthInverseRandom             = BoxLengthInverse                     ! Inverse of box length                 (initial value)
  BoxVolumeRandom                    = BoxVolume                            ! Box volume                            (initial value)
  MaxTranslationalDisplacement       = maxTranslationalDisplcRandom         ! Maximum translational displacement    (initial value)
  MaxAngularDisplacement             = maxRotationalDisplcRandom            ! Maximum rotational displacement       (initial value)
  MaxIsoVolumetricDisplacement       = MaxIsoVolumetricDisplacementRandom   ! Maximum isovolumetric displacement    (initial value)
  MaxAnisoVolumetricDisplacement     = MaxAnisoVolumetricDisplacementRandom ! Maximum anisovolumetric displacement  (initial value)
  DensityRandom                      = DBLE( nParticles ) / BoxVolumeRandom

  ! Initial total energy calculation
  CALL Calculation_Total_Energy( BoxLengthRandom, BoxLengthInverseRandom, tPotentialEnergy, VirialContribution, tAssociationTerm, &
  &                              tMieSWTerm, tSteeleTerm, .FALSE., Overlap )
  IF( Overlap ) THEN
    WRITE( *, "(G0)" ) "Error: Overlapping particles detected in the initial configuration for NPT simulation."
    STOP
  END IF

  ! Open initial configuration file
  OPEN( Unit= 55, File= "initial_configuration.xyz", Status= "REPLACE" )

  ! Isothermal-isobaric Monte Carlo simulation
  SimulationNPT: DO

    ! Choose between displacement of molecules or volume scaling
    CALL Random_Number( RandomNumber )
    IF( RandomNumber < MovementProbability ) THEN
      MovementTranslationLogical       = .TRUE.  ! Enable translation
      MovementRotationLogical          = .TRUE.  ! Enable rotation
      MovementIsoVolumeChangeLogical   = .FALSE. ! Disable volume scaling
      MovementAnisoVolumeChangeLogical = .FALSE. ! Disable volume scaling
    ELSE IF( RandomNumber >= MovementProbability ) THEN
      MovementTranslationLogical       = .FALSE. ! Disable translation
      MovementRotationLogical          = .FALSE. ! Disable rotation
      MovementIsoVolumeChangeLogical   = .TRUE.  ! Enable volume scaling
      MovementAnisoVolumeChangeLogical = .TRUE.  ! Enable volume scaling
    END IF

    ! Movement (translation or rotation)
    IF( MovementTranslationLogical .OR. MovementRotationLogical ) THEN

      ! Pseudorandom number generator (uniform distribution)
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

      ! Pseudorandom number generator (uniform distribution)
      CALL Random_Number( RandomNumber )
      ! Random selection of particles of component cComponent
      iParticle = INT( RandomNumber * DBLE( nParticles ) ) + 1

      ! Assignment of previous configuration (microstate m)
      iOldPosition(:)          = pPositionRandom(:,iParticle)         ! Old position
      iOldQuaternion(:)        = pQuaternionRandom(:,iParticle)       ! Old quaternion
      iOldOrientation(:)       = pOrientationRandom(:,iParticle)      ! Old orientation
      iOldSitePosition(:,:)    = sitePositionRandom(:,:,iParticle)    ! Position (sites)
      iOldSiteOrientation(:,:) = siteOrientationRandom(:,:,iParticle) ! Orientation (sites)
      iOldSiteQuaternion(:,:)  = siteQuaternionRandom(:,:,iParticle)  ! Quaternion (sites)

      ! Translational movement
      IF( MovementTranslationLogical ) THEN
        ! Random translation along x-axis
        IF( .NOT. xyTranslationIndependency ) THEN
          CALL Random_Number( RandomNumber )
          iNewPosition(1) = iOldPosition(1) + ( ( 2.D0 * RandomNumber ) - 1.D0 ) * MaxTranslationalDisplacement ! Range [-drmax,drmax]
          ! Random translation along y-axis
          CALL Random_Number( RandomNumber )
          iNewPosition(2) = iOldPosition(2) + ( ( 2.D0 * RandomNumber ) - 1.D0 ) * MaxTranslationalDisplacement ! Range [-drmax,drmax]
          ! Random translation along z-axis
          CALL Random_Number( RandomNumber )
          iNewPosition(3) = iOldPosition(3) + ( ( 2.D0 * RandomNumber ) - 1.D0 ) * MaxTranslationalDisplacement ! Range [-drmax,drmax]
        ELSE
          CALL Random_Number( RandomNumber )
          ! Random translation along the XY-plane
          IF( RandomNumber < TranslationalProbabilityXY ) THEN
            CALL Random_Number( RandomNumber )
            iNewPosition(1) = iOldPosition(1) + ( ( 2.D0 * RandomNumber ) - 1.D0 ) * MaxTranslationalDisplacement ! Range [-drmax,drmax]
            ! Random translation along y-axis
            CALL Random_Number( RandomNumber )
            iNewPosition(2) = iOldPosition(2) + ( ( 2.D0 * RandomNumber ) - 1.D0 ) * MaxTranslationalDisplacement ! Range [-drmax,drmax]
            ! Z is the same
            iNewPosition(3) = iOldPosition(3)
          ! Random translation along Z
          ELSE
            ! X and Y are the same
            iNewPosition(1) = iOldPosition(1)
            iNewPosition(2) = iOldPosition(2)
            ! Random translation along z-axis
            CALL Random_Number( RandomNumber )
            iNewPosition(3) = iOldPosition(3) + ( ( 2.D0 * RandomNumber ) - 1.D0 ) * MaxTranslationalDisplacement ! Range [-drmax,drmax]
          END IF
        END IF
        ! Prevent molecules from approximating too close to the pore walls
        IF( Confinement ) THEN
          IF( (iNewPosition(3) + 0.5D0 * poreWidth) < (0.5D0 * SigmaFluid) ) THEN
            iNewPosition(3) = 0.D0
          ELSE IF( (iNewPosition(3) + 0.5D0 * poreWidth) > (poreWidth - 0.5D0 * SigmaFluid) ) THEN
            iNewPosition(3) = 0.D0
          END IF
        END IF
        ! Minimum image convention
        CALL MatrixVectorMultiplication( BoxLengthInverseRandom, iNewPosition, ScalingDistanceUnitBox ) ! Spatial transformation
        IF( Confinement ) THEN
          ScalingDistanceUnitBox(1:2) = ScalingDistanceUnitBox(1:2) - ANINT( ScalingDistanceUnitBox(1:2) ) ! Periodic boundary conditions (system is confined along z-axis)
        ELSE
          ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
        END IF
        CALL MatrixVectorMultiplication( BoxLengthRandom, ScalingDistanceUnitBox, iNewPosition ) ! Spatial transformation
      ! Disable translation
      ELSE IF( .NOT. MovementTranslationLogical ) THEN
        iNewPosition = iOldPosition
      END IF

      ! Rotational movement
      IF( MovementRotationLogical ) THEN
        ! Random composed unit quaternion
        CALL QuaternionCombination( iOldQuaternion, iNewQuaternion, MaxAngularDisplacement )
        ! Active transformation
        CALL VectorRotation( zAxis, iNewQuaternion, iNewOrientation )
      ! Disable rotation
      ELSE IF( .NOT. MovementRotationLogical ) THEN
        iNewQuaternion  = iOldQuaternion
        iNewOrientation = iOldOrientation
      END IF

      ! Site positions and orientations (after translation or rotation)
      DO iSite = 1, 4
        CALL VectorRotation( ReferencePositions(:,iSite ), iNewQuaternion, RotatedPosition(:) )
        iNewSitePosition(:,iSite) = iNewPosition(:) + RotatedPosition(:)
        iNewSiteOrientation(:,iSite) = RotatedPosition(:) / SQRT( DOT_PRODUCT( RotatedPosition(:), RotatedPosition(:) ) )
        CALL QuaternionMultiplication( iNewQuaternion, ReferenceQuaternions(:,iSite), iNewSiteQuaternion(:,iSite) )
      END DO

      ! Calculate old energy
      CALL Particle_Energy( iParticle, iOldPosition, iOldSitePosition, iOldPotential, iOldVirial, iOldAssociationTerm, &
      &                     iOldMieSWTerm, iOldSteeleTerm, BoxLengthRandom, BoxLengthInverseRandom, .TRUE., .FALSE., .FALSE., &
      &                     .FALSE., Overlap )
      IF( Overlap ) STOP "ERROR: Overlap detected in old configuration during displacement."

      ! Calculate new energy
      CALL Particle_Energy( iParticle, iNewPosition, iNewSitePosition, iNewPotential, iNewVirial, iNewAssociationTerm, &
      &                     iNewMieSWTerm, iNewSteeleTerm, BoxLengthRandom, BoxLengthInverseRandom, .TRUE., .TRUE., .FALSE., &
      &                     .FALSE., Overlap )

      ! Energy differences
      PotentialEnergyDifference = iNewPotential - iOldPotential
      AssociationEnergyDifference = iNewAssociationTerm - iOldAssociationTerm
      MieEnergyDifference = iNewMieSWTerm - iOldMieSWTerm
      SteeleEnergyDifference = iNewSteeleTerm - iOldSteeleTerm

      ! Virial differences
      VirialDifference = iNewVirial - iOldVirial

      ! Metropolis criterion
      IF( .NOT. Overlap ) THEN
        IF( PotentialEnergyDifference <= 75.D0 ) THEN
          IF( PotentialEnergyDifference <= 0.D0 ) THEN
            ! Update position
            pPositionRandom(:,iParticle) = iNewPosition(:)
            ! Update quaternion
            pQuaternionRandom(:,iParticle) = iNewQuaternion(:)
            ! Update orientation
            pOrientationRandom(:,iParticle) = iNewOrientation(:)
            ! Update site positions and orientations
            sitePositionRandom(:,:,iParticle)    = iNewSitePosition(:,:)
            siteOrientationRandom(:,:,iParticle) = iNewSiteOrientation(:,:)
            siteQuaternionRandom(:,:,iParticle)  = iNewSiteQuaternion(:,:)
            ! Update total energy
            tPotentialEnergy = tPotentialEnergy + PotentialEnergyDifference
            tAssociationTerm = tAssociationTerm + AssociationEnergyDifference
            tMieSWTerm = tMieSWTerm + MieEnergyDifference
            tSteeleTerm = tSteeleTerm + SteeleEnergyDifference
            ! Update virial
            VirialContribution = VirialContribution + VirialDifference
            ! Update displacement counter
            IF( MovementTranslationLogical ) THEN
              nAcceptanceTranslation = nAcceptanceTranslation + 1 ! Translational move counter
            ELSE IF( MovementRotationLogical ) THEN
              nAcceptanceRotation = nAcceptanceRotation + 1 ! Rotational move counter
            END IF
          ELSE
            CALL Random_Number( RandomNumber ) 
            IF( EXP( - PotentialEnergyDifference ) >= RandomNumber ) THEN
              ! Update position
              pPositionRandom(:,iParticle) = iNewPosition(:)
              ! Update quaternion
              pQuaternionRandom(:,iParticle) = iNewQuaternion(:)
              ! Update orientation
              pOrientationRandom(:,iParticle) = iNewOrientation(:)
              ! Update site positions and orientations
              sitePositionRandom(:,:,iParticle)    = iNewSitePosition(:,:)
              siteOrientationRandom(:,:,iParticle) = iNewSiteOrientation(:,:)
              siteQuaternionRandom(:,:,iParticle)  = iNewSiteQuaternion(:,:)
              ! Update total energy
              tPotentialEnergy = tPotentialEnergy + PotentialEnergyDifference
              tAssociationTerm = tAssociationTerm + AssociationEnergyDifference
              tMieSWTerm = tMieSWTerm + MieEnergyDifference
              tSteeleTerm = tSteeleTerm + SteeleEnergyDifference
              ! Update virial
              VirialContribution = VirialContribution + VirialDifference
              ! Update displacement counter
              IF( MovementTranslationLogical ) THEN
                nAcceptanceTranslation = nAcceptanceTranslation + 1 ! Translational move counter
              ELSE IF( MovementRotationLogical ) THEN
                nAcceptanceRotation = nAcceptanceRotation + 1 ! Rotational move counter
              END IF
            END IF
          END IF
        END IF
      END IF

    ! Volume scaling (isotropic or anisotropic)
    ELSE IF( MovementIsoVolumeChangeLogical .OR. MovementAnisoVolumeChangeLogical ) THEN

      ! Assignment of previous configuration (microstate m)
      OldBoxLength          = BoxLengthRandom        ! Box length
      OldBoxLengthInverse   = BoxLengthInverseRandom ! Box length (inverse)
      OldBoxVolume          = BoxVolumeRandom        ! Box volume
      OldPotentialEnergy    = tPotentialEnergy
      OldAssociationTerm    = tAssociationTerm
      OldMieTerm            = tMieSWTerm
      OldSteeleTerm         = tSteeleTerm
      OldVirialContribution = VirialContribution

      ! Expansion/compression type
      CALL Random_Number( RandomNumber )

      ! Isotropic volume scaling
      IF( RandomNumber < IsotropicProbability ) THEN
        ! Random scaling factor
        CALL Random_Number( RandomNumber )
        ! Random walk on the logarithm of the volume
        IF( Confinement ) THEN
          VolumeScalingFactor = DLOG( OldBoxLength(1) * OldBoxLength(5) ) + (2.D0 * RandomNumber - 1.D0) * &
          &                     MaxIsoVolumetricDisplacement
          VolumeScalingFactor = DEXP( VolumeScalingFactor )
          VolumeScalingFactor = SQRT( VolumeScalingFactor / ( OldBoxLength(1) * OldBoxLength(5) ) )
          ! Proportional box length
          NewBoxLength = OldBoxLength
          NewBoxLength(1) = OldBoxLength(1) * VolumeScalingFactor
          NewBoxLength(5) = OldBoxLength(5) * VolumeScalingFactor
        ELSE
          VolumeScalingFactor = DLOG( OldBoxVolume ) + (2.D0 * RandomNumber - 1.D0) * MaxIsoVolumetricDisplacement
          VolumeScalingFactor = DEXP( VolumeScalingFactor )
          VolumeScalingFactor = (VolumeScalingFactor / OldBoxVolume) ** (1.D0 / 3.D0)
          ! Proportional box length
          NewBoxLength = OldBoxLength * VolumeScalingFactor
        END IF
        ! New volume
        CALL InverseMatrixCofactorVec( NewBoxLength, NewBoxLengthInverse, NewBoxVolume )
        ! Movement counter
        nMovementIsoVolumeChangeCounter = nMovementIsoVolumeChangeCounter + 1
        ! Movement type
        MovementIsoVolumeChangeLogical   = .TRUE.  ! Enable isotropic volume scaling
        MovementAnisoVolumeChangeLogical = .FALSE. ! Disable anisotropic volume scaling
      ! Anisotropic volume scaling
      ELSE IF( RandomNumber >= IsotropicProbability ) THEN
        ! Random box component
        CALL Random_Number( RandomNumber )
        IF( Confinement ) THEN
          BoxMatrixComponent = INT( RandomNumber * 2.D0 ) + 1
        ELSE
          BoxMatrixComponent = INT( RandomNumber * 3.D0 ) + 1
        END IF
        IF( BoxMatrixComponent == 1 ) THEN
          BoxMatrixComponent = 1 ! XX component
        ELSE IF( BoxMatrixComponent == 2 ) THEN
          BoxMatrixComponent = 5 ! YY component
        ELSE IF( BoxMatrixComponent == 3 ) THEN
          BoxMatrixComponent = 9 ! YY component
        ELSE
          WRITE( *, "(G0)" ) "Error: wrong box matrix component selected during anisotropic volume scaling."
        END IF
        NewBoxLength = OldBoxLength
        ! Random factor
        CALL Random_Number( RandomNumber )
        NewBoxLength(BoxMatrixComponent) = OldBoxLength(BoxMatrixComponent) + MaxAnisoVolumetricDisplacement * &
        &                                  (2.D0 * RandomNumber - 1.D0)
        ! Calculate the new reciprocal box basis vectors and the volume of the system
        CALL InverseMatrixCofactorVec( NewBoxLength, NewBoxLengthInverse, NewBoxVolume )
        ! Movement counter
        nMovementAnisoVolumeChangeCounter = nMovementAnisoVolumeChangeCounter + 1
        ! Movement type
        MovementIsoVolumeChangeLogical   = .FALSE. ! Disable isotropic volume scaling
        MovementAnisoVolumeChangeLogical = .TRUE.  ! Enable anisotropic volume scaling
      END IF

      ! Condition of anisotropic volume scaling (box distortion)
      CheckBoxDistortion = .FALSE.
      IF( MovementAnisoVolumeChangeLogical ) THEN
        ! Box length
        BoxEdgeLength(1) = DSQRT( DOT_PRODUCT( NewBoxLength(1:3), NewBoxLength(1:3) ) )
        BoxEdgeLength(2) = DSQRT( DOT_PRODUCT( NewBoxLength(4:6), NewBoxLength(4:6) ) )
        BoxEdgeLength(3) = DSQRT( DOT_PRODUCT( NewBoxLength(7:9), NewBoxLength(7:9) ) )
        ! Length ratio
        BoxEdgeRatio(1) = BoxEdgeLength(1) / BoxEdgeLength(2)
        BoxEdgeRatio(2) = BoxEdgeLength(1) / BoxEdgeLength(3)
        BoxEdgeRatio(3) = BoxEdgeLength(2) / BoxEdgeLength(3)
        ! Avoid big distortions of the simulation box
        DO bEdge = 1, 1
          ! Length distortion
          IF( BoxEdgeRatio(bEdge) > BoxEdgeMaxRatio .OR. BoxEdgeRatio(bEdge) < 1.D0 / BoxEdgeMaxRatio ) THEN
            BoxLengthRandom        = OldBoxLength
            BoxLengthInverseRandom = OldBoxLengthInverse
            BoxVolumeRandom        = OldBoxVolume
            CheckBoxDistortion     = .TRUE.
            EXIT
          END IF
        END DO
      END IF

      ! Box not too distorted
      IF( .NOT. CheckBoxDistortion ) THEN

        ! Enthalpy change (weighing function)
        IF(  MovementIsoVolumeChangeLogical ) THEN
          EnthalpyChange = ( ReducedTargetPressure * ( NewBoxVolume - OldBoxVolume ) ) - &
          &                ( DBLE( nParticles + 1 ) * DLOG( NewBoxVolume / OldBoxVolume ) )
        ELSE IF( MovementAnisoVolumeChangeLogical ) THEN
          EnthalpyChange = ( ReducedTargetPressure * ( NewBoxVolume - OldBoxVolume ) ) - &
          &                ( DBLE( nParticles ) * DLOG( NewBoxVolume / OldBoxVolume ) )
        END IF

        ! System configuration update
        pPositionRandomSave = pPositionRandom
        sitePositionRandomSave = sitePositionRandom

        ! Isotropic volume scaling
        IF( MovementIsoVolumeChangeLogical ) THEN
          ! Rescale positions of all particles accordingly
          DO pParticle = 1, nParticles
            pPositionRandom(:,pParticle) = pPositionRandom(:,pParticle) * VolumeScalingFactor
            ! Site positions and orientations (after translation or rotation)
            DO iSite = 1, 4
              CALL VectorRotation( ReferencePositions(:,iSite), pQuaternionRandom(:,pParticle), RotatedPosition(:) )
              sitePositionRandom(:,iSite,pParticle) = pPositionRandom(:,pParticle) + RotatedPosition(:)
            END DO
          END DO
        ! Anisotropic volume scaling
        ELSE IF( MovementAnisoVolumeChangeLogical ) THEN
          ! Rescale positions of all particles accordingly
          DO pParticle = 1, nParticles
            ! Scaling coordinates using the old box length
            CALL MatrixVectorMultiplication( OldBoxLengthInverse, pPositionRandom(:,pParticle), ScalingDistanceUnitBox )
            ! New real coordinates using the new box length
            CALL MatrixVectorMultiplication( NewBoxLength, ScalingDistanceUnitBox, pPositionRandom(:,pParticle) )
            ! Site positions and orientations (after translation or rotation)
            DO iSite = 1, 4
              CALL VectorRotation( ReferencePositions(:,iSite), pQuaternionRandom(:,pParticle), RotatedPosition(:) )
              sitePositionRandom(:,iSite,pParticle) = pPositionRandom(:,pParticle) + RotatedPosition(:)
            END DO
          END DO
        END IF

        ! Potential
        CALL Calculation_Total_Energy( NewBoxLength, NewBoxLengthInverse, NewPotentialEnergy, NewVirialContribution, &
        &                              NewAssociationTerm, NewMieSWTerm, NewSteeleTerm, .TRUE., Overlap )

        ! Energy change
        EnergyChange = NewPotentialEnergy - OldPotentialEnergy

        ! Update enthalpy
        EnthalpyChange = EnthalpyChange + EnergyChange

        ! Metropolis criterion
        IF( .NOT. Overlap ) THEN
          IF( EnthalpyChange <= 75.D0 ) THEN
            IF( EnthalpyChange <= 0.D0 ) THEN
              ! Assigns the simulation box properties of a trial volume scaling to the system configuration
              BoxVolumeRandom        = NewBoxVolume        ! Update box volume
              BoxLengthRandom        = NewBoxLength        ! Update box length
              BoxLengthInverseRandom = NewBoxLengthInverse ! Update box length (inverse)
              ! Update total energy
              tPotentialEnergy = NewPotentialEnergy
              tAssociationTerm = NewAssociationTerm
              tMieSWTerm = NewMieSWTerm
              tSteeleTerm = NewSteeleTerm
              ! Update virial
              VirialContribution = NewVirialContribution
              ! Displacement counter update
              IF( MovementIsoVolumeChangeLogical ) THEN
                nAcceptanceIsotropicVolumeChange = nAcceptanceIsotropicVolumeChange + 1 ! Isotropic move counter
              ELSE IF( MovementAnisoVolumeChangeLogical ) THEN
                nAcceptanceAnisotropicVolumeChange = nAcceptanceAnisotropicVolumeChange + 1 ! Anisotropic move counter
              END IF
              ! Update packing fraction and reduced number density
              DensityRandom = ( DBLE( nParticles ) / BoxVolumeRandom )
              ! Re-initialization
              CheckBoxDistortion = .FALSE.
            ELSE
              CALL Random_Number( RandomNumber ) 
              IF( EXP( - EnthalpyChange ) >= RandomNumber ) THEN
                ! Assigns the simulation box properties of a trial volume scaling to the system configuration
                BoxVolumeRandom        = NewBoxVolume        ! Update box volume
                BoxLengthRandom        = NewBoxLength        ! Update box length
                BoxLengthInverseRandom = NewBoxLengthInverse ! Update box length (inverse)
                ! Update total energy
                tPotentialEnergy = NewPotentialEnergy
                tAssociationTerm = NewAssociationTerm
                tMieSWTerm = NewMieSWTerm
                tSteeleTerm = NewSteeleTerm
                ! Update virial
                VirialContribution = NewVirialContribution
                ! Displacement counter update
                IF( MovementIsoVolumeChangeLogical ) THEN
                  nAcceptanceIsotropicVolumeChange = nAcceptanceIsotropicVolumeChange + 1 ! Isotropic move counter
                ELSE IF( MovementAnisoVolumeChangeLogical ) THEN
                  nAcceptanceAnisotropicVolumeChange = nAcceptanceAnisotropicVolumeChange + 1 ! Anisotropic move counter
                END IF
                ! Update packing fraction and reduced number density
                DensityRandom = ( DBLE( nParticles ) / BoxVolumeRandom )
                ! Re-initialization
                CheckBoxDistortion = .FALSE.
              ELSE
                ! Retrieve old properties of the simulation box
                BoxVolumeRandom        = OldBoxVolume           ! Retrieve old box volume
                BoxLengthRandom        = OldBoxLength           ! Retrieve old box length
                BoxLengthInverseRandom = OldBoxLengthInverse    ! Retrieve old box length (inverse)
                pPositionRandom        = pPositionRandomSave    ! Retrieve old position of particles
                sitePositionRandom     = sitePositionRandomSave ! Retrieve old site positions
              END IF
            END IF
          ELSE
            ! Retrieve old properties of the simulation box
            BoxVolumeRandom        = OldBoxVolume           ! Retrieve old box volume
            BoxLengthRandom        = OldBoxLength           ! Retrieve old box length
            BoxLengthInverseRandom = OldBoxLengthInverse    ! Retrieve old box length (inverse)
            pPositionRandom        = pPositionRandomSave    ! Retrieve old position of particles
            sitePositionRandom     = sitePositionRandomSave ! Retrieve old site positions
          END IF
        ELSE
          ! Retrieve old properties of the simulation box
          BoxVolumeRandom        = OldBoxVolume           ! Retrieve old box volume
          BoxLengthRandom        = OldBoxLength           ! Retrieve old box length
          BoxLengthInverseRandom = OldBoxLengthInverse    ! Retrieve old box length (inverse)
          pPositionRandom        = pPositionRandomSave    ! Retrieve old position of particles
          sitePositionRandom     = sitePositionRandomSave ! Retrieve old site positions
        END IF

      END IF ! Box distortion criterion

    END IF

    ! Iteration
    nAttempts = nAttempts + 1

    ! Adjustment of maximum displacement (translation and rotation)
    IF( MOD( nAttempts, nAdjustment ) == 0 ) THEN

      ! Translational adjustment
      IF( nMovementTranslationCounter > 0 ) THEN
        ! Acceptance ratio (non-overlapping microstates over sampled microstates)
        Ratio = DBLE( nAcceptanceTranslation ) / DBLE( nMovementTranslationCounter )
        ! Translational adjustment
        IF( Ratio <= AcceptanceRatioTranslation ) THEN
          MaxTranslationalDisplacement = 0.95D0 * MaxTranslationalDisplacement
        ELSE
          MaxTranslationalDisplacement = 1.05D0 * MaxTranslationalDisplacement
        END IF
        ! Reset counter
        nAcceptanceTranslation      = 0
        nMovementTranslationCounter = 0
      END IF

      ! Avoid multiple turns (arbitrary)
      BoxEdgeLength(1) = DSQRT( DOT_PRODUCT( BoxLengthRandom(1:3), BoxLengthRandom(1:3) ) )
      BoxEdgeLength(2) = DSQRT( DOT_PRODUCT( BoxLengthRandom(4:6), BoxLengthRandom(4:6) ) )
      BoxEdgeLength(3) = DSQRT( DOT_PRODUCT( BoxLengthRandom(7:9), BoxLengthRandom(7:9) ) )
      IF( MaxTranslationalDisplacement >= 2.D0 * MAXVAL( BoxEdgeLength ) ) THEN
        MaxTranslationalDisplacement = MaxTranslationalDisplacement - MAXVAL( BoxEdgeLength )
      END IF

      ! Rotational adjustment
      IF( nMovementRotationCounter > 0 ) THEN
        ! Acceptance ratio (non-overlapping microstates over sampled microstates)
        Ratio = DBLE( nAcceptanceRotation ) / DBLE( nMovementRotationCounter )
        ! Rotational adjustment
        IF( Ratio <= AcceptanceRatioRotation ) THEN
          MaxAngularDisplacement = 0.95D0 * MaxAngularDisplacement
        ELSE
          MaxAngularDisplacement = 1.05D0 * MaxAngularDisplacement
        END IF
        ! Reset counter
        nAcceptanceRotation      = 0
        nMovementRotationCounter = 0
      END IF

      ! Avoid multiple turns (arbitrary)
      IF( MaxAngularDisplacement >= 4.D0 * cPi ) THEN
        MaxAngularDisplacement = MaxAngularDisplacement - 2.D0 * cPi
      END IF

    END IF

    ! Adjustment of maximum displacement (volume scaling)
    IF( MOD( nAttempts, nAdjustment ) == 0 ) THEN

      ! Volumetric adjustment (isotropic)
      IF( nMovementIsoVolumeChangeCounter > 0 ) THEN
        ! Acceptance ratio (non-overlapping microstates over sampled microstates)
        Ratio = DBLE( nAcceptanceIsotropicVolumeChange ) / DBLE( nMovementIsoVolumeChangeCounter )
        ! Volumetric adjustment
        IF( Ratio <= AcceptanceRatioIsotropic ) THEN
          MaxIsoVolumetricDisplacement = 0.95D0 * MaxIsoVolumetricDisplacement
        ELSE
          MaxIsoVolumetricDisplacement = 1.05D0 * MaxIsoVolumetricDisplacement
        END IF
        ! Reset counter
        nAcceptanceIsotropicVolumeChange = 0
        nMovementIsoVolumeChangeCounter  = 0
      END IF

      ! Volumetric adjustment (anisotropic)
      IF( nMovementAnisoVolumeChangeCounter > 0 ) THEN
        ! Acceptance ratio (non-overlapping microstates over sampled microstates)
        Ratio = DBLE( nAcceptanceAnisotropicVolumeChange ) / DBLE( nMovementAnisoVolumeChangeCounter )
        ! Volumetric adjustment
        IF( Ratio <= AcceptanceRatioAnisotropic ) THEN
          MaxAnisoVolumetricDisplacement = 0.95D0 * MaxAnisoVolumetricDisplacement
        ELSE
          MaxAnisoVolumetricDisplacement = 1.05D0 * MaxAnisoVolumetricDisplacement
        END IF
        ! Reset counter
        nAcceptanceAnisotropicVolumeChange = 0
        nMovementAnisoVolumeChangeCounter  = 0
      END IF

    END IF

    ! Summary
    CALL ProgressBarRandomConfigNPT( nAttempts, ReducedTargetPressure, DensityRandom )
 
    ! Initial configuration (partial)
    REWIND(55)
    WRITE( 55, "(G0)" ) nParticles * 5
    DescriptorString = "(G0,8(G0.9,1X),G0.9,G0,2(G0.9,1X),G0.9,2G0)"
    WRITE( 55, DescriptorString ) 'Lattice="', BoxLengthRandom(1:9), '" Origin="', -0.5D0 * ( BoxLengthRandom(1) + &
    &                             BoxLengthRandom(4) + BoxLengthRandom(7) ), -0.5D0 * ( BoxLengthRandom(2) + BoxLengthRandom(5) + &
    &                             BoxLengthRandom(8) ), -0.5D0 * ( BoxLengthRandom(3) + BoxLengthRandom(6) + BoxLengthRandom(9) &
    &                             ), '" ', "Properties=id:I:1:species:S:1:pos:R:3:orientation:R:4:aspherical_shape:R:3"
    DO iParticle = 1, nParticles
      WRITE( 55, "(12(G0,1X))" ) iParticle, "C", pPositionRandom(1:3,iParticle), pQuaternionRandom(1:3,iParticle), &
      &                          pQuaternionRandom(0,iParticle), 0.5D0 * SigmaFluid, 0.5D0 * SigmaFluid, 0.5D0 * SigmaFluid
      DO iSite = 1, 4
        WRITE( 55, "(12(G0,1X))" ) iParticle, SiteNames(iSite), sitePositionRandom(1:3,iSite,iParticle), &
        &                          siteQuaternionRandom(1:3,iSite,iParticle), siteQuaternionRandom(0,iSite,iParticle), &
        &                          0.1D0 * SigmaFluid, 0.1D0 * SigmaFluid, 0.4D0 * SigmaFluid
      END DO
    END DO
    FLUSH( 55 )

    ! Target packing fraction
    IF( nAttempts >= maxAttemptsNPT ) THEN
      EXIT SimulationNPT
    END IF

  END DO SimulationNPT

  ! Close output unit
  CLOSE( 55 )

  ! Assign final configuration
  BoxLength = BoxLengthRandom
  BoxLengthInverse = BoxLengthInverseRandom
  BoxVolume = BoxVolumeRandom
  pPosition = pPositionRandom
  pOrientation = pOrientationRandom
  pQuaternion = pQuaternionRandom
  sitePosition = sitePositionRandom
  siteOrientation = siteOrientationRandom
  siteQuaternion = siteQuaternionRandom

END IF

RETURN

END SUBROUTINE RandomConfiguration

! ----------------------------------------------------------------------------------------- !
! Subroutine to print the progress bar of the hit-and-miss algorithm                        !
! ----------------------------------------------------------------------------------------- !
SUBROUTINE ProgressBarHitAndMiss( iParticle, OverlappingParticles )

! Use Output_Unit
USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: Output_Unit
USE GlobalVariables

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: iParticle            ! Counter (particle)
INTEGER( Kind= Int64 ) :: OverlappingParticles ! Counter of overlapping particles
INTEGER( Kind= Int64 ) :: AuxiliarInt1         ! Auxiliar
INTEGER( Kind= Int64 ) :: AuxiliarInt2         ! Auxiliar

! *********************************************************************************************** !
! CHARACTER STRINGS                                                                               !
! *********************************************************************************************** !
CHARACTER( LEN= 59 ) :: ProgressBar ! Progress bar
CHARACTER( LEN= 02 ) :: StringSize1 ! String size
CHARACTER( LEN= 08 ) :: StringSize2 ! String size

! Initialization
AuxiliarInt1 = 0
AuxiliarInt2 = 0

! Progress bar (FORMAT)
IF( iParticle < 10 ) THEN
  AuxiliarInt1 = 0
  ProgressBar(1:(12+AuxiliarInt1)) = "Attempts: ? "
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt1)), Fmt= "(I1)" ) iParticle
ELSE IF( iParticle < 100 ) THEN
  AuxiliarInt1 = 1
  ProgressBar(1:(12+AuxiliarInt1)) = "Attempts: ?? "
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt1)), Fmt= "(I2)" ) iParticle
ELSE IF( iParticle < 1000 ) THEN
  AuxiliarInt1 = 2
  ProgressBar(1:(12+AuxiliarInt1)) = "Attempts: ??? "
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt1)), Fmt= "(I3)" ) iParticle
ELSE IF( iParticle < 10000 ) THEN
  AuxiliarInt1 = 3
  ProgressBar(1:(12+AuxiliarInt1)) = "Attempts: ???? "
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt1)), Fmt= "(I4)" ) iParticle
ELSE IF( iParticle < 100000 ) THEN
  AuxiliarInt1 = 4
  ProgressBar(1:(12+AuxiliarInt1)) = "Attempts: ????? "
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt1)), Fmt= "(I5)" ) iParticle
ELSE IF( iParticle < 1000000 ) THEN
  AuxiliarInt1 = 5
  ProgressBar(1:(12+AuxiliarInt1)) = "Attempts: ?????? "
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt1)), Fmt= "(I6)" ) iParticle
ELSE IF( iParticle < 10000000 ) THEN
  AuxiliarInt1 = 6
  ProgressBar(1:(12+AuxiliarInt1)) = "Attempts: ??????? "
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt1)), Fmt= "(I7)" ) iParticle
ELSE IF( iParticle < 100000000 ) THEN
  AuxiliarInt1 = 7
  ProgressBar(1:(12+AuxiliarInt1)) = "Attempts: ???????? "
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt1)), Fmt= "(I8)" ) iParticle
ELSE IF( iParticle < 1000000000 ) THEN
  AuxiliarInt1 = 8
  ProgressBar(1:(12+AuxiliarInt1)) = "Attempts: ????????? "
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt1)), Fmt= "(I9)" ) iParticle
ELSE IF( iParticle >= 1000000000 ) THEN
  AuxiliarInt1 = 10
  ProgressBar(1:(12+AuxiliarInt1)) = "Attempts: > 1 billion "
END IF
IF( OverlappingParticles < 10 ) THEN
  AuxiliarInt2 = 0
  ProgressBar((13+AuxiliarInt1):(38+AuxiliarInt1+AuxiliarInt2)) = "| Overlapping Particles: ?"
  WRITE( Unit= ProgressBar((38+AuxiliarInt1):(38+AuxiliarInt1+AuxiliarInt2)), Fmt= "(I1)" ) OverlappingParticles
ELSE IF( OverlappingParticles < 100 ) THEN
  AuxiliarInt2 = 1
  ProgressBar((13+AuxiliarInt1):(38+AuxiliarInt1+AuxiliarInt2)) = "| Overlapping Particles: ??"
  WRITE( Unit= ProgressBar((38+AuxiliarInt1):(38+AuxiliarInt1+AuxiliarInt2)), Fmt= "(I2)" ) OverlappingParticles
ELSE IF( OverlappingParticles < 1000 ) THEN
  AuxiliarInt2 = 2
  ProgressBar((13+AuxiliarInt1):(38+AuxiliarInt1+AuxiliarInt2)) = "| Overlapping Particles: ???"
  WRITE( Unit= ProgressBar((38+AuxiliarInt1):(38+AuxiliarInt1+AuxiliarInt2)), Fmt= "(I3)" ) OverlappingParticles
ELSE IF( OverlappingParticles < 10000 ) THEN
  AuxiliarInt2 = 3
  ProgressBar((13+AuxiliarInt1):(38+AuxiliarInt1+AuxiliarInt2)) = "| Overlapping Particles: ????"
  WRITE( Unit= ProgressBar((38+AuxiliarInt1):(38+AuxiliarInt1+AuxiliarInt2)), Fmt= "(I4)" ) OverlappingParticles
ELSE IF( OverlappingParticles < 100000 ) THEN
  AuxiliarInt2 = 4
  ProgressBar((13+AuxiliarInt1):(38+AuxiliarInt1+AuxiliarInt2)) = "| Overlapping Particles: ?????"
  WRITE( Unit= ProgressBar((38+AuxiliarInt1):(38+AuxiliarInt1+AuxiliarInt2)), Fmt= "(I5)" ) OverlappingParticles
ELSE IF( OverlappingParticles < 1000000 ) THEN
  AuxiliarInt2 = 5
  ProgressBar((13+AuxiliarInt1):(38+AuxiliarInt1+AuxiliarInt2)) = "| Overlapping Particles: ??????"
  WRITE( Unit= ProgressBar((38+AuxiliarInt1):(38+AuxiliarInt1+AuxiliarInt2)), Fmt= "(I6)" ) OverlappingParticles
ELSE IF( OverlappingParticles >= 1000000 ) THEN
  AuxiliarInt2 = 10
  ProgressBar((13+AuxiliarInt1):(38+AuxiliarInt1+AuxiliarInt2)) = "| Overlapping Particles: > 1 million"
END IF
ProgressBar((39+AuxiliarInt1+AuxiliarInt2):59) = REPEAT( " ", ( (20 - AuxiliarInt1 - AuxiliarInt2) + 1 ) )

! Print progress bar
WRITE( StringSize1, "(I0.2)" ) 38 + AuxiliarInt1 + AuxiliarInt2 + 1
StringSize2 = "(A1,A"//TRIM( StringSize1 )//")"
WRITE( Unit= Output_Unit, Fmt= StringSize2, Advance= "NO" ) CHAR(13), ProgressBar(1:(38+AuxiliarInt1+AuxiliarInt2+1))

! Flush standard output unit
FLUSH( Unit= Output_Unit )

RETURN

END SUBROUTINE ProgressBarHitAndMiss

! ----------------------------------------------------------------------------------------- !
! Subroutine to print the progress bar of the equilibration stage of the hit-and-miss       !
! algorithm                                                                                 !
! ----------------------------------------------------------------------------------------- !
SUBROUTINE ProgressBarEquilibration( Counter, Ensemble )

! Use Output_Unit
USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: Output_Unit
USE GlobalVariables

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: Counter      ! Counter
INTEGER( Kind= Int64 ) :: AuxiliarInt1 ! Auxiliar

! *********************************************************************************************** !
! CHARACTER STRINGS                                                                               !
! *********************************************************************************************** !
CHARACTER( LEN= 23 ) :: ProgressBar ! Progress bar
CHARACTER( LEN= 03 ) :: Ensemble    ! Ensemble type
CHARACTER( LEN= 02 ) :: StringSize1 ! String size
CHARACTER( LEN= 08 ) :: StringSize2 ! String size

! Initialization
AuxiliarInt1 = 0

! Progress bar (FORMAT)
IF( ( DBLE( Counter ) / DBLE( nEquilRandom ) ) * 100.D0 < 10.D0 ) THEN
  AuxiliarInt1 = 0
  ProgressBar(1:(20+AuxiliarInt1)) = "Progress("//TRIM( Ensemble )//"): ????%"
  WRITE( Unit= ProgressBar(16:(19+AuxiliarInt1)), Fmt= "(RD,F4.2)" ) ( DBLE( Counter ) / DBLE( nEquilRandom ) ) * 100.D0
ELSE IF( ( DBLE( Counter ) / DBLE( nEquilRandom ) ) * 100.D0 < 100.D0 ) THEN
  AuxiliarInt1 = 1
  ProgressBar(1:(20+AuxiliarInt1)) = "Progress("//TRIM( Ensemble )//"): ?????%"
  WRITE( Unit= ProgressBar(16:(19+AuxiliarInt1)), Fmt= "(RD,F5.2)" ) ( DBLE( Counter ) / DBLE( nEquilRandom ) ) * 100.D0
ELSE
  AuxiliarInt1 = 2
  ProgressBar(1:(20+AuxiliarInt1)) = "Progress("//TRIM( Ensemble )//"): ??????%"
  WRITE( Unit= ProgressBar(16:(19+AuxiliarInt1)), Fmt= "(RD,F6.2)" ) ( DBLE( Counter ) / DBLE( nEquilRandom ) ) * 100.D0
END IF
ProgressBar((21+AuxiliarInt1):23) = REPEAT( " ", ( (2 - AuxiliarInt1) + 1 ) )

! Print progress bar
WRITE( StringSize1, "(I0.2)" ) 20 + AuxiliarInt1 + 1
StringSize2 = "(A1,A"//TRIM( StringSize1 )//")"
WRITE( Unit= Output_Unit, Fmt= StringSize2, Advance= "NO" ) CHAR(13), ProgressBar(1:(20+AuxiliarInt1+1))

! Flush standard output unit
FLUSH( Unit= Output_Unit )

RETURN

END SUBROUTINE ProgressBarEquilibration

! *********************************************************************************************** !
!                                          PROGRESS BAR                                           !
! *********************************************************************************************** !
!       This subroutine generates a progress bar for the algorithm that compress/expand the       !
!              volume of the simulation box to the desirable target packing fraction              !
! *********************************************************************************************** !
SUBROUTINE ProgressBarRandomConfigNPT( iParticle, PressureRandom, DensityRandom )

! Use Output_Unit
USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: Output_Unit
USE GlobalVariables

IMPLICIT NONE

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: iParticle    ! Counter
INTEGER( Kind= Int64 ) :: AuxiliarInt1 ! Auxiliar

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 ) :: PressureRandom, DensityRandom ! Packing fraction

! *********************************************************************************************** !
! CHARACTER STRINGS                                                                               !
! *********************************************************************************************** !
CHARACTER( LEN= 78 ) :: ProgressBar ! Progress bar
CHARACTER( LEN= 02 ) :: StringSize1 ! String size
CHARACTER( LEN= 08 ) :: StringSize2 ! String size

! Initialization
AuxiliarInt1 = 0

! Progress bar (FORMAT)
IF( iParticle < 10 ) THEN
  AuxiliarInt1 = 0
  ProgressBar(1:(12+AuxiliarInt1)) = "Attempts: ? "
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt1)), Fmt= "(I1)" ) iParticle
ELSE IF( iParticle < 100 ) THEN
  AuxiliarInt1 = 1
  ProgressBar(1:(12+AuxiliarInt1)) = "Attempts: ?? "
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt1)), Fmt= "(I2)" ) iParticle
ELSE IF( iParticle < 1000 ) THEN
  AuxiliarInt1 = 2
  ProgressBar(1:(12+AuxiliarInt1)) = "Attempts: ??? "
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt1)), Fmt= "(I3)" ) iParticle
ELSE IF( iParticle < 10000 ) THEN
  AuxiliarInt1 = 3
  ProgressBar(1:(12+AuxiliarInt1)) = "Attempts: ???? "
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt1)), Fmt= "(I4)" ) iParticle
ELSE IF( iParticle < 100000 ) THEN
  AuxiliarInt1 = 4
  ProgressBar(1:(12+AuxiliarInt1)) = "Attempts: ????? "
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt1)), Fmt= "(I5)" ) iParticle
ELSE IF( iParticle < 1000000 ) THEN
  AuxiliarInt1 = 5
  ProgressBar(1:(12+AuxiliarInt1)) = "Attempts: ?????? "
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt1)), Fmt= "(I6)" ) iParticle
ELSE IF( iParticle < 10000000 ) THEN
  AuxiliarInt1 = 6
  ProgressBar(1:(12+AuxiliarInt1)) = "Attempts: ??????? "
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt1)), Fmt= "(I7)" ) iParticle
ELSE IF( iParticle < 100000000 ) THEN
  AuxiliarInt1 = 7
  ProgressBar(1:(12+AuxiliarInt1)) = "Attempts: ???????? "
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt1)), Fmt= "(I8)" ) iParticle
ELSE IF( iParticle < 1000000000 ) THEN
  AuxiliarInt1 = 8
  ProgressBar(1:(12+AuxiliarInt1)) = "Attempts: ????????? "
  WRITE( Unit= ProgressBar(11:(11+AuxiliarInt1)), Fmt= "(I9)" ) iParticle
ELSE IF( iParticle >= 1000000000 ) THEN
  AuxiliarInt1 = 10
  ProgressBar(1:(12+AuxiliarInt1)) = "Attempts: > 1 billion "
END IF
ProgressBar((13+AuxiliarInt1):(36+AuxiliarInt1)) = "| Density: ???????????? "
WRITE( Unit= ProgressBar((24+AuxiliarInt1):(35+AuxiliarInt1)), Fmt= "(E12.6)" ) DensityRandom
ProgressBar((37+AuxiliarInt1):(58+AuxiliarInt1)) = "(TARGET: ????????????)"
WRITE( Unit= ProgressBar((46+AuxiliarInt1):(57+AuxiliarInt1)), Fmt= "(E12.6)" ) PressureRandom
ProgressBar((59+AuxiliarInt1):78) = REPEAT( " ", ( (10 - AuxiliarInt1) + 1 ) )

! Print progress bar
WRITE( StringSize1, "(I0.2)" ) 58 + AuxiliarInt1 + 1
StringSize2 = "(A1,A"//TRIM( StringSize1 )//")"
WRITE( Unit= Output_Unit, Fmt= StringSize2, Advance= "NO" ) CHAR(13), ProgressBar(1:(58+AuxiliarInt1+1))

! Flush standard output unit
FLUSH( Unit= Output_Unit )

RETURN

END SUBROUTINE ProgressBarRandomConfigNPT