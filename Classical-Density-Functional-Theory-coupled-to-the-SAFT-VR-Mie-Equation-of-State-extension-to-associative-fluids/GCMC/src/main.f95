! ----------------------------------------------------------------------- !
! Main Program                                                            !
! ----------------------------------------------------------------------- !
PROGRAM MonteCarlo_Mie

USE GlobalVariables
USE QuaternionOperations
USE VectorOperations

IMPLICIT NONE

! INTEGER VARIABLES -*- THIS IS SINGLE PRECISION -*-
INTEGER :: SeedSize ! Seed array size

! INTEGER VARIABLES (ALLOCATABLE) -*- THIS IS SINGLE PRECISION -*-
INTEGER, DIMENSION( : ), ALLOCATABLE :: Seed ! Random seed

INTEGER( Kind= Int64 ) :: iParticle, rParticle, iCycle, iSite, iSlab, bin, iOldParticles, counterWidom
INTEGER( Kind= Int64 ) :: index1, index2, nValidConfigs, iDestroyedParticle, iCreatedParticle
INTEGER( Kind= Int64 ) :: nAcceptanceTranslation, nAcceptanceRotation, nAcceptanceCreation, nAcceptanceDestruction
INTEGER( Kind= Int64 ) :: nMovementTranslationCounter, nMovementRotationCounter, nCreationCounter, nDestructionCounter

INTEGER( Kind= Int64 ), ALLOCATABLE :: AssociatedSitesSave( :, :, : ), AssociatedSitesTemp( :, :, : ) ! Indexes of the associated site

REAL( Kind= Real64 ) :: RotatedPosition(3), ScalingDistanceUnitBox(3), binWidth
REAL( Kind= Real64 ) :: iOldPosition(3), iOldQuaternion(0:3), iOldOrientation(3), iOldSitePosition(3,4), iOldSiteOrientation(3,4), &
&                       iOldSiteQuaternion(0:3,4), iOldPotential, iOldVirial, iDestroyedPosition(3), iDestroyedQuaternion(0:3), &
&                       iDestroyedOrientation(3), iDestroyedSitePosition(3,4), iDestroyedSiteOrientation(3,4), &
&                       iDestroyedSiteQuaternion(0:3,4), iDestroyedPotential, iDestroyedVirial, iDestroyedAssociationTerm, &
&                       iDestroyedMieSWTerm, iDestroyedSteeleTerm
REAL( Kind= Real64 ) :: iNewPosition(3), iNewQuaternion(0:3), iNewOrientation(3), iNewSitePosition(3,4), iNewSiteOrientation(3,4), &
&                       iNewSiteQuaternion(0:3,4), iNewPotential, iNewVirial, iCreatedPosition(3), iCreatedQuaternion(0:3), &
&                       iCreatedOrientation(3), iCreatedSitePosition(3,4), iCreatedSiteOrientation(3,4), &
&                       iCreatedSiteQuaternion(0:3,4), iCreatedPotential, iCreatedVirial, iCreatedAssociationTerm, &
&                       iCreatedMieSWTerm, iCreatedSteeleTerm
REAL( Kind= Real64 ) :: PotentialEnergyDifference, VirialDifference, Ratio, AssociationEnergyDifference, SteeleEnergyDifference 
REAL( Kind= Real64 ) :: tPotentialEnergy, VirialContribution, tAssociationTerm, tMieSWTerm, tSteeleTerm, iOldAssociationTerm
REAL( Kind= Real64 ) :: pPotential, iOldMieSWTerm, iOldSteeleTerm, iNewAssociationTerm, iNewMieSWTerm, iNewSteeleTerm
REAL( Kind= Real64 ) :: MieEnergyDifference, GrandCanonicalEnergyDifference, WidomEnergy, ExcessWidomChemicalPotential

REAL( Kind= Real64 ), ALLOCATABLE :: avgHistogram(:)
REAL( Kind= Real64 ), ALLOCATABLE :: Histogram(:)
REAL( Kind= Real64 ), ALLOCATABLE :: pPositionTemp(:,:), pQuaternionTemp(:,:), pOrientationTemp(:,:)
REAL( Kind= Real64 ), ALLOCATABLE :: sitePositionTemp(:,:,:), siteOrientationTemp(:,:,:), siteQuaternionTemp(:,:,:)

CHARACTER( LEN= 4 )   :: EnsembleType
CHARACTER( LEN= 400 ) :: Lattice
CHARACTER( LEN= 140 ) :: DescriptorString ! Descriptor for strings

LOGICAL :: MovementTranslationLogical, MovementRotationLogical, MovementCreationLogical, MovementDestructionLogical, Overlap

LOGICAL, ALLOCATABLE :: isSiteAssociatedSave( :, : ), isSiteAssociatedTemp( :, : ) ! Checks if a site is associated

 
! Get the maximum number of threads
#ifdef _OPENMP
  nThreads = OMP_GET_MAX_THREADS(  )
#else
  nThreads = 1
#endif

WRITE( *, "(A,I0)" ) "Number of threads: ", nThreads

! Initialize variables
CALL Variable_Initialization(  )

! Construct sites
CALL Site_Initialization(  )

! Initialization
EnsembleType = "NVT"
IF( MovementProbability < 1.D0 ) THEN
  EnsembleType = "μVT"
END IF

! Histogram
ALLOCATE( Histogram(nSlabs), avgHistogram(nSlabs) )
nValidConfigs = 0
Histogram = 0.D0
avgHistogram = 0.D0

! Bin width
binWidth = poreWidth / DBLE( nSlabs )

! Fixed seed
IF( .NOT. RandomSeed ) THEN
  CALL RANDOM_SEED( Size= SeedSize )
  ALLOCATE( Seed(SeedSize) )
  Seed = 123456789
  CALL RANDOM_SEED( Put= Seed )
END IF

! Generate random configuration
IF( RandomConf ) THEN
  CALL RandomConfiguration(  )
END IF

! Initial configuration
IF( .NOT. RandomConf ) THEN
  OPEN( Unit= 15, File= "initial_configuration.xyz", Action= "READ" )
  READ( 15, * ) Dummy
  READ( 15, "(A)" ) Lattice
  index1 = INDEX( TRIM( Lattice ), '"' )
  index2 = INDEX( TRIM( Lattice(index1+1:) ), '"' )
  READ( Lattice(index1+1:index2+index1-1), * ) BoxLength
  DO iParticle = 1, nParticles
    ! Positions and quaternions
    READ( 15, * ) Dummy, Dummy, pPosition(1:3,iParticle), pQuaternion(1:3,iParticle), pQuaternion(0,iParticle)
    ! Generate orientation
    CALL VectorRotation( zAxis, pQuaternion(:,iParticle), pOrientation(:,iParticle) )
    ! Sites
    DO iSite = 1, 4
      READ( 15, * ) Dummy, Dummy, sitePosition(1:3,iSite,iParticle), siteQuaternion(1:3,iSite,iParticle), &
      &             siteQuaternion(0,iSite,iParticle)
      ! Generate orientation
      CALL VectorRotation( ReferencePositions(:,iSite), pQuaternion(:,iParticle), RotatedPosition(:) )
      siteOrientation(:,iSite,iParticle) = RotatedPosition(:) / SQRT( DOT_PRODUCT( RotatedPosition(:), RotatedPosition(:) ) )
    END DO
  END DO
  CLOSE( Unit= 15 )
END IF

! Check associated sites
WRITE( *, "(A)", Advance= "NO" ) "Now checking for associated sites... "

! Attribute sites
CALL Attribute_Sites( BoxLength, BoxLengthInverse )

! Number of associated sites
nAssocSites = COUNT( isSiteAssociated ) / 2 ! Removes redundancy

! Fraction of associated sites
fractionAssociatedSites = DBLE( nAssocSites ) / DBLE( tAssocSites )
fractionFreeSites = 1.D0 - fractionAssociatedSites

! Status
WRITE( *, "(A,A,G0.6,A,G0,A,G0,A)", Advance= "YES" ) CHAR(13), "Now checking associated sites. Fraction of free sites: ", &
&                                                    fractionFreeSites, " [", nAssocSites, " associated site(s) in ", &
&                                                    tAssocSites, " site(s)]"
WRITE( *, "(A)" ) " "

! Initial total energy calculation
CALL Calculation_Total_Energy( BoxLength, BoxLengthInverse, tPotentialEnergy, VirialContribution, tAssociationTerm, tMieSWTerm, &
&                              tSteeleTerm, .FALSE., Overlap )
IF( Overlap ) THEN
  WRITE( *, "(A)" ) "Initial configuration has overlapping particles. Exiting..."
  STOP
END IF
WRITE( *, "(2G0)" ) "Total potential energy: ", tPotentialEnergy
WRITE( *, "(2G0)" ) "Total Mie/SW energy: ", tMieSWTerm
WRITE( *, "(2G0)" ) "Total association energy: ", tAssociationTerm
WRITE( *, "(2G0)" ) "Total Steele energy: ", tSteeleTerm
WRITE( *, "(G0)" ) " "

! Potential energy per particle
pPotential = tPotentialEnergy / DBLE( nParticles )
WRITE( *, "(2G0)" ) 'Potential energy per particle: ', pPotential
WRITE( *, "(2G0)" ) 'Mie/SW energy per particle: ', tMieSWTerm / DBLE( nParticles )
WRITE( *, "(2G0)" ) 'Association energy per particle: ', tAssociationTerm / DBLE( nParticles )
WRITE( *, "(2G0)" ) 'Steele energy per particle: ', tSteeleTerm / DBLE( nParticles )
WRITE( *, "(G0)" ) " "

! File for equilibration data
OPEN( Unit= 20, File= "ratio_translation.dat" )
OPEN( Unit= 25, File= "ratio_creation.dat" )
OPEN( Unit= 30, File= "ratio_rotation.dat" )
OPEN( Unit= 35, File= "ratio_destruction.dat" )
OPEN( Unit= 40, File= "potential.dat" )
OPEN( Unit= 50, File= "trajectory.xyz" )
OPEN( Unit= 55, File= "fraction_free_sites.dat" )
OPEN( Unit= 60, File= "density.dat" )
IF( isWidomEnabled ) OPEN( Unit= 65, File= "widom.dat" )

! Initialization
nMovementTranslationCounter = 0
nMovementRotationCounter = 0
nCreationCounter = 0
nDestructionCounter = 0
nAcceptanceTranslation = 0
nAcceptanceRotation = 0
nAcceptanceCreation = 0
nAcceptanceDestruction = 0
counterWidom = 0
WidomEnergy = 0.D0

! Allocation
ALLOCATE( AssociatedSitesSave( 2, 4, nParticles ), isSiteAssociatedSave( 4, nParticles ) )

! Status
CALL Progress_Bar( 1_INT64, nCycles, EnsembleType )

! Monte Carlo algorithm
DO iCycle = 1, nCycles

  ! Progress bar
  IF( MOD( iCycle, nSave ) == 0 ) THEN
    CALL Progress_Bar( iCycle, nCycles, EnsembleType )
  END IF

  ! Random sampling (molecular movement or particle creation/destruction)
  CALL Random_Number( RandomNumber )

  ! Perform particle translation/rotation
  IF( RandomNumber < MovementProbability ) THEN

    ! Disable particle creation/destruction
    MovementCreationLogical = .FALSE.
    MovementDestructionLogical = .FALSE.

    ! Random particle
    CALL Random_Number( RandomNumber )
    rParticle = INT( RandomNumber * DBLE( nParticles ) ) + 1

    ! Assignment of previous configuration
    iOldPosition(:)          = pPosition(:,rParticle)         ! Position
    iOldQuaternion(:)        = pQuaternion(:,rParticle)       ! Quaternion
    iOldOrientation(:)       = pOrientation(:,rParticle)      ! Orientation
    iOldSitePosition(:,:)    = sitePosition(:,:,rParticle)    ! Position (sites)
    iOldSiteOrientation(:,:) = siteOrientation(:,:,rParticle) ! Orientation (sites)
    iOldSiteQuaternion(:,:)  = siteQuaternion(:,:,rParticle)  ! Quaternion (sites)
    AssociatedSitesSave      = AssociatedSites
    isSiteAssociatedSave     = isSiteAssociated
    
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

    ! Translation move
    IF( MovementTranslationLogical ) THEN
      ! Random translation along x-, y-, and z-axis
      IF( .NOT. xyTranslationIndependency ) THEN
        CALL Random_Number( RandomNumber )
        iNewPosition(1) = iOldPosition(1) + ( ( 2.D0 * RandomNumber ) - 1.D0 ) * maxTranslationalDisplc ! Range [-drmax,drmax]
        ! Random translation along y-axis
        CALL Random_Number( RandomNumber )
        iNewPosition(2) = iOldPosition(2) + ( ( 2.D0 * RandomNumber ) - 1.D0 ) * maxTranslationalDisplc ! Range [-drmax,drmax]
        ! Random translation along z-axis
        CALL Random_Number( RandomNumber )
        iNewPosition(3) = iOldPosition(3) + ( ( 2.D0 * RandomNumber ) - 1.D0 ) * maxTranslationalDisplc ! Range [-drmax,drmax]
      ! Random translation along the XY-plane or along z-axis
      ELSE
        CALL Random_Number( RandomNumber )
        ! Random translation along the XY-plane
        IF( RandomNumber < TranslationalProbabilityXY ) THEN
          CALL Random_Number( RandomNumber )
          iNewPosition(1) = iOldPosition(1) + ( ( 2.D0 * RandomNumber ) - 1.D0 ) * maxTranslationalDisplc ! Range [-drmax,drmax]
          ! Random translation along y-axis
          CALL Random_Number( RandomNumber )
          iNewPosition(2) = iOldPosition(2) + ( ( 2.D0 * RandomNumber ) - 1.D0 ) * maxTranslationalDisplc ! Range [-drmax,drmax]
          ! Z is the same
          iNewPosition(3) = iOldPosition(3)
        ! Random translation along Z
        ELSE
          ! X and Y are the same
          iNewPosition(1) = iOldPosition(1)
          iNewPosition(2) = iOldPosition(2)
          ! Random translation along z-axis
          CALL Random_Number( RandomNumber )
          iNewPosition(3) = iOldPosition(3) + ( ( 2.D0 * RandomNumber ) - 1.D0 ) * maxTranslationalDisplc ! Range [-drmax,drmax]
        END IF
      END IF
      ! Apply periodic boundary conditions
      CALL MatrixVectorMultiplication( BoxLengthInverse, iNewPosition, ScalingDistanceUnitBox )
      IF( Confinement ) THEN
        ScalingDistanceUnitBox(1:2) = ScalingDistanceUnitBox(1:2) - ANINT( ScalingDistanceUnitBox(1:2) ) ! Periodic boundary conditions (system is confined along z-axis)
      ELSE
        ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
      END IF
      CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, iNewPosition )
    ELSE IF( .NOT. MovementTranslationLogical ) THEN
      iNewPosition(:) = iOldPosition(:)
    END IF

    ! Rotation move
    IF( MovementRotationLogical ) THEN
      ! Random Composed Unit Quaternion
      CALL QuaternionCombination( iOldQuaternion, iNewQuaternion, maxRotationalDisplc )
      ! Active transformation (rotation)
      CALL VectorRotation( zAxis, iNewQuaternion, iNewOrientation )
    ELSE IF( .NOT. MovementRotationLogical ) THEN
      iNewQuaternion(:) = iOldQuaternion(:)
      iNewOrientation(:) = iOldOrientation(:)
    END IF

    ! Site positions and orientations (after translation or rotation)
    DO iSite = 1, 4
      CALL VectorRotation( ReferencePositions(:,iSite ), iNewQuaternion, RotatedPosition(:) )
      iNewSitePosition(:,iSite) = iNewPosition(:) + RotatedPosition(:)
      iNewSiteOrientation(:,iSite) = RotatedPosition(:) / SQRT( DOT_PRODUCT( RotatedPosition(:), RotatedPosition(:) ) )
      CALL QuaternionMultiplication( iNewQuaternion, ReferenceQuaternions(:,iSite), iNewSiteQuaternion(:,iSite) )
    END DO

    ! Calculate old energy
    CALL Particle_Energy( rParticle, iOldPosition, iOldSitePosition, iOldPotential, iOldVirial, iOldAssociationTerm, &
    &                     iOldMieSWTerm, iOldSteeleTerm, BoxLength, BoxLengthInverse, .FALSE., .FALSE., .FALSE., .FALSE., Overlap ) ! LOGICAL: RandomConfig, Displaced, Destroyed, Created
    IF( Overlap ) STOP "ERROR: Overlap detected in old configuration during displacement."

    ! Calculate new energy
    CALL Particle_Energy( rParticle, iNewPosition, iNewSitePosition, iNewPotential, iNewVirial, iNewAssociationTerm, &
    &                     iNewMieSWTerm, iNewSteeleTerm, BoxLength, BoxLengthInverse, .FALSE., .TRUE., .FALSE., .FALSE., Overlap ) ! LOGICAL: RandomConfig, Displaced, Destroyed, Created

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
          pPosition(:,rParticle) = iNewPosition(:)
          ! Update quaternion
          pQuaternion(:,rParticle) = iNewQuaternion(:)
          ! Update orientation
          pOrientation(:,rParticle) = iNewOrientation(:)
          ! Update site positions and orientations
          sitePosition(:,:,rParticle)    = iNewSitePosition(:,:)
          siteOrientation(:,:,rParticle) = iNewSiteOrientation(:,:)
          siteQuaternion(:,:,rParticle)  = iNewSiteQuaternion(:,:)
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
            pPosition(:,rParticle) = iNewPosition(:)
            ! Update quaternion
            pQuaternion(:,rParticle) = iNewQuaternion(:)
            ! Update orientation
            pOrientation(:,rParticle) = iNewOrientation(:)
            ! Update site positions and orientations
            sitePosition(:,:,rParticle)    = iNewSitePosition(:,:)
            siteOrientation(:,:,rParticle) = iNewSiteOrientation(:,:)
            siteQuaternion(:,:,rParticle)  = iNewSiteQuaternion(:,:)
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
            ! Recover associations
            isSiteAssociated = isSiteAssociatedSave
            AssociatedSites = AssociatedSitesSave
          END IF
        END IF
      ELSE
        ! Recover associations
        isSiteAssociated = isSiteAssociatedSave
        AssociatedSites = AssociatedSitesSave
      END IF
    ELSE
      ! Recover associations
      isSiteAssociated = isSiteAssociatedSave
      AssociatedSites = AssociatedSitesSave
    END IF

  ! Particle creation/destruction
  ELSE IF( RandomNumber >= MovementProbability ) THEN

    ! Disable translation and rotation
    MovementTranslationLogical = .FALSE.
    MovementRotationLogical    = .FALSE.

    ! Assignment of previous configuration
    iOldParticles = nParticles

    ! Creation/destruction type
    CALL Random_Number( RandomNumber )

    ! Creation
    IF( RandomNumber < CreationProbability ) THEN

      ! Update logicals
      MovementCreationLogical = .TRUE.
      MovementDestructionLogical = .FALSE. ! Disable destruction
      nCreationCounter = nCreationCounter + 1

      ! New particle
      iCreatedParticle = nParticles + 1

      ! Old properties
      isSiteAssociatedSave = isSiteAssociated
      AssociatedSitesSave = AssociatedSites

      ! Create random position within unit box
      CALL Random_Number( RandomNumber )
      ScalingDistanceUnitBox(1) = RandomNumber - 0.5D0 ! X coordinate
      CALL Random_Number( RandomNumber )
      ScalingDistanceUnitBox(2) = RandomNumber - 0.5D0 ! Y coordinate
      CALL Random_Number( RandomNumber )
      ScalingDistanceUnitBox(3) = RandomNumber - 0.5D0 ! Z coordinate

      ! Put new particle into the simulation box
      CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, iCreatedPosition )

      ! Create random quaternion
      CALL RandomQuaternionGenerator( iCreatedQuaternion, cPi ) ! Angle in [-π, π]

      ! Create random orientation
      CALL VectorRotation( zAxis, iCreatedQuaternion, iCreatedOrientation )

      ! Random sites
      DO iSite = 1, 4
        CALL VectorRotation( ReferencePositions(:,iSite), iCreatedQuaternion, RotatedPosition(:) )
        ! Random positions
        iCreatedSitePosition(:,iSite) = iCreatedPosition(:) + RotatedPosition(:)
        ! Random orientations
        iCreatedSiteOrientation(:,iSite) = RotatedPosition(:) / SQRT( DOT_PRODUCT( RotatedPosition(:), RotatedPosition(:) ) )
        ! Random quaternions
        CALL QuaternionMultiplication( iCreatedQuaternion, ReferenceQuaternions(:,iSite), iCreatedSiteQuaternion(:,iSite) )
      END DO

      ! New number of particles
      nParticles = nParticles + 1

      ! Temporary arrays
      pPositionTemp = pPosition
      pQuaternionTemp = pQuaternion
      pOrientationTemp = pOrientation
      sitePositionTemp = sitePosition
      siteOrientationTemp = siteOrientation
      siteQuaternionTemp = siteQuaternion
      isSiteAssociatedTemp = isSiteAssociated
      AssociatedSitesTemp = AssociatedSites

      ! Deallocate and reallocate arrays
      DEALLOCATE( pPosition, pQuaternion, pOrientation, sitePosition, siteOrientation, siteQuaternion, &
      &           isSiteAssociated, AssociatedSites )
      ALLOCATE( pPosition(3,nParticles), pQuaternion(0:3,nParticles), pOrientation(3,nParticles), &
      &         sitePosition(3,4,nParticles), siteOrientation(3,4,nParticles), siteQuaternion(0:3,4,nParticles), &
      &         isSiteAssociated(4,nParticles), AssociatedSites(2,4,nParticles) )

      ! Associate arrays
      pPosition(:,1:iOldParticles) = pPositionTemp(:,1:iOldParticles)
      pPosition(:,iCreatedParticle) = iCreatedPosition
      pOrientation(:,1:iOldParticles) = pOrientationTemp(:,1:iOldParticles)
      pOrientation(:,iCreatedParticle) = iCreatedOrientation
      pQuaternion(:,1:iOldParticles) = pQuaternionTemp(:,1:iOldParticles)
      pQuaternion(:,iCreatedParticle) = iCreatedQuaternion
      sitePosition(:,:,1:iOldParticles) = sitePositionTemp(:,:,1:iOldParticles)
      sitePosition(:,:,iCreatedParticle) = iCreatedSitePosition
      siteOrientation(:,:,1:iOldParticles) = siteOrientationTemp(:,:,1:iOldParticles)
      siteOrientation(:,:,iCreatedParticle) = iCreatedSiteOrientation
      siteQuaternion(:,:,1:iOldParticles) = siteQuaternionTemp(:,:,1:iOldParticles)
      siteQuaternion(:,:,iCreatedParticle) = iCreatedSiteQuaternion
      isSiteAssociated(:,1:iOldParticles) = isSiteAssociatedTemp(:,1:iOldParticles)
      isSiteAssociated(:,iCreatedParticle) = .FALSE.
      AssociatedSites(:,:,1:iOldParticles) = AssociatedSitesTemp(:,:,1:iOldParticles)
      AssociatedSites(:,:,iCreatedParticle) = 0

      ! Calculate new energy of inserted particle
      CALL Particle_Energy( iCreatedParticle, iCreatedPosition, iCreatedSitePosition, iCreatedPotential, &
      &                     iCreatedVirial, iCreatedAssociationTerm, iCreatedMieSWTerm, iCreatedSteeleTerm, BoxLength, &
      &                     BoxLengthInverse, .FALSE., .FALSE., .FALSE., .TRUE., Overlap ) ! LOGICAL: RandomConfig, Displaced, Destroyed, Created

      ! Energy differences
      PotentialEnergyDifference = iCreatedPotential
      AssociationEnergyDifference = iCreatedAssociationTerm
      MieEnergyDifference = iCreatedMieSWTerm
      SteeleEnergyDifference = iCreatedSteeleTerm

      ! Virial differences
      VirialDifference = iCreatedVirial

      ! Compute the grand-canonical energy difference
      GrandCanonicalEnergyDifference = PotentialEnergyDifference - ReducedChemicalPotential - GrandCanonicalConstant + &
      &                                LOG( DBLE( iOldParticles + 1 ) )

    ! Destruction
    ELSE IF( RandomNumber >= CreationProbability ) THEN

      ! Update logicals
      MovementDestructionLogical = .TRUE.
      MovementCreationLogical = .FALSE. ! Disable creation
      nDestructionCounter = nDestructionCounter + 1

      ! Random particle to be destroyed
      CALL Random_Number( RandomNumber )
      iDestroyedParticle = INT( RandomNumber * DBLE( iOldParticles ) ) + 1

      ! Assignment of previous configuration
      iDestroyedPosition(:)          = pPosition(:,iDestroyedParticle)         ! Position
      iDestroyedQuaternion(:)        = pQuaternion(:,iDestroyedParticle)       ! Quaternion
      iDestroyedOrientation(:)       = pOrientation(:,iDestroyedParticle)      ! Orientation
      iDestroyedSitePosition(:,:)    = sitePosition(:,:,iDestroyedParticle)    ! Position (sites)
      iDestroyedSiteOrientation(:,:) = siteOrientation(:,:,iDestroyedParticle) ! Orientation (sites)
      iDestroyedSiteQuaternion(:,:)  = siteQuaternion(:,:,iDestroyedParticle)  ! Quaternion (sites)
      isSiteAssociatedSave           = isSiteAssociated
      AssociatedSitesSave            = AssociatedSites

      ! Calculate old energy of destroyed particle
      CALL Particle_Energy( iDestroyedParticle, iDestroyedPosition, iDestroyedSitePosition, iDestroyedPotential, &
      &                     iDestroyedVirial, iDestroyedAssociationTerm, iDestroyedMieSWTerm, iDestroyedSteeleTerm, BoxLength, &
      &                     BoxLengthInverse, .FALSE., .FALSE., .TRUE., .FALSE., Overlap ) ! LOGICAL: RandomConfig, Displaced, Destroyed, Created
      IF( Overlap ) STOP "ERROR: Overlap detected in old configuration during displacement."

      ! New number of particles
      nParticles = nParticles - 1

      ! Temporary arrays
      pPositionTemp = pPosition
      pQuaternionTemp = pQuaternion
      pOrientationTemp = pOrientation
      sitePositionTemp = sitePosition
      siteOrientationTemp = siteOrientation
      siteQuaternionTemp = siteQuaternion
      AssociatedSitesTemp = AssociatedSites
      isSiteAssociatedTemp = isSiteAssociated

      ! Deallocate and reallocate arrays
      DEALLOCATE( pPosition, pQuaternion, pOrientation, sitePosition, siteOrientation, siteQuaternion, &
      &           AssociatedSites, isSiteAssociated )
      ALLOCATE( pPosition(3,nParticles), pQuaternion(0:3,nParticles), pOrientation(3,nParticles), &
      &         sitePosition(3,4,nParticles), siteOrientation(3,4,nParticles), siteQuaternion(0:3,4,nParticles), &
      &         isSiteAssociated(4,nParticles), AssociatedSites(2,4,nParticles) )

      ! Associate arrays
      IF( iDestroyedParticle == 1 ) THEN
        pPosition(:,iDestroyedParticle:nParticles) = pPositionTemp(:,iDestroyedParticle+1:iOldParticles)
        pOrientation(:,iDestroyedParticle:nParticles) = pOrientationTemp(:,iDestroyedParticle+1:iOldParticles)
        pQuaternion(:,iDestroyedParticle:nParticles) = pQuaternionTemp(:,iDestroyedParticle+1:iOldParticles)
        sitePosition(:,:,iDestroyedParticle:nParticles) = sitePositionTemp(:,:,iDestroyedParticle+1:iOldParticles)
        siteOrientation(:,:,iDestroyedParticle:nParticles) = siteOrientationTemp(:,:,iDestroyedParticle+1:iOldParticles)
        siteQuaternion(:,:,iDestroyedParticle:nParticles) = siteQuaternionTemp(:,:,iDestroyedParticle+1:iOldParticles)
        isSiteAssociated(:,iDestroyedParticle:nParticles) = isSiteAssociatedTemp(:,iDestroyedParticle+1:iOldParticles)
        IF( ANY( AssociatedSitesTemp(2,:,:) == iDestroyedParticle ) ) STOP "ERROR: Associated site found after destruction!"
        AssociatedSites(:,:,iDestroyedParticle:nParticles) = AssociatedSitesTemp(:,:,iDestroyedParticle+1:iOldParticles)
        WHERE (AssociatedSites(2,:,:) > iDestroyedParticle)
          AssociatedSites(2,:,:) = AssociatedSites(2,:,:) - 1
        END WHERE
      ELSE IF( iDestroyedParticle == iOldParticles ) THEN
        pPosition(:,1:iDestroyedParticle-1) = pPositionTemp(:,1:iDestroyedParticle-1)
        pOrientation(:,1:iDestroyedParticle-1) = pOrientationTemp(:,1:iDestroyedParticle-1)
        pQuaternion(:,1:iDestroyedParticle-1) = pQuaternionTemp(:,1:iDestroyedParticle-1)
        sitePosition(:,:,1:iDestroyedParticle-1) = sitePositionTemp(:,:,1:iDestroyedParticle-1)
        siteOrientation(:,:,1:iDestroyedParticle-1) = siteOrientationTemp(:,:,1:iDestroyedParticle-1)
        siteQuaternion(:,:,1:iDestroyedParticle-1) = siteQuaternionTemp(:,:,1:iDestroyedParticle-1)
        isSiteAssociated(:,1:iDestroyedParticle-1) = isSiteAssociatedTemp(:,1:iDestroyedParticle-1)
        IF( ANY( AssociatedSitesTemp(2,:,:) == iDestroyedParticle ) ) STOP "ERROR: Associated site found after destruction!"
        AssociatedSites(:,:,1:iDestroyedParticle-1) = AssociatedSitesTemp(:,:,1:iDestroyedParticle-1)
      ELSE
        pPosition(:,1:iDestroyedParticle-1) = pPositionTemp(:,1:iDestroyedParticle-1)
        pPosition(:,iDestroyedParticle:nParticles) = pPositionTemp(:,iDestroyedParticle+1:iOldParticles)
        pOrientation(:,1:iDestroyedParticle-1) = pOrientationTemp(:,1:iDestroyedParticle-1)
        pOrientation(:,iDestroyedParticle:nParticles) = pOrientationTemp(:,iDestroyedParticle+1:iOldParticles)
        pQuaternion(:,1:iDestroyedParticle-1) = pQuaternionTemp(:,1:iDestroyedParticle-1)
        pQuaternion(:,iDestroyedParticle:nParticles) = pQuaternionTemp(:,iDestroyedParticle+1:iOldParticles)
        sitePosition(:,:,1:iDestroyedParticle-1) = sitePositionTemp(:,:,1:iDestroyedParticle-1)
        sitePosition(:,:,iDestroyedParticle:nParticles) = sitePositionTemp(:,:,iDestroyedParticle+1:iOldParticles)
        siteOrientation(:,:,1:iDestroyedParticle-1) = siteOrientationTemp(:,:,1:iDestroyedParticle-1)
        siteOrientation(:,:,iDestroyedParticle:nParticles) = siteOrientationTemp(:,:,iDestroyedParticle+1:iOldParticles)
        siteQuaternion(:,:,1:iDestroyedParticle-1) = siteQuaternionTemp(:,:,1:iDestroyedParticle-1)
        siteQuaternion(:,:,iDestroyedParticle:nParticles) = siteQuaternionTemp(:,:,iDestroyedParticle+1:iOldParticles)
        isSiteAssociated(:,1:iDestroyedParticle-1) = isSiteAssociatedTemp(:,1:iDestroyedParticle-1)
        isSiteAssociated(:,iDestroyedParticle:nParticles) = isSiteAssociatedTemp(:,iDestroyedParticle+1:iOldParticles)
        IF( ANY( AssociatedSitesTemp(2,:,:) == iDestroyedParticle ) ) STOP "ERROR: Associated site found after destruction!"
        AssociatedSites(:,:,1:iDestroyedParticle-1) = AssociatedSitesTemp(:,:,1:iDestroyedParticle-1)
        AssociatedSites(:,:,iDestroyedParticle:nParticles) = AssociatedSitesTemp(:,:,iDestroyedParticle+1:iOldParticles)
        WHERE (AssociatedSites(2,:,:) > iDestroyedParticle)
          AssociatedSites(2,:,:) = AssociatedSites(2,:,:) - 1
        END WHERE
      END IF

      ! Energy differences
      PotentialEnergyDifference = - iDestroyedPotential
      AssociationEnergyDifference = - iDestroyedAssociationTerm
      MieEnergyDifference = - iDestroyedMieSWTerm
      SteeleEnergyDifference = - iDestroyedSteeleTerm

      ! Virial differences
      VirialDifference = - iDestroyedVirial

      ! Compute the grand-canonical energy difference
      GrandCanonicalEnergyDifference = PotentialEnergyDifference + ReducedChemicalPotential + GrandCanonicalConstant - &
      &                                LOG( DBLE( iOldParticles ) )

    END IF

    ! Metropolis criterion
    IF( .NOT. Overlap ) THEN
      IF( GrandCanonicalEnergyDifference <= 75.D0 ) THEN
        IF( GrandCanonicalEnergyDifference <= 0.D0 ) THEN
          ! Update total energy
          tPotentialEnergy = tPotentialEnergy + PotentialEnergyDifference
          tAssociationTerm = tAssociationTerm + AssociationEnergyDifference
          tMieSWTerm = tMieSWTerm + MieEnergyDifference
          tSteeleTerm = tSteeleTerm + SteeleEnergyDifference
          ! Update virial
          VirialContribution = VirialContribution + VirialDifference
          ! Update displacement counter
          IF( MovementCreationLogical ) THEN
            nAcceptanceCreation = nAcceptanceCreation + 1 ! Creation move counter
          ELSE IF( MovementDestructionLogical ) THEN
            nAcceptanceDestruction = nAcceptanceDestruction + 1 ! Destruction move counter
          END IF
          ! Update total number of association sites
          tAssocSites = 4 * nParticles / 2
        ELSE
          CALL Random_Number( RandomNumber ) 
          IF( EXP( - GrandCanonicalEnergyDifference ) >= RandomNumber ) THEN
            ! Update total energy
            tPotentialEnergy = tPotentialEnergy + PotentialEnergyDifference
            tAssociationTerm = tAssociationTerm + AssociationEnergyDifference
            tMieSWTerm = tMieSWTerm + MieEnergyDifference
            tSteeleTerm = tSteeleTerm + SteeleEnergyDifference
            ! Update virial
            VirialContribution = VirialContribution + VirialDifference
            ! Update displacement counter
            IF( MovementCreationLogical ) THEN
              nAcceptanceCreation = nAcceptanceCreation + 1 ! Translational move counter
            ELSE IF( MovementDestructionLogical ) THEN
              nAcceptanceDestruction = nAcceptanceDestruction + 1 ! Rotational move counter
            END IF
            ! Update total number of association sites
            tAssocSites = 4 * nParticles / 2
          ELSE
            ! Recover number of particles
            nParticles = iOldParticles
            ! Recover arrays
            DEALLOCATE( pPosition, pQuaternion, pOrientation, sitePosition, siteOrientation, siteQuaternion, &
            &           isSiteAssociated, AssociatedSites )
            pPosition = pPositionTemp
            pQuaternion = pQuaternionTemp
            pOrientation = pOrientationTemp
            sitePosition = sitePositionTemp
            siteOrientation = siteOrientationTemp
            siteQuaternion = siteQuaternionTemp
            isSiteAssociated = isSiteAssociatedSave
            AssociatedSites = AssociatedSitesSave
            ! Deallocate temporary arrays
            DEALLOCATE( pPositionTemp, pQuaternionTemp, pOrientationTemp, sitePositionTemp, siteOrientationTemp, &
            &           siteQuaternionTemp, isSiteAssociatedTemp, AssociatedSitesTemp )
          END IF
        END IF
      ELSE
        ! Recover number of particles
        nParticles = iOldParticles
        ! Recover arrays
        DEALLOCATE( pPosition, pQuaternion, pOrientation, sitePosition, siteOrientation, siteQuaternion, &
        &           isSiteAssociated, AssociatedSites )
        pPosition = pPositionTemp
        pQuaternion = pQuaternionTemp
        pOrientation = pOrientationTemp
        sitePosition = sitePositionTemp
        siteOrientation = siteOrientationTemp
        siteQuaternion = siteQuaternionTemp
        isSiteAssociated = isSiteAssociatedSave
        AssociatedSites = AssociatedSitesSave
        ! Deallocate temporary arrays
        DEALLOCATE( pPositionTemp, pQuaternionTemp, pOrientationTemp, sitePositionTemp, siteOrientationTemp, &
        &           siteQuaternionTemp, isSiteAssociatedTemp, AssociatedSitesTemp )
      END IF
    ELSE
      ! Recover number of particles
      nParticles = iOldParticles
      ! Recover arrays
      DEALLOCATE( pPosition, pQuaternion, pOrientation, sitePosition, siteOrientation, siteQuaternion, &
      &           isSiteAssociated, AssociatedSites )
      pPosition = pPositionTemp
      pQuaternion = pQuaternionTemp
      pOrientation = pOrientationTemp
      sitePosition = sitePositionTemp
      siteOrientation = siteOrientationTemp
      siteQuaternion = siteQuaternionTemp
      isSiteAssociated = isSiteAssociatedSave
      AssociatedSites = AssociatedSitesSave
      ! Deallocate temporary arrays
      DEALLOCATE( pPositionTemp, pQuaternionTemp, pOrientationTemp, sitePositionTemp, siteOrientationTemp, &
      &           siteQuaternionTemp, isSiteAssociatedTemp, AssociatedSitesTemp )
    END IF

  END IF

  ! Widom insertion
  IF( isWidomEnabled .AND. MOD( iCycle, nWidomFrequency ) == 0 ) THEN
    ! Create random position within unit box
    CALL Random_Number( RandomNumber )
    ScalingDistanceUnitBox(1) = RandomNumber - 0.5D0 ! X coordinate
    CALL Random_Number( RandomNumber )
    ScalingDistanceUnitBox(2) = RandomNumber - 0.5D0 ! Y coordinate
    CALL Random_Number( RandomNumber )
    ScalingDistanceUnitBox(3) = RandomNumber - 0.5D0 ! Z coordinate
    ! Put new particle into the simulation box
    CALL MatrixVectorMultiplication( BoxLength, ScalingDistanceUnitBox, iCreatedPosition )
    ! Create random quaternion
    CALL RandomQuaternionGenerator( iCreatedQuaternion, cPi ) ! Angle in [-π, π]
    ! Create random orientation
    CALL VectorRotation( zAxis, iCreatedQuaternion, iCreatedOrientation )
    ! Random sites
    DO iSite = 1, 4
      CALL VectorRotation( ReferencePositions(:,iSite), iCreatedQuaternion, RotatedPosition(:) )
      ! Random positions
      iCreatedSitePosition(:,iSite) = iCreatedPosition(:) + RotatedPosition(:)
      ! Random orientations
      iCreatedSiteOrientation(:,iSite) = RotatedPosition(:) / SQRT( DOT_PRODUCT( RotatedPosition(:), RotatedPosition(:) ) )
      ! Random quaternions
      CALL QuaternionMultiplication( iCreatedQuaternion, ReferenceQuaternions(:,iSite), iCreatedSiteQuaternion(:,iSite) )
    END DO
    ! Calculate new energy of inserted particle
    iCreatedPotential = 0.D0
    CALL Widom_Energy( iCreatedPosition, iCreatedSitePosition, iCreatedPotential, iCreatedVirial, iCreatedAssociationTerm, &
    &                  iCreatedMieSWTerm, iCreatedSteeleTerm, BoxLength, BoxLengthInverse, Overlap )
    ! Evaluate energy of insertion
    IF( .NOT. Overlap ) THEN
      WidomEnergy = WidomEnergy + EXP( - iCreatedPotential )
    END IF
    ! Increment counter
    counterWidom = counterWidom + 1
    ! Widom data
    WRITE( 65, "(G0,1X,G0,1X,G0.6,1X,G0.6,1X,G0.6)" ) iCycle, counterWidom, WidomEnergy, &
    &                                                 LOG( nParticles * cbBroglieWavelength / BoxVolume ), &
    &                                                 - LOG( WidomEnergy / DBLE( counterWidom ) )
  END IF

  ! Adjustment of maximum displacements
  IF( iCycle <= nEquilibration ) THEN ! During equilibration only
    ! Translational adjustment
    IF( MOD( iCycle, nAdjustment ) == 0 ) THEN
      IF( nMovementTranslationCounter > 0 ) THEN
        ! Acceptance ratio (non-overlapping microstates over sampled microstates)
        Ratio = DBLE( nAcceptanceTranslation ) / DBLE( nMovementTranslationCounter )
        ! Translational adjustment
        IF( Ratio <= AcceptanceRatioTranslation ) THEN
          maxTranslationalDisplc  = 0.95D0 * maxTranslationalDisplc
        ELSE
          maxTranslationalDisplc  = 1.05D0 * maxTranslationalDisplc
        END IF
        ! Ratio data
        WRITE( 20, "(G0,3(' ',G0.6))" ) iCycle, Ratio, maxTranslationalDisplc, AcceptanceRatioTranslation
        FLUSH( 20 )
        ! Reset counter
        nAcceptanceTranslation = 0
        nMovementTranslationCounter = 0
        ! Avoid multiple turns
        IF( maxTranslationalDisplc >= MAXVAL( BoxLength ) ) THEN
          maxTranslationalDisplc = MAXVAL( BoxLength )
        END IF
      END IF
    END IF
    ! Rotational adjustment
    IF( MOD( iCycle, nAdjustment ) == 0 ) THEN
      IF( nMovementRotationCounter > 0 ) THEN
        ! Acceptance ratio (non-overlapping microstates over sampled microstates)
        Ratio = DBLE( nAcceptanceRotation ) / DBLE( nMovementRotationCounter )
        ! Rotation adjustment
        IF( Ratio <= AcceptanceRatioRotation ) THEN
          maxRotationalDisplc = 0.95D0 * maxRotationalDisplc
        ELSE
          maxRotationalDisplc = 1.05D0 * maxRotationalDisplc
        END IF
        ! Ratio data
        WRITE( 30, "(G0,3(' ',G0.6))" ) iCycle, Ratio, maxRotationalDisplc, AcceptanceRatioRotation
        FLUSH( 30 )
        ! Reset counter
        nAcceptanceRotation = 0
        nMovementRotationCounter  = 0
        ! Avoid multiple turns
        IF( maxRotationalDisplc >= cPi ) THEN
          maxRotationalDisplc = cPi
        END IF
      END IF
    END IF
    ! Creation adjustment
    IF( MOD( iCycle, nAdjustment ) == 0 ) THEN
      IF( nCreationCounter > 0 ) THEN
        ! Acceptance ratio (non-overlapping microstates over sampled microstates)
        Ratio = DBLE( nAcceptanceCreation ) / DBLE( nCreationCounter )
        ! Ratio data
        WRITE( 25, "(G0,(' ',G0.6))" ) iCycle, Ratio
        FLUSH( 25 )
        ! Reset counter
        nAcceptanceCreation = 0
        nCreationCounter = 0
      END IF
    END IF
    ! Destruction adjustment
    IF( MOD( iCycle, nAdjustment ) == 0 ) THEN
      IF( nDestructionCounter > 0 ) THEN
        ! Acceptance ratio (non-overlapping microstates over sampled microstates)
        Ratio = DBLE( nAcceptanceDestruction ) / DBLE( nDestructionCounter )
        ! Ratio data
        WRITE( 35, "(G0,(' ',G0.6))" ) iCycle, Ratio
        FLUSH( 35 )
        ! Reset counter
        nAcceptanceDestruction = 0
        nDestructionCounter = 0
      END IF
    END IF
  END IF

  ! Histogram
  IF( MOD( iCycle, nSaveConfig ) == 0 ) THEN
    IF( iCycle >= nEquilibration ) THEN  
      DO iParticle = 1, nParticles
        bin = FLOOR( ( pPosition(3,iParticle) + ( 0.5D0 * poreWidth ) ) / binWidth ) + 1
        Histogram(bin) = Histogram(bin) + 1.D0
      END DO
      nValidConfigs = nValidConfigs + 1
    END IF
  END IF

  ! Save potential energy data
  IF( MOD( iCycle, nSave ) == 0 ) THEN
    pPotential = tPotentialEnergy / DBLE( nParticles )
    WRITE( 40, "(5(G0,1X))" ) iCycle, pPotential, tAssociationTerm / DBLE( nParticles ), tMieSWTerm / DBLE( nParticles ), &
    &                         tSteeleTerm / DBLE( nParticles )
    FLUSH( 40 )
  END IF

  ! Save fraction of free sites data
  IF( MOD( iCycle, nSave ) == 0 ) THEN
    nAssocSites = COUNT( isSiteAssociated ) / 2
    fractionAssociatedSites = DBLE( nAssocSites ) / DBLE( tAssocSites )
    fractionFreeSites = 1.D0 - fractionAssociatedSites
    WRITE( 55, "(G0,1X,G0,1X,G0,1X,G0)" ) iCycle, nAssocSites, tAssocSites, fractionFreeSites
    FLUSH( 55 )
  END IF

  ! Save number of particles data
  IF( MOD( iCycle, nSave ) == 0 ) THEN
    calculatedReducedChemPot = excessReducedChemPot + LOG( DBLE( nParticles ) ) - GrandCanonicalConstant
    WRITE( 60, "(G0,1X,G0,1X,G0,1X,G0)" ) iCycle, DBLE( nParticles ) / BoxVolume, nParticles, calculatedReducedChemPot
    FLUSH( 60 )
  END IF

  ! Save trajectory configuration
  IF( MOD( iCycle, nSaveConfig ) == 0 ) THEN
    WRITE( 50, "(G0)" ) nParticles * 5
    DescriptorString = "(G0,8(G0.9,1X),G0.9,G0,2(G0.9,1X),G0.9,2G0)"
    WRITE( 50, DescriptorString ) 'Lattice="', BoxLength(1:9), '" Origin="', -0.5D0 * ( BoxLength(1) + BoxLength(4) + &
    &                             BoxLength(7) ), -0.5D0 * ( BoxLength(2) + BoxLength(5) + BoxLength(8) ), -0.5D0 * &
    &                             ( BoxLength(3) + BoxLength(6) + BoxLength(9) ), '" ', &
    &                             "Properties=id:I:1:species:S:1:pos:R:3:orientation:R:4:aspherical_shape:R:3"
    DO iParticle = 1, nParticles
      WRITE( 50, "(12(G0,1X))" ) iParticle, "C", pPosition(1:3,iParticle), pQuaternion(1:3,iParticle), pQuaternion(0,iParticle), &
      &                          0.5D0 * SigmaFluid, 0.5D0 * SigmaFluid, 0.5D0 * SigmaFluid
      DO iSite = 1, 4
        WRITE( 50, "(12(G0,1X))" ) iParticle, SiteNames(iSite), sitePosition(1:3,iSite,iParticle), &
        &                          siteQuaternion(1:3,iSite,iParticle), siteQuaternion(0,iSite,iParticle), &
        &                          0.1D0 * SigmaFluid, 0.1D0 * SigmaFluid, 0.4D0 * SigmaFluid
      END DO
    END DO
    FLUSH( 50 )
  END IF

END DO

! Close files
CLOSE( Unit= 20 )
CLOSE( Unit= 25 )
CLOSE( Unit= 30 )
CLOSE( Unit= 35 )
CLOSE( Unit= 40 )
CLOSE( Unit= 50 )
CLOSE( Unit= 55 )
CLOSE( Unit= 60 )
IF( isWidomEnabled ) CLOSE( Unit= 65 )

! Widom's chemical potential
IF( isWidomEnabled ) THEN
  WidomEnergy = WidomEnergy / DBLE( counterWidom )
  ExcessWidomChemicalPotential = - LOG( WidomEnergy )
  WRITE( *, "(G0,1X,G0.6)" ) "Widom Excess Chemical Potential (reduced units):", ExcessWidomChemicalPotential
END IF

! Write radial distribution function data
OPEN( Unit= 60, File= "rdf.dat" )
DO iSlab = 1 , nslabs
  avgHistogram(iSlab) = Histogram(iSlab) / DBLE(nValidConfigs)
  WRITE( 60, "(G0,1X,G0)" ) 0.5D0 * binWidth + DBLE( iSlab - 1 ) * binWidth, avgHistogram(iSlab)
END DO
CLOSE( Unit= 60 )

! Deallocate arrays
DEALLOCATE( Histogram, avgHistogram )

END PROGRAM MonteCarlo_Mie

! ************************************************************************************* !
!                                     PROGRESS BAR                                      !
! ************************************************************************************* !
!          This subroutine generates a progress bar for the simplex algorithm.          !
! ************************************************************************************* !
SUBROUTINE Progress_Bar( CurrentCycle, TotalCycles, EnsembleType )

! Uses one module: global variables
USE GlobalVariables

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 ) :: CurrentCycle ! Counter (cycles)
INTEGER( Kind= Int64 ) :: TotalCycles  ! Total/Maximum number of cycles
INTEGER( Kind= Int64 ) :: AuxiliarInt1 ! Auxiliary integer variable

! CHARACTER STRINGS
CHARACTER( Len= 22 ) :: ProgressBarMC ! Progress bar
CHARACTER( Len= 04 ) :: EnsembleType  ! Ensemble type

! Progress bar (FORMAT)
IF( ( DBLE(CurrentCycle) / DBLE(TotalCycles) ) * 100.D0 < 10.D0 ) THEN
  ProgressBarMC = "Progress("//TRIM( EnsembleType )//"): ???%"
ELSE IF( ( DBLE(CurrentCycle) / DBLE(TotalCycles) ) * 100.D0 < 100.D0 ) THEN
  ProgressBarMC = "Progress("//TRIM( EnsembleType )//"): ????%"
ELSE IF( ( DBLE(CurrentCycle) / DBLE(TotalCycles) ) * 100.D0 < 1000.D0 ) THEN
  ProgressBarMC = "Progress("//TRIM( EnsembleType )//"): ?????%"
END IF

! Auxiliar integer variable
AuxiliarInt1 = LEN_TRIM( EnsembleType ) - 3

! Progress bar (replace character positions)
IF( ( DBLE( CurrentCycle ) / DBLE( TotalCycles ) ) * 100.D0 < 10.D0 ) THEN
  WRITE( Unit= ProgressBarMC(16+AuxiliarInt1:18+AuxiliarInt1), Fmt= "(F3.1)" ) ( DBLE( CurrentCycle ) / DBLE( TotalCycles ) ) * &
  &                                                                            100.D0 - 0.05D0
  IF( ( DBLE( CurrentCycle ) / DBLE( TotalCycles ) ) * 100.D0 - 0.05D0 < 0.D0 ) THEN
    WRITE( Unit= ProgressBarMC(16+AuxiliarInt1:18+AuxiliarInt1), Fmt= "(F3.1)" ) 0.D0
  END IF
ELSE IF( ( DBLE( CurrentCycle ) / DBLE( TotalCycles ) ) * 100.D0 < 100.D0 ) THEN
  WRITE( Unit= ProgressBarMC(16+AuxiliarInt1:19+AuxiliarInt1), Fmt= "(F4.1)" ) ( DBLE( CurrentCycle ) / DBLE( TotalCycles ) ) * &
  &                                                                            100.D0 - 0.05D0
ELSE IF( ( DBLE( CurrentCycle ) / DBLE( TotalCycles ) ) * 100.D0 < 1000.D0 ) THEN
  WRITE( Unit= ProgressBarMC(16+AuxiliarInt1:20+AuxiliarInt1), Fmt= "(F5.1)" ) ( DBLE( CurrentCycle ) / DBLE( TotalCycles ) ) * &
  &                                                                            100.D0
END IF

! Print progress bar
IF( ( DBLE( CurrentCycle ) / DBLE( TotalCycles ) ) * 100.D0 < 10.D0 ) THEN
  IF( AuxiliarInt1 == 0 ) THEN
    WRITE( Unit= Output_Unit, Fmt= "(A1,A19)" , Advance= "NO" ) CHAR(13), ProgressBarMC
  ELSE
    WRITE( Unit= Output_Unit, Fmt= "(A1,A20)" , Advance= "NO" ) CHAR(13), ProgressBarMC
  END IF
ELSE IF( ( DBLE( CurrentCycle ) / DBLE( TotalCycles ) ) * 100.D0 < 100.D0 ) THEN
  IF( AuxiliarInt1 == 0 ) THEN
    WRITE( Unit= Output_Unit, Fmt= "(A1,A20)" , Advance= "NO" ) CHAR(13), ProgressBarMC
  ELSE
    WRITE( Unit= Output_Unit, Fmt= "(A1,A21)" , Advance= "NO" ) CHAR(13), ProgressBarMC
  END IF
ELSE IF( ( DBLE( CurrentCycle ) / DBLE( TotalCycles ) ) * 100.D0 < 1000.D0 ) THEN
  IF( AuxiliarInt1 == 0 ) THEN
    WRITE( Unit= Output_Unit, Fmt= "(A1,A21)" , Advance= "NO" ) CHAR(13), ProgressBarMC
  ELSE
    WRITE( Unit= Output_Unit, Fmt= "(A1,A22)" , Advance= "NO" ) CHAR(13), ProgressBarMC
  END IF
END IF

! Flush standard output unit
FLUSH( Unit= Output_Unit )

RETURN


END SUBROUTINE Progress_Bar
