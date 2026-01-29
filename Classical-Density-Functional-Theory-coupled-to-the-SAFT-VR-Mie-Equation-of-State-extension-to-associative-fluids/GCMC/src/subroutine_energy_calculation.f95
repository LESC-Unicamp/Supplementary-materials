! ----------------------------------------------------------------------- !
! Total Energy Calculation Subroutine                                     !
! ----------------------------------------------------------------------- !
SUBROUTINE Calculation_Total_Energy( BoxLengthEnergy, BoxLengthInverseEnergy, tPotentialEnergy, VirialContribution, &
&                                    AssociationTerm, MieSWTerm, SteeleTerm, RandomConfig, Overlap )

USE GlobalVariables
USE VectorOperations

IMPLICIT NONE

INTEGER( Kind= Int64 ) :: iParticle, jParticle, iSite, jSite

REAL( Kind= Real64 ) :: tPotentialEnergy, AccumulateEnergy, VirialContribution, AccumulateVirial
REAL( Kind= Real64 ) :: AssociationTerm, MieSWTerm, SteeleTerm, BoxLengthInverseEnergy(9), BoxLengthEnergy(9)
REAL( Kind= Real64 ) :: iPosition(3), jPosition(3), VectorDistance(3), ScalingDistanceUnitBox(3), SquaredDistance
REAL( Kind= Real64 ) :: iSitePosition(3), jSitePosition(3)

LOGICAL :: RandomConfig, Overlap

! Initialization
tPotentialEnergy = 0.D0
VirialContribution = 0.D0
AssociationTerm = 0.D0
MieSWTerm = 0.D0
SteeleTerm = 0.D0
Overlap = .FALSE.

! Loop through particles
DO iParticle = 1, nParticles
  DO jParticle = 1, nParticles
    ! Cycle if index j <= i (prevents double counting)
    IF( jParticle <= iParticle ) CYCLE
    ! Position of particle i
    iPosition(:) = pPosition(:,iParticle)
    IF( RandomConfig ) iPosition(:) = pPositionRandom(:,iParticle)
    ! Position of particle j
    jPosition(:) = pPosition(:,jParticle)
    IF( RandomConfig ) jPosition(:) = pPositionRandom(:,jParticle)
    ! Vector distance
    VectorDistance = jPosition - iPosition
    ! Apply periodic boundary conditions
    CALL MatrixVectorMultiplication( BoxLengthInverseEnergy, VectorDistance, ScalingDistanceUnitBox )
    IF( Confinement ) THEN
      ScalingDistanceUnitBox(1:2) = ScalingDistanceUnitBox(1:2) - ANINT( ScalingDistanceUnitBox(1:2) ) ! Periodic boundary conditions (system is confined along z-axis)
    ELSE
      ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
    END IF
    CALL MatrixVectorMultiplication( BoxLengthEnergy, ScalingDistanceUnitBox, VectorDistance )
    ! Squared distance
    SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
    ! Check cutoff (Mie)
    IF( isMiePotentialEnabled .AND. SquaredDistance < cutoffSquaredRadiusMie ) THEN
      ! Long-range dispersive effect
      CALL Mie_Potential( SquaredDistance, AccumulateEnergy, AccumulateVirial )
      ! Increment energy and virial
      tPotentialEnergy = tPotentialEnergy + AccumulateEnergy
      VirialContribution = VirialContribution + AccumulateVirial
      ! Set Mie contribution
      MieSWTerm = MieSWTerm + AccumulateEnergy
    END IF
    ! Check cutoff (SW)
    IF( isSWPotentialEnabled .AND. ( SquaredDistance <= sqSWRange .AND. SquaredDistance > squaredSigmaFluid ) ) THEN
      ! Long-range square-well effect
      CALL SW_Potential( AccumulateEnergy )
      ! Increment energy
      tPotentialEnergy = tPotentialEnergy + AccumulateEnergy
      ! Set SW contribution
      MieSWTerm = MieSWTerm + AccumulateEnergy
    END IF
    ! Check cutoff
    ! Check overlap
    IF( isHardSpherePotentialEnabled .AND. SquaredDistance <= squaredSigmaFluid ) THEN
      Overlap = .TRUE.
      RETURN ! Exit subroutine if overlap detected
    END IF
    ! Association effect
    IF( isAssociativePotentialEnabled ) THEN
      DO iSite = 1, 4
        ! Conditions
        IF( .NOT. isSiteAssociated(iSite,iParticle) ) CYCLE ! Site not associated
        IF( AssociatedSites(2,iSite,iParticle) /= jParticle ) CYCLE ! Site associated but not with particle j 
        DO jSite = 1, 4
          ! Conditions
          IF( .NOT. isSiteAssociated(jSite,jParticle) ) CYCLE ! Site not associated
          IF( AssociatedSites(2,jSite,jParticle) /= iParticle ) CYCLE ! Site associated but not with particle i
          IF( AssociatedSites(1,iSite,iParticle) /= jSite ) CYCLE ! Site A of particle i not associated with site B of particle j
          IF( iSite == jSite ) CYCLE ! Skip interaction O-O or H-H
          IF( iSite == 1 .AND. jSite == 2 ) CYCLE ! Skip interaction O-O
          IF( iSite == 2 .AND. jSite == 1 ) CYCLE ! Skip interaction O-O
          IF( iSite == 3 .AND. jSite == 4 ) CYCLE ! Skip interaction H-H
          IF( iSite == 4 .AND. jSite == 3 ) CYCLE ! Skip interaction H-H
          ! Site position of particle i
          iSitePosition(:) = sitePosition(:,iSite,iParticle)
          IF( RandomConfig ) iSitePosition(:) = sitePositionRandom(:,iSite,iParticle)
          ! Site position of particle j
          jSitePosition(:) = sitePosition(:,jSite,jParticle)
          IF( RandomConfig ) jSitePosition(:) = sitePositionRandom(:,jSite,jParticle)
          ! Vector distance
          VectorDistance = jSitePosition - iSitePosition
          ! Apply periodic boundary conditions
          CALL MatrixVectorMultiplication( BoxLengthInverseEnergy, VectorDistance, ScalingDistanceUnitBox )
          IF( Confinement ) THEN
            ScalingDistanceUnitBox(1:2) = ScalingDistanceUnitBox(1:2) - ANINT( ScalingDistanceUnitBox(1:2) ) ! Periodic boundary conditions (system is confined along z-axis)
          ELSE
            ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
          END IF
          CALL MatrixVectorMultiplication( BoxLengthEnergy, ScalingDistanceUnitBox, VectorDistance )
          ! Squared distance
          SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
          ! Four-site model potential
          CALL FourSite_Model( SquaredDistance, AccumulateEnergy )
          ! Increment energy
          tPotentialEnergy = tPotentialEnergy + AccumulateEnergy
          ! Set association contribution
          AssociationTerm = AssociationTerm + AccumulateEnergy
        END DO
      END DO
    END IF
  END DO
END DO

! Interaction with hard walls
IF( Confinement .AND. isSolidPotentialEnabled ) THEN
  DO iParticle = 1, nParticles
    ! Position of particle i
    iPosition(:) = pPosition(:,iParticle)
    IF( RandomConfig ) iPosition(:) = pPositionRandom(:,iParticle)
    ! Steele potential
    CALL Steele_Potential( iPosition(3), AccumulateEnergy, Overlap )
    IF( Overlap ) RETURN
    ! Increment energy
    tPotentialEnergy = tPotentialEnergy + AccumulateEnergy
    ! Set Steele contribution
    SteeleTerm = SteeleTerm + AccumulateEnergy
  END DO
END IF

RETURN

END SUBROUTINE Calculation_Total_Energy

! ----------------------------------------------------------------------- !
! Particle Energy Calculation Subroutine                                  !
! ----------------------------------------------------------------------- !
SUBROUTINE Particle_Energy( iParticle, iPosition, iSitePosition, iPotentialEnergy, iVirial, iAssociationTerm, iMieSWTerm, &
&                           iSteeleTerm, BoxLengthEnergy, BoxLengthInverseEnergy, RandomConfig, isParticleDisplaced, &
&                           isParticleDestroyed, isParticleCreated, Overlap )

USE GlobalVariables
USE VectorOperations

IMPLICIT NONE

INTEGER( Kind= Int64 ) :: iParticle, jParticle, kParticle, iSite, jSite, kSite, SaveParticle(4), SaveSite(4), Counter, AssocCounter

REAL( Kind= Real64 ) :: iPotentialEnergy, iVirial, AccumulateEnergy, AccumulateVirial, VectorDistance(3), ScalingDistanceUnitBox(3)
REAL( Kind= Real64 ) :: iAssociationTerm, iMieSWTerm, iSteeleTerm, BoxLengthEnergy(9), BoxLengthInverseEnergy(9)
REAL( Kind= Real64 ) :: iPosition(3), jPosition(3), iSitePosition(3,4), jSitePosition(3,4), kSitePosition(3,4), SquaredDistance

LOGICAL :: RandomConfig, isParticleDisplaced, isParticleDestroyed, isParticleCreated, Overlap

! Initialization
iPotentialEnergy = 0.D0
iVirial = 0.D0
iAssociationTerm = 0.D0
iMieSWTerm = 0.D0
iSteeleTerm = 0.D0
Overlap = .FALSE.

! Check for associations before displacement and correct them if necessary
IF( isAssociativePotentialEnabled .AND. isParticleDisplaced .AND. ANY( isSiteAssociated(:,iParticle) ) ) THEN
  ! Initialization
  AssocCounter = 0
  SaveSite = 0
  SaveParticle = 0
  ! Particle i has associated sites
  DO iSite = 1, 4
    ! That is not the associated site
    IF( .NOT. isSiteAssociated(iSite,iParticle) ) CYCLE
    ! Get associated particle j
    jParticle = AssociatedSites(2,iSite,iParticle)
    ! Site positions of particle j
    jSitePosition(:,:) = sitePosition(:,:,jParticle)
    ! Get associated site B
    jSite = AssociatedSites(1,iSite,iParticle)
    ! Vector distance
    VectorDistance = jSitePosition(:,jSite) - iSitePosition(:,iSite)
    ! Apply periodic boundary conditions
    CALL MatrixVectorMultiplication( BoxLengthInverseEnergy, VectorDistance, ScalingDistanceUnitBox )
    IF( Confinement ) THEN
      ScalingDistanceUnitBox(1:2) = ScalingDistanceUnitBox(1:2) - ANINT( ScalingDistanceUnitBox(1:2) ) ! Periodic boundary conditions (system is confined along z-axis)
    ELSE
      ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
    END IF
    CALL MatrixVectorMultiplication( BoxLengthEnergy, ScalingDistanceUnitBox, VectorDistance )
    ! Squared distance
    SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
    ! Four-site model potential
    CALL FourSite_Model( SquaredDistance, AccumulateEnergy )
    ! Increment energy
    iPotentialEnergy = iPotentialEnergy + AccumulateEnergy
    ! Set association contribution
    iAssociationTerm = iAssociationTerm + AccumulateEnergy
    ! Sites no longer associated
    IF( SquaredDistance >= cutoffSquaredRadiusWater ) THEN
      ! Disassociate particles
      isSiteAssociated(iSite,iParticle) = .FALSE.
      isSiteAssociated(jSite,jParticle) = .FALSE.
      AssociatedSites(1,iSite,iParticle) = 0
      AssociatedSites(2,iSite,iParticle) = 0
      AssociatedSites(1,jSite,jParticle) = 0
      AssociatedSites(2,jSite,jParticle) = 0
      ! Save particles and sites
      AssocCounter = AssocCounter + 1
      IF( AssocCounter > 4 ) STOP "ERROR: Too many disassociated sites!"
      SaveParticle(AssocCounter) = jParticle
      SaveSite(AssocCounter) = jSite
    END IF
  END DO
  ! Check disassociated j particles (ignore particle i for now)
  IF( AssocCounter > 0 ) THEN
    DO Counter = 1, AssocCounter
      ! Get particle and site indexes
      jParticle = SaveParticle(Counter)
      jSite = SaveSite(Counter)
      ! Site positions of particle j
      jSitePosition(:,:) = sitePosition(:,:,jParticle)
      ! Loop though k particles (ignore particle i for now)
      DO kParticle = 1, nParticles
        ! Skip itself
        IF( kParticle == jParticle ) CYCLE
        ! Ignore particle i
        IF( kParticle == iParticle ) CYCLE
        ! Site positions of particle k
        kSitePosition(:,:) = sitePosition(:,:,kParticle)
        ! Loop through sites of k
        DO kSite = 1, 4
          IF( isSiteAssociated(kSite,kParticle) ) CYCLE ! Consider only free sites of k
          IF( isSiteAssociated(jSite,jParticle) ) CYCLE ! Consider only free sites of j (safe guard)
          ! Conditions
          IF( jSite == kSite ) CYCLE ! Skip interaction O-O or H-H
          IF( jSite == 1 .AND. kSite == 2 ) CYCLE ! Skip interaction O-O
          IF( jSite == 2 .AND. kSite == 1 ) CYCLE ! Skip interaction O-O
          IF( jSite == 3 .AND. kSite == 4 ) CYCLE ! Skip interaction H-H
          IF( jSite == 4 .AND. kSite == 3 ) CYCLE ! Skip interaction H-H
          ! Vector distance
          VectorDistance = kSitePosition(:,kSite) - jSitePosition(:,jSite)
          ! Apply periodic boundary conditions
          CALL MatrixVectorMultiplication( BoxLengthInverseEnergy, VectorDistance, ScalingDistanceUnitBox )
          IF( Confinement ) THEN
            ScalingDistanceUnitBox(1:2) = ScalingDistanceUnitBox(1:2) - ANINT( ScalingDistanceUnitBox(1:2) ) ! Periodic boundary conditions (system is confined along z-axis)
          ELSE
            ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
          END IF
          CALL MatrixVectorMultiplication( BoxLengthEnergy, ScalingDistanceUnitBox, VectorDistance )
          ! Squared distance
          SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
          ! Four-site model potential
          CALL FourSite_Model( SquaredDistance, AccumulateEnergy )
          ! Increment energy
          iPotentialEnergy = iPotentialEnergy + AccumulateEnergy
          ! Set association contribution
          iAssociationTerm = iAssociationTerm + AccumulateEnergy
          ! Associate particles
          IF( SquaredDistance < cutoffSquaredRadiusWater ) THEN
            isSiteAssociated(jSite,jParticle) = .TRUE.
            isSiteAssociated(kSite,kParticle) = .TRUE.
            AssociatedSites(1,jSite,jParticle) = kSite
            AssociatedSites(2,jSite,jParticle) = kParticle
            AssociatedSites(1,kSite,kParticle) = jSite
            AssociatedSites(2,kSite,kParticle) = jParticle
          END IF
        END DO
      END DO
    END DO
  END IF
END IF

! Check for associations before destruction and correct them if necessary
IF( isAssociativePotentialEnabled .AND. isParticleDestroyed .AND. ANY( isSiteAssociated(:,iParticle) ) ) THEN
  ! Initialization
  AssocCounter = 0
  SaveSite = 0
  SaveParticle = 0
  ! Particle i has associated sites
  DO iSite = 1, 4
    ! That is not the associated site
    IF( .NOT. isSiteAssociated(iSite,iParticle) ) CYCLE
    ! Get associated particle j
    jParticle = AssociatedSites(2,iSite,iParticle)
    ! Site positions of particle j
    jSitePosition(:,:) = sitePosition(:,:,jParticle)
    ! Get associated site B
    jSite = AssociatedSites(1,iSite,iParticle)
    ! Vector distance
    VectorDistance = jSitePosition(:,jSite) - iSitePosition(:,iSite)
    ! Apply periodic boundary conditions
    CALL MatrixVectorMultiplication( BoxLengthInverseEnergy, VectorDistance, ScalingDistanceUnitBox )
    IF( Confinement ) THEN
      ScalingDistanceUnitBox(1:2) = ScalingDistanceUnitBox(1:2) - ANINT( ScalingDistanceUnitBox(1:2) ) ! Periodic boundary conditions (system is confined along z-axis)
    ELSE
      ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
    END IF
    CALL MatrixVectorMultiplication( BoxLengthEnergy, ScalingDistanceUnitBox, VectorDistance )
    ! Squared distance
    SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
    ! Four-site model potential
    CALL FourSite_Model( SquaredDistance, AccumulateEnergy )
    ! Increment energy
    iPotentialEnergy = iPotentialEnergy + AccumulateEnergy
    ! Set association contribution
    iAssociationTerm = iAssociationTerm + AccumulateEnergy
    ! Disassociate particles regardless of the distance (particle is destroyed)
    isSiteAssociated(iSite,iParticle) = .FALSE.
    isSiteAssociated(jSite,jParticle) = .FALSE.
    AssociatedSites(1,iSite,iParticle) = 0
    AssociatedSites(2,iSite,iParticle) = 0
    AssociatedSites(1,jSite,jParticle) = 0
    AssociatedSites(2,jSite,jParticle) = 0
    ! Save particles and sites
    AssocCounter = AssocCounter + 1
    IF( AssocCounter > 4 ) STOP "ERROR: Too many disassociated sites!"
    SaveParticle(AssocCounter) = jParticle
    SaveSite(AssocCounter) = jSite
  END DO
  ! Check disassociated j particles
  IF( AssocCounter > 0 ) THEN
    DO Counter = 1, AssocCounter
      ! Get particle and site indexes
      jParticle = SaveParticle(Counter)
      jSite = SaveSite(Counter)
      ! Site positions of particle j
      jSitePosition(:,:) = sitePosition(:,:,jParticle)
      ! Loop though k particles (ignore particle i for now)
      DO kParticle = 1, nParticles
        ! Skip itself
        IF( kParticle == jParticle ) CYCLE
        ! Ignore particle i (destroyed)
        IF( kParticle == iParticle ) CYCLE
        ! Site positions of particle k
        kSitePosition(:,:) = sitePosition(:,:,kParticle)
        ! Loop though sites of k
        DO kSite = 1, 4
          IF( isSiteAssociated(kSite,kParticle) ) CYCLE ! Consider only free sites of k
          IF( isSiteAssociated(jSite,jParticle) ) CYCLE ! Consider only free sites of j (safe guard)
          ! Conditions
          IF( jSite == kSite ) CYCLE ! Skip interaction O-O or H-H
          IF( jSite == 1 .AND. kSite == 2 ) CYCLE ! Skip interaction O-O
          IF( jSite == 2 .AND. kSite == 1 ) CYCLE ! Skip interaction O-O
          IF( jSite == 3 .AND. kSite == 4 ) CYCLE ! Skip interaction H-H
          IF( jSite == 4 .AND. kSite == 3 ) CYCLE ! Skip interaction H-H
          ! Vector distance
          VectorDistance = kSitePosition(:,kSite) - jSitePosition(:,jSite)
          ! Apply periodic boundary conditions
          CALL MatrixVectorMultiplication( BoxLengthInverseEnergy, VectorDistance, ScalingDistanceUnitBox )
          IF( Confinement ) THEN
            ScalingDistanceUnitBox(1:2) = ScalingDistanceUnitBox(1:2) - ANINT( ScalingDistanceUnitBox(1:2) ) ! Periodic boundary conditions (system is confined along z-axis)
          ELSE
            ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
          END IF
          CALL MatrixVectorMultiplication( BoxLengthEnergy, ScalingDistanceUnitBox, VectorDistance )
          ! Squared distance
          SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
          ! Four-site model potential
          CALL FourSite_Model( SquaredDistance, AccumulateEnergy )
          ! Increment energy
          iPotentialEnergy = iPotentialEnergy + AccumulateEnergy
          ! Set association contribution
          iAssociationTerm = iAssociationTerm + AccumulateEnergy
          ! Associate particles
          IF( SquaredDistance < cutoffSquaredRadiusWater ) THEN
            isSiteAssociated(jSite,jParticle) = .TRUE.
            isSiteAssociated(kSite,kParticle) = .TRUE.
            AssociatedSites(1,jSite,jParticle) = kSite
            AssociatedSites(2,jSite,jParticle) = kParticle
            AssociatedSites(1,kSite,kParticle) = jSite
            AssociatedSites(2,kSite,kParticle) = jParticle
          END IF
        END DO
      END DO
    END DO
  END IF
END IF

! Loop through j particles
DO jParticle = 1, nParticles
  ! Skip itself
  IF( iParticle == jParticle ) CYCLE
  ! Position of particle j
  jPosition(:) = pPosition(:,jParticle)
  IF( RandomConfig ) jPosition(:) = pPositionRandom(:,jParticle)
  ! Vector distance
  VectorDistance = jPosition - iPosition
  ! Apply periodic boundary conditions
  CALL MatrixVectorMultiplication( BoxLengthInverseEnergy, VectorDistance, ScalingDistanceUnitBox )
  IF( Confinement ) THEN
    ScalingDistanceUnitBox(1:2) = ScalingDistanceUnitBox(1:2) - ANINT( ScalingDistanceUnitBox(1:2) ) ! Periodic boundary conditions (system is confined along z-axis)
  ELSE
    ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
  END IF
  CALL MatrixVectorMultiplication( BoxLengthEnergy, ScalingDistanceUnitBox, VectorDistance )
  ! Squared distance
  SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
  ! Check cutoff
  IF( isMiePotentialEnabled .AND. SquaredDistance < cutoffSquaredRadiusMie ) THEN
    ! Long-range dispersive effect
    CALL Mie_Potential( SquaredDistance, AccumulateEnergy, AccumulateVirial )
    ! Increment energy and virial
    iPotentialEnergy = iPotentialEnergy + AccumulateEnergy
    iVirial = iVirial + AccumulateVirial
    ! Set Mie contribution
    iMieSWTerm = iMieSWTerm + AccumulateEnergy
  END IF
  ! Check cutoff (SW)
  IF( isSWPotentialEnabled .AND. ( SquaredDistance <= sqSWRange .AND. SquaredDistance > squaredSigmaFluid ) ) THEN
    ! Long-range square-well effect
    CALL SW_Potential( AccumulateEnergy )
    ! Increment energy
    iPotentialEnergy = iPotentialEnergy + AccumulateEnergy
    ! Set SW contribution
    iMieSWTerm = iMieSWTerm + AccumulateEnergy
  END IF
  ! Check overlap
  IF( isHardSpherePotentialEnabled .AND. SquaredDistance <= squaredSigmaFluid ) THEN
    Overlap = .TRUE.
    RETURN ! Exit subroutine if overlap detected
  END IF
  ! Association effect
  IF( isAssociativePotentialEnabled ) THEN
    ! Site positions of particle j
    jSitePosition(:,:) = sitePosition(:,:,jParticle)
    IF( RandomConfig ) jSitePosition(:,:) = sitePositionRandom(:,:,jParticle)
    ! Association effect
    DO iSite = 1, 4
      DO jSite = 1, 4
        ! Conditions
        IF( iSite == jSite ) CYCLE ! Skip interaction O-O or H-H
        IF( iSite == 1 .AND. jSite == 2 ) CYCLE ! Skip interaction O-O
        IF( iSite == 2 .AND. jSite == 1 ) CYCLE ! Skip interaction O-O
        IF( iSite == 3 .AND. jSite == 4 ) CYCLE ! Skip interaction H-H
        IF( iSite == 4 .AND. jSite == 3 ) CYCLE ! Skip interaction H-H
        ! Different strategies
        IF( isParticleDisplaced ) THEN ! Check if displaced particle now interacts with any other particle
          ! Skip if site is still associated
          IF( isSiteAssociated(jSite,jParticle) ) CYCLE ! We are checking if particle i can associate with others after displacement
          IF( isSiteAssociated(iSite,iParticle) ) CYCLE ! Only free sites of particle i (bounded sites already accounted for)
          ! Vector distance
          VectorDistance = jSitePosition(:,jSite) - iSitePosition(:,iSite)
          ! Apply periodic boundary conditions
          CALL MatrixVectorMultiplication( BoxLengthInverseEnergy, VectorDistance, ScalingDistanceUnitBox )
          IF( Confinement ) THEN
            ScalingDistanceUnitBox(1:2) = ScalingDistanceUnitBox(1:2) - ANINT( ScalingDistanceUnitBox(1:2) ) ! Periodic boundary conditions (system is confined along z-axis)
          ELSE
            ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
          END IF
          CALL MatrixVectorMultiplication( BoxLengthEnergy, ScalingDistanceUnitBox, VectorDistance )
          ! Squared distance
          SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
          ! Four-site model potential
          CALL FourSite_Model( SquaredDistance, AccumulateEnergy )
          ! Increment energy
          iPotentialEnergy = iPotentialEnergy + AccumulateEnergy
          ! Set association contribution
          iAssociationTerm = iAssociationTerm + AccumulateEnergy
          ! Associate sites
          IF( SquaredDistance < cutoffSquaredRadiusWater ) THEN
            isSiteAssociated(iSite,iParticle) = .TRUE.
            isSiteAssociated(jSite,jParticle) = .TRUE.
            AssociatedSites(1,iSite,iParticle) = jSite
            AssociatedSites(2,iSite,iParticle) = jParticle
            AssociatedSites(1,jSite,jParticle) = iSite
            AssociatedSites(2,jSite,jParticle) = iParticle
          END IF
        ELSE IF( isParticleCreated ) THEN ! Check if created particle interacts with any other particle
          ! Skip if associated
          IF( isSiteAssociated(jSite,jParticle) ) CYCLE ! Only free sites
          IF( isSiteAssociated(iSite,iParticle) ) CYCLE ! Only free sites (safe guard)
          ! Vector distance
          VectorDistance = jSitePosition(:,jSite) - iSitePosition(:,iSite)
          ! Apply periodic boundary conditions
          CALL MatrixVectorMultiplication( BoxLengthInverseEnergy, VectorDistance, ScalingDistanceUnitBox )
          IF( Confinement ) THEN
            ScalingDistanceUnitBox(1:2) = ScalingDistanceUnitBox(1:2) - ANINT( ScalingDistanceUnitBox(1:2) ) ! Periodic boundary conditions (system is confined along z-axis)
          ELSE
            ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
          END IF
          CALL MatrixVectorMultiplication( BoxLengthEnergy, ScalingDistanceUnitBox, VectorDistance )
          ! Squared distance
          SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
          ! Four-site model potential
          CALL FourSite_Model( SquaredDistance, AccumulateEnergy )
          ! Increment energy
          iPotentialEnergy = iPotentialEnergy + AccumulateEnergy
          ! Set association contribution
          iAssociationTerm = iAssociationTerm + AccumulateEnergy
          ! Associate sites
          IF( SquaredDistance < cutoffSquaredRadiusWater ) THEN
            isSiteAssociated(iSite,iParticle) = .TRUE.
            isSiteAssociated(jSite,jParticle) = .TRUE.
            AssociatedSites(1,iSite,iParticle) = jSite
            AssociatedSites(2,iSite,iParticle) = jParticle
            AssociatedSites(1,jSite,jParticle) = iSite
            AssociatedSites(2,jSite,jParticle) = iParticle
          END IF
        ELSE ! If it is a calculation of the old potential before displacement, then check the logical array
          ! Condition
          IF( .NOT. isSiteAssociated(iSite,iParticle) ) CYCLE ! Site A of particle i is not associated
          IF( .NOT. isSiteAssociated(jSite,jParticle) ) CYCLE ! Site B of particle j is not associated
          IF( AssociatedSites(2,iSite,iParticle) /= jParticle ) CYCLE ! Particle i not associated to j
          IF( AssociatedSites(2,jSite,jParticle) /= iParticle ) CYCLE ! Particle j not associated to i
          IF( AssociatedSites(1,iSite,iParticle) /= jSite ) CYCLE ! Site A of particle i not associated to site B of particle j
          IF( AssociatedSites(1,jSite,jParticle) /= iSite ) CYCLE ! Site B of particle j not associated to site A of particle i
          ! Vector distance
          VectorDistance = jSitePosition(:,jSite) - iSitePosition(:,iSite)
          ! Apply periodic boundary conditions
          CALL MatrixVectorMultiplication( BoxLengthInverseEnergy, VectorDistance, ScalingDistanceUnitBox )
          IF( Confinement ) THEN
            ScalingDistanceUnitBox(1:2) = ScalingDistanceUnitBox(1:2) - ANINT( ScalingDistanceUnitBox(1:2) ) ! Periodic boundary conditions (system is confined along z-axis)
          ELSE
            ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
          END IF
          CALL MatrixVectorMultiplication( BoxLengthEnergy, ScalingDistanceUnitBox, VectorDistance )
          ! Squared distance
          SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
          ! Four-site model potential
          CALL FourSite_Model( SquaredDistance, AccumulateEnergy )
          ! Increment energy
          iPotentialEnergy = iPotentialEnergy + AccumulateEnergy
          ! Set association contribution
          iAssociationTerm = iAssociationTerm + AccumulateEnergy
        END IF
      END DO
    END DO
  END IF
END DO

! Interaction with hard walls
IF( Confinement .AND. isSolidPotentialEnabled ) THEN
  CALL Steele_Potential( iPosition(3), AccumulateEnergy, Overlap )
  IF( Overlap ) RETURN
  ! Increment energy
  iPotentialEnergy = iPotentialEnergy + AccumulateEnergy
  ! Set Steele contribution
  iSteeleTerm = iSteeleTerm + AccumulateEnergy
END IF

RETURN

END SUBROUTINE Particle_Energy