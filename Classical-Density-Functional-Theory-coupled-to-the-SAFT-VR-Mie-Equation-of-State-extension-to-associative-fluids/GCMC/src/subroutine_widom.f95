! ----------------------------------------------------------------------- !
! Particle Energy Calculation Subroutine                                  !
! ----------------------------------------------------------------------- !
SUBROUTINE Widom_Energy( iPosition, iSitePosition, iPotentialEnergy, iVirial, iAssociationTerm, iMieSWTerm, &
&                        iSteeleTerm, BoxLengthEnergy, BoxLengthInverseEnergy, Overlap )

USE GlobalVariables
USE VectorOperations

IMPLICIT NONE

INTEGER( Kind= Int64 ) :: jParticle, iSite, jSite

REAL( Kind= Real64 ) :: iPotentialEnergy, iVirial, AccumulateEnergy, AccumulateVirial, VectorDistance(3), ScalingDistanceUnitBox(3)
REAL( Kind= Real64 ) :: iAssociationTerm, iMieSWTerm, iSteeleTerm, BoxLengthEnergy(9), BoxLengthInverseEnergy(9)
REAL( Kind= Real64 ) :: iPosition(3), jPosition(3), iSitePosition(3,4), jSitePosition(3,4), SquaredDistance

LOGICAL :: Overlap, iWidomAssociated(4), jWidomAssociated(4)

! Initialization
iPotentialEnergy = 0.D0
iVirial = 0.D0
iAssociationTerm = 0.D0
iMieSWTerm = 0.D0
iSteeleTerm = 0.D0
Overlap = .FALSE.
iWidomAssociated = .FALSE.

! Loop through j particles
DO jParticle = 1, nParticles
  ! Position of particle j
  jPosition(:) = pPosition(:,jParticle)
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
    ! Associated sites
    jWidomAssociated(:) = isSiteAssociated(:,jParticle)
    ! Association effect
    DO iSite = 1, 4
      DO jSite = 1, 4
        ! Conditions
        IF( iSite == jSite ) CYCLE ! Skip interaction O-O or H-H
        IF( iSite == 1 .AND. jSite == 2 ) CYCLE ! Skip interaction O-O
        IF( iSite == 2 .AND. jSite == 1 ) CYCLE ! Skip interaction O-O
        IF( iSite == 3 .AND. jSite == 4 ) CYCLE ! Skip interaction H-H
        IF( iSite == 4 .AND. jSite == 3 ) CYCLE ! Skip interaction H-H
        ! Skip if associated
        IF( jWidomAssociated(jSite) ) CYCLE ! Only free sites
        IF( iWidomAssociated(iSite) ) CYCLE ! Only free sites (safe guard)
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
          iWidomAssociated(iSite) = .TRUE.
          jWidomAssociated(jSite) = .TRUE.
        END IF
      END DO
    END DO
  END IF
END DO

! Interaction with hard walls
IF( Confinement .AND. isSolidPotentialEnabled ) THEN
  CALL Steele_Potential( iPosition(3), AccumulateEnergy, Overlap )
  Overlap = .TRUE.
  ! Increment energy
  iPotentialEnergy = iPotentialEnergy + AccumulateEnergy
  ! Set Steele contribution
  iSteeleTerm = iSteeleTerm + AccumulateEnergy
END IF

RETURN

END SUBROUTINE Widom_Energy