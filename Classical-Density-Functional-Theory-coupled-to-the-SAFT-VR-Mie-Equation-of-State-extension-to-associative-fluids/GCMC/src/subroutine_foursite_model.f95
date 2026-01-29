! ----------------------------------------------------------------------- !
! Site Initialization Subroutine                                          !
! ----------------------------------------------------------------------- !
SUBROUTINE Site_Initialization(  )

USE GlobalVariables
USE QuaternionOperations, ONLY: VectorRotation

IMPLICIT NONE

REAL( Kind= Real64 ) :: wQuaternion(0:3), FakePosition(3), RotationAngle

! Fake position vector
FakePosition = [ 0.D0, 0.D0, SiteRadius ]

! First site
RotationAngle = 0.5D0 * AngleWater
wQuaternion(0) = COS( RotationAngle * 0.5D0 )
wQuaternion(1) = SIN( RotationAngle * 0.5D0 ) * yAxis(1)
wQuaternion(2) = SIN( RotationAngle * 0.5D0 ) * yAxis(2)
wQuaternion(3) = SIN( RotationAngle * 0.5D0 ) * yAxis(3)
ReferenceQuaternions(0,1) = wQuaternion(0)
ReferenceQuaternions(1,1) = wQuaternion(1)
ReferenceQuaternions(2,1) = wQuaternion(2)
ReferenceQuaternions(3,1) = wQuaternion(3)
CALL VectorRotation( FakePosition, wQuaternion, ReferencePositions( :, 1 ) )

! Second site
RotationAngle = - 0.5D0 * AngleWater
wQuaternion(0) = COS( RotationAngle * 0.5D0 )
wQuaternion(1) = SIN( RotationAngle * 0.5D0 ) * yAxis(1)
wQuaternion(2) = SIN( RotationAngle * 0.5D0 ) * yAxis(2)
wQuaternion(3) = SIN( RotationAngle * 0.5D0 ) * yAxis(3)
ReferenceQuaternions(0,2) = wQuaternion(0)
ReferenceQuaternions(1,2) = wQuaternion(1)
ReferenceQuaternions(2,2) = wQuaternion(2)
ReferenceQuaternions(3,2) = wQuaternion(3)
CALL VectorRotation( FakePosition, wQuaternion, ReferencePositions( :, 2 ) )

! Third site
RotationAngle = cPi - 0.5D0 * AngleWater
wQuaternion(0) = COS( RotationAngle * 0.5D0 )
wQuaternion(1) = SIN( RotationAngle * 0.5D0 ) * xAxis(1)
wQuaternion(2) = SIN( RotationAngle * 0.5D0 ) * xAxis(2)
wQuaternion(3) = SIN( RotationAngle * 0.5D0 ) * xAxis(3)
ReferenceQuaternions(0,3) = wQuaternion(0)
ReferenceQuaternions(1,3) = wQuaternion(1)
ReferenceQuaternions(2,3) = wQuaternion(2)
ReferenceQuaternions(3,3) = wQuaternion(3)
CALL VectorRotation( FakePosition, wQuaternion, ReferencePositions( :, 3 ) )

! Fourth site
RotationAngle = - cPi + 0.5D0 * AngleWater
wQuaternion(0) = COS( RotationAngle * 0.5D0 )
wQuaternion(1) = SIN( RotationAngle * 0.5D0 ) * xAxis(1)
wQuaternion(2) = SIN( RotationAngle * 0.5D0 ) * xAxis(2)
wQuaternion(3) = SIN( RotationAngle * 0.5D0 ) * xAxis(3)
ReferenceQuaternions(0,4) = wQuaternion(0)
ReferenceQuaternions(1,4) = wQuaternion(1)
ReferenceQuaternions(2,4) = wQuaternion(2)
ReferenceQuaternions(3,4) = wQuaternion(3)
CALL VectorRotation( FakePosition, wQuaternion, ReferencePositions( :, 4 ) )

! Site names
SiteNames(1) = "O"
SiteNames(2) = "O"
SiteNames(3) = "H"
SiteNames(4) = "H"

RETURN

END SUBROUTINE Site_Initialization

! ----------------------------------------------------------------------- !
! Attribute Associated Sites Subroutine                                   !
! ----------------------------------------------------------------------- !
SUBROUTINE Attribute_Sites( BoxLengthAssociation, BoxLengthInverseAssociation )

USE GlobalVariables
USE VectorOperations

IMPLICIT NONE

INTEGER( Kind= Int64 ) :: iParticle, jParticle, iSite, jSite

REAL( Kind= Real64 ) :: BoxLengthInverseAssociation(9), BoxLengthAssociation(9)
REAL( Kind= Real64 ) :: iPosition(3), jPosition(3), VectorDistance(3), ScalingDistanceUnitBox(3), SquaredDistance
REAL( Kind= Real64 ) :: iSitePosition(3), jSitePosition(3)

! Initialization
isSiteAssociated = .FALSE.
AssociatedSites = 0

! Loop through all particles
IF( isAssociativePotentialEnabled ) THEN
  DO iParticle = 1, nParticles
    DO jParticle = 1, nParticles
      ! Cycle if index j <= i (prevents double counting)
      IF( jParticle <= iParticle ) CYCLE
      ! Position of particle i
      iPosition(:) = pPosition(:,iParticle)
      ! Position of particle j
      jPosition(:) = pPosition(:,jParticle)
      ! Vector distance
      VectorDistance = jPosition - iPosition
      ! Apply periodic boundary conditions
      CALL MatrixVectorMultiplication( BoxLengthInverseAssociation, VectorDistance, ScalingDistanceUnitBox )
      IF( Confinement ) THEN
        ScalingDistanceUnitBox(1:2) = ScalingDistanceUnitBox(1:2) - ANINT( ScalingDistanceUnitBox(1:2) ) ! Periodic boundary conditions (system is confined along z-axis)
      ELSE
        ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
      END IF
      CALL MatrixVectorMultiplication( BoxLengthAssociation, ScalingDistanceUnitBox, VectorDistance )
      ! Squared distance
      SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
      ! Association effect
      IF( SquaredDistance < cutoffSquaredRadiusMie ) THEN
        ! Sites of particle i
        DO iSite = 1, 4
          ! Sites of particle j
          DO jSite = 1, 4
            ! Do not evaluate site A of particle i if it is already associated
            IF( isSiteAssociated(iSite,iParticle) ) CYCLE
            ! Do not evaluate site B of particle j if it is already associated
            IF( isSiteAssociated(jSite,jParticle) ) CYCLE
            ! Conditions for non-bonding interactions
            IF( iSite == jSite ) CYCLE ! Skip interaction O-O or H-H
            IF( iSite == 1 .AND. jSite == 2 ) CYCLE ! Skip interaction O-O
            IF( iSite == 2 .AND. jSite == 1 ) CYCLE ! Skip interaction O-O
            IF( iSite == 3 .AND. jSite == 4 ) CYCLE ! Skip interaction H-H
            IF( iSite == 4 .AND. jSite == 3 ) CYCLE ! Skip interaction H-H
            ! Site position of particle i
            iSitePosition(:) = sitePosition(:,iSite,iParticle)
            ! Site position of particle j
            jSitePosition(:) = sitePosition(:,jSite,jParticle)
            ! Vector distance
            VectorDistance = jSitePosition - iSitePosition
            ! Apply periodic boundary conditions
            CALL MatrixVectorMultiplication( BoxLengthInverseAssociation, VectorDistance, ScalingDistanceUnitBox )
            IF( Confinement ) THEN
              ScalingDistanceUnitBox(1:2) = ScalingDistanceUnitBox(1:2) - ANINT( ScalingDistanceUnitBox(1:2) ) ! Periodic boundary conditions (system is confined along z-axis)
            ELSE
              ScalingDistanceUnitBox = ScalingDistanceUnitBox - ANINT( ScalingDistanceUnitBox )
            END IF
            CALL MatrixVectorMultiplication( BoxLengthAssociation, ScalingDistanceUnitBox, VectorDistance )
            ! Squared distance
            SquaredDistance = DOT_PRODUCT( VectorDistance, VectorDistance )
            ! Four-site model potential
            IF( SquaredDistance < cutoffSquaredRadiusWater ) THEN
              ! Mutual association
              isSiteAssociated(iSite,iParticle) = .TRUE.
              isSiteAssociated(jSite,jParticle) = .TRUE.
              ! Associated sites
              AssociatedSites(1,iSite,iParticle) = jSite     ! 1 for associated site
              AssociatedSites(2,iSite,iParticle) = jParticle ! 2 for associated particle
              AssociatedSites(1,jSite,jParticle) = iSite     ! 1 for associated site
              AssociatedSites(2,jSite,jParticle) = iParticle ! 2 for associated particle
            END IF
          END DO
        END DO
      END IF
    END DO
  END DO
END IF

RETURN

END SUBROUTINE Attribute_Sites