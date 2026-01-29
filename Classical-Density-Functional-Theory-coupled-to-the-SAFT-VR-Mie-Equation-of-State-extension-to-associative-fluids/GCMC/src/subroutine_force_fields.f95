! ----------------------------------------------------------------------- !
! Mie Potential Subroutine                                                !
! ----------------------------------------------------------------------- !
SUBROUTINE Mie_Potential( sqDistance, PotentialEnergy, Virial )

USE GlobalVariables

IMPLICIT NONE

REAL( Kind= Real64 ) :: sqDistance, PotentialEnergy, Virial, distanceRatio, distanceRatioAttractive, distanceRatioRepulsive

! Initialization
PotentialEnergy = 0.D0
Virial = 0.D0

! Mie potential calculation
distanceRatio = SQRT( squaredSigmaFluid / sqDistance )
distanceRatioAttractive = distanceRatio ** ( attractiveRange )
distanceRatioRepulsive = distanceRatio ** ( repulsiveRange )
PotentialEnergy = distanceRatioRepulsive - distanceRatioAttractive
PotentialEnergy = PotentialEnergy * reducedEpsilonFluid * NormalisationCoefficient
Virial = 2D-2

RETURN

END SUBROUTINE Mie_Potential

! ----------------------------------------------------------------------- !
! SW Potential Subroutine                                                 !
! ----------------------------------------------------------------------- !
SUBROUTINE SW_Potential( PotentialEnergy )

USE GlobalVariables

IMPLICIT NONE

REAL( Kind= Real64 ) :: PotentialEnergy

! Initialization
PotentialEnergy = - reducedEpsilonFluid

RETURN

END SUBROUTINE SW_Potential

! ----------------------------------------------------------------------- !
! Associative Potential Subroutine                                                !
! ----------------------------------------------------------------------- !
SUBROUTINE FourSite_Model( sqDistance, PotentialEnergy )

USE GlobalVariables

IMPLICIT NONE

REAL( Kind= Real64 ) :: sqDistance, PotentialEnergy

! Initialization
PotentialEnergy = 0.D0

! Association potential
IF( sqDistance < cutoffSquaredRadiusWater ) THEN
  PotentialEnergy = - reducedEpsilonWater
END IF

RETURN

END SUBROUTINE FourSite_Model

! ----------------------------------------------------------------------- !
! Steele Potential Subroutine                                             !
! ----------------------------------------------------------------------- !
SUBROUTINE Steele_Potential( zPosition, PotentialEnergy, Overlap )

USE GlobalVariables

IMPLICIT NONE

REAL( Kind= Real64 ) :: zPosition, adjustedPosition
REAL( Kind= Real64 ) :: PotentialEnergy

LOGICAL :: Overlap

! Initialization
PotentialEnergy = 0.D0
Overlap = .FALSE.

! Adjusted position
adjustedPosition = zPosition + 0.5D0 * poreWidth

! Steele potential
PotentialEnergy = 2.D0 * cPi * SolidDensity * (EpsilonSolid_Fluid / Temperature) * (SigmaSolid_Fluid * SigmaSolid_Fluid) * &
&                 InterlayerSpacing * ( 0.4D0 * ( SigmaSolid_Fluid / adjustedPosition ) ** 10 - ( SigmaSolid_Fluid / &
&                 adjustedPosition ) ** 4 - ( SigmaSolid_Fluid ** 4 ) / ( 3.D0 * InterlayerSpacing * ( adjustedPosition + &
&                 0.61D0 * InterlayerSpacing ) ** 3 ) + 0.4D0 * ( SigmaSolid_Fluid / ( poreWidth - adjustedPosition ) ) ** 10 - &
&                 ( SigmaSolid_Fluid / ( poreWidth - adjustedPosition ) ) ** 4 - ( SigmaSolid_Fluid ** 4 ) / ( 3.D0 * &
&                 InterlayerSpacing * ( poreWidth - adjustedPosition + 0.61D0 * InterlayerSpacing ) ** 3 ) )

! Too close to the wall
IF( adjustedPosition < ( 0.5D0 * SigmaFluid ) .OR. adjustedPosition > ( poreWidth - 0.5D0 * SigmaFluid ) ) THEN
  Overlap = .TRUE.
  RETURN
END IF

RETURN

END SUBROUTINE Steele_Potential