! ************************************************************************************************ !
!                                   Combining and Mixing Rules                                     !
! ************************************************************************************************ !
!           This subroutine is used to calculate the combining and mixing rules for the            !
!             molecular parameters (e.g., diameter, well depth, and potential range).              !
! ************************************************************************************************ !
!     => AUTHOR:     Nathan Barros de Souza                                                        !
!     => E-MAIL:     n264179@dac.unicamp.br                                                        !
!     => SUPERVISOR: Luís Fernando Mercier Franco                                                  !
! ************************************************************************************************ !
!                       University of Campinas, Campinas - São Paulo, Brazil                       !
!                                         March 7th, 2024                                          !
! ************************************************************************************************ !
! MAIN REFERENCES:     A. Galindo, L. A. Davies, A. Gil-Villegas, G. Jackson                       !
!                                 J. Mol. Phys., 93, 241-252, 1998                                 !
! ************************************************************************************************ !
! Version 1.0.0                                                                                    !
! ************************************************************************************************ !

! ************************************************************************************************ !
! Combining Rules                                                                                  !
! ************************************************************************************************ !
SUBROUTINE Combining_Rules(  )

! Uses two modules: global variables and substances
USE GlobalVariables
USE Substances

IMPLICIT NONE

! ************************************************************************************************ !
! INTEGER VARIABLES                                                                                !
! ************************************************************************************************ !
INTEGER( Kind= Int64 ) :: iComponent, jComponent ! Component indexes

! Combining rule for the diameter
DO iComponent = 1, nComponents
  DO jComponent = 1, nComponents
    ijDiameter( iComponent, jComponent ) = 0.5D0 * ( cDiameter( iComponent ) + cDiameter( jComponent ) )
  END DO
END DO

! Combining rule for the length
DO iComponent = 1, nComponents
  DO jComponent = 1, nComponents
    ijLength( iComponent, jComponent ) = 0.5D0 * ( cLength( iComponent ) + cLength( jComponent ) )
  END DO
END DO

! Combining rule for the aspect ratio
DO iComponent = 1, nComponents
  DO jComponent = 1, nComponents
    ijAspectRatio( iComponent, jComponent ) = ijLength( iComponent, jComponent ) / ijDiameter( iComponent, jComponent )
  END DO
END DO

! Combining rule for the diameter of a sphere with the same volume of the non-spherical geometry
DO iComponent = 1, nComponents
  DO jComponent = 1, nComponents
    ijDiameterSphere( iComponent, jComponent ) = 0.5D0 * ( cDiameterSphere( iComponent ) + cDiameterSphere( jComponent ) )
  END DO
END DO

! Combining rule for the well depth
DO iComponent = 1, nComponents
  DO jComponent = 1, nComponents
    ijWellDepth( iComponent, jComponent ) = DSQRT( cWellDepth( iComponent ) * cWellDepth( jComponent ) )
    ! Correction for the well depth
    IF( iComponent /= jComponent ) THEN
      ijWellDepth( iComponent, jComponent ) = ijWellDepth( iComponent, jComponent ) * ijWellDepthCorrection
    END IF
  END DO
END DO

! Combining rule for the potential range
DO iComponent = 1, nComponents
  DO jComponent = 1, nComponents
    ijPotentialRange( iComponent, jComponent ) = ( cDiameter( iComponent ) * cPotentialRange( iComponent ) + &
    &     cDiameter( jComponent ) * cPotentialRange( jComponent ) ) / ( cDiameter( iComponent ) + cDiameter( jComponent ) )
    ! Correction for the potential range
    IF( iComponent /= jComponent ) THEN
      ijPotentialRange( iComponent, jComponent ) = ijPotentialRange( iComponent, jComponent ) * ijPotentialRangeCorrection
    END IF
  END DO
END DO

RETURN

END SUBROUTINE Combining_Rules

! ************************************************************************************************ !
! Mixing Rules                                                                                     !
! ************************************************************************************************ !
SUBROUTINE Mixing_Rules( mFraction, cMatrix, mProperty )

! Uses two modules: global variables and substances
USE GlobalVariables
USE Substances

IMPLICIT NONE

! ************************************************************************************************ !
! INTEGER VARIABLES                                                                                !
! ************************************************************************************************ !
INTEGER( Kind= Int64 ) :: iComponent, jComponent ! Component indexes

! ************************************************************************************************ !
! REAL VARIABLES                                                                                   !
! ************************************************************************************************ !
REAL( Kind= Real64 )                                        :: mProperty ! Mixed property
REAL( Kind= Real64 ), DIMENSION( nComponents )              :: mFraction ! Molar fraction
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents ) :: cMatrix   ! Matrix of the combined properties

! Mixing rule
mProperty = 0.D0
DO iComponent = 1, nComponents
  DO jComponent = 1, nComponents
    mProperty = mProperty + mFraction(iComponent) * mFraction(jComponent) * cMatrix(iComponent, jComponent)
  END DO
END DO

RETURN

END SUBROUTINE Mixing_Rules

! ************************************************************************************************ !
! Mixing Rules                                                                                     !
! ************************************************************************************************ !
SUBROUTINE Mixing_Rules_Single( mFraction, cMatrix, mProperty )

! Uses two modules: global variables and substances
USE GlobalVariables
USE Substances

IMPLICIT NONE

! ************************************************************************************************ !
! INTEGER VARIABLES                                                                                !
! ************************************************************************************************ !
INTEGER( Kind= Int64 ) :: iComponent ! Component index

! ************************************************************************************************ !
! REAL VARIABLES                                                                                   !
! ************************************************************************************************ !
REAL( Kind= Real64 )                                        :: mProperty ! Mixed property
REAL( Kind= Real64 ), DIMENSION( nComponents )              :: mFraction ! Molar fraction
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents ) :: cMatrix   ! Matrix of the combined properties

! Mixing rule
mProperty = 0.D0
DO iComponent = 1, nComponents
  mProperty = mProperty + mFraction(iComponent) * cMatrix(iComponent, iComponent)
END DO

RETURN

END SUBROUTINE Mixing_Rules_Single

! ************************************************************************************************ !
! Derivative of mixing rules with respect to the number of particles                               !
! ************************************************************************************************ !
SUBROUTINE Derivative_Mixing_Rules_dN( cComponent, mFraction, cMatrix, dFraction, dMatrix, mProperty )

! Uses two modules: global variables and substances
USE GlobalVariables
USE Substances

IMPLICIT NONE

! ************************************************************************************************ !
! INTEGER VARIABLES                                                                                !
! ************************************************************************************************ !
INTEGER( Kind= Int64 ) :: iComponent, jComponent, cComponent ! Component indexes

! ************************************************************************************************ !
! REAL VARIABLES                                                                                   !
! ************************************************************************************************ !
REAL( Kind= Real64 )                                                     :: mProperty ! Mixed property
REAL( Kind= Real64 ), DIMENSION( nComponents )                           :: mFraction ! Molar fraction
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )              :: cMatrix   ! Matrix of the properties
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )              :: dFraction ! Matrix of the derivatives of the molar fractions
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents, nComponents ) :: dMatrix   ! Matrix of the derivatives of the coefficients

! Derivative of the mixing rules
mProperty = 0.D0
DO iComponent = 1, nComponents
  DO jComponent = 1, nComponents
    mProperty = mProperty + dFraction(cComponent,iComponent) * mFraction(jComponent) * cMatrix(iComponent,jComponent) + &
    &     mFraction(iComponent) * dFraction(cComponent,jComponent) * cMatrix(iComponent,jComponent) + mFraction(iComponent) * &
    &     mFraction(jComponent) * dMatrix(cComponent,iComponent,jComponent)
  END DO
END DO

RETURN

END SUBROUTINE Derivative_Mixing_Rules_dN