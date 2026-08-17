! ############################################################################################### !
!                          FREE-ENERGY METHODS FOR NON-CONVEX MOLECULES                           !
!             This code contains the symmetry-dependent orientational field functions             !
!                                    used in the main program.                                    !
!                                                                                                 !
! Version number: 1.0.0                                                                           !
! ############################################################################################### !
!                                University of Campinas (Unicamp)                                 !
!                                 School of Chemical Engineering                                  !
!                                     Nathan Barros de Souza                                      !
!                             --------------------------------------                              !
!                             Supervisor: Luís Fernando Mercier Franco                            !
!                             --------------------------------------                              !
!                                        August 13th, 2026                                        !
! ############################################################################################### !
! Main References:                  M. P. Allen, D. J. Tildesley                                  !
!                           Oxford University Press, 2nd Edition (2017)                           !
!                             DOI: 10.1093/oso/9780198803195.001.0001                             !
!                             --------------------------------------                              !
!                                          G. Marsaglia                                           !
!                            Ann. Math. Statist. 43(2), 645-646 (1972)                            !
!                                  DOI: 10.1214/aoms/1177692644                                   !
!                             ---------------------------------------                             !
!                     J. Graaf, L. Filion, M. Marechal, R. Roij, M. Dijkstra                      !
!                                J. Chem. Phys. 137, 214101 (2012)                                !
!                                     DOI: 10.1063/1.4767529                                      !
! ############################################################################################### !
! Disclaimer note: Authors assume no responsibility or liability for the use of this code.        !
! ############################################################################################### !
MODULE SymmetryOrientationalFields

! Uses one module: global variables
USE GlobalVariables

CONTAINS

! *********************************************************************************************** !
!                                    Multidimensional Function                                    !
!             This is the function to be integrated through Monte Carlo integration.              !
! *********************************************************************************************** !
DOUBLE PRECISION FUNCTION SymmetryRotationFunction( IntegralDimension, IndependentVariable )

! *********************************************************************************************** !
! INTEGER VARIABLES                                                                               !
! *********************************************************************************************** !
INTEGER( Kind= Int64 ) :: IntegralDimension ! Dimension of the integral

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 ), DIMENSION( IntegralDimension ) :: IndependentVariable ! Independent variables

! Function
IF( IntegralDimension == 1 ) THEN
  SymmetryRotationFunction = IndependentVariable(1) ! Ignore
ELSE IF( IntegralDimension == 2 ) THEN
  SymmetryRotationFunction = IndependentVariable(1) * IndependentVariable(2) ! Ignore
ELSE IF( IntegralDimension == 3 ) THEN
  ! Point group symmetry
  SymmetryRotationFunction = PointGroup_Dnh( IndependentVariable(2), IndependentVariable(3) )
  SymmetryRotationFunction = DEXP( - rSpringConstant * SymmetryRotationFunction ) * DSIN( IndependentVariable(1) )
END IF

RETURN

END FUNCTION SymmetryRotationFunction

! *********************************************************************************************** !
!                 Potential energy for a molecule with a point group of type Dn,h                 !
! *********************************************************************************************** !
DOUBLE PRECISION FUNCTION PointGroup_Dnh( PsiAngleA, PsiAngleB )

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 ) :: PsiAngleA, PsiAngleB ! Lattice-restricted angles

! Function
PointGroup_Dnh = DSIN( 2.D0 * PsiAngleA ) * DSIN( 2.D0 * PsiAngleA ) + DSIN( PsiAngleB ) * DSIN( PsiAngleB )

RETURN

END FUNCTION PointGroup_Dnh

END MODULE SymmetryOrientationalFields