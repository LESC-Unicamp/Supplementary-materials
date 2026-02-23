MODULE VectorOperations

! Uses one module: global variables
USE GlobalVariables

IMPLICIT NONE

CONTAINS

! ********************************************************************************************************************************* !
! This subroutine calculates the cross product of two vectors <VectorA> and <VectorB>, resuling in <OutputVector>.                  !
! ********************************************************************************************************************************* !
SUBROUTINE Cross_Product( VectorA, VectorB, OutputVector )

IMPLICIT NONE

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( : ), INTENT( IN )  :: VectorA      ! Vector A
REAL( Kind= Real64 ), DIMENSION( : ), INTENT( IN )  :: VectorB      ! Vector B
REAL( Kind= Real64 ), DIMENSION( : ), INTENT( OUT ) :: OutputVector ! Output vector

! Cross product
OutputVector(1) = VectorA(2) * VectorB(3) - VectorA(3) * VectorB(2)
OutputVector(2) = VectorA(3) * VectorB(1) - VectorA(1) * VectorB(3)
OutputVector(3) = VectorA(1) * VectorB(2) - VectorA(2) * VectorB(1)

RETURN

END SUBROUTINE Cross_Product

! ********************************************************************************************************************************* !
! This subroutine calculates the inverse matrix (<InverseMatrix>) using cofactors.                                                  !
! ********************************************************************************************************************************* !
SUBROUTINE InverseMatrixCofactorVec( InputMatrix, InverseMatrix, Determinant )

IMPLICIT NONE

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 )              :: DeterminantInverse ! Inverse of determinant
REAL( Kind= Real64 ), INTENT(OUT) :: Determinant        ! Determinant

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 9 )              :: tCofactorMatrix ! Inverse matrix
REAL( Kind= Real64 ), DIMENSION( : ), INTENT(OUT) :: InverseMatrix   ! Matrix (inverse)
REAL( Kind= Real64 ), DIMENSION( : ), INTENT(IN)  :: InputMatrix     ! Matrix (input)

! Transpose matrix of the matrix of cofactors
tCofactorMatrix(1) = InputMatrix(5) * InputMatrix(9) - InputMatrix(6) * InputMatrix(8)
tCofactorMatrix(2) = InputMatrix(3) * InputMatrix(8) - InputMatrix(2) * InputMatrix(9)
tCofactorMatrix(3) = InputMatrix(2) * InputMatrix(6) - InputMatrix(3) * InputMatrix(5)
tCofactorMatrix(4) = InputMatrix(6) * InputMatrix(7) - InputMatrix(4) * InputMatrix(9)
tCofactorMatrix(5) = InputMatrix(1) * InputMatrix(9) - InputMatrix(3) * InputMatrix(7)
tCofactorMatrix(6) = InputMatrix(3) * InputMatrix(4) - InputMatrix(1) * InputMatrix(6)
tCofactorMatrix(7) = InputMatrix(4) * InputMatrix(8) - InputMatrix(5) * InputMatrix(7)
tCofactorMatrix(8) = InputMatrix(2) * InputMatrix(7) - InputMatrix(1) * InputMatrix(8)
tCofactorMatrix(9) = InputMatrix(1) * InputMatrix(5) - InputMatrix(2) * InputMatrix(4)

! Determinant of matrix
Determinant = InputMatrix(1) * tCofactorMatrix(1) + InputMatrix(4) * tCofactorMatrix(2) + InputMatrix(7) * tCofactorMatrix(3)

! Inverse of the determinant of matrix
DeterminantInverse = 0.D0
IF( DABS( Determinant ) > 0.D0 ) THEN
  DeterminantInverse = 1.0D0 / Determinant
END IF

! Inverse matrix
InverseMatrix = DeterminantInverse * tCofactorMatrix

RETURN

END SUBROUTINE InverseMatrixCofactorVec

! ********************************************************************************************************************************* !
! This subroutine multiplies a matrix (<InputMatrix>) and a vector (<InputVector>), resulting in another vector (<OutputVector>).   !
! ********************************************************************************************************************************* !
SUBROUTINE MatrixVectorMultiplication( InputMatrix, InputVector, OutputVector )

IMPLICIT NONE

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( : ), INTENT(IN)  :: InputMatrix  ! Input matrix
REAL( Kind= Real64 ), DIMENSION( : ), INTENT(IN)  :: InputVector  ! Input vector
REAL( Kind= Real64 ), DIMENSION( : ), INTENT(OUT) :: OutputVector ! Output vector

! Multiplication of a matrix and a vector
OutputVector(1) = InputMatrix(1) * InputVector(1) + InputMatrix(4) * InputVector(2) + InputMatrix(7) * InputVector(3)
OutputVector(2) = InputMatrix(2) * InputVector(1) + InputMatrix(5) * InputVector(2) + InputMatrix(8) * InputVector(3)
OutputVector(3) = InputMatrix(3) * InputVector(1) + InputMatrix(6) * InputVector(2) + InputMatrix(9) * InputVector(3)

RETURN

END SUBROUTINE MatrixVectorMultiplication

! ********************************************************************************************************************************* !
! This subroutine converts a vector (<CartesianVector>) in Cartesian coordinates to spherical coordinates.                          !
! ********************************************************************************************************************************* !
SUBROUTINE ConvertCartesianToSpherical( CartesianVector, RadialDistance, PolarAngle, AzimuthalAngle )

IMPLICIT NONE

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 ), INTENT( OUT ) :: RadialDistance ! Radial distance
REAL( Kind= Real64 ), INTENT( OUT ) :: PolarAngle     ! Polar angle
REAL( Kind= Real64 ), INTENT( OUT ) :: AzimuthalAngle ! Azimuthal angle

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( : ), INTENT( IN ) :: CartesianVector ! Cartesian vector

! Radial distance
RadialDistance = DSQRT( DOT_PRODUCT( CartesianVector, CartesianVector ) )

! Polar angle
PolarAngle = DACOS( CartesianVector(3) / RadialDistance ) ! Range: [0, π]

! Azimuthal angle
AzimuthalAngle = DATAN2( CartesianVector(2), CartesianVector(1) ) ! Range: [-π, π]
IF( AzimuthalAngle < 0.D0 ) AzimuthalAngle = AzimuthalAngle + 2.D0 * cPi ! Range: [0, 2π] (to be consistent with the notation in Grey and Gubbins)

RETURN

END SUBROUTINE ConvertCartesianToSpherical

END MODULE VectorOperations