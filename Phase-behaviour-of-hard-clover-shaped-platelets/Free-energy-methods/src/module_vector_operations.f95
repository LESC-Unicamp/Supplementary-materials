! ############################################################################################### !
!                          FREE-ENERGY METHODS FOR NON-CONVEX MOLECULES                           !
!       This code contains vector operation routines subroutines used in the main program.        !
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
MODULE VectorOperations

! Uses one module: global variables
USE GlobalVariables

CONTAINS

! *********************************************************************************************** !
!               This subroutine calculates the inverse of a matrix using cofactors                !
! *********************************************************************************************** !
SUBROUTINE InverseMatrixCofactorVec( InputMatrix, InverseMatrix, Determinant )

IMPLICIT NONE

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 )                 :: Determinant        ! Determinant
REAL( Kind= Real64 )                 :: DeterminantInverse ! Inverse of determinant
REAL( Kind= Real64 ), DIMENSION( 9 ) :: InverseMatrix      ! Matrix (inverse)
REAL( Kind= Real64 ), DIMENSION( 9 ) :: InputMatrix        ! Matrix (input)
REAL( Kind= Real64 ), DIMENSION( 9 ) :: tCofactorMatrix    ! Cofactor matrix (transpose)

! Tranpose matrix of the matrix of cofactors
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

! *********************************************************************************************** !
!                             This subroutine multiplies two matrices                             !
! *********************************************************************************************** !
SUBROUTINE MatrixVectorMultiplication( InputMatrix, InputVector, OutputVector )

IMPLICIT NONE

! *********************************************************************************************** !
! REAL VARIABLES                                                                                  !
! *********************************************************************************************** !
REAL( Kind= Real64 ), DIMENSION( 9 ) :: InputMatrix  ! Input matrix
REAL( Kind= Real64 ), DIMENSION( 3 ) :: InputVector  ! Vector (input)
REAL( Kind= Real64 ), DIMENSION( 3 ) :: OutputVector ! Vector (output)

! Multiplication of a matrix and a vector
OutputVector(1) = InputMatrix(1) * InputVector(1) + InputMatrix(4) * InputVector(2) + InputMatrix(7) * InputVector(3)
OutputVector(2) = InputMatrix(2) * InputVector(1) + InputMatrix(5) * InputVector(2) + InputMatrix(8) * InputVector(3)
OutputVector(3) = InputMatrix(3) * InputVector(1) + InputMatrix(6) * InputVector(2) + InputMatrix(9) * InputVector(3)

RETURN

END SUBROUTINE MatrixVectorMultiplication

END MODULE VectorOperations