MODULE Quaternions

! Uses one module: global variables
USE GlobalVariables

IMPLICIT NONE

CONTAINS

! ********************************************************************************************************************************* !
! This subroutine takes a body-fixed orientation/position (<BodyOrientationAxis>) and a rotation quaternion (<RotationQuaternion>)  !
! and generates a space-fixed orientation/position (<SpaceFixedAxis>) using a 3D-rotation matrix.                                   !
! ********************************************************************************************************************************* !
! For more information, please read: Allen and Tildesley, 2nd Edition (2017), pages 106-111.                                        !
! ********************************************************************************************************************************* !
SUBROUTINE ActiveTransformation( BodyOrientationAxis, RotationQuaternion, SpaceFixedAxis )

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 ) :: iRow, jCol ! Counters

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 3 ), INTENT(IN)   :: BodyOrientationAxis ! Body-fixed orientation/position [XYZ]
REAL( Kind= Real64 ), DIMENSION( 3 ), INTENT(OUT)  :: SpaceFixedAxis      ! Space-fixed orientation/position [XYZ]
REAL( Kind= Real64 ), DIMENSION( 0:3 ), INTENT(IN) :: RotationQuaternion  ! Rotation quaternion [WXYZ]
REAL( Kind= Real64 ), DIMENSION( 3, 3 )            :: RotationMatrix      ! Rotation matrix
REAL( Kind= Real64 ), DIMENSION( 3, 3 )            :: tRotationMatrix     ! Transpose of rotation matrix

! Rotation Matrix - Allen and Tildesley, 2nd edition, page 110
RotationMatrix(1,1) = ( RotationQuaternion(0) * RotationQuaternion(0) ) + ( RotationQuaternion(1) * RotationQuaternion(1) ) - &
&                     ( RotationQuaternion(2) * RotationQuaternion(2) ) - ( RotationQuaternion(3) * RotationQuaternion(3) )
RotationMatrix(1,2) = 2.D0 * ( RotationQuaternion(1) * RotationQuaternion(2) + RotationQuaternion(0) * RotationQuaternion(3) )
RotationMatrix(1,3) = 2.D0 * ( RotationQuaternion(1) * RotationQuaternion(3) - RotationQuaternion(0) * RotationQuaternion(2) )
RotationMatrix(2,1) = 2.D0 * ( RotationQuaternion(1) * RotationQuaternion(2) - RotationQuaternion(0) * RotationQuaternion(3) )
RotationMatrix(2,2) = ( RotationQuaternion(0) * RotationQuaternion(0) ) - ( RotationQuaternion(1) * RotationQuaternion(1) ) + &
&                     ( RotationQuaternion(2) * RotationQuaternion(2) ) - ( RotationQuaternion(3) * RotationQuaternion(3) )
RotationMatrix(2,3) = 2.D0 * ( RotationQuaternion(2) * RotationQuaternion(3) + RotationQuaternion(0) * RotationQuaternion(1) )
RotationMatrix(3,1) = 2.D0 * ( RotationQuaternion(1) * RotationQuaternion(3) + RotationQuaternion(0) * RotationQuaternion(2) )
RotationMatrix(3,2) = 2.D0 * ( RotationQuaternion(2) * RotationQuaternion(3) - RotationQuaternion(0) * RotationQuaternion(1) )
RotationMatrix(3,3) = ( RotationQuaternion(0) * RotationQuaternion(0) ) - ( RotationQuaternion(1) * RotationQuaternion(1) ) - &
&                     ( RotationQuaternion(2) * RotationQuaternion(2) ) + ( RotationQuaternion(3) * RotationQuaternion(3) )

! Transpose of rotation matrix
DO iRow = 1, 3
  DO jCol = 1, 3
    tRotationMatrix(iRow,jCol) = RotationMatrix(jCol,iRow)
  END DO
END DO

! Active tranformation (dot product of body-fixed vector and transpose of rotation matrix)
SpaceFixedAxis(1) = ( BodyOrientationAxis(1) * tRotationMatrix(1,1) ) + ( BodyOrientationAxis(2) * tRotationMatrix(1,2) ) + &
&                    ( BodyOrientationAxis(3) * tRotationMatrix(1,3) )
SpaceFixedAxis(2) = ( BodyOrientationAxis(1) * tRotationMatrix(2,1) ) + ( BodyOrientationAxis(2) * tRotationMatrix(2,2) ) + &
&                    ( BodyOrientationAxis(3) * tRotationMatrix(2,3) )
SpaceFixedAxis(3) = ( BodyOrientationAxis(1) * tRotationMatrix(3,1) ) + ( BodyOrientationAxis(2) * tRotationMatrix(3,2) ) + &
&                    ( BodyOrientationAxis(3) * tRotationMatrix(3,3) )

RETURN

END SUBROUTINE ActiveTransformation

! ********************************************************************************************************************************* !
! This subroutine rotates a vector (<PointVector>) using Hamilton's quaternion product.                                             !
! ********************************************************************************************************************************* !
SUBROUTINE VectorRotation( PointVector, RotationQuaternion, RotatedVector )

IMPLICIT NONE

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 0:3 )             :: ConjugateQuaternion ! Conjugate quaternion [WXYZ]
REAL( Kind= Real64 ), DIMENSION( 0:3 )             :: RotatedQuaternion   ! Rotated quaternion [WXYZ]
REAL( Kind= Real64 ), DIMENSION( 0:3 )             :: TempQuaternion      ! Temporary quaternion [WXYZ]
REAL( Kind= Real64 ), DIMENSION( 3 ), INTENT(IN)   :: PointVector         ! Point vector [XYZ]
REAL( Kind= Real64 ), DIMENSION( 3 ), INTENT(OUT)  :: RotatedVector       ! Rotated vector [XYZ]
REAL( Kind= Real64 ), DIMENSION( 0:3 ), INTENT(IN) :: RotationQuaternion  ! Rotation quaternion [WXYZ]

! Multiplication of quaternions to find the rotated vector
CALL QuaternionMultiplication( RotationQuaternion, [ 0.D0, PointVector ], TempQuaternion )

! Conjugate quaternion
ConjugateQuaternion = [ RotationQuaternion(0), -RotationQuaternion(1), -RotationQuaternion(2), -RotationQuaternion(3) ]

! Multiplication of quaternions to find the rotated vector
CALL QuaternionMultiplication( TempQuaternion, ConjugateQuaternion, RotatedQuaternion )

! Rotated vector
RotatedVector(1:3) = RotatedQuaternion(1:3)

RETURN

END SUBROUTINE VectorRotation

! ********************************************************************************************************************************* !
! This subroutine creates a composed rotation quaternion (<ComposedUnitQuaternion>) by multiplying a reference quaternion           !
! (<ReferenceUnitQuaternion>) with a random quaternion (<RandomUnitQuaternion>).                                                    !
! ********************************************************************************************************************************* !
SUBROUTINE QuaternionMultiplication( RandomUnitQuaternion, ReferenceUnitQuaternion, ComposedUnitQuaternion )

IMPLICIT NONE

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 0:3 ), INTENT(IN)  :: RandomUnitQuaternion    ! Random quaternion [WXYZ]
REAL( Kind= Real64 ), DIMENSION( 0:3 ), INTENT(IN)  :: ReferenceUnitQuaternion ! Reference quaternion [WXYZ]
REAL( Kind= Real64 ), DIMENSION( 0:3 ), INTENT(OUT) :: ComposedUnitQuaternion  ! Composed quaternion [WXYZ]

! Cross product of quaternions (qr × qm)
ComposedUnitQuaternion(0) = ( RandomUnitQuaternion(0) * ReferenceUnitQuaternion(0) ) - &
&                           ( RandomUnitQuaternion(1) * ReferenceUnitQuaternion(1) ) - &
&                           ( RandomUnitQuaternion(2) * ReferenceUnitQuaternion(2) ) - &
&                           ( RandomUnitQuaternion(3) * ReferenceUnitQuaternion(3) )
ComposedUnitQuaternion(1) = ( RandomUnitQuaternion(1) * ReferenceUnitQuaternion(0) ) + &
&                           ( RandomUnitQuaternion(0) * ReferenceUnitQuaternion(1) ) - &
&                           ( RandomUnitQuaternion(3) * ReferenceUnitQuaternion(2) ) + &
&                           ( RandomUnitQuaternion(2) * ReferenceUnitQuaternion(3) )
ComposedUnitQuaternion(2) = ( RandomUnitQuaternion(2) * ReferenceUnitQuaternion(0) ) + &
&                           ( RandomUnitQuaternion(3) * ReferenceUnitQuaternion(1) ) + &
&                           ( RandomUnitQuaternion(0) * ReferenceUnitQuaternion(2) ) - &
&                           ( RandomUnitQuaternion(1) * ReferenceUnitQuaternion(3) )
ComposedUnitQuaternion(3) = ( RandomUnitQuaternion(3) * ReferenceUnitQuaternion(0) ) - &
&                           ( RandomUnitQuaternion(2) * ReferenceUnitQuaternion(1) ) + &
&                           ( RandomUnitQuaternion(1) * ReferenceUnitQuaternion(2) ) + &
&                           ( RandomUnitQuaternion(0) * ReferenceUnitQuaternion(3) )

RETURN

END SUBROUTINE QuaternionMultiplication

! ********************************************************************************************************************************* !
! This subroutine generates a random quaternion (<RandomUnitQuaternion>) from a random angle (<Angle>) and random unit axis         !
! (<RandomUnitVector>).                                                                                                             !
! ********************************************************************************************************************************* !
SUBROUTINE RandomQuaternionGenerator( RandomUnitQuaternion, Angle )

IMPLICIT NONE

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 )               :: RandomAngle = 0.D0 ! Random angle
REAL( Kind= Real64 ), INTENT( IN ) :: Angle              ! Maximum angular displacement

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 3 )                  :: RandomUnitVector = 0.D0 ! Random vector
REAL( Kind= Real64 ), DIMENSION( 0:3 ), INTENT( OUT ) :: RandomUnitQuaternion    ! Random unit quaternion

! Random angle
CALL Random_Number( RandomNumber )
RandomAngle = ( (2.D0 * RandomNumber) - 1.D0 ) * Angle ! Range [-angmax,angmax]

! Random unit vector
CALL RandomVectorGenerator( RandomUnitVector )

! Quaternion algebra
RandomUnitQuaternion(0) = DCOS( RandomAngle * 0.5D0 )                       ! Real part
RandomUnitQuaternion(1) = DSIN( RandomAngle * 0.5D0 ) * RandomUnitVector(1) ! Imaginary part (Vector)
RandomUnitQuaternion(2) = DSIN( RandomAngle * 0.5D0 ) * RandomUnitVector(2) ! Imaginary part (Vector)
RandomUnitQuaternion(3) = DSIN( RandomAngle * 0.5D0 ) * RandomUnitVector(3) ! Imaginary part (Vector)

RETURN

END SUBROUTINE RandomQuaternionGenerator

! ********************************************************************************************************************************* !
! This subroutine generates a random unit vector (<RandomUnitVector>) on the surface of a unit sphere.                              !
! ********************************************************************************************************************************* !
! For more information, please read: Allen and Tildesley, 2nd Edition (2017), page 514.                                             !
! ********************************************************************************************************************************* !   
SUBROUTINE RandomVectorGenerator( RandomUnitVector )

IMPLICIT NONE

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 ) :: Alpha, Beta, Zeta ! Random numbers

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( : ), INTENT( OUT ) :: RandomUnitVector ! Random vector

! Marsaglia's routine
UnitVectorCriterion: DO
  ! Uniform random number, α
  CALL Random_Number( RandomNumber )
  Alpha = (2.D0 * RandomNumber) - 1.D0
  ! Uniform random number, β
  CALL Random_Number( RandomNumber )
  Beta = (2.D0 * RandomNumber) - 1.D0
  ! Sum of squares, ζ
  Zeta  = (Alpha * Alpha) + (Beta * Beta)
  ! Marseglia's criterion
  IF( Zeta < 1.D0 ) THEN
    EXIT UnitVectorCriterion
  END IF
END DO UnitVectorCriterion

! Random vector
RandomUnitVector(1) = 2.D0 * Alpha * DSQRT( 1.D0 - Zeta )
RandomUnitVector(2) = 2.D0 * Beta * DSQRT( 1.D0 - Zeta )
RandomUnitVector(3) = 1.D0 - (2.D0 * Zeta)

RETURN

END SUBROUTINE RandomVectorGenerator

END MODULE Quaternions