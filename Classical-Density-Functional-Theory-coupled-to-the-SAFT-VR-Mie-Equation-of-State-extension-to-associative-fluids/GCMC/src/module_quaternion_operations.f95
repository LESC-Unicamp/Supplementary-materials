! ----------------------------------------------------------------------- !
! Quaternion Operations Module                                            !
! ----------------------------------------------------------------------- !
MODULE QuaternionOperations

! Uses one module: global variables
USE GlobalVariables

CONTAINS

! *********************************************************************************************** !
!              This subroutine rotates a vector using Hamilton's quaternion product.              !
! *********************************************************************************************** !
SUBROUTINE VectorRotation( PointVector, RotationQuaternion, RotatedVector )

IMPLICIT NONE

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 3 )   :: PointVector         ! Point vector [XYZ]
REAL( Kind= Real64 ), DIMENSION( 3 )   :: RotatedVector       ! Rotated vector [XYZ]
REAL( Kind= Real64 ), DIMENSION( 0:3 ) :: RotationQuaternion  ! Rotation quaternion [WXYZ]
REAL( Kind= Real64 ), DIMENSION( 0:3 ) :: ConjugateQuaternion ! Conjugate quaternion [WXYZ]
REAL( Kind= Real64 ), DIMENSION( 0:3 ) :: RotatedQuaternion   ! Rotated quaternion [WXYZ]
REAL( Kind= Real64 ), DIMENSION( 0:3 ) :: TempQuaternion      ! Temporary quaternion [WXYZ]

! Multiplication of quaternions to find the rotated quaternion
CALL QuaternionMultiplication( RotationQuaternion, [ 0.D0, PointVector ], TempQuaternion )

! Conjugate quaternion
ConjugateQuaternion = [ RotationQuaternion(0), -RotationQuaternion(1), -RotationQuaternion(2), -RotationQuaternion(3) ]

! Multiplication of the rotated quaternion and the conjugate of the rotation quaternion to find the rotated vector
CALL QuaternionMultiplication( TempQuaternion, ConjugateQuaternion, RotatedQuaternion )

! Rotated vector
RotatedVector(1:3) = RotatedQuaternion(1:3)

RETURN

END SUBROUTINE VectorRotation

! *********************************************************************************************** !
!        This subroutine takes a rotation quaternion (qm) and combines it with a randomly         !
!        generated quaternion (qr) through quaternion multiplication, creating a randomly         !
!                                composed rotation quaternion (qn)                                !
! *********************************************************************************************** !
SUBROUTINE QuaternionCombination( ReferenceUnitQuaternion, ComposedUnitQuaternion, Angle )

IMPLICIT NONE

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 ) :: Angle ! Maximum angular displacement

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 0:3 ) :: ReferenceUnitQuaternion ! Reference rotation quaternion
REAL( Kind= Real64 ), DIMENSION( 0:3 ) :: RandomUnitQuaternion    ! Random quaternion
REAL( Kind= Real64 ), DIMENSION( 0:3 ) :: ComposedUnitQuaternion  ! Composed rotation quaternion

! Random quaternion generator
CALL RandomQuaternionGenerator( RandomUnitQuaternion, Angle )

! Quaternion multiplication (composed rotation)
CALL QuaternionMultiplication( RandomUnitQuaternion, ReferenceUnitQuaternion, ComposedUnitQuaternion )

RETURN

END SUBROUTINE QuaternionCombination

! *********************************************************************************************** !
!        This subroutine generates a random quaternion from a random angle and random axis        !
! *********************************************************************************************** !
SUBROUTINE RandomQuaternionGenerator( RandomUnitQuaternion, Angle )

IMPLICIT NONE

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 ) :: Angle       ! Maximum angular displacement
REAL( Kind= Real64 ) :: RandomAngle ! Random angle

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 3 )   :: RandomUnitVector     ! Random vector
REAL( Kind= Real64 ), DIMENSION( 0:3 ) :: RandomUnitQuaternion ! Random unit quaternion

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

! *********************************************************************************************** !
!            This subroutine generates a random vector on the surface of a unit sphere            !
!           See Allen and Tildesley, 2nd Edition (2017), page 514 for more information.           !
!        (Routine 'maths_module.f90' of Code A.1. in Marsaglia, Ann. Math. Statist., 1972)        !
! *********************************************************************************************** !   
SUBROUTINE RandomVectorGenerator( RandomUnitVector )

IMPLICIT NONE

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 ) :: Alpha, Beta, Zeta ! Random numbers

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 3 ) :: RandomUnitVector ! Random vector

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

! *********************************************************************************************** !
!            This subroutine creates a composed rotation quaternion by multiplying the            !
!                         reference quaternion with a random quaternion                           !
! *********************************************************************************************** !
SUBROUTINE QuaternionMultiplication( RandomUnitQuaternion, ReferenceUnitQuaternion, ComposedUnitQuaternion )

IMPLICIT NONE

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( 0:3 ) :: RandomUnitQuaternion    ! Random quaternion
REAL( Kind= Real64 ), DIMENSION( 0:3 ) :: ReferenceUnitQuaternion ! Reference quaternion
REAL( Kind= Real64 ), DIMENSION( 0:3 ) :: ComposedUnitQuaternion  ! Composed quaternion

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

END MODULE QuaternionOperations