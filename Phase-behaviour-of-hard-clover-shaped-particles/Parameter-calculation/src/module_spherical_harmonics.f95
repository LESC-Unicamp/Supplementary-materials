MODULE Spherical_Harmonics_LM

! Uses one module: global variables
USE GlobalVariables

CONTAINS

! ********************************************************************************************************************************* !
! This subroutine calculates the spherical harmonics (Yℓm) of ℓth-degree (<lQuantumNumber>) and mth-order (<mQuantumNumber>).       !
!                                                                                                                                   !
! ● The spherical harmonics are calculated as a function of the polar angle (<PolarAngle>) and azimuthal angle (<AzimuthalAngle>)   !
!   in a spherical coordinate system. The real part of the spherical harmonics is stored in <RealTerm>, and the imaginary part is   !
!   stored in <ImaginaryTerm>.                                                                                                      !
!                                                                                                                                   !
! ● The spherical harmonics are given by the formula:                                                                               !
!                                                                                                                                   !
!       Yℓm(θ,φ) = (-1)ᵐ × √[ (2ℓ + 1) / (4π) × (ℓ - m)! / (ℓ + m)! ] × Pℓm(cos(θ)) × eⁱᵐᶲ,                                         !
!                                                                                                                                   !
!   where Pℓm(x) is the associated Legendre polynomial of ℓth-degree and mth-order, and eⁱᵐᶲ is expressed through the Euler's       !
!   formula: eⁱᵐᶲ = cos(mφ) + i × sin(mφ).                                                                                          !
! ********************************************************************************************************************************* !
SUBROUTINE Spherical_Harmonics( lQuantumNumber, mQuantumNumber, PolarAngle, AzimuthalAngle, RealTerm, ImaginaryTerm )

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 ), INTENT( IN ) :: mQuantumNumber ! Magnetic quantum number (m)
INTEGER( Kind= Int64 ), INTENT( IN ) :: lQuantumNumber ! Azimuthal quantum number (ℓ)

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 )                :: FactorialTerm      ! Factorial term
REAL( Kind= Real64 )                :: Prefactor          ! Prefactor of the spherical harmonics
REAL( Kind= Real64 )                :: SphericalHarmonics ! Spherical harmonics of ℓth-degree and mth-order
REAL( Kind= Real64 ), INTENT( OUT ) :: RealTerm           ! Real part of the spherical harmonics
REAL( Kind= Real64 ), INTENT( OUT ) :: ImaginaryTerm      ! Imaginary part of the spherical harmonics
REAL( Kind= Real64 ), INTENT( IN )  :: PolarAngle         ! Polar angle
REAL( Kind= Real64 ), INTENT( IN )  :: AzimuthalAngle     ! Azimuthal angle

! Check if the degree is non-negative
IF( lQuantumNumber < 0 ) THEN
  WRITE( *, "(G0)" ) " "
  WRITE( *, "(G0)" ) " "
  WRITE( *, "(A,I0,A)" ) "The azimuthal quantum number (", lQuantumNumber, ") must be non-negative!"
  STOP
END IF

! Check if the order is within the degree
IF( ABS( mQuantumNumber ) > lQuantumNumber ) THEN
  WRITE( *, "(G0)" ) " "
  WRITE( *, "(G0)" ) " "
  WRITE( *, "(A,I0,A,I0,A)" ) "The magnetic quantum number (", mQuantumNumber, &
  &                           ") must be within the range of the azimuthal quantum number (±", lQuantumNumber, ")!"
  STOP
END IF

! Check if degree is less than six
IF( lQuantumNumber > 12 ) THEN
  WRITE( *, "(G0)" ) " "
  WRITE( *, "(G0)" ) " "
  WRITE( *, "(A,I0,A)" ) "The azimuthal quantum number (", lQuantumNumber, ") of the spherical harmonics cannot exceed six!"
  STOP
END IF

! Prefactor
FactorialTerm = Factorial( lQuantumNumber - mQuantumNumber ) / Factorial( lQuantumNumber + mQuantumNumber )
Prefactor = DSQRT( FactorialTerm * ( 2.D0 * DBLE( lQuantumNumber ) + 1.D0 ) / ( 4.D0 * cPi ) )
Prefactor = ( - 1.D0 ) ** ( ABS( mQuantumNumber ) ) * Prefactor

! Spherical harmonics of ℓth-degree and mth-order
SphericalHarmonics = Prefactor * AssociatedLegendrePolynomials( lQuantumNumber, mQuantumNumber, DCOS( PolarAngle ) )

! Real and imaginary parts of the spherical harmonics
RealTerm = SphericalHarmonics * DCOS( mQuantumNumber * AzimuthalAngle )
ImaginaryTerm = SphericalHarmonics * DSIN( mQuantumNumber * AzimuthalAngle )

RETURN

END SUBROUTINE Spherical_Harmonics

! ********************************************************************************************************************************* !
! This subroutine calculates the associated Legendre polynomials of ℓth-degree (<pDegree>) and mth-order (<pOrder>) of the cossine  !
! of the polar angle (<Argument>).                                                                                                  !
!                                                                                                                                   !
! ● The associated Legendre polynomials are calculated as:                                                                          !
!                                                                                                                                   !
!       Pℓm(x) = √(1 - x²)ᵐ × dᵐ/dxᵐ [ Pℓ(x) ],                                                                                      !
!                                                                                                                                   !
!   where Pℓ(x) is the Legendre polynomial of ℓth-degree, and dᵐ/dxᵐ [ Pℓ(x) ] is the mth derivative of the Legendre polynomial.    !
!                                                                                                                                   !
!   If the order is negative, the associated Legendre polynomial is calculated as:                                                  !
!                                                                                                                                   !
!       Pℓm(x) = (-1)ᵐ × (ℓ - |m|)! / (ℓ + |m|)! × √(1 - x²)ᵐ × dᵐ/dxᵐ [ Pℓ(x) ]                                                     !
! ********************************************************************************************************************************* !
DOUBLE PRECISION FUNCTION AssociatedLegendrePolynomials( pDegree, pOrder, Argument )

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 ), INTENT( IN ) :: pDegree ! Degree (l) of the associated Legendre polynomial
INTEGER( Kind= Int64 ), INTENT( IN ) :: pOrder  ! Order (m) of the associated Legendre polynomial

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 )               :: dLegendrePolynomial ! Derivative of mth-order of the Legendre polynomial of ℓth-degree
REAL( Kind= Real64 )               :: FactorialTerm       ! Term of the factorial
REAL( Kind= Real64 ), INTENT( IN ) :: Argument            ! Argument (x) of the associated Legendre polynomial

! Derivatives of mth-order of the Legendre polynomial of ℓth-degree
IF( pDegree == 0 ) THEN ! Zeroth degree
  ! Zeroth derivative
  dLegendrePolynomial = 1.D0
ELSE IF( pDegree == 1 ) THEN ! First degree
  IF( ABS( pOrder ) == 0 ) THEN
    ! Zeroth derivative
    dLegendrePolynomial = Argument
  ELSE IF( ABS( pOrder ) == 1 ) THEN
    ! First derivative
    dLegendrePolynomial = 1.D0
  END IF
ELSE IF( pDegree == 2 ) THEN ! Second degree
  IF( ABS( pOrder ) == 0 ) THEN
    ! Zeroth derivative
    dLegendrePolynomial = 0.5D0 * ( 3.D0 * Argument ** 2 - 1.D0 )
  ELSE IF( ABS( pOrder ) == 1 ) THEN
    ! First derivative
    dLegendrePolynomial = 3.D0 * Argument
  ELSE IF( ABS( pOrder ) == 2 ) THEN
    ! Second derivative
    dLegendrePolynomial = 3.D0
  END IF
ELSE IF( pDegree == 3 ) THEN ! Third degree
  IF( ABS( pOrder ) == 0 ) THEN
    ! Zeroth derivative
    dLegendrePolynomial = 0.5D0 * ( 5.D0 * Argument ** 3 - 3.D0 * Argument )
  ELSE IF( ABS( pOrder ) == 1 ) THEN
    ! First derivative
    dLegendrePolynomial = 1.5D0 * ( 5.D0 * Argument ** 2 - 1.D0 )
  ELSE IF( ABS( pOrder ) == 2 ) THEN
    ! Second derivative
    dLegendrePolynomial = 15.D0 * Argument
  ELSE IF( ABS( pOrder ) == 3 ) THEN
    dLegendrePolynomial = 15.D0
  END IF
ELSE IF( pDegree == 4 ) THEN ! Fourth degree
  IF( ABS( pOrder ) == 0 ) THEN
    ! Zeroth derivative
    dLegendrePolynomial = 0.125D0 * ( 35.D0 * Argument ** 4 - 30.D0 * Argument ** 2 + 3.D0 )
  ELSE IF( ABS( pOrder ) == 1 ) THEN
    ! First derivative
    dLegendrePolynomial = 0.5D0 * ( 35.D0 * Argument ** 3 - 15.D0 * Argument )
  ELSE IF( ABS( pOrder ) == 2 ) THEN
    ! Second derivative
    dLegendrePolynomial = 0.5D0 * ( 105.D0 * Argument ** 2 - 15.D0 )
  ELSE IF( ABS( pOrder ) == 3 ) THEN
    ! Third derivative
    dLegendrePolynomial = 105.D0 * Argument
  ELSE IF( ABS( pOrder ) == 4 ) THEN
    ! Fourth derivative
    dLegendrePolynomial = 105.D0
  END IF
ELSE IF( pDegree == 5 ) THEN ! Fifth degree
  IF( ABS( pOrder ) == 0 ) THEN
    ! Zeroth derivative
    dLegendrePolynomial = 0.125D0 * ( 63.D0 * Argument ** 5 - 70.D0 * Argument ** 3 + 15.D0 * Argument )
  ELSE IF( ABS( pOrder ) == 1 ) THEN
    ! First derivative
    dLegendrePolynomial = 0.625D0 * ( 63.D0 * Argument ** 4 - 42.D0 * Argument ** 2 + 3.D0 )
  ELSE IF( ABS( pOrder ) == 2 ) THEN
    ! Second derivative
    dLegendrePolynomial = 2.5D0 * ( 63.D0 * Argument ** 3 - 21.D0 * Argument )
  ELSE IF( ABS( pOrder ) == 3 ) THEN
    ! Third derivative
    dLegendrePolynomial = 2.5D0 * ( 189.D0 * Argument ** 2 - 21.D0 )
  ELSE IF( ABS( pOrder ) == 4 ) THEN
    ! Fourth derivative
    dLegendrePolynomial = 945.D0 * Argument
  ELSE IF( ABS( pOrder ) == 5 ) THEN
    ! Fifth derivative
    dLegendrePolynomial = 945.D0
  END IF
ELSE IF( pDegree == 6 ) THEN ! Sixth degree
  IF( ABS( pOrder ) == 0 ) THEN
    ! Zeroth derivative
    dLegendrePolynomial = 0.0625D0 * ( 231.D0 * Argument ** 6 - 315.D0 * Argument ** 4 + 105.D0 * Argument ** 2 - 5.D0 )
  ELSE IF( ABS( pOrder ) == 1 ) THEN
    ! First derivative
    dLegendrePolynomial = 0.375D0 * ( 231.D0 * Argument ** 5 - 210.D0 * Argument ** 3 + 35.D0 * Argument )
  ELSE IF( ABS( pOrder ) == 2 ) THEN
    ! Second derivative
    dLegendrePolynomial = 1.875D0 * ( 231.D0 * Argument ** 4 - 126.D0 * Argument ** 2 + 7.D0 )
  ELSE IF( ABS( pOrder ) == 3 ) THEN
    ! Third derivative
    dLegendrePolynomial = 7.5D0 * ( 231.D0 * Argument ** 3 - 63.D0 * Argument )
  ELSE IF( ABS( pOrder ) == 4 ) THEN
    ! Fourth derivative
    dLegendrePolynomial = 22.5D0 * ( 231.D0 * Argument ** 2 - 21.D0 )
  ELSE IF( ABS( pOrder ) == 5 ) THEN
    ! Fifth derivative
    dLegendrePolynomial = 10395.D0 * Argument
  ELSE IF( ABS( pOrder ) == 6 ) THEN
    ! Sixth derivative
    dLegendrePolynomial = 10395.D0
  END IF
ELSE IF( pDegree == 7 ) THEN ! Seventh degree
  IF( ABS( pOrder ) == 0 ) THEN
    ! Zeroth derivative
    dLegendrePolynomial = 0.0625D0 * ( 429.D0 * Argument ** 7 - 693.D0 * Argument ** 5 + 315.D0 * Argument ** 3 - 35.D0 * Argument )
  ELSE IF( ABS( pOrder ) == 1 ) THEN
    ! First derivative
    dLegendrePolynomial = 0.4375D0 * ( 429.D0 * Argument ** 6 - 495.D0 * Argument ** 4 + 135.D0 * Argument ** 2 - 5.D0 )
  ELSE IF( ABS( pOrder ) == 2 ) THEN
    ! Second derivative
    dLegendrePolynomial = 2.625D0 * ( 429.D0 * Argument ** 5 - 330.D0 * Argument ** 3 + 45.D0 * Argument )
  ELSE IF( ABS( pOrder ) == 3 ) THEN
    ! Third derivative
    dLegendrePolynomial = 13.125D0 * ( 429.D0 * Argument ** 4 - 198.D0 * Argument ** 2 + 9.D0 )
  ELSE IF( ABS( pOrder ) == 4 ) THEN
    ! Fourth derivative
    dLegendrePolynomial = 52.5D0 * ( 429.D0 * Argument ** 3 - 99.D0 * Argument )
  ELSE IF( ABS( pOrder ) == 5 ) THEN
    ! Fifth derivative
    dLegendrePolynomial = 157.5D0 * ( 429.D0 * Argument ** 2 - 33.D0 )
  ELSE IF( ABS( pOrder ) == 6 ) THEN
    ! Sixth derivative
    dLegendrePolynomial = 135135.D0 * Argument
  ELSE IF( ABS( pOrder ) == 7 ) THEN
    ! Seventh derivative
    dLegendrePolynomial = 135135.D0
  END IF
ELSE IF( pDegree == 8 ) THEN ! Eighth degree
  IF( ABS( pOrder ) == 0 ) THEN
    ! Zeroth derivative
    dLegendrePolynomial = 0.0078125D0 * ( 6435.D0 * Argument ** 8 - 12012.D0 * Argument ** 6 + 6930.D0 * Argument ** 4 - 1260.D0 * &
    &                     Argument ** 2 + 35.D0 )
  ELSE IF( ABS( pOrder ) == 1 ) THEN
    ! First derivative
    dLegendrePolynomial = 0.0625D0 * ( 6435.D0 * Argument ** 7 - 9009.D0 * Argument ** 5 + 3465.D0 * Argument ** 3 - 315.D0 * &
    &                     Argument )
  ELSE IF( ABS( pOrder ) == 2 ) THEN
    ! Second derivative
    dLegendrePolynomial = 0.4375D0 * ( 6435.D0 * Argument ** 6 - 6435.D0 * Argument ** 4 + 1485.D0 * Argument ** 2 - 45.D0 )
  ELSE IF( ABS( pOrder ) == 3 ) THEN
    ! Third derivative
    dLegendrePolynomial = 2.625D0 * ( 6435.D0 * Argument ** 5 - 4290.D0 * Argument ** 3 + 495.D0 * Argument )
  ELSE IF( ABS( pOrder ) == 4 ) THEN
    ! Fourth derivative
    dLegendrePolynomial = 13.125D0 * ( 6435.D0 * Argument ** 4 - 2574.D0 * Argument ** 2 + 99.D0 )
  ELSE IF( ABS( pOrder ) == 5 ) THEN
    ! Fifth derivative
    dLegendrePolynomial = 52.5D0 * ( 6435.D0 * Argument ** 3 - 1287.D0 * Argument )
  ELSE IF( ABS( pOrder ) == 6 ) THEN
    ! Sixth derivative
    dLegendrePolynomial = 157.5D0 * ( 6435.D0 * Argument ** 2 - 429.D0 )
  ELSE IF( ABS( pOrder ) == 7 ) THEN
    ! Seventh derivative
    dLegendrePolynomial = 315.D0 * ( 6435.D0 * Argument )
  ELSE IF( ABS( pOrder ) == 8 ) THEN
    ! Eighth derivative
    dLegendrePolynomial = 2027025.D0
  END IF
ELSE IF( pDegree == 9 ) THEN ! Ninth degree
  IF( ABS( pOrder ) == 0 ) THEN
    ! Zeroth derivative
    dLegendrePolynomial = 0.0078125D0 * ( 12155.D0 * Argument ** 9 - 25740.D0 * Argument ** 7 + 18018.D0 * Argument ** 5 - &
    &                     4620.D0 * Argument ** 3 + 315.D0 * Argument )
  ELSE IF( ABS( pOrder ) == 1 ) THEN
    ! First derivative
    dLegendrePolynomial = 0.0703125D0 * ( 12155.D0 * Argument ** 8 - 20020.D0 * Argument ** 6 + 10010.D0 * Argument ** 4 - &
    &                     1540.D0 * Argument ** 2 + 35.D0 )
  ELSE IF( ABS( pOrder ) == 2 ) THEN
    ! Second derivative
    dLegendrePolynomial = 0.5625D0 * ( 12155.D0 * Argument ** 7 - 15015.D0 * Argument ** 5 + 5005.D0 * Argument ** 3 - &
    &                     385.D0 * Argument )
  ELSE IF( ABS( pOrder ) == 3 ) THEN
    ! Third derivative
    dLegendrePolynomial = 3.9375D0 * ( 12155.D0 * Argument ** 6 - 10725.D0 * Argument ** 4 + 2145.D0 * Argument ** 2 - &
    &                     55.D0 )
  ELSE IF( ABS( pOrder ) == 4 ) THEN
    ! Fourth derivative
    dLegendrePolynomial = 23.625D0 * ( 12155.D0 * Argument ** 5 - 7150.D0 * Argument ** 3 + 715.D0 * Argument )
  ELSE IF( ABS( pOrder ) == 5 ) THEN
    ! Fifth derivative
    dLegendrePolynomial = 118.125D0 * ( 12155.D0 * Argument ** 4 - 4290.D0 * Argument ** 2 + 143.D0 )
  ELSE IF( ABS( pOrder ) == 6 ) THEN
    ! Sixth derivative
    dLegendrePolynomial = 472.5D0 * ( 12155.D0 * Argument ** 3 - 2145.D0 * Argument )
  ELSE IF( ABS( pOrder ) == 7 ) THEN
    ! Seventh derivative
    dLegendrePolynomial = 1417.5D0 * ( 12155.D0 * Argument ** 2 - 715.D0 )
  ELSE IF( ABS( pOrder ) == 8 ) THEN
    ! Eighth derivative
    dLegendrePolynomial = 2835.D0 * ( 12155.D0 * Argument )
  ELSE IF( ABS( pOrder ) == 9 ) THEN
    ! Ninth derivative
    dLegendrePolynomial = 34459425.D0
  END IF
ELSE IF( pDegree == 10 ) THEN ! Tenth degree
  IF( ABS( pOrder ) == 0 ) THEN
    ! Zeroth derivative
    dLegendrePolynomial = 0.00390625D0 * ( 46189.D0 * Argument ** 10 - 109395.D0 * Argument ** 8 + 90090.D0 * Argument ** 6 - &
    &                     30030.D0 * Argument ** 4 + 3465.D0 * Argument ** 2 - 63.D0 )
  ELSE IF( ABS( pOrder ) == 1 ) THEN
    ! First derivative
    dLegendrePolynomial = 0.0390625D0 * ( 46189.D0 * Argument ** 9 - 87516.D0 * Argument ** 7 + 54054.D0 * Argument ** 5 - &
    &                     12012.D0 * Argument ** 3 + 693.D0 * Argument )
  ELSE IF( ABS( pOrder ) == 2 ) THEN
    ! Second derivative
    dLegendrePolynomial = 0.3515625D0 * ( 46189.D0 * Argument ** 8 - 68068.D0 * Argument ** 6 + 30030.D0 * Argument ** 4 - &
    &                     4004.D0 * Argument ** 2 + 77.D0 )
  ELSE IF( ABS( pOrder ) == 3 ) THEN
    ! Third derivative
    dLegendrePolynomial = 2.8125D0 * ( 46189.D0 * Argument ** 7 - 51051.D0 * Argument ** 5 + 15015.D0 * Argument ** 3 - &
    &                     1001.D0 * Argument )
  ELSE IF( ABS( pOrder ) == 4 ) THEN
    ! Fourth derivative
    dLegendrePolynomial = 19.6875D0 * ( 46189.D0 * Argument ** 6 - 36465.D0 * Argument ** 4 + 6435.D0 * Argument ** 2 - &
    &                     143.D0 )
  ELSE IF( ABS( pOrder ) == 5 ) THEN
    ! Fifth derivative
    dLegendrePolynomial = 118.125D0 * ( 46189.D0 * Argument ** 5 - 24310.D0 * Argument ** 3 + 2145.D0 * Argument )
  ELSE IF( ABS( pOrder ) == 6 ) THEN
    ! Sixth derivative
    dLegendrePolynomial = 590.625D0 * ( 46189.D0 * Argument ** 4 - 14586.D0 * Argument ** 2 + 429.D0 )
  ELSE IF( ABS( pOrder ) == 7 ) THEN
    ! Seventh derivative
    dLegendrePolynomial = 2362.5D0 * ( 46189.D0 * Argument ** 3 - 7293.D0 * Argument )
  ELSE IF( ABS( pOrder ) == 8 ) THEN
    ! Eighth derivative
    dLegendrePolynomial = 7087.5D0 * ( 46189.D0 * Argument ** 2 - 2431.D0 )
  ELSE IF( ABS( pOrder ) == 9 ) THEN
    ! Ninth derivative
    dLegendrePolynomial = 14175.D0 * ( 46189.D0 * Argument )
  ELSE IF( ABS( pOrder ) == 10 ) THEN
    ! Tenth derivative
    dLegendrePolynomial = 654729075.D0
  END IF
ELSE IF( pDegree == 11 ) THEN ! Eleventh degree
  IF( ABS( pOrder ) == 0 ) THEN
    ! Zeroth derivative
    dLegendrePolynomial = 0.00390625D0 * ( 88179.D0 * Argument ** 11 - 230945.D0 * Argument ** 9 + 218790.D0 * Argument ** 7 - &
    &                     90090.D0 * Argument ** 5 + 15015.D0 * Argument ** 3 - 693.D0 * Argument )
  ELSE IF( ABS( pOrder ) == 1 ) THEN
    ! First derivative
    dLegendrePolynomial = 0.04296875D0 * ( 88179.D0 * Argument ** 10 - 188955.D0 * Argument ** 8 + 139230.D0 * Argument ** 6 - &
    &                     40950.D0 * Argument ** 4 + 4095.D0 * Argument ** 2 - 63.D0 )
  ELSE IF( ABS( pOrder ) == 2 ) THEN
    ! Second derivative
    dLegendrePolynomial = 0.4296875D0 * ( 88179.D0 * Argument ** 9 - 151164.D0 * Argument ** 7 + 83538.D0 * Argument ** 5 - &
    &                     16380.D0 * Argument ** 3 + 819.D0 * Argument )
  ELSE IF( ABS( pOrder ) == 3 ) THEN
    ! Third derivative
    dLegendrePolynomial = 3.8671875D0 * ( 88179.D0 * Argument ** 8 - 117572.D0 * Argument ** 6 + 46410.D0 * Argument ** 4 - &
    &                     5460.D0 * Argument ** 2 + 91.D0 )
  ELSE IF( ABS( pOrder ) == 4 ) THEN
    ! Fourth derivative
    dLegendrePolynomial = 30.9375D0 * ( 88179.D0 * Argument ** 7 - 88179.D0 * Argument ** 5 + 23205.D0 * Argument ** 3 - &
    &                     1365.D0 * Argument )
  ELSE IF( ABS( pOrder ) == 5 ) THEN
    ! Fifth derivative
    dLegendrePolynomial = 216.5625D0 * ( 88179.D0 * Argument ** 6 - 62985.D0 * Argument ** 4 + 9945.D0 * Argument ** 2 - &
    &                     195.D0 )
  ELSE IF( ABS( pOrder ) == 6 ) THEN
    ! Sixth derivative
    dLegendrePolynomial = 1299.375D0 * ( 88179.D0 * Argument ** 5 - 41990.D0 * Argument ** 3 + 3315.D0 * Argument )
  ELSE IF( ABS( pOrder ) == 7 ) THEN
    ! Seventh derivative
    dLegendrePolynomial = 6496.875D0 * ( 88179.D0 * Argument ** 4 - 25194.D0 * Argument ** 2 + 663.D0 )
  ELSE IF( ABS( pOrder ) == 8 ) THEN
    ! Eighth derivative
    dLegendrePolynomial = 25987.5D0 * ( 88179.D0 * Argument ** 3 - 12597.D0 * Argument )
  ELSE IF( ABS( pOrder ) == 9 ) THEN
    ! Ninth derivative
    dLegendrePolynomial = 77962.5D0 * ( 88179.D0 * Argument ** 2 - 4199.D0 )
  ELSE IF( ABS( pOrder ) == 10 ) THEN
    ! Tenth derivative
    dLegendrePolynomial = 155925.D0 * ( 88179.D0 * Argument )
  ELSE IF( ABS( pOrder ) == 11 ) THEN
    ! Eleventh derivative
    dLegendrePolynomial = 13749310575.D0
  END IF
ELSE IF( pDegree == 12 ) THEN ! Twelfth degree
  IF( ABS( pOrder ) == 0 ) THEN
    ! Zeroth derivative
    dLegendrePolynomial = 0.0009765625D0 * ( 676039.D0 * Argument ** 12 - 1939938.D0 * Argument ** 10 + 2078505.D0 * Argument ** 8 &
    &                     - 1021020.D0 * Argument ** 6 + 225225.D0 * Argument ** 4 - 18018.D0 * Argument ** 2 + 231.D0 )
  ELSE IF( ABS( pOrder ) == 1 ) THEN
    ! First derivative
    dLegendrePolynomial = 0.01171875D0 * ( 676039.D0 * Argument ** 11 - 1616615.D0 * Argument ** 9 + 1385670.D0 * Argument ** 7 - &
    &                     510510.D0 * Argument ** 5 + 75075.D0 * Argument ** 3 - 3003.D0 * Argument )
  ELSE IF( ABS( pOrder ) == 2 ) THEN
    ! Second derivative
    dLegendrePolynomial = 0.12890625D0 * ( 676039.D0 * Argument ** 10 - 1322685.D0 * Argument ** 8 + 881790.D0 * Argument ** 6 - &
    &                     232050.D0 * Argument ** 4 + 20475.D0 * Argument ** 2 - 273.D0 )
  ELSE IF( ABS( pOrder ) == 3 ) THEN
    ! Third derivative
    dLegendrePolynomial = 1.2890625D0 * ( 676039.D0 * Argument ** 9 - 1058148.D0 * Argument ** 7 + 529074.D0 * Argument ** 5 - &
    &                     92820.D0 * Argument ** 3 + 4095.D0 * Argument )
  ELSE IF( ABS( pOrder ) == 4 ) THEN
    ! Fourth derivative
    dLegendrePolynomial = 11.6015625D0 * ( 676039.D0 * Argument ** 8 - 823004.D0 * Argument ** 6 + 293930.D0 * Argument ** 4 - &
    &                     30940.D0 * Argument ** 2 + 455.D0 )
  ELSE IF( ABS( pOrder ) == 5 ) THEN
    ! Fifth derivative
    dLegendrePolynomial = 92.8125D0 * ( 676039.D0 * Argument ** 7 - 617253.0D0 * Argument ** 5 + 146965.D0 * Argument ** 3 - &
    &                     7735.D0 * Argument )
  ELSE IF( ABS( pOrder ) == 6 ) THEN
    ! Sixth derivative
    dLegendrePolynomial = 649.6875D0 * ( 676039.D0 * Argument ** 6 - 440895.D0 * Argument ** 4 + 62985.D0 * Argument ** 2 - &
    &                     1105.D0 )
  ELSE IF( ABS( pOrder ) == 7 ) THEN
    ! Seventh derivative
    dLegendrePolynomial = 3898.125D0 * ( 676039.D0 * Argument ** 5 - 293930.D0 * Argument ** 3 + 20995.D0 * Argument )
  ELSE IF( ABS( pOrder ) == 8 ) THEN
    ! Eighth derivative
    dLegendrePolynomial = 19490.625D0 * ( 676039.D0 * Argument ** 4 - 176358.D0 * Argument ** 2 + 4199.D0 )
  ELSE IF( ABS( pOrder ) == 9 ) THEN
    ! Ninth derivative
    dLegendrePolynomial = 77962.5D0 * ( 676039.D0 * Argument ** 3 - 88179.D0 * Argument )
  ELSE IF( ABS( pOrder ) == 10 ) THEN
    ! Tenth derivative
    dLegendrePolynomial = 233887.5D0 * ( 676039.D0 * Argument ** 2 - 29393.D0 )
  ELSE IF( ABS( pOrder ) == 11 ) THEN
    ! Eleventh derivative
    dLegendrePolynomial = 467775.D0 * ( 676039.D0 * Argument )
  ELSE IF( ABS( pOrder ) == 12 ) THEN
    ! Twelfth derivative
    dLegendrePolynomial = 316234143225.D0
  END IF
END IF

! Associated Legendre polynomial
AssociatedLegendrePolynomials = DSQRT( 1.D0 - Argument * Argument ) ** ( ABS( pOrder ) ) * dLegendrePolynomial
IF( pOrder < 0 ) THEN ! Correct if the magnetic quantum number is negative
  ! Factorial term
  FactorialTerm = Factorial( pDegree - ABS( pOrder ) ) / Factorial( pDegree + ABS( pOrder ) )
  ! Associated Legendre polynomial
  AssociatedLegendrePolynomials = ( - 1.D0 ) ** ( ABS( pOrder ) ) * FactorialTerm * AssociatedLegendrePolynomials
END IF

RETURN

END FUNCTION AssociatedLegendrePolynomials

! ********************************************************************************************************************************* !
! This subroutine calculates the factorial of a number (<Argument>).                                                                !
! ********************************************************************************************************************************* !
DOUBLE PRECISION FUNCTION Factorial( Argument )

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 )               :: Counter  ! Counter
INTEGER( Kind= Int64 ), INTENT( IN ) :: Argument ! Argument of the factorial

! Calculate the factorial
Factorial = 1.D0
DO Counter = 1, Argument
  Factorial = Factorial * DBLE( Counter )
END DO

RETURN

END FUNCTION Factorial

END MODULE Spherical_Harmonics_LM