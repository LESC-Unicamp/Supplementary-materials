! ************************************************************************************************ !
!                                      GOLDEN-SECTION SEARCH                                       !
! ************************************************************************************************ !
!            This subroutine uses the golden-section search determine the local minimum            !
!        (critical point) of a function f(x) in the interval [a,b]. It can also be used to         !
!        calculate the local maximum in the same interval by making f(x) = -f(x). Finally,         !
!          it can also be used to calculate the inflection point in the interval [a,b] by          !
!       checking the critical points of the first derivative of f with respect to x, f'(x).        !
! ************************************************************************************************ !
! => AUTHOR:     Nathan Barros de Souza                                                            !
! => E-MAIL:     n264179@dac.unicamp.br                                                            !
! => SUPERVISOR: Luís Fernando Mercier Franco                                                      !
! ************************************************************************************************ !
!                       University of Campinas, Campinas - São Paulo, Brazil                       !
!                                         March 7th, 2024                                          !
! ************************************************************************************************ !
! Version 1.0.0                                                                                    !
! ************************************************************************************************ !

! ************************************************************************************************ !
!                             GOLDEN-SECTION SEARCH METHOD (MIXTURES)                              !
! ************************************************************************************************ !
SUBROUTINE Golden_Section_Search( mFraction, Temperature, LowerBound, UpperBound, OptimizationOption, CriticalPoint )

! Uses two modules: global variables and substances
USE GlobalVariables
USE Substances

IMPLICIT NONE
    
! ************************************************************************************************ !
! INTEGER VARIABLES                                                                                !
! ************************************************************************************************ !
INTEGER( Kind= Int64 ) :: OptimizationOption ! Optimization option (1 = inflexion point; 2 = local minimum; 3 = local maximum)

! ************************************************************************************************ !
! REAL VARIABLES                                                                                   !
! ************************************************************************************************ !
REAL( Kind= Real64 )                           :: LowerBound, UpperBound               ! Extrema points of interval
REAL( Kind= Real64 )                           :: aMidpointInterval, bMidpointInterval ! Midpoints of interval
REAL( Kind= Real64 )                           :: aFunctionValue, bFunctionValue       ! Function values at midpoints
REAL( Kind= Real64 )                           :: CriticalPoint                        ! Critical point
REAL( Kind= Real64 )                           :: Pressure                             ! Pressure
REAL( Kind= Real64 )                           :: Temperature                          ! Temperature
REAL( Kind= Real64 )                           :: IsothermalCompressibility            ! Isothermal compressibility
REAL( Kind= Real64 )                           :: dPressure_dDensity                   ! First derivative of pressure with respect to the molar density
REAL( Kind= Real64 )                           :: ThermalExpansionCoefficient          ! Thermal expansion coefficient
REAL( Kind= Real64 )                           :: CompressibilityFactor                ! Compressibility factor
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: mFraction                            ! Molar fraction of a component

! ************************************************************************************************ !
! REAL VARIABLES (PARAMETERS)                                                                      !
! ************************************************************************************************ !
REAL( Kind= Real64 ), PARAMETER :: GoldenRatio = 0.5D0 * ( DSQRT( 5.D0 ) + 1.D0 ) ! Golden ratio
REAL( Kind= Real64 ), PARAMETER :: Tolerance   = 1.D-10                           ! Numerical tolerance

! Probe point selection
aMidpointInterval = UpperBound - (UpperBound - LowerBound) / GoldenRatio
CALL Calculate_Pressure( mFraction, aMidpointInterval, Temperature, Pressure, IsothermalCompressibility, &
&                        ThermalExpansionCoefficient, CompressibilityFactor )
dPressure_dDensity = aMidpointInterval / IsothermalCompressibility

! Optimization option (inflection point, local minimum, or local maximum)
IF( OptimizationOption == 1 ) THEN
  ! Inflexion point
  aFunctionValue = dPressure_dDensity
ELSE IF( OptimizationOption == 2 ) THEN
  ! Local minimum
  aFunctionValue = Pressure
ELSE IF( OptimizationOption == 3 ) THEN
  ! Local maximum
  aFunctionValue = - Pressure
END IF

! Probe point selection
bMidpointInterval = LowerBound + (UpperBound - LowerBound) / GoldenRatio
CALL Calculate_Pressure( mFraction, bMidpointInterval, Temperature, Pressure, IsothermalCompressibility, &
&                        ThermalExpansionCoefficient, CompressibilityFactor )
dPressure_dDensity = bMidpointInterval / IsothermalCompressibility

! Optimization option (inflection point, local minimum, or local maximum)
IF( OptimizationOption == 1 ) THEN
  ! Inflexion point
  bFunctionValue = dPressure_dDensity
ELSE IF( OptimizationOption == 2 ) THEN
  ! Local minimum
  bFunctionValue = Pressure
ELSE IF( OptimizationOption == 3 ) THEN
  ! Local maximum
  bFunctionValue = - Pressure
END IF

! Criterion
DO WHILE( DABS( (aMidpointInterval - bMidpointInterval) / aMidpointInterval ) >= Tolerance )

  IF( aFunctionValue < bFunctionValue ) THEN

    ! Swap intervals
    UpperBound = bMidpointInterval
    bMidpointInterval = aMidpointInterval
    bFunctionValue = aFunctionValue

    ! Probe point selection
    aMidpointInterval = UpperBound - (UpperBound - LowerBound) / GoldenRatio
    CALL Calculate_Pressure( mFraction, aMidpointInterval, Temperature, Pressure, IsothermalCompressibility, &
    &                        ThermalExpansionCoefficient, CompressibilityFactor )
    dPressure_dDensity = aMidpointInterval / IsothermalCompressibility

    ! Optimization option (inflection point, local minimum, or local maximum)
    IF( OptimizationOption == 1 ) THEN
      ! Inflexion point
      aFunctionValue = dPressure_dDensity
    ELSE IF( OptimizationOption == 2 ) THEN
      ! Local minimum
      aFunctionValue = Pressure
    ELSE IF( OptimizationOption == 3 ) THEN
      ! Local maximum
      aFunctionValue = - Pressure
    END IF

  ELSE IF( aFunctionValue >= bFunctionValue ) THEN

    ! Swap intervals
    LowerBound = aMidpointInterval
    aMidpointInterval = bMidpointInterval
    aFunctionValue = bFunctionValue

    ! Probe point selection
    bMidpointInterval = LowerBound + (UpperBound - LowerBound) / GoldenRatio
    CALL Calculate_Pressure( mFraction, bMidpointInterval, Temperature, Pressure, IsothermalCompressibility, &
    &                        ThermalExpansionCoefficient, CompressibilityFactor )
    dPressure_dDensity = bMidpointInterval / IsothermalCompressibility

    ! Optimization option (inflection point, local minimum, or local maximum)
    IF( OptimizationOption == 1 ) THEN
      ! Inflexion point
      bFunctionValue = dPressure_dDensity
    ELSE IF( OptimizationOption == 2 ) THEN
      ! Local minimum
      bFunctionValue = Pressure
    ELSE IF( OptimizationOption == 3 ) THEN
      ! Local maximum
      bFunctionValue = - Pressure
    END IF

  END IF

END DO

! Density at critical point
CriticalPoint = 0.5D0 * (LowerBound + UpperBound)

RETURN
 
END SUBROUTINE Golden_Section_Search

! ************************************************************************************************ !
!                          GOLDEN-SECTION SEARCH METHOD (PURE COMPONENTS)                          !
! ************************************************************************************************ !
SUBROUTINE Golden_Section_Single_Component( cComponent, Temperature, LowerBound, UpperBound, OptimizationOption, CriticalPoint )

! Uses two modules: global variables and substances
USE GlobalVariables
USE Substances

IMPLICIT NONE
    
! ************************************************************************************************ !
! INTEGER VARIABLES                                                                                !
! ************************************************************************************************ !
INTEGER( Kind= Int64 ) :: cComponent         ! Component index
INTEGER( Kind= Int64 ) :: OptimizationOption ! Optimization option (1 = inflexion point; 2 = local minimum; 3 = local maximum)

! ************************************************************************************************ !
! REAL VARIABLES                                                                                   !
! ************************************************************************************************ !
REAL( Kind= Real64 ) :: LowerBound, UpperBound               ! Extrema points of interval
REAL( Kind= Real64 ) :: aMidpointInterval, bMidpointInterval ! Midpoints of interval
REAL( Kind= Real64 ) :: aFunctionValue, bFunctionValue       ! Function values at midpoints
REAL( Kind= Real64 ) :: CriticalPoint                        ! Critical point
REAL( Kind= Real64 ) :: Pressure                             ! Pressure
REAL( Kind= Real64 ) :: Temperature                          ! Temperature
REAL( Kind= Real64 ) :: IsothermalCompressibility            ! Isothermal compressibility
REAL( Kind= Real64 ) :: dPressure_dDensity                   ! First derivative of pressure with respect to the molar density
REAL( Kind= Real64 ) :: ThermalExpansionCoefficient          ! Thermal expansion coefficient
REAL( Kind= Real64 ) :: CompressibilityFactor                ! Compressibility factor

! ************************************************************************************************ !
! REAL VARIABLES (PARAMETERS)                                                                      !
! ************************************************************************************************ !
REAL( Kind= Real64 ), PARAMETER :: GoldenRatio = 0.5D0 * ( DSQRT( 5.D0 ) + 1.D0 ) ! Golden ratio
REAL( Kind= Real64 ), PARAMETER :: Tolerance   = 1.D-10                           ! Numerical tolerance

! Probe point selection
aMidpointInterval = UpperBound - ( UpperBound - LowerBound ) / GoldenRatio
CALL Calculate_Pressure_Single_Component( cComponent, aMidpointInterval, Temperature, Pressure, IsothermalCompressibility, &
&                                         ThermalExpansionCoefficient, CompressibilityFactor )
dPressure_dDensity = aMidpointInterval / IsothermalCompressibility

! Optimization option (inflection point, local minimum, or local maximum)
IF( OptimizationOption == 1 ) THEN
  ! Inflexion point
  aFunctionValue = dPressure_dDensity
ELSE IF( OptimizationOption == 2 ) THEN
  ! Local minimum
  aFunctionValue = Pressure
ELSE IF( OptimizationOption == 3 ) THEN
  ! Local maximum
  aFunctionValue = - Pressure
END IF

! Probe point selection
bMidpointInterval = LowerBound + (UpperBound - LowerBound) / GoldenRatio
CALL Calculate_Pressure_Single_Component( cComponent, bMidpointInterval, Temperature, Pressure, IsothermalCompressibility, &
&                                         ThermalExpansionCoefficient, CompressibilityFactor )
dPressure_dDensity = bMidpointInterval / IsothermalCompressibility

! Optimization option (inflection point, local minimum, or local maximum)
IF( OptimizationOption == 1 ) THEN
  ! Inflexion point
  bFunctionValue = dPressure_dDensity
ELSE IF( OptimizationOption == 2 ) THEN
  ! Local minimum
  bFunctionValue = Pressure
ELSE IF( OptimizationOption == 3 ) THEN
  ! Local maximum
  bFunctionValue = - Pressure
END IF

! Criterion
DO WHILE( DABS( (aMidpointInterval - bMidpointInterval) / aMidpointInterval ) >= Tolerance )

  IF( aFunctionValue < bFunctionValue ) THEN

    ! Swap intervals
    UpperBound = bMidpointInterval
    bMidpointInterval = aMidpointInterval
    bFunctionValue = aFunctionValue

    ! Probe point selection
    aMidpointInterval = UpperBound - (UpperBound - LowerBound) / GoldenRatio
    CALL Calculate_Pressure_Single_Component( cComponent, aMidpointInterval, Temperature, Pressure, IsothermalCompressibility, &
    &                                         ThermalExpansionCoefficient, CompressibilityFactor )
    dPressure_dDensity = aMidpointInterval / IsothermalCompressibility

    ! Optimization option (inflection point, local minimum, or local maximum)
    IF( OptimizationOption == 1 ) THEN
      ! Inflexion point
      aFunctionValue = dPressure_dDensity
    ELSE IF( OptimizationOption == 2 ) THEN
      ! Local minimum
      aFunctionValue = Pressure
    ELSE IF( OptimizationOption == 3 ) THEN
      ! Local maximum
      aFunctionValue = - Pressure
    END IF

  ELSE IF( aFunctionValue >= bFunctionValue ) THEN

    ! Swap intervals
    LowerBound = aMidpointInterval
    aMidpointInterval = bMidpointInterval
    aFunctionValue = bFunctionValue

    ! Probe point selection
    bMidpointInterval = LowerBound + (UpperBound - LowerBound) / GoldenRatio
    CALL Calculate_Pressure_Single_Component( cComponent, bMidpointInterval, Temperature, Pressure, IsothermalCompressibility, &
    &                                         ThermalExpansionCoefficient, CompressibilityFactor )
    dPressure_dDensity = bMidpointInterval / IsothermalCompressibility

    ! Optimization option (inflection point, local minimum, or local maximum)
    IF( OptimizationOption == 1 ) THEN
      ! Inflexion point
      bFunctionValue = dPressure_dDensity
    ELSE IF( OptimizationOption == 2 ) THEN
      ! Local minimum
      bFunctionValue = Pressure
    ELSE IF( OptimizationOption == 3 ) THEN
      ! Local maximum
      bFunctionValue = - Pressure
    END IF

  END IF

END DO

! Final result
CriticalPoint = 0.5D0 * (LowerBound + UpperBound)

RETURN

END SUBROUTINE Golden_Section_Single_Component