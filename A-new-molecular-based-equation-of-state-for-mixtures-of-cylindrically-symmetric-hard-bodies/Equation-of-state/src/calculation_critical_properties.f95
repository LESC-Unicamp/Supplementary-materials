! ************************************************************************************************ !
!                                       CRITICAL PROPERTIES                                        !
! ************************************************************************************************ !
!           This subroutine is used to estimate the critical properties of a substance.            !
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
SUBROUTINE Critical_Properties_Pure_Components( cComponent, cTemperature, cDensity, cPressure, cVolume )

! Uses two modules: global variables and substances
USE GlobalVariables
USE Substances

IMPLICIT NONE

! ************************************************************************************************ !
! INTEGER VARIABLES                                                                                !
! ************************************************************************************************ !
INTEGER( Kind= Int64 ) :: cComponent         ! Component index
INTEGER( Kind= Int64 ) :: OptimizationOption ! Optimization option of the golden-section search (1 = inflexion point; 2 = local minimum; 3 = local maximum)
INTEGER( Kind= Int64 ) :: Counter            ! Counter

! ************************************************************************************************ !
! REAL VARIABLES                                                                                   !
! ************************************************************************************************ !
REAL( Kind= Real64 ) :: LowerBoundVolume, UpperBoundVolume           ! Molar volume bounds
REAL( Kind= Real64 ) :: LowerBoundTemperature, UpperBoundTemperature ! Temperature bounds (subcritical and supercritical)
REAL( Kind= Real64 ) :: aInterval, bInterval                         ! Golden-section search interval (lower and upper bounds)
REAL( Kind= Real64 ) :: mVolumeInflection                            ! Molar volume at inflection point
REAL( Kind= Real64 ) :: cTemperature                                 ! Critical temperature
REAL( Kind= Real64 ) :: cPressure                                    ! Critical pressure
REAL( Kind= Real64 ) :: cVolume                                      ! Critical molar volume
REAL( Kind= Real64 ) :: cDensity                                     ! Critical molar density
REAL( Kind= Real64 ) :: Pressure                                     ! Pressure
REAL( Kind= Real64 ) :: IsothermalCompressibility                    ! Isothermal compressibility
REAL( Kind= Real64 ) :: ThermalExpansionCoefficient                  ! Thermal expansion coefficient
REAL( Kind= Real64 ) :: CompressibilityFactor                        ! Compressibility factor
REAL( Kind= Real64 ) :: TempVar                                      ! Temporary variable
REAL( Kind= Real64 ) :: PreviousPreviousTemperature                  ! Temperature (root guess)
REAL( Kind= Real64 ) :: PreviousTemperature                          ! Temperature (root guess)
REAL( Kind= Real64 ) :: PreviousFunction                             ! Derivative of pressure with respect to the molar density at the inflection point (guess)
REAL( Kind= Real64 ) :: MidpointTemperature                          ! Temperature (midpoint)
REAL( Kind= Real64 ) :: dPressure_dDensityInflection_SubT            ! First derivative of pressure with respect to the molar density at the inflection point (subcritical isotherm)
REAL( Kind= Real64 ) :: dPressure_dDensityInflection_SuperT          ! First derivative of pressure with respect to the molar density at the inflection point (super critical isotherm)
REAL( Kind= Real64 ) :: dPressure_dDensityInflection_CriticalT       ! First derivative of pressure with respect to the molar density at the inflection point (critical isotherm)

! ************************************************************************************************ !
! REAL VARIABLES (PARAMETER)                                                                       !
! ************************************************************************************************ !
REAL( Kind= Real64 ), PARAMETER :: Tolerance = 1.D-8 ! Numerical tolerance

! ************************************************************************************************ !
! LOGICAL VARIABLES                                                                                !
! ************************************************************************************************ !
LOGICAL :: BisectionLogical ! Checks whether the root from the bisection method will be used or not (Brent's method)

! Initialization
Counter = 0
OptimizationOption = 1

! Molar volume bounds [m³ / mol]
LowerBoundVolume = cPi * ijDiameterSphereCubic(cComponent,cComponent) * cAvogadro / 6.D0
UpperBoundVolume = 1.D10

! Temperature bounds (subcritical and supercritical)
LowerBoundTemperature = 1.D0 ! A very low temperature (subcritical isotherms)
UpperBoundTemperature = 1.D6 ! A very high temperature (supercritical isotherms)

! Molar volume at inflection point (subcritical temperature)
aInterval = LowerBoundVolume
bInterval = UpperBoundVolume
CALL Golden_Section_Single_Component( cComponent, LowerBoundTemperature, aInterval, bInterval, OptimizationOption, &
&                                     mVolumeInflection )

! Derivative of pressure with respect to the molar density at the inflection point (subcritical temperature)
CALL Calculate_Pressure_Single_Component( cComponent, mVolumeInflection, LowerBoundTemperature, Pressure, &
&                                         IsothermalCompressibility, ThermalExpansionCoefficient, CompressibilityFactor )
dPressure_dDensityInflection_subT = mVolumeInflection / IsothermalCompressibility

! Molar volume at inflection point (supercritical temperature)
aInterval = LowerBoundVolume
bInterval = UpperBoundVolume
CALL Golden_Section_Single_Component( cComponent, UpperBoundTemperature, aInterval, bInterval, OptimizationOption, &
&                                     mVolumeInflection )

! Derivative of pressure with respect to the molar density at the inflection point (supercritical temperature)
CALL Calculate_Pressure_Single_Component( cComponent, mVolumeInflection, UpperBoundTemperature, Pressure, &
&                                         IsothermalCompressibility, ThermalExpansionCoefficient, CompressibilityFactor )
dPressure_dDensityInflection_superT = mVolumeInflection / IsothermalCompressibility

! ************************************************************************************************ !
! Brent's method                                                                                   !
! ************************************************************************************************ !
!  The critical isotherm has a first derivative of 0 at the inflection point.                      !
! ************************************************************************************************ !
    
! Initialization
PreviousPreviousTemperature = 0.D0    ! Third-to-last root guess
BisectionLogical            = .FALSE. ! Bisection method flag
Counter                     = 0       ! Iteration counter

! Brent's condition
IF( ( dPressure_dDensityInflection_subT * dPressure_dDensityInflection_superT ) < 0.D0 ) THEN
  ! Swap condition
  IF( DABS( dPressure_dDensityInflection_subT ) < DABS( dPressure_dDensityInflection_superT ) ) THEN
    ! Swap temperature bounds
    TempVar = LowerBoundTemperature
    LowerBoundTemperature = UpperBoundTemperature
    UpperBoundTemperature = TempVar
    ! Swap function values
    TempVar = dPressure_dDensityInflection_subT
    dPressure_dDensityInflection_subT = dPressure_dDensityInflection_superT
    dPressure_dDensityInflection_superT = TempVar
  END IF
  ! Initialization of previous bound
  PreviousTemperature = LowerBoundTemperature
  PreviousFunction = dPressure_dDensityInflection_subT
  ! Initialize value of objective function (derivative of pressure with respect to the molar density at the inflection point)
  dPressure_dDensityInflection_criticalT = 1.D0
  ! Initialize 'BISECTION' flag as TRUE since the last iteration used was the bisection method
  BisectionLogical = .TRUE.
  ! Stop criterion
  DO WHILE( ( DABS( dPressure_dDensityInflection_superT ) >= Tolerance .OR. DABS( dPressure_dDensityInflection_subT ) >= &
    &       Tolerance ) .AND. DABS( dPressure_dDensityInflection_criticalT ) >= Tolerance .AND. DABS( ( LowerBoundTemperature - &
    &       UpperBoundTemperature ) / LowerBoundTemperature ) >= Tolerance )
    ! Initialize root of the function (critical temperature)
    IF( DABS( dPressure_dDensityInflection_subT - PreviousFunction ) >= EPSILON( 1.D0 ) .AND. &
    &   DABS( dPressure_dDensityInflection_superT - PreviousFunction ) >= EPSILON( 1.D0 ) ) THEN
      ! Inverse quadratic interpolation root-finding procedure
      MidpointTemperature = ( LowerBoundTemperature * dPressure_dDensityInflection_superT * PreviousFunction ) / &
      &                     ( ( dPressure_dDensityInflection_subT - dPressure_dDensityInflection_superT ) * &
      &                     ( dPressure_dDensityInflection_subT - PreviousFunction ) ) + ( UpperBoundTemperature * &
      &                     dPressure_dDensityInflection_subT * PreviousFunction ) / ( ( dPressure_dDensityInflection_superT - &
      &                     dPressure_dDensityInflection_subT ) * ( dPressure_dDensityInflection_superT - PreviousFunction ) ) + &
      &                     ( PreviousTemperature * dPressure_dDensityInflection_subT * dPressure_dDensityInflection_superT ) / &
      &                     ( ( PreviousFunction - dPressure_dDensityInflection_subT ) * ( PreviousFunction - &
      &                     dPressure_dDensityInflection_superT ) )
    ELSE
      ! False position formula
      MidpointTemperature = UpperBoundTemperature - ( dPressure_dDensityInflection_superT * ( UpperBoundTemperature - &
      &                     LowerBoundTemperature ) ) / ( dPressure_dDensityInflection_superT - dPressure_dDensityInflection_subT )
    END IF
    ! Check whether the root obtained from the interpolation method or false position formula will be used; otherwise, use the midpoint from bisection method
    IF( ( ( MidpointTemperature - ( 3.D0 * LowerBoundTemperature + UpperBoundTemperature ) / 4.D0 ) * ( MidpointTemperature - &
    &   UpperBoundTemperature ) >= 0.D0 ) .OR. ( BisectionLogical .AND. ( DABS( MidpointTemperature - UpperBoundTemperature ) >= &
    &   DABS( ( UpperBoundTemperature - PreviousTemperature ) / 2.D0 ) ) ) .OR. ( .NOT. BisectionLogical .AND. &
    &   ( DABS( MidpointTemperature - UpperBoundTemperature ) >= DABS( ( PreviousTemperature - PreviousPreviousTemperature ) / &
    &   2.D0 ) ) ) .OR. ( BisectionLogical .AND. ( DABS( MidpointTemperature - PreviousTemperature ) < Tolerance ) ) .OR. ( .NOT. &
    &   BisectionLogical .AND. ( DABS( PreviousTemperature - PreviousPreviousTemperature ) < Tolerance ) ) ) THEN
      ! Root from bisection method
      BisectionLogical    = .TRUE.
      MidpointTemperature = 0.5D0 * ( LowerBoundTemperature + UpperBoundTemperature )
    ELSE
      ! Root from inverse quadratic interpolation or false position formula
      BisectionLogical = .FALSE.
    END IF
    ! Molar volume at inflection point (midpoint/critical temperature)
    aInterval = LowerBoundVolume
    bInterval = UpperBoundVolume
    CALL Golden_Section_Single_Component( cComponent, MidpointTemperature, aInterval, bInterval, OptimizationOption, &
    &                                     mVolumeInflection )
    ! Derivative of pressure with respect to the molar density at the inflection point (midpoint temperature)
    CALL Calculate_Pressure_Single_Component( cComponent, mVolumeInflection, MidpointTemperature, Pressure, &
    &                                         IsothermalCompressibility, ThermalExpansionCoefficient, CompressibilityFactor )
    dPressure_dDensityInflection_criticalT = mVolumeInflection / IsothermalCompressibility
    ! Set third-to-last root guess
    PreviousPreviousTemperature = PreviousTemperature
    ! Set second-to-last root guess
    PreviousTemperature = UpperBoundTemperature
    PreviousFunction = dPressure_dDensityInflection_superT
    ! Check signs of functions
    IF( dPressure_dDensityInflection_subT * dPressure_dDensityInflection_criticalT < 0.D0 ) THEN
      UpperBoundTemperature = MidpointTemperature
      dPressure_dDensityInflection_superT = dPressure_dDensityInflection_criticalT
    ELSE
      LowerBoundTemperature = MidpointTemperature
      dPressure_dDensityInflection_subT = dPressure_dDensityInflection_criticalT
    END IF
    ! Swap condition
    IF( DABS( dPressure_dDensityInflection_subT ) < DABS( dPressure_dDensityInflection_superT ) ) THEN
      ! Swap bounds
      TempVar = LowerBoundTemperature
      LowerBoundTemperature = UpperBoundTemperature
      UpperBoundTemperature = TempVar
      ! Swap function values
      TempVar = dPressure_dDensityInflection_subT
      dPressure_dDensityInflection_subT = dPressure_dDensityInflection_superT
      dPressure_dDensityInflection_superT = TempVar
    END IF
    ! Iteration
    Counter = Counter + 1
    IF( Counter > 50 ) THEN
      WRITE( *, "(4G0)" ) "Brent's method failed to converge after ", Counter, " iterations during the routine procedure that ", &
      &                   "searches for the critical properties of pure components. Exiting..."
      CALL Exit(  )
    END IF
  END DO
  ! Root (iterative process or initial guess)
  IF( DABS( dPressure_dDensityInflection_criticalT ) < Tolerance ) THEN
    cTemperature = MidpointTemperature
  ELSE IF( DABS( dPressure_dDensityInflection_superT ) < Tolerance .AND. DABS( ( LowerBoundTemperature - UpperBoundTemperature ) / &
    &      LowerBoundTemperature ) < Tolerance ) THEN
    cTemperature = UpperBoundTemperature
  ELSE IF( DABS( dPressure_dDensityInflection_subT ) < Tolerance .AND. DABS( ( LowerBoundTemperature - UpperBoundTemperature ) / &
    &      LowerBoundTemperature ) < Tolerance ) THEN
    cTemperature = LowerBoundTemperature
  ELSE
    cTemperature = 0.5D0 * ( LowerBoundTemperature + UpperBoundTemperature )
  END IF
! Algorithm should not reach this point
ELSE
  ! Brent's condition not satisfied
  WRITE( *, "(2G0)" ) "Brent's condition not satisfied during the routine procedure that searches for the critical properties ", &
  &                   "of pure components. Exiting..."
  CALL Exit(  )
END IF

! ************************************************************************************************ !
! Critical Properties                                                                              !
! ************************************************************************************************ !

! Critical temperature [K]
cTemperature = MidpointTemperature

! Critical molar density [mol / m³]
cDensity = 1.D0 / mVolumeInflection

! Critical molar volume [m³ / mol]
cVolume = mVolumeInflection

! Critical pressure [Pa]
CALL Calculate_Pressure_Single_Component( cComponent, cVolume, cTemperature, cPressure, IsothermalCompressibility, &
&                                         ThermalExpansionCoefficient, CompressibilityFactor )

RETURN

END SUBROUTINE Critical_Properties_Pure_Components