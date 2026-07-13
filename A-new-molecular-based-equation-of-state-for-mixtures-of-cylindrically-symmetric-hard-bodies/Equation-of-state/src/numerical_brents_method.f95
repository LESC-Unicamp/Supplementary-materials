! ************************************************************************************************ !
!                                          BRENT'S METHOD                                          !
! ************************************************************************************************ !
!           This subroutine applies the Brent's method to find the roots of a function.            !
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
!                                  BRENT'S METHOD (DENSITY ROOT)                                   !
! ************************************************************************************************ !
! This subroutine is used to determine the density roots through the Brent's method.               !
! ************************************************************************************************ !
SUBROUTINE Brent_Method_Density_Root( cComponent, mFraction, Temperature, SpecifiedPressure, LowerBound, UpperBound, mVolumeRoot, &
&                                     PureComponent )

! Uses two modules: global variables and substances
USE GlobalVariables
USE Substances

IMPLICIT NONE

! ************************************************************************************************ !
! INTEGER VARIABLES                                                                                !
! ************************************************************************************************ !
INTEGER( Kind= Int64 ) :: cComponent ! Component index
INTEGER( Kind= Int64 ) :: Counter    ! Iteration counter

! ************************************************************************************************ !
! REAL VARIABLES (PARAMETERS)                                                                      !
! ************************************************************************************************ !
REAL( Kind= Real64 ) :: Tolerance = 1.D-10 ! Numerical tolerance

! ************************************************************************************************ !
! REAL VARIABLES                                                                                   !
! ************************************************************************************************ !
REAL( Kind= Real64 )                           :: LowerBound, UpperBound      ! Molar volume bounds
REAL( Kind= Real64 )                           :: SpecifiedPressure           ! Pressure (specified)
REAL( Kind= Real64 )                           :: Pressure                    ! Pressure (calculated)
REAL( Kind= Real64 )                           :: Temperature                 ! Temperature
REAL( Kind= Real64 )                           :: mVolumeRoot                 ! Molar volume (root)
REAL( Kind= Real64 )                           :: IsothermalCompressibility   ! Isothermal compressibility
REAL( Kind= Real64 )                           :: ThermalExpansionCoefficient ! Thermal expansion coefficient
REAL( Kind= Real64 )                           :: CompressibilityFactor       ! Compressibility factor
REAL( Kind= Real64 )                           :: TempVar                     ! Temporary variable
REAL( Kind= Real64 )                           :: LowerPressureDifference     ! Lower pressure difference (objective function)
REAL( Kind= Real64 )                           :: UpperPressureDifference     ! Upper pressure difference (objective function)
REAL( Kind= Real64 )                           :: PreviousPreviousVolume      ! Volume (root guess)
REAL( Kind= Real64 )                           :: PreviousVolume              ! Volume (root guess)
REAL( Kind= Real64 )                           :: PreviousPDifference         ! Pressure difference (guess)
REAL( Kind= Real64 )                           :: MidpointPressureDifference  ! Pressure difference (objective function)
REAL( Kind= Real64 )                           :: MidpointVolume              ! Molar volume (midpoint)
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: mFraction                   ! Molar fraction of a component

! ************************************************************************************************ !
! LOGICAL VARIABLES                                                                                !
! ************************************************************************************************ !
LOGICAL :: PureComponent    ! Checks whether the system is a pure component or a mixture
LOGICAL :: BisectionLogical ! Checks whether or not the maximum iterations will be increased

! Initialization
Counter = 0

! Initialize pressure difference (lower bound)
IF( PureComponent ) THEN
  CALL Calculate_Pressure_Single_Component( cComponent, LowerBound, Temperature, Pressure, IsothermalCompressibility, &
  &                                         ThermalExpansionCoefficient, CompressibilityFactor )
ELSE
  CALL Calculate_Pressure( mFraction, LowerBound, Temperature, Pressure, IsothermalCompressibility, ThermalExpansionCoefficient, &
  &                        CompressibilityFactor )
END IF
LowerPressureDifference = Pressure - SpecifiedPressure

! Initialize pressure difference (Upper bound)
IF( PureComponent ) THEN
  CALL Calculate_Pressure_Single_Component( cComponent, UpperBound, Temperature, Pressure, IsothermalCompressibility, &
  &                                         ThermalExpansionCoefficient, CompressibilityFactor )
ELSE
  CALL Calculate_Pressure( mFraction, UpperBound, Temperature, Pressure, IsothermalCompressibility, ThermalExpansionCoefficient, &
  &                        CompressibilityFactor )
END IF
UpperPressureDifference = Pressure - SpecifiedPressure

! ************************************************************************************************ !
! Brent's method                                                                                   !
! ************************************************************************************************ !
    
! Initialization
PreviousPreviousVolume = 0.D0    ! Third-to-last root guess
BisectionLogical       = .FALSE. ! Bisection method flag
Counter                = 0       ! Iteration counter

! Brent's condition
IF( ( LowerPressureDifference * UpperPressureDifference ) < 0.D0 ) THEN
  ! Swap condition
  IF( DABS( LowerPressureDifference ) < DABS( UpperPressureDifference ) ) THEN
    ! Swap volume bounds
    TempVar = LowerBound
    LowerBound = UpperBound
    UpperBound = TempVar
    ! Swap function values
    TempVar = LowerPressureDifference
    LowerPressureDifference = UpperPressureDifference
    UpperPressureDifference = TempVar
  END IF
  ! Initialization of previous bound
  PreviousVolume = LowerBound
  PreviousPDifference = LowerPressureDifference
  ! Initialize value of objective function (derivative of pressure with respect to the molar density at the inflection point)
  MidpointPressureDifference = 1.D0
  ! Initialize 'BISECTION' flag as TRUE since the last iteration used was the bisection method
  BisectionLogical = .TRUE.
  ! Stop criterion
  DO WHILE( ( DABS( UpperPressureDifference ) >= Tolerance .OR. DABS( LowerPressureDifference ) >= Tolerance ) .AND. &
    &       DABS( MidpointPressureDifference ) >= Tolerance .AND. DABS( ( LowerBound - UpperBound ) / LowerBound ) >= Tolerance )
    ! Initialize root of the function (critical temperature)
    IF( DABS( LowerPressureDifference - PreviousPDifference ) >= EPSILON( 1.D0 ) .AND. DABS( UpperPressureDifference - &
    &   PreviousPDifference ) >= EPSILON( 1.D0 ) ) THEN
      ! Inverse quadratic interpolation root-finding procedure
      MidpointVolume = ( LowerBound * UpperPressureDifference * PreviousPDifference ) / ( ( LowerPressureDifference - &
      &                UpperPressureDifference ) * ( LowerPressureDifference - PreviousPDifference ) ) + ( UpperBound * &
      &                LowerPressureDifference * PreviousPDifference ) / ( ( UpperPressureDifference - LowerPressureDifference ) * &
      &                ( UpperPressureDifference - PreviousPDifference ) ) + ( PreviousVolume * LowerPressureDifference * &
      &                UpperPressureDifference ) / ( ( PreviousPDifference - LowerPressureDifference ) * ( PreviousPDifference - &
      &                UpperPressureDifference ) )
    ELSE
      ! False position formula
      MidpointVolume = UpperBound - ( UpperPressureDifference * ( UpperBound - LowerBound ) ) / ( UpperPressureDifference - &
      &                LowerPressureDifference )
    END IF
    ! Check whether the root obtained from the interpolation method or false position formula will be used; otherwise, use the midpoint from bisection method
    IF( ( ( MidpointVolume - ( 3.D0 * LowerBound + UpperBound ) / 4.D0 ) * ( MidpointVolume - UpperBound ) >= 0.D0 ) .OR. &
    &   ( BisectionLogical .AND. ( DABS( MidpointVolume - UpperBound ) >= DABS( ( UpperBound - PreviousVolume ) / 2.D0 ) ) ) .OR. &
    &   ( .NOT. BisectionLogical .AND. ( DABS( MidpointVolume - UpperBound ) >= DABS( ( PreviousVolume - PreviousPreviousVolume ) &
    &   / 2.D0 ) ) ) .OR. ( BisectionLogical .AND. ( DABS( MidpointVolume - PreviousVolume ) < Tolerance ) ) .OR. &
    &   ( .NOT. BisectionLogical .AND. ( DABS( PreviousVolume - PreviousPreviousVolume ) < Tolerance ) ) ) THEN
      ! Root from bisection method
      BisectionLogical = .TRUE.
      MidpointVolume   = 0.5D0 * ( LowerBound + UpperBound )
    ELSE
      ! Root from inverse quadratic interpolation or false position formula
      BisectionLogical = .FALSE.
    END IF
    ! Pressure difference at midpoint volume
    IF( PureComponent ) THEN
      CALL Calculate_Pressure_Single_Component( cComponent, MidpointVolume, Temperature, Pressure, IsothermalCompressibility, &
      &                                         ThermalExpansionCoefficient, CompressibilityFactor )
    ELSE
      CALL Calculate_Pressure( mFraction, MidpointVolume, Temperature, Pressure, IsothermalCompressibility, &
      &                        ThermalExpansionCoefficient, CompressibilityFactor )
    END IF
    MidpointPressureDifference = Pressure - SpecifiedPressure
    ! Set third-to-last root guess
    PreviousPreviousVolume = PreviousVolume
    ! Set second-to-last root guess
    PreviousVolume = UpperBound
    PreviousPDifference = UpperPressureDifference
    ! Check signs of functions
    IF( LowerPressureDifference * MidpointPressureDifference < 0.D0 ) THEN
      UpperBound = MidpointVolume
      UpperPressureDifference = MidpointPressureDifference
    ELSE
      LowerBound = MidpointVolume
      LowerPressureDifference = MidpointPressureDifference
    END IF
    ! Swap condition
    IF( DABS( LowerPressureDifference ) < DABS( UpperPressureDifference ) ) THEN
      ! Swap bounds
      TempVar = LowerBound
      LowerBound = UpperBound
      UpperBound = TempVar
      ! Swap function values
      TempVar = LowerPressureDifference
      LowerPressureDifference = UpperPressureDifference
      UpperPressureDifference = TempVar
    END IF
    ! Iteration
    Counter = Counter + 1
    IF( Counter > 200 ) THEN
      WRITE( *, "(4G0)" ) "Brent's method failed to converge after ", Counter, " iterations during the routine procedure that ", &
      &                   "calculates the volume roots. Exiting..."
      CALL Exit(  )
    END IF
  END DO
  ! Root (iterative process or initial guess)
  IF( DABS( MidpointPressureDifference ) < Tolerance ) THEN
    mVolumeRoot = MidpointVolume
  ELSE IF( DABS( UpperPressureDifference ) < Tolerance .AND. DABS( ( LowerBound - UpperBound ) / &
  &        LowerBound ) < Tolerance ) THEN
    mVolumeRoot = UpperBound
  ELSE IF( DABS( LowerPressureDifference ) < Tolerance .AND. DABS( ( LowerBound - UpperBound ) / &
  &        LowerBound ) < Tolerance ) THEN
    mVolumeRoot = LowerBound
  ELSE
    mVolumeRoot = 0.5D0 * ( LowerBound + UpperBound )
  END IF
! Algorithm should not reach this point
ELSE
  ! Brent's condition not satisfied
  WRITE( *, "(G0)" ) "Brent's condition not satisfied during the routine procedure that calculates the volume roots. Exiting..."
  CALL Exit(  )
END IF

RETURN

END SUBROUTINE Brent_Method_Density_Root

! ************************************************************************************************ !
!                                  BRENT'S METHOD (DENSITY ROOT)                                   !
! ************************************************************************************************ !
! This subroutine is used to determine the density roots through the Brent's method.               !
! ************************************************************************************************ !
SUBROUTINE Brent_Method_Rachford_Rice( LowerBound, UpperBound, mFraction, EquilibriumFactors, VaporPhaseFraction )

! Uses two modules: global variables and substances
USE GlobalVariables
USE Substances

IMPLICIT NONE

! ************************************************************************************************ !
! INTEGER VARIABLES                                                                                !
! ************************************************************************************************ !
INTEGER( Kind= Int64 ) :: Counter ! Iteration counter

! ************************************************************************************************ !
! REAL VARIABLES (PARAMETERS)                                                                      !
! ************************************************************************************************ !
REAL( Kind= Real64 ) :: Tolerance = 1.D-10 ! Numerical tolerance

! ************************************************************************************************ !
! REAL VARIABLES                                                                                   !
! ************************************************************************************************ !
REAL( Kind= Real64 )                           :: LowerBound, UpperBound ! Vapor phase fraction bounds
REAL( Kind= Real64 )                           :: VaporPhaseFraction     ! Molar fraction of the vapor phase
REAL( Kind= Real64 )                           :: TempVar                ! Temporary variable
REAL( Kind= Real64 )                           :: LowerBoundRachfordRice ! Lower Rachford-Rice function (objective function)
REAL( Kind= Real64 )                           :: UpperBoundRachfordRice ! Upper Rachford-Rice function (objective function)
REAL( Kind= Real64 )                           :: PreviousPreviousBeta   ! Vapor phase fraction (root guess)
REAL( Kind= Real64 )                           :: PreviousBeta           ! Vapor phase fraction (root guess)
REAL( Kind= Real64 )                           :: PreviousRachfordRice   ! Rachford-Rice function (guess)
REAL( Kind= Real64 )                           :: MidpointRachfordRice   ! Rachford-Rice function (objective function)
REAL( Kind= Real64 )                           :: MidpointBeta           ! Vapor phase fraction (midpoint)
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: EquilibriumFactors     ! Equilibrium factors
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: mFraction              ! Molar fraction of a component

! ************************************************************************************************ !
! LOGICAL VARIABLES                                                                                !
! ************************************************************************************************ !
LOGICAL :: BisectionLogical ! Checks whether or not the maximum iterations will be increased

! Initialization
Counter = 0

! Rachford-Rice function at the lower bound
LowerBoundRachfordRice = SUM( mFraction * ( EquilibriumFactors - 1.D0 ) / ( 1.D0 - LowerBound + LowerBound * EquilibriumFactors ) )

! Rachford-Rice function at the upper bound
UpperBoundRachfordRice = SUM( mFraction * ( EquilibriumFactors - 1.D0 ) / ( 1.D0 - UpperBound + UpperBound * EquilibriumFactors ) )

! ************************************************************************************************ !
! Brent's method                                                                                   !
! ************************************************************************************************ !
    
! Initialization
PreviousPreviousBeta = 0.D0    ! Third-to-last root guess
BisectionLogical     = .FALSE. ! Bisection method flag
Counter              = 0       ! Iteration counter

! Brent's condition
IF( ( LowerBoundRachfordRice * UpperBoundRachfordRice ) < 0.D0 ) THEN
  ! Swap condition
  IF( DABS( LowerBoundRachfordRice ) < DABS( UpperBoundRachfordRice ) ) THEN
    ! Swap volume bounds
    TempVar = LowerBound
    LowerBound = UpperBound
    UpperBound = TempVar
    ! Swap function values
    TempVar = LowerBoundRachfordRice
    LowerBoundRachfordRice = UpperBoundRachfordRice
    UpperBoundRachfordRice = TempVar
  END IF
  ! Initialization of previous bound
  PreviousBeta = LowerBound
  PreviousRachfordRice = LowerBoundRachfordRice
  ! Initialize value of objective function (Rachford-Rice function)
  MidpointRachfordRice = 1.D0
  ! Initialize 'BISECTION' flag as TRUE since the last iteration used was the bisection method
  BisectionLogical = .TRUE.
  ! Stop criterion
  DO WHILE( ( DABS( UpperBoundRachfordRice ) >= Tolerance .OR. DABS( LowerBoundRachfordRice ) >= Tolerance ) .AND. &
    &       DABS( MidpointRachfordRice ) >= Tolerance .AND. DABS( ( LowerBound - UpperBound ) / LowerBound ) >= Tolerance )
    ! Initialize root of the function (vapor phase fraction)
    IF( DABS( LowerBoundRachfordRice - PreviousRachfordRice ) >= EPSILON( 1.D0 ) .AND. DABS( UpperBoundRachfordRice - &
    &   PreviousRachfordRice ) >= EPSILON( 1.D0 ) ) THEN
      ! Inverse quadratic interpolation root-finding procedure
      MidpointBeta = ( LowerBound * UpperBoundRachfordRice * PreviousRachfordRice ) / ( ( LowerBoundRachfordRice - &
      &              UpperBoundRachfordRice ) * ( LowerBoundRachfordRice - PreviousRachfordRice ) ) + ( UpperBound * &
      &              LowerBoundRachfordRice * PreviousRachfordRice ) / ( ( UpperBoundRachfordRice - LowerBoundRachfordRice ) * &
      &              ( UpperBoundRachfordRice - PreviousRachfordRice ) ) + ( PreviousBeta * LowerBoundRachfordRice * &
      &              UpperBoundRachfordRice ) / ( ( PreviousRachfordRice - LowerBoundRachfordRice ) * ( PreviousRachfordRice - &
      &              UpperBoundRachfordRice ) )
    ELSE
      ! False position formula
      MidpointBeta = UpperBound - ( UpperBoundRachfordRice * ( UpperBound - LowerBound ) ) / ( UpperBoundRachfordRice - &
      &              LowerBoundRachfordRice )
    END IF
    ! Check whether the root obtained from the interpolation method or false position formula will be used; otherwise, use the midpoint from bisection method
    IF( ( ( MidpointBeta - ( 3.D0 * LowerBound + UpperBound ) / 4.D0 ) * ( MidpointBeta - UpperBound ) >= 0.D0 ) .OR. &
    &   ( BisectionLogical .AND. ( DABS( MidpointBeta - UpperBound ) >= DABS( ( UpperBound - PreviousBeta ) / 2.D0 ) ) ) .OR. &
    &   ( .NOT. BisectionLogical .AND. ( DABS( MidpointBeta - UpperBound ) >= DABS( ( PreviousBeta - PreviousPreviousBeta ) &
    &   / 2.D0 ) ) ) .OR. ( BisectionLogical .AND. ( DABS( MidpointBeta - PreviousBeta ) < Tolerance ) ) .OR. &
    &   ( .NOT. BisectionLogical .AND. ( DABS( PreviousBeta - PreviousPreviousBeta ) < Tolerance ) ) ) THEN
      ! Root from bisection method
      BisectionLogical = .TRUE.
      MidpointBeta = 0.5D0 * ( LowerBound + UpperBound )
    ELSE
      ! Root from inverse quadratic interpolation or false position formula
      BisectionLogical = .FALSE.
    END IF
    ! Rachford-Rice function at midpoint vapor phase fraction
    MidpointRachfordRice = SUM( mFraction * ( EquilibriumFactors - 1.D0 ) / (1.D0 - MidpointBeta + MidpointBeta * &
    &                      EquilibriumFactors ) )
    ! Set third-to-last root guess
    PreviousPreviousBeta = PreviousBeta
    ! Set second-to-last root guess
    PreviousBeta = UpperBound
    PreviousRachfordRice = UpperBoundRachfordRice
    ! Check signs of functions
    IF( LowerBoundRachfordRice * MidpointRachfordRice < 0.D0 ) THEN
      UpperBound = MidpointBeta
      UpperBoundRachfordRice = MidpointRachfordRice
    ELSE
      LowerBound = MidpointBeta
      LowerBoundRachfordRice = MidpointRachfordRice
    END IF
    ! Swap condition
    IF( DABS( LowerBoundRachfordRice ) < DABS( UpperBoundRachfordRice ) ) THEN
      ! Swap bounds
      TempVar = LowerBound
      LowerBound = UpperBound
      UpperBound = TempVar
      ! Swap function values
      TempVar = LowerBoundRachfordRice
      LowerBoundRachfordRice = UpperBoundRachfordRice
      UpperBoundRachfordRice = TempVar
    END IF
    ! Iteration
    Counter = Counter + 1
    IF( Counter > 200 ) THEN
      WRITE( *, "(4G0)" ) "Brent's method failed to converge after ", Counter, " iterations during the routine procedure that ", &
      &                   "calculates the vapor phase fractions by minimizing the Rachford-Rice function. Exiting..."
      CALL Exit(  )
    END IF
  END DO
  ! Root (iterative process or initial guess)
  IF( DABS( MidpointRachfordRice ) < Tolerance ) THEN
    VaporPhaseFraction = MidpointBeta
  ELSE IF( DABS( UpperBoundRachfordRice ) < Tolerance .AND. DABS( ( LowerBound - UpperBound ) / LowerBound ) < Tolerance ) THEN
    VaporPhaseFraction = UpperBound
  ELSE IF( DABS( LowerBoundRachfordRice ) < Tolerance .AND. DABS( ( LowerBound - UpperBound ) / LowerBound ) < Tolerance ) THEN
    VaporPhaseFraction = LowerBound
  ELSE
    VaporPhaseFraction = 0.5D0 * ( LowerBound + UpperBound )
  END IF
! Algorithm should not reach this point
ELSE
  ! Brent's condition not satisfied
  WRITE( *, "(G0)" ) " "
  WRITE( *, "(G0)" ) " "
  WRITE( *, "(2G0)" ) "Brent's condition not satisfied during the routine procedure that calculates the vapor phase fractions ", &
  &                   "by minimizing the Rachford-Rice function. Exiting..."
  CALL Exit(  )
END IF

RETURN

END SUBROUTINE Brent_Method_Rachford_Rice