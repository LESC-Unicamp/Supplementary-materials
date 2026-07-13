! ************************************************************************************************ !
!                                     SUCCESSIVE SUBSTITUTION                                      !
! ************************************************************************************************ !
!           This subroutine used the method of successive substitution to calculate the            !
!            vapor-liquid equilibrium of binary mixtures (isothermal two-phase flash).             !
! ************************************************************************************************ !
! => ORIGINAL DEVELOPER: Luís Fernando Mercier Franco                                              !
! => AUTHOR:             Nathan Barros de Souza                                                    !
! => E-MAIL:             n264179@dac.unicamp.br                                                    !
! => SUPERVISOR:         Luís Fernando Mercier Franco                                              !
! ************************************************************************************************ !
! Main References:                M. L. Michelsen, J. M. Mollerup                                  !
!                   "Themodynamic Models", 2nd Ed., Tie-Line Publications (2007)                   !
!                                                                                                  !
! ************************************************************************************************ !
!                       University of Campinas, Campinas - São Paulo, Brazil                       !
!                                         March 7th, 2024                                          !
! ************************************************************************************************ !
! Version 1.0.0                                                                                    !
! ************************************************************************************************ !

! ************************************************************************************************ !
!               This subroutine tests if the feed composition is a subcooled liquid                !
! ************************************************************************************************ !
SUBROUTINE Subcooled_Liquid_Test( Temperature, Pressure, FeedComposition, EquilibriumFactor, LiquidComposition, VaporComposition, &
&                                 Subcooled_Liquid )

! Uses two modules: global variables and substances
USE GlobalVariables
USE Substances

IMPLICIT NONE

! ************************************************************************************************ !
! INTEGER VARIABLES                                                                                !
! ************************************************************************************************ !
INTEGER( Kind= Int64 ) :: Counter ! Counter

! ************************************************************************************************ !
! REAL PARAMETERS                                                                                  !
! ************************************************************************************************ !
REAL( Kind= Real64 ), PARAMETER :: Tolerance = 1.D-6 ! Numerical tolerance

! ************************************************************************************************ !
! REAL VARIABLES                                                                                   !
! ************************************************************************************************ !
REAL( Kind= Real64 )                           :: LowerBoundRachfordRice ! Lower bound (Rachford-Rice equation)
REAL( Kind= Real64 )                           :: Temperature            ! Temperature
REAL( Kind= Real64 )                           :: Pressure               ! Pressure
REAL( Kind= Real64 )                           :: ErrorDifference        ! Difference between old and new equilibrium factors
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: FeedComposition        ! Feed composition
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: LiquidComposition      ! Liquid composition
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: VaporComposition       ! Vapor composition
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: EquilibriumFactor      ! Equilibrium factor

! ************************************************************************************************ !
! LOGICAL VARIABLES                                                                                !
! ************************************************************************************************ !
LOGICAL :: Subcooled_Liquid ! Subcooled liquid flag

! Initialization
Subcooled_Liquid = .FALSE.
Counter = 0

! Rachford-Rice equation at β = 0
LowerBoundRachfordRice = SUM( FeedComposition * EquilibriumFactor ) - 1.D0

! Violation of a condition for the Rachford-Rice equation be monotonically decreasing
IF( LowerBoundRachfordRice <= 0.D0 ) THEN
  LiquidComposition = FeedComposition
  ErrorDifference = 1.D0
  DO WHILE( DABS( ErrorDifference ) >= Tolerance .AND. Counter <= 20 )
    VaporComposition = LiquidComposition * EquilibriumFactor / SUM( FeedComposition * EquilibriumFactor )
    CALL Equilibrium_Factor_Calculation( Temperature, Pressure, VaporComposition, LiquidComposition, EquilibriumFactor, &
    &                                    ErrorDifference )
    Counter = Counter + 1
  END DO
  LowerBoundRachfordRice = SUM( FeedComposition * EquilibriumFactor ) - 1.D0
  IF( LowerBoundRachfordRice <= 0.D0 ) THEN
    Subcooled_Liquid = .TRUE.
  END IF
END IF

RETURN

END SUBROUTINE Subcooled_Liquid_Test

! ************************************************************************************************ !
!               This subroutine tests if the feed composition is a superheated vapor               !
! ************************************************************************************************ !
SUBROUTINE Superheated_Vapor_Test( Temperature, Pressure, FeedComposition, EquilibriumFactor, LiquidComposition, VaporComposition, &
&                                  Superheated_Vapor )

! Uses two modules: global variables and substances
USE GlobalVariables
USE Substances
  
IMPLICIT NONE

! ************************************************************************************************ !
! INTEGER VARIABLES                                                                                !
! ************************************************************************************************ !
INTEGER( Kind= Int64 ) :: Counter ! Counter

! ************************************************************************************************ !
! REAL PARAMETERS                                                                                  !
! ************************************************************************************************ !
REAL( Kind= Real64 ), PARAMETER :: Tolerance = 1.D-6 ! Numerical tolerance

! ************************************************************************************************ !
! REAL VARIABLES                                                                                   !
! ************************************************************************************************ !
REAL( Kind= Real64 )                           :: UpperBoundRachfordRice ! Lower bound (Rachford-Rice equation)
REAL( Kind= Real64 )                           :: Temperature            ! Temperature
REAL( Kind= Real64 )                           :: Pressure               ! Pressure
REAL( Kind= Real64 )                           :: ErrorDifference        ! Difference between old and new equilibrium factors
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: FeedComposition        ! Feed composition
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: LiquidComposition      ! Liquid composition
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: VaporComposition       ! Vapor composition
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: EquilibriumFactor      ! Equilibrium factor

! ************************************************************************************************ !
! LOGICAL VARIABLES                                                                                !
! ************************************************************************************************ !
LOGICAL :: Superheated_Vapor ! Superheated vapor flag

! Initialization
Superheated_Vapor = .FALSE.
Counter = 0

! Rachford-Rice equation at β = 1
UpperBoundRachfordRice = 1.D0 - SUM( FeedComposition / EquilibriumFactor )

! Violation of a condition for the Rachford-Rice equation be monotonically decreasing
IF( UpperBoundRachfordRice >= 0.D0 ) THEN
  VaporComposition = FeedComposition
  ErrorDifference = 1.D0
  DO WHILE( DABS( ErrorDifference ) >= Tolerance .AND. Counter <= 20 )
    LiquidComposition = VaporComposition / EquilibriumFactor / SUM( FeedComposition / EquilibriumFactor )
    CALL Equilibrium_Factor_Calculation( Temperature, Pressure, VaporComposition, LiquidComposition, EquilibriumFactor, &
    &                                    ErrorDifference )
    Counter = Counter + 1
  END DO
  UpperBoundRachfordRice = 1.D0 - SUM( FeedComposition / EquilibriumFactor )
  IF( UpperBoundRachfordRice >= 0.D0 ) THEN
    Superheated_Vapor = .TRUE.
  END IF
END IF

RETURN

END SUBROUTINE Superheated_Vapor_Test

! ************************************************************************************************ !
!                  This subroutine calculates the equilibrium factor of a mixture                  !
! ************************************************************************************************ !
SUBROUTINE Equilibrium_Factor_Calculation( Temperature, Pressure, VaporComposition, LiquidComposition, OldEquilibriumFactor, Error )

! Uses two modules: global variables and substances
USE GlobalVariables
USE Substances
  
IMPLICIT NONE

! ************************************************************************************************ !
! REAL VARIABLES                                                                                   !
! ************************************************************************************************ !
REAL( Kind= Real64 )                           :: IsothermalCompressibility   ! Isothermal compressibility
REAL( Kind= Real64 )                           :: ThermalExpansionCoefficient ! Thermal expansion coefficient
REAL( Kind= Real64 )                           :: Temperature                 ! Temperature
REAL( Kind= Real64 )                           :: Pressure                    ! Pressure
REAL( Kind= Real64 )                           :: mVolumeLiquidPhase          ! Molar volume of the liquid phase
REAL( Kind= Real64 )                           :: mVolumeVaporPhase           ! Molar volume of the vapor phase
REAL( Kind= Real64 )                           :: CompressibilityFactor       ! Compressibility factor
REAL( Kind= Real64 )                           :: Error                       ! Difference between old and new equilibrium factors
REAL( Kind= Real64 )                           :: Anynumber                   ! Any number (dummy)
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: VaporComposition            ! Vapor composition
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: LiquidComposition           ! Liquid composition
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: OldEquilibriumFactor        ! Old equilibrium factors
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: NewEquilibriumFactor        ! New equilibrium factors
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: TotalChemicalPotential      ! Total chemical potential
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: ResidualChemicalPotential   ! Residual chemical potential
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: cFugacityVaporPhase         ! Fugacity of species i in the vapor phase
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: cFugacityLiquidPhase        ! Fugacity of species i in the liquid phase

! ************************************************************************************************ !
! CHARACTER STRINGS                                                                                !
! ************************************************************************************************ !
CHARACTER( Len= 01 ) :: CurveType ! Isotherm types (A, B, or C)

! ************************************************************************************************ !
! LOGICAL VARIABLES                                                                                !
! ************************************************************************************************ !
LOGICAL, DIMENSION( 4 ) :: FluidPhase, PhaseTest ! Phase type (1: Supercritical fluid; 2: Critical fluid; 3: Subcritical liquid phase; 4: Subcritical vapor phase)

! ************************************************************************************************ !
! Calculating the chemical potential of species i in the vapor phase                              !
! ************************************************************************************************ !

! Initialization
FluidPhase    = .FALSE.
FluidPhase(4) = .TRUE. ! Subcritical vapor phase
PhaseTest     = FluidPhase

! Calculate the molar volume of the vapor phase
CALL Topliss_Algorithm( 1_Int64, VaporComposition, Temperature, Pressure, mVolumeVaporPhase, FluidPhase, CurveType, .FALSE. )
IF( ALL( PhaseTest ) .NEQV. ALL( FluidPhase ) ) THEN ! Should not fail
  WRITE( *, "(G0)" ) "Error: 'PhaseTest' and 'FluidPhase' logical arrays are not equal. Exiting..."
  CALL EXIT(  )
END IF

! Calculate compressibility factor of the vapor phase
CALL Calculate_Pressure( VaporComposition, mVolumeVaporPhase, Temperature, Anynumber, IsothermalCompressibility, &
&                        ThermalExpansionCoefficient, CompressibilityFactor )

! Calculate residual chemical potential of species i in the vapor phase
CALL Calculate_CPotential( VaporComposition, mVolumeVaporPhase, Temperature, TotalChemicalPotential, ResidualChemicalPotential )

! Calculate fugacity of species i in the vapor phase
cFugacityVaporPhase = ( ResidualChemicalPotential / ( cUniversalGas * Temperature ) ) - DLOG( CompressibilityFactor )

! ************************************************************************************************ !
! Calculating the chemical potential of species i in the liquid phase                              !
! ************************************************************************************************ !

! Initialization
FluidPhase    = .FALSE.
FluidPhase(3) = .TRUE. ! Subcritical liquid phase
PhaseTest     = FluidPhase

! Calculate the molar volume of the liquid phase
CALL Topliss_Algorithm( 1_Int64, LiquidComposition, Temperature, Pressure, mVolumeLiquidPhase, FluidPhase, CurveType, .FALSE. )
IF( ALL( PhaseTest ) .NEQV. ALL( FluidPhase ) ) THEN ! Should not fail
  WRITE( *, "(G0)" ) "Error: 'PhaseTest' and 'FluidPhase' logical arrays are not equal. Exiting..."
  CALL EXIT(  )
END IF

! Calculate compressibility factor of the liquid phase
CALL Calculate_Pressure( LiquidComposition, mVolumeLiquidPhase, Temperature ,Anynumber, IsothermalCompressibility, &
&                        ThermalExpansionCoefficient, CompressibilityFactor )

! Calculate residual chemical potential of species i in the liquid phase
CALL Calculate_CPotential( LiquidComposition, mVolumeLiquidPhase, Temperature, TotalChemicalPotential, ResidualChemicalPotential )

! Calculate fugacity of species i in the liquid phase
cFugacityLiquidPhase = ( ResidualChemicalPotential / ( cUniversalGas * Temperature ) ) - DLOG( CompressibilityFactor )

! ************************************************************************************************ !
! Calculating the difference between old and new equilibrium factors                               !
! ************************************************************************************************ !

! Calculating the new equilibrium factor
NewEquilibriumFactor = DEXP( cFugacityLiquidPhase - cFugacityVaporPhase )

! Evaluating the error between the old and new equilibrium factors
Error = SUM( DABS( NewEquilibriumFactor - OldEquilibriumFactor ) )

! Updating the equilibrium factors
OldEquilibriumFactor = NewEquilibriumFactor

RETURN

END SUBROUTINE Equilibrium_Factor_Calculation