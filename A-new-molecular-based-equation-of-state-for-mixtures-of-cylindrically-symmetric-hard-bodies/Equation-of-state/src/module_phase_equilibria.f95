! ************************************************************************************************ !
!                                     VAPOR-LIQUID EQUILIBRIUM                                     !
! ************************************************************************************************ !
!               This subroutine calculates the coexistence curve of a pure component               !
!                            or the vapor-liquid envelope of a mixture.                            !
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
! Coexistence curve                                                                                !
! ************************************************************************************************ !
SUBROUTINE Coexistence_Curve(  )

USE GlobalVariables
USE Substances

IMPLICIT NONE

! ************************************************************************************************ !
! INTEGER VARIABLES                                                                                !
! ************************************************************************************************ !
INTEGER( Kind= Int64 ) :: cComponent ! Counter (components)
INTEGER( Kind= Int64 ) :: Counter    ! Counter (generic)
INTEGER( Kind= Int64 ) :: dCounter   ! Counter (damping)

! ************************************************************************************************ !
! REAL PARAMETERS                                                                                  !
! ************************************************************************************************ !
REAL( Kind= Real64 ), PARAMETER :: Tolerance = 1.D-8 ! Tolerance

! ************************************************************************************************ !
! REAL VARIABLES                                                                                   !
! ************************************************************************************************ !
REAL( Kind= Real64 )                           :: InitialTemperature       ! Temperature (initial)
REAL( Kind= Real64 )                           :: FinalTemperature         ! Temperature (final)
REAL( Kind= Real64 )                           :: MaxTemperature           ! Temperature (maximum)
REAL( Kind= Real64 )                           :: StepTemperature          ! Temperature (increment)
REAL( Kind= Real64 )                           :: Temperature              ! Temperature
REAL( Kind= Real64 )                           :: MinimumPressure          ! Pressure (minimum)
REAL( Kind= Real64 )                           :: MaximumPressure          ! Pressure (maximum)
REAL( Kind= Real64 )                           :: MinimumDensity           ! Molar density (minimum)
REAL( Kind= Real64 )                           :: MaximumDensity           ! Molar density (maximum)
REAL( Kind= Real64 )                           :: aPressure, bPressure     ! Pressure interval
REAL( Kind= Real64 )                           :: MidpointPressure         ! Pressure (midpoint)
REAL( Kind= Real64 )                           :: LiquidVolume             ! Molar volume (liquid phase)
REAL( Kind= Real64 )                           :: VaporVolume              ! Molar volume (vapor phase)
REAL( Kind= Real64 )                           :: LiquidChemicalPotential  ! Chemical potential (liquid phase)
REAL( Kind= Real64 )                           :: VaporChemicalPotential   ! Chemical potential (vapor phase)
REAL( Kind= Real64 )                           :: rLiquidChemicalPotential ! Residual chemical potential (liquid phase)
REAL( Kind= Real64 )                           :: rVaporChemicalPotential  ! Residual chemical potential (vapor phase)
REAL( Kind= Real64 )                           :: FugacityRatio            ! Fugacity ratio
REAL( Kind= Real64 )                           :: FugacityError            ! Fugacity convergence error
REAL( Kind= Real64 )                           :: DampingFactor            ! Damping factor
REAL( Kind= Real64 )                           :: DampingFactorIteration   ! Damping factor (iteration)
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: mFraction                ! Molar fraction of a component
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: cCriticalVolume          ! Molar volume (critical)
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: cCriticalDensity         ! Molar density (critical)
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: cCriticalPressure        ! Pressure (critical)
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: cCriticalTemperature     ! Temperature (critical)

! ************************************************************************************************ !
! CHARACTER STRINGS                                                                                !
! ************************************************************************************************ !
CHARACTER( Len= 99 ) :: Filename  ! File name
CHARACTER( Len= 01 ) :: CurveType ! Isotherm types (A, B, or C)
CHARACTER( Len= 01 ) :: Dummy     ! Dummy

! ************************************************************************************************ !
! LOGICAL VARIABLES                                                                                !
! ************************************************************************************************ !
LOGICAL                 :: FolderExists ! Checks whether a folder exists
LOGICAL, DIMENSION( 4 ) :: FluidPhase   ! Phase type (1: Supercritical fluid; 2: Critical fluid; 3: Subcritical liquid phase; 4: Subcritical vapor phase)

! Status
WRITE( *, "(G0)" ) "Coexistence Curve Route"
WRITE( *, "(G0)" ) " "

! Check parent subfolders
INQUIRE( File= "ROUTE_02_COEXISTENCE_CURVE/", Exist= FolderExists )
IF( .NOT. FolderExists ) THEN
  CALL SYSTEM( "mkdir ROUTE_02_COEXISTENCE_CURVE/" )
END IF
INQUIRE( File= "ROUTE_02_COEXISTENCE_CURVE/"//TRIM( DescriptorDate )//"/", Exist= FolderExists )
IF( .NOT. FolderExists ) THEN
  CALL SYSTEM( "mkdir ROUTE_02_COEXISTENCE_CURVE/"//TRIM( DescriptorDate )//"/" )
END IF

! File name
Filename = "ROUTE_02_COEXISTENCE_CURVE/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )//"_"//&
&          TRIM( cMoleculeName(1) )//".dat"

! Status
WRITE( *, "(G0)" ) "Finding critical properties of "//TRIM( cMoleculeName(1) )//"..."
WRITE( *, "(G0)" ) " "

! Find the critical properties of components
cComponent = 1_INT64
CALL Critical_Properties_Pure_Components( cComponent, cCriticalTemperature(cComponent), cCriticalDensity(cComponent), &
&                                         cCriticalPressure(cComponent), cCriticalVolume(cComponent) )

! Create output file
OPEN( Unit= 20, File= Filename )
WRITE( 20, "(G0)" ) "Geometry: "//TRIM( cGeometryName(cComponent) )//""
WRITE( *, "(G0)" ) "Geometry: "//TRIM( cGeometryName(cComponent) )//""
WRITE( 20, "(G0)" ) " "
WRITE( *, "(G0)" ) " "
IF( cGeometry(cComponent) == 1 ) THEN
  WRITE( 20, "(G0)" ) "EOR model: Isihara-Hadwiger Theorem"
  WRITE( *, "(G0)" ) "EOR model: Isihara-Hadwiger Theorem"
ELSE IF( cGeometry(cComponent) == 2 ) THEN
  WRITE( 20, "(G0)" ) "SPC model: Isihara-Hadwiger Theorem"
  WRITE( *, "(G0)" ) "SPC model: Isihara-Hadwiger Theorem"
ELSE IF( cGeometry(cComponent) == 3 ) THEN
  WRITE( 20, "(G0)" ) "CYL model: Isihara-Hadwiger Theorem"
  WRITE( *, "(G0)" ) "CYL model: Isihara-Hadwiger Theorem"
END IF
WRITE( 20, "(G0)" ) " "
WRITE( *, "(G0)" ) " "
IF( ZhangCorrectionLogical ) WRITE( 20, "(G0)" ) "Zhang Correction Applied"
IF( ZhangCorrectionLogical ) WRITE( *, "(G0)" ) "Zhang Correction Applied"
IF( .NOT. ZhangCorrectionLogical ) WRITE( 20, "(G0)" ) "Zhang Correction NOT Applied"
IF( .NOT. ZhangCorrectionLogical ) WRITE( *, "(G0)" ) "Zhang Correction NOT Applied"
WRITE( 20, "(G0)" ) " "
WRITE( *, "(G0)" ) " "
IF( PotentialTypeLogical(4) ) THEN
  IF( PYHCBCorrectionLogical ) WRITE( 20, "(G0)" ) "Percus-Yevick Correction Applied"
  IF( PYHCBCorrectionLogical ) WRITE( *, "(G0)" ) "Percus-Yevick Correction Applied"
  IF( .NOT. PYHCBCorrectionLogical ) WRITE( 20, "(G0)" ) "Percus-Yevick Correction NOT Applied"
  IF( .NOT. PYHCBCorrectionLogical ) WRITE( *, "(G0)" ) "Percus-Yevick Correction NOT Applied"
END IF
WRITE( 20, "(G0)" ) " "
WRITE( *, "(G0)" ) " "
IF( PotentialTypeLogical(4) ) THEN
  IF( UseA1ForA2Logical ) WRITE( 20, "(G0)" ) "Overriding effective packing fraction coefficients for A2"
  IF( UseA1ForA2Logical ) WRITE( *, "(G0)" ) "Overriding effective packing fraction coefficients for A2"
  IF( .NOT. UseA1ForA2Logical ) WRITE( 20, "(G0)" ) "Keeping effective packing fraction coefficients for A2"
  IF( .NOT. UseA1ForA2Logical ) WRITE( *, "(G0)" ) "Keeping effective packing fraction coefficients for A2"
END IF
WRITE( 20, "(G0)" ) " "
WRITE( *, "(G0)" ) " "
WRITE( 20, "(G0,G0.10,G0)" ) "Critical Temperature: ", cCriticalTemperature(cComponent), " K"
WRITE( 20, "(G0,G0.10,G0)" ) "Critical Pressure: ", cCriticalPressure(cComponent) / 1.D6, " MPa"
WRITE( 20, "(G0,G0.10,G0)" ) "Critical Molar Density: ", cCriticalDensity(cComponent), " mol/m³"
WRITE( 20, "(G0,G0.10,G0)" ) "Critical Molar Volume: ", cCriticalVolume(cComponent), " m³/mol"
WRITE( 20, "(G0)" ) " "
WRITE( 20, "(6(G0,','))", Advance= "NO" ) "'Temperature [K]'", "'Minimum Pressure [MPa]'", "'Saturation Pressure [MPa]'", &
&                                         "'Maximum Pressure [MPa]'", "'Liq. Molar Density [mol/m³]'", &
&                                         "'Vap. Molar Density [mol/m³]'"
WRITE( 20, "(G0,',',G0,',',G0)", Advance= "YES" ) "'Fugacity Ratio'", "'Damping Factor'", "'Iterations'"
FLUSH( 20 )

! Critical properties
WRITE( *, "(G0,G0.5,G0)" ) "Critical Temperature: ", cCriticalTemperature(cComponent), "K"
WRITE( *, "(G0,G0.5,G0)" ) "Critical Pressure: ", cCriticalPressure(cComponent) / 1.D6, "MPa"
WRITE( *, "(G0,G0.5,G0)" ) "Critical Density: ", cCriticalDensity(cComponent), "mol/m³"
WRITE( *, "(G0,G0.5,G0)" ) "Critical Volume: ", cCriticalVolume(cComponent), "m³/mol"
WRITE( *, "(G0)" ) " "

! Initial temperature
OPEN( Unit= 10, File= "process_2_specs.ini", Action= "READ" )
READ( 10, * ) Dummy, Temperature
IF( Temperature <= 0.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "Error: the initial temperature [", Temperature,"K] must be greater than zero. Exiting..."
  CLOSE( 10 )
  CALL EXIT(  )
ELSE IF( Temperature >= cCriticalTemperature(cComponent) ) THEN
  WRITE( *, "(G0,G0.5,2G0,G0.5,G0)" ) "Error: the initial temperature [", Temperature,"K] must be less than the critical ", &
  &                                   "temperature [", cCriticalTemperature(cComponent),"K]. Exiting..."
  CLOSE( 10 )
  CALL EXIT(  )
END IF
READ( 10, * ) Dummy, StepTemperature
IF( StepTemperature <= 0.D0 ) THEN
  WRITE( *, "(G0,G0.5,G0)" ) "Error: the temperature increment [", StepTemperature,"K] must be greater than zero. Exiting..."
  CLOSE( 10 )
  CALL EXIT(  )
END IF
READ( 10, * ) Dummy, DampingFactorIteration
IF( DampingFactorIteration <= 0.D0 .OR. DampingFactorIteration > 1.D0 ) THEN
  WRITE( *, "(G0,G0.5,2G0)" ) "Error: the damping factor increment [", DampingFactorIteration,"] must be greater than zero and ", &
  &                           "less than one. Exiting..."
  CLOSE( 10 )
  CALL EXIT(  )
END IF
CLOSE( 10 )

! Initial temperature
InitialTemperature = Temperature

WRITE( *, "(G0)" ) "Running the coexistence curve route..."
WRITE( *, "(G0)" ) " "

! Maximum temperature
MaxTemperature = cCriticalTemperature(cComponent)

! Initialization
mFraction = 1.D0

! Coexistence curve
ConvergenceCheck: DO
  DO WHILE( Temperature <= MaxTemperature )
    ! Initialization
    FinalTemperature = Temperature
    DampingFactor    = 1.D0
    Counter          = 0
    dCounter         = 0
    ! Status
    CALL Progress_Bar_Path_02( Temperature )
    ! Find pressure interval
    CALL Find_Pressure_Interval( cComponent, mFraction, Temperature, MinimumPressure, MinimumDensity, MaximumPressure, &
    &                            MaximumDensity, .TRUE. )
    ! Pressure interval
    aPressure = MinimumPressure
    bPressure = MaximumPressure
    ! Change minimum pressure to 0 when minimum pressure is negative
    IF( aPressure < 0.D0 ) aPressure = 0.D0
    ! Midpoint pressure
    MidpointPressure = 0.5D0 * (aPressure + bPressure)
    ! Liquid Phase
    FluidPhase(:) = .FALSE.
    FluidPhase(3) = .TRUE.
    ! Find volume of the liquid phase
    CALL Topliss_Algorithm( 1_Int64, mFraction, Temperature, MidpointPressure, LiquidVolume, FluidPhase, CurveType, .TRUE. )
    ! Find chemical potential of the liquid phase
    CALL Calculate_CPotential_Single_Component( 1_Int64, LiquidVolume, Temperature, LiquidChemicalPotential, &
    &                                           rLiquidChemicalPotential )
    ! Vapor Phase
    FluidPhase(:) = .FALSE.
    FluidPhase(4) = .TRUE.
    ! Find volume of the vapor phase
    CALL Topliss_Algorithm( 1_Int64, mFraction, Temperature, MidpointPressure, VaporVolume, FluidPhase, CurveType, .TRUE. )
    ! Find chemical potential of the liquid phase
    CALL Calculate_CPotential_Single_Component( 1_Int64, VaporVolume, Temperature, VaporChemicalPotential, rVaporChemicalPotential )
    ! Fugacity ratio
    FugacityRatio = DEXP( ( LiquidChemicalPotential - VaporChemicalPotential ) / ( cUniversalGas * Temperature ) )
    ! Convergence criterion
    FugacityError = FugacityRatio - 1.D0
    DO WHILE( DABS( FugacityError ) >= Tolerance )
      ! Update midpoint pressure
      MidpointPressure = MidpointPressure * FugacityRatio
      ! Liquid Phase
      FluidPhase(:) = .FALSE.
      FluidPhase(3) = .TRUE.
      ! Find volume of the liquid phase
      CALL Topliss_Algorithm( 1_Int64, mFraction, Temperature, MidpointPressure, LiquidVolume, FluidPhase, CurveType, .TRUE. )
      ! Find chemical potential of the liquid phase
      CALL Calculate_CPotential_Single_Component( 1_Int64, LiquidVolume, Temperature, LiquidChemicalPotential, &
      &                                           rLiquidChemicalPotential )
      ! Vapor Phase
      FluidPhase(:) = .FALSE.
      FluidPhase(4) = .TRUE.
      ! Find volume of the vapor phase
      CALL Topliss_Algorithm( 1_Int64, mFraction, Temperature, MidpointPressure, VaporVolume, FluidPhase, CurveType, .TRUE. )
      ! Find chemical potential of the liquid phase
      CALL Calculate_CPotential_Single_Component( 1_Int64, VaporVolume, Temperature, VaporChemicalPotential, &
      &                                           rVaporChemicalPotential )
      ! Fugacity ratio
      FugacityRatio = DEXP( ( LiquidChemicalPotential - VaporChemicalPotential ) / ( cUniversalGas * Temperature ) ) ** &
      &               DampingFactor
      ! Convergence criterion
      FugacityError = FugacityRatio - 1.D0
      ! Update counter
      Counter  = Counter + 1
      dCounter = dCounter + 1
      ! Update damping factor
      IF( dCounter > 10 ) THEN
        dCounter = 0
        DampingFactor = DampingFactor * DampingFactorIteration
      END IF
      ! Check damping factor
      IF( DampingFactor <= 1.D-6 ) THEN
        WRITE( *, "(G0)" ) " "
        WRITE( *, "(G0)" ) " "
        WRITE( *, "(G0,G0.5,G0)" ) "Error: the damping factor is becoming too small. Convergence will probably fail. Finishing..."
        EXIT ConvergenceCheck
      END IF
      ! Check iterations
      IF( Counter >= 20000 ) THEN
        WRITE( *, "(G0)" ) " "
        WRITE( *, "(G0)" ) " "
        WRITE( *, "(G0,G0.5,G0)" ) "Error: the number of iterations is too high. Convergence will probably fail. Finishing..."
        EXIT ConvergenceCheck
      END IF
    END DO
    ! Write results
    WRITE( 20, "(6(G0,','))", Advance= "NO" ) Temperature, aPressure / 1.D6, MidpointPressure / 1.D6, bPressure / 1.D6, &
    &                                         1.D0 / LiquidVolume, 1.D0 / VaporVolume
    WRITE( 20, "(G0,',',G0,',',G0)", Advance= "YES" ) FugacityRatio, DampingFactor, Counter
    FLUSH( 20 )
    ! Update temperature
    Temperature = Temperature + StepTemperature
  END DO
  EXIT ConvergenceCheck
END DO ConvergenceCheck

! Include critical point
IF( FinalTemperature < MaxTemperature ) THEN
  WRITE( 20, "(6(G0,','))", Advance= "NO" ) MaxTemperature, cCriticalPressure(cComponent) / 1.D6, cCriticalPressure(cComponent) / &
  &                                         1.D6, cCriticalPressure(cComponent) / 1.D6, cCriticalDensity(cComponent), &
  &                                         cCriticalDensity(cComponent)
  WRITE( 20, "(G0,',',G0,',',G0)", Advance= "YES" ) 1.D0, 1.D0, 0_Int64
  FLUSH( 20 )
END IF

! Status
WRITE( *, "(G0)" ) " "
WRITE( *, "(G0)" ) " "

! Close opened unit
CLOSE( 20 )

RETURN

END SUBROUTINE Coexistence_Curve

! ************************************************************************************************ !
! Vapor-liquid envelope                                                                            !
! ************************************************************************************************ !
SUBROUTINE Binary_Vapor_Liquid_Equilibrium(  )

USE GlobalVariables
USE Substances
USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: Output_Unit

IMPLICIT NONE

! ************************************************************************************************ !
! INTEGER VARIABLES                                                                                !
! ************************************************************************************************ !
INTEGER( Kind= Int64 ) :: cComponent, nCounter ! Counter

! ************************************************************************************************ !
! REAL CONSTANTS                                                                                   !
! ************************************************************************************************ !
REAL( Kind= Real64 ), PARAMETER :: Tolerance = 1.D-9 ! Numerical tolerance

! ************************************************************************************************ !
! REAL VARIABLES                                                                                   !
! ************************************************************************************************ !
REAL( Kind= Real64 )                           :: Temperature          ! Temperature
REAL( Kind= Real64 )                           :: UserTemperature      ! Temperature (user-defined)
REAL( Kind= Real64 )                           :: mVolume              ! Molar volume
REAL( Kind= Real64 )                           :: Pressure             ! Pressure
REAL( Kind= Real64 )                           :: VaporPressure        ! Vapor pressure
REAL( Kind= Real64 )                           :: MinimumPressure      ! Pressure (minimum)
REAL( Kind= Real64 )                           :: MaximumPressure      ! Pressure (maximum)
REAL( Kind= Real64 )                           :: MinimumDensity       ! Molar density (minimum)
REAL( Kind= Real64 )                           :: MaximumDensity       ! Molar density (maximum)
REAL( Kind= Real64 )                           :: mVolumeLiquidPhase   ! Molar volume (liquid phase)
REAL( Kind= Real64 )                           :: mVolumeVaporPhase    ! Molar volume (vapor phase)
REAL( Kind= Real64 )                           :: cPotentialLiquid     ! Chemical potential (liquid phase)
REAL( Kind= Real64 )                           :: cPotentialVapor      ! Chemical potential (vapor phase)
REAL( Kind= Real64 )                           :: FugacityRatio        ! Fugacity ratio
REAL( Kind= Real64 )                           :: ErrorFugacity        ! Error
REAL( Kind= Real64 )                           :: StoreError           ! Old error
REAL( Kind= Real64 )                           :: UpdateError          ! New error
REAL( Kind= Real64 )                           :: CalculationPressure  ! System pressure
REAL( Kind= Real64 )                           :: FinalPressure        ! Final pressure (stop criterion)
REAL( Kind= Real64 )                           :: MinVaporFraction     ! Vapor fraction (minimum)
REAL( Kind= Real64 )                           :: MaxVaporFraction     ! Vapor fraction (maximum)
REAL( Kind= Real64 )                           :: VaporFraction        ! Vapor fraction
REAL( Kind= Real64 )                           :: PressureIncrement    ! Pressure increments
REAL( Kind= Real64 )                           :: Error                ! Error
REAL( Kind= Real64 )                           :: DummyNumber          ! Dummy
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: cVaporPressureAcc    ! Critical pressure (accentric factor calculation)
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: cVaporPressure       ! Critical pressure
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: cCriticalTemperature ! Critical temperature
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: cCriticalPressure    ! Critical pressure
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: cCriticalDensity     ! Critical density
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: cCriticalVolume      ! Critical volume
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: cAccentricFactor     ! Accentric factor
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: cEquilibriumFactors  ! Equilibrium factor
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: FeedComposition      ! Feed composition
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: VaporComposition     ! Vapor composition
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: LiquidComposition    ! Liquid composition

! ************************************************************************************************ !
! CHARACTER STRINGS                                                                                !
! ************************************************************************************************ !
CHARACTER( Len= 150 ) :: Filename ! File name
CHARACTER( Len= 01 )  :: Dummy    ! Dummy string

! ************************************************************************************************ !
! LOGICAL VARIABLES                                                                                !
! ************************************************************************************************ !
LOGICAL                 :: SubcooledFlag, SuperheatedFlag ! Flags for subcooled liquid and superheated vapor
LOGICAL                 :: EquilibriumFactorAboveUnity    ! Checks if any equilibrium factors are above unity
LOGICAL                 :: FolderExists                   ! Checks whether a folder exists
LOGICAL, DIMENSION( 4 ) :: FluidPhase                     ! Phase type
LOGICAL, DIMENSION( 4 ) :: PhaseTest                      ! Phase type

! Status
WRITE( *, "(G0)" ) "Vapor-Liquid Equilibrium Route"
WRITE( *, "(G0)" ) " "

! Check parent subfolders
INQUIRE( File= "ROUTE_02_COEXISTENCE_CURVE/", Exist= FolderExists )
IF( .NOT. FolderExists ) THEN
  CALL SYSTEM( "mkdir ROUTE_02_COEXISTENCE_CURVE/" )
END IF
INQUIRE( File= "ROUTE_02_COEXISTENCE_CURVE/"//TRIM( DescriptorDate )//"/", Exist= FolderExists )
IF( .NOT. FolderExists ) THEN
  CALL SYSTEM( "mkdir ROUTE_02_COEXISTENCE_CURVE/"//TRIM( DescriptorDate )//"/" )
END IF

! Critical properties
DO cComponent = 1, nComponents
  CALL Critical_Properties_Pure_Components( cComponent, cCriticalTemperature(cComponent), cCriticalDensity(cComponent), &
  &                                         cCriticalPressure(cComponent), cCriticalVolume(cComponent) )
END DO

! File name
Filename = "ROUTE_02_COEXISTENCE_CURVE/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )
! Add component names to the file name
DO cComponent = 1, nComponents
  Filename = TRIM( Filename )//"_"//TRIM( cMoleculeName(cComponent) )
END DO
! Add file extension
Filename = TRIM( Filename )//".dat"

! Create output file
OPEN( Unit= 20, File= Filename )

! Write component-specific information
DO cComponent = 1, nComponents
  WRITE( 20, "(G0)" ) "Component: "//TRIM( cMoleculeName(cComponent) )
  WRITE( *,  "(G0)" ) "Component: "//TRIM( cMoleculeName(cComponent) )
  WRITE( 20, "(G0)" ) "Geometry: "//TRIM( cGeometryName(cComponent) )
  WRITE( *,  "(G0)" ) "Geometry: "//TRIM( cGeometryName(cComponent) )
  IF( cGeometry(cComponent) == 1 ) THEN
    WRITE( 20, "(G0)" ) "EOR model: Isihara-Hadwiger Theorem"
    WRITE( *,  "(G0)" ) "EOR model: Isihara-Hadwiger Theorem"
  ELSE IF( cGeometry(cComponent) == 2 ) THEN
    WRITE( 20, "(G0)" ) "SPC model: Isihara-Hadwiger Theorem"
    WRITE( *,  "(G0)" ) "SPC model: Isihara-Hadwiger Theorem"
  ELSE IF( cGeometry(cComponent) == 3 ) THEN
    WRITE( 20, "(G0)" ) "CYL model: Isihara-Hadwiger Theorem"
    WRITE( *,  "(G0)" ) "CYL model: Isihara-Hadwiger Theorem"
  END IF
  WRITE( 20, "(G0,G0.10,G0)" ) "Critical Temperature: ", cCriticalTemperature(cComponent), " K"
  WRITE( 20, "(G0,G0.10,G0)" ) "Critical Pressure: ", cCriticalPressure(cComponent) / 1.D6, " MPa"
  WRITE( 20, "(G0,G0.10,G0)" ) "Critical Molar Density: ", cCriticalDensity(cComponent), " mol/m³"
  WRITE( 20, "(G0,G0.10,G0)" ) "Critical Molar Volume: ", cCriticalVolume(cComponent), " m³/mol"
  WRITE( 20, "(G0)" ) " "
  WRITE( *,  "(G0)" ) " "
END DO

! Write model information
IF( ZhangCorrectionLogical ) THEN
  WRITE( 20, "(G0)" ) "Zhang Correction Applied"
  WRITE( *,  "(G0)" ) "Zhang Correction Applied"
ELSE
  WRITE( 20, "(G0)" ) "Zhang Correction NOT Applied"
  WRITE( *,  "(G0)" ) "Zhang Correction NOT Applied"
END IF

WRITE( 20, "(G0)" ) " "
WRITE( *,  "(G0)" ) " "

IF( PotentialTypeLogical(4) ) THEN
  IF( PYHCBCorrectionLogical ) THEN
    WRITE( 20, "(G0)" ) "Percus-Yevick Correction Applied"
    WRITE( *,  "(G0)" ) "Percus-Yevick Correction Applied"
  ELSE
    WRITE( 20, "(G0)" ) "Percus-Yevick Correction NOT Applied"
    WRITE( *,  "(G0)" ) "Percus-Yevick Correction NOT Applied"
  END IF
  WRITE( 20, "(G0)" ) " "
  WRITE( *,  "(G0)" ) " "
  IF( UseA1ForA2Logical ) THEN
    WRITE( 20, "(G0)" ) "Overriding effective packing fraction coefficients for A2"
    WRITE( *,  "(G0)" ) "Overriding effective packing fraction coefficients for A2"
  ELSE
    WRITE( 20, "(G0)" ) "Keeping effective packing fraction coefficients for A2"
    WRITE( *,  "(G0)" ) "Keeping effective packing fraction coefficients for A2"
  END IF
END IF

WRITE( 20, "(G0)" ) " "
WRITE( *,  "(G0)" ) " "

! Output-file header
WRITE( 20, "(7(G0,:,','))" ) "'Temperature [K]'",                    &
&                            "'Pressure [MPa]'",                     &
&                            "'Vapor Composition of Volatile'",      &
&                            "'Liquid Composition of Volatile'",     &
&                            "'Equilibrium Factor of Volatile'",     &
&                            "'Error'",                              &
&                            "'Iterations'"
FLUSH( 20 )

! Temperature
WRITE( *, "(G0)" ) "Enter temperature of the mixture (in K):"
READ( *, * ) UserTemperature

! Calculate the saturation pressure of each component at the user-specified temperature
DO cComponent = 1, nComponents
  Temperature = UserTemperature
  IF( Temperature > cCriticalTemperature(cComponent) ) THEN
    cVaporPressure(cComponent) = cCriticalPressure(cComponent)
    CYCLE
  END IF
  CALL Find_Pressure_Interval( cComponent, [0.5D0, 0.5D0], Temperature, MinimumPressure, MinimumDensity, MaximumPressure, &
  &                            MaximumDensity, .TRUE. )
  ! Pressure (initial guess)
  VaporPressure = 0.5D0 * ( MinimumPressure + MaximumPressure )
  IF( VaporPressure < 0.D0 ) THEN
    MinimumPressure = 0.D0
    VaporPressure = 0.5D0 * ( MinimumPressure + MaximumPressure )
  END IF
  FluidPhase    = .FALSE.
  FluidPhase(3) = .TRUE.
  PhaseTest     = FluidPhase
  CALL Topliss_Algorithm( cComponent, [0.5D0, 0.5D0], Temperature, VaporPressure, mVolume, FluidPhase, Dummy, .TRUE. )
  IF( ALL( PhaseTest ) .NEQV. ALL( FluidPhase ) ) THEN
    STOP "ERROR: PhaseTest and FluidPhase are not equal."
  END IF
  mVolumeLiquidPhase = mVolume
  FluidPhase    = .FALSE.
  FluidPhase(4) = .TRUE.
  PhaseTest     = FluidPhase
  CALL Topliss_Algorithm( cComponent, [0.5D0, 0.5D0], Temperature, VaporPressure, mVolume, FluidPhase, Dummy, .TRUE. )
  IF( ALL( PhaseTest ) .NEQV. ALL( FluidPhase ) ) THEN
    STOP "ERROR: PhaseTest and FluidPhase are not equal."
  END IF
  mVolumeVaporPhase = mVolume
  ! Chemical potential (liquid phase)
  CALL Calculate_CPotential_Single_Component( cComponent, mVolumeLiquidPhase, Temperature, cPotentialLiquid, DummyNumber )
  CALL Calculate_CPotential_Single_Component( cComponent, mVolumeVaporPhase, Temperature, cPotentialVapor, DummyNumber )
  FugacityRatio = DEXP( ( cPotentialLiquid - cPotentialVapor ) / ( cUniversalGas * Temperature ) )
  ErrorFugacity = FugacityRatio - 1.D0
  ! Initialization
  nCounter = 0
  ! Saturation pressure
  DO WHILE( DABS( ErrorFugacity ) >= Tolerance  .AND. nCounter < 100 )
    ! Store old error
    StoreError = ErrorFugacity
    ! Change pressure based on the fugacity ratio
    VaporPressure = VaporPressure * FugacityRatio
    FluidPhase    = .FALSE.
    FluidPhase(3) = .TRUE.
    PhaseTest     = FluidPhase
    CALL Topliss_Algorithm( cComponent, [0.5D0, 0.5D0], Temperature, VaporPressure, mVolume, FluidPhase, Dummy, .TRUE. )
    IF( ALL( PhaseTest ) .NEQV. ALL( FluidPhase ) ) THEN
      STOP "ERROR: PhaseTest and FluidPhase are not equal."
    END IF
    mVolumeLiquidPhase = mVolume
    FluidPhase    = .FALSE.
    FluidPhase(4) = .TRUE.
    PhaseTest     = FluidPhase
    CALL Topliss_Algorithm( cComponent, [0.5D0, 0.5D0], Temperature, VaporPressure, mVolume, FluidPhase, Dummy, .TRUE. )
    IF( ALL( PhaseTest ) .NEQV. ALL( FluidPhase ) ) THEN
      STOP "ERROR: PhaseTest and FluidPhase are not equal."
    END IF
    mVolumeVaporPhase = mVolume
    ! Chemical potential (liquid phase)
    CALL Calculate_CPotential_Single_Component( cComponent, mVolumeLiquidPhase, Temperature, cPotentialLiquid, DummyNumber )
    CALL Calculate_CPotential_Single_Component( cComponent, mVolumeVaporPhase, Temperature, cPotentialVapor, DummyNumber )
    ! New fugacity ratio
    FugacityRatio = DEXP( (cPotentialLiquid - cPotentialVapor) / (cUniversalGas * Temperature) )
    ! Update error
    ErrorFugacity = FugacityRatio - 1.D0
    ! Update counter
    nCounter = nCounter + 1
    ! Store new error
    UpdateError = ErrorFugacity
  END DO
  ! Vapor pressure
  cVaporPressure(cComponent) = VaporPressure
END DO

! Calculate the saturation pressure at Tr = 0.7 for the acentric-factor calculation
DO cComponent = 1, nComponents
  Temperature = 0.7D0 * cCriticalTemperature(cComponent)
  Pressure = cCriticalPressure(cComponent)
  mVolume = cCriticalVolume(cComponent)
  CALL Find_Pressure_Interval( cComponent, [0.5D0, 0.5D0], Temperature, MinimumPressure, MinimumDensity, MaximumPressure, &
  &                            MaximumDensity, .TRUE. )
  ! Pressure (initial guess)
  VaporPressure = 0.5D0 * ( MinimumPressure + MaximumPressure )
  IF( VaporPressure < 0.D0 ) THEN
    MinimumPressure = 0.D0
    VaporPressure = 0.5D0 * ( MinimumPressure + MaximumPressure )
  END IF
  FluidPhase    = .FALSE.
  FluidPhase(3) = .TRUE.
  PhaseTest     = FluidPhase
  CALL Topliss_Algorithm( cComponent, [0.5D0, 0.5D0], Temperature, VaporPressure, mVolume, FluidPhase, Dummy, .TRUE. )
  IF( ALL( PhaseTest ) .NEQV. ALL( FluidPhase ) ) THEN
   STOP "ERROR: PhaseTest and FluidPhase are not equal."
  END IF
  mVolumeLiquidPhase = mVolume
  FluidPhase    = .FALSE.
  FluidPhase(4) = .TRUE.
  PhaseTest     = FluidPhase
  CALL Topliss_Algorithm( cComponent, [0.5D0, 0.5D0], Temperature, VaporPressure, mVolume, FluidPhase, Dummy, .TRUE. )
  IF( ALL( PhaseTest ) .NEQV. ALL( FluidPhase ) ) THEN
   STOP "ERROR: PhaseTest and FluidPhase are not equal."
  END IF
  mVolumeVaporPhase = mVolume
  ! Chemical potential (liquid phase)
  CALL Calculate_CPotential_Single_Component( cComponent, mVolumeLiquidPhase, Temperature, cPotentialLiquid, DummyNumber )
  CALL Calculate_CPotential_Single_Component( cComponent, mVolumeVaporPhase, Temperature, cPotentialVapor, DummyNumber )
  FugacityRatio = DEXP( ( cPotentialLiquid - cPotentialVapor ) / ( cUniversalGas * Temperature ) )
  ErrorFugacity = FugacityRatio - 1.D0
  ! Initialization
  nCounter = 0
  ! Saturation pressure
  DO WHILE( DABS( ErrorFugacity ) >= Tolerance .AND. nCounter < 100 )
    ! Store old error
    StoreError = ErrorFugacity
    ! Change pressure based on the fugacity ratio
    VaporPressure = VaporPressure * FugacityRatio
    FluidPhase    = .FALSE.
    FluidPhase(3) = .TRUE.
    PhaseTest     = FluidPhase
    CALL Topliss_Algorithm( cComponent, [0.5D0, 0.5D0], Temperature, VaporPressure, mVolume, FluidPhase, Dummy, .TRUE. )
    IF( ALL( PhaseTest ) .NEQV. ALL( FluidPhase ) ) THEN
      STOP "ERROR: PhaseTest and FluidPhase are not equal."
    END IF
    mVolumeLiquidPhase = mVolume
    FluidPhase    = .FALSE.
    FluidPhase(4) = .TRUE.
    PhaseTest     = FluidPhase
    CALL Topliss_Algorithm( cComponent, [0.5D0, 0.5D0], Temperature, VaporPressure, mVolume, FluidPhase, Dummy, .TRUE. )
    IF( ALL( PhaseTest ) .NEQV. ALL( FluidPhase ) ) THEN
      STOP "ERROR: PhaseTest and FluidPhase are not equal."
    END IF
    mVolumeVaporPhase = mVolume
    ! Chemical potential (liquid phase)
    CALL Calculate_CPotential_Single_Component( cComponent, mVolumeLiquidPhase, Temperature, cPotentialLiquid, DummyNumber )
    CALL Calculate_CPotential_Single_Component( cComponent, mVolumeVaporPhase, Temperature, cPotentialVapor, DummyNumber )
    ! New fugacity ratio
    FugacityRatio = DEXP( (cPotentialLiquid - cPotentialVapor) / (cUniversalGas * Temperature) )
    ! Update error
    ErrorFugacity = FugacityRatio - 1.D0
    ! Update counter
    nCounter = nCounter + 1
    ! Store new error
    UpdateError = ErrorFugacity
  END DO
  ! Vapor pressure (accentric factor calculation)
  cVaporPressureAcc(cComponent) = VaporPressure
END DO

! Sort in vapor pressures ascending order
CALL Bubble_Sort_VLE( nComponents, cVaporPressure )

! Accentric factor
cAccentricFactor = - 1.D0 - DLOG10( cVaporPressureAcc / cCriticalPressure )

! Initialization
Temperature = UserTemperature
CalculationPressure = cVaporPressure(1) * 1.001D0

! Feed composition
WRITE( *, * ) ' '
WRITE( *, "(G0)" ) "Enter initial feed composition of the most volatile component: "
READ( *, * ) FeedComposition(1)
FeedComposition(2) = 1.D0 - FeedComposition(1)

! Final pressure
WRITE( *, * ) ' '
WRITE( *, "(G0)" ) "Enter upper pressure for stopping calculation (in MPa): "
READ( *, * ) FinalPressure

! Pressure increment
WRITE( *, * ) ' '
WRITE( *, "(G0)" ) "Enter a pressure increment (in MPa): "
READ( *, * ) PressureIncrement

! Equilibrium factors (Wilson's correlation)
cEquilibriumFactors(:) = DLOG( cCriticalPressure(:) / CalculationPressure ) + 5.373D0 * ( 1.D0 + cAccentricFactor(:) ) * ( 1.D0 - &
&                        cCriticalTemperature(:) / Temperature )
cEquilibriumFactors = DEXP( cEquilibriumFactors )

! Screen header
WRITE( Output_Unit, "(A)" ) "  T [K]       P [MPa]  y          x           K          Error"

! Pressure loop
DO
  ! Initialization
  Error    = 1.D0
  nCounter = 0
  ! Successive substitution method
  DO WHILE( ( DABS( Error ) >= 1.D-5 ) .AND. ( nCounter <= 200 ) )
    ! Testing for a subcooled liquid
    CALL Subcooled_Liquid_Test( Temperature, CalculationPressure, FeedComposition, cEquilibriumFactors, VaporComposition, &
    &                           LiquidComposition, SubcooledFlag )
    ! Not subcooled liquid
    IF( .NOT. SubcooledFlag ) THEN
      ! Testing for a superheated vapor
      CALL Superheated_Vapor_Test( Temperature, CalculationPressure, FeedComposition, cEquilibriumFactors, VaporComposition, &
      &                            LiquidComposition, SuperheatedFlag )
      ! Not superheated vapor
      IF( .NOT. SuperheatedFlag ) THEN
        ! Set the initial vapor-fraction interval
        MinVaporFraction = 0.D0
        MaxVaporFraction = 1.D0
        ! Rachford-Rice equation
        CALL Brent_Method_Rachford_Rice( MinVaporFraction, MaxVaporFraction, FeedComposition, cEquilibriumFactors, VaporFraction )
        ! Calculate the liquid- and vapor-phase compositions
        LiquidComposition(:) = FeedComposition(:) / ( 1.D0 - VaporFraction + VaporFraction * cEquilibriumFactors(:) )
        VaporComposition(:) = cEquilibriumFactors(:) * LiquidComposition(:)
      END IF
    END IF
    ! Calculate the equilibrium factors and convergence error
    CALL Equilibrium_Factor_Calculation( Temperature, CalculationPressure, VaporComposition, LiquidComposition, &
    &                                    cEquilibriumFactors, Error )
    ! Update the iteration counter
    nCounter = nCounter + 1
  END DO
  ! Write results
  WRITE( 20, "(6(G0,','))", Advance= "NO" ) Temperature, CalculationPressure / 1.D6, VaporComposition(1), &
  &                                         LiquidComposition(1), cEquilibriumFactors(1), Error
  WRITE( 20, "(I0)", Advance= "YES" ) nCounter
  FLUSH( 20 )
  WRITE( Output_Unit, "(A,F9.3,2X,F9.4,2X,F9.6,2X,F9.6,2X,F9.5,2X,ES11.3)", Advance= "NO" ) CHAR(13), Temperature, &
  &     CalculationPressure / 1.D6, VaporComposition(1), LiquidComposition(1), cEquilibriumFactors(1), Error
  FLUSH( Output_Unit )
  ! Increment pressure
  CalculationPressure = CalculationPressure + PressureIncrement * 1.D6
  ! Recalculate feed composition
  FeedComposition(1) = VaporComposition(1) * 0.999D0
  FeedComposition(2) = 1.D0 - FeedComposition(1)
  ! Stop condition
  IF( CalculationPressure / 1.D6 >= FinalPressure ) EXIT
END DO

! Status
WRITE( *, * ) ' '
WRITE( *, * ) ' '

RETURN

END SUBROUTINE Binary_Vapor_Liquid_Equilibrium