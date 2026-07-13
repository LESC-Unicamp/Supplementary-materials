! ************************************************************************************************ !
!                                     THERMODYNAMIC PROCESSES                                      !
! ************************************************************************************************ !
!          This subroutine calculates the thermodynamic properties of the system through           !
!            three possible processes: [1] ISOTHERMAL, [2] ISOBARIC, and [3] ISOCHORIC             !
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
SUBROUTINE Thermodynamic_Properties( PureComponent )

USE GlobalVariables
USE Substances

IMPLICIT NONE

! ************************************************************************************************ !
! INTEGER VARIABLES                                                                                !
! ************************************************************************************************ !
INTEGER( Kind= Int64 ) :: iComponent           ! Counter (components)
INTEGER( Kind= Int64 ) :: cComponent           ! Counter (components)
INTEGER( Kind= Int64 ) :: Counter              ! Counter
INTEGER( Kind= Int64 ) :: PhaseCounter         ! Counter (phases)
INTEGER( Kind= Int64 ) :: ThermodynamicProcess ! Thermodynamic process

! ************************************************************************************************ !
! REAL VARIABLES                                                                                   !
! ************************************************************************************************ !
REAL( Kind= Real64 )                           :: LowerBoundTemperature, UpperBoundTemperature ! Temperature range (initial and final)
REAL( Kind= Real64 )                           :: StepTemperature                              ! Temperature increment
REAL( Kind= Real64 )                           :: LowerBoundPressure, UpperBoundPressure       ! Pressure range (initial and final)
REAL( Kind= Real64 )                           :: StepPressure                                 ! Pressure increment
REAL( Kind= Real64 )                           :: LowerBoundDensity, UpperBoundDensity         ! Molar density range (initial and final)
REAL( Kind= Real64 )                           :: StepDensity                                  ! Molar density increment
REAL( Kind= Real64 )                           :: Temperature                                  ! Temperature
REAL( Kind= Real64 )                           :: Pressure                                     ! Pressure
REAL( Kind= Real64 )                           :: Density                                      ! Molar density
REAL( Kind= Real64 )                           :: Volume                                       ! Molar volume
REAL( Kind= Real64 )                           :: IsothermalCompressibility                    ! Isothermal compressibility
REAL( Kind= Real64 )                           :: ThermalExpansionCoefficient                  ! Thermal expansion coefficient
REAL( Kind= Real64 )                           :: CompressibilityFactor                        ! Compressibility factor
REAL( Kind= Real64 )                           :: Entropy                                      ! Entropy
REAL( Kind= Real64 )                           :: vSpecificHeat                                ! Specific heat at constant volume
REAL( Kind= Real64 )                           :: pSpecificHeat                                ! Specific heat at constant pressure
REAL( Kind= Real64 )                           :: MassDensity                                  ! Mass density
REAL( Kind= Real64 )                           :: SoundSpeed                                   ! Speed of sound
REAL( Kind= Real64 )                           :: JTCoefficient                                ! Joule-Thompson coefficient
REAL( Kind= Real64 )                           :: Anynumber                                    ! Any number (dummy)
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: mFraction                                    ! Molar fraction of a component
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: tChemicalPotential                           ! Total chemical potential
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: rChemicalPotential                           ! Residual chemical potential

! ************************************************************************************************ !
! CHARACTER STRINGS                                                                                !
! ************************************************************************************************ !
CHARACTER( Len= 200 ) :: Filename   ! File name
CHARACTER( Len= 01  ) :: CurveType  ! Isotherm types (A, B, or C)
CHARACTER( Len= 01  ) :: Dummy      ! Dummy

! ************************************************************************************************ !
! LOGICAL VARIABLES                                                                                !
! ************************************************************************************************ !
LOGICAL                 :: PureComponent            ! Checks whether the system is a pure component or a mixture
LOGICAL                 :: FolderExists             ! Checks whether a folder exists
LOGICAL                 :: LeaveLoop01, LeaveLoop02 ! Checks when a loop must be stopped
LOGICAL, DIMENSION( 4 ) :: FluidPhase               ! Phase type (1: Supercritical fluid; 2: Critical fluid; 3: Subcritical liquid phase; 4: Supercritical vapor phase)
LOGICAL, DIMENSION( 4 ) :: PhaseTest                ! Phase type (1: Supercritical fluid; 2: Critical fluid; 3: Subcritical liquid phase; 4: Supercritical vapor phase)

! Read mixture specification
IF( nComponents > 1 ) THEN
  ! Read input file
  CALL Read_Mixture_Specification( mFraction )
  ! Normalize the mole fractions
  mFraction = mFraction / SUM( mFraction )
ELSE
  mFraction = 1.D0
END IF

! Read thermodynamic route
OPEN( Unit= 10, File= "path.ini", Action= "READ" )
READ( 10, * ) Dummy, Dummy
READ( 10, * ) Dummy, ThermodynamicProcess
IF( ThermodynamicProcess < 1 .OR. ThermodynamicProcess > 3 ) THEN
  WRITE( *, "(2G0)" ) "Error: the thermodynamic process must be 1 (ISOTHERMAL), 2 (ISOBARIC), or 3 (ISOCHORIC). Exiting..."
  CALL EXIT(  )
END IF
CLOSE( 10 )

! Isothermal route
IF( ThermodynamicProcess == 1 ) THEN
  ! Status
  WRITE( *, "(G0)" ) "Isothermal Route"
  WRITE( *, "(G0)" ) " "
  ! Read specifications
  OPEN( Unit= 10, File= "process_1_specs.ini", Action= "READ" )
  READ( 10, * ) Dummy, LowerBoundTemperature
  IF( LowerBoundTemperature < 0.D0 ) THEN
    WRITE( *, "(G0)" ) "Error: the lower bound temperature must be greater than zero. Exiting..."
    CALL EXIT(  )
  END IF
  READ( 10, * ) Dummy, UpperBoundTemperature
  IF( UpperBoundTemperature < LowerBoundTemperature ) THEN
    WRITE( *, "(G0)" ) "Error: the upper bound temperature must be greater than or equal to the lower bound temperature. Exiting..."
    CALL EXIT(  )
  END IF
  READ( 10, * ) Dummy, StepTemperature
  IF( StepTemperature < 0.D0 ) THEN
    WRITE( *, "(G0)" ) "Error: the temperature step must be greater than zero. Exiting..."
    CALL EXIT(  )
  END IF
  READ( 10, * ) Dummy, LowerBoundPressure
  IF( LowerBoundPressure < 0.D0 ) THEN
    WRITE( *, "(G0)" ) "Error: the lower bound pressure must be greater than zero. Exiting..."
    CALL EXIT(  )
  END IF
  READ( 10, * ) Dummy, UpperBoundPressure
  IF( UpperBoundPressure < LowerBoundPressure ) THEN
    WRITE( *, "(G0)" ) "Error: the upper bound pressure must be greater than or equal to the lower bound pressure. Exiting..."
    CALL EXIT(  )
  END IF
  READ( 10, * ) Dummy, StepPressure
  IF( StepPressure < 0.D0 ) THEN
    WRITE( *, "(G0)" ) "Error: the pressure step must be greater than zero. Exiting..."
    CALL EXIT(  )
  END IF
  ! Status
  WRITE( *, "(G0,G0.5,G0,G0.5,G0)" ) "Temperature Range: [", LowerBoundTemperature, ", ", UpperBoundTemperature, "] K"
  WRITE( *, "(G0,G0.5,G0,G0.5,G0)" ) "Pressure Range: [", LowerBoundPressure, ", ", UpperBoundPressure, "] MPa"
  WRITE( *, "(G0)" ) " "
  ! Conversion from MPa to Pa
  LowerBoundPressure = 1.D6 * LowerBoundPressure
  UpperBoundPressure = 1.D6 * UpperBoundPressure
  StepPressure = 1.D6 * StepPressure
  CLOSE( 10 )
! Isobaric route
ELSE IF( ThermodynamicProcess == 2 ) THEN
  ! Status
  WRITE( *, "(G0)" ) "Isobaric Route"
  WRITE( *, "(G0)" ) " "
  ! Read specifications
  OPEN( Unit= 10, File= "process_1_specs.ini", Action= "READ" )
  ! Skip lines
  DO Counter = 1, 6
    READ( 10, * ) Dummy, Dummy
  END DO
  READ( 10, * ) Dummy, LowerBoundPressure
  IF( LowerBoundPressure < 0.D0 ) THEN
    WRITE( *, "(G0)" ) "Error: the lower bound pressure must be greater than zero. Exiting..."
    CALL EXIT(  )
  END IF
  READ( 10, * ) Dummy, UpperBoundPressure
  IF( UpperBoundPressure < LowerBoundPressure ) THEN
    WRITE( *, "(G0)" ) "Error: the upper bound pressure must be greater than or equal to the lower bound pressure. Exiting..."
    CALL EXIT(  )
  END IF
  READ( 10, * ) Dummy, StepPressure
  IF( StepPressure < 0.D0 ) THEN
    WRITE( *, "(G0)" ) "Error: the pressure step must be greater than zero. Exiting..."
    CALL EXIT(  )
  END IF
  READ( 10, * ) Dummy, LowerBoundTemperature
  IF( LowerBoundTemperature < 0.D0 ) THEN
    WRITE( *, "(G0)" ) "Error: the lower bound temperature must be greater than zero. Exiting..."
    CALL EXIT(  )
  END IF
  READ( 10, * ) Dummy, UpperBoundTemperature
  IF( UpperBoundTemperature < LowerBoundTemperature ) THEN
    WRITE( *, "(G0)" ) "Error: the upper bound temperature must be greater than or equal to the lower bound temperature. Exiting..."
    CALL EXIT(  )
  END IF
  READ( 10, * ) Dummy, StepTemperature
  IF( StepTemperature < 0.D0 ) THEN
    WRITE( *, "(G0)" ) "Error: the temperature step must be greater than zero. Exiting..."
    CALL EXIT(  )
  END IF
  ! Status
  WRITE( *, "(G0,G0.5,G0,G0.5,G0)" ) "Pressure Range: [", LowerBoundPressure, ", ", UpperBoundPressure, "] MPa"
  WRITE( *, "(G0,G0.5,G0,G0.5,G0)" ) "Temperature Range: [", LowerBoundTemperature, ", ", UpperBoundTemperature, "] K"
  WRITE( *, "(G0)" ) " "
  ! Conversion from MPa to Pa
  LowerBoundPressure = 1.D6 * LowerBoundPressure
  UpperBoundPressure = 1.D6 * UpperBoundPressure
  StepPressure = 1.D6 * StepPressure
  CLOSE( 10 )
! Isochoric route
ELSE IF( ThermodynamicProcess == 3 ) THEN
  ! Status
  WRITE( *, "(G0)" ) "Isochoric Route"
  WRITE( *, "(G0)" ) " "
  ! Read specifications
  OPEN( Unit= 10, File= "process_1_specs.ini", Action= "READ" )
  ! Skip lines
  DO Counter = 1, 12
    READ( 10, * ) Dummy, Dummy
  END DO
  READ( 10, * ) Dummy, LowerBoundDensity
  IF( LowerBoundDensity < 0.D0 ) THEN
    WRITE( *, "(G0)" ) "Error: the lower bound density must be greater than zero. Exiting..."
    CALL EXIT(  )
  END IF
  READ( 10, * ) Dummy, UpperBoundDensity
  IF( UpperBoundDensity < LowerBoundDensity ) THEN
    WRITE( *, "(G0)" ) "Error: the upper bound density must be greater than or equal to the lower bound density. Exiting..."
    CALL EXIT(  )
  END IF
  READ( 10, * ) Dummy, StepDensity
  IF( StepDensity < 0.D0 ) THEN
    WRITE( *, "(G0)" ) "Error: the density step must be greater than zero. Exiting..."
    CALL EXIT(  )
  END IF
  READ( 10, * ) Dummy, LowerBoundTemperature
  IF( LowerBoundTemperature < 0.D0 ) THEN
    WRITE( *, "(G0)" ) "Error: the lower bound temperature must be greater than zero. Exiting..."
    CALL EXIT(  )
  END IF
  READ( 10, * ) Dummy, UpperBoundTemperature
  IF( UpperBoundTemperature < LowerBoundTemperature ) THEN
    WRITE( *, "(G0)" ) "Error: the upper bound temperature must be greater than or equal to the lower bound temperature. Exiting..."
    CALL EXIT(  )
  END IF
  READ( 10, * ) Dummy, StepTemperature
  IF( StepTemperature < 0.D0 ) THEN
    WRITE( *, "(G0)" ) "Error: the temperature step must be greater than zero. Exiting..."
    CALL EXIT(  )
  END IF
  ! Status
  WRITE( *, "(G0,G0.8,G0,G0.8,G0)" ) "Density Range: [", LowerBoundDensity, ", ", UpperBoundDensity, "] mol/m³"
  WRITE( *, "(G0,G0.5,G0,G0.5,G0)" ) "Temperature Range: [", LowerBoundTemperature, ", ", UpperBoundTemperature, "] K"
  WRITE( *, "(G0)" ) " "
  CLOSE( 10 )
END IF

! Check parent folders
INQUIRE( File= "ROUTE_01_THERMODYNAMIC_PROCESS/", Exist= FolderExists )
IF( .NOT. FolderExists ) THEN
  CALL SYSTEM( "mkdir ROUTE_01_THERMODYNAMIC_PROCESS/" )
END IF

! Isothermal process
IF( ThermodynamicProcess == 1 ) THEN
  ! Check parent subfolders
  INQUIRE( File= "ROUTE_01_THERMODYNAMIC_PROCESS/ISOTHERMAL/", Exist= FolderExists )
  IF( .NOT. FolderExists ) THEN
    CALL SYSTEM( "mkdir ROUTE_01_THERMODYNAMIC_PROCESS/ISOTHERMAL/" )
  END IF
  INQUIRE( File= "ROUTE_01_THERMODYNAMIC_PROCESS/ISOTHERMAL/"//TRIM( DescriptorDate )//"/", Exist= FolderExists )
  IF( .NOT. FolderExists ) THEN
    CALL SYSTEM( "mkdir ROUTE_01_THERMODYNAMIC_PROCESS/ISOTHERMAL/"//TRIM( DescriptorDate )//"/" )
  END IF
  ! File name
  IF( PureComponent ) THEN
    Filename = "ROUTE_01_THERMODYNAMIC_PROCESS/ISOTHERMAL/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )//"_"//&
    &          TRIM( cMoleculeName(1) )//".dat"
  ELSE
    Filename = "ROUTE_01_THERMODYNAMIC_PROCESS/ISOTHERMAL/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )//"_"//&
    &          TRIM( cMoleculeName(1) )
    DO iComponent = 2, nComponents - 1
      Filename = TRIM( Filename )//"_"//TRIM( cMoleculeName(iComponent) )
    END DO
    Filename = TRIM( Filename )//"_"//TRIM( cMoleculeName(nComponents) )//".dat"
  END IF
  ! Create output file
  OPEN( Unit= 20, File= Filename )
  WRITE( 20, "(G0)" ) "Isothermal Route"
  WRITE( 20, "(G0)" ) " "
  WRITE( 20, "(G0,G0.5,G0,G0.5,G0)" ) "Temperature Range: [", LowerBoundTemperature, ", ", UpperBoundTemperature, "] K"
  WRITE( 20, "(G0,G0.5,G0,G0.5,G0)" ) "Pressure Range: [", LowerBoundPressure / 1.D6, ", ", UpperBoundPressure / 1.D6, "] MPa"
  WRITE( 20, "(G0)" ) " "
  IF( nComponents == 1 ) THEN
    cComponent = 1_Int64
    WRITE( 20, "(G0)" ) "Molecule Formula: "//TRIM( cFormulaName(cComponent) )//""
    WRITE( *, "(G0)" ) "Molecule Formula: "//TRIM( cFormulaName(cComponent) )//""
    WRITE( 20, "(G0)" ) " "
    WRITE( *, "(G0)" ) " "
    WRITE( 20, "(G0)" ) "Geometry: "//TRIM( cGeometryName(cComponent) )//""
    WRITE( *, "(G0)" ) "Geometry: "//TRIM( cGeometryName(cComponent) )//""
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
  ELSE IF( nComponents > 1 ) THEN
    WRITE( 20, "(G0)", Advance= "NO" ) "Mixture: "
    WRITE( *, "(G0)", Advance= "NO" ) "Mixture: "
    DO cComponent = 1, nComponents - 1
      WRITE( 20, "(3G0)", Advance= "NO" ) ""//TRIM( cFormulaName(cComponent) )//"(", cComponent, ") + "
      WRITE( *, "(3G0)", Advance= "NO" ) ""//TRIM( cFormulaName(cComponent) )//"(", cComponent, ") + "
    END DO
    WRITE( 20, "(3G0)", Advance= "YES" ) ""//TRIM( cFormulaName(cComponent) )//"(", nComponents, ")"
    WRITE( *, "(3G0)", Advance= "YES" ) ""//TRIM( cFormulaName(cComponent) )//"(", nComponents, ")"
    WRITE( 20, "(G0)" ) " "
    WRITE( *, "(G0)" ) " "
    WRITE( 20, "(G0)", Advance= "NO" ) "Molar Fraction: "
    WRITE( *, "(G0)", Advance= "NO" ) "Molar Fraction: "
    DO cComponent = 1, nComponents - 1
      WRITE( 20, "(G0.5,3G0)", Advance= "NO" ) mFraction(cComponent), "(", cComponent, ") + "
      WRITE( *, "(G0.5,3G0)", Advance= "NO" ) mFraction(cComponent), "(", cComponent, ") + "
    END DO
    WRITE( 20, "(G0.5,3G0)", Advance= "YES" ) mFraction(cComponent), "(", nComponents, ")"
    WRITE( *, "(G0.5,3G0)", Advance= "YES" ) mFraction(cComponent), "(", nComponents, ")"
    WRITE( 20, "(G0)" ) " "
    WRITE( *, "(G0)" ) " "
    WRITE( 20, "(G0)", Advance= "NO" ) "Geometry: "
    WRITE( *, "(G0)", Advance= "NO" ) "Geometry: "
    DO cComponent = 1, nComponents - 1
      WRITE( 20, "(3G0)", Advance= "NO" ) ""//TRIM( cGeometryName(cComponent) )//"(", cComponent, ") + "
      WRITE( *, "(3G0)", Advance= "NO" ) ""//TRIM( cGeometryName(cComponent) )//"(", cComponent, ") + "
    END DO
    WRITE( 20, "(3G0)", Advance= "YES" ) ""//TRIM( cGeometryName(cComponent) )//"(", nComponents, ")"
    WRITE( *, "(3G0)", Advance= "YES" ) ""//TRIM( cGeometryName(cComponent) )//"(", nComponents, ")"
    WRITE( 20, "(G0)" ) " "
    WRITE( *, "(G0)" ) " "
    IF( ANY( cGeometry == 1 ) ) THEN
      WRITE( 20, "(G0)" ) "EOR model: Isihara-Hadwiger Theorem"
      WRITE( *, "(G0)" ) "EOR model: Isihara-Hadwiger Theorem"
      WRITE( 20, "(G0)" ) " "
      WRITE( *, "(G0)" ) " "
    END IF
    IF( ANY( cGeometry == 2 ) ) THEN
      WRITE( 20, "(G0)" ) "SPC model: Isihara-Hadwiger Theorem"
      WRITE( *, "(G0)" ) "SPC model: Isihara-Hadwiger Theorem"
      WRITE( 20, "(G0)" ) " "
      WRITE( *, "(G0)" ) " "
    END IF
    IF( ANY( cGeometry == 3 ) ) THEN
      WRITE( 20, "(G0)" ) "CYL model: Isihara-Hadwiger Theorem"
      WRITE( *, "(G0)" ) "CYL model: Isihara-Hadwiger Theorem"
      WRITE( 20, "(G0)" ) " "
      WRITE( *, "(G0)" ) " "
    END IF
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
    IF( NonSphericalMixingRule == 1 ) THEN
      WRITE( 20, "(G0)" ) "Non-Spherical Mixing Rule: Aspect Ratio Rule"
      WRITE( *, "(G0)" ) "Non-Spherical Mixing Rule: Aspect Ratio Rule"
    ELSE IF( NonSphericalMixingRule == 2 ) THEN
      WRITE( 20, "(G0)" ) "Non-Spherical Mixing Rule: Angle Average Rule"
      WRITE( *, "(G0)" ) "Non-Spherical Mixing Rule: Angle Average Rule"
    ELSE IF( NonSphericalMixingRule == 3 ) THEN
      WRITE( 20, "(G0)" ) "Non-Spherical Mixing Rule: Ratio of Second Virial Coefficients"
      WRITE( *, "(G0)" ) "Non-Spherical Mixing Rule: Ratio of Second Virial Coefficients"
    ELSE IF( NonSphericalMixingRule == 4 ) THEN
      WRITE( 20, "(G0)" ) "Non-Spherical Mixing Rule: Second Virial Coefficients"
      WRITE( *, "(G0)" ) "Non-Spherical Mixing Rule: Second Virial Coefficients"
    END IF
    WRITE( 20, "(G0)" ) " "
    WRITE( *, "(G0)" ) " "
    WRITE( 20, "(G0)" ) "Zhang Factor of the Mixture: Mixed Spherical Diameter"
    WRITE( *, "(G0)" ) "Zhang Factor of the Mixture: Mixed Spherical Diameter"
    WRITE( 20, "(G0)" ) " "
    WRITE( *, "(G0)" ) " "
    IF( EffPFractionMixingRule == 1 ) THEN
      WRITE( 20, "(G0)" ) "Effective Packing Fraction Mixing Rule: Overall Reduced Density Mixing Rule"
      WRITE( *, "(G0)" ) "Effective Packing Fraction Mixing Rule: Overall Reduced Density Mixing Rule"
    ELSE IF( EffPFractionMixingRule == 2 ) THEN
      WRITE( 20, "(G0)" ) "Effective Packing Fraction Mixing Rule: One-Fluid van der Waals Mixing Rule"
      WRITE( *, "(G0)" ) "Effective Packing Fraction Mixing Rule: One-Fluid van der Waals Mixing Rule"
    END IF
    WRITE( 20, "(G0)" ) " "
    WRITE( *, "(G0)" ) " "
  END IF
  WRITE( 20, "(13(G0,','))", Advance= "NO" ) "'Temperature [K]'", "'Pressure [MPa]'", "'Molar Density [mol/m³]'", &
  &     "'Mass Density [kg/m³]'", "'Molar Volume [m³/mol]'", "'Isothermal Compressibility [1/MPa]'", &
  &     "'Thermal Expansion Coefficient [1/K]'", "'Compressibility Factor'", "'Entropy [J/(mol.K)]'", &
  &     "'Specific Heat at Constant Volume [J/(mol.K)]'", "'Specific Heat at Constant Pressure [J/(mol.K)]'", &
  &     "'Speed of Sound [m/s]'", "'Joule-Thompson Coefficient [K/MPa]'"
  DO iComponent = 1, nComponents - 1
    WRITE( 20, "(2(G0,','))", Advance= "NO" ) "'Total Chemical Potential [J/mol] of "//TRIM( cMoleculeName(iComponent) )//"'", &
    &     "'Residual Chemical Potential [J/mol] of "//TRIM( cMoleculeName(iComponent) )//"'"
  END DO
  WRITE( 20, "(G0,',',G0)", Advance= "YES" ) "'Total Chemical Potential [J/mol] of "//TRIM( cMoleculeName(nComponents) )//"'", &
  &     "'Residual Chemical Potential [J/mol] of "//TRIM( cMoleculeName(nComponents) )//"'"
  FLUSH( 20 )
  ! Initialization
  LeaveLoop01 = .FALSE.
  Temperature = LowerBoundTemperature
  ! Temperature loop
  TemperatureLoop01: DO
    ! Initialization
    LeaveLoop02 = .FALSE.
    Pressure = LowerBoundPressure
    ! Pressure loop
    PressureLoop02: DO
      ! Progress bar
      CALL Progress_Bar_Path_01( ThermodynamicProcess, Temperature, Pressure )
      ! Phase loop
      DO PhaseCounter = 1, 4
        ! Initialization
        FluidPhase = .FALSE.
        PhaseTest  = FluidPhase
        ! Phase type
        FluidPhase(PhaseCounter) = .TRUE.
        IF( PureComponent ) mFraction = 1.D0
        ! Find volume for the given temperature, pressure, and phase
        CALL Topliss_Algorithm( 1_Int64, mFraction, Temperature, Pressure, Volume, FluidPhase, CurveType, PureComponent )
        ! Cycle phase type if the selected phase is not the same as the calculated phase
        IF( ALL( PhaseTest ) .NEQV. ALL( FluidPhase ) ) THEN
          CYCLE
        END IF
      END DO
      ! Calculate isothermal compressibility [1 / Pa], thermal expansion coefficient [1 / K], and compressibility factor
      IF( PureComponent ) THEN
        CALL Calculate_Pressure_Single_Component( 1_Int64, Volume, Temperature, Anynumber, IsothermalCompressibility, &
        &                                         ThermalExpansionCoefficient, CompressibilityFactor )
      ELSE
        CALL Calculate_Pressure( mFraction, Volume, Temperature, Anynumber, IsothermalCompressibility, &
        &                        ThermalExpansionCoefficient, CompressibilityFactor )
      END IF
      ! Calculate entropy and specific heat at constant volume [J / (mol . K)]
      IF( PureComponent ) THEN
        CALL Calculate_Entropy_Single_Component( 1_Int64, Volume, Temperature, Entropy, vSpecificHeat )
      ELSE
        CALL Calculate_Entropy( mFraction, Volume, Temperature, Entropy, vSpecificHeat )
      END IF
      ! Calculate total chemical potential and residual chemical potential [Joule / mol]
      IF( PureComponent ) THEN
        CALL Calculate_CPotential_Single_Component( 1_Int64, Volume, Temperature, tChemicalPotential(1), rChemicalPotential(1) )
      ELSE
        CALL Calculate_CPotential( mFraction, Volume, Temperature, tChemicalPotential, rChemicalPotential )
      END IF
      ! Calculate molar density [mol / m³]
      Density = 1.D0 / Volume
      ! Calculate mass density [kg / m³]
      IF( PureComponent ) THEN
        MassDensity = Density * cMolarMass(1)
      ELSE
        MassDensity = SUM( mFraction * Density * cMolarMass )
      END IF
      ! Calculate specific heat at constant pressure using Mayer's relation [J / (mol . K)]
      pSpecificHeat = vSpecificHeat + Volume * Temperature * ThermalExpansionCoefficient * ThermalExpansionCoefficient / &
      &               IsothermalCompressibility
      ! Calculate the speed of sound using Newton-Laplace equation [m / s]
      SoundSpeed = DSQRT( ( pSpecificHeat / vSpecificHeat ) * ( 1.D0 / ( IsothermalCompressibility * MassDensity ) ) )
      ! Calculate the Joule-Thompson coefficient [K / Pa]
      JTCoefficient = ( Volume / pSpecificHeat ) * ( Temperature * ThermalExpansionCoefficient - 1.D0 )
      ! Write results
      WRITE( 20, "(13(G0,','))", Advance= "NO" ) Temperature, Pressure / 1.D6, Density, MassDensity, Volume, &
      &     IsothermalCompressibility * 1.D6, ThermalExpansionCoefficient, CompressibilityFactor, Entropy, vSpecificHeat, &
      &     pSpecificHeat, SoundSpeed, JTCoefficient
      DO iComponent = 1, nComponents - 1
        WRITE( 20, "(2(G0,','))", Advance= "NO" ) tChemicalPotential(iComponent), rChemicalPotential(iComponent)
      END DO
      WRITE( 20, "(G0,',',G0)", Advance= "YES" ) tChemicalPotential(nComponents), rChemicalPotential(nComponents)
      FLUSH( 20 )
      ! Break loop condition
      IF( UpperBoundPressure == LowerBoundPressure ) EXIT PressureLoop02
      ! Increment pressure
      Pressure = Pressure + StepPressure
      ! Pressure limit
      IF( .NOT. LeaveLoop02 ) THEN
        ! Break loop condition
        IF( Pressure >= UpperBoundPressure ) THEN
          Pressure = UpperBoundPressure
          LeaveLoop02 = .TRUE.
        END IF
      ELSE
        EXIT PressureLoop02
      END IF
    END DO PressureLoop02
    ! Break loop condition
    IF( UpperBoundTemperature == LowerBoundTemperature ) EXIT TemperatureLoop01
    ! Increment temperature
    Temperature = Temperature + StepTemperature
    ! Temperature limit
    IF( .NOT. LeaveLoop01 ) THEN
      ! Break loop condition
      IF( Temperature >= UpperBoundTemperature ) THEN
        Temperature = UpperBoundTemperature
        LeaveLoop01 = .TRUE.
      END IF
    ELSE
      EXIT TemperatureLoop01
    END IF
  END DO TemperatureLoop01
  CLOSE( 10 )
END IF

! Isobaric process
IF( ThermodynamicProcess == 2 ) THEN
  ! Check parent subfolders
  INQUIRE( File= "ROUTE_01_THERMODYNAMIC_PROCESS/ISOBARIC/", Exist= FolderExists )
  IF( .NOT. FolderExists ) THEN
    CALL SYSTEM( "mkdir ROUTE_01_THERMODYNAMIC_PROCESS/ISOBARIC/" )
  END IF
  INQUIRE( File= "ROUTE_01_THERMODYNAMIC_PROCESS/ISOBARIC/"//TRIM( DescriptorDate )//"/", Exist= FolderExists )
  IF( .NOT. FolderExists ) THEN
    CALL SYSTEM( "mkdir ROUTE_01_THERMODYNAMIC_PROCESS/ISOBARIC/"//TRIM( DescriptorDate )//"/" )
  END IF
  ! File name
  IF( PureComponent ) THEN
    Filename = "ROUTE_01_THERMODYNAMIC_PROCESS/ISOBARIC/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )//"_"//&
    &          TRIM( cMoleculeName(1) )//".dat"
  ELSE
    Filename = "ROUTE_01_THERMODYNAMIC_PROCESS/ISOBARIC/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )//"_"//&
    &          TRIM( cMoleculeName(1) )
    DO iComponent = 2, nComponents - 1
      Filename = TRIM( Filename )//"_"//TRIM( cMoleculeName(iComponent) )
    END DO
    Filename = TRIM( Filename )//"_"//TRIM( cMoleculeName(nComponents) )//".dat"
  END IF
  ! Create output file
  OPEN( Unit= 20, File= Filename )
  WRITE( 20, "(G0)" ) "Isobaric Route"
  WRITE( 20, "(G0)" ) " "
  WRITE( 20, "(G0,G0.5,G0,G0.5,G0)" ) "Pressure Range: [", LowerBoundPressure / 1.D6, ", ", UpperBoundPressure / 1.D6, "] MPa"
  WRITE( 20, "(G0,G0.5,G0,G0.5,G0)" ) "Temperature Range: [", LowerBoundTemperature, ", ", UpperBoundTemperature, "] K"
  WRITE( 20, "(G0)" ) " "
  IF( nComponents == 1 ) THEN
    cComponent = 1_Int64
    WRITE( 20, "(G0)" ) "Molecule Formula: "//TRIM( cFormulaName(cComponent) )//""
    WRITE( *, "(G0)" ) "Molecule Formula: "//TRIM( cFormulaName(cComponent) )//""
    WRITE( 20, "(G0)" ) " "
    WRITE( *, "(G0)" ) " "
    WRITE( 20, "(G0)" ) "Geometry: "//TRIM( cGeometryName(cComponent) )//""
    WRITE( *, "(G0)" ) "Geometry: "//TRIM( cGeometryName(cComponent) )//""
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
  ELSE IF( nComponents > 1 ) THEN
    WRITE( 20, "(G0)", Advance= "NO" ) "Mixture: "
    WRITE( *, "(G0)", Advance= "NO" ) "Mixture: "
    DO cComponent = 1, nComponents - 1
      WRITE( 20, "(3G0)", Advance= "NO" ) ""//TRIM( cFormulaName(cComponent) )//"(", cComponent, ") + "
      WRITE( *, "(3G0)", Advance= "NO" ) ""//TRIM( cFormulaName(cComponent) )//"(", cComponent, ") + "
    END DO
    WRITE( 20, "(3G0)", Advance= "YES" ) ""//TRIM( cFormulaName(cComponent) )//"(", nComponents, ")"
    WRITE( *, "(3G0)", Advance= "YES" ) ""//TRIM( cFormulaName(cComponent) )//"(", nComponents, ")"
    WRITE( 20, "(G0)" ) " "
    WRITE( *, "(G0)" ) " "
    WRITE( 20, "(G0)", Advance= "NO" ) "Molar Fraction: "
    WRITE( *, "(G0)", Advance= "NO" ) "Molar Fraction: "
    DO cComponent = 1, nComponents - 1
      WRITE( 20, "(G0.5,3G0)", Advance= "NO" ) mFraction(cComponent), "(", cComponent, ") + "
      WRITE( *, "(G0.5,3G0)", Advance= "NO" ) mFraction(cComponent), "(", cComponent, ") + "
    END DO
    WRITE( 20, "(G0.5,3G0)", Advance= "YES" ) mFraction(cComponent), "(", nComponents, ")"
    WRITE( *, "(G0.5,3G0)", Advance= "YES" ) mFraction(cComponent), "(", nComponents, ")"
    WRITE( 20, "(G0)" ) " "
    WRITE( *, "(G0)" ) " "
    WRITE( 20, "(G0)", Advance= "NO" ) "Geometry: "
    WRITE( *, "(G0)", Advance= "NO" ) "Geometry: "
    DO cComponent = 1, nComponents - 1
      WRITE( 20, "(3G0)", Advance= "NO" ) ""//TRIM( cGeometryName(cComponent) )//"(", cComponent, ") + "
      WRITE( *, "(3G0)", Advance= "NO" ) ""//TRIM( cGeometryName(cComponent) )//"(", cComponent, ") + "
    END DO
    WRITE( 20, "(3G0)", Advance= "YES" ) ""//TRIM( cGeometryName(cComponent) )//"(", nComponents, ")"
    WRITE( *, "(3G0)", Advance= "YES" ) ""//TRIM( cGeometryName(cComponent) )//"(", nComponents, ")"
    WRITE( 20, "(G0)" ) " "
    WRITE( *, "(G0)" ) " "
    IF( ANY( cGeometry == 1 ) ) THEN
      WRITE( 20, "(G0)" ) "EOR model: Isihara-Hadwiger Theorem"
      WRITE( *, "(G0)" ) "EOR model: Isihara-Hadwiger Theorem"
      WRITE( 20, "(G0)" ) " "
      WRITE( *, "(G0)" ) " "
    END IF
    IF( ANY( cGeometry == 2 ) ) THEN
      WRITE( 20, "(G0)" ) "SPC model: Isihara-Hadwiger Theorem"
      WRITE( *, "(G0)" ) "SPC model: Isihara-Hadwiger Theorem"
      WRITE( 20, "(G0)" ) " "
      WRITE( *, "(G0)" ) " "
    END IF
    IF( ANY( cGeometry == 3 ) ) THEN
      WRITE( 20, "(G0)" ) "CYL model: Isihara-Hadwiger Theorem"
      WRITE( *, "(G0)" ) "CYL model: Isihara-Hadwiger Theorem"
      WRITE( 20, "(G0)" ) " "
      WRITE( *, "(G0)" ) " "
    END IF
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
    IF( NonSphericalMixingRule == 1 ) THEN
      WRITE( 20, "(G0)" ) "Non-Spherical Mixing Rule: Aspect Ratio Rule"
      WRITE( *, "(G0)" ) "Non-Spherical Mixing Rule: Aspect Ratio Rule"
    ELSE IF( NonSphericalMixingRule == 2 ) THEN
      WRITE( 20, "(G0)" ) "Non-Spherical Mixing Rule: Angle Average Rule"
      WRITE( *, "(G0)" ) "Non-Spherical Mixing Rule: Angle Average Rule"
    ELSE IF( NonSphericalMixingRule == 3 ) THEN
      WRITE( 20, "(G0)" ) "Non-Spherical Mixing Rule: Ratio of Second Virial Coefficients"
      WRITE( *, "(G0)" ) "Non-Spherical Mixing Rule: Ratio of Second Virial Coefficients"
    ELSE IF( NonSphericalMixingRule == 4 ) THEN
      WRITE( 20, "(G0)" ) "Non-Spherical Mixing Rule: Second Virial Coefficients"
      WRITE( *, "(G0)" ) "Non-Spherical Mixing Rule: Second Virial Coefficients"
    END IF
    WRITE( 20, "(G0)" ) " "
    WRITE( *, "(G0)" ) " "
    WRITE( 20, "(G0)" ) "Zhang Factor of the Mixture: Mixed Spherical Diameter"
    WRITE( *, "(G0)" ) "Zhang Factor of the Mixture: Mixed Spherical Diameter"
    WRITE( 20, "(G0)" ) " "
    WRITE( *, "(G0)" ) " "
    IF( EffPFractionMixingRule == 1 ) THEN
      WRITE( 20, "(G0)" ) "Effective Packing Fraction Mixing Rule: Overall Reduced Density Mixing Rule"
      WRITE( *, "(G0)" ) "Effective Packing Fraction Mixing Rule: Overall Reduced Density Mixing Rule"
    ELSE IF( EffPFractionMixingRule == 2 ) THEN
      WRITE( 20, "(G0)" ) "Effective Packing Fraction Mixing Rule: One-Fluid van der Waals Mixing Rule"
      WRITE( *, "(G0)" ) "Effective Packing Fraction Mixing Rule: One-Fluid van der Waals Mixing Rule"
    END IF
    WRITE( 20, "(G0)" ) " "
    WRITE( *, "(G0)" ) " "
  END IF
  WRITE( 20, "(13(G0,','))", Advance= "NO" ) "'Temperature [K]'", "'Pressure [MPa]'", "'Molar Density [mol/m³]'", &
  &     "'Mass Density [kg/m³]'", "'Molar Volume [m³/mol]'", "'Isothermal Compressibility [1/MPa]'", &
  &     "'Thermal Expansion Coefficient [1/K]'", "'Compressibility Factor'", "'Entropy [J/(mol.K)]'", &
  &     "'Specific Heat at Constant Volume [J/(mol.K)]'", "'Specific Heat at Constant Pressure [J/(mol.K)]'", &
  &     "'Speed of Sound [m/s]'", "'Joule-Thompson Coefficient [K/MPa]'"
  DO iComponent = 1, nComponents - 1
    WRITE( 20, "(2(G0,','))", Advance= "NO" ) "'Total Chemical Potential [J/mol] of "//TRIM( cMoleculeName(iComponent) )//"'", &
    &     "'Residual Chemical Potential [J/mol] of "//TRIM( cMoleculeName(iComponent) )//"'"
  END DO
  WRITE( 20, "(G0,',',G0)", Advance= "YES" ) "'Total Chemical Potential [J/mol] of "//TRIM( cMoleculeName(nComponents) )//"'", &
  &     "'Residual Chemical Potential [J/mol] of "//TRIM( cMoleculeName(nComponents) )//"'"
  FLUSH( 20 )
  ! Initialization
  LeaveLoop01 = .FALSE.
  Pressure = LowerBoundPressure
  ! Pressure loop
  PressureLoop01: DO
    ! Initialization
    LeaveLoop02 = .FALSE.
    Temperature = LowerBoundTemperature
    ! Temperature loop
    TemperatureLoop02A: DO
      ! Progress bar
      CALL Progress_Bar_Path_01( ThermodynamicProcess, Pressure, Temperature )
      ! Phase loop
      DO PhaseCounter = 1, 4
        ! Initialization
        FluidPhase = .FALSE.
        PhaseTest  = FluidPhase
        ! Phase type
        FluidPhase(PhaseCounter) = .TRUE.
        IF( PureComponent ) mFraction = 1.D0
        ! Find volume for the given temperature, pressure, and phase
        CALL Topliss_Algorithm( 1_Int64, mFraction, Temperature, Pressure, Volume, FluidPhase, CurveType, PureComponent )
        ! Cycle phase type if the selected phase is not the same as the calculated phase
        IF( ALL( PhaseTest ) .NEQV. ALL( FluidPhase ) ) THEN
          CYCLE
        END IF
      END DO
      ! Calculate isothermal compressibility [1 / Pa], thermal expansion coefficient [1 / K], and compressibility factor
      IF( PureComponent ) THEN
        CALL Calculate_Pressure_Single_Component( 1_Int64, Volume, Temperature, Anynumber, IsothermalCompressibility, &
        &                                         ThermalExpansionCoefficient, CompressibilityFactor )
      ELSE
        CALL Calculate_Pressure( mFraction, Volume, Temperature, Anynumber, IsothermalCompressibility, &
        &                        ThermalExpansionCoefficient, CompressibilityFactor )
      END IF
      ! Calculate entropy and specific heat at constant volume [J / (mol . K)]
      IF( PureComponent ) THEN
        CALL Calculate_Entropy_Single_Component( 1_Int64, Volume, Temperature, Entropy, vSpecificHeat )
      ELSE
        CALL Calculate_Entropy( mFraction, Volume, Temperature, Entropy, vSpecificHeat )
      END IF
      ! Calculate total chemical potential and residual chemical potential [Joule / mol]
      IF( PureComponent ) THEN
        CALL Calculate_CPotential_Single_Component( 1_Int64, Volume, Temperature, tChemicalPotential(1), rChemicalPotential(1) )
      ELSE
        CALL Calculate_CPotential( mFraction, Volume, Temperature, tChemicalPotential, rChemicalPotential )
      END IF
      ! Calculate molar density [mol / m³]
      Density = 1.D0 / Volume
      ! Calculate mass density [kg / m³]
      IF( PureComponent ) THEN
        MassDensity = Density * cMolarMass(1)
      ELSE
        MassDensity = SUM( mFraction * Density * cMolarMass )
      END IF
      ! Calculate specific heat at constant pressure using Mayer's relation [J / (mol . K)]
      pSpecificHeat = vSpecificHeat + Volume * Temperature * ThermalExpansionCoefficient * ThermalExpansionCoefficient / &
      &               IsothermalCompressibility
      ! Calculate the speed of sound using Newton-Laplace equation [m / s]
      SoundSpeed = DSQRT( ( pSpecificHeat / vSpecificHeat ) * ( 1.D0 / ( IsothermalCompressibility * MassDensity ) ) )
      ! Calculate the Joule-Thompson coefficient [K / Pa]
      JTCoefficient = ( Volume / pSpecificHeat ) * ( Temperature * ThermalExpansionCoefficient - 1.D0 )
      ! Write results
      WRITE( 20, "(13(G0,','))", Advance= "NO" ) Temperature, Pressure / 1.D6, Density, MassDensity, Volume, &
      &     IsothermalCompressibility * 1.D6, ThermalExpansionCoefficient, CompressibilityFactor, Entropy, vSpecificHeat, &
      &     pSpecificHeat, SoundSpeed, JTCoefficient
      DO iComponent = 1, nComponents - 1
        WRITE( 20, "(2(G0,','))", Advance= "NO" ) tChemicalPotential(iComponent), rChemicalPotential(iComponent)
      END DO
      WRITE( 20, "(G0,',',G0)", Advance= "YES" ) tChemicalPotential(nComponents), rChemicalPotential(nComponents)
      FLUSH( 20 )
      ! Break loop condition
      IF( UpperBoundTemperature == LowerBoundTemperature ) EXIT TemperatureLoop02A
      ! Increment temperature
      Temperature = Temperature + StepTemperature
      ! Temperature limit
      IF( .NOT. LeaveLoop02 ) THEN
        ! Break loop condition
        IF( Temperature >= UpperBoundTemperature ) THEN
          Temperature = UpperBoundTemperature
          LeaveLoop02 = .TRUE.
        END IF
      ELSE
        EXIT TemperatureLoop02A
      END IF
    END DO TemperatureLoop02A
    ! Break loop condition
    IF( UpperBoundPressure == LowerBoundPressure ) EXIT PressureLoop01
    ! Increment pressure
    Pressure = Pressure + StepPressure
    ! Pressure limit
    IF( .NOT. LeaveLoop01 ) THEN
      ! Break loop condition
      IF( Pressure >= UpperBoundPressure ) THEN
        Pressure = UpperBoundPressure
        LeaveLoop01 = .TRUE.
      END IF
    ELSE
      EXIT PressureLoop01
    END IF
  END DO PressureLoop01
END IF

! Isochoric process
IF( ThermodynamicProcess == 3 ) THEN
  ! Check parent subfolders
  INQUIRE( File= "ROUTE_01_THERMODYNAMIC_PROCESS/ISOCHORIC/", Exist= FolderExists )
  IF( .NOT. FolderExists ) THEN
    CALL SYSTEM( "mkdir ROUTE_01_THERMODYNAMIC_PROCESS/ISOCHORIC/" )
  END IF
  INQUIRE( File= "ROUTE_01_THERMODYNAMIC_PROCESS/ISOCHORIC/"//TRIM( DescriptorDate )//"/", Exist= FolderExists )
  IF( .NOT. FolderExists ) THEN
    CALL SYSTEM( "mkdir ROUTE_01_THERMODYNAMIC_PROCESS/ISOCHORIC/"//TRIM( DescriptorDate )//"/" )
  END IF
  ! File name
  IF( PureComponent ) THEN
    Filename = "ROUTE_01_THERMODYNAMIC_PROCESS/ISOCHORIC/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )//"_"//&
    &          TRIM( cMoleculeName(1) )//".dat"
  ELSE
    Filename = "ROUTE_01_THERMODYNAMIC_PROCESS/ISOCHORIC/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )//"_"//&
    &          TRIM( cMoleculeName(1) )
    DO iComponent = 2, nComponents - 1
      Filename = TRIM( Filename )//"_"//TRIM( cMoleculeName(iComponent) )
    END DO
    Filename = TRIM( Filename )//"_"//TRIM( cMoleculeName(nComponents) )//".dat"
  END IF
  ! Create output file
  OPEN( Unit= 20, File= Filename )
  WRITE( 20, "(G0)" ) "Isochoric Route"
  WRITE( 20, "(G0)" ) " "
  WRITE( 20, "(G0,G0.7,G0,G0.7,G0)" ) "Density Range: [", LowerBoundDensity, ", ", UpperBoundDensity, "] mol/m³"
  WRITE( 20, "(G0,G0.5,G0,G0.5,G0)" ) "Temperature Range: [", LowerBoundTemperature, ", ", UpperBoundTemperature, "] K"
  WRITE( 20, "(G0)" ) " "
  IF( nComponents == 1 ) THEN
    cComponent = 1_Int64
    WRITE( 20, "(G0)" ) "Molecule Formula: "//TRIM( cFormulaName(cComponent) )//""
    WRITE( *, "(G0)" ) "Molecule Formula: "//TRIM( cFormulaName(cComponent) )//""
    WRITE( 20, "(G0)" ) " "
    WRITE( *, "(G0)" ) " "
    WRITE( 20, "(G0)" ) "Geometry: "//TRIM( cGeometryName(cComponent) )//""
    WRITE( *, "(G0)" ) "Geometry: "//TRIM( cGeometryName(cComponent) )//""
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
  ELSE IF( nComponents > 1 ) THEN
    WRITE( 20, "(G0)", Advance= "NO" ) "Mixture: "
    WRITE( *, "(G0)", Advance= "NO" ) "Mixture: "
    DO cComponent = 1, nComponents - 1
      WRITE( 20, "(3G0)", Advance= "NO" ) ""//TRIM( cFormulaName(cComponent) )//"(", cComponent, ") + "
      WRITE( *, "(3G0)", Advance= "NO" ) ""//TRIM( cFormulaName(cComponent) )//"(", cComponent, ") + "
    END DO
    WRITE( 20, "(3G0)", Advance= "YES" ) ""//TRIM( cFormulaName(cComponent) )//"(", nComponents, ")"
    WRITE( *, "(3G0)", Advance= "YES" ) ""//TRIM( cFormulaName(cComponent) )//"(", nComponents, ")"
    WRITE( 20, "(G0)" ) " "
    WRITE( *, "(G0)" ) " "
    WRITE( 20, "(G0)", Advance= "NO" ) "Molar Fraction: "
    WRITE( *, "(G0)", Advance= "NO" ) "Molar Fraction: "
    DO cComponent = 1, nComponents - 1
      WRITE( 20, "(G0.5,3G0)", Advance= "NO" ) mFraction(cComponent), "(", cComponent, ") + "
      WRITE( *, "(G0.5,3G0)", Advance= "NO" ) mFraction(cComponent), "(", cComponent, ") + "
    END DO
    WRITE( 20, "(G0.5,3G0)", Advance= "YES" ) mFraction(cComponent), "(", nComponents, ")"
    WRITE( *, "(G0.5,3G0)", Advance= "YES" ) mFraction(cComponent), "(", nComponents, ")"
    WRITE( 20, "(G0)" ) " "
    WRITE( *, "(G0)" ) " "
    WRITE( 20, "(G0)", Advance= "NO" ) "Geometry: "
    WRITE( *, "(G0)", Advance= "NO" ) "Geometry: "
    DO cComponent = 1, nComponents - 1
      WRITE( 20, "(3G0)", Advance= "NO" ) ""//TRIM( cGeometryName(cComponent) )//"(", cComponent, ") + "
      WRITE( *, "(3G0)", Advance= "NO" ) ""//TRIM( cGeometryName(cComponent) )//"(", cComponent, ") + "
    END DO
    WRITE( 20, "(3G0)", Advance= "YES" ) ""//TRIM( cGeometryName(cComponent) )//"(", nComponents, ")"
    WRITE( *, "(3G0)", Advance= "YES" ) ""//TRIM( cGeometryName(cComponent) )//"(", nComponents, ")"
    WRITE( 20, "(G0)" ) " "
    WRITE( *, "(G0)" ) " "
    IF( ANY( cGeometry == 1 ) ) THEN
      WRITE( 20, "(G0)" ) "EOR model: Isihara-Hadwiger Theorem"
      WRITE( *, "(G0)" ) "EOR model: Isihara-Hadwiger Theorem"
      WRITE( 20, "(G0)" ) " "
      WRITE( *, "(G0)" ) " "
    END IF
    IF( ANY( cGeometry == 2 ) ) THEN
      WRITE( 20, "(G0)" ) "SPC model: Isihara-Hadwiger Theorem"
      WRITE( *, "(G0)" ) "SPC model: Isihara-Hadwiger Theorem"
      WRITE( 20, "(G0)" ) " "
      WRITE( *, "(G0)" ) " "
    END IF
    IF( ANY( cGeometry == 3 ) ) THEN
      WRITE( 20, "(G0)" ) "CYL model: Isihara-Hadwiger Theorem"
      WRITE( *, "(G0)" ) "CYL model: Isihara-Hadwiger Theorem"
      WRITE( 20, "(G0)" ) " "
      WRITE( *, "(G0)" ) " "
    END IF
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
    IF( NonSphericalMixingRule == 1 ) THEN
      WRITE( 20, "(G0)" ) "Non-Spherical Mixing Rule: Aspect Ratio Rule"
      WRITE( *, "(G0)" ) "Non-Spherical Mixing Rule: Aspect Ratio Rule"
    ELSE IF( NonSphericalMixingRule == 2 ) THEN
      WRITE( 20, "(G0)" ) "Non-Spherical Mixing Rule: Angle Average Rule"
      WRITE( *, "(G0)" ) "Non-Spherical Mixing Rule: Angle Average Rule"
    ELSE IF( NonSphericalMixingRule == 3 ) THEN
      WRITE( 20, "(G0)" ) "Non-Spherical Mixing Rule: Ratio of Second Virial Coefficients"
      WRITE( *, "(G0)" ) "Non-Spherical Mixing Rule: Ratio of Second Virial Coefficients"
    ELSE IF( NonSphericalMixingRule == 4 ) THEN
      WRITE( 20, "(G0)" ) "Non-Spherical Mixing Rule: Second Virial Coefficients"
      WRITE( *, "(G0)" ) "Non-Spherical Mixing Rule: Second Virial Coefficients"
    END IF
    WRITE( 20, "(G0)" ) " "
    WRITE( *, "(G0)" ) " "
    WRITE( 20, "(G0)" ) "Zhang Factor of the Mixture: Mixed Spherical Diameter"
    WRITE( *, "(G0)" ) "Zhang Factor of the Mixture: Mixed Spherical Diameter"
    WRITE( 20, "(G0)" ) " "
    WRITE( *, "(G0)" ) " "
    IF( EffPFractionMixingRule == 1 ) THEN
      WRITE( 20, "(G0)" ) "Effective Packing Fraction Mixing Rule: Overall Reduced Density Mixing Rule"
      WRITE( *, "(G0)" ) "Effective Packing Fraction Mixing Rule: Overall Reduced Density Mixing Rule"
    ELSE IF( EffPFractionMixingRule == 2 ) THEN
      WRITE( 20, "(G0)" ) "Effective Packing Fraction Mixing Rule: One-Fluid van der Waals Mixing Rule"
      WRITE( *, "(G0)" ) "Effective Packing Fraction Mixing Rule: One-Fluid van der Waals Mixing Rule"
    END IF
    WRITE( 20, "(G0)" ) " "
    WRITE( *, "(G0)" ) " "
  END IF
  WRITE( 20, "(13(G0,','))", Advance= "NO" ) "'Temperature [K]'", "'Pressure [MPa]'", "'Molar Density [mol/m³]'", &
  &     "'Mass Density [kg/m³]'", "'Molar Volume [m³/mol]'", "'Isothermal Compressibility [1/MPa]'", &
  &     "'Thermal Expansion Coefficient [1/K]'", "'Compressibility Factor'", "'Entropy [J/(mol.K)]'", &
  &     "'Specific Heat at Constant Volume [J/(mol.K)]'", "'Specific Heat at Constant Pressure [J/(mol.K)]'", &
  &     "'Speed of Sound [m/s]'", "'Joule-Thompson Coefficient [K/MPa]'"
  DO iComponent = 1, nComponents - 1
    WRITE( 20, "(2(G0,','))", Advance= "NO" ) "'Total Chemical Potential [J/mol] of "//TRIM( cMoleculeName(iComponent) )//"'", &
    &     "'Residual Chemical Potential [J/mol] of "//TRIM( cMoleculeName(iComponent) )//"'"
  END DO
  WRITE( 20, "(G0,',',G0)", Advance= "YES" ) "'Total Chemical Potential [J/mol] of "//TRIM( cMoleculeName(nComponents) )//"'", &
  &     "'Residual Chemical Potential [J/mol] of "//TRIM( cMoleculeName(nComponents) )//"'"
  FLUSH( 20 )
  ! Initialization
  LeaveLoop01 = .FALSE.
  Density = LowerBoundDensity
  ! Density loop
  DensityLoop01: DO
    ! Initialization
    LeaveLoop02 = .FALSE.
    Temperature = LowerBoundTemperature
    ! Temperature loop
    TemperatureLoop02B: DO
      ! Progress bar
      CALL Progress_Bar_Path_01( ThermodynamicProcess, Density, Temperature )
      ! Calculate molar volume [m³ / mol]
      Volume = 1.D0 / Density
      ! Calculate isothermal compressibility [1 / Pa], thermal expansion coefficient [1 / K], and compressibility factor
      IF( PureComponent ) THEN
        CALL Calculate_Pressure_Single_Component( 1_Int64, Volume, Temperature, Pressure, IsothermalCompressibility, &
        &                                         ThermalExpansionCoefficient, CompressibilityFactor )
      ELSE
        CALL Calculate_Pressure( mFraction, Volume, Temperature, Pressure, IsothermalCompressibility, &
        &                        ThermalExpansionCoefficient, CompressibilityFactor )
      END IF
      ! Phase loop
      DO PhaseCounter = 1, 4
        ! Initialization
        FluidPhase = .FALSE.
        PhaseTest  = FluidPhase
        ! Phase type
        FluidPhase(PhaseCounter) = .TRUE.
        IF( PureComponent ) mFraction = 1.D0
        ! Find phase and curve type for the given temperature and pressure (density already defined)
        CALL Topliss_Algorithm( 1_Int64, mFraction, Temperature, Pressure, Anynumber, FluidPhase, CurveType, PureComponent )
        ! Cycle phase type if the selected phase is not the same as the calculated phase
        IF( ALL( PhaseTest ) .NEQV. ALL( FluidPhase ) ) THEN
          CYCLE
        END IF
      END DO
      ! Calculate entropy and specific heat at constant volume [J / (mol . K)]
      IF( PureComponent ) THEN
        CALL Calculate_Entropy_Single_Component( 1_Int64, Volume, Temperature, Entropy, vSpecificHeat )
      ELSE
        CALL Calculate_Entropy( mFraction, Volume, Temperature, Entropy, vSpecificHeat )
      END IF
      ! Calculate total chemical potential and residual chemical potential [Joule / mol]
      IF( PureComponent ) THEN
        CALL Calculate_CPotential_Single_Component( 1_Int64, Volume, Temperature, tChemicalPotential(1), rChemicalPotential(1) )
      ELSE
        CALL Calculate_CPotential( mFraction, Volume, Temperature, tChemicalPotential, rChemicalPotential )
      END IF
      ! Calculate mass density [kg / m³]
      IF( PureComponent ) THEN
        MassDensity = Density * cMolarMass(1)
      ELSE
        MassDensity = SUM( mFraction * Density * cMolarMass )
      END IF
      ! Calculate specific heat at constant pressure using Mayer's relation [J / (mol . K)]
      pSpecificHeat = vSpecificHeat + Volume * Temperature * ThermalExpansionCoefficient * ThermalExpansionCoefficient / &
      &               IsothermalCompressibility
      ! Calculate the speed of sound using Newton-Laplace equation [m / s]
      SoundSpeed = DSQRT( ( pSpecificHeat / vSpecificHeat ) * ( 1.D0 / ( IsothermalCompressibility * MassDensity ) ) )
      ! Calculate the Joule-Thompson coefficient [K / Pa]
      JTCoefficient = ( Volume / pSpecificHeat ) * ( Temperature * ThermalExpansionCoefficient - 1.D0 )
      ! Write results
      WRITE( 20, "(13(G0,','))", Advance= "NO" ) Temperature, Pressure / 1.D6, Density, MassDensity, Volume, &
      &     IsothermalCompressibility * 1.D6, ThermalExpansionCoefficient, CompressibilityFactor, Entropy, vSpecificHeat, &
      &     pSpecificHeat, SoundSpeed, JTCoefficient
      DO iComponent = 1, nComponents - 1
        WRITE( 20, "(2(G0,','))", Advance= "NO" ) tChemicalPotential(iComponent), rChemicalPotential(iComponent)
      END DO
      WRITE( 20, "(G0,',',G0)", Advance= "YES" ) tChemicalPotential(nComponents), rChemicalPotential(nComponents)
      FLUSH( 20 )
      ! Break loop condition
      IF( UpperBoundTemperature == LowerBoundTemperature ) EXIT TemperatureLoop02B
      ! Increment temperature
      Temperature = Temperature + StepTemperature
      ! Temperature limit
      IF( .NOT. LeaveLoop02 ) THEN
        ! Break loop condition
        IF( Temperature >= UpperBoundTemperature ) THEN
          Temperature = UpperBoundTemperature
          LeaveLoop02 = .TRUE.
        END IF
      ELSE
        EXIT TemperatureLoop02B
      END IF
    END DO TemperatureLoop02B
    ! Break loop condition
    IF( UpperBoundDensity == LowerBoundDensity ) EXIT DensityLoop01
    ! Increment density
    Density = Density + StepDensity
    ! Density limit
    IF( .NOT. LeaveLoop01 ) THEN
      ! Break loop condition
      IF( Density >= UpperBoundDensity ) THEN
        Density = UpperBoundDensity
        LeaveLoop01 = .TRUE.
      END IF
    ELSE
      EXIT DensityLoop01
    END IF
  END DO DensityLoop01
END IF

! Status
WRITE( *, "(G0)" ) " "

RETURN

END SUBROUTINE Thermodynamic_Properties