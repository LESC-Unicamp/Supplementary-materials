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
SUBROUTINE Simplex_Optimization( PureComponents )

! Uses two modules: global variables and substances
USE GlobalVariables
USE Substances

IMPLICIT NONE

! ************************************************************************************************ !
! INTEGER VARIABLES                                                                                !
! ************************************************************************************************ !
INTEGER( Kind= Int64 ) :: cComponent      ! Counter (components)
INTEGER( Kind= Int64 ) :: iPoint, nPoints ! Number of points in file
INTEGER( Kind= Int64 ) :: nTotalSkips     ! Total number of skips due to temperature bounds
INTEGER( Kind= Int64 ) :: nSkipsBelow     ! Number of skips due to upper temperature bound
INTEGER( Kind= Int64 ) :: nSkipsAbove     ! Number of skips due to lower temperature bound
INTEGER( Kind= Int64 ) :: nParameters     ! Number of optimization parameters
INTEGER( Kind= Int64 ) :: kLine           ! Counter (lines)
INTEGER( Kind= Int64 ) :: cPos            ! Counter (array position)
INTEGER( Kind= Int64 ) :: cSkips          ! Counter (skips)
INTEGER( Kind= Int64 ) :: nLineSkips      ! Number of line skips
INTEGER( Kind= Int64 ) :: Counter         ! Counter (generic)
INTEGER( Kind= Int64 ) :: dCounter        ! Counter (damping)
INTEGER( Kind= Int64 ) :: cGlobalSimplex  ! Number of global simplex cycles

! ************************************************************************************************ !
! REAL PARAMETERS                                                                                  !
! ************************************************************************************************ !
REAL( Kind= Real64 ), PARAMETER :: Tolerance = 1.D-9 ! Tolerance

! ************************************************************************************************ !
! REAL VARIABLES                                                                                   !
! ************************************************************************************************ !
REAL( Kind= Real64 )                           :: Temperature              ! Temperature
REAL( Kind= Real64 )                           :: TemperatureMixture       ! Temperature of the mixture
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
REAL( Kind= Real64 )                           :: SkipAboveTripleTemp      ! Skip above triple point (Kelvin)
REAL( Kind= Real64 )                           :: SkipBelowCritTemp        ! Skip below critical point (Kelvin)
REAL( Kind= Real64 )                           :: TriplePointTemperature   ! Temperature of the triple point
REAL( Kind= Real64 )                           :: CriticalPointTemperature ! Temperature of the critical point
REAL( Kind= Real64 )                           :: TemperatureAbove         ! Temperature (above triple point)
REAL( Kind= Real64 )                           :: TemperatureBelow         ! Temperature (below critical point)
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: cAccentricFactor         ! Accentric factor
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: cCriticalPressure        ! Critical pressure
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: cCriticalTemperature     ! Critical temperature
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: cCriticalDensity         ! Critical density
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: cCriticalVolume          ! Critical volume
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: cVaporPressure           ! Vapor pressure
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: mFraction                ! Mole fraction

! ************************************************************************************************ !
! REAL VARIABLES (ALLOCATABLE)                                                                     !
! ************************************************************************************************ !
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE :: EqLiquidComposition       ! Liquid composition (equilibrium)
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE :: EqVaporComposition        ! Vapor composition (equilibrium)
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE :: EqPressure                ! Pressure (equilibrium)
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE :: EqSaturatedVaporDens      ! Saturated vapor density (equilibrium)
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE :: EqSaturatedLiquidDens     ! Saturated liquid density (equilibrium)
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE :: EqVaporPressure           ! Vapor pressure (equilibrium)
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE :: EqTemperature             ! Temperature (equilibrium)
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE :: EqSaturatedVaporDensTemp  ! Saturated vapor density (equilibrium)
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE :: EqSaturatedLiquidDensTemp ! Saturated liquid density (equilibrium)
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE :: EqVaporPressureTemp       ! Vapor pressure (equilibrium)
REAL( Kind= Real64 ), DIMENSION( : ), ALLOCATABLE :: EqTemperatureTemp         ! Temperature (equilibrium)

! ************************************************************************************************ !
! CHARACTER STRINGS                                                                                !
! ************************************************************************************************ !
CHARACTER( LEN= 100 ) :: FileName       ! File name
CHARACTER( LEN= 15 )  :: TemperatureStr ! Temperature string
CHARACTER( Len= 01 )  :: CurveType      ! Isotherm types (A, B, or C)
CHARACTER( LEN= 01 )  :: Dummy          ! Dummy string

! ************************************************************************************************ !
! LOGICAL VARIABLES                                                                                !
! ************************************************************************************************ !
LOGICAL                 :: PureComponents ! Flag to indicate if the components are pure
LOGICAL                 :: FileExists     ! Flag to indicate if the file exists
LOGICAL, DIMENSION( 4 ) :: FluidPhase     ! Phase type (1: Supercritical fluid; 2: Critical fluid; 3: Subcritical liquid phase; 4: Subcritical vapor phase)

! Pure components
IF( PureComponents ) THEN
  ! File name
  FileName = "../src/VLE_DATA/PURE/"//TRIM( cFormulaName(1) )//"/data.csv"
  ! Inquire if file exists
  INQUIRE( File= FileName, Exist= FileExists )
  IF( .NOT. FileExists ) THEN
    WRITE( *, "(G0)" ) "File not found! Exiting..."
    CALL EXIT(  )
  END IF
  ! Read the data from the file
  OPEN( Unit= 35, File= FileName, Action= "READ" )
  READ( 35, * ) Dummy, Dummy, Dummy, Dummy
  READ( 35, * ) nPoints
  IF( nPoints <= 1 ) THEN
    WRITE( *, "(G0)" ) "Error: the number of points must be a positive integer greater than 1. Exiting..."
    CALL EXIT(  )
  END IF
  CLOSE( 35 )
  ! Read parameters
  OPEN( Unit= 35, File= "process_3_specs.ini", Action= "READ" )
  READ( 35, * ) Dummy, Dummy
  READ( 35, * ) Dummy, DampingFactorIteration
  IF( DampingFactorIteration < 0.D0 .OR. DampingFactorIteration > 1.D0 ) THEN
    WRITE( *, "(G0,G0.5,G0)" ) "Error: the damping factor iteration is out of range [0,1]. Exiting..."
    CALL EXIT(  )
  END IF
  READ( 35, * ) Dummy, Dummy
  READ( 35, * ) Dummy, nLineSkips
  IF( nLineSkips < 0 ) THEN
    WRITE( *, "(G0)" ) "Error: the number of lines must be a positive integer. Exiting..."
    CALL EXIT(  )
  END IF
  READ( 35, * ) Dummy, SkipAboveTripleTemp
  IF( SkipAboveTripleTemp < 0.D0 ) THEN
    WRITE( *, "(G0)" ) "Error: the temperature skip above the triple point must be a positive value. Exiting..."
    CALL EXIT(  )
  END IF
  READ( 35, * ) Dummy, SkipBelowCritTemp
  IF( SkipBelowCritTemp < 0.D0 ) THEN
    WRITE( *, "(G0)" ) "Error: the temperature skip below the critical point must be a positive value. Exiting..."
    CALL EXIT(  )
  END IF
  READ( 35, * ) Dummy, cGlobalSimplex
  IF( cGlobalSimplex <= 0 ) THEN
    WRITE( *, "(G0)" ) "Error: the number of global simplex cycles must be a non-zero positive integer. Exiting..."
    CALL EXIT(  )
  END IF
  CLOSE( 35 )
  ! Find triple and critical temperature points
  OPEN( Unit= 35, File= FileName, Action= "READ" )
  READ( 35, * ) Dummy, Dummy, Dummy, Dummy
  READ( 35, * ) Dummy
  READ( 35, * ) TriplePointTemperature
  DO iPoint = 2, nPoints - 1
    READ( 35, * ) Dummy
  END DO
  READ( 35, * ) CriticalPointTemperature
  CLOSE( 35 )
  ! Temperature limits
  TemperatureAbove = TriplePointTemperature + SkipAboveTripleTemp
  TemperatureBelow = CriticalPointTemperature - SkipBelowCritTemp
  IF( ( TemperatureBelow - TemperatureAbove ) < 0.D0 ) THEN
    WRITE( *, "(2(G0,G0.5),G0)" ) "Error: the lower bound [", TemperatureAbove, "K] is greater than the upper bound [", &
    &                             TemperatureBelow, "K]. Exiting..."
    CALL EXIT(  )
  END IF
  ! Initialize
  nSkipsAbove = 0
  nSkipsBelow = 0
  ! Find number of skips due to temperature bounds
  OPEN( Unit= 35, File= FileName, Action= "READ" )
  READ( 35, * ) Dummy, Dummy, Dummy, Dummy
  READ( 35, * ) Dummy
  DO iPoint = 1, nPoints
    READ( 35, * ) Temperature
    IF( Temperature < TemperatureAbove ) nSkipsAbove = nSkipsAbove + 1
    IF( Temperature > TemperatureBelow ) nSkipsBelow = nSkipsBelow + 1
    nTotalSkips = nSkipsAbove + nSkipsBelow
  END DO
  CLOSE( 35 )
  ! Correct number of points based on skips
  nPoints = nPoints - nTotalSkips
  ! Correct number of skips
  IF( nLineSkips >= (nPoints - 1) ) nLineSkips = nPoints - 2
  ! Deallocation of arrays
  IF( ALLOCATED( EqSaturatedLiquidDens ) ) DEALLOCATE( EqSaturatedLiquidDens )
  IF( ALLOCATED( EqSaturatedVaporDens ) ) DEALLOCATE( EqSaturatedVaporDens )
  IF( ALLOCATED( EqVaporPressure ) ) DEALLOCATE( EqVaporPressure )
  IF( ALLOCATED( EqTemperature ) ) DEALLOCATE( EqTemperature )
  ! Initial allocation of arrays
  ALLOCATE( EqSaturatedLiquidDens(1), EqSaturatedVaporDens(1), EqVaporPressure(1), EqTemperature(1) )
  ! Experimental data (saturated data points)
  OPEN( Unit= 35, File= FileName, Action= "READ" )
  READ( 35, * ) Dummy, Dummy, Dummy, Dummy
  READ( 35, * ) Dummy
  ! Skip points above the triple point
  DO iPoint = 1, nSkipsAbove
    READ( 35, * ) Dummy, Dummy, Dummy, Dummy
  END DO
  ! First real line
  READ( 35, * ) EqTemperature(1), EqVaporPressure(1), EqSaturatedVaporDens(1), EqSaturatedLiquidDens(1)
  ! Initialization
  cSkips = 0
  kLine  = 2
  cPos   = 2
  ! Dynamic allocation of arrays
  DO WHILE( kLine < nPoints )
    ! Line skips
    IF( cSkips < nLineSkips ) THEN
      ! Skip line
      READ( 35, * ) Dummy, Dummy, Dummy, Dummy
      ! Increment counter
      cSkips = cSkips + 1
    ELSE
      ! Deallocation of temporary arrays
      IF( ALLOCATED( EqSaturatedLiquidDensTemp ) ) DEALLOCATE( EqSaturatedLiquidDensTemp )
      IF( ALLOCATED( EqSaturatedVaporDensTemp ) ) DEALLOCATE( EqSaturatedVaporDensTemp )
      IF( ALLOCATED( EqVaporPressureTemp ) ) DEALLOCATE( EqVaporPressureTemp )
      IF( ALLOCATED( EqTemperatureTemp ) ) DEALLOCATE( EqTemperatureTemp )
      ! Re-allocation of temporary arrays
      ALLOCATE( EqSaturatedLiquidDensTemp(cPos), EqSaturatedVaporDensTemp(cPos), EqVaporPressureTemp(cPos), &
      &         EqTemperatureTemp(cPos) )
      ! Transfer data
      EqSaturatedLiquidDensTemp(1:(cPos-1)) = EqSaturatedLiquidDens(1:(cPos-1))
      EqSaturatedVaporDensTemp(1:(cPos-1))  = EqSaturatedVaporDens(1:(cPos-1))
      EqVaporPressureTemp(1:(cPos-1))       = EqVaporPressure(1:(cPos-1))
      EqTemperatureTemp(1:(cPos-1))         = EqTemperature(1:(cPos-1))
      ! Read data
      READ( 35, * ) EqTemperatureTemp(cPos), EqVaporPressureTemp(cPos), EqSaturatedVaporDensTemp(cPos), &
      &             EqSaturatedLiquidDensTemp(cPos)
      ! Deallocation of arrays
      IF( ALLOCATED( EqSaturatedLiquidDens ) ) DEALLOCATE( EqSaturatedLiquidDens )
      IF( ALLOCATED( EqSaturatedVaporDens ) ) DEALLOCATE( EqSaturatedVaporDens )
      IF( ALLOCATED( EqVaporPressure ) ) DEALLOCATE( EqVaporPressure )
      IF( ALLOCATED( EqTemperature ) ) DEALLOCATE( EqTemperature )
      ! Re-allocation of arrays
      ALLOCATE( EqSaturatedLiquidDens(cPos), EqSaturatedVaporDens(cPos), EqVaporPressure(cPos), EqTemperature(cPos) )
      ! Transfer data
      EqSaturatedLiquidDens(:) = EqSaturatedLiquidDensTemp(:)
      EqSaturatedVaporDens(:)  = EqSaturatedVaporDensTemp(:)
      EqVaporPressure(:)       = EqVaporPressureTemp(:)
      EqTemperature(:)         = EqTemperatureTemp(:)
      ! Increment counter
      cPos = cPos + 1
      ! Reset
      cSkips = 0
    END IF
    ! Increment counter
    kLine = kLine + 1
  END DO
  ! Deallocation of temporary arrays
  IF( ALLOCATED( EqSaturatedLiquidDensTemp ) ) DEALLOCATE( EqSaturatedLiquidDensTemp )
  IF( ALLOCATED( EqSaturatedVaporDensTemp ) ) DEALLOCATE( EqSaturatedVaporDensTemp )
  IF( ALLOCATED( EqVaporPressureTemp ) ) DEALLOCATE( EqVaporPressureTemp )
  IF( ALLOCATED( EqTemperatureTemp ) ) DEALLOCATE( EqTemperatureTemp )
  ! Re-allocation of temporary arrays
  ALLOCATE( EqSaturatedLiquidDensTemp(cPos), EqSaturatedVaporDensTemp(cPos), EqVaporPressureTemp(cPos), &
  &         EqTemperatureTemp(cPos) )
  ! Transfer data
  EqSaturatedLiquidDensTemp(1:(cPos-1)) = EqSaturatedLiquidDens(1:(cPos-1))
  EqSaturatedVaporDensTemp(1:(cPos-1))  = EqSaturatedVaporDens(1:(cPos-1))
  EqVaporPressureTemp(1:(cPos-1))       = EqVaporPressure(1:(cPos-1))
  EqTemperatureTemp(1:(cPos-1))         = EqTemperature(1:(cPos-1))
  ! Read data
  READ( 35, * ) EqTemperatureTemp(cPos), EqVaporPressureTemp(cPos), EqSaturatedVaporDensTemp(cPos), &
  &             EqSaturatedLiquidDensTemp(cPos)
  CLOSE( 35 )
  ! Deallocation of arrays
  IF( ALLOCATED( EqSaturatedLiquidDens ) ) DEALLOCATE( EqSaturatedLiquidDens )
  IF( ALLOCATED( EqSaturatedVaporDens ) ) DEALLOCATE( EqSaturatedVaporDens )
  IF( ALLOCATED( EqVaporPressure ) ) DEALLOCATE( EqVaporPressure )
  IF( ALLOCATED( EqTemperature ) ) DEALLOCATE( EqTemperature )
  ! Re-allocation of arrays
  ALLOCATE( EqSaturatedLiquidDens(cPos), EqSaturatedVaporDens(cPos), EqVaporPressure(cPos), EqTemperature(cPos) )
  ! Transfer data of last real line
  EqSaturatedLiquidDens(:) = EqSaturatedLiquidDensTemp(:)
  EqSaturatedVaporDens(:)  = EqSaturatedVaporDensTemp(:)
  EqVaporPressure(:)       = EqVaporPressureTemp(:)
  EqTemperature(:)         = EqTemperatureTemp(:)
  ! Status
  WRITE( *, "(2G0)" ) "Total number of points: ", cPos
  WRITE( *, "(2G0)" ) "Total number of skipped points: ", nPoints + nTotalSkips - cPos
  WRITE( *, "(G0)" ) " "
  WRITE( *, "(G0,G0.5,G0)" ) "Triple point temperature: ", TriplePointTemperature, "K"
  WRITE( *, "(G0,G0.5,G0)" ) "Critical point temperature: ", CriticalPointTemperature, "K"
  WRITE( *, "(G0)" ) " "
  WRITE( *, "(G0,G0.5,G0)" ) "Lower bound temperature: ", EqTemperature(1), "K"
  WRITE( *, "(G0,G0.5,G0)" ) "Upper bound temperature: ", EqTemperature(cPos), "K"
  WRITE( *, "(G0)" ) " "
  ! Simplex algorithm
  nParameters = 4
  CALL Simplex_Algorithm_Pure( nParameters, cPos, cGlobalSimplex, DampingFactorIteration, EqSaturatedLiquidDens, EqVaporPressure, &
  &                            EqTemperature )
! Mixture
ELSE
  ! Initialization
  OPEN( Unit= 35, File= "process_3_specs.ini", Action= "READ" )
  READ( 35, * ) Dummy, TemperatureMixture
  IF( TemperatureMixture < 0.D0 ) THEN
    WRITE( *, "(G0,G0.5,G0)" ) "Error: the temperature of the mixture is negative. Exiting..."
    CALL EXIT(  )
  END IF
  READ( 35, * ) Dummy, DampingFactorIteration
  IF( DampingFactorIteration < 0.D0 .OR. DampingFactorIteration > 1.D0 ) THEN
    WRITE( *, "(G0,G0.5,G0)" ) "Error: the damping factor iteration is out of range [0,1]. Exiting..."
    CALL EXIT(  )
  END IF
  CLOSE( 35 )
  ! Critical properties
  DO cComponent = 1, nComponents
    CALL Critical_Properties_Pure_Components( cComponent, cCriticalTemperature(cComponent), cCriticalDensity(cComponent), &
    &                                         cCriticalPressure(cComponent), cCriticalVolume(cComponent) )
  END DO
  ! Accentric factor
  DO cComponent = 1, nComponents
    ! Initialization
    DampingFactor = 1.D0
    Temperature   = 0.7D0 * cCriticalTemperature(cComponent)
    mFraction     = 1.D0
    Counter       = 0
    dCounter      = 0
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
    CALL Topliss_Algorithm( cComponent, mFraction, Temperature, MidpointPressure, LiquidVolume, FluidPhase, CurveType, .TRUE. )
    ! Find chemical potential of the liquid phase
    CALL Calculate_CPotential_Single_Component( cComponent, LiquidVolume, Temperature, LiquidChemicalPotential, &
    &                                           rLiquidChemicalPotential )
    ! Vapor Phase
    FluidPhase(:) = .FALSE.
    FluidPhase(4) = .TRUE.
    ! Find volume of the vapor phase
    CALL Topliss_Algorithm( cComponent, mFraction, Temperature, MidpointPressure, VaporVolume, FluidPhase, CurveType, .TRUE. )
    ! Find chemical potential of the liquid phase
    CALL Calculate_CPotential_Single_Component( cComponent, VaporVolume, Temperature, VaporChemicalPotential, &
    &                                           rVaporChemicalPotential )
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
      CALL Topliss_Algorithm( cComponent, mFraction, Temperature, MidpointPressure, LiquidVolume, FluidPhase, CurveType, .TRUE. )
      ! Find chemical potential of the liquid phase
      CALL Calculate_CPotential_Single_Component( cComponent, LiquidVolume, Temperature, LiquidChemicalPotential, &
      &                                           rLiquidChemicalPotential )
      ! Vapor Phase
      FluidPhase(:) = .FALSE.
      FluidPhase(4) = .TRUE.
      ! Find volume of the vapor phase
      CALL Topliss_Algorithm( cComponent, mFraction, Temperature, MidpointPressure, VaporVolume, FluidPhase, CurveType, .TRUE. )
      ! Find chemical potential of the liquid phase
      CALL Calculate_CPotential_Single_Component( cComponent, VaporVolume, Temperature, VaporChemicalPotential, &
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
        WRITE( *, "(G0,G0.5,G0)" ) "Error: the damping factor is becoming too small. Convergence will probably fail. Exiting..."
        CALL EXIT(  )
      END IF
      ! Check iterations
      IF( Counter >= 20000 ) THEN
        WRITE( *, "(G0)" ) " "
        WRITE( *, "(G0)" ) " "
        WRITE( *, "(G0,G0.5,G0)" ) "Error: the number of iterations is too high. Convergence will probably fail. Exiting..."
        CALL EXIT(  )
      END IF
    END DO
    ! Vapor pressure
    cVaporPressure(cComponent) = MidpointPressure
    ! Accentric factor
    cAccentricFactor(cComponent) = - 1.D0 - DLOG10( cVaporPressure(cComponent) / cCriticalPressure(cComponent) )
  END DO
  ! Temperature string
  WRITE( TemperatureStr, "(F0.5)" ) TemperatureMixture
  ! File name
  FileName = "../src/VLE_DATA/MIXTURES/"//TRIM( cFormulaName(1) )//"+"//TRIM( cFormulaName(2) )//"_"//TRIM( TemperatureStr )//&
  &          "K/data.csv"
  ! Inquire if file exists
  INQUIRE( File= FileName, Exist= FileExists )
  IF( .NOT. FileExists ) THEN
    FileName = "../src/VLE_DATA/MIXTURES/"//TRIM( cFormulaName(2) )//"+"//TRIM( cFormulaName(1) )//"_"//TRIM( TemperatureStr )//&
    &          "K/data.csv"
    INQUIRE( File= FileName, Exist= FileExists )
    IF( .NOT. FileExists ) THEN
      print *,FileName
      WRITE( *, "(G0)" ) "File not found! Exiting..."
      CALL EXIT(  )
    END IF
  END IF
  ! Read the data from the file
  OPEN( Unit= 35, File= FileName, Action= "READ" )
  READ( 35, * )
  READ( 35, * ) nPoints
  IF( ALLOCATED( EqLiquidComposition ) ) DEALLOCATE( EqLiquidComposition )
  IF( ALLOCATED( EqVaporComposition ) ) DEALLOCATE( EqVaporComposition )
  IF( ALLOCATED( EqPressure ) ) DEALLOCATE( EqPressure )
  ALLOCATE( EqLiquidComposition( nPoints ), EqVaporComposition( nPoints ), EqPressure( nPoints ) )
  DO iPoint = 1, nPoints
    READ( 35, * ) EqLiquidComposition(iPoint), EqVaporComposition(iPoint), EqPressure(iPoint)
  END DO
  CLOSE( 35 )
  ! Read number of optmiization parameters
  OPEN( Unit= 35, File= "process_3_specs.ini", Action= "READ" )
  READ( 35, * )
  READ( 35, * )
  READ( 35, * ) Dummy, nParameters
  IF( nParameters <= 0 .OR. nParameters > 2 ) THEN
    WRITE( *, "(G0)" ) "Error: the number of optimization parameters must be 1 or 2. Exiting..."
    CALL EXIT(  )
  END IF
  READ( 35, * )
  READ( 35, * )
  READ( 35, * )
  READ( 35, * ) Dummy, cGlobalSimplex
  IF( cGlobalSimplex <= 0 ) THEN
    WRITE( *, "(G0)" ) "Error: the number of global simplex cycles must be a non-zero positive integer. Exiting..."
    CALL EXIT(  )
  END IF
  CLOSE( 35 )
  ! Simplex algorithm
  CALL Simplex_Algorithm_Mixture( nParameters, nPoints, cGlobalSimplex, TemperatureMixture, cAccentricFactor, cCriticalPressure, &
  &                               cCriticalTemperature, EqLiquidComposition, EqVaporComposition, EqPressure )
END IF

RETURN

END SUBROUTINE Simplex_Optimization

! ************************************************************************************************ !
!                          Simplex algorithm to search zero of functions                           !
!                Reference: Nelder, J. A.; Mead, R. A., Comp. J., 7, 308-313, 1965                 !
! ************************************************************************************************ !
SUBROUTINE Simplex_Algorithm_Pure( nParameters, nPoints, cGlobalSimplex, DampingFactorIteration, EqLiquidDens, EqVaporPressure, &
&                                  EqTemperature)

! Uses two modules: global variables and substances
USE GlobalVariables
USE Substances

IMPLICIT NONE

! ************************************************************************************************ !
! INTEGER VARIABLES                                                                                !
! ************************************************************************************************ !
INTEGER( Kind= Int64 ) :: iGuess         ! Counter (guesses)
INTEGER( Kind= Int64 ) :: iSimplex       ! Counter (simplex)
INTEGER( Kind= Int64 ) :: iParameter     ! Counter (parameters)
INTEGER( Kind= Int64 ) :: jParameter     ! Counter (parameters)
INTEGER( Kind= Int64 ) :: cComponent     ! Counter (components)
INTEGER( Kind= Int64 ) :: nPoints        ! Number of points
INTEGER( Kind= Int64 ) :: nParameters    ! Number of optimization parameters
INTEGER( Kind= Int64 ) :: CaseSelection  ! Case selection (1: Reflection or Expansion; 2: Contraction or Reduction/Shrink)
INTEGER( Kind= Int64 ) :: SimplexType    ! Simplex type (-1: Initial guesses; 0: Midpoint; 1: Reflection; 2: Expansion; 3: Contraction; 4: Reduction)
INTEGER( Kind= Int64 ) :: cGlobalSimplex ! Number of global simplex cycles

! ************************************************************************************************ !
! REAL PARAMETERS                                                                                  !
! ************************************************************************************************ !
REAL( Kind= Real64 ), PARAMETER :: Tolerance = 1.D-8 ! Tolerance
REAL( Kind= Real64 ), PARAMETER :: Alpha     = 1.D0  ! Reflection coefficient
REAL( Kind= Real64 ), PARAMETER :: Beta      = 0.5D0 ! Contraction coefficient
REAL( Kind= Real64 ), PARAMETER :: Gamma     = 2.D0  ! Expansion coefficient
REAL( Kind= Real64 ), PARAMETER :: MinDelta  = 1.D-3 ! Minimum variation
REAL( Kind= Real64 ), PARAMETER :: MaxDelta  = 1.D-2 ! Maximum variation

! ************************************************************************************************ !
! REAL VARIABLES                                                                                   !
! ************************************************************************************************ !
REAL( Kind= Real64 )                                          :: RandomNumber                  ! Random number
REAL( Kind= Real64 )                                          :: Delta                         ! Variation
REAL( Kind= Real64 )                                          :: DampingFactorIteration        ! Damping factor (iteration)
REAL( Kind= Real64 )                                          :: ObjectiveFunctionMidpoint     ! Objective function (midpoint)
REAL( Kind= Real64 )                                          :: ObjectiveFunctionReflection   ! Objective function (reflection)
REAL( Kind= Real64 )                                          :: ObjectiveFunctionExpansion    ! Objective function (expansion)
REAL( Kind= Real64 )                                          :: ObjectiveFunctionContraction  ! Objective function (contraction)
REAL( Kind= Real64 )                                          :: ObjectiveFunctionReduction    ! Objective function (reduction)
REAL( Kind= Real64 )                                          :: ObjectiveFunctionAverage      ! Objective function (average)
REAL( Kind= Real64 )                                          :: ConvergenceCriterion          ! Convergence criterion
REAL( Kind= Real64 )                                          :: SecondLastObjectiveFunction   ! Second last value of the objective function
REAL( Kind= Real64 )                                          :: AuxSum                        ! Auxiliary sum
REAL( Kind= Real64 )                                          :: Dummy                         ! Dummy variable
REAL( Kind= Real64 ), DIMENSION( nPoints )                    :: EqLiquidDens                  ! Liquid composition (equilibrium)
REAL( Kind= Real64 ), DIMENSION( nPoints )                    :: EqVaporPressure               ! Pressure (equilibrium)
REAL( Kind= Real64 ), DIMENSION( nPoints )                    :: EqTemperature                 ! Temperature (equilibrium)
REAL( Kind= Real64 ), DIMENSION( nComponents )                :: cCriticalVolume               ! Molar volume (critical)
REAL( Kind= Real64 ), DIMENSION( nComponents )                :: cCriticalDensity              ! Molar density (critical)
REAL( Kind= Real64 ), DIMENSION( nComponents )                :: cCriticalPressure             ! Pressure (critical)
REAL( Kind= Real64 ), DIMENSION( nComponents )                :: cCriticalTemperature          ! Temperature (critical)
REAL( Kind= Real64 ), DIMENSION( nParameters )                :: LowerBoundOptimization        ! Lower bound
REAL( Kind= Real64 ), DIMENSION( nParameters )                :: UpperBoundOptimization        ! Upper bound
REAL( Kind= Real64 ), DIMENSION( nParameters )                :: gLowerBoundOptimization       ! Upper bound (global)
REAL( Kind= Real64 ), DIMENSION( nParameters )                :: gUpperBoundOptimization       ! Upper bound (global)
REAL( Kind= Real64 ), DIMENSION( nParameters )                :: MidpointSet                   ! Set of parameters (midpoint)
REAL( Kind= Real64 ), DIMENSION( nParameters )                :: ReflectionSet                 ! Set of parameters (reflection)
REAL( Kind= Real64 ), DIMENSION( nParameters )                :: ExpansionSet                  ! Set of parameters (expansion)
REAL( Kind= Real64 ), DIMENSION( nParameters )                :: ContractionSet                ! Set of parameters (contraction)
REAL( Kind= Real64 ), DIMENSION( nParameters )                :: ReductionSet                  ! Set of parameters (reduction)
REAL( Kind= Real64 ), DIMENSION( nParameters+1 )              :: ObjectiveFunction             ! Objective function
REAL( Kind= Real64 ), DIMENSION( nParameters+1, nParameters ) :: ParameterSets                 ! Parameter sets
REAL( Kind= Real64 ), DIMENSION( nParameters+1, nParameters ) :: Guesses                       ! Guesses

! ************************************************************************************************ !
! CHARACTER STRINGS                                                                                !
! ************************************************************************************************ !
CHARACTER( Len= 100 ) :: FileName ! File name

! ************************************************************************************************ !
! LOGICAL VARIABLES                                                                                !
! ************************************************************************************************ !
LOGICAL :: FolderExists ! Flag to indicate if the folder exists

! Initialization
cComponent = 1
nSimplexCycles = 1
LastObjectiveFunction = -1.D0
SecondLastObjectiveFunction = LastObjectiveFunction

! Initial guess for the optimization parameters (original values)
Guesses(1,1) = cWellDepth(cComponent)      ! Well depth (J)
Guesses(1,2) = cPotentialRange(cComponent) ! Attractive range
Guesses(1,3) = cDiameter(cComponent)       ! Minor axis diameter (Å)
Guesses(1,4) = cAspectRatio(cComponent)    ! Aspect ratio

! Global bounds
IF( PotentialTypeLogical(1) ) gUpperBoundOptimization(2) = 1.8D0  ! Square-well potential
IF( PotentialTypeLogical(2) ) gUpperBoundOptimization(2) = 12.0D0 ! Sutherland potential
IF( PotentialTypeLogical(3) ) gUpperBoundOptimization(2) = 4.0D0  ! Yukawa potential
IF( PotentialTypeLogical(4) ) gUpperBoundOptimization(2) = 0.8D0  ! Convex square-well potential
IF( PotentialTypeLogical(1) ) gLowerBoundOptimization(2) = 1.1D0  ! Square-well potential
IF( PotentialTypeLogical(2) ) gLowerBoundOptimization(2) = 3.0D0  ! Sutherland potential
IF( PotentialTypeLogical(3) ) gLowerBoundOptimization(2) = 1.1D0  ! Yukawa potential
IF( PotentialTypeLogical(4) ) gLowerBoundOptimization(2) = 0.1D0  ! Convex square-well potential
IF( cGeometry(cComponent) == 1 ) gUpperBoundOptimization(4) = 6.D0
IF( cGeometry(cComponent) == 2 ) gUpperBoundOptimization(4) = 6.D0
IF( cGeometry(cComponent) == 3 ) gUpperBoundOptimization(4) = 6.D0
IF( cGeometry(cComponent) == 1 .AND. PotentialTypeLogical(4) ) gLowerBoundOptimization(4) = 0.2D0
IF( cGeometry(cComponent) == 2 .AND. PotentialTypeLogical(4) ) gLowerBoundOptimization(4) = 0.0D0
IF( cGeometry(cComponent) == 3 .AND. PotentialTypeLogical(4) ) gLowerBoundOptimization(4) = 0.2D0
IF( cGeometry(cComponent) == 1 .AND. .NOT. PotentialTypeLogical(4) ) gLowerBoundOptimization(4) = 0.0D0
IF( cGeometry(cComponent) == 2 .AND. .NOT. PotentialTypeLogical(4) ) gLowerBoundOptimization(4) = 0.0D0
IF( cGeometry(cComponent) == 3 .AND. .NOT. PotentialTypeLogical(4) ) gLowerBoundOptimization(4) = 0.0D0

! Upper bound (initial guesses)
UpperBoundOptimization(1) = 2.0D0 * cWellDepth(cComponent)
IF( PotentialTypeLogical(1) ) UpperBoundOptimization(2) = 1.2D0 * cPotentialRange(cComponent) ! Square-well potential
IF( PotentialTypeLogical(2) ) UpperBoundOptimization(2) = 2.0D0 * cPotentialRange(cComponent) ! Sutherland potential
IF( PotentialTypeLogical(3) ) UpperBoundOptimization(2) = 1.5D0 * cPotentialRange(cComponent) ! Yukawa potential
IF( PotentialTypeLogical(4) ) UpperBoundOptimization(2) = 1.2D0 * cPotentialRange(cComponent) ! Convex square-well potential
UpperBoundOptimization(3) = 2.0D0 * cDiameter(cComponent)
UpperBoundOptimization(4) = 2.0D0 * cAspectRatio(cComponent)
IF( UpperBoundOptimization(2) > gUpperBoundOptimization(2) ) UpperBoundOptimization(2) = gUpperBoundOptimization(2)
IF( UpperBoundOptimization(4) > gUpperBoundOptimization(4) ) UpperBoundOptimization(4) = gUpperBoundOptimization(4)

! Lower bound (initial guesses)
LowerBoundOptimization(1) = 0.5D0 * cWellDepth(cComponent)
IF( PotentialTypeLogical(1) ) LowerBoundOptimization(2) = 0.8D0 * cPotentialRange(cComponent) ! Square-well potential
IF( PotentialTypeLogical(2) ) LowerBoundOptimization(2) = 0.5D0 * cPotentialRange(cComponent) ! Sutherland potential
IF( PotentialTypeLogical(3) ) LowerBoundOptimization(2) = 0.7D0 * cPotentialRange(cComponent) ! Yukawa potential
IF( PotentialTypeLogical(4) ) LowerBoundOptimization(2) = 0.8D0 * cPotentialRange(cComponent) ! Convex square-well potential
LowerBoundOptimization(3) = 0.5D0 * cDiameter(cComponent)
LowerBoundOptimization(4) = 0.5D0 * cAspectRatio(cComponent) 
IF( LowerBoundOptimization(2) < gLowerBoundOptimization(2) ) LowerBoundOptimization(2) = gLowerBoundOptimization(2)
IF( LowerBoundOptimization(4) < gLowerBoundOptimization(4) ) LowerBoundOptimization(4) = gLowerBoundOptimization(4)

! Check parent subfolders
INQUIRE( File= "ROUTE_03_OPTIMIZATION/", Exist= FolderExists )
IF( .NOT. FolderExists ) THEN
  CALL SYSTEM( "mkdir ROUTE_03_OPTIMIZATION/" )
END IF
INQUIRE( File= "ROUTE_03_OPTIMIZATION/"//TRIM( DescriptorDate )//"/", Exist= FolderExists )
IF( .NOT. FolderExists ) THEN
  CALL SYSTEM( "mkdir ROUTE_03_OPTIMIZATION/"//TRIM( DescriptorDate )//"/" )
END IF
FileName = "ROUTE_03_OPTIMIZATION/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )//"_"//&
&          TRIM( cMoleculeName(cComponent) )//"_Simplex_Optimization.dat"
OPEN( Unit= 335, File= FileName, Action= "WRITE" )
WRITE( 335, "(G0)" ) "Geometry: "//TRIM( cGeometryName(cComponent) )//""
WRITE( *, "(G0)" ) "Geometry: "//TRIM( cGeometryName(cComponent) )//""
WRITE( 335, "(G0)" ) " "
WRITE( *, "(G0)" ) " "
IF( cGeometry(cComponent) == 1 ) THEN
  WRITE( 335, "(G0)" ) "EOR model: Isihara-Hadwiger Theorem"
  WRITE( *, "(G0)" ) "EOR model: Isihara-Hadwiger Theorem"
ELSE IF( cGeometry(cComponent) == 2 ) THEN
  WRITE( 335, "(G0)" ) "SPC model: Isihara-Hadwiger Theorem"
  WRITE( *, "(G0)" ) "SPC model: Isihara-Hadwiger Theorem"
ELSE IF( cGeometry(cComponent) == 3 ) THEN
  WRITE( 335, "(G0)" ) "CYL model: Isihara-Hadwiger Theorem"
  WRITE( *, "(G0)" ) "CYL model: Isihara-Hadwiger Theorem"
END IF
WRITE( 335, "(G0)" ) " "
WRITE( *, "(G0)" ) " "
IF( ZhangCorrectionLogical ) WRITE( 335, "(G0)" ) "Zhang Correction Applied"
IF( ZhangCorrectionLogical ) WRITE( *, "(G0)" ) "Zhang Correction Applied"
IF( .NOT. ZhangCorrectionLogical ) WRITE( 335, "(G0)" ) "Zhang Correction NOT Applied"
IF( .NOT. ZhangCorrectionLogical ) WRITE( *, "(G0)" ) "Zhang Correction NOT Applied"
WRITE( 335, "(G0)" ) " "
WRITE( *, "(G0)" ) " "
IF( PotentialTypeLogical(4) ) THEN
  IF( PYHCBCorrectionLogical ) WRITE( 335, "(G0)" ) "Percus-Yevick Correction Applied"
  IF( PYHCBCorrectionLogical ) WRITE( *, "(G0)" ) "Percus-Yevick Correction Applied"
  IF( .NOT. PYHCBCorrectionLogical ) WRITE( 335, "(G0)" ) "Percus-Yevick Correction NOT Applied"
  IF( .NOT. PYHCBCorrectionLogical ) WRITE( *, "(G0)" ) "Percus-Yevick Correction NOT Applied"
END IF
WRITE( 335, "(G0)" ) " "
WRITE( *, "(G0)" ) " "
IF( PotentialTypeLogical(4) ) THEN
  IF( UseA1ForA2Logical ) WRITE( 335, "(G0)" ) "Overriding effective packing fraction coefficients for A2"
  IF( UseA1ForA2Logical ) WRITE( *, "(G0)" ) "Overriding effective packing fraction coefficients for A2"
  IF( .NOT. UseA1ForA2Logical ) WRITE( 335, "(G0)" ) "Keeping effective packing fraction coefficients for A2"
  IF( .NOT. UseA1ForA2Logical ) WRITE( *, "(G0)" ) "Keeping effective packing fraction coefficients for A2"
END IF
WRITE( 335, "(G0)" ) " "
WRITE( *, "(G0)" ) " "
IF( PotentialTypeLogical(4) ) THEN
  IF( ReferenceBoublikLogical ) WRITE( 335, "(G0)" ) "Using Boublik's reference term"
  IF( ReferenceBoublikLogical ) WRITE( *, "(G0)" ) "Using Boublik's reference term"
  IF( .NOT. ReferenceBoublikLogical ) WRITE( 335, "(G0)" ) "Using Parsons-Lee's reference term"
  IF( .NOT. ReferenceBoublikLogical ) WRITE( *, "(G0)" ) "Using Parsons-Lee's reference term"
END IF
WRITE( 335, "(G0)" ) " "
WRITE( *, "(G0)" ) " "
WRITE( 335, "(13G0)" ) "'# Simplex Cycles'", ",", "'Objective Function'", ",", "'Well Depth [K]'", ",", "'Well Range'", ",", &
&                      "'Diameter [Å]'", ",", "'Aspect Ratio'", ",", "'Convergence Criterion'"
FLUSH( 335 )
FileName = "ROUTE_03_OPTIMIZATION/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )//"_"//&
&          TRIM( cMoleculeName(cComponent) )//"_Critical_Properties.dat"
OPEN( Unit= 336, File= FileName, Action= "WRITE" )
WRITE( 336, "(G0)" ) "Geometry: "//TRIM( cGeometryName(cComponent) )//""
WRITE( 336, "(G0)" ) " "
IF( cGeometry(cComponent) == 1 ) THEN
  WRITE( 336, "(G0)" ) "EOR model: Isihara-Hadwiger Theorem"
ELSE IF( cGeometry(cComponent) == 2 ) THEN
  WRITE( 336, "(G0)" ) "SPC model: Isihara-Hadwiger Theorem"
ELSE IF( cGeometry(cComponent) == 3 ) THEN
  WRITE( 336, "(G0)" ) "CYL model: Isihara-Hadwiger Theorem"
END IF
WRITE( 336, "(G0)" ) " "
IF( ZhangCorrectionLogical ) WRITE( 336, "(G0)" ) "Zhang Correction Applied"
IF( .NOT. ZhangCorrectionLogical ) WRITE( 336, "(G0)" ) "Zhang Correction NOT Applied"
WRITE( 336, "(G0)" ) " "
IF( PotentialTypeLogical(4) ) THEN
  IF( PYHCBCorrectionLogical ) WRITE( 336, "(G0)" ) "Percus-Yevick Correction Applied"
  IF( .NOT. PYHCBCorrectionLogical ) WRITE( 336, "(G0)" ) "Percus-Yevick Correction NOT Applied"
END IF
WRITE( 336, "(G0)" ) " "
IF( PotentialTypeLogical(4) ) THEN
  IF( UseA1ForA2Logical ) WRITE( 336, "(G0)" ) "Overriding effective packing fraction coefficients for A2"
  IF( .NOT. UseA1ForA2Logical ) WRITE( 336, "(G0)" ) "Keeping effective packing fraction coefficients for A2"
END IF
WRITE( 336, "(G0)" ) " "
WRITE( 336, "(9G0)" ) "'# Simplex Cycles'", ",", "'Critical Temperature [K]'", ",", "'Critical Pressure [MPa]'", ",", &
&                     "'Critical Density [mol/m³]'", ",", "'Critical Volume [m³/mol]'"
FLUSH( 336 )

! Global simplex cycles
DO iSimplex = 1, cGlobalSimplex

  ! Other guesses
  CALL Random_Seed(  )
  DO iParameter = 2, nParameters + 1
    DO jParameter = 1, nParameters
      ! Random number
      CALL Random_Number( RandomNumber )
      ! Variation
      Delta = MinDelta + RandomNumber * MaxDelta
      ! Random number
      CALL Random_Number( RandomNumber )
      ! Change the sign of the variation
      IF( RandomNumber >= 0.D0 .AND. RandomNumber < 0.5D0 ) THEN
        Delta = - Delta
      END IF
      ! Second guess for the optimization parameters (original values towards the positive direction)
      Guesses(iParameter,jParameter) = Guesses(1,jParameter) + Delta * ( UpperBoundOptimization(jParameter) - &
      &                                LowerBoundOptimization(jParameter) )
      IF( Guesses(iParameter,jParameter) < LowerBoundOptimization(jParameter) ) Guesses(iParameter,jParameter) = &
      &   LowerBoundOptimization(jParameter)
      IF( Guesses(iParameter,jParameter) > UpperBoundOptimization(jParameter) ) Guesses(iParameter,jParameter) = &
      &   UpperBoundOptimization(jParameter)
    END DO
  END DO

  ! Status
  WRITE( *, "(5G0)" ) "Initializing Nelder-Mead Simplex Nonlinear Optimization... (", iSimplex, "/", cGlobalSimplex, ")"
  WRITE( *, "(G0)" ) " "

  ! Objective function (initial guesses)
  SimplexType = -1
  DO iGuess = 1, nParameters + 1
    CALL Objective_Function_Pure( nParameters, nPoints, DampingFactorIteration, cCriticalTemperature(cComponent), &
    &    cCriticalDensity(cComponent), cCriticalPressure(cComponent), cCriticalVolume(cComponent), EqLiquidDens, EqVaporPressure, &
    &    EqTemperature, Guesses(iGuess,:), ObjectiveFunction(iGuess), SimplexType, iGuess, 0.D0, .FALSE. )
  END DO

  ! Status
  WRITE( *, "(2G0)" ) CHAR(13), "Calculating Initial Guesses... Done!"
  WRITE( *, "(G0)" ) " "

  ! Sort the objective function values and the parameter sets
  CALL Bubble_Sort( nParameters + 1, ObjectiveFunction, Guesses )

  ! Parameter sets
  ParameterSets = Guesses

  ! Objective function (average)
  ObjectiveFunctionAverage = SUM( ObjectiveFunction ) / ( nParameters + 1 )

  ! Criteria for convergence
  ConvergenceCriterion = 0
  DO iParameter = 1, nParameters + 1
    ConvergenceCriterion = ConvergenceCriterion + ( ObjectiveFunction(iParameter) - ObjectiveFunctionAverage ) * &
    &                      ( ObjectiveFunction(iParameter) - ObjectiveFunctionAverage )
  END DO
  ConvergenceCriterion = ConvergenceCriterion / nParameters
  ConvergenceCriterion = DSQRT( ConvergenceCriterion )

  ! Simplex algorithm
  DO WHILE( ConvergenceCriterion >= Tolerance )

    ! Midpoint (centroid) between sets (excluding the worst set)
    DO iParameter = 1, nParameters
      AuxSum = SUM( ParameterSets(1:nParameters,iParameter), Dim= 1 )
      MidpointSet(iParameter) = AuxSum / nParameters
    END DO
    ! Objective function (midpoint)
    SimplexType = 0
    CALL Objective_Function_Pure( nParameters, nPoints, DampingFactorIteration, cCriticalTemperature(cComponent), &
    &    cCriticalDensity(cComponent), cCriticalPressure(cComponent), cCriticalVolume(cComponent), EqLiquidDens, EqVaporPressure, &
    &    EqTemperature, MidpointSet, ObjectiveFunctionMidpoint, SimplexType, 0_Int64, ConvergenceCriterion, .FALSE. )

    ! Reflection of the worst set
    ReflectionSet(:) = ( 1.D0 + Alpha ) * MidpointSet(:) - Alpha * ParameterSets(nParameters+1,:)
    DO jParameter = 1, nParameters
      IF( ReflectionSet(jParameter) < LowerBoundOptimization(jParameter) ) ReflectionSet(jParameter) = &
      &   LowerBoundOptimization(jParameter)
      IF( ReflectionSet(jParameter) > UpperBoundOptimization(jParameter) ) ReflectionSet(jParameter) = &
      &   UpperBoundOptimization(jParameter)
    END DO
    ! Objective function (reflection)
    SimplexType = 1
    CALL Objective_Function_Pure( nParameters, nPoints, DampingFactorIteration, cCriticalTemperature(cComponent), &
    &    cCriticalDensity(cComponent), cCriticalPressure(cComponent), cCriticalVolume(cComponent), EqLiquidDens, EqVaporPressure, &
    &    EqTemperature, ReflectionSet, ObjectiveFunctionReflection, SimplexType, 0_Int64, ConvergenceCriterion, .FALSE. )

    ! Case selection
    IF( ObjectiveFunctionReflection <  ObjectiveFunction(nParameters) ) CaseSelection = 1 ! Reflection or Expansion
    IF( ObjectiveFunctionReflection >= ObjectiveFunction(nParameters) ) CaseSelection = 2 ! Contraction or Reduction/Shrink

    ! Cases
    SELECT CASE( CaseSelection )

    CASE( 1 ) ! Reflection or Expansion

      ! Apply REFLECTION
      IF( ObjectiveFunction(1) < ObjectiveFunctionReflection ) THEN
        ParameterSets(nParameters+1,:) = ReflectionSet(:)
        ObjectiveFunction(nParameters+1) = ObjectiveFunctionReflection
      ! Attempt EXPANSION
      ELSE
        ! Expansion of the reflection set
        ExpansionSet(:) = Gamma * ReflectionSet(:) + ( 1.D0 - Gamma ) * MidpointSet(:)
        DO jParameter = 1, nParameters
          IF( ExpansionSet(jParameter) < LowerBoundOptimization(jParameter) ) ExpansionSet(jParameter) = &
          &   LowerBoundOptimization(jParameter)
          IF( ExpansionSet(jParameter) > UpperBoundOptimization(jParameter) ) ExpansionSet(jParameter) = &
          &   UpperBoundOptimization(jParameter)
        END DO
        ! Objective function (expansion)
        SimplexType = 2
        CALL Objective_Function_Pure( nParameters, nPoints, DampingFactorIteration, cCriticalTemperature(cComponent), &
        &    cCriticalDensity(cComponent), cCriticalPressure(cComponent), cCriticalVolume(cComponent), EqLiquidDens, &
        &    EqVaporPressure, EqTemperature, ExpansionSet, ObjectiveFunctionExpansion, SimplexType, 0_Int64, ConvergenceCriterion, &
        &    .FALSE. )
        IF( ObjectiveFunctionExpansion < ObjectiveFunction(1) ) THEN ! Apply EXPANSION
          ParameterSets(nParameters+1,:) = ExpansionSet(:)
          ObjectiveFunction(nParameters+1) = ObjectiveFunctionExpansion
        ELSE ! Apply REFLECTION
          ParameterSets(nParameters+1,:) = ReflectionSet(:)
          ObjectiveFunction(nParameters+1) = ObjectiveFunctionReflection
        END IF
      END IF

    CASE( 2 ) ! Contraction or Reduction/Shrink

      ! Replace WORST SET with REFLECTION
      IF( ObjectiveFunctionReflection < ObjectiveFunction(nParameters+1) ) THEN
        ParameterSets(nParameters+1,:) = ReflectionSet(:)
        ObjectiveFunction(nParameters+1) = ObjectiveFunctionReflection
      END IF

      ! Attempt CONTRACTION
      ContractionSet(:) = ( 1.D0 - Beta ) * MidpointSet(:) + ParameterSets(nParameters+1,:) * Beta
      DO jParameter = 1, nParameters
        IF( ContractionSet(jParameter) < LowerBoundOptimization(jParameter) ) ContractionSet(jParameter) = &
        &   LowerBoundOptimization(jParameter)
        IF( ContractionSet(jParameter) > UpperBoundOptimization(jParameter) ) ContractionSet(jParameter) = &
        &   UpperBoundOptimization(jParameter)
      END DO
      ! Objective function (contraction)
      SimplexType = 3
      CALL Objective_Function_Pure( nParameters, nPoints, DampingFactorIteration, cCriticalTemperature(cComponent), &
      &    cCriticalDensity(cComponent), cCriticalPressure(cComponent), cCriticalVolume(cComponent), EqLiquidDens, &
      &    EqVaporPressure, EqTemperature, ContractionSet, ObjectiveFunctionContraction, SimplexType, 0_Int64, &
      &    ConvergenceCriterion, .FALSE. )

      ! Apply CONTRACTION
      IF( ObjectiveFunctionContraction < ObjectiveFunction(nParameters+1) ) THEN
        ParameterSets(nParameters+1,:) = ContractionSet(:)
        ObjectiveFunction(nParameters+1) = ObjectiveFunctionContraction
      ! Apply REDUCTION/SHRINK
      ELSE
        ! Reduction/Shrink of all sets except the best one
        DO jParameter = 2, nParameters + 1
          ReductionSet(:) = 0.5D0 * ( ParameterSets(1,:) + ParameterSets(jParameter,:) ) ! Midpoint between best and other sets
          DO iParameter = 1, nParameters
            IF( ReductionSet(iParameter) < LowerBoundOptimization(iParameter) ) ReductionSet(iParameter) = &
            &   LowerBoundOptimization(iParameter)
            IF( ReductionSet(iParameter) > UpperBoundOptimization(iParameter) ) ReductionSet(iParameter) = &
            &   UpperBoundOptimization(iParameter)
          END DO
          ! Objective function (reduction)
          SimplexType = 4
          CALL Objective_Function_Pure( nParameters, nPoints, DampingFactorIteration, cCriticalTemperature(cComponent), &
          &    cCriticalDensity(cComponent), cCriticalPressure(cComponent), cCriticalVolume(cComponent), EqLiquidDens, &
          &    EqVaporPressure, EqTemperature, ReductionSet, ObjectiveFunctionReduction, SimplexType, 0_Int64, &
          &    ConvergenceCriterion, .FALSE. )
          ParameterSets(jParameter,:) = ReductionSet(:)
          ObjectiveFunction(jParameter) = ObjectiveFunctionReduction
        END DO
      END IF

    END SELECT

    ! Sort the objective function values and the parameter sets
    CALL Bubble_Sort( nParameters + 1, ObjectiveFunction, ParameterSets )

    ! Last (best) objective function
    LastObjectiveFunction = ObjectiveFunction(1)

    ! Objective function (average)
    ObjectiveFunctionAverage = SUM( ObjectiveFunction ) / ( nParameters + 1 )

    ! Criteria for convergence
    ConvergenceCriterion = 0
    DO iParameter = 1, nParameters + 1
      ConvergenceCriterion = ConvergenceCriterion + ( ObjectiveFunction(iParameter) - ObjectiveFunctionAverage ) * &
      &                      ( ObjectiveFunction(iParameter) - ObjectiveFunctionAverage )
    END DO
    ConvergenceCriterion = ConvergenceCriterion / nParameters
    ConvergenceCriterion = DSQRT( ConvergenceCriterion )

    ! Update counter
    nSimplexCycles = nSimplexCycles + 1

    ! Output file
    WRITE( 335, "(I0,5(',',G0.10),',',G0.10)" ) nSimplexCycles - 1, LastObjectiveFunction, ParameterSets(1,1) / cBoltzmann, &
    &                                           ParameterSets(1,2), ParameterSets(1,3) / 1.D-10, ParameterSets(1,4), &
    &                                           ConvergenceCriterion
    FLUSH( 335 )

    ! Critical properties
    IF( SecondLastObjectiveFunction /= LastObjectiveFunction ) THEN
      SecondLastObjectiveFunction = LastObjectiveFunction
      SimplexType = 5
      CALL Objective_Function_Pure( nParameters, nPoints, DampingFactorIteration, cCriticalTemperature(cComponent), &
      &    cCriticalDensity(cComponent), cCriticalPressure(cComponent), cCriticalVolume(cComponent), EqLiquidDens, &
      &    EqVaporPressure, EqTemperature, ParameterSets(1,:), Dummy, SimplexType, 0_Int64, ConvergenceCriterion, .TRUE. )
      WRITE( 336, "(G0,4(',',G0.10))" ) nSimplexCycles - 1, cCriticalTemperature(cComponent), cCriticalPressure(cComponent) / &
      &                                 1.D6, cCriticalDensity(cComponent), cCriticalVolume(cComponent)
      FLUSH( 336 )
    END IF

  END DO

  ! Status
  WRITE( *, "(G0)" ) " "
  WRITE( *, "(G0)" ) " "

  ! Re-initialize guesses with best result from previous simplex
  Guesses(1,1) = ParameterSets(1,1)
  Guesses(1,2) = ParameterSets(1,2)
  Guesses(1,3) = ParameterSets(1,3)
  Guesses(1,4) = ParameterSets(1,4)

  ! Upper bound (initial guesses)
  UpperBoundOptimization(1) = 2.0D0 * Guesses(1,1)
  IF( PotentialTypeLogical(1) ) UpperBoundOptimization(2) = 1.2D0 * Guesses(1,2) ! Square-well potential
  IF( PotentialTypeLogical(2) ) UpperBoundOptimization(2) = 2.0D0 * Guesses(1,2) ! Sutherland potential
  IF( PotentialTypeLogical(3) ) UpperBoundOptimization(2) = 1.5D0 * Guesses(1,2) ! Yukawa potential
  IF( PotentialTypeLogical(4) ) UpperBoundOptimization(2) = 1.2D0 * Guesses(1,2) ! Convex square-well potential
  UpperBoundOptimization(3) = 2.0D0 * Guesses(1,3)
  UpperBoundOptimization(4) = 2.0D0 * Guesses(1,4)
  IF( UpperBoundOptimization(2) > gUpperBoundOptimization(2) ) UpperBoundOptimization(2) = gUpperBoundOptimization(2)
  IF( UpperBoundOptimization(4) > gUpperBoundOptimization(4) ) UpperBoundOptimization(4) = gUpperBoundOptimization(4)

  ! Lower bound (initial guesses)
  LowerBoundOptimization(1) = 0.5D0 * Guesses(1,1)
  IF( PotentialTypeLogical(1) ) LowerBoundOptimization(2) = 0.8D0 * Guesses(1,2) ! Square-well potential
  IF( PotentialTypeLogical(2) ) LowerBoundOptimization(2) = 0.5D0 * Guesses(1,2) ! Sutherland potential
  IF( PotentialTypeLogical(3) ) LowerBoundOptimization(2) = 0.7D0 * Guesses(1,2) ! Yukawa potential
  IF( PotentialTypeLogical(4) ) LowerBoundOptimization(2) = 0.8D0 * Guesses(1,2) ! Convex square-well potential
  LowerBoundOptimization(3) = 0.5D0 * Guesses(1,3)
  LowerBoundOptimization(4) = 0.5D0 * Guesses(1,4)
  IF( LowerBoundOptimization(2) < gLowerBoundOptimization(2) ) LowerBoundOptimization(2) = gLowerBoundOptimization(2)
  IF( LowerBoundOptimization(4) < gLowerBoundOptimization(4) ) LowerBoundOptimization(4) = gLowerBoundOptimization(4)

END DO

! Close files
CLOSE( 335 )
CLOSE( 336 )

! Status
WRITE( *, "(G0)" ) "Finishing Nelder-Mead Simplex Nonlinear Optimization..."
WRITE( *, "(G0)" ) " "

! Last parameters
WRITE( *, "(2G0)" ) "Number of Simplex Cycles: ", nSimplexCycles - 1
WRITE( *, "(G0,G0.8)" ) "Objective Function Minimized: ", LastObjectiveFunction
WRITE( *, "(G0,G0.8)" ) "Optimized Well Depth [K]: ", ParameterSets(1,1) / cBoltzmann
WRITE( *, "(G0,G0.8)" ) "Optimized Well Range: ", ParameterSets(1,2)
WRITE( *, "(G0,G0.8)" ) "Optimized Diameter [Å]: ", ParameterSets(1,3) / 1.D-10
WRITE( *, "(G0,G0.8)" ) "Optimized Aspect Ratio: ", ParameterSets(1,4)
WRITE( *, "(G0,G0.8)" ) "Convergence Criterion: ", ConvergenceCriterion
WRITE( *, "(G0)" ) " "

! Critical properties
SimplexType = 5
CALL Objective_Function_Pure( nParameters, nPoints, DampingFactorIteration, cCriticalTemperature(cComponent), &
&    cCriticalDensity(cComponent), cCriticalPressure(cComponent), cCriticalVolume(cComponent), EqLiquidDens, EqVaporPressure, &
&    EqTemperature, ParameterSets(1,:), Dummy, SimplexType, 0_Int64, ConvergenceCriterion, .TRUE. )
WRITE( *, "(G0,G0.8)" ) "Optimized Critical Temperature [K]: ", cCriticalTemperature(cComponent)
WRITE( *, "(G0,G0.8)" ) "Optimized Critical Pressure [MPa]: ", cCriticalPressure(cComponent) / 1.D6
WRITE( *, "(G0,G0.8)" ) "Optimized Critical Density [mol/m³]: ", cCriticalDensity(cComponent)
WRITE( *, "(G0,G0.8)" ) "Optimized Critical Volume [m³/mol]: ", cCriticalVolume(cComponent)
WRITE( *, "(G0)" ) " "

RETURN

END SUBROUTINE Simplex_Algorithm_Pure

! ************************************************************************************************ !
!                          Simplex algorithm to search zero of functions                           !
!                Reference: Nelder, J. A.; Mead, R. A., Comp. J., 7, 308-313, 1965                 !
! ************************************************************************************************ !
SUBROUTINE Simplex_Algorithm_Mixture( nParameters, nPoints, cGlobalSimplex, Temperature, cAccentricFactor, cCriticalPressure, &
&                                     cCriticalTemperature, EqLiquidComposition, EqVaporComposition, EqPressure )

USE GlobalVariables
USE Substances

IMPLICIT NONE

! ************************************************************************************************ !
! INTEGER VARIABLES                                                                                !
! ************************************************************************************************ !
INTEGER( Kind= Int64 ) :: cComponent     ! Counter (components)
INTEGER( Kind= Int64 ) :: iGuess         ! Counter (guesses)
INTEGER( Kind= Int64 ) :: iSimplex       ! Counter (simplex)
INTEGER( Kind= Int64 ) :: iParameter     ! Counter (parameters)
INTEGER( Kind= Int64 ) :: jParameter     ! Counter (parameters)
INTEGER( Kind= Int64 ) :: nPoints        ! Number of points
INTEGER( Kind= Int64 ) :: nParameters    ! Number of optimization parameters
INTEGER( Kind= Int64 ) :: CaseSelection  ! Case selection (1: Reflection or Expansion; 2: Contraction or Reduction/Shrink)
INTEGER( Kind= Int64 ) :: SimplexType    ! Simplex type (-1: Initial guesses; 0: Midpoint; 1: Reflection; 2: Expansion; 3: Contraction; 4: Reduction)
INTEGER( Kind= Int64 ) :: cGlobalSimplex ! Number of global simplex cycles

! ************************************************************************************************ !
! REAL PARAMETERS                                                                                  !
! ************************************************************************************************ !
REAL( Kind= Real64 ), PARAMETER :: Tolerance = 1.D-8 ! Expansion coefficient
REAL( Kind= Real64 ), PARAMETER :: Alpha     = 1.D0  ! Reflection coefficient
REAL( Kind= Real64 ), PARAMETER :: Beta      = 0.5D0 ! Contraction coefficient
REAL( Kind= Real64 ), PARAMETER :: Gamma     = 2.D0  ! Expansion coefficient
REAL( Kind= Real64 ), PARAMETER :: MinDelta  = 1.D-5 ! Minimum variation
REAL( Kind= Real64 ), PARAMETER :: MaxDelta  = 1.D-4 ! Maximum variation

! ************************************************************************************************ !
! REAL VARIABLES                                                                                   !
! ************************************************************************************************ !
REAL( Kind= Real64 )                                          :: Temperature                   ! Temperature
REAL( Kind= Real64 )                                          :: RandomNumber                  ! Random number
REAL( Kind= Real64 )                                          :: Delta                         ! Variation
REAL( Kind= Real64 )                                          :: ObjectiveFunctionMidpoint     ! Objective function (midpoint)
REAL( Kind= Real64 )                                          :: ObjectiveFunctionReflection   ! Objective function (reflection)
REAL( Kind= Real64 )                                          :: ObjectiveFunctionExpansion    ! Objective function (expansion)
REAL( Kind= Real64 )                                          :: ObjectiveFunctionContraction  ! Objective function (contraction)
REAL( Kind= Real64 )                                          :: ObjectiveFunctionReduction    ! Objective function (reduction)
REAL( Kind= Real64 )                                          :: ObjectiveFunctionAverage      ! Objective function (average)
REAL( Kind= Real64 )                                          :: ConvergenceCriterion          ! Convergence criterion
REAL( Kind= Real64 )                                          :: SecondLastObjectiveFunction   ! Second last value of the objective function
REAL( Kind= Real64 )                                          :: AuxSum                        ! Auxiliary sum
REAL( Kind= Real64 ), DIMENSION( nPoints )                    :: EqLiquidComposition           ! Liquid composition (equilibrium)
REAL( Kind= Real64 ), DIMENSION( nPoints )                    :: EqVaporComposition            ! Vapor composition (equilibrium)
REAL( Kind= Real64 ), DIMENSION( nPoints )                    :: EqPressure                    ! Pressure (equilibrium)
REAL( Kind= Real64 ), DIMENSION( nComponents )                :: cAccentricFactor              ! Accentric factor
REAL( Kind= Real64 ), DIMENSION( nComponents )                :: cCriticalPressure             ! Pressure (critical)
REAL( Kind= Real64 ), DIMENSION( nComponents )                :: cCriticalTemperature          ! Temperature (critical)
REAL( Kind= Real64 ), DIMENSION( nParameters )                :: LowerBoundOptimization        ! Lower bound
REAL( Kind= Real64 ), DIMENSION( nParameters )                :: UpperBoundOptimization        ! Upper bound
REAL( Kind= Real64 ), DIMENSION( nParameters )                :: MidpointSet                   ! Set of parameters (midpoint)
REAL( Kind= Real64 ), DIMENSION( nParameters )                :: ReflectionSet                 ! Set of parameters (reflection)
REAL( Kind= Real64 ), DIMENSION( nParameters )                :: ExpansionSet                  ! Set of parameters (expansion)
REAL( Kind= Real64 ), DIMENSION( nParameters )                :: ContractionSet                ! Set of parameters (contraction)
REAL( Kind= Real64 ), DIMENSION( nParameters )                :: ReductionSet                  ! Set of parameters (reduction)
REAL( Kind= Real64 ), DIMENSION( nParameters+1 )              :: ObjectiveFunction             ! Objective function
REAL( Kind= Real64 ), DIMENSION( nParameters+1, nParameters ) :: ParameterSets                 ! Parameter sets
REAL( Kind= Real64 ), DIMENSION( nParameters+1, nParameters ) :: Guesses                       ! Guesses

! ************************************************************************************************ !
! CHARACTER STRINGS                                                                                !
! ************************************************************************************************ !
CHARACTER( Len= 100 ) :: FileName ! File name

! ************************************************************************************************ !
! LOGICAL VARIABLES                                                                                !
! ************************************************************************************************ !
LOGICAL :: FolderExists ! Flag to indicate if the folder exists

! Initialization
nSimplexCycles = 1
LastObjectiveFunction = -1.D0
SecondLastObjectiveFunction = LastObjectiveFunction

! Initial guess for the optimization parameters (original values)
Guesses(1,1) = ijWellDepthCorrection ! Well depth correction
IF( nParameters == 2 ) Guesses(1,nParameters) = ijPotentialRangeCorrection ! Potential range correction

! Upper bound (initial guesses)
UpperBoundOptimization(1) = 2.0D0
IF( nParameters == 2 ) UpperBoundOptimization(nParameters) = 2.0D0

! Lower bound (initial guesses)
LowerBoundOptimization(1) = 0.0D0
IF( nParameters == 2 ) LowerBoundOptimization(nParameters) = 0.0D0

! Check parent subfolders
INQUIRE( File= "ROUTE_03_OPTIMIZATION/", Exist= FolderExists )
IF( .NOT. FolderExists ) THEN
  CALL SYSTEM( "mkdir ROUTE_03_OPTIMIZATION/" )
END IF
INQUIRE( File= "ROUTE_03_OPTIMIZATION/"//TRIM( DescriptorDate )//"/", Exist= FolderExists )
IF( .NOT. FolderExists ) THEN
  CALL SYSTEM( "mkdir ROUTE_03_OPTIMIZATION/"//TRIM( DescriptorDate )//"/" )
END IF
FileName = "ROUTE_03_OPTIMIZATION/"//TRIM( DescriptorDate )//"/"//TRIM( DescriptorHour )//"_"//TRIM( cFormulaName(1) )//"+"//&
&          TRIM( cFormulaName(2) )//"_Simplex_Optimization.dat"
OPEN( Unit= 335, File= FileName, Action= "WRITE" )
WRITE( 335, "(G0)", Advance= "NO" ) "Mixture: "
WRITE( *, "(G0)", Advance= "NO" ) "Mixture: "
DO cComponent = 1, nComponents - 1
  WRITE( 335, "(3G0)", Advance= "NO" ) ""//TRIM( cFormulaName(cComponent) )//"(", cComponent, ") + "
  WRITE( *, "(3G0)", Advance= "NO" ) ""//TRIM( cFormulaName(cComponent) )//"(", cComponent, ") + "
END DO
WRITE( 335, "(3G0)", Advance= "YES" ) ""//TRIM( cFormulaName(cComponent) )//"(", nComponents, ")"
WRITE( *, "(3G0)", Advance= "YES" ) ""//TRIM( cFormulaName(cComponent) )//"(", nComponents, ")"
WRITE( 335, "(G0)" ) " "
WRITE( *, "(G0)" ) " "
WRITE( 335, "(G0)", Advance= "NO" ) "Geometry: "
WRITE( *, "(G0)", Advance= "NO" ) "Geometry: "
DO cComponent = 1, nComponents - 1
  WRITE( 335, "(3G0)", Advance= "NO" ) ""//TRIM( cGeometryName(cComponent) )//"(", cComponent, ") + "
  WRITE( *, "(3G0)", Advance= "NO" ) ""//TRIM( cGeometryName(cComponent) )//"(", cComponent, ") + "
END DO
WRITE( 335, "(3G0)", Advance= "YES" ) ""//TRIM( cGeometryName(cComponent) )//"(", nComponents, ")"
WRITE( *, "(3G0)", Advance= "YES" ) ""//TRIM( cGeometryName(cComponent) )//"(", nComponents, ")"
WRITE( 335, "(G0)" ) " "
WRITE( *, "(G0)" ) " "
IF( ANY( cGeometry == 1 ) ) THEN
  WRITE( 335, "(G0)" ) "EOR model: Isihara-Hadwiger Theorem"
  WRITE( *, "(G0)" ) "EOR model: Isihara-Hadwiger Theorem"
  WRITE( 335, "(G0)" ) " "
  WRITE( *, "(G0)" ) " "
END IF
IF( ANY( cGeometry == 2 ) ) THEN
  WRITE( 335, "(G0)" ) "SPC model: Isihara-Hadwiger Theorem"
  WRITE( *, "(G0)" ) "SPC model: Isihara-Hadwiger Theorem"
  WRITE( 335, "(G0)" ) " "
  WRITE( *, "(G0)" ) " "
END IF
IF( ANY( cGeometry == 3 ) ) THEN
  WRITE( 335, "(G0)" ) "CYL model: Isihara-Hadwiger Theorem"
  WRITE( *, "(G0)" ) "CYL model: Isihara-Hadwiger Theorem"
  WRITE( 335, "(G0)" ) " "
  WRITE( *, "(G0)" ) " "
END IF
IF( ZhangCorrectionLogical ) WRITE( 335, "(G0)" ) "Zhang Correction Applied"
IF( ZhangCorrectionLogical ) WRITE( *, "(G0)" ) "Zhang Correction Applied"
IF( .NOT. ZhangCorrectionLogical ) WRITE( 335, "(G0)" ) "Zhang Correction NOT Applied"
IF( .NOT. ZhangCorrectionLogical ) WRITE( *, "(G0)" ) "Zhang Correction NOT Applied"
IF( PotentialTypeLogical(4) ) THEN
  IF( PYHCBCorrectionLogical ) WRITE( 335, "(G0)" ) "Percus-Yevick Correction Applied"
  IF( PYHCBCorrectionLogical ) WRITE( *, "(G0)" ) "Percus-Yevick Correction Applied"
  IF( .NOT. PYHCBCorrectionLogical ) WRITE( 335, "(G0)" ) "Percus-Yevick Correction NOT Applied"
  IF( .NOT. PYHCBCorrectionLogical ) WRITE( *, "(G0)" ) "Percus-Yevick Correction NOT Applied"
END IF
WRITE( 335, "(G0)" ) " "
WRITE( *, "(G0)" ) " "
IF( PotentialTypeLogical(4) ) THEN
  IF( UseA1ForA2Logical ) WRITE( 335, "(G0)" ) "Overriding effective packing fraction coefficients for A2"
  IF( UseA1ForA2Logical ) WRITE( *, "(G0)" ) "Overriding effective packing fraction coefficients for A2"
  IF( .NOT. UseA1ForA2Logical ) WRITE( 335, "(G0)" ) "Keeping effective packing fraction coefficients for A2"
  IF( .NOT. UseA1ForA2Logical ) WRITE( *, "(G0)" ) "Keeping effective packing fraction coefficients for A2"
END IF
WRITE( 335, "(G0)" ) " "
WRITE( *, "(G0)" ) " "
IF( NonSphericalMixingRule == 1 ) THEN
  WRITE( 335, "(G0)" ) "Non-Spherical Mixing Rule: Aspect Ratio Rule"
  WRITE( *, "(G0)" ) "Non-Spherical Mixing Rule: Aspect Ratio Rule"
ELSE IF( NonSphericalMixingRule == 2 ) THEN
  WRITE( 335, "(G0)" ) "Non-Spherical Mixing Rule: Angle Average Rule"
  WRITE( *, "(G0)" ) "Non-Spherical Mixing Rule: Angle Average Rule"
ELSE IF( NonSphericalMixingRule == 3 ) THEN
  WRITE( 335, "(G0)" ) "Non-Spherical Mixing Rule: Ratio of Second Virial Coefficients"
  WRITE( *, "(G0)" ) "Non-Spherical Mixing Rule: Ratio of Second Virial Coefficients"
ELSE IF( NonSphericalMixingRule == 4 ) THEN
  WRITE( 335, "(G0)" ) "Non-Spherical Mixing Rule: Second Virial Coefficients"
  WRITE( *, "(G0)" ) "Non-Spherical Mixing Rule: Second Virial Coefficients"
END IF
WRITE( 335, "(G0)" ) " "
WRITE( *, "(G0)" ) " "
WRITE( 335, "(G0)" ) "Zhang Factor of the Mixture: Mixed Spherical Diameter"
WRITE( *, "(G0)" ) "Zhang Factor of the Mixture: Mixed Spherical Diameter"
WRITE( 335, "(G0)" ) " "
WRITE( *, "(G0)" ) " "
IF( EffPFractionMixingRule == 1 ) THEN
  WRITE( 335, "(G0)" ) "Effective Packing Fraction Mixing Rule: Overall Reduced Density Mixing Rule"
  WRITE( *, "(G0)" ) "Effective Packing Fraction Mixing Rule: Overall Reduced Density Mixing Rule"
ELSE IF( EffPFractionMixingRule == 2 ) THEN
  WRITE( 335, "(G0)" ) "Effective Packing Fraction Mixing Rule: One-Fluid van der Waals Mixing Rule"
  WRITE( *, "(G0)" ) "Effective Packing Fraction Mixing Rule: One-Fluid van der Waals Mixing Rule"
END IF
WRITE( 335, "(G0)" ) " "
WRITE( *, "(G0)" ) " "
IF( nParameters == 1 ) THEN
  WRITE( 335, "(7G0)" ) "'# Simplex Cycles'", ",", "'Objective Function'", ",", "'Well Depth Correction'", ",", &
  &                     "'Convergence Criterion'"
ELSE IF( nParameters == 2 ) THEN
  WRITE( 335, "(9G0)" ) "'# Simplex Cycles'", ",", "'Objective Function'", ",", "'Well Depth Correction'", ",", &
  &                     "'Well Range Correction'", ",", "'Convergence Criterion'"
END IF
FLUSH( 335 )

! Global simplex cycles
DO iSimplex = 1, cGlobalSimplex

  ! Other guesses
  CALL Random_Seed(  )
  DO iParameter = 2, nParameters + 1
    DO jParameter = 1, nParameters
      ! Random number
      CALL Random_Number( RandomNumber )
      ! Variation
      Delta = MinDelta + RandomNumber * ( MaxDelta - MinDelta )
      ! Random number
      CALL Random_Number( RandomNumber )
      ! Change the sign of the variation
      IF( RandomNumber >= 0.D0 .AND. RandomNumber < 0.5D0 ) THEN
        Delta = - Delta
      END IF
      ! Second guess for the optimization parameters (original values towards the positive direction)
      Guesses(iParameter,jParameter) = Guesses(1,jParameter) + Delta * ( UpperBoundOptimization(jParameter) - &
      &                                LowerBoundOptimization(jParameter) )
    END DO
  END DO

  ! Status
  WRITE( *, "(5G0)" ) "Initializing Nelder-Mead Simplex Nonlinear Optimization... (", iSimplex, "/", cGlobalSimplex, ")"
  WRITE( *, "(G0)" ) " "

  ! Objective function (initial guesses)
  SimplexType = -1
  DO iGuess = 1, nParameters + 1
    CALL Objective_Function_Mixture( nParameters, nPoints, Temperature, cAccentricFactor, cCriticalPressure, cCriticalTemperature, &
    &    EqLiquidComposition, EqVaporComposition, EqPressure, Guesses(iGuess,:), ObjectiveFunction(iGuess), SimplexType, iGuess, &
    &    0.D0 )
  END DO

  ! Status
  WRITE( *, "(2G0)" ) CHAR(13), "Calculating Initial Guesses... Done!"
  WRITE( *, "(G0)" ) " "

  ! Sort the objective function values and the parameter sets
  CALL Bubble_Sort( nParameters + 1, ObjectiveFunction, Guesses )

  ! Parameter sets
  ParameterSets = Guesses

  ! Objective function (average)
  ObjectiveFunctionAverage = SUM( ObjectiveFunction ) / ( nParameters + 1 )

  ! Criteria for convergence
  ConvergenceCriterion = 0
  DO iParameter = 1, nParameters + 1
    ConvergenceCriterion = ConvergenceCriterion + ( ObjectiveFunction(iParameter) - ObjectiveFunctionAverage ) * &
    &                      ( ObjectiveFunction(iParameter) - ObjectiveFunctionAverage )
  END DO
  ConvergenceCriterion = ConvergenceCriterion / nParameters
  ConvergenceCriterion = DSQRT( ConvergenceCriterion )

  ! Simplex algorithm
  DO WHILE( ConvergenceCriterion >= Tolerance )

    ! Midpoint (centroid) between sets (excluding the worst set)
    DO iParameter = 1, nParameters
      AuxSum = SUM( ParameterSets(1:nParameters,iParameter), Dim= 1 )
      MidpointSet(iParameter) = AuxSum / nParameters
    END DO
    ! Objective function (midpoint)
    SimplexType = 0
    CALL Objective_Function_Mixture( nParameters, nPoints, Temperature, cAccentricFactor, cCriticalPressure, cCriticalTemperature, &
    &    EqLiquidComposition, EqVaporComposition, EqPressure, MidpointSet, ObjectiveFunctionMidpoint, SimplexType, 0_Int64, &
    &    ConvergenceCriterion )

    ! Reflection of the worst set
    ReflectionSet(:) = ( 1.D0 + Alpha ) * MidpointSet(:) - Alpha * ParameterSets(nParameters+1,:)
    ! Objective function (reflection)
    SimplexType = 1
    CALL Objective_Function_Mixture( nParameters, nPoints, Temperature, cAccentricFactor, cCriticalPressure, cCriticalTemperature, &
    &    EqLiquidComposition, EqVaporComposition, EqPressure, ReflectionSet, ObjectiveFunctionReflection, SimplexType, 0_Int64, &
    &    ConvergenceCriterion )

    ! Case selection
    IF( ObjectiveFunctionReflection <  ObjectiveFunction(nParameters) ) CaseSelection = 1 ! Reflection or Expansion
    IF( ObjectiveFunctionReflection >= ObjectiveFunction(nParameters) ) CaseSelection = 2 ! Contraction or Reduction/Shrink

    ! Cases
    SELECT CASE( CaseSelection )

    CASE( 1 ) ! Reflection or Expansion

      ! Apply REFLECTION
      IF( ObjectiveFunction(1) < ObjectiveFunctionReflection ) THEN
        ParameterSets(nParameters+1,:) = ReflectionSet(:)
        ObjectiveFunction(nParameters+1) = ObjectiveFunctionReflection
      ! Attempt EXPANSION
      ELSE
        ! Expansion of the reflection set
        ExpansionSet(:) = Gamma * ReflectionSet(:) + ( 1.D0 - Gamma ) * MidpointSet(:)
        ! Objective function (expansion)
        SimplexType = 2
        CALL Objective_Function_Mixture( nParameters, nPoints, Temperature, cAccentricFactor, cCriticalPressure, &
        &    cCriticalTemperature, EqLiquidComposition, EqVaporComposition, EqPressure, ExpansionSet, ObjectiveFunctionExpansion, &
        &    SimplexType, 0_Int64, ConvergenceCriterion )
        IF( ObjectiveFunctionExpansion < ObjectiveFunction(1) ) THEN ! Apply EXPANSION
          ParameterSets(nParameters+1,:) = ExpansionSet(:)
          ObjectiveFunction(nParameters+1) = ObjectiveFunctionExpansion
        ELSE ! Apply REFLECTION
          ParameterSets(nParameters+1,:) = ReflectionSet(:)
          ObjectiveFunction(nParameters+1) = ObjectiveFunctionReflection
        END IF
      END IF

    CASE( 2 ) ! Contraction or Reduction/Shrink

      ! Replace WORST SET with REFLECTION
      IF( ObjectiveFunctionReflection < ObjectiveFunction(nParameters+1) ) THEN
        ParameterSets(nParameters+1,:) = ReflectionSet(:)
        ObjectiveFunction(nParameters+1) = ObjectiveFunctionReflection
      END IF

      ! Attempt CONTRACTION
      ContractionSet(:) = ( 1.D0 - Beta ) * MidpointSet(:) + ParameterSets(nParameters+1,:) * Beta
      ! Objective function (contraction)
      SimplexType = 3
      CALL Objective_Function_Mixture( nParameters, nPoints, Temperature, cAccentricFactor, cCriticalPressure, &
      &    cCriticalTemperature, EqLiquidComposition, EqVaporComposition, EqPressure, ContractionSet, &
      &    ObjectiveFunctionContraction, SimplexType, 0_Int64, ConvergenceCriterion )

      ! Apply CONTRACTION
      IF( ObjectiveFunctionContraction < ObjectiveFunction(nParameters+1) ) THEN
        ParameterSets(nParameters+1,:) = ContractionSet(:)
        ObjectiveFunction(nParameters+1) = ObjectiveFunctionContraction
      ! Apply REDUCTION/SHRINK
      ELSE
        ! Reduction/Shrink of all sets except the best one
        DO jParameter = 2, nParameters + 1
          ReductionSet(:) = 0.5D0 * ( ParameterSets(1,:) + ParameterSets(jParameter,:) ) ! Midpoint between best and other sets
          ! Objective function (reduction)
          SimplexType = 4
          CALL Objective_Function_Mixture( nParameters, nPoints, Temperature, cAccentricFactor, cCriticalPressure, &
          &    cCriticalTemperature, EqLiquidComposition, EqVaporComposition, EqPressure, ReductionSet, &
          &    ObjectiveFunctionReduction, SimplexType, 0_Int64, ConvergenceCriterion )
          ParameterSets(jParameter,:) = ReductionSet(:)
          ObjectiveFunction(jParameter) = ObjectiveFunctionReduction
        END DO
      END IF

    END SELECT

    ! Sort the objective function values and the parameter sets
    CALL Bubble_Sort( nParameters + 1, ObjectiveFunction, ParameterSets )

    ! Last (best) objective function
    LastObjectiveFunction = ObjectiveFunction(1)

    ! Objective function (average)
    ObjectiveFunctionAverage = SUM( ObjectiveFunction ) / ( nParameters + 1 )

    ! Criteria for convergence
    ConvergenceCriterion = 0
    DO iParameter = 1, nParameters + 1
      ConvergenceCriterion = ConvergenceCriterion + ( ObjectiveFunction(iParameter) - ObjectiveFunctionAverage ) * &
      &                      ( ObjectiveFunction(iParameter) - ObjectiveFunctionAverage )
    END DO
    ConvergenceCriterion = ConvergenceCriterion / nParameters
    ConvergenceCriterion = DSQRT( ConvergenceCriterion )

    ! Update counter
    nSimplexCycles = nSimplexCycles + 1

    ! Output file
    IF( nParameters == 1 ) THEN
      WRITE( 335, "(I0,2(',',G0.10),',',G0.10)" ) nSimplexCycles - 1, LastObjectiveFunction, ParameterSets(1,1), &
      &                                           ConvergenceCriterion
    ELSE IF( nParameters == 2 ) THEN
      WRITE( 335, "(I0,3(',',G0.10),',',G0.10)" ) nSimplexCycles - 1, LastObjectiveFunction, ParameterSets(1,1), &
      &                                           ParameterSets(1,2), ConvergenceCriterion
    END IF
    FLUSH( 335 )

    ! Critical properties
    IF( SecondLastObjectiveFunction /= LastObjectiveFunction ) THEN
      SecondLastObjectiveFunction = LastObjectiveFunction
    END IF

  END DO

  ! Status
  WRITE( *, "(G0)" ) " "
  WRITE( *, "(G0)" ) " "

  ! Re-initialize guesses with best result from previous simplex
  Guesses(1,1) = ParameterSets(1,1)
  IF( nParameters == 2 ) Guesses(1,nParameters) = ParameterSets(1,2)

  ! Upper bound (initial guesses)
  UpperBoundOptimization(1) = 2.0D0
  IF( nParameters == 2 ) UpperBoundOptimization(nParameters) = 2.0D0

  ! Lower bound (initial guesses)
  LowerBoundOptimization(1) = 0.0D0
  IF( nParameters == 2 ) LowerBoundOptimization(nParameters) = 0.0D0

END DO

RETURN

END SUBROUTINE Simplex_Algorithm_Mixture

! ************************************************************************************************ !
!                          Simplex algorithm to search zero of functions                           !
!                Reference: Nelder, J. A.; Mead, R. A., Comp. J., 7, 308-313, 1965                 !
! ************************************************************************************************ !
SUBROUTINE Simplex_test( nParameters, cGlobalSimplex, Temperature, Pressure, INIVAPOR, INILIQUID,ERROR )

USE GlobalVariables
USE Substances

IMPLICIT NONE

! ************************************************************************************************ !
! INTEGER VARIABLES                                                                                !
! ************************************************************************************************ !
INTEGER( Kind= Int64 ) :: cComponent     ! Counter (components)
INTEGER( Kind= Int64 ) :: iGuess         ! Counter (guesses)
INTEGER( Kind= Int64 ) :: iSimplex       ! Counter (simplex)
INTEGER( Kind= Int64 ) :: iParameter     ! Counter (parameters)
INTEGER( Kind= Int64 ) :: jParameter     ! Counter (parameters)
INTEGER( Kind= Int64 ) :: nPoints        ! Number of points
INTEGER( Kind= Int64 ) :: nParameters    ! Number of optimization parameters
INTEGER( Kind= Int64 ) :: CaseSelection  ! Case selection (1: Reflection or Expansion; 2: Contraction or Reduction/Shrink)
INTEGER( Kind= Int64 ) :: SimplexType    ! Simplex type (-1: Initial guesses; 0: Midpoint; 1: Reflection; 2: Expansion; 3: Contraction; 4: Reduction)
INTEGER( Kind= Int64 ) :: cGlobalSimplex ! Number of global simplex cycles

! ************************************************************************************************ !
! REAL PARAMETERS                                                                                  !
! ************************************************************************************************ !
REAL( Kind= Real64 ), PARAMETER :: Tolerance = 1.D-6 ! Expansion coefficient
REAL( Kind= Real64 ), PARAMETER :: Alpha     = 1.D0  ! Reflection coefficient
REAL( Kind= Real64 ), PARAMETER :: Beta      = 0.5D0 ! Contraction coefficient
REAL( Kind= Real64 ), PARAMETER :: Gamma     = 2.D0  ! Expansion coefficient
REAL( Kind= Real64 ), PARAMETER :: MinDelta  = 1.D-5 ! Minimum variation
REAL( Kind= Real64 ), PARAMETER :: MaxDelta  = 1.D-4 ! Maximum variation

! ************************************************************************************************ !
! REAL VARIABLES                                                                                   !
! ************************************************************************************************ !
REAL( Kind= Real64 )                                          :: Temperature, Pressure                   ! Temperature
REAL( Kind= Real64 )                                          :: RandomNumber                  ! Random number
REAL( Kind= Real64 )                           :: mVolumeLiquidPhase          ! Molar volume of the liquid phase
REAL( Kind= Real64 )                           :: mVolumeVaporPhase           ! Molar volume of the vapor phase
REAL( Kind= Real64 )                                          :: Delta                         ! Variation
REAL( Kind= Real64 )                           :: CompressibilityFactor       ! Compressibility factor
REAL( Kind= Real64 )                           :: IsothermalCompressibility   ! Isothermal compressibility
REAL( Kind= Real64 )                           :: ThermalExpansionCoefficient ! Thermal expansion coefficient
REAL( Kind= Real64 )                           :: Anynumber,INILIQUID,INIVAPOR                   ! Any number (dummy)
REAL( Kind= Real64 )                           :: Error                       ! Difference between old and new equilibrium factors
REAL( Kind= Real64 )                                          :: ObjectiveFunctionMidpoint     ! Objective function (midpoint)
REAL( Kind= Real64 )                                          :: ObjectiveFunctionReflection   ! Objective function (reflection)
REAL( Kind= Real64 )                                          :: ObjectiveFunctionExpansion    ! Objective function (expansion)
REAL( Kind= Real64 )                                          :: ObjectiveFunctionContraction  ! Objective function (contraction)
REAL( Kind= Real64 )                                          :: ObjectiveFunctionReduction    ! Objective function (reduction)
REAL( Kind= Real64 )                                          :: ObjectiveFunctionAverage      ! Objective function (average)
REAL( Kind= Real64 )                                          :: ConvergenceCriterion          ! Convergence criterion
REAL( Kind= Real64 )                                          :: SecondLastObjectiveFunction   ! Second last value of the objective function
REAL( Kind= Real64 )                                          :: AuxSum                        ! Auxiliary sum
REAL( Kind= Real64 ), DIMENSION( nComponents )                :: cAccentricFactor              ! Accentric factor
REAL( Kind= Real64 ), DIMENSION( nComponents )                :: cCriticalPressure             ! Pressure (critical)
REAL( Kind= Real64 ), DIMENSION( nComponents )                :: cCriticalTemperature          ! Temperature (critical)
REAL( Kind= Real64 ), DIMENSION( nParameters )                :: LowerBoundOptimization        ! Lower bound
REAL( Kind= Real64 ), DIMENSION( nParameters )                :: UpperBoundOptimization        ! Upper bound
REAL( Kind= Real64 ), DIMENSION( nParameters )                :: MidpointSet                   ! Set of parameters (midpoint)
REAL( Kind= Real64 ), DIMENSION( nParameters )                :: ReflectionSet                 ! Set of parameters (reflection)
REAL( Kind= Real64 ), DIMENSION( nParameters )                :: ExpansionSet                  ! Set of parameters (expansion)
REAL( Kind= Real64 ), DIMENSION( nParameters )                :: ContractionSet                ! Set of parameters (contraction)
REAL( Kind= Real64 ), DIMENSION( nParameters )                :: ReductionSet                  ! Set of parameters (reduction)
REAL( Kind= Real64 ), DIMENSION( nParameters+1 )              :: ObjectiveFunction             ! Objective function
REAL( Kind= Real64 ), DIMENSION( nParameters+1, nParameters ) :: ParameterSets                 ! Parameter sets
REAL( Kind= Real64 ), DIMENSION( nParameters+1, nParameters ) :: Guesses                       ! Guesses
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: TotalChemicalPotential      ! Total chemical potential
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: ResidualChemicalPotential   ! Residual chemical potential
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: liqpot, vappot        ! New equilibrium factors
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: VaporComposition            ! Vapor composition
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: LiquidComposition           ! Liquid composition

! ************************************************************************************************ !
! CHARACTER STRINGS                                                                                !
! ************************************************************************************************ !
CHARACTER( Len= 100 ) :: FileName ! File name
CHARACTER( Len= 01 ) :: CurveType ! Isotherm types (A, B, or C)

! ************************************************************************************************ !
! LOGICAL VARIABLES                                                                                !
! ************************************************************************************************ !
LOGICAL :: FolderExists ! Flag to indicate if the folder exists
LOGICAL, DIMENSION( 4 ) :: FluidPhase, PhaseTest ! Phase type (1: Supercritical fluid; 2: Critical fluid; 3: Subcritical liquid phase; 4: Subcritical vapor phase)

! Initialization
nSimplexCycles = 1
LastObjectiveFunction = -1.D0
SecondLastObjectiveFunction = LastObjectiveFunction

! Initial guess for the optimization parameters (original values)
Guesses(1,1) = INIVAPOR
Guesses(1,2) = INILIQUID

! Upper bound (initial guesses)
UpperBoundOptimization(1) = 0.9999D0
UpperBoundOptimization(2) = 0.9999D0

! Lower bound (initial guesses)
LowerBoundOptimization(1) = 0.1D-5
LowerBoundOptimization(2) = 0.1D-5

! Global simplex cycles
DO iSimplex = 1, cGlobalSimplex

  ! Other guesses
  CALL Random_Seed(  )
  DO iParameter = 2, nParameters + 1
    DO jParameter = 1, nParameters
      ! Random number
      CALL Random_Number( RandomNumber )
      ! Variation
      Delta = MinDelta + RandomNumber * ( MaxDelta - MinDelta )
      ! Random number
      CALL Random_Number( RandomNumber )
      ! Change the sign of the variation
      IF( RandomNumber >= 0.D0 .AND. RandomNumber < 0.5D0 ) THEN
        Delta = - Delta
      END IF
      ! Second guess for the optimization parameters (original values towards the positive direction)
      Guesses(iParameter,jParameter) = Guesses(1,jParameter) + Delta * ( UpperBoundOptimization(jParameter) - &
      &                                LowerBoundOptimization(jParameter) )
      IF( Guesses(iParameter,jParameter) < LowerBoundOptimization(jParameter) ) Guesses(iParameter,jParameter) = &
      &   LowerBoundOptimization(jParameter)
      IF( Guesses(iParameter,jParameter) > UpperBoundOptimization(jParameter) ) Guesses(iParameter,jParameter) = &
      &   UpperBoundOptimization(jParameter)
    END DO
  END DO

  ! Status
  !WRITE( *, "(5G0)" ) "Initializing Nelder-Mead Simplex Nonlinear Optimization... (", iSimplex, "/", cGlobalSimplex, ")"
  !WRITE( *, "(G0)" ) " "

  ! Objective function (initial guesses)
  DO iGuess = 1, nParameters + 1
    ! Initialization
    FluidPhase    = .FALSE.
    FluidPhase(4) = .TRUE. ! Subcritical vapor phase
    PhaseTest     = FluidPhase
    ! Calculate the molar volume of the vapor phase
    VaporComposition(1) = Guesses(iGuess,1)
    VaporComposition(2) = 1.D0 - VaporComposition(1)
    CALL Topliss_Algorithm( 1_Int64, VaporComposition, Temperature, Pressure, mVolumeVaporPhase, FluidPhase, CurveType, &
    &                       .FALSE. )
    IF( ALL( PhaseTest ) .NEQV. ALL( FluidPhase ) ) THEN ! Should not fail
      WRITE( *, "(G0)" ) "Error: 'PhaseTest' and 'FluidPhase' logical arrays are not equal. Exiting..."
      CALL EXIT(  )
    END IF
    ! Calculate compressibility factor of the vapor phase
    CALL Calculate_Pressure( VaporComposition, mVolumeVaporPhase, Temperature, Anynumber, IsothermalCompressibility, &
    &                        ThermalExpansionCoefficient, CompressibilityFactor )
    ! Calculate residual chemical potential of species i in the vapor phase
    CALL Calculate_CPotential( VaporComposition, mVolumeVaporPhase, Temperature, TotalChemicalPotential, &
    &                          ResidualChemicalPotential )
    vappot = TotalChemicalPotential
    ! Initialization
    FluidPhase    = .FALSE.
    FluidPhase(3) = .TRUE. ! Subcritical liquid phase
    PhaseTest     = FluidPhase
    ! Calculate the molar volume of the liquid phase
    LiquidComposition(1) = Guesses(iGuess,2)
    LiquidComposition(2) = 1.D0 - LiquidComposition(1)
    CALL Topliss_Algorithm( 1_Int64, LiquidComposition, Temperature, Pressure, mVolumeLiquidPhase, FluidPhase, CurveType, &
    &                       .FALSE. )
    IF( ALL( PhaseTest ) .NEQV. ALL( FluidPhase ) ) THEN ! Should not fail
      WRITE( *, "(G0)" ) "Error: 'PhaseTest' and 'FluidPhase' logical arrays are not equal. Exiting..."
      CALL EXIT(  )
    END IF
    ! Calculate compressibility factor of the liquid phase
    CALL Calculate_Pressure( LiquidComposition, mVolumeLiquidPhase, Temperature, Anynumber, IsothermalCompressibility, &
    &                        ThermalExpansionCoefficient, CompressibilityFactor )
    ! Calculate residual chemical potential of species i in the liquid phase
    CALL Calculate_CPotential( LiquidComposition, mVolumeLiquidPhase, Temperature, TotalChemicalPotential, &
    &                          ResidualChemicalPotential )
    liqpot = TotalChemicalPotential
    ! Objective function
    ObjectiveFunction(iGuess) = sqrt((vappot(1)-liqpot(1))*(vappot(1)-liqpot(1))+(vappot(2)-liqpot(2))*(vappot(2)-liqpot(2)))
  END DO

  ! Status
  !WRITE( *, "(2G0)" ) CHAR(13), "Calculating Initial Guesses... Done!"
  !WRITE( *, "(G0)" ) " "

  ! Sort the objective function values and the parameter sets
  CALL Bubble_Sort( nParameters + 1, ObjectiveFunction, Guesses )

  ! Parameter sets
  ParameterSets = Guesses

  ! Objective function (average)
  ObjectiveFunctionAverage = SUM( ObjectiveFunction ) / ( nParameters + 1 )

  ! Criteria for convergence
  ConvergenceCriterion = 0
  DO iParameter = 1, nParameters + 1
    ConvergenceCriterion = ConvergenceCriterion + ( ObjectiveFunction(iParameter) - ObjectiveFunctionAverage ) * &
    &                      ( ObjectiveFunction(iParameter) - ObjectiveFunctionAverage )
  END DO
  ConvergenceCriterion = ConvergenceCriterion / nParameters
  ConvergenceCriterion = DSQRT( ConvergenceCriterion )

  ! Simplex algorithm
  DO WHILE( ConvergenceCriterion >= Tolerance )

    ! Midpoint (centroid) between sets (excluding the worst set)
    DO iParameter = 1, nParameters
      AuxSum = SUM( ParameterSets(1:nParameters,iParameter), Dim= 1 )
      MidpointSet(iParameter) = AuxSum / nParameters
      IF( MidpointSet(iParameter) < LowerBoundOptimization(iParameter) ) MidpointSet(iParameter) = &
        & LowerBoundOptimization(iParameter)
      IF( MidpointSet(iParameter) > UpperBoundOptimization(iParameter) ) MidpointSet(iParameter) = &
        & UpperBoundOptimization(iParameter)
    END DO

    ! Initialization
    FluidPhase    = .FALSE.
    FluidPhase(4) = .TRUE. ! Subcritical vapor phase
    PhaseTest     = FluidPhase
    ! Calculate the molar volume of the vapor phase
    VaporComposition(1) = MidpointSet(1)
    VaporComposition(2) = 1.D0 - VaporComposition(1)    
    CALL Topliss_Algorithm( 1_Int64, VaporComposition, Temperature, Pressure, mVolumeVaporPhase, FluidPhase, CurveType, &
    &                       .FALSE. )
    IF( ALL( PhaseTest ) .NEQV. ALL( FluidPhase ) ) THEN ! Should not fail
      WRITE( *, "(G0)" ) "Error: 'PhaseTest' and 'FluidPhase' logical arrays are not equal. Exiting..."
      CALL EXIT(  )
    END IF
    ! Calculate compressibility factor of the vapor phase
    CALL Calculate_Pressure( VaporComposition, mVolumeVaporPhase, Temperature, Anynumber, IsothermalCompressibility, &
    &                        ThermalExpansionCoefficient, CompressibilityFactor )
    ! Calculate residual chemical potential of species i in the vapor phase
    CALL Calculate_CPotential( VaporComposition, mVolumeVaporPhase, Temperature, TotalChemicalPotential, &
    &                          ResidualChemicalPotential )
    vappot = TotalChemicalPotential
    ! Initialization
    FluidPhase    = .FALSE.
    FluidPhase(3) = .TRUE. ! Subcritical liquid phase
    PhaseTest     = FluidPhase
    ! Calculate the molar volume of the liquid phase
    LiquidComposition(1) = MidpointSet(2)
    LiquidComposition(2) = 1.D0 - LiquidComposition(1)
    CALL Topliss_Algorithm( 1_Int64, LiquidComposition, Temperature, Pressure, mVolumeLiquidPhase, FluidPhase, CurveType, &
    &                       .FALSE. )
    IF( ALL( PhaseTest ) .NEQV. ALL( FluidPhase ) ) THEN ! Should not fail
      WRITE( *, "(G0)" ) "Error: 'PhaseTest' and 'FluidPhase' logical arrays are not equal. Exiting..."
      CALL EXIT(  )
    END IF
    ! Calculate compressibility factor of the liquid phase
    CALL Calculate_Pressure( LiquidComposition, mVolumeLiquidPhase, Temperature, Anynumber, IsothermalCompressibility, &
    &                        ThermalExpansionCoefficient, CompressibilityFactor )
    ! Calculate residual chemical potential of species i in the liquid phase
    CALL Calculate_CPotential( LiquidComposition, mVolumeLiquidPhase, Temperature, TotalChemicalPotential, &
    &                          ResidualChemicalPotential )
    liqpot = TotalChemicalPotential
    ! Objective function
    ObjectiveFunctionMidpoint = sqrt((vappot(1)-liqpot(1))*(vappot(1)-liqpot(1))+(vappot(2)-liqpot(2))*(vappot(2)-liqpot(2)))

    ! Reflection of the worst set
    ReflectionSet(:) = ( 1.D0 + Alpha ) * MidpointSet(:) - Alpha * ParameterSets(nParameters+1,:)
    DO iParameter = 1, nParameters
      IF( ReflectionSet(iParameter) < LowerBoundOptimization(iParameter) ) ReflectionSet(iParameter) = &
      &   LowerBoundOptimization(iParameter)
      IF( ReflectionSet(iParameter) > UpperBoundOptimization(iParameter) ) ReflectionSet(iParameter) = &
      &   UpperBoundOptimization(iParameter)
    END DO

    ! Initialization
    FluidPhase    = .FALSE.
    FluidPhase(4) = .TRUE. ! Subcritical vapor phase
    PhaseTest     = FluidPhase
    ! Calculate the molar volume of the vapor phase
    VaporComposition(1) = ReflectionSet(1)
    VaporComposition(2) = 1.D0 - VaporComposition(1)
    CALL Topliss_Algorithm( 1_Int64, VaporComposition, Temperature, Pressure, mVolumeVaporPhase, FluidPhase, CurveType, &
    &                       .FALSE. )
    IF( ALL( PhaseTest ) .NEQV. ALL( FluidPhase ) ) THEN ! Should not fail
      WRITE( *, "(G0)" ) "Error: 'PhaseTest' and 'FluidPhase' logical arrays are not equal. Exiting..."
      CALL EXIT(  )
    END IF
    ! Calculate compressibility factor of the vapor phase
    CALL Calculate_Pressure( VaporComposition, mVolumeVaporPhase, Temperature, Anynumber, IsothermalCompressibility, &
    &                        ThermalExpansionCoefficient, CompressibilityFactor )
    ! Calculate residual chemical potential of species i in the vapor phase
    CALL Calculate_CPotential( VaporComposition, mVolumeVaporPhase, Temperature, TotalChemicalPotential, &
    &                          ResidualChemicalPotential )
    vappot = TotalChemicalPotential
    ! Initialization
    FluidPhase    = .FALSE.
    FluidPhase(3) = .TRUE. ! Subcritical liquid phase
    PhaseTest     = FluidPhase
    ! Calculate the molar volume of the liquid phase
    LiquidComposition(1) = ReflectionSet(2)
    LiquidComposition(2) = 1.D0 - LiquidComposition(1)
    CALL Topliss_Algorithm( 1_Int64, LiquidComposition, Temperature, Pressure, mVolumeLiquidPhase, FluidPhase, CurveType, &
    &                       .FALSE. )
    IF( ALL( PhaseTest ) .NEQV. ALL( FluidPhase ) ) THEN ! Should not fail
      WRITE( *, "(G0)" ) "Error: 'PhaseTest' and 'FluidPhase' logical arrays are not equal. Exiting..."
      CALL EXIT(  )
    END IF
    ! Calculate compressibility factor of the liquid phase
    CALL Calculate_Pressure( LiquidComposition, mVolumeLiquidPhase, Temperature, Anynumber, IsothermalCompressibility, &
    &                        ThermalExpansionCoefficient, CompressibilityFactor )
    ! Calculate residual chemical potential of species i in the liquid phase
    CALL Calculate_CPotential( LiquidComposition, mVolumeLiquidPhase, Temperature, TotalChemicalPotential, &
    &                          ResidualChemicalPotential )
    liqpot = TotalChemicalPotential
    ! Objective function
    ObjectiveFunctionReflection = sqrt((vappot(1)-liqpot(1))*(vappot(1)-liqpot(1))+(vappot(2)-liqpot(2))*(vappot(2)-liqpot(2)))

    ! Case selection
    IF( ObjectiveFunctionReflection <  ObjectiveFunction(nParameters) ) CaseSelection = 1 ! Reflection or Expansion
    IF( ObjectiveFunctionReflection >= ObjectiveFunction(nParameters) ) CaseSelection = 2 ! Contraction or Reduction/Shrink

    ! Cases
    SELECT CASE( CaseSelection )

    CASE( 1 ) ! Reflection or Expansion

      ! Apply REFLECTION
      IF( ObjectiveFunction(1) < ObjectiveFunctionReflection ) THEN
        ParameterSets(nParameters+1,:) = ReflectionSet(:)
        ObjectiveFunction(nParameters+1) = ObjectiveFunctionReflection
      ! Attempt EXPANSION
      ELSE
        ! Expansion of the reflection set
        ExpansionSet(:) = Gamma * ReflectionSet(:) + ( 1.D0 - Gamma ) * MidpointSet(:)
        DO iParameter = 1, nParameters
          IF( ExpansionSet(iParameter) < LowerBoundOptimization(iParameter) ) ExpansionSet(iParameter) = &
          &   LowerBoundOptimization(iParameter)
          IF( ExpansionSet(iParameter) > UpperBoundOptimization(iParameter) ) ExpansionSet(iParameter) = &
          &   UpperBoundOptimization(iParameter)
        END DO
        ! Initialization
        FluidPhase    = .FALSE.
        FluidPhase(4) = .TRUE. ! Subcritical vapor phase
        PhaseTest     = FluidPhase
        ! Calculate the molar volume of the vapor phase
        VaporComposition(1) = ExpansionSet(1)
        VaporComposition(2) = 1.D0 - VaporComposition(1)
        CALL Topliss_Algorithm( 1_Int64, VaporComposition, Temperature, Pressure, mVolumeVaporPhase, FluidPhase, CurveType, &
        &                       .FALSE. )
        IF( ALL( PhaseTest ) .NEQV. ALL( FluidPhase ) ) THEN ! Should not fail
          WRITE( *, "(G0)" ) "Error: 'PhaseTest' and 'FluidPhase' logical arrays are not equal. Exiting..."
          CALL EXIT(  )
        END IF
        ! Calculate compressibility factor of the vapor phase
        CALL Calculate_Pressure( VaporComposition, mVolumeVaporPhase, Temperature, Anynumber, IsothermalCompressibility, &
        &                        ThermalExpansionCoefficient, CompressibilityFactor )
        ! Calculate residual chemical potential of species i in the vapor phase
        CALL Calculate_CPotential( VaporComposition, mVolumeVaporPhase, Temperature, TotalChemicalPotential, &
        &                          ResidualChemicalPotential )
        vappot = TotalChemicalPotential
        ! Initialization
        FluidPhase    = .FALSE.
        FluidPhase(3) = .TRUE. ! Subcritical liquid phase
        PhaseTest     = FluidPhase
        ! Calculate the molar volume of the liquid phase
        LiquidComposition(1) = ExpansionSet(2)
        LiquidComposition(2) = 1.D0 - LiquidComposition(1)
        CALL Topliss_Algorithm( 1_Int64, LiquidComposition, Temperature, Pressure, mVolumeLiquidPhase, FluidPhase, CurveType, &
        &                       .FALSE. )
        IF( ALL( PhaseTest ) .NEQV. ALL( FluidPhase ) ) THEN ! Should not fail
          WRITE( *, "(G0)" ) "Error: 'PhaseTest' and 'FluidPhase' logical arrays are not equal. Exiting..."
          CALL EXIT(  )
        END IF
        ! Calculate compressibility factor of the liquid phase
        CALL Calculate_Pressure( LiquidComposition, mVolumeLiquidPhase, Temperature, Anynumber, IsothermalCompressibility, &
        &                        ThermalExpansionCoefficient, CompressibilityFactor )
        ! Calculate residual chemical potential of species i in the liquid phase
        CALL Calculate_CPotential( LiquidComposition, mVolumeLiquidPhase, Temperature, TotalChemicalPotential, &
        &                          ResidualChemicalPotential )
        liqpot = TotalChemicalPotential
        ! Objective function
        ObjectiveFunctionExpansion = sqrt((vappot(1)-liqpot(1))*(vappot(1)-liqpot(1))+(vappot(2)-liqpot(2))*(vappot(2)-liqpot(2)))
        IF( ObjectiveFunctionExpansion < ObjectiveFunction(1) ) THEN ! Apply EXPANSION
          ParameterSets(nParameters+1,:) = ExpansionSet(:)
          ObjectiveFunction(nParameters+1) = ObjectiveFunctionExpansion
        ELSE ! Apply REFLECTION
          ParameterSets(nParameters+1,:) = ReflectionSet(:)
          ObjectiveFunction(nParameters+1) = ObjectiveFunctionReflection
        END IF
      END IF

    CASE( 2 ) ! Contraction or Reduction/Shrink

      ! Replace WORST SET with REFLECTION
      IF( ObjectiveFunctionReflection < ObjectiveFunction(nParameters+1) ) THEN
        ParameterSets(nParameters+1,:) = ReflectionSet(:)
        ObjectiveFunction(nParameters+1) = ObjectiveFunctionReflection
      END IF

      ! Attempt CONTRACTION
      ContractionSet(:) = ( 1.D0 - Beta ) * MidpointSet(:) + ParameterSets(nParameters+1,:) * Beta
      DO iParameter = 1, nParameters
        IF( ContractionSet(iParameter) < LowerBoundOptimization(iParameter) ) ContractionSet(iParameter) = &
        &   LowerBoundOptimization(iParameter)
        IF( ContractionSet(iParameter) > UpperBoundOptimization(iParameter) ) ContractionSet(iParameter) = &
        &   UpperBoundOptimization(iParameter)
      END DO
      ! Initialization
      FluidPhase    = .FALSE.
      FluidPhase(4) = .TRUE. ! Subcritical vapor phase
      PhaseTest     = FluidPhase
      ! Calculate the molar volume of the vapor phase
      VaporComposition(1) = ContractionSet(1)
      VaporComposition(2) = 1.D0 - VaporComposition(1)
      CALL Topliss_Algorithm( 1_Int64, VaporComposition, Temperature, Pressure, mVolumeVaporPhase, FluidPhase, CurveType, &
      &                       .FALSE. )
      IF( ALL( PhaseTest ) .NEQV. ALL( FluidPhase ) ) THEN ! Should not fail
        WRITE( *, "(G0)" ) "Error: 'PhaseTest' and 'FluidPhase' logical arrays are not equal. Exiting..."
        CALL EXIT(  )
      END IF
      ! Calculate compressibility factor of the vapor phase
      CALL Calculate_Pressure( VaporComposition, mVolumeVaporPhase, Temperature, Anynumber, IsothermalCompressibility, &
      &                        ThermalExpansionCoefficient, CompressibilityFactor )
      ! Calculate residual chemical potential of species i in the vapor phase
      CALL Calculate_CPotential( VaporComposition, mVolumeVaporPhase, Temperature, TotalChemicalPotential, &
      &                          ResidualChemicalPotential )
      vappot = TotalChemicalPotential
      ! Initialization
      FluidPhase    = .FALSE.
      FluidPhase(3) = .TRUE. ! Subcritical liquid phase
      PhaseTest     = FluidPhase
      ! Calculate the molar volume of the liquid phase
      LiquidComposition(1) = ContractionSet(2)
      LiquidComposition(2) = 1.D0 - LiquidComposition(1)
      CALL Topliss_Algorithm( 1_Int64, LiquidComposition, Temperature, Pressure, mVolumeLiquidPhase, FluidPhase, CurveType, &
      &                       .FALSE. )
      IF( ALL( PhaseTest ) .NEQV. ALL( FluidPhase ) ) THEN ! Should not fail
        WRITE( *, "(G0)" ) "Error: 'PhaseTest' and 'FluidPhase' logical arrays are not equal. Exiting..."
        CALL EXIT(  )
      END IF
      ! Calculate compressibility factor of the liquid phase
      CALL Calculate_Pressure( LiquidComposition, mVolumeLiquidPhase, Temperature, Anynumber, IsothermalCompressibility, &
      &                        ThermalExpansionCoefficient, CompressibilityFactor )
      ! Calculate residual chemical potential of species i in the liquid phase
      CALL Calculate_CPotential( LiquidComposition, mVolumeLiquidPhase, Temperature, TotalChemicalPotential, &
      &                          ResidualChemicalPotential )
      liqpot = TotalChemicalPotential
      ! Objective function
      ObjectiveFunctionContraction = sqrt((vappot(1)-liqpot(1))*(vappot(1)-liqpot(1))+(vappot(2)-liqpot(2))*(vappot(2)-liqpot(2)))
      ! Apply CONTRACTION
      IF( ObjectiveFunctionContraction < ObjectiveFunction(nParameters+1) ) THEN
        ParameterSets(nParameters+1,:) = ContractionSet(:)
        ObjectiveFunction(nParameters+1) = ObjectiveFunctionContraction
      ! Apply REDUCTION/SHRINK
      ELSE
        ! Reduction/Shrink of all sets except the best one
        DO jParameter = 2, nParameters + 1
          ReductionSet(:) = 0.5D0 * ( ParameterSets(1,:) + ParameterSets(jParameter,:) ) ! Midpoint between best and other sets
          DO iParameter = 1, nParameters
            IF( ReductionSet(iParameter) < LowerBoundOptimization(iParameter) ) ReductionSet(iParameter) = &
            &   LowerBoundOptimization(iParameter)
            IF( ReductionSet(iParameter) > UpperBoundOptimization(iParameter) ) ReductionSet(iParameter) = &
            &   UpperBoundOptimization(iParameter)
          END DO
          ! Initialization
          FluidPhase    = .FALSE.
          FluidPhase(4) = .TRUE. ! Subcritical vapor phase
          PhaseTest     = FluidPhase
          ! Calculate the molar volume of the vapor phase
          VaporComposition(1) = ReductionSet(1)
          VaporComposition(2) = 1.D0 - VaporComposition(1)
          CALL Topliss_Algorithm( 1_Int64, VaporComposition, Temperature, Pressure, mVolumeVaporPhase, FluidPhase, CurveType, &
          &                       .FALSE. )
          IF( ALL( PhaseTest ) .NEQV. ALL( FluidPhase ) ) THEN ! Should not fail
            WRITE( *, "(G0)" ) "Error: 'PhaseTest' and 'FluidPhase' logical arrays are not equal. Exiting..."
            CALL EXIT(  )
          END IF
          ! Calculate compressibility factor of the vapor phase
          CALL Calculate_Pressure( VaporComposition, mVolumeVaporPhase, Temperature, Anynumber, IsothermalCompressibility, &
          &                        ThermalExpansionCoefficient, CompressibilityFactor )
          ! Calculate residual chemical potential of species i in the vapor phase
          CALL Calculate_CPotential( VaporComposition, mVolumeVaporPhase, Temperature, TotalChemicalPotential, &
          &                          ResidualChemicalPotential )
          vappot = TotalChemicalPotential
          ! Initialization
          FluidPhase    = .FALSE.
          FluidPhase(3) = .TRUE. ! Subcritical liquid phase
          PhaseTest     = FluidPhase
          ! Calculate the molar volume of the liquid phase
          LiquidComposition(1) = ReductionSet(2)
          LiquidComposition(2) = 1.D0 - LiquidComposition(1)
          CALL Topliss_Algorithm( 1_Int64, LiquidComposition, Temperature, Pressure, mVolumeLiquidPhase, FluidPhase, CurveType, &
          &                       .FALSE. )
          IF( ALL( PhaseTest ) .NEQV. ALL( FluidPhase ) ) THEN ! Should not fail
            WRITE( *, "(G0)" ) "Error: 'PhaseTest' and 'FluidPhase' logical arrays are not equal. Exiting..."
            CALL EXIT(  )
          END IF
          ! Calculate compressibility factor of the liquid phase
          CALL Calculate_Pressure( LiquidComposition, mVolumeLiquidPhase, Temperature, Anynumber, IsothermalCompressibility, &
          &                        ThermalExpansionCoefficient, CompressibilityFactor )
          ! Calculate residual chemical potential of species i in the liquid phase
          CALL Calculate_CPotential( LiquidComposition, mVolumeLiquidPhase, Temperature, TotalChemicalPotential, &
          &                          ResidualChemicalPotential )
          liqpot = TotalChemicalPotential
          ! Objective function
          ObjectiveFunctionReduction = sqrt((vappot(1)-liqpot(1))*(vappot(1)-liqpot(1))+(vappot(2)-liqpot(2))*(vappot(2)-&
          &   liqpot(2)))
          ParameterSets(jParameter,:) = ReductionSet(:)
          ObjectiveFunction(jParameter) = ObjectiveFunctionReduction
        END DO
      END IF

    END SELECT

    ! Sort the objective function values and the parameter sets
    CALL Bubble_Sort( nParameters + 1, ObjectiveFunction, ParameterSets )

    ! Last (best) objective function
    LastObjectiveFunction = ObjectiveFunction(1)

    ! Objective function (average)
    ObjectiveFunctionAverage = SUM( ObjectiveFunction ) / ( nParameters + 1 )

    ! Criteria for convergence
    ConvergenceCriterion = 0
    DO iParameter = 1, nParameters + 1
      ConvergenceCriterion = ConvergenceCriterion + ( ObjectiveFunction(iParameter) - ObjectiveFunctionAverage ) * &
      &                      ( ObjectiveFunction(iParameter) - ObjectiveFunctionAverage )
    END DO
    ConvergenceCriterion = ConvergenceCriterion / nParameters
    ConvergenceCriterion = DSQRT( ConvergenceCriterion )

    ! Update counter
    nSimplexCycles = nSimplexCycles + 1

    ! Critical properties
    IF( SecondLastObjectiveFunction /= LastObjectiveFunction ) THEN
      SecondLastObjectiveFunction = LastObjectiveFunction
    END IF

  END DO

  ! Status
  !WRITE( *, "(G0)" ) " "
  !WRITE( *, "(G0)" ) " "

  ! Re-initialize guesses with best result from previous simplex
  Guesses = ParameterSets

  !write(*,*) Guesses(1,1), Guesses(1,2),LastObjectiveFunction

  inivapor = Guesses(1,1)
  iniliquid = Guesses(1,2)
  ERROR = ConvergenceCriterion

END DO

RETURN

END SUBROUTINE Simplex_test