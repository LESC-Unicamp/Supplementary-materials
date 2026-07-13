! ************************************************************************************* !
!                           OBJECTIVE FUNCTION (OPTIMIZATION)                           !
! ************************************************************************************* !
!     This subroutine corresponds to the objective function to be minimized by the      !
!      simplex algorithm, which is the average absolute relative deviation (AARD)       !
!                       between calculated and experimental data.                       !
! ************************************************************************************* !
!     => AUTHOR:     Nathan Barros de Souza                                             !
!     => E-MAIL:     n264179@dac.unicamp.br                                             !
!     => SUPERVISOR: Luís Fernando Mercier Franco                                       !
!************************************************************************************** !
!                 University of Campinas, Campinas - São Paulo, Brazil                  !
!                                   August 29th, 2022                                   !
!                                         v1.0                                          !
! ************************************************************************************* !
SUBROUTINE Objective_Function_Pure( nParameters, nPoints, DampingFactorIteration, cCriticalTemperature, cCriticalDensity, &
&                                   cCriticalPressure, cCriticalVolume, EqLiquidDens, EqVaporPressure, EqTemperature, Set, &
&                                   ObjectiveFunction, SimplexType, iGuess, ConvergenceCriterion, cPropCheck )

! Uses two modules: global variables and substances
USE GlobalVariables
USE Substances

IMPLICIT NONE

! ************************************************************************************************ !
! INTEGER VARIABLES                                                                                !
! ************************************************************************************************ !
INTEGER( Kind= Int64 ) :: cComponent      ! Counter (components)
INTEGER( Kind= Int64 ) :: Counter         ! Counter (generic)
INTEGER( Kind= Int64 ) :: dCounter        ! Counter (damping)
INTEGER( Kind= Int64 ) :: iPoint, nPoints ! Number of points in file
INTEGER( Kind= Int64 ) :: nParameters     ! Number of optimization parameters
INTEGER( Kind= Int64 ) :: IgnoredPoints   ! Number of ignored data points
INTEGER( Kind= Int64 ) :: SimplexType     ! Simplex type (-1: Initial guesses; 0: Midpoint; 1: Reflection; 2: Expansion; 3: Contraction; 4: Reduction)
INTEGER( Kind= Int64 ) :: iGuess          ! Guess index

! ************************************************************************************************ !
! REAL PARAMETERS                                                                                  !
! ************************************************************************************************ !
REAL( Kind= Real64 ), PARAMETER :: Tolerance = 1.D-8 ! Tolerance

! ************************************************************************************************ !
! REAL VARIABLES                                                                                   !
! ************************************************************************************************ !
REAL( Kind= Real64 )                           :: CalcLiquidDens           ! Saturated liquid density (calculated)
REAL( Kind= Real64 )                           :: DampingFactor            ! Damping factor
REAL( Kind= Real64 )                           :: DampingFactorIteration   ! Damping factor (iteration)
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
REAL( Kind= Real64 )                           :: MinimumPressure          ! Pressure (minimum)
REAL( Kind= Real64 )                           :: MaximumPressure          ! Pressure (maximum)
REAL( Kind= Real64 )                           :: MinimumDensity           ! Molar density (minimum)
REAL( Kind= Real64 )                           :: MaximumDensity           ! Molar density (maximum)
REAL( Kind= Real64 )                           :: ObjectiveFunction        ! Objective function
REAL( Kind= Real64 )                           :: ConvergenceCriterion     ! Convergence criterion
REAL( Kind= Real64 ), DIMENSION( nPoints )     :: EqLiquidDens             ! Saturated liquid density (experimental)
REAL( Kind= Real64 ), DIMENSION( nPoints )     :: EqVaporPressure          ! Vapor pressure (experimental)
REAL( Kind= Real64 ), DIMENSION( nPoints )     :: EqTemperature            ! Temperature (experimental)
REAL( Kind= Real64 ), DIMENSION( nParameters ) :: Set                      ! Set of parameters
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: mFraction                ! Molar fraction of a component
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: cCriticalVolume          ! Molar volume (critical)
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: cCriticalDensity         ! Molar density (critical)
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: cCriticalPressure        ! Pressure (critical)
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: cCriticalTemperature     ! Temperature (critical)

! ************************************************************************************************ !
! CHARACTER STRINGS                                                                                !
! ************************************************************************************************ !
CHARACTER( Len= 01 ) :: CurveType ! Isotherm types (A, B, or C)

! ************************************************************************************************ !
! LOGICAL VARIABLES                                                                                !
! ************************************************************************************************ !
LOGICAL                 :: cPropCheck ! Check critical properties only
LOGICAL, DIMENSION( 4 ) :: FluidPhase ! Phase type (1: Supercritical fluid; 2: Critical fluid; 3: Subcritical liquid phase; 4: Subcritical vapor phase)

! Initialization
cComponent = 1
ObjectiveFunction = 0.D0

! Parameters
cWellDepth(cComponent)      = Set(1) ! Joule
cPotentialRange(cComponent) = Set(2)
cDiameter(cComponent)       = Set(3) ! Meter
cAspectRatio(cComponent)    = Set(4)
cLength(cComponent)         = Set(4) * Set(3) ! Meter

! Diameter of the equivalent sphere (recalculated)
IF( GeometrySpecification(cComponent,1) ) THEN ! Ellipsoids-of-revolution
  cDiameterSphere(cComponent) = cDiameter(cComponent) * ( cAspectRatio(cComponent) ) ** ( 1.D0 / 3.D0 ) ! Meter
ELSE IF( GeometrySpecification(cComponent,2) ) THEN ! Spherocylinders
  cDiameterSphere(cComponent) = cDiameter(cComponent) * ( 1.D0 + 1.5D0 * cAspectRatio(cComponent) ) ** ( 1.D0 / 3.D0 ) ! Meter
ELSE IF( GeometrySpecification(cComponent,3) ) THEN ! Cylinders
  cDiameterSphere(cComponent) = cDiameter(cComponent) * ( 1.5D0 * cAspectRatio(cComponent) ) ** ( 1.D0 / 3.D0 ) ! Meter
END IF

! Field properties (recalculated)
fDiameter(cComponent) = cDiameter(cComponent) + cPotentialRange(cComponent) * cDiameterSphere(cComponent)
IF( GeometrySpecification(cComponent,1) ) THEN ! Ellipsoids-of-revolution
  fLength(cComponent) = cLength(cComponent) + cPotentialRange(cComponent) * cDiameterSphere(cComponent)
ELSE IF( GeometrySpecification(cComponent,2) ) THEN ! Spherocylinders
  fLength(cComponent) = cLength(cComponent)
ELSE IF( GeometrySpecification(cComponent,3) ) THEN ! Cylinders
  fLength(cComponent) = cLength(cComponent)
END IF
fAspectRatio(cComponent) = fLength(cComponent) / fDiameter(cComponent)

! Recalculate combining rules
CALL Combining_Rules(  )

! Recalculate combined properties
CALL Combined_Properties(  )

! Recalculate Isihara-Hadwiger theorem
CALL Isihara_Hadwiger_Theorem(  )

! Recalculate critical properties
CALL Critical_Properties_Pure_Components( cComponent, cCriticalTemperature(cComponent), cCriticalDensity(cComponent), &
&                                         cCriticalPressure(cComponent), cCriticalVolume(cComponent) )

! Check critical properties only
IF( cPropCheck ) RETURN

! Initialize
IgnoredPoints = 0
mFraction = 1.D0

DataLoop: DO iPoint = 1, nPoints
  ! Check temperature and ignore points above the critical temperature
  IF( EqTemperature(iPoint) > cCriticalTemperature(cComponent) ) THEN
    IgnoredPoints = IgnoredPoints + 1
    CYCLE DataLoop
  END IF
  ! Initialization
  DampingFactor = 1.D0
  Counter       = 0
  dCounter      = 0
  ! Find pressure interval
  CALL Find_Pressure_Interval( cComponent, mFraction, EqTemperature(iPoint), MinimumPressure, MinimumDensity, MaximumPressure, &
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
  CALL Topliss_Algorithm( 1_Int64, mFraction, EqTemperature(iPoint), MidpointPressure, LiquidVolume, FluidPhase, CurveType, .TRUE. )
  ! Find chemical potential of the liquid phase
  CALL Calculate_CPotential_Single_Component( 1_Int64, LiquidVolume, EqTemperature(iPoint), LiquidChemicalPotential, &
  &                                           rLiquidChemicalPotential )
  ! Vapor Phase
  FluidPhase(:) = .FALSE.
  FluidPhase(4) = .TRUE.
  ! Find volume of the vapor phase
  CALL Topliss_Algorithm( 1_Int64, mFraction, EqTemperature(iPoint), MidpointPressure, VaporVolume, FluidPhase, CurveType, .TRUE. )
  ! Find chemical potential of the liquid phase
  CALL Calculate_CPotential_Single_Component( 1_Int64, VaporVolume, EqTemperature(iPoint), VaporChemicalPotential, &
  &                                           rVaporChemicalPotential )
  ! Fugacity ratio
  FugacityRatio = DEXP( ( LiquidChemicalPotential - VaporChemicalPotential ) / ( cUniversalGas * EqTemperature(iPoint) ) )
  ! Convergence criterion
  FugacityError = FugacityRatio - 1.D0
  ! Progress
  CALL Progress_Bar_Path_03( iPoint, nPoints, iGuess, SimplexType, LastObjectiveFunction, ConvergenceCriterion )
  ! Convergence loop
  DO WHILE( DABS( FugacityError ) >= Tolerance )
    ! Update midpoint pressure
    MidpointPressure = MidpointPressure * FugacityRatio
    ! Liquid Phase
    FluidPhase(:) = .FALSE.
    FluidPhase(3) = .TRUE.
    ! Find volume of the liquid phase
    CALL Topliss_Algorithm( 1_Int64, mFraction, EqTemperature(iPoint), MidpointPressure, LiquidVolume, FluidPhase, CurveType, &
    &                       .TRUE. )
    ! Find chemical potential of the liquid phase
    CALL Calculate_CPotential_Single_Component( 1_Int64, LiquidVolume, EqTemperature(iPoint), LiquidChemicalPotential, &
    &                                           rLiquidChemicalPotential )
    ! Vapor Phase
    FluidPhase(:) = .FALSE.
    FluidPhase(4) = .TRUE.
    ! Find volume of the vapor phase
    CALL Topliss_Algorithm( 1_Int64, mFraction, EqTemperature(iPoint), MidpointPressure, VaporVolume, FluidPhase, CurveType, &
    &                       .TRUE. )
    ! Find chemical potential of the liquid phase
    CALL Calculate_CPotential_Single_Component( 1_Int64, VaporVolume, EqTemperature(iPoint), VaporChemicalPotential, &
    &                                           rVaporChemicalPotential )
    ! Fugacity ratio
    FugacityRatio = DEXP( ( LiquidChemicalPotential - VaporChemicalPotential ) / ( cUniversalGas * EqTemperature(iPoint) ) ) ** &
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
      IgnoredPoints = IgnoredPoints + 1
      CYCLE DataLoop
    END IF
    ! Check iterations
    IF( Counter >= 2000 ) THEN
      IgnoredPoints = IgnoredPoints + 1
      CYCLE DataLoop
    END IF
  END DO

  ! Liquid density
  CalcLiquidDens = 1.D0 / LiquidVolume

  ! Objective function (minimization)
  ObjectiveFunction = ObjectiveFunction + DABS( ( CalcLiquidDens - EqLiquidDens(iPoint) ) / EqLiquidDens(iPoint) ) + &
  &                   DABS( ( ( MidpointPressure / 1.D6 ) - EqVaporPressure(iPoint) ) / EqVaporPressure(iPoint) )

END DO DataLoop

! Check number of data points
IF( (nPoints - IgnoredPoints) <= 0 ) THEN
  WRITE( *, "(G0)" ) " "
  WRITE( *, "(G0,G0.5,G0)" ) "Error: the number of data points is less than one. Finishing..."
  CALL EXIT(  )
END IF

! Objective function
ObjectiveFunction = ObjectiveFunction / DBLE( nPoints - IgnoredPoints ) * 100.D0

RETURN

END SUBROUTINE Objective_Function_Pure

SUBROUTINE Objective_Function_Mixture( nParameters, nPoints, Temperature, cAccentricFactor, cCriticalPressure, &
&                                      cCriticalTemperature, EqLiquidComposition, EqVaporComposition, EqPressure, &
&                                      Set, ObjectiveFunction, SimplexType, iGuess, ConvergenceCriterion )

! Uses two modules: global variables and substances
USE GlobalVariables
USE Substances

IMPLICIT NONE

! ************************************************************************************************ !
! INTEGER VARIABLES                                                                                !
! ************************************************************************************************ !
INTEGER( Kind= Int64 ) :: cComponent  ! Counter (components)
INTEGER( Kind= Int64 ) :: iPoint      ! Counter (points)
INTEGER( Kind= Int64 ) :: cCounter    ! Counter
INTEGER( Kind= Int64 ) :: nParameters ! Number of parameters
INTEGER( Kind= Int64 ) :: nPoints     ! Number of data points
INTEGER( Kind= Int64 ) :: SimplexType ! Simplex type (-1: Initial guesses; 0: Midpoint; 1: Reflection; 2: Expansion; 3: Contraction; 4: Reduction)
INTEGER( Kind= Int64 ) :: iGuess      ! Guess index

! ************************************************************************************************ !
! REAL VARIABLES                                                                                   !
! ************************************************************************************************ !
REAL( Kind= Real64 ) :: Temperature ! Temperature
REAL( Kind= Real64 ) :: Pressure                         ! Pressure
REAL( Kind= Real64 ) :: Error        ! Difference between old and new equilibrium factors
REAL( Kind= Real64 ) :: beta_min, beta_max           ! Dummy (number)
REAL( Kind= Real64 )            :: beta_initial           ! Dummy (number)
REAL( Kind= Real64 )            :: g                  ! Speed of sound
REAL( Kind= Real64 )            :: beta_frac           ! Dummy (number)
REAL( Kind= Real64 )                           :: ConvergenceCriterion     ! Convergence criterion
REAL( Kind= Real64 )                 :: ObjectiveFunction
REAL( Kind= Real64 ), DIMENSION(nPoints)  :: EqLiquidComposition, EqVaporComposition, EqPressure
REAL( Kind= Real64 ), DIMENSION(nParameters) :: Set
REAL( Kind= Real64 ), DIMENSION(2) :: FeedComposition, VaporComposition, LiquidComposition, test_beta_max, test_beta_min
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: cAccentricFactor
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: cCriticalPressure
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: cEquilibriumFactors
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: cCriticalTemperature ! Critical temperature

! ************************************************************************************************ !
! LOGICAL VARIABLES                                                                                !
! ************************************************************************************************ !
LOGICAL :: subcooled_l   ! Various uses
LOGICAL :: superheated_v ! Various uses
LOGICAL :: any_ki_1      ! Various uses

! Parameters
ijWellDepthCorrection = Set(1)
IF( nParameters == 2 ) ijPotentialRangeCorrection = Set(2)

! Recalculate combining rules
CALL Combining_Rules(  )

! Recalculate combined properties
ijaWellDepth = ijWellDepth / cBoltzmann ! K

! Combined potential range
IF( nParameters == 2 ) ijPotentialRangeSquared = ijPotentialRange * ijPotentialRange
IF( nParameters == 2 ) ijPotentialRangeCubic = ijPotentialRange * ijPotentialRange * ijPotentialRange

! ************************************************************************************* !
! Minimization function (initial value)                                                 !
! ************************************************************************************* !
ObjectiveFunction = 0.D0

DO iPoint = 1, nPoints

  Pressure = EqPressure(iPoint)
  FeedComposition(2) = EqVaporComposition(iPoint)
  FeedComposition(1) = 1.D0 - FeedComposition(2)

  cEquilibriumFactors(:) = DLOG( cCriticalPressure(:) / Pressure ) + 5.373D0 * ( 1.D0 + cAccentricFactor(:) ) * ( 1.D0 - &
  & cCriticalTemperature(:) / Temperature )
  cEquilibriumFactors = DEXP( cEquilibriumFactors )

  Error = 1.D0
  cCounter=0

  CALL Progress_Bar_Path_03( iPoint, nPoints, iGuess, SimplexType, LastObjectiveFunction, ConvergenceCriterion )

  do while ((dabs(Error) >= 1d-5) .and. (cCounter <= 40))
  
    ! Testing for a subcooled liquid
    call subcooled_liquid_test(Temperature,Pressure,FeedComposition,cEquilibriumFactors,VaporComposition,LiquidComposition,&
    & subcooled_l)
        
    if (.not. subcooled_l) then
  
        ! Testing for a superheated vapor
      call superheated_vapor_test(Temperature,Pressure,FeedComposition,cEquilibriumFactors,VaporComposition,LiquidComposition,&
      & superheated_v)
  
      if (.not. superheated_v) then
  
        ! Setting beta interval between 0 and 1
        beta_min = 0.d0
        beta_max = 1.d0
  
        ! Testing condition expressed by Equation (7),
        test_beta_min(:) = 0.d0
        do cComponent=1,nComponents
          if (cEquilibriumFactors(cComponent) > 1.d0) then
            test_beta_min(cComponent) = (cEquilibriumFactors(cComponent)*FeedComposition(cComponent)-1.d0) &
            &                  /(cEquilibriumFactors(cComponent)-1.d0)
            any_ki_1         = .true.
          end if
        end do
        if (any_ki_1) then
          beta_min = maxval(test_beta_min)
        end if
  
        any_ki_1 = .false.

        ! Testing condition expressed by Equation (8),
        test_beta_max(:) = 1.d0
        do cComponent=1,nComponents
          if (cEquilibriumFactors(cComponent) < 1.d0) then
            test_beta_max(cComponent) = (1.d0-FeedComposition(cComponent))       &
            &                  /(1.d0-cEquilibriumFactors(cComponent))
            any_ki_1         = .true.
          end if
        end do
        if (any_ki_1) then
          beta_max = minval(test_beta_max)
        end if
        any_ki_1 = .false.
  
        ! Computing initial estimate for beta
        ! (fraction of vapor phase)
        beta_initial = 0.5d0*(beta_min+beta_max)
  
        ! Rachford-Rice equation calculated with the initial
        ! estimate of beta (fraction of vapor phase)
        ! Equation (2), Chapter 10, p. 252
        g = 0.d0
        do cComponent=1,nComponents
          g = g+FeedComposition(cComponent)*(cEquilibriumFactors(cComponent)-1.d0)                 &
          &   /(1.d0-beta_initial+beta_initial*cEquilibriumFactors(cComponent))
        end do
  
        ! Updating boundary values of beta
        !beta_min = 0.d0
        !beta_max = 1.d0
        if (g > 0.d0) then
            beta_min = beta_initial
        else 
            beta_max = beta_initial
        end if
  
        beta_min = 0.0d0
        beta_max = 1.d0

        call Brent_Method_Rachford_Rice( beta_min, beta_max, FeedComposition, cEquilibriumFactors, beta_frac )
  
        ! Evaluating new liquid and vapor mole fractions
        LiquidComposition(:) = FeedComposition(:)/(1.d0-beta_frac+beta_frac*cEquilibriumFactors(:))
        VaporComposition(:)  = cEquilibriumFactors(:)*LiquidComposition(:)
  
      end if
  
    end if
  
    call Equilibrium_Factor_Calculation( Temperature, Pressure, VaporComposition, LiquidComposition, cEquilibriumFactors, Error )

    cCounter = cCounter+1

  end do

  ObjectiveFunction = ObjectiveFunction + DABS(VaporComposition(1)-EqVaporComposition(iPoint)) + &
  & DABS(LiquidComposition(1)-EqLiquidComposition(iPoint))

END DO

ObjectiveFunction = ObjectiveFunction / DBLE( nPoints ) * 100.D0

!PRINT *, ObjectiveFunction,ijWellDepthCorrection,ijPotentialRangeCorrection

RETURN

END SUBROUTINE Objective_Function_Mixture