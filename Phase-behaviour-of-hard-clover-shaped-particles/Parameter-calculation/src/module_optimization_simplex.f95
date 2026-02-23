MODULE SimplexAlgorithm

! Uses one module: global variables
USE GlobalVariables

CONTAINS

! ********************************************************************************************************************************* !
! Simplex algorithm for optimization problems.                                                                                      !
!                                                                                                                                   !
! ● Because the simplex involves a single parameter, two initial guesses are needed. The first guress is user-defined               !
!   (<InitialLayerSpacing>), while the second guess is randomly generated in the range [<MinLayerSpacing>, <MaxLayerSpacing>].      !
!                                                                                                                                   !
! ● After defining the two initial guesses, the algorithm calculates the objective function for each guess, as well as the          !
!   convergence criterion (<ConvergenceCriterion>).                                                                                 !
!                                                                                                                                   !
! ● The Nelder-Mead simplex algorithm is based on the following steps:                                                              !
!                                                                                                                                   !
!     1. Calculate the midpoint (centroid) between sets (excluding the worst set).                                                  !
!     2. Calculate the reflection of the worst set.                                                                                 !
!     3. Depending on the objective function value of the reflected set, apply:                                                     !
!        a. Reflection or Expansion, if the reflected set is better than the worst set.                                             !
!        b. Contraction or Reduction/Shrink, if the reflected set is worse than the worst set.                                      !
!     4. The expansion operation attempts to expand the reflection set, speeding up the convergence. If it fails,                   !
!        the reflection set is used.                                                                                                !
!     5. The contraction operation attempts to contract the reflection set, improving the value of the worst set. If it fails,      !
!        the reduction/shrink operation is applied to all sets except the best one.                                                 !
!     6. The algorithm sorts the objective function values and the parameter sets.                                                  !
!     7. Calculate the last (best) objective function and the objective function (average).                                         !
!     8. Re-calculate the convergence criterion.                                                                                    !
!     9. Loops until the convergence criterion is less than the tolerance.                                                          !
!                                                                                                                                   !
! ● OBS.: Due to the cyclial nature of the objective function, the algorithm may fall into a local minimum. To avoid this, the      !
!         algorithm is repeated a number of times (<cGlobalSimplex>). This however does not guarantee that the global minimum will  !
!         be found, but it increases the chances of finding a better solution. The post-processing of the results is essential to   !
!         verify the quality of the optimization.                                                                                   ! 
! ********************************************************************************************************************************* !
! For more information, please read: Nelder, J. A.; Mead, R. A., Comp. J., 7, 308-313, 1965                                         !
! ********************************************************************************************************************************* !
SUBROUTINE Nelder_Mead( nParameters, CurrentFrame, TotalFrames, PhaseDirector, pCurrentPosition, SmecticOrder, OptimalLayerSpacing )

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 )               :: iGuess                       ! Counter (guesses)
INTEGER( Kind= Int64 )               :: iSimplex                     ! Counter (simplex)
INTEGER( Kind= Int64 )               :: iParameter                   ! Counter (parameters)
INTEGER( Kind= Int64 )               :: jParameter                   ! Counter (parameters)
INTEGER( Kind= Int64 )               :: CaseSelection                ! Case selection (1: Reflection or Expansion; 2: Contraction or Reduction/Shrink)
INTEGER( Kind= Int64 )               :: nSimplexCycles               ! Number of simplex cycles
INTEGER( Kind= Int64 )               :: OldLenMessage, MaxLenMessage ! Length of the print message
INTEGER( Kind= Int64 ), INTENT( IN ) :: nParameters                  ! Number of optimization parameters
INTEGER( Kind= Int64 ), INTENT( IN ) :: CurrentFrame, TotalFrames    ! Frames of the configuration file

! REAL PARAMETERS (SCALAR,CONSTANTS)
REAL( Kind= Real64 ), PARAMETER :: Tolerance = 1.D-6 ! Tolerance
REAL( Kind= Real64 ), PARAMETER :: Alpha     = 1.D0  ! Reflection coefficient
REAL( Kind= Real64 ), PARAMETER :: Beta      = 0.5D0 ! Contraction coefficient
REAL( Kind= Real64 ), PARAMETER :: Gamma     = 2.D0  ! Expansion coefficient

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 )                :: RandomNumber                 ! Random number
REAL( Kind= Real64 )                :: ObjectiveFunctionMidpoint    ! Objective function (midpoint)
REAL( Kind= Real64 )                :: ObjectiveFunctionReflection  ! Objective function (reflection)
REAL( Kind= Real64 )                :: ObjectiveFunctionExpansion   ! Objective function (expansion)
REAL( Kind= Real64 )                :: ObjectiveFunctionContraction ! Objective function (contraction)
REAL( Kind= Real64 )                :: ObjectiveFunctionReduction   ! Objective function (reduction)
REAL( Kind= Real64 )                :: ObjectiveFunctionAverage     ! Objective function (average)
REAL( Kind= Real64 )                :: ConvergenceCriterion         ! Convergence criterion
REAL( Kind= Real64 )                :: LastObjectiveFunction        ! Last value of the objective function
REAL( Kind= Real64 )                :: SecondLastObjectiveFunction  ! Second last value of the objective function
REAL( Kind= Real64 )                :: AuxSum                       ! Auxiliary sum
REAL( Kind= Real64 ), INTENT( OUT ) :: SmecticOrder                 ! Smectic order parameter
REAL( Kind= Real64 ), INTENT( OUT ) :: OptimalLayerSpacing          ! Optmized layer spacing

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( nParameters )                :: LowerBoundOptimization ! Lower bound
REAL( Kind= Real64 ), DIMENSION( nParameters )                :: UpperBoundOptimization ! Upper bound
REAL( Kind= Real64 ), DIMENSION( nParameters )                :: MidpointSet            ! Set of parameters (midpoint)
REAL( Kind= Real64 ), DIMENSION( nParameters )                :: ReflectionSet          ! Set of parameters (reflection)
REAL( Kind= Real64 ), DIMENSION( nParameters )                :: ExpansionSet           ! Set of parameters (expansion)
REAL( Kind= Real64 ), DIMENSION( nParameters )                :: ContractionSet         ! Set of parameters (contraction)
REAL( Kind= Real64 ), DIMENSION( nParameters )                :: ReductionSet           ! Set of parameters (reduction)
REAL( Kind= Real64 ), DIMENSION( nParameters+1 )              :: ObjectiveFunction      ! Objective function
REAL( Kind= Real64 ), DIMENSION( nParameters+1, nParameters ) :: ParameterSets          ! Parameter sets
REAL( Kind= Real64 ), DIMENSION( nParameters+1, nParameters ) :: Guesses                ! Guesses
REAL( Kind= Real64 ), DIMENSION( : ), INTENT( IN )            :: PhaseDirector          ! Phase director of the nematic phase
REAL( Kind= Real64 ), DIMENSION( :, : ), INTENT( IN )         :: pCurrentPosition       ! Position of particles

! CHARACTER STRINGS
CHARACTER( Len= 1000 ) :: PrintStatus ! Print status

! Initialization
nSimplexCycles = 1
LastObjectiveFunction = -1.D0
SecondLastObjectiveFunction = LastObjectiveFunction

! Initial guess for the optimization parameter (original values)
DO iParameter = 1, nParameters
  Guesses(1,iParameter) = InitialLayerSpacing ! Layer spacing (smectic phase)
END DO

! Upper bound and lower bound (initial guesses)
DO iParameter = 1, nParameters
  UpperBoundOptimization(iParameter) = MaxLayerSpacing
  LowerBoundOptimization(iParameter) = MinLayerSpacing
END DO

! Initialization
MaxLenMessage = 0

! Global simplex cycles
DO iSimplex = 1, cGlobalSimplex

  ! Status
  WRITE( PrintStatus, "(9G0)" ) "Reading frame: ", CurrentFrame, " of ", TotalFrames, " | Nelder-Mead simplex repetition ", &
  &                             iSimplex, " of ", cGlobalSimplex, "..."
  OldLenMessage = LEN( TRIM( PrintStatus ) )
  IF( OldLenMessage > MaxLenMessage ) MaxLenMessage = OldLenMessage
  WRITE( *, "(2G0)", Advance= "No" ) CHAR(13), TRIM( PrintStatus )

  ! Other guesses
  CALL Random_Seed(  )
  DO iParameter = 2, nParameters + 1
    DO jParameter = 1, nParameters
      ! Random number [0, 1[
      CALL Random_Number( RandomNumber )
      ! Second guess for the optimization parameters (original values towards the positive direction)
      Guesses(iParameter,jParameter) = LowerBoundOptimization(jParameter) + RandomNumber * ( UpperBoundOptimization(jParameter) - &
      &                                LowerBoundOptimization(jParameter) )
    END DO
  END DO

  ! Objective function (initial guesses)
  DO iGuess = 1, nParameters + 1
    CALL Objective_Function_SmecticOrder( Guesses(iGuess,:), ObjectiveFunction(iGuess), PhaseDirector, pCurrentPosition, &
    &                                     SmecticOrder )
  END DO

  ! Sort the objective function values and the parameter sets
  CALL Bubble_Sort( nParameters + 1, ObjectiveFunction, Guesses )

  ! Parameter sets
  ParameterSets = Guesses

  ! Objective function (average)
  ObjectiveFunctionAverage = SUM( ObjectiveFunction ) / ( nParameters + 1 )

  ! Criteria for convergence
  ConvergenceCriterion = 0.D0
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
    CALL Objective_Function_SmecticOrder( MidpointSet, ObjectiveFunctionMidpoint, PhaseDirector, pCurrentPosition, SmecticOrder )

    ! Reflection of the worst set
    ReflectionSet(:) = ( 1.D0 + Alpha ) * MidpointSet(:) - Alpha * ParameterSets(nParameters+1,:)
    ! Objective function (reflection)
    CALL Objective_Function_SmecticOrder( ReflectionSet, ObjectiveFunctionReflection, PhaseDirector, pCurrentPosition, &
    &                                     SmecticOrder )

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
        CALL Objective_Function_SmecticOrder( ExpansionSet, ObjectiveFunctionExpansion, PhaseDirector, pCurrentPosition, &
        &                                     SmecticOrder )
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
      CALL Objective_Function_SmecticOrder( ContractionSet, ObjectiveFunctionContraction, PhaseDirector, pCurrentPosition, &
      &                                     SmecticOrder )

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
          CALL Objective_Function_SmecticOrder( ReductionSet, ObjectiveFunctionReduction, PhaseDirector, pCurrentPosition, &
          &                                     SmecticOrder )
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

  END DO

  ! Re-initialize guesses with best result from previous simplex
  DO iParameter = 1, nParameters
    Guesses(1,iParameter) = ParameterSets(1,iParameter)
  END DO

END DO

! Status
WRITE( PrintStatus, "(A,I0,A,I0)" ) "Reading frame: ", CurrentFrame, " of ", TotalFrames
OldLenMessage = LEN( TRIM( PrintStatus ) )
WRITE( *, "(6G0)", Advance= "No" ) CHAR(13), "Reading frame: ", CurrentFrame, " of ", TotalFrames, &
&                                  REPEAT( " ", MaxLenMessage - OldLenMessage )

! Optimal layer spacing
OptimalLayerSpacing = DABS( ParameterSets(1,1) )

! Smectic order parameter
SmecticOrder = 1.D0 - LastObjectiveFunction

RETURN

END SUBROUTINE Nelder_Mead

! ********************************************************************************************************************************* !
! Objective function for the smectic order parameter.                                                                               !
!                                                                                                                                   !
! ● <Set> refers to the layer spacing of the smectic phase that is being optimized to maximize the smectic order parameter. In this !
!   case, since the maximum value of the smectic parameter is 1, we define the function 1 - <SmecticOrder> and try to minimize it   !
!   instead.                                                                                                                        !
!                                                                                                                                   !
! ● The smectic order parameter is defined as:                                                                                      !
!                                                                                                                                   !
!     τ = | Σⱼ exp[ 2π × i × (n · rⱼ) / d ] | / N                                                                                   !
!                                                                                                                                   !
!   where n is the phase director, rⱼ is the position of the particles, d is the layer spacing, and N is the number of particles.   !
!                                                                                                                                   !
! ● The term | … | can be represented as √| … |², which is the absolute squared magnitude of a complex number, given as:            !
!                                                                                                                                   !
!     √| a + i × b |² = √[ ( a + i × b ) × ( a - i × b ) ] = √(a² + b²).                                                            !
! ********************************************************************************************************************************* !
SUBROUTINE Objective_Function_SmecticOrder( Set, ObjectiveFunction, PhaseDirector, pCurrentPosition, SmecticOrder )

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 ) :: iParticle = 0 ! Counter (particle)

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 )                :: TrigonometricArgument ! Argument of the trigonometric functions
REAL( Kind= Real64 )                :: Cossine               ! Cossine term
REAL( Kind= Real64 )                :: Sine                  ! Sine term
REAL( Kind= Real64 )                :: CossineTerm           ! Sum of cossine terms
REAL( Kind= Real64 )                :: SineTerm              ! Sum of sine terms
REAL( Kind= Real64 ), INTENT( OUT ) :: ObjectiveFunction     ! Objective function
REAL( Kind= Real64 ), INTENT( OUT ) :: SmecticOrder          ! Smectic order parameter

! REAL PARAMETERS (SCALAR,CONSTANTS)
REAL( Kind= Real64 ), PARAMETER :: Tolerance = 1.D-8 ! Tolerance

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( : ), INTENT( IN )    :: Set              ! Set of parameters
REAL( Kind= Real64 ), DIMENSION( : ), INTENT( IN )    :: PhaseDirector    ! Phase director of the nematic phase
REAL( Kind= Real64 ), DIMENSION( :, : ), INTENT( IN ) :: pCurrentPosition ! Position of particles

! Initialization
ObjectiveFunction = 0.D0
CossineTerm = 0.D0
SineTerm = 0.D0

! Loop over particles
DO iParticle = 1, nParticles
  TrigonometricArgument = 2.D0 * cPi * DOT_PRODUCT( PhaseDirector, pCurrentPosition(iParticle,1:3) ) / Set(1) ! Smectic ordering
  Cossine = DCOS( TrigonometricArgument ) ! Cossine term
  Sine = DSIN( TrigonometricArgument )  ! Sine term
  CossineTerm = CossineTerm + Cossine ! Sum of cossine terms
  SineTerm = SineTerm + Sine ! Sum of sine terms
END DO

! Smectic order parameter
SmecticOrder = DSQRT( CossineTerm * CossineTerm + SineTerm * SineTerm ) / DBLE( nParticles ) ! Square root of the squared absolute magnitude of the complex number

! Objective function (minimization)
ObjectiveFunction = 1.D0 - SmecticOrder

RETURN

END SUBROUTINE Objective_Function_SmecticOrder

! ********************************************************************************************************************************* !
! This subroutine sorts an array in ascending order using the bubble sort.                                                          !
! ********************************************************************************************************************************* !
SUBROUTINE Bubble_Sort( ArraySize, Array, ArrayAux )

IMPLICIT NONE

! INTEGER VARIABLES (SCALAR)
INTEGER( Kind= Int64 )               :: iArray, jArray ! Counters
INTEGER( Kind= Int64 ), INTENT( IN ) :: ArraySize      ! Array size

! REAL VARIABLES (SCALAR)
REAL( Kind= Real64 ) :: ArrayTemp ! Temporary array

! REAL VARIABLES (ARRAY)
REAL( Kind= Real64 ), DIMENSION( ArraySize )              :: Array        ! Array
REAL( Kind= Real64 ), DIMENSION( ArraySize-1 )            :: ArrayAuxTemp ! Temporary array (auxiliary)
REAL( Kind= Real64 ), DIMENSION( ArraySize, ArraySize-1 ) :: ArrayAux     ! Array (auxiliary)

! Sort array
DO iArray = 1, ArraySize-1
  DO jArray = 1, ArraySize-iArray
    IF (Array(jArray) > Array(jArray+1)) THEN
      ArrayTemp = Array(jArray)
      ArrayAuxTemp(:) = ArrayAux(jArray,:)
      Array(jArray) = Array(jArray+1)
      Array(jArray+1) = ArrayTemp
      ArrayAux(jArray,:) = ArrayAux(jArray+1,:)
      ArrayAux(jArray+1,:) = ArrayAuxTemp(:)
    END IF
  END DO
END DO

RETURN

END SUBROUTINE Bubble_Sort

END MODULE SimplexAlgorithm