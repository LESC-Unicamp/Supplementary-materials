! ************************************************************************************************ !
!                         DENSITY-FINDING PROCEDURE (TOPLISS'S ALGORITHM)                          !
! ************************************************************************************************ !
!        This subroutine is used to determine the density roots at a specified pressure and        !
!           specified temperature following the Topliss's Algorithm (See R. J. Topliss,            !
!            D. Dimitrelis, and J. M. Prausnitz, Computers & Chemical Engineering, 12,             !
!                              483-489, 1988., for more information).                              !
! ************************************************************************************************ !
!           Although the SAFT EoS is non-cubic, it behaves as a cubic EoS within certain           !
!                                         density regions.                                         !
! ************************************************************************************************ !
! The density roots can occur on three curve types (isotherms):                                    !
!                                                                                                  !
! # Type C Isotherm: only one root above the critical isotherm                                     !
!                    (curve has no inflection points and the derivative of pressure with respect   ! 
!                    to the density is always positive)                                            !
! # Type B Isotherm: only one root at the critical isotherm                                        !
!                    (curve has one inflection point and the derivative of pressure with respect   ! 
!                    to the density at the inflection point is positive)                           !
! # Type A Isotherm: three roots below the critical isotherm (binodal curve)                       !
!                    (curve has one inflection point and the derivative of pressure with respect   ! 
!                    to the density at the inflection point is negative)                           !
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
SUBROUTINE Topliss_Algorithm( cComponent, mFraction, Temperature, SpecifiedPressure, mVolumeRoot, FluidPhase, CurveType, &
&                             PureComponent )

! Uses two modules: global variables and substances
USE GlobalVariables
USE Substances

IMPLICIT NONE

! ************************************************************************************************ !
! INTEGER VARIABLES                                                                                !
! ************************************************************************************************ !
INTEGER( Kind= Int64 ) :: OptimizationOption ! Optimization option of the golden-section search (1 = inflexion point; 2 = local minimum; 3 = local maximum)
INTEGER( Kind= Int64 ) :: cComponent         ! Component index

! ************************************************************************************************ !
! REAL VARIABLES                                                                                   !
! ************************************************************************************************ !
REAL( Kind= Real64 )                           :: Temperature                              ! Temperature
REAL( Kind= Real64 )                           :: mVolumeRoot                              ! Molar volume (root)
REAL( Kind= Real64 )                           :: mVolumeInflection                        ! Molar volume at the inflection point
REAL( Kind= Real64 )                           :: LowerLimitingVolume, UpperLimitingVolume ! Molar volume limits
REAL( Kind= Real64 )                           :: LowerBoundVolume, UpperBoundVolume       ! Molar volume bounds
REAL( Kind= Real64 )                           :: Pressure                                 ! Pressure (calculated)
REAL( Kind= Real64 )                           :: SpecifiedPressure                        ! Pressure (specified)
REAL( Kind= Real64 )                           :: MinimumPressure, MaximumPressure         ! Pressure (critical points)
REAL( Kind= Real64 )                           :: MinimumVolume, MaximumVolume             ! Molar volume (critical points)
REAL( Kind= Real64 )                           :: IsothermalCompressibility                ! Isothermal compressibility
REAL( Kind= Real64 )                           :: ThermalExpansionCoefficient              ! Thermal expansion coefficient
REAL( Kind= Real64 )                           :: CompressibilityFactor                    ! Compressibility factor
REAL( Kind= Real64 )                           :: dPressure_dDensity                       ! First derivative of pressure with respect to the molar density
REAL( Kind= Real64 )                           :: d2Pressure_d2Density                     ! Second derivative of pressure with respect to the molar density
REAL( Kind= Real64 )                           :: dPressure_dDensityInflection             ! First derivative of pressure with respect to the molar density at the inflection point
REAL( Kind= Real64 )                           :: aInterval, bInterval                     ! Golden-section search interval
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: mFraction                                ! Molar fraction of a component

! ************************************************************************************************ !
! CHARACTER STRINGS                                                                                !
! ************************************************************************************************ !
CHARACTER( Len= 01 ) :: CurveType ! Isotherm types (A, B, or C)

! ************************************************************************************************ !
! LOGICAL VARIABLES                                                                                !
! ************************************************************************************************ !
LOGICAL                 :: PureComponent ! Checks whether the system is a pure component or a mixture
LOGICAL, DIMENSION( 4 ) :: FluidPhase    ! Phase type (1: Supercritical fluid; 2: Critical fluid; 3: Subcritical liquid phase; 4: Subcritical vapor phase)

! Molar volume bounds
IF( PureComponent ) THEN
  LowerLimitingVolume = cPi * ijDiameterSphereCubic(cComponent,cComponent) * cAvogadro / 6.D0 / 0.999D0
ELSE
  LowerLimitingVolume = cPi * cAvogadro * SUM( mFraction * cDiameterSphere * cDiameterSphere * cDiameterSphere ) / 6.D0 / 0.999D0
END IF

! Molar density bounds
LowerBoundVolume = LowerLimitingVolume / 4.D-1
UpperBoundVolume = 1.D10

! Upper bound density / Lower bound volume
IF( PureComponent ) THEN
  CALL Calculate_Pressure_Single_Component( cComponent, LowerBoundVolume, Temperature, Pressure, IsothermalCompressibility, &
  &                                         ThermalExpansionCoefficient, CompressibilityFactor )
ELSE
  CALL Calculate_Pressure( mFraction, LowerBoundVolume, Temperature, Pressure, IsothermalCompressibility, &
  &                        ThermalExpansionCoefficient, CompressibilityFactor )
END IF
dPressure_dDensity = LowerBoundVolume / IsothermalCompressibility
! Conditions
DO WHILE( ( Pressure <= SpecifiedPressure ) .OR. ( dPressure_dDensity <= 0.D0 ) )
  LowerBoundVolume = LowerBoundVolume + ( 5.D-5 * ( LowerLimitingVolume - LowerBoundVolume ) )
  IF( PureComponent ) THEN
    CALL Calculate_Pressure_Single_Component( cComponent, LowerBoundVolume, Temperature, Pressure, IsothermalCompressibility, &
    &                                         ThermalExpansionCoefficient, CompressibilityFactor )
  ELSE
    CALL Calculate_Pressure( mFraction, LowerBoundVolume, Temperature, Pressure, IsothermalCompressibility, &
    &                        ThermalExpansionCoefficient, CompressibilityFactor )
  END IF
  dPressure_dDensity = LowerBoundVolume / IsothermalCompressibility
END DO

! Numerical approximation for the d²p/dρ² at the low-density limit (ρ = 0 / V → ∞)
UpperLimitingVolume = UpperBoundVolume
IF( PureComponent ) THEN
  CALL Calculate_Pressure_Single_Component( cComponent, UpperLimitingVolume, Temperature, Pressure, IsothermalCompressibility, &
  &                                         ThermalExpansionCoefficient, CompressibilityFactor )
ELSE
  CALL Calculate_Pressure( mFraction, UpperLimitingVolume, Temperature, Pressure, IsothermalCompressibility, &
  &                        ThermalExpansionCoefficient, CompressibilityFactor )
END IF
dPressure_dDensity = UpperLimitingVolume / IsothermalCompressibility
d2Pressure_d2Density = ( dPressure_dDensity - ( Pressure * UpperLimitingVolume ) ) * UpperLimitingVolume

! Type-C Isotherm (supercritical fluid)
IF( d2Pressure_d2Density > 0.D0 ) THEN
  ! Density root for the supercritical fluid
  CALL Brent_Method_Density_Root( cComponent, mFraction, Temperature, SpecifiedPressure, LowerBoundVolume, UpperBoundVolume, &
  &                               mVolumeRoot, PureComponent )
  ! Curve type
  CurveType = "C"
  ! Supercritical fluid
  FluidPhase    = .FALSE.
  FluidPhase(1) = .TRUE.
  RETURN
! Inflection Point
ELSE
  ! Inflection category
  OptimizationOption = 1
  aInterval = LowerBoundVolume
  bInterval = UpperLimitingVolume
  ! Molar volume at the inflection point
  IF( PureComponent ) THEN
    CALL Golden_Section_Single_Component( cComponent, Temperature, aInterval, bInterval, OptimizationOption, mVolumeInflection )
    CALL Calculate_Pressure_Single_Component( cComponent, mVolumeInflection, Temperature, Pressure, IsothermalCompressibility, &
    &                                         ThermalExpansionCoefficient, CompressibilityFactor )
  ELSE
    CALL Golden_Section_Search( mFraction, Temperature, aInterval, bInterval, OptimizationOption, mVolumeInflection )
    CALL Calculate_Pressure( mFraction, mVolumeInflection, Temperature, Pressure, IsothermalCompressibility, &
    &                        ThermalExpansionCoefficient, CompressibilityFactor )
  END IF
  ! First derivative of pressure with respect to the molar density at the inflection point
  dPressure_dDensityInflection = mVolumeInflection / IsothermalCompressibility
  ! Type-B Isotherm (critical fluid)
  IF( dPressure_dDensityInflection >= 0.D0 ) THEN
    ! Density root for the supercritical fluid
    CALL Brent_Method_Density_Root( cComponent, mFraction, Temperature, SpecifiedPressure, LowerBoundVolume, UpperBoundVolume, &
    &                               mVolumeRoot, PureComponent )
    ! Curve type
    CurveType = "B"
    ! Critical fluid
    FluidPhase    = .FALSE.
    FluidPhase(2) = .TRUE.
    RETURN
  ! Type-A Isotherm (subcritical fluid)
  ELSE
    ! Local minimum in P(ρ)
    IF( FluidPhase(3) ) THEN
      ! Local minimum category
      OptimizationOption = 2
      aInterval = LowerBoundVolume
      bInterval = mVolumeInflection
      ! Pressure at local minimum
      IF( PureComponent ) THEN
        CALL Golden_Section_Single_Component( cComponent, Temperature, aInterval, bInterval, OptimizationOption, MinimumVolume )
        CALL Calculate_Pressure_Single_Component( cComponent, MinimumVolume, Temperature, Pressure, IsothermalCompressibility, &
        &                                         ThermalExpansionCoefficient, CompressibilityFactor )
      ELSE
        CALL Golden_Section_Search( mFraction, Temperature, aInterval, bInterval, OptimizationOption, MinimumVolume )
        CALL Calculate_Pressure( mFraction, MinimumVolume, Temperature, Pressure, IsothermalCompressibility, &
        &                        ThermalExpansionCoefficient, CompressibilityFactor )
      END IF
      ! Pressure of the local extremum (minimum)
      MinimumPressure = Pressure
      ! Density root for the liquid phase
      IF( MinimumPressure < SpecifiedPressure ) THEN
        CALL Brent_Method_Density_Root( cComponent, mFraction, Temperature, SpecifiedPressure, LowerBoundVolume, MinimumVolume, &
        &                               mVolumeRoot, PureComponent )
        ! Curve type
        CurveType = "A"
        ! Subcritical liquid
        FluidPhase    = .FALSE.
        FluidPhase(3) = .TRUE.
        RETURN
      ELSE
        ! Curve type
        CurveType = "A"
        ! Subcritical vapor
        FluidPhase    = .FALSE.
        FluidPhase(4) = .TRUE.
        RETURN
      END IF
    ! Local maximum in P(ρ)
    ELSE IF( FluidPhase(4) ) THEN
      ! Local maximum category
      OptimizationOption = 3
      aInterval = mVolumeInflection
      bInterval = UpperBoundVolume
      ! Pressure at local maximum
      IF( PureComponent ) THEN
        CALL Golden_Section_Single_Component( cComponent, Temperature, aInterval, bInterval, OptimizationOption, MaximumVolume )
        CALL Calculate_Pressure_Single_Component( cComponent, MaximumVolume, Temperature, Pressure, IsothermalCompressibility, &
        &                                         ThermalExpansionCoefficient, CompressibilityFactor )
      ELSE
        CALL Golden_Section_Search( mFraction, Temperature, aInterval, bInterval, OptimizationOption, MaximumVolume )
        CALL Calculate_Pressure( mFraction, MaximumVolume, Temperature, Pressure, IsothermalCompressibility, &
        &                        ThermalExpansionCoefficient, CompressibilityFactor )
      END IF
      ! Pressure of the local extremum (maximum)
      MaximumPressure = Pressure
      ! Density root for the vapor phase
      IF( MaximumPressure > SpecifiedPressure ) THEN
        CALL Brent_Method_Density_Root( cComponent, mFraction, Temperature, SpecifiedPressure, MaximumVolume, UpperBoundVolume, &
        &                               mVolumeRoot, PureComponent )
        ! Curve type
        CurveType = "A"
        ! Subcritical vapor
        FluidPhase    = .FALSE.
        FluidPhase(4) = .TRUE.
        RETURN
      ELSE
        ! Curve type
        CurveType = "A"
        ! Subcritical liquid
        FluidPhase    = .FALSE.
        FluidPhase(3) = .TRUE.
        RETURN
      END IF
    END IF
  END IF
END IF

RETURN

END SUBROUTINE Topliss_Algorithm