! ************************************************************************************************ !
!                                      FIND PRESSURE INTERVAL                                      !
! ************************************************************************************************ !
!       This subroutine is used to find the maximum and minimum pressures (critical points)        !
!                                      of a type-A isotherm.                                       !
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
SUBROUTINE Find_Pressure_Interval( cComponent, mFraction, Temperature, MinimumPressure, MinimumDensity, MaximumPressure, &
&                                  MaximumDensity, PureComponent )

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
REAL( Kind= Real64 )                           :: Temperature                        ! Temperature
REAL( Kind= Real64 )                           :: MinimumDensity, MaximumDensity     ! Molar densities of the critical points on a type-A curve
REAL( Kind= Real64 )                           :: UpperBoundVolume, LowerBoundVolume ! Molar density bounds
REAL( Kind= Real64 )                           :: mVolumeInflection                  ! Molar density at the inflection point
REAL( Kind= Real64 )                           :: MinimumVolume, MaximumVolume       ! Molar volumes of the critical points on a type-A curve
REAL( Kind= Real64 )                           :: LimitingVolume                     ! Molar volume (limit)
REAL( Kind= Real64 )                           :: Pressure                           ! Pressure (calculated)
REAL( Kind= Real64 )                           :: MinimumPressure, MaximumPressure   ! Pressures of the critical points on a type-A curve
REAL( Kind= Real64 )                           :: IsothermalCompressibility          ! Isothermal compressibility
REAL( Kind= Real64 )                           :: dPressure_dDensityInflection       ! First derivative of pressure with respect to the molar density at the inflection point
REAL( Kind= Real64 )                           :: aInterval, bInterval               ! Golden-section search interval (upper and lower bounds)
REAL( Kind= Real64 )                           :: ThermalExpansionCoefficient        ! Thermal expansion coefficient
REAL( Kind= Real64 )                           :: CompressibilityFactor              ! Compressibility factor
REAL( Kind= Real64 ), DIMENSION( nComponents ) :: mFraction                          ! Molar fraction of a component

! ************************************************************************************************ !
! LOGICAL VARIABLES                                                                                !
! ************************************************************************************************ !
LOGICAL :: PureComponent ! Checks whether the system is a pure component or a mixture

! Limiting volume (hard spheres)
IF( PureComponent ) THEN
  LimitingVolume = cPi * ijDiameterSphereCubic(cComponent,cComponent) * cAvogadro / 6.D0 / 0.999D0
ELSE
  LimitingVolume = cPi * cAvogadro * SUM( mFraction * cDiameterSphere * cDiameterSphere * cDiameterSphere ) / 6.D0 / 0.999D0
END IF

! Molar density bounds
LowerBoundVolume = LimitingVolume / 4.D-1
UpperBoundVolume = 1.D10

! ************************************************************************************************ !
! Inflection Category                                                                              !
! ************************************************************************************************ !
OptimizationOption = 1
aInterval = LowerBoundVolume
bInterval = UpperBoundVolume

! Find the inflection point with golden-section search method
IF( PureComponent ) THEN
  CALL Golden_Section_Single_Component( cComponent, Temperature, aInterval, bInterval, OptimizationOption, mVolumeInflection )
  CALL Calculate_Pressure_Single_Component( cComponent, mVolumeInflection, Temperature, Pressure, IsothermalCompressibility, &
  &                                         ThermalExpansionCoefficient, CompressibilityFactor )
ELSE
  CALL Golden_Section_Search( Temperature, aInterval, bInterval, OptimizationOption, mVolumeInflection )
  CALL Calculate_Pressure( mFraction, mVolumeInflection, Temperature, Pressure, IsothermalCompressibility, &
  &                        ThermalExpansionCoefficient, CompressibilityFactor )
END IF

! Density at the inflection point
dPressure_dDensityInflection = mVolumeInflection / IsothermalCompressibility

! ************************************************************************************************ !
! Local Minimum Category                                                                           !
! ************************************************************************************************ !
OptimizationOption = 2
aInterval = LowerBoundVolume
bInterval = mVolumeInflection

! Find the local minimum with golden-section search method
IF( PureComponent ) THEN
  CALL Golden_Section_Single_Component( cComponent, Temperature, aInterval, bInterval, OptimizationOption, MinimumVolume )
  CALL Calculate_Pressure_Single_Component( cComponent, MinimumVolume, Temperature, Pressure, IsothermalCompressibility, &
  &                                         ThermalExpansionCoefficient, CompressibilityFactor )
ELSE
  CALL Golden_Section_Search( Temperature, aInterval, bInterval, OptimizationOption, MinimumVolume )
  CALL Calculate_Pressure( mFraction, MinimumVolume, Temperature, Pressure, IsothermalCompressibility, &
  &                        ThermalExpansionCoefficient, CompressibilityFactor )
END IF

! Pressure and molar density at the local extremum (minimum)
MinimumPressure = Pressure
MinimumDensity  = 1.D0 / MinimumVolume

! ************************************************************************************************ !
! Local Maximum Category                                                                           !
! ************************************************************************************************ !
OptimizationOption = 3
aInterval = mVolumeInflection
bInterval = UpperBoundVolume

! Find the local maximum with golden-section search method
IF( PureComponent ) THEN
  CALL Golden_Section_Single_Component( cComponent, Temperature, aInterval, bInterval, OptimizationOption, MaximumVolume )
  CALL Calculate_Pressure_Single_Component( cComponent, MaximumVolume, Temperature, Pressure, IsothermalCompressibility, &
  &                                         ThermalExpansionCoefficient, CompressibilityFactor )
ELSE
  CALL Golden_Section_Search( Temperature, aInterval, bInterval, OptimizationOption, MaximumVolume )
  CALL Calculate_Pressure( mFraction, MaximumVolume, Temperature, Pressure, IsothermalCompressibility, &
  &                        ThermalExpansionCoefficient, CompressibilityFactor )
END IF

! Pressure and density at the local extremum (maximum)
MaximumPressure = Pressure
MaximumDensity  = 1.D0 / MaximumVolume

RETURN

END SUBROUTINE Find_Pressure_Interval