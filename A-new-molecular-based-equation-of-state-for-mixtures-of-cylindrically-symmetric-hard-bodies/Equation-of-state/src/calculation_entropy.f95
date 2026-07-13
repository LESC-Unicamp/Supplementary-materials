! ************************************************************************************************ !
!                                       ENTROPY CALCULATION                                        !
! ************************************************************************************************ !
!             This subroutine is used to calculate the entropy and the derivatives of              !
!                    the Helmholtz free energy with respect to the temperature.                    !
! ************************************************************************************************ !
! => AUTHOR:     Nathan Barros de Souza                                                            !
! => E-MAIL:     n264179@dac.unicamp.br                                                            !
! => SUPERVISOR: Luís Fernando Mercier Franco                                                      !
! ************************************************************************************************ !
! Main References:                          B.-J. Zhang                                            !
!                             Fluid Phase Equilibria, 154, 1-10 (1999)                             !
!                        DOI: https://doi.org/10.1016/S0378-3812(98)00431-2                        !
!                             --------------------------------------                               !
!                                 P. Gurin, S. Varga, G. Odriozola                                 !
!                              Molecular Liquids, 360, 119528 (2022)                               !
!                        DOI: https://doi.org/10.1016/j.molliq.2022.119528                         !
!                             --------------------------------------                               !
!       A. Gil-Villegas, A. Galindo, P. J. Whitehead, S. J. Mills, G. Jackson, A. N. Burgess       !
!                               J. Chem. Phys. 106, 4168–4186 (1997)                               !
!                              DOI: https://doi.org/10.1063/1.473101                               !
!                             --------------------------------------                               !
!                      A. Galindo, L. A. Davies, A. Gil-Villegas, G. Jackson                       !
!                                 J. Mol. Phys. 93, 241-252 (1998)                                 !
!                           DOI: https://doi.org/10.1080/002689798169249                           !
!                             --------------------------------------                               !
!                                   J. T. Lopes, L. F. M. Franco                                   !
!                              Molecular Liquids, 330, 115676 (2021)                               !
!                           DOI: https://doi.org/10.1080/002689798169249                           !
!                             --------------------------------------                               !
!                              A. K. Singh, U. N. Singh, S. K. Sinha                               !
!                                   J. Phys., 28, 343-353 (1987)                                   !
!                             DOI: https://doi.org/10.1007/BF02847095                              !
!                             --------------------------------------                               !
!                                      A. Galindo, F. J. Blas                                      !
!                             J. Phys. Chem. B, 106, 4503-4515 (2002)                              !
!                              DOI: https://doi.org/10.1021/jp013402h                              !
!                             --------------------------------------                               !
!                                          D. Constantin                                           !
!                                 Eur. Phys. J. E, 38, 116 (2015)                                  !
!                         DOI: https://doi.org/10.1140/epje/i2015-15116-2                          !
!                                                                                                  !
! ************************************************************************************************ !
!                       University of Campinas, Campinas - São Paulo, Brazil                       !
!                                         March 7th, 2024                                          !
! ************************************************************************************************ !
! Version 1.0.0                                                                                    !
! ************************************************************************************************ !

! ************************************************************************************************ !
! Calculates the entropy and the derivatives of the Helmholtz free energy with respect to the      !
! temperature for mixtures.                                                                        !
! ************************************************************************************************ !
SUBROUTINE Calculate_Entropy( mFraction, mVolume, Temperature, Entropy, vSpecificHeat )

! Uses two modules: global variables and substances
USE GlobalVariables
USE Substances

IMPLICIT NONE

! ************************************************************************************************ !
! INTEGER VARIABLES                                                                                !
! ************************************************************************************************ !
INTEGER( Kind= Int64 ) :: cComponent, iComponent, jComponent ! Component indexes
INTEGER( Kind= Int64 ) :: nOrder                             ! Counter of TPT coefficients

! ************************************************************************************************ !
! REAL PARAMETERS                                                                                  !
! ************************************************************************************************ !
REAL( Kind= Real64 ), PARAMETER                    :: CondensationPFraction = 0.493D0                                          ! Condensation packing fraction
REAL( Kind= Real64 ), PARAMETER                    :: ZhangCorrection = 1.D0 / (CondensationPFraction * CondensationPFraction) ! Zhang's correction
REAL( Kind= Real64 ), DIMENSION( 3, 3 ), PARAMETER :: SWCoefficientsMatrix = RESHAPE( [ [2.258550D0, -0.669270D0, 1.015760D1], &
&                                                                                       [-1.50349D0, 1.4004900D0, -1.50427D1], &
&                                                                                       [0.249434D0, -0.827739D0, 5.308270D0]  &
&                                                                                     ], [3, 3] )                              ! Matrix of coefficients of the effective packing fraction
REAL( Kind= Real64 ), DIMENSION( 2, 4 ), PARAMETER :: SLCoefficientsMatrix = RESHAPE( [ [-0.943973D0, 0.370942D0], &
&                                                                                       [0.422543D0, -0.173333D0], &
&                                                                                       [-3.71763D-2, 1.75599D-2], &
&                                                                                       [1.16901D-3, -5.72729D-4]  &
&                                                                                     ], [2, 4] )                              ! Matrix of coefficients of the effective packing fraction
REAL( Kind= Real64 ), DIMENSION( 2, 3 ), PARAMETER :: YK1CoefficientsMatrix = RESHAPE( [ [0.900678D0, -0.314300D0], &
&                                                                                        [-1.50051D0, 0.2571010D0], &
&                                                                                        [0.776577D0, -4.31566D-2]  &
&                                                                                      ], [2, 3] )                             ! Matrix of coefficients of the effective packing fraction
REAL( Kind= Real64 ), DIMENSION( 2, 5 ), PARAMETER :: YK2CoefficientsMatrix = RESHAPE( [ [0.989601D0, -1.19152D-2], &
&                                                                                        [-0.872203D0, -1.24029D0], &
&                                                                                        [0.320808D0, 2.41636D0], &
&                                                                                        [0.00000D0, -2.01922D0], &
&                                                                                        [0.00000D0, 0.647565D0]  &
&                                                                                      ], [2, 5] )                             ! Matrix of coefficients of the effective packing fraction

! ************************************************************************************************ !
! REAL VARIABLES                                                                                   !
! ************************************************************************************************ !
REAL( Kind= Real64 )                                           :: mVolume                                    ! Molar volume
REAL( Kind= Real64 )                                           :: Temperature                                ! Temperature
REAL( Kind= Real64 )                                           :: rNumberDensity                             ! Reduced number density
REAL( Kind= Real64 )                                           :: AngleAverageMixture                        ! Angle average of the mixture
REAL( Kind= Real64 )                                           :: HSSecondVirialCoefficientMixture           ! Second virial coefficient of a mixture of hard spheres
REAL( Kind= Real64 )                                           :: SecondVirialCoefficientMixtureSingle       ! Second virial coefficient of a mixture of non-spherical rigid bodies
REAL( Kind= Real64 )                                           :: SecondVirialCoefficientMixture             ! Second virial coefficient of a mixture of non-spherical rigid bodies
REAL( Kind= Real64 )                                           :: HigherOrderTerms                           ! Higher-order perturbation terms
REAL( Kind= Real64 )                                           :: Factorial                                  ! Factorial function
REAL( Kind= Real64 )                                           :: mHigherOrderHelmholtzFreeEnergy            ! Higher-order contribution to the Helmholtz free energy (molar basis)
REAL( Kind= Real64 )                                           :: mIdealHelmholtzFreeEnergy                  ! Ideal contribution to the Helmholtz free energy (molar basis)
REAL( Kind= Real64 )                                           :: mExcludedVolumeHelmholtzFreeEnergy         ! Excluded-volume contribution to the Helmholtz free energy (molar basis)
REAL( Kind= Real64 )                                           :: mMeanAttractiveHelmholtzFreeEnergy         ! Mean-attractive energy contribution to the Helmholtz free energy (molar basis)
REAL( Kind= Real64 )                                           :: mMeanAttFluctuationHelmholtzFreeEnergy     ! Mean-attractive energy fluctuation contribution to the Helmholtz free energy (molar basis)
REAL( Kind= Real64 )                                           :: mMonomerMonomerHelmholtzFreeEnergy         ! Monomer-monomer contribution to the Helmholtz free energy (molar basis)
REAL( Kind= Real64 )                                           :: mTotalHelmholtzFreeEnergy                  ! Total Helmholtz free energy (molar basis)
REAL( Kind= Real64 )                                           :: HSBoublikMansoori                          ! Boublik-Mansoori hard-sphere mixture term
REAL( Kind= Real64 )                                           :: HCBBoublik                                 ! Boublik hard convex body term
REAL( Kind= Real64 )                                           :: AuxDiameterRelationship                    ! Auxiliary factor (diameter relationship)
REAL( Kind= Real64 )                                           :: HSIsothermalCompressibility                ! Hard-sphere isothermal compressibility (Percus-Yevick expression)
REAL( Kind= Real64 )                                           :: HCBIsothermalCompressibility               ! Hard convex-body isothermal compressibility (Boublik expression)
REAL( Kind= Real64 )                                           :: pFirstOrderCoefficient                     ! First-order perturbation coefficient
REAL( Kind= Real64 )                                           :: pSecondOrderCoefficient                    ! Second-order perturbation coefficient
REAL( Kind= Real64 )                                           :: Entropy                                    ! Entropy
REAL( Kind= Real64 )                                           :: vSpecificHeat                              ! Specific heat capacity at constant volume
REAL( Kind= Real64 )                                           :: dIdealHelmholtzFreeEnergy_dTemperature     ! First derivative of the ideal Helmholtz free energy with respect to the temperature
REAL( Kind= Real64 )                                           :: dExcludedVolumeFreeEnergy_dTemperature     ! First derivative of the excluded-volume contribution to the Helmholtz free energy with respect to the temperature
REAL( Kind= Real64 )                                           :: dMeanAttFluctuationFEnergy_dTemperature    ! First derivative of the mean-attractive energy fluctuation contribution to the Helmholtz free energy with respect to the temperature
REAL( Kind= Real64 )                                           :: dHigherOrderTerms_dTemperature             ! First derivative of the higher-order perturbation terms with respect to the temperature
REAL( Kind= Real64 )                                           :: dHigherOrderFEnergy_dTemperature           ! First derivative of the higher-order perturbation contribution to the Helmholtz free energy with respect to the temperature
REAL( Kind= Real64 )                                           :: dMonomerMonomerFEnergy_dTemperature        ! First derivative of the monomer-monomer contribution to the Helmholtz free energy with respect to the temperature
REAL( Kind= Real64 )                                           :: dTotalHelmholtzFEnergy_dTemperature        ! First derivative of the total Helmholtz free energy with respect to the temperature
REAL( Kind= Real64 )                                           :: d2IdealHelmholtzFEnergy_d2Temperature      ! Second derivative of the ideal Helmholtz free energy with respect to the temperature
REAL( Kind= Real64 )                                           :: d2MeanAttFluctFEnergy_d2Temperature        ! Second derivative of the mean-attractive energy fluctuation contribution to the Helmholtz free energy with respect to the temperature
REAL( Kind= Real64 )                                           :: d2HigherOrderTerms_d2Temperature           ! Second derivative of the higher-order perturbation terms with respect to the temperature
REAL( Kind= Real64 )                                           :: d2HigherOrderFEnergy_d2Temperature         ! Second derivative of the higher-order perturbation contribution to the Helmholtz free energy with respect to the temperature
REAL( Kind= Real64 )                                           :: d2MonomerMonomerFEnergy_d2Temperature      ! Second derivative of the monomer-monomer contribution to the Helmholtz free energy with respect to the temperature
REAL( Kind= Real64 )                                           :: d2TotalHelmholtzFEnergy_d2Temperature      ! Second derivative of the total Helmholtz free energy with respect to the temperature
REAL( Kind= Real64 )                                           :: AspectRatioMixture                         ! Aspect ratio of the mixture
REAL( Kind= Real64 )                                           :: sumVolume, sumRadius, sumArea, sumSqRadius ! Morphological descriptors (sum)
REAL( Kind= Real64 )                                           :: sumVolumeSurface                           ! Morphological descriptors (sum)
REAL( Kind= Real64 )                                           :: geoAux, dAux1, dAux2                       ! Auxiliars
REAL( Kind= Real64 )                                           :: mNonSphericityMixture                      ! Non-sphericity parameter of the mixture
REAL( Kind= Real64 )                                           :: DiameterSphereMixture                      ! Diameter of the representative sphere of the mixture
REAL( Kind= Real64 )                                           :: DiameterSphereMixtureCubic                 ! Diameter of the representative sphere of the mixture (cubic)
REAL( Kind= Real64 )                                           :: rDensityMixture                            ! Reduced density of the mixture
REAL( Kind= Real64 )                                           :: ZhangFactor                                ! Zhang's correction factor
REAL( Kind= Real64 )                                           :: gAux1, gAux2                               ! Auxiliars (nonsphericity expressions
REAL( Kind= Real64 ), DIMENSION( 0:3 )                         :: rDensity                                   ! Reduced densities
REAL( Kind= Real64 ), DIMENSION( 0:3 )                         :: rDensityConstants                          ! Reduced density constants
REAL( Kind= Real64 ), DIMENSION( 2, 9 )                        :: CSWAlphaCoefficients                       ! Coefficients of the nonsphericity (Convex square-well potential)
REAL( Kind= Real64 ), DIMENSION( 3, 9 )                        :: CSW1CoefficientsMatrix                     ! Matrix of coefficients of the effective packing fraction (first-order coefficient)
REAL( Kind= Real64 ), DIMENSION( 3, 9 )                        :: CSW2CoefficientsMatrix                     ! Matrix of coefficients of the effective packing fraction (second-order coefficient)
REAL( Kind= Real64 ), DIMENSION( nComponents )                 :: MomentInertia                              ! Moment of inertia (non-spherical rigid bodies)
REAL( Kind= Real64 ), DIMENSION( nComponents )                 :: DeBroglie_Wavelength_Translational         ! Thermal de Broglie wavelength (translational contribution)
REAL( Kind= Real64 ), DIMENSION( nComponents )                 :: DeBroglie_Wavelength_Translational_Cb      ! Thermal de Broglie wavelength (translational contribution, cubic)
REAL( Kind= Real64 ), DIMENSION( nComponents )                 :: DeBroglie_Wavelength_Rotational            ! Thermal de Broglie wavelength (rotational contribution)
REAL( Kind= Real64 ), DIMENSION( nComponents )                 :: vIdealSpecificHeat                         ! Specific heat capacity at constant volume (ideal contribution)
REAL( Kind= Real64 ), DIMENSION( nComponents )                 :: mFraction                                  ! Molar fraction of a component
REAL( Kind= Real64 ), DIMENSION( nComponents )                 :: cMolecularVolume                           ! Molecular volume (component) [Å³]
REAL( Kind= Real64 ), DIMENSION( nComponents )                 :: cSurfaceArea                               ! Surface area (component) [Å²]
REAL( Kind= Real64 ), DIMENSION( nComponents )                 :: cCurvatureRadius                           ! Mean radius of curvature (component) [Å]
REAL( Kind= Real64 ), DIMENSION( nComponents )                 :: cVolumeSurface                             ! Ratio between volume and surface area (component) [Å]
REAL( Kind= Real64 ), DIMENSION( nComponents )                 :: cRadius                                    ! Radius (component) [Å]
REAL( Kind= Real64 ), DIMENSION( nComponents )                 :: fMolecularVolume                           ! Molecular volume (potential) [Å³]
REAL( Kind= Real64 ), DIMENSION( nComponents )                 :: fSurfaceArea                               ! Surface area (potential) [Å²]
REAL( Kind= Real64 ), DIMENSION( nComponents )                 :: fCurvatureRadius                           ! Mean radius of curvature (potential) [Å]
REAL( Kind= Real64 ), DIMENSION( nComponents )                 :: fRadius                                    ! Radius (potential) [Å]
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: ijShapeAnisotropy                          ! Shape anisotropy of ellipsoids-of-revolution
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: ijAngleAverage                             ! Angle average of the excluded volume of a pair of non-spherical rigid bodies
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: dEffRadialDistributionFunct_dDensity       ! First derivative of the radial distribution function with respect to the density
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: dEffRadialDistributionFunctYK_dDensity     ! First derivative of the radial distribution function with respect to the density (Yukawa potential)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: dEffRadialDistributionFunctSL_dDensity     ! First derivative of the radial distribution function with respect to the density (Sutherland potential)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: dEffRadialDistributionFunctCSW_dDensity    ! First derivative of the radial distribution function with respect to the density (Convex square-well potential)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: dMeanAttractiveEnergy_dDensity             ! First derivative of the mean-attractive energy with respect to the density
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: cAverageDiameterRelationship               ! Diameter relationship between components (the product of diameters divided by the sum of diameters)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: cRadialDistributionFunction                ! Contact radial distribution function
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: cEffectiveRadialDistributionFunction       ! Contact radial distribution function for an effective packing fraction
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: cEffectiveRadialDistributionFunctionYK     ! Contact radial distribution function for an effective packing fraction (Yukawa potential)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: cEffectiveRadialDistributionFunctionSL     ! Contact radial distribution function for an effective packing fraction (Sutherland potential)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: cEffectiveRadialDistributionFunctionCSW    ! Contact radial distribution function for an effective packing fraction (Convex square-well potential)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: cMeanAttractiveEnergy                      ! Mean-attractive energy between components
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: cMeanAttractiveEnergySutherland            ! Mean-attractive energy between components (Sutherland potential)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: cMeanAttractiveEnergyYukawa                ! Mean-attractive energy between components (Yukawa potential)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: cMeanAttractiveEnergyCSW                   ! Mean-attractive energy between components (Convex square-well potential)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: cMeanAttractiveEnergyFluctuations          ! Mean-attractive energy fluctuations between components
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: EffPackingFraction                         ! Effective packing fraction
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: EffPackingFractionYukawa                   ! Effective packing fraction (Yukawa potential)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: EffPackingFractionSutherland               ! Effective packing fraction (Sutherland potential)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: EffPackingFractionCSW                      ! Effective packing fraction (Convex square-well potential)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: dEffPackingFraction_dDensity               ! First derivative of the effective packing fraction with respect to the density
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: dEffPackingFractionYukawa_dDensity         ! First derivative of the effective packing fraction with respect to the density (Yukawa potential)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: dEffPackingFractionSutherland_dDensity     ! First derivative of the effective packing fraction with respect to the density (Sutherland potential)
REAL( Kind= Real64 ), DIMENSION( nComponents, nComponents )    :: dEffPackingFractionCSW_dDensity            ! First derivative of the effective packing fraction with respect to the density (Convex square-well potential)
REAL( Kind= Real64 ), DIMENSION( 3, nComponents, nComponents ) :: EffPackingFractionCoefficients             ! Coefficients of the effective packing fraction
REAL( Kind= Real64 ), DIMENSION( 3, nComponents, nComponents ) :: EffPackingFractionCoefficientsYukawa       ! Coefficients of the effective packing fraction (Yukawa potential)
REAL( Kind= Real64 ), DIMENSION( 3, nComponents, nComponents ) :: EffPackingFractionCoefficientsSutherland   ! Coefficients of the effective packing fraction (Sutherland potential)
REAL( Kind= Real64 ), DIMENSION( 3, nComponents, nComponents ) :: EffPackingFractionCoefficientsCSW          ! Coefficients of the effective packing fraction (convex square-well potential)

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !
!                                          INITIALIZATION                                          !
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !

! Reduced number density [1 / m³]
rNumberDensity = cAvogadro / mVolume

! Reduced number density [1 / Å³]
rNumberDensity = rNumberDensity * 1.D-30

! Mixing rule for the ratio of second virial coefficients [Å³]
CALL Mixing_Rules_Single( mFraction, ijSecondVirialCoefficient, SecondVirialCoefficientMixtureSingle )

! Zhang's correction [unitless]
IF( ZhangCorrectionLogical ) THEN
  IF( .NOT. PotentialTypeLogical(4) ) THEN
    ZhangFactor = 1.D0 + ( 2.D0 * ZhangCorrection * ( cPi * rNumberDensity * SUM( mFraction * aDiameterSphere * aDiameterSphere * &
    &     aDiameterSphere ) / 6.D0 ) * ( cPi * rNumberDensity * SUM( mFraction * aDiameterSphere * aDiameterSphere * &
    &     aDiameterSphere ) / 6.D0 ) )
  ELSE
    ZhangFactor = 1.D0 + ( 0.125D0 * ZhangCorrection * ( rNumberDensity * SecondVirialCoefficientMixtureSingle ) * ( &
    &     rNumberDensity * SecondVirialCoefficientMixtureSingle ) )
  END IF
END IF

! Morphological descriptors
DO cComponent = 1, nComponents
  ! Radius
  cRadius(cComponent) = 0.5D0 * aDiameter(cComponent)
  IF( GeometrySpecification(cComponent,1) ) THEN ! Ellipsoids-of-revolution
    ! Molecular volume [Å³]
    cMolecularVolume(cComponent) = ( cPi / 6.D0 ) * aDiameter(cComponent) * aDiameter(cComponent) * aDiameter(cComponent) * &
    &     cAspectRatio(cComponent)
    ! Prolate ellipsoids-of-revolution
    IF( cAspectRatio(cComponent) > 1.D0 ) THEN
      ! Radius of curvature [Å]
      cCurvatureRadius(cComponent) = 0.5D0 * cAspectRatio(cComponent) * cRadius(cComponent) + ( 0.5D0 * cRadius(cComponent) / &
      &     DSQRT( cAspectRatio(cComponent) * cAspectRatio(cComponent) - 1.D0 ) ) * DATANH( DSQRT( cAspectRatio(cComponent) * &
      &     cAspectRatio(cComponent) - 1.D0 ) / cAspectRatio(cComponent) )
      ! Surface area [Å²]
      cSurfaceArea(cComponent) = ( 2.D0 * cPi * cRadius(cComponent) * cRadius(cComponent) / DSQRT( cAspectRatio(cComponent) * &
      &     cAspectRatio(cComponent) - 1.D0 ) ) * ( ( cAspectRatio(cComponent) * cAspectRatio(cComponent) * DASIN( DSQRT( &
      &     cAspectRatio(cComponent) * cAspectRatio(cComponent) - 1.D0 ) / cAspectRatio(cComponent) ) ) + ( DSQRT( &
      &     cAspectRatio(cComponent) * cAspectRatio(cComponent) - 1.D0 ) ) )
    ! Oblate ellipsoids-of-revolution
    ELSE IF( cAspectRatio(cComponent) < 1.D0 ) THEN
      ! Radius of curvature [Å]
      cCurvatureRadius(cComponent) = 0.5D0 * cAspectRatio(cComponent) * cRadius(cComponent) + ( 0.5D0 * cRadius(cComponent) / &
      &     DSQRT( 1.D0 - cAspectRatio(cComponent) * cAspectRatio(cComponent) ) ) * DATAN( DSQRT( 1.D0 - cAspectRatio(cComponent) &
      &     * cAspectRatio(cComponent) ) / cAspectRatio(cComponent) )
      ! Surface area [Å²]
      cSurfaceArea(cComponent) = ( 2.D0 * cPi * cRadius(cComponent) * cRadius(cComponent) / DSQRT( 1.D0 - cAspectRatio(cComponent) &
      &     * cAspectRatio(cComponent) ) ) * ( ( cAspectRatio(cComponent) * cAspectRatio(cComponent) * DLOG( ( 1.D0 + DSQRT( 1.D0 &
      &     - cAspectRatio(cComponent) * cAspectRatio(cComponent) ) ) / cAspectRatio(cComponent) ) ) + ( DSQRT( 1.D0 - &
      &     cAspectRatio(cComponent) * cAspectRatio(cComponent) ) ) )
    ! Spheres
    ELSE
      ! Radius of curvature [Å]
      cCurvatureRadius(cComponent) = cRadius(cComponent)
      ! Surface area [Å²]
      cSurfaceArea(cComponent) = 4.D0 * cPi * cRadius(cComponent) * cRadius(cComponent)
    END IF
  ELSE IF( GeometrySpecification(cComponent,2) ) THEN ! Spherocylinders
    ! Molecular volume [Å³]
    cMolecularVolume(cComponent) = ( cPi / 4.D0 ) * aDiameter(cComponent) * aDiameter(cComponent) * aDiameter(cComponent) * &
    &     ( cAspectRatio(cComponent) + (2.D0 / 3.D0) )
    ! Prolate spherocylinders
    IF( cAspectRatio(cComponent) > 0.D0 ) THEN
      ! Radius of curvature [Å]
      cCurvatureRadius(cComponent) = 0.5D0 * ( cAspectRatio(cComponent) + 2.D0 ) * cRadius(cComponent)
      ! Surface area [Å²]
      cSurfaceArea(cComponent) = 4.D0 * cPi * cRadius(cComponent) * cRadius(cComponent) * ( cAspectRatio(cComponent) + 1.D0 )
    ! Spheres
    ELSE
      ! Radius of curvature [Å]
      cCurvatureRadius(cComponent) = cRadius(cComponent)
      ! Surface area [Å²]
      cSurfaceArea(cComponent) = 4.D0 * cPi * cRadius(cComponent) * cRadius(cComponent)
    END IF
  ELSE IF( GeometrySpecification(cComponent,3) ) THEN ! Cylinders
    ! Molecular volume [Å³]
    cMolecularVolume(cComponent) = ( cPi / 4.D0 ) * aDiameter(cComponent) * aDiameter(cComponent) * aDiameter(cComponent) * &
    &     cAspectRatio(cComponent)
    ! Radius of curvature [Å]
    cCurvatureRadius(cComponent) = 0.5D0 * ( cAspectRatio(cComponent) + ( 0.5D0 * cPi ) ) * cRadius(cComponent)
    ! Surface area [Å²]
    cSurfaceArea(cComponent) = 2.D0 * cPi * cRadius(cComponent) * cRadius(cComponent) * ( 1.D0 + 2.D0 * cAspectRatio(cComponent) )
  END IF
  cVolumeSurface(cComponent) = 4.D0 * cPi * cCurvatureRadius(cComponent) * cCurvatureRadius(cComponent) * &
  &     cCurvatureRadius(cComponent) / cSurfaceArea(cComponent)
END DO

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !
!                                HELMHOLTZ FREE ENERGY CALCULATION                                 !
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !

! ************************************************************************************************ !
! IDEAL CONTRIBUTION                                                                               !
! ************************************************************************************************ !

! Ideal contribution to the Helmholtz free energy (molar basis) [J / mol] - Translational de Broglie wavelength contribution
DeBroglie_Wavelength_Translational(:) = 1.D0 / DSQRT( 2.D0 * cPi * cMolarMass(:) * cUniversalGas * Temperature )
DeBroglie_Wavelength_Translational_Cb = DeBroglie_Wavelength_Translational * DeBroglie_Wavelength_Translational * &
&     DeBroglie_Wavelength_Translational

! Ideal contribution to the Helmholtz free energy (molar basis) [J / mol] - Rotational de Broglie wavelength contribution
DO cComponent = 1, nComponents
  IF( GeometrySpecification(cComponent,1) ) THEN ! Ellipsoids of revolution
    MomentInertia(cComponent) = ( 1.D0 / 20.D0 ) * cMolarMass(cComponent) * ( aDiameter(cComponent) * aDiameter(cComponent) + &
    &     aLength(cComponent) * aLength(cComponent) ) ! [Kg . Å² / mol]
    DeBroglie_Wavelength_Rotational(cComponent) = (cDeBroglie * cDeBroglie) / (cPi * MomentInertia(cComponent) * cUniversalGas * &
    &     Temperature) ! [m² / Å²]
    DeBroglie_Wavelength_Rotational(cComponent) = DeBroglie_Wavelength_Rotational(cComponent) * 1.D20 ! [unitless]
  ELSE IF( GeometrySpecification(cComponent,2) ) THEN ! Spherocylinders
    MomentInertia(cComponent) = (3.D0 / 8.D0) * ( ( cMolarMass(cComponent) / ( 3.D0 * cAspectRatio(cComponent) + 2.D0 ) ) * &
    &     ( aDiameter(cComponent) * aDiameter(cComponent) ) ) * ( ( ( cAspectRatio(cComponent) / 6.D0 ) * ( 3.D0 + ( 4.D0 * &
    &     cAspectRatio(cComponent) * cAspectRatio(cComponent) ) ) ) + ( (4.D0 / 3.D0) * ( (83.D0 / 320.D0) + &
    &     ( cAspectRatio(cComponent) + (3.D0 / 8.D0) ) * ( cAspectRatio(cComponent) + (3.D0 / 8.D0) ) ) ) ) ! [Kg . Å² / mol]
    DeBroglie_Wavelength_Rotational(cComponent) = (cDeBroglie * cDeBroglie) / (cPi * MomentInertia(cComponent) * cUniversalGas * &
    &     Temperature) ! [m² / Å²]
    DeBroglie_Wavelength_Rotational(cComponent) = DeBroglie_Wavelength_Rotational(cComponent) * 1.D20 ! [unitless]
  ELSE IF( GeometrySpecification(cComponent,3) ) THEN ! Cylinders
    MomentInertia(cComponent) = (1.D0 / 48.D0) * cMolarMass(cComponent) * ( ( 4.D0 * aLength(cComponent) * aLength(cComponent) ) + &
    &     ( 3.D0 * aDiameter(cComponent) * aDiameter(cComponent) ) ) ! [Kg . Å² / mol]
    DeBroglie_Wavelength_Rotational(cComponent) = (cDeBroglie * cDeBroglie) / (cPi * MomentInertia(cComponent) * cUniversalGas * &
    &     Temperature) ! [m² / Å²]
    DeBroglie_Wavelength_Rotational(cComponent) = DeBroglie_Wavelength_Rotational(cComponent) * 1.D20 ! [unitless]
  END IF
END DO

! Ideal contribution to the Helmholtz free energy (molar basis) [J / mol]
mIdealHelmholtzFreeEnergy = 0.D0
DO iComponent = 1, nComponents
  ! Ideal contribution to the Helmholtz free energy per component [unitless]
  mIdealHelmholtzFreeEnergy = mIdealHelmholtzFreeEnergy + mFraction(iComponent) * DLOG( cDeBroglieIdeal * mFraction(iComponent) * &
  &     DeBroglie_Wavelength_Translational_Cb(iComponent) * DeBroglie_Wavelength_Rotational(iComponent) / mVolume )
END DO
mIdealHelmholtzFreeEnergy = mIdealHelmholtzFreeEnergy - 1.D0
mIdealHelmholtzFreeEnergy = cUniversalGas * Temperature * mIdealHelmholtzFreeEnergy ! Proven units

! ************************************************************************************************ !
! EXCLUDED-VOLUME CONTRIBUTION                                                                     !
! ************************************************************************************************ !

! Reduced density constants
rDensityConstants(0) = (cPi / 6.D0) ! [unitless]
rDensityConstants(1) = (cPi / 6.D0) * SUM( mFraction * aDiameterSphere ) ! [Å]
rDensityConstants(2) = (cPi / 6.D0) * SUM( mFraction * aDiameterSphere * aDiameterSphere ) ! [Å²]
rDensityConstants(3) = (cPi / 6.D0) * SUM( mFraction * aDiameterSphere * aDiameterSphere * aDiameterSphere ) ! [Å³]

! Reduced densities
rDensity(0) = rDensityConstants(0) * rNumberDensity ! [1 / Å³]
rDensity(1) = rDensityConstants(1) * rNumberDensity ! [1 / Å²]
rDensity(2) = rDensityConstants(2) * rNumberDensity ! [1 / Å]
rDensity(3) = rDensityConstants(3) * rNumberDensity ! [unitless]

! Boublik-Mansoori hard-sphere mixture term [1 / Å³]
HSBoublikMansoori = ( ( rDensity(2) * rDensity(2) * rDensity(2) / ( rDensity(3) * rDensity(3) ) ) - rDensity(0) ) * DLOG( 1.D0 - &
&     rDensity(3) ) + ( 3.D0 * rDensity(1) * rDensity(2) / ( 1.D0 - rDensity(3) ) ) + ( rDensity(2) * rDensity(2) * rDensity(2) / &
&     ( rDensity(3) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) )

! Boublik-Mansoori hard-sphere mixture term [unitless]
HSBoublikMansoori = 6.D0 * HSBoublikMansoori / cPi / rNumberDensity

! Mixing rule for the angle average of a mixture of non-spherical rigid bodies [unitless]
IF( NonSphericalMixingRule == 1 ) THEN ! [1] Rule applied directly to the aspect ratio of the mixture
  ! Mixing rule for the aspect ratio [unitless]
  CALL Mixing_Rules( mFraction, ijAspectRatio, AspectRatioMixture )
  ! Angle average of the excluded volume of a pair of non-spherical rigid bodies [unitless]
  IF( ALL( GeometrySpecification(:,1) ) ) THEN ! Ellipsoids-of-revolution
    ! Analytical expression of the angle average of the excluded volume of a pair of non-spherical rigid bodies [unitless]
    IF( AspectRatioMixture /= 1.D0 ) THEN
      IF( AspectRatioMixture > 1.D0 ) THEN
        AngleAverageMixture = 2.D0 + 1.5D0 * ( 1.D0 + DATANH( DSQRT( ( AspectRatioMixture * AspectRatioMixture ) - 1.D0 ) / &
        &     AspectRatioMixture ) / ( AspectRatioMixture * DSQRT( ( AspectRatioMixture * AspectRatioMixture ) - 1.D0 ) ) ) * &
        &     ( 1.D0 + ( ( AspectRatioMixture * AspectRatioMixture ) / DSQRT( ( AspectRatioMixture * AspectRatioMixture ) - &
        &     1.D0 ) ) * DATAN( DSQRT( ( AspectRatioMixture * AspectRatioMixture ) - 1.D0 ) ) )
      ELSE IF( AspectRatioMixture < 1.D0 ) THEN
        AngleAverageMixture = 2.D0 + 1.5D0 * ( 1.D0 + DATANH( DSQRT( ( ( 1.D0 / AspectRatioMixture ) * ( 1.D0 / &
        &     AspectRatioMixture ) ) - 1.D0 ) / ( 1.D0 / AspectRatioMixture ) ) / ( ( 1.D0 / AspectRatioMixture ) * &
        &     DSQRT( ( ( 1.D0 / AspectRatioMixture ) * ( 1.D0 / AspectRatioMixture ) ) - 1.D0 ) ) ) * ( 1.D0 + ( ( ( 1.D0 / &
        &     AspectRatioMixture ) * ( 1.D0 / AspectRatioMixture ) ) / DSQRT( ( ( 1.D0 / AspectRatioMixture ) * ( 1.D0 &
        &     / AspectRatioMixture ) ) - 1.D0 ) ) * DATAN( DSQRT( ( ( 1.D0 / AspectRatioMixture ) * ( 1.D0 / &
        &     AspectRatioMixture ) ) - 1.D0 ) ) )
      END IF
      AngleAverageMixture = 0.125D0 * AngleAverageMixture
    ELSE
      AngleAverageMixture = 1.D0
    END IF
  ELSE IF( ALL( GeometrySpecification(:,2) ) ) THEN ! Spherocylinders
    AngleAverageMixture = ( 12.D0 / ( ( 3.D0 * AspectRatioMixture ) + 2.D0 ) ) * ( (4.D0 / 3.D0) + ( 2.D0 * AspectRatioMixture ) + &
    &     ( ( AspectRatioMixture * AspectRatioMixture ) / 2.D0 ) )
    AngleAverageMixture = 0.125D0 * AngleAverageMixture
  ELSE IF( ALL( GeometrySpecification(:,3) ) ) THEN ! Cylinders
    AngleAverageMixture = 2.D0 * AspectRatioMixture + cPi + 3.D0 + ( cPi / ( 2.D0 * AspectRatioMixture ) )
    AngleAverageMixture = 0.125D0 * AngleAverageMixture
  END IF
ELSE IF( NonSphericalMixingRule == 2 ) THEN ! [2] Rule applied directly to the angle average of the excluded volume of a mixture of non-spherical rigid bodies
  ! Angle average of the excluded volume of a pair of non-spherical rigid bodies [unitless]
  IF( ALL( GeometrySpecification(:,1) ) ) THEN ! Ellipsoids-of-revolution
    DO iComponent = 1, nComponents
      DO jComponent = 1, nComponents
        ijShapeAnisotropy(iComponent,jComponent) = ( ( ijAspectRatio(iComponent,jComponent) * &
        &     ijAspectRatio(iComponent,jComponent) ) - 1.D0 ) / ( ( ijAspectRatio(iComponent,jComponent) * &
        &     ijAspectRatio(iComponent,jComponent) ) + 1.D0 )
        ! Analytical expression of the angle average of the excluded volume of a pair of non-spherical rigid bodies [unitless]
        IF( ijAspectRatio(iComponent,jComponent) /= 1.D0 ) THEN
          IF( ijAspectRatio(iComponent,jComponent) > 1.D0 ) THEN
            ijAngleAverage(iComponent,jComponent) = 2.D0 + 1.5D0 * ( 1.D0 + DATANH( DSQRT( ( &
            &     ijAspectRatio(iComponent,jComponent) * ijAspectRatio(iComponent,jComponent) ) - 1.D0 ) / &
            &     ijAspectRatio(iComponent,jComponent) ) / ( ijAspectRatio(iComponent,jComponent) * DSQRT( ( &
            &     ijAspectRatio(iComponent,jComponent) * ijAspectRatio(iComponent,jComponent) ) - 1.D0 ) ) ) * ( 1.D0 + ( ( &
            &     ijAspectRatio(iComponent,jComponent) * ijAspectRatio(iComponent,jComponent) ) / DSQRT( ( &
            &     ijAspectRatio(iComponent,jComponent) * ijAspectRatio(iComponent,jComponent) ) - 1.D0 ) ) * DATAN( DSQRT( ( &
            &     ijAspectRatio(iComponent,jComponent) * ijAspectRatio(iComponent,jComponent) ) - 1.D0 ) ) )
          ELSE IF( ijAspectRatio(iComponent,jComponent) < 1.D0 ) THEN
            ijAngleAverage(iComponent,jComponent) = 2.D0 + 1.5D0 * ( 1.D0 + DATANH( DSQRT( ( ( 1.D0 / &
            &     ijAspectRatio(iComponent,jComponent) ) * ( 1.D0 / ijAspectRatio(iComponent,jComponent) ) ) - 1.D0 ) / &
            &     ( 1.D0 / ijAspectRatio(iComponent,jComponent) ) ) / ( ( 1.D0 / ijAspectRatio(iComponent,jComponent) ) * &
            &     DSQRT( ( ( 1.D0 / ijAspectRatio(iComponent,jComponent) ) * ( 1.D0 / ijAspectRatio(iComponent,jComponent) ) ) &
            &     - 1.D0 ) ) ) * ( 1.D0 + ( ( ( 1.D0 / ijAspectRatio(iComponent,jComponent) ) * ( 1.D0 / &
            &     ijAspectRatio(iComponent,jComponent) ) ) / DSQRT( ( ( 1.D0 / ijAspectRatio(iComponent,jComponent) ) * ( 1.D0 &
            &     / ijAspectRatio(iComponent,jComponent) ) ) - 1.D0 ) ) * DATAN( DSQRT( ( ( 1.D0 / &
            &     ijAspectRatio(iComponent,jComponent) ) * ( 1.D0 / ijAspectRatio(iComponent,jComponent) ) ) - 1.D0 ) ) )
          END IF
          ijAngleAverage(iComponent,jComponent) = 0.125D0 * ijAngleAverage(iComponent,jComponent)
        ELSE
          ijAngleAverage(iComponent,jComponent) = 1.D0
        END IF
      END DO
    END DO
  ELSE IF( ALL( GeometrySpecification(:,2) ) ) THEN ! Spherocylinders
    DO iComponent = 1, nComponents
      DO jComponent = 1, nComponents
        ijAngleAverage(iComponent,jComponent) = ( 12.D0 / ( ( 3.D0 * ijAspectRatio(iComponent,jComponent) ) + 2.D0 ) ) * ( (4.D0 / &
        &     3.D0) + ( 2.D0 * ijAspectRatio(iComponent,jComponent) ) + ( ( ijAspectRatio(iComponent,jComponent) * &
        &     ijAspectRatio(iComponent,jComponent) ) / 2.D0 ) )
        ijAngleAverage(iComponent,jComponent) = 0.125D0 * ijAngleAverage(iComponent,jComponent)
      END DO
    END DO
  ELSE IF( ALL( GeometrySpecification(:,3) ) ) THEN ! Cylinders
    DO iComponent = 1, nComponents
      DO jComponent = 1, nComponents
        ijAngleAverage(iComponent,jComponent) = 2.D0 * ijAspectRatio(iComponent,jComponent) + cPi + 3.D0 + ( cPi / ( 2.D0 * &
        &     ijAspectRatio(iComponent,jComponent) ) )
        ijAngleAverage(iComponent,jComponent) = 0.125D0 * ijAngleAverage(iComponent,jComponent)
      END DO
    END DO
  END IF
  ! Mixing rule for the aspect ratio [unitless]
  CALL Mixing_Rules( mFraction, ijAngleAverage, AngleAverageMixture )
ELSE IF( NonSphericalMixingRule == 3 ) THEN ! [3] Rule applied to the ratio of the second virial coefficients (Isihara-Hadwiger theorem)
  ! Mixing rule for the ratio of second virial coefficients [unitless]
  CALL Mixing_Rules( mFraction, ijRatioSecondVirialCoefficient, AngleAverageMixture )
ELSE IF( NonSphericalMixingRule == 4 ) THEN ! [4] Rule applied to the second virial coefficients (Isihara-Hadwiger theorem)
  ! Mixing rule for the ratio of second virial coefficients [Å³]
  CALL Mixing_Rules( mFraction, ijSecondVirialCoefficient, SecondVirialCoefficientMixture )
  ! Mixing rule for the ratio of second virial coefficients of hard spheres [Å³]
  CALL Mixing_Rules( mFraction, ijHSSecondVirialCoefficient, HSSecondVirialCoefficientMixture )
  ! Mixing rule for the angle average of the excluded volume of a mixture [unitless]
  AngleAverageMixture = SecondVirialCoefficientMixture / HSSecondVirialCoefficientMixture
END IF

! Excluded-volume contribution to the Helmholtz free energy (molar basis) [J / mol]
mExcludedVolumeHelmholtzFreeEnergy = HSBoublikMansoori * cUniversalGas * AngleAverageMixture * Temperature ! Proven units

! Boublik's reference term (hard convex bodies)
IF( ReferenceBoublikLogical ) THEN
  ! Morphological descriptors
  sumVolume = SUM( mFraction * cMolecularVolume )
  sumRadius = SUM( mFraction * cCurvatureRadius )
  sumArea = SUM( mFraction * cSurfaceArea )
  sumSqRadius = SUM( mFraction * cCurvatureRadius * cCurvatureRadius )
  ! Reference term
  HCBBoublik = ( ( ( sumSqRadius * sumArea * sumArea / ( 9.D0 * sumVolume * sumVolume ) ) - 1.D0 ) * DLOG( 1.D0 - rDensity(3) ) ) &
  &     + ( ( sumRadius * sumArea * rDensity(3) ) / ( sumVolume * ( 1.D0 - rDensity(3) ) ) ) + ( ( sumSqRadius * sumArea * &
  &     sumArea * rDensity(3) ) / ( 9.D0 * sumVolume * sumVolume * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) )
  ! Excluded-volume contribution to the Helmholtz free energy (molar basis) [J / mol]
  mExcludedVolumeHelmholtzFreeEnergy = HCBBoublik * cUniversalGas * Temperature ! Proven units
END IF

! ************************************************************************************************ !
! FIRST-ORDER PERTURBATION CONTRIBUTION                                                            !
! ************************************************************************************************ !

! Auxiliary diameter relationship [1 / Å]
AuxDiameterRelationship = SUM( mFraction * aDiameterSphere * aDiameterSphere ) / SUM( mFraction * aDiameterSphere * &
&     aDiameterSphere * aDiameterSphere )

! Mixing rule for the diameter of the mixture [Å³]
CALL Mixing_Rules( mFraction, ijaDiameterSphereCubic, DiameterSphereMixtureCubic )

! Diameter of mixture [Å]
DiameterSphereMixture = DiameterSphereMixtureCubic ** (1.D0 / 3.D0)

! Reduced density of the mixture [unitless]
rDensityMixture = (cPi / 6.D0) * rNumberDensity * DiameterSphereMixtureCubic

! Non-sphericity parameter of the mixture [unitless]
mNonSphericityMixture = SUM( mFraction * cCurvatureRadius * cSurfaceArea ) / ( 3.D0 * SUM( mFraction * cMolecularVolume ) )

! Morphological descriptors (field)
DO cComponent = 1, nComponents
  IF( GeometrySpecification(cComponent,1) ) THEN ! Ellipsoids-of-revolution
    ! Molecular volume [Å³]
    fMolecularVolume(cComponent) = ( cPi / 6.D0 ) * aDiameterField(cComponent) * aDiameterField(cComponent) * &
    &     aDiameterField(cComponent) * fAspectRatio(cComponent)
    ! Prolate ellipsoids-of-revolution
    IF( fAspectRatio(cComponent) > 1.D0 ) THEN
      ! Radius of curvature [Å]
      fCurvatureRadius(cComponent) = 0.5D0 * fAspectRatio(cComponent) * fRadius(cComponent) + ( 0.5D0 * fRadius(cComponent) / &
      &     DSQRT( fAspectRatio(cComponent) * fAspectRatio(cComponent) - 1.D0 ) ) * DATANH( DSQRT( fAspectRatio(cComponent) * &
      &     fAspectRatio(cComponent) - 1.D0 ) / fAspectRatio(cComponent) )
      ! Surface area [Å²]
      fSurfaceArea(cComponent) = ( 2.D0 * cPi * fRadius(cComponent) * fRadius(cComponent) / DSQRT( fAspectRatio(cComponent) * &
      &     fAspectRatio(cComponent) - 1.D0 ) ) * ( ( fAspectRatio(cComponent) * fAspectRatio(cComponent) * DASIN( DSQRT( &
      &     fAspectRatio(cComponent) * fAspectRatio(cComponent) - 1.D0 ) / fAspectRatio(cComponent) ) ) + ( DSQRT( &
      &     fAspectRatio(cComponent) * fAspectRatio(cComponent) - 1.D0 ) ) )
    ! Oblate ellipsoids-of-revolution
    ELSE IF( fAspectRatio(cComponent) < 1.D0 ) THEN
      ! Radius of curvature [Å]
      fCurvatureRadius(cComponent) = 0.5D0 * fAspectRatio(cComponent) * fRadius(cComponent) + ( 0.5D0 * fRadius(cComponent) / &
      &     DSQRT( 1.D0 - fAspectRatio(cComponent) * fAspectRatio(cComponent) ) ) * DATAN( DSQRT( 1.D0 - fAspectRatio(cComponent) &
      &     * fAspectRatio(cComponent) ) / fAspectRatio(cComponent) )
      ! Surface area [Å²]
      fSurfaceArea(cComponent) = ( 2.D0 * cPi * fRadius(cComponent) * fRadius(cComponent) / DSQRT( 1.D0 - fAspectRatio(cComponent) &
      &     * fAspectRatio(cComponent) ) ) * ( ( fAspectRatio(cComponent) * fAspectRatio(cComponent) * DLOG( ( 1.D0 + DSQRT( 1.D0 &
      &     - fAspectRatio(cComponent) * fAspectRatio(cComponent) ) ) / fAspectRatio(cComponent) ) ) + ( DSQRT( 1.D0 - &
      &     fAspectRatio(cComponent) * fAspectRatio(cComponent) ) ) )
    ! Spheres
    ELSE
      ! Radius of curvature [Å]
      fCurvatureRadius(cComponent) = fRadius(cComponent)
      ! Surface area [Å²]
      fSurfaceArea(cComponent) = 4.D0 * cPi * fRadius(cComponent) * fRadius(cComponent)
    END IF
  ELSE IF( GeometrySpecification(cComponent,2) .OR. GeometrySpecification(cComponent,3) ) THEN ! Spherocylinders and cylinders
    ! Molecular volume [Å³]
    fMolecularVolume(cComponent) = ( cPi / 4.D0 ) * aDiameterField(cComponent) * aDiameterField(cComponent) * &
    &     aDiameterField(cComponent) * ( fAspectRatio(cComponent) + (2.D0 / 3.D0) )
    ! Prolate spherocylinders
    IF( fAspectRatio(cComponent) > 0.D0 ) THEN
      ! Radius of curvature [Å]
      fCurvatureRadius(cComponent) = 0.5D0 * ( fAspectRatio(cComponent) + 2.D0 ) * fRadius(cComponent)
      ! Surface area [Å²]
      fSurfaceArea(cComponent) = 4.D0 * cPi * fRadius(cComponent) * fRadius(cComponent) * ( fAspectRatio(cComponent) + 1.D0 )
    ! Spheres
    ELSE
      ! Radius of curvature [Å]
      fCurvatureRadius(cComponent) = fRadius(cComponent)
      ! Surface area [Å²]
      fSurfaceArea(cComponent) = 4.D0 * cPi * fRadius(cComponent) * fRadius(cComponent)
    END IF
  END IF
END DO

! Morphological descriptors
sumVolume = SUM( mFraction * cMolecularVolume )
sumRadius = SUM( mFraction * cCurvatureRadius )
sumArea = SUM( mFraction * cSurfaceArea )
sumSqRadius = SUM( mFraction * cCurvatureRadius * cCurvatureRadius )
sumVolumeSurface = SUM( mFraction * cVolumeSurface )

! Mean-attractive energy between components
DO iComponent = 1, nComponents
  DO jComponent = 1, nComponents
    ! Diameter relationship between components (the product of diameters divided by the sum of diameters) [Å]
    cAverageDiameterRelationship(iComponent,jComponent) = aDiameterSphere(iComponent) * aDiameterSphere(jComponent) / &
    &     ( aDiameterSphere(iComponent) + aDiameterSphere(jComponent) )
    ! Effective packing fraction coefficients
    IF( PotentialTypeLogical(4) ) THEN ! Convex square-well potential
      CSW1CoefficientsMatrix = CSWAllCoefficientA1
      CSW2CoefficientsMatrix = CSWAllCoefficientA2
      IF( PYHCBCorrectionLogical )CSW2CoefficientsMatrix = CSWAllCoefficientA2
      IF( UseA1ForA2Logical ) CSW2CoefficientsMatrix = CSWAllCoefficientA1
    END IF
    ! Parametrization of the effective packing fraction [unitless]
    IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
      EffPackingFractionCoefficients(1,iComponent,jComponent) = SWCoefficientsMatrix(1,1) + SWCoefficientsMatrix(1,2) * &
      &     ijPotentialRange(iComponent,jComponent) + SWCoefficientsMatrix(1,3) * ijPotentialRange(iComponent,jComponent) * &
      &     ijPotentialRange(iComponent,jComponent)
      EffPackingFractionCoefficients(2,iComponent,jComponent) = SWCoefficientsMatrix(2,1) + SWCoefficientsMatrix(2,2) * &
      &     ijPotentialRange(iComponent,jComponent) + SWCoefficientsMatrix(2,3) * ijPotentialRange(iComponent,jComponent) * &
      &     ijPotentialRange(iComponent,jComponent)
      EffPackingFractionCoefficients(3,iComponent,jComponent) = SWCoefficientsMatrix(3,1) + SWCoefficientsMatrix(3,2) * &
      &     ijPotentialRange(iComponent,jComponent) + SWCoefficientsMatrix(3,3) * ijPotentialRange(iComponent,jComponent) * &
      &     ijPotentialRange(iComponent,jComponent)
    ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
      EffPackingFractionCoefficients(1,iComponent,jComponent) = SLCoefficientsMatrix(1,1) + SLCoefficientsMatrix(1,2) * &
      &     ijPotentialRange(iComponent,jComponent) + SLCoefficientsMatrix(1,3) * ijPotentialRange(iComponent,jComponent) * &
      &     ijPotentialRange(iComponent,jComponent) + SLCoefficientsMatrix(1,4) * ijPotentialRange(iComponent,jComponent) * &
      &     ijPotentialRange(iComponent,jComponent) * ijPotentialRange(iComponent,jComponent)
      EffPackingFractionCoefficients(2,iComponent,jComponent) = SLCoefficientsMatrix(2,1) + SLCoefficientsMatrix(2,2) * &
      &     ijPotentialRange(iComponent,jComponent) + SLCoefficientsMatrix(2,3) * ijPotentialRange(iComponent,jComponent) * &
      &     ijPotentialRange(iComponent,jComponent) + SLCoefficientsMatrix(2,4) * ijPotentialRange(iComponent,jComponent) * &
      &     ijPotentialRange(iComponent,jComponent) * ijPotentialRange(iComponent,jComponent)
      EffPackingFractionCoefficients(3,iComponent,jComponent) = 0.D0
      EffPackingFractionCoefficientsSutherland(1,iComponent,jComponent) = SLCoefficientsMatrix(1,1) + 2.D0 * &
      &     SLCoefficientsMatrix(1,2) * ijPotentialRange(iComponent,jComponent) + 4.D0 * SLCoefficientsMatrix(1,3) * &
      &     ijPotentialRange(iComponent,jComponent) * ijPotentialRange(iComponent,jComponent) + 8.D0 * SLCoefficientsMatrix(1,4) * &
      &     ijPotentialRange(iComponent,jComponent) * ijPotentialRange(iComponent,jComponent) * &
      &     ijPotentialRange(iComponent,jComponent)
      EffPackingFractionCoefficientsSutherland(2,iComponent,jComponent) = SLCoefficientsMatrix(2,1) + 2.D0 * &
      &     SLCoefficientsMatrix(2,2) * ijPotentialRange(iComponent,jComponent) + 4.D0 * SLCoefficientsMatrix(2,3) * &
      &     ijPotentialRange(iComponent,jComponent) * ijPotentialRange(iComponent,jComponent) + 8.D0 * SLCoefficientsMatrix(2,4) * &
      &     ijPotentialRange(iComponent,jComponent) * ijPotentialRange(iComponent,jComponent) * &
      &     ijPotentialRange(iComponent,jComponent)
      EffPackingFractionCoefficientsSutherland(3,iComponent,jComponent) = 0.D0
    ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
      EffPackingFractionCoefficients(1,iComponent,jComponent) = YK1CoefficientsMatrix(1,1) + YK1CoefficientsMatrix(1,2) * &
      &     ( 1.D0 / ijPotentialRange(iComponent,jComponent) ) + YK1CoefficientsMatrix(1,3) * ( 1.D0 / &
      &     ( ijPotentialRange(iComponent,jComponent) * ijPotentialRange(iComponent,jComponent) ) )
      EffPackingFractionCoefficients(2,iComponent,jComponent) = YK1CoefficientsMatrix(2,1) + YK1CoefficientsMatrix(2,2) * &
      &     ( 1.D0 / ijPotentialRange(iComponent,jComponent) ) + YK1CoefficientsMatrix(2,3) * ( 1.D0 / &
      &     ( ijPotentialRange(iComponent,jComponent) * ijPotentialRange(iComponent,jComponent) ) )
      EffPackingFractionCoefficients(3,iComponent,jComponent) = 0.D0
      EffPackingFractionCoefficientsYukawa(1,iComponent,jComponent) = YK2CoefficientsMatrix(1,1) + YK2CoefficientsMatrix(1,2) * &
      &     ( 1.D0 / ijPotentialRange(iComponent,jComponent) ) + YK2CoefficientsMatrix(1,3) * ( 1.D0 / &
      &     ( ijPotentialRange(iComponent,jComponent) * ijPotentialRange(iComponent,jComponent) ) ) + YK2CoefficientsMatrix(1,4) * &
      &     ( 1.D0 / ( ijPotentialRange(iComponent,jComponent) * ijPotentialRange(iComponent,jComponent) * &
      &     ijPotentialRange(iComponent,jComponent) ) ) + YK2CoefficientsMatrix(1,5) * ( 1.D0 / &
      &     ( ijPotentialRange(iComponent,jComponent) * ijPotentialRange(iComponent,jComponent) * &
      &     ijPotentialRange(iComponent,jComponent) * ijPotentialRange(iComponent,jComponent) ) )
      EffPackingFractionCoefficientsYukawa(2,iComponent,jComponent) = YK2CoefficientsMatrix(2,1) + YK2CoefficientsMatrix(2,2) * &
      &     ( 1.D0 / ijPotentialRange(iComponent,jComponent) ) + YK2CoefficientsMatrix(2,3) * ( 1.D0 / &
      &     ( ijPotentialRange(iComponent,jComponent) * ijPotentialRange(iComponent,jComponent) ) ) + YK2CoefficientsMatrix(2,4) * &
      &     ( 1.D0 / ( ijPotentialRange(iComponent,jComponent) * ijPotentialRange(iComponent,jComponent) * &
      &     ijPotentialRange(iComponent,jComponent) ) ) + YK2CoefficientsMatrix(2,5) * ( 1.D0 / &
      &     ( ijPotentialRange(iComponent,jComponent) * ijPotentialRange(iComponent,jComponent) * &
      &     ijPotentialRange(iComponent,jComponent) * ijPotentialRange(iComponent,jComponent) ) )
      EffPackingFractionCoefficientsYukawa(3,iComponent,jComponent) = 0.D0
    ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square-well potential
      IF( EffPFractionMixingRule == 1 ) THEN ! Reduced density 3 mixing rule
        CSWAlphaCoefficients(1,1) = CSW1CoefficientsMatrix(1,1) + CSW1CoefficientsMatrix(2,1) * ( & 
        &     ijNonSphericity(iComponent,jComponent) - 1.D0 ) + CSW1CoefficientsMatrix(3,1) * ( &
        &     ijNonSphericity(iComponent,jComponent) - 1.D0 ) * ( ijNonSphericity(iComponent,jComponent) - 1.D0 )
        CSWAlphaCoefficients(1,2) = CSW1CoefficientsMatrix(1,2) + CSW1CoefficientsMatrix(2,2) * ( & 
        &     ijNonSphericity(iComponent,jComponent) - 1.D0 ) + CSW1CoefficientsMatrix(3,2) * ( &
        &     ijNonSphericity(iComponent,jComponent) - 1.D0 ) * ( ijNonSphericity(iComponent,jComponent) - 1.D0 )
        CSWAlphaCoefficients(1,3) = CSW1CoefficientsMatrix(1,3) + CSW1CoefficientsMatrix(2,3) * ( & 
        &     ijNonSphericity(iComponent,jComponent) - 1.D0 ) + CSW1CoefficientsMatrix(3,3) * ( &
        &     ijNonSphericity(iComponent,jComponent) - 1.D0 ) * ( ijNonSphericity(iComponent,jComponent) - 1.D0 )
        CSWAlphaCoefficients(1,4) = CSW1CoefficientsMatrix(1,4) + CSW1CoefficientsMatrix(2,4) * ( & 
        &     ijNonSphericity(iComponent,jComponent) - 1.D0 ) + CSW1CoefficientsMatrix(3,4) * ( &
        &     ijNonSphericity(iComponent,jComponent) - 1.D0 ) * ( ijNonSphericity(iComponent,jComponent) - 1.D0 )
        CSWAlphaCoefficients(1,5) = CSW1CoefficientsMatrix(1,5) + CSW1CoefficientsMatrix(2,5) * ( & 
        &     ijNonSphericity(iComponent,jComponent) - 1.D0 ) + CSW1CoefficientsMatrix(3,5) * ( &
        &     ijNonSphericity(iComponent,jComponent) - 1.D0 ) * ( ijNonSphericity(iComponent,jComponent) - 1.D0 )
        CSWAlphaCoefficients(1,6) = CSW1CoefficientsMatrix(1,6) + CSW1CoefficientsMatrix(2,6) * ( & 
        &     ijNonSphericity(iComponent,jComponent) - 1.D0 ) + CSW1CoefficientsMatrix(3,6) * ( &
        &     ijNonSphericity(iComponent,jComponent) - 1.D0 ) * ( ijNonSphericity(iComponent,jComponent) - 1.D0 )
        CSWAlphaCoefficients(1,7) = CSW1CoefficientsMatrix(1,7) + CSW1CoefficientsMatrix(2,7) * ( & 
        &     ijNonSphericity(iComponent,jComponent) - 1.D0 ) + CSW1CoefficientsMatrix(3,7) * ( &
        &     ijNonSphericity(iComponent,jComponent) - 1.D0 ) * ( ijNonSphericity(iComponent,jComponent) - 1.D0 )
        CSWAlphaCoefficients(1,8) = CSW1CoefficientsMatrix(1,8) + CSW1CoefficientsMatrix(2,8) * ( & 
        &     ijNonSphericity(iComponent,jComponent) - 1.D0 ) + CSW1CoefficientsMatrix(3,8) * ( &
        &     ijNonSphericity(iComponent,jComponent) - 1.D0 ) * ( ijNonSphericity(iComponent,jComponent) - 1.D0 )
        CSWAlphaCoefficients(1,9) = CSW1CoefficientsMatrix(1,9) + CSW1CoefficientsMatrix(2,9) * ( & 
        &     ijNonSphericity(iComponent,jComponent) - 1.D0 ) + CSW1CoefficientsMatrix(3,9) * ( &
        &     ijNonSphericity(iComponent,jComponent) - 1.D0 ) * ( ijNonSphericity(iComponent,jComponent) - 1.D0 )
      ELSE IF( EffPFractionMixingRule == 2 ) THEN ! One-fluid van der Waals mixing rule
        CSWAlphaCoefficients(1,1) = CSW1CoefficientsMatrix(1,1) + CSW1CoefficientsMatrix(2,1) * ( mNonSphericityMixture - 1.D0 ) + &
        &     CSW1CoefficientsMatrix(3,1) * ( mNonSphericityMixture - 1.D0 ) * ( mNonSphericityMixture - 1.D0 )
        CSWAlphaCoefficients(1,2) = CSW1CoefficientsMatrix(1,2) + CSW1CoefficientsMatrix(2,2) * ( mNonSphericityMixture - 1.D0 ) + &
        &     CSW1CoefficientsMatrix(3,2) * ( mNonSphericityMixture - 1.D0 ) * ( mNonSphericityMixture - 1.D0 )
        CSWAlphaCoefficients(1,3) = CSW1CoefficientsMatrix(1,3) + CSW1CoefficientsMatrix(2,3) * ( mNonSphericityMixture - 1.D0 ) + &
        &     CSW1CoefficientsMatrix(3,3) * ( mNonSphericityMixture - 1.D0 ) * ( mNonSphericityMixture - 1.D0 )
        CSWAlphaCoefficients(1,4) = CSW1CoefficientsMatrix(1,4) + CSW1CoefficientsMatrix(2,4) * ( mNonSphericityMixture - 1.D0 ) + &
        &     CSW1CoefficientsMatrix(3,4) * ( mNonSphericityMixture - 1.D0 ) * ( mNonSphericityMixture - 1.D0 )
        CSWAlphaCoefficients(1,5) = CSW1CoefficientsMatrix(1,5) + CSW1CoefficientsMatrix(2,5) * ( mNonSphericityMixture - 1.D0 ) + &
        &     CSW1CoefficientsMatrix(3,5) * ( mNonSphericityMixture - 1.D0 ) * ( mNonSphericityMixture - 1.D0 )
        CSWAlphaCoefficients(1,6) = CSW1CoefficientsMatrix(1,6) + CSW1CoefficientsMatrix(2,6) * ( mNonSphericityMixture - 1.D0 ) + &
        &     CSW1CoefficientsMatrix(3,6) * ( mNonSphericityMixture - 1.D0 ) * ( mNonSphericityMixture - 1.D0 )
        CSWAlphaCoefficients(1,7) = CSW1CoefficientsMatrix(1,7) + CSW1CoefficientsMatrix(2,7) * ( mNonSphericityMixture - 1.D0 ) + &
        &     CSW1CoefficientsMatrix(3,7) * ( mNonSphericityMixture - 1.D0 ) * ( mNonSphericityMixture - 1.D0 )
        CSWAlphaCoefficients(1,8) = CSW1CoefficientsMatrix(1,8) + CSW1CoefficientsMatrix(2,8) * ( mNonSphericityMixture - 1.D0 ) + &
        &     CSW1CoefficientsMatrix(3,8) * ( mNonSphericityMixture - 1.D0 ) * ( mNonSphericityMixture - 1.D0 )
        CSWAlphaCoefficients(1,9) = CSW1CoefficientsMatrix(1,9) + CSW1CoefficientsMatrix(2,9) * ( mNonSphericityMixture - 1.D0 ) + &
        &     CSW1CoefficientsMatrix(3,9) * ( mNonSphericityMixture - 1.D0 ) * ( mNonSphericityMixture - 1.D0 )
      END IF
      EffPackingFractionCoefficients(1,iComponent,jComponent) = CSWAlphaCoefficients(1,1) + CSWAlphaCoefficients(1,2) * ( 1.D0 + &
      &     ijPotentialRange(iComponent,jComponent) ) + CSWAlphaCoefficients(1,3) * ( 1.D0 + &
      &     ijPotentialRange(iComponent,jComponent) ) * ( 1.D0 + ijPotentialRange(iComponent,jComponent) )
      EffPackingFractionCoefficients(2,iComponent,jComponent) = CSWAlphaCoefficients(1,4) + CSWAlphaCoefficients(1,5) * ( 1.D0 + &
      &     ijPotentialRange(iComponent,jComponent) ) + CSWAlphaCoefficients(1,6) * ( 1.D0 + &
      &     ijPotentialRange(iComponent,jComponent) ) * ( 1.D0 + ijPotentialRange(iComponent,jComponent) )
      EffPackingFractionCoefficients(3,iComponent,jComponent) = CSWAlphaCoefficients(1,7) + CSWAlphaCoefficients(1,8) * ( 1.D0 + &
      &     ijPotentialRange(iComponent,jComponent) ) + CSWAlphaCoefficients(1,9) * ( 1.D0 + &
      &     ijPotentialRange(iComponent,jComponent) ) * ( 1.D0 + ijPotentialRange(iComponent,jComponent) )
      IF( EffPFractionMixingRule == 1 ) THEN ! Reduced density 3 mixing rule
        CSWAlphaCoefficients(2,1) = CSW2CoefficientsMatrix(1,1) + CSW2CoefficientsMatrix(2,1) * ( &
        &     ijNonSphericity(iComponent,jComponent) - 1.D0 ) + CSW2CoefficientsMatrix(3,1) * ( &
        &     ijNonSphericity(iComponent,jComponent) - 1.D0 ) * ( ijNonSphericity(iComponent,jComponent) - 1.D0 )
        CSWAlphaCoefficients(2,2) = CSW2CoefficientsMatrix(1,2) + CSW2CoefficientsMatrix(2,2) * ( &
        &     ijNonSphericity(iComponent,jComponent) - 1.D0 ) + CSW2CoefficientsMatrix(3,2) * ( &
        &     ijNonSphericity(iComponent,jComponent) - 1.D0 ) * ( ijNonSphericity(iComponent,jComponent) - 1.D0 )
        CSWAlphaCoefficients(2,3) = CSW2CoefficientsMatrix(1,3) + CSW2CoefficientsMatrix(2,3) * ( &
        &     ijNonSphericity(iComponent,jComponent) - 1.D0 ) + CSW2CoefficientsMatrix(3,3) * ( &
        &     ijNonSphericity(iComponent,jComponent) - 1.D0 ) * ( ijNonSphericity(iComponent,jComponent) - 1.D0 )
        CSWAlphaCoefficients(2,4) = CSW2CoefficientsMatrix(1,4) + CSW2CoefficientsMatrix(2,4) * ( &
        &     ijNonSphericity(iComponent,jComponent) - 1.D0 ) + CSW2CoefficientsMatrix(3,4) * ( &
        &     ijNonSphericity(iComponent,jComponent) - 1.D0 ) * ( ijNonSphericity(iComponent,jComponent) - 1.D0 )
        CSWAlphaCoefficients(2,5) = CSW2CoefficientsMatrix(1,5) + CSW2CoefficientsMatrix(2,5) * ( &
        &     ijNonSphericity(iComponent,jComponent) - 1.D0 ) + CSW2CoefficientsMatrix(3,5) * ( &
        &     ijNonSphericity(iComponent,jComponent) - 1.D0 ) * ( ijNonSphericity(iComponent,jComponent) - 1.D0 )
        CSWAlphaCoefficients(2,6) = CSW2CoefficientsMatrix(1,6) + CSW2CoefficientsMatrix(2,6) * ( &
        &     ijNonSphericity(iComponent,jComponent) - 1.D0 ) + CSW2CoefficientsMatrix(3,6) * ( &
        &     ijNonSphericity(iComponent,jComponent) - 1.D0 ) * ( ijNonSphericity(iComponent,jComponent) - 1.D0 )
        CSWAlphaCoefficients(2,7) = CSW2CoefficientsMatrix(1,7) + CSW2CoefficientsMatrix(2,7) * ( &
        &     ijNonSphericity(iComponent,jComponent) - 1.D0 ) + CSW2CoefficientsMatrix(3,7) * ( &
        &     ijNonSphericity(iComponent,jComponent) - 1.D0 ) * ( ijNonSphericity(iComponent,jComponent) - 1.D0 )
        CSWAlphaCoefficients(2,8) = CSW2CoefficientsMatrix(1,8) + CSW2CoefficientsMatrix(2,8) * ( &
        &     ijNonSphericity(iComponent,jComponent) - 1.D0 ) + CSW2CoefficientsMatrix(3,8) * ( &
        &     ijNonSphericity(iComponent,jComponent) - 1.D0 ) * ( ijNonSphericity(iComponent,jComponent) - 1.D0 )
        CSWAlphaCoefficients(2,9) = CSW2CoefficientsMatrix(1,9) + CSW2CoefficientsMatrix(2,9) * ( &
        &     ijNonSphericity(iComponent,jComponent) - 1.D0 ) + CSW2CoefficientsMatrix(3,9) * ( &
        &     ijNonSphericity(iComponent,jComponent) - 1.D0 ) * ( ijNonSphericity(iComponent,jComponent) - 1.D0 )
      ELSE IF( EffPFractionMixingRule == 2 ) THEN ! One-fluid van der Waals mixing rule
        CSWAlphaCoefficients(2,1) = CSW2CoefficientsMatrix(1,1) + CSW2CoefficientsMatrix(2,1) * ( mNonSphericityMixture - 1.D0 ) + &
        &     CSW2CoefficientsMatrix(3,1) * ( mNonSphericityMixture - 1.D0 ) * ( mNonSphericityMixture - 1.D0 )
        CSWAlphaCoefficients(2,2) = CSW2CoefficientsMatrix(1,2) + CSW2CoefficientsMatrix(2,2) * ( mNonSphericityMixture - 1.D0 ) + &
        &     CSW2CoefficientsMatrix(3,2) * ( mNonSphericityMixture - 1.D0 ) * ( mNonSphericityMixture - 1.D0 )
        CSWAlphaCoefficients(2,3) = CSW2CoefficientsMatrix(1,3) + CSW2CoefficientsMatrix(2,3) * ( mNonSphericityMixture - 1.D0 ) + &
        &     CSW2CoefficientsMatrix(3,3) * ( mNonSphericityMixture - 1.D0 ) * ( mNonSphericityMixture - 1.D0 )
        CSWAlphaCoefficients(2,4) = CSW2CoefficientsMatrix(1,4) + CSW2CoefficientsMatrix(2,4) * ( mNonSphericityMixture - 1.D0 ) + &
        &     CSW2CoefficientsMatrix(3,4) * ( mNonSphericityMixture - 1.D0 ) * ( mNonSphericityMixture - 1.D0 )
        CSWAlphaCoefficients(2,5) = CSW2CoefficientsMatrix(1,5) + CSW2CoefficientsMatrix(2,5) * ( mNonSphericityMixture - 1.D0 ) + &
        &     CSW2CoefficientsMatrix(3,5) * ( mNonSphericityMixture - 1.D0 ) * ( mNonSphericityMixture - 1.D0 )
        CSWAlphaCoefficients(2,6) = CSW2CoefficientsMatrix(1,6) + CSW2CoefficientsMatrix(2,6) * ( mNonSphericityMixture - 1.D0 ) + &
        &     CSW2CoefficientsMatrix(3,6) * ( mNonSphericityMixture - 1.D0 ) * ( mNonSphericityMixture - 1.D0 )
        CSWAlphaCoefficients(2,7) = CSW2CoefficientsMatrix(1,7) + CSW2CoefficientsMatrix(2,7) * ( mNonSphericityMixture - 1.D0 ) + &
        &     CSW2CoefficientsMatrix(3,7) * ( mNonSphericityMixture - 1.D0 ) * ( mNonSphericityMixture - 1.D0 )
        CSWAlphaCoefficients(2,8) = CSW2CoefficientsMatrix(1,8) + CSW2CoefficientsMatrix(2,8) * ( mNonSphericityMixture - 1.D0 ) + &
        &     CSW2CoefficientsMatrix(3,8) * ( mNonSphericityMixture - 1.D0 ) * ( mNonSphericityMixture - 1.D0 )
        CSWAlphaCoefficients(2,9) = CSW2CoefficientsMatrix(1,9) + CSW2CoefficientsMatrix(2,9) * ( mNonSphericityMixture - 1.D0 ) + &
        &     CSW2CoefficientsMatrix(3,9) * ( mNonSphericityMixture - 1.D0 ) * ( mNonSphericityMixture - 1.D0 )
      END IF
      EffPackingFractionCoefficientsCSW(1,iComponent,jComponent) = CSWAlphaCoefficients(2,1) + CSWAlphaCoefficients(2,2) * &
      &     ( 1.D0 + ijPotentialRange(iComponent,jComponent) ) + CSWAlphaCoefficients(2,3) * ( 1.D0 + &
      &     ijPotentialRange(iComponent,jComponent) ) * ( 1.D0 + ijPotentialRange(iComponent,jComponent) )
      EffPackingFractionCoefficientsCSW(2,iComponent,jComponent) = CSWAlphaCoefficients(2,4) + CSWAlphaCoefficients(2,5) * &
      &     ( 1.D0 + ijPotentialRange(iComponent,jComponent) ) + CSWAlphaCoefficients(2,6) * ( 1.D0 + &
      &     ijPotentialRange(iComponent,jComponent) ) * ( 1.D0 + ijPotentialRange(iComponent,jComponent) )
      EffPackingFractionCoefficientsCSW(3,iComponent,jComponent) = CSWAlphaCoefficients(2,7) + CSWAlphaCoefficients(2,8) * &
      &     ( 1.D0 + ijPotentialRange(iComponent,jComponent) ) + CSWAlphaCoefficients(2,9) * ( 1.D0 + &
      &     ijPotentialRange(iComponent,jComponent) ) * ( 1.D0 + ijPotentialRange(iComponent,jComponent) )
    END IF
    ! Effective packing fraction [unitless]
    IF( EffPFractionMixingRule == 1 ) THEN ! Reduced density 3 mixing rule
      IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
        EffPackingFraction(iComponent,jComponent) = EffPackingFractionCoefficients(1,iComponent,jComponent) * rDensity(3) + &
        &     EffPackingFractionCoefficients(2,iComponent,jComponent) * rDensity(3) * rDensity(3) + &
        &     EffPackingFractionCoefficients(3,iComponent,jComponent) * rDensity(3) * rDensity(3) * rDensity(3)
      ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
        EffPackingFraction(iComponent,jComponent) = EffPackingFractionCoefficients(1,iComponent,jComponent) * rDensity(3) + &
        &     EffPackingFractionCoefficients(2,iComponent,jComponent) * rDensity(3) * rDensity(3)
        EffPackingFractionSutherland(iComponent,jComponent) = EffPackingFractionCoefficientsSutherland(1,iComponent,jComponent) * &
        &     rDensity(3) + EffPackingFractionCoefficientsSutherland(2,iComponent,jComponent) * rDensity(3) * rDensity(3)
      ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
        EffPackingFraction(iComponent,jComponent) = EffPackingFractionCoefficients(1,iComponent,jComponent) * rDensity(3) + &
        &     EffPackingFractionCoefficients(2,iComponent,jComponent) * rDensity(3) * rDensity(3)
        EffPackingFractionYukawa(iComponent,jComponent) = EffPackingFractionCoefficientsYukawa(1,iComponent,jComponent) * &
        &     rDensity(3) + EffPackingFractionCoefficientsYukawa(2,iComponent,jComponent) * rDensity(3) * rDensity(3)
      ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square-well potential
        EffPackingFraction(iComponent,jComponent) = EffPackingFractionCoefficients(1,iComponent,jComponent) * rDensity(3) + &
        &     EffPackingFractionCoefficients(2,iComponent,jComponent) * rDensity(3) * rDensity(3) + &
        &     EffPackingFractionCoefficients(3,iComponent,jComponent) * rDensity(3) * rDensity(3) * rDensity(3)
        EffPackingFractionCSW(iComponent,jComponent) = EffPackingFractionCoefficientsCSW(1,iComponent,jComponent) * rDensity(3) + &
        &     EffPackingFractionCoefficientsCSW(2,iComponent,jComponent) * rDensity(3) * rDensity(3) + &
        &     EffPackingFractionCoefficientsCSW(3,iComponent,jComponent) * rDensity(3) * rDensity(3) * rDensity(3)
      END IF
    ELSE IF( EffPFractionMixingRule == 2 ) THEN ! One-fluid van der Waals mixing rule
      IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
        EffPackingFraction(iComponent,jComponent) = EffPackingFractionCoefficients(1,iComponent,jComponent) * rDensityMixture + &
        &     EffPackingFractionCoefficients(2,iComponent,jComponent) * rDensityMixture * rDensityMixture + &
        &     EffPackingFractionCoefficients(3,iComponent,jComponent) * rDensityMixture * rDensityMixture * rDensityMixture
      ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
        EffPackingFraction(iComponent,jComponent) = EffPackingFractionCoefficients(1,iComponent,jComponent) * rDensityMixture + &
        &     EffPackingFractionCoefficients(2,iComponent,jComponent) * rDensityMixture * rDensityMixture
        EffPackingFractionSutherland(iComponent,jComponent) = EffPackingFractionCoefficientsSutherland(1,iComponent,jComponent) * &
        &     rDensityMixture + EffPackingFractionCoefficientsSutherland(2,iComponent,jComponent) * rDensityMixture * &
        &     rDensityMixture
      ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
        EffPackingFraction(iComponent,jComponent) = EffPackingFractionCoefficients(1,iComponent,jComponent) * rDensityMixture + &
        &     EffPackingFractionCoefficients(2,iComponent,jComponent) * rDensityMixture * rDensityMixture
        EffPackingFractionYukawa(iComponent,jComponent) = EffPackingFractionCoefficientsYukawa(1,iComponent,jComponent) * &
        &     rDensityMixture + EffPackingFractionCoefficientsYukawa(2,iComponent,jComponent) * rDensityMixture * rDensityMixture
      ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square-well potential mNonSphericityMixture
        EffPackingFraction(iComponent,jComponent) = EffPackingFractionCoefficients(1,iComponent,jComponent) * rDensityMixture + &
        &     EffPackingFractionCoefficients(2,iComponent,jComponent) * rDensityMixture * rDensityMixture + &
        &     EffPackingFractionCoefficients(3,iComponent,jComponent) * rDensityMixture * rDensityMixture * rDensityMixture
        EffPackingFractionCSW(iComponent,jComponent) = EffPackingFractionCoefficientsCSW(1,iComponent,jComponent) * &
        &     rDensityMixture + EffPackingFractionCoefficientsCSW(2,iComponent,jComponent) * rDensityMixture * rDensityMixture + &
        &     EffPackingFractionCoefficientsCSW(3,iComponent,jComponent) * rDensityMixture * rDensityMixture * rDensityMixture
      END IF
    END IF
    ! Contact radial distribution function [unitless]
    cRadialDistributionFunction(iComponent,jComponent) = ( 1.D0 / ( 1.D0 - rDensity(3) ) ) + ( ( 3.D0 * rDensity(2) * &
    &     cAverageDiameterRelationship(iComponent,jComponent) ) / ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) ) + &
    &     ( ( 2.D0 * cAverageDiameterRelationship(iComponent,jComponent) * cAverageDiameterRelationship(iComponent,jComponent) * &
    &     rDensity(2) * rDensity(2) ) / ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) )
    ! Contact radial distribution function for an effective packing fraction [unitless]
    IF( EffPFractionMixingRule == 1 ) THEN ! Reduced density 3 mixing rule
      IF( .NOT. PotentialTypeLogical(4) ) THEN
        cEffectiveRadialDistributionFunction(iComponent,jComponent) = ( 1.D0 / ( 1.D0 - EffPackingFraction(iComponent,jComponent) &
        &     ) ) + ( ( 3.D0 * cAverageDiameterRelationship(iComponent,jComponent) * EffPackingFraction(iComponent,jComponent) * &
        &     AuxDiameterRelationship ) / ( ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) ) ) + ( ( 2.D0 * cAverageDiameterRelationship(iComponent,jComponent) * &
        &     cAverageDiameterRelationship(iComponent,jComponent) * EffPackingFraction(iComponent,jComponent) * &
        &     AuxDiameterRelationship * EffPackingFraction(iComponent,jComponent) * AuxDiameterRelationship ) / ( ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) ) )
      ELSE
        geoAux = 3.D0 * cMolecularVolume(iComponent) + 2.D0 * cSurfaceArea(iComponent) * cCurvatureRadius(jComponent) + &
        &     cCurvatureRadius(iComponent) * cSurfaceArea(jComponent)
        cEffectiveRadialDistributionFunction(iComponent,jComponent) = 1.D0 / ( 1.D0 - EffPackingFraction(iComponent,jComponent) &
        &     ) + ( 3.D0 * cMolecularVolume(iComponent) * cSurfaceArea(jComponent) * EffPackingFraction(iComponent,jComponent) * &
        &     ( sumRadius - sumVolumeSurface ) + sumArea * EffPackingFraction(iComponent,jComponent) * ( 3.D0 * &
        &     cMolecularVolume(iComponent) * cVolumeSurface(jComponent) + cCurvatureRadius(iComponent) * &
        &     cCurvatureRadius(iComponent) * cSurfaceArea(jComponent) ) ) / ( ( 1.D0 - EffPackingFraction(iComponent,jComponent) &
        &     ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * geoAux * sumVolume ) + 6.D0 * &
        &     cMolecularVolume(iComponent) * cSurfaceArea(jComponent) * sumSqRadius * EffPackingFraction(iComponent,jComponent) &
        &     * sumArea * EffPackingFraction(iComponent,jComponent) / ( 9.D0 * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * geoAux * sumVolume * sumVolume )
      END IF
      IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
        cEffectiveRadialDistributionFunctionSL(iComponent,jComponent) = ( 1.D0 / ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) ) + ( ( 3.D0 * &
        &     cAverageDiameterRelationship(iComponent,jComponent) * EffPackingFractionSutherland(iComponent,jComponent) * &
        &     AuxDiameterRelationship ) / ( ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) ) ) + ( ( 2.D0 * &
        &     cAverageDiameterRelationship(iComponent,jComponent) * cAverageDiameterRelationship(iComponent,jComponent) * &
        &     EffPackingFractionSutherland(iComponent,jComponent) * AuxDiameterRelationship * &
        &     EffPackingFractionSutherland(iComponent,jComponent) * AuxDiameterRelationship ) / ( ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) &
        &     ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) ) )
      ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
        cEffectiveRadialDistributionFunctionYK(iComponent,jComponent) = ( 1.D0 / ( 1.D0 - &
        &     EffPackingFractionYukawa(iComponent,jComponent) ) ) + ( ( 3.D0 * cAverageDiameterRelationship(iComponent,jComponent) &
        &     * EffPackingFractionYukawa(iComponent,jComponent) * AuxDiameterRelationship ) / ( ( 1.D0 - &
        &     EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) ) ) + &
        &     ( ( 2.D0 * cAverageDiameterRelationship(iComponent,jComponent) * cAverageDiameterRelationship(iComponent,jComponent) &
        &     * EffPackingFractionYukawa(iComponent,jComponent) * AuxDiameterRelationship * &
        &     EffPackingFractionYukawa(iComponent,jComponent) * AuxDiameterRelationship ) / ( ( 1.D0 - &
        &     EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * &
        &     ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) ) )
      ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square-well potential
        geoAux = 3.D0 * cMolecularVolume(iComponent) + 2.D0 * cSurfaceArea(iComponent) * cCurvatureRadius(jComponent) + &
        &     cCurvatureRadius(iComponent) * cSurfaceArea(jComponent)
        cEffectiveRadialDistributionFunctionCSW(iComponent,jComponent) = 1.D0 / ( 1.D0 - &
        &     EffPackingFractionCSW(iComponent,jComponent) ) + ( 3.D0 * cMolecularVolume(iComponent) * cSurfaceArea(jComponent) &
        &     * EffPackingFractionCSW(iComponent,jComponent) * ( sumRadius - sumVolumeSurface ) + sumArea * &
        &     EffPackingFractionCSW(iComponent,jComponent) * ( 3.D0 * cMolecularVolume(iComponent) * cVolumeSurface(jComponent) &
        &     + cCurvatureRadius(iComponent) * cCurvatureRadius(iComponent) * cSurfaceArea(jComponent) ) ) / ( ( 1.D0 - &
        &     EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * geoAux &
        &     * sumVolume ) + 6.D0 * cMolecularVolume(iComponent) * cSurfaceArea(jComponent) * sumSqRadius * &
        &     EffPackingFractionCSW(iComponent,jComponent) * sumArea * EffPackingFractionCSW(iComponent,jComponent) / ( 9.D0 * &
        &     ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) &
        &     * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * geoAux * sumVolume * sumVolume )
      END IF
    ELSE IF( EffPFractionMixingRule == 2 ) THEN ! One-fluid van der Waals mixing rule
      IF( .NOT. PotentialTypeLogical(4) ) THEN
        cEffectiveRadialDistributionFunction(iComponent,jComponent) = ( 1.D0 - 0.5D0 * EffPackingFraction(iComponent,jComponent) ) &
        &     / ( ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( &
        &     1.D0 - EffPackingFraction(iComponent,jComponent) ) )
      ELSE
        cEffectiveRadialDistributionFunction(iComponent,jComponent) = 1.D0 / ( 1.D0 - EffPackingFraction(iComponent,jComponent) &
        &     ) + 3.D0 * EffPackingFraction(iComponent,jComponent) * mNonSphericityMixture * ( 1.D0 + mNonSphericityMixture ) / &
        &     ( ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( &
        &     1.D0 + 3.D0 * mNonSphericityMixture ) ) + 2.D0 * EffPackingFraction(iComponent,jComponent) * &
        &     EffPackingFraction(iComponent,jComponent) *  mNonSphericityMixture * mNonSphericityMixture / ( ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 + 3.D0 * mNonSphericityMixture ) )
      END IF
      IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
        cEffectiveRadialDistributionFunctionSL(iComponent,jComponent) = ( 1.D0 - 0.5D0 * &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) / ( ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) &
        &     ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) )
      ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
        cEffectiveRadialDistributionFunctionYK(iComponent,jComponent) = ( 1.D0 - 0.5D0 * &
        &     EffPackingFractionYukawa(iComponent,jComponent) ) / ( ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * &
        &     ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionYukawa(iComponent,jComponent) ) )
      ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square-well potential
        cEffectiveRadialDistributionFunctionCSW(iComponent,jComponent) = 1.D0 / ( 1.D0 - &
        &     EffPackingFractionCSW(iComponent,jComponent) ) + 3.D0 * EffPackingFractionCSW(iComponent,jComponent) * &
        &     mNonSphericityMixture * ( 1.D0 + mNonSphericityMixture ) / ( ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) &
        &     ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 + 3.D0 * mNonSphericityMixture ) ) + 2.D0 * &
        &     EffPackingFractionCSW(iComponent,jComponent) * EffPackingFractionCSW(iComponent,jComponent) * &
        &     mNonSphericityMixture * mNonSphericityMixture / ( ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 &
        &     - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * &
        &     ( 1.D0 + 3.D0 * mNonSphericityMixture ) )
      END IF
    END IF
    ! Mean-attractive energy between components [K]
    IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
      cMeanAttractiveEnergy(iComponent,jComponent) = - 2.D0 * cPi * rNumberDensity * ijaWellDepth(iComponent,jComponent) * &
      &     ( ijaDiameterSphereCubic(iComponent,jComponent) / 3.D0 ) * ( ijPotentialRangeCubic(iComponent,jComponent) - 1.D0 ) * &
      &     cEffectiveRadialDistributionFunction(iComponent,jComponent)
    ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
      cMeanAttractiveEnergy(iComponent,jComponent) = - 2.D0 * cPi * rNumberDensity * ijaWellDepth(iComponent,jComponent) * &
      &     ijaDiameterSphereCubic(iComponent,jComponent) * ( 1.D0 / ( ijPotentialRange(iComponent,jComponent) - 3.D0 ) ) * &
      &     cEffectiveRadialDistributionFunction(iComponent,jComponent)
      cMeanAttractiveEnergySutherland(iComponent,jComponent) = - 2.D0 * cPi * rNumberDensity * ijaWellDepth(iComponent,jComponent) &
      &     * ijaDiameterSphereCubic(iComponent,jComponent) * ( 1.D0 / ( 2.D0 * ijPotentialRange(iComponent,jComponent) - 3.D0 ) ) &
      &     * cEffectiveRadialDistributionFunctionSL(iComponent,jComponent)
    ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
      cMeanAttractiveEnergy(iComponent,jComponent) = - 2.D0 * cPi * rNumberDensity * ijaWellDepth(iComponent,jComponent) * &
      &     ijaDiameterSphereCubic(iComponent,jComponent) * ( ( 1.D0 / ijPotentialRange(iComponent,jComponent) ) + ( 1.D0 / &
      &     ( ijPotentialRange(iComponent,jComponent) * ijPotentialRange(iComponent,jComponent) ) ) ) * &
      &     cEffectiveRadialDistributionFunction(iComponent,jComponent)
      cMeanAttractiveEnergyYukawa(iComponent,jComponent) = - cPi * rNumberDensity * ijaWellDepth(iComponent,jComponent) * &
      &     ijaDiameterSphereCubic(iComponent,jComponent) * ( 1.D0 / ijPotentialRange(iComponent,jComponent) ) * &
      &     cEffectiveRadialDistributionFunctionYK(iComponent,jComponent)
    ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square-well potential
      cMeanAttractiveEnergy(iComponent,jComponent) = - rNumberDensity * ijaWellDepth(iComponent,jComponent) * &
      &     cEffectiveRadialDistributionFunction(iComponent,jComponent) * ( ijSecondVirialCoefficientField(iComponent,jComponent) &
      &     - ijSecondVirialCoefficient(iComponent,jComponent) )
      cMeanAttractiveEnergyCSW(iComponent,jComponent) = - rNumberDensity * ijaWellDepth(iComponent,jComponent) * &
      &     cEffectiveRadialDistributionFunctionCSW(iComponent,jComponent) * &
      &     ( ijSecondVirialCoefficientField(iComponent,jComponent) - ijSecondVirialCoefficient(iComponent,jComponent) )
    END IF
  END DO
END DO

! Mixing rule for the first-order perturbation coefficient [K]
CALL Mixing_Rules( mFraction, cMeanAttractiveEnergy, pFirstOrderCoefficient )

! First-order perturbation contribution to the Helmholtz free energy (molar basis) [J / mol]
mMeanAttractiveHelmholtzFreeEnergy = pFirstOrderCoefficient * cUniversalGas ! Proven units

! ************************************************************************************************ !
! SECOND-ORDER PERTURBATION CONTRIBUTION                                                           !
! ************************************************************************************************ !

! Hard-sphere isothermal compressibility given by the Percus-Yevick expression [unitless]
HSIsothermalCompressibility = ( rDensity(0) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - &
&     rDensity(3) ) ) / ( ( rDensity(0) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) + ( 6.D0 * rDensity(1) * rDensity(2) * &
&     ( 1.D0 - rDensity(3) ) ) + ( 9.D0 * rDensity(2) * rDensity(2) * rDensity(2) ) )

! Hard convex-body isothermal compressibility given by the Boublík expression [unitless]
IF( PYHCBCorrectionLogical ) THEN
  HCBIsothermalCompressibility = ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) &
  &     / ( 1.D0 + 2.D0 * rDensity(3) * ( sumRadius * sumArea / sumVolume - 1.D0 ) + rDensity(3) * rDensity(3) * ( 1.D0 - 2.D0 * &
  &     sumRadius * sumArea / sumVolume + sumSqRadius * sumArea * sumArea / ( sumVolume * sumVolume ) ) )
ELSE
  HCBIsothermalCompressibility = ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) &
  &     / ( 1.D0 + 2.D0 * rDensity(3) * ( sumRadius * sumArea / sumVolume - 1.D0 ) + rDensity(3) * rDensity(3) * ( 1.D0 - 2.D0 * &
  &     sumRadius * sumArea / sumVolume + sumSqRadius * sumArea * sumArea / ( sumVolume * sumVolume ) ) + ( sumSqRadius * sumArea &
  &     * sumArea / ( 9.D0 * sumVolume * sumVolume ) ) * ( rDensity(3) * rDensity(3) * rDensity(3) * rDensity(3) - 4.D0 * &
  &     rDensity(3) * rDensity(3) * rDensity(3) ) )
END IF

! Mean-attractive energy fluctuation between components
DO iComponent = 1, nComponents
  DO jComponent = 1, nComponents
    ! First derivative of the effective packing fraction with respect to the density [Å³]
    IF( EffPFractionMixingRule == 1 ) THEN ! Reduced density 3 mixing rule
      IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
        dEffPackingFraction_dDensity(iComponent,jComponent) = EffPackingFractionCoefficients(1,iComponent,jComponent) * &
        &    rDensity(3) + 2.D0 * EffPackingFractionCoefficients(2,iComponent,jComponent) * rDensity(3) * rDensity(3) + 3.D0 * &
        &    rDensity(3) * rDensity(3) * rDensity(3) * EffPackingFractionCoefficients(3,iComponent,jComponent)
      ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
        dEffPackingFraction_dDensity(iComponent,jComponent) = EffPackingFractionCoefficients(1,iComponent,jComponent) * &
        &    rDensity(3) + 2.D0 * EffPackingFractionCoefficients(2,iComponent,jComponent) * rDensity(3) * rDensity(3)
        dEffPackingFractionSutherland_dDensity(iComponent,jComponent) = rDensity(3) * &
        &    EffPackingFractionCoefficientsSutherland(1,iComponent,jComponent) + 2.D0 * &
        &    EffPackingFractionCoefficientsSutherland(2,iComponent,jComponent) * rDensity(3) * rDensity(3)
      ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
        dEffPackingFraction_dDensity(iComponent,jComponent) = EffPackingFractionCoefficients(1,iComponent,jComponent) * &
        &    rDensity(3) + 2.D0 * EffPackingFractionCoefficients(2,iComponent,jComponent) * rDensity(3) * rDensity(3)
        dEffPackingFractionYukawa_dDensity(iComponent,jComponent) = EffPackingFractionCoefficientsYukawa(1,iComponent,jComponent) &
        &    * rDensity(3) + 2.D0 * EffPackingFractionCoefficientsYukawa(2,iComponent,jComponent) * rDensity(3) * rDensity(3)
      ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square-well potential
        dEffPackingFraction_dDensity(iComponent,jComponent) = EffPackingFractionCoefficients(1,iComponent,jComponent) * &
        &    rDensity(3) + 2.D0 * EffPackingFractionCoefficients(2,iComponent,jComponent) * rDensity(3) * rDensity(3) + 3.D0 * &
        &    rDensity(3) * rDensity(3) * rDensity(3) * EffPackingFractionCoefficients(3,iComponent,jComponent)
        dEffPackingFractionCSW_dDensity(iComponent,jComponent) = EffPackingFractionCoefficientsCSW(1,iComponent,jComponent) * &
        &    rDensity(3) + 2.D0 * rDensity(3) * rDensity(3) * EffPackingFractionCoefficientsCSW(2,iComponent,jComponent) + 3.D0 * &
        &    EffPackingFractionCoefficientsCSW(3,iComponent,jComponent) * rDensity(3) * rDensity(3) * rDensity(3)
      END IF
    ELSE IF( EffPFractionMixingRule == 2 ) THEN ! One-fluid van der Waals mixing rule
      IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
        dEffPackingFraction_dDensity(iComponent,jComponent) = EffPackingFractionCoefficients(1,iComponent,jComponent) * &
        &     rDensityMixture + 2.D0 * EffPackingFractionCoefficients(2,iComponent,jComponent) * rDensityMixture * rDensityMixture &
        &     + 3.D0 * rDensityMixture * rDensityMixture * rDensityMixture * EffPackingFractionCoefficients(3,iComponent,jComponent)
      ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
        dEffPackingFraction_dDensity(iComponent,jComponent) = EffPackingFractionCoefficients(1,iComponent,jComponent) * &
        &     rDensityMixture + 2.D0 * EffPackingFractionCoefficients(2,iComponent,jComponent) * rDensityMixture * rDensityMixture
        dEffPackingFractionSutherland_dDensity(iComponent,jComponent) = rDensityMixture * &
        &     EffPackingFractionCoefficientsSutherland(1,iComponent,jComponent) + 2.D0 * &
        &     EffPackingFractionCoefficientsSutherland(2,iComponent,jComponent) * rDensityMixture * rDensityMixture
      ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
        dEffPackingFraction_dDensity(iComponent,jComponent) = EffPackingFractionCoefficients(1,iComponent,jComponent) * &
        &     rDensityMixture + 2.D0 * EffPackingFractionCoefficients(2,iComponent,jComponent) * rDensityMixture * rDensityMixture
        dEffPackingFractionYukawa_dDensity(iComponent,jComponent) = EffPackingFractionCoefficientsYukawa(1,iComponent,jComponent) &
        &     * rDensityMixture + 2.D0 * EffPackingFractionCoefficientsYukawa(2,iComponent,jComponent) * rDensityMixture * &
        &     rDensityMixture
      ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square-well potential
        dEffPackingFraction_dDensity(iComponent,jComponent) = EffPackingFractionCoefficients(1,iComponent,jComponent) * &
        &    rDensityMixture + 2.D0 * EffPackingFractionCoefficients(2,iComponent,jComponent) * rDensityMixture * rDensityMixture &
        &    + 3.D0 * rDensityMixture * rDensityMixture * rDensityMixture * EffPackingFractionCoefficients(3,iComponent,jComponent)
        dEffPackingFractionCSW_dDensity(iComponent,jComponent) = EffPackingFractionCoefficientsCSW(1,iComponent,jComponent) * &
        &    rDensityMixture + 2.D0 * EffPackingFractionCoefficientsCSW(2,iComponent,jComponent) * rDensityMixture * &
        &    rDensityMixture + 3.D0 * EffPackingFractionCoefficientsCSW(3,iComponent,jComponent) * rDensityMixture * &
        &    rDensityMixture * rDensityMixture
      END IF
    END IF
    dEffPackingFraction_dDensity(iComponent,jComponent) = dEffPackingFraction_dDensity(iComponent,jComponent) / rNumberDensity
    IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
      dEffPackingFractionSutherland_dDensity(iComponent,jComponent) = &
      &     dEffPackingFractionSutherland_dDensity(iComponent,jComponent) / rNumberDensity
    ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
      dEffPackingFractionYukawa_dDensity(iComponent,jComponent) = dEffPackingFractionYukawa_dDensity(iComponent,jComponent) / &
      &     rNumberDensity
    ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square-well potential
      dEffPackingFractionCSW_dDensity(iComponent,jComponent) = dEffPackingFractionCSW_dDensity(iComponent,jComponent) / &
      &     rNumberDensity
    END IF
    ! First derivative of the radial distribution function with respect to the density [Å³]
    IF( EffPFractionMixingRule == 1 ) THEN ! Reduced density 3 mixing rule
      IF( .NOT. PotentialTypeLogical(4) ) THEN
        dEffRadialDistributionFunct_dDensity(iComponent,jComponent) = ( dEffPackingFraction_dDensity(iComponent,jComponent) / ( ( &
        &     1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) ) ) + ( ( &
        &     3.D0 * cAverageDiameterRelationship(iComponent,jComponent) * ( dEffPackingFraction_dDensity(iComponent,jComponent) * &
        &     ( 1.D0 + EffPackingFraction(iComponent,jComponent) ) ) * AuxDiameterRelationship ) / ( ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) ) ) + ( ( 2.D0 * cAverageDiameterRelationship(iComponent,jComponent) * &
        &     cAverageDiameterRelationship(iComponent,jComponent) * AuxDiameterRelationship * AuxDiameterRelationship * ( ( 2.D0 + &
        &     EffPackingFraction(iComponent,jComponent) ) * EffPackingFraction(iComponent,jComponent) * &
        &     dEffPackingFraction_dDensity(iComponent,jComponent) ) ) / ( ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * &
        &     ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 &
        &     - EffPackingFraction(iComponent,jComponent) ) ) )
      ELSE
        geoAux = 3.D0 * cMolecularVolume(iComponent) + 2.D0 * cSurfaceArea(iComponent) * cCurvatureRadius(jComponent) + &
        &     cCurvatureRadius(iComponent) * cSurfaceArea(jComponent)
        dAux1 = 3.D0 * cMolecularVolume(iComponent) * cSurfaceArea(jComponent) * ( sumRadius - sumVolumeSurface ) + sumArea * &
        &     ( 3.D0 * cMolecularVolume(iComponent) * cVolumeSurface(jComponent) + cCurvatureRadius(iComponent) * &
        &     cCurvatureRadius(iComponent) * cSurfaceArea(jComponent) ) / ( sumVolume * geoAux )
        dAux2 = 2.D0 * cMolecularVolume(iComponent) * cSurfaceArea(jComponent) * sumSqRadius * sumArea / ( 3.D0 * sumVolume * &
        &     sumVolume * geoAux )
        dEffRadialDistributionFunct_dDensity(iComponent,jComponent) = dEffPackingFraction_dDensity(iComponent,jComponent) * &
        &     ( 1.D0 + dAux1 + 2.D0 * EffPackingFraction(iComponent,jComponent) * ( dAux2 - 1.D0 ) + &
        &     EffPackingFraction(iComponent,jComponent) * EffPackingFraction(iComponent,jComponent) * ( 1.D0 - dAux1 + dAux2 ) ) &
        &     / ( ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * &
        &     ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) )
      END IF
      IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
        dEffRadialDistributionFunctSL_dDensity(iComponent,jComponent) = &
        &     ( dEffPackingFractionSutherland_dDensity(iComponent,jComponent) / ( ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) &
        &     ) ) ) + ( ( 3.D0 * cAverageDiameterRelationship(iComponent,jComponent) * ( &
        &     dEffPackingFractionSutherland_dDensity(iComponent,jComponent) * ( 1.D0 + &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) ) * AuxDiameterRelationship ) / ( ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) &
        &     ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) ) ) + ( ( 2.D0 * &
        &     cAverageDiameterRelationship(iComponent,jComponent) * cAverageDiameterRelationship(iComponent,jComponent) * &
        &     AuxDiameterRelationship * AuxDiameterRelationship * ( ( 2.D0 + EffPackingFractionSutherland(iComponent,jComponent) ) &
        &     * dEffPackingFractionSutherland_dDensity(iComponent,jComponent) * &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) ) / ( ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) &
        &     ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) ) )
      ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
        dEffRadialDistributionFunctYK_dDensity(iComponent,jComponent) = &
        &     ( dEffPackingFractionYukawa_dDensity(iComponent,jComponent) / ( ( 1.D0 - &
        &     EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) ) ) + &
        &     ( ( 3.D0 * cAverageDiameterRelationship(iComponent,jComponent) * ( &
        &     dEffPackingFractionYukawa_dDensity(iComponent,jComponent) * ( 1.D0 + EffPackingFractionYukawa(iComponent,jComponent) &
        &     ) ) * AuxDiameterRelationship ) / ( ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) ) ) + &
        &     ( ( 2.D0 * cAverageDiameterRelationship(iComponent,jComponent) * cAverageDiameterRelationship(iComponent,jComponent) &
        &     * AuxDiameterRelationship * AuxDiameterRelationship * ( ( 2.D0 + EffPackingFractionYukawa(iComponent,jComponent) ) * &
        &     EffPackingFractionYukawa(iComponent,jComponent) * dEffPackingFractionYukawa_dDensity(iComponent,jComponent) ) ) / &
        &     ( ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * &
        &     ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) ) )
      ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square-well potential
        geoAux = 3.D0 * cMolecularVolume(iComponent) + 2.D0 * cSurfaceArea(iComponent) * cCurvatureRadius(jComponent) + &
        &     cCurvatureRadius(iComponent) * cSurfaceArea(jComponent)
        dAux1 = 3.D0 * cMolecularVolume(iComponent) * cSurfaceArea(jComponent) * ( sumRadius - sumVolumeSurface ) + sumArea * &
        &     ( 3.D0 * cMolecularVolume(iComponent) * cVolumeSurface(jComponent) + cCurvatureRadius(iComponent) * &
        &     cCurvatureRadius(iComponent) * cSurfaceArea(jComponent) ) / ( sumVolume * geoAux )
        dAux2 = 2.D0 * cMolecularVolume(iComponent) * cSurfaceArea(jComponent) * sumSqRadius * sumArea / ( 3.D0 * sumVolume * &
        &     sumVolume * geoAux )
        dEffRadialDistributionFunctCSW_dDensity(iComponent,jComponent) = dEffPackingFraction_dDensity(iComponent,jComponent) * &
        &     ( 1.D0 + dAux1 + 2.D0 * EffPackingFractionCSW(iComponent,jComponent) * ( dAux2 - 1.D0 ) + &
        &     EffPackingFractionCSW(iComponent,jComponent) * EffPackingFractionCSW(iComponent,jComponent) * ( 1.D0 - dAux1 + &
        &     dAux2 ) ) / ( ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * &
        &     ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) )
      END IF
    ELSE IF( EffPFractionMixingRule == 2 ) THEN ! One-fluid van der Waals mixing rule
      IF( .NOT. PotentialTypeLogical(4) ) THEN
        dEffRadialDistributionFunct_dDensity(iComponent,jComponent) = - ( 0.5D0 * &
        &     dEffPackingFraction_dDensity(iComponent,jComponent) / ( ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) ) ) + ( ( 3.D0 * ( 1.D0 - ( EffPackingFraction(iComponent,jComponent) / &
        &     2.D0 ) ) ) / ( ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) ) ) * dEffPackingFraction_dDensity(iComponent,jComponent)
      ELSE
        gAux1 = mNonSphericityMixture * ( 1.D0 + mNonSphericityMixture ) / ( 1.D0 + 3.D0 * mNonSphericityMixture )
        gAux2 = mNonSphericityMixture * mNonSphericityMixture / ( 1.D0 + 3.D0 * mNonSphericityMixture )
        dEffRadialDistributionFunct_dDensity(iComponent,jComponent) = ( dEffPackingFraction_dDensity(iComponent,jComponent) / ( &
        &     ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) ) ) + &
        &     ( 3.D0 * gAux1 * ( dEffPackingFraction_dDensity(iComponent,jComponent) * ( 1.D0 + &
        &     EffPackingFraction(iComponent,jComponent) ) ) / ( ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) ) ) + ( ( 2.D0 &
        &     * gAux2 * EffPackingFraction(iComponent,jComponent) * ( ( 2.D0 + EffPackingFraction(iComponent,jComponent) ) * &
        &     dEffPackingFraction_dDensity(iComponent,jComponent) ) ) / ( ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * &
        &     ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) * &
        &     ( 1.D0 - EffPackingFraction(iComponent,jComponent) ) ) )
      END IF
      IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
        dEffRadialDistributionFunctSL_dDensity(iComponent,jComponent) = - ( 0.5D0 * &
        &     dEffPackingFractionSutherland_dDensity(iComponent,jComponent) / ( ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) &
        &     ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) ) ) + ( ( 3.D0 * ( 1.D0 - ( &
        &     EffPackingFractionSutherland(iComponent,jComponent) / 2.D0 ) ) ) / ( ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) &
        &     ) * ( 1.D0 - EffPackingFractionSutherland(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionSutherland(iComponent,jComponent) ) ) ) * &
        &     dEffPackingFractionSutherland_dDensity(iComponent,jComponent)
      ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
        dEffRadialDistributionFunctYK_dDensity(iComponent,jComponent) = - ( 0.5D0 * &
        &     dEffPackingFractionYukawa_dDensity(iComponent,jComponent) / ( ( 1.D0 - &
        &     EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * &
        &     ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) ) ) + ( ( 3.D0 * ( 1.D0 - ( &
        &     EffPackingFractionYukawa(iComponent,jComponent) / 2.D0 ) ) ) / ( ( 1.D0 - &
        &     EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * &
        &     ( 1.D0 - EffPackingFractionYukawa(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionYukawa(iComponent,jComponent) ) ) ) * dEffPackingFractionYukawa_dDensity(iComponent,jComponent)
      ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square-well potential
        gAux1 = mNonSphericityMixture * ( 1.D0 + mNonSphericityMixture ) / ( 1.D0 + 3.D0 * mNonSphericityMixture )
        gAux2 = mNonSphericityMixture * mNonSphericityMixture / ( 1.D0 + 3.D0 * mNonSphericityMixture )
        dEffRadialDistributionFunctCSW_dDensity(iComponent,jComponent) = ( dEffPackingFraction_dDensity(iComponent,jComponent) / &
        &     ( ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) &
        &     ) ) ) + ( 3.D0 * gAux1 * ( dEffPackingFraction_dDensity(iComponent,jComponent) * ( 1.D0 + &
        &     EffPackingFractionCSW(iComponent,jComponent) ) ) / ( ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( &
        &     1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) ) &
        &     ) + ( ( 2.D0 * gAux2 * EffPackingFractionCSW(iComponent,jComponent) * ( ( 2.D0 + &
        &     EffPackingFractionCSW(iComponent,jComponent) ) * dEffPackingFraction_dDensity(iComponent,jComponent) ) ) / ( ( &
        &     1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * &
        &     ( 1.D0 - EffPackingFractionCSW(iComponent,jComponent) ) * ( 1.D0 - &
        &     EffPackingFractionCSW(iComponent,jComponent) ) ) )
      END IF
    END IF
    ! First derivative of the mean-attractive energy with respect to the density [K . Å³]
    IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
      dMeanAttractiveEnergy_dDensity(iComponent,jComponent) = cMeanAttractiveEnergy(iComponent,jComponent) / rNumberDensity &
      &     + cMeanAttractiveEnergy(iComponent,jComponent) / cEffectiveRadialDistributionFunction(iComponent,jComponent) * &
      &     dEffRadialDistributionFunct_dDensity(iComponent,jComponent)
    ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
      dMeanAttractiveEnergy_dDensity(iComponent,jComponent) = cMeanAttractiveEnergySutherland(iComponent,jComponent) / &
      &     rNumberDensity + cMeanAttractiveEnergySutherland(iComponent,jComponent) / &
      &     cEffectiveRadialDistributionFunctionSL(iComponent,jComponent) * &
      &     dEffRadialDistributionFunctSL_dDensity(iComponent,jComponent)
    ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
      dMeanAttractiveEnergy_dDensity(iComponent,jComponent) = cMeanAttractiveEnergyYukawa(iComponent,jComponent) / &
      &     rNumberDensity + cMeanAttractiveEnergyYukawa(iComponent,jComponent) / &
      &     cEffectiveRadialDistributionFunctionYK(iComponent,jComponent) * &
      &     dEffRadialDistributionFunctYK_dDensity(iComponent,jComponent)
    ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square-well potential
      dMeanAttractiveEnergy_dDensity(iComponent,jComponent) = cMeanAttractiveEnergyCSW(iComponent,jComponent) / rNumberDensity + &
      &     cMeanAttractiveEnergyCSW(iComponent,jComponent) / cEffectiveRadialDistributionFunctionCSW(iComponent,jComponent) * &
      &     dEffRadialDistributionFunctCSW_dDensity(iComponent,jComponent)
    END IF
    ! Mean-attractive energy fluctuation between components [K²]
    IF( .NOT. PotentialTypeLogical(4) ) THEN
      cMeanAttractiveEnergyFluctuations(iComponent,jComponent) = 0.5D0 * dMeanAttractiveEnergy_dDensity(iComponent,jComponent) * &
      &     ijaWellDepth(iComponent,jComponent) * rNumberDensity * HSIsothermalCompressibility
    ELSE
      cMeanAttractiveEnergyFluctuations(iComponent,jComponent) = 0.5D0 * dMeanAttractiveEnergy_dDensity(iComponent,jComponent) * &
      &     ijaWellDepth(iComponent,jComponent) * rNumberDensity * HCBIsothermalCompressibility
    END IF
  END DO
END DO

! Zhang's correction [K²]
IF( ZhangCorrectionLogical ) THEN
  cMeanAttractiveEnergyFluctuations(:,:) = cMeanAttractiveEnergyFluctuations(:,:) * ZhangFactor
END IF

! Mixing rule for the second-order perturbation coefficient [K²]
CALL Mixing_Rules( mFraction, cMeanAttractiveEnergyFluctuations, pSecondOrderCoefficient )

! Second-order perturbation contribution to the Helmholtz free energy (molar basis) [J / mol]
mMeanAttFluctuationHelmholtzFreeEnergy = pSecondOrderCoefficient * cUniversalGas / Temperature ! Proven units

! ************************************************************************************************ !
! HIGHER-ORDER PERTURBATION CONTRIBUTION                                                           !
! ************************************************************************************************ !

! Higher order perturbation contribution to the Helmholtz free energy (molar basis) [J / mol]
IF( HigherOrderTPTLogical ) THEN
  ! Barker-Henderson's higher-order perturbation theory (approximation)
  HigherOrderTerms = 0.D0
  DO nOrder = 3, nHigherOrder
    HigherOrderTerms = HigherOrderTerms + ( mMeanAttractiveHelmholtzFreeEnergy * ( ( (2.D0 * &
    &     mMeanAttFluctuationHelmholtzFreeEnergy) / mMeanAttractiveHelmholtzFreeEnergy ) ** ( nOrder - 1 ) ) / &
    &     Factorial( nOrder ) )
  END DO
  ! Higher order perturbation contribution to the Helmholtz free energy (molar basis) [J / mol]
  mHigherOrderHelmholtzFreeEnergy = HigherOrderTerms ! Proven units
ELSE
  ! Higher order perturbation contribution to the Helmholtz free energy (molar basis) [J / mol]
  mHigherOrderHelmholtzFreeEnergy = 0.D0 ! Proven units
END IF

! ************************************************************************************************ !
! MONOMER-MONOMER TOTAL CONTRIBUTION                                                               !
! ************************************************************************************************ !

! Monomer-monomer contribution to the Helmholtz free energy (molar basis) [J / mol]
mMonomerMonomerHelmholtzFreeEnergy = mExcludedVolumeHelmholtzFreeEnergy + mMeanAttractiveHelmholtzFreeEnergy + &
&     mMeanAttFluctuationHelmholtzFreeEnergy + mHigherOrderHelmholtzFreeEnergy ! Proven units

! ************************************************************************************************ !
! TOTAL HELMHOLTZ FREE ENERGY                                                                      !
! ************************************************************************************************ !

! Total Helmholtz free energy (molar basis) [J / mol]
mTotalHelmholtzFreeEnergy = mIdealHelmholtzFreeEnergy + mMonomerMonomerHelmholtzFreeEnergy ! Proven units

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !
!          FIRST DERIVATIVES OF THE HELMHOLTZ FREE ENERGY WITH RESPECT TO THE TEMPERATURE          !
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !

! ************************************************************************************************ !
! IDEAL CONTRIBUTION                                                                               !
! ************************************************************************************************ !

! First derivative of the ideal Helmholtz free energy with respect to the temperature [J / (mol . K)]
dIdealHelmholtzFreeEnergy_dTemperature = 0.D0
DO iComponent = 1, nComponents
  ! Ideal contribution to the Helmholtz free energy per component [unitless]
  dIdealHelmholtzFreeEnergy_dTemperature = dIdealHelmholtzFreeEnergy_dTemperature + mFraction(iComponent) * DLOG( cDeBroglieIdeal &
  &     * mFraction(iComponent) * DeBroglie_Wavelength_Translational_Cb(iComponent) * DeBroglie_Wavelength_Rotational(iComponent) &
  &     / mVolume )
END DO
dIdealHelmholtzFreeEnergy_dTemperature = dIdealHelmholtzFreeEnergy_dTemperature - 3.5D0
dIdealHelmholtzFreeEnergy_dTemperature = cUniversalGas * dIdealHelmholtzFreeEnergy_dTemperature ! Proven units

! ************************************************************************************************ !
! EXCLUDED-VOLUME CONTRIBUTION                                                                     !
! ************************************************************************************************ !

! First derivative of the excluded-volume contribution to the Helmholtz free energy with respect to the temperature [J / (mol . K)]
dExcludedVolumeFreeEnergy_dTemperature = HSBoublikMansoori * cUniversalGas * AngleAverageMixture ! Proven units
IF( ReferenceBoublikLogical ) THEN
  dExcludedVolumeFreeEnergy_dTemperature = HCBBoublik * cUniversalGas
END IF

! ************************************************************************************************ !
! SECOND-ORDER PERTURBATION CONTRIBUTION                                                           !
! ************************************************************************************************ !

! First derivative of the second-order perturbation contribution to the Helmholtz free energy with respect to the temperature [J / (mol . K)]
dMeanAttFluctuationFEnergy_dTemperature = - pSecondOrderCoefficient * cUniversalGas / ( Temperature * Temperature ) ! Proven units

! ************************************************************************************************ !
! HIGHER-ORDER PERTURBATION CONTRIBUTION                                                           !
! ************************************************************************************************ !

! First derivative of the higher order perturbation contribution to the Helmholtz free energy with respect to the temperature [J / (mol . K)]
IF( HigherOrderTPTLogical ) THEN
  ! Barker-Henderson's higher-order perturbation theory (approximation)
  dHigherOrderTerms_dTemperature = 0.D0
  DO nOrder = 3, nHigherOrder
    dHigherOrderTerms_dTemperature = dHigherOrderTerms_dTemperature + ( 2.D0 * ( DBLE( nOrder - 1 ) / Factorial( nOrder ) ) ) * &
    &     ( ( (2.D0 * mMeanAttFluctuationHelmholtzFreeEnergy) / mMeanAttractiveHelmholtzFreeEnergy ) ** ( nOrder - 2 ) ) * &
    &     ( dMeanAttFluctuationFEnergy_dTemperature )
  END DO
  ! First derivative of the higher order perturbation contribution to the Helmholtz free energy with respect to the temperature [J / (mol . K)]
  dHigherOrderFEnergy_dTemperature = dHigherOrderTerms_dTemperature ! Proven units
ELSE
  ! First derivative of the higher order perturbation contribution to the Helmholtz free energy with respect to the temperature [J / (mol . K)]
  dHigherOrderFEnergy_dTemperature = 0.D0 ! Proven units
END IF

! ************************************************************************************************ !
! MONOMER-MONOMER TOTAL CONTRIBUTION                                                               !
! ************************************************************************************************ !

! First derivative of the monomer-monomer contribution to the Helmholtz free energy with respect to the temperature [J / (mol . K)]
dMonomerMonomerFEnergy_dTemperature = dExcludedVolumeFreeEnergy_dTemperature + dMeanAttFluctuationFEnergy_dTemperature + &
&     dHigherOrderFEnergy_dTemperature ! Proven units

! ************************************************************************************************ !
! TOTAL HELMHOLTZ FREE ENERGY                                                                      !
! ************************************************************************************************ !

! First derivative of the total Helmholtz free energy with respect to the temperature [J / (mol . K)]
dTotalHelmholtzFEnergy_dTemperature = dIdealHelmholtzFreeEnergy_dTemperature + dMonomerMonomerFEnergy_dTemperature ! Proven units

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !
!         SECOND DERIVATIVES OF THE HELMHOLTZ FREE ENERGY WITH RESPECT TO THE TEMPERATURE          !
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !

! ************************************************************************************************ !
! IDEAL CONTRIBUTION                                                                               !
! ************************************************************************************************ !

! Second derivative of the ideal Helmholtz free energy with respect to the temperature [J / (mol . K²)]
d2IdealHelmholtzFEnergy_d2Temperature = - 2.5D0 * cUniversalGas / Temperature ! Proven units

! ************************************************************************************************ !
! SECOND-ORDER PERTURBATION CONTRIBUTION                                                           !
! ************************************************************************************************ !

! Second derivative of the second-order perturbation contribution to the Helmholtz free energy with respect to the temperature [J / (mol . K²)]
d2MeanAttFluctFEnergy_d2Temperature = 2.D0 * pSecondOrderCoefficient * cUniversalGas / ( Temperature * Temperature * Temperature ) ! Proven units

! ************************************************************************************************ !
! HIGHER-ORDER PERTURBATION CONTRIBUTION                                                           !
! ************************************************************************************************ !

! Second derivative of the higher order perturbation contribution to the Helmholtz free energy with respect to the temperature [J / (mol . K²)]
IF( HigherOrderTPTLogical ) THEN
  ! Barker-Henderson's higher-order perturbation theory (approximation)
  d2HigherOrderTerms_d2Temperature = 0.D0
  DO nOrder = 3, nHigherOrder
    d2HigherOrderTerms_d2Temperature = d2HigherOrderTerms_d2Temperature + ( ( 4.D0 * ( DBLE( nOrder - 1 ) / &
    &     Factorial( nOrder ) ) ) * ( ( (2.D0 * mMeanAttFluctuationHelmholtzFreeEnergy) / mMeanAttractiveHelmholtzFreeEnergy ) ** &
    &     ( nOrder - 3 ) ) / mMeanAttractiveHelmholtzFreeEnergy ) * ( ( d2MeanAttFluctFEnergy_d2Temperature * &
    &     mMeanAttFluctuationHelmholtzFreeEnergy ) + ( DBLE( nOrder - 2 ) * dMeanAttFluctuationFEnergy_dTemperature * &
    &     dMeanAttFluctuationFEnergy_dTemperature ) )
  END DO
  ! Second derivative of the higher order perturbation contribution to the Helmholtz free energy with respect to the temperature [J / (mol . K²)]
  d2HigherOrderFEnergy_d2Temperature = d2HigherOrderTerms_d2Temperature ! Proven units
ELSE
  ! Second derivative of the higher order perturbation contribution to the Helmholtz free energy with respect to the temperature [J / (mol . K²)]
  d2HigherOrderFEnergy_d2Temperature = 0.D0 ! Proven units
END IF

! ************************************************************************************************ !
! MONOMER-MONOMER TOTAL CONTRIBUTION                                                               !
! ************************************************************************************************ !

! Second derivative of the monomer-monomer contribution to the Helmholtz free energy with respect to the temperature [J / (mol . K²)]
d2MonomerMonomerFEnergy_d2Temperature = d2MeanAttFluctFEnergy_d2Temperature + d2HigherOrderFEnergy_d2Temperature ! Proven units

! ************************************************************************************************ !
! TOTAL HELMHOLTZ FREE ENERGY                                                                      !
! ************************************************************************************************ !

! Second derivative of the total Helmholtz free energy with respect to the temperature [J / (mol . K²)]
d2TotalHelmholtzFEnergy_d2Temperature = d2IdealHelmholtzFEnergy_d2Temperature + d2MonomerMonomerFEnergy_d2Temperature ! Proven units

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !
!                                     THERMODYNAMIC PROPERTIES                                     !
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !

! Specific heat at constant volume (ideal contribution) [J / (mol . K)]
DO iComponent = 1, nComponents
  IF( SpecificHeatReference == "TRC" ) THEN
    vIdealSpecificHeat(iComponent) = mFraction(iComponent) * cUniversalGas * ( TemperatureParameter(iComponent,1) + &
    &     TemperatureParameter(iComponent,2) * Temperature + TemperatureParameter(iComponent,3) * Temperature * Temperature + &
    &     TemperatureParameter(iComponent,4) * Temperature * Temperature * Temperature + TemperatureParameter(iComponent,5) * &
    &     Temperature * Temperature * Temperature * Temperature - 1.D0 )
  ELSE IF( SpecificHeatReference == "NIST" ) THEN
    vIdealSpecificHeat(iComponent) = mFraction(iComponent) * ( TemperatureParameter(iComponent,1) + &
    &     TemperatureParameter(iComponent,2) * ( Temperature / 1.D3 ) + TemperatureParameter(iComponent,3) * ( Temperature * &
    &     Temperature / 1.D6 ) + TemperatureParameter(iComponent,4) * ( Temperature * Temperature * Temperature / 1.D9 ) + &
    &     TemperatureParameter(iComponent,5) / ( Temperature * Temperature / 1.D6 ) - cUniversalGas )
  END IF
END DO

! Specific heat at constant volume (ideal contribution + residual contribution) [J / (mol . K)]
vSpecificHeat = SUM( vIdealSpecificHeat ) - (Temperature * d2MonomerMonomerFEnergy_d2Temperature)

! Entropy [J / (mol . K)]
Entropy = - dTotalHelmholtzFEnergy_dTemperature

RETURN

END SUBROUTINE Calculate_Entropy

! ************************************************************************************************ !
! Calculates the entropy and the derivatives of the Helmholtz free energy with respect to the      !
! temperature for pure components.                                                                 !
! ************************************************************************************************ !
SUBROUTINE Calculate_Entropy_Single_Component( cComponent, mVolume, Temperature, Entropy, vSpecificHeat )

! Uses two modules: global variables and substances
USE GlobalVariables
USE Substances

IMPLICIT NONE

! ************************************************************************************************ !
! INTEGER VARIABLES                                                                                !
! ************************************************************************************************ !
INTEGER( Kind= Int64 ) :: cComponent ! Component indexes
INTEGER( Kind= Int64 ) :: nOrder     ! Counter of TPT coefficients

! ************************************************************************************************ !
! REAL PARAMETERS                                                                                  !
! ************************************************************************************************ !
REAL( Kind= Real64 ), PARAMETER                    :: CondensationPFraction = 0.493D0                                          ! Condensation packing fraction
REAL( Kind= Real64 ), PARAMETER                    :: ZhangCorrection = 1.D0 / (CondensationPFraction * CondensationPFraction) ! Zhang's correction
REAL( Kind= Real64 ), DIMENSION( 3, 3 ), PARAMETER :: SWCoefficientsMatrix = RESHAPE( [ [2.258550D0, -0.669270D0, 1.015760D1], &
&                                                                                       [-1.50349D0, 1.4004900D0, -1.50427D1], &
&                                                                                       [0.249434D0, -0.827739D0, 5.308270D0]  &
&                                                                                     ], [3, 3] )                              ! Matrix of coefficients of the effective packing fraction
REAL( Kind= Real64 ), DIMENSION( 2, 4 ), PARAMETER :: SLCoefficientsMatrix = RESHAPE( [ [-0.943973D0, 0.370942D0], &
&                                                                                       [0.422543D0, -0.173333D0], &
&                                                                                       [-3.71763D-2, 1.75599D-2], &
&                                                                                       [1.16901D-3, -5.72729D-4]  &
&                                                                                     ], [2, 4] )                              ! Matrix of coefficients of the effective packing fraction
REAL( Kind= Real64 ), DIMENSION( 2, 3 ), PARAMETER :: YK1CoefficientsMatrix = RESHAPE( [ [0.900678D0, -0.314300D0], &
&                                                                                        [-1.50051D0, 0.2571010D0], &
&                                                                                        [0.776577D0, -4.31566D-2]  &
&                                                                                      ], [2, 3] )                             ! Matrix of coefficients of the effective packing fraction
REAL( Kind= Real64 ), DIMENSION( 2, 5 ), PARAMETER :: YK2CoefficientsMatrix = RESHAPE( [ [0.989601D0, -1.19152D-2], &
&                                                                                        [-0.872203D0, -1.24029D0], &
&                                                                                        [0.320808D0, 2.41636D0], &
&                                                                                        [0.00000D0, -2.01922D0], &
&                                                                                        [0.00000D0, 0.647565D0]  &
&                                                                                      ], [2, 5] )                             ! Matrix of coefficients of the effective packing fraction

! ************************************************************************************************ !
! REAL VARIABLES                                                                                   !
! ************************************************************************************************ !
REAL( Kind= Real64 )                    :: mVolume                                  ! Molar volume
REAL( Kind= Real64 )                    :: cMolecularVolume                         ! Molecular volume (component) [Å³]
REAL( Kind= Real64 )                    :: cSurfaceArea                             ! Surface area (component) [Å²]
REAL( Kind= Real64 )                    :: cCurvatureRadius                         ! Mean radius of curvature (component) [Å]
REAL( Kind= Real64 )                    :: cRadius                                  ! Radius (component) [Å]
REAL( Kind= Real64 )                    :: cNonSphericity                           ! Nonsphericity parameter (component)
REAL( Kind= Real64 )                    :: cSecondVirialCoefficient                 ! Second virial coefficient of the component
REAL( Kind= Real64 )                    :: fMolecularVolume                         ! Molecular volume (potential) [Å³]
REAL( Kind= Real64 )                    :: fSurfaceArea                             ! Surface area (potential) [Å²]
REAL( Kind= Real64 )                    :: fCurvatureRadius                         ! Mean radius of curvature (potential) [Å]
REAL( Kind= Real64 )                    :: fRadius                                  ! Radius (potential) [Å]
REAL( Kind= Real64 )                    :: fNonSphericity                           ! Nonsphericity parameter (potential)
REAL( Kind= Real64 )                    :: fSecondVirialCoefficient                 ! Second virial coefficient of the field
REAL( Kind= Real64 )                    :: MomentInertia                            ! Moment of inertia (non-spherical rigid bodies)
REAL( Kind= Real64 )                    :: Temperature                              ! Temperature
REAL( Kind= Real64 )                    :: rNumberDensity                           ! Reduced number density
REAL( Kind= Real64 )                    :: HigherOrderTerms                         ! Higher-order perturbation terms
REAL( Kind= Real64 )                    :: Factorial                                ! Factorial function
REAL( Kind= Real64 )                    :: mHigherOrderHelmholtzFreeEnergy          ! Higher-order contribution to the Helmholtz free energy (molar basis)
REAL( Kind= Real64 )                    :: mIdealHelmholtzFreeEnergy                ! Ideal contribution to the Helmholtz free energy (molar basis)
REAL( Kind= Real64 )                    :: mExcludedVolumeHelmholtzFreeEnergy       ! Excluded-volume contribution to the Helmholtz free energy (molar basis)
REAL( Kind= Real64 )                    :: mMeanAttractiveHelmholtzFreeEnergy       ! Mean-attractive energy contribution to the Helmholtz free energy (molar basis)
REAL( Kind= Real64 )                    :: mMeanAttFluctuationHelmholtzFreeEnergy   ! Mean-attractive energy fluctuation contribution to the Helmholtz free energy (molar basis)
REAL( Kind= Real64 )                    :: mMonomerMonomerHelmholtzFreeEnergy       ! Monomer-monomer contribution to the Helmholtz free energy (molar basis)
REAL( Kind= Real64 )                    :: mTotalHelmholtzFreeEnergy                ! Total Helmholtz free energy (molar basis)
REAL( Kind= Real64 )                    :: HSBoublikMansoori                        ! Boublik-Mansoori hard-sphere mixture term
REAL( Kind= Real64 )                    :: AngleAverage                             ! Angle average of the excluded volume of a pair of non-spherical rigid bodies
REAL( Kind= Real64 )                    :: ZhangFactor                              ! Zhang's correction factor
REAL( Kind= Real64 )                    :: AuxDiameterRelationship                  ! Auxiliary factor (diameter relationship)
REAL( Kind= Real64 )                    :: HSIsothermalCompressibility              ! Hard-sphere isothermal compressibility (Percus-Yevick expression)
REAL( Kind= Real64 )                    :: HCBIsothermalCompressibility             ! Hard convex-body isothermal compressibility (Boublik expression)
REAL( Kind= Real64 )                    :: Entropy                                  ! Entropy
REAL( Kind= Real64 )                    :: vSpecificHeat                            ! Specific heat capacity at constant volume
REAL( Kind= Real64 )                    :: dIdealHelmholtzFreeEnergy_dTemperature   ! First derivative of the ideal Helmholtz free energy with respect to the temperature
REAL( Kind= Real64 )                    :: dExcludedVolumeFreeEnergy_dTemperature   ! First derivative of the excluded-volume contribution to the Helmholtz free energy with respect to the temperature
REAL( Kind= Real64 )                    :: dMeanAttFluctuationFEnergy_dTemperature  ! First derivative of the mean-attractive energy fluctuation contribution to the Helmholtz free energy with respect to the temperature
REAL( Kind= Real64 )                    :: dHigherOrderTerms_dTemperature           ! First derivative of the higher-order perturbation terms with respect to the temperature
REAL( Kind= Real64 )                    :: dHigherOrderFEnergy_dTemperature         ! First derivative of the higher-order perturbation contribution to the Helmholtz free energy with respect to the temperature
REAL( Kind= Real64 )                    :: dMonomerMonomerFEnergy_dTemperature      ! First derivative of the monomer-monomer contribution to the Helmholtz free energy with respect to the temperature
REAL( Kind= Real64 )                    :: dTotalHelmholtzFEnergy_dTemperature      ! First derivative of the total Helmholtz free energy with respect to the temperature
REAL( Kind= Real64 )                    :: d2IdealHelmholtzFEnergy_d2Temperature    ! Second derivative of the ideal Helmholtz free energy with respect to the temperature
REAL( Kind= Real64 )                    :: d2MeanAttFluctFEnergy_d2Temperature      ! Second derivative of the mean-attractive energy fluctuation contribution to the Helmholtz free energy with respect to the temperature
REAL( Kind= Real64 )                    :: d2HigherOrderTerms_d2Temperature         ! Second derivative of the higher-order perturbation terms with respect to the temperature
REAL( Kind= Real64 )                    :: d2HigherOrderFEnergy_d2Temperature       ! Second derivative of the higher-order perturbation contribution to the Helmholtz free energy with respect to the temperature
REAL( Kind= Real64 )                    :: d2MonomerMonomerFEnergy_d2Temperature    ! Second derivative of the monomer-monomer contribution to the Helmholtz free energy with respect to the temperature
REAL( Kind= Real64 )                    :: d2TotalHelmholtzFEnergy_d2Temperature    ! Second derivative of the total Helmholtz free energy with respect to the temperature
REAL( Kind= Real64 )                    :: vIdealSpecificHeat                       ! Specific heat capacity at constant volume (ideal contribution)
REAL( Kind= Real64 )                    :: DeBroglie_Wavelength_Translational       ! Thermal de Broglie wavelength
REAL( Kind= Real64 )                    :: DeBroglie_Wavelength_Translational_Cb    ! Thermal de Broglie wavelength (cubic)
REAL( Kind= Real64 )                    :: DeBroglie_Wavelength_Rotational          ! Thermal de Broglie wavelength (rotational contribution)
REAL( Kind= Real64 )                    :: dEffRadialDistributionFunct_dDensity     ! First derivative of the radial distribution function with respect to the density
REAL( Kind= Real64 )                    :: dEffRadialDistributionFunctYK_dDensity   ! First derivative of the radial distribution function with respect to the density (Yukawa potential)
REAL( Kind= Real64 )                    :: dEffRadialDistributionFunctSL_dDensity   ! First derivative of the radial distribution function with respect to the density (Sutherland potential)
REAL( Kind= Real64 )                    :: dEffRadialDistributionFunctCSW_dDensity  ! First derivative of the radial distribution function with respect to the density (Convex square-well potential)
REAL( Kind= Real64 )                    :: dMeanAttractiveEnergy_dDensity           ! First derivative of the mean-attractive energy with respect to the density
REAL( Kind= Real64 )                    :: cAverageDiameterRelationship             ! Diameter relationship between components (the product of diameters divided by the sum of diameters)
REAL( Kind= Real64 )                    :: cRadialDistributionFunction              ! Contact radial distribution function
REAL( Kind= Real64 )                    :: cEffectiveRadialDistributionFunction     ! Contact radial distribution function for an effective packing fraction
REAL( Kind= Real64 )                    :: cEffectiveRadialDistributionFunctionYK   ! Contact radial distribution function for an effective packing fraction (Yukawa potential)
REAL( Kind= Real64 )                    :: cEffectiveRadialDistributionFunctionSL   ! Contact radial distribution function for an effective packing fraction (Sutherland potential)
REAL( Kind= Real64 )                    :: cEffectiveRadialDistributionFunctionCSW  ! Contact radial distribution function for an effective packing fraction (Convex square-well potential)
REAL( Kind= Real64 )                    :: cMeanAttractiveEnergy                    ! Mean-attractive energy between components
REAL( Kind= Real64 )                    :: cMeanAttractiveEnergySutherland          ! Mean-attractive energy between components (Sutherland potential)
REAL( Kind= Real64 )                    :: cMeanAttractiveEnergyYukawa              ! Mean-attractive energy between components (Yukawa potential)
REAL( Kind= Real64 )                    :: cMeanAttractiveEnergyCSW                 ! Mean-attractive energy between components (Convex square-well potential)
REAL( Kind= Real64 )                    :: cMeanAttractiveEnergyFluctuations        ! Mean-attractive energy fluctuations between components
REAL( Kind= Real64 )                    :: gAux1, gAux2                             ! Auxiliars (nonsphericity expressions)
REAL( Kind= Real64 )                    :: EffPackingFraction                       ! Effective packing fraction
REAL( Kind= Real64 )                    :: EffPackingFractionYukawa                 ! Effective packing fraction (Yukawa potential)
REAL( Kind= Real64 )                    :: EffPackingFractionSutherland             ! Effective packing fraction (Sutherland potential)
REAL( Kind= Real64 )                    :: EffPackingFractionCSW                    ! Effective packing fraction (Convex square-well potential)
REAL( Kind= Real64 )                    :: dEffPackingFraction_dDensity             ! First derivative of the effective packing fraction with respect to the density
REAL( Kind= Real64 )                    :: dEffPackingFractionYukawa_dDensity       ! First derivative of the effective packing fraction with respect to the density (Yukawa potential)
REAL( Kind= Real64 )                    :: dEffPackingFractionSutherland_dDensity   ! First derivative of the effective packing fraction with respect to the density (Sutherland potential)
REAL( Kind= Real64 )                    :: dEffPackingFractionCSW_dDensity          ! First derivative of the effective packing fraction with respect to the density (Convex square-well potential)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: EffPackingFractionCoefficients           ! Coefficients of the effective packing fraction
REAL( Kind= Real64 ), DIMENSION( 3 )    :: EffPackingFractionCoefficientsYukawa     ! Coefficients of the effective packing fraction (Yukawa potential)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: EffPackingFractionCoefficientsSutherland ! Coefficients of the effective packing fraction (Sutherland potential)
REAL( Kind= Real64 ), DIMENSION( 3 )    :: EffPackingFractionCoefficientsCSW        ! Coefficients of the effective packing fraction (Convex square-well potential)
REAL( Kind= Real64 ), DIMENSION( 0:3 )  :: rDensity                                 ! Reduced densities
REAL( Kind= Real64 ), DIMENSION( 0:3 )  :: rDensityConstants                        ! Reduced density constants
REAL( Kind= Real64 ), DIMENSION( 2, 9 ) :: CSWAlphaCoefficients                     ! Coefficients of the nonsphericity (Convex square-well potential)
REAL( Kind= Real64 ), DIMENSION( 3, 9 ) :: CSW1CoefficientsMatrix                   ! Matrix of coefficients of the effective packing fraction (first-order coefficient)
REAL( Kind= Real64 ), DIMENSION( 3, 9 ) :: CSW2CoefficientsMatrix                   ! Matrix of coefficients of the effective packing fraction (second-order coefficient)

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !
!                                          INITIALIZATION                                          !
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !

! Reduced number density [1 / m³]
rNumberDensity = cAvogadro / mVolume

! Reduced number density [1 / Å³]
rNumberDensity = rNumberDensity * 1.D-30

! Radius
cRadius = 0.5D0 * aDiameter(cComponent)

! Morphological descriptors
IF( GeometrySpecification(cComponent,1) ) THEN ! Ellipsoids-of-revolution
  ! Molecular volume [Å³]
  cMolecularVolume = ( cPi / 6.D0 ) * aDiameter(cComponent) * aDiameter(cComponent) * aDiameter(cComponent) * &
  &     cAspectRatio(cComponent)
  ! Prolate ellipsoids-of-revolution
  IF( cAspectRatio(cComponent) > 1.D0 ) THEN
    ! Radius of curvature [Å]
    cCurvatureRadius = 0.5D0 * cAspectRatio(cComponent) * cRadius + ( 0.5D0 * cRadius / DSQRT( cAspectRatio(cComponent) * &
    &     cAspectRatio(cComponent) - 1.D0 ) ) * DATANH( DSQRT( cAspectRatio(cComponent) * cAspectRatio(cComponent) - 1.D0 ) / &
    &     cAspectRatio(cComponent) )
    ! Surface area [Å²]
    cSurfaceArea = ( 2.D0 * cPi * cRadius * cRadius / DSQRT( cAspectRatio(cComponent) * cAspectRatio(cComponent) - 1.D0 ) ) * ( ( &
    &     cAspectRatio(cComponent) * cAspectRatio(cComponent) * DASIN( DSQRT( cAspectRatio(cComponent) * &
    &     cAspectRatio(cComponent) - 1.D0 ) / cAspectRatio(cComponent) ) ) + ( DSQRT( cAspectRatio(cComponent) * &
    &     cAspectRatio(cComponent) - 1.D0 ) ) )
  ! Oblate ellipsoids-of-revolution
  ELSE IF( cAspectRatio(cComponent) < 1.D0 ) THEN
    ! Radius of curvature [Å]
    cCurvatureRadius = 0.5D0 * cAspectRatio(cComponent) * cRadius + ( 0.5D0 * cRadius / DSQRT( 1.D0 - cAspectRatio(cComponent) * &
    &     cAspectRatio(cComponent) ) ) * DATAN( DSQRT( 1.D0 - cAspectRatio(cComponent) * cAspectRatio(cComponent) ) / &
    &     cAspectRatio(cComponent) )
    ! Surface area [Å²]
    cSurfaceArea = ( 2.D0 * cPi * cRadius * cRadius / DSQRT( 1.D0 - cAspectRatio(cComponent) * cAspectRatio(cComponent) ) ) * ( ( &
    &     cAspectRatio(cComponent) * cAspectRatio(cComponent) * DLOG( ( 1.D0 + DSQRT( 1.D0 - cAspectRatio(cComponent) * &
    &     cAspectRatio(cComponent) ) ) / cAspectRatio(cComponent) ) ) + ( DSQRT( 1.D0 - cAspectRatio(cComponent) * &
    &     cAspectRatio(cComponent) ) ) )
  ! Spheres
  ELSE
    ! Radius of curvature [Å]
    cCurvatureRadius = cRadius
    ! Surface area [Å²]
    cSurfaceArea = 4.D0 * cPi * cRadius * cRadius
  END IF
ELSE IF( GeometrySpecification(cComponent,2) ) THEN ! Spherocylinders
  ! Molecular volume [Å³]
  cMolecularVolume = ( cPi / 4.D0 ) * aDiameter(cComponent) * aDiameter(cComponent) * aDiameter(cComponent) * &
  &     ( cAspectRatio(cComponent) + (2.D0 / 3.D0) )
  ! Prolate spherocylinders
  IF( cAspectRatio(cComponent) > 0.D0 ) THEN
    ! Radius of curvature [Å]
    cCurvatureRadius = 0.5D0 * ( cAspectRatio(cComponent) + 2.D0 ) * cRadius
    ! Surface area [Å²]
    cSurfaceArea = 4.D0 * cPi * cRadius * cRadius * ( cAspectRatio(cComponent) + 1.D0 )
  ! Spheres
  ELSE
    ! Radius of curvature [Å]
    cCurvatureRadius = cRadius
    ! Surface area [Å²]
    cSurfaceArea = 4.D0 * cPi * cRadius * cRadius
  END IF
ELSE IF( GeometrySpecification(cComponent,3) ) THEN ! Cylinders
  ! Molecular volume [Å³]
  cMolecularVolume = ( cPi / 4.D0 ) * aDiameter(cComponent) * aDiameter(cComponent) * aDiameter(cComponent) * &
  &     cAspectRatio(cComponent)
  ! Radius of curvature [Å]
  cCurvatureRadius = 0.5D0 * ( cAspectRatio(cComponent) + (0.5D0 * cPi) ) * cRadius
  ! Surface area [Å²]
  cSurfaceArea = 2.D0 * cPi * cRadius * cRadius * ( 1.D0 + 2.D0 * cAspectRatio(cComponent) )
END IF

! Non-sphericity parameter (Boublík, J. Chem. Phys., 63, 1975)
cNonSphericity = cCurvatureRadius * cSurfaceArea / 3.D0 / cMolecularVolume

! Second virial coefficient of the component (Isihara-Hadwiger theorem)
cSecondVirialCoefficient = cMolecularVolume + cCurvatureRadius * cSurfaceArea

! Zhang's correction factor [unitless]
IF( ZhangCorrectionLogical ) THEN
  IF( .NOT. PotentialTypeLogical(4) ) THEN
    ZhangFactor = 1.D0 + ( 2.D0 * ZhangCorrection * ( cPi * rNumberDensity * aDiameterSphere(cComponent) * &
    &     aDiameterSphere(cComponent) * aDiameterSphere(cComponent) / 6.D0 ) * ( cPi * rNumberDensity * &
    &     aDiameterSphere(cComponent) * aDiameterSphere(cComponent) * aDiameterSphere(cComponent) / 6.D0 ) )
  ELSE
    ZhangFactor = 1.D0 + ( 0.125D0 * ZhangCorrection * ( rNumberDensity * ijSecondVirialCoefficient(cComponent,cComponent) ) * ( &
    &     rNumberDensity * ijSecondVirialCoefficient(cComponent,cComponent) ) )
  END IF
END IF

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !
!                                HELMHOLTZ FREE ENERGY CALCULATION                                 !
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !

! ************************************************************************************************ !
! IDEAL CONTRIBUTION                                                                               !
! ************************************************************************************************ !

! Ideal contribution to the Helmholtz free energy (molar basis) [J / mol] - Translational de Broglie wavelength contribution
DeBroglie_Wavelength_Translational = 1.D0 / DSQRT( 2.D0 * cPi * cMolarMass(cComponent) * cUniversalGas * Temperature )
DeBroglie_Wavelength_Translational_Cb = DeBroglie_Wavelength_Translational * DeBroglie_Wavelength_Translational * &
&     DeBroglie_Wavelength_Translational

! Ideal contribution to the Helmholtz free energy (molar basis) [J / mol] - Rotational de Broglie wavelength contribution
IF( GeometrySpecification(cComponent,1) ) THEN ! Ellipsoids of revolution
  MomentInertia = ( 1.D0 / 20.D0 ) * cMolarMass(cComponent) * ( aDiameter(cComponent) * aDiameter(cComponent) + &
  &     aLength(cComponent) * aLength(cComponent) ) ! [Kg . Å² / mol]
  DeBroglie_Wavelength_Rotational = (cDeBroglie * cDeBroglie) / (cPi * MomentInertia * cUniversalGas * Temperature) ! [m² / Å²]
  DeBroglie_Wavelength_Rotational = DeBroglie_Wavelength_Rotational * 1.D20 ! [unitless]
ELSE IF( GeometrySpecification(cComponent,2) ) THEN ! Spherocylinders
  MomentInertia = (3.D0 / 8.D0) * ( ( cMolarMass(cComponent) / ( 3.D0 * cAspectRatio(cComponent) + 2.D0 ) ) * &
  &     ( aDiameter(cComponent) * aDiameter(cComponent) ) ) * ( ( ( cAspectRatio(cComponent) / 6.D0 ) * ( 3.D0 + ( 4.D0 * &
  &     cAspectRatio(cComponent) * cAspectRatio(cComponent) ) ) ) + ( (4.D0 / 3.D0) * ( (83.D0 / 320.D0) + &
  &     ( cAspectRatio(cComponent) + (3.D0 / 8.D0) ) * ( cAspectRatio(cComponent) + (3.D0 / 8.D0) ) ) ) ) ! [Kg . Å² / mol]
  DeBroglie_Wavelength_Rotational = (cDeBroglie * cDeBroglie) / (cPi * MomentInertia * cUniversalGas * Temperature) ! [m² / Å²]
  DeBroglie_Wavelength_Rotational = DeBroglie_Wavelength_Rotational * 1.D20 ! [unitless]
ELSE IF( GeometrySpecification(cComponent,3) ) THEN ! Cylinders
  MomentInertia = (1.D0 / 48.D0) * cMolarMass(cComponent) * ( ( 4.D0 * aLength(cComponent) * aLength(cComponent) ) + &
  &     ( 3.D0 * aDiameter(cComponent) * aDiameter(cComponent) ) ) ! [Kg . Å² / mol]
  DeBroglie_Wavelength_Rotational = (cDeBroglie * cDeBroglie) / (cPi * MomentInertia * cUniversalGas * Temperature) ! [m² / Å²]
  DeBroglie_Wavelength_Rotational = DeBroglie_Wavelength_Rotational * 1.D20 ! [unitless]
END IF

! Ideal contribution to the Helmholtz free energy (molar basis) [J / mol]
mIdealHelmholtzFreeEnergy = DLOG( cDeBroglieIdeal * DeBroglie_Wavelength_Translational_Cb * DeBroglie_Wavelength_Rotational / &
&     mVolume )
mIdealHelmholtzFreeEnergy = mIdealHelmholtzFreeEnergy - 1.D0
mIdealHelmholtzFreeEnergy = cUniversalGas * Temperature * mIdealHelmholtzFreeEnergy ! Proven units

! ************************************************************************************************ !
! EXCLUDED-VOLUME CONTRIBUTION                                                                     !
! ************************************************************************************************ !

! Reduced density constants
rDensityConstants(0) = (cPi / 6.D0) ! [unitless]
rDensityConstants(1) = (cPi / 6.D0) * aDiameterSphere(cComponent) ! [Å]
rDensityConstants(2) = (cPi / 6.D0) * aDiameterSphere(cComponent) * aDiameterSphere(cComponent) ! [Å²]
rDensityConstants(3) = (cPi / 6.D0) * aDiameterSphere(cComponent) * aDiameterSphere(cComponent) * aDiameterSphere(cComponent) ! [Å³]

! Reduced densities
rDensity(0) = rDensityConstants(0) * rNumberDensity ! [1 / Å³]
rDensity(1) = rDensityConstants(1) * rNumberDensity ! [1 / Å²]
rDensity(2) = rDensityConstants(2) * rNumberDensity ! [1 / Å]
rDensity(3) = rDensityConstants(3) * rNumberDensity ! [unitless]

! Boublik-Mansoori hard-sphere mixture term [1 / Å³]
HSBoublikMansoori = ( ( rDensity(2) * rDensity(2) * rDensity(2) / ( rDensity(3) * rDensity(3) ) ) - rDensity(0) ) * DLOG( 1.D0 - &
&     rDensity(3) ) + ( 3.D0 * rDensity(1) * rDensity(2) / ( 1.D0 - rDensity(3) ) ) + ( rDensity(2) * rDensity(2) * rDensity(2) / &
&     ( rDensity(3) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) )

! Boublik-Mansoori hard-sphere mixture term [unitless]
HSBoublikMansoori = 6.D0 * HSBoublikMansoori / cPi / rNumberDensity

! Angle average of the excluded volume of a pair of non-spherical rigid bodies [unitless]
IF( GeometrySpecification(cComponent,1) ) THEN ! Ellipsoids-of-revolution
  ! Analytical expression of the angle average of the excluded volume of a pair of non-spherical rigid bodies [unitless]
  IF( cAspectRatio(cComponent) /= 1.D0 ) THEN
    IF( cAspectRatio(cComponent) > 1.D0 ) THEN
      AngleAverage = 2.D0 + 1.5D0 * ( 1.D0 + DATANH( DSQRT( ( cAspectRatio(cComponent) * cAspectRatio(cComponent) ) - 1.D0 ) / &
      &     cAspectRatio(cComponent) ) / ( cAspectRatio(cComponent) * DSQRT( ( cAspectRatio(cComponent) * &
      &     cAspectRatio(cComponent) ) - 1.D0 ) ) ) * ( 1.D0 + ( ( cAspectRatio(cComponent) * cAspectRatio(cComponent) ) / &
      &     DSQRT( ( cAspectRatio(cComponent) * cAspectRatio(cComponent) ) - 1.D0 ) ) * DATAN( DSQRT( ( &
      &     cAspectRatio(cComponent) * cAspectRatio(cComponent) ) - 1.D0 ) ) )
    ELSE IF( cAspectRatio(cComponent) < 1.D0 ) THEN
      AngleAverage = 2.D0 + 1.5D0 * ( 1.D0 + DATANH( DSQRT( ( ( 1.D0 / cAspectRatio(cComponent) ) * ( 1.D0 / &
      &     cAspectRatio(cComponent) ) ) - 1.D0 ) / ( 1.D0 / cAspectRatio(cComponent) ) ) / ( ( 1.D0 / &
      &     cAspectRatio(cComponent) ) * DSQRT( ( ( 1.D0 / cAspectRatio(cComponent) ) * ( 1.D0 / cAspectRatio(cComponent) ) ) &
      &     - 1.D0 ) ) ) * ( 1.D0 + ( ( ( 1.D0 / cAspectRatio(cComponent) ) * ( 1.D0 / cAspectRatio(cComponent) ) ) / &
      &     DSQRT( ( ( 1.D0 / cAspectRatio(cComponent) ) * ( 1.D0 / cAspectRatio(cComponent) ) ) - 1.D0 ) ) * &
      &     DATAN( DSQRT( ( ( 1.D0 / cAspectRatio(cComponent) ) * ( 1.D0 / cAspectRatio(cComponent) ) ) - 1.D0 ) ) )
    END IF
    AngleAverage = 0.125D0 * AngleAverage
  ELSE
    AngleAverage = 1.D0
  END IF
ELSE IF( GeometrySpecification(cComponent,2) ) THEN ! Spherocylinders
  AngleAverage = ( 12.D0 / ( ( 3.D0 * cAspectRatio(cComponent) ) + 2.D0 ) ) * ( (4.D0 / 3.D0) + ( 2.D0 * &
  &     cAspectRatio(cComponent) ) + ( ( cAspectRatio(cComponent) * cAspectRatio(cComponent) ) / 2.D0 ) )
  AngleAverage = 0.125D0 * AngleAverage
ELSE IF( GeometrySpecification(cComponent,3) ) THEN ! Cylinders
  AngleAverage = 2.D0 * cAspectRatio(cComponent) + cPi + 3.D0 + ( cPi / ( 2.D0 * cAspectRatio(cComponent) ) )
  AngleAverage = 0.125D0 * AngleAverage
END IF

! Excluded-volume contribution to the Helmholtz free energy (molar basis) [J / mol]
mExcludedVolumeHelmholtzFreeEnergy = HSBoublikMansoori * cUniversalGas * Temperature * AngleAverage ! Proven units
IF( ReferenceBoublikLogical ) THEN
  mExcludedVolumeHelmholtzFreeEnergy = ( cNonSphericity * rDensity(3) * ( cNonSphericity + 3.D0 - 3.D0 * rDensity(3) ) / &
  &     ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) ) - ( ( 1.D0 - cNonSphericity * cNonSphericity ) * &
  &     DLOG( 1.D0 - rDensity(3) ) )
  mExcludedVolumeHelmholtzFreeEnergy = mExcludedVolumeHelmholtzFreeEnergy * cUniversalGas * Temperature
END IF

! ************************************************************************************************ !
! FIRST-ORDER PERTURBATION CONTRIBUTION                                                            !
! ************************************************************************************************ !

! Radius
fRadius = 0.5D0 * aDiameterField(cComponent)

! Morphological descriptors (field)
IF( GeometrySpecification(cComponent,1) ) THEN ! Ellipsoids-of-revolution
  ! Molecular volume [Å³]
  fMolecularVolume = ( cPi / 6.D0 ) * aDiameterField(cComponent) * aDiameterField(cComponent) * aDiameterField(cComponent) * &
  &     fAspectRatio(cComponent)
  ! Prolate ellipsoids-of-revolution
  IF( fAspectRatio(cComponent) > 1.D0 ) THEN
    ! Radius of curvature [Å]
    fCurvatureRadius = 0.5D0 * fAspectRatio(cComponent) * fRadius + ( 0.5D0 * fRadius / DSQRT( fAspectRatio(cComponent) * &
    &     fAspectRatio(cComponent) - 1.D0 ) ) * DATANH( DSQRT( fAspectRatio(cComponent) * fAspectRatio(cComponent) - 1.D0 ) / &
    &     fAspectRatio(cComponent) )
    ! Surface area [Å²]
    fSurfaceArea = ( 2.D0 * cPi * fRadius * fRadius / DSQRT( fAspectRatio(cComponent) * fAspectRatio(cComponent) - 1.D0 ) ) * ( ( &
    &     fAspectRatio(cComponent) * fAspectRatio(cComponent) * DASIN( DSQRT( fAspectRatio(cComponent) * &
    &     fAspectRatio(cComponent) - 1.D0 ) / fAspectRatio(cComponent) ) ) + ( DSQRT( fAspectRatio(cComponent) * &
    &     fAspectRatio(cComponent) - 1.D0 ) ) )
  ! Oblate ellipsoids-of-revolution
  ELSE IF( fAspectRatio(cComponent) < 1.D0 ) THEN
    ! Radius of curvature [Å]
    fCurvatureRadius = 0.5D0 * fAspectRatio(cComponent) * fRadius + ( 0.5D0 * fRadius / DSQRT( 1.D0 - fAspectRatio(cComponent) * &
    &     fAspectRatio(cComponent) ) ) * DATAN( DSQRT( 1.D0 - fAspectRatio(cComponent) * fAspectRatio(cComponent) ) / &
    &     fAspectRatio(cComponent) )
    ! Surface area [Å²]
    fSurfaceArea = ( 2.D0 * cPi * fRadius * fRadius / DSQRT( 1.D0 - fAspectRatio(cComponent) * fAspectRatio(cComponent) ) ) * ( ( &
    &     fAspectRatio(cComponent) * fAspectRatio(cComponent) * DLOG( ( 1.D0 + DSQRT( 1.D0 - fAspectRatio(cComponent) * &
    &     fAspectRatio(cComponent) ) ) / fAspectRatio(cComponent) ) ) + ( DSQRT( 1.D0 - fAspectRatio(cComponent) * &
    &     fAspectRatio(cComponent) ) ) )
  ! Spheres
  ELSE
    ! Radius of curvature [Å]
    fCurvatureRadius = fRadius
    ! Surface area [Å²]
    fSurfaceArea = 4.D0 * cPi * fRadius * fRadius
  END IF
ELSE IF( GeometrySpecification(cComponent,2) .OR. GeometrySpecification(cComponent,3) ) THEN ! Spherocylinders and cylinders
  ! Molecular volume [Å³]
  fMolecularVolume = ( cPi / 4.D0 ) * aDiameterField(cComponent) * aDiameterField(cComponent) * aDiameterField(cComponent) * &
  &     ( fAspectRatio(cComponent) + (2.D0 / 3.D0) )
  ! Prolate spherocylinders
  IF( fAspectRatio(cComponent) > 0.D0 ) THEN
    ! Radius of curvature [Å]
    fCurvatureRadius = 0.5D0 * ( fAspectRatio(cComponent) + 2.D0 ) * fRadius
    ! Surface area [Å²]
    fSurfaceArea = 4.D0 * cPi * fRadius * fRadius * ( fAspectRatio(cComponent) + 1.D0 )
  ! Spheres
  ELSE
    ! Radius of curvature [Å]
    fCurvatureRadius = fRadius
    ! Surface area [Å²]
    fSurfaceArea = 4.D0 * cPi * fRadius * fRadius
  END IF
END IF

! Non-sphericity parameter (Boublík, J. Chem. Phys., 63, 1975)
fNonSphericity = fCurvatureRadius * fSurfaceArea / 3.D0 / fMolecularVolume

! Second virial coefficient of the force field (Isihara-Hadwiger theorem)
fSecondVirialCoefficient = fMolecularVolume + fCurvatureRadius * fSurfaceArea

! Auxiliary diameter relationship [1 / Å]
AuxDiameterRelationship = 1.D0 / aDiameterSphere(cComponent)

! Diameter relationship between components (the product of diameters divided by the sum of diameters) [Å]
cAverageDiameterRelationship = 0.5D0 * aDiameterSphere(cComponent)

! Effective packing fraction coefficients
IF( PotentialTypeLogical(4) ) THEN ! Convex square-well potential
  IF( GeometrySpecification(cComponent,1) ) THEN ! Ellipsoids-of-revolution
    CSW1CoefficientsMatrix = CSWEORCoefficientA1
    IF( PYHCBCorrectionLogical ) THEN
      CSW2CoefficientsMatrix = CSWEORCoefficientA2PY
    ELSE
      CSW2CoefficientsMatrix = CSWEORCoefficientA2
    END IF
  ELSE IF( GeometrySpecification(cComponent,2) ) THEN ! Spherocylinders
    CSW1CoefficientsMatrix = CSWSPCCoefficientA1
    IF( PYHCBCorrectionLogical ) THEN
      CSW2CoefficientsMatrix = CSWSPCCoefficientA2PY
    ELSE
      CSW2CoefficientsMatrix = CSWSPCCoefficientA2
    END IF
  ELSE IF( GeometrySpecification(cComponent,3) ) THEN ! Cylinders
    CSW1CoefficientsMatrix = CSWCYLCoefficientA1
    IF( PYHCBCorrectionLogical ) THEN
      CSW2CoefficientsMatrix = CSWCYLCoefficientA2PY
    ELSE
      CSW2CoefficientsMatrix = CSWCYLCoefficientA2
    END IF
  END IF
  ! Overwrite the first-order coefficient matrix with a matrix for all geometries (Convex square-well potential)
  IF( UseA1AllGeometries ) CSW1CoefficientsMatrix = CSWAllCoefficientA1
  ! Overwrite the second-order coefficient matrix with a matrix for all geometries (Convex square-well potential)
  IF( UseA2AllGeometries ) CSW2CoefficientsMatrix = CSWAllCoefficientA2
  ! Overwrite the second-order coefficient matrix with the first-order coefficient matrix (Convex square-well potential)
  IF( UseA1ForA2Logical ) CSW2CoefficientsMatrix = CSW1CoefficientsMatrix
END IF

! Parametrization of the effective packing fraction [unitless]
IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
  EffPackingFractionCoefficients(1) = SWCoefficientsMatrix(1,1) + SWCoefficientsMatrix(1,2) * &
  &     ijPotentialRange(cComponent,cComponent) + SWCoefficientsMatrix(1,3) * ijPotentialRange(cComponent,cComponent) * &
  &     ijPotentialRange(cComponent,cComponent)
  EffPackingFractionCoefficients(2) = SWCoefficientsMatrix(2,1) + SWCoefficientsMatrix(2,2) * &
  &     ijPotentialRange(cComponent,cComponent) + SWCoefficientsMatrix(2,3) * ijPotentialRange(cComponent,cComponent) * &
  &     ijPotentialRange(cComponent,cComponent)
  EffPackingFractionCoefficients(3) = SWCoefficientsMatrix(3,1) + SWCoefficientsMatrix(3,2) * &
  &     ijPotentialRange(cComponent,cComponent) + SWCoefficientsMatrix(3,3) * ijPotentialRange(cComponent,cComponent) * &
  &     ijPotentialRange(cComponent,cComponent)
ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
  EffPackingFractionCoefficients(1) = SLCoefficientsMatrix(1,1) + SLCoefficientsMatrix(1,2) * &
  &     ijPotentialRange(cComponent,cComponent) + SLCoefficientsMatrix(1,3) * ijPotentialRange(cComponent,cComponent) * &
  &     ijPotentialRange(cComponent,cComponent) + SLCoefficientsMatrix(1,4) * ijPotentialRange(cComponent,cComponent) * &
  &     ijPotentialRange(cComponent,cComponent) * ijPotentialRange(cComponent,cComponent)
  EffPackingFractionCoefficients(2) = SLCoefficientsMatrix(2,1) + SLCoefficientsMatrix(2,2) * &
  &     ijPotentialRange(cComponent,cComponent) + SLCoefficientsMatrix(2,3) * ijPotentialRange(cComponent,cComponent) * &
  &     ijPotentialRange(cComponent,cComponent) + SLCoefficientsMatrix(2,4) * ijPotentialRange(cComponent,cComponent) * &
  &     ijPotentialRange(cComponent,cComponent) * ijPotentialRange(cComponent,cComponent)
  EffPackingFractionCoefficients(3) = 0.D0
  EffPackingFractionCoefficientsSutherland(1) = SLCoefficientsMatrix(1,1) + 2.D0 * SLCoefficientsMatrix(1,2) * &
  &     ijPotentialRange(cComponent,cComponent) + 4.D0 * SLCoefficientsMatrix(1,3) * ijPotentialRange(cComponent,cComponent) * &
  &     ijPotentialRange(cComponent,cComponent) + 8.D0 * SLCoefficientsMatrix(1,4) * ijPotentialRange(cComponent,cComponent) * &
  &     ijPotentialRange(cComponent,cComponent) * ijPotentialRange(cComponent,cComponent)
  EffPackingFractionCoefficientsSutherland(2) = SLCoefficientsMatrix(2,1) + 2.D0 * SLCoefficientsMatrix(2,2) * &
  &     ijPotentialRange(cComponent,cComponent) + 4.D0 * SLCoefficientsMatrix(2,3) * ijPotentialRange(cComponent,cComponent) * &
  &     ijPotentialRange(cComponent,cComponent) + 8.D0 * SLCoefficientsMatrix(2,4) * ijPotentialRange(cComponent,cComponent) * &
  &     ijPotentialRange(cComponent,cComponent) * ijPotentialRange(cComponent,cComponent)
  EffPackingFractionCoefficientsSutherland(3) = 0.D0
ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
  EffPackingFractionCoefficients(1) = YK1CoefficientsMatrix(1,1) + YK1CoefficientsMatrix(1,2) * ( 1.D0 / &
  &     ijPotentialRange(cComponent,cComponent) ) + YK1CoefficientsMatrix(1,3) * ( 1.D0 / &
  &     ( ijPotentialRange(cComponent,cComponent) * ijPotentialRange(cComponent,cComponent) ) )
  EffPackingFractionCoefficients(2) = YK1CoefficientsMatrix(2,1) + YK1CoefficientsMatrix(2,2) * ( 1.D0 / &
  &     ijPotentialRange(cComponent,cComponent) ) + YK1CoefficientsMatrix(2,3) * ( 1.D0 / &
  &     ( ijPotentialRange(cComponent,cComponent) * ijPotentialRange(cComponent,cComponent) ) )
  EffPackingFractionCoefficients(3) = 0.D0
  EffPackingFractionCoefficientsYukawa(1) = YK2CoefficientsMatrix(1,1) + YK2CoefficientsMatrix(1,2) * ( 1.D0 / &
  &     ijPotentialRange(cComponent,cComponent) ) + YK2CoefficientsMatrix(1,3) * ( 1.D0 / &
  &     ( ijPotentialRange(cComponent,cComponent) * ijPotentialRange(cComponent,cComponent) ) ) + YK2CoefficientsMatrix(1,4) * &
  &     ( 1.D0 / ( ijPotentialRange(cComponent,cComponent) * ijPotentialRange(cComponent,cComponent) * &
  &     ijPotentialRange(cComponent,cComponent) ) ) + YK2CoefficientsMatrix(1,5) * ( 1.D0 / &
  &     ( ijPotentialRange(cComponent,cComponent) * ijPotentialRange(cComponent,cComponent) * &
  &     ijPotentialRange(cComponent,cComponent) * ijPotentialRange(cComponent,cComponent) ) )
  EffPackingFractionCoefficientsYukawa(2) = YK2CoefficientsMatrix(2,1) + YK2CoefficientsMatrix(2,2) * ( 1.D0 / &
  &     ijPotentialRange(cComponent,cComponent) ) + YK2CoefficientsMatrix(2,3) * ( 1.D0 / &
  &     ( ijPotentialRange(cComponent,cComponent) * ijPotentialRange(cComponent,cComponent) ) ) + YK2CoefficientsMatrix(2,4) * &
  &     ( 1.D0 / ( ijPotentialRange(cComponent,cComponent) * ijPotentialRange(cComponent,cComponent) * &
  &     ijPotentialRange(cComponent,cComponent) ) ) + YK2CoefficientsMatrix(2,5) * ( 1.D0 / &
  &     ( ijPotentialRange(cComponent,cComponent) * ijPotentialRange(cComponent,cComponent) * &
  &     ijPotentialRange(cComponent,cComponent) * ijPotentialRange(cComponent,cComponent) ) )
  EffPackingFractionCoefficientsYukawa(3) = 0.D0
ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square well
  CSWAlphaCoefficients(1,1) = CSW1CoefficientsMatrix(1,1) + CSW1CoefficientsMatrix(2,1) * ( cNonSphericity - 1.D0 ) + &
  &     CSW1CoefficientsMatrix(3,1) * ( cNonSphericity - 1.D0 ) * ( cNonSphericity - 1.D0 )
  CSWAlphaCoefficients(1,2) = CSW1CoefficientsMatrix(1,2) + CSW1CoefficientsMatrix(2,2) * ( cNonSphericity - 1.D0 ) + &
  &     CSW1CoefficientsMatrix(3,2) * ( cNonSphericity - 1.D0 ) * ( cNonSphericity - 1.D0 )
  CSWAlphaCoefficients(1,3) = CSW1CoefficientsMatrix(1,3) + CSW1CoefficientsMatrix(2,3) * ( cNonSphericity - 1.D0 ) + &
  &     CSW1CoefficientsMatrix(3,3) * ( cNonSphericity - 1.D0 ) * ( cNonSphericity - 1.D0 )
  CSWAlphaCoefficients(1,4) = CSW1CoefficientsMatrix(1,4) + CSW1CoefficientsMatrix(2,4) * ( cNonSphericity - 1.D0 ) + &
  &     CSW1CoefficientsMatrix(3,4) * ( cNonSphericity - 1.D0 ) * ( cNonSphericity - 1.D0 )
  CSWAlphaCoefficients(1,5) = CSW1CoefficientsMatrix(1,5) + CSW1CoefficientsMatrix(2,5) * ( cNonSphericity - 1.D0 ) + &
  &     CSW1CoefficientsMatrix(3,5) * ( cNonSphericity - 1.D0 ) * ( cNonSphericity - 1.D0 )
  CSWAlphaCoefficients(1,6) = CSW1CoefficientsMatrix(1,6) + CSW1CoefficientsMatrix(2,6) * ( cNonSphericity - 1.D0 ) + &
  &     CSW1CoefficientsMatrix(3,6) * ( cNonSphericity - 1.D0 ) * ( cNonSphericity - 1.D0 )
  CSWAlphaCoefficients(1,7) = CSW1CoefficientsMatrix(1,7) + CSW1CoefficientsMatrix(2,7) * ( cNonSphericity - 1.D0 ) + &
  &     CSW1CoefficientsMatrix(3,7) * ( cNonSphericity - 1.D0 ) * ( cNonSphericity - 1.D0 )
  CSWAlphaCoefficients(1,8) = CSW1CoefficientsMatrix(1,8) + CSW1CoefficientsMatrix(2,8) * ( cNonSphericity - 1.D0 ) + &
  &     CSW1CoefficientsMatrix(3,8) * ( cNonSphericity - 1.D0 ) * ( cNonSphericity - 1.D0 )
  CSWAlphaCoefficients(1,9) = CSW1CoefficientsMatrix(1,9) + CSW1CoefficientsMatrix(2,9) * ( cNonSphericity - 1.D0 ) + &
  &     CSW1CoefficientsMatrix(3,9) * ( cNonSphericity - 1.D0 ) * ( cNonSphericity - 1.D0 )
  EffPackingFractionCoefficients(1) = CSWAlphaCoefficients(1,1) + CSWAlphaCoefficients(1,2) * ( 1.D0 + &
  &     ijPotentialRange(cComponent,cComponent) ) + CSWAlphaCoefficients(1,3) * ( 1.D0 + &
  &     ijPotentialRange(cComponent,cComponent) ) * ( 1.D0 + ijPotentialRange(cComponent,cComponent) )
  EffPackingFractionCoefficients(2) = CSWAlphaCoefficients(1,4) + CSWAlphaCoefficients(1,5) * ( 1.D0 + &
  &     ijPotentialRange(cComponent,cComponent) ) + CSWAlphaCoefficients(1,6) * ( 1.D0 + &
  &     ijPotentialRange(cComponent,cComponent) ) * ( 1.D0 + ijPotentialRange(cComponent,cComponent) )
  EffPackingFractionCoefficients(3) = CSWAlphaCoefficients(1,7) + CSWAlphaCoefficients(1,8) * ( 1.D0 + &
  &     ijPotentialRange(cComponent,cComponent) ) + CSWAlphaCoefficients(1,9) * ( 1.D0 + &
  &     ijPotentialRange(cComponent,cComponent) ) * ( 1.D0 + ijPotentialRange(cComponent,cComponent) )
  CSWAlphaCoefficients(2,1) = CSW2CoefficientsMatrix(1,1) + CSW2CoefficientsMatrix(2,1) * ( cNonSphericity - 1.D0 ) + &
  &     CSW2CoefficientsMatrix(3,1) * ( cNonSphericity - 1.D0 ) * ( cNonSphericity - 1.D0 )
  CSWAlphaCoefficients(2,2) = CSW2CoefficientsMatrix(1,2) + CSW2CoefficientsMatrix(2,2) * ( cNonSphericity - 1.D0 ) + &
  &     CSW2CoefficientsMatrix(3,2) * ( cNonSphericity - 1.D0 ) * ( cNonSphericity - 1.D0 )
  CSWAlphaCoefficients(2,3) = CSW2CoefficientsMatrix(1,3) + CSW2CoefficientsMatrix(2,3) * ( cNonSphericity - 1.D0 ) + &
  &     CSW2CoefficientsMatrix(3,3) * ( cNonSphericity - 1.D0 ) * ( cNonSphericity - 1.D0 )
  CSWAlphaCoefficients(2,4) = CSW2CoefficientsMatrix(1,4) + CSW2CoefficientsMatrix(2,4) * ( cNonSphericity - 1.D0 ) + &
  &     CSW2CoefficientsMatrix(3,4) * ( cNonSphericity - 1.D0 ) * ( cNonSphericity - 1.D0 )
  CSWAlphaCoefficients(2,5) = CSW2CoefficientsMatrix(1,5) + CSW2CoefficientsMatrix(2,5) * ( cNonSphericity - 1.D0 ) + &
  &     CSW2CoefficientsMatrix(3,5) * ( cNonSphericity - 1.D0 ) * ( cNonSphericity - 1.D0 )
  CSWAlphaCoefficients(2,6) = CSW2CoefficientsMatrix(1,6) + CSW2CoefficientsMatrix(2,6) * ( cNonSphericity - 1.D0 ) + &
  &     CSW2CoefficientsMatrix(3,6) * ( cNonSphericity - 1.D0 ) * ( cNonSphericity - 1.D0 )
  CSWAlphaCoefficients(2,7) = CSW2CoefficientsMatrix(1,7) + CSW2CoefficientsMatrix(2,7) * ( cNonSphericity - 1.D0 ) + &
  &     CSW2CoefficientsMatrix(3,7) * ( cNonSphericity - 1.D0 ) * ( cNonSphericity - 1.D0 )
  CSWAlphaCoefficients(2,8) = CSW2CoefficientsMatrix(1,8) + CSW2CoefficientsMatrix(2,8) * ( cNonSphericity - 1.D0 ) + &
  &     CSW2CoefficientsMatrix(3,8) * ( cNonSphericity - 1.D0 ) * ( cNonSphericity - 1.D0 )
  CSWAlphaCoefficients(2,9) = CSW2CoefficientsMatrix(1,9) + CSW2CoefficientsMatrix(2,9) * ( cNonSphericity - 1.D0 ) + &
  &     CSW2CoefficientsMatrix(3,9) * ( cNonSphericity - 1.D0 ) * ( cNonSphericity - 1.D0 )
  EffPackingFractionCoefficientsCSW(1) = CSWAlphaCoefficients(2,1) + CSWAlphaCoefficients(2,2) * ( 1.D0 + &
  &     ijPotentialRange(cComponent,cComponent) ) + CSWAlphaCoefficients(2,3) * ( 1.D0 + &
  &     ijPotentialRange(cComponent,cComponent) ) * ( 1.D0 + ijPotentialRange(cComponent,cComponent) )
  EffPackingFractionCoefficientsCSW(2) = CSWAlphaCoefficients(2,4) + CSWAlphaCoefficients(2,5) * ( 1.D0 + &
  &     ijPotentialRange(cComponent,cComponent) ) + CSWAlphaCoefficients(2,6) * ( 1.D0 + &
  &     ijPotentialRange(cComponent,cComponent) ) * ( 1.D0 + ijPotentialRange(cComponent,cComponent) )
  EffPackingFractionCoefficientsCSW(3) = CSWAlphaCoefficients(2,7) + CSWAlphaCoefficients(2,8) * ( 1.D0 + &
  &     ijPotentialRange(cComponent,cComponent) ) + CSWAlphaCoefficients(2,9) * ( 1.D0 + &
  &     ijPotentialRange(cComponent,cComponent) ) * ( 1.D0 + ijPotentialRange(cComponent,cComponent) )
END IF

! Effective packing fraction [unitless]
IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
  EffPackingFraction = EffPackingFractionCoefficients(1) * rDensity(3) + EffPackingFractionCoefficients(2) * rDensity(3) * &
  &     rDensity(3) + EffPackingFractionCoefficients(3) * rDensity(3) * rDensity(3) * rDensity(3)
ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
  EffPackingFraction = EffPackingFractionCoefficients(1) * rDensity(3) + EffPackingFractionCoefficients(2) * rDensity(3) * &
  &     rDensity(3)
  EffPackingFractionSutherland = EffPackingFractionCoefficientsSutherland(1) * rDensity(3) + rDensity(3) * rDensity(3) * &
  &     EffPackingFractionCoefficientsSutherland(2)
ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
  EffPackingFraction = EffPackingFractionCoefficients(1) * rDensity(3) + EffPackingFractionCoefficients(2) * rDensity(3) * &
  &     rDensity(3)
  EffPackingFractionYukawa = EffPackingFractionCoefficientsYukawa(1) * rDensity(3) + EffPackingFractionCoefficientsYukawa(2) * &
  &     rDensity(3) * rDensity(3)
ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square well
  EffPackingFraction = EffPackingFractionCoefficients(1) * rDensity(3) + EffPackingFractionCoefficients(2) * rDensity(3) * &
  &     rDensity(3) + EffPackingFractionCoefficients(3) * rDensity(3) * rDensity(3) * rDensity(3)
  EffPackingFractionCSW = EffPackingFractionCoefficientsCSW(1) * rDensity(3) + EffPackingFractionCoefficientsCSW(2) * &
  &     rDensity(3) * rDensity(3) + EffPackingFractionCoefficientsCSW(3) * rDensity(3) * rDensity(3) * rDensity(3)
END IF

! Contact radial distribution function [unitless]
cRadialDistributionFunction = ( 1.D0 / ( 1.D0 - rDensity(3) ) ) + ( ( 3.D0 * rDensity(2) * cAverageDiameterRelationship ) / &
&     ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) ) + ( ( 2.D0 * cAverageDiameterRelationship * &
&     cAverageDiameterRelationship * rDensity(2) * rDensity(2) ) / ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * &
&     ( 1.D0 - rDensity(3) ) ) )

! Contact radial distribution function for an effective packing fraction [unitless]
IF( .NOT. PotentialTypeLogical(4) ) THEN
  cEffectiveRadialDistributionFunction = ( 1.D0 / ( 1.D0 - EffPackingFraction ) ) + ( ( 3.D0 * cAverageDiameterRelationship * &
  &     EffPackingFraction * AuxDiameterRelationship ) / ( ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) ) ) + &
  &     ( ( 2.D0 * cAverageDiameterRelationship * cAverageDiameterRelationship * EffPackingFraction * AuxDiameterRelationship * &
  &     EffPackingFraction * AuxDiameterRelationship ) / ( ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) * ( 1.D0 &
  &     - EffPackingFraction ) ) )
ELSE
  cEffectiveRadialDistributionFunction = 1.D0 / ( 1.D0 - EffPackingFraction ) + 3.D0 * EffPackingFraction * cNonSphericity * &
  &     ( 1.D0 + cNonSphericity ) / ( ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) * ( 1.D0 + 3.D0 * &
  &     cNonSphericity ) ) + 2.D0 * EffPackingFraction * EffPackingFraction * cNonSphericity * cNonSphericity / ( ( 1.D0 - &
  &     EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) * ( 1.D0 + 3.D0 * cNonSphericity ) )
END IF
IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
  cEffectiveRadialDistributionFunctionSL = ( 1.D0 / ( 1.D0 - EffPackingFractionSutherland ) ) + ( ( 3.D0 * &
  &     cAverageDiameterRelationship * EffPackingFractionSutherland * AuxDiameterRelationship ) / ( ( 1.D0 - &
  &     EffPackingFractionSutherland ) * ( 1.D0 - EffPackingFractionSutherland ) ) ) + ( ( 2.D0 * cAverageDiameterRelationship * &
  &     cAverageDiameterRelationship * EffPackingFractionSutherland * AuxDiameterRelationship * EffPackingFractionSutherland * &
  &     AuxDiameterRelationship ) / ( ( 1.D0 - EffPackingFractionSutherland ) * ( 1.D0 - EffPackingFractionSutherland ) * ( 1.D0 - &
  &     EffPackingFractionSutherland ) ) )
ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
  cEffectiveRadialDistributionFunctionYK = ( 1.D0 / ( 1.D0 - EffPackingFractionYukawa ) ) + ( ( 3.D0 * &
  &     cAverageDiameterRelationship * EffPackingFractionYukawa * AuxDiameterRelationship ) / ( ( 1.D0 - EffPackingFractionYukawa &
  &     ) * ( 1.D0 - EffPackingFractionYukawa ) ) ) + ( ( 2.D0 * cAverageDiameterRelationship * cAverageDiameterRelationship * &
  &     EffPackingFractionYukawa * AuxDiameterRelationship * EffPackingFractionYukawa * AuxDiameterRelationship ) / ( ( 1.D0 - &
  &     EffPackingFractionYukawa ) * ( 1.D0 - EffPackingFractionYukawa ) * ( 1.D0 - EffPackingFractionYukawa ) ) )
ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square well
  cEffectiveRadialDistributionFunctionCSW = 1.D0 / ( 1.D0 - EffPackingFractionCSW ) + 3.D0 * EffPackingFractionCSW * &
  &     cNonSphericity * ( 1.D0 + cNonSphericity ) / ( ( 1.D0 - EffPackingFractionCSW ) * ( 1.D0 - EffPackingFractionCSW ) * &
  &     ( 1.D0 + 3.D0 * cNonSphericity ) ) + 2.D0 * EffPackingFractionCSW * EffPackingFractionCSW * cNonSphericity * &
  &     cNonSphericity / ( ( 1.D0 - EffPackingFractionCSW ) * ( 1.D0 - EffPackingFractionCSW ) * ( 1.D0 - EffPackingFractionCSW ) &
  &     * ( 1.D0 + 3.D0 * cNonSphericity ) )
END IF

! Mean-attractive energy between components [K]
IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
  cMeanAttractiveEnergy = - 2.D0 * cPi * rNumberDensity * ijaWellDepth(cComponent,cComponent) * &
  &     ( ijaDiameterSphereCubic(cComponent,cComponent) / 3.D0 ) * ( ijPotentialRangeCubic(cComponent,cComponent) - 1.D0 ) * &
  &     cEffectiveRadialDistributionFunction
ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
  cMeanAttractiveEnergy = - 2.D0 * cPi * rNumberDensity * ijaWellDepth(cComponent,cComponent) * &
  &     ijaDiameterSphereCubic(cComponent,cComponent) * ( 1.D0 / ( ijPotentialRange(cComponent,cComponent) - 3.D0 ) ) * &
  &     cEffectiveRadialDistributionFunction
  cMeanAttractiveEnergySutherland = - 2.D0 * cPi * rNumberDensity * ijaWellDepth(cComponent,cComponent) &
  &     * ijaDiameterSphereCubic(cComponent,cComponent) * ( 1.D0 / ( 2.D0 * ijPotentialRange(cComponent,cComponent) - 3.D0 ) ) &
  &     * cEffectiveRadialDistributionFunctionSL
ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
  cMeanAttractiveEnergy = - 2.D0 * cPi * rNumberDensity * ijaWellDepth(cComponent,cComponent) * &
  &     ijaDiameterSphereCubic(cComponent,cComponent) * ( ( 1.D0 / ijPotentialRange(cComponent,cComponent) ) + ( 1.D0 / &
  &     ( ijPotentialRange(cComponent,cComponent) * ijPotentialRange(cComponent,cComponent) ) ) ) * &
  &     cEffectiveRadialDistributionFunction
  cMeanAttractiveEnergyYukawa = - cPi * rNumberDensity * ijaWellDepth(cComponent,cComponent) * &
  &     ijaDiameterSphereCubic(cComponent,cComponent) * ( 1.D0 / ijPotentialRange(cComponent,cComponent) ) * &
  &     cEffectiveRadialDistributionFunctionYK
ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square well
  cMeanAttractiveEnergy = - rNumberDensity * cEffectiveRadialDistributionFunction * ijaWellDepth(cComponent,cComponent) * &
  &     ( fMolecularVolume * ( 1.D0 + 3.D0 * fNonSphericity ) - cMolecularVolume * ( 1.D0 + 3.D0 * cNonSphericity ) )
  cMeanAttractiveEnergyCSW = - rNumberDensity * cEffectiveRadialDistributionFunctionCSW * ijaWellDepth(cComponent,cComponent) * &
  &     ( fMolecularVolume * ( 1.D0 + 3.D0 * fNonSphericity ) - cMolecularVolume * ( 1.D0 + 3.D0 * cNonSphericity ) )
END IF

! First-order perturbation contribution to the Helmholtz free energy (molar basis) [J / mol]
mMeanAttractiveHelmholtzFreeEnergy = cMeanAttractiveEnergy * cUniversalGas ! Proven units

! ************************************************************************************************ !
! SECOND-ORDER PERTURBATION CONTRIBUTION                                                           !
! ************************************************************************************************ !

! Hard-sphere isothermal compressibility given by the Percus-Yevick expression [unitless]
HSIsothermalCompressibility = ( rDensity(0) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - &
&     rDensity(3) ) ) / ( ( rDensity(0) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) + ( 6.D0 * rDensity(1) * rDensity(2) * &
&     ( 1.D0 - rDensity(3) ) ) + ( 9.D0 * rDensity(2) * rDensity(2) * rDensity(2) ) )

! Hard convex-body isothermal compressibility given by the Boublík expression [unitless]
IF( PYHCBCorrectionLogical ) THEN
  HCBIsothermalCompressibility = ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) &
  &     / ( 1.D0 + 2.D0 * rDensity(3) * ( 3.D0 * cNonSphericity - 1.D0 ) + rDensity(3) * rDensity(3) * ( 9.D0 * cNonSphericity * &
  &     cNonSphericity - 6.D0 * cNonSphericity + 1.D0 ) )
ELSE
  HCBIsothermalCompressibility = ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) &
  &     / ( 1.D0 + 2.D0 * rDensity(3) * ( 3.D0 * cNonSphericity - 1.D0 ) + rDensity(3) * rDensity(3) * ( 9.D0 * cNonSphericity * &
  &     cNonSphericity - 6.D0 * cNonSphericity + 1.D0 ) - 4.D0 * rDensity(3) * rDensity(3) * rDensity(3) * cNonSphericity * &
  &     cNonSphericity + rDensity(3) * rDensity(3) * rDensity(3) * rDensity(3) * cNonSphericity * cNonSphericity )
END IF

! First derivative of the effective packing fraction with respect to the density [Å³]
IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
  dEffPackingFraction_dDensity = EffPackingFractionCoefficients(1) * rDensity(3) + 2.D0 * EffPackingFractionCoefficients(2) * &
  &     rDensity(3) * rDensity(3) + 3.D0 * rDensity(3) * rDensity(3) * rDensity(3) * EffPackingFractionCoefficients(3)
ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
  dEffPackingFraction_dDensity = EffPackingFractionCoefficients(1) * rDensity(3) + 2.D0 * EffPackingFractionCoefficients(2) * &
  &     rDensity(3) * rDensity(3)
  dEffPackingFractionSutherland_dDensity = EffPackingFractionCoefficientsSutherland(1) * rDensity(3) + 2.D0 * &
  &     EffPackingFractionCoefficientsSutherland(2) * rDensity(3) * rDensity(3)
ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
  dEffPackingFraction_dDensity = EffPackingFractionCoefficients(1) * rDensity(3) + 2.D0 * EffPackingFractionCoefficients(2) * &
  &     rDensity(3) * rDensity(3)
  dEffPackingFractionYukawa_dDensity = EffPackingFractionCoefficientsYukawa(1) * rDensity(3) + 2.D0 * &
  &     EffPackingFractionCoefficientsYukawa(2) * rDensity(3) * rDensity(3)
ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square well
  dEffPackingFraction_dDensity = EffPackingFractionCoefficients(1) * rDensity(3) + 2.D0 * EffPackingFractionCoefficients(2) * &
  &     rDensity(3) * rDensity(3) + 3.D0 * rDensity(3) * rDensity(3) * rDensity(3) * EffPackingFractionCoefficients(3)
  dEffPackingFractionCSW_dDensity = EffPackingFractionCoefficientsCSW(1) * rDensity(3) + 2.D0 * rDensity(3) * rDensity(3) * &
  &     EffPackingFractionCoefficientsCSW(2) + 3.D0 * EffPackingFractionCoefficientsCSW(3) * rDensity(3) * rDensity(3) * &
  &     rDensity(3)
END IF
dEffPackingFraction_dDensity = dEffPackingFraction_dDensity / rNumberDensity
IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
  dEffPackingFractionSutherland_dDensity = dEffPackingFractionSutherland_dDensity / rNumberDensity
ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
  dEffPackingFractionYukawa_dDensity = dEffPackingFractionYukawa_dDensity / rNumberDensity
ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square well
  dEffPackingFractionCSW_dDensity = dEffPackingFractionCSW_dDensity / rNumberDensity
END IF

! First derivative of the radial distribution function with respect to the density [Å³]
IF( .NOT. PotentialTypeLogical(4) ) THEN
  dEffRadialDistributionFunct_dDensity = ( dEffPackingFraction_dDensity / ( ( 1.D0 - EffPackingFraction ) * ( 1.D0 - &
  &     EffPackingFraction ) ) ) + ( ( 3.D0 * cAverageDiameterRelationship * ( dEffPackingFraction_dDensity * ( 1.D0 + &
  &     EffPackingFraction ) ) * AuxDiameterRelationship ) / ( ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) * &
  &     ( 1.D0 - EffPackingFraction ) ) ) + ( ( 2.D0 * cAverageDiameterRelationship * cAverageDiameterRelationship * &
  &     AuxDiameterRelationship * AuxDiameterRelationship * ( ( 2.D0 + EffPackingFraction ) * EffPackingFraction * &
  &     dEffPackingFraction_dDensity ) ) / ( ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) * ( 1.D0 - &
  &     EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) ) )
ELSE
  gAux1 = cNonSphericity * ( 1.D0 + cNonSphericity ) / ( 1.D0 + 3.D0 * cNonSphericity )
  gAux2 = cNonSphericity * cNonSphericity / ( 1.D0 + 3.D0 * cNonSphericity )
  dEffRadialDistributionFunct_dDensity = ( dEffPackingFraction_dDensity / ( ( 1.D0 - EffPackingFraction ) * ( 1.D0 - &
  &     EffPackingFraction ) ) ) + ( 3.D0 * gAux1 * ( dEffPackingFraction_dDensity * ( 1.D0 + EffPackingFraction ) ) / ( ( 1.D0 - &
  &     EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) ) ) + ( ( 2.D0 * gAux2 * &
  &     EffPackingFraction * ( ( 2.D0 + EffPackingFraction ) * dEffPackingFraction_dDensity ) ) / ( ( 1.D0 - EffPackingFraction ) &
  &     * ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) * ( 1.D0 - EffPackingFraction ) ) )
END IF
IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
  dEffRadialDistributionFunctSL_dDensity = ( dEffPackingFractionSutherland_dDensity / ( ( 1.D0 - EffPackingFractionSutherland ) * &
  &     ( 1.D0 - EffPackingFractionSutherland ) ) ) + ( ( 3.D0 * cAverageDiameterRelationship * ( &
  &     dEffPackingFractionSutherland_dDensity * ( 1.D0 + EffPackingFractionSutherland ) ) * AuxDiameterRelationship ) / ( ( 1.D0 &
  &     - EffPackingFractionSutherland ) * ( 1.D0 - EffPackingFractionSutherland ) * ( 1.D0 - EffPackingFractionSutherland ) ) ) + &
  &     ( ( 2.D0 * cAverageDiameterRelationship * cAverageDiameterRelationship * AuxDiameterRelationship * AuxDiameterRelationship &
  &     * ( ( 2.D0 + EffPackingFractionSutherland ) * EffPackingFractionSutherland * dEffPackingFractionSutherland_dDensity ) ) / &
  &     ( ( 1.D0 - EffPackingFractionSutherland ) * ( 1.D0 - EffPackingFractionSutherland ) * ( 1.D0 - &
  &     EffPackingFractionSutherland ) * ( 1.D0 - EffPackingFractionSutherland ) ) )
ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
  dEffRadialDistributionFunctYK_dDensity = ( dEffPackingFractionYukawa_dDensity / ( ( 1.D0 - EffPackingFractionYukawa ) * &
  &     ( 1.D0 - EffPackingFractionYukawa ) ) ) + ( ( 3.D0 * cAverageDiameterRelationship * ( dEffPackingFractionYukawa_dDensity * &
  &     ( 1.D0 + EffPackingFractionYukawa ) ) * AuxDiameterRelationship ) / ( ( 1.D0 - EffPackingFractionYukawa ) * ( 1.D0 - &
  &     EffPackingFractionYukawa ) * ( 1.D0 - EffPackingFractionYukawa ) ) ) + ( ( 2.D0 * cAverageDiameterRelationship * &
  &     cAverageDiameterRelationship * AuxDiameterRelationship * AuxDiameterRelationship * ( ( 2.D0 + EffPackingFractionYukawa ) * &
  &     EffPackingFractionYukawa * dEffPackingFractionYukawa_dDensity ) ) / ( ( 1.D0 - EffPackingFractionYukawa ) * ( 1.D0 - &
  &     EffPackingFractionYukawa ) * ( 1.D0 - EffPackingFractionYukawa ) * ( 1.D0 - EffPackingFractionYukawa ) ) )
ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square well
  dEffRadialDistributionFunctCSW_dDensity = ( dEffPackingFractionCSW_dDensity / ( ( 1.D0 - EffPackingFractionCSW ) * ( 1.D0 - &
  &     EffPackingFractionCSW ) ) ) + ( 3.D0 * gAux1 * ( dEffPackingFractionCSW_dDensity * ( 1.D0 + EffPackingFractionCSW ) ) / ( &
  &     ( 1.D0 - EffPackingFractionCSW ) * ( 1.D0 - EffPackingFractionCSW ) * ( 1.D0 - EffPackingFractionCSW ) ) ) + ( ( 2.D0 * &
  &     gAux2 * EffPackingFractionCSW * ( ( 2.D0 + EffPackingFractionCSW ) * dEffPackingFractionCSW_dDensity ) ) / ( ( 1.D0 - &
  &     EffPackingFractionCSW ) * ( 1.D0 - EffPackingFractionCSW ) * ( 1.D0 - EffPackingFractionCSW ) * ( 1.D0 - &
  &     EffPackingFractionCSW ) ) )
END IF

! First derivative of the mean-attractive energy with respect to the density [K . Å³]
IF( PotentialTypeLogical(1) ) THEN ! Square-well potential
  dMeanAttractiveEnergy_dDensity = cMeanAttractiveEnergy / rNumberDensity + cMeanAttractiveEnergy / &
  &     cEffectiveRadialDistributionFunction * dEffRadialDistributionFunct_dDensity
ELSE IF( PotentialTypeLogical(2) ) THEN ! Sutherland potential
  dMeanAttractiveEnergy_dDensity = cMeanAttractiveEnergySutherland / rNumberDensity + cMeanAttractiveEnergySutherland / &
  &     cEffectiveRadialDistributionFunctionSL * dEffRadialDistributionFunctSL_dDensity
ELSE IF( PotentialTypeLogical(3) ) THEN ! Yukawa potential
  dMeanAttractiveEnergy_dDensity = cMeanAttractiveEnergyYukawa / rNumberDensity + cMeanAttractiveEnergyYukawa / &
  &     cEffectiveRadialDistributionFunctionYK * dEffRadialDistributionFunctYK_dDensity
ELSE IF( PotentialTypeLogical(4) ) THEN ! Convex square well
  dMeanAttractiveEnergy_dDensity = cMeanAttractiveEnergyCSW / rNumberDensity + cMeanAttractiveEnergyCSW / &
  &     cEffectiveRadialDistributionFunctionCSW * dEffRadialDistributionFunctCSW_dDensity
END IF

! Mean-attractive energy fluctuation between components [K²]
IF( .NOT. PotentialTypeLogical(4) ) THEN
  cMeanAttractiveEnergyFluctuations = 0.5D0 * dMeanAttractiveEnergy_dDensity * ijaWellDepth(cComponent,cComponent) * &
  &     rNumberDensity * HSIsothermalCompressibility
ELSE
  cMeanAttractiveEnergyFluctuations = 0.5D0 * dMeanAttractiveEnergy_dDensity * ijaWellDepth(cComponent,cComponent) * &
  &     rNumberDensity * HCBIsothermalCompressibility
END IF

! Zhang's correction [K²]
IF( ZhangCorrectionLogical ) THEN
  cMeanAttractiveEnergyFluctuations = cMeanAttractiveEnergyFluctuations * ZhangFactor
END IF

! Second-order perturbation contribution to the Helmholtz free energy (molar basis) [J / mol]
mMeanAttFluctuationHelmholtzFreeEnergy = cMeanAttractiveEnergyFluctuations * cUniversalGas / Temperature ! Proven units

! ************************************************************************************************ !
! HIGHER-ORDER PERTURBATION CONTRIBUTION                                                           !
! ************************************************************************************************ !

! Higher order perturbation contribution to the Helmholtz free energy (molar basis) [J / mol]
IF( HigherOrderTPTLogical ) THEN
  ! Barker-Henderson's higher-order perturbation theory (approximation)
  HigherOrderTerms = 0.D0
  DO nOrder = 3, nHigherOrder
    HigherOrderTerms = HigherOrderTerms + ( mMeanAttractiveHelmholtzFreeEnergy * ( ( (2.D0 * &
    &     mMeanAttFluctuationHelmholtzFreeEnergy) / mMeanAttractiveHelmholtzFreeEnergy ) ** ( nOrder - 1 ) ) / &
    &     Factorial( nOrder ) )
  END DO
  ! Higher order perturbation contribution to the Helmholtz free energy (molar basis) [J / mol]
  mHigherOrderHelmholtzFreeEnergy = HigherOrderTerms ! Proven units
ELSE
  ! Higher order perturbation contribution to the Helmholtz free energy (molar basis) [J / mol]
  mHigherOrderHelmholtzFreeEnergy = 0.D0 ! Proven units
END IF

! ************************************************************************************************ !
! MONOMER-MONOMER TOTAL CONTRIBUTION                                                               !
! ************************************************************************************************ !

! Monomer-monomer contribution to the Helmholtz free energy (molar basis) [J / mol]
mMonomerMonomerHelmholtzFreeEnergy = mExcludedVolumeHelmholtzFreeEnergy + mMeanAttractiveHelmholtzFreeEnergy + &
&     mMeanAttFluctuationHelmholtzFreeEnergy + mHigherOrderHelmholtzFreeEnergy ! Proven units

! ************************************************************************************************ !
! TOTAL HELMHOLTZ FREE ENERGY                                                                      !
! ************************************************************************************************ !

! Total Helmholtz free energy (molar basis) [J / mol]
mTotalHelmholtzFreeEnergy = mIdealHelmholtzFreeEnergy + mMonomerMonomerHelmholtzFreeEnergy ! Proven units

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !
!          FIRST DERIVATIVES OF THE HELMHOLTZ FREE ENERGY WITH RESPECT TO THE TEMPERATURE          !
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !

! ************************************************************************************************ !
! IDEAL CONTRIBUTION                                                                               !
! ************************************************************************************************ !

! First derivative of the ideal Helmholtz free energy with respect to the temperature [J / (mol . K)]
dIdealHelmholtzFreeEnergy_dTemperature = DLOG( cDeBroglieIdeal * DeBroglie_Wavelength_Translational_Cb * &
&     DeBroglie_Wavelength_Rotational / mVolume )
dIdealHelmholtzFreeEnergy_dTemperature = dIdealHelmholtzFreeEnergy_dTemperature - 3.5D0
dIdealHelmholtzFreeEnergy_dTemperature = cUniversalGas * dIdealHelmholtzFreeEnergy_dTemperature ! Proven units

! ************************************************************************************************ !
! EXCLUDED-VOLUME CONTRIBUTION                                                                     !
! ************************************************************************************************ !

! First derivative of the excluded-volume contribution to the Helmholtz free energy with respect to the temperature [J / (mol . K)]
dExcludedVolumeFreeEnergy_dTemperature = HSBoublikMansoori * cUniversalGas * AngleAverage ! Proven units
IF( ReferenceBoublikLogical ) THEN
  dExcludedVolumeFreeEnergy_dTemperature = ( cNonSphericity * rDensity(3) * ( cNonSphericity + 3.D0 - 3.D0 * rDensity(3) ) / &
  &     ( ( 1.D0 - rDensity(3) ) * ( 1.D0 - rDensity(3) ) ) ) - ( ( 1.D0 - cNonSphericity * cNonSphericity ) * &
  &     DLOG( 1.D0 - rDensity(3) ) )
  dExcludedVolumeFreeEnergy_dTemperature = dExcludedVolumeFreeEnergy_dTemperature * cUniversalGas
END IF

! ************************************************************************************************ !
! SECOND-ORDER PERTURBATION CONTRIBUTION                                                           !
! ************************************************************************************************ !

! First derivative of the second-order perturbation contribution to the Helmholtz free energy with respect to the temperature [J / (mol . K)]
dMeanAttFluctuationFEnergy_dTemperature = - cMeanAttractiveEnergyFluctuations * cUniversalGas / ( Temperature * Temperature ) ! Proven units

! ************************************************************************************************ !
! HIGHER-ORDER PERTURBATION CONTRIBUTION                                                           !
! ************************************************************************************************ !

! First derivative of the higher order perturbation contribution to the Helmholtz free energy with respect to the temperature [J / (mol . K)]
IF( HigherOrderTPTLogical ) THEN
  ! Barker-Henderson's higher-order perturbation theory (approximation)
  dHigherOrderTerms_dTemperature = 0.D0
  DO nOrder = 3, nHigherOrder
    dHigherOrderTerms_dTemperature = dHigherOrderTerms_dTemperature + ( 2.D0 * ( DBLE( nOrder - 1 ) / Factorial( nOrder ) ) ) * &
    &     ( ( (2.D0 * mMeanAttFluctuationHelmholtzFreeEnergy) / mMeanAttractiveHelmholtzFreeEnergy ) ** ( nOrder - 2 ) ) * &
    &     ( dMeanAttFluctuationFEnergy_dTemperature )
  END DO
  ! First derivative of the higher order perturbation contribution to the Helmholtz free energy with respect to the temperature [J / (mol . K)]
  dHigherOrderFEnergy_dTemperature = dHigherOrderTerms_dTemperature ! Proven units
ELSE
  ! First derivative of the higher order perturbation contribution to the Helmholtz free energy with respect to the temperature [J / (mol . K)]
  dHigherOrderFEnergy_dTemperature = 0.D0 ! Proven units
END IF

! ************************************************************************************************ !
! MONOMER-MONOMER TOTAL CONTRIBUTION                                                               !
! ************************************************************************************************ !

! First derivative of the monomer-monomer contribution to the Helmholtz free energy with respect to the temperature [J / (mol . K)]
dMonomerMonomerFEnergy_dTemperature = dExcludedVolumeFreeEnergy_dTemperature + dMeanAttFluctuationFEnergy_dTemperature + &
&     dHigherOrderFEnergy_dTemperature ! Proven units

! ************************************************************************************************ !
! TOTAL HELMHOLTZ FREE ENERGY                                                                      !
! ************************************************************************************************ !

! First derivative of the total Helmholtz free energy with respect to the temperature [J / (mol . K)]
dTotalHelmholtzFEnergy_dTemperature = dIdealHelmholtzFreeEnergy_dTemperature + dMonomerMonomerFEnergy_dTemperature ! Proven units

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !
!         SECOND DERIVATIVES OF THE HELMHOLTZ FREE ENERGY WITH RESPECT TO THE TEMPERATURE          !
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !

! ************************************************************************************************ !
! IDEAL CONTRIBUTION                                                                               !
! ************************************************************************************************ !

! Second derivative of the ideal Helmholtz free energy with respect to the temperature [J / (mol . K²)]
d2IdealHelmholtzFEnergy_d2Temperature = - 2.5D0 * cUniversalGas / Temperature ! Proven units

! ************************************************************************************************ !
! SECOND-ORDER PERTURBATION CONTRIBUTION                                                           !
! ************************************************************************************************ !

! Second derivative of the second-order perturbation contribution to the Helmholtz free energy with respect to the temperature [J / (mol . K²)]
d2MeanAttFluctFEnergy_d2Temperature = 2.D0 * cMeanAttractiveEnergyFluctuations * cUniversalGas / ( Temperature * Temperature * &
&     Temperature ) ! Proven units

! ************************************************************************************************ !
! HIGHER-ORDER PERTURBATION CONTRIBUTION                                                           !
! ************************************************************************************************ !

! Second derivative of the higher order perturbation contribution to the Helmholtz free energy with respect to the temperature [J / (mol . K²)]
IF( HigherOrderTPTLogical ) THEN
  ! Barker-Henderson's higher-order perturbation theory (approximation)
  d2HigherOrderTerms_d2Temperature = 0.D0
  DO nOrder = 3, nHigherOrder
    d2HigherOrderTerms_d2Temperature = d2HigherOrderTerms_d2Temperature + ( ( 4.D0 * ( DBLE( nOrder - 1 ) / &
    &     Factorial( nOrder ) ) ) * ( ( (2.D0 * mMeanAttFluctuationHelmholtzFreeEnergy) / mMeanAttractiveHelmholtzFreeEnergy ) ** &
    &     ( nOrder - 3 ) ) / mMeanAttractiveHelmholtzFreeEnergy ) * ( ( d2MeanAttFluctFEnergy_d2Temperature * &
    &     mMeanAttFluctuationHelmholtzFreeEnergy ) + ( DBLE( nOrder - 2 ) * dMeanAttFluctuationFEnergy_dTemperature * &
    &     dMeanAttFluctuationFEnergy_dTemperature ) )
  END DO
  ! Second derivative of the higher order perturbation contribution to the Helmholtz free energy with respect to the temperature [J / (mol . K²)]
  d2HigherOrderFEnergy_d2Temperature = d2HigherOrderTerms_d2Temperature ! Proven units
ELSE
  ! Second derivative of the higher order perturbation contribution to the Helmholtz free energy with respect to the temperature [J / (mol . K²)]
  d2HigherOrderFEnergy_d2Temperature = 0.D0 ! Proven units
END IF

! ************************************************************************************************ !
! MONOMER-MONOMER TOTAL CONTRIBUTION                                                               !
! ************************************************************************************************ !

! Second derivative of the monomer-monomer contribution to the Helmholtz free energy with respect to the temperature [J / (mol . K²)]
d2MonomerMonomerFEnergy_d2Temperature = d2MeanAttFluctFEnergy_d2Temperature + d2HigherOrderFEnergy_d2Temperature ! Proven units

! ************************************************************************************************ !
! TOTAL HELMHOLTZ FREE ENERGY                                                                      !
! ************************************************************************************************ !

! Second derivative of the total Helmholtz free energy with respect to the temperature [J / (mol . K²)]
d2TotalHelmholtzFEnergy_d2Temperature = d2IdealHelmholtzFEnergy_d2Temperature + d2MonomerMonomerFEnergy_d2Temperature ! Proven units

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !
!                                     THERMODYNAMIC PROPERTIES                                     !
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !

! Specific heat at constant volume (ideal contribution) [J / (mol . K)]
IF( SpecificHeatReference == "TRC" ) THEN
  vIdealSpecificHeat = cUniversalGas * ( TemperatureParameter(cComponent,1) + TemperatureParameter(cComponent,2) * Temperature &
  &     + TemperatureParameter(cComponent,3) * Temperature * Temperature + TemperatureParameter(cComponent,4) * Temperature * &
  &     Temperature * Temperature + TemperatureParameter(cComponent,5) * Temperature * Temperature * Temperature * Temperature - &
  &     1.D0 )
ELSE IF( SpecificHeatReference == "NIST" ) THEN
  vIdealSpecificHeat = ( TemperatureParameter(cComponent,1) + TemperatureParameter(cComponent,2) * ( Temperature / 1.D3 ) + &
  &     TemperatureParameter(cComponent,3) * ( Temperature * Temperature / 1.D6 ) + TemperatureParameter(cComponent,4) * &
  &     ( Temperature * Temperature * Temperature / 1.D9 ) + TemperatureParameter(cComponent,5) / ( Temperature * Temperature / &
  &     1.D6 ) - cUniversalGas )
END IF

! Specific heat at constant volume (ideal contribution + residual contribution) [J / (mol . K)]
vSpecificHeat = vIdealSpecificHeat - ( Temperature * d2MonomerMonomerFEnergy_d2Temperature )

! Entropy [J / (mol . K)]
Entropy = - dTotalHelmholtzFEnergy_dTemperature

RETURN

END SUBROUTINE Calculate_Entropy_Single_Component